"""E16 -- four further external genes: records, coordinates and labels; E16b -- a
fifth, SBDS, and two checks that change nothing.

CTCF, SFPQ, RAD51D and XRCC2 have HAP1 saturation genome editing deposits in MaveDB
(CC0), scored for cell fitness like the seven genes and DDX3X. They are added as
external genes, for a purely prospective test: nothing here is fitted, and no
threshold, score or likelihood ratio is computed at this stage. SBDS (E16b, run with
--gene SBDS) is a fifth deposit of the same platform and licence; its files carry the
prefix sbds_ where E16's carry newgenes_, and a run without the argument writes E16's
files as before.

What it does, offline, from files scripts/stage_E16/fetch_mavedb.py,
fetch_reference.py and (E16b) fetch_igvf.py put on disk:

  records   newgenes_record_check.csv and docs/newgenes-record-check.md: each score
            set's URN, title, licence, publication links, dates and versions, target
            transcript and row count, read from the metadata MaveDB served, and for
            the two URNs given in advance checked against the author's record.
  variants  <gene>_variants.parquet: every intron-side SNV of the score table
            (1 <= |offset|; |offset| 1-2 and beyond 50 kept and flagged), with GRCh38
            coordinates, signed offset, side, stratum, score, range class and label.
  labels    the investigator-provided score ranges of the deposit, applied
            mechanically: 'abnormal' -> 1, 'normal' -> 0, anything else NA (the seven
            genes' rule for intermediate classes, config.ASSAY_LABELS); a score in no
            range, or on a point two ranges share, is 'unassigned' (NA). The second
            case is the author's rule for CTCF, decided 2026-10-20: its Indeterminate
            and Functionally Normal ranges both include their common end point
            (docs/newgenes-record-check.md gives it). SBDS has no such ranges: its
            classes are copied from the score table's per-variant class column,
            mapped by wording through CLASS_NAME_MAP, by the order of sources the
            author fixed on 2026-10-20 before any predictor result (label_source).
  counts    newgenes_counts_by_offset.csv, newgenes_label_summary.csv.
  checks    newgenes_coordinate_check.csv, newgenes_offset_sample.csv,
            newgenes_controls.csv.
  E16b      also sbds_igvf_crosscheck.csv (the MaveDB classes against the IGVF copy
            of the same experiment, which never gives a label) and
            label_version_check.csv, label_version_diff.csv (the BARD1 and PALB2
            label files against every archived version of their score sets). Both
            are checks only: no label, threshold or analysis-set row changes.

Coordinates. MaveDB holds no GRCh38 postMapped record for any intron-side SNV of
these five score sets (their mapping state is 'incomplete'), so they are made the
way the DDX3X deposit's were (evid_external.prepare_ddx3x): the atlas mapping
module maps each c. description through the Mutalyzer transcript model, anchored to
GRCh38 by the Ensembl lookup of the MANE Select partner, and checks the stated
reference base against an Ensembl sequence slice. All three are read from the
tracked reference cache; this stage never opens a network connection. One check is
independent of the mapping: each intron-side coordinate is recomputed from MaveDB's
own postMapped coordinate of the exonic base its c. position names, plus the signed
offset. Another only checks the mapping's transcript model against itself: the
distance to the exon boundary, computed from that model's exon structure, must equal
|offset|, on the nearer side. It cannot catch a wrongly anchored exon; the MaveDB
check and the reference bases can. scripts/stage_E16/check_reference.py and the
tests check every reference base against the Ensembl release-112 chromosome. Where a
score table has the authors' own chrom/pos/ref/alt columns (SBDS), or MaveDB a
postMapped allele for an intron-side SNV, each is compared too
(extra_coordinate_checks).

Run (PYTHONPATH=phase1, from phase1/):
    python -m src.evid_newgenes               # E16
    python -m src.evid_newgenes --gene SBDS   # E16b
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from . import config as C
from . import evid_common as K
from . import evid_external as E
from .evid_build_set import _norm_hgvs

GENES = ["CTCF", "SFPQ", "RAD51D", "XRCC2"]
EXT_DIR = Path("data/evidence/external")
REF_CACHE = Path("data/evidence/reference_cache")
REPORT_DIR = Path("reports/evidence")
DOC_PATH = Path("../docs/newgenes-record-check.md")
# E16 (no argument) and E16b (--gene SBDS): the genes, the prefix of the report files,
# the record page, and whether the summary carries a total row over the genes
RUNS = {
    None: {"genes": GENES, "prefix": "newgenes", "doc": DOC_PATH, "total": "ALL4"},
    "SBDS": {"genes": ["SBDS"], "prefix": "sbds", "doc": Path("../docs/sbds-record-check.md"),
             "total": None},
}
RAW = {g: EXT_DIR / f"{g.lower()}_raw" for g in GENES + ["SBDS"]}

# the DDX3X pattern, on the transcripts these deposits write hgvs_nt against (RefSeq,
# or for SBDS the MANE Select Ensembl transcript)
TX = r"(?:N[MR]_|ENST)[\d.]+"
HGVS_NT = E.intronic_snv_regex(TX)
C_POS = re.compile(r"^c\.([-*]?\d+)([+-]\d+)[ACGT]>[ACGT]$")
EXONIC_C = re.compile(rf"^({TX}):c\.([-*]?\d+)[ACGT]>[ACGT]$")
G_HGVS = re.compile(r"^(NC_\d+\.\d+):g\.(\d+)([ACGT])>([ACGT])$")
DEEP = "out_of_scope_deep"
MAX_DEPTH = 50          # the analysis window, 1 <= |offset| <= 50 (evid_build_set)
LABEL_OF = {"abnormal": 1.0, "normal": 0.0}       # anything else: NA, as config.ASSAY_LABELS
UNASSIGNED = "unassigned"
ON_SHARED_POINT = "score on a point two ranges share (author's rule, 2026-10-20)"
N_SAMPLE = 10
SNV_NT = re.compile(rf"^{TX}:c\.[-*]?\d+([+-]\d+)?[ACGT]>[ACGT]$")
# The author's record for the URNs given in advance (E16.1, E16b.1): URN, licence,
# publication status (public; no paper or preprint linked to the score set or its
# experiment) and the month of first publication. The stage stops if the metadata
# MaveDB served disagrees with any of them.
_NO_PAPER = {"private": False, "score_set_publications": "not_provided",
             "experiment_publications": "not_provided"}
AUTHOR_RECORD = {
    "CTCF": {"urn": "urn:mavedb:00001262-a-3", "license": "CC0", **_NO_PAPER,
             "first_published": "2026-01"},
    "SFPQ": {"urn": "urn:mavedb:00001265-a-2", "license": "CC0", **_NO_PAPER,
             "first_published": "2026-01"},
    "SBDS": {"urn": "urn:mavedb:00001284-a-1", "license": "CC0", **_NO_PAPER,
             "first_published": "2026-08"},
}
AUTHOR_RECORD_DATE = {"CTCF": "2026-10-19", "SFPQ": "2026-10-19", "SBDS": "2026-10-20"}
# Where a deposit's classes come from (E16b.0, fixed by the author on 2026-10-20,
# before any predictor result): (a) exactly one investigator-provided set of ranges,
# checked variant by variant against a class column if the score table has one;
# (b) else a per-variant class column of the score table, copied as it stands;
# (c) else stop. The class names of a column map by wording only, through this table;
# any other value is 'unassigned' and listed.
CLASS_COLUMNS = ("functional_consequence",)
CLASS_NAME_MAP = {"functionally_abnormal": "abnormal", "functionally_normal": "normal",
                  "indeterminate": "not_specified"}


# ---------------------------------------------------------------------------
# inputs
# ---------------------------------------------------------------------------
def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _here(repo_path: str) -> Path:
    """A repository-relative path, as the manifests record them, from phase1/."""
    return Path(repo_path).relative_to("phase1")


def load_raw(gene: str) -> dict:
    """The raw manifest, with every file it lists checked against its sha256."""
    man = json.loads((RAW[gene] / "manifest.json").read_text())
    for entry in man["files"]:
        path = _here(entry["path"])
        if not path.is_file() or _sha256(path) != entry["sha256"]:
            raise SystemExit(f"[E16] {entry['path']} is missing or differs from its "
                             "manifest; rerun scripts/stage_E16/fetch_mavedb.py")
    for key, item in man["transcript_model"]["cache_files"].items():
        path = _here(item["path"])
        if not path.is_file() or _sha256(path) != item["sha256"]:
            raise SystemExit(f"[E16] {gene} {key} {item['path']} is missing or differs; "
                             "rerun scripts/stage_E16/fetch_reference.py")
    return man


def _matches(x: float, c: dict) -> bool:
    lo, hi = c["range"]
    ok_lo = lo is None or (x >= lo if c["inclusiveLowerBound"] else x > lo)
    ok_hi = hi is None or (x <= hi if c["inclusiveUpperBound"] else x < hi)
    return ok_lo and ok_hi


def shared_points(classes: list[dict]) -> list[dict]:
    """Points where two ranges both include a shared end, and any real overlap."""
    out = []
    for i, a in enumerate(classes):
        for b in classes[i + 1:]:
            a_lo = -math.inf if a["range"][0] is None else a["range"][0]
            a_hi = math.inf if a["range"][1] is None else a["range"][1]
            b_lo = -math.inf if b["range"][0] is None else b["range"][0]
            b_hi = math.inf if b["range"][1] is None else b["range"][1]
            lo, hi = max(a_lo, b_lo), min(a_hi, b_hi)
            if lo < hi:
                out.append({"classes": [a["label"], b["label"]], "interval": [lo, hi]})
            elif lo == hi and _matches(lo, a) and _matches(lo, b):
                out.append({"classes": [a["label"], b["label"]], "point": lo})
    return out


def require_mechanical(gene: str, classes: list[dict]) -> None:
    """Stop (E16.10) when two ranges share more than an end point, or an open end lies
    inside the order: the author's rule of 2026-10-20 covers a shared point only."""
    overlap = [p for p in shared_points(classes) if "interval" in p]
    if overlap:
        raise SystemExit(f"[E16] {gene}: investigator-provided ranges overlap {overlap}; "
                         "they cannot be applied mechanically (E16.10)")


def classify(score: str, classes: list[dict]) -> tuple[str, str]:
    """(range_class, why) for one score as the deposit wrote it."""
    if score in (None, "", "NA"):
        return UNASSIGNED, "no score"
    x = float(score)
    hit = [c for c in classes if _matches(x, c)]
    if len(hit) == 1:
        return hit[0]["label"], ""
    if not hit:
        return UNASSIGNED, "score outside every range"
    return UNASSIGNED, ON_SHARED_POINT


def column_class(value: str | None) -> tuple[str, str]:
    """(range_class, why) for a value of a per-variant class column: the deposit's own
    word if the class-name table knows it, else 'unassigned'."""
    if value in (None, "", "NA"):
        return UNASSIGNED, "no class in the column"
    if value in CLASS_NAME_MAP:
        return value, ""
    return UNASSIGNED, f"class value not in the class-name table: {value}"


def label_source(gene: str, table: list[dict], ranges_blob: dict) -> dict:
    """The deposit's classes by the rule fixed before any predictor result (E16b.0):
    ranges if exactly one investigator-provided set exists (checked against a class
    column variant by variant), else the class column, else stop."""
    inv = ranges_blob["investigator_provided"]
    col = next((c for c in CLASS_COLUMNS if table and c in table[0]), None)
    if len(inv) > 1:
        raise SystemExit(f"[E16] {gene}: {len(inv)} investigator-provided calibrations; "
                         "exactly one is required (E16.4)")
    if len(inv) == 1:
        classes = inv[0]["functionalClassifications"]
        require_mechanical(gene, classes)
        fc = {c["label"]: c["functionalClassification"] for c in classes}
        if col:
            differ = []
            for r in table:
                by_range = fc.get(classify(r["score"], classes)[0])
                by_col = CLASS_NAME_MAP.get(r.get(col))
                if by_range != by_col:
                    differ.append((r["accession"], r["score"], r.get(col)))
            if differ:
                raise SystemExit(f"[E16] {gene}: the investigator-provided ranges and the "
                                 f"score-table column {col} disagree on {len(differ)} variants, "
                                 f"e.g. {differ[:3]} (E16b.11)")
        return {"kind": "ranges", "calibration": inv[0], "classes": classes, "fc": fc,
                "column": col}
    if col:
        classes = [{"label": v, "functionalClassification": f} for v, f in CLASS_NAME_MAP.items()]
        return {"kind": "column", "calibration": None, "classes": classes,
                "fc": dict(CLASS_NAME_MAP), "column": col}
    raise SystemExit(f"[E16] {gene}: no investigator-provided ranges and no per-variant "
                     "class column; no labels can be taken (E16b.11)")


def require_known_classes(gene: str, df: pd.DataFrame) -> None:
    """Stop (E16b.11) when an intron-side SNV within 50 nt carries a class-column value
    outside the class-name table: it would be unassigned, and a count would move."""
    bad = df[(df["range_class"] == UNASSIGNED) & (df["offset"].abs() <= MAX_DEPTH)
             & df["range_class_note"].str.startswith("class value not in")]
    if len(bad):
        raise SystemExit(f"[E16b] {gene}: {len(bad)} intron-side SNVs carry a class outside "
                         f"the class-name table, e.g. {bad['range_class_note'].iloc[0]} (E16b.11)")


# ---------------------------------------------------------------------------
# coordinates
# ---------------------------------------------------------------------------
def _atlas_mapping():
    sys.path.insert(0, str(K.ATLAS_REPO / "src"))
    from atlas import mapping as M
    M.REF_CACHE = REF_CACHE

    def offline(url, **_kw):
        raise SystemExit(f"[E16] the stage tried the network ({url}); run "
                         "scripts/stage_E16/fetch_reference.py first")
    M.http_get = offline
    return M


def anchors_from_mavedb(gene: str, man: dict) -> tuple[dict, dict]:
    """c. position of each exonic SNV -> its GRCh38 position, from MaveDB's own
    postMapped records; and the reference accessions those records name."""
    acc = man["transcript_model"]["mane_select"]["grch38_chr_accession"]
    pos, refs = {}, {}
    with (RAW[gene] / "variant_details.ndjson").open() as f:
        for line in f:
            if not line.strip():
                continue
            r = json.loads(line)
            post = r.get("postMapped")
            m = EXONIC_C.match(r.get("targetHgvs") or "")
            # current records with a postMapped allele only, as scripts/13_probe_postmapped
            if r.get("isCurrent") is not True or not post or not m:
                continue
            for e in post.get("expressions", []):
                g = G_HGVS.match(e.get("value", "")) if e.get("syntax") == "hgvs.g" else None
                if g:
                    refs[g.group(1)] = refs.get(g.group(1), 0) + 1
                    if g.group(1) == acc:
                        pos.setdefault(f"c.{m.group(2)}", set()).add(int(g.group(2)))
    return pos, refs


def exon_distances(pos: int, exons: list, strand: int) -> tuple[int, int] | None:
    """(distance to the donor-side boundary, to the acceptor-side boundary) of an
    intronic position, transcript orientation, from the exon structure."""
    ex = sorted(exons, key=lambda e: e.g_start)
    for a, b in zip(ex, ex[1:]):
        if a.g_end < pos < b.g_start:
            left, right = pos - a.g_end, b.g_start - pos
            return (left, right) if strand == 1 else (right, left)
    return None


def build_gene(gene: str, M) -> tuple[pd.DataFrame, dict]:
    man = load_raw(gene)
    meta = json.loads((RAW[gene] / "scoreset_meta.json").read_text())
    target = meta["targetGenes"][0]["targetAccession"]["accession"]
    csv_entry = next(e for e in man["files"] if e["path"].endswith(".csv"))
    table = list(csv.DictReader(_here(csv_entry["path"]).open(newline="")))
    src = label_source(gene, table, json.loads((RAW[gene] / "score_ranges.json").read_text()))
    cal, classes, fc = src["calibration"], src["classes"], src["fc"]

    tm = man["transcript_model"]
    mane = M.load_mane_records([gene])[gene]
    tmap, prov = M.build_transcript_map(target, mane)
    lo = min(e.g_start for e in tmap.exons) - M.REGION_PAD
    hi = max(e.g_end for e in tmap.exons) + M.REGION_PAD
    if [lo, hi] != [tm["region"]["start"], tm["region"]["end"]]:
        raise SystemExit(f"[E16] {gene}: region {lo}-{hi} is not the one fetched")
    s0, _, seq = M.fetch_ensembl_region(mane.chr_name, lo, hi)
    rseq = M.RegionSeq(s0, seq)
    strand = tmap.strand_sign
    anchors, anchor_refs = anchors_from_mavedb(gene, man)

    rna = "rna_score" in (table[0] if table else {})
    rows, failures = [], {}
    for r in table:
        m = HGVS_NT.search(r["hgvs_nt"] or "")
        if not m:
            continue
        hgvs, c = r["hgvs_nt"], m.group(2)
        depth = int(m.group(3))
        sign = hgvs[m.start(3) - 1]
        offset = depth if sign == "+" else -depth
        try:
            pos, ref, alt = M.map_variant(tmap, rseq, c)
        except M.MappingFailure as e:
            failures[f"{e.category}: {hgvs}"] = str(e)
            continue
        range_class, why = (classify(r["score"], classes) if src["kind"] == "ranges"
                            else column_class(r.get(src["column"])))
        anchor_c = "c." + C_POS.match(c).group(1)
        a = anchors.get(anchor_c)
        anchor_pos = next(iter(a)) if a and len(a) == 1 else None
        expected = None if anchor_pos is None else anchor_pos + strand * offset
        d = exon_distances(pos, tmap.exons, strand)
        stratum = E.stratum_of(depth) if depth <= MAX_DEPTH else DEEP
        in_scope = stratum in K.IN_SCOPE_STRATA
        rows.append({
            "gene": gene, "urn": man["urn"], "variant_urn": r["accession"],
            "hgvs_nt": hgvs, "hgvs_nt_norm": _norm_hgvs(hgvs),
            "chrom": mane.chr_name, "pos": int(pos), "ref": ref, "alt": alt,
            "variant_id": f"{mane.chr_name}-{pos}-{ref}-{alt}",
            "transcript": m.group(1), "offset": offset,
            "side": "donor" if sign == "+" else "acceptor",
            "stratum": stratum, "in_scope": in_scope,
            "scope": "pvs1_pm12" if stratum == "pm12" else ("in_scope" if in_scope else DEEP),
            "score": float(r["score"]) if r["score"] not in ("", "NA", None) else np.nan,
            "range_class": range_class, "range_class_note": why,
            "label": LABEL_OF.get(fc.get(range_class), np.nan),
            "anchor_c": anchor_c,
            "anchor_pos_mavedb": anchor_pos,
            "pos_from_anchor": expected,
            "anchor_agrees": None if expected is None else bool(expected == pos),
            "exon_distance": None if d is None else (d[0] if sign == "+" else d[1]),
            "exon_distance_other_side": None if d is None else (d[1] if sign == "+" else d[0]),
            "rna_score_present": (rna and r.get("rna_score") not in ("", "NA", None)),
        })
    if failures:
        raise SystemExit(f"[E16] {gene}: {len(failures)} intron-side SNVs failed to map "
                         f"(E16.10): {list(failures.items())[:3]}")
    df = pd.DataFrame(rows)
    n_before = len(df)
    dup = df.duplicated(["chrom", "pos", "ref", "alt"], keep="first")
    df = df[~dup].sort_values(["pos", "alt"], kind="mergesort").reset_index(drop=True)
    for col in ("anchor_pos_mavedb", "pos_from_anchor", "exon_distance",
                "exon_distance_other_side"):
        df[col] = df[col].astype("Int64")
    df["anchor_agrees"] = df["anchor_agrees"].astype("boolean")

    info = {
        "man": man, "meta": meta, "calibration": cal, "classes": classes,
        "shared_points": shared_points(classes) if src["kind"] == "ranges" else [],
        "target": target, "mane": mane,
        "provenance": prov, "strand": strand, "anchor_reference_accessions": anchor_refs,
        "n_rows": len(table), "n_intronic": n_before, "n_duplicates": int(dup.sum()),
        "rna_column": rna, "table": table, "fc": fc, "label_source": src,
        "extra_coordinate_checks": extra_coordinate_checks(gene, man, table, df),
    }
    return df, info


def extra_coordinate_checks(gene: str, man: dict, table: list[dict], df: pd.DataFrame) -> dict:
    """Two further sources for the coordinates, where a deposit has them (E16b.3): the
    authors' own chrom/pos/ref/alt columns of the score table, and MaveDB's postMapped
    GRCh38 allele of an intron-side SNV. Each is compared, variant by variant, with the
    mapped coordinate; any difference stops the stage. Neither replaces it."""
    ours = {r.hgvs_nt: (str(r.chrom), int(r.pos), r.ref, r.alt) for r in df.itertuples()}
    out = {}
    if table and all(k in table[0] for k in ("chrom", "pos", "ref", "alt")):
        n, bad = 0, []
        for r in table:
            if r["hgvs_nt"] in ours:
                n += 1
                theirs = (r["chrom"].removeprefix("chr"), int(r["pos"]), r["ref"], r["alt"])
                if theirs != ours[r["hgvs_nt"]]:
                    bad.append({"hgvs_nt": r["hgvs_nt"], "ours": ours[r["hgvs_nt"]], "authors": theirs})
        out["author_columns"] = {"n_compared": n, "n_agree": n - len(bad), "differ": bad[:20]}
        if bad:
            raise SystemExit(f"[E16] {gene}: {len(bad)} intron-side SNVs differ from the "
                             f"authors' chrom/pos/ref/alt, e.g. {bad[:3]} (E16b.11)")
    acc = man["transcript_model"]["mane_select"]["grch38_chr_accession"]
    n, bad = 0, []
    with (RAW[gene] / "variant_details.ndjson").open() as f:
        for line in f:
            if not line.strip():
                continue
            r = json.loads(line)
            th = r.get("targetHgvs") or ""
            if r.get("isCurrent") is not True or th not in ours or not r.get("postMapped"):
                continue
            for e in r["postMapped"].get("expressions", []):
                g = G_HGVS.match(e.get("value", "")) if e.get("syntax") == "hgvs.g" else None
                if g and g.group(1) == acc:
                    n += 1
                    theirs = (ours[th][0], int(g.group(2)), g.group(3), g.group(4))
                    if theirs != ours[th]:
                        bad.append({"hgvs_nt": th, "ours": ours[th], "postMapped": theirs})
    out["intronic_postmapped"] = {"n_compared": n, "n_agree": n - len(bad), "differ": bad[:20]}
    if bad:
        raise SystemExit(f"[E16] {gene}: {len(bad)} intron-side SNVs differ from MaveDB's "
                         f"postMapped coordinate, e.g. {bad[:3]} (E16b.11)")
    return out


# ---------------------------------------------------------------------------
# outputs
# ---------------------------------------------------------------------------
def _chain(gene: str, meta: dict, man: dict) -> list[dict]:
    """The score set's versions, oldest first: this one and those it supersedes."""
    old = [json.loads((RAW[gene] / f"scoreset_meta_{urn.split(':')[-1]}.json").read_text())
           for urn in man["superseded_versions"]]
    return sorted(old, key=lambda m: m["publishedDate"] or "") + [meta]


def _first_published(gene: str, meta: dict, man: dict) -> str:
    return min(m["publishedDate"] for m in _chain(gene, meta, man) if m["publishedDate"])


def _pubs(obj: dict) -> str:
    ids = [f"{p.get('dbName', '?')}:{p.get('identifier')}"
           for k in ("primaryPublicationIdentifiers", "secondaryPublicationIdentifiers")
           for p in obj.get(k) or []]
    return "; ".join(ids) if ids else "not_provided"


def record_row(gene: str, info: dict) -> dict:
    man, meta = info["man"], info["meta"]
    entry = next(e for e in man["files"] if e["path"].endswith("/scoreset_meta.json"))
    csv_entry = next(e for e in man["files"] if e["path"].endswith(".csv"))
    exp = json.loads((RAW[gene] / "experiment_meta.json").read_text())
    acc = meta["targetGenes"][0]["targetAccession"]
    lic = man["license"]
    return {
        "gene": gene, "urn": meta["urn"], "title": meta["title"],
        "short_description": meta["shortDescription"],
        "license": f"{lic['shortName']} {lic['version']} ({lic['link']})",
        "publication_identifiers": _pubs(meta),
        "created_date": meta["creationDate"],
        "first_published_date": _first_published(gene, meta, man),
        "current_version_date": meta["publishedDate"],
        "target_gene": meta["targetGenes"][0]["name"],
        "target_accession": acc["accession"], "target_assembly": acc.get("assembly", "not_provided"),
        "n_variants_total": meta["numVariants"],
        "score_ranges_present": info["calibration"] is not None,
        "fetched_utc": entry["fetched_utc"], "source_url": entry["source_url"],
        "raw_sha256": entry["sha256"],
        # beyond the fields E16.1 asks for
        "experiment_urn": exp["urn"],
        "experiment_publication_identifiers": _pubs(exp),
        "superseded_versions": "; ".join(man["superseded_versions"]) or "none",
        "mapping_state": meta.get("mappingState"),
        "n_score_rows": info["n_rows"], "scores_csv_sha256": csv_entry["sha256"],
        "calibration_title": (info["calibration"]["title"] if info["calibration"] is not None
                              else f"not_provided (classes from the score-table column "
                                   f"{info['label_source']['column']})"),
        "api_version": entry["api_version"],
        "mavedb_private": meta["private"], "processing_state": meta["processingState"],
        "selection": man["selection"]["how"],
        "author_record_check": author_check(gene, meta, exp, man),
    }


def author_check(gene: str, meta: dict, exp: dict, man: dict) -> str:
    """The fetched record against the author's of 2026-10-19, where there is one."""
    want = AUTHOR_RECORD.get(gene)
    if want is None:
        return "no author record: URN found by search (the manifest's selection)"
    got = {"urn": meta["urn"], "license": meta["license"]["shortName"],
           "private": meta["private"], "score_set_publications": _pubs(meta),
           "experiment_publications": _pubs(exp),
           "first_published": _first_published(gene, meta, man)[:7]}
    bad = [f"{k} {got[k]!r}, the author's record {want[k]!r}" for k in want if got[k] != want[k]]
    if bad:
        raise SystemExit(f"[E16] {gene}: {'; '.join(bad)} (E16.10)")
    return "agrees: " + ", ".join(f"{k} {v}" for k, v in want.items())


def _readout(info: dict, df: pd.DataFrame) -> str:
    if not info["rna_column"]:
        return "fitness only (the score table has no rna_score column)"
    n = int(df["rna_score_present"].sum())
    if n == 0:
        return ("fitness only (the score table has an rna_score column, but no "
                "intron-side row carries a value)")
    return f"fitness; rna_score present on {n} intron-side rows"


def _bins_per_nt(depth: int) -> str:
    return str(depth) if depth <= MAX_DEPTH else f">{MAX_DEPTH}"


def counts_tables(frames: dict, infos: dict, total: str | None = "ALL4"
                  ) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows, summ = [], []
    per_nt = [str(i) for i in range(1, MAX_DEPTH + 1)] + [f">{MAX_DEPTH}"]
    strata = ["pm12", "s3_10", "s11_50", "s3_50", DEEP]
    for gene, df in frames.items():
        classes = [c["label"] for c in infos[gene]["classes"]] + [UNASSIGNED]
        depth = df["offset"].abs()
        for b in per_nt:
            sub = df[depth.map(_bins_per_nt) == b]
            for cl in classes:
                rows.append({"gene": gene, "scheme": "per_nt", "bin": b,
                             "range_class": cl, "n": int((sub["range_class"] == cl).sum())})
        for st in strata:
            sub = K.stratum_frame(df[df["stratum"] != DEEP], st) if st != DEEP \
                else df[df["stratum"] == DEEP]
            for cl in classes:
                rows.append({"gene": gene, "scheme": "stratum", "bin": st,
                             "range_class": cl, "n": int((sub["range_class"] == cl).sum())})
    by_offset = pd.DataFrame(rows)

    def summarise(gene: str, st: str, sub: pd.DataFrame) -> dict:
        lab = sub["label"]
        n_dam, n_norm = int((lab == 1).sum()), int((lab == 0).sum())
        n_un = int((sub["range_class"] == UNASSIGNED).sum())
        n_unc = int(len(sub) - n_dam - n_norm - n_un)
        labelled = n_dam + n_norm
        return {"gene": gene, "stratum": st, "n_total": int(len(sub)),
                "n_damaging": n_dam, "n_normal": n_norm, "n_uncertain": n_unc,
                "n_unassigned": n_un,
                "damaging_rate": n_dam / labelled if labelled else np.nan,
                "evaluable_min10": bool(n_dam >= K.MIN_POS and n_norm >= K.MIN_NEG)}
    allf = pd.concat(frames.values(), ignore_index=True)
    for gene, df in list(frames.items()) + ([(total, allf)] if total else []):
        for st in strata:
            sub = K.stratum_frame(df[df["stratum"] != DEEP], st) if st != DEEP \
                else df[df["stratum"] == DEEP]
            summ.append(summarise(gene, st, sub))
    return by_offset, pd.DataFrame(summ)


def coordinate_tables(frames: dict, infos: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    genes = list(frames)
    rows, samples = [], []
    for gene, df in frames.items():
        info = infos[gene]
        ach = df["anchor_agrees"]
        same_side = df["exon_distance"] <= df["exon_distance_other_side"]
        rows.append({
            "gene": gene, "chrom": info["mane"].chr_name,
            "strand": "+" if info["strand"] == 1 else "-",
            "n_intronic_snv": info["n_intronic"], "n_mapped": int(len(df)),
            "n_duplicates_dropped": info["n_duplicates"],
            "n_anchor_checked": int(ach.notna().sum()),
            "n_anchor_agree": int(ach.fillna(False).sum()),
            "n_no_mavedb_anchor": int(ach.isna().sum()),
            "n_exon_distance_checked": int(df["exon_distance"].notna().sum()),
            "n_exon_distance_equals_offset": int(
                (df["exon_distance"] == df["offset"].abs()).fillna(False).sum()),
            "n_offset_side_is_nearest": int(same_side.fillna(False).sum()),
            "anchor_reference_accessions": "; ".join(
                f"{k} ({v})" for k, v in sorted(info["anchor_reference_accessions"].items())),
            "grch38_chr_accession": info["mane"].grch38_chr_accession,
        })
        rng = np.random.default_rng([C.RANDOM_SEED, genes.index(gene)])
        idx = np.sort(rng.choice(len(df), size=min(N_SAMPLE, len(df)), replace=False))
        s = df.iloc[idx]
        samples += [{
            "gene": gene, "hgvs_nt": v.hgvs_nt, "chrom": v.chrom, "pos": int(v.pos),
            "ref": v.ref, "alt": v.alt, "offset": int(v.offset), "side": v.side,
            "exon_distance": None if pd.isna(v.exon_distance) else int(v.exon_distance),
            "exon_distance_other_side": (None if pd.isna(v.exon_distance_other_side)
                                         else int(v.exon_distance_other_side)),
            "anchor_c": v.anchor_c,
            "anchor_pos_mavedb": None if pd.isna(v.anchor_pos_mavedb) else int(v.anchor_pos_mavedb),
            "pos_from_anchor": None if pd.isna(v.pos_from_anchor) else int(v.pos_from_anchor),
            "anchor_agrees": None if pd.isna(v.anchor_agrees) else bool(v.anchor_agrees),
        } for v in s.itertuples()]
    return pd.DataFrame(rows), pd.DataFrame(samples)


def _variant_type(hgvs_nt: str) -> str:
    if SNV_NT.match(hgvs_nt or ""):
        return "snv"
    return "deletion" if (hgvs_nt or "").endswith("del") else "other"


def controls_table(infos: dict) -> pd.DataFrame:
    """Synonymous and nonsense variants of each deposit under its ranges, found with
    the DDX3X control patterns. These score tables leave hgvs_pro empty, so the protein
    consequence is read from MaveDB's own variant-details record (mapped.protein, its
    projection of each current mapped variant); the score is the score table's. Each
    count is split by the score table's nucleotide change: a 3-bp deletion can project
    to a stop as an SNV does, and only the SNV rows compare with the intron-side SNVs."""
    rows = []
    for gene, info in infos.items():
        score = {r["accession"]: r["score"] for r in info["table"]}
        vtype = {r["accession"]: _variant_type(r["hgvs_nt"]) for r in info["table"]}
        n_pro = sum(1 for r in info["table"] if (r.get("hgvs_pro") or "NA") != "NA")
        src = info["label_source"]
        by_column = ({r["accession"]: column_class(r.get(src["column"]))[0] for r in info["table"]}
                     if src["kind"] == "column" else None)
        counts, n_proj = {}, 0
        with (RAW[gene] / "variant_details.ndjson").open() as f:
            for line in f:
                if not line.strip():
                    continue
                r = json.loads(line)
                pro = (((r.get("mapped") or {}).get("protein") or {}).get("hgvs") or "")
                if r.get("isCurrent") is not True or not pro:
                    continue
                n_proj += 1
                kind = ("synonymous" if E._PRO_SYN.search(pro) else
                        "nonsense" if E._PRO_NON.search(pro) else None)
                if kind is None:
                    continue
                cl = (by_column.get(r["urn"], UNASSIGNED) if by_column is not None
                      else classify(score.get(r["urn"]), info["classes"])[0])
                key = (kind, vtype.get(r["urn"], "other"), cl)
                counts[key] = counts.get(key, 0) + 1
        for kind in ("synonymous", "nonsense"):
            for vt in ("snv", "deletion", "other"):
                for cl in [c["label"] for c in info["classes"]] + [UNASSIGNED]:
                    rows.append({"gene": gene, "control": kind, "variant_type": vt,
                                 "range_class": cl, "n": counts.get((kind, vt, cl), 0),
                                 "score_table_rows_with_hgvs_pro": n_pro,
                                 "variant_details_with_protein_projection": n_proj})
    return pd.DataFrame(rows)


def _src(gene: str, man: dict, *names: str) -> str:
    """'[`<gene>_raw/<file>`, fetched <time>; ...]' for files the raw manifest lists."""
    entry = {Path(e["path"]).name: e for e in man["files"]}
    parts = []
    for name in names:
        e, g = entry[name], gene.lower()
        if "derived_from" in e:
            src = entry[Path(e["derived_from"]).name]
            parts.append(f"`{g}_raw/{name}`, extracted from `{g}_raw/"
                         f"{Path(e['derived_from']).name}`, fetched {src['fetched_utc']}")
        else:
            parts.append(f"`{g}_raw/{name}`, fetched {e['fetched_utc']}")
    return "[" + "; ".join(parts) + "]"


def _interval(c: dict) -> str:
    lo, hi = c["range"]
    return (("(-inf" if lo is None else ("[" if c["inclusiveLowerBound"] else "(") + repr(lo))
            + ", " + ("inf)" if hi is None else repr(hi) + ("]" if c["inclusiveUpperBound"] else ")")))


def write_doc(records: pd.DataFrame, frames: dict, infos: dict, doc_path: Path = DOC_PATH,
              title: str = "Four further external genes: MaveDB records", stage: str = "E16",
              preface: list[str] | None = None) -> None:
    lines = [f"# {title}",
             "",
             f"Written by `phase1/src/evid_newgenes.py` (stage {stage}) from the files "
             "`scripts/stage_E16/fetch_mavedb.py` saved; do not edit by hand. Each fact "
             "is followed by the file it was read from (under "
             "`phase1/data/evidence/external/`) and the time that file was fetched.",
             ""] + (preface or [])
    for r in records.itertuples():
        info = infos[r.gene]
        man, meta, df = info["man"], info["meta"], frames[r.gene]
        meta_src = _src(r.gene, man, "scoreset_meta.json")
        chain = _chain(r.gene, meta, man)
        chain_src = _src(r.gene, man, *[
            "scoreset_meta.json" if m is meta else f"scoreset_meta_{m['urn'].split(':')[-1]}.json"
            for m in chain])
        versions = "; ".join(
            f"{m['urn']}{' (this one)' if m is meta else ''}, published {m['publishedDate']}, "
            f"{m['numVariants']:,} variants, mapping {m['mappingState']}, investigator-provided "
            "classes " + "/".join(str(len(c.get("functionalClassifications") or []))
                                  for c in m.get("scoreCalibrations") or []
                                  if c.get("investigatorProvided"))
            for m in chain)
        src = info["label_source"]
        classes = ("; ".join(f"{c['label']} ({c['functionalClassification']}) {_interval(c)}"
                             for c in info["classes"]) if src["kind"] == "ranges" else "")
        shared = "; ".join(
            f"{' and '.join(p['classes'])} share "
            + (f"the point {p['point']!r}" if "point" in p else f"the interval {p['interval']}")
            for p in info["shared_points"]) or "none"
        on_shared = int((df["range_class_note"] == ON_SHARED_POINT).sum())
        csv_name = f"urn_mavedb_{r.urn.split(':')[-1]}.csv"
        lines += [f"## {r.gene}", "",
                  f"- Score set: {r.urn}, \"{r.title}\" {meta_src}",
                  f"- Description: {r.short_description} {meta_src}",
                  f"- Status: {'private' if r.mavedb_private else 'public'} (private = "
                  f"{str(r.mavedb_private).lower()}), processing "
                  f"state {r.processing_state}, mapping state {r.mapping_state} {meta_src}",
                  f"- Licence: {r.license} {meta_src}",
                  f"- Publication links: score set {r.publication_identifiers} {meta_src}; "
                  f"experiment {r.experiment_urn} {r.experiment_publication_identifiers} "
                  f"{_src(r.gene, man, 'experiment_meta.json')}",
                  f"- Dates: this version created {r.created_date}, published "
                  f"{r.current_version_date} {meta_src}; first published "
                  f"{r.first_published_date}, the earliest publication date of its versions "
                  f"{chain_src}",
                  f"- Versions: {versions} {chain_src}"]
        if r.gene in AUTHOR_RECORD:
            lines.append(
                f"- Checked against the author's record of {AUTHOR_RECORD_DATE[r.gene]} (URN, "
                f"licence, publication status, first published "
                f"{AUTHOR_RECORD[r.gene]['first_published']}, the current version dated "
                f"separately): {r.author_record_check}; "
                f"first published {r.first_published_date}, current version published "
                f"{r.current_version_date} {chain_src}")
        else:
            lines.append(
                f"- URN found by MaveDB's search, not given in advance: {r.selection}; the "
                f"requests, the candidates and the reason each was kept or excluded are in "
                f"`{man['selection']['record']}`")
        lines += [f"- Target: {r.target_gene}, {r.target_accession} ({r.target_assembly}) "
                  f"{meta_src}",
                  f"- Variants: {r.n_variants_total:,} in the metadata {meta_src}; "
                  f"{r.n_score_rows:,} rows in the score table {_src(r.gene, man, csv_name)}",
                  ]
        if src["kind"] == "ranges":
            lines += [f"- Investigator-provided ranges, \"{r.calibration_title}\": {classes} "
                      f"{_src(r.gene, man, 'score_ranges.json')}",
                      f"- Shared range ends: {shared}; intron-side SNVs scored there "
                      f"({UNASSIGNED!r} by the author's rule of 2026-10-20): {on_shared} "
                      f"{_src(r.gene, man, 'score_ranges.json', csv_name)}"]
        else:
            seen = sorted({row.get(src["column"]) or "" for row in info["table"]})
            lines += [f"- Investigator-provided ranges: none (no scoreCalibrations entry has "
                      f"investigatorProvided true) {_src(r.gene, man, 'score_ranges.json')}",
                      f"- Classes: copied from the score-table column `{src['column']}`, mapped "
                      "by wording only: " + "; ".join(f"{k} -> {v}" for k, v in CLASS_NAME_MAP.items())
                      + f" (label 1 / 0 / NA); values in the column: {', '.join(seen)} "
                      f"{_src(r.gene, man, csv_name)}"]
        lines.append("")
    doc_path.write_text("\n".join(lines))


CANON = {"abnormal": "abnormal", "normal": "normal", "not_specified": "intermediate"}
IGVF_DIR = EXT_DIR / "igvf_sbds"
LABEL_DIR = EXT_DIR / "label_versions"
LABEL_GENES = ["BARD1", "PALB2"]


def _check_manifest(path: Path) -> dict:
    man = json.loads(path.read_text())
    for e in man["files"]:
        f = _here(e["path"])
        if not f.is_file() or _sha256(f) != e["sha256"]:
            raise SystemExit(f"[E16b] {e['path']} is missing or differs from {path}")
    return man


def igvf_crosscheck(df: pd.DataFrame, info: dict) -> pd.DataFrame:
    """MaveDB against the IGVF copy of the SBDS scores (E16b.5): aligned on the
    normalised c. HGVS; counted as same class, different class, or one side only, over
    all rows and over the intron-side SNVs at 1-50 nt; every difference listed. Labels
    are not touched: they come from MaveDB only."""
    import gzip
    man = _check_manifest(IGVF_DIR / "manifest.json")
    gz = next(_here(e["path"]) for e in man["files"] if e["path"].endswith(".tsv.gz"))
    igvf = {}
    for r in csv.DictReader(gzip.open(gz, "rt", newline=""), delimiter="\t"):
        igvf[_norm_hgvs(r["hgvs_c"])] = r
    col = info["label_source"]["column"]
    mave = {_norm_hgvs(r["hgvs_nt"]): r for r in info["table"]}
    intronic = {_norm_hgvs(h): (int(o), s) for h, o, s in zip(df["hgvs_nt"], df["offset"], df["stratum"])
                if abs(int(o)) <= MAX_DEPTH}
    rows = []
    for scope, keys in (("all rows", set(mave) | set(igvf)),
                        ("intron-side SNVs 1-50 nt", set(intronic))):
        cat = {"same class": 0, "different class": 0, "MaveDB only": 0, "IGVF only": 0}
        for k in sorted(keys):
            a, b = mave.get(k), igvf.get(k)
            if a and b:
                same = (a.get(col) or "") == (b.get("functional_consequence") or "")
                cat["same class" if same else "different class"] += 1
                if not same and scope == "all rows":
                    off, st = intronic.get(k, (None, "exonic or deeper than 50 nt"))
                    rows.append({"row_type": "variant", "scope": scope, "category": "different class",
                                 "n": 1, "hgvs_nt": a["hgvs_nt"], "offset": off, "stratum": st,
                                 "mavedb_score": a["score"], "mavedb_class": a.get(col),
                                 "igvf_score": b["score"], "igvf_class": b.get("functional_consequence")})
            else:
                cat["MaveDB only" if a else "IGVF only"] += 1
        rows[0:0] = [{"row_type": "count", "scope": scope, "category": c, "n": n} for c, n in cat.items()]
    out = pd.DataFrame(rows, columns=["row_type", "scope", "category", "n", "hgvs_nt", "offset",
                                      "stratum", "mavedb_score", "mavedb_class", "igvf_score",
                                      "igvf_class"])
    order = {"count": 0, "variant": 1}
    out = out.sort_values(["row_type", "scope"], key=lambda s: s.map(order) if s.name == "row_type" else s,
                          kind="mergesort").reset_index(drop=True)
    out["offset"] = out["offset"].astype("Int64")
    return out


def _version_classes(table: list[dict], ranges_blob: dict) -> tuple[str, dict]:
    """Each variant's class in one archived version, by the E16b.0 order: ranges if
    exactly one investigator set, else the class column, else none."""
    inv = ranges_blob["investigator_provided"]
    col = next((c for c in CLASS_COLUMNS if table and c in table[0]), None)
    out = {}
    if len(inv) == 1:
        classes = inv[0]["functionalClassifications"]
        fc = {c["label"]: c["functionalClassification"] for c in classes}
        for r in table:
            cl = classify(r["score"], classes)[0]
            out[_norm_hgvs(r["hgvs_nt"])] = (CANON.get(fc.get(cl), UNASSIGNED), r["score"])
        return ("ranges" + (f" (the table also has the column {col})" if col else "")), out
    if col:
        for r in table:
            cl = column_class(r.get(col))[0]
            out[_norm_hgvs(r["hgvs_nt"])] = (CANON.get(CLASS_NAME_MAP.get(cl), UNASSIGNED), r["score"])
        return f"class column {col}", out
    return "no classes", out


def label_version_check() -> tuple[pd.DataFrame, pd.DataFrame]:
    """The training genes' label files against every archived version (E16b.6), on the
    analysis set's intron-side SNVs (1-50 nt). A version written in c. notation is
    aligned on the normalised c. HGVS; one written only in genomic notation (the -a-1
    versions) on the normalised g. HGVS built from the analysis set's GRCh38 allele.
    A check only: no label, threshold or analysis-set row changes."""
    aset = pd.read_parquet(K.SET_PATH)
    checks, diffs = [], []
    for gene in LABEL_GENES:
        spec = C.ASSAY_LABELS[gene]
        lab = pd.read_csv(Path("data/assay_labels") / spec["file"],
                          sep="\t" if spec["file"].endswith(".tsv") else ",", dtype=str)
        ours_class = {_norm_hgvs(h): {1.0: "abnormal", 0.0: "normal", None: "intermediate"}
                      [spec["map"][c]] for h, c in zip(lab[spec["key"]], lab[spec["class_col"]])}
        v = aset[(aset["gene"] == gene) & (aset["intron_offset_abs"].between(1, MAX_DEPTH))]
        ours = [{"c": _norm_hgvs(r.hgvs_c), "g": f"g.{r.pos}{r.ref}>{r.alt}".lower(),
                 "hgvs": r.hgvs_c, "offset": int(r.intron_offset_abs), "stratum": r.stratum,
                 "cls": ours_class.get(_norm_hgvs(r.hgvs_c))} for r in v.itertuples()]
        man = _check_manifest(LABEL_DIR / f"{gene.lower()}_raw" / "manifest.json")
        for ver in man["versions"]:
            vid = ver["urn"].split(":")[-1]
            table = list(csv.DictReader((LABEL_DIR / f"{gene.lower()}_raw" / f"urn_mavedb_{vid}.csv")
                                        .open(newline="")))
            source, arch = _version_classes(table, json.loads(
                (LABEL_DIR / f"{gene.lower()}_raw" / f"score_ranges_{vid}.json").read_text()))
            genomic = all((r["hgvs_nt"] or "").startswith("NC_") for r in table[:50])
            key = "g" if genomic else "c"
            members = {_norm_hgvs(r["hgvs_nt"]) for r in table}
            found = agree = disagree = unclassified = 0
            for o in ours:
                if o[key] not in members:
                    continue
                found += 1
                hit = arch.get(o[key])
                if hit is None or hit[0] == UNASSIGNED or o["cls"] is None:
                    unclassified += 1
                elif hit[0] == o["cls"]:
                    agree += 1
                else:
                    disagree += 1
                    diffs.append({"gene": gene, "version": ver["urn"], "hgvs_nt": o["hgvs"],
                                  "offset": o["offset"], "stratum": o["stratum"],
                                  "our_class": o["cls"], "archive_class": hit[0],
                                  "archive_score": hit[1]})
            if key == "c":
                arch_intronic = set()
                for r in table:
                    m = HGVS_NT.search(r["hgvs_nt"] or "")
                    if m and int(m.group(3)) <= MAX_DEPTH:
                        arch_intronic.add(_norm_hgvs(r["hgvs_nt"]))
                only_archive = len(arch_intronic - {o["c"] for o in ours})
            else:
                only_archive = None              # genomic notation: intron-side not identifiable
            checks.append({"gene": gene, "version": ver["urn"], "current": ver["superseded_by"] is None,
                           "class_source": source, "align_key": "c. HGVS" if key == "c" else
                           "g. HGVS from the analysis set's GRCh38 allele",
                           "n_ours": len(ours), "n_found": found, "n_compared": agree + disagree,
                           "n_agree": agree, "n_disagree": disagree,
                           "n_found_without_class": unclassified,
                           "n_only_ours": len(ours) - found, "n_only_archive": only_archive})
    check = pd.DataFrame(checks)
    check["n_only_archive"] = check["n_only_archive"].astype("Int64")
    diff = pd.DataFrame(diffs, columns=["gene", "version", "hgvs_nt", "offset", "stratum",
                                        "our_class", "archive_class", "archive_score"])
    return check, diff


PRESPECIFIED = [
    "## Pre-specified label source (预先规定)",
    "",
    "Fixed by the author on 2026-10-20, before any predictor was run on SBDS, and not to "
    "be changed afterwards:",
    "",
    "- SBDS labels come from MaveDB only, as for every external gene so far (DDX3X, TP53 "
    "and the four genes of E16). The IGVF copy of the same experiment "
    "(`igvf_sbds/IGVFFI9790OCYB.tsv.gz`) is used only for a cross-check "
    "(`phase1/reports/evidence/sbds_igvf_crosscheck.csv`) and never to build a label.",
    "- Why the source had to be fixed in advance: in the author's SGE survey (2026-10-20) "
    "the 3-10 nt band holds 10 functionally abnormal SNVs by MaveDB and 9 by IGVF, either "
    "side of the evaluability threshold (at least 10 damaging), so the choice of source "
    "could decide whether the band is evaluated at all.",
    "- Order of class sources: (a) exactly one investigator-provided set of ranges, checked "
    "variant by variant against a class column if the score table has one; (b) otherwise a "
    "per-variant class column of the score table, copied as it stands; (c) otherwise stop. "
    "SBDS has no investigator-provided ranges, so (b) applies: the column "
    "`functional_consequence`.",
    "- Class names map by wording only: functionally_abnormal -> 1; functionally_normal -> "
    "0; indeterminate -> NA (row kept); any other value -> unassigned, listed.",
    "",
]


def main(argv: list[str] | None = None) -> None:
    import argparse
    ap = argparse.ArgumentParser(description="E16 (no argument) or E16b (--gene SBDS)")
    ap.add_argument("--gene", choices=[k for k in RUNS if k])
    args = ap.parse_args(argv)
    run = RUNS[args.gene]
    genes, pre, stage = run["genes"], run["prefix"], "E16b" if args.gene else "E16"
    M = _atlas_mapping()
    frames, infos = {}, {}
    for gene in genes:
        df, info = build_gene(gene, M)
        frames[gene], infos[gene] = df, info
        if info["label_source"]["kind"] == "column":
            require_known_classes(gene, df)

    records = pd.DataFrame([record_row(g, infos[g]) for g in genes])
    by_offset, summary = counts_tables(frames, infos, run["total"])
    coords, samples = coordinate_tables(frames, infos)
    controls = controls_table(infos)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    records.to_csv(REPORT_DIR / f"{pre}_record_check.csv", index=False)
    by_offset.to_csv(REPORT_DIR / f"{pre}_counts_by_offset.csv", index=False)
    summary.to_csv(REPORT_DIR / f"{pre}_label_summary.csv", index=False)
    coords.to_csv(REPORT_DIR / f"{pre}_coordinate_check.csv", index=False)
    samples.to_csv(REPORT_DIR / f"{pre}_offset_sample.csv", index=False)
    controls.to_csv(REPORT_DIR / f"{pre}_controls.csv", index=False)
    if args.gene:
        write_doc(records, frames, infos, run["doc"],
                  title=f"{args.gene}: MaveDB record (fifth new external gene)", stage=stage,
                  preface=PRESPECIFIED)
        igvf_crosscheck(frames["SBDS"], infos["SBDS"]).to_csv(
            REPORT_DIR / "sbds_igvf_crosscheck.csv", index=False)
        check, diff = label_version_check()
        check.to_csv(REPORT_DIR / "label_version_check.csv", index=False)
        diff.to_csv(REPORT_DIR / "label_version_diff.csv", index=False)
    else:
        write_doc(records, frames, infos)

    for gene in genes:
        df, info = frames[gene], infos[gene]
        out = EXT_DIR / f"{gene.lower()}_variants.parquet"
        df.drop(columns=["rna_score_present"]).to_parquet(out, index=False)
        man = info["man"]
        src = info["label_source"]
        prov = {
            "gene": gene, "urn": man["urn"], "source": "MaveDB",
            "license": man["license"], "readout": _readout(info, df),
            "inputs": {e["path"]: e["sha256"] for e in man["files"]},
            "coordinates": {
                "source": man["transcript_model"]["route"],
                "mane_select": man["transcript_model"]["mane_select"],
                "transcript_map": info["provenance"],
                "region": man["transcript_model"]["region"],
                "cache_files": man["transcript_model"]["cache_files"],
                "mavedb_postMapped": ("none for any intron-side SNV of this score set; "
                                      "used only as the exonic anchors of the check"),
            },
        }
        if src["kind"] == "ranges":
            prov |= {
                "score_ranges": info["calibration"],
                "range_class": ("the label of the one range the score falls in, in the "
                                "deposit's own wording; 'unassigned' for no score, a score "
                                "outside every range, or a score on a point two ranges share "
                                "(author's rule, 2026-10-20)"),
                "shared_points": info["shared_points"],
                "label_map": {"functionalClassification abnormal": 1,
                              "functionalClassification normal": 0,
                              "any other class, and unassigned": None},
            }
        else:
            extra = info["extra_coordinate_checks"]
            pm = extra["intronic_postmapped"]
            prov["coordinates"]["mavedb_postMapped"] = (
                f"{pm['n_compared']} intron-side SNVs carry a GRCh38 postMapped allele, "
                f"{pm['n_agree']} equal the mapped one; used as an extra check only")
            prov["coordinates"]["extra_checks"] = extra
            values = sorted({r.get(src["column"]) or "" for r in info["table"]})
            prov |= {
                "score_ranges": "none: no scoreCalibrations entry has investigatorProvided true",
                "label_source": (f"the score-table column {src['column']}, not ranges (E16b.0 "
                                 "rule 2b, fixed by the author on 2026-10-20 before any "
                                 "predictor result); MaveDB only, IGVF cross-checked only"),
                "class_name_map": {k: {"class": v, "label": LABEL_OF.get(v)}
                                   for k, v in CLASS_NAME_MAP.items()},
                "column_values": values,
                "unmapped_values": [v for v in values if v not in CLASS_NAME_MAP],
                "range_class": ("the deposit's own word in the column; 'unassigned' for an "
                                "empty value or a value outside the class-name table"),
            }
        prov |= {
            "scope": {"pvs1_pm12": "|offset| 1-2, kept and flagged",
                      "in_scope": "|offset| 3-50", DEEP: "|offset| > 50, kept and flagged"},
            "n_score_rows": info["n_rows"], "n_intronic_snv": info["n_intronic"],
            "n_duplicates_dropped": info["n_duplicates"], "n_rows_written": int(len(df)),
            "deepest_offset": int(df["offset"].abs().max()),
            "sha256": _sha256(out),
            "built_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        }
        Path(str(out) + ".provenance.json").write_text(
            json.dumps(prov, indent=1, ensure_ascii=False, default=str) + "\n")

    s = summary.set_index(["gene", "stratum"])
    print(f"[{stage}] " + ("four external genes, labelled by the deposits' own ranges" if not args.gene
                          else f"{args.gene}, labelled by the deposit's own class column"))
    for gene in genes:
        c = coords.set_index("gene").loc[gene]
        print(f"  {gene:7s} {c.n_mapped:>5,} intron-side SNVs mapped "
              f"(chr{c.chrom} {c.strand}); anchor check {c.n_anchor_agree}/"
              f"{c.n_anchor_checked}; 3-10 nt {s.loc[(gene, 's3_10'), 'n_damaging']} "
              f"damaging, 11-50 nt {s.loc[(gene, 's11_50'), 'n_damaging']}")
    if run["total"]:
        print(f"  distal band (11-50 nt), four genes: "
              f"{s.loc[(run['total'], 's11_50'), 'n_damaging']} damaging variants added")
    print(f"  wrote {REPORT_DIR}/{pre}_*.csv, {EXT_DIR}/<gene>_variants.parquet, {run['doc']}")


if __name__ == "__main__":
    main()
