"""E16 -- four further external genes: records, coordinates and labels.

CTCF, SFPQ, RAD51D and XRCC2 have HAP1 saturation genome editing deposits in MaveDB
(CC0), scored for cell fitness like the seven genes and DDX3X. They are added as
external genes, for a purely prospective test: nothing here is fitted, and no
threshold, score or likelihood ratio is computed at this stage.

What it does, offline, from files scripts/stage_E16/fetch_mavedb.py and
fetch_reference.py put on disk:

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
            (docs/newgenes-record-check.md gives it).
  counts    newgenes_counts_by_offset.csv, newgenes_label_summary.csv.
  checks    newgenes_coordinate_check.csv, newgenes_offset_sample.csv,
            newgenes_controls.csv.

Coordinates. MaveDB holds no GRCh38 postMapped record for any intron-side SNV of
these four score sets (their mapping state is 'incomplete'), so they are made the
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
tests check every reference base against the Ensembl release-112 chromosome.

Run (PYTHONPATH=phase1, from phase1/):
    python -m src.evid_newgenes
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
RAW = {g: EXT_DIR / f"{g.lower()}_raw" for g in GENES}

# the DDX3X pattern, on the RefSeq transcripts these deposits write hgvs_nt against
HGVS_NT = E.intronic_snv_regex(r"N[MR]_[\d.]+")
C_POS = re.compile(r"^c\.([-*]?\d+)([+-]\d+)[ACGT]>[ACGT]$")
EXONIC_C = re.compile(r"^(N[MR]_[\d.]+):c\.([-*]?\d+)[ACGT]>[ACGT]$")
G_HGVS = re.compile(r"^(NC_\d+\.\d+):g\.(\d+)([ACGT])>([ACGT])$")
DEEP = "out_of_scope_deep"
MAX_DEPTH = 50          # the analysis window, 1 <= |offset| <= 50 (evid_build_set)
LABEL_OF = {"abnormal": 1.0, "normal": 0.0}       # anything else: NA, as config.ASSAY_LABELS
UNASSIGNED = "unassigned"
ON_SHARED_POINT = "score on a point two ranges share (author's rule, 2026-10-20)"
N_SAMPLE = 10
SNV_NT = re.compile(r"^N[MR]_[\d.]+:c\.[-*]?\d+([+-]\d+)?[ACGT]>[ACGT]$")
# The author's record of 2026-10-19 for the two URNs given in advance (E16.1): URN,
# licence, publication status (public; no paper or preprint linked to the score set
# or its experiment) and the month of first publication. The stage stops if the
# metadata MaveDB served disagrees with any of them.
_NO_PAPER = {"private": False, "score_set_publications": "not_provided",
             "experiment_publications": "not_provided"}
AUTHOR_RECORD = {
    "CTCF": {"urn": "urn:mavedb:00001262-a-3", "license": "CC0", **_NO_PAPER,
             "first_published": "2026-01"},
    "SFPQ": {"urn": "urn:mavedb:00001265-a-2", "license": "CC0", **_NO_PAPER,
             "first_published": "2026-01"},
}


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


def investigator_ranges(gene: str) -> dict:
    blob = json.loads((RAW[gene] / "score_ranges.json").read_text())
    inv = blob["investigator_provided"]
    if len(inv) != 1:
        raise SystemExit(f"[E16] {gene}: {len(inv)} investigator-provided calibrations; "
                         "exactly one is required (E16.4)")
    return inv[0]


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
    cal = investigator_ranges(gene)
    classes = cal["functionalClassifications"]
    fc = {c["label"]: c["functionalClassification"] for c in classes}
    require_mechanical(gene, classes)

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

    csv_entry = next(e for e in man["files"] if e["path"].endswith(".csv"))
    table = list(csv.DictReader(_here(csv_entry["path"]).open(newline="")))
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
        range_class, why = classify(r["score"], classes)
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
        "shared_points": shared_points(classes), "target": target, "mane": mane,
        "provenance": prov, "strand": strand, "anchor_reference_accessions": anchor_refs,
        "n_rows": len(table), "n_intronic": n_before, "n_duplicates": int(dup.sum()),
        "rna_column": rna, "table": table, "fc": fc,
    }
    return df, info


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
        "target_accession": acc["accession"], "target_assembly": acc["assembly"],
        "n_variants_total": meta["numVariants"],
        "score_ranges_present": bool(info["classes"]),
        "fetched_utc": entry["fetched_utc"], "source_url": entry["source_url"],
        "raw_sha256": entry["sha256"],
        # beyond the fields E16.1 asks for
        "experiment_urn": exp["urn"],
        "experiment_publication_identifiers": _pubs(exp),
        "superseded_versions": "; ".join(man["superseded_versions"]) or "none",
        "mapping_state": meta.get("mappingState"),
        "n_score_rows": info["n_rows"], "scores_csv_sha256": csv_entry["sha256"],
        "calibration_title": info["calibration"]["title"],
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


def counts_tables(frames: dict, infos: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
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
    for gene, df in list(frames.items()) + [("ALL4", allf)]:
        for st in strata:
            sub = K.stratum_frame(df[df["stratum"] != DEEP], st) if st != DEEP \
                else df[df["stratum"] == DEEP]
            summ.append(summarise(gene, st, sub))
    return by_offset, pd.DataFrame(summ)


def coordinate_tables(frames: dict, infos: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
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
        rng = np.random.default_rng([C.RANDOM_SEED, GENES.index(gene)])
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
                cl, _ = classify(score.get(r["urn"]), info["classes"])
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


def write_doc(records: pd.DataFrame, frames: dict, infos: dict) -> None:
    lines = ["# Four further external genes: MaveDB records",
             "",
             "Written by `phase1/src/evid_newgenes.py` (stage E16) from the files "
             "`scripts/stage_E16/fetch_mavedb.py` saved; do not edit by hand. Each fact "
             "is followed by the file it was read from (under "
             "`phase1/data/evidence/external/`) and the time that file was fetched.",
             ""]
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
        classes = "; ".join(f"{c['label']} ({c['functionalClassification']}) {_interval(c)}"
                            for c in info["classes"])
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
                f"- Checked against the author's record of 2026-10-19 (URN, licence, "
                f"publication status, first published 2026-01, the current version dated "
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
                  f"- Investigator-provided ranges, \"{r.calibration_title}\": {classes} "
                  f"{_src(r.gene, man, 'score_ranges.json')}",
                  f"- Shared range ends: {shared}; intron-side SNVs scored there "
                  f"({UNASSIGNED!r} by the author's rule of 2026-10-20): {on_shared} "
                  f"{_src(r.gene, man, 'score_ranges.json', csv_name)}",
                  ""]
    DOC_PATH.write_text("\n".join(lines))


def main() -> None:
    M = _atlas_mapping()
    frames, infos = {}, {}
    for gene in GENES:
        df, info = build_gene(gene, M)
        frames[gene], infos[gene] = df, info

    records = pd.DataFrame([record_row(g, infos[g]) for g in GENES])
    by_offset, summary = counts_tables(frames, infos)
    coords, samples = coordinate_tables(frames, infos)
    controls = controls_table(infos)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    records.to_csv(REPORT_DIR / "newgenes_record_check.csv", index=False)
    by_offset.to_csv(REPORT_DIR / "newgenes_counts_by_offset.csv", index=False)
    summary.to_csv(REPORT_DIR / "newgenes_label_summary.csv", index=False)
    coords.to_csv(REPORT_DIR / "newgenes_coordinate_check.csv", index=False)
    samples.to_csv(REPORT_DIR / "newgenes_offset_sample.csv", index=False)
    controls.to_csv(REPORT_DIR / "newgenes_controls.csv", index=False)
    write_doc(records, frames, infos)

    for gene in GENES:
        df, info = frames[gene], infos[gene]
        out = EXT_DIR / f"{gene.lower()}_variants.parquet"
        df.drop(columns=["rna_score_present"]).to_parquet(out, index=False)
        man = info["man"]
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
            "score_ranges": info["calibration"],
            "range_class": ("the label of the one range the score falls in, in the "
                            "deposit's own wording; 'unassigned' for no score, a score "
                            "outside every range, or a score on a point two ranges share "
                            "(author's rule, 2026-10-20)"),
            "shared_points": info["shared_points"],
            "label_map": {"functionalClassification abnormal": 1,
                          "functionalClassification normal": 0,
                          "any other class, and unassigned": None},
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
    print("[E16] four external genes, labelled by the deposits' own ranges")
    for gene in GENES:
        c = coords.set_index("gene").loc[gene]
        print(f"  {gene:7s} {c.n_mapped:>5,} intron-side SNVs mapped "
              f"(chr{c.chrom} {c.strand}); anchor check {c.n_anchor_agree}/"
              f"{c.n_anchor_checked}; 3-10 nt {s.loc[(gene, 's3_10'), 'n_damaging']} "
              f"damaging, 11-50 nt {s.loc[(gene, 's11_50'), 'n_damaging']}")
    print(f"  distal band (11-50 nt), four genes: "
          f"{s.loc[('ALL4', 's11_50'), 'n_damaging']} damaging variants added")
    print(f"  wrote {REPORT_DIR}/newgenes_*.csv, {EXT_DIR}/<gene>_variants.parquet, {DOC_PATH}")


if __name__ == "__main__":
    main()
