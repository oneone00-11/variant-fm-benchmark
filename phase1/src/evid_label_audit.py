"""E16c -- the seven training genes' label files against the sources they cite: a
read-only audit.

Each label file in data/assay_labels/ (config.ASSAY_LABELS) is compared, variant by
variant, with what the score set it cites (supplementary/Table_S1.tsv) gives under
each candidate rule of E16c.2:

  R1  each investigator-provided calibration of the cited version: a score placed in
      its ranges (point estimate), or, where its classes are classes rather than
      ranges, MaveDB's own per-variant membership of them;
  R2  each per-variant class column of the cited version's score table;
  R3  the investigator-provided ranges with the score table's 95% confidence
      interval: a class only when the whole interval lies in one range, intermediate
      when it crosses a bound (only where the table has the interval);
  R4  the rule the cited paper states, where it states one a program can apply
      (label_audit/<gene>_raw/paper_rule.json: the rule quoted verbatim, with its
      place in the paper);
  R5  R1 and R2 on every other version of the score set's chain;
  R6  a paper supplementary table's per-variant classes, where the author made one a
      label source (BRCA2's Supplementary Table 3, 2026-10-21), grouped as the paper's
      own text groups them (paper_rule.json, quoted).

Every rule gives a row, whether it reproduces the file or not. Nothing is fitted: a
rule is a field of the archive or a sentence of the paper, never a cut-off found by
trying.

Classes are compared as abnormal / normal / intermediate. The label file's class goes
through config.ASSAY_LABELS (1 / 0 / None); a calibration's through MaveDB's own
functionalClassification (abnormal / normal / not_specified); a class column's through
the functionalClassification an investigator calibration of the same version gives
the same class, else through COLUMN_WORDS, a wording table printed in the record page.
A variant the source leaves without a class ('unassigned') differs. Variants are
aligned on the normalised c. HGVS (evid_build_set._norm_hgvs). A version written only
in genomic notation is aligned on GRCh38 (chrom, pos, ref, alt) instead, the label
file's c. HGVS mapped through the atlas mapping module as E16 maps it (MANE transcript
model, tracked Ensembl slice, reference base checked).

Three large calibration memberships and the BRCA2 table are not committed (the
author's decision of 2026-10-21): the audit reads the compact extracts committed in
their place (evid_label_sources), which scripts/reproduce_evidence.py --fetch-inputs
checks against the files as served.

Writes reports/evidence/label_audit_summary.csv, label_audit_diff.csv,
label_audit_impact.csv, label_audit_unclassified.csv and docs/label-source-audit.md,
and nothing else: the label files, the configuration and the analysis set are only
read.

Run (PYTHONPATH=phase1, from phase1/):
    python -m src.evid_label_audit
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd

from . import config as C
from . import evid_common as K
from . import evid_label_sources as LS
from . import evid_newgenes as N
from .evid_build_set import _norm_hgvs

AUDIT_DIR = Path("data/evidence/label_audit")
TABLE_S1 = Path("../supplementary/Table_S1.tsv")
REPORT_DIR = K.REPORT_DIR
OUT_SUMMARY = REPORT_DIR / "label_audit_summary.csv"
OUT_DIFF = REPORT_DIR / "label_audit_diff.csv"
OUT_IMPACT = REPORT_DIR / "label_audit_impact.csv"
OUT_UNCLASSIFIED = REPORT_DIR / "label_audit_unclassified.csv"
DOC_PATH = Path("../docs/label-source-audit.md")
DEPTH_COUNTS = REPORT_DIR / "depth_counts.csv"
SET_COUNTS = REPORT_DIR / "set_counts.csv"
# The seven label files entered the repository together, in commit 10d94bb on
# 2026-07-08 (git log of phase1/data/assay_labels/; a test checks it where git is
# here). A calibration modified after that day would mean the archive changed its
# classes since (E16c.7).
LABEL_FILES_COMMITTED = "2026-07-08"
CLASS_COLUMNS = ("functional_consequence", "functional_classification", "class",
                 "coarse_class", "func.class")
CI_LOWER, CI_UPPER = "95_ci_lower", "95_ci_upper"
ABN, NOR, INT, UNA = "abnormal", "normal", "intermediate", "unassigned"
FC_CANON = {"abnormal": ABN, "normal": NOR, "not_specified": INT}
# A class column's words where no investigator calibration of the same version names
# the class: wording only, as E16b maps SBDS's column (evid_newgenes.CLASS_NAME_MAP)
COLUMN_WORDS = {k: FC_CANON[v] for k, v in N.CLASS_NAME_MAP.items()}
WINDOW = (3, N.MAX_DEPTH)          # the in-scope window, 3 <= |offset| <= 50
STOP_PCT = 1.0                     # E16c.7: no rule reproduces and every one differs by more
STRATA = ["pm12", "s3_10", "s11_50", "s3_50"]
DEPTH_BIN = {"s3_10": "3-10", "s11_50": "11-50", "s3_50": "3-50"}
OFFSET = re.compile(r"c\.[-*]?\d+([+-]\d+)")
UNCLASSIFIED_STRATA = ["pm12", "s3_10", "s11_50", "deep", "exonic_or_utr"]
# The author's decisions of 2026-10-21 on what E16c found (E16d), printed with the
# findings they answer
AUTHOR_DECISIONS = {
    "BRCA2": "the paper's Supplementary Table 3 is BRCA2's label source (rule R6)",
    "BARD1": ("relabel by R1, the cited a-2's investigator-provided ranges, and refit "
              "(stage E16d)"),
    "PALB2": ("relabel by R1, the cited a-2's investigator-provided ranges, and refit "
              "(stage E16d)"),
    "VHL": ("the score-table variants in no class of the investigator calibration stay "
            "unclassified, as in the label file, and the Methods say so"),
    "BAP1": ("the sub-class names strongly / weakly depleted are accepted: they split the "
             "depleted variants at the median the paper's Figure 4e legend prints, "
             "-0.1260642; c.1132del scores exactly that median, the label file calls it "
             "strongly depleted and the paper does not say on which side it falls; both "
             "sub-classes are damaging, so the analysis does not depend on it"),
}


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _fold(name: str) -> str:
    """A class name with case and separators folded, for the class-name check only."""
    return " ".join(str(name).lower().replace("_", " ").replace("-", " ").split())


def _offset(hgvs: str) -> int | None:
    m = OFFSET.search(str(hgvs).split(":")[-1])
    return int(m.group(1)) if m else None


def _depth_stratum(offset: int | None) -> str:
    if offset is None:
        return "exonic_or_utr"
    d = abs(offset)
    return "deep" if d > N.MAX_DEPTH else N.E.stratum_of(d)


# ---------------------------------------------------------------------------
# inputs
# ---------------------------------------------------------------------------
def table_s1() -> dict:
    with open(TABLE_S1, newline="") as f:
        return {r["gene"]: r for r in csv.DictReader(f, delimiter="\t")}


def cited_paper(s1_row: dict) -> str | None:
    """The paper a Table S1 row cites, as 'Author et al. year' (no venue), or None when
    the row cites the MaveDB score set alone."""
    ref = s1_row["reference_pmid_source"]
    if ref.startswith("MaveDB score set"):
        return None
    m = re.match(r"^(.+? et al\. \d{4})", ref)
    return m.group(1) if m else ref


def our_labels(gene: str) -> dict:
    """The label file as config.ASSAY_LABELS reads it: one class per normalised HGVS."""
    spec = C.ASSAY_LABELS[gene]
    path = C.ASSAY_LABEL_DIR / spec["file"]
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t" if spec["file"].endswith(".tsv") else ","))
    variants, dup = {}, 0
    for r in rows:
        key = _norm_hgvs(r[spec["key"]])
        cls = r[spec["class_col"]]
        canon = ({1.0: ABN, 0.0: NOR, None: INT}[spec["map"][cls]] if cls in spec["map"]
                 else UNA)
        if key in variants:
            dup += 1
            continue
        variants[key] = {"hgvs": r[spec["key"]], "cls": cls, "canon": canon}
    return {"path": path, "sha256": _sha256(path), "rows": len(rows), "duplicates": dup,
            "variants": variants, "classes": sorted({v["cls"] for v in variants.values()})}


def load_manifest(gene: str) -> dict:
    """The audit's raw manifest, every file it lists checked against its sha256."""
    path = AUDIT_DIR / f"{gene.lower()}_raw" / "manifest.json"
    man = json.loads(path.read_text())
    for e in man["files"]:
        if e.get("cached"):        # not committed; --fetch-inputs and the tests check it
            continue
        f = N._here(e["path"])
        if not f.is_file() or _sha256(f) != e["sha256"]:
            raise SystemExit(f"[E16c] {e['path']} is missing or differs from {path}; rerun "
                             "scripts/stage_E16c/fetch_sources.py")
    return man


def _file(man: dict, name: str) -> Path | None:
    hit = [N._here(e["path"]) for e in man["files"] if Path(e["path"]).name == name]
    return hit[0] if hit else None


def version_data(man: dict, urn: str) -> dict:
    vid = urn.split(":")[-1]
    meta_path = _file(man, f"scoreset_meta_{vid}.json")
    table_path = _file(man, f"urn_mavedb_{vid}.csv")
    with open(table_path, newline="") as f:
        table = list(csv.DictReader(f))
    meta = json.loads(meta_path.read_text())
    members = {e["calibration_urn"]: N._here(e["path"]) for e in man["files"]
               if e.get("score_set") == urn and "compact_of" in e}
    genomic = bool(table) and all((r["hgvs_nt"] or "").startswith("NC_") for r in table[:100])
    return {"urn": urn, "meta": meta, "table": table, "table_path": table_path,
            "table_sha256": _sha256(table_path), "members": members, "genomic": genomic,
            "calibrations": [c for c in meta.get("scoreCalibrations") or []
                             if c.get("investigatorProvided") is True]}


def _gkey(chrom: str, pos, ref: str, alt: str) -> tuple:
    return (str(chrom).removeprefix("chr"), int(pos), ref, alt)


def genomic_index(gene: str, ours: dict, M) -> tuple[dict, dict]:
    """The label file's variants on GRCh38 (chrom, pos, ref, alt), mapped as E16 maps:
    the MANE Select transcript model of the file's transcript, the tracked Ensembl
    slice, the reference base checked. Returns (genomic key -> our key, failures)."""
    txs = {v["hgvs"].split(":")[0] for v in ours["variants"].values() if ":" in v["hgvs"]}
    if len(txs) != 1:
        raise SystemExit(f"[E16c] {gene}: the label file names {sorted(txs)}, not one transcript")
    mane = M.load_mane_records([gene])[gene]
    tmap, _prov = M.build_transcript_map(txs.pop(), mane)
    lo = min(e.g_start for e in tmap.exons) - M.REGION_PAD
    hi = max(e.g_end for e in tmap.exons) + M.REGION_PAD
    s0, _, seq = M.fetch_ensembl_region(mane.chr_name, lo, hi)
    rseq = M.RegionSeq(s0, seq)
    index, failures = {}, {}
    for key, v in ours["variants"].items():
        try:
            pos, ref, alt = M.map_variant(tmap, rseq, v["hgvs"].split(":")[-1])
        except M.MappingFailure as e:
            failures[v["hgvs"]] = e.category
            continue
        index[_gkey(mane.chr_name, pos, ref, alt)] = key
    return index, failures


def analysis_window(aset: pd.DataFrame, gene: str) -> dict:
    """Normalised HGVS -> (|offset|, stratum) of the gene's analysis-set rows."""
    v = aset[aset["gene"] == gene]
    return {_norm_hgvs(h): (int(o), s)
            for h, o, s in zip(v["hgvs_c"], v["intron_offset_abs"], v["stratum"])}


# ---------------------------------------------------------------------------
# the rules: each gives {source key -> (class as written, canonical class, score)}
# ---------------------------------------------------------------------------
def _ranges(cal: dict) -> bool:
    fcs = cal.get("functionalClassifications") or []
    return bool(fcs) and all(fc.get("range") is not None for fc in fcs)


def _row_key(r: dict, genomic: bool):
    if not genomic:
        return _norm_hgvs(r["hgvs_nt"])
    if not str(r.get("pos") or "").isdigit():
        return ("no position", r["accession"])
    return _gkey(r["chrom"], r["pos"], r["ref"], r["alt"])


def rule_ranges(vd: dict, cal: dict) -> dict:
    classes = cal["functionalClassifications"]
    fc = {c["label"]: FC_CANON.get(c["functionalClassification"], UNA) for c in classes}
    out = {}
    for r in vd["table"]:
        label, why = N.classify(r["score"], classes)
        out[_row_key(r, vd["genomic"])] = (label if label != N.UNASSIGNED else f"{UNA} ({why})",
                                           fc.get(label, UNA), r["score"])
    return out


def rule_members(vd: dict, cal: dict) -> dict:
    """A class-based calibration: each variant's class is the class whose membership
    MaveDB serves it in; a score-table row in no class is unassigned."""
    by_id = {str(c["id"]): c for c in cal["functionalClassifications"]}
    member = {r["variant_urn"]: by_id[r["functional_classification_id"]]
              for r in LS.read_tsv(vd["members"][cal["urn"]])}
    out = {}
    for r in vd["table"]:
        c = member.get(r["accession"])
        out[_row_key(r, vd["genomic"])] = (
            (c["label"], FC_CANON.get(c["functionalClassification"], UNA), r["score"]) if c
            else (f"{UNA} (in no class of the calibration)", UNA, r["score"]))
    return out


def class_words(vd: dict) -> dict:
    """A class-column word -> (canonical class, where the mapping comes from)."""
    words = {k: (v, "wording table COLUMN_WORDS") for k, v in COLUMN_WORDS.items()}
    for cal in vd["calibrations"]:                      # the first calibration that names
        for c in cal["functionalClassifications"]:      # a class gives its mapping
            have = words.get(c.get("class"))
            if c.get("class") is not None and (have is None or have[1].startswith("wording")):
                words[c["class"]] = (FC_CANON.get(c["functionalClassification"], UNA),
                                     f"calibration '{cal.get('title')}'")
    return words


def rule_column(vd: dict, col: str) -> tuple[dict, dict]:
    words = class_words(vd)
    out, used = {}, {}
    for r in vd["table"]:
        w = r.get(col) or ""
        canon, src = words.get(w, (UNA, "not in the table"))
        used[w] = (canon, src)
        out[_row_key(r, vd["genomic"])] = (w if canon != UNA else f"{UNA} ({w or 'empty'})",
                                           canon, r["score"])
    return out, used


def rule_interval(vd: dict, cal: dict) -> dict:
    """A class only when the whole 95% interval lies in one range; an interval across
    a bound is intermediate; a row without an interval is unassigned."""
    classes = cal["functionalClassifications"]
    fc = {c["label"]: FC_CANON.get(c["functionalClassification"], UNA) for c in classes}
    out = {}
    for r in vd["table"]:
        lo, hi = r.get(CI_LOWER), r.get(CI_UPPER)
        key = _row_key(r, vd["genomic"])
        if lo in (None, "", "NA") or hi in (None, "", "NA"):
            out[key] = (f"{UNA} (no interval)", UNA, r["score"])
            continue
        both = [c["label"] for c in classes
                if N._matches(float(lo), c) and N._matches(float(hi), c)]
        out[key] = ((both[0], fc[both[0]], r["score"]) if len(both) == 1
                    else ("intermediate (interval crosses a bound)", INT, r["score"]))
    return out


def rule_paper(vd: dict, rule: dict) -> dict:
    """The paper's rule as paper_rule.json writes it: ordered classes, each a list of
    conditions on named score-table columns, the first class whose conditions all hold
    wins, else the default class."""
    ops = {"<": lambda a, b: a < b, "<=": lambda a, b: a <= b,
           ">": lambda a, b: a > b, ">=": lambda a, b: a >= b}
    out = {}
    for r in vd["table"]:
        key = _row_key(r, vd["genomic"])
        hit = None
        try:
            for c in rule["classes"]:
                if all(ops[cond["op"]](float(r[cond["column"]]), float(cond["value"]))
                       for cond in c["conditions"]):
                    hit = c
                    break
        except (KeyError, ValueError):
            out[key] = (f"{UNA} (a column the rule reads is empty)", UNA, r["score"])
            continue
        hit = hit or rule["default"]
        out[key] = (hit["class"], hit["canonical"], r["score"])
    return out


def supplementary_data(entry: dict) -> dict:
    """A supplementary table's compact extract, as a source with the fields run() reads."""
    spec = entry["supplementary"]
    return {"urn": f"{spec['pmcid']} {spec['member']} sheet '{spec['sheet']}'",
            "table": LS.read_tsv(N._here(entry["path"])), "genomic": False}


def rule_supplementary(sd: dict, rule: dict) -> dict:
    """R6: the table's own category of each variant, grouped as rule['groups'] (the
    paper's text, quoted in paper_rule.json); a category outside the groups is
    unassigned."""
    return {_norm_hgvs(r["hgvs_nt"]): (r["category"], rule["groups"].get(r["category"], UNA),
                                       r["score"]) for r in sd["table"]}


# ---------------------------------------------------------------------------
# comparison
# ---------------------------------------------------------------------------
def compare(gene: str, ours: dict, source: dict, window: dict, row: dict,
            genomic_to_ours: dict | None = None) -> tuple[dict, list[dict]]:
    """One rule's source classes against the label file, on the variants both have."""
    n_source = len(source)
    if genomic_to_ours is not None:
        source = {genomic_to_ours[k]: v for k, v in source.items() if k in genomic_to_ours}
    both = [k for k in ours["variants"] if k in source]
    diffs, n_same, n_written, n_win, n_win_diff = [], 0, 0, 0, 0
    for k in both:
        o, (s_cls, s_canon, s_score) = ours["variants"][k], source[k]
        off = _offset(o["hgvs"])
        in_win = k in window and WINDOW[0] <= window[k][0] <= WINDOW[1]
        n_win += in_win
        n_written += _fold(o["cls"]) == _fold(s_cls)
        if o["canon"] == s_canon:
            n_same += 1
            continue
        n_win_diff += in_win
        diffs.append({"gene": gene, "rule": row["rule"], "version": row["version"],
                      "rule_detail": row["rule_detail"], "hgvs_nt": o["hgvs"],
                      "offset": off, "stratum": window[k][1] if k in window else _depth_stratum(off),
                      "our_class": o["cls"], "source_class": s_cls, "source_score": s_score,
                      "in_analysis_set": in_win, "our_canonical": o["canon"],
                      "source_canonical": s_canon})
    n = len(both)
    row |= {"n_compared": n, "n_identical": n_same, "n_differ": n - n_same,
            "pct_differ_all": 100 * (n - n_same) / n if n else np.nan,
            "n_compared_3_50nt": n_win, "n_differ_3_50nt": n_win_diff,
            "pct_differ_3_50nt": 100 * n_win_diff / n_win if n_win else np.nan,
            "exact_match": ("yes" if n and n == n_same else "no") if n else "not applicable",
            "n_only_ours": len(ours["variants"]) - n,
            "n_only_source": n_source - n,
            "n_identical_as_written": n_written}
    return row, diffs


def _na_row(base: dict, rule: str, version: str, detail: str, sha: str = "") -> dict:
    return base | {"rule": rule, "version": version, "rule_detail": detail,
                   "n_source": np.nan, "n_compared": 0, "n_identical": 0, "n_differ": 0,
                   "pct_differ_all": np.nan, "n_compared_3_50nt": 0, "n_differ_3_50nt": 0,
                   "pct_differ_3_50nt": np.nan, "exact_match": "not applicable",
                   "n_only_ours": np.nan, "n_only_source": np.nan,
                   "n_identical_as_written": np.nan, "source_sha256": sha}


def audit_gene(gene: str, s1: dict, aset: pd.DataFrame, M) -> dict:
    ours = our_labels(gene)
    man = load_manifest(gene)
    cited = s1[gene]["mavedb_urn"]
    if man["cited"]["urn"] != cited:
        raise SystemExit(f"[E16c] {gene}: the manifest cites {man['cited']['urn']}, "
                         f"Table S1 {cited}")
    window = analysis_window(aset, gene)
    base = {"gene": gene, "label_file": C.ASSAY_LABELS[gene]["file"], "urn": cited,
            "n_our_file": ours["rows"]}
    served = [v["urn"] for v in man["versions"]]
    data = {u: version_data(man, u) for u in served}
    gindex, gfail = None, {}
    if any(d["genomic"] for d in data.values()):
        gindex, gfail = genomic_index(gene, ours, M)
    rows, diffs, notes = [], [], {"column_words": {}, "genomic_mapping_failures": gfail}

    def run(rule: str, vd: dict, detail: str, source: dict, sha: str) -> None:
        if vd["genomic"]:
            detail += (f"; aligned on GRCh38 (the version writes only g. HGVS; "
                       f"{len(gfail)} label-file variants unmapped)")
        row = base | {"rule": rule, "version": vd["urn"], "rule_detail": detail,
                      "n_source": len(vd["table"]), "source_sha256": sha}
        row, d = compare(gene, ours, source, window, row, gindex if vd["genomic"] else None)
        rows.append(row)
        diffs.extend(d)

    def r1_r2(vd: dict, rule: str, prefix: str = "") -> int:
        n = 0
        for cal in vd["calibrations"]:
            if _ranges(cal):
                run(rule, vd, f"{prefix}ranges of calibration '{cal.get('title')}' "
                    f"({len(cal['functionalClassifications'])} classes), point estimate",
                    rule_ranges(vd, cal), vd["table_sha256"])
            elif cal["urn"] in vd["members"]:
                run(rule, vd, f"{prefix}membership of class-based calibration "
                    f"'{cal.get('title')}' ({len(cal['functionalClassifications'])} classes) "
                    f"served by MaveDB", rule_members(vd, cal),
                    _sha256(vd["members"][cal["urn"]]))
            else:
                rows.append(_na_row(base, rule, vd["urn"], f"{prefix}calibration "
                                    f"'{cal.get('title')}' has neither ranges nor a "
                                    "served membership"))
            n += 1
        if not vd["calibrations"]:
            rows.append(_na_row(base, rule, vd["urn"], f"{prefix}no investigator-provided "
                                "calibration"))
        return n

    def r2(vd: dict, rule: str, prefix: str = "") -> None:
        cols = [c for c in CLASS_COLUMNS if vd["table"] and c in vd["table"][0]]
        for col in cols:
            src, used = rule_column(vd, col)
            notes["column_words"][f"{vd['urn']} {col}"] = used
            run(rule, vd, f"{prefix}score-table column '{col}'", src, vd["table_sha256"])
        if not cols:
            rows.append(_na_row(base, rule, vd["urn"], f"{prefix}no class column in the "
                                "score table"))

    vd = data[cited]
    r1_r2(vd, "R1")
    r2(vd, "R2")
    has_ci = bool(vd["table"]) and CI_LOWER in vd["table"][0] and CI_UPPER in vd["table"][0]
    range_cals = [c for c in vd["calibrations"] if _ranges(c)]
    if has_ci and range_cals:
        for cal in range_cals:
            run("R3", vd, f"ranges of calibration '{cal.get('title')}' with the score "
                f"table's {CI_LOWER} / {CI_UPPER}", rule_interval(vd, cal), vd["table_sha256"])
    else:
        rows.append(_na_row(base, "R3", cited, "no 95% interval columns in the score table"
                            if not has_ci else "no investigator-provided ranges"))
    paper = AUDIT_DIR / f"{gene.lower()}_raw" / "paper_rule.json"
    if cited_paper(s1[gene]) is None:
        rows.append(_na_row(base, "R4", cited, "Table S1 cites no paper: not applicable"))
    elif not paper.is_file():
        raise SystemExit(f"[E16c] {gene}: Table S1 cites {cited_paper(s1[gene])}, and "
                         f"{paper} is missing")
    else:
        pr = json.loads(paper.read_text())
        check = json.loads((paper.parent / "paper_quote_check.json").read_text())
        notes["paper_rule"], notes["paper_quote_check"] = pr, check
        ids = pr["source"]["citation_ids"]
        sha = _sha256(paper)
        if check["paper_rule_sha256"] != sha or not all(q["found"] for q in check["quotes"]):
            raise SystemExit(f"[E16c] {gene}: {paper} is not the file the quote check read, or "
                             "a quote was not found in the paper; rerun "
                             "scripts/stage_E16c/check_paper_quotes.py")
        if pr["status"] != "stated":
            rows.append(_na_row(base, "R4", cited, f"paper rule {pr['status']} ({ids}): "
                                f"{pr['not_applicable_because']}", sha))
        for rule in pr["rules"] if pr["status"] == "stated" else []:
            need = sorted({cond["column"] for c in rule["classes"] for cond in c["conditions"]})
            lack = [c for c in need if c not in vd["table"][0]]
            detail = f"paper rule '{rule['name']}' ({ids}; {rule['where']})"
            if lack:
                rows.append(_na_row(base, "R4", cited, f"{detail} reads columns the score table "
                                    f"does not have: {', '.join(lack)}", sha))
            else:
                run("R4", vd, detail, rule_paper(vd, rule), vd["table_sha256"] + ";" + sha)
    sup = next((e for e in man["files"] if "supplementary" in e), None)
    if sup is None:
        rows.append(_na_row(base, "R6", cited, "a paper supplementary table is a candidate "
                            "for BRCA2 only (the author's decision of 2026-10-21)"))
    else:
        rule6 = notes["paper_rule"]["supplementary_rule"]
        sd = supplementary_data(sup)
        run("R6", sd, f"supplementary table: {rule6['name']} ({rule6['quote_where']})",
            rule_supplementary(sd, rule6), sup["sha256"] + ";" + _sha256(paper))
    others = [u for u in served if u != cited]
    for u in others:
        r1_r2(data[u], "R5", "R1 on another version: ")
        r2(data[u], "R5", "R2 on another version: ")
    for miss in man.get("versions_not_served") or []:
        rows.append(_na_row(base, "R5", miss["urn"], f"version not served by MaveDB "
                            f"(HTTP {miss['http_status']})"))
    if not others and not man.get("versions_not_served"):
        rows.append(_na_row(base, "R5", cited, "no other version in the chain"))
    return {"gene": gene, "ours": ours, "manifest": man, "data": data, "rows": rows,
            "diffs": diffs, "notes": notes, "window": window, "gindex": gindex}


# ---------------------------------------------------------------------------
# impact, class names, flags
# ---------------------------------------------------------------------------
def current_counts(gene: str) -> dict:
    """(stratum) -> counts now: pm12 from set_counts.csv (summed over the arms), the
    in-scope strata from depth_counts.csv."""
    out = {}
    with open(DEPTH_COUNTS, newline="") as f:
        for r in csv.DictReader(f):
            for st, b in DEPTH_BIN.items():
                if r["gene"] == gene and r["depth_nt"] == b:
                    out[st] = {"n": int(r["n"]), "labelled": int(r["labelled"]),
                               "damaging": int(r["damaging"]), "normal": int(r["normal"])}
    with open(SET_COUNTS, newline="") as f:
        pm = [r for r in csv.DictReader(f) if r["gene"] == gene and r["stratum"] == "pm12"]
    out["pm12"] = {"n": sum(int(r["n"]) for r in pm),
                   "labelled": sum(int(r["n_labelled"]) for r in pm),
                   "damaging": sum(int(r["n_damaging"]) for r in pm),
                   "normal": sum(int(r["n_normal"]) for r in pm)}
    return out


def closest(rows: list[dict]) -> dict | None:
    live = [r for r in rows if r["exact_match"] == "no"]
    if not live:
        return None
    return min(live, key=lambda r: (r["pct_differ_3_50nt"] if r["pct_differ_3_50nt"] == r["pct_differ_3_50nt"]
                                    else math.inf, r["pct_differ_all"]))


def impact(res: dict, aset: pd.DataFrame) -> list[dict]:
    gene = res["gene"]
    now = current_counts(gene)
    v = aset[aset["gene"] == gene]
    for st in STRATA:                                 # the counts files against the set
        sub = K.stratum_frame(v, st)
        lab = sub["y_assay"]
        got = {"n": len(sub), "labelled": int(lab.notna().sum()),
               "damaging": int((lab == 1).sum()), "normal": int((lab == 0).sum())}
        if got != now[st]:
            raise SystemExit(f"[E16c] {gene} {st}: the counts files say {now[st]}, the "
                             f"analysis set {got}")
    best = closest(res["rows"])
    out = []
    for st in STRATA:
        sub = K.stratum_frame(v, st)
        row = {"gene": gene, "stratum": st, "closest_rule": best["rule"] if best else "none",
               "version": best["version"] if best else "",
               "rule_detail": best["rule_detail"] if best else "no archive rule applies",
               "pct_differ_3_50nt": best["pct_differ_3_50nt"] if best else np.nan,
               "n_variants": now[st]["n"],
               "n_labelled_now": now[st]["labelled"], "n_damaging_now": now[st]["damaging"],
               "n_normal_now": now[st]["normal"],
               "n_excluded_now": now[st]["n"] - now[st]["labelled"]}
        if best is None:
            row |= {k: np.nan for k in ("n_labelled_rule", "n_damaging_rule", "n_normal_rule",
                                        "n_excluded_rule")}
        else:
            src = source_of(res, best)
            canon = []
            for h, chrom, pos, ref, alt in zip(sub["hgvs_c"], sub["chrom"], sub["pos"],
                                               sub["ref"], sub["alt"]):
                key = (_gkey(chrom, pos, ref, alt) if res["data"][best["version"]]["genomic"]
                       else _norm_hgvs(h))
                canon.append(src[key][1] if key in src else UNA)
            dmg, nor = canon.count(ABN), canon.count(NOR)
            row |= {"n_labelled_rule": dmg + nor, "n_damaging_rule": dmg, "n_normal_rule": nor,
                    "n_excluded_rule": len(canon) - dmg - nor}
        out.append(row)
    return out


def source_of(res: dict, row: dict) -> dict:
    """The source classes behind one summary row, recomputed (for the impact counts)."""
    detail = row["rule_detail"]
    if detail.startswith("supplementary table"):
        sup = next(e for e in res["manifest"]["files"] if "supplementary" in e)
        return rule_supplementary(supplementary_data(sup),
                                  res["notes"]["paper_rule"]["supplementary_rule"])
    vd = res["data"][row["version"]]
    for cal in vd["calibrations"]:
        title = f"'{cal.get('title')}'"
        if title in detail:
            if "95_ci" in detail:
                return rule_interval(vd, cal)
            return rule_ranges(vd, cal) if _ranges(cal) else rule_members(vd, cal)
    for col in CLASS_COLUMNS:
        if f"column '{col}'" in detail:
            return rule_column(vd, col)[0]
    if detail.startswith("paper rule"):
        for rule in res["notes"]["paper_rule"]["rules"]:
            if f"'{rule['name']}'" in detail:
                return rule_paper(vd, rule)
    raise SystemExit(f"[E16c] no source for {row['rule']} {detail}")


def unclassified(res: dict) -> list[dict]:
    """The cited score table's variants that the label file does not carry, by stratum
    (the author asked for VHL's on 2026-10-21; any gene with such variants is listed)."""
    vd = res["data"][res["manifest"]["cited"]["urn"]]
    missing = [r for r in vd["table"] if _row_key(r, vd["genomic"]) not in res["ours"]["variants"]]
    if not missing:
        return []
    n, n_set = {}, {}
    for r in missing:
        k = _norm_hgvs(r["hgvs_nt"])
        st = res["window"][k][1] if k in res["window"] else _depth_stratum(_offset(r["hgvs_nt"]))
        n[st] = n.get(st, 0) + 1
        n_set[st] = n_set.get(st, 0) + (k in res["window"])
    return [{"gene": res["gene"], "version": vd["urn"], "stratum": st, "n": n.get(st, 0),
             "n_in_analysis_set": n_set.get(st, 0)}
            for st in UNCLASSIFIED_STRATA + sorted(set(n) - set(UNCLASSIFIED_STRATA))]


def archive_names(res: dict) -> dict:
    """Every class name the archive writes for the gene, in any version: calibration
    labels and classes, functionalClassification values, class-column values."""
    names = {}
    for u, vd in res["data"].items():
        for cal in vd["meta"].get("scoreCalibrations") or []:
            for c in cal.get("functionalClassifications") or []:
                for field in ("label", "class", "functionalClassification"):
                    if c.get(field):
                        names.setdefault(c[field], set()).add(f"{u} calibration "
                                                              f"'{cal.get('title')}' {field}")
        for col in CLASS_COLUMNS:
            if vd["table"] and col in vd["table"][0]:
                for r in vd["table"]:
                    if r[col]:
                        names.setdefault(r[col], set()).add(f"{u} column {col}")
    return names


def class_name_check(res: dict) -> list[dict]:
    names = archive_names(res)
    folded = {}
    for n, where in names.items():
        folded.setdefault(_fold(n), set()).update(where)
    out = []
    for cls in res["ours"]["classes"]:
        out.append({"class": cls, "as_written": cls in names,
                    "case_and_separators_folded": _fold(cls) in folded,
                    "where": sorted(names.get(cls) or folded.get(_fold(cls)) or [])[:4]})
    return out


def flags(res: dict, s1: dict) -> list[str]:
    """The conditions of E16c.7, each with its evidence."""
    gene, man, out = res["gene"], res["manifest"], []
    cited = next(v for v in man["versions"] if v["role"] == "cited")
    if cited["num_variants_metadata"] != res["ours"]["rows"] or \
            cited["n_rows_score_table"] != res["ours"]["rows"]:
        out.append(f"variant count: the cited {cited['urn']} has {cited['num_variants_metadata']:,} "
                   f"variants in its metadata and {cited['n_rows_score_table']:,} score-table rows; "
                   f"the label file has {res['ours']['rows']:,} rows")
    if not any(r["exact_match"] == "yes" for r in res["rows"]):
        best = closest(res["rows"])
        if best is None:
            out.append("no rule reproduces the file, and no rule applies at all")
        elif not best["pct_differ_3_50nt"] <= STOP_PCT:
            out.append(f"no rule reproduces the file; the closest ({best['rule']}, "
                       f"{best['rule_detail']}) differs on {best['pct_differ_3_50nt']:.2f}% at 3-50 nt")
    for c in class_name_check(res):
        if not c["case_and_separators_folded"]:
            out.append(f"class name '{c['class']}' appears nowhere in the archive's records "
                       "of this score set (any version, any calibration or class column)")
    for cal in cited["score_calibrations"]:
        if cal["investigator_provided"] is True and cal["modified"] != "not_provided" and \
                cal["modified"] > LABEL_FILES_COMMITTED:
            out.append(f"calibration '{cal['title']}' of the cited version was modified "
                       f"{cal['modified']}, after the label files were committed "
                       f"({LABEL_FILES_COMMITTED})")
    lic = cited["license"]["shortName"]
    if lic != s1[gene]["license"]:
        out.append(f"licence: the archive says {lic}, Table S1 {s1[gene]['license']}")
    if cited.get("private") is True:
        out.append("the cited version is private")
    return out


# ---------------------------------------------------------------------------
# outputs
# ---------------------------------------------------------------------------
SUMMARY_COLUMNS = ["gene", "label_file", "urn", "version", "rule", "rule_detail",
                   "n_our_file", "n_source", "n_compared", "n_identical", "n_differ",
                   "pct_differ_all", "pct_differ_3_50nt", "exact_match", "source_sha256",
                   "n_compared_3_50nt", "n_differ_3_50nt", "n_only_ours", "n_only_source",
                   "n_identical_as_written"]
DIFF_COLUMNS = ["gene", "rule", "version", "rule_detail", "hgvs_nt", "offset", "stratum",
                "our_class", "source_class", "source_score", "in_analysis_set",
                "our_canonical", "source_canonical"]
IMPACT_COLUMNS = ["gene", "stratum", "closest_rule", "version", "rule_detail",
                  "pct_differ_3_50nt", "n_variants", "n_labelled_now", "n_labelled_rule",
                  "n_damaging_now", "n_damaging_rule", "n_normal_now", "n_normal_rule",
                  "n_excluded_now", "n_excluded_rule"]


def _fmt_pct(x) -> str:
    return "n/a" if x != x else f"{x:.2f}%"


def write_doc(results: list[dict], summary: pd.DataFrame, impact_df: pd.DataFrame,
              s1: dict, all_flags: dict) -> None:
    L = ["# Label source audit: the seven training genes",
         "",
         "Written by `phase1/src/evid_label_audit.py` (stage E16c) from the files "
         "`scripts/stage_E16c/fetch_sources.py` saved under `phase1/data/evidence/label_audit/`; "
         "do not edit by hand. A read-only audit: no label file, configuration entry, "
         "threshold or analysis-set row is changed. Each label file in "
         "`phase1/data/assay_labels/` is compared, variant by variant, with what the score set "
         "Supplementary Table S1 cites gives under each candidate rule (R1 investigator-provided "
         "calibration; R2 score-table class column; R3 ranges with the 95% interval; R4 the "
         "cited paper's stated rule; R5 R1 and R2 on the other versions of the chain; R6 a "
         "paper supplementary table the author made a source, for BRCA2 only). Classes "
         "are compared as abnormal / normal / intermediate; a variant the source leaves without "
         "a class counts as different. Where no rule reproduces a file, its source is written "
         "as not established; no reason is offered.",
         "",
         "The full rows are in `phase1/reports/evidence/label_audit_summary.csv`, every "
         "differing variant in `label_audit_diff.csv`, and the counts each closest rule would "
         "give the analysis set in `label_audit_impact.csv`.",
         ""]
    for res in results:
        gene, man, ours = res["gene"], res["manifest"], res["ours"]
        rows = summary[summary["gene"] == gene]
        exact = rows[rows["exact_match"] == "yes"]
        L += [f"## {gene}", "",
              f"- Label file: `phase1/data/assay_labels/{ours['path'].name}`, {ours['rows']:,} rows, "
              f"classes {', '.join(ours['classes'])} (sha256 {ours['sha256'][:16]}...)",
              f"- Cited: {man['cited']['urn']} (Table S1, licence there "
              f"{man['cited']['license_in_table_s1']}); paper: "
              f"{cited_paper(s1[gene]) or 'none cited (the score set alone)'}"]
        for v in man["versions"]:
            cals = "; ".join(f"'{c['title']}' (investigator {c['investigator_provided']}, created "
                             f"{c['created']}, modified {c['modified']})"
                             for c in v["score_calibrations"]) or "none"
            L.append(f"- {v['urn']} ({v['role']}): score set created {v['created']}, modified "
                     f"{v['modified']}, published {v['published']}; {v['num_variants_metadata']:,} "
                     f"variants, {v['n_rows_score_table']:,} score-table rows; licence "
                     f"{v['license']['shortName']}; calibrations: {cals}")
        for miss in man.get("versions_not_served") or []:
            L.append(f"- {miss['urn']}: linked from the chain, not served by MaveDB "
                     f"(HTTP {miss['http_status']}, {miss['at_utc']})")
        L += ["", "| rule | version | detail | compared | differ | differ, all | differ, 3-50 nt | "
              "only in the file | only in the source | reproduced |",
              "|---|---|---|---|---|---|---|---|---|---|"]
        for r in rows.itertuples():
            only = ("n/a | n/a" if pd.isna(r.n_only_ours) else
                    f"{int(r.n_only_ours):,} | {int(r.n_only_source):,}")
            L.append(f"| {r.rule} | {r.version.split(':')[-1]} | {r.rule_detail} | "
                     f"{int(r.n_compared):,} | {int(r.n_differ):,} | {_fmt_pct(r.pct_differ_all)} | "
                     f"{_fmt_pct(r.pct_differ_3_50nt)} | {only} | {r.exact_match} |")
        L.append("")
        if len(exact):
            L.append("- Reproduced by: " + "; ".join(
                f"{r.rule} ({r.rule_detail}, {r.version}; class names identical as written, "
                f"case and separators folded: {int(r.n_identical_as_written):,} of "
                f"{int(r.n_compared):,})" for r in exact.itertuples()))
        else:
            L.append("- Reproduced by: no rule. Source not established.")
        names = class_name_check(res)
        L.append("- Class names of the label file in the archive (case and separators "
                 "folded): " + "; ".join(
                     f"'{c['class']}' {'yes' if c['case_and_separators_folded'] else 'NO'}"
                     for c in names))
        words = res["notes"]["column_words"]
        for k, used in words.items():
            L.append(f"- Class column {k}: " + "; ".join(
                f"{w or '(empty)'} -> {canon} ({src})" for w, (canon, src) in sorted(used.items())))
        if res["notes"]["genomic_mapping_failures"]:
            L.append(f"- Label-file variants the GRCh38 mapping could not place (excluded from "
                     f"the genomic alignment): {len(res['notes']['genomic_mapping_failures'])}")
        pr = res["notes"].get("paper_rule")
        if pr:
            L.append(f"- Paper ({pr['source']['citation_ids']}): {pr['status']}; "
                     + "; ".join(f"\"{q['text']}\" ({q['where']})" for q in pr.get("quotes", [])))
        for f in all_flags[gene]:
            L.append(f"- **Condition for the author (E16c.7):** {f}")
        if gene in AUTHOR_DECISIONS:
            L.append(f"- **Decided by the author (2026-10-21):** {AUTHOR_DECISIONS[gene]}")
        unc = unclassified(res)
        if unc:
            L.append("- Variants of the cited score table that the label file does not carry, "
                     "by stratum (n, of which in the analysis set): " + "; ".join(
                         f"{u['stratum']} {u['n']} ({u['n_in_analysis_set']})" for u in unc))
        sub = impact_df[impact_df["gene"] == gene]
        if len(sub):
            L.append("- Counts in the analysis set under the closest rule (labelled / damaging / "
                     "normal / excluded, now -> rule): " + "; ".join(
                         f"{r.stratum} {r.n_labelled_now}/{r.n_damaging_now}/{r.n_normal_now}/"
                         f"{r.n_excluded_now} -> " + (
                             "n/a" if pd.isna(r.n_labelled_rule) else
                             f"{int(r.n_labelled_rule)}/{int(r.n_damaging_rule)}/"
                             f"{int(r.n_normal_rule)}/{int(r.n_excluded_rule)}")
                         for r in sub.itertuples()))
        L.append("")
    DOC_PATH.write_text("\n".join(L))


def main() -> None:
    s1 = table_s1()
    aset = pd.read_parquet(K.SET_PATH)
    M = N._atlas_mapping()
    results, all_flags = [], {}
    for gene in C.ASSAY_LABELS:
        res = audit_gene(gene, s1, aset, M)
        results.append(res)
        all_flags[gene] = flags(res, s1)
    summary = pd.DataFrame([r for res in results for r in res["rows"]], columns=SUMMARY_COLUMNS)
    diff = pd.DataFrame([d for res in results for d in res["diffs"]], columns=DIFF_COLUMNS)
    diff["offset"] = diff["offset"].astype("Int64")
    impact_rows = [r for res in results
                   if not any(x["exact_match"] == "yes" for x in res["rows"])
                   for r in impact(res, aset)]
    impact_df = pd.DataFrame(impact_rows, columns=IMPACT_COLUMNS)
    for col in ("n_labelled_rule", "n_damaging_rule", "n_normal_rule", "n_excluded_rule"):
        impact_df[col] = impact_df[col].astype("Int64")
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    summary.to_csv(OUT_SUMMARY, index=False)
    diff.to_csv(OUT_DIFF, index=False)
    impact_df.to_csv(OUT_IMPACT, index=False)
    pd.DataFrame([u for res in results for u in unclassified(res)],
                 columns=["gene", "version", "stratum", "n", "n_in_analysis_set"]
                 ).to_csv(OUT_UNCLASSIFIED, index=False)
    write_doc(results, summary, impact_df, s1, all_flags)
    print("[E16c] label source audit (read-only)")
    for res in results:
        rows = summary[summary["gene"] == res["gene"]]
        ok = rows[rows["exact_match"] == "yes"]
        best = closest(res["rows"])
        print(f"  {res['gene']:7s} " + (f"reproduced by {', '.join(ok['rule'])}" if len(ok) else
              "not reproduced" + (f"; closest {best['rule']} differs on "
                                  f"{_fmt_pct(best['pct_differ_3_50nt'])} at 3-50 nt" if best
                                  else "; no rule applies")))
        for f in all_flags[res["gene"]]:
            print(f"          E16c.7: {f}")
    print(f"  wrote {OUT_SUMMARY}, {OUT_DIFF}, {OUT_IMPACT}, {OUT_UNCLASSIFIED}, {DOC_PATH}")


if __name__ == "__main__":
    main()
