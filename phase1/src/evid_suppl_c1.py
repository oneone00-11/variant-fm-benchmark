"""E19 C1 -- the predictors against ClinVar, as a concordance check.

ClinVar is not the standard here: it is the evaluation the evidence feeds, and this
stage only asks how the frozen thresholds and the published cut points read against
it. No likelihood ratio is formed from ClinVar labels and no evidence tier is read
(E19.0 rule 3): sensitivity, specificity, AUROC and counts only.

The reference. The analysis set's germline classification, `clnsig_class`, the
project's collapse of the VCF's CLNSIG ("Aggregate germline classification for this
single variant", ClinVar VCF of 2026-06-15; evid_build_set.annotate_clinvar): positive
P/LP, negative B/LB; Conflicting, VUS, Other and unrecorded variants carry no
reference label. Stars come from CLNREVSTAT ("ClinVar review status of germline
classification for the Variation ID") by the conversion ClinVar's documentation gives,
read from the page as saved (phase1/data/evidence/e19/clinvar_review_status/) and
written out beside the results; a review status the page does not grade stops the
stage. The main analysis keeps assertions of at least one star, the sensitivity
analysis at least two.

The bands: 3-10, 11-50 and 3-50 nt. A band with fewer than MIN_POS P/LP or MIN_NEG
B/LB variants scored by a predictor is not evaluable there, and its row gives counts
only.

For every predictor of E4 (the thirteen columns and the fusion, rebuilt leave one gene
out as E3, E4 and E15 rebuild it), in every band:
  * at two kinds of threshold, the counts, sensitivity and specificity with E15's
    Wilson intervals, by E15's own functions and names: Walker's published cut points
    for both SpliceAI columns (PP3 >= 0.2, BP4 <= 0.1), and the in-sample thresholds of
    Table 3(b), read at full precision from the frozen evidence_thresholds.csv (PP3 at
    Supporting, Moderate and Strong, BP4 at Supporting and Moderate, where E3 reached
    the tier in-sample). On the BP4 side E15's names hold: bp4_sensitivity is the share
    of negatives at or below the threshold, bp4_specificity the share of positives
    above it.
  * AUROC with its interval, by E4's own estimator: per gene, at least MIN_POS and
    MIN_NEG in the gene, pooled across genes on the logit scale (Hanley-McNeil
    variances, DerSimonian-Laird), as the companion atlas pools it. Each gene's AUROC
    is written beside the pooled one. Where a gene's AUROC is 0 or 1 its Hanley-McNeil
    variance is 0 and the logit pooling is undefined; the row says so ("undefined")
    and nothing else is put in its place.

The cross-tabulation: ClinVar (P/LP, B/LB, unrecorded) x SGE (damaging, normal) x the
predictor's side of each threshold, as counts. In the P/LP group the share of SGE
damaging variants at or above each PP3 threshold (E15's `sensitivity`); in the B/LB
group the share of SGE normal variants at or below each BP4 threshold (E15's
`bp4_sensitivity`; the author's decision of 2026-10-25); in the unrecorded group both,
for comparison. No likelihood ratio.

Outputs, in reports/evidence/suppl_analyses/:
  c1_clinvar_star_conversion.csv     ClinVar's review statuses and their stars, as the
                                     saved page gives them, with the analysis set's
                                     count of each
  c1_clinvar_reference_metrics.csv   the metrics above, per star level, predictor, band
                                     and threshold, and the AUROC rows
  c1_clinvar_sge_crosstab.csv        the cross-tabulation

Run (PYTHONPATH=phase1):  python -m src.evid_suppl_c1
"""
from __future__ import annotations

import html
import re
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from . import evid_common as K
from . import evid_practical_metrics as P
from . import evid_territory_metrics as T

REPORT_DIR = Path("reports/evidence")
OUT_DIR = REPORT_DIR / "suppl_analyses"
OUT_STARS = OUT_DIR / "c1_clinvar_star_conversion.csv"
OUT_METRICS = OUT_DIR / "c1_clinvar_reference_metrics.csv"
OUT_CROSSTAB = OUT_DIR / "c1_clinvar_sge_crosstab.csv"
REVIEW_STATUS_PAGE = Path("data/evidence/e19/clinvar_review_status/review_status.html")
# the page's section for the review status of a germline classification on aggregate
# records (VCV and RCV), which is CLNREVSTAT's, and how it writes the number of stars
TABLE_HEADINGS = ("Review status on aggregate records (VCV and RCV)",
                  "Germline classification and oncogenicity")
STAR_WORDS = {"four": 4, "three": 3, "two": 2, "one": 1, "none": 0}
CONFIG_PATH = Path("config/walker2023.yaml")

BANDS = ["s3_10", "s11_50", "s3_50"]
STAR_LEVELS = [("stars_ge_1", 1, "main"), ("stars_ge_2", 2, "sensitivity")]
PP3_TIERS = ["supporting", "moderate", "strong"]
BP4_TIERS = ["supporting", "moderate"]
SPLICEAI = ("spliceai", "spliceai_walker")
GROUPS = [("P/LP", "plp"), ("B/LB", "blb"), ("unrecorded", "unrecorded")]
# E15's names for the counts and the two shares a row keeps; its band ratios, tiers,
# posteriors, coverage and predictive values are not carried (E19.0 rule 3)
PP3_KEEP = ["n_pos_above", "n_neg_above", "n_pos_below", "n_neg_below",
            "sensitivity", "sensitivity_lo", "sensitivity_hi",
            "specificity", "specificity_lo", "specificity_hi"]
BP4_KEEP = ["bp4_n_neg_at_or_below", "bp4_n_pos_at_or_below", "bp4_n_neg_above",
            "bp4_n_pos_above", "bp4_sensitivity", "bp4_sensitivity_lo", "bp4_sensitivity_hi",
            "bp4_specificity", "bp4_specificity_lo", "bp4_specificity_hi"]


# ---------------------------------------------------------------------------
# ClinVar's stars, from its documentation as saved
# ---------------------------------------------------------------------------
def _norm(status: str) -> str:
    """A review status as ClinVar's pages write it: the VCF's underscores as spaces,
    runs of white space as one, lower case."""
    return " ".join(status.replace("_", " ").split()).lower()


def _plain(s: str) -> str:
    return " ".join(html.unescape(re.sub(r"<[^>]+>", " ", s)).split())


def star_table(page: Path = REVIEW_STATUS_PAGE) -> pd.DataFrame:
    """(review status, stars) as the saved page gives them for the review status
    CLNREVSTAT carries, that of the germline classification of a Variation ID: the
    table under the page's headings TABLE_HEADINGS (its other two star tables are for
    submitted records and for somatic classifications). Its rows are read as they
    stand; a header other than the page's, or a number of stars written other than in
    STAR_WORDS, stops the stage."""
    t = page.read_text(encoding="utf-8")
    heads = [(m.start(), m.end(), _plain(m.group(0)))
             for m in re.finditer(r"<h([1-6])[^>]*>.*?</h\1>", t, re.S | re.I)]
    try:
        i = next(k for k, h in enumerate(heads) if h[2] == TABLE_HEADINGS[0])
        j = next(k for k, h in enumerate(heads) if k > i and h[2] == TABLE_HEADINGS[1])
    except StopIteration:
        raise SystemExit(f"[C1] {page}: no headings {TABLE_HEADINGS}: stopped") from None
    seg = t[heads[j][1]:heads[j + 1][0] if j + 1 < len(heads) else len(t)]
    m = re.search(r"<table.*?</table>", seg, re.S | re.I)
    if m is None:
        raise SystemExit(f"[C1] {page}: no table under {TABLE_HEADINGS}: stopped")
    cells = [[_plain(c) for c in re.findall(r"<t[dh][^>]*>(.*?)</t[dh]>", r, re.S | re.I)]
             for r in re.findall(r"<tr[^>]*>(.*?)</tr>", m.group(0), re.S | re.I)]
    if cells[0] != ["Number of gold stars", "Review status", "Description"]:
        raise SystemExit(f"[C1] {page}: the table's header is {cells[0]}: stopped")
    rows = []
    for c in cells[1:]:
        if len(c) != 3 or c[0] not in STAR_WORDS:
            raise SystemExit(f"[C1] {page}: a row of the star table reads {c}: stopped")
        rows.append({"review_status": c[1], "stars_as_written": c[0], "stars": STAR_WORDS[c[0]],
                     "description": c[2], "table": " / ".join(TABLE_HEADINGS)})
    return pd.DataFrame(rows)


def stars_of(df: pd.DataFrame, table: pd.DataFrame) -> pd.Series:
    lut = dict(zip(table["review_status"].map(_norm), table["stars"]))
    seen = {_norm(s) for s in df["clnrevstat"].dropna() if s}
    missing = sorted(seen - set(lut))
    if missing:
        raise SystemExit(f"[C1] review statuses the saved page does not grade: {missing}")
    return df["clnrevstat"].map(lambda s: lut[_norm(s)] if isinstance(s, str) and s else np.nan)


# ---------------------------------------------------------------------------
# the thresholds
# ---------------------------------------------------------------------------
def thresholds(tool: str, stratum: str, ev: pd.DataFrame, cfg: dict) -> list[tuple]:
    """(side, threshold_kind, threshold_source, threshold) in a fixed order."""
    out = []
    if tool in SPLICEAI:
        out += [("pp3", "walker_pp3", "published", float(cfg["thresholds"]["pp3"]["value"])),
                ("bp4", "walker_bp4", "published", float(cfg["thresholds"]["bp4"]["value"]))]
    e = ev[(ev.tool == tool) & (ev.stratum == stratum) & (ev.status == "ok")].set_index("tier")
    for t in PP3_TIERS:
        if t in e.index and bool(e.loc[t, "pp3_threshold_reachable"]):
            out.append(("pp3", t, "in_sample", float(e.loc[t, "pp3_threshold_insample"])))
    for t in BP4_TIERS:
        if t in e.index and bool(e.loc[t, "bp4_threshold_reachable"]):
            out.append(("bp4", f"bp4_{t}", "in_sample", float(e.loc[t, "bp4_threshold_insample"])))
    return out


def _per_gene_auroc(sub: pd.DataFrame, tool: str, ce) -> str:
    """Each gene's AUROC on the genes E4's estimator keeps (T.per_gene_and_pool's rule),
    with its P/LP and B/LB counts."""
    out = []
    for gene, g in sub.groupby("gene", observed=True):
        y = g["y_assay"].to_numpy(dtype=float)
        s = g[tool].to_numpy(dtype=float)
        m = np.isfinite(y) & np.isfinite(s)
        y, s = y[m], s[m]
        n_pos, n_neg = int((y == 1).sum()), int((y == 0).sum())
        if n_pos >= T.MIN_POS and n_neg >= T.MIN_NEG:
            out.append(f"{gene} {ce.roc_auc(y, s)!r} ({n_pos} P/LP, {n_neg} B/LB)")
    return "; ".join(out)


def _side_values(c: dict, side: str, bands: tuple[dict, dict]) -> dict:
    v = P._values(c, side, bands)
    return {k: v[k] for k in (PP3_KEEP if side == "pp3" else BP4_KEEP)}


# ---------------------------------------------------------------------------
# the stage
# ---------------------------------------------------------------------------
def reference_metrics(df: pd.DataFrame, tools: list[str], ev: pd.DataFrame, cfg: dict,
                      ce, levels: list[tuple] = STAR_LEVELS) -> pd.DataFrame:
    """The metrics on the variants whose y_clinvar is set and whose stars reach each
    level; the low-circularity subset (evid_suppl_c1_low) passes the main level only."""
    bands = K.acmg_bands(cfg)
    rows = []
    for level, min_stars, role in levels:
        ref = df[df["y_clinvar"].notna() & (df["clinvar_stars"] >= min_stars)]
        for tool in tools:
            for stratum in BANDS:
                frame = K.stratum_frame(ref, stratum)
                c0 = P._cell(frame, tool, "y_clinvar", np.inf)
                base = {"star_level": level, "analysis": role, "tool": tool,
                        "stratum": stratum, "n_scored": c0["n_scored"],
                        "n_plp": c0["n_pos"], "n_blb": c0["n_neg"]}
                ok = c0["n_pos"] >= K.MIN_POS and c0["n_neg"] >= K.MIN_NEG
                status = "ok" if ok else "not evaluable"
                note = "" if ok else (f"fewer than {K.MIN_POS} P/LP or {K.MIN_NEG} B/LB "
                                      "variants scored: counts only")
                # AUROC, E4's estimator on the ClinVar label
                auc = {"status": status, "note": note}
                if ok:
                    sub = frame.dropna(subset=["y_clinvar", tool]).assign(
                        y_assay=lambda d: d["y_clinvar"])
                    p = T.per_gene_and_pool(sub, tool, ce)
                    per_gene = _per_gene_auroc(sub, tool, ce)
                    if p["status"] == "ok" and not np.isfinite(p.get("auroc", np.nan)):
                        p = p | {"status": "undefined", "reason": (
                            "a gene's AUROC is 0 or 1, so its Hanley-McNeil variance is 0 and "
                            "E4's logit pooling is undefined")}
                    auc = {"status": p["status"], "note": p.get("reason", ""),
                           "auroc_per_gene": per_gene,
                           "auroc": p.get("auroc", np.nan), "auroc_lo": p.get("auroc_lo", np.nan),
                           "auroc_hi": p.get("auroc_hi", np.nan),
                           "auroc_k_genes": p.get("k_genes", np.nan),
                           "auroc_n_plp_used": p.get("n_pos", np.nan),
                           "auroc_n_blb_used": p.get("n_neg", np.nan)}
                rows.append(base | {"measure": "auroc", "side": "", "threshold_kind": "",
                                    "threshold_source": ""} | auc)
                for side, kind, source, thr in thresholds(tool, stratum, ev, cfg):
                    row = base | {"measure": "threshold", "side": side,
                                  "threshold_kind": kind, "threshold_source": source,
                                  "threshold": thr, "status": status, "note": note}
                    c = P._cell(frame, tool, "y_clinvar", thr)
                    vals = _side_values(c, side, bands)
                    if not ok:
                        # counts only
                        vals = {k: v for k, v in vals.items() if "_n_" in f"_{k}" or k.startswith("n_")}
                    rows.append(row | vals)
    lead = ["star_level", "analysis", "tool", "stratum", "measure", "side", "threshold_kind",
            "threshold_source", "threshold", "status", "note", "n_scored", "n_plp", "n_blb"]
    auc_cols = ["auroc", "auroc_lo", "auroc_hi", "auroc_k_genes", "auroc_n_plp_used",
                "auroc_n_blb_used", "auroc_per_gene"]
    return pd.DataFrame(rows).reindex(columns=lead + PP3_KEEP + BP4_KEEP + auc_cols)


def crosstab(df: pd.DataFrame, tools: list[str], ev: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    bands = K.acmg_bands(cfg)
    rows = []
    for level, min_stars, role in STAR_LEVELS:
        member = {"plp": (df["clnsig_class"] == "P/LP") & (df["clinvar_stars"] >= min_stars),
                  "blb": (df["clnsig_class"] == "B/LB") & (df["clinvar_stars"] >= min_stars),
                  "unrecorded": df["clinvar_arm"] == "unrecorded"}
        for tool in tools:
            for stratum in BANDS:
                for side, kind, source, thr in thresholds(tool, stratum, ev, cfg):
                    for group, short in GROUPS:
                        frame = K.stratum_frame(df[member[short]], stratum)
                        c = P._cell(frame, tool, "y_assay", thr)
                        row = {"star_level": level, "analysis": role, "tool": tool,
                               "stratum": stratum, "side": side, "threshold_kind": kind,
                               "threshold_source": source, "threshold": thr,
                               "clinvar_group": group, "n_scored": c["n_scored"],
                               "n_sge_damaging": c["n_pos"], "n_sge_normal": c["n_neg"]}
                        v = P._values(c, side, bands)
                        if side == "pp3":
                            row |= {k: v[k] for k in ("n_pos_above", "n_pos_below",
                                                      "n_neg_above", "n_neg_below")}
                            if group in ("P/LP", "unrecorded"):
                                row |= {k: v[k] for k in ("sensitivity", "sensitivity_lo",
                                                          "sensitivity_hi")}
                        else:
                            row |= {k: v[k] for k in ("bp4_n_pos_at_or_below", "bp4_n_pos_above",
                                                      "bp4_n_neg_at_or_below", "bp4_n_neg_above")}
                            if group in ("B/LB", "unrecorded"):
                                row |= {k: v[k] for k in ("bp4_sensitivity", "bp4_sensitivity_lo",
                                                          "bp4_sensitivity_hi")}
                        rows.append(row)
    cols = ["star_level", "analysis", "tool", "stratum", "side", "threshold_kind",
            "threshold_source", "threshold", "clinvar_group", "n_scored", "n_sge_damaging",
            "n_sge_normal", "n_pos_above", "n_pos_below", "n_neg_above", "n_neg_below",
            "sensitivity", "sensitivity_lo", "sensitivity_hi",
            "bp4_n_pos_at_or_below", "bp4_n_pos_above", "bp4_n_neg_at_or_below",
            "bp4_n_neg_above", "bp4_sensitivity", "bp4_sensitivity_lo", "bp4_sensitivity_hi"]
    return pd.DataFrame(rows).reindex(columns=cols)


def main() -> None:
    cfg = yaml.safe_load(CONFIG_PATH.read_text())
    df = K.load_set()
    tools = K.panel_of(df)
    df[K.FUSION] = K.logo_fusion(df, [t for t in K.FUSION_FEATURES if t in df.columns])
    tools = tools + [K.FUSION]
    table = star_table()
    df["clinvar_stars"] = stars_of(df, table)
    df["y_clinvar"] = df["clnsig_class"].map({"P/LP": 1.0, "B/LB": 0.0})
    ev = K.read_back(REPORT_DIR / "evidence_thresholds.csv")
    ce = T._atlas_clinical_evidence()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    counts = df["clnrevstat"].fillna("").map(_norm).value_counts()
    table.assign(n_variants_in_set=table["review_status"].map(_norm).map(counts).fillna(0)
                 .astype(int)).to_csv(OUT_STARS, index=False)
    m = reference_metrics(df, tools, ev, cfg, ce)
    m.to_csv(OUT_METRICS, index=False)
    x = crosstab(df, tools, ev, cfg)
    x.to_csv(OUT_CROSSTAB, index=False)
    print(f"[C1] wrote {OUT_STARS}, {OUT_METRICS} ({len(m)} rows) and {OUT_CROSSTAB} "
          f"({len(x)} rows)")
    a = m[(m.measure == "auroc") & (m.star_level == "stars_ge_1")]
    print(a.pivot_table(index="tool", columns="stratum", values="auroc").round(3).to_string())


if __name__ == "__main__":
    main()
