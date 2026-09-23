"""The four main tables, as printed.

Each table is derived from outputs the pipeline has already written and is saved
as a CSV of the strings that appear in print, so that the table in the manuscript
and the file in the repository cannot drift apart. The counts that the text
states in words ("in eleven of the thirteen columns") are computed here as well,
into the same files, rather than being counted by hand.

Table 1  composition of the analysis set, in-scope strata only
Table 2  what each predictor was trained on, and how that relates to the readout
Table 3  the evidence tier each predictor reaches, in-sample and held out
Table 4  evidence by ClinVar record status

Two conventions the tables share with the figures: the canonical dinucleotides
are outside the PP3/BP4 recommendation and appear only where a pool containing
them is needed as a comparison, labelled out of scope; and SpliceAI's published
cut points are applied to the two SpliceAI columns only, because they are defined
on SpliceAI's score and on nothing else. Every other column is compared at its
own fitted threshold.

Run (PYTHONPATH=phase1):  python -m src.evid_tables
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .evid_figures import COLUMNS, NAME

REPORT_DIR = Path("reports/evidence")
OUT = REPORT_DIR / "tables"

STRATA = ["s3_10", "s11_50", "s3_50"]
STRATUM_LABEL = {"s3_10": "3–10 bp", "s11_50": "11–50 bp", "s3_50": "3–50 bp",
                 "all_1_50": "1–50 bp including the canonical dinucleotides (out of scope)"}
ARMS = ["classified", "recorded_unclassified", "unrecorded"]
# one label set for the arms across text, tables and figures
ARM_LABEL = {"classified": "Classified", "recorded_unclassified": "Recorded, unclassified",
             "unrecorded": "Unrecorded"}
# the combined Atlas score's model was selected on the BRCA1 and RAD51C assays
# (predictor_training_provenance.csv, row avi), so its held-out counts are given
# over the genes its training did not see
AVI_SEEN_IN_TRAINING = ("BRCA1", "RAD51C")


def _it(gene: str) -> str:
    """Gene symbols are italic in print; the assembler renders *...* as italic."""
    return f"*{gene}*"
TIER_SHORT = {"supporting": "Supporting", "moderate": "Moderate", "strong": "Strong",
              "very_strong": "Very strong"}


def _pct(k: int, n: int) -> str:
    return f"{k:,} ({100 * k / n:.1f}%)" if n else "–"


# ---------------------------------------------------------------------------
# Table 1
# ---------------------------------------------------------------------------
def table1() -> tuple[pd.DataFrame, pd.DataFrame]:
    sc = pd.read_csv(REPORT_DIR / "set_counts.csv")
    ins = sc[sc.stratum != "pm12"].copy()
    ins["s3_50"] = True
    genes = sorted(ins.gene.unique())
    rows = []
    for g in genes + ["Total"]:
        sub = ins if g == "Total" else ins[ins.gene == g]
        r = {"Gene": g if g == "Total" else _it(g), "Variants, 3–50 bp": f"{int(sub.n.sum()):,}"}
        for st in STRATA:
            s2 = sub if st == "s3_50" else sub[sub.stratum == st]
            lab, dmg = int(s2.n_labelled.sum()), int(s2.n_damaging.sum())
            r[f"{STRATUM_LABEL[st]}: labelled"] = f"{lab:,}" if lab else "0"
            r[f"{STRATUM_LABEL[st]}: damaging (%)"] = _pct(dmg, lab)
        rows.append(r)
    a = pd.DataFrame(rows)

    rows = []
    strict_order = [("classified", "classified", "Classified pathogenic or benign"),
                    ("recorded_unclassified", None, "Recorded without such a classification"),
                    ("recorded_unclassified", "recorded_no_assertion",
                     "   of which no clinical assertion"),
                    ("unrecorded", "unrecorded", "Unrecorded (not in the ClinVar release)")]
    for arm, strict, label in strict_order:
        sub = ins[ins.clinvar_arm == arm]
        if strict is not None:
            sub = sub[sub.clinvar_arm_strict == strict]
        lab, dmg = int(sub.n_labelled.sum()), int(sub.n_damaging.sum())
        rows.append({"ClinVar record status (3–50 bp)": label,
                     "Variants": f"{int(sub.n.sum()):,}", "Labelled": f"{lab:,}",
                     "Damaging (%)": _pct(dmg, lab),
                     "Genes": str(sub[sub.n > 0].gene.nunique())})
    b = pd.DataFrame(rows)
    return a, b


# ---------------------------------------------------------------------------
# Table 2
# ---------------------------------------------------------------------------
READOUT_SHORT = {
    "Splicing (annotated junctions)":
        "Indirect: predicts a splice event, while the assay measures gene function",
    "Splicing (measured site usage)":
        "Indirect: predicts a splice event, while the assay measures gene function",
    "Functional genomics tracks, including splicing":
        "Indirect: predicts molecular effects including splicing",
    "Proxy contrast of simulated against fixed derived variants":
        "Shared cause: the training contrast is shaped by selection",
    "Cross-species sequence constraint":
        "Shared cause: constraint and functional damage both reflect selection",
    "Cross-species sequence constraint (self-supervised)":
        "Shared cause: constraint and functional damage both reflect selection",
    "Supervised on population allele frequency, with AlphaGenome, AlphaMissense "
    "and conservation inputs":
        "Direct for *BRCA1*, *RAD51C* and *DDX3X*, whose assays selected the model; "
        "shared cause (allele frequency, constraint) otherwise",
    "Genomic sequence (self-supervised)": "None documented",
    "This study's own functional labels, refitted per fold":
        "Direct: fitted on the functional standard, so evaluated leave-one-gene-out",
}


# Short printed form of each column's terms of use (full text in the provenance
# file and Supplementary Table S3). A row group must share one form.
TERMS_SHORT = {
    "spliceai": "Non-commercial trained models (CC BY-NC 4.0); commercial use needs a licence",
    "spliceai_walker": "Non-commercial trained models (CC BY-NC 4.0); commercial use needs a licence",
    "pangolin": "Open software (GPL-3.0)",
    "alphagenome": "Research only: no clinical decision-making, no training of other models",
    "avi": "Research only: no clinical decision-making, no training of other models",
    "avi_splice_sites": "Research only: no clinical decision-making, no training of other models",
    "avi_splice_site_usage": "Research only: no clinical decision-making, no training of other models",
    "avi_splice_junctions": "Research only: no clinical decision-making, no training of other models",
    "cadd": "Non-commercial; commercial licence available",
    "phylop": "Free, including commercial use", "phastcons": "Free, including commercial use",
    "gpn_msa": "Open (MIT)", "nt": "Non-commercial (CC BY-NC-SA 4.0)",
    "fusion_enet": "Terms of its inputs",
}


def table2() -> pd.DataFrame:
    t = pd.read_csv(REPORT_DIR / "predictor_training_provenance.csv")
    names = dict(NAME, fusion_enet="Elastic-net combination (Supplementary Table S9)")
    rows = []
    for cls, g in t.groupby("training_signal_class", sort=False):
        mave = g.documented_mave_or_sge_in_training
        if (mave == "None documented").all():
            mave_cell = "None documented"
        elif (g.score_column == "fusion_enet").all():
            mave_cell = "By construction (this study's labels)"
        elif (g.score_column == "avi").all():
            mave_cell = ("Yes: four saturation genome editing assays, including the "
                         "*BRCA1*, *RAD51C* and *DDX3X* assays used here, selected "
                         "the model")
        else:
            raise SystemExit(f"[tables] Table 2: no printed form for {list(g.score_column)}")
        rows.append({
            "Score columns": "; ".join(names.get(c, c) for c in
                                       sorted(g.score_column, key=_display_order)),
            "Training signal": cls,
            "Clinical labels in training":
                "No" if (g.contains_clinical_labels == "No").all() else "Yes",
            "Multiplexed-assay measurements in training": mave_cell,
            "Relation to the functional readout": READOUT_SHORT[cls],
            "Terms of use": _one({TERMS_SHORT[c] for c in g.score_column}, cls),
        })
    return pd.DataFrame(rows)


def _display_order(c: str) -> int:
    return COLUMNS.index(c) if c in COLUMNS else len(COLUMNS)


def _one(values: set, what: str) -> str:
    if len(values) != 1:
        raise SystemExit(f"[tables] Table 2 row '{what}' mixes terms of use: {values}")
    return next(iter(values))


# ---------------------------------------------------------------------------
# Table 3
# ---------------------------------------------------------------------------
def _highest_tier(ev: pd.DataFrame, tool: str, stratum: str) -> str:
    e = ev[(ev.tool == tool) & (ev.stratum == stratum) & (ev.status == "ok")
           & ev.pp3_threshold_reachable.astype(bool)]
    order = ["supporting", "moderate", "strong", "very_strong"]
    reached = [t for t in order if t in set(e.tier)]
    return TIER_SHORT[reached[-1]] if reached else "None"


def _folds(t: pd.DataFrame, tool: str, stratum: str, tier: str) -> str:
    r = t[(t.tool == tool) & (t.stratum == stratum) & (t.tier == tier)]
    if not len(r) or r.status.iloc[0] != "ok" or not bool(r.in_sample_tier_reached.iloc[0]):
        return "–"
    r = r.iloc[0]
    n, ev_, k = int(r.n_folds), int(r.folds_with_heldout_lr), int(r.folds_heldout_lr_above_cut)
    # "k/n (m)": k held-out genes clear the cut out of n, m folds not evaluable
    return f"{k}/{n}" + (f" ({n - ev_})" if ev_ < n else "")


def _folds_excluding(t: pd.DataFrame, folds: pd.DataFrame, tool: str, stratum: str,
                     tier: str, exclude: tuple[str, ...]) -> str:
    """k/n (m) over the held-out genes not in `exclude`, from the per-fold file.

    The same rule as tier_logo_table: a fold clears when its held-out ratio is at
    or above the tier's boundary. Checked against tier_logo_table on all genes."""
    r = t[(t.tool == tool) & (t.stratum == stratum) & (t.tier == tier)]
    if not len(r) or r.status.iloc[0] != "ok" or not bool(r.in_sample_tier_reached.iloc[0]):
        return "–"
    cut = float(r.tier_lr_cut.iloc[0])
    f = folds[(folds.tool == tool) & (folds.stratum == stratum) & (folds.tier == tier)
              & (folds.side == "pp3")]
    def count(ff):
        lr = ff.heldout_lr
        return len(ff), int(lr.notna().sum()), int((lr >= cut).sum())
    n_all, ev_all, k_all = count(f)
    if (n_all, ev_all, k_all) != (int(r.n_folds.iloc[0]), int(r.folds_with_heldout_lr.iloc[0]),
                                  int(r.folds_heldout_lr_above_cut.iloc[0])):
        raise SystemExit(f"[tables] {tool} {stratum} {tier}: per-fold file disagrees "
                         "with tier_logo_table")
    n, ev_, k = count(f[~f.heldout_gene.isin(exclude)])
    return f"{k}/{n}" + (f" ({n - ev_})" if ev_ < n else "")


def table3() -> pd.DataFrame:
    ev = pd.read_csv(REPORT_DIR / "evidence_thresholds.csv")
    t = pd.read_csv(REPORT_DIR / "tier_logo_table.csv")
    folds = pd.read_csv(REPORT_DIR / "evidence_thresholds_logo_folds.csv")
    rows = []
    for c in COLUMNS:
        r = {"Predictor": NAME[c] + (" †" if c == "avi" else "")}
        for st in STRATA:
            lab = STRATUM_LABEL[st]
            r[f"{lab}: In-sample"] = _highest_tier(ev, c, st)
            for tier, head in (("moderate", "Moderate, held out"), ("strong", "Strong, held out")):
                r[f"{lab}: {head}"] = (
                    _folds_excluding(t, folds, c, st, tier, AVI_SEEN_IN_TRAINING)
                    if c == "avi" else _folds(t, c, st, tier))
        rows.append(r)
    return pd.DataFrame(rows)


def table3_counts() -> pd.DataFrame:
    """Counts the text states in words, derived from the same file as Table 3."""
    ev = pd.read_csv(REPORT_DIR / "evidence_thresholds.csv")
    t = pd.read_csv(REPORT_DIR / "tier_logo_table.csv")
    rows = []
    for st in STRATA:
        for tier in ("moderate", "strong"):
            e = ev[(ev.stratum == st) & (ev.tier == tier) & ev.tool.isin(COLUMNS)
                   & (ev.status == "ok")]
            reached = sorted(e[e.pp3_threshold_reachable.astype(bool)].tool)
            rows.append({"stratum": st, "tier": tier,
                         "n_columns_reaching_in_sample": len(reached),
                         "columns_reaching_in_sample": "; ".join(reached)})
            tt = t[(t.stratum == st) & (t.tier == tier) & t.tool.isin(COLUMNS)
                   & (t.status == "ok")]
            for r in tt.itertuples():
                rows.append({"stratum": st, "tier": tier, "tool": r.tool,
                             "n_folds": r.n_folds, "folds_evaluable": r.folds_with_heldout_lr,
                             "folds_clearing": r.folds_heldout_lr_above_cut})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Table 4
# ---------------------------------------------------------------------------
def table4() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    at = pd.read_csv(REPORT_DIR / "arms_at_tool_threshold.csv")
    ok = at[at.status == "ok"]

    # (a) SpliceAI, published basis, at the published cut point
    rows = []
    pub = ok[(ok.tool == "spliceai_walker") & (ok.threshold_basis == "published cut point")]
    for gs, gl in (("all_genes", "Seven genes"), ("no_BRCA1", "Without *BRCA1*")):
        for st in STRATA:
            s = pub[(pub.gene_set == gs) & (pub.stratum == st)]
            r = {"Genes": gl, "Band": STRATUM_LABEL[st]}
            for arm in ARMS:
                v = s[s.clinvar_arm == arm]
                r[ARM_LABEL[arm]] = f"{v.band_lr.iloc[0]:.1f}" if len(v) and np.isfinite(v.band_lr.iloc[0]) else "n.e."
            rows.append(r)
    a = pd.DataFrame(rows)

    # (b) every column at its own fitted threshold: which arm is highest
    fitted = ok[ok.threshold_basis.str.startswith("fitted") & ok.tool.isin(COLUMNS)]
    rows, detail = [], []
    for gs, gl in (("all_genes", "Seven genes"), ("no_BRCA1", "Without *BRCA1*")):
        for st in STRATA:
            s = fitted[(fitted.gene_set == gs) & (fitted.stratum == st)
                       & fitted.clinvar_arm.isin(ARMS)]
            p = s.pivot_table(index="tool", columns="clinvar_arm", values="band_lr")
            p = p.reindex(columns=ARMS).dropna()
            top = p.idxmax(axis=1)
            r = {"Genes": gl, "Band": STRATUM_LABEL[st],
                 "Columns with all three arms estimable": str(len(p))}
            for arm in ARMS:
                r[f"Highest: {ARM_LABEL[arm]}"] = str(int((top == arm).sum()))
            rows.append(r)
            for tool, arm in top.items():
                detail.append({"gene_set": gs, "stratum": st, "tool": tool,
                               "highest_arm": arm, **{f"lr_{a2}": p.loc[tool, a2] for a2 in ARMS}})
    b = pd.DataFrame(rows)

    # (c) AUROC, threshold-free, every column
    am = pd.read_csv(REPORT_DIR / "territory_metrics_arms_noBRCA1.csv")
    am = am[(am.status == "ok") & am.tool.isin(COLUMNS)]
    rows = []
    for st in ("s3_50", "all_1_50"):
        for gs, gl in (("all_genes", "Seven genes"), ("no_BRCA1", "Without *BRCA1*")):
            s = am[(am.stratum == st) & (am.gene_set == gs)]
            p = s.pivot_table(index="tool", columns="clinvar_arm", values="auroc")
            p = p.reindex(columns=ARMS).dropna()
            top, bottom = p.idxmax(axis=1), p.idxmin(axis=1)
            r = {"Pool": STRATUM_LABEL[st], "Genes": gl, "Columns": str(len(p))}
            for arm in ARMS:
                r[f"Highest AUROC: {ARM_LABEL[arm]}"] = str(int((top == arm).sum()))
            r["Lowest AUROC: Classified"] = str(int((bottom == "classified").sum()))
            rows.append(r)
    c = pd.DataFrame(rows)
    return a, b, c, pd.DataFrame(detail)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    t1a, t1b = table1()
    t1a.to_csv(OUT / "table1a_composition_by_gene.csv", index=False)
    t1b.to_csv(OUT / "table1b_composition_by_clinvar_status.csv", index=False)
    table2().to_csv(OUT / "table2_training_signal.csv", index=False)
    table3().to_csv(OUT / "table3_evidence_tiers.csv", index=False)
    table3_counts().to_csv(OUT / "table3_counts.csv", index=False)
    a, b, c, d = table4()
    a.to_csv(OUT / "table4a_spliceai_published_cut_by_arm.csv", index=False)
    b.to_csv(OUT / "table4b_fitted_threshold_highest_arm.csv", index=False)
    c.to_csv(OUT / "table4c_auroc_highest_arm.csv", index=False)
    d.to_csv(OUT / "table4_detail_fitted_threshold_by_arm.csv", index=False)
    for f in sorted(OUT.glob("*.csv")):
        print(f"[tables] wrote {f.name} ({len(pd.read_csv(f))} rows)")


if __name__ == "__main__":
    main()
