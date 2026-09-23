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
                 "all_1_50": "1–50 bp incl. ±1,2 (out of scope)"}
ARMS = ["classified", "recorded_unclassified", "unrecorded"]
ARM_LABEL = {"classified": "Classified", "recorded_unclassified": "Recorded, unclassified",
             "unrecorded": "Not in ClinVar"}
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
        r = {"Gene": g, "Variants, 3–50 bp": f"{int(sub.n.sum()):,}"}
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
                    ("unrecorded", "unrecorded", "Not in the ClinVar release")]
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
    "Genomic sequence (self-supervised)": "None documented",
    "This study's own functional labels, refitted per fold":
        "Direct: fitted on the functional standard, so evaluated leave-one-gene-out",
}


def table2() -> pd.DataFrame:
    t = pd.read_csv(REPORT_DIR / "predictor_training_provenance.csv")
    names = dict(NAME, fusion_enet="Elastic-net combination (Supplementary Table S9)")
    rows = []
    for cls, g in t.groupby("training_signal_class", sort=False):
        rows.append({
            "Score columns": "; ".join(names.get(c, c) for c in g.score_column),
            "Training signal": cls,
            "Clinical classifications in training":
                "None" if (g.contains_clinical_labels == "No").all() else "Yes",
            "Multiplexed-assay measurements in training":
                "None documented" if (g.documented_mave_or_sge_in_training
                                      == "None documented").all() else "This study's",
            "Relation to the functional readout": READOUT_SHORT[cls],
        })
    return pd.DataFrame(rows)


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


def table3() -> pd.DataFrame:
    ev = pd.read_csv(REPORT_DIR / "evidence_thresholds.csv")
    t = pd.read_csv(REPORT_DIR / "tier_logo_table.csv")
    rows = []
    for c in COLUMNS:
        r = {"Predictor": NAME[c]}
        for st in STRATA:
            lab = STRATUM_LABEL[st]
            r[f"{lab}: In-sample"] = _highest_tier(ev, c, st)
            r[f"{lab}: Moderate, held out"] = _folds(t, c, st, "moderate")
            r[f"{lab}: Strong, held out"] = _folds(t, c, st, "strong")
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
    for gs, gl in (("all_genes", "Seven genes"), ("no_BRCA1", "Without BRCA1")):
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
    for gs, gl in (("all_genes", "Seven genes"), ("no_BRCA1", "Without BRCA1")):
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
        for gs, gl in (("all_genes", "Seven genes"), ("no_BRCA1", "Without BRCA1")):
            s = am[(am.stratum == st) & (am.gene_set == gs)]
            p = s.pivot_table(index="tool", columns="clinvar_arm", values="auroc")
            p = p.reindex(columns=ARMS).dropna()
            top = p.idxmax(axis=1)
            r = {"Pool": STRATUM_LABEL[st], "Genes": gl, "Columns": str(len(p))}
            for arm in ARMS:
                r[f"Highest AUROC: {ARM_LABEL[arm]}"] = str(int((top == arm).sum()))
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
