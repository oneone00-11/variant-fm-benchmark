"""E19 C8, C9, C10 -- read from existing outputs; nothing is fitted (E19.9).

C8, the distal band without BAP1.
  * composition: the labelled variants of the 11-50 nt band by gene (the analysis set),
    each gene's damaging variants and their share of the band's damaging variants;
  * thresholds: from E11's threshold_dilution.csv, the in-sample threshold and whether
    the tier is reached for the subset of every gene but BAP1 at 11-50 and 3-50 nt
    (five genes at 11-50 nt, where BRCA2 has none; six at 3-50), and that subset's
    held-out result on BAP1, beside the fit on every gene (E11's four columns, the PP3
    side, Moderate and Strong).

C9, family and function held out.
  (a) from threshold_dilution.csv, two splits: the homologous recombination genes
      (BRCA1, BRCA2, BARD1, PALB2, RAD51C; at 11-50 nt the four of them with variants
      there) as the training subset, VHL and BAP1 held out; and VHL and BAP1 as the
      training subset, each of those genes held out. E11 fitted every subset, so both
      splits are rows it already has; exactly one row per training subset and held-out
      gene, or the stage stops.
  (b) the external genes, descriptively: RAD51C's paralogues RAD51D and XRCC2 against
      the other five external genes (CTCF, SFPQ, SBDS, DDX3X, TP53), each column, band,
      side and tier, counted from Supplementary Table S12 by the supplement's own
      external_summary (k/n (m), Table 3's way, the cells that are no external test left
      out as there), and beside them the seven external genes' row of
      external_summary.csv as it stands.

C10, confusion matrices, MCC and F1, from the counts in
practical_metrics_by_threshold.csv, for its in-sample and its held-out (leave-one-gene-
out) rows. The formulas are written in the file's header:
  PP3 side: damaging is the positive class, a score at or above the threshold the call:
      TP = n_pos_above, FP = n_neg_above, FN = n_pos_below, TN = n_neg_below;
  BP4 side: normal is the positive class, a score at or below the threshold the call:
      TP = bp4_n_neg_at_or_below, FP = bp4_n_pos_at_or_below, FN = bp4_n_neg_above,
      TN = bp4_n_pos_above;
  precision = TP / (TP + FP); recall = TP / (TP + FN); F1 = 2 TP / (2 TP + FP + FN);
  MCC = (TP TN - FP FN) / sqrt((TP + FP)(TP + FN)(TN + FP)(TN + FN)), undefined (empty)
  where a factor under the root is zero.
For the supplement only.

Outputs, in reports/evidence/suppl_analyses/:
  c8_without_bap1.csv       part 'composition' and part 'thresholds'
  c9_family_holdout.csv     part '(a) family and function held out' and part
                            '(b) external genes'
  c10_confusion_mcc_f1.csv

Run (PYTHONPATH=phase1):  python -m src.evid_suppl_c8_c10
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

from . import evid_common as K
from . import evid_supp_tables as S

OUT_DIR = Path("reports/evidence/suppl_analyses")
OUT_C8 = OUT_DIR / "c8_without_bap1.csv"
OUT_C9 = OUT_DIR / "c9_family_holdout.csv"
OUT_C10 = OUT_DIR / "c10_confusion_mcc_f1.csv"
DILUTION = Path("reports/evidence/threshold_dilution.csv")
DILUTION_SUMMARY = Path("reports/evidence/threshold_dilution_summary.csv")
PRACTICAL = Path("reports/evidence/practical_metrics_by_threshold.csv")
S12 = Path("reports/evidence/supplement/tableS12_external_genes.csv")
EXT_SUMMARY = Path("reports/evidence/external_summary.csv")
HR = ["BRCA1", "BRCA2", "BARD1", "PALB2", "RAD51C"]
FUNCTION = ["VHL", "BAP1"]
PARALOGUES = ["RAD51D", "XRCC2"]
OTHER_EXTERNAL = [g for g in S.S12_GENES if g not in PARALOGUES]
HELDOUT_COLS = ["heldout_gene", "heldout_n_damaging", "heldout_n_normal", "heldout_n_in_band",
                "heldout_n_damaging_in_band", "heldout_n_normal_in_band", "heldout_lr",
                "evaluable", "cleared_cut", "not_evaluable_reason"]
C10_HEADER = """E19 C10: confusion matrices, MCC and F1 at the study's thresholds (src/evid_suppl_c8_c10.py).
Counts from practical_metrics_by_threshold.csv, in-sample and held-out (leave-one-gene-out) rows.
PP3 side: damaging is the positive class, score at or above the threshold the call:
  TP = n_pos_above, FP = n_neg_above, FN = n_pos_below, TN = n_neg_below
BP4 side: normal is the positive class, score at or below the threshold the call:
  TP = bp4_n_neg_at_or_below, FP = bp4_n_pos_at_or_below, FN = bp4_n_neg_above, TN = bp4_n_pos_above
precision = TP / (TP + FP); recall = TP / (TP + FN); F1 = 2 TP / (2 TP + FP + FN)
MCC = (TP x TN - FP x FN) / sqrt((TP + FP)(TP + FN)(TN + FP)(TN + FN)); empty where a factor is zero
For the supplement only."""


def _genes_of(g: pd.DataFrame) -> set[str]:
    genes = set(g["heldout_gene"])
    for t in g["training_genes"]:
        genes |= set(t.split(";"))
    return genes


def _one(g: pd.DataFrame, train: set[str], held: str, what: str) -> pd.Series:
    hit = g[g["training_genes"].map(lambda t: set(t.split(";")) == train)
            & (g["heldout_gene"] == held)]
    if len(hit) != 1:
        raise SystemExit(f"[C8-C10] {what}: threshold_dilution.csv has {len(hit)} rows for "
                         f"training {sorted(train)}, held out {held}: stopped")
    return hit.iloc[0]


def c8(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    d = df[(df["stratum"] == "s11_50") & df["y_assay"].notna()]
    n_dam = int((d["y_assay"] == 1).sum())
    for gene in sorted(d["gene"].unique()):
        g = d[d["gene"] == gene]
        k = int((g["y_assay"] == 1).sum())
        rows.append({"part": "composition", "stratum": "s11_50", "gene": gene,
                     "n_labelled": len(g), "n_damaging": k, "n_normal": int((g["y_assay"] == 0).sum()),
                     "share_of_band_damaging": k / n_dam})
    dil = K.read_back(DILUTION)
    summ = K.read_back(DILUTION_SUMMARY)
    for (tool, st, tier), g in dil[dil["stratum"].isin(["s11_50", "s3_50"])].groupby(
            ["tool", "stratum", "tier"], sort=False):
        genes = _genes_of(g)
        ref = float(summ[(summ.tool == tool) & (summ.stratum == st)
                         & (summ.tier == tier)]["threshold_all_genes"].iloc[0])
        base = {"part": "thresholds", "tool": tool, "stratum": st, "tier": tier,
                "path_cut_lr": float(g["path_cut_lr"].iloc[0]),
                "threshold_every_gene": ref, "tier_reached_every_gene": bool(np.isfinite(ref))}
        rows.append(base | {"subset": "every gene", "n_genes": len(genes),
                            "training_genes": ";".join(sorted(genes)), "threshold": ref,
                            "tier_reached": bool(np.isfinite(ref))})
        train = genes - {"BAP1"}
        r = _one(g, train, "BAP1", f"C8 {tool} {st} {tier}")
        rows.append(base | {"subset": "without BAP1", "n_genes": len(train),
                            "training_genes": ";".join(sorted(train)),
                            "threshold": float(r.threshold), "tier_reached": bool(r.tier_reached)}
                    | {c: r[c] for c in HELDOUT_COLS})
    return pd.DataFrame(rows)


def c9a() -> pd.DataFrame:
    dil = K.read_back(DILUTION)
    summ = K.read_back(DILUTION_SUMMARY)
    rows = []
    for (tool, st, tier), g in dil.groupby(["tool", "stratum", "tier"], sort=False):
        genes = _genes_of(g)
        hr = set(HR) & genes
        ref = float(summ[(summ.tool == tool) & (summ.stratum == st)
                         & (summ.tier == tier)]["threshold_all_genes"].iloc[0])
        for split, train, held in (
                ("homologous recombination genes in training, VHL and BAP1 held out", hr,
                 FUNCTION),
                ("VHL and BAP1 in training, homologous recombination genes held out",
                 set(FUNCTION), sorted(hr))):
            for h in held:
                r = _one(g, train, h, f"C9 {tool} {st} {tier}")
                rows.append({"part": "(a) family and function held out", "split": split,
                             "tool": tool, "predictor": S._tool(tool), "stratum": st,
                             "side": "PP3", "tier": tier,
                             "path_cut_lr": float(r.path_cut_lr),
                             "training_genes": ";".join(sorted(train)),
                             "threshold": float(r.threshold), "tier_reached": bool(r.tier_reached),
                             "threshold_every_gene": ref}
                            | {c: r[c] for c in HELDOUT_COLS})
    return pd.DataFrame(rows)


def c9b() -> pd.DataFrame:
    s12 = K.read_back(S12, dtype=str, keep_default_na=False)
    tools = {S._tool(t): t for t in S.TOOLS}
    sets = [("RAD51C paralogues: RAD51D, XRCC2", PARALOGUES),
            ("the other five external genes: " + ", ".join(OTHER_EXTERNAL), OTHER_EXTERNAL)]
    parts = [S.external_summary(s12, name, genes) for name, genes in sets]
    seven = K.read_back(EXT_SUMMARY, dtype=str, keep_default_na=False)
    parts.append(seven[seven["Gene set"] == "Seven external genes"])
    t = pd.concat(parts, ignore_index=True)
    return pd.DataFrame({
        "part": "(b) external genes", "split": t["Gene set"], "tool": t["Tool"].map(tools),
        "predictor": t["Tool"], "stratum_label": t["Stratum (intronic offset, bp)"],
        "side": t["Side"], "tier": t["Tier"].str.lower(),
        "genes_reaching_the_tier": t["External genes reaching the tier, k/n (m)"],
        "note": t["Note"]})


def _metrics(tp, fp, fn, tn) -> dict:
    if any(pd.isna(x) for x in (tp, fp, fn, tn)):
        return {k: np.nan for k in ("precision", "recall", "f1", "mcc")}
    tp, fp, fn, tn = (int(x) for x in (tp, fp, fn, tn))
    den = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
    return {"precision": tp / (tp + fp) if tp + fp else np.nan,
            "recall": tp / (tp + fn) if tp + fn else np.nan,
            "f1": 2 * tp / (2 * tp + fp + fn) if 2 * tp + fp + fn else np.nan,
            "mcc": (tp * tn - fp * fn) / math.sqrt(den) if den else np.nan}


def c10() -> pd.DataFrame:
    p = K.read_back(PRACTICAL)
    p = p[p["threshold_source"].isin(["in_sample", "logo"])]
    rows = []
    for r in p.itertuples():
        if r.side == "pp3":
            cm = (r.n_pos_above, r.n_neg_above, r.n_pos_below, r.n_neg_below)
            pos, call = "damaging", "score at or above the threshold"
        else:
            cm = (r.bp4_n_neg_at_or_below, r.bp4_n_pos_at_or_below, r.bp4_n_neg_above,
                  r.bp4_n_pos_above)
            pos, call = "normal", "score at or below the threshold"
        rows.append({"tool": r.tool, "stratum": r.stratum, "side": r.side,
                     "threshold_kind": r.threshold_kind, "threshold_source": r.threshold_source,
                     "threshold": r.threshold, "status": r.status,
                     "positive_class": pos, "call": call,
                     "tp": cm[0], "fp": cm[1], "fn": cm[2], "tn": cm[3]} | _metrics(*cm))
    return pd.DataFrame(rows)


def main() -> None:
    df = K.load_set()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    t8 = c8(df)
    t8.to_csv(OUT_C8, index=False)
    t9 = pd.concat([c9a(), c9b()], ignore_index=True)
    t9.to_csv(OUT_C9, index=False)
    t10 = c10()
    with open(OUT_C10, "w", newline="") as f:
        f.write("".join(f"# {x}\n" for x in C10_HEADER.splitlines()))
        t10.to_csv(f, index=False)
    w = t8[(t8.part == "thresholds") & (t8.subset == "without BAP1")]
    print(f"[C8] without BAP1: {int((w.tier_reached != w.tier_reached_every_gene).sum())} of "
          f"{len(w)} cells change whether the tier is reached")
    a = t9[t9.part.str.startswith("(a)")]
    print(f"[C9] (a) {len(a)} held-out rows, {int(a.evaluable.sum())} evaluable, "
          f"{int((a.cleared_cut == True).sum())} clearing the cut")  # noqa: E712
    print(f"[C10] {len(t10)} rows, {int(t10.mcc.notna().sum())} with an MCC")


if __name__ == "__main__":
    main()
