"""M3-B: sample-size matrix + a-priori precision estimates.

Spearman 95% CI half-width: Fisher z, SE_z = 1/sqrt(N-3), back-transformed,
  at assumed rho = 0.3 / 0.5 / 0.7.
AUROC 95% CI half-width: Hanley-McNeil SE with the subset's actual clean_PB
  class balance, at assumed AUROC = 0.75 / 0.85 / 0.90.

These are PRIOR precision estimates under ASSUMED effect sizes, NOT real results
(real values require model scoring, which this milestone does not do).

Denominators (documented):
  - Spearman (continuous, predictor vs functional_score): needs only a functional
    score -> primary N = full functional-SNV set in the region (N_func). We also
    report N on the merged analysis-ready set (N_merged).
  - AUROC (binary, predictor vs ClinVar clean_PB): needs ClinVar P/B labels ->
    n_pos = P/LP, n_neg = B/LB from the merged set.
"""
import numpy as np
import pandas as pd

from config import DATA_PROCESSED, RESULTS_TABLES

FINAL_GENES = ["BRCA1", "BRCA2", "BARD1", "PALB2", "RAD51C"]

ar = pd.read_csv(DATA_PROCESSED / "analysis_ready.tsv", sep="\t",
                 dtype={"chrom": str}, low_memory=False)
fm = pd.read_csv(DATA_PROCESSED / "functional_scores_master.tsv", sep="\t",
                 dtype={"chrom": str}, low_memory=False)
fm = fm[(fm.is_snv == True) & (fm.gene.isin(FINAL_GENES))].copy()


def func_gap_region(off):
    if pd.isna(off):
        return None, None
    a = abs(off)
    if a <= 2:
        return "splice", "splice_core"
    if a <= 8:
        return "splice", "splice_region"
    return "intronic", None


fm[["fregion", "fsplice"]] = fm["intron_offset"].apply(
    lambda o: pd.Series(func_gap_region(o)))

# --- subset definitions (mask functions on a dataframe) ---
SUBSETS = {
    "whole_gap":     lambda d: d.region_class.isin(["splice", "intronic"]),
    "splice":        lambda d: d.region_class == "splice",
    "splice_core":   lambda d: d.splice_class == "splice_core",
    "splice_region": lambda d: d.splice_class == "splice_region",
    "intronic":      lambda d: d.region_class == "intronic",
    "coding":        lambda d: d.region_class.isin(["missense", "synonymous", "nonsense"]),
    "missense":      lambda d: d.region_class == "missense",
    "utr":           lambda d: d.region_class == "utr",
}
# functional-only ceiling masks (gap subsets only)
FSUBSETS = {
    "whole_gap":     lambda d: d.fregion.isin(["splice", "intronic"]),
    "splice":        lambda d: d.fregion == "splice",
    "splice_core":   lambda d: d.fsplice == "splice_core",
    "splice_region": lambda d: d.fsplice == "splice_region",
    "intronic":      lambda d: d.fregion == "intronic",
}


def spearman_hw(N, rho):
    if N is None or N <= 4:
        return np.nan
    se = 1.0 / np.sqrt(N - 3)
    z = np.arctanh(rho)
    lo, hi = np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)
    return (hi - lo) / 2.0


def auroc_hw(auc, n_pos, n_neg):
    if n_pos < 1 or n_neg < 1:
        return np.nan
    q1 = auc / (2 - auc)
    q2 = 2 * auc ** 2 / (1 + auc)
    var = (auc * (1 - auc) + (n_pos - 1) * (q1 - auc ** 2)
           + (n_neg - 1) * (q2 - auc ** 2)) / (n_pos * n_neg)
    return 1.96 * np.sqrt(max(var, 0))


def verdict_sp(hw):
    if np.isnan(hw):
        return "N_too_small"
    return "sufficient" if hw <= 0.10 else ("marginal" if hw <= 0.15 else "insufficient")


def verdict_au(hw):
    if np.isnan(hw):
        return "no_PB"
    return "sufficient" if hw <= 0.05 else ("marginal" if hw <= 0.08 else "insufficient")


groups = [(g, ar[ar.gene == g], fm[fm.gene == g]) for g in FINAL_GENES]
groups.append(("POOLED", ar[ar.gene.isin(FINAL_GENES)], fm))

ss_rows, pr_rows = [], []
for gene, ard, fmd in groups:
    for sub, mask in SUBSETS.items():
        d = ard[mask(ard)]
        N_merged = len(d)
        n_plp = int((d.clnsig_class == "P/LP").sum())
        n_blb = int((d.clnsig_class == "B/LB").sum())
        n_vus = int((d.label_set == "vus").sum())
        n_conf = int((d.label_set == "conflicting").sum())
        n_other = int((d.label_set == "other").sum())
        N_func = int(FSUBSETS[sub](fmd).sum()) if sub in FSUBSETS else np.nan
        ss_rows.append(dict(gene=gene, region_subset=sub, N_merged=N_merged,
                            n_PLP=n_plp, n_BLB=n_blb, n_clean_PB=n_plp + n_blb,
                            n_VUS=n_vus, n_conflicting=n_conf, n_other=n_other,
                            N_func_total=N_func))
        # precision: prefer functional N for Spearman where available
        N_sp = N_func if (sub in FSUBSETS and not np.isnan(N_func) and N_func > 0) else N_merged
        pr_rows.append({
            "gene": gene, "region_subset": sub, "N_spearman": int(N_sp),
            "N_merged": N_merged,
            "sp_hw_rho0.3": round(spearman_hw(N_sp, 0.3), 4),
            "sp_hw_rho0.5": round(spearman_hw(N_sp, 0.5), 4),
            "sp_hw_rho0.7": round(spearman_hw(N_sp, 0.7), 4),
            "n_pos": n_plp, "n_neg": n_blb,
            "auroc_hw_0.75": round(auroc_hw(0.75, n_plp, n_blb), 4),
            "auroc_hw_0.85": round(auroc_hw(0.85, n_plp, n_blb), 4),
            "auroc_hw_0.90": round(auroc_hw(0.90, n_plp, n_blb), 4),
            "verdict_spearman": verdict_sp(spearman_hw(N_sp, 0.5)),
            "verdict_auroc": verdict_au(auroc_hw(0.85, n_plp, n_blb)),
        })

ss = pd.DataFrame(ss_rows)
pr = pd.DataFrame(pr_rows)
ss.to_csv(RESULTS_TABLES / "sample_size_matrix.tsv", sep="\t", index=False)
pr.to_csv(RESULTS_TABLES / "precision_estimate.tsv", sep="\t", index=False)

pd.set_option("display.width", 220); pd.set_option("display.max_columns", 40)
print("=== sample_size_matrix (clean_PB = P/LP vs B/LB; VUS separate) ===")
print(ss.to_string(index=False))
print("\n=== precision_estimate (PRIOR, assumed effect sizes — not real results) ===")
show = ["gene", "region_subset", "N_spearman", "sp_hw_rho0.5", "n_pos", "n_neg",
        "auroc_hw_0.85", "verdict_spearman", "verdict_auroc"]
print(pr[show].to_string(index=False))
