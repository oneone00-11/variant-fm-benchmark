"""Eval 2 (AUXILIARY): AUROC + AUPRC vs ClinVar clean_PB (P/LP vs B/LB), pooled
across genes, stratified bootstrap 95% CI, in two versions (incl/excl BRCA1 for
the circularity sensitivity). Regions: all, splice, coding (splice_core skipped
- no benign). gnomAD-absent kept as missing (not AF=0).
"""
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score

from config import DATA_PROCESSED, RESULTS_TABLES
from phase3_lib import MODELS, oriented, region_mask, stratified_boot_indices

df = pd.read_csv(DATA_PROCESSED / "score_matrix_final.tsv", sep="\t",
                 dtype={"chrom": str}, low_memory=False)
df = df[df.label_set == "clean_PB"].copy()
df["y"] = df["path_binary"].astype(float)
for m in MODELS:
    df["_o_" + m] = oriented(df, m)

N_BOOT = 2000
REGIONS = ["all", "splice", "coding"]
rng = np.random.default_rng(20260619)
rows = []

for version, drop in [("incl_BRCA1", None), ("excl_BRCA1", "BRCA1")]:
    base = df if drop is None else df[df.gene != drop]
    for region in REGIONS:
        sub0 = base[region_mask(base, region)]
        for model in MODELS:
            sub = sub0[["_o_" + model, "y", "gene"]].dropna()
            y = sub["y"].values
            s = sub["_o_" + model].values
            npos, nneg = int((y == 1).sum()), int((y == 0).sum())
            if npos < 3 or nneg < 3:
                rows.append({"model": model, "region": region, "version": version,
                             "auroc": "", "auroc_lo": "", "auroc_hi": "",
                             "auprc": "", "auprc_lo": "", "auprc_hi": "",
                             "n_pos": npos, "n_neg": nneg, "note": "too few"})
                continue
            auroc = roc_auc_score(y, s)
            auprc = average_precision_score(y, s)
            genes = sub["gene"].values
            br, bp = [], []
            for _ in range(N_BOOT):
                idx = stratified_boot_indices(genes, rng)
                yb, sb = y[idx], s[idx]
                if len(np.unique(yb)) < 2:
                    continue
                br.append(roc_auc_score(yb, sb))
                bp.append(average_precision_score(yb, sb))
            ro = (np.percentile(br, 2.5), np.percentile(br, 97.5)) if br else (np.nan, np.nan)
            pr = (np.percentile(bp, 2.5), np.percentile(bp, 97.5)) if bp else (np.nan, np.nan)
            rows.append({"model": model, "region": region, "version": version,
                         "auroc": round(auroc, 4), "auroc_lo": round(ro[0], 4),
                         "auroc_hi": round(ro[1], 4), "auprc": round(auprc, 4),
                         "auprc_lo": round(pr[0], 4), "auprc_hi": round(pr[1], 4),
                         "n_pos": npos, "n_neg": nneg, "note": ""})

res = pd.DataFrame(rows)
res.to_csv(RESULTS_TABLES / "auroc_clinvar.tsv", sep="\t", index=False)

pd.set_option("display.width", 180)
print("=== AUROC vs ClinVar clean_PB (excl_BRCA1 = main; region=all) ===")
a = res[(res.version == "excl_BRCA1") & (res.region == "all")]
for r in a.itertuples():
    if r.auroc != "":
        print(f"  {r.model:20s} AUROC={r.auroc} [{r.auroc_lo},{r.auroc_hi}]  "
              f"AUPRC={r.auprc}  n_pos={r.n_pos} n_neg={r.n_neg}")
print("\n=== splice region, both versions (circularity sensitivity) ===")
for r in res[res.region == "splice"].itertuples():
    if r.auroc != "":
        print(f"  {r.model:20s} {r.version:11s} AUROC={r.auroc} [{r.auroc_lo},{r.auroc_hi}] "
              f"n_pos={r.n_pos} n_neg={r.n_neg}")
