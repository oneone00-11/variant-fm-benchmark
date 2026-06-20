"""Eval 3 (clinical, exploratory) + pairwise splice comparison.
- pairwise_splice.tsv: paired bootstrap of Delta(mean per-gene rho) between
  DNA/splice models on the splice subset (same variants).
- vus_reclass.tsv: exploratory - per-gene functional threshold (Youden J on
  clean_PB) labels VUS functionally abnormal/normal; each model's AUROC to
  separate them. Small n -> CI reported, exploratory only.
- sens_at_95spec.tsv: sensitivity at 95% specificity, clean_PB excl BRCA1.
"""
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score, roc_curve

from config import DATA_PROCESSED, RESULTS_TABLES
from phase3_lib import (MODELS, GENES, DNA_SPLICE, oriented, func_path,
                        region_mask, stratified_boot_indices)

df = pd.read_csv(DATA_PROCESSED / "score_matrix_final.tsv", sep="\t",
                 dtype={"chrom": str}, low_memory=False)
df["fp"] = func_path(df)
for m in MODELS:
    df["_o_" + m] = oriented(df, m)
rng = np.random.default_rng(20260619)
N_BOOT = 2000

# ---------- pairwise splice: Delta(mean per-gene rho) ----------
PAIR_MODELS = DNA_SPLICE + ["nucleotide_transformer", "cadd_phred"]
sp = df[region_mask(df, "splice")].copy()


def mean_gene_rho(frame, model, idx=None):
    f = frame if idx is None else frame.iloc[idx]
    rs = []
    for g in GENES:
        gg = f[f.gene == g]
        x = gg["_o_" + model].values
        y = gg["fp"].values
        m = np.isfinite(x) & np.isfinite(y)
        if m.sum() >= 4 and np.std(x[m]) > 0:
            r = stats.spearmanr(x[m], y[m]).correlation
            if np.isfinite(r):
                rs.append(r)
    return np.mean(rs) if rs else np.nan


sp_idx_genes = sp["gene"].values
pair_rows = []
point = {m: mean_gene_rho(sp, m) for m in PAIR_MODELS}
for i, a in enumerate(PAIR_MODELS):
    for b in PAIR_MODELS[i + 1:]:
        diffs = []
        for _ in range(N_BOOT):
            idx = stratified_boot_indices(sp_idx_genes, rng)
            diffs.append(mean_gene_rho(sp, a, idx) - mean_gene_rho(sp, b, idx))
        diffs = [d for d in diffs if np.isfinite(d)]
        lo, hi = np.percentile(diffs, [2.5, 97.5])
        pair_rows.append({"model_a": a, "model_b": b,
                          "rho_a": round(point[a], 4), "rho_b": round(point[b], 4),
                          "delta": round(point[a] - point[b], 4),
                          "ci_lo": round(lo, 4), "ci_hi": round(hi, 4),
                          "sig": "yes" if (lo > 0 or hi < 0) else "no"})
pd.DataFrame(pair_rows).to_csv(RESULTS_TABLES / "pairwise_splice.tsv", sep="\t", index=False)
print("=== pairwise splice Delta(mean per-gene rho) ===")
for r in pair_rows:
    print(f"  {r['model_a']:18s} - {r['model_b']:18s} d={r['delta']:+.3f} "
          f"[{r['ci_lo']:+.3f},{r['ci_hi']:+.3f}] sig={r['sig']}")

# ---------- VUS reclassification (exploratory) ----------
# per-gene functional threshold from clean_PB via Youden J (damaging = low fp? fp
# is higher=pathogenic, so threshold on fp; abnormal = fp >= thr)
pb = df[df.label_set == "clean_PB"]
vus = df[df.label_set == "vus"].copy()
thr = {}
for g in GENES:
    c = pb[pb.gene == g]
    y = c["path_binary"].values.astype(float)
    f = c["fp"].values
    m = np.isfinite(f)
    if m.sum() < 10 or len(np.unique(y[m])) < 2:
        continue
    fpr, tpr, t = roc_curve(y[m], f[m])
    thr[g] = t[np.argmax(tpr - fpr)]
vus["func_abnormal"] = np.nan
for g, t in thr.items():
    sel = vus.gene == g
    vus.loc[sel, "func_abnormal"] = (vus.loc[sel, "fp"] >= t).astype(float)
vrows = []
vv = vus[vus.func_abnormal.notna()]
for model in MODELS:
    s = vv["_o_" + model].values
    y = vv["func_abnormal"].values
    m = np.isfinite(s) & np.isfinite(y)
    npos, nneg = int((y[m] == 1).sum()), int((y[m] == 0).sum())
    if npos < 5 or nneg < 5:
        vrows.append({"model": model, "auroc": "", "ci_lo": "", "ci_hi": "",
                      "n_abnormal": npos, "n_normal": nneg, "note": "too few"})
        continue
    auroc = roc_auc_score(y[m], s[m])
    genes = vv["gene"].values[m]
    yb_, sb_ = y[m], s[m]
    boots = []
    for _ in range(N_BOOT):
        idx = stratified_boot_indices(genes, rng)
        if len(np.unique(yb_[idx])) < 2:
            continue
        boots.append(roc_auc_score(yb_[idx], sb_[idx]))
    lo, hi = (np.percentile(boots, [2.5, 97.5]) if boots else (np.nan, np.nan))
    vrows.append({"model": model, "auroc": round(auroc, 4),
                  "ci_lo": round(lo, 4), "ci_hi": round(hi, 4),
                  "n_abnormal": npos, "n_normal": nneg, "note": "EXPLORATORY"})
pd.DataFrame(vrows).to_csv(RESULTS_TABLES / "vus_reclass.tsv", sep="\t", index=False)
print(f"\n=== VUS reclassification (EXPLORATORY; n VUS with func label) ===")
for r in vrows:
    if r["auroc"] != "":
        print(f"  {r['model']:20s} AUROC={r['auroc']} [{r['ci_lo']},{r['ci_hi']}] "
              f"n_abn={r['n_abnormal']} n_norm={r['n_normal']}")

# ---------- sensitivity at 95% specificity (clean_PB excl BRCA1) ----------
pbx = df[(df.label_set == "clean_PB") & (df.gene != "BRCA1")]
srows = []
for model in MODELS:
    sub = pbx[["_o_" + model, "path_binary", "gene"]].dropna()
    y = sub["path_binary"].values.astype(float)
    s = sub["_o_" + model].values
    if (y == 1).sum() < 10 or (y == 0).sum() < 10:
        srows.append({"model": model, "sens_at_95spec": "", "ci_lo": "", "ci_hi": "",
                      "n_pos": int((y == 1).sum()), "n_neg": int((y == 0).sum())})
        continue

    def sens95(yy, ss):
        thr = np.quantile(ss[yy == 0], 0.95)  # 95% specificity
        return (ss[yy == 1] >= thr).mean()

    point_s = sens95(y, s)
    genes = sub["gene"].values
    boots = []
    for _ in range(N_BOOT):
        idx = stratified_boot_indices(genes, rng)
        if (y[idx] == 0).sum() > 5 and (y[idx] == 1).sum() > 5:
            boots.append(sens95(y[idx], s[idx]))
    lo, hi = (np.percentile(boots, [2.5, 97.5]) if boots else (np.nan, np.nan))
    srows.append({"model": model, "sens_at_95spec": round(point_s, 4),
                  "ci_lo": round(lo, 4), "ci_hi": round(hi, 4),
                  "n_pos": int((y == 1).sum()), "n_neg": int((y == 0).sum())})
pd.DataFrame(srows).to_csv(RESULTS_TABLES / "sens_at_95spec.tsv", sep="\t", index=False)
print("\n=== sensitivity @ 95% specificity (clean_PB, excl BRCA1) ===")
for r in srows:
    if r["sens_at_95spec"] != "":
        print(f"  {r['model']:20s} sens={r['sens_at_95spec']} [{r['ci_lo']},{r['ci_hi']}] "
              f"n_pos={r['n_pos']} n_neg={r['n_neg']}")
