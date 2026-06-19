"""Eval 1 (MAIN): Spearman rho of each model vs the continuous functional gold
standard, per gene (each gene's own functional scale), then DL random-effects
meta across genes. Stratified by region. Positive rho = model agrees with
functional pathogenicity.
"""
import pandas as pd

from config import DATA_PROCESSED, RESULTS_TABLES
from phase3_lib import (MODELS, GENES, REGION_SUBSETS, oriented, func_path,
                        region_mask, spearman_ci, dl_meta)

df = pd.read_csv(DATA_PROCESSED / "score_matrix_final.tsv", sep="\t",
                 dtype={"chrom": str}, low_memory=False)
df["fp"] = func_path(df)
for m in MODELS:
    df["_o_" + m] = oriented(df, m)

by_gene_rows, meta_rows = [], []
for model in MODELS:
    for subset in REGION_SUBSETS:
        sub = df[region_mask(df, subset)]
        rhos, ns = [], []
        for gene in GENES:
            g = sub[sub.gene == gene]
            rho, n, lo, hi = spearman_ci(g["_o_" + model].values, g["fp"].values)
            by_gene_rows.append({
                "model": model, "gene": gene, "region": subset,
                "rho": round(rho, 4) if pd.notna(rho) else "",
                "n": n, "ci_lo": round(lo, 4) if pd.notna(lo) else "",
                "ci_hi": round(hi, 4) if pd.notna(hi) else "",
            })
            rhos.append(rho); ns.append(n)
        meta = dl_meta(rhos, ns)
        meta_rows.append({
            "model": model, "region": subset,
            "meta_rho": round(meta["rho"], 4) if pd.notna(meta["rho"]) else "",
            "ci_lo": round(meta["lo"], 4) if pd.notna(meta["lo"]) else "",
            "ci_hi": round(meta["hi"], 4) if pd.notna(meta["hi"]) else "",
            "I2_pct": round(meta["I2"], 1), "tau2": round(meta["tau2"], 4),
            "k_genes": meta["k"],
            "total_n": int(sub[["_o_" + model, "fp"]].dropna().shape[0]),
        })

pd.DataFrame(by_gene_rows).to_csv(RESULTS_TABLES / "spearman_by_gene.tsv",
                                  sep="\t", index=False)
meta = pd.DataFrame(meta_rows)
meta.to_csv(RESULTS_TABLES / "spearman_meta.tsv", sep="\t", index=False)

pd.set_option("display.width", 170)
print("=== meta Spearman rho (model x region), main result ===")
piv = meta.pivot(index="model", columns="region", values="meta_rho")
print(piv[["all", "splice", "intronic", "coding", "missense"]].to_string())
print("\n=== SPLICE subset: meta rho [CI], k genes, total n ===")
sp = meta[meta.region == "splice"].copy()
for r in sp.itertuples():
    print(f"  {r.model:20s} rho={r.meta_rho}  CI[{r.ci_lo},{r.ci_hi}]  "
          f"I2={r.I2_pct}%  k={r.k_genes}  n={r.total_n}")
print("\n=== AlphaMissense coverage in splice (the gap) ===")
sm = df[region_mask(df, "splice")]
print(f"  splice variants: {len(sm)}; AlphaMissense non-NA: "
      f"{sm['alphamissense'].notna().sum()}")
