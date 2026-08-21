"""Phase 3 figures (simple, reproducible):
 1) splice-subset model performance (meta rho + CI) with coverage-gap annotation
 2) forest plot: per-gene rho + meta for the splice subset (SpliceAI)
 3) region-stratified heatmap of meta rho (model x region)
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from config import DATA_PROCESSED, RESULTS_TABLES, RESULTS_FIGURES
from phase3_lib import MODELS, GENES, DNA_SPLICE, region_mask

# Publication settings: 300 dpi, opaque white background (no alpha).
DPI = 300
SAVE_KW = dict(dpi=DPI, facecolor="white")

# Clean display names for predictors (no raw column keys in figures).
NAME = {
    "cadd_phred": "CADD", "alphamissense": "AlphaMissense",
    "phylop100way": "phyloP", "phastcons100way": "phastCons",
    "spliceai_ds": "SpliceAI", "pangolin_score": "Pangolin",
    "gpn_msa_score": "GPN-MSA", "alphagenome_splice": "AlphaGenome",
    "gnomad_af_global": "gnomAD AF (global)", "gnomad_af_popmax": "gnomAD AF (popmax)",
    "nucleotide_transformer": "Nucleotide Transformer",
}
def disp(m):
    return NAME.get(m, m)

meta = pd.read_csv(RESULTS_TABLES / "spearman_meta.tsv", sep="\t")
bygene = pd.read_csv(RESULTS_TABLES / "spearman_by_gene.tsv", sep="\t")
sm = pd.read_csv(DATA_PROCESSED / "score_matrix_final.tsv", sep="\t",
                 dtype={"chrom": str}, low_memory=False)


def num(x):
    return pd.to_numeric(x, errors="coerce")


# ---- Fig 1: splice subset performance + coverage gap ----
sp = meta[meta.region == "splice"].copy()
sp["meta_rho"] = num(sp["meta_rho"]); sp = sp.dropna(subset=["meta_rho"])
sp = sp.sort_values("meta_rho", ascending=True)
splice_n = region_mask(sm, "splice").sum()
covn = {m: int(sm.loc[region_mask(sm, "splice"), m].notna().sum()) for m in MODELS}

fig, ax = plt.subplots(figsize=(8, 5))
colors = ["#d62728" if m in DNA_SPLICE else ("#7f7f7f" if m == "alphamissense"
          else "#1f77b4") for m in sp.model]
lo = num(sp.ci_lo).values; hi = num(sp.ci_hi).values
rho = sp.meta_rho.values
ax.barh(range(len(sp)), rho, color=colors,
        xerr=[rho - lo, hi - rho], capsize=3)
ax.set_yticks(range(len(sp)))
ax.set_yticklabels([f"{disp(m)}  (n={covn.get(m,0)})" for m in sp.model])
ax.set_xlabel("meta Spearman rho vs functional gold standard (splice subset)")
ax.set_title(f"Splice variants (N={splice_n}): DNA/splice (red) vs protein "
             f"AlphaMissense (grey, n={covn.get('alphamissense',0)})")
ax.axvline(0, color="k", lw=0.5)
plt.tight_layout()
plt.savefig(RESULTS_FIGURES / "fig1_splice_performance.png", **SAVE_KW)
plt.close()

# ---- Fig 2: forest plot, SpliceAI splice, per gene + meta ----
fg = bygene[(bygene.model == "spliceai_ds") & (bygene.region == "splice")].copy()
fg["rho"] = num(fg["rho"]); fg = fg.dropna(subset=["rho"])
fig, ax = plt.subplots(figsize=(7, 5))
y = range(len(fg))
ax.errorbar(fg.rho, y, xerr=[fg.rho - num(fg.ci_lo), num(fg.ci_hi) - fg.rho],
            fmt="o", color="#1f77b4", capsize=3)
ax.set_yticks(list(y)); ax.set_yticklabels([f"{g} (n={n})" for g, n in zip(fg.gene, fg.n)])
mrow = meta[(meta.model == "spliceai_ds") & (meta.region == "splice")].iloc[0]
mr = float(mrow.meta_rho)
ax.axvline(mr, color="#d62728", ls="--", label=f"meta rho={mr:.3f}")
ax.set_xlabel("Spearman rho (SpliceAI vs functional, splice)")
ax.set_title("Per-gene rho + random-effects meta (SpliceAI, splice)")
ax.legend()
plt.tight_layout()
plt.savefig(RESULTS_FIGURES / "fig2_forest_spliceai_splice.png", **SAVE_KW)
plt.close()

# ---- Fig 3: region heatmap of meta rho ----
regions = ["all", "splice", "intronic", "coding", "missense"]
mat = meta.copy(); mat["meta_rho"] = num(mat["meta_rho"])
H = mat.pivot(index="model", columns="region", values="meta_rho").reindex(
    index=MODELS, columns=regions)
fig, ax = plt.subplots(figsize=(7, 6))
im = ax.imshow(H.values, cmap="RdBu_r", vmin=-0.4, vmax=0.8, aspect="auto")
ax.set_xticks(range(len(regions))); ax.set_xticklabels(regions, rotation=30)
ax.set_yticks(range(len(MODELS))); ax.set_yticklabels([disp(m) for m in MODELS])
for i in range(len(MODELS)):
    for j in range(len(regions)):
        v = H.values[i, j]
        if np.isfinite(v):
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7,
                    color="white" if abs(v) > 0.5 else "black")
fig.colorbar(im, label="meta Spearman rho")
ax.set_title("Meta rho vs functional gold standard (model x region)")
plt.tight_layout()
plt.savefig(RESULTS_FIGURES / "fig3_region_heatmap.png", **SAVE_KW)
plt.close()

print("wrote fig1_splice_performance.png, fig2_forest_spliceai_splice.png, "
      "fig3_region_heatmap.png")
