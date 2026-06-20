"""Phase 3 shared utilities: orientation, Fisher-z CIs, DerSimonian-Laird
random-effects meta-analysis, and stratified bootstrap.

Orientation: every model oriented so HIGHER = more pathogenic.
Functional pathogenicity = -functional_score (SGE: more negative = more damaging),
so a POSITIVE Spearman(model_oriented, func_path) = model agrees with pathogenicity.
"""
import numpy as np
import pandas as pd
from scipy import stats

# +1 keep, -1 negate, so higher = more pathogenic
ORIENT = {
    "cadd_phred": 1, "alphamissense": 1, "phylop100way": 1, "phastcons100way": 1,
    "spliceai_ds": 1, "pangolin_score": 1, "alphagenome_splice": 1,
    "gpn_msa_score": -1,        # lower = pathogenic -> negate
    "gnomad_af_global": -1, "gnomad_af_popmax": -1,  # higher = benign -> negate
    "nucleotide_transformer": 1,  # M8b: higher = pathogenic (already aligned)
}
MODELS = list(ORIENT)
# DNA/splice-reading vs protein/missense (for the banner comparison)
DNA_SPLICE = ["spliceai_ds", "pangolin_score", "gpn_msa_score", "alphagenome_splice"]
DNA_LM = ["gpn_msa_score", "nucleotide_transformer"]  # the two DNA language models
PROTEIN = ["alphamissense"]

GENES = ["BRCA1", "BRCA2", "BARD1", "PALB2", "RAD51C", "VHL", "BAP1"]


def oriented(df, model):
    """Return model score oriented higher=pathogenic (NaN where uncovered)."""
    return ORIENT[model] * pd.to_numeric(df[model], errors="coerce")


def func_path(df):
    """Functional pathogenicity, higher = more pathogenic."""
    return -pd.to_numeric(df["functional_score"], errors="coerce")


def region_mask(df, subset):
    rc = df["region_class"]
    if subset == "all":
        return pd.Series(True, index=df.index)
    if subset == "splice":
        return rc == "splice"
    if subset == "splice_core":
        return df["splice_class"] == "splice_core"
    if subset == "splice_region":
        return df["splice_class"] == "splice_region"
    if subset == "intronic":
        return rc == "intronic"
    if subset == "coding":
        return rc.isin(["missense", "synonymous", "nonsense"])
    if subset == "missense":
        return rc == "missense"
    raise ValueError(subset)


REGION_SUBSETS = ["all", "splice", "splice_core", "splice_region",
                  "intronic", "coding", "missense"]


def spearman_ci(x, y):
    """Spearman rho + n + 95% Fisher-z CI. Returns (rho, n, lo, hi)."""
    m = np.isfinite(x) & np.isfinite(y)
    x, y = np.asarray(x)[m], np.asarray(y)[m]
    n = len(x)
    if n < 4 or np.std(x) == 0 or np.std(y) == 0:
        return (np.nan, n, np.nan, np.nan)
    rho = stats.spearmanr(x, y).correlation
    if not np.isfinite(rho):
        return (np.nan, n, np.nan, np.nan)
    z = np.arctanh(np.clip(rho, -0.999999, 0.999999))
    se = 1.0 / np.sqrt(n - 3)
    return (rho, n, np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se))


def dl_meta(rhos, ns):
    """DerSimonian-Laird random-effects meta over Fisher-z of per-gene rho.
    Returns dict: pooled rho, CI, I2, tau2, k."""
    rhos = np.asarray(rhos, float)
    ns = np.asarray(ns, float)
    ok = np.isfinite(rhos) & (ns > 3)
    rhos, ns = rhos[ok], ns[ok]
    k = len(rhos)
    if k == 0:
        return dict(rho=np.nan, lo=np.nan, hi=np.nan, I2=np.nan, tau2=np.nan, k=0)
    z = np.arctanh(np.clip(rhos, -0.999999, 0.999999))
    v = 1.0 / (ns - 3)
    w = 1.0 / v
    z_fixed = np.sum(w * z) / np.sum(w)
    Q = np.sum(w * (z - z_fixed) ** 2)
    df = k - 1
    C = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0.0, (Q - df) / C) if C > 0 and k > 1 else 0.0
    wster = 1.0 / (v + tau2)
    z_re = np.sum(wster * z) / np.sum(wster)
    se_re = np.sqrt(1.0 / np.sum(wster))
    I2 = max(0.0, (Q - df) / Q) * 100 if Q > 0 and k > 1 else 0.0
    return dict(rho=np.tanh(z_re), lo=np.tanh(z_re - 1.96 * se_re),
                hi=np.tanh(z_re + 1.96 * se_re), I2=I2, tau2=tau2, k=k)


def stratified_boot_indices(genes, rng):
    """Resample row indices with replacement WITHIN each gene (stratified)."""
    idx = []
    g = pd.Series(genes).reset_index(drop=True)
    for gene, sub in g.groupby(g):
        ii = sub.index.values
        idx.append(rng.choice(ii, size=len(ii), replace=True))
    return np.concatenate(idx)
