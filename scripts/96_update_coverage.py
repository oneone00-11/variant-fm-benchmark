"""Regenerate full_model_coverage.tsv from the NT-augmented score_matrix_final.
(Does NOT rebuild the matrix; just recomputes coverage incl. NT.)"""
import pandas as pd
from config import DATA_PROCESSED, RESULTS_TABLES
from phase3_lib import MODELS, REGION_SUBSETS

df = pd.read_csv(DATA_PROCESSED / "score_matrix_final.tsv", sep="\t",
                 dtype={"chrom": str}, low_memory=False)
ALL = MODELS + ["evo2", "esm"]  # evo2/esm still NA (pending)
rows = []
for model in ALL:
    col = df[model] if model in df.columns else pd.Series([pd.NA] * len(df))
    s = pd.to_numeric(col, errors="coerce")
    for reg in REGION_SUBSETS:
        from phase3_lib import region_mask
        g = s[region_mask(df, reg)]
        n_sc = int(g.notna().sum())
        n_low = int((g.notna() & (g.abs() < 0.01)).sum())
        rows.append({"model": model, "region": reg, "n_total": int(region_mask(df, reg).sum()),
                     "n_scored": n_sc, "n_scored_low_abs<0.01": n_low,
                     "n_NA": int(region_mask(df, reg).sum()) - n_sc})
cov = pd.DataFrame(rows)
cov.to_csv(RESULTS_TABLES / "full_model_coverage.tsv", sep="\t", index=False)
print("=== coverage (region=all) ===")
a = cov[cov.region == "all"]
for r in a.itertuples():
    print(f"  {r.model:24s} scored {r.n_scored:6d}/{r.n_total}")
