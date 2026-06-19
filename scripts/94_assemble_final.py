"""M7: merge GPN-MSA + AlphaGenome into score_matrix_v2 -> score_matrix_final.tsv.
NT/Evo2/ESM are deferred (no GPU) -> added as NA columns. Builds the full
model x region coverage table and spot-checks known pathogenic splice variants.
Aligned on exact GRCh38 (chrom,pos,ref,alt) (ref verified by construction).
"""
import pandas as pd
from config import DATA_PROCESSED, DATA_RAW, RESULTS_TABLES

KEY = ["chrom", "pos", "ref", "alt"]
GPN = DATA_RAW / "scores" / "gpn" / "gpn.tsv"
AG = DATA_RAW / "scores" / "alphagenome" / "ag_cache.tsv"

m = pd.read_csv(DATA_PROCESSED / "score_matrix_v2.tsv", sep="\t",
                dtype={"chrom": str}, low_memory=False)
m["pos"] = m["pos"].astype(int)

# GPN-MSA (directionality: LOWER/more-negative = more deleterious)
gpn = pd.read_csv(GPN, sep="\t", dtype={"chrom": str})
gpn["pos"] = gpn["pos"].astype(int)
m = m.merge(gpn[KEY + ["gpn_msa_score"]], on=KEY, how="left")

# AlphaGenome (directionality: HIGHER = more splice-altering / pathogenic)
if AG.exists():
    ag = pd.read_csv(AG, sep="\t", dtype={"chrom": str})
    ag["pos"] = ag["pos"].astype(int)
    ag = ag.drop_duplicates(KEY)
    m = m.merge(ag[KEY + ["alphagenome_splice"]], on=KEY, how="left",
                suffixes=("", "_new"))
    if "alphagenome_splice_new" in m.columns:  # replace old NA placeholder
        m["alphagenome_splice"] = m["alphagenome_splice_new"]
        m = m.drop(columns="alphagenome_splice_new")
    print(f"AlphaGenome merged: {m['alphagenome_splice'].notna().sum()} scored")
else:
    print("AlphaGenome cache missing!")

# Deferred DNA-LMs (no GPU) -> explicit NA columns
for col in ["nucleotide_transformer", "evo2", "esm"]:
    m[col] = pd.NA

# directionality registry
DIRECTION = {
    "cadd_phred": "higher=pathogenic", "alphamissense": "higher=pathogenic",
    "phylop100way": "higher=pathogenic", "phastcons100way": "higher=pathogenic",
    "gnomad_af_global": "higher=benign", "gnomad_af_popmax": "higher=benign",
    "spliceai_ds": "higher=pathogenic", "pangolin_score": "higher=pathogenic",
    "gpn_msa_score": "LOWER=pathogenic", "alphagenome_splice": "higher=pathogenic",
    "nucleotide_transformer": "deferred", "evo2": "deferred", "esm": "deferred",
}
MODELS = list(DIRECTION)
for col in MODELS:
    if col not in m.columns:
        m[col] = pd.NA
    m[col + "_covered"] = m[col].notna()

m.to_csv(DATA_PROCESSED / "score_matrix_final.tsv", sep="\t", index=False)
print(f"Wrote score_matrix_final.tsv: {len(m)} rows, models: {MODELS}")

# coverage: model x region (scored / scored-low(|v|<0.01) / NA)
rows = []
for model in MODELS:
    for reg, g in m.groupby("region_class"):
        s = g[model]
        n_sc = int(s.notna().sum())
        n_low = int((s.notna() & (pd.to_numeric(s, errors="coerce").abs() < 0.01)).sum())
        rows.append({"model": model, "direction": DIRECTION[model],
                     "region_class": reg, "n_total": len(g),
                     "n_scored": n_sc, "n_scored_low_abs<0.01": n_low,
                     "n_NA": len(g) - n_sc})
cov = pd.DataFrame(rows)
cov.to_csv(RESULTS_TABLES / "full_model_coverage.tsv", sep="\t", index=False)

print("\n=== overall coverage per model (scored / total) ===")
for model in MODELS:
    print(f"  {model:24s} {int(m[model].notna().sum()):6d}/{len(m)}   {DIRECTION[model]}")

print("\n=== SPOT-CHECK known pathogenic splice (GPN-MSA low=pathogenic; AG high) ===")
chk = m[m.hgvs_nt.isin([
    "NM_007294.3:c.5194-1G>C", "NM_007294.3:c.4986+1G>A",
    "ENST00000256474.3:c.340+1G>T"]) |
    ((m.gene == "VHL") & (m.splice_class == "splice_core") & (m.clnsig_class == "P/LP"))]
print(chk[["gene", "hgvs_nt", "splice_class", "functional_score",
           "gpn_msa_score", "spliceai_ds", "alphagenome_splice", "clnsig_class"]]
      .head(10).to_string(index=False))
