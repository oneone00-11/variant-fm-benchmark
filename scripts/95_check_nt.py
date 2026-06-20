"""Verify the NT-augmented score matrix before re-running Phase 3."""
import pandas as pd
from pathlib import Path

SRC = Path(r"C:\Users\张宁一\Desktop\variant-fm-benchmark\score_matrix_final_with_nt.tsv")
df = pd.read_csv(SRC, sep="\t", dtype={"chrom": str}, low_memory=False)
print("rows:", len(df))
print("has nucleotide_transformer col:", "nucleotide_transformer" in df.columns)
nt = pd.to_numeric(df["nucleotide_transformer"], errors="coerce")
print(f"NT non-NA: {nt.notna().sum()} / {len(df)}")
print(f"NT range: {nt.min():.4f} .. {nt.max():.4f}")
# original key cols present?
need = ["gene", "region_class", "splice_class", "functional_score", "clnsig_class",
        "path_binary", "label_set", "circularity_flag", "spliceai_ds", "gpn_msa_score",
        "alphagenome_splice", "alphamissense", "cadd_phred"]
miss = [c for c in need if c not in df.columns]
print("missing expected cols:", miss if miss else "none")
# direction spot-check: pathogenic splice (high NT?) vs benign deep-intronic
print("\nNT by clnsig on splice (mean):")
sp = df[df.region_class == "splice"]
print(sp.groupby("clnsig_class")["nucleotide_transformer"].apply(
    lambda s: round(pd.to_numeric(s, errors="coerce").mean(), 4)).to_string())
print("\nknown pathogenic splice NT scores:")
chk = df[df.hgvs_nt.isin(["NM_007294.3:c.5194-1G>C", "NM_007294.3:c.4986+1G>A",
                          "ENST00000256474.3:c.340+1G>T"])]
print(chk[["gene", "hgvs_nt", "nucleotide_transformer", "spliceai_ds",
           "gpn_msa_score", "functional_score", "clnsig_class"]].to_string(index=False))
