"""Inspect both master tables before the M3 clean/merge: clnvc distribution,
mc_terms / consequence format, functional columns, intron_offset, is_snv."""
import pandas as pd
from collections import Counter
from config import DATA_PROCESSED, TARGET_GENES

cv = pd.read_csv(DATA_PROCESSED / "clinvar_target_variants.tsv", sep="\t",
                 dtype={"chrom": str}, low_memory=False)
fm = pd.read_csv(DATA_PROCESSED / "functional_scores_master.tsv", sep="\t",
                 dtype={"chrom": str}, low_memory=False)

print("=== ClinVar master ===")
print("cols:", list(cv.columns))
print("rows:", len(cv))
print("\nclnvc value counts:")
print(cv.clnvc.value_counts().to_string())
print("\nsample mc_terms (10):")
print(cv[cv.mc_terms.notna()].mc_terms.head(10).tolist())
print("\nconsequence_class counts:")
print(cv.consequence_class.value_counts().to_string())

# all distinct SO terms across mc_terms
so = Counter()
for s in cv.mc_terms.dropna():
    for t in str(s).split(";"):
        if t:
            so[t] += 1
print("\nAll SO terms in mc_terms (sorted by freq):")
for t, n in so.most_common():
    print(f"  {t}: {n}")

print("\n=== Functional master ===")
print("cols:", list(fm.columns))
print("rows:", len(fm), "| is_snv True:", int(fm.is_snv.sum()))
print("region_class counts:")
print(fm.region_class.value_counts().to_string())
print("\nintron_offset non-null:", fm.intron_offset.notna().sum(),
      "| range:", fm.intron_offset.min(), "to", fm.intron_offset.max())
print("functional_call non-empty:",
      (fm.functional_call.astype(str).str.strip() != "").sum())
print("assay_type counts:", fm.assay_type.value_counts().to_dict())
print("source_urn counts:", fm.source_urn.value_counts().to_dict())
