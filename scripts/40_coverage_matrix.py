"""Intersect functional SNVs with the ClinVar target SNV set and build the
gold-standard coverage matrix (Milestone 2, subtask 4 / deliverable 2).

coverage matrix: gene x region_class with
  (a) n_functional       -- functional-scored variants in that region
  (b) n_functional_snv   -- of those, with a usable GRCh38 SNV coordinate
  (c) n_overlap_clinvar  -- (b) that are also a ClinVar target SNV (any star)
  (d) n_overlap_ge1star  -- (c) restricted to ClinVar review stars >= 1
  (e) n_overlap_vus      -- (c) that are ClinVar VUS (the reclassification target)

Also writes a ClinVar-consequence view of the overlap (caliber-consistent with
Milestone 1's SO-based labels) for the missense/splice/noncoding split.
"""
import pandas as pd

from config import DATA_PROCESSED, RESULTS_TABLES, TARGET_GENES


def norm_chrom(x):
    s = str(x)
    return s[:-2] if s.endswith(".0") else s


def key_series(df):
    return (df["chrom"].map(norm_chrom) + ":" + df["pos"].astype("Int64").astype(str)
            + ":" + df["ref"].astype(str) + ":" + df["alt"].astype(str))


fm = pd.read_csv(DATA_PROCESSED / "functional_scores_master.tsv", sep="\t",
                 dtype={"chrom": str}, low_memory=False)
cv = pd.read_csv(DATA_PROCESSED / "clinvar_target_variants.tsv", sep="\t",
                 dtype={"chrom": str}, low_memory=False)
cv = cv[cv["clnvc"] == "single_nucleotide_variant"].copy()

# Keys
fm_snv = fm[fm["is_snv"] == True].copy()
fm_snv["key"] = key_series(fm_snv)
cv["key"] = key_series(cv)

cv_idx = cv.set_index("key")
cv_keys = set(cv.key)
cv_ge1 = set(cv[cv.stars >= 1].key)
cv_vus = set(cv[(cv.stars >= 1) & (cv.clnsig_class == "VUS")].key)

fm_snv["in_clinvar"] = fm_snv.key.isin(cv_keys)
fm_snv["in_clinvar_ge1"] = fm_snv.key.isin(cv_ge1)
fm_snv["in_clinvar_vus"] = fm_snv.key.isin(cv_vus)

# --- coverage matrix: gene x region_class ---
rows = []
for gene in TARGET_GENES:
    g_all = fm[fm.gene == gene]
    g_snv = fm_snv[fm_snv.gene == gene]
    regions = sorted(set(g_all.region_class))
    for reg in regions:
        a = (g_all.region_class == reg).sum()
        sub = g_snv[g_snv.region_class == reg]
        rows.append({
            "gene": gene, "region_class": reg,
            "n_functional": int(a),
            "n_functional_snv": int(len(sub)),
            "n_overlap_clinvar": int(sub.in_clinvar.sum()),
            "n_overlap_ge1star": int(sub.in_clinvar_ge1.sum()),
            "n_overlap_vus": int(sub.in_clinvar_vus.sum()),
        })
mat = pd.DataFrame(rows)
mat.to_csv(RESULTS_TABLES / "goldstandard_coverage_matrix.tsv", sep="\t", index=False)

pd.set_option("display.width", 160)
print("=== Gold-standard coverage matrix (gene x region) ===")
print(mat.to_string(index=False))

# --- headline: splice/non-coding coverage ---
gap = mat[mat.region_class.isin(["splice", "intronic", "noncoding"])]
print("\n=== Splice/intronic/non-coding (the project's gap) ===")
print(gap.groupby("gene")[["n_functional_snv", "n_overlap_clinvar",
                           "n_overlap_ge1star", "n_overlap_vus"]].sum().to_string())

# --- ClinVar-consequence view of the overlap (SO labels, Milestone-1 caliber) ---
ov = fm_snv[fm_snv.in_clinvar].copy()
ov = ov.merge(cv_idx[["consequence_class", "clnsig_class", "stars"]],
              left_on="key", right_index=True, how="left")
print("\n=== Overlap by ClinVar SO consequence_class (>=1 star) ===")
ov1 = ov[ov.stars >= 1]
print(ov1.pivot_table(index="consequence_class", columns="gene",
                      values="key", aggfunc="count", fill_value=0).to_string())

print("\n=== Overlap clnsig (>=1 star) per gene ===")
print(ov1.pivot_table(index="gene", columns="clnsig_class",
                      values="key", aggfunc="count", fill_value=0)
      .reindex(TARGET_GENES).fillna(0).astype(int).to_string())

# totals
print(f"\nTotal functional SNVs: {len(fm_snv)}")
print(f"  overlapping ClinVar target SNVs (any star): {fm_snv.in_clinvar.sum()}")
print(f"  overlapping ClinVar >=1 star: {fm_snv.in_clinvar_ge1.sum()}")
print(f"  overlapping ClinVar VUS (>=1 star): {fm_snv.in_clinvar_vus.sum()}")
