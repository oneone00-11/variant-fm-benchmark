"""M3-A: clean ClinVar to SNV-only, classify region precisely, inner-join with
the functional gold standard -> analysis_ready.tsv.

Region rule (unified, documented):
  - intron-side variants (functional intron_offset present):
      |off| <= 2          -> splice  (splice_class=splice_core)
      3 <= |off| <= 8     -> splice  (splice_class=splice_region)
      |off| > 8           -> intronic (intronic_depth: proximal <=50, else deep)
  - exonic / UTR (no offset): from ClinVar mc_terms SO:
      splice_donor/acceptor -> splice/splice_core (safety)
      nonsense              -> nonsense
      missense_variant      -> missense
      synonymous_variant    -> synonymous
      5'/3'_UTR_variant     -> utr
      frameshift_variant    -> frameshift_SNV (≈0; SNV can't frameshift)
      intron_variant only   -> intronic (depth unknown)
      else                  -> other
NB: ClinVar MC lacks a splice_region term, so intron-side splice_region comes
from intron_offset; exon-side last-1-3 splice-region is NOT detectable here.
"""
import numpy as np
import pandas as pd

from config import DATA_PROCESSED, TARGET_GENES

FINAL_GENES = ["BRCA1", "BRCA2", "BARD1", "PALB2", "RAD51C"]
CIRCULARITY_TRUE = {"BRCA1"}  # Findlay SGE deposited to ClinVar; verified by ~100% overlap


def norm_chrom(x):
    s = str(x)
    return s[:-2] if s.endswith(".0") else s


def region_of(off, mc):
    terms = set(str(mc).split(";")) if mc and str(mc) != "nan" else set()
    if pd.notna(off):
        a = abs(off)
        if a <= 2:
            return "splice", "splice_core", None
        if a <= 8:
            return "splice", "splice_region", None
        return "intronic", None, ("proximal" if a <= 50 else "deep")
    if {"splice_donor_variant", "splice_acceptor_variant"} & terms:
        return "splice", "splice_core", None
    if "nonsense" in terms or "stop_gained" in terms:
        return "nonsense", None, None
    if "missense_variant" in terms:
        return "missense", None, None
    if "synonymous_variant" in terms:
        return "synonymous", None, None
    if {"5_prime_UTR_variant", "3_prime_UTR_variant"} & terms:
        return "utr", None, None
    if "frameshift_variant" in terms:
        return "frameshift_SNV", None, None
    if "intron_variant" in terms:
        return "intronic", None, "unknown"
    return "other", None, None


def label_set_of(clnsig):
    if clnsig in ("P/LP", "B/LB"):
        return "clean_PB"
    if clnsig == "Conflicting":
        return "conflicting"
    if clnsig == "VUS":
        return "vus"
    return "other"


# --- load ---
cv = pd.read_csv(DATA_PROCESSED / "clinvar_target_variants.tsv", sep="\t",
                 dtype={"chrom": str}, low_memory=False)
fm = pd.read_csv(DATA_PROCESSED / "functional_scores_master.tsv", sep="\t",
                 dtype={"chrom": str}, low_memory=False)

# --- clean ClinVar: SNV vs non-SNV ---
cv["chrom"] = cv["chrom"].map(norm_chrom)
is_snv = cv["clnvc"] == "single_nucleotide_variant"
cv_excluded = cv[~is_snv].copy()
cv_excluded["exclude_reason"] = "ClinVar clnvc=" + cv_excluded["clnvc"]
cv_excluded.to_csv(DATA_PROCESSED / "excluded_nonSNV.tsv", sep="\t", index=False)
cv_snv = cv[is_snv].copy()
print("=== ClinVar cleaning ===")
print(f"  SNV kept: {len(cv_snv):,} | non-SNV excluded: {len(cv_excluded):,}")
print(cv_excluded["clnvc"].value_counts().to_string())

# --- functional: SNV vs non-SNV (codon delins/MNV) ---
fm["chrom"] = fm["chrom"].map(lambda x: norm_chrom(x) if pd.notna(x) else x)
fm_snv = fm[fm["is_snv"] == True].copy()
fm_non = fm[fm["is_snv"] != True].copy()
fm_non.to_csv(DATA_PROCESSED / "functional_nonSNV.tsv", sep="\t", index=False)
print(f"\n=== Functional cleaning ===\n  SNV-coord: {len(fm_snv):,} | "
      f"non-SNV/unmapped held out: {len(fm_non):,}")

# per-gene overlap fraction (for circularity reasoning)
print("\n=== per-gene functional-SNV overlap with ClinVar (circularity signal) ===")
cv_keys = set(cv_snv.chrom + ":" + cv_snv.pos.astype(str) + ":" + cv_snv.ref + ":" + cv_snv.alt)
for g in FINAL_GENES:
    sub = fm_snv[fm_snv.gene == g]
    k = sub.chrom + ":" + sub.pos.astype("Int64").astype(str) + ":" + sub.ref + ":" + sub.alt
    ov = k.isin(cv_keys).sum()
    print(f"  {g}: {ov}/{len(sub)} = {ov/max(len(sub),1):.0%} in ClinVar")

# --- merge (inner) on chrom,pos,ref,alt ---
fm_snv["pos"] = fm_snv["pos"].astype("Int64")
key_cols = ["chrom", "pos", "ref", "alt"]
fcols = key_cols + ["gene", "intron_offset", "functional_score",
                    "functional_call", "source_urn", "assay_type", "hgvs_nt"]
ccols = key_cols + ["variation_id", "clnsig_class", "stars", "mc_terms",
                    "consequence_class"]
merged = fm_snv[fcols].merge(cv_snv[ccols], on=key_cols, how="inner",
                             suffixes=("", "_cv"))
print(f"\n=== Inner join: {len(merged):,} variants with BOTH functional + ClinVar ===")

# --- region classification ---
reg = merged.apply(lambda r: region_of(r["intron_offset"], r["mc_terms"]),
                   axis=1, result_type="expand")
merged["region_class"] = reg[0]
merged["splice_class"] = reg[1]
merged["intronic_depth"] = reg[2]

# --- flags ---
merged["circularity_flag"] = merged["gene"].isin(CIRCULARITY_TRUE)
merged["label_set"] = merged["clnsig_class"].map(label_set_of)
merged["path_binary"] = np.where(merged["clnsig_class"] == "P/LP", 1,
                          np.where(merged["clnsig_class"] == "B/LB", 0, np.nan))

out_cols = ["gene", "chrom", "pos", "ref", "alt", "region_class", "splice_class",
            "intronic_depth", "intron_offset", "functional_score", "functional_call",
            "clnsig_class", "stars", "path_binary", "label_set", "circularity_flag",
            "source_urn", "assay_type", "variation_id", "mc_terms", "hgvs_nt"]
merged[out_cols].to_csv(DATA_PROCESSED / "analysis_ready.tsv", sep="\t", index=False)
print(f"Wrote analysis_ready.tsv ({len(merged):,} rows)")

print("\n=== region_class x gene (analysis_ready) ===")
print(merged.pivot_table(index="region_class", columns="gene", values="pos",
                         aggfunc="count", fill_value=0).to_string())
print("\n=== splice_class breakdown ===")
print(merged[merged.region_class == "splice"]
      .groupby(["gene", "splice_class"]).size().unstack(fill_value=0).to_string())
print("\n=== label_set counts ===")
print(merged.label_set.value_counts().to_string())

# --- spot-check: 5 known BRCA1 splice-core variants should be P/LP + damaging ---
print("\n=== SPOT-CHECK: BRCA1 splice_core, P/LP, lowest functional scores ===")
sc = merged[(merged.gene == "BRCA1") & (merged.splice_class == "splice_core")
            & (merged.clnsig_class == "P/LP")].nsmallest(5, "functional_score")
print(sc[["hgvs_nt", "chrom", "pos", "ref", "alt", "intron_offset",
          "functional_score", "clnsig_class", "stars", "variation_id"]].to_string(index=False))
print("\n=== SPOT-CHECK: BRCA1 deep-intronic B/LB (should be near-zero score) ===")
bi = merged[(merged.gene == "BRCA1") & (merged.region_class == "intronic")
            & (merged.clnsig_class == "B/LB")].head(5)
print(bi[["hgvs_nt", "pos", "intron_offset", "functional_score", "clnsig_class"]].to_string(index=False))
