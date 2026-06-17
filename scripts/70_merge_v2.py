"""M5-A: merge BAP1 + VHL into the 7-gene analysis_ready_v2.tsv, using the EXACT
same column structure and region rule as analysis_ready (script 51).
"""
import numpy as np
import pandas as pd

from config import DATA_PROCESSED

CIRCULARITY_TRUE = {"BRCA1"}  # M4: BAP1=49%, VHL=49/32% overlap -> FALSE
CAND = ["BAP1", "VHL"]


def norm_chrom(x):
    s = str(x)
    return s[:-2] if s.endswith(".0") else s


def region_of(off, mc):
    terms = set(str(mc).split(";")) if mc and str(mc) != "nan" else set()
    if pd.notna(off):
        a = abs(off)
        if a <= 2: return "splice", "splice_core", None
        if a <= 8: return "splice", "splice_region", None
        return "intronic", None, ("proximal" if a <= 50 else "deep")
    if {"splice_donor_variant", "splice_acceptor_variant"} & terms:
        return "splice", "splice_core", None
    if "nonsense" in terms or "stop_gained" in terms: return "nonsense", None, None
    if "missense_variant" in terms: return "missense", None, None
    if "synonymous_variant" in terms: return "synonymous", None, None
    if {"5_prime_UTR_variant", "3_prime_UTR_variant"} & terms: return "utr", None, None
    if "frameshift_variant" in terms: return "frameshift_SNV", None, None
    if "intron_variant" in terms: return "intronic", None, "unknown"
    return "other", None, None


def label_set_of(c):
    if c in ("P/LP", "B/LB"): return "clean_PB"
    if c == "Conflicting": return "conflicting"
    if c == "VUS": return "vus"
    return "other"


# --- load existing 5-gene set + candidate sources ---
ar = pd.read_csv(DATA_PROCESSED / "analysis_ready.tsv", sep="\t",
                 dtype={"chrom": str}, low_memory=False)
fm = pd.read_csv(DATA_PROCESSED / "candidate_functional_master.tsv", sep="\t",
                 dtype={"chrom": str}, low_memory=False)
cv = pd.read_csv(DATA_PROCESSED / "candidate_clinvar_variants.tsv", sep="\t",
                 dtype={"chrom": str}, low_memory=False)

fm = fm[(fm.is_snv == True) & (fm.gene.isin(CAND))].copy()
fm["chrom"] = fm["chrom"].map(norm_chrom)
fm["pos"] = fm["pos"].astype("Int64")
cv = cv[cv.clnvc == "single_nucleotide_variant"].copy()
cv["chrom"] = cv["chrom"].map(norm_chrom)

key = ["chrom", "pos", "ref", "alt"]
fcols = key + ["gene", "intron_offset", "functional_score", "functional_call",
               "source_urn", "assay_type", "hgvs_nt"]
ccols = key + ["variation_id", "clnsig_class", "stars", "mc_terms", "consequence_class"]
merged = fm[fcols].merge(cv[ccols], on=key, how="inner")
print(f"Candidate inner-join (BAP1+VHL): {len(merged):,} variants")

reg = merged.apply(lambda r: region_of(r["intron_offset"], r["mc_terms"]),
                   axis=1, result_type="expand")
merged["region_class"], merged["splice_class"], merged["intronic_depth"] = reg[0], reg[1], reg[2]
merged["circularity_flag"] = merged["gene"].isin(CIRCULARITY_TRUE)
merged["label_set"] = merged["clnsig_class"].map(label_set_of)
merged["path_binary"] = np.where(merged.clnsig_class == "P/LP", 1,
                          np.where(merged.clnsig_class == "B/LB", 0, np.nan))

merged = merged[ar.columns]  # exact same column order
v2 = pd.concat([ar, merged], ignore_index=True)
v2.to_csv(DATA_PROCESSED / "analysis_ready_v2.tsv", sep="\t", index=False)
print(f"Wrote analysis_ready_v2.tsv: {len(v2):,} rows, {v2.gene.nunique()} genes")

# --- completeness check ---
print("\n=== rows per gene ===")
print(v2.gene.value_counts().reindex(
    ["BRCA1","BRCA2","BARD1","PALB2","RAD51C","BAP1","VHL"]).to_string())
print("\n=== region_class x gene ===")
print(v2.pivot_table(index="region_class", columns="gene", values="pos",
                     aggfunc="count", fill_value=0).to_string())
print("\n=== label_set per gene ===")
print(v2.pivot_table(index="gene", columns="label_set", values="pos",
                     aggfunc="count", fill_value=0).to_string())
print("\n=== circularity_flag per gene ===")
print(v2.groupby("gene").circularity_flag.first().to_string())
print(f"\nexisting 5-gene rows unchanged: {len(ar):,} (expect 17,272)")

# --- spot-check 5 known candidate variants ---
print("\n=== SPOT-CHECK: VHL & BAP1 splice_core P/LP, lowest functional_score ===")
sc = merged[(merged.splice_class == "splice_core") & (merged.clnsig_class == "P/LP")]
sc = sc.sort_values("functional_score").groupby("gene").head(3)
print(sc[["gene","hgvs_nt","chrom","pos","ref","alt","intron_offset",
          "functional_score","clnsig_class","stars","variation_id"]].to_string(index=False))
