"""M4: validate candidate mapping, intersect with ClinVar, compute circularity,
and build the coverage + splice-N comparison tables.  No model scoring.
"""
import numpy as np
import pandas as pd

from config import DATA_PROCESSED, RESULTS_TABLES
from mutalyzer_map import map_one, ctoken_from_hgvs

EXISTING = ["BRCA1", "BRCA2", "BARD1", "PALB2", "RAD51C"]
CAND = ["VHL", "BAP1"]


def norm_chrom(x):
    s = str(x)
    return s[:-2] if s.endswith(".0") else s


def offset_region(off):
    if pd.isna(off):
        return None
    a = abs(off)
    if a <= 2: return "splice_core"
    if a <= 8: return "splice_region"
    return "intronic"


def func_region(row):
    """Region for a functional-only variant: offset-based for gap, else coarse."""
    r = offset_region(row["intron_offset"])
    if r: return r
    rc = row["region_class"]
    if rc == "coding_other": return "coding"
    if rc == "noncoding": return "utr"
    return rc or "other"


# --- load ---
cand = pd.read_csv(DATA_PROCESSED / "candidate_functional_master.tsv", sep="\t",
                   dtype={"chrom": str}, low_memory=False)
fm = pd.read_csv(DATA_PROCESSED / "functional_scores_master.tsv", sep="\t",
                 dtype={"chrom": str}, low_memory=False)
# VHL/BAP1 ClinVar SNVs (extracted from full VCF by script 62 — they are NOT in
# the original 6-gene master).
cv = pd.read_csv(DATA_PROCESSED / "candidate_clinvar_variants.tsv", sep="\t",
                 dtype={"chrom": str}, low_memory=False)
cv = cv[cv.clnvc == "single_nucleotide_variant"].copy()
cv["chrom"] = cv["chrom"].map(norm_chrom)

print("=== VALIDATION A: ENST->NM coding cross-check vs MaveDB post-mapping ===")
for gene in CAND:
    sub = cand[(cand.gene == gene) & (cand.is_snv) &
               (cand.mapping_source == "mavedb_postmapped") &
               (cand.region_class == "coding_other")].head(12)
    ok = bad = 0
    for _, row in sub.iterrows():
        try:
            res = map_one(gene, ctoken_from_hgvs(row["hgvs_nt"]))
        except Exception as e:
            print(f"  {gene} {row['hgvs_nt']}: ERR {e}"); continue
        truth = (norm_chrom(row["chrom"]), int(row["pos"]), row["ref"], row["alt"])
        if res == truth: ok += 1
        else:
            bad += 1
            print(f"  MISMATCH {gene} {row['hgvs_nt']}: mut={res} mavedb={truth}")
    print(f"  {gene}: {ok}/{len(sub)} coding coords match")

print("\n=== give-up / unmapped rate per candidate ===")
for gene in CAND:
    g = cand[cand.gene == gene]
    single_sub = g.hgvs_nt.str.contains(r":c\.[\*\-]?\d+(?:[+\-]\d+)?[ACGT]>[ACGT]$",
                                        regex=True, na=False)
    n_snv = int(g.is_snv.sum())
    n_unmapped_sub = int((single_sub & ~g.mapped).sum())
    print(f"  {gene}: {n_snv} SNV-coord; {n_unmapped_sub} single-sub still unmapped "
          f"of {len(g)} total variants")

# --- ClinVar overlap + circularity ---
cv_ge1 = set((cv[cv.stars >= 1].chrom + ":" + cv[cv.stars >= 1].pos.astype(str)
              + ":" + cv[cv.stars >= 1].ref + ":" + cv[cv.stars >= 1].alt))
cv_vus = cv[(cv.stars >= 1) & (cv.clnsig_class == "VUS")]
cv_vus_keys = set(cv_vus.chrom + ":" + cv_vus.pos.astype(str) + ":" + cv_vus.ref + ":" + cv_vus.alt)
cv_all = set(cv.chrom + ":" + cv.pos.astype(str) + ":" + cv.ref + ":" + cv.alt)

print("\n=== VALIDATION B: candidate splice -> ClinVar P/LP, score direction ===")
for gene in CAND:
    g = cand[(cand.gene == gene) & cand.is_snv].copy()
    g["key"] = g.chrom + ":" + g.pos.astype("Int64").astype(str) + ":" + g.ref + ":" + g.alt
    g["freg"] = g.apply(func_region, axis=1)
    sp = g[g.freg.isin(["splice_core", "splice_region"])].copy()
    # join clinsig
    cvk = cv.copy()
    cvk["key"] = cvk.chrom + ":" + cvk.pos.astype(str) + ":" + cvk.ref + ":" + cvk.alt
    cvk = cvk[cvk.stars >= 1].set_index("key")[["clnsig_class", "stars"]]
    sp = sp.join(cvk, on="key")
    plp = sp[sp.clnsig_class == "P/LP"].nsmallest(5, "functional_score")
    print(f"  {gene}: {len(plp)} P/LP splice shown (lowest score = most damaging)")
    print(plp[["hgvs_nt", "freg", "functional_score", "clnsig_class", "stars"]]
          .to_string(index=False))

# --- coverage table (candidates) + circularity ---
rows = []
circ = {}
for gene in CAND:
    g = cand[(cand.gene == gene) & cand.is_snv].copy()
    g["key"] = g.chrom + ":" + g.pos.astype("Int64").astype(str) + ":" + g.ref + ":" + g.alt
    g["freg"] = g.apply(func_region, axis=1)
    circ[gene] = (g.key.isin(cv_all).sum(), len(g))
    for reg, sub in g.groupby("freg"):
        rows.append({
            "gene": gene, "region": reg,
            "n_functional_snv": len(sub),
            "n_overlap_clinvar_ge1": int(sub.key.isin(cv_ge1).sum()),
            "n_overlap_vus": int(sub.key.isin(cv_vus_keys).sum()),
        })
cov = pd.DataFrame(rows)
cov.to_csv(RESULTS_TABLES / "candidate_gene_coverage.tsv", sep="\t", index=False)
print("\n=== candidate_gene_coverage.tsv ===")
print(cov.to_string(index=False))
print("\n=== circularity ratio (functional SNV already in ClinVar) ===")
for gene in CAND:
    inn, tot = circ[gene]
    print(f"  {gene}: {inn}/{tot} = {inn/max(tot,1):.0%}  "
          f"-> circularity_flag={'TRUE' if inn/max(tot,1) >= 0.9 else 'FALSE'}")

# --- splice-N comparison: existing 5 vs candidates ---
print("\n=== SPLICE (core+region) evaluable-N comparison ===")
allfm = pd.concat([fm[fm.gene.isin(EXISTING)], cand], ignore_index=True)
allfm = allfm[allfm.is_snv].copy()
allfm["freg"] = allfm.apply(func_region, axis=1)
comp = []
for gene in EXISTING + CAND:
    g = allfm[allfm.gene == gene]
    sc = int((g.freg == "splice_core").sum())
    sr = int((g.freg == "splice_region").sum())
    comp.append({"gene": gene, "splice_core": sc, "splice_region": sr,
                 "splice_total": sc + sr,
                 "group": "existing" if gene in EXISTING else "CANDIDATE"})
comp = pd.DataFrame(comp)
pooled_existing = comp[comp.group == "existing"].splice_total.sum()
print(comp.to_string(index=False))
print(f"\nPooled splice N (existing 5 genes): {pooled_existing}")
for gene in CAND:
    add = comp[comp.gene == gene].splice_total.values[0]
    print(f"  + {gene}: +{add} -> pooled {pooled_existing + add} "
          f"(+{add/pooled_existing:.0%})")
both = comp[comp.group == 'CANDIDATE'].splice_total.sum()
print(f"  + both VHL&BAP1: +{both} -> pooled {pooled_existing + both} (+{both/pooled_existing:.0%})")
comp.to_csv(RESULTS_TABLES / "candidate_splice_comparison.tsv", sep="\t", index=False)
