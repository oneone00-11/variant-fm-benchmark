"""Validate the Mutalyzer mapping BEFORE the full run (subtask 7 + de-risk):

(A) ENST->NM MANE equivalence: take CODING BRCA2 & RAD51C variants that MaveDB
    already post-mapped (genomic truth), re-map their c. via the RefSeq NM with
    Mutalyzer, and confirm identical GRCh38 coords. If they match, mapping
    intronic variants of these genes through the NM equivalent is justified.

(B) Known BRCA1 SGE pathogenic splice variants: map a couple and confirm the
    genomic coords land in the ClinVar master with a pathogenic label.
"""
import pandas as pd

from config import DATA_PROCESSED
from mutalyzer_map import map_one, ctoken_from_hgvs

fm = pd.read_csv(DATA_PROCESSED / "functional_scores_master.tsv", sep="\t")
cv = pd.read_csv(DATA_PROCESSED / "clinvar_target_variants.tsv", sep="\t",
                 dtype={"chrom": str})

print("=== (A) ENST->NM coding cross-check vs MaveDB post-mapped ===")
for gene in ["BRCA2", "RAD51C"]:
    sub = fm[(fm.gene == gene) & (fm.is_snv) & (fm.region_class == "coding_other")].head(6)
    ok = bad = 0
    for _, row in sub.iterrows():
        ct = ctoken_from_hgvs(row["hgvs_nt"])
        try:
            res = map_one(gene, ct)
        except Exception as e:
            print(f"  {gene} {row['hgvs_nt']}: ERROR {e}")
            continue
        truth = (str(row["chrom"]), int(row["pos"]), row["ref"], row["alt"])
        if res == truth:
            ok += 1
        else:
            bad += 1
            print(f"  MISMATCH {gene} {row['hgvs_nt']}: mutalyzer={res} mavedb={truth}")
    print(f"  {gene}: {ok} match, {bad} mismatch (of {len(sub)})")

print("\n=== (B) BRCA1 splice/intronic -> ClinVar lookup ===")
b1 = fm[(fm.gene == "BRCA1") & (fm.region_class.isin(["splice", "intronic"]))
        & (~fm.mapped)].head(8)
cv_snv = cv[cv.clnvc == "single_nucleotide_variant"].copy()
cv_snv["key"] = (cv_snv.chrom.astype(str) + ":" + cv_snv.pos.astype(str)
                 + ":" + cv_snv.ref + ":" + cv_snv.alt)
cvmap = cv_snv.set_index("key")[["clnsig_class", "consequence_class", "stars"]].to_dict("index")
for _, row in b1.iterrows():
    ct = ctoken_from_hgvs(row["hgvs_nt"])
    try:
        res = map_one("BRCA1", ct)
    except Exception as e:
        print(f"  {row['hgvs_nt']}: ERROR {e}")
        continue
    if res is None:
        print(f"  {row['hgvs_nt']}: unmappable")
        continue
    key = f"{res[0]}:{res[1]}:{res[2]}:{res[3]}"
    hit = cvmap.get(key)
    print(f"  {row['hgvs_nt']} [{row['region_class']}] -> {key}  "
          f"funcScore={row['functional_score']:.3f}  ClinVar={hit}")
