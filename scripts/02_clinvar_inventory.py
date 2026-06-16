"""ClinVar inventory for the target genes (Milestone 1).

Streams the GRCh38 ClinVar VCF, keeps variants whose GENEINFO names one of the
target genes, classifies each by clinical significance + review stars + coarse
molecular consequence, and writes:

  data/processed/clinvar_target_variants.tsv   -- the filtered variant master seed
  results/tables/clinvar_inventory_by_gene.tsv -- P/B/VUS counts per gene (>=1 star)
  results/tables/clinvar_inventory_by_class.tsv-- gene x consequence breakdown (>=1 star)

Prints a readable summary to stdout.

Run:  python scripts/02_clinvar_inventory.py
"""
import gzip
import sys
from collections import defaultdict

import pandas as pd

from config import (
    CLINVAR_VCF_GZ, DATA_PROCESSED, RESULTS_TABLES, TARGET_GENES,
    classify_clnsig, classify_consequence, revstat_to_stars,
)

TARGET_SET = set(TARGET_GENES)


def parse_info(info_str):
    d = {}
    for field in info_str.split(";"):
        if "=" in field:
            k, v = field.split("=", 1)
            d[k] = v
        else:
            d[field] = True
    return d


def genes_from_geneinfo(geneinfo):
    """GENEINFO='BRCA1:672|RPL21P4:...' -> {'BRCA1', 'RPL21P4'}."""
    if not geneinfo:
        return set()
    return {tok.split(":", 1)[0] for tok in geneinfo.split("|") if tok}


def mc_terms(mc_field):
    """MC='SO:0001583|missense_variant,SO:...|...' -> ['missense_variant', ...]."""
    if not mc_field:
        return []
    terms = []
    for item in mc_field.split(","):
        parts = item.split("|")
        if len(parts) == 2:
            terms.append(parts[1])
    return terms


def iter_variants(vcf_gz):
    with gzip.open(vcf_gz, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            chrom, pos, vid, ref, alt, _qual, _filt, info = cols[:8]
            yield chrom, pos, vid, ref, alt, parse_info(info)


def main():
    if not CLINVAR_VCF_GZ.exists():
        sys.exit(f"Missing {CLINVAR_VCF_GZ}. Run 01_download_clinvar.py first.")

    rows = []
    n_total = 0
    for chrom, pos, vid, ref, alt, info in iter_variants(CLINVAR_VCF_GZ):
        n_total += 1
        genes = genes_from_geneinfo(info.get("GENEINFO", ""))
        hit = genes & TARGET_SET
        if not hit:
            continue
        # If multiple target genes overlap (rare), keep them joined; pick the
        # first target gene as the primary label for grouping.
        primary_gene = sorted(hit)[0]
        clnsig = info.get("CLNSIG", "")
        revstat = info.get("CLNREVSTAT", "")
        terms = mc_terms(info.get("MC", ""))
        rows.append({
            "variation_id": vid,
            "chrom": chrom,
            "pos": int(pos),
            "ref": ref,
            "alt": alt,
            "gene": primary_gene,
            "all_target_genes": "|".join(sorted(hit)),
            "clnsig_raw": clnsig,
            "clnsig_class": classify_clnsig(clnsig),
            "revstat": revstat,
            "stars": revstat_to_stars(revstat),
            "clnvc": info.get("CLNVC", ""),
            "mc_terms": ";".join(terms),
            "consequence_class": classify_consequence(terms),
            "rs": info.get("RS", ""),
            "allele_id": info.get("ALLELEID", ""),
        })

    if not rows:
        sys.exit("No target-gene variants found — check GENEINFO parsing.")

    df = pd.DataFrame(rows)
    print(f"Scanned {n_total:,} ClinVar records; "
          f"{len(df):,} on target genes {TARGET_GENES}")

    # Save the full filtered table (all stars) as the variant master seed.
    out_master = DATA_PROCESSED / "clinvar_target_variants.tsv"
    df.to_csv(out_master, sep="\t", index=False)
    print(f"Wrote variant master seed -> {out_master}")

    # --- Primary inventory: >= 1 star --------------------------------------
    hi = df[df["stars"] >= 1].copy()

    by_gene = (
        hi.pivot_table(index="gene", columns="clnsig_class",
                       values="variation_id", aggfunc="count", fill_value=0)
        .reindex(TARGET_GENES)
    )
    # ensure all class columns exist
    for c in ["P/LP", "B/LB", "VUS", "Conflicting", "Other"]:
        if c not in by_gene.columns:
            by_gene[c] = 0
    by_gene = by_gene[["P/LP", "B/LB", "VUS", "Conflicting", "Other"]]
    by_gene["TOTAL_ge1star"] = by_gene.sum(axis=1)
    by_gene.loc["TOTAL"] = by_gene.sum(axis=0)
    by_gene.to_csv(RESULTS_TABLES / "clinvar_inventory_by_gene.tsv", sep="\t")

    # --- Stratified: gene x consequence class (>=1 star) -------------------
    by_class = (
        hi.pivot_table(index="gene", columns="consequence_class",
                       values="variation_id", aggfunc="count", fill_value=0)
        .reindex(TARGET_GENES)
    )
    by_class.to_csv(RESULTS_TABLES / "clinvar_inventory_by_class.tsv", sep="\t")

    # --- Console summary ---------------------------------------------------
    pd.set_option("display.width", 140)
    pd.set_option("display.max_columns", 30)
    print("\n=== ClinVar inventory by gene (review stars >= 1) ===")
    print(by_gene.to_string())

    print("\n=== Variant consequence classes (review stars >= 1) ===")
    print(by_class.fillna(0).astype(int).to_string())

    # Splice + non-coding focus (the project's gap of interest)
    gap = hi[hi["consequence_class"].isin(["splice", "noncoding"])]
    print("\n=== Splice / non-coding subset (>=1 star) — the project's gap ===")
    print(gap.groupby(["gene", "consequence_class"]).size()
          .unstack(fill_value=0).reindex(TARGET_GENES).fillna(0).astype(int)
          .to_string())
    print(f"\nTotal splice+noncoding (>=1 star): {len(gap)}")
    print(f"  of which P/LP: {(gap['clnsig_class']=='P/LP').sum()}, "
          f"B/LB: {(gap['clnsig_class']=='B/LB').sum()}, "
          f"VUS: {(gap['clnsig_class']=='VUS').sum()}")


if __name__ == "__main__":
    main()
