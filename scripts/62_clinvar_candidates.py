"""M4: extract ClinVar SNVs for candidate genes (VHL, BAP1) from the full VCF,
same parsing/caliber as Milestone 1 -> data/processed/candidate_clinvar_variants.tsv.
"""
import gzip
import pandas as pd

from config import (CLINVAR_VCF_GZ, DATA_PROCESSED, classify_clnsig,
                    classify_consequence, revstat_to_stars)

CAND_GENES = {"VHL", "BAP1"}


def parse_info(s):
    d = {}
    for f in s.split(";"):
        if "=" in f:
            k, v = f.split("=", 1); d[k] = v
        else:
            d[f] = True
    return d


def genes_from_geneinfo(gi):
    return {t.split(":", 1)[0] for t in gi.split("|") if t} if gi else set()


def mc_terms(mc):
    out = []
    if mc:
        for it in mc.split(","):
            p = it.split("|")
            if len(p) == 2: out.append(p[1])
    return out


rows = []
n = 0
with gzip.open(CLINVAR_VCF_GZ, "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 8:
            continue
        n += 1
        info = parse_info(c[7])
        hit = genes_from_geneinfo(info.get("GENEINFO", "")) & CAND_GENES
        if not hit:
            continue
        terms = mc_terms(info.get("MC", ""))
        rows.append({
            "variation_id": c[2], "chrom": c[0], "pos": int(c[1]),
            "ref": c[3], "alt": c[4], "gene": sorted(hit)[0],
            "clnsig_raw": info.get("CLNSIG", ""),
            "clnsig_class": classify_clnsig(info.get("CLNSIG", "")),
            "revstat": info.get("CLNREVSTAT", ""),
            "stars": revstat_to_stars(info.get("CLNREVSTAT", "")),
            "clnvc": info.get("CLNVC", ""),
            "mc_terms": ";".join(terms),
            "consequence_class": classify_consequence(terms),
        })

df = pd.DataFrame(rows)
out = DATA_PROCESSED / "candidate_clinvar_variants.tsv"
df.to_csv(out, sep="\t", index=False)
print(f"Scanned {n:,} records; {len(df):,} on VHL/BAP1")
print(f"Wrote {out}")
print("\nSNV counts (>=1 star) by gene x clnsig:")
snv = df[(df.clnvc == "single_nucleotide_variant") & (df.stars >= 1)]
print(snv.pivot_table(index="gene", columns="clnsig_class", values="variation_id",
                      aggfunc="count", fill_value=0).to_string())
