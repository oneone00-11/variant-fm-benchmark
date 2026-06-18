"""M6: emit a GRCh38 VCF of all unique SNVs in analysis_ready_v2 for SpliceAI /
Pangolin input. Chrom naming = no 'chr' prefix (Ensembl style, matching our
table's '17'/'3' and the Ensembl reference FASTA/GTF). Sorted by chrom,pos.
"""
import pandas as pd
from config import DATA_PROCESSED

CHR_LEN = {  # GRCh38 lengths for the chromosomes our genes live on
    "2": 242193529, "3": 198295559, "13": 114364328,
    "16": 90338345, "17": 83257441,
}

df = pd.read_csv(DATA_PROCESSED / "analysis_ready_v2.tsv", sep="\t",
                 dtype={"chrom": str}, low_memory=False)
df["pos"] = df["pos"].astype(int)
v = df[["chrom", "pos", "ref", "alt"]].drop_duplicates().copy()
v["c_int"] = v["chrom"].map(lambda x: int(x) if x.isdigit() else 99)
v = v.sort_values(["c_int", "pos"])
print(f"{len(v)} unique SNVs on chroms {sorted(v.chrom.unique(), key=lambda x:int(x))}")

out = DATA_PROCESSED / "variants_grch38.vcf"
with open(out, "w", encoding="utf-8", newline="\n") as f:
    f.write("##fileformat=VCFv4.2\n")
    for c, ln in CHR_LEN.items():
        f.write(f"##contig=<ID={c},length={ln}>\n")
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for r in v.itertuples():
        f.write(f"{r.chrom}\t{r.pos}\t.\t{r.ref}\t{r.alt}\t.\t.\t.\n")
print(f"Wrote {out}")
