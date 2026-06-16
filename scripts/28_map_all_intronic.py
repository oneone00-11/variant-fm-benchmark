"""Map all currently-unmapped single-nucleotide functional variants (splice /
intronic / UTR that MaveDB couldn't place) to GRCh38 via Mutalyzer, then write
the completed functional_scores_master.tsv. Cached + restartable.
"""
import re
import sys
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

from config import DATA_PROCESSED
from mutalyzer_map import (build_query, ctoken_from_hgvs, load_cache, map_one,
                           save_cache, GENE_REF)

MASTER = DATA_PROCESSED / "functional_scores_master.tsv"
SINGLE_SUB = re.compile(r":c\.[\*\-]?\d+(?:[+\-]\d+)?[ACGT]>[ACGT]$")

df = pd.read_csv(MASTER, sep="\t", dtype={"chrom": str},
                 low_memory=False)

todo = df[(~df["mapped"]) & (df["gene"].isin(GENE_REF))
          & (df["hgvs_nt"].str.contains(SINGLE_SUB, regex=True, na=False))]
work = todo[["gene", "hgvs_nt"]].drop_duplicates()
print(f"{len(todo)} unmapped single-sub rows; {len(work)} unique to map", flush=True)

cache = load_cache()
lock = threading.Lock()
done = [0]
errors = [0]


def task(gene, hgvs_nt):
    ct = ctoken_from_hgvs(hgvs_nt)
    key = build_query(gene, ct)
    if key in cache:
        return hgvs_nt, cache[key]
    try:
        res = map_one(gene, ct)
    except Exception:
        with lock:
            errors[0] += 1
        return hgvs_nt, None
    val = list(res) if res else None
    with lock:
        cache[key] = val
        done[0] += 1
        if done[0] % 250 == 0:
            save_cache(cache)
            print(f"  mapped {done[0]} new (errors {errors[0]})", flush=True)
    return hgvs_nt, val


results = {}
with ThreadPoolExecutor(max_workers=6) as ex:
    futs = [ex.submit(task, r.gene, r.hgvs_nt) for r in work.itertuples()]
    for f in as_completed(futs):
        h, v = f.result()
        if v:
            results[h] = v
save_cache(cache)
# also fold in any cache hits not in results
for r in work.itertuples():
    if r.hgvs_nt not in results:
        key = build_query(r.gene, ctoken_from_hgvs(r.hgvs_nt))
        if cache.get(key):
            results[r.hgvs_nt] = cache[key]

print(f"newly mapped: {len(results)} (give-up errors: {errors[0]})", flush=True)

# Apply to the master
def apply_row(row):
    if row["mapped"] or row["hgvs_nt"] not in results:
        return row
    chrom, pos, ref, alt = results[row["hgvs_nt"]]
    row["chrom"], row["pos"], row["ref"], row["alt"] = chrom, pos, ref, alt
    row["is_snv"], row["mapped"], row["mapping_source"] = True, True, "mutalyzer"
    return row

if "mapping_source" not in df.columns:
    df["mapping_source"] = df["mapped"].map(lambda m: "mavedb_postmapped" if m else "")
df = df.apply(apply_row, axis=1)
df["chrom"] = df["chrom"].astype("string")
df.to_csv(MASTER, sep="\t", index=False)

n_snv = int(df["is_snv"].sum())
print(f"master now has {n_snv} SNV-coord variants of {len(df)} total", flush=True)
print(df.groupby(["gene"]).agg(n=("variant_urn", "count"),
                               mapped_snv=("is_snv", "sum")).to_string())
