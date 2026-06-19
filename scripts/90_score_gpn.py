"""GPN-MSA scoring: remote-query precomputed scores for our 7 gene regions.
Robust: fresh TabixFile + retry (with CDN-URL re-resolve) per gene region.
Run in WSL (spliceai env, pysam). argv[1] = initial resolved direct CDN bgz URL.
"""
import csv
import os
import subprocess
import sys
import time

import pysam

RESOLVE = ("https://huggingface.co/datasets/songlab/gpn-msa-hg38-scores/"
           "resolve/main/scores.tsv.bgz")
TBI = "/mnt/d/variant-fm-benchmark/refs/gpn_scores.tsv.bgz.tbi"
AR = "/mnt/d/variant-fm-benchmark/data/processed/analysis_ready_v2.tsv"
OUT = "/mnt/d/variant-fm-benchmark/data/raw/scores/gpn/gpn.tsv"


def resolve_url():
    return subprocess.check_output(
        ["curl", "-s", "-o", "/dev/null", "-w", "%{url_effective}", "-I", "-L", RESOLVE]
    ).decode().strip()


url = sys.argv[1] if len(sys.argv) > 1 else resolve_url()

keys, order, gene_region = set(), [], {}
with open(AR) as f:
    for d in csv.DictReader(f, delimiter="\t"):
        c, p = d["chrom"], int(float(d["pos"]))
        ref, alt, gene = d["ref"], d["alt"], d["gene"]
        keys.add((c, p, ref, alt)); order.append((c, p, ref, alt))
        if gene not in gene_region:
            gene_region[gene] = [c, p, p]
        gene_region[gene][1] = min(gene_region[gene][1], p)
        gene_region[gene][2] = max(gene_region[gene][2], p)

score = {}
for g, (c, lo, hi) in gene_region.items():
    for attempt in range(5):
        try:
            tbx = pysam.TabixFile(url, index=TBI)
            n = 0
            for line in tbx.fetch(c, lo - 1, hi + 1):
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5:
                    continue
                k = (parts[0], int(parts[1]), parts[2], parts[3])
                if k in keys:
                    score[k] = parts[4]; n += 1
            tbx.close()
            print(f"  {g} {c}:{lo}-{hi}: matched {n}", flush=True)
            break
        except Exception as e:
            print(f"  {g} attempt {attempt} failed: {str(e)[:60]}; re-resolving", flush=True)
            time.sleep(3)
            try:
                url = resolve_url()
            except Exception:
                pass
    else:
        print(f"  {g}: GAVE UP after retries", flush=True)

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as o:
    o.write("chrom\tpos\tref\talt\tgpn_msa_score\n")
    for (c, p, ref, alt) in order:
        o.write(f"{c}\t{p}\t{ref}\t{alt}\t{score.get((c, p, ref, alt), '')}\n")
print(f"GPN-MSA scored {len(score)} of {len(order)} -> {OUT}", flush=True)
