"""AlphaGenome variant scoring for all 21,410 SNVs (all regions).
- Key from env only (never printed/cached/logged).
- Splice scorers (SPLICE_SITES, SPLICE_SITE_USAGE, SPLICE_JUNCTIONS), 1MB interval.
- Aggregate = max |raw_score| across splice tracks (REF/ALT effect; higher = more
  splice-altering = more likely pathogenic).
- Concurrent (gentle), exponential-backoff retry, incremental cache, resumable
  (cached variants skipped). Run in WSL `alphagenome` env.
"""
import csv
import os
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

if not os.environ.get("ALPHAGENOME_API_KEY"):
    sys.exit("ALPHAGENOME_API_KEY not set - aborting (key must come from env).")

from alphagenome.models import dna_client, variant_scorers
from alphagenome.data import genome

# Repo-relative paths (portable). API key is read from the environment only.
ROOT = Path(__file__).resolve().parents[1]
AR = str(ROOT / "data" / "processed" / "analysis_ready_v2.tsv")
OUTDIR = str(ROOT / "data" / "raw" / "scores" / "alphagenome")
CACHE = os.path.join(OUTDIR, "ag_cache.tsv")
os.makedirs(OUTDIR, exist_ok=True)
N_WORKERS = 6

model = dna_client.create(os.environ["ALPHAGENOME_API_KEY"])
recs = variant_scorers.RECOMMENDED_VARIANT_SCORERS
SPLICE = [recs[n] for n in recs if "SPLICE" in str(n).upper()]

# load variants
variants = []
with open(AR) as f:
    for d in csv.DictReader(f, delimiter="\t"):
        variants.append((d["chrom"], int(float(d["pos"])), d["ref"], d["alt"]))
variants = list(dict.fromkeys(variants))  # unique, keep order

# resume: load cached keys
done = {}
if os.path.exists(CACHE):
    with open(CACHE) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) == 5 and p[0] != "chrom":
                done[(p[0], int(p[1]), p[2], p[3])] = p[4]
todo = [v for v in variants if v not in done]
print(f"{len(variants)} variants; {len(done)} cached; {len(todo)} to score", flush=True)

lock = threading.Lock()
cnt = [0]
if not os.path.exists(CACHE):
    with open(CACHE, "w") as f:
        f.write("chrom\tpos\tref\talt\talphagenome_splice\n")


def score_one(v):
    chrom, pos, ref, alt = v
    var = genome.Variant(chromosome="chr" + chrom, position=pos,
                         reference_bases=ref, alternate_bases=alt)
    iv = var.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
    last = None
    for attempt in range(5):
        try:
            out = model.score_variant(interval=iv, variant=var, variant_scorers=SPLICE)
            df = variant_scorers.tidy_scores(out)
            if df is None or len(df) == 0 or "raw_score" not in df.columns:
                return v, ""  # scored but no splice rows -> empty (rare)
            return v, f"{df['raw_score'].abs().max():.6f}"
        except Exception as e:
            last = str(e)[:80]
            time.sleep(2 * (attempt + 1))
    return v, None  # give up -> NA (leave uncached so a rerun retries)


with ThreadPoolExecutor(max_workers=N_WORKERS) as ex:
    futs = {ex.submit(score_one, v): v for v in todo}
    for fut in as_completed(futs):
        v, val = fut.result()
        if val is None:
            continue  # give-up: don't cache, rerun will retry
        with lock:
            with open(CACHE, "a") as f:
                f.write(f"{v[0]}\t{v[1]}\t{v[2]}\t{v[3]}\t{val}\n")
            cnt[0] += 1
            if cnt[0] % 200 == 0:
                print(f"  scored {cnt[0]}/{len(todo)} new", flush=True)

print(f"AlphaGenome done: cache now has "
      f"{sum(1 for _ in open(CACHE)) - 1} records", flush=True)
