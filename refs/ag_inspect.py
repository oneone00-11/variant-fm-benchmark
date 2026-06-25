"""Inspect AlphaGenome variant scoring output + time one call, to fix the
per-variant aggregation before the full run. Key from env only."""
import os, time, sys
import numpy as np
from alphagenome.models import dna_client, variant_scorers
from alphagenome.data import genome

model = dna_client.create(os.environ["ALPHAGENOME_API_KEY"])
recs = variant_scorers.RECOMMENDED_VARIANT_SCORERS
splice = [recs[n] for n in recs if "SPLICE" in str(n).upper()]
print("n splice scorers:", len(splice))

v = genome.Variant(chromosome="chr17", position=43057136,
                   reference_bases="C", alternate_bases="G")
for L, name in [(dna_client.SEQUENCE_LENGTH_100KB, "100KB"),
                (dna_client.SEQUENCE_LENGTH_1MB, "1MB")]:
    iv = v.reference_interval.resize(L)
    t = time.time()
    out = model.score_variant(interval=iv, variant=v, variant_scorers=splice)
    dt = time.time() - t
    df = variant_scorers.tidy_scores(out)
    print(f"\n=== interval {name}: {dt:.1f}s, tidy rows={len(df)} ===")
    print("columns:", list(df.columns))
    num = [c for c in df.columns if df[c].dtype.kind in "fc"]
    print("numeric cols:", num)
    sc = "raw_score" if "raw_score" in df.columns else (num[-1] if num else None)
    if sc:
        print(f"max |{sc}| = {df[sc].abs().max():.4f}")
    print(df.head(8).to_string())
