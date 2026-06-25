"""AlphaGenome Step-0 gate: minimal 1-variant scoring test.
Key from env only; never printed. Run inside the alphagenome conda env.
"""
import os
import sys

key = os.environ.get("ALPHAGENOME_API_KEY")
if not key:
    sys.exit("ALPHAGENOME_API_KEY not set in environment - aborting (do not pass key as arg).")

from alphagenome.models import dna_client, variant_scorers
from alphagenome.data import genome

model = dna_client.create(key)
print("client created OK")

# list recommended variant scorers (look for splice-related)
recs = variant_scorers.RECOMMENDED_VARIANT_SCORERS
print("n recommended scorers:", len(recs))
for name in list(recs)[:40]:
    print("  scorer:", name)

# minimal test: BRCA1 c.5194-1G>C  -> chr17:43057136 C>G (AlphaGenome uses 'chr')
variant = genome.Variant(chromosome="chr17", position=43057136,
                         reference_bases="C", alternate_bases="G")
interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
# pick splice-related recommended scorers
splice_scorers = [s for n, s in recs.items()
                  if "SPLICE" in str(n).upper()] or list(recs.values())[:2]
print(f"\nscoring 1 variant with {len(splice_scorers)} scorer(s)...")
out = model.score_variant(interval=interval, variant=variant,
                          variant_scorers=splice_scorers)
print("returned", len(out), "score objects")
for sd in out:
    df = sd.values if hasattr(sd, "values") else None
    print("  scorer output type:", type(sd).__name__,
          "| metadata:", getattr(sd, "name", "?"))
print("TEST OK")
