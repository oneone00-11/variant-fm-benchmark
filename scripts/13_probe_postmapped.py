"""Find a real (current, non-placeholder) postMapped VRS allele and dump it,
so we learn how to extract GRCh38 chrom/pos/ref/alt."""
import json
import requests

S = requests.Session()
S.trust_env = False
URN = "urn:mavedb:00000097-0-2"
data = S.get(f"https://api.mavedb.org/api/v1/score-sets/{URN}/mapped-variants",
             timeout=120).json()
print("total mapped records:", len(data))
current = [m for m in data if m.get("current")]
print("current=True records:", len(current))
with_post = [m for m in current if m.get("postMapped")]
print("current & postMapped not null:", len(with_post))

if with_post:
    m = with_post[0]
    print("\n=== variantUrn:", m["variantUrn"], "clingen:", m.get("clingenAlleleId"))
    print("=== postMapped full ===")
    print(json.dumps(m["postMapped"], indent=2, ensure_ascii=False)[:3000])
    print("\n=== preMapped full ===")
    print(json.dumps(m.get("preMapped"), indent=2, ensure_ascii=False)[:1500])
