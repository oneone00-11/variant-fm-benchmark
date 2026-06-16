"""Inspect the shape of the MaveDB score-set search response."""
import json
import requests

S = requests.Session()
S.trust_env = False

r = S.post("https://api.mavedb.org/api/v1/score-sets/search",
           json={"text": "BRCA1"}, timeout=60)
data = r.json()
print("top-level type:", type(data).__name__)
if isinstance(data, dict):
    print("top-level keys:", list(data.keys()))
    # find the list of results
    for k, v in data.items():
        print(f"  {k}: {type(v).__name__}"
              + (f" len={len(v)}" if hasattr(v, '__len__') else ""))
    items = data.get("scoreSets") or data.get("data") or data.get("items")
elif isinstance(data, list):
    items = data
    print("list len:", len(items))
else:
    items = None

if items:
    print("\n=== first item keys ===")
    print(list(items[0].keys()))
    print("\n=== full first item ===")
    print(json.dumps(items[0], indent=2, ensure_ascii=False)[:4000])
    # Show how target gene names appear across all items
    print("\n=== targetGenes summary across results ===")
    for it in items[:62]:
        tgs = it.get("targetGenes") or it.get("targetGene") or []
        names = [tg.get("name") for tg in tgs] if isinstance(tgs, list) else tgs
        print(f"  {it.get('urn')}  targets={names}  nVariants={it.get('numVariants')}")
