"""Probe a single score set's detail + scores + mapped-variants endpoints.

Goal: learn (a) how targetGenes encodes the gene + transcript + reference genome,
(b) how to pull the per-variant scores, (c) whether GRCh38 post-mapped genomic
coordinates are available and in what shape.

Uses BRCA1 Findlay SGE function-score set urn:mavedb:00000097-0-2.
"""
import json
import requests

S = requests.Session()
S.trust_env = False
URN = "urn:mavedb:00000097-0-2"
BASE = "https://api.mavedb.org/api/v1"


def show(name, r):
    ct = r.headers.get("content-type", "")
    print(f"[{r.status_code}] {name}  len={len(r.content)}  ctype={ct[:30]}")
    return r


# 1) Score set detail
r = show("GET score-sets/{urn}", S.get(f"{BASE}/score-sets/{URN}", timeout=60))
if r.status_code == 200:
    d = r.json()
    print("  top keys:", list(d.keys()))
    tgs = d.get("targetGenes", [])
    print("\n  === targetGenes[0] full ===")
    if tgs:
        print(json.dumps(tgs[0], indent=2, ensure_ascii=False)[:2500])

# 2) Scores CSV
r = show("\nGET score-sets/{urn}/scores",
         S.get(f"{BASE}/score-sets/{URN}/scores", timeout=120))
if r.status_code == 200:
    text = r.text
    lines = text.splitlines()
    print("  scores CSV header:", lines[0] if lines else "(empty)")
    for ln in lines[1:4]:
        print("   ", ln)
    print(f"  total score rows: {len(lines)-1}")

# 3) Mapped variants (genomic coordinates)
for ep in ("mapped-variants", "variants", "mapped-variants/"):
    r = S.get(f"{BASE}/score-sets/{URN}/{ep}", timeout=120)
    print(f"\n[{r.status_code}] GET score-sets/{{urn}}/{ep}  len={len(r.content)}")
    if r.status_code == 200:
        try:
            data = r.json()
            if isinstance(data, list) and data:
                print("  list len:", len(data), "| first item keys:", list(data[0].keys()))
                print(json.dumps(data[0], indent=2, ensure_ascii=False)[:2500])
            elif isinstance(data, dict):
                print("  dict keys:", list(data.keys()))
                print(json.dumps(data, indent=2, ensure_ascii=False)[:2000])
        except Exception as e:
            print("  (non-json)", str(e)[:80], r.text[:200])
        break
