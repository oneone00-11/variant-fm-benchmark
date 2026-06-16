"""Test Mutalyzer with chromosomal-reference + transcript-selector form, which
is the documented way to map transcript c. (incl. intronic) -> genomic g."""
import json
import requests

S = requests.Session(); S.trust_env = False

tests = [
    "NC_000017.11(NM_007294.3):c.5467+20C>A",   # BRCA1 intronic
    "NC_000013.11(ENST00000380152.8):c.7436-10T>A",  # BRCA2 intronic (Ensembl)
]
for hgvs in tests:
    print(f"\n===== {hgvs} =====")
    try:
        r = S.get(f"https://mutalyzer.nl/api/normalize/{hgvs}", timeout=90)
        print(f"[{r.status_code}] len={len(r.content)}")
        d = r.json()
        print("top keys:", list(d.keys()))
        for k in ("normalized_description", "genomic_description"):
            if d.get(k):
                print(f"  {k}: {d[k]}")
        eq = d.get("equivalent_descriptions") or d.get("chromosomal_descriptions")
        if eq:
            print("  equivalents:", json.dumps(eq)[:400])
        if "errors" in d:
            print("  errors:", json.dumps(d["errors"])[:300])
        if "messages" in d:
            print("  messages:", json.dumps(d["messages"])[:300])
    except Exception as e:
        print(f"[ERR] {type(e).__name__}: {str(e)[:150]}")
