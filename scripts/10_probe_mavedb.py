"""Probe current MaveDB API + ProteinGym access (Milestone 2, subtask 0).

Does NOT assume old endpoints — just checks reachability and shape of responses.
Run:  python scripts/10_probe_mavedb.py
"""
import json
import requests

S = requests.Session()
S.trust_env = False  # bypass broken system proxy
S.headers.update({"Accept": "application/json"})


def probe(name, url, method="GET", body=None):
    try:
        if method == "GET":
            r = S.get(url, timeout=40)
        else:
            r = S.post(url, json=body, timeout=40)
        ctype = r.headers.get("content-type", "")[:40]
        print(f"[{r.status_code}] {name}  len={len(r.content)}  ctype={ctype}")
        return r
    except Exception as e:
        print(f"[ERR] {name}: {type(e).__name__}: {str(e)[:140]}")
        return None


print("=== Reachability ===")
probe("mavedb_api_root", "https://api.mavedb.org/api/v1/")
probe("mavedb_site", "https://www.mavedb.org/")
r_openapi = probe("mavedb_openapi", "https://api.mavedb.org/api/v1/openapi.json")
probe("proteingym_site", "https://proteingym.org/")

# Enumerate available API paths from the OpenAPI spec, if we got it.
if r_openapi is not None and r_openapi.status_code == 200:
    try:
        spec = r_openapi.json()
        paths = sorted(spec.get("paths", {}).keys())
        print(f"\n=== MaveDB API paths ({len(paths)}) ===")
        for p in paths:
            # highlight the ones we care about
            if any(k in p for k in ("score", "target", "search", "gene", "experiment", "mapped", "variant")):
                methods = ",".join(m.upper() for m in spec["paths"][p].keys())
                print(f"  {p}   [{methods}]")
    except Exception as e:
        print(f"openapi parse error: {e}")

# Try a score-set search for BRCA1 to learn the response shape.
print("\n=== Try score-set search (BRCA1) ===")
for url, body in [
    ("https://api.mavedb.org/api/v1/score-sets/search", {"text": "BRCA1"}),
    ("https://api.mavedb.org/api/v1/score-sets/search", {"targets": ["BRCA1"]}),
]:
    r = probe(f"POST {url} body={body}", url, method="POST", body=body)
    if r is not None and r.status_code == 200:
        try:
            data = r.json()
            print(f"   -> returned {len(data)} items; first keys: "
                  f"{list(data[0].keys()) if data else 'none'}")
            if data:
                print("   sample item (trimmed):")
                print(json.dumps(data[0], indent=2)[:1500])
        except Exception as e:
            print(f"   json parse error: {e}")
        break
