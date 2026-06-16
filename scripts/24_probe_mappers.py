"""Test transcript-HGVS -> GRCh38 mapping tools on real intronic/splice variants.
Pick the one that works reliably from here, then batch it."""
import json

import pandas as pd
import requests
from config import DATA_PROCESSED

S = requests.Session(); S.trust_env = False

df = pd.read_csv(DATA_PROCESSED / "functional_scores_master.tsv", sep="\t")
ex = (df[(df.gene == "BRCA1") & (df.region_class.isin(["splice", "intronic"])) & (~df.mapped)]
      ["hgvs_nt"].head(4).tolist())
print("example intronic/splice HGVS:", ex)

# --- Ensembl VEP REST (RefSeq HGVS) ---
print("\n=== Ensembl VEP REST /vep/human/hgvs ===")
try:
    r = S.post("https://rest.ensembl.org/vep/human/hgvs",
               headers={"Content-Type": "application/json", "Accept": "application/json"},
               json={"hgvs_notations": ex}, timeout=60)
    print(f"[{r.status_code}] len={len(r.content)}")
    if r.status_code == 200:
        d = r.json()
        for item in d[:2]:
            print("  input:", item.get("input"),
                  "| seq:", item.get("seq_region_name"),
                  "start:", item.get("start"),
                  "allele:", item.get("allele_string"),
                  "| conseq:", item.get("most_severe_consequence"))
    else:
        print("  body:", r.text[:300])
except Exception as e:
    print(f"[ERR] {type(e).__name__}: {str(e)[:150]}")

# --- Mutalyzer 3 normalize ---
print("\n=== Mutalyzer normalize ===")
for hgvs in ex[:1]:
    try:
        r = S.get(f"https://mutalyzer.nl/api/normalize/{hgvs}", timeout=60)
        print(f"[{r.status_code}] {hgvs}  len={len(r.content)}")
        if r.status_code == 200:
            d = r.json()
            print("  keys:", list(d.keys()))
            # look for chromosomal / genomic equivalents
            for k in ("chromosomal_descriptions", "equivalent_descriptions",
                      "genomic_descriptions", "normalized_description"):
                if k in d:
                    print(f"  {k}: {json.dumps(d[k])[:300]}")
        else:
            print("  body:", r.text[:300])
    except Exception as e:
        print(f"[ERR] {type(e).__name__}: {str(e)[:150]}")
