"""For unmapped (no postMapped) intronic/splice variants, check whether the
current MappedVariant record still carries a clingenAlleleId we could resolve
to GRCh38 via the ClinGen Allele Registry."""
import json
from collections import Counter

import requests
from config import DATA_RAW

CACHE = DATA_RAW / "mavedb" / "mapped"
S = requests.Session(); S.trust_env = False

for urn, gene in [("urn_mavedb_00000097-0-1", "BRCA1?"),
                  ("urn_mavedb_00000097-0-2", "BRCA1"),
                  ("urn_mavedb_00000673-0-1", "RAD51C")]:
    p = CACHE / f"{urn}.json"
    if not p.exists():
        continue
    data = json.loads(p.read_text(encoding="utf-8"))
    cur = [m for m in data if m.get("current")]
    no_post = [m for m in cur if not m.get("postMapped")]
    with_caid = [m for m in no_post if m.get("clingenAlleleId")]
    print(f"\n=== {gene} {urn} ===")
    print(f"  current: {len(cur)}  | current w/o postMapped: {len(no_post)}"
          f"  | of those WITH clingenAlleleId: {len(with_caid)}")
    # sample a few unmapped to see variantUrn + caid + errorMessage
    for m in no_post[:4]:
        print(f"    {m['variantUrn']}  caid={m.get('clingenAlleleId')}"
              f"  err={(m.get('errorMessage') or '')[:50]}")

# Probe ClinGen Allele Registry for one CAID to confirm coords are retrievable
print("\n=== ClinGen Allele Registry probe ===")
for caid in ["CA500142878"]:
    for url in [f"https://reg.clinicalgenome.org/allele/{caid}"]:
        try:
            r = S.get(url, timeout=40)
            print(f"[{r.status_code}] {url}  len={len(r.content)}")
            if r.status_code == 200:
                d = r.json()
                # find GRCh38 genomic HGVS
                ga = d.get("genomicAlleles", [])
                for g in ga:
                    refg = g.get("referenceGenome")
                    hg = g.get("hgvs", [])
                    if refg == "GRCh38":
                        print("   GRCh38 hgvs:", hg)
        except Exception as e:
            print(f"[ERR] {url}: {type(e).__name__}: {str(e)[:100]}")
