"""Test NCBI Variation Services HGVS->SPDI (genomic) on intronic RefSeq HGVS.
Same provider as ClinVar -> coordinates will match our ClinVar master exactly."""
import urllib.parse
import requests

S = requests.Session(); S.trust_env = False
BASE = "https://api.ncbi.nlm.nih.gov/variation/v0"

tests = [
    "NM_007294.3:c.5467+20C>A",   # BRCA1 intronic
    "NM_007294.3:c.5333-1G>A",    # BRCA1 splice acceptor-ish
    "NM_000059.4:c.7436-10T>A",   # BRCA2 (RefSeq equiv of ENST MANE)
    "NM_058216.3:c.706-2A>G",     # RAD51C (RefSeq equiv)
]
for hgvs in tests:
    enc = urllib.parse.quote(hgvs, safe="")
    url = f"{BASE}/hgvs/{enc}/contextuals"
    try:
        r = S.get(url, timeout=60)
        print(f"\n[{r.status_code}] {hgvs}")
        if r.status_code == 200:
            d = r.json()
            spdis = d.get("data", {}).get("spdi", [])
            for s in spdis:
                print(f"   SPDI: seq={s['seq_id']} pos0={s['position']} "
                      f"del='{s['deleted_sequence']}' ins='{s['inserted_sequence']}'")
        else:
            print("   body:", r.text[:200])
    except Exception as e:
        print(f"[ERR] {type(e).__name__}: {str(e)[:120]}")
