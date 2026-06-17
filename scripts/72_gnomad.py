"""M5-B: gnomAD allele frequency (global + popmax) per gene via the GraphQL API.
Naive 'common = likely benign' baseline. Directionality: higher AF = more likely
BENIGN (inverse of pathogenicity). gnomAD v4 (GRCh38).
"""
import json
import time

import pandas as pd
import requests

from config import DATA_PROCESSED, DATA_RAW

GENES = ["BRCA1", "BRCA2", "BARD1", "PALB2", "RAD51C", "BAP1", "VHL"]
DATASET = "gnomad_r4"
SCORES = DATA_RAW / "scores"; SCORES.mkdir(parents=True, exist_ok=True)
S = requests.Session(); S.trust_env = False

QUERY = """
query($symbol: String!) {
  gene(gene_symbol: $symbol, reference_genome: GRCh38) {
    variants(dataset: %s) {
      variant_id
      genome { af populations { id ac an } }
      exome  { af populations { id ac an } }
    }
  }
}
""" % DATASET

# continental population ids for popmax (exclude small/bottleneck groups)
POPMAX_IDS = {"afr", "amr", "eas", "nfe", "sas", "fin", "mid", "asj"}


def fetch_gene(symbol, retries=3):
    for i in range(retries):
        r = S.post("https://gnomad.broadinstitute.org/api",
                   json={"query": QUERY, "variables": {"symbol": symbol}},
                   headers={"Content-Type": "application/json"}, timeout=90)
        if r.status_code == 200:
            j = r.json()
            if j.get("data", {}).get("gene"):
                return j["data"]["gene"]["variants"]
            if j.get("errors"):
                raise RuntimeError(j["errors"][0].get("message", "gql error")[:120])
        time.sleep(2 * (i + 1))
    raise RuntimeError(f"gnomAD fetch failed for {symbol}")


def popmax(seq_obj):
    if not seq_obj:
        return None
    best = None
    for p in (seq_obj.get("populations") or []):
        if p["id"].lower() in POPMAX_IDS and p.get("an"):
            af = p["ac"] / p["an"] if p["an"] else None
            if af is not None:
                best = af if best is None else max(best, af)
    return best


def main():
    rows = []
    for g in GENES:
        vs = fetch_gene(g)
        print(f"{g}: {len(vs)} gnomAD variants", flush=True)
        for v in vs:
            c, p, ref, alt = v["variant_id"].split("-")
            gen, exo = v.get("genome"), v.get("exome")
            afs = [x.get("af") for x in (gen, exo) if x and x.get("af") is not None]
            af_global = max(afs) if afs else None
            pms = [x for x in (popmax(gen), popmax(exo)) if x is not None]
            af_popmax = max(pms) if pms else None
            rows.append({"chrom": c, "pos": int(p), "ref": ref, "alt": alt,
                         "gnomad_af_global": af_global, "gnomad_af_popmax": af_popmax})
        time.sleep(0.5)
    df = pd.DataFrame(rows).drop_duplicates(["chrom", "pos", "ref", "alt"])
    df.to_csv(SCORES / "gnomad_af.tsv", sep="\t", index=False)
    (SCORES / "gnomad_source.txt").write_text(
        f"dataset={DATASET}\napi=https://gnomad.broadinstitute.org/api\n", encoding="utf-8")
    print(f"Wrote {SCORES/'gnomad_af.tsv'}: {len(df)} gnomAD variants "
          f"({df.gnomad_af_global.notna().sum()} with global AF)")


if __name__ == "__main__":
    main()
