import json, requests
S = requests.Session(); S.trust_env = False
Q = """
query($symbol: String!) {
  gene(gene_symbol: $symbol, reference_genome: GRCh38) {
    variants(dataset: gnomad_r4) {
      variant_id
      genome { af populations { id af } }
      exome  { af populations { id af } }
    }
  }
}"""
r = S.post("https://gnomad.broadinstitute.org/api",
           json={"query": Q, "variables": {"symbol": "BRCA1"}},
           headers={"Content-Type": "application/json"}, timeout=90)
print("status", r.status_code)
j = r.json()
print(json.dumps(j, indent=2)[:1500])
