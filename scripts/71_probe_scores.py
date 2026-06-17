"""M5-B: probe reachability + response shape of each scoring source from THIS
machine (China network, flaky proxy). Decide which are feasible before scoring.
"""
import json
import requests

S = requests.Session(); S.trust_env = False
S.headers.update({"User-Agent": "variant-fm-benchmark/0.1"})


def probe(name, method, url, **kw):
    try:
        r = S.request(method, url, timeout=40, **kw)
        body = r.text[:200].replace("\n", " ")
        print(f"[{r.status_code}] {name}  len={len(r.content)}  {body}")
        return r
    except Exception as e:
        print(f"[ERR] {name}: {type(e).__name__}: {str(e)[:120]}")
        return None


print("=== gnomAD GraphQL (per-gene AF) ===")
q = {"query": '{gene(gene_symbol:"VHL",reference_genome:GRCh38){gene_id chrom start stop}}'}
probe("gnomad_api", "POST", "https://gnomad.broadinstitute.org/api",
      json=q, headers={"Content-Type": "application/json"})

print("\n=== Broad SpliceAI / Pangolin lookup API ===")
probe("spliceai_lookup", "GET",
      "https://spliceailookup-api.broadinstitute.org/spliceai/?hg=38&distance=50&variant=17-43045705-T-A")
probe("pangolin_lookup", "GET",
      "https://spliceailookup-api.broadinstitute.org/pangolin/?hg=38&distance=50&variant=17-43045705-T-A")

print("\n=== CADD API ===")
for u in ["https://cadd.gs.washington.edu/api/v1.0/v1.7/17:43045705_T_A",
          "https://cadd.gs.washington.edu/api/v1.0/GRCh38-v1.7/17:43045705_T_A"]:
    probe(u.split("/api")[1], "GET", u)

print("\n=== UCSC track API (phyloP / GERP, per-region) ===")
probe("ucsc_phyloP100way", "GET",
      "https://api.genome.ucsc.edu/getData/track?genome=hg38;track=phyloP100way;chrom=chr17;start=43045700;end=43045710")
probe("ucsc_list_tracks", "GET",
      "https://api.genome.ucsc.edu/list/tracks?genome=hg38;trackLeavesOnly=1")

print("\n=== AlphaMissense (Zenodo hg38 tsv.gz) HEAD ===")
probe("alphamissense_zenodo", "HEAD",
      "https://zenodo.org/records/8208688/files/AlphaMissense_hg38.tsv.gz")

print("\n=== AlphaGenome API base ===")
probe("alphagenome", "GET", "https://www.alphagenomedocs.com/")
