"""Catalog MaveDB score sets targeting the 6 project genes (Milestone 2, subtask 1).

Strategy:
  - text-search each gene (fuzzy, returns off-target hits like BAP1/CHEK2/TP53)
  - fetch full detail for every unique candidate URN (cached to data/raw/mavedb/details)
  - KEEP only score sets where a targetGene's targetAccession.gene / mappedHgncName /
    name exactly matches one of our genes  -> drops off-target text matches
  - record URN, gene, transcript accession, assembly, assay type (flag SGE),
    numVariants, primary PMID + title, license

Output: results/tables/mavedb_catalog.tsv  (also prints a summary)
"""
import json
import re
import time

import pandas as pd
import requests

from config import TARGET_GENES, RESULTS_TABLES, DATA_RAW

BASE = "https://api.mavedb.org/api/v1"
CACHE = DATA_RAW / "mavedb" / "details"
CACHE.mkdir(parents=True, exist_ok=True)
TARGET_SET = set(TARGET_GENES)

S = requests.Session()
S.trust_env = False
S.headers.update({"Accept": "application/json"})


def search_gene(gene):
    r = S.post(f"{BASE}/score-sets/search", json={"text": gene}, timeout=60)
    r.raise_for_status()
    return r.json().get("scoreSets", [])


def get_detail(urn):
    cached = CACHE / f"{urn.replace(':', '_')}.json"
    if cached.exists():
        return json.loads(cached.read_text(encoding="utf-8"))
    r = S.get(f"{BASE}/score-sets/{urn}", timeout=60)
    r.raise_for_status()
    d = r.json()
    cached.write_text(json.dumps(d, ensure_ascii=False), encoding="utf-8")
    time.sleep(0.1)  # be polite
    return d


def matched_genes(detail):
    """Return the set of our target genes this score set actually targets."""
    hits = set()
    for tg in detail.get("targetGenes", []):
        names = {tg.get("name"), tg.get("mappedHgncName")}
        acc = tg.get("targetAccession") or {}
        names.add(acc.get("gene"))
        # also handle names like "BRCA1 RING domain"
        for nm in list(names):
            if nm:
                first = nm.split()[0]
                names.add(first)
        hits |= (names & TARGET_SET)
    return hits


def assay_type(detail):
    """Heuristic assay label + SGE flag from title/abstract/method text."""
    txt = " ".join([
        detail.get("title", ""),
        detail.get("experiment", {}).get("abstractText", "") or "",
        detail.get("experiment", {}).get("methodText", "") or "",
        detail.get("methodText", "") or "",
    ]).lower()
    is_sge = ("saturation genome editing" in txt
              or re.search(r"\bsge\b", txt) is not None)
    if is_sge:
        label = "SGE"
    elif "vamp-seq" in txt or "vampseq" in txt:
        label = "VAMP-seq"
    elif "yeast" in txt and ("two-hybrid" in txt or "complementation" in txt):
        label = "yeast-functional"
    elif "homology" in txt and "repair" in txt:
        label = "HDR"
    else:
        label = "other-MAVE"
    return label, is_sge


def main():
    candidates = {}  # urn -> summary
    for gene in TARGET_GENES:
        sets = search_gene(gene)
        print(f"text-search '{gene}': {len(sets)} candidate score sets")
        for ss in sets:
            candidates[ss["urn"]] = ss
    print(f"\n{len(candidates)} unique candidate URNs; fetching details + filtering...")

    rows = []
    for urn in sorted(candidates):
        try:
            d = get_detail(urn)
        except Exception as e:
            print(f"  detail failed {urn}: {e}")
            continue
        hits = matched_genes(d)
        if not hits:
            continue  # off-target text match
        label, is_sge = assay_type(d)
        accs = []
        assemblies = set()
        for tg in d.get("targetGenes", []):
            acc = tg.get("targetAccession") or {}
            if acc.get("accession"):
                accs.append(acc["accession"])
            if acc.get("assembly"):
                assemblies.add(acc["assembly"])
        pubs = d.get("primaryPublicationIdentifiers", []) or []
        pmid = pubs[0].get("identifier") if pubs else ""
        ptitle = pubs[0].get("title", "") if pubs else ""
        lic = (d.get("license") or {}).get("shortName", "")
        rows.append({
            "urn": urn,
            "gene": "|".join(sorted(hits)),
            "title": d.get("title", ""),
            "assay_type": label,
            "is_sge": is_sge,
            "num_variants": d.get("numVariants", 0),
            "transcripts": "|".join(sorted(set(accs))),
            "assembly": "|".join(sorted(assemblies)),
            "pmid": pmid,
            "pub_title": ptitle[:90],
            "license": lic,
            "mapping_state": d.get("mappingState", ""),
        })

    df = pd.DataFrame(rows).sort_values(["gene", "is_sge", "num_variants"],
                                        ascending=[True, False, False])
    out = RESULTS_TABLES / "mavedb_catalog.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"\nWrote {out}  ({len(df)} on-target score sets)")

    pd.set_option("display.width", 200)
    pd.set_option("display.max_colwidth", 45)
    print("\n=== On-target MaveDB score sets ===")
    print(df[["urn", "gene", "assay_type", "is_sge", "num_variants",
              "transcripts", "assembly", "pmid", "license"]].to_string(index=False))

    print("\n=== Counts per gene ===")
    print(df.groupby("gene").agg(
        n_scoresets=("urn", "count"),
        n_sge=("is_sge", "sum"),
        max_variants=("num_variants", "max"),
    ).to_string())


if __name__ == "__main__":
    main()
