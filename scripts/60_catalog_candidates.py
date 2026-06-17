"""M4-1: catalog MaveDB score sets for candidate panel-expansion genes.
Reuses the Milestone-2 catalog logic (off-target filtering via targetAccession.gene).
"""
import json
import re
import time

import pandas as pd
import requests

from config import RESULTS_TABLES, DATA_RAW

CANDIDATE_GENES = ["VHL", "BAP1", "MSH2", "TP53", "PTEN", "MLH1", "NF1"]
BASE = "https://api.mavedb.org/api/v1"
CACHE = DATA_RAW / "mavedb" / "details"
CACHE.mkdir(parents=True, exist_ok=True)
TARGET_SET = set(CANDIDATE_GENES)

S = requests.Session(); S.trust_env = False
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
    time.sleep(0.1)
    return d


def matched_genes(detail):
    hits = set()
    for tg in detail.get("targetGenes", []):
        names = {tg.get("name"), tg.get("mappedHgncName")}
        acc = tg.get("targetAccession") or {}
        names.add(acc.get("gene"))
        for nm in list(names):
            if nm:
                names.add(nm.split()[0])
        hits |= (names & TARGET_SET)
    return hits


def assay_type(detail):
    txt = " ".join([
        detail.get("title", ""),
        detail.get("experiment", {}).get("abstractText", "") or "",
        detail.get("experiment", {}).get("methodText", "") or "",
        detail.get("methodText", "") or "",
    ]).lower()
    is_sge = ("saturation genome editing" in txt or re.search(r"\bsge\b", txt) is not None)
    label = "SGE" if is_sge else ("VAMP-seq" if "vamp" in txt else "other-MAVE")
    return label, is_sge


def main():
    candidates = {}
    for gene in CANDIDATE_GENES:
        sets = search_gene(gene)
        print(f"text-search '{gene}': {len(sets)} candidates")
        for ss in sets:
            candidates[ss["urn"]] = ss
    print(f"\n{len(candidates)} unique URNs; filtering by targetAccession.gene...")

    rows = []
    for urn in sorted(candidates):
        try:
            d = get_detail(urn)
        except Exception as e:
            print(f"  detail failed {urn}: {e}")
            continue
        hits = matched_genes(d)
        if not hits:
            continue
        label, is_sge = assay_type(d)
        accs, assemblies = [], set()
        for tg in d.get("targetGenes", []):
            acc = tg.get("targetAccession") or {}
            if acc.get("accession"):
                accs.append(acc["accession"])
            if acc.get("assembly"):
                assemblies.add(acc["assembly"])
        pubs = d.get("primaryPublicationIdentifiers", []) or []
        rows.append({
            "urn": urn, "gene": "|".join(sorted(hits)), "title": d.get("title", ""),
            "assay_type": label, "is_sge": is_sge,
            "num_variants": d.get("numVariants", 0),
            "transcripts": "|".join(sorted(set(accs))),
            "assembly": "|".join(sorted(assemblies)),
            "pmid": pubs[0].get("identifier") if pubs else "",
            "license": (d.get("license") or {}).get("shortName", ""),
        })

    df = pd.DataFrame(rows).sort_values(["gene", "is_sge", "num_variants"],
                                        ascending=[True, False, False])
    out = RESULTS_TABLES / "candidate_mavedb_catalog.tsv"
    df.to_csv(out, sep="\t", index=False)
    pd.set_option("display.width", 200); pd.set_option("display.max_colwidth", 40)
    print(f"\nWrote {out} ({len(df)} on-target score sets)\n")
    print(df[["urn", "gene", "assay_type", "is_sge", "num_variants",
              "transcripts", "assembly", "pmid", "license"]].to_string(index=False))
    print("\n=== SGE availability per candidate ===")
    for g in CANDIDATE_GENES:
        sge = df[(df.gene == g) & (df.is_sge)]
        n = len(sge); mx = sge.num_variants.max() if n else 0
        print(f"  {g}: {n} SGE score set(s), max {mx} variants")


if __name__ == "__main__":
    main()
