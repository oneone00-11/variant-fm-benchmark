"""Download the primary SGE score sets, extract GRCh38 coordinates from MaveDB
post-mapped VRS, classify region, and build functional_scores_master.tsv.

(Milestone 2, subtasks 2-3.)

Coordinates come straight from MaveDB's post-mapped genomic HGVS expression
(e.g. 'NC_000017.11:g.43045705T>A') -> no manual coordinate math. Variants
without a current post-mapped record are kept (score + region) but with empty
coords and mapped=False, so we can report how many (and which regions) lack
genomic placement.
"""
import csv
import io
import json
import re
import time

import pandas as pd
import requests

from config import DATA_RAW, DATA_PROCESSED
from hgvs_region import classify_region

BASE = "https://api.mavedb.org/api/v1"
CACHE = DATA_RAW / "mavedb"
(CACHE / "scores").mkdir(parents=True, exist_ok=True)
(CACHE / "mapped").mkdir(parents=True, exist_ok=True)

S = requests.Session()
S.trust_env = False

# Primary SGE function-score set per gene (ATM has none in MaveDB).
SELECTED = {
    "BRCA1":  "urn:mavedb:00000097-0-2",   # Findlay 2018 SGE
    "BARD1":  "urn:mavedb:00001250-a-2",   # BARD1 SGE
    "BRCA2":  "urn:mavedb:00001225-a-1",   # BRCA2 SGE (2025)
    "PALB2":  "urn:mavedb:00001259-a-2",   # PALB2 SGE
    "RAD51C": "urn:mavedb:00000673-0-1",   # RAD51C SGE
}
ASSAY = {g: "SGE" for g in SELECTED}

# RefSeq NC accession (GRCh38) -> ClinVar-style chrom label
NC_CHROM = {
    **{f"NC_0000{n:02d}": str(i) for i, n in zip(range(1, 23), range(1, 23))},
    "NC_000023": "X", "NC_000024": "Y", "NC_012920": "MT",
}
_GHGVS = re.compile(r"(NC_\d+)\.\d+:g\.(\d+)([ACGTN]+)>([ACGTN]+)")


def fetch_text(url, cache_path):
    if cache_path.exists():
        return cache_path.read_text(encoding="utf-8")
    r = S.get(url, timeout=180)
    r.raise_for_status()
    cache_path.write_text(r.text, encoding="utf-8")
    time.sleep(0.2)
    return r.text


def fetch_json(url, cache_path):
    if cache_path.exists():
        return json.loads(cache_path.read_text(encoding="utf-8"))
    r = S.get(url, timeout=240)
    r.raise_for_status()
    cache_path.write_text(r.text, encoding="utf-8")
    time.sleep(0.2)
    return r.json()


def parse_genomic(post_mapped):
    """Return (chrom, pos, ref, alt, is_snv) from a postMapped VRS allele,
    using its hgvs.g expression. Non-substitutions -> coords from location,
    ref/alt blank, is_snv=False."""
    if not post_mapped:
        return None
    for expr in post_mapped.get("expressions", []):
        if expr.get("syntax") == "hgvs.g":
            m = _GHGVS.match(expr["value"])
            if m:
                nc, pos, ref, alt = m.groups()
                chrom = NC_CHROM.get(nc)
                if chrom and len(ref) == 1 and len(alt) == 1:
                    return chrom, int(pos), ref, alt, True
                if chrom:  # indel/MNV substitution-like: record locus only
                    return chrom, int(pos), "", "", False
    # fallback: location only
    loc = post_mapped.get("location", {})
    sr = (loc.get("sequenceReference") or {}).get("label", "")
    nc = sr.split(".")[0]
    chrom = NC_CHROM.get(nc)
    if chrom and loc.get("start") is not None:
        return chrom, int(loc["start"]) + 1, "", "", False
    return None


def load_mapped(urn):
    """variant_urn -> (chrom,pos,ref,alt,is_snv) for current post-mapped records."""
    data = fetch_json(f"{BASE}/score-sets/{urn}/mapped-variants",
                      CACHE / "mapped" / f"{urn.replace(':', '_')}.json")
    out = {}
    n_current = n_post = 0
    for m in data:
        if not m.get("current"):
            continue
        n_current += 1
        pm = m.get("postMapped")
        if not pm:
            continue
        g = parse_genomic(pm)
        if g:
            n_post += 1
            out[m["variantUrn"]] = g
    return out, n_current, n_post


def main():
    rows = []
    stats = []
    for gene, urn in SELECTED.items():
        print(f"\n=== {gene}  {urn} ===")
        scores_csv = fetch_text(f"{BASE}/score-sets/{urn}/scores",
                                CACHE / "scores" / f"{urn.replace(':', '_')}.csv")
        mapped, n_cur, n_post = load_mapped(urn)
        reader = csv.DictReader(io.StringIO(scores_csv))
        n = n_mapped_snv = 0
        for r in reader:
            n += 1
            vurn = r["accession"]
            hgvs_nt = r.get("hgvs_nt", "") or ""
            hgvs_pro = r.get("hgvs_pro", "") or ""
            hgvs_splice = r.get("hgvs_splice", "") or ""
            region, off = classify_region(hgvs_nt, hgvs_pro, hgvs_splice)
            try:
                score = float(r["score"]) if r.get("score") not in (None, "", "NA") else None
            except ValueError:
                score = None
            g = mapped.get(vurn)
            if g:
                chrom, pos, ref, alt, is_snv = g
            else:
                chrom = pos = ref = alt = None
                is_snv = False
            if is_snv:
                n_mapped_snv += 1
            rows.append({
                "gene": gene, "source_urn": urn, "assay_type": ASSAY[gene],
                "variant_urn": vurn,
                "hgvs_nt": hgvs_nt, "hgvs_pro": hgvs_pro,
                "region_class": region, "intron_offset": off,
                "functional_score": score, "functional_call": "",
                "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                "is_snv": is_snv, "mapped": g is not None,
            })
        print(f"  scores rows: {n}  | current-mapped: {n_cur}  "
              f"post-mapped: {n_post}  | usable SNV coords: {n_mapped_snv}")
        stats.append({"gene": gene, "urn": urn, "n_variants": n,
                      "n_mapped_snv": n_mapped_snv})

    df = pd.DataFrame(rows)
    out = DATA_PROCESSED / "functional_scores_master.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"\nWrote {out}  ({len(df)} functional variants)")

    print("\n=== region_class x gene (all functional variants) ===")
    print(df.pivot_table(index="region_class", columns="gene",
                         values="variant_urn", aggfunc="count", fill_value=0)
          .to_string())

    print("\n=== mapping completeness (SNV coords) ===")
    print(pd.DataFrame(stats).to_string(index=False))


if __name__ == "__main__":
    main()
