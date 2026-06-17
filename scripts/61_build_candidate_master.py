"""M4: build candidate functional master (VHL, BAP1) — download SGE scores +
post-mapped coords, then Mutalyzer-map the splice/intronic SNVs. Same pipeline
as Milestone 2 (scripts 21+28). Output: data/processed/candidate_functional_master.tsv
(separate from analysis_ready; nothing existing is touched).
"""
import csv
import io
import json
import re
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
import requests

from config import DATA_RAW, DATA_PROCESSED
from hgvs_region import classify_region
from mutalyzer_map import (build_query, ctoken_from_hgvs, load_cache, map_one,
                           save_cache, NC_CHROM)

BASE = "https://api.mavedb.org/api/v1"
CACHE = DATA_RAW / "mavedb"
(CACHE / "scores").mkdir(parents=True, exist_ok=True)
(CACHE / "mapped").mkdir(parents=True, exist_ok=True)
S = requests.Session(); S.trust_env = False

SELECTED = {
    "VHL":  "urn:mavedb:00000675-a-1",
    "BAP1": "urn:mavedb:00000662-0-1",
}
_GHGVS = re.compile(r"(NC_\d+)\.\d+:g\.(\d+)([ACGTN]+)>([ACGTN]+)")
SINGLE_SUB = re.compile(r":c\.[\*\-]?\d+(?:[+\-]\d+)?[ACGT]>[ACGT]$")


def fetch_text(url, path):
    if path.exists():
        return path.read_text(encoding="utf-8")
    r = S.get(url, timeout=180); r.raise_for_status()
    path.write_text(r.text, encoding="utf-8"); time.sleep(0.2)
    return r.text


def fetch_json(url, path):
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    r = S.get(url, timeout=240); r.raise_for_status()
    path.write_text(r.text, encoding="utf-8"); time.sleep(0.2)
    return r.json()


def parse_genomic(pm):
    if not pm:
        return None
    for e in pm.get("expressions", []):
        if e.get("syntax") == "hgvs.g":
            m = _GHGVS.match(e["value"])
            if m:
                nc, pos, ref, alt = m.groups()
                chrom = NC_CHROM.get(nc)
                if chrom and len(ref) == 1 and len(alt) == 1:
                    return chrom, int(pos), ref, alt, True
                if chrom:
                    return chrom, int(pos), "", "", False
    return None


def load_mapped(urn):
    data = fetch_json(f"{BASE}/score-sets/{urn}/mapped-variants",
                      CACHE / "mapped" / f"{urn.replace(':', '_')}.json")
    out = {}
    for m in data:
        if m.get("current") and m.get("postMapped"):
            g = parse_genomic(m["postMapped"])
            if g:
                out[m["variantUrn"]] = g
    return out


# --- 1. build rows from scores + postMapped ---
rows = []
for gene, urn in SELECTED.items():
    print(f"=== {gene} {urn} ===", flush=True)
    scores = fetch_text(f"{BASE}/score-sets/{urn}/scores",
                        CACHE / "scores" / f"{urn.replace(':', '_')}.csv")
    mapped = load_mapped(urn)
    n = npost = 0
    for r in csv.DictReader(io.StringIO(scores)):
        n += 1
        vurn = r["accession"]
        hgvs_nt = r.get("hgvs_nt", "") or ""
        region, off = classify_region(hgvs_nt, r.get("hgvs_pro", "") or "",
                                      r.get("hgvs_splice", "") or "")
        try:
            score = float(r["score"]) if r.get("score") not in (None, "", "NA") else None
        except ValueError:
            score = None
        g = mapped.get(vurn)
        if g:
            chrom, pos, ref, alt, is_snv = g; npost += 1
        else:
            chrom = pos = ref = alt = None; is_snv = False
        rows.append({"gene": gene, "source_urn": urn, "assay_type": "SGE",
                     "variant_urn": vurn, "hgvs_nt": hgvs_nt,
                     "hgvs_pro": r.get("hgvs_pro", ""), "region_class": region,
                     "intron_offset": off, "functional_score": score,
                     "functional_call": "", "chrom": chrom, "pos": pos,
                     "ref": ref, "alt": alt, "is_snv": is_snv,
                     "mapped": g is not None,
                     "mapping_source": "mavedb_postmapped" if g else ""})
    print(f"  scores {n} | postMapped {npost}", flush=True)

df = pd.DataFrame(rows)

# --- 2. Mutalyzer-map unmapped single-sub (splice/intronic) ---
todo = df[(~df["mapped"]) & (df["hgvs_nt"].str.contains(SINGLE_SUB, regex=True, na=False))]
work = todo[["gene", "hgvs_nt"]].drop_duplicates()
print(f"\nMutalyzer-mapping {len(work)} unique splice/intronic variants...", flush=True)
cache = load_cache()
lock = threading.Lock(); done = [0]; err = [0]


def task(gene, hgvs_nt):
    key = build_query(gene, ctoken_from_hgvs(hgvs_nt))
    if key in cache:
        return hgvs_nt, cache[key]
    try:
        res = map_one(gene, ctoken_from_hgvs(hgvs_nt))
    except Exception:
        with lock: err[0] += 1
        return hgvs_nt, None
    val = list(res) if res else None
    with lock:
        cache[key] = val; done[0] += 1
        if done[0] % 250 == 0:
            save_cache(cache); print(f"  mapped {done[0]} (err {err[0]})", flush=True)
    return hgvs_nt, val


results = {}
with ThreadPoolExecutor(max_workers=6) as ex:
    futs = [ex.submit(task, r.gene, r.hgvs_nt) for r in work.itertuples()]
    for f in as_completed(futs):
        h, v = f.result()
        if v: results[h] = v
save_cache(cache)
for r in work.itertuples():
    if r.hgvs_nt not in results:
        k = build_query(r.gene, ctoken_from_hgvs(r.hgvs_nt))
        if cache.get(k): results[r.hgvs_nt] = cache[k]
print(f"newly mapped {len(results)} (give-up {err[0]})", flush=True)


def apply_row(row):
    if (not row["mapped"]) and row["hgvs_nt"] in results:
        c, p, rf, al = results[row["hgvs_nt"]]
        row["chrom"], row["pos"], row["ref"], row["alt"] = c, p, rf, al
        row["is_snv"], row["mapped"], row["mapping_source"] = True, True, "mutalyzer"
    return row


df = df.apply(apply_row, axis=1)
df["chrom"] = df["chrom"].astype("string")
out = DATA_PROCESSED / "candidate_functional_master.tsv"
df.to_csv(out, sep="\t", index=False)
print(f"\nWrote {out} ({len(df)} variants)", flush=True)
print(df.groupby("gene").agg(n=("variant_urn", "count"),
                             mapped_snv=("is_snv", "sum")).to_string())
