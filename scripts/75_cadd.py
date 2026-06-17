"""M5-B: CADD PHRED scores for all SNVs in analysis_ready_v2 via the CADD API
(GRCh38-v1.7). All regions. Cached + restartable. Directionality: higher = more
deleterious.
"""
import json
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
import requests

from config import DATA_PROCESSED, DATA_RAW

SCORES = DATA_RAW / "scores"; SCORES.mkdir(parents=True, exist_ok=True)
CACHE = SCORES / "cadd_cache.json"
API = "https://cadd.gs.washington.edu/api/v1.0/GRCh38-v1.7/{c}:{p}_{r}_{a}"

_local = threading.local()
lock = threading.Lock()


def sess():
    if not hasattr(_local, "s"):
        s = requests.Session(); s.trust_env = False
        _local.s = s
    return _local.s


def fetch(c, p, r, a, retries=3):
    url = API.format(c=c, p=p, r=r, a=a)
    for i in range(retries):
        try:
            resp = sess().get(url, timeout=40)
            if resp.status_code == 200:
                data = resp.json()
                if isinstance(data, list) and data:
                    # may return multiple annotations; take max PHRED
                    phred = max(float(d["PHRED"]) for d in data if d.get("PHRED"))
                    raw = max(float(d["RawScore"]) for d in data if d.get("RawScore"))
                    return phred, raw
                return None  # empty = not scored
            if resp.status_code in (429, 500, 502, 503):
                import time; time.sleep(1.5 * (i + 1)); continue
            return None
        except Exception:
            import time; time.sleep(1.0 * (i + 1))
    return "ERR"


def main():
    df = pd.read_csv(DATA_PROCESSED / "analysis_ready_v2.tsv", sep="\t",
                     dtype={"chrom": str}, low_memory=False)
    variants = df[["chrom", "pos", "ref", "alt"]].drop_duplicates()
    variants["pos"] = variants["pos"].astype(int)
    cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}
    todo = [(str(r.chrom), int(r.pos), r.ref, r.alt) for r in variants.itertuples()
            if f"{r.chrom}:{r.pos}:{r.ref}:{r.alt}" not in cache]
    print(f"{len(variants)} unique variants; {len(todo)} to fetch, "
          f"{len(cache)} cached", flush=True)
    done = [0]; err = [0]

    def task(v):
        c, p, r, a = v
        res = fetch(c, p, r, a)
        key = f"{c}:{p}:{r}:{a}"
        with lock:
            if res == "ERR":
                err[0] += 1
            else:
                cache[key] = res  # (phred,raw) or None
            done[0] += 1
            if done[0] % 500 == 0:
                CACHE.write_text(json.dumps(cache))
                print(f"  {done[0]}/{len(todo)} (err {err[0]})", flush=True)

    with ThreadPoolExecutor(max_workers=5) as ex:
        futs = [ex.submit(task, v) for v in todo]
        for f in as_completed(futs):
            f.result()
    CACHE.write_text(json.dumps(cache))

    rows = []
    for k, v in cache.items():
        c, p, r, a = k.split(":")
        phred, raw = (v if v else (None, None))
        rows.append({"chrom": c, "pos": int(p), "ref": r, "alt": a,
                     "cadd_phred": phred, "cadd_raw": raw})
    out = pd.DataFrame(rows)
    out.to_csv(SCORES / "cadd.tsv", sep="\t", index=False)
    print(f"Wrote {SCORES/'cadd.tsv'}: {out.cadd_phred.notna().sum()} scored "
          f"of {len(out)} (give-up errors {err[0]})", flush=True)


if __name__ == "__main__":
    main()
