"""M5-B: conservation scores (phyloP100way, phastCons100way) via the UCSC track
API, per gene region. Naive baselines. Directionality: higher = more conserved =
more likely pathogenic. GERP is not natively served for hg38 by the UCSC API
(documented in the report; deferred to an external bigWig if needed).
"""
import time

import pandas as pd
import requests

from config import DATA_PROCESSED, DATA_RAW

SCORES = DATA_RAW / "scores"; SCORES.mkdir(parents=True, exist_ok=True)
S = requests.Session(); S.trust_env = False
TRACKS = ["phyloP100way", "phastCons100way"]
API = "https://api.genome.ucsc.edu/getData/track"


def fetch_track(track, chrom, start, end, retries=3):
    """Return list of (start,end,value) intervals (0-based half-open)."""
    for i in range(retries):
        try:
            r = S.get(API, params={"genome": "hg38", "track": track,
                                   "chrom": f"chr{chrom}", "start": start, "end": end},
                      timeout=120)
            if r.status_code == 200:
                j = r.json()
                data = j.get(track)
                if isinstance(data, dict):
                    data = data.get(f"chr{chrom}", [])
                return data or []
            time.sleep(2 * (i + 1))
        except Exception:
            time.sleep(2 * (i + 1))
    return []


def main():
    df = pd.read_csv(DATA_PROCESSED / "analysis_ready_v2.tsv", sep="\t",
                     dtype={"chrom": str}, low_memory=False)
    df["pos"] = df["pos"].astype(int)
    out = df[["chrom", "pos", "ref", "alt"]].drop_duplicates().copy()

    for track in TRACKS:
        col = "phylop100way" if "phyloP" in track else "phastcons100way"
        pos2val = {}
        # group by GENE so each query is one gene's compact span (genes on the
        # same chromosome can be megabases apart -> chrom-wide span overflows the API)
        for gene, g in df.groupby("gene"):
            chrom = g.chrom.iloc[0]
            lo, hi = int(g.pos.min()), int(g.pos.max())
            intervals = fetch_track(track, chrom, lo - 1, hi)
            n = 0
            for it in intervals:
                v = it.get("value")
                s, e = it.get("start"), it.get("end")
                if v is None or s is None:
                    continue
                for p in range(int(s) + 1, int(e) + 1):  # 0-based half-open -> 1-based
                    pos2val[(chrom, p)] = v; n += 1
            print(f"{track} {gene} (chr{chrom} {lo}-{hi}): "
                  f"{len(intervals)} intervals -> {n} positions", flush=True)
            time.sleep(0.3)
        out[col] = out.apply(lambda r: pos2val.get((r["chrom"], r["pos"])), axis=1)
        print(f"  {col}: {out[col].notna().sum()}/{len(out)} positions scored", flush=True)

    out.to_csv(SCORES / "conservation.tsv", sep="\t", index=False)
    (SCORES / "conservation_source.txt").write_text(
        "tracks=phyloP100way,phastCons100way\napi=https://api.genome.ucsc.edu/getData/track\n"
        "genome=hg38\nGERP=not served natively for hg38 by UCSC API (deferred)\n",
        encoding="utf-8")
    print(f"Wrote {SCORES/'conservation.tsv'}")


if __name__ == "__main__":
    main()
