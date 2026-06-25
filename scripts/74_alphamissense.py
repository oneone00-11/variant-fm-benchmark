"""M5-B: AlphaMissense (hg38) pathogenicity for MISSENSE variants only.
Downloads the Zenodo tsv.gz to data/raw/scores and stream-filters to the 7
gene coordinate ranges, then extracts am_pathogenicity for our (chrom,pos,ref,alt).
Non-missense variants are NOT scored (left absent — no extrapolation).
Directionality: higher am_pathogenicity = more likely pathogenic.
"""
import gzip
import sys

import pandas as pd
import requests
from tqdm import tqdm

from config import DATA_PROCESSED, DATA_RAW

URL = "https://zenodo.org/records/8208688/files/AlphaMissense_hg38.tsv.gz"
SCORES = DATA_RAW / "scores"; SCORES.mkdir(parents=True, exist_ok=True)
GZ = SCORES / "AlphaMissense_hg38.tsv.gz"


def download():
    if GZ.exists() and GZ.stat().st_size > 5e8:
        print(f"already have {GZ} ({GZ.stat().st_size/1e6:.0f} MB)"); return
    s = requests.Session(); s.trust_env = False
    with s.get(URL, stream=True, timeout=120) as r:
        r.raise_for_status()
        total = int(r.headers.get("Content-Length", 0))
        tmp = GZ.with_suffix(".part")
        with open(tmp, "wb") as f, tqdm(total=total, unit="B", unit_scale=True) as bar:
            for b in r.iter_content(1 << 20):
                f.write(b); bar.update(len(b))
        tmp.replace(GZ)


def main():
    df = pd.read_csv(DATA_PROCESSED / "analysis_ready_v2.tsv", sep="\t",
                     dtype={"chrom": str}, low_memory=False)
    df["pos"] = df["pos"].astype(int)
    # build per-chrom coordinate window + exact key set
    chroms = {f"chr{c}" for c in df.chrom.unique()}
    bounds = {f"chr{c}": (g.pos.min(), g.pos.max()) for c, g in df.groupby("chrom")}
    want = set(df.chrom + ":" + df.pos.astype(str) + ":" + df.ref + ":" + df.alt)
    print(f"target chroms {chroms}; {len(want)} variant keys", flush=True)

    download()
    print("streaming AlphaMissense (filtering to gene windows)...", flush=True)
    rows = []
    with gzip.open(GZ, "rt") as f:
        for line in f:
            if line.startswith("#") or line.startswith("CHROM"):
                continue
            p = line.rstrip("\n").split("\t")
            # CHROM POS REF ALT genome uniprot transcript protein_variant am_path am_class
            chrom = p[0]
            if chrom not in chroms:
                continue
            pos = int(p[1])
            lo, hi = bounds[chrom]
            if pos < lo or pos > hi:
                continue
            c = chrom[3:]
            key = f"{c}:{pos}:{p[2]}:{p[3]}"
            if key in want:
                rows.append({"chrom": c, "pos": pos, "ref": p[2], "alt": p[3],
                             "alphamissense": float(p[8]), "am_class": p[9]})
    out = pd.DataFrame(rows).drop_duplicates(["chrom", "pos", "ref", "alt"])
    out.to_csv(SCORES / "alphamissense.tsv", sep="\t", index=False)
    print(f"Wrote {SCORES/'alphamissense.tsv'}: {len(out)} missense scored", flush=True)


if __name__ == "__main__":
    sys.exit(main())
