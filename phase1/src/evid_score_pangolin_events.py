"""E6 support -- Pangolin's gain and loss scores AND their positions.

Same gap as for SpliceAI: the atlas keeps the aggregate and discards where the
predicted change is, and E6 needs where. Pangolin's own output carries it --
each transcript record is `gene|<pos>:<increase>|<pos>:<decrease>|Warnings`, with
the positions relative to the variant -- so this runs the same pinned CLI the atlas
runs, with a parser that keeps all four numbers instead of two.

Only the seven frozen genes are in scope. Pangolin needs a gffutils database of the
annotation, and the one the atlas built covers chromosomes 2, 3, 13, 16 and 17;
DDX3X is on chrX and would need a new database, which is a separate job.

Runs in the atlas's pinned Pangolin environment:
    <atlas>/models/pangolin/.venv/bin/python phase1/src/evid_score_pangolin_events.py \
        --variants <parquet> --out <parquet> --run
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ATLAS = Path(os.environ.get("evidence-strength_ATLAS_REPO",
                            "/Users/cliffzhang/work/functional-standard-atlas"))
REF_FASTA = ATLAS / "data" / "refs" / "grch38_subset.fa"
GTF_DB = ATLAS / "data" / "refs" / "grch38_subset.gtf.db"
DISTANCE = 50           # the atlas column's setting; Pangolin has no Walker basis

# gene|<gain_pos>:<gain>|<loss_pos>:<loss>|Warnings
PANG_RE = re.compile(r"\|(-?\d+):(-?[0-9.]+)\|(-?\d+):(-?[0-9.]+)\|Warnings")


def emit_csv(df: pd.DataFrame, path: Path) -> int:
    snv = df[(df["ref"].str.len() == 1) & (df["alt"].str.len() == 1)].copy()
    snv["chrom_s"] = snv["chrom"].astype(str)
    snv = snv.sort_values(["chrom_s", "pos"])
    with open(path, "w") as fh:
        fh.write("CHROM,POS,REF,ALT\n")
        for r in snv.itertuples():
            fh.write(f"{r.chrom_s},{int(r.pos)},{r.ref},{r.alt}\n")
    return len(snv)


def parse_with_positions(path: Path) -> dict:
    """(chrom,pos,ref,alt) -> (score, gain, gain_pos, loss, loss_pos).

    The winning transcript record is the one with the largest max(gain, -loss),
    matching the atlas's aggregate, so the positions describe the same record the
    score comes from.
    """
    out: dict = {}
    with open(path, encoding="utf-8") as fh:
        fh.readline()
        for line in fh:
            parts = line.rstrip("\n").split(",", 4)
            if len(parts) < 5:
                continue
            chrom, pos, ref, alt, val = parts
            best = None
            rec = None
            for g_pos, gain, l_pos, loss in PANG_RE.findall(val):
                gn, ls = float(gain), float(loss)
                agg = max(gn, -ls)
                if best is None or agg > best:
                    best = agg
                    rec = (agg, gn, int(g_pos), ls, int(l_pos))
            if rec is not None:
                out[(chrom, int(pos), ref, alt)] = rec
    return out


def run_chunks(csv_path: Path, chunks_dir: Path, n_chunks: int) -> dict:
    chunks_dir.mkdir(parents=True, exist_ok=True)
    lines = csv_path.read_text().splitlines()
    header, body = lines[0], lines[1:]
    size = math.ceil(len(body) / n_chunks)
    chunk_csvs = []
    for i in range(n_chunks):
        part = body[i * size:(i + 1) * size]
        if not part:
            continue
        cc = chunks_dir / f"chunk_{i:02d}.csv"
        cc.write_text("\n".join([header] + part) + "\n")
        chunk_csvs.append((i, cc))

    exe = Path(sys.executable).parent / "pangolin"
    env = dict(os.environ)
    env.update({"OMP_NUM_THREADS": "2", "MKL_NUM_THREADS": "2"})

    def one(item) -> Path:
        i, cc = item
        out = cc.with_suffix(".out.csv")
        marker = cc.with_suffix(".done")
        want = hashlib.md5(cc.read_bytes()).hexdigest()
        if marker.exists() and out.exists() and marker.read_text() == want:
            return out
        marker.unlink(missing_ok=True)
        time.sleep((i % 5) * 10)      # stagger torch model loads
        cmd = [str(exe), "-d", str(DISTANCE), str(cc), str(REF_FASTA),
               str(GTF_DB), str(cc.with_suffix(".out"))]
        proc = subprocess.run(cmd, capture_output=True, text=True, env=env)
        if proc.returncode != 0:
            raise RuntimeError(f"pangolin failed on {cc.name}: {proc.stderr[-400:]}")
        marker.write_text(want)
        print(f"  chunk {i} done", flush=True)
        return out

    with ThreadPoolExecutor(max_workers=min(len(chunk_csvs), 5)) as ex:
        outs = list(ex.map(one, chunk_csvs))
    scores: dict = {}
    for o in outs:
        scores.update(parse_with_positions(o))
    return scores


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--variants", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--work", default="data/evidence/pangolin_events_work")
    ap.add_argument("--chunks", type=int, default=5)
    ap.add_argument("--run", action="store_true")
    args = ap.parse_args()
    if not args.run:
        ap.error("pass --run")

    for p in (REF_FASTA, GTF_DB):
        if not p.exists():
            raise SystemExit(f"missing reference: {p}")

    df = pd.read_parquet(args.variants)
    work = Path(args.work)
    work.mkdir(parents=True, exist_ok=True)
    csv_path = work / "variants.csv"
    n = emit_csv(df, csv_path)
    print(f"[pangolin-events] {n} SNVs, {args.chunks} chunks, -d {DISTANCE}", flush=True)
    t0 = time.time()
    scores = run_chunks(csv_path, work / "chunks", args.chunks)

    key = list(zip(df["chrom"].astype(str), df["pos"].astype(int), df["ref"], df["alt"]))
    hit = [scores.get(k) for k in key]
    out = df[["variant_id", "gene", "chrom", "pos", "ref", "alt"]].copy()
    out["pangolin_event_max"] = [h[0] if h else np.nan for h in hit]
    out["pangolin_gain"] = [h[1] if h else np.nan for h in hit]
    out["pangolin_gain_pos"] = [h[2] if h else np.nan for h in hit]
    out["pangolin_loss"] = [h[3] if h else np.nan for h in hit]
    out["pangolin_loss_pos"] = [h[4] if h else np.nan for h in hit]

    path = Path(args.out)
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(path, index=False)
    prov = {
        "product": "Pangolin gain/loss scores with their positions",
        "distance": DISTANCE,
        "reference_fasta": str(REF_FASTA), "gtf_db": str(GTF_DB),
        "n_variants": int(len(out)),
        "n_scored": int(out["pangolin_event_max"].notna().sum()),
        "output_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "scored_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "runtime_minutes": round((time.time() - t0) / 60, 1),
    }
    Path(str(path) + ".provenance.json").write_text(json.dumps(prov, indent=2) + "\n")
    print(json.dumps(prov, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
