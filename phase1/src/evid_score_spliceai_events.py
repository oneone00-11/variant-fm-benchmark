"""E6 support -- SpliceAI's four delta scores AND their four positions.

The atlas persisted only the maximum delta score, and the E0 re-score kept the four
components but not the positions. E6 needs the positions: whether a predicted
splicing change is in-frame depends on WHERE the lost or gained site is, not only on
how strong the prediction is. Rather than re-score everything, this scores the
subset E6 asks about -- variants the tool calls and the assay does not -- and emits
the full SpliceAI record.

Output columns follow SpliceAI's own field order:
  DS_AG DS_AL DS_DG DS_DL  (delta scores)
  DP_AG DP_AL DP_DG DP_DL  (positions, relative to the variant)

Runs in the atlas's pinned SpliceAI environment:
    <atlas>/models/spliceai/.venv/bin/python phase1/src/evid_score_spliceai_events.py \
        --variants <parquet> --distance 4999 --out <parquet> --run
"""
from __future__ import annotations

import argparse
import os
import hashlib
import inspect
import json
import re
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

# the companion atlas: EVID_ATLAS_REPO, or a checkout beside this repository,
# the rule evid_common uses (this script runs in the atlas SpliceAI environment
# and cannot import it)
ATLAS = Path(os.environ.get(
    "EVID_ATLAS_REPO",
    str(Path(__file__).resolve().parents[2].parent / "functional-standard-atlas")))
DEFAULT_FASTA = ATLAS / "data" / "refs" / "grch38_subset.fa"
ANNOTATION = "grch38"
MASK = 0

SCORE_FIELDS = ["ds_ag", "ds_al", "ds_dg", "ds_dl"]
POS_FIELDS = ["dp_ag", "dp_al", "dp_dg", "dp_dl"]


def build_fullprec():
    from spliceai import utils as U

    src = inspect.getsource(U.get_delta_scores)
    patched, n = re.subn(r"\{:\.2f\}", "{:.17g}", src)
    if n != 4:
        raise RuntimeError(f"expected 4 '{{:.2f}}' fields to patch, found {n}")
    patched = patched.replace("def get_delta_scores(", "def get_delta_scores_fullprec(", 1)
    ns = dict(vars(U))
    exec(compile(patched, "<spliceai-fullprec>", "exec"), ns)
    return ns["get_delta_scores_fullprec"]


class Rec:
    def __init__(self, chrom, pos, ref, alt):
        self.chrom, self.pos, self.ref, self.alts = chrom, pos, ref, [alt]


def parse(fields: list[str]):
    """Full record at the gene whose maximum delta is the largest."""
    best = np.nan
    out = (np.nan,) * 8
    symbol = ""
    for f in fields:
        p = f.split("|")
        if len(p) < 10 or p[2] == ".":
            continue
        ds = [float(x) for x in p[2:6]]
        dp = [float(x) for x in p[6:10]]
        m = max(ds)
        if np.isnan(best) or m > best:
            best, out, symbol = m, tuple(ds) + tuple(dp), p[1]
    return best, out, symbol


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--variants", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--fasta", default=str(DEFAULT_FASTA))
    ap.add_argument("--distance", type=int, default=4999)
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--report-every", type=int, default=100)
    args = ap.parse_args()
    if not args.run:
        ap.error("pass --run")

    from spliceai.utils import Annotator
    import spliceai

    variants = pd.read_parquet(args.variants)
    variants["chrom_s"] = variants["chrom"].astype(str)
    ann = Annotator(args.fasta, ANNOTATION)
    fullprec = build_fullprec()

    rows, symbols, maxes, t0 = [], [], [], time.time()
    for i, r in enumerate(variants.itertuples(), 1):
        rec = Rec(r.chrom_s, int(r.pos), r.ref, r.alt)
        try:
            m, rec_out, sym = parse(fullprec(rec, ann, args.distance, MASK))
        except Exception:
            m, rec_out, sym = np.nan, (np.nan,) * 8, ""
        rows.append(rec_out)
        maxes.append(m)
        symbols.append(sym)
        if i % args.report_every == 0:
            el = time.time() - t0
            print(f"{i:>6}/{len(variants)}  {el/60:5.1f} min elapsed, "
                  f"{(len(variants)-i)*el/i/60:5.1f} min remaining", flush=True)

    arr = np.asarray(rows, dtype=float)
    out = variants[["variant_id", "gene", "chrom", "pos", "ref", "alt"]].copy()
    out["spliceai_event_max"] = maxes
    out["spliceai_event_symbol"] = symbols
    for j, f in enumerate(SCORE_FIELDS + POS_FIELDS):
        out[f"spliceai_{f}"] = arr[:, j]

    path = Path(args.out)
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(path, index=False)
    prov = {
        "product": "SpliceAI full record (four deltas and four positions)",
        "package_version": getattr(spliceai, "__version__", "1.3.1"),
        "distance": args.distance, "mask": MASK,
        "reference_fasta": args.fasta, "annotation": ANNOTATION,
        "n_variants": int(len(out)),
        "n_scored": int(out["spliceai_event_max"].notna().sum()),
        "output_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "scored_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "runtime_minutes": round((time.time() - t0) / 60, 1),
    }
    Path(str(path) + ".provenance.json").write_text(json.dumps(prov, indent=2) + "\n")
    print(json.dumps(prov, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
