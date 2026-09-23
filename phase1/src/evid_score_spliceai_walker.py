"""Score SpliceAI under Walker et al. 2023's own settings, for the evidence-strength analysis set.

phase1/config/walker2023.yaml records why this column exists: the PP3/BP4 cut
points of 0.2 and 0.1 were calibrated on the maximum RAW delta score at a maximum
distance of 10,000 nt (+/-4,999), and the atlas column this project otherwise uses
was produced at -D 50. Version (1.3.1), statistic (max of the four deltas) and
masking (raw, -M 0) already agree; only the distance does not, so only the
distance changes here.

The scoring function is taken verbatim from the installed package with the four
``{:.2f}`` output fields rewritten to ``{:.17g}``, exactly as the atlas does in
``models/spliceai/score_fullprec.py``; ``--validate`` re-rounds the patched output
and requires it to reproduce the stock function field for field.

The four delta components are kept alongside the maximum, because E6 needs to know
which event type a call came from and the atlas persisted only the maximum.

Scoring 8,853 variants at -D 4999 takes about two and a half hours, so it
checkpoints: every --checkpoint-every variants the scores so far are written to
``<out>.partial.parquet``, and a re-run loads that file and scores only what is
missing. Killing the process costs at most one checkpoint interval, and the command
to resume is the same command.

Runs in the atlas's pinned SpliceAI environment, not this repo's .venv:

    <atlas>/models/spliceai/.venv/bin/python phase1/src/evid_score_spliceai_walker.py \
        --variants phase1/data/evidence/analysis_set_variants.parquet --validate
    ... --run
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
REF_FASTA = ATLAS / "data" / "refs" / "grch38_subset.fa"
ANNOTATION = "grch38"

# Walker et al. 2023, Methods ("Splicing prediction analysis"); see
# phase1/config/walker2023.yaml -> score_definition.
DISTANCE = 4999
MASK = 0

FIELDS = ["ds_ag", "ds_al", "ds_dg", "ds_dl"]


def build_fullprec():
    """`get_delta_scores` with only the output precision changed."""
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
    """Minimal stand-in for the pysam VCF record get_delta_scores expects."""

    def __init__(self, chrom, pos, ref, alt):
        self.chrom, self.pos, self.ref, self.alts = chrom, pos, ref, [alt]


def aggregate(fields: list[str]) -> tuple[float, float, float, float, float]:
    """(max, ds_ag, ds_al, ds_dg, ds_dl) over the annotated genes at this variant.

    The maximum is the atlas's and Walker's definition. The four components are
    reported at the gene whose maximum is the largest, so the component vector and
    the maximum always describe the same predicted event.
    """
    best = np.nan
    comp = (np.nan,) * 4
    for f in fields:
        parts = f.split("|")
        if len(parts) < 6 or parts[2] == ".":
            continue
        vals = [float(x) for x in parts[2:6]]
        m = max(vals)
        if np.isnan(best) or m > best:
            best, comp = m, tuple(vals)
    return (best,) + comp


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--variants", required=True,
                    help="parquet with variant_id, chrom, pos, ref, alt")
    ap.add_argument("--out", default="phase1/data/evidence/spliceai_walker.parquet")
    ap.add_argument("--fasta", default=str(REF_FASTA),
                    help="GRCh38 fasta; the atlas subset covers chr 2/3/13/16/17 only, "
                         "so an external gene on another chromosome needs its own")
    ap.add_argument("--distance", type=int, default=DISTANCE,
                    help="SpliceAI -D; 4999 is Walker's setting, 50 is the atlas column's")
    ap.add_argument("--column", default=None,
                    help="output column name; defaults to spliceai_walker for -D 4999")
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--report-every", type=int, default=250)
    ap.add_argument("--checkpoint-every", type=int, default=250,
                    help="write partial results this often; a re-run resumes from them")
    ap.add_argument("--restart", action="store_true",
                    help="ignore any checkpoint and score every variant again")
    args = ap.parse_args()

    from spliceai.utils import Annotator, get_delta_scores
    import spliceai

    variants = pd.read_parquet(args.variants)
    if "gene" not in variants.columns:
        variants["gene"] = ""
    variants["chrom_s"] = variants["chrom"].astype(str)
    variants = variants.sort_values(["chrom_s", "pos"]).reset_index(drop=True)

    distance = args.distance
    column = args.column or ("spliceai_walker" if distance == DISTANCE
                             else f"spliceai_d{distance}")
    ann = Annotator(args.fasta, ANNOTATION)
    fullprec = build_fullprec()

    if args.validate:
        sample = variants.sample(n=min(120, len(variants)), random_state=20260917)
        checked = mismatched = 0
        for r in sample.itertuples():
            rec = Rec(r.chrom_s, int(r.pos), r.ref, r.alt)
            stock = get_delta_scores(rec, ann, distance, MASK)
            fine = fullprec(rec, ann, distance, MASK)
            if len(stock) != len(fine):
                mismatched += 1
                continue
            for a, b in zip(stock, fine):
                pa, pb = a.split("|"), b.split("|")
                if pa[2] == "." or pb[2] == ".":
                    continue
                checked += 1
                if [f"{float(x):.2f}" for x in pb[2:6]] != pa[2:6]:
                    mismatched += 1
                    print("MISMATCH", pa[2:6], pb[2:6])
        print(f"validated {checked} score fields on {len(sample)} variants at "
              f"distance {distance}, mask {MASK}; {mismatched} mismatches")
        return 1 if mismatched else 0

    if not args.run:
        ap.error("pass --validate or --run")

    if args.limit:
        variants = variants.head(args.limit)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    ckpt_path = Path(str(out_path) + ".partial.parquet")
    cols = [column] + [f"{column}_{f}" for f in FIELDS]

    done: dict[str, tuple] = {}
    if ckpt_path.exists() and not args.restart:
        prev = pd.read_parquet(ckpt_path)
        if all(c in prev.columns for c in cols):
            for rec in prev[["variant_id"] + cols].itertuples(index=False):
                done[rec[0]] = tuple(rec[1:])
        print(f"[resume] {len(done):,} of {len(variants):,} already scored in "
              f"{ckpt_path.name}", flush=True)

    def flush(scored: dict) -> None:
        frame = variants[["variant_id", "gene", "chrom", "pos", "ref", "alt"]].copy()
        arr = np.asarray([scored.get(v, (np.nan,) * 5) for v in frame["variant_id"]],
                         dtype=float)
        for j, c in enumerate(cols):
            frame[c] = arr[:, j]
        frame[frame["variant_id"].isin(scored)].to_parquet(ckpt_path, index=False)

    todo = [r for r in variants.itertuples() if r.variant_id not in done]
    t0 = time.time()
    for i, r in enumerate(todo, 1):
        rec = Rec(r.chrom_s, int(r.pos), r.ref, r.alt)
        try:
            done[r.variant_id] = aggregate(fullprec(rec, ann, distance, MASK))
        except Exception:
            done[r.variant_id] = (np.nan,) * 5
        if i % args.checkpoint_every == 0:
            flush(done)
        if i % args.report_every == 0:
            el = time.time() - t0
            print(f"{len(done):>6}/{len(variants)}  {el/60:6.1f} min elapsed, "
                  f"{(len(todo)-i)*el/i/60:6.1f} min remaining", flush=True)
    flush(done)

    out = variants[["variant_id", "gene", "chrom", "pos", "ref", "alt"]].copy()
    arr = np.asarray([done.get(v, (np.nan,) * 5) for v in out["variant_id"]], dtype=float)
    for j, c in enumerate(cols):
        out[c] = arr[:, j]
    out.to_parquet(out_path, index=False)
    ckpt_path.unlink(missing_ok=True)

    sha = hashlib.sha256(out_path.read_bytes()).hexdigest()
    prov = {
        "column": column,
        "tool": "SpliceAI",
        "package_version": getattr(spliceai, "__version__", "1.3.1"),
        "distance": distance,
        "mask": MASK,
        "statistic": "max over ds_ag, ds_al, ds_dg, ds_dl; full float precision",
        "reference_fasta": args.fasta,
        "annotation": ANNOTATION,
        "n_variants": int(len(out)),
        "n_scored": int(out["spliceai_walker"].notna().sum()),
        "basis": "Walker et al. 2023 Methods; phase1/config/walker2023.yaml",
        "output_sha256": sha,
        "scored_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "runtime_minutes_this_run": round((time.time() - t0) / 60, 1),
        "resumed": bool(len(todo) < len(variants)),
    }
    Path(str(out_path) + ".provenance.json").write_text(json.dumps(prov, indent=2) + "\n")
    print(json.dumps(prov, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
