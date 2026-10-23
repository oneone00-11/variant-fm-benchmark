"""E18.2: Pangolin for the external genes, by the analysis set's own scorer.

The atlas's models/pangolin/score_fullprec.py (Pangolin at git 5cf94b8, -d 50, the
default masking, max(gain, |loss|) over transcript records, full float precision), in
its pinned environment (models/pangolin/.venv of the atlas checkout that ran it: torch
2.13.0, CPU), on one thread. The thread count matters to the last bits: on one thread
the scorer reproduces all 39 anchors bit for bit, and so on 3 and 5, which give the same
values as one on external variants too; on 2, 4, 6 and 8 it does not (E18 report,
section 4). So the run requires OMP_NUM_THREADS = MKL_NUM_THREADS = 1 and sets
torch.set_num_threads(1). The scorer is loaded from that checkout, where its gene
database lies, and must be byte for byte the release's file.

Each worker (--worker k --workers n) scores the variants at positions k, k + n, ... of
e18_score_input.parquet, the anchors among them, so every worker's anchors are scored in
its own run; workers are separate processes on one thread each, so a variant's value
does not depend on n. Reference and gene database: the atlas's own (grch38_subset.fa,
grch38_subset.gtf.db) on chromosomes 2, 3, 13, 16 and 17; on 1, 7 and X, the Ensembl
release 112 chromosome (phase1/data/evidence/refs) and the database
scripts/stage_E18/build_pangolin_refs.py built by the atlas's recipe. A variant the
scorer cannot score is NaN, as in the atlas's run.

Cache: phase1/data/e18_cache/pangolin_w<k>.tsv (gitignored), appended as scored. --trial
scores the anchors and the first 50 external variants, as one worker, into a cache of
its own, and projects the full run's time on --workers processes.

Run from the repository root, one process per worker:
    OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 <atlas>/models/pangolin/.venv/bin/python \\
        scripts/stage_E18/run_pangolin.py --worker 0 --workers 6
"""
from __future__ import annotations

import argparse
import importlib.util
import os
import platform
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel, utc_now  # noqa: E402
import _e18 as E  # noqa: E402

# the checkout whose models/pangolin/.venv scored the atlas's column
PANGOLIN_ATLAS = Path(os.environ.get(
    "EVID_PANGOLIN_ATLAS", "/Users/cliffzhang/Documents/kimi/workspace/functional-standard-atlas"))
SCORER = PANGOLIN_ATLAS / "models/pangolin/score_fullprec.py"
REFS = REPO / "phase1/data/evidence/refs"
DBS = E.CACHE / "refs"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--worker", type=int, default=0)
    ap.add_argument("--workers", type=int, required=True)
    ap.add_argument("--trial", action="store_true",
                    help=f"the anchors and the first {E.TRIAL_N} external variants, one worker")
    args = ap.parse_args()
    if (os.environ.get("OMP_NUM_THREADS"), os.environ.get("MKL_NUM_THREADS")) != ("1", "1"):
        raise SystemExit("set OMP_NUM_THREADS=1 MKL_NUM_THREADS=1: the anchors hold on one thread")
    import numpy as np
    import torch
    torch.set_num_threads(1)
    sha = E.release_checked(SCORER, "models/pangolin/score_fullprec.py")
    spec = importlib.util.spec_from_file_location("pangolin_fullprec", SCORER)
    pf = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pf)
    import gffutils
    import pangolin.pangolin as PP
    with RunLog(__file__) as log:
        log.input(SCORER)
        log.input(E.INPUT)
        df = E.load_input(args.trial)
        if not args.trial:
            df = df.iloc[args.worker::args.workers]
        E.CACHE.mkdir(parents=True, exist_ok=True)
        name = "pangolin" if args.trial else f"pangolin_w{args.worker}"
        cache = E.cache_path(name, args.trial, ".tsv")
        done = set()
        if cache.exists():
            done = {line.split("\t", 1)[0] for line in cache.read_text().splitlines()[1:] if line}
        else:
            cache.write_text("variant_id\tvalue\tscored_utc\n")
        todo = df[~df["variant_id"].isin(done)]
        print(f"[{name}{' trial' if args.trial else ''}] {len(df)} variants, {len(done)} cached, "
              f"{len(todo)} to score; torch threads {torch.get_num_threads()}", flush=True)
        gtf, models = pf.setup(PP)
        base = pf.Args(pf.DISTANCE)
        pick = {}
        for c in sorted(set(df["chrom"])):
            if c in {"2", "3", "13", "16", "17"}:
                pick[c] = (gtf, base, pf.REF_FASTA, pf.GTF_DB)
            else:
                a = pf.Args(pf.DISTANCE)
                a.reference_file = str(REFS / f"chr{c}.fa")
                db = DBS / f"chr{c}.gtf.db"
                pick[c] = (gffutils.FeatureDB(str(db)), a, Path(a.reference_file), db)
        for _, _, fa, db in pick.values():
            log.input(fa)
            log.input(db)
        fullprec = pf.build_fullprec(PP)
        start, clock, failed = utc_now(), E.Clock(), 0
        for i, r in enumerate(todo.itertuples(index=False), 1):
            g, a, _, _ = pick[r.chrom]
            try:
                sc = fullprec(0, str(r.chrom), int(r.pos), r.ref, r.alt, g, models, a)
                val = pf.aggregate(sc) if sc != -1 else np.nan
            except Exception:                                   # noqa: BLE001
                val, failed = np.nan, failed + 1
            with open(cache, "a") as f:
                f.write(f"{r.variant_id}\t{float(val)!r}\t{utc_now()}\n")
            if i % 100 == 0:
                el = clock.seconds()
                print(f"  {name}: {i}/{len(todo)}  {el / 60:.1f} min, {el / i:.2f} s each",
                      flush=True)
        seconds = clock.seconds()
        scored = {line.split("\t")[0]: float(line.split("\t")[1])
                  for line in cache.read_text().splitlines()[1:] if line}
        tp = (E.throughput(seconds, len(todo), len(E.load_input(False)), args.workers)
              if args.trial and len(todo) else None)
        E.report(log, E.anchor_check("pangolin", scored), tp)
        freeze = subprocess.run([sys.executable, "-m", "pip", "freeze"], check=True,
                                capture_output=True, text=True).stdout.split()
        env = {"scorer": "the atlas's models/pangolin/score_fullprec.py", "scorer_sha256": sha,
               "distance": pf.DISTANCE, "mask": base.mask, "torch_threads": torch.get_num_threads(),
               "OMP_NUM_THREADS": os.environ["OMP_NUM_THREADS"],
               "MKL_NUM_THREADS": os.environ["MKL_NUM_THREADS"],
               "worker": args.worker, "workers": args.workers,
               "references": {c: [rel(fa), rel(db)] for c, (_, _, fa, db) in pick.items()},
               "interpreter": sys.executable, "python": platform.python_version(),
               "platform": platform.platform(), "pip_freeze": freeze,
               "run_started_utc": start, "run_ended_utc": utc_now(),
               "seconds": round(seconds, 1), "n_scored_this_run": len(todo),
               "n_exceptions": failed}
        envf = E.cache_path(name, args.trial, ".environment.json")
        E.write_json(envf, env)
        log.output(cache)
        log.output(envf)
        log.note("n_scored_this_run", len(todo))
        log.note("n_exceptions", failed)
        print(f"[{name}] done: {rel(cache)}")


if __name__ == "__main__":
    main()
