"""E18.2: SpliceAI for the external genes, by the analysis set's own scorers.

Two columns, two scorers, both run in the atlas's pinned SpliceAI environment
(<atlas>/models/spliceai/.venv: SpliceAI 1.3.1, tensorflow 2.19.1 CPU):
  --mode d50     spliceai: the atlas's models/spliceai/score_fullprec.py (-D 50, mask 0,
                 full float precision), the scorer of the analysis set's column; it is
                 loaded from the atlas checkout, byte for byte the release's file;
  --mode walker  spliceai_walker: phase1/src/evid_score_spliceai_walker.py (-D 4999,
                 mask 0, full float precision, the four deltas beside the maximum), the
                 scorer of the analysis set's column.
Each scores the 39 anchors and the 6,667 external variants of e18_score_input.parquet in
one run. The reference is the one the analysis set's column was scored on for the
chromosomes it covers (the atlas's grch38_subset.fa: 2, 3, 13, 16, 17), and the Ensembl
release 112 chromosome the E16 reference check fetched for 1, 7 and X; SpliceAI's own
bundled annotation (grch38) throughout. Every value is cached as Python's repr in
phase1/data/e18_cache/spliceai_<mode>.tsv (gitignored), appended as it is scored, so a
stopped run resumes; scripts/stage_E18/build_columns.py checks the anchors again and
writes the committed files from the cache. --trial scores the anchors and the first 50
external variants into a cache of its own and projects the full run's time.

Run from the repository root:
    <atlas>/models/spliceai/.venv/bin/python scripts/stage_E18/run_spliceai.py --mode d50
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

ATLAS = Path(os.environ.get("EVID_ATLAS_REPO", str(REPO.parent / "functional-standard-atlas")))
SUBSET = ATLAS / "data/refs/grch38_subset.fa"
REFS = REPO / "phase1/data/evidence/refs"
# mode -> (scorer, its path in the atlas release or None, what it is, the column)
SCORERS = {
    "d50": (ATLAS / "models/spliceai/score_fullprec.py", "models/spliceai/score_fullprec.py",
            "the atlas's models/spliceai/score_fullprec.py", "spliceai"),
    "walker": (REPO / "phase1/src/evid_score_spliceai_walker.py", None,
               "phase1/src/evid_score_spliceai_walker.py", "spliceai_walker"),
}


def fasta_for(chrom: str) -> Path:
    return SUBSET if chrom in {"2", "3", "13", "16", "17"} else REFS / f"chr{chrom}.fa"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=sorted(SCORERS), required=True)
    ap.add_argument("--trial", action="store_true",
                    help=f"the anchors and the first {E.TRIAL_N} external variants only")
    args = ap.parse_args()
    from spliceai.utils import Annotator
    path, release_rel, what, column = SCORERS[args.mode]
    code_sha = E.release_checked(path, release_rel) if release_rel else E.sha(path)
    spec = importlib.util.spec_from_file_location(f"scorer_{args.mode}", path)
    S = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(S)
    with RunLog(__file__) as log:
        log.input(path)
        log.input(E.INPUT)
        df = E.load_input(args.trial)
        E.CACHE.mkdir(parents=True, exist_ok=True)
        cache = E.cache_path(f"spliceai_{args.mode}", args.trial, ".tsv")
        done = set()
        if cache.exists():
            done = {line.split("\t", 1)[0] for line in cache.read_text().splitlines()[1:] if line}
        else:
            cache.write_text("variant_id\t" + ("value" if args.mode == "d50" else
                             "value\tds_ag\tds_al\tds_dg\tds_dl") + "\tscored_utc\n")
        todo = df[~df["variant_id"].isin(done)]
        print(f"[{args.mode}{' trial' if args.trial else ''}] {len(df)} variants, "
              f"{len(done)} cached, {len(todo)} to score", flush=True)
        annotators, fullprec = {}, S.build_fullprec()
        refs = sorted({str(fasta_for(c)) for c in df["chrom"]})
        for r in refs:
            log.input(Path(r))
        start, clock = utc_now(), E.Clock()
        for i, r in enumerate(todo.itertuples(index=False), 1):
            ref = str(fasta_for(r.chrom))
            if ref not in annotators:
                annotators[ref] = Annotator(ref, S.ANNOTATION)
            rec = type("Rec", (), {})()
            rec.chrom, rec.pos, rec.ref, rec.alts = str(r.chrom), int(r.pos), r.ref, [r.alt]
            fields = fullprec(rec, annotators[ref], S.DISTANCE, S.MASK)
            vals = S.aggregate(fields)
            vals = vals if isinstance(vals, tuple) else (vals,)
            with open(cache, "a") as f:
                f.write("\t".join([r.variant_id] + [repr(float(v)) for v in vals]
                                  + [utc_now()]) + "\n")
            if i % 200 == 0:
                el = clock.seconds()
                print(f"  {i}/{len(todo)}  {el / 60:.1f} min, {el / i:.2f} s each", flush=True)
        seconds = clock.seconds()
        scored = {line.split("\t")[0]: float(line.split("\t")[1])
                  for line in cache.read_text().splitlines()[1:] if line}
        tp = (E.throughput(seconds, len(todo), len(E.load_input(False)))
              if args.trial and len(todo) else None)
        E.report(log, E.anchor_check(column, scored), tp)
        freeze = subprocess.run([sys.executable, "-m", "pip", "freeze"], check=True,
                                capture_output=True, text=True).stdout.split()
        env = {"scorer": what, "scorer_sha256": code_sha,
               "distance": S.DISTANCE, "mask": S.MASK, "annotation": S.ANNOTATION,
               "references": {r: rel(Path(r)) for r in refs},
               "interpreter": sys.executable, "python": platform.python_version(),
               "platform": platform.platform(), "pip_freeze": freeze,
               "run_started_utc": start, "run_ended_utc": utc_now(),
               "seconds": round(seconds, 1), "n_scored_this_run": len(todo)}
        envf = E.cache_path(f"spliceai_{args.mode}", args.trial, ".environment.json")
        E.write_json(envf, env)
        log.output(cache)
        log.output(envf)
        log.note("n_scored_this_run", len(todo))
        print(f"[{args.mode}] done: {rel(cache)}")


if __name__ == "__main__":
    main()
