"""E20.1 (C12): the true positives' event records, by E6's own scorers, unchanged.

Runs phase1/src/evid_score_spliceai_events.py or phase1/src/evid_score_pangolin_events.py
on one score input, with the arguments E6 and E2.7 used (E20.0 rule 3), in the model
environments, from the repository root:

  spliceai  the seven genes' true positives and anchor 2's SpliceAI variants
  pangolin  the seven genes' true positives and anchor 2's Pangolin variants
  ddx3x     DDX3X's true positives (SpliceAI; chromosome X, so E2.7's --fasta)
  tp53      TP53's true positives (SpliceAI)

SpliceAI: the atlas's models/spliceai/.venv (SpliceAI 1.3.1), --distance 4999, the
scorer's own mask (0, no masking) and reference (the atlas's grch38_subset.fa; DDX3X
phase1/data/evidence/refs/chrX.fa), as E6's and E2.7's provenance records give them.
Pangolin: models/pangolin/.venv of the atlas checkout whose environment scored the
atlas's column (scripts/stage_E18/run_pangolin.py; the only Pangolin environment on this
machine; E6's provenance does not name its interpreter): Pangolin at git 5cf94b8, torch
2.13.0, the scorer's own settings (-d 50, the CLI's default masking, five chunks, two
threads per chunk), its own work folder, phase1/data/evidence/pangolin_events_work_tp.
No thread variable is set for the scorer beyond what the scorer sets itself.

The scorer writes <out>.provenance.json; this writes <out>.environment.json beside it
(the command, the interpreter and its pip freeze --all, the scorer's sha256, the reference
files' sizes and sha256, the run's start, end and duration) and a run record in run_log/.

    .venv/bin/python scripts/stage_E20/run_c12_events.py spliceai
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel, sha256_of, utc_now  # noqa: E402

HOME = Path.home()
SPLICEAI_ATLAS = HOME / "work/functional-standard-atlas"
PANGOLIN_ATLAS = HOME / "Documents/kimi/workspace/functional-standard-atlas"
SPLICEAI_PY = SPLICEAI_ATLAS / "models/spliceai/.venv/bin/python"
PANGOLIN_PY = PANGOLIN_ATLAS / "models/pangolin/.venv/bin/python"
SPLICEAI_SCORER = "phase1/src/evid_score_spliceai_events.py"
PANGOLIN_SCORER = "phase1/src/evid_score_pangolin_events.py"
# the scorers' default references (EVID_ATLAS_REPO unset: the checkout beside this one)
ATLAS_REFS = SPLICEAI_ATLAS / "data/refs"
EV = "phase1/data/evidence"
JOBS = {
    "spliceai": {"python": SPLICEAI_PY, "scorer": SPLICEAI_SCORER,
                 "input": f"{EV}/inframe_score_input_tp.parquet",
                 "out": f"{EV}/inframe_events_tp.parquet",
                 "extra": ["--distance", "4999"],
                 "references": [ATLAS_REFS / "grch38_subset.fa"]},
    "pangolin": {"python": PANGOLIN_PY, "scorer": PANGOLIN_SCORER,
                 "input": f"{EV}/inframe_score_input_pangolin_tp.parquet",
                 "out": f"{EV}/inframe_events_pangolin_tp.parquet",
                 "extra": ["--work", f"{EV}/pangolin_events_work_tp"],
                 "references": [ATLAS_REFS / "grch38_subset.fa",
                                ATLAS_REFS / "grch38_subset.gtf.db"]},
    "ddx3x": {"python": SPLICEAI_PY, "scorer": SPLICEAI_SCORER,
              "input": f"{EV}/external/ddx3x_inframe_subset_tp.parquet",
              "out": f"{EV}/external/ddx3x_inframe_events_tp.parquet",
              "extra": ["--distance", "4999", "--fasta", f"{EV}/refs/chrX.fa"],
              "references": [REPO / f"{EV}/refs/chrX.fa"]},
    "tp53": {"python": SPLICEAI_PY, "scorer": SPLICEAI_SCORER,
             "input": f"{EV}/external/tp53_inframe_subset_12nt_tp.parquet",
             "out": f"{EV}/external/tp53_inframe_events_12nt_tp.parquet",
             "extra": ["--distance", "4999"],
             "references": [ATLAS_REFS / "grch38_subset.fa"]},
}
THREAD_VARS = ["OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
               "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "TF_NUM_INTRAOP_THREADS",
               "TF_NUM_INTEROP_THREADS"]


def as_given(path: Path | str) -> str:
    """The path as passed, symlinks not followed (rel() resolves them: a venv's python
    would read as its base interpreter)."""
    p = Path(path)
    if not p.is_absolute():
        return str(p)
    for base, prefix in ((REPO, ""), (HOME, "~/")):
        if p.is_relative_to(base):
            return prefix + str(p.relative_to(base))
    return str(p)


def pip_freeze(python: Path) -> list[str]:
    out = subprocess.run([str(python), "-m", "pip", "freeze", "--all"], check=True,
                         capture_output=True, text=True).stdout
    return [line for line in out.splitlines() if line.strip()]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("job", choices=sorted(JOBS))
    job = JOBS[ap.parse_args().job]
    python = job["python"]
    if not python.exists():
        raise SystemExit(f"{python} not found: stopped")
    out = REPO / job["out"]
    if out.exists():
        raise SystemExit(f"{job['out']} exists already: stopped")
    cmd = [str(python), job["scorer"], "--variants", job["input"]] + job["extra"] + \
          ["--out", job["out"], "--run"]
    env = {k: v for k, v in os.environ.items() if k not in THREAD_VARS + ["EVID_ATLAS_REPO"]}
    with RunLog(__file__) as log:
        log.input(REPO / job["scorer"])
        log.input(REPO / job["input"])
        refs = [{"path": as_given(p), "resolved": rel(Path(p).resolve()),
                 "bytes": Path(p).stat().st_size, "sha256": sha256_of(Path(p))}
                for p in job["references"]]
        freeze = pip_freeze(python)
        python_version = subprocess.run([str(python), "-c", "import sys; print(sys.version.split()[0])"],
                                        check=True, capture_output=True, text=True).stdout.strip()
        log.note("command", [as_given(cmd[0])] + cmd[1:])
        started, t0 = utc_now(), time.time()
        proc = subprocess.run(cmd, cwd=REPO, env=env, capture_output=True, text=True)
        seconds = round(time.time() - t0, 1)
        ended = utc_now()
        log.note("returncode", proc.returncode)
        log.note("stdout_tail", proc.stdout.splitlines()[-40:])
        log.note("stderr_tail", proc.stderr.splitlines()[-20:])
        if proc.returncode != 0:
            raise SystemExit(f"the scorer exited {proc.returncode}: stopped")
        environment = {
            "file": job["out"],
            "command": [as_given(cmd[0])] + cmd[1:],
            "working_directory": "the repository root",
            "interpreter": as_given(python),
            "python": python_version,
            "platform": platform.platform(),
            "pip_freeze": freeze,
            "scorer": job["scorer"],
            "scorer_sha256": sha256_of(REPO / job["scorer"]),
            "references": refs,
            "thread_variables_set_by_this_run": "none (the scorer's own settings only)",
            "EVID_ATLAS_REPO": "unset (the scorer's default: the checkout beside this one)",
            "run_started_utc": started,
            "run_ended_utc": ended,
            "seconds": seconds,
        }
        if job["scorer"] == PANGOLIN_SCORER:
            environment["pangolin_cli"] = as_given(python.parent / "pangolin")
            environment["pangolin_settings"] = (
                "as evid_score_pangolin_events.py runs the CLI: -d 50, the CLI's default mask "
                "(True), five chunks, OMP_NUM_THREADS=2 and MKL_NUM_THREADS=2 per chunk process")
        else:
            environment["spliceai_settings"] = (
                "as evid_score_spliceai_events.py runs it: --distance 4999, mask 0, annotation "
                "grch38, the four delta scores at full precision")
        env_path = Path(str(out) + ".environment.json")
        env_path.write_text(json.dumps(environment, indent=1, ensure_ascii=False) + "\n")
        for p in (out, Path(str(out) + ".provenance.json"), env_path):
            log.output(p)
        if job["scorer"] == PANGOLIN_SCORER:
            work = REPO / f"{EV}/pangolin_events_work_tp"
            for p in sorted(work.rglob("*")):
                if p.is_file():
                    log.output(p)
        print(f"{job['out']}: {seconds / 60:.1f} min")


if __name__ == "__main__":
    main()
