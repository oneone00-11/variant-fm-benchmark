"""E20b.1 (commit 1): E6's event records under the current labels, by E6's own scorers,
unchanged, on the whole of each new subset, in E20.1's environments with E20.1's settings.

  spliceai  inframe_subset.parquet (615): the atlas's models/spliceai/.venv (SpliceAI
            1.3.1), --distance 4999, the scorer's mask (0) and reference (the atlas's
            grch38_subset.fa), as E6's provenance records them
  pangolin  inframe_subset_pangolin.parquet (531): models/pangolin/.venv of the atlas
            checkout E20.1 used (Pangolin at git 5cf94b8, torch 2.13.0), the scorer's own
            settings (-d 50, the CLI's default masking, five chunks, two threads each)

Each run writes first to phase1/data/e20b_staging/ (the events, the scorer's provenance,
Pangolin's work folder); check_e6_rescore.py holds them against the committed records and
the score columns, and only then puts them in place. Beside each events file this writes
<out>.environment.json, which E6 did not record: the command, the interpreter and its pip
freeze --all, the scorer's sha256, the references' sha256, the run's start, end and length.

    .venv/bin/python scripts/stage_E20b/run_e6_events.py spliceai
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
STAGING = "phase1/data/e20b_staging"
JOBS = {
    "spliceai": {"python": SPLICEAI_PY, "scorer": SPLICEAI_SCORER,
                 "input": f"{EV}/inframe_subset.parquet",
                 "out": f"{STAGING}/inframe_events.parquet",
                 "final": f"{EV}/inframe_events.parquet",
                 "extra": ["--distance", "4999"],
                 "references": [ATLAS_REFS / "grch38_subset.fa"]},
    "pangolin": {"python": PANGOLIN_PY, "scorer": PANGOLIN_SCORER,
                 "input": f"{EV}/inframe_subset_pangolin.parquet",
                 "out": f"{STAGING}/inframe_events_pangolin.parquet",
                 "final": f"{EV}/inframe_events_pangolin.parquet",
                 "extra": ["--work", f"{STAGING}/pangolin_events_work"],
                 "references": [ATLAS_REFS / "grch38_subset.fa",
                                ATLAS_REFS / "grch38_subset.gtf.db"]},
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
            "file": job["final"],
            "written_to": job["out"] + " (moved by check_e6_rescore.py once both anchors held)",
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
            work = REPO / f"{STAGING}/pangolin_events_work"
            for p in sorted(work.rglob("*")):
                if p.is_file():
                    log.output(p)
        print(f"{job['out']}: {seconds / 60:.1f} min")


if __name__ == "__main__":
    main()
