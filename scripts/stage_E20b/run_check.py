"""E20b.4 (commit 4): `bash reproduce.sh --check` on this machine, its result in run_log.

The quick path reads every recorded checksum and runs the test suite; after E20b.4 a file
kept only in a gitignored cache is read at its manifest's cached_at (checked when the
cache holds it, waited for when it does not). This runs the command as a reader would and
notes the lines that carry the result: the checksum count, each file found missing or
different, each input waited for, the test summary and the verdict, with the output's
last lines.

Run from the repository root (no network, no stage run):
    .venv/bin/python scripts/stage_E20b/run_check.py
"""
from __future__ import annotations

import re
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, utc_now  # noqa: E402

RESULT = re.compile(r"checksums over|^\s+(missing|differs)\s|is not here, so its checksums wait|"
                    r"\d+ passed|\d+ failed|^CHECKED|^CHECK FAILED")


def main() -> None:
    with RunLog(__file__) as log:
        start = time.monotonic()
        run = subprocess.run(["bash", "reproduce.sh", "--check"], cwd=REPO, capture_output=True, text=True)
        out = (run.stdout + run.stderr).splitlines()
        lines = [x.rstrip() for x in out if RESULT.search(x)]
        log.note("command", "bash reproduce.sh --check")
        log.note("finished_utc", utc_now())
        log.note("minutes", round((time.monotonic() - start) / 60, 1))
        log.note("returncode", run.returncode)
        log.note("result_lines", lines)
        log.note("output_tail", out[-40:])
        print("\n".join(lines))
        print(f"returncode {run.returncode}")
        if run.returncode != 0:
            raise SystemExit("reproduce.sh --check failed: see the run record")


if __name__ == "__main__":
    main()
