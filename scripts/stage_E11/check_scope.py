"""E11: scripts/stage_E18/check_scope.py, run as it stands, its run record written here.

The two commits made after release 3.1.0 was published (E11) are each checked against the
commit before it the way the commits of E5 to E10 were, each scope declared in
scripts/stage_E11/scope_e11_c<N>.json. Neither changes an output: commit 1 adds the run
records of the binding table (v7) built and checked on the submitted text with the release,
its commit and its DOI filled in, the table itself unchanged; commit 2 corrects a wording.
So each scope declares nothing and every output must stay byte for byte the same. This runs
E18's script unchanged; only the run record goes to this folder's run_log/, because this
folder's _runlog is imported first and E18's script then finds it already loaded.

Run from the repository root:
    .venv/bin/python scripts/stage_E11/check_scope.py --scope scripts/stage_E11/scope_e11_c1.json --base 29bbbe5
"""
from __future__ import annotations

import runpy
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import _runlog  # noqa: E402,F401  this folder's run record, found first by E18's script

if __name__ == "__main__":
    runpy.run_path(str(HERE.parent / "stage_E18" / "check_scope.py"), run_name="__main__")
