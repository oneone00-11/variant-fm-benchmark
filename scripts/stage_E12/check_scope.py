"""E12: scripts/stage_E18/check_scope.py, run as it stands, its run record written here.

The two commits of E12 are each checked against the commit before it the way the commits
of E5 to E11 were, each scope declared in scripts/stage_E12/scope_e12_c<N>.json. Neither
changes an output: commit 1 lets the front page's author line of the supplementary PDF be
given on the command line (the PDF is not in the repository); commit 2 adds the run records
of the binding table (v7) built and checked again on the submitted text, the table itself
unchanged. So each scope declares nothing and every output must stay byte for byte the
same. This runs E18's script unchanged; only the run record goes to this folder's run_log/,
because this folder's _runlog is imported first and E18's script then finds it already
loaded.

Run from the repository root:
    .venv/bin/python scripts/stage_E12/check_scope.py --scope scripts/stage_E12/scope_e12_c1.json --base aba04f9
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
