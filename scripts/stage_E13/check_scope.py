"""E13: scripts/stage_E18/check_scope.py, run as it stands, its run record written here.

The commit of E13 is checked against the commit before it the way the commits of E5 to E12
were, its scope declared in scripts/stage_E13/scope_e13_c1.json. It changes no output: the
binding table of the submitted text (v7, docs/), its build script, its tests and its run
records follow the abstract written anew on 2026-11-17. So the scope declares nothing and
every output must stay byte for byte the same. This runs E18's script unchanged; only the
run record goes to this folder's run_log/, because this folder's _runlog is imported first
and E18's script then finds it already loaded.

Run from the repository root:
    .venv/bin/python scripts/stage_E13/check_scope.py --scope scripts/stage_E13/scope_e13_c1.json --base 10128f9
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
