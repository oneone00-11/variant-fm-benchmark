"""E5: scripts/stage_E18/check_scope.py, run as it stands, its run record written here.

E5 checks every commit's outputs against the previous commit the way E18 did. Its three
commits change no output: commits 1 and 2 change code, tests and documents only, so their
scope declares no file and every output must stay byte for byte the same; commit 3 adds
the layout stage's files and declares exactly those (E5's rules, E19.0 rule 1). This runs
E18's script unchanged; only the run record goes to this folder's run_log/, because this
folder's _runlog is imported first and E18's script then finds it already loaded.

Run from the repository root:
    .venv/bin/python scripts/stage_E5/check_scope.py --scope scripts/stage_E5/scope_gate.json --base <commit>
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
