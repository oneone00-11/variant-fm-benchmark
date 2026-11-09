"""E8: scripts/stage_E18/check_scope.py, run as it stands, its run record written here.

Each E8 commit is checked against the one before it the way E5, E6 and E7 were, declared
to the row, field or line that the check can name (scripts/stage_E8/scope_e8_c<N>.json);
a file the check can only name whole (the workbook, contents.csv and manifest.csv after
the renumbering, a figure) is declared as a file with its reason and checked cell by
cell, row by row or pixel by pixel by check_e8.py. Every other output must stay byte for
byte the same. This runs E18's script unchanged; only the run record goes to this
folder's run_log/, because this folder's _runlog is imported first and E18's script then
finds it already loaded.

Run from the repository root:
    .venv/bin/python scripts/stage_E8/check_scope.py --scope scripts/stage_E8/scope_e8_c1.json --base d6f7bc5
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
