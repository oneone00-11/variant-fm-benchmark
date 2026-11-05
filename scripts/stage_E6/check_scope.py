"""E6: scripts/stage_E18/check_scope.py, run as it stands, its run record written here.

E6 checks its commit's outputs against the previous commit the way E5 did. The commit
adds (c) to Figure 1 (the author's decision of 2026-11-04) and declares exactly the
four outputs that change (scope_figure1c.json): Figure 1's PDF and PNG; contents.csv,
of which only Figure 1's line; and manifest.csv, of which only the rows of those three
files. Every other output, the supplementary-table workbook among them, must stay byte
for byte the same. This runs E18's script unchanged; only the run record goes to this
folder's run_log/, because this folder's _runlog is imported first and E18's script then
finds it already loaded.

Run from the repository root:
    .venv/bin/python scripts/stage_E6/check_scope.py --scope scripts/stage_E6/scope_figure1c.json --base 6ce66eb
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
