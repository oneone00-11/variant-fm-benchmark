"""E7: scripts/stage_E18/check_scope.py, run as it stands, its run record written here.

E7 checks its commit's outputs against the previous commit the way E5 and E6 did,
declared to the cell or row that the check can name (scope_e7.json): Supplementary Table
S32's rows each without the two threshold-change columns; S6's nine rows of the combined
Atlas score, in their count, summary and per-fold fields; S16 (a) with one column
renamed and every cell as it was; contents.csv in the lines of S6, S7 and S32;
manifest.csv in the rows of the files that change; predictor_training_provenance.csv in
the fusion row's training data; the supplement's manifest in S3's hash of that file.
Figure 2, Supplementary Figure S1 and the workbook are declared as files and checked
pixel by pixel and cell by cell by check_e7.py. Every other output must stay byte for
byte the same. This runs E18's script unchanged; only the run record goes to this
folder's run_log/, because this folder's _runlog is imported first and E18's script then
finds it already loaded.

Run from the repository root:
    .venv/bin/python scripts/stage_E7/check_scope.py --scope scripts/stage_E7/scope_e7.json --base 9b8eae4
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
