"""E9: scripts/stage_E18/check_scope.py, run as it stands, its run record written here.

Each E9 commit is checked against the one before it the way E5 to E8 were, its scope
declared in scripts/stage_E9/scope_e9_c<N>.json. Commits 1 to 3 (the fetching, the
workflows, the binding table of the submitted text) change no output, so their scopes
declare nothing and every output must stay byte for byte the same; commit 4 (stage 42)
declares the files it adds and the rows that register them. This runs E18's script
unchanged; only the run record goes to this folder's run_log/, because this folder's
_runlog is imported first and E18's script then finds it already loaded.

Run from the repository root:
    .venv/bin/python scripts/stage_E9/check_scope.py --scope scripts/stage_E9/scope_e9_c1.json --base 028dbe5
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
