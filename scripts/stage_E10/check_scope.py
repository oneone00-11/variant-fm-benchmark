"""E10: scripts/stage_E18/check_scope.py, run as it stands, its run record written here.

The release commit of 3.1.0 (E10, commit 2) is checked against the commit before it the
way the commits of E5 to E9 were, its scope declared in
scripts/stage_E10/scope_e10_c2.json. It changes no output: the version, its DOI, the
addresses of the release's record, comments, the README and tests; so the scope declares
nothing and every output must stay byte for byte the same. This runs E18's script
unchanged; only the run record goes to this folder's run_log/, because this folder's
_runlog is imported first and E18's script then finds it already loaded.

Run from the repository root:
    .venv/bin/python scripts/stage_E10/check_scope.py --scope scripts/stage_E10/scope_e10_c2.json --base 020683b
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
