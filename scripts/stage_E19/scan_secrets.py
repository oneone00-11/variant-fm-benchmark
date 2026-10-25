"""E19 (E17.0 rule 3): tests/test_no_secrets.py's patterns, unchanged, over the files
this stage writes, by scripts/stage_E18/scan_secrets.py as it stands: E18's folders,
and this stage's scripts and run logs (scripts/stage_E19), its outputs
(phase1/reports/evidence/suppl_analyses) and its fetched records and extracts
(phase1/data/evidence/e19). The ClinVar archive in the gitignored cache
(phase1/data/e19_cache) is downloaded data, kept as served and read by no code that
reads the key, and is not scanned, as E18's scanner skips GPN-MSA's local copy.

    .venv/bin/python scripts/stage_E19/scan_secrets.py [extra folder ...]
"""
from __future__ import annotations

import runpy
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
EXTRA = ["scripts/stage_E19", "phase1/reports/evidence/suppl_analyses", "phase1/data/evidence/e19"]

if __name__ == "__main__":
    sys.argv = [sys.argv[0]] + [str(REPO / e) for e in EXTRA] + sys.argv[1:]
    runpy.run_path(str(HERE.parent / "stage_E18" / "scan_secrets.py"), run_name="__main__")
