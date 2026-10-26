"""E20b.1 (commit 1): E6's two subsets under the current labels, by `python -m
src.evid_inframe --subset` itself, the code unchanged.

The committed subsets are read first; after the module has rewritten them from HEAD's
analysis set, each new subset must be the committed one, row for row and value for value,
with exactly the variants the pre-check found added (precheck_label_inputs.py's run
record): 6 for SpliceAI, 7 for Pangolin, none removed. Otherwise this stops, and the
committed files are restored from what was read.

Run in the reference image (修改单 §0 rule 9), offline:
    docker run ... IMAGE bash scripts/stage_E20b/in_image.sh scripts/stage_E20b/write_e6_subsets.py
"""
from __future__ import annotations

import glob
import json
import os
import runpy
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel  # noqa: E402

PHASE1 = REPO / "phase1"
sys.path.insert(0, str(PHASE1))
os.chdir(PHASE1)

EV = PHASE1 / "data/evidence"
FILES = {"spliceai": EV / "inframe_subset.parquet", "pangolin": EV / "inframe_subset_pangolin.parquet"}
PRECHECK = sorted(glob.glob(str(Path(__file__).resolve().parent / "run_log/precheck_label_inputs_*.json")))


def main() -> None:
    if not PRECHECK:
        raise SystemExit("no pre-check run record: run precheck_label_inputs.py first")
    pre = json.loads(Path(PRECHECK[-1]).read_text())
    expected = {r["file"]: r["only_derived"] for r in pre["notes"]["compared"]}
    with RunLog(__file__) as log:
        log.input(PRECHECK[-1])
        old = {}
        for tool, p in FILES.items():
            log.input(p)
            old[tool] = (pd.read_parquet(p), p.read_bytes())
        sys.argv = ["src.evid_inframe", "--subset"]
        runpy.run_module("src.evid_inframe", run_name="__main__")
        checks = {}
        try:
            for tool, p in FILES.items():
                o, n = old[tool][0], pd.read_parquet(p)
                added = n[~n["variant_id"].isin(o["variant_id"])]
                kept = n[n["variant_id"].isin(o["variant_id"])].reset_index(drop=True)
                same_old = (len(kept) == len(o) and list(kept.columns) == list(o.columns)
                            and all(((kept[c].isna() & o[c].isna()) | (kept[c] == o[c])).all()
                                    for c in o.columns))
                want = sorted(expected[p.name])
                checks[tool] = {"committed": len(o), "new": len(n), "added": sorted(added["variant_id"]),
                                "removed": sorted(set(o["variant_id"]) - set(n["variant_id"])),
                                "committed_rows_unchanged": bool(same_old),
                                "added_as_the_precheck_found": sorted(added["variant_id"]) == want,
                                "by_stratum": n.groupby("stratum").size().to_dict()}
                if not (same_old and checks[tool]["added_as_the_precheck_found"]
                        and not checks[tool]["removed"]):
                    raise SystemExit(f"{p.name}: not the committed subset plus the pre-check's "
                                     "variants: stopped")
        except SystemExit:
            for tool, p in FILES.items():
                p.write_bytes(old[tool][1])
            log.note("checks", checks)
            log.note("restored", [rel(p) for p in FILES.values()])
            raise
        for p in FILES.values():
            log.output(p)
        log.note("checks", checks)
        print(json.dumps(checks, indent=1))


if __name__ == "__main__":
    main()
