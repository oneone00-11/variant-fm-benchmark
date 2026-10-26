"""E20b.3 (commit 3): C13 after the corrections, against ebd6b50 and against C6.

  * the 24 band rows: every count, ratio, interval and tier as at ebd6b50, cell for cell
    (and the rows' identity: gene, rule, band, column, status); the text cells that
    changed are listed (the reason of a band without labelled variants, the score
    column's basis for ENIGMA's rules after E20b.2);
  * c13 at HEAD is ebd6b50's byte for byte (commits 1 and 2 left it alone);
  * the BP7 rows: ENIGMA's equal c6_bp4_bp7.csv's condition (a), cell for cell; every
    BP7 row is listed;
  * the header: the lines that changed.
Any difference in the band rows' numbers stops.

Run in the reference image (修改单 §0 rule 9), offline, after the stage has run:
    docker run ... IMAGE bash scripts/stage_E20b/in_image.sh scripts/stage_E20b/check_c13.py
"""
from __future__ import annotations

import io
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog  # noqa: E402

sys.path.insert(0, str(REPO / "phase1"))
from src import evid_common as K  # noqa: E402

PATH = "phase1/reports/evidence/suppl_analyses/c13_vcep_band_lr.csv"
C6 = REPO / "phase1/reports/evidence/suppl_analyses/c6_bp4_bp7.csv"
BASE = "ebd6b50"
KEY = ["gene", "criterion", "stratum"]
SAME = ["part", "vcep", "cspec_id", "version", "rule", "strength_allowed", "score_column", "status",
        "n_damaging", "n_normal", "n_in_band", "n_damaging_in_band", "n_normal_in_band",
        "lr", "lr_lo", "lr_hi", "lr_for_tier", "tier_at_point", "bound_for_tier", "tier_at_bound"]
TEXT = ["score_column_basis", "not_evaluated", "reason"]
C6_FIELDS = ["n_at_positions", "n_damaging_at_positions", "n_normal_at_positions", "n_meeting",
             "n_damaging", "n_normal", "share_damaging"]


def at(rev: str) -> bytes:
    return subprocess.run(["git", "-C", str(REPO), "show", f"{rev}:{PATH}"], capture_output=True,
                          check=True).stdout


def equal(a, b) -> bool:
    if isinstance(a, float) and isinstance(b, float):
        return a == b or (np.isnan(a) and np.isnan(b))
    return (pd.isna(a) and pd.isna(b)) or a == b


def main() -> None:
    failures, notes = [], {}
    with RunLog(__file__) as log:
        log.input(REPO / PATH)
        log.input(C6)
        base = at(BASE)
        notes["head_is_base"] = at("HEAD") == base
        if not notes["head_is_base"]:
            failures.append(f"c13 at HEAD is not {BASE}'s")
        old = K.read_back(io.BytesIO(base), comment="#")
        new = K.read_back(REPO / PATH, comment="#")
        o = old[old["part"] == "band_lr"].set_index(KEY)
        n = new[new["part"] == "band_lr"].set_index(KEY)
        notes["band_rows"] = {"base": len(o), "new": len(n)}
        if list(o.index) != list(n.index) or len(o) != 24:
            failures.append("the band rows are not ebd6b50's 24, in its order")
        cells, text = 0, []
        for i in o.index:
            for c in SAME:
                cells += 1
                if not equal(o.loc[i, c], n.loc[i, c]):
                    failures.append(f"band row {' '.join(map(str, i))} {c}: {o.loc[i, c]!r} -> {n.loc[i, c]!r}")
            for c in TEXT:
                if not equal(o.loc[i, c], n.loc[i, c]):
                    text.append({"row": " ".join(map(str, i)), "column": c, "base": o.loc[i, c],
                                 "new": n.loc[i, c]})
        notes["band_cells_compared"] = cells
        notes["band_text_changed"] = text
        # BP7
        bp7 = new[new["part"] == "bp7"]
        c6 = K.read_back(C6, comment="#")
        c6 = c6[c6["condition"].str.startswith("(a) ")].set_index(["gene", "side"])
        enigma = bp7[bp7["cspec_id"].isin(["GN092", "GN097"])]
        mism = [f"{r.gene} {r.side} {f}" for r in enigma.itertuples() for f in C6_FIELDS
                if not equal(float(getattr(r, f)), float(c6.loc[(r.gene, r.side), f]))]
        notes["enigma_bp7_against_c6_a"] = {"rows": len(enigma), "cells": len(enigma) * len(C6_FIELDS),
                                            "differing": mism}
        if mism or len(enigma) != 6:
            failures.append(f"ENIGMA's BP7 rows are not c6 (a)'s: {mism}")
        cols = ["gene", "side", "stratum", "bp7_condition", "n_at_positions", "n_damaging_at_positions",
                "n_normal_at_positions", "n_phylop_missing", "n_meeting", "n_damaging", "n_normal",
                "share_damaging", "reason"]
        notes["bp7_rows"] = json.loads(bp7[[c for c in cols if c in bp7.columns]].to_json(orient="records"))
        notes["parts"] = new["part"].value_counts().to_dict()
        notes["columns_removed"] = [c for c in old.columns if c not in new.columns]
        notes["columns_added"] = [c for c in new.columns if c not in old.columns]
        hl = [x for x in base.decode().splitlines() if x.startswith("#")]
        nl = [x for x in (REPO / PATH).read_text().splitlines() if x.startswith("#")]
        notes["header"] = {"removed": [x for x in hl if x not in nl], "added": [x for x in nl if x not in hl]}
        log.note("checks", notes)
        log.note("failures", failures)
        print(json.dumps({k: v for k, v in notes.items() if k != "bp7_rows"}, indent=1, ensure_ascii=False)[:5000])
        if failures:
            for f in failures:
                print("FAIL", f)
            raise SystemExit("C13 differs from ebd6b50 beyond the corrections: stopped")
        print(f"band rows: {cells} cells as at {BASE}; ENIGMA's BP7 rows are c6 (a)'s")


if __name__ == "__main__":
    main()
