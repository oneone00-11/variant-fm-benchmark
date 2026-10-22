"""E17, gate 1 review R6b: the notebook's own words around the code it scores with.

NOTEBOOK_SOURCE.json quotes the cells that compute the score. Three markdown cells of
the same notebook say what the score is for, why its window is 1 Mb, and how far it
reaches: cell 1 (the approach the AlphaGenome paper used to score ClinVar variants for
missplicing, and the recommended method), cell 10 (the 1 Mb window) and cell 16 (the
score's range, and that a recommended cutoff is still being worked on). This fetches the
notebook from the URL NOTEBOOK_SOURCE.json records, at the commit it pins, checks it
against the sha256 recorded there, and quotes those three cells verbatim, each with the
sha256 of its text, in phase1/data/evidence/alphagenome_official/NOTEBOOK_CONTEXT.json.
Cells are numbered from 0, as in NOTEBOOK_SOURCE.json, which is read and not rewritten.

Run from the repository root:
    .venv-ag-official/bin/python scripts/stage_E17/notebook_context.py
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel, sha256_of, utc_now  # noqa: E402
from notebook_source import _get                         # noqa: E402

DIR = REPO / "phase1" / "data" / "evidence" / "alphagenome_official"
SOURCE = DIR / "NOTEBOOK_SOURCE.json"
OUT = DIR / "NOTEBOOK_CONTEXT.json"
QUOTED = {1: "what the score is for", 10: "why the window is 1 Mb",
          16: "the score's range, and the cutoff"}


def main() -> None:
    with RunLog(__file__) as log:
        src = json.loads(SOURCE.read_text())
        log.input(SOURCE)
        blob = _get(src["raw_url"])
        got = hashlib.sha256(blob).hexdigest()
        if got != src["sha256"]:
            raise SystemExit(f"{src['path']} at {src['commit'][:12]} hashes to {got}, "
                             f"not {src['sha256']}")
        cells = json.loads(blob)["cells"]
        quoted = []
        for i, shows in QUOTED.items():
            if cells[i]["cell_type"] != "markdown":
                raise SystemExit(f"cell {i} is {cells[i]['cell_type']}, not markdown")
            text = "".join(cells[i].get("source", []))
            quoted.append({"cell": i, "type": "markdown", "shows": shows, "source": text,
                           "sha256": hashlib.sha256(text.encode()).hexdigest()})
        record = {
            "written_by": rel(Path(__file__).resolve()),
            "repository": src["repository"], "commit": src["commit"], "path": src["path"],
            "raw_url": src["raw_url"], "sha256": got, "bytes": len(blob),
            "code_record": {"file": rel(SOURCE), "sha256": sha256_of(SOURCE)},
            "accessed_utc": utc_now(),
            "cell_numbering": "from 0, as in NOTEBOOK_SOURCE.json",
            "quoted_cells": quoted,
        }
        OUT.write_text(json.dumps(record, indent=1, ensure_ascii=False) + "\n")
        log.output(OUT)
        log.note("sha256", got)
        print(f"{src['path']} @ {src['commit'][:12]}: sha256 ok; cells "
              f"{', '.join(str(i) for i in QUOTED)} -> {rel(OUT)}")


if __name__ == "__main__":
    main()
