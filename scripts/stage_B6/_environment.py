"""container/ENVIRONMENT.json: the reference image's environment, one section per
script that measured it. B6.0 rule 1: every fact in it comes from a command's output,
written by the script that ran the command; nothing in it is typed by hand. Each
section names the script that wrote it and when; that script's run log
(run_log/<script>_<UTC timestamp>.json) holds the rest.
"""
from __future__ import annotations

import json
from pathlib import Path

from _runlog import REPO, rel, utc_now

PATH = REPO / "container" / "ENVIRONMENT.json"


def read() -> dict:
    return json.loads(PATH.read_text()) if PATH.exists() else {}


def update(section: str, data: dict, script: str, log) -> None:
    """Replace one section, keeping the others."""
    doc = read()
    doc[section] = {"written_by": rel(Path(script).resolve()), "written_utc": utc_now(),
                    **data}
    PATH.parent.mkdir(parents=True, exist_ok=True)
    PATH.write_text(json.dumps(doc, indent=1, ensure_ascii=False) + "\n")
    log.output(PATH)
