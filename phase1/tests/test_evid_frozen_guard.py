"""E18.0 rule 1: the frozen thresholds and everything they rest on are as tagged.

phase1/data/evidence/FROZEN_AT_TAG.json records, at the tag frozen-thresholds-v1, the
sha256 of the analysis set, the seven genes' AlphaGenome scores, the thresholds, Table 3,
Supplementary Tables S6 and S7, every interval curve and every output of the fusion
(scripts/stage_E18/freeze_guard.py). HEAD must hold each of them byte for byte as
recorded; in a release archive, which has no .git, the tree as unpacked must.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

PHASE1 = Path(__file__).resolve().parents[1]
REPO = PHASE1.parent
RECORD = PHASE1 / "data/evidence/FROZEN_AT_TAG.json"


def _committed(paths: list[str]) -> dict[str, bytes]:
    """Each path as HEAD has it, through one `git cat-file --batch`; the tree as it is
    where there is no .git."""
    if not (REPO / ".git").exists():
        return {p: (REPO / p).read_bytes() for p in paths}
    out = subprocess.run(["git", "-C", str(REPO), "cat-file", "--batch"],
                         input="".join(f"HEAD:{p}\n" for p in paths).encode(),
                         check=True, capture_output=True).stdout
    blobs, pos = {}, 0
    for p in paths:
        nl = out.index(b"\n", pos)
        head = out[pos:nl].split()
        assert head[-1] != b"missing", f"HEAD has no {p}"
        size = int(head[2])
        blobs[p] = out[nl + 1: nl + 1 + size]
        pos = nl + 1 + size + 1
    return blobs


def test_the_record_is_the_tags():
    rec = json.loads(RECORD.read_text())
    assert rec["tag"] == "frozen-thresholds-v1"
    assert rec["commit"] == "1d030192ff7c91f8eed63c95154122437c819dba"
    assert rec["tagged_utc"] == "2026-10-22T16:34:25Z"
    groups = {g: set(v) for g, v in rec["groups"].items()}
    assert set().union(*groups.values()) == set(rec["files"])
    assert len(groups["the interval likelihood-ratio curves"]) == 70
    assert "phase1/data/evidence/alphagenome_official/all_scores.tsv" in rec["files"]


def test_head_holds_every_frozen_file_as_tagged():
    rec = json.loads(RECORD.read_text())["files"]
    blobs = _committed(sorted(rec))
    moved = [p for p, e in rec.items()
             if hashlib.sha256(blobs[p]).hexdigest() != e["sha256"]]
    assert not moved, f"frozen files changed since the tag: {moved}"
