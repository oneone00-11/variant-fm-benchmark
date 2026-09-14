"""The prose pass on Methods 2.1 and 2.2 and on Limitations 4.1 is in the manuscript.

That pass split sentences and paragraphs in three sections, changed no number, and added
one clause to 4.1 ("and the ranking difference does not"). Three revision rounds skipped
it without anyone noticing: it changes no number, so nothing in the number checks could
see whether it had been applied. This reads the manuscript for five sentences only the
rewrite contains, and for the absence of the run-on sentence it replaced.
"""
from __future__ import annotations

import os
from pathlib import Path

import pytest

SIGNATURES = [
    "1,781 variants in all",
    "A variant qualifies if it carries",
    "The Nucleotide Transformer is scored as a masked-token",
    "and the ranking difference does not",
    "That is not a design choice",
]
REPLACED = ("which carry no intron offset — the companion classifier files them as coding/UTR — "
            "and are grouped with the core — 1,781 in all")


def _manuscript() -> Path | None:
    """$VARIANT_FM_MANUSCRIPT, else the newest calibration draft on the Desktop, as the other tests find it."""
    env = os.environ.get("VARIANT_FM_MANUSCRIPT")
    if env:
        p = Path(env).expanduser()
        return p if p.exists() else None
    desk = Path.home() / "Desktop"
    hits = sorted([*desk.glob("calibration_draft*.docx"),
                   *(desk / "calibration").glob("calibration_draft*.docx"),
                   *desk.glob("draft_reframed_*.docx")],
                  key=lambda q: q.stat().st_mtime, reverse=True)
    return hits[0] if hits else None


def test_the_prose_pass_is_applied():
    manuscript = _manuscript()
    if manuscript is None:
        pytest.skip("manuscript not available in this checkout")
    from docx import Document

    text = "\n".join(p.text for p in Document(str(manuscript)).paragraphs)
    missing = [s for s in SIGNATURES if s not in text]
    assert not missing, f"sentences of the prose pass missing from the manuscript: {missing}"
    assert REPLACED not in text, "the run-on sentence the prose pass replaced is still in the manuscript"
