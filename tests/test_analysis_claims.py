"""Every analysis the manuscript claims must have code and output behind it.

A value check asks whether a number is right. It cannot ask whether the analysis
a sentence describes exists at all. The sampling-frame paragraph reported a
whole-gene enrichment of 0.713 versus 0.578; both numbers were plausible, both
were near real values elsewhere in the pipeline, and the comparison had never
been computed -- ipw_frame_balance.csv is restricted to the splice window.

Claims bind by content anchor, never by paragraph number: rewriting the abstract
shifted every later paragraph by two and left all six `para` values pointing at
the wrong text, with nothing to notice.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "phase1"))

MANIFEST = REPO / "phase1" / "config" / "analysis_claims.json"
def _manuscript() -> Path | None:
    """The .docx these claims bind to, if it is on this machine.

    Resolved from $VARIANT_FM_MANUSCRIPT, else the newest draft_reframed_*.docx
    on the Desktop. The filename is not hard-coded: it carries a working title
    that is nobody's business but the author's, and the tests only need some
    draft to check the anchors against.
    """
    env = os.environ.get("VARIANT_FM_MANUSCRIPT")
    if env:
        p = Path(env).expanduser()
        return p if p.exists() else None
    hits = sorted((Path.home() / "Desktop").glob("draft_reframed_*.docx"),
                  key=lambda q: q.stat().st_mtime, reverse=True)
    return hits[0] if hits else None


MANUSCRIPT = _manuscript()
VALID_STATUS = {"verified", "narrative", "gap"}


def _claims() -> list[dict]:
    assert MANIFEST.exists(), f"{MANIFEST.relative_to(REPO)} is missing"
    return json.loads(MANIFEST.read_text())["claims"]


def _outputs(c) -> list[str]:
    out = c.get("output")
    return [] if not out else (out if isinstance(out, list) else [out])


def test_manifest_is_well_formed():
    claims = _claims()
    assert claims, "no claims recorded"
    for c in claims:
        assert c.get("claim", "").strip(), c
        assert c.get("status") in VALID_STATUS, c
        assert c.get("anchor"), f"no content anchor: {c['claim']}"


def test_every_named_script_and_output_exists():
    missing = [f"{c['claim']}: {rel}" for c in _claims()
               for rel in ([c["script"]] if c.get("script") else []) + _outputs(c)
               if not (REPO / rel).exists()]
    assert not missing, "claims naming files that do not exist:\n  " + "\n  ".join(missing)


def test_every_anchor_hits_exactly_once():
    """Zero hits means the sentence was reworded and the row now describes
    nothing; more than one means the anchor could bind to the wrong sentence.
    Both are invisible under paragraph numbers, which always resolve to *some*
    paragraph however far the text has moved."""
    if MANUSCRIPT is None:
        pytest.skip("manuscript not available in this checkout")
    from src.check_manuscript_numbers import _blocks

    blocks = _blocks(MANUSCRIPT)
    bad = [f"{c['para']}: {sum(t.count(c['anchor']) for _, t in blocks)} hits "
           f"for {c['anchor']!r}" for c in _claims()
           if sum(t.count(c["anchor"]) for _, t in blocks) != 1]
    assert not bad, "anchors that do not hit exactly once:\n  " + "\n  ".join(bad)


def test_paragraph_numbers_are_not_load_bearing(tmp_path):
    """Regression for the abstract rewrite: with every `para` deliberately
    wrong, resolution must be unaffected, because nothing matches on it."""
    if MANUSCRIPT is None:
        pytest.skip("manuscript not available in this checkout")
    import src.check_manuscript_numbers as cmn

    before = cmn.resolve_anchors(MANUSCRIPT)
    spec = json.loads(MANIFEST.read_text())
    for c in spec["claims"]:
        c["para"] = f"P{int(c['para'][1:]) + 23}"
    shifted = tmp_path / "shifted.json"
    shifted.write_text(json.dumps(spec))
    assert before == cmn.resolve_anchors(MANUSCRIPT, shifted) and before
