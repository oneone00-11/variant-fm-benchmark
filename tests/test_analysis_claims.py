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


# --- the scoped rule, after the global-fallback escape hatch was removed ------

def _backup_with_the_defect():
    """The draft as it stood before 74.5 was corrected, if it is on this machine."""
    root = Path.home() / "Desktop"
    hits = sorted(root.glob("**/draft_reframed_*_backup_*preMerge*.docx"))
    return hits[-1] if hits else None


def test_a_number_absent_from_its_own_declared_source_is_reported(tmp_path, monkeypatch):
    """The regression this rule exists for.

    74.5 was absent from the yield table its paragraph declares, matched exactly
    one unrelated value in the global pool, and passed under the old rule while
    being the wrong number for its own sentence. Synthetic so it runs anywhere.
    """
    import docx as docxlib
    import src.check_manuscript_numbers as cmn

    d = docxlib.Document()
    d.add_paragraph("padding")
    d.add_paragraph("The share moved from 17.7% to 74.5% under the new definition.")
    doc = tmp_path / "m.docx"
    d.save(doc)

    results = tmp_path / "results"
    results.mkdir()
    # declared source holds 17.7 but not 74.5; a second file holds 74.5 once, so
    # the old rule would have called it "a specific source, just not this one"
    (results / "declared.csv").write_text("a\n0.177\n")
    (results / "unrelated.csv").write_text("b\n74.5\n")
    claims = tmp_path / "claims.json"
    claims.write_text(json.dumps({"claims": [
        {"para": "P1", "claim": "share under the new definition",
         "anchor": "The share moved from", "script": None,
         "output": "results/declared.csv", "status": "verified"}]}))

    # scoped paths resolve against the repository root, so point it here
    monkeypatch.setattr(cmn, "REPO", tmp_path)
    r = cmn.classify(doc, results=results, whitelist=tmp_path / "none.json", claims=claims)
    flagged = {t["token"] for t in r["scoped_mismatch"]}
    assert "74.5" in flagged, r["scoped_mismatch"]
    assert "17.7" not in flagged, "a value present in the declared source must pass"


def test_the_rule_reports_it_on_the_actual_pre_fix_draft():
    backup = _backup_with_the_defect()
    if backup is None:
        pytest.skip("the pre-fix draft is not on this machine")
    import src.check_manuscript_numbers as cmn

    r = cmn.classify(backup)
    tokens = {t["token"] for t in r["scoped_mismatch"]}
    assert "74.5" in tokens, sorted(tokens)


def test_a_declared_derivation_licenses_exactly_one_value(tmp_path, monkeypatch):
    """Route (b): a summary of stored values carries its arithmetic. It must
    admit that value and nothing else, so it cannot become a paragraph pass."""
    import docx as docxlib
    import src.check_manuscript_numbers as cmn

    d = docxlib.Document()
    d.add_paragraph("padding")
    d.add_paragraph("There were 30 disagreements and 41 agreements.")
    doc = tmp_path / "m.docx"
    d.save(doc)
    results = tmp_path / "results"
    results.mkdir()
    (results / "declared.csv").write_text("a\n18\n12\n")
    claims = tmp_path / "claims.json"
    claims.write_text(json.dumps({"claims": [
        {"para": "P1", "claim": "disagreements", "anchor": "There were 30 disagreements",
         "script": None, "output": "results/declared.csv", "status": "verified",
         "derived": [{"value": 30, "how": "18 + 12"}]}]}))

    monkeypatch.setattr(cmn, "REPO", tmp_path)
    r = cmn.classify(doc, results=results, whitelist=tmp_path / "none.json", claims=claims)
    flagged = {t["token"] for t in r["scoped_mismatch"]}
    assert "30" not in flagged          # licensed by its declared arithmetic
    assert "41" in flagged              # not licensed, still reported
