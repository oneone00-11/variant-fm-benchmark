"""E20.3 (C13): the VCEP specifications' records and the rules file read from them.

Every record is the file its manifest names; the registry's specifications for the seven
training genes are the ones the listing names, with their status; every quote of the rules
file stands verbatim in the field of the record it cites, or, for ENIGMA's Appendix J
(E20b.2), in both copies of the specification file the quote check read, whose sha256 is the
manifest's; the strengths a rule allows are the record's; and the cut point read stands in
the quote.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

PHASE1 = Path(__file__).resolve().parents[1]
REPO = PHASE1.parent
VCEP = PHASE1 / "data/evidence/vcep"
RAW = VCEP / "raw"


def _json(p: Path):
    return json.loads(p.read_text())


def test_every_record_is_the_one_its_manifest_names():
    man = _json(RAW / "manifest.json")
    for f in man["files"]:
        data = (REPO / f["file"]).read_bytes()
        assert hashlib.sha256(data).hexdigest() == f["sha256"] and len(data) == f["bytes"]
        assert f["http_status"] == 200


def test_the_registry_gives_these_specifications():
    man = _json(RAW / "manifest.json")
    listing = _json(RAW / "cspec_svis.json")["data"]
    for gene, entries in man["found"].items():
        named = [e["@id"].rsplit("/", 1)[1] for e in listing
                 if gene in {g.get("label") for rs in e.get("ruleSets", []) for g in rs.get("genes", [])}]
        assert sorted(named) == sorted(e["id"] for e in entries)
    rules = _json(VCEP / "vcep_rules.json")
    current = {g: v["specification"] for g, v in rules["genes"].items()}
    assert current == {"BAP1": None, "BARD1": None, "BRCA1": "GN092", "BRCA2": "GN097",
                       "PALB2": "GN077", "RAD51C": None, "VHL": "GN078"}
    assert [(e["id"], e["status"]) for e in man["found"]["RAD51C"]] == [("GN096", "Classification Rules In Prep")]


def _appendix_j_found() -> set:
    check = _json(RAW / "enigma_quote_check.json")
    return {q["quote"] for q in check["quotes"] if all(v["found"] for v in q["found"].values())}


def test_every_quote_is_verbatim_in_the_record_it_cites():
    rules = _json(VCEP / "vcep_rules.json")["rules"]
    assert {(r["gene"], r["criterion"]) for r in rules} == {
        (g, c) for g in ("BRCA1", "BRCA2", "PALB2", "VHL") for c in ("PP3", "BP4", "BP7")}
    appendix = _appendix_j_found()
    for r in rules:
        rec = _json(REPO / r["source"]["record"])["data"]
        crit = [c["entContent"] for c in rec["ld"]["CriteriaCode"]
                if c["entContent"].get("label") == r["criterion"]]
        assert len(crit) == 1
        c = crit[0]
        applicable = {s["strength"]: s.get("text") or "" for s in c["strengthDescriptor"]
                      if s.get("applicability") == "Applicable"}
        assert r["strength_allowed"] == list(applicable)
        assert r["quote"] in applicable[r["strength_quoted"]]
        for q in r["further_quotes"]:
            field = q["section"].split(", ", 1)[1]
            if field == "Appendix J":
                assert r["gene"] in ("BRCA1", "BRCA2") and q["quote"] in appendix
                continue
            assert q["quote"] in (c.get(field) or "")
        assert r["version"] == rec["entContent"]["version"]


def test_enigmas_setting_was_read_from_the_files_the_manifest_names():
    """E20b.2: the specification files stay in the gitignored cache (licence unclear); the
    quote check records their sha256, which must be the manifest's, and finds every quote
    in both copies of Appendix J. The rules file records the setting for every ENIGMA rule."""
    spec = _json(RAW / "manifest.json")["enigma_specification"]
    check = _json(RAW / "enigma_quote_check.json")
    cached = {f"{gn} {a['name']}": a["sha256"] for gn, s in spec["specifications"].items()
              for a in s["attachments"] if a.get("cached_at")}
    assert len(cached) == 4 and {k: v["sha256"] for k, v in check["files"].items()} == cached
    assert check["all_found"] and len(_appendix_j_found()) == len(check["quotes"]) == 7
    assert all(f["caption_found"] and len(f["images"]) == 2 for f in check["figure_1a"].values())
    for r in _json(VCEP / "vcep_rules.json")["rules"]:
        stated = r["score_definition_stated"]
        if r["gene"] in ("BRCA1", "BRCA2"):
            assert "stated in the specification's Appendix J" in stated and "Max Distance 10 000" in stated
            assert any(q["section"] == f"{r['criterion']}, Appendix J" for q in r["further_quotes"])
        else:
            assert "Appendix J" not in stated


def test_the_cut_point_read_stands_in_the_quote():
    sign = {">=": "≥", "<=": "≤", ">": ">"}
    for r in _json(VCEP / "vcep_rules.json")["rules"]:
        if r["criterion"] in ("PP3", "BP4"):
            assert f"{sign[r['direction']]}{r['cut']}" in r["quote"].replace("**", "").replace(" ", "")
            assert r["strength_allowed"] == ["Supporting"]
