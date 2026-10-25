"""E19 C3, BRCA2: the second set's records, label rule and mapping (before any predictor score).

Every record is the file its manifest names; the screening chooses exactly one score set
and gives each other its failing conditions; the paper's rule is quoted and checked and
does not apply (R4), so the classes are Supplementary Table 5's (R6), whose extract is the
table's and whose counts, grouped as the paper counts them, are the paper's; every row of
the table is checked against its own g., with exactly the 28 miswritten c. and the 4
dropped rows the author allowed; and the coverage record adds up both ways.
"""
from __future__ import annotations

import csv
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

import pytest

PHASE1 = Path(__file__).resolve().parents[1]
REPO = PHASE1.parent
RAW = PHASE1 / "data/evidence/second_set/brca2_raw"

if str(PHASE1) not in sys.path:
    sys.path.insert(0, str(PHASE1))

from src import evid_label_sources as LS  # noqa: E402
from src import evid_second_set as SS  # noqa: E402


def _json(name: str) -> dict:
    return json.loads((RAW / name).read_text())


def test_every_record_is_the_one_its_manifest_names():
    for f in _json("manifest.json")["files"]:
        data = (REPO / f["file"]).read_bytes()
        assert hashlib.sha256(data).hexdigest() == f.get("sha256_gz", f["sha256"])
    man = _json("label_source_manifest.json")
    assert [e for e in man["files"] if e.get("cached")][0]["path"].startswith("phase1/data/e19_cache/")
    differ, absent = LS.extract_differences(man, REPO)
    assert not differ


def test_the_screening_chooses_one_score_set_and_says_why_not_the_others():
    with open(RAW / "screening.csv", newline="") as f:
        rows = {r["urn"]: r for r in csv.DictReader(f)}
    assert len(rows) == 10
    meets = [u for u, r in rows.items() if r["meets_all"] == "True"]
    assert meets == ["urn:mavedb:00001242-a-1"]
    assert rows["urn:mavedb:00001242-a-1"]["d_rule_used"] == "R6"
    fails = {u: r["fails"] for u, r in rows.items()}
    assert "(c)" in fails["urn:mavedb:00001263-a-2"]
    assert all("(c)" in fails[f"urn:mavedb:00001268-{x}"] for x in ("0-1", "a-1", "b-1", "c-1"))
    assert "(a)" in fails["urn:mavedb:00001224-a-1"] and "(a)" in fails["urn:mavedb:00001223-a-1"]
    assert "(b)" in fails["urn:mavedb:00001279-a-1"] and "(b)" in fails["urn:mavedb:00001225-a-1"]
    d = _json("decision.json")
    assert d["outcome"] == "use urn:mavedb:00001242-a-1"
    assert d["item_7"]["investigator_provided_calibration"] is False


def test_the_paper_rule_is_quoted_checked_and_does_not_apply():
    rule, check = _json("paper_rule.json"), _json("paper_quote_check.json")
    assert rule["status"] == "not_stated" and rule["supplementary_rule"]["rule"] == "R6"
    assert check["paper_rule_sha256"] == hashlib.sha256((RAW / "paper_rule.json").read_bytes()).hexdigest()
    assert len(check["quotes"]) == len(rule["quotes"]) and all(q["found"] for q in check["quotes"])
    assert check["pdf"]["sha256"] == rule["source"]["pdf"]["sha256"]
    assert "ACMG.class.AVENGERS" in rule["supplementary_rule"]["table"]["columns_not_used"]


def test_the_classes_grouped_as_the_paper_counts_them_are_the_papers():
    rule = _json("paper_rule.json")
    sup = rule["supplementary_rule"]
    rows = LS.read_tsv(RAW / "supplementary_table5_classes.tsv")
    assert set(r["category"] for r in rows) == set(sup["groups"])
    assert dict(Counter(sup["paper_groups"][r["category"]] for r in rows)) == \
        rule["reported_counts"]["counts"]
    assert len(rows) == rule["reported_counts"]["total"]
    assert set(sup["groups"].values()) == {"abnormal", "normal", "intermediate"}


def test_every_row_is_checked_against_its_own_g(monkeypatch):
    monkeypatch.chdir(PHASE1)
    kept, check = SS.table5_rows()
    assert (check["rows"], check["kept"]) == (6551, 6547)
    assert len(check["c_miswritten"]) == SS.N_MISWRITTEN == 28
    assert sorted(x["hgvs_nt"].partition(":")[2] for x in check["dropped"]) == sorted(SS.DROPPED)
    written = {x["hgvs_nt"].partition(":")[2] for x in check["c_miswritten"]}
    assert {"c.9505A>C", "c.9505A>G", "c.9505A>T", "c.8633-2-1A>T"} <= written
    assert len({r["key"] for r in kept}) == len(kept)


def test_the_coverage_record_adds_up_both_ways():
    with open(RAW / "coverage_3_10.csv", newline="") as f:
        cov = list(csv.DictReader(f))
    total = lambda c: sum(int(r[c]) for r in cov)  # noqa: E731
    assert total("n_ours") == 551
    assert (total("n_ours_in_table5"), total("n_ours_in_mavedb_score_set")) == (201, 198)
    plus4 = next(r for r in cov if r["offset"] == "+4")
    assert int(plus4["n_ours_in_table5"]) - int(plus4["n_ours_in_mavedb_score_set"]) == 3
    with open(RAW / "screening.csv", newline="") as f:
        scr = {r["urn"]: r for r in csv.DictReader(f)}
    assert int(scr["urn:mavedb:00001242-a-1"]["c_shared_3_10"]) == 198
