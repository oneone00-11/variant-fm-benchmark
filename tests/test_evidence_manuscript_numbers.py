"""The manuscript numbers gate (scripts/check_evidence_manuscript_numbers.py).

E5.0, the check before the revision, found three holes in it: its self-test passed
whenever the text had a number it could not find, whether or not the planted numbers
were caught; it compared every number with stored proportions times one hundred, so a
planted 12.34 matched a proportion of 0.12336; and it declared every "95" the confidence
level, "at least 95% of the cases" included. These tests hold the repairs, the reasons the
context rules give, and the cells mode, which holds each number to the one cell a binding
table names: on small texts and tables made here, on the v5 binding table against the
outputs of 72dbcb1 (the commit E5.0 read), and on the v5 texts themselves when they are
named in EVID_V5_MAIN and EVID_V5_SUPPLEMENT (they are not in the repository).
"""
from __future__ import annotations

import csv
import os
import subprocess
import sys
from decimal import Decimal
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import check_evidence_manuscript_numbers as gate  # noqa: E402

V5_TABLE = REPO / "docs" / "manuscript_bindings_v5.tsv"
V5_EXPECTED = REPO / "tests" / "fixtures" / "manuscript_bindings_v5_expected.tsv"
E5_0_COMMIT = "72dbcb1d7418b168e84e5311b3f7ebc21711606f"     # the outputs E5.0 read


def _md(body: str) -> str:
    return f"# Title page\n\nNot read 999.\n\n# Abstract\n\n{body}\n\n# References\n\nNot read 888.\n"


def _reasons(sentence: str) -> dict[str, str | None]:
    return {t.printed: t.reason for t in gate.read_tokens(_md(sentence))}


# ---------------------------------------------------------------- declared tokens
@pytest.mark.parametrize("sentence, printed, reason", [
    ("A ratio of 8.3 (95% CI 6.3 to 11.0) was seen.", "95%", "confidence level"),
    ("A ratio of 8.3 (95% confidence interval 6.3–11.0) was seen.", "95%", "confidence level"),
    ("The one-sided 95% lower bound clears the boundary.", "95%", "confidence level"),
    ("Its lower edge, the one-sided 95% bound, defines each threshold.", "95%", "confidence level"),
    ("Confidence intervals (CI) are 95% and come from resamples.", "95%", "confidence level"),
    ("We ask (1) what and (2) whether.", "1", "inline enumeration, (1) (2) (3)"),
    ("The set was superseded on 11 September 2026 by a new version.", "11",
     "day of the date 11 September 2026"),
    ("Variants at 3–10 bp were kept.", "10", "stratum bound, the upper end of 3-10 nt"),
    ("Variants at 11–50 bp were kept.", "50", "stratum bound"),
    ("It was scored at a distance of 50 nucleotides.", "50",
     "SpliceAI scoring distance of 50 nt (the atlas column, the fusion input)"),
    ("The fusion takes SpliceAI at distance 50 and Pangolin.", "50",
     "SpliceAI scoring distance of 50 nt (the atlas column, the fusion input)"),
    ("Variants above a filtering allele frequency of 0.1% were set apart.", "0.1%",
     "gnomAD filtering allele frequency of 0.1% (the Atlas model's labels)"),
    ("BP4 applies at a score of 0.1 or below.", "0.1", "BP4 cut point, Walker et al. Table 1 (config/walker2023.yaml)"),
])
def test_a_constant_is_declared_with_the_reason_of_its_context(sentence, printed, reason):
    assert _reasons(sentence)[printed] == reason


def test_the_score_scale_anchors_are_not_stratum_bounds():
    toks = gate.read_tokens(_md("The score is anchored at −1 for synonymous and +1 for nonsense controls."))
    assert [(t.printed, t.reason) for t in toks] == \
        [("1", "anchor of the score scale, -1 synonymous and +1 nonsense controls")] * 2
    assert _reasons("Variants at ±1–2 nt are canonical.")["1"] == "stratum bound"


def test_95_measured_is_a_measurement():
    """v5 line 89: "in at least 95% of evaluable cases" is a measured lower bound."""
    assert _reasons("The genes cleared the boundary in at least 95% of evaluable cases.")["95%"] is None
    assert _reasons("Of the variants, 95 were damaging.")["95"] is None


# ---------------------------------------------------------------- percentages
def test_a_proportion_times_one_hundred_matches_only_a_printed_percentage():
    pairs = [(0.12336, "walker_thresholds.csv:posterior_bp4_prior0.2")]
    assert gate._match_value("12.34", pairs, percent=False) is None
    assert gate._match_value("12.34", pairs, percent=True) is not None
    toks = gate.read_tokens(_md("It was 12.34% and 12.34 per cent and 12.34 here."))
    assert [t.context.percent for t in toks] == [True, True, False]


def _toy_dirs(monkeypatch, tmp_path: Path, csv_text: str) -> None:
    (tmp_path / "toy.csv").write_text(csv_text)
    monkeypatch.setattr(gate, "SEARCH_DIRS", [tmp_path])


def test_the_bound_mode_no_longer_matches_a_plain_number_to_a_proportion(monkeypatch, tmp_path):
    _toy_dirs(monkeypatch, tmp_path, "tool,share\nspliceai,0.12336\n")
    res = gate.bound_results(_md("The share was 12.34 [toy.csv]. The share was 12.34% [toy.csv]."))
    assert [t.printed for t in res.fails] == ["12.34"]
    assert [t.printed for t in res.bound] == ["12.34%"]


# ---------------------------------------------------------------- self-test
def _selftest_md(tmp_path: Path, body: str) -> Path:
    p = tmp_path / "ms.md"
    p.write_text(_md(body))
    return p


def test_the_selftest_counts_each_planted_number(monkeypatch, tmp_path, capsys):
    _toy_dirs(monkeypatch, tmp_path, "tool,share\nspliceai,0.12336\nsplice,0.5\n")
    # the text's own number fails, which no longer decides the verdict
    md = _selftest_md(tmp_path, "The share was 0.5 [toy.csv]. A count of 4,321 [toy.csv].")
    assert gate.selftest(md) == 0
    out = capsys.readouterr().out
    assert "SELFTEST caught 6 / 6: PASSED" in out and "MISSED" not in out


def test_the_selftest_refuses_a_number_the_named_file_holds(monkeypatch, tmp_path, capsys):
    _toy_dirs(monkeypatch, tmp_path, "tool,share\nspliceai,0.3712\n")
    md = _selftest_md(tmp_path, "The share was 0.3712 [toy.csv].")
    assert gate.selftest(md) == 2
    assert "cannot run" in capsys.readouterr().out


def test_the_selftest_refuses_a_declared_constant(monkeypatch, tmp_path, capsys):
    _toy_dirs(monkeypatch, tmp_path, "tool,share\nspliceai,0.5\n")
    md = _selftest_md(tmp_path, "The share was 0.5 [toy.csv].")
    assert gate.selftest(md, planted=["0.3712", "0.2"]) == 2
    assert "declared constants 0.2" in capsys.readouterr().out


def test_the_selftest_fails_when_a_planted_number_is_missed(monkeypatch, tmp_path, capsys):
    """A missed number fails it even when the text's own numbers fail too: the old verdict
    was PASSED whenever anything at all was reported."""
    _toy_dirs(monkeypatch, tmp_path, "tool,share\nspliceai,0.5\n")
    md = _selftest_md(tmp_path, "The share was 0.5 [toy.csv]. A count of 4,321 [toy.csv].")
    real = gate.bound_results

    def hides_one(text):
        res = real(text)
        res.fails = [t for t in res.fails if t.printed != "987"]
        return res

    monkeypatch.setattr(gate, "bound_results", hides_one)
    assert gate.selftest(md) == 1
    out = capsys.readouterr().out
    assert "MISSED 987" in out and "SELFTEST caught 5 / 6: FAILED" in out


def test_the_selftest_needs_a_sentence_naming_a_table(tmp_path):
    md = _selftest_md(tmp_path, "No file is named here, 0.5.")
    assert gate.selftest(md) == 2


# ---------------------------------------------------------------- cells mode
TOY_TEXT = ("The ratio was 8.3 [t.csv]. It was 8.3 again in 12% of 1,234 [t.csv]. "
            "The two ends were 0.5 and 0.9 [t.csv]. All round to 2.2 [t.csv]. "
            "There are 3 rows and 2 set [t.csv]. The record says 16,384 [r.json]. "
            "The sets give 71 in all and at least 0.99 [r.json].")


def _toy_root(tmp_path: Path) -> Path:
    root = tmp_path / "repo"
    (root / "out").mkdir(parents=True)
    (root / "out" / "t.csv").write_text(
        "tool,stratum,lr,share,n,lo\n"
        "a,s3_10,8.3,0.12,1234,2.21\n"
        "b,s3_10,8.26,0.5,10,2.24\n"
        "a,s11_50,7.0,0.9,0,\n")
    (root / "out" / "r.json").write_text(
        '{"code": "resize(SEQUENCE_LENGTH_1MB)", "sets": {"urn:x:q-1": {"n": 40, "auc": 0.51}, '
        '"urn:x:p-1": {"n": 30, "auc": 0.6}, "urn:x:r-1": {"n": 1, "auc": 0.993}}}')
    return root


def _binding(line=7, printed="8.3", occurrence=1, excerpt="", kind="measured", reason="",
             file="out/t.csv", rows="", column="lr", derive="value", scale=1) -> gate.Binding:
    return gate.Binding("main", line, occurrence, printed, excerpt, kind, reason, file, rows,
                        column, derive, scale)


def _toy_bindings() -> list[gate.Binding]:
    return [
        _binding(printed="8.3", rows="tool=a; stratum=s3_10", excerpt="The ratio was 8.3"),
        _binding(printed="8.3", occurrence=2, rows="tool=b"),                      # 8.26 -> 8.3
        _binding(printed="12%", rows="tool=a; stratum=s3_10", column="share", scale=100),
        _binding(printed="1,234", column="n", derive="max"),
        _binding(printed="0.5", rows="lr>=8", column="share", derive="max"),     # max of 0.12, 0.5
        _binding(printed="0.9", rows="stratum not in s3_10", column="share"),
        _binding(printed="2.2", rows="stratum in s3_10|s3_50", column="lo", derive="all"),
        _binding(printed="3", derive="count", column=""),
        _binding(printed="2", column="lo", derive="count_nonnull"),
        _binding(printed="16,384", file="out/r.json", column="code", derive="text"),
        _binding(printed="71", file="out/r.json", column="sets.*.n", derive="sum"),
        _binding(printed="0.99", file="out/r.json", column="sets.*.auc", rows="key not in urn:x:q-1|urn:x:p-1",
                 derive="min"),
    ]


def test_the_cells_mode_reads_each_kind_of_binding(tmp_path):
    root = _toy_root(tmp_path)
    res = gate.cells_results(_md(TOY_TEXT), _toy_bindings(), root)
    got = {(r.binding.printed, r.binding.occurrence): (r.status, r.now) for r in res}
    assert got == {
        ("8.3", 1): ("match", "8.3"), ("8.3", 2): ("match", "8.3"), ("12%", 1): ("match", "12%"),
        ("1,234", 1): ("match", "1,234"), ("0.5", 1): ("match", "0.5"), ("0.9", 1): ("match", "0.9"),
        ("2.2", 1): ("match", "2.2"), ("3", 1): ("match", "3"), ("2", 1): ("match", "2"),
        ("16,384", 1): ("mismatch", "(text) resize(SEQUENCE_LENGTH_1MB)"),
        ("71", 1): ("match", "71"), ("0.99", 1): ("match", "0.99"),
    }


def test_a_cell_is_compared_at_the_printed_precision(tmp_path):
    root = _toy_root(tmp_path)
    b = _binding(printed="8.30", rows="tool=b")                                  # 8.26 is not 8.30
    t = [x for x in gate.read_tokens(_md("It was 8.30 [t.csv].")) if x.printed == "8.30"]
    b.line = t[0].line
    r = gate.cells_results(_md("It was 8.30 [t.csv]."), [b], root)[0]
    assert (r.status, r.now) == ("mismatch", "8.26")
    assert gate.rounds_to(Decimal("0.125"), "0.12") and gate.rounds_to(Decimal("0.125"), "0.13")
    assert not gate.rounds_to(Decimal("0.1251"), "0.12")


@pytest.mark.parametrize("change, why", [
    (dict(file="out/none.csv"), "no file out/none.csv"),
    (dict(column="nope"), "no column 'nope'"),
    (dict(rows="nope=a"), "no column 'nope'"),
    (dict(rows="stratum=s3_10"), "the rows select 2 cells, not one"),
    (dict(rows="tool=c"), "the rows select 0 cells, not one"),
    (dict(rows="tool=a; stratum=s3_10", column="share", scale=100), "scale 100 on a number not printed"),
    (dict(excerpt="not this sentence"), "the excerpt is not in the sentence"),
    (dict(kind="maybe"), "is neither constant nor measured"),
    (dict(derive="mean"), "is not one of"),
    (dict(file="out/r.json", rows="", column='sets["urn:x:z-1"].n'), "no key 'urn:x:z-1'"),
    (dict(file="out/r.json", rows="tool=a", column='sets.*.n'), "the rows can only test the key"),
    (dict(file="out/r.json", rows="key in urn:x:q-1", column='sets["urn:x:q-1"].n'), "the path has no '*'"),
])
def test_a_row_that_cannot_be_resolved_is_reported(tmp_path, change, why):
    root = _toy_root(tmp_path)
    b = _binding(rows="tool=a; stratum=s3_10")
    for k, v in change.items():
        setattr(b, k, v)
    r = gate.cells_results(_md("The ratio was 8.3 [t.csv]."), [b], root)[0]
    assert r.status == "cannot_bind" and why in r.why


def test_a_condition_splits_at_its_first_operator(tmp_path):
    """A value may hold another operator's word or sign: "note = present in ClinVar"
    tests the column note."""
    assert gate._conditions("note = present in ClinVar; tool in a|b; label startswith x = y") == [
        ("note", "=", "present in ClinVar"), ("tool", "in", "a|b"), ("label", "startswith", "x = y")]
    assert gate._conditions("n>=3; kind not in a|b; lo=") == [
        ("n", ">=", "3"), ("kind", "not in", "a|b"), ("lo", "=", "")]
    for bad in ("= a", "tool a", "tool in"):
        with pytest.raises(gate.BindError, match="cannot read the row condition"):
            gate._conditions(bad)
    root = tmp_path / "repo"
    (root / "out").mkdir(parents=True)
    (root / "out" / "t.csv").write_text("note,lr\npresent in ClinVar,8.3\nabsent,7.0\n")
    r = gate.cells_results(_md("The ratio was 8.3 [t.csv]."),
                           [_binding(rows="note = present in ClinVar")], root)[0]
    assert (r.status, r.now) == ("match", "8.3")


def test_numbers_without_a_row_and_rows_without_a_number_are_reported(tmp_path):
    root = _toy_root(tmp_path)
    stale = _binding(line=40, printed="7.0")
    res = gate.cells_results(_md("The ratio was 8.3 [t.csv]."), [stale], root)
    assert sorted(r.status for r in res) == ["cannot_bind", "not_in_table"]
    assert "no number 7.0 (occurrence 1) on line 40" in [r.why for r in res if r.binding][0]


def test_the_cells_mode_fails_on_each_list(tmp_path, capsys):
    root = _toy_root(tmp_path)
    md = tmp_path / "ms.md"
    md.write_text(_md("The ratio was 8.3 [t.csv]."))
    table = tmp_path / "b.tsv"
    good = dict(document="main", line=7, occurrence=1, printed="8.3", excerpt="ratio was 8.3", kind="measured",
                reason="", file="out/t.csv", rows="tool=a; stratum=s3_10", column="lr", derive="value", scale=1)
    for rows, rc in (([good], 0),
                     ([dict(good, rows="tool=b")], 0),                     # 8.26 prints as 8.3
                     ([dict(good, column="share", rows="tool=b")], 1),     # 0.5
                     ([dict(good, file="out/x.csv")], 1),
                     ([dict(good, kind="constant", reason="a parameter", file="")], 0),
                     ([], 1)):
        gate.write_bindings(table, rows)
        assert gate.cells_check(md, table, "main", root=root) == rc, rows
    out = capsys.readouterr().out
    assert "MISMATCH line 7: 8.3, now 0.5" in out and "CANNOT BIND line 7" in out and "NOT IN TABLE line 7" in out


# ---------------------------------------------------------------- the v5 table
def _expected() -> dict:
    with open(V5_EXPECTED, newline="") as fh:
        return {(r["document"], int(r["line"]), r["printed"], int(r["occurrence"])): r
                for r in csv.DictReader(fh, delimiter="\t")}


def test_the_v5_table_has_a_valid_row_for_each_number_e5_0_read():
    rows = gate.read_bindings(V5_TABLE)
    exp = _expected()
    assert len(rows) == len(exp) == 182
    assert {(b.document, b.line, b.printed, b.occurrence) for b in rows} == set(exp)
    for b in rows:
        assert b.kind in ("constant", "measured")
        if b.kind == "constant":
            assert b.reason and not b.file
        else:
            assert b.file and b.derive in gate.DERIVES and b.scale in (1, 100)
            assert b.scale == 1 or b.printed.endswith("%")
    by_doc = {d: [b for b in rows if b.document == d] for d in ("main", "supplement")}
    assert (len(by_doc["main"]), len(by_doc["supplement"])) == (120, 62)


def _v5_now(res: gate.CellResult, want: str) -> bool:
    if res.now == want:
        return True
    try:
        nums = [Decimal(x.strip().replace(",", "").replace("%", "")) for x in want.split("–")]
    except Exception:  # noqa: BLE001 - E5.0 wrote words: the status is what is compared
        return True
    got = [min(res.values), max(res.values)] if len(nums) == 2 else res.values[:1]
    places = [Decimal(1).scaleb(-(len(str(n).split(".")[1]) if "." in str(n) else 0)) for n in nums]
    return all((g * res.binding.scale).quantize(q) == n for g, n, q in zip(got, nums, places))


def test_the_v5_table_reproduces_e5_0_on_the_outputs_it_read(tmp_path):
    """Every binding of v5, read from 72dbcb1's files: the verdict E5.0 §2 reached on each
    number, and E5.0's value where it does not match (main text 19 mismatches, E5.0's 14
    and the 5 its bound mode could not find; supplementary notes 8, E5.0's 3 and the 5 of
    its 10 not found that have changed)."""
    if subprocess.run(["git", "-C", str(REPO), "cat-file", "-e", f"{E5_0_COMMIT}^{{commit}}"],
                      capture_output=True).returncode != 0:
        pytest.skip(f"E5.0's commit {E5_0_COMMIT[:7]} is not in this clone")
    rows = gate.read_bindings(V5_TABLE)
    for f in sorted({b.file for b in rows if b.kind == "measured"}):
        blob = subprocess.run(["git", "-C", str(REPO), "show", f"{E5_0_COMMIT}:{f}"],
                              capture_output=True, check=True).stdout
        (tmp_path / f).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / f).write_bytes(blob)
    exp = _expected()
    tally: dict[tuple, int] = {}
    for b in rows:
        r = gate.evaluate(b, None, tmp_path)
        e = exp[(b.document, b.line, b.printed, b.occurrence)]
        assert r.status == e["expected_status"], (b, r.why)
        if r.status == "mismatch":
            assert _v5_now(r, e["expected_now"]), (b, r.now, e["expected_now"])
        tally[(b.document, r.status)] = tally.get((b.document, r.status), 0) + 1
    assert tally == {("main", "constant"): 53, ("main", "match"): 48, ("main", "mismatch"): 19,
                     ("supplement", "constant"): 23, ("supplement", "match"): 31,
                     ("supplement", "mismatch"): 8}


@pytest.mark.parametrize("doc, env", [("main", "EVID_V5_MAIN"), ("supplement", "EVID_V5_SUPPLEMENT")])
def test_the_v5_texts_read_number_by_number(doc, env):
    """With the v5 texts at hand, every number of each has its row, its excerpt and its
    verdict: nothing NOT IN TABLE and nothing that CANNOT BIND."""
    path = os.environ.get(env)
    if not path or not Path(path).expanduser().is_file():
        pytest.skip(f"set {env} to the v5 {doc} text to run this")
    if subprocess.run(["git", "-C", str(REPO), "cat-file", "-e", f"{E5_0_COMMIT}^{{commit}}"],
                      capture_output=True).returncode != 0:
        pytest.skip(f"E5.0's commit {E5_0_COMMIT[:7]} is not in this clone")
    rows = [b for b in gate.read_bindings(V5_TABLE) if b.document == doc]
    text = Path(path).expanduser().read_text()
    head = subprocess.run(["git", "-C", str(REPO), "diff", "--quiet", E5_0_COMMIT, "--",
                           *sorted({b.file for b in rows if b.kind == "measured"})]).returncode
    if head != 0:
        pytest.skip("the bound files differ from 72dbcb1's here")
    res = gate.cells_results(text, rows)
    exp = _expected()
    assert not [r for r in res if r.status in ("cannot_bind", "not_in_table")]
    for r in res:
        b = r.binding
        assert r.status == exp[(doc, b.line, b.printed, b.occurrence)]["expected_status"]


def test_the_pool_mode_compares_percentages_only_when_printed(monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(gate, "pool_values", lambda: {0.12336: ["toy.csv:share"]})
    md = tmp_path / "ms.md"
    md.write_text(_md("The share was 12.34. The share was 12.34%."))
    assert gate.pool_check(md) == 1
    out = capsys.readouterr().out
    assert "tokens: 1 matched, 0 declared, 1 unmatched" in out and "UNMATCHED 12.34" in out
