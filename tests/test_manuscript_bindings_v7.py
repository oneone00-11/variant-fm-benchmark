"""E9.3: the binding table of the submitted text (docs/manuscript_bindings_v7.tsv) and the
scripts that build and check it (scripts/stage_E9/build_bindings_v7.py,
check_bindings_v7.py). The submitted files are not in the repository, so what is held here
needs none of them: the table against the outputs (every measured number is its one cell at
the printed precision, every constant has its reason), its difference from v6's, and the
scripts' own steps on small texts made here."""
from __future__ import annotations

import csv
import sys
from collections import Counter
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts" / "stage_E9"))
sys.path.insert(0, str(REPO / "scripts"))
import build_bindings_v7 as B  # noqa: E402
import check_bindings_v7 as C  # noqa: E402
import check_evidence_manuscript_numbers as gate  # noqa: E402

V6, V7 = REPO / "docs" / "manuscript_bindings_v6.tsv", REPO / "docs" / "manuscript_bindings_v7.tsv"


def _holds(b) -> bool:
    """A measured row against its cell without the text: the gate's own comparison
    (gate.evaluate), which reads from the sentence only whether the number is printed with
    a minus sign, so either sign is taken here."""
    cells = gate.resolve(b, REPO)
    if b.derive == "text":
        return gate._number_in_text(b.printed, cells[0])
    values = cells if b.derive == "all" else cells[:1]
    return all(any(gate.rounds_to(v, b.printed, b.scale, minus) for minus in (False, True)) for v in values)


def test_every_measured_number_is_its_cell_and_every_constant_has_a_reason():
    bindings = gate.read_bindings(V7)
    kinds = Counter((b.document, b.kind) for b in bindings)
    assert kinds == {("main", "measured"): 113, ("main", "constant"): 90,
                     ("supplement", "measured"): 7, ("supplement", "constant"): 29}
    for b in bindings:
        if b.kind == "measured":
            assert _holds(b), (b.line, b.printed, b.where())
        else:
            assert b.reason and not b.file, (b.line, b.printed)
    # a row for each number, in reading order: lines never go back within a document
    for doc in ("main", "supplement"):
        lines = [b.line for b in bindings if b.document == doc]
        assert lines == sorted(lines)


def test_v7_differs_from_v6_where_the_submitted_text_does():
    """The legends' constants that went when the legends were shortened; the title page's
    counts and the statements' numbers, read for the first time; and no number of the
    abstract, the main text or the notes with another value or another cell."""
    d = C.difference(C.table_rows(V6), C.table_rows(V7))
    assert (d["numbers_before"], d["numbers_after"]) == (235, 239)
    assert sorted(x["printed"] for x in d["vanished"]) == ["0.2", "0.2", "18.7", "2,000", "4.33"]
    assert {x["kind"] for x in d["vanished"]} == {"constant"}
    assert sorted(x["printed"] for x in d["new"]) == ["1.2", "11", "11", "15", "247", "3", "3,934", "4", "42"]
    assert {x["kind"] for x in d["new"]} == {"constant"} and {x["document"] for x in d["new"]} == {"main"}
    assert sum("counted from the manuscript" in x["entry"] for x in d["new"]) == 5
    assert sorted(x["before"]["printed"] for x in d["reworded"]) == ["0.1", "0.2", "10", "3"]
    assert all(x["before"]["entry"] == x["after"]["entry"] for x in d["reworded"])
    assert d["kept_in_place"] + len(d["moved"]) + len(d["reworded"]) + len(d["vanished"]) == 235
    assert not [x for x in d["moved"] if x["document"] == "supplement"]
    # every measured number of v6 is in v7 with the cell it had
    cell = lambda r: tuple(r[k] for k in ("document", "printed", "file", "rows", "column", "derive", "scale"))   # noqa: E731
    before = Counter(cell(r) for r in C.table_rows(V6) if r["kind"] == "measured")
    after = Counter(cell(r) for r in C.table_rows(V7) if r["kind"] == "measured")
    assert before == after and sum(after.values()) == 120


def _rows(*items) -> list[dict]:
    blank = {"reason": "", "file": "", "rows": "", "column": "", "derive": "", "scale": ""}
    return [{"document": "main", "line": str(line), "occurrence": "1", "printed": printed, "excerpt": excerpt,
             "kind": kind} | blank | entry for line, printed, excerpt, kind, entry in items]


def test_the_difference_of_two_tables_names_what_vanished_came_was_reworded_or_moved():
    reason, cell = {"reason": "a bound"}, {"file": "f.csv", "column": "n", "derive": "value", "scale": "1"}
    before = _rows((3, "10", "within 10 bp of", "constant", reason), (5, "42", "of 42 genes", "measured", cell),
                   (7, "0.2", "a cut point of 0.2 on", "constant", reason), (9, "7", "in 7 folds", "measured", cell),
                   (11, "2,000", "of 2,000 draws", "constant", reason))
    after = _rows((3, "10", "within 10 bp of", "constant", reason), (5, "7", "in 7 folds", "measured", cell),
                  (7, "42", "of 42 genes", "measured", cell), (9, "0.2", "the cut point, 0.2, on", "constant", reason),
                  (13, "5", "5 figures", "constant", {"reason": "counted"}))
    d = C.difference(before, after)
    assert [x["printed"] for x in d["vanished"]] == ["2,000"] and [x["printed"] for x in d["new"]] == ["5"]
    assert [(x["before"]["excerpt"], x["after"]["excerpt"]) for x in d["reworded"]] == \
        [("a cut point of 0.2 on", "the cut point, 0.2, on")]
    assert len(d["moved"]) == 1 and d["moved"][0]["printed"] in ("42", "7") and d["kept_in_place"] == 2
    same = C.difference(before, before)
    assert (same["kept_in_place"], same["moved"], same["reworded"], same["vanished"], same["new"]) == (5, [], [], [], [])


STRUCTURE = '''Document
  H2 (block):
    "Supplementary Notes"
  H3 (block):
    "Supplementary Note 1. Scores"
  P (block):
     /Placement /Block
    Span (inline)
      "SpliceAI was run at version 1.3.1, a splice-"
      "region predictor, "
      "within 50 bp."
  L (block):
    LI (block):
      Lbl (block):
        "1."
      LBody (block):
        P (block):
          "The first item of 3. "
    LI (block):
      Lbl (block):
        "2."
      LBody (block):
        P (block):
          "The second item."
  H3 (block):
    "References cited in the Supplementary Notes"
  P (block):
    "21. Zhang N. A reference of 2026;12:34–56."
  P (block):
    "Supplementary Figure S1. A caption."
'''


def test_the_notes_are_read_from_the_pdfs_structure_without_the_list_labels(tmp_path, monkeypatch):
    blocks = B.struct_blocks(STRUCTURE)
    assert [t for t, _ in blocks] == ["H2", "H3", "P", "Lbl", "P", "Lbl", "P", "H3", "P", "P"]
    assert blocks[2][1] == "SpliceAI was run at version 1.3.1, a splice-region predictor, within 50 bp."
    structure = tmp_path / "structure.txt"
    structure.write_text(STRUCTURE)
    monkeypatch.setattr(B, "confirm_in_pdf", lambda pdf, text: {"list labels between them": 2})
    text, notes = B.supplement_text(tmp_path / "any.pdf", structure)
    assert text.endswith("\n\n") and text.rstrip("\n").split("\n\n") == [
        "## Supplementary Notes", "### Supplementary Note 1. Scores",
        "SpliceAI was run at version 1.3.1, a splice-region predictor, within 50 bp.",
        "The first item of 3.", "The second item.", "# References cited in the Supplementary Notes"]
    assert notes["list labels left out"] == 2 and notes["entries of the notes' reference list"] == 1
    assert notes["numbers in the notes' reference list"] == 5
    # the gate reads the notes' numbers (a version such as 1.3.1 is not one) and stops at their references
    assert [t.printed for t in gate.read_tokens(text)] == ["50", "3"]
    # a count of labels in the pages that is not the structure's stops the run
    monkeypatch.setattr(B, "confirm_in_pdf", lambda pdf, text: {"list labels between them": 1})
    with pytest.raises(SystemExit, match="2 list labels in the structure, 1 in the pages"):
        B.supplement_text(tmp_path / "any.pdf", structure)


def test_the_title_pages_counts_are_held_to_the_manuscript():
    line = "Word count: Main text 3,934 words; abstract 247 words; 3 figures; 4 tables; 42 references"
    counts = {"main text words": 3934, "abstract words": 247, "figures": 3, "tables": 4, "references": 42}
    got = B.counted_reasons(f"# Title page\n\n{line}\n", counts)
    assert [(doc, printed) for doc, _, printed, _, _ in got] == [("main", p) for p in ("3,934", "247", "3", "4", "42")]
    assert all(piece in line and reason.startswith(B.COUNTED) for _, piece, _, _, reason in got)
    with pytest.raises(SystemExit, match="prints 3,934 for the main text words, the manuscript has 3935"):
        B.counted_reasons(f"{line}\n", counts | {"main text words": 3935})
    with pytest.raises(SystemExit, match="not in the form expected"):
        B.counted_reasons("Word count: about four thousand\n", counts)


def test_a_number_without_an_entry_and_an_entry_without_a_number_stop_the_build():
    text = "# Abstract\n\nThe set holds 8,853 variants and 77 more.\n\n# References\n"
    bound = [("main", "The set holds", "8,853", 1, dict(file="phase1/reports/evidence/set_counts.csv", rows="", column="n",
                                                         derive="sum", scale=1))]
    with pytest.raises(SystemExit, match="77 is measured and has no binding"):
        B.rows_for("main", text, bound, [])
    rows, used = B.rows_for("main", text, bound, [("main", "The set holds", "77", 1, "a number made for this test")])
    assert [(r["printed"], r["kind"]) for r in rows] == [("8,853", "measured"), ("77", "constant")] and len(used) == 2
    with pytest.raises(SystemExit, match="has 2 entries"):
        B.rows_for("main", text, bound, [("main", "The set holds", "8,853", 1, "declared as well")])


def test_the_table_is_written_as_the_gate_reads_it():
    with open(V7, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    assert list(rows[0]) == gate.BINDING_COLUMNS and len(rows) == 239
    assert all(r["excerpt"] and r["printed"] in r["excerpt"] for r in rows)
