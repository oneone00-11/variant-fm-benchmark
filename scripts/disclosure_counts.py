"""The counts Methods 2.6 quotes, written from the record they were measured into.

Methods 2.6 states two kinds of count that no report table holds: how many tests pass
from an archive extract of the release, with both documents present and in a clone,
and how many measured numbers the number check binds. Both are measured and recorded in
phase1/config/pipeline_facts.json: the test counts from the runs, the binding counts by
`python -m src.check_manuscript_numbers <docx> --record-facts`. The test count went stale
once while the sentence was maintained by hand. This renders every such count from the
record and, with --apply, writes it into the manuscript in place, so no count in the
sentence is typed. tests/test_reported_numbers.py compares the sentence with the record,
and the binding counts with a fresh measurement.

    python scripts/disclosure_counts.py <manuscript.docx>           # list differences; exit 1 if any
    python scripts/disclosure_counts.py <manuscript.docx> --apply   # write the recorded counts
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
FACTS = REPO / "phase1" / "config" / "pipeline_facts.json"

UNITS = ["zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten", "eleven",
         "twelve", "thirteen", "fourteen", "fifteen", "sixteen", "seventeen", "eighteen", "nineteen"]
TENS = {2: "twenty", 3: "thirty", 4: "forty", 5: "fifty", 6: "sixty", 7: "seventy", 8: "eighty", 9: "ninety"}
WORD = r"\b([A-Za-z]+(?:-[a-z]+)?)"
NUM = r"(\d[\d,]*)"


def words(n: int) -> str:
    if n < 20:
        return UNITS[n]
    if n < 100:
        t, u = divmod(n, 10)
        return TENS[t] + (f"-{UNITS[u]}" if u else "")
    raise ValueError(f"{n}: the sentence spells counts out, and this table stops at ninety-nine")


def expected(facts: dict) -> list[tuple[str, str]]:
    """(pattern whose one group is the count, the text that group should read)."""
    a = facts["archive_breakdown"]
    docs = a["with_manuscript_visible"]
    b = facts["manuscript_binding"]
    return [
        (WORD + r" automated tests pass from an extract", words(a["passed"]).capitalize()),
        (WORD + r" are collected: [a-z-]+ read the manuscript", words(a["collected"]).capitalize()),
        # the skips that disappear once both documents are visible are the tests that read them
        (r"are collected: " + WORD + r" read the manuscript", words(a["skipped"] - docs["skipped"])),
        (r"With both documents present " + WORD + r" pass", words(docs["passed"])),
        (r"in a clone of the same commit all " + WORD, words(a["in_a_clone_with_manuscript_visible"]["passed"])),
        (r"All " + NUM + r" measured numbers are bound", f"{b['measured']:,}"),
        (r"are bound: " + NUM + r" to a single cell", f"{b['single_cell']:,}"),
        (r"to a single cell and " + NUM + r" to one of several", f"{b['several_cells']:,}"),
        (r"The other " + NUM + r" numeric tokens", f"{b['non_measurement']:,}"),
    ]


def _replace_span(p, start: int, end: int, new: str) -> None:
    """Replace characters [start, end) of a paragraph inside the run that holds them."""
    if "".join(r.text for r in p.runs) != p.text:
        raise ValueError("paragraph text is not the concatenation of its runs; offsets would be wrong")
    pos = 0
    for r in p.runs:
        if pos <= start and end <= pos + len(r.text):
            r.text = r.text[:start - pos] + new + r.text[end - pos:]
            return
        pos += len(r.text)
    raise ValueError("the count spans more than one run")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("docx")
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args()
    import docx

    d = docx.Document(a.docx)
    changes = []
    for pattern, want in expected(json.loads(FACTS.read_text())):
        rx = re.compile(pattern)
        hits = [p for p in d.paragraphs if rx.search(p.text)]
        if len(hits) != 1:
            print(f"{len(hits)} paragraphs match {pattern!r}; expected exactly one")
            return 2
        m = rx.search(hits[0].text)
        if m.group(1) != want:
            changes.append((m.group(1), want, pattern))
            if a.apply:
                _replace_span(hits[0], m.start(1), m.end(1), want)
    for got, want, pattern in changes:
        print(f"  {got!r} -> {want!r}   ({pattern})")
    if a.apply and changes:
        d.save(a.docx)
    print(f"{len(changes)} count(s) {'written' if a.apply else 'differ from the record'}")
    return 0 if a.apply or not changes else 1


if __name__ == "__main__":
    sys.exit(main())
