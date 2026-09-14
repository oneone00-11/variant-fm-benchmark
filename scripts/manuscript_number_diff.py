"""Every numeric token that changed between two versions of the manuscript.

    python scripts/manuscript_number_diff.py OLD.docx NEW.docx [--out docs/manuscript-number-diff.md]

Tokens are taken with the same tokeniser as the forward number check
(phase1/src/check_manuscript_numbers.manuscript_tokens), so the list is the
position / old / new record of what the reverse pass and the sentence rewrites
did to the numbers, independent of how the edits were made. Paragraphs are
aligned by id (P<n> for body paragraphs, T<t>r<r>c<c> for table cells) and the
tokens inside each block by difflib, so an inserted or deleted number shows up
as old-only or new-only rather than shifting everything after it.
"""
from __future__ import annotations

import argparse
import difflib
import sys
from collections import OrderedDict
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "phase1"))
from src.check_manuscript_numbers import manuscript_tokens  # noqa: E402


def blocks(path):
    """Ordered (uid, text, tokens) for body paragraphs, and {uid: tokens} for table cells."""
    import docx
    toks = OrderedDict()
    for t in manuscript_tokens(Path(path)):
        toks.setdefault(t["uid"], []).append(t["token"])
    d = docx.Document(str(path))
    paras = [(f"P{i}", p.text, toks.get(f"P{i}", [])) for i, p in enumerate(d.paragraphs)]
    cells = OrderedDict((u, v) for u, v in toks.items() if not u.startswith("P"))
    return paras, cells


def diff_tokens(old, new, rows):
    """old/new: lists of (uid, token). Appends (uid, old, new) for every change."""
    sm = difflib.SequenceMatcher(a=[t for _, t in old], b=[t for _, t in new], autojunk=False)
    for tag, i1, i2, j1, j2 in sm.get_opcodes():
        if tag == "equal":
            continue
        if tag == "replace" and (i2 - i1) == (j2 - j1):
            for (_, x), (u, y) in zip(old[i1:i2], new[j1:j2]):
                rows.append((u, x, y))
            continue
        for u, x in old[i1:i2]:
            rows.append((u, x, "(removed)"))
        for u, y in new[j1:j2]:
            rows.append((u, "(new)", y))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("old"); ap.add_argument("new")
    ap.add_argument("--out", default="docs/manuscript-number-diff.md")
    a = ap.parse_args()
    (po, co), (pn, cn) = blocks(a.old), blocks(a.new)
    rows = []
    # paragraphs are aligned by their text, not their index: a paragraph split or
    # inserted early in the document renumbers everything after it, and aligning by
    # index then reports every later number as removed and re-added
    sm = difflib.SequenceMatcher(a=[t for _, t, _ in po], b=[t for _, t, _ in pn], autojunk=False)
    for tag, i1, i2, j1, j2 in sm.get_opcodes():
        old = [(u, x) for u, _, ts in po[i1:i2] for x in ts]
        new = [(u, y) for u, _, ts in pn[j1:j2] for y in ts]
        if tag == "equal" and [x for _, x in old] == [y for _, y in new]:
            continue
        diff_tokens(old, new, rows)
    for uid in list(co) + [u for u in cn if u not in co]:
        diff_tokens([(uid, x) for x in co.get(uid, [])], [(uid, y) for y in cn.get(uid, [])], rows)
    n_old = sum(len(ts) for _, _, ts in po) + sum(map(len, co.values()))
    n_new = sum(len(ts) for _, _, ts in pn) + sum(map(len, cn.values()))
    changed = sum(1 for r in rows if r[1] != "(new)" and r[2] != "(removed)")
    added = sum(1 for r in rows if r[1] == "(new)")
    removed = sum(1 for r in rows if r[2] == "(removed)")
    lines = [f"# Numeric tokens: `{Path(a.old).name}` → `{Path(a.new).name}`", "",
             f"Tokens: {n_old} before, {n_new} after. Changed in place: {changed}; new: {added}; removed: {removed}. "
             "Paragraphs are aligned by text, so a split or inserted paragraph does not renumber the rest; the block "
             "column is the paragraph id (P<n>) or table cell (T<table>r<row>c<col>) in the NEW file, or in the old "
             "file for a removed number. Rewritten sentences show their numbers as removed and new.", "",
             "| block | old | new |", "|---|---|---|"]
    lines += [f"| {u} | {x} | {y} |" for u, x, y in rows]
    Path(a.out).write_text("\n".join(lines) + "\n")
    print(f"{n_old} -> {n_new} tokens; {changed} changed, {added} new, {removed} removed -> {a.out}")


if __name__ == "__main__":
    main()
