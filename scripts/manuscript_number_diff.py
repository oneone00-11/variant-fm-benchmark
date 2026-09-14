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


def by_block(path):
    out = OrderedDict()
    for t in manuscript_tokens(Path(path)):
        out.setdefault(t["uid"], []).append(t["token"])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("old"); ap.add_argument("new")
    ap.add_argument("--out", default="docs/manuscript-number-diff.md")
    a = ap.parse_args()
    old, new = by_block(a.old), by_block(a.new)
    rows, n_old, n_new = [], 0, 0
    for uid in list(old) + [u for u in new if u not in old]:
        o, n = old.get(uid, []), new.get(uid, [])
        n_old += len(o); n_new += len(n)
        sm = difflib.SequenceMatcher(a=o, b=n, autojunk=False)
        for tag, i1, i2, j1, j2 in sm.get_opcodes():
            if tag == "equal":
                continue
            if tag == "replace" and (i2 - i1) == (j2 - j1):
                for x, y in zip(o[i1:i2], n[j1:j2]):
                    rows.append((uid, x, y))
            else:
                for x in o[i1:i2]:
                    rows.append((uid, x, "(removed)"))
                for y in n[j1:j2]:
                    rows.append((uid, "(new)", y))
    changed = sum(1 for r in rows if r[1] not in ("(new)",) and r[2] != "(removed)")
    added = sum(1 for r in rows if r[1] == "(new)")
    removed = sum(1 for r in rows if r[2] == "(removed)")
    lines = [f"# Numeric tokens: `{Path(a.old).name}` → `{Path(a.new).name}`", "",
             f"Tokens: {n_old} before, {n_new} after. Changed in place: {changed}; new: {added}; removed: {removed}. "
             "Blocks are paragraph ids (P<n>) or table cells (T<table>r<row>c<col>) of the NEW file where the block "
             "exists there, of the old file otherwise.", "",
             "| block | old | new |", "|---|---|---|"]
    lines += [f"| {u} | {x} | {y} |" for u, x, y in rows]
    Path(a.out).write_text("\n".join(lines) + "\n")
    print(f"{n_old} -> {n_new} tokens; {changed} changed, {added} new, {removed} removed -> {a.out}")


if __name__ == "__main__":
    main()
