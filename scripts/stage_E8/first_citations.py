"""E8.1: the order in which the manuscript first cites each supplementary table, figure
and note, recomputed from the manuscript itself (read only) and set beside the numbering
config/layout.yaml now gives (E8.1's table of new and old numbers).

The rule (E8.1): the text first, a figure's legend counted where it stands in the text;
then the tables' titles and notes; a list of several ("Supplementary Tables S10, S11 and
S15") in the order it is written, a range in its order; what only the supplementary
notes cite follows, in the notes' order. The manuscript still carries the numbers it was
written with, the layout's before E8 (--base), so each citation is read as an item by
that numbering, and the item's place in the order is its new number.

Only numbers and paragraph indices are printed and recorded, never the manuscript's text.

    .venv/bin/python scripts/stage_E8/first_citations.py ~/Desktop/calibration/CALInewversion.docx --base d6f7bc5
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog  # noqa: E402

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
KEY = re.compile(r"Supplementary (Table|Figure|Note)s?\s+")
ITEM = re.compile(r"S?(\d+)(?:\([a-z](?:[,–-] ?[a-z])*\))?")
SEP = re.compile(r"(,\s*and\s+|\s+and\s+|,\s*|\s*[–-]\s*|\s+to\s+)")
HEADINGS = {"Introduction": "text", "Materials and methods": "text", "Results": "text",
            "Discussion": "text", "References": "references", "Tables": "tables"}


def blocks(docx: Path) -> list[tuple[str, str, list[str]]]:
    """The body in order: (kind, style, texts), a paragraph's text or a table's cells."""
    root = ET.fromstring(zipfile.ZipFile(docx).read("word/document.xml"))

    def text(p) -> str:
        return "".join(t.text or "" for t in p.iter(W + "t"))

    out = []
    for el in root.find(W + "body"):
        if el.tag == W + "p":
            st = el.find(f"{W}pPr/{W}pStyle")
            out.append(("p", st.get(W + "val") if st is not None else "", [text(el)]))
        elif el.tag == W + "tbl":
            out.append(("tbl", "", [" ".join(text(p) for p in tc.iter(W + "p"))
                                    for tc in el.iter(W + "tc")]))
    return out


def cites(text: str) -> list[tuple[str, int]]:
    """(kind, number) in writing order; a range expands in its order."""
    out = []
    for m in KEY.finditer(text):
        kind, i, prev, rng = m.group(1), m.end(), None, False
        while True:
            im = ITEM.match(text, i)
            if not im:
                break
            n = int(im.group(1))
            out += ([(kind, k) for k in range(prev + 1, n + 1)] if rng and prev is not None
                    else [(kind, n)])
            prev, i = n, im.end()
            sm = SEP.match(text, i)
            if not sm or not ITEM.match(text, sm.end()) or (kind != "Note" and text[sm.end()] != "S"):
                break
            rng = sm.group(1).strip() in ("–", "-", "to")
            i = sm.end()
    return out


def sections(bl) -> list[str]:
    """Each block's place: text, references, tables, or a supplementary note's number."""
    out, cur = [], "front"
    for kind, style, texts in bl:
        t = texts[0].strip() if kind == "p" else ""
        if kind == "p" and t in HEADINGS:
            cur = HEADINGS[t]
        elif kind == "p" and t.startswith("Supplementary Notes"):
            cur = "notes"
        m = re.match(r"Supplementary Note (\d+)\.", t) if kind == "p" else None
        if cur.startswith("note") and m:
            cur = f"note {m.group(1)}"
        out.append(cur)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("docx", type=Path)
    ap.add_argument("--base", default="d6f7bc5", help="the commit whose numbers the text carries")
    args = ap.parse_args()
    with RunLog(__file__) as log:
        log.input(args.docx)
        bl = blocks(args.docx)
        where = sections(bl)
        old = yaml.safe_load(subprocess.run(["git", "-C", str(REPO), "show",
                                             f"{args.base}:phase1/config/layout.yaml"],
                                            check=True, capture_output=True, text=True).stdout)
        new = yaml.safe_load((REPO / "phase1/config/layout.yaml").read_text())
        old_items = {it["number"]: it["content"] for it in old["supplementary"]}
        new_number = {it["content"]: it["number"] for it in new["supplementary"]}
        # the text first, then the tables, then the notes in the order of their new numbers
        notes_new = {1: 1, 2: 2, 4: 3, 3: 4, 6: 5, 5: 6}           # E8.1: new <- old, inverted
        rank = {"text": 0, "tables": 1}
        seen, first = set(), []
        order = sorted(range(len(bl)), key=lambda i: (
            rank.get(where[i], 2), notes_new.get(int(where[i].split()[1]), 0)
            if where[i].startswith("note ") else 0, i))
        for i in order:
            if where[i] in ("front", "references"):
                continue
            for t in bl[i][2]:
                for kind, n in cites(t):
                    if (kind, n) not in seen:
                        seen.add((kind, n))
                        first.append({"kind": kind, "old": n, "block": i, "place": where[i]})
        ok = True
        for kind, word in (("Table", "Supplementary Table"), ("Figure", "Supplementary Figure")):
            got = [f for f in first if f["kind"] == kind]
            print(f"{word}s, in the order first cited (old number -> place -> number now):")
            for k, f in enumerate(got, 1):
                content = old_items[f"{word} S{f['old']}"]
                now = new_number[content]
                good = now == f"{word} S{k}"
                ok &= good
                print(f"  {'ok  ' if good else 'FAIL'} S{f['old']:<2} -> S{k:<2} {content:32s} "
                      f"now {now}  (block {f['block']}, {f['place']})")
            ok &= len(got) == sum(1 for it in new["supplementary"] if it["kind"] == kind.lower())
            log.note(f"{kind.lower()}s", [{"old": f["old"], "place": k, "block": f["block"],
                                           "where": f["place"]} for k, f in enumerate(got, 1)])
        got = [f for f in first if f["kind"] == "Note"]
        print("Supplementary Notes, in the order first cited (old number -> place):")
        for k, f in enumerate(got, 1):
            good = notes_new[f["old"]] == k
            ok &= good
            print(f"  {'ok  ' if good else 'FAIL'} {f['old']} -> {k}  (block {f['block']}, {f['place']})")
        log.note("notes", [{"old": f["old"], "place": k, "block": f["block"]} for k, f in enumerate(got, 1)])
        log.note("matches_layout", bool(ok))
        print("the first-citation order matches config/layout.yaml's numbering" if ok else
              "MISMATCH: the first-citation order differs from config/layout.yaml's numbering")
        if not ok:
            raise SystemExit(1)


if __name__ == "__main__":
    main()
