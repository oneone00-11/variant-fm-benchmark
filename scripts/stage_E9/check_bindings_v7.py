"""E9.3: the number gate (scripts/check_evidence_manuscript_numbers.py) on the text as
submitted, in each mode v6 was checked in (scripts/stage_E8/check_bindings_v6.py, which
stays as it is), the cells mode on the four main tables' files, and what changed between
the binding tables of v6 and v7.

Before the modes, three things are held, and the run stops at the first that fails:
- the texts in --texts are the ones build_bindings_v7.py writes from the files given
  (the manuscript, the four tables' files, the supplementary PDF), and
  docs/manuscript_bindings_v7.tsv is the table it writes for them: both are built again
  here, in memory, and compared byte for byte. The build holds the title page's counts to
  the manuscript and the notes' paragraphs to the PDF's pages, as it always does;
- each table file given is, byte for byte, one of the layout's four
  (phase1/reports/evidence/layout/main/<content>.docx), each once.

Then:
- pool mode, the weak check: on the manuscript as it stands (.docx, which that mode reads
  from its first line to the heading "References") and on the two texts;
- bound mode: each sentence's numbers against the files its source tag names (the keyed
  texts), and the self-test, which plants numbers the named files do not hold and passes
  only if the bound mode reports each;
- cells mode on the two texts with the binding table: every measured number against its
  one cell, every constant with its reason;
- cells mode on each of the four table files: each data cell a paragraph, bound to the
  same cell of the layout's CSV (scripts/stage_E8/check_e8_part2.cells_mode).

Recorded whole: each mode's verdict; a tally of every number written in the texts by what
the gate does with it (read, or left out: headings, identifiers, citations, labels, digits
within a word or a symbol), for each text and for each part of the main text; the numbers
the texts leave out whole (the two reference lists, the notes' list labels); and the
difference between the two binding tables, read from the tables alone:
- vanished: a number of v6 that v7 does not hold;
- new: a number of v7 that v6 does not hold;
- reworded: the same number with the same entry (its cell, or its reason), the words
  beside it changed, so that its excerpt differs;
- moved: a number that both hold, excerpt and entry alike, whose place among the others
  changed.

Exits 1 if the cells mode, on the texts or on a table, finds a number that does not match
its cell or that cannot be bound, or if the self-test fails.

    .venv/bin/python scripts/stage_E9/check_bindings_v7.py \\
        --manuscript ~/Desktop/submission/Manuscript.docx \\
        --tables ~/Desktop/submission/Table_1.docx ~/Desktop/submission/Table_2.docx \\
                 ~/Desktop/submission/Table_3.docx ~/Desktop/submission/Table_4.docx \\
        --supplement ~/Desktop/submission/Supplementary_Information.pdf \\
        --texts ~/Desktop/整理文件/E9_工作文件/v7_texts

(--structure as for build_bindings_v7.py, where poppler is not installed.)
"""
from __future__ import annotations

import argparse
import csv
import difflib
import hashlib
import sys
import tempfile
from collections import Counter
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from _runlog import REPO, RunLog  # noqa: E402
import build_bindings_v7 as B  # noqa: E402

sys.path.insert(0, str(REPO / "scripts" / "stage_E8"))
import check_bindings_v6 as C6  # noqa: E402  its tally of the numbers, its table parts
from check_e8_part2 import cells_mode, docx_blocks  # noqa: E402

sys.path.insert(0, str(REPO / "scripts"))
import check_evidence_manuscript_numbers as gate  # noqa: E402

TABLE, BEFORE = B.OUT, REPO / "docs" / "manuscript_bindings_v6.tsv"
LAYOUT = REPO / "phase1" / "reports" / "evidence" / "layout"
ENTRY = ("kind", "reason", "file", "rows", "column", "derive", "scale")
TEXTS = ("v7_main.md", "v7_supplement.md", "v7_main_keyed.md", "v7_supplement_keyed.md")


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def rebuilt(ms: Path, tables: list[Path], pdf: Path, structure: Path | None) -> tuple[dict[str, str], bytes]:
    """The four texts and the binding table as build_bindings_v7.py writes them for these
    files, built in memory."""
    main_md, counts = B.main_text(ms, tables)
    supp_md, _ = B.supplement_text(pdf, structure)
    reasons = B.V6.REASONS + B.NEW_REASONS + B.counted_reasons(main_md, counts)
    texts, rows = {"v7_main.md": main_md, "v7_supplement.md": supp_md}, []
    for doc, text in (("main", main_md), ("supplement", supp_md)):
        r, _ = B.rows_for(doc, text, B.V6.BINDINGS, reasons)
        rows += r
        texts[f"v7_{doc}_keyed.md"] = B.V6.keyed(text, r)
    with tempfile.TemporaryDirectory() as tmp:
        gate.write_bindings(Path(tmp) / "t.tsv", rows)
        return texts, (Path(tmp) / "t.tsv").read_bytes()


def parts_of(text: str) -> dict[str, str]:
    """A text's parts by their headings of the first level, up to where the gate stops."""
    start, end = gate.body_span(text)
    out: dict[str, list[str]] = {}
    here = None
    for line in text[start:end].split("\n"):
        if line.startswith("# "):
            here = line[2:]
            out[here] = []
        elif here is not None:
            out[here].append(line)
    return {k: "\n".join(v) for k, v in out.items()}


def table_rows(path: Path) -> list[dict]:
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def difference(before: list[dict], after: list[dict]) -> dict:
    """What changed between two binding tables (see the module's docstring)."""
    full = lambda r: (r["document"], r["printed"], r["excerpt"]) + tuple(r[k] for k in ENTRY)      # noqa: E731
    entry = lambda r: (r["document"], r["printed"]) + tuple(r[k] for k in ENTRY)                   # noqa: E731
    show = lambda r: {k: r[k] for k in ("document", "line", "printed", "excerpt", "kind")} | {     # noqa: E731
        "entry": r["reason"] if r["kind"] == "constant" else f"{r['file']} [{r['rows']}] {r['column']} ({r['derive']})"}

    def numbered(rows):                 # each row's key with its rank among the rows keyed alike
        seen: Counter = Counter()
        out = []
        for r in rows:
            seen[full(r)] += 1
            out.append((full(r), seen[full(r)]))
        return out

    a, b = numbered(before), numbered(after)
    both = set(a) & set(b)
    only_a = [r for r, k in zip(before, a) if k not in both]
    only_b = [r for r, k in zip(after, b) if k not in both]
    reworded, left_b = [], list(only_b)
    vanished = []
    for r in only_a:
        hit = next((x for x in left_b if entry(x) == entry(r)), None)
        if hit is None:
            vanished.append(show(r))
        else:
            left_b.remove(hit)
            reworded.append({"before": show(r), "after": show(hit)})
    kept_a, kept_b = [k for k in a if k in both], [k for k in b if k in both]
    blocks = difflib.SequenceMatcher(a=kept_a, b=kept_b, autojunk=False).get_matching_blocks()
    in_place = {kept_a[m.a + i] for m in blocks for i in range(m.size)}
    line_b = {k: r["line"] for r, k in zip(after, b)}
    moved = [show(r) | {"line_after": line_b[k]} for r, k in zip(before, a) if k in both and k not in in_place]
    return {"numbers_before": len(before), "numbers_after": len(after), "kept_in_place": len(in_place),
            "moved": moved, "reworded": reworded, "vanished": vanished, "new": [show(r) for r in left_b]}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manuscript", required=True, type=Path)
    ap.add_argument("--tables", required=True, type=Path, nargs=4)
    ap.add_argument("--supplement", required=True, type=Path)
    ap.add_argument("--structure", type=Path)
    ap.add_argument("--texts", required=True, type=Path)
    args = ap.parse_args()
    ms, pdf = args.manuscript.expanduser().resolve(), args.supplement.expanduser().resolve()
    tables = [p.expanduser().resolve() for p in args.tables]
    structure = args.structure.expanduser().resolve() if args.structure else None
    folder = args.texts.expanduser().resolve()
    failed = []
    with RunLog(__file__) as log:
        log.note("manuscript_sha256", sha(ms.read_bytes()))
        log.note("supplement_sha256", sha(pdf.read_bytes()))
        log.note("bindings_sha256", sha(TABLE.read_bytes()))
        for name in TEXTS:
            log.note(f"{name}_sha256", sha((folder / name).read_bytes()))

        # the texts and the table are the ones these files give
        texts, table = rebuilt(ms, tables, pdf, structure)
        differ = [n for n in TEXTS if texts[n].encode() != (folder / n).read_bytes()]
        if differ or table != TABLE.read_bytes():
            raise SystemExit(f"not the builder's output for the files given: {differ + ([TABLE.name] if table != TABLE.read_bytes() else [])}")
        log.note("texts_and_table_are_the_builders_for_these_files", True)
        print("== the four texts and the binding table are the builder's for the files given")

        # the table files are the layout's
        layout_files = {sha((LAYOUT / "main" / f"{c}.docx").read_bytes()): name for name, c in C6.FILES.items()}
        given = {}
        for p in tables:
            name = layout_files.get(sha(p.read_bytes()))
            if name is None or name in given:
                raise SystemExit(f"{p.name}: not one of the layout's four table files, or given twice")
            given[name] = p
            log.note(f"{p.name}_sha256", sha(p.read_bytes()))
        log.note("table_files", {p.name: f"layout/main/{C6.FILES[name]}.docx ({name}), byte for byte" for name, p in given.items()})
        print("== the four table files are the layout's, byte for byte: "
              + "; ".join(f"{p.name} = {C6.FILES[name]}.docx" for name, p in given.items()))

        # pool mode
        for label, target in (("pool_manuscript", ms), ("pool_main", folder / "v7_main.md"),
                              ("pool_supplement", folder / "v7_supplement.md")):
            rc, out = C6.run_gate(str(target), *(["--pool"] if target.suffix == ".md" else []))
            log.note(label, {"exit": rc, "output": out})
            print(f"== {label}: exit {rc}\n" + "\n".join(
                ln for ln in out.splitlines() if ln.startswith(("tokens:", "  UNMATCHED"))))
        # bound mode and the self-test
        for doc in ("main", "supplement"):
            rc, out = C6.run_gate(str(folder / f"v7_{doc}_keyed.md"))
            log.note(f"bound_{doc}", {"exit": rc, "output": out})
            print(f"== bound mode, {doc}: exit {rc}\n" + "\n".join(
                ln for ln in out.splitlines() if ln.startswith(("bound mode", "  NOT IN", "  UNBOUND"))))
        rc, out = C6.run_gate(str(folder / "v7_main_keyed.md"), "--selftest")
        log.note("selftest", {"exit": rc, "output": out})
        print(f"== self-test: exit {rc}\n" + "\n".join(out.strip().splitlines()[-2:]))
        if rc:
            failed.append("selftest")

        # cells mode on the texts, and the tally of the numbers written
        bindings = gate.read_bindings(TABLE)
        _, counts = B.main_text(ms, tables)
        _, notes = B.supplement_text(pdf, structure)
        whole = {"main references": counts["numbers in the reference list"],
                 "notes references": notes["numbers in the notes' reference list"],
                 "notes list numbers": notes["list labels left out"]}
        log.note("counted_in_the_manuscript", counts)
        log.note("excluded_whole", whole)
        print(f"== left out whole: {whole}")
        for doc in ("main", "supplement"):
            text = texts[f"v7_{doc}.md"]
            tally = C6.excluded_in_text(text)
            log.note(f"numbers_{doc}", tally)
            print(f"== numbers written in the {doc} text: {tally}")
            if doc == "main":
                by_part = {name: C6.excluded_in_text(part) for name, part in parts_of(text).items()}
                log.note("numbers_main_by_part", by_part)
                for name, t in by_part.items():
                    print(f"     {name}: " + ", ".join(f"{k} {v}" for k, v in t.items() if v))
            res = gate.cells_results(text, [b for b in bindings if b.document == doc])
            n = {k: sum(1 for r in res if r.status == k)
                 for k in ("constant", "match", "mismatch", "cannot_bind", "not_in_table")}
            bad = [{"line": (r.binding or r.token).line, "printed": (r.binding or r.token).printed,
                    "status": r.status, "now": r.now, "why": r.why,
                    "where": r.binding.where() if r.binding else "",
                    "excerpt": r.binding.excerpt if r.binding else r.token.excerpt_space[:120]}
                   for r in res if r.status in ("mismatch", "cannot_bind", "not_in_table")]
            log.note(f"cells_{doc}", n | {"numbers": sum(1 for r in res if r.token is not None), "problems": bad})
            print(f"== cells mode, {doc}: {n}")
            for x in bad:
                print(f"   {x['status'].upper()} line {x['line']}: {x['printed']} ({x['why'] or x['now']})  "
                      f"[{x['where']}]  ...{x['excerpt'][:80]}...")
            if bad:
                failed.append(f"cells_{doc}")

        # cells mode on each table file
        with tempfile.TemporaryDirectory() as tmp:
            for name, parts in C6.PARTS.items():
                p = given[name]
                found = [b[1] for b in docx_blocks(p.read_bytes())[0] if b[0] == "tbl"]
                if len(found) != len(parts):
                    failed.append(f"parts_{p.name}")
                    log.note(f"cells_{p.name}", {"parts_found": len(found), "parts": len(parts)})
                    continue
                g = cells_mode(dict(zip(parts, found)), LAYOUT, Path(tmp), p.stem)
                log.note(f"cells_{p.name}", g | {"parts": parts})
                print(f"== cells mode, {p.name} ({name}): {g['numbers']} numbers, {g['constant']} constants, "
                      f"{g['match']} match, {g['mismatch']} mismatch, {g['cannot_bind']} cannot bind, "
                      f"{g['not_in_table']} not in the table")
                if g["mismatch"] or g["cannot_bind"] or g["not_in_table"]:
                    failed.append(f"cells_{p.name}")

        # v6 to v7
        log.note("before_sha256", sha(BEFORE.read_bytes()))
        d = difference(table_rows(BEFORE), table_rows(TABLE))
        log.note("v6_to_v7", d)
        print(f"== v6 to v7: {d['numbers_before']} numbers before, {d['numbers_after']} after; "
              f"{d['kept_in_place']} kept in place, {len(d['moved'])} moved, {len(d['reworded'])} reworded, "
              f"{len(d['vanished'])} vanished, {len(d['new'])} new")
        for key in ("vanished", "new"):
            for x in d[key]:
                print(f"   {key.upper()} {x['document']} line {x['line']}: {x['printed']} ({x['kind']}: {x['entry'][:70]})  ...{x['excerpt']}...")
        for x in d["reworded"]:
            print(f"   REWORDED {x['before']['printed']}: ...{x['before']['excerpt']}...  ->  ...{x['after']['excerpt']}...")
        for x in d["moved"]:
            print(f"   MOVED {x['document']} line {x['line']} -> {x['line_after']}: {x['printed']}  ...{x['excerpt']}...")

        log.note("failed", failed)
        print(f"{len(failed)} check(s) failed: {failed}" if failed else "0 check(s) failed")
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
