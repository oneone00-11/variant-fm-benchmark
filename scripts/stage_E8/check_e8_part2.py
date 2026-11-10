"""E8 part two (E8.10-E8.13): what each commit changed in the layout's outputs, cell by
cell, row by row and pixel by pixel, against the commit before it (--base), and whether
what the layout prints is the author's approved text word for word. It reads the
outputs, and for the numbering check draws the whole layout again into a temporary
folder with every item numbered in another scheme; it writes only its run record.
Part one's checks are scripts/stage_E8/check_e8.py (commits 1-3); this imports its
readers of a CSV and of the printed texts.

  --commit 4  E8.10. The approved texts (phase1/tests/fixtures/approved_texts_20261109.yaml,
              the author's file of 2026-11-09 byte for byte): every supplementary table's
              title, part titles and notes, every supplementary figure's title and legend,
              and the conventions, as the configuration prints them and where the stage
              prints them (contents.csv, both workbooks), character for character. Every
              CSV against the base cell by cell: only the seven reworded cells of the
              approved cell_fixes differ, each from the text it printed to the text the
              author gave. contents.csv field by field (titles, notes, the new legend
              column, the one-page estimate of the tables whose part notes became the
              table's); manifest.csv row by row (the files that differ, and one new row).
              The repository's workbook sheet by sheet: its title rows, notes and part
              titles are the approved texts, the part notes moved to their table gone from
              above their parts, every other row (each source's comment lines, each header
              and row) cell for cell as before (text and style: a number's format, bold),
              the seven cells aside; README's titles. The workbook for submission against
              the repository's: its contents sheet, and each sheet the same rows, cell for
              cell (the cell's XML: value, number format, italic runs), less the sources'
              comment lines and nothing else. Supplementary Figure S4 differs from the base
              only within its x-axis label's box; every other figure byte for byte. What is
              printed holds no "provisional", revision code, "v5" or "Table 3(b)"; and in
              the layout drawn with every item, the notes among them, renumbered (Q1, Q2,
              ...), no printed text names an item of this paper by its number.

  --commit 5  E8.11. The manuscript (read only, --manuscript): its paragraph of each of
              Tables 1-4 is "Table N. " and the approved title, a space and the approved
              notes, and each approved part title stands on a line of its own. The main
              tables' titles, part titles and notes and the main figures' titles, as the
              configuration prints them, are the approved texts. Against the base, only
              contents.csv (the main items' titles and notes, the tables' files) and
              manifest.csv (the four new files' rows and contents.csv's) change; every other
              output byte for byte. Each main table's file: "Table N." in bold and the
              title, the part titles, the notes, the approved texts; its tables the CSVs'
              cells, each cell; A4 at 2.5 cm, Tables 2 and 4 turned; ruled across at the
              top, below the header and at the foot, each rule at least 1 pt, no rule down
              and no shading; Arial; its members dated as the workbook's. The number gate's
              cells mode on each file (each data cell a paragraph, bound to the same cell of
              the CSV): every number matches. What is printed, the files' texts among it,
              holds no residue, and no number of this paper is left fixed when renumbered.

Run from the repository root:
    .venv/bin/python scripts/stage_E8/check_e8_part2.py --commit 4 [--base 1495bd9]
or in the reference image:
    ... IMAGE bash scripts/stage_E8/in_image.sh scripts/stage_E8/check_e8_part2.py --commit 4
(--commit 5 reads the manuscript, --manuscript PATH: in the image, a copy under
phase1/results/scratch/, which git ignores, removed after the run.)
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import os
import re
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog  # noqa: E402
from check_e8 import RESIDUES, SOURCE_PAPERS, _git_ls, _git_show, csv_cells, printed_texts  # noqa: E402

PHASE1 = REPO / "phase1"
OUT = Path("reports/evidence/layout")
APPROVED = PHASE1 / "tests" / "fixtures" / "approved_texts_20261109.yaml"
APPROVED_SHA256 = "654a1cc95d01ec84f872182ca289b8bd4af44e450e0864c7112f888a7c31b2c8"
# a number of this paper or of a source paper, as printed; the notes too (E8.10)
NUMBERED = re.compile(r"(?:Supplementary\s+)?(?:Table|Figure|Note|Fig\.)s?\s+S?\d+[a-z]?(?:\([a-z]\))?")
OURS = re.compile(r"(Supplementary )?(Table|Figure|Note)s? S?\d")
# a source paper's number the approved texts add: the second BRCA2 dataset's table
SOURCE_PAPERS_TOO = re.compile(r"Supplementary Table 5 of that publication")
# the parts whose notes the configuration moved to their table (E8.10: a table's notes are
# the table's), and how many each had
MOVED_PART_NOTES = {("label_sources", "effect"): 1, ("practical_metrics", "by_heldout_gene"): 1,
                    ("practical_metrics", "curve_areas"): 1, ("heldout_intervals", "intervals"): 1}
BOLD = 1                       # the workbook's style of a bold cell (evid_layout._styles_xml)


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _unxml(s: str) -> str:
    return s.replace("&lt;", "<").replace("&gt;", ">").replace("&quot;", '"').replace("&amp;", "&")


def _xml_text(inner: str) -> str:
    """A cell's text from what its XML holds: a string's runs joined, each italic run
    marked *GENE* again as the CSV marks it; a number as written."""
    runs = re.findall(r"<r><rPr>(.*?)</rPr><t[^>]*>(.*?)</t></r>", inner)
    if runs:
        return "".join(f"*{_unxml(t)}*" if "<i/>" in p else _unxml(t) for p, t in runs)
    v = re.search(r"<v>(.*?)</v>", inner)
    return v.group(1) if v else _unxml("".join(re.findall(r"<t[^>]*>(.*?)</t>", inner)))


def book(data: bytes) -> dict[str, list[dict]]:
    """A workbook's sheets, each its rows in order (an empty row {}), each row
    {column: (text, style)}, the style the cell's index (0 none, 1 bold, from 2 a number's
    format, evid_layout._styles_xml)."""
    z = zipfile.ZipFile(io.BytesIO(data))
    out = {}
    for name, i in re.findall(r'<sheet name="([^"]+)" sheetId="(\d+)"', z.read("xl/workbook.xml").decode()):
        rows: list[dict] = []
        for r, body in re.findall(r'<row r="(\d+)">(.*?)</row>', z.read(f"xl/worksheets/sheet{i}.xml").decode()):
            while len(rows) < int(r) - 1:
                rows.append({})
            row = {}
            for col, attrs, inner in re.findall(r'<c r="([A-Z]+)\d+"([^>]*)>(.*?)</c>', body):
                s = re.search(r's="(\d+)"', attrs)
                row[col] = (_xml_text(inner), int(s.group(1)) if s else 0)
            rows.append(row)
        out[_unxml(name)] = rows
    return out


def book_xml(data: bytes) -> dict[str, list[list[str]]]:
    """A workbook's sheets, each its rows, each its cells' XML without the row number in
    their references: what a cell holds, its style and its runs, so that two rows compare
    cell for cell wherever they stand."""
    z = zipfile.ZipFile(io.BytesIO(data))
    out = {}
    for name, i in re.findall(r'<sheet name="([^"]+)" sheetId="(\d+)"', z.read("xl/workbook.xml").decode()):
        xml = z.read(f"xl/worksheets/sheet{i}.xml").decode()
        out[_unxml(name)] = [re.findall(r'<c r="[A-Z]+\d+"(.*?)</c>', body)
                             for body in re.findall(r'<row r="\d+">(.*?)</row>', xml)]
    return out


def _bold(row: dict) -> bool:
    return bool(row) and all(s == BOLD for _, s in row.values())


def _texts(row: dict) -> list[str]:
    return [t for _, (t, _) in sorted(row.items(), key=lambda kv: (len(kv[0]), kv[0]))]


def cut(rows: list[dict], n_parts: int, lengths: list[int]) -> dict:
    """A table's sheet as the workbook writes it (evid_layout._sheet_rows): its title row,
    its notes, and for each part, after an empty row, its title (for a table of several
    parts), the plain one-cell rows above its header (its sources' comment lines and any
    notes of its own), its bold header and its rows. Fails loudly if the sheet is not so."""
    k, notes = 1, []
    while k < len(rows) and rows[k]:
        notes.append(rows[k])
        k += 1
    parts = []
    for j in range(n_parts):
        if k >= len(rows) or rows[k]:
            raise ValueError(f"no empty row before part {j} (row {k + 1})")
        k += 1
        title = None
        if n_parts > 1:
            title = rows[k]
            k += 1
        above = []
        while k < len(rows) and len(rows[k]) == 1 and not _bold(rows[k]):
            above.append(rows[k])
            k += 1
        if k >= len(rows) or not _bold(rows[k]):
            raise ValueError(f"part {j}: row {k + 1} is not a header")
        parts.append({"title": title, "above": above, "header": rows[k],
                      "body": rows[k + 1:k + 1 + lengths[j]]})
        k += 1 + lengths[j]
    if k != len(rows):
        raise ValueError(f"{len(rows) - k} rows after the last part")
    return {"title": rows[0], "notes": notes, "parts": parts}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--commit", type=int, required=True, choices=[4, 5])
    ap.add_argument("--base", default=None, help="the commit before (default: HEAD)")
    ap.add_argument("--manuscript", default=None, help="commit 5: the manuscript, read only")
    args = ap.parse_args()
    base = args.base or subprocess.run(["git", "-C", str(REPO), "rev-parse", "--short", "HEAD"],
                                       check=True, capture_output=True, text=True).stdout.strip()
    if args.manuscript:                         # from the folder it is run in, before the next line
        args.manuscript = str(Path(args.manuscript).resolve())
    os.chdir(PHASE1)
    sys.path.insert(0, str(PHASE1))
    import numpy as np
    import yaml
    from PIL import Image
    from src import evid_layout as L

    failed: list[str] = []

    def check(name: str, ok: bool, value=None) -> None:
        log.note(name, {"ok": bool(ok), "value": value} if value is not None else bool(ok))
        shown = "" if value is None else f": {value}"
        print(f"  {'ok  ' if ok else 'FAIL'} {name}{shown[:400]}")
        if not ok:
            failed.append(name)

    with RunLog(__file__) as log:
        log.note("base", base)
        log.note("commit", args.commit)
        cfg = L.load_config()
        old_cfg = yaml.safe_load(_git_show(base, "phase1/config/layout.yaml"))
        if args.commit == 4:
            commit4(cfg, old_cfg, base, check, log, np, Image, L, yaml)
        if args.commit == 5:
            commit5(cfg, old_cfg, base, check, log, L, yaml, args.manuscript)
        if failed:
            print(f"{len(failed)} check(s) failed: {failed}")
            raise SystemExit(1)
        print("0 check(s) failed")


# ------------------------------------------------------------------------- commit 4
def _notes_of(entry: dict) -> list[str]:
    n = entry.get("notes") or []
    return [n] if isinstance(n, str) else list(n)


def commit4(cfg, old_cfg, base, check, log, np, Image, L, yaml) -> None:
    A = yaml.safe_load(APPROVED.read_text())
    numbers = L._numbers(cfg)
    check("approved_texts_file_is_the_authors", _sha(APPROVED.read_bytes()) == APPROVED_SHA256)

    def printed(t: str) -> str:
        return L.grouped(L._resolve(t, numbers))

    items = {it["content"]: it for it in cfg["supplementary"]}
    tables = [it for it in cfg["supplementary"] if it["kind"] == "table"]

    def part_file(it, p) -> str:
        return "supplementary/" + L._part_path(Path("."), it, p).name

    # ------------------------------------------------------------- the configuration's texts
    print("E8.10 the configuration's texts against the approved ones")
    bad, compared = [], 0

    def same(where, got, want):
        nonlocal compared
        compared += 1
        if got != want:
            bad.append((where, got[:160], want[:160]))

    same("conventions", printed(cfg["supplementary_format"]["conventions"]), A["conventions"])
    for c, a in A["supplementary_figures"].items():
        same(f"{c} title", printed(items[c]["title"]), a["title"])
        same(f"{c} legend", printed(items[c]["legend"]), a["legend"])
        same(f"{c} notes", str(items[c].get("notes", [])), "[]")
    for c, a in A["supplementary_tables"].items():
        it = items[c]
        same(f"{c} title", printed(it["title"]), a["title"])
        same(f"{c} notes", str([printed(n) for n in it.get("notes", [])]), str(_notes_of(a)))
        want_parts = list((a.get("parts") or {}).values())
        got_parts = [printed(p["title"]) for p in it["parts"]] if want_parts else []
        same(f"{c} part titles", str(got_parts), str(want_parts))
        same(f"{c} part notes", str([p.get("notes") for p in it["parts"] if p.get("notes")]), "[]")
    check("every_item_of_the_approved_file_is_in_the_configuration",
          set(A["supplementary_tables"]) | set(A["supplementary_figures"]) == set(items))
    check("configuration_prints_the_approved_texts", not bad, bad[:6])
    log.note("approved_texts_compared_in_the_configuration", compared)
    print(f"    {compared} texts compared")

    # ------------------------------------------------------------- every CSV, cell by cell
    print("E8.10 every CSV against the base, cell by cell")
    fixes = A["cell_fixes"]
    # each fix as the base printed it and as the author gave it, by file
    rename = {
        "supplementary/label_sources.audit.csv": [("Table S1 cites no paper", printed(fixes[0]["print"]))],
        "supplementary/bp4_at_bp7_positions.csv": [
            (fixes[1]["printed_now"], printed(fixes[1]["print"])),
            (fixes[2]["printed_now"], fixes[2]["print"])],
        "supplementary/clinvar_coverage.csv": [(fixes[3]["printed_now"], fixes[3]["print"])],
        "supplementary/new_external_genes.counts_by_offset.csv":
            list(zip(fixes[4]["printed_now"].split("; "), fixes[4]["print"].split("; ")))
            + list(zip(fixes[5]["printed_now"].split("; "), fixes[5]["print"].split("; "))),
    }
    log.note("the_seven_fixes_as_checked", rename)

    def explained(rel: str, x: str, y: str) -> bool:
        for a_, b_ in rename.get(rel, []):
            if x == a_ and y == b_:
                return True
            if rel.endswith("label_sources.audit.csv") and a_ in x and x.replace(a_, b_) == y:
                return True                        # "...: not applicable" after the text
        return False

    base_files = {p.split("layout/", 1)[1] for p in _git_ls(base, f"phase1/{OUT}")}
    here_files = {str(p.relative_to(OUT)) for p in OUT.rglob("*") if p.is_file()}
    check("files_are_the_bases_and_the_workbook_for_submission",
          here_files - base_files == {"submission/supplementary_tables.xlsx"} and base_files <= here_files,
          sorted(base_files ^ here_files))
    changed_cells, unexplained, fixed = {}, [], {}
    for rel in sorted(f for f in base_files if f.endswith(".csv") and f not in ("contents.csv", "manifest.csv")):
        a = csv_cells(_git_show(base, f"phase1/{OUT}/{rel}"))
        b = csv_cells((OUT / rel).read_bytes())
        if [len(a[0]), len(a[1]), len(a[2])] != [len(b[0]), len(b[1]), len(b[2])]:
            unexplained.append((rel, "shape"))
            continue
        n = 0
        for x, y in zip(a[0] + a[1] + [c for r in a[2] for c in r], b[0] + b[1] + [c for r in b[2] for c in r]):
            if x != y:
                n += 1
                if not explained(rel, x, y):
                    unexplained.append((rel, x[:120], y[:120]))
                else:
                    fixed[(rel, x, y)] = fixed.get((rel, x, y), 0) + 1
        if n:
            changed_cells[rel] = n
    check("csv_changes_are_the_seven_fixes", not unexplained, unexplained[:8])
    check("csv_changes_where_the_fixes_are", set(changed_cells) == set(rename), sorted(changed_cells))
    check("label_audit_two_cells_and_conditions_24_each",
          changed_cells.get("supplementary/label_sources.audit.csv") == 2
          and sorted(v for (r, _, _), v in fixed.items() if r.endswith("bp4_at_bp7_positions.csv")) == [24, 24]
          and changed_cells.get("supplementary/clinvar_coverage.csv") == 1)
    log.note("csv_changed_cells", changed_cells)
    log.note("csv_fixes_made", [[r, x, y, v] for (r, x, y), v in sorted(fixed.items())])
    print("    changed cells:", changed_cells)
    off = csv_cells((OUT / "supplementary/new_external_genes.counts_by_offset.csv").read_bytes())
    bins = {r[2] for r in off[2] if r[1] == "Per position"}
    check("per_position_bins_as_they_were", all(re.fullmatch(r"\d+|>50", v) for v in bins), sorted(bins)[:5])

    # ------------------------------------------------------------- contents.csv
    print("E8.10 contents.csv, field by field")
    text_a = _git_show(base, f"phase1/{OUT}/contents.csv").decode()
    text_b = (OUT / "contents.csv").read_text()
    a, b = list(csv.DictReader(io.StringIO(text_a))), list(csv.DictReader(io.StringIO(text_b)))
    head_a, head_b = next(csv.reader(io.StringIO(text_a))), next(csv.reader(io.StringIO(text_b)))
    at = head_a.index("title") + 1
    check("contents_has_a_legend_column_after_the_title", head_b == head_a[:at] + ["legend"] + head_a[at:], head_b)
    key = lambda r: r["content"] or r["sources"]                              # noqa: E731
    ka, kb = {key(r): r for r in a}, {key(r): r for r in b}
    check("contents_same_rows_in_the_same_order", [key(r) for r in a] == [key(r) for r in b])
    bad, per_field = [], {}
    estimate = re.compile(r"^(yes|no) \(([\d,]+) lines, ([\d,]+) mm wide\)$")
    for k, rb in kb.items():
        ra, it = ka[k], items.get(k)
        for f in head_b:
            x, y = ra.get(f, ""), rb[f]
            if x != y:
                per_field[f] = per_field.get(f, 0) + 1
            if it is not None and f in ("title", "legend", "notes"):
                entry = (A["supplementary_figures"] if it["kind"] == "figure" else A["supplementary_tables"])[k]
                want = {"title": entry["title"], "legend": entry.get("legend", ""),
                        "notes": " | ".join(_notes_of(entry))}[f]
            elif it is not None and f == "fits_one_page_estimate" and x != y:
                moved = sum(v for (c, _), v in MOVED_PART_NOTES.items() if c == k)
                mx, my = estimate.match(x), estimate.match(y)
                ok = bool(mx and my and moved and mx.group(3) == my.group(3) and
                          int(mx.group(2).replace(",", "")) - moved == int(my.group(2).replace(",", "")))
                want = y if ok else f"{x}, less {moved} line(s)"
            else:
                want = x
            if want != y:
                bad.append((k, f, x[:90], y[:90]))
    check("contents_changes_explained", not bad, bad[:8])
    log.note("contents_changed_fields", per_field)
    print("    changed fields:", per_field)

    # ------------------------------------------------------------- manifest.csv
    print("E8.10 manifest.csv, row by row")
    a = list(csv.DictReader(io.StringIO(_git_show(base, f"phase1/{OUT}/manifest.csv").decode())))
    b = list(csv.DictReader(io.StringIO((OUT / "manifest.csv").read_text())))
    ka, kb = {r["file"]: r for r in a}, {r["file"]: r for r in b}
    new_row = "phase1/reports/evidence/layout/submission/supplementary_tables.xlsx"
    book_row = "phase1/reports/evidence/layout/supplementary_tables.xlsx"
    order = [r["file"] for r in b]
    check("manifest_rows_are_the_bases_and_one_new", set(kb) - set(ka) == {new_row} and set(ka) <= set(kb)
          and [f for f in order if f != new_row] == [r["file"] for r in a])
    check("manifest_new_row_follows_the_workbooks", order.index(new_row) == order.index(book_row) + 1)
    differs = {f for f in base_files if _git_show(base, f"phase1/{OUT}/{f}") != (OUT / f).read_bytes()}
    bad, per_field = [], {}
    for k, ra in ka.items():
        rb = kb[k]
        for f in ra:
            if ra[f] != rb[f]:
                per_field[f] = per_field.get(f, 0) + 1
                if not (f in ("bytes", "sha256") and k.split("layout/", 1)[1] in differs):
                    bad.append((k, f, ra[f][:60], rb[f][:60]))
    check("manifest_changes_are_sizes_and_hashes_of_the_files_that_differ", not bad, bad[:6])
    check("manifest_new_row_has_the_workbooks_sources",
          (kb[new_row]["sources"], kb[new_row]["sources_sha256"]) == (kb[book_row]["sources"], kb[book_row]["sources_sha256"]))
    check("manifest_hashes_true", all(_sha((REPO / r["file"]).read_bytes()) == r["sha256"]
                                      for r in b if not r["file"].endswith("manifest.csv")))
    log.note("manifest_changed_fields", per_field)
    log.note("files_that_differ_from_the_base", sorted(differs))
    print("    changed fields:", per_field, " files that differ:", len(differs))

    # ------------------------------------------------------------- the repository's workbook
    print("E8.10 the repository's workbook, sheet by sheet")
    old_numbers = L._numbers(old_cfg)
    old_items = {it["content"]: it for it in old_cfg["supplementary"]}
    sa = book(_git_show(base, f"phase1/{OUT}/supplementary_tables.xlsx"))
    sb = book((OUT / "supplementary_tables.xlsx").read_bytes())
    check("sheets_as_before", list(sb) == list(sa))
    bad, per_sheet, compared_book = [], {}, 0
    for it in tables:
        name = L._sheet_name(it["number"])
        lengths = [len(csv_cells((OUT / part_file(it, p)).read_bytes())[2]) for p in it["parts"]]
        try:
            ga, gb = cut(sa[name], len(it["parts"]), lengths), cut(sb[name], len(it["parts"]), lengths)
        except ValueError as e:
            bad.append((name, f"cannot be cut: {e}"))
            continue
        entry = A["supplementary_tables"][it["content"]]
        if _texts(gb["title"]) != [f"{it['number']}. {entry['title']}"] or not _bold(gb["title"]):
            bad.append((name, "title row", _texts(gb["title"])))
        if [_texts(r) for r in gb["notes"]] != [[x] for x in _notes_of(entry)]:
            bad.append((name, "notes", [_texts(r) for r in gb["notes"]][:2]))
        compared_book += 1 + len(_notes_of(entry))
        n = (ga["title"] != gb["title"]) + (ga["notes"] != gb["notes"])
        old_part_notes = {p["part"]: [L.grouped(L._resolve(x, old_numbers)) for x in p.get("notes", [])]
                          for p in old_items[it["content"]]["parts"]}
        for j, (pa, pb) in enumerate(zip(ga["parts"], gb["parts"])):
            part = it["parts"][j]
            if len(it["parts"]) > 1:
                want = [f"({chr(97 + j)}) {list(entry['parts'].values())[j]}"]
                if _texts(pb["title"]) != want or not _bold(pb["title"]):
                    bad.append((name, f"part {j} title", _texts(pb["title"])))
                compared_book += 1
                n += pa["title"] != pb["title"]
            moved = old_part_notes.get(part["part"], [])
            if len(moved) != MOVED_PART_NOTES.get((it["content"], part["part"]), 0):
                bad.append((name, f"part {j}: notes on the part in the base", moved))
            above_a, above_b = [_texts(r) for r in pa["above"]], [_texts(r) for r in pb["above"]]
            # the source's comment lines as before; the part's own notes, last above the
            # header in the base, gone
            if above_a != above_b + [[m] for m in moved]:
                bad.append((name, f"part {j}: rows above the header", above_a[-2:], above_b[-2:]))
            n += len(above_a) - len(above_b)
            rel = part_file(it, part)
            for ra_, rb_ in zip([pa["header"]] + pa["body"], [pb["header"]] + pb["body"]):
                if sorted(ra_) != sorted(rb_):
                    bad.append((name, rel, "row shape"))
                    continue
                for col in ra_:
                    x, y = ra_[col], rb_[col]
                    if x != y:
                        n += 1
                        if not (x[1] == y[1] and explained(rel, x[0], y[0])):
                            bad.append((name, rel, x, y))
        if n:
            per_sheet[f"{name} ({it['content']})"] = n
    # README: the title column alone changes, to the approved titles
    ra_, rb_ = sa["README"], sb["README"]
    check("readme_same_shape", [sorted(r) for r in ra_] == [sorted(r) for r in rb_])
    readme_bad, n_readme = [], 0
    for x, y in zip(ra_, rb_):
        for col in x:
            if x[col] != y.get(col):
                n_readme += 1
                c = x.get("C", ("", 0))[0]
                ok = (col == "D" and c in A["supplementary_tables"] and y[col][1] == x[col][1]
                      and y[col][0] == A["supplementary_tables"][c]["title"])
                if not ok:
                    readme_bad.append((x.get("A"), col, x[col][0][:60], str(y.get(col))[:60]))
    check("readme_changes_are_the_titles", not readme_bad and n_readme == len(tables), readme_bad[:6] or n_readme)
    per_sheet["README"] = n_readme
    check("workbook_changes_explained", not bad, bad[:8])
    log.note("workbook_changed_rows_or_cells", per_sheet)
    log.note("approved_texts_compared_in_the_workbook", compared_book)
    print("    changed rows or cells:", per_sheet)

    # ------------------------------------------------------------- the workbook for submission
    print("E8.10 the workbook for submission against the repository's")
    sub_x = book_xml((OUT / "submission" / "supplementary_tables.xlsx").read_bytes())
    book_x = book_xml((OUT / "supplementary_tables.xlsx").read_bytes())
    sub = book((OUT / "submission" / "supplementary_tables.xlsx").read_bytes())
    excel = [it for it in tables if it["destination"] == "excel"]
    pdf = [it for it in tables if it["destination"] == "pdf"]
    check("submission_sheets", list(sub_x) == ["Contents"] + [L._sheet_name(it["number"]) for it in excel],
          list(sub_x)[:8])
    check("submission_holds_29_tables_and_names_5_for_the_pdf", (len(excel), len(pdf)) == (29, 5))
    bad, dropped_total = [], 0
    for it in excel:
        name = L._sheet_name(it["number"])
        lengths = [len(csv_cells((OUT / part_file(it, p)).read_bytes())[2]) for p in it["parts"]]
        comments = [_texts(r) for p in cut(sb[name], len(it["parts"]), lengths)["parts"] for r in p["above"]]
        a_rows, b_rows = book_x[name], sub_x.get(name, [])
        gone, j = [], 0
        for row in a_rows:
            if j < len(b_rows) and row == b_rows[j]:
                j += 1
            else:
                gone.append(row)
        texts = [[_xml_text(re.sub(r"^[^>]*>", "", c, count=1)) for c in r] for r in gone]
        if j != len(b_rows) or texts != comments:
            bad.append((name, j, len(b_rows), texts[:2], comments[:2]))
        dropped_total += len(gone)
    check("submission_sheets_are_the_workbooks_less_the_comment_lines", not bad, bad[:5])
    log.note("comment_lines_left_out", dropped_total)
    sub_fmt = cfg["supplementary_format"]["submission_workbook"]
    want = ([[sub_fmt["heading"]], [], [A["conventions"]], []]
            + [[it["number"], printed(it["title"])] for it in excel]
            + [[], [f"Supplementary Tables S5, S6, S12, S25 and S28 {sub_fmt['pdf_tables']}"]])
    got = [_texts(r) for r in sub["Contents"]]
    check("submission_contents_sheet", got == want,
          [(k, x, y) for k, (x, y) in enumerate(zip(got, want)) if x != y][:3] or (len(got), len(want)))
    check("submission_contents_heading_bold_the_rest_plain",
          _bold(sub["Contents"][0]) and not any(_bold(r) for r in sub["Contents"][1:]))
    check("submission_titles_are_the_approved",
          [r[1] for r in got[4:4 + len(excel)]] == [A["supplementary_tables"][it["content"]]["title"] for it in excel])
    log.note("submission_contents_rows", len(got))

    # ------------------------------------------------------------- figures
    print("E8.10 figures: S4 only, within its x-axis label")
    drawn = sorted(f for f in base_files if f.endswith((".png", ".pdf")))
    moved = [f for f in drawn if _git_show(base, f"phase1/{OUT}/{f}") != (OUT / f).read_bytes()]
    check("only_s4_differs_byte_for_byte", moved == ["supplementary/clinvar_coverage_by_position.pdf",
                                                     "supplementary/clinvar_coverage_by_position.png"], moved)
    old = np.asarray(Image.open(io.BytesIO(_git_show(
        base, f"phase1/{OUT}/supplementary/clinvar_coverage_by_position.png"))).convert("RGBA"))
    new = np.asarray(Image.open(OUT / "supplementary/clinvar_coverage_by_position.png").convert("RGBA"))
    box = s4_label_box(L, cfg, np, old_label=fixes[6]["printed_now"])
    same_shape = old.shape == new.shape
    diff = np.argwhere((old != new).any(axis=2)) if same_shape else np.zeros((0, 2))
    inside = same_shape and len(diff) > 0 and bool(
        ((diff[:, 0] >= box[1]) & (diff[:, 0] <= box[3]) & (diff[:, 1] >= box[0]) & (diff[:, 1] <= box[2])).all())
    check("s4_differs_only_within_its_x_label", inside,
          {"box_px_x0_y0_x1_y1": box, "pixels_differing": int(len(diff)), "same_size": same_shape})
    log.note("s4_label_box_px", box)

    # ------------------------------------------------------------- residues, numbering
    print("E8.10 residues in what is printed")
    texts = printed_texts(OUT)
    for name, pat in RESIDUES.items():
        hits = [(w, m.group(0)) for w, t in texts for m in pat.finditer(t)]
        check(f"no_{name.replace(' ', '_').replace('(', '').replace(')', '')}", not hits, hits[:8])
    print("E8.10 every number of this paper printed through the configuration")
    with tempfile.TemporaryDirectory() as tmp:
        q_texts = q_render(Path(tmp), L)
    left = sorted({(w.split(" ")[0].split("!")[0], m.group(0), t) for w, t in q_texts
                   for m in NUMBERED.finditer(t)})
    ours = [x[:2] + (x[2][:140],) for x in left
            if OURS.match(x[1]) and not SOURCE_PAPERS.search(x[2]) and not SOURCE_PAPERS_TOO.search(x[2])]
    check("no_number_of_this_paper_left_fixed", not ours, ours[:10])
    log.note("numbers_left_when_renumbered", [[x[0], x[1], x[2][:200]] for x in left])
    print("    left in the renumbered layout (the source papers' numbers):")
    for x in left:
        print(f"      {x[0]}: «{x[1]}»  {x[2][:110]}")
    q_notes = sorted({m.group(0) for _, t in q_texts for m in re.finditer(r"Supplementary Notes? Q\d+", t)})
    check("the_notes_are_numbered_by_the_configuration", len(q_notes) >= 5, q_notes)


# ------------------------------------------------------------------------- commit 5
W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"


def docx_blocks(data: bytes):
    """A main table's file: its body's blocks in order, each ("p", text, its bold runs'
    texts) or ("tbl", rows of cell texts, the element), an italic run as *GENE*; and its
    section's properties (page size, orientation, margins)."""
    from xml.etree import ElementTree as ET
    root = ET.fromstring(zipfile.ZipFile(io.BytesIO(data)).read("word/document.xml"))

    def runs(el):
        out, bold = [], []
        for r in el.iter(f"{W}r"):
            t = "".join(x.text or "" for x in r.iter(f"{W}t"))
            rpr = r.find(f"{W}rPr")
            out.append(f"*{t}*" if rpr is not None and rpr.find(f"{W}i") is not None else t)
            if rpr is not None and rpr.find(f"{W}b") is not None:
                bold.append(t)
        return "".join(out), bold

    body = root.find(f"{W}body")
    blocks = []
    for el in body:
        if el.tag == f"{W}p":
            blocks.append(("p",) + runs(el))
        elif el.tag == f"{W}tbl":
            blocks.append(("tbl", [[runs(c)[0] for c in tr.findall(f"{W}tc")] for tr in el.findall(f"{W}tr")], el))
    sect = body.find(f"{W}sectPr")
    size, mar = sect.find(f"{W}pgSz"), sect.find(f"{W}pgMar")
    page = {"w": int(size.get(f"{W}w")), "h": int(size.get(f"{W}h")), "orient": size.get(f"{W}orient"),
            "margins": sorted({int(mar.get(f"{W}{k}")) for k in ("top", "right", "bottom", "left")})}
    return blocks, page


def docx_texts(root: Path) -> list[tuple[str, str]]:
    """Every text the main tables' files print under `root`: (where, text) for each
    paragraph and cell."""
    out = []
    for p in sorted(root.rglob("*.docx")):
        blocks, _ = docx_blocks(p.read_bytes())
        rel = str(p.relative_to(root))
        for k, b in enumerate(blocks):
            if b[0] == "p":
                out.append((f"{rel} paragraph {k}", b[1]))
            else:
                out += [(f"{rel} table {k} row {i} col {j}", c) for i, r in enumerate(b[1]) for j, c in enumerate(r)]
    return out


def cells_mode(tables: dict[str, list[list[str]]], csv_root: Path, tmp: Path, name: str) -> dict:
    """The number gate's cells mode (scripts/check_evidence_manuscript_numbers.py) on a
    printed table: each data cell written as a paragraph of its own, and each number in
    it bound to the same cell of the layout's CSV (the row by its first columns, enough
    to tell the rows apart; the column by its header): a number cell's value, or a text
    cell that holds the number. Header rows are not read. Returns the gate's counts and
    its last lines; writes only into `tmp`."""
    sys.path.insert(0, str(REPO / "scripts"))
    import check_evidence_manuscript_numbers as G
    lines, where = [], {}
    for rel, rows in tables.items():
        lay = list(csv.reader(io.StringIO((csv_root / rel).read_text())))
        nkey = 1
        while len({tuple(r[:nkey]) for r in lay[1:]}) < len(lay) - 1:
            nkey += 1
        lines += [f"# {name} {rel}", ""]
        for i, row in enumerate(rows[1:], 1):
            for j, cell in enumerate(row):
                where[len(lines) + 1] = (rel, lay[0], lay[i], lay[1:], j, nkey)
                lines += [cell.strip(), ""]
    text = "\n".join(lines) + "\n"
    md = tmp / f"{name}.md"
    md.write_text(text)
    out = []
    for t in G.read_tokens(text):
        rel, head, row, rows, j, nkey = where[t.line]
        if any(row[k] != row[k].strip() for k in range(nkey)):
            # a label the condition reader would strip: the row as the one no other is
            conds = f"{head[0]} not in {'|'.join(sorted({r[0] for r in rows if r[0] != row[0]}))}"
        else:
            conds = "; ".join(f"{head[k]}={row[k]}" for k in range(nkey))
        measured = not t.reason
        out.append({"document": "main", "line": t.line, "occurrence": t.occurrence, "printed": t.printed,
                    "excerpt": t.excerpt_space, "kind": "measured" if measured else "constant",
                    "reason": t.reason or "", "file": f"phase1/reports/evidence/layout/{rel}" if measured else "",
                    "rows": conds if measured else "", "column": head[j] if measured else "",
                    "derive": ("value" if G._decimal(row[j]) is not None else "text") if measured else "",
                    "scale": 1 if measured else ""})
    tsv = tmp / f"{name}_bindings.tsv"
    G.write_bindings(tsv, out)
    res = G.cells_results(text, G.read_bindings(tsv))
    n = {k: sum(1 for r in res if r.status == k)
         for k in ("constant", "match", "mismatch", "cannot_bind", "not_in_table")}
    bad = [(r.binding.printed if r.binding else r.token.printed, r.status, r.why) for r in res
           if r.status in ("mismatch", "cannot_bind", "not_in_table")]
    return {"numbers": sum(1 for r in res if r.token is not None), **n, "problems": bad[:5]}


def commit5(cfg, old_cfg, base, check, log, L, yaml, manuscript) -> None:
    A = yaml.safe_load(APPROVED.read_text())
    numbers = L._numbers(cfg)
    check("approved_texts_file_is_the_authors", _sha(APPROVED.read_bytes()) == APPROVED_SHA256)

    def printed(t: str) -> str:
        return L.grouped(L._resolve(t, numbers))

    main = {it["content"]: it for it in cfg["main"]}
    tables = [it for it in cfg["main"] if it["kind"] == "table"]

    # ------------------------------------------------------------- the manuscript, read only
    print("E8.11 the manuscript's paragraphs of Tables 1-4 against the approved texts")
    if not manuscript:
        check("manuscript_given", False, "--manuscript")
    else:
        import docx
        doc = docx.Document(manuscript)
        log.note("manuscript_sha256", _sha(Path(manuscript).read_bytes()))

        def with_italics(p) -> str:
            return re.sub(r"\*\*", "", "".join(f"*{r.text}*" if r.italic and r.text.strip() else r.text
                                               for r in p.runs))
        paras = doc.paragraphs
        found = {}
        for c, a in A["main_tables"].items():
            want = f"{a['number']}. {a['title']} {a['notes']}"
            hits = [p for p in paras if p.text.startswith(a["number"] + ".")]
            lines = [f"({chr(97 + k)}) {t}" for k, t in enumerate(a["parts"].values()) if t]
            found[c] = {"paragraphs": len(hits),
                        "equal": len(hits) == 1 and with_italics(hits[0]) == want,
                        "bold": len(hits) == 1 and [r.text.strip() for r in hits[0].runs if r.bold] == [f"{a['number']}."],
                        "part_lines": [sum(1 for p in paras if p.text.strip() == x) for x in lines]}
        check("manuscript_table_paragraphs_are_the_approved_texts",
              all(v["equal"] and v["bold"] and all(n == 1 for n in v["part_lines"]) for v in found.values()), found)

    # ------------------------------------------------------------- the configuration's texts
    print("E8.11 the configuration's texts against the approved ones")
    bad = []
    for c, a in A["main_tables"].items():
        it = main[c]
        titles = [t for t in a["parts"].values() if t]
        got = (it["number"], printed(it["title"]), [printed(n) for n in it.get("notes", [])],
               [printed(p["title"]) for p in it["parts"] if "title" in p])
        if got != (a["number"], a["title"], [a["notes"]], titles):
            bad.append((c, got))
    for c, a in A["main_figures"].items():
        if (main[c]["number"], printed(main[c]["title"])) != (a["number"], a["title"]) or "legend" in main[c]:
            bad.append((c, main[c]["title"][:80]))
    check("configuration_prints_the_approved_main_texts", not bad, bad[:4])

    # ------------------------------------------------------------- files against the base
    print("E8.11 every output against the base")
    new_files = {f"main/{it['content']}.docx" for it in tables}
    base_files = {p.split("layout/", 1)[1] for p in _git_ls(base, f"phase1/{OUT}")}
    here_files = {str(p.relative_to(OUT)) for p in OUT.rglob("*") if p.is_file()}
    check("files_are_the_bases_and_the_four_main_table_files",
          here_files - base_files == new_files and base_files <= here_files, sorted(base_files ^ here_files))
    differs = sorted(f for f in base_files if _git_show(base, f"phase1/{OUT}/{f}") != (OUT / f).read_bytes())
    check("only_contents_and_manifest_differ", differs == ["contents.csv", "manifest.csv"], differs)

    # contents.csv: the main items' titles and notes, the tables' files
    text_a = _git_show(base, f"phase1/{OUT}/contents.csv").decode()
    a_rows = list(csv.DictReader(io.StringIO(text_a)))
    b_rows = list(csv.DictReader(io.StringIO((OUT / "contents.csv").read_text())))
    check("contents_same_rows_and_columns", [(r["content"], list(r)) for r in a_rows]
          == [(r["content"], list(r)) for r in b_rows])
    bad, per_field = [], {}
    for ra, rb in zip(a_rows, b_rows):
        it = main.get(ra["content"])
        for f in ra:
            x, y = ra[f], rb[f]
            if x != y:
                per_field[f] = per_field.get(f, 0) + 1
            want = x
            if it is not None and f == "title":
                want = (A["main_tables"] if it["kind"] == "table" else A["main_figures"])[it["content"]]["title"]
            elif it is not None and it["kind"] == "table" and f == "notes":
                want = A["main_tables"][it["content"]]["notes"]
            elif it is not None and it["kind"] == "table" and f == "files":
                want = f"{x}; phase1/reports/evidence/layout/main/{it['content']}.docx"
            if want != y:
                bad.append((ra["content"], f, x[:80], y[:80]))
    check("contents_changes_explained", not bad, bad[:6])
    log.note("contents_changed_fields", per_field)
    print("    contents.csv changed fields:", per_field)

    # manifest.csv: a row for each new file, right after its table's CSVs, and contents.csv's
    a_man = list(csv.DictReader(io.StringIO(_git_show(base, f"phase1/{OUT}/manifest.csv").decode())))
    b_man = list(csv.DictReader(io.StringIO((OUT / "manifest.csv").read_text())))
    want_order = []
    for r in a_man:
        want_order.append(r["file"])
        for it in tables:
            last = [x["file"] for x in a_man if x["content"] == it["content"]][-1]
            if r["file"] == last:
                want_order.append(f"phase1/reports/evidence/layout/main/{it['content']}.docx")
    check("manifest_rows_the_bases_and_the_four_files_after_their_tables", [r["file"] for r in b_man] == want_order)
    ka = {r["file"]: r for r in a_man}
    bad = [(r["file"], f) for r in b_man if r["file"] in ka for f in r
           if r[f] != ka[r["file"]][f] and not (r["file"].endswith("contents.csv") and f in ("bytes", "sha256"))]
    new_rows = [r for r in b_man if r["file"] not in ka]
    for r in new_rows:
        it = main[r["content"]]
        srcs = [L._rel(s) for s in L._sources(it)]
        if (r["number"], r["sources"], r["sources_sha256"].split("; ")) != (
                it["number"], "; ".join(srcs), [_sha((PHASE1 / s).read_bytes()) for s in L._sources(it)]):
            bad.append((r["file"], "sources"))
    check("manifest_changes_explained", not bad, bad[:6])
    check("manifest_hashes_true", all(_sha((REPO / r["file"]).read_bytes()) == r["sha256"]
                                      for r in b_man if not r["file"].endswith("manifest.csv")))

    # ------------------------------------------------------------- the four files
    print("E8.11 the main tables' files")
    turned = {"evidence_tiers_and_thresholds", "practical_performance"}
    margin = round(cfg["docx"]["margin_cm"] * 10 * 1440 / 25.4)
    bad, sizes, gate = [], {}, {}
    with tempfile.TemporaryDirectory() as tmp:
        for it in tables:
            a = A["main_tables"][it["content"]]
            data = (OUT / "main" / f"{it['content']}.docx").read_bytes()
            blocks, page = docx_blocks(data)
            parts = [t for t in a["parts"].values() if t]
            kinds = [b[0] for b in blocks]
            want_kinds = ["p"] + (["p", "tbl"] * len(parts) if parts else ["tbl"]) + ["p"]
            if kinds != want_kinds:
                bad.append((it["content"], "blocks", kinds))
                continue
            if blocks[0][1] != f"{it['number']}. {a['title']}" or blocks[0][2] != [f"{it['number']}."]:
                bad.append((it["content"], "title paragraph", blocks[0][1][:80]))
            if [b[1] for b in blocks[1:-1] if b[0] == "p"] != [f"({chr(97 + k)}) {t}" for k, t in enumerate(parts)]:
                bad.append((it["content"], "part titles"))
            if blocks[-1][1] != a["notes"] or blocks[-1][2]:
                bad.append((it["content"], "notes"))
            got_tables = [b for b in blocks if b[0] == "tbl"]
            csvs = {}
            for part, (_, rows, el) in zip(it["parts"], got_tables):
                rel = "main/" + L._part_path(Path("."), it, part).name
                if rows != L.written_rows(OUT / rel):
                    bad.append((it["content"], rel, "cells"))
                csvs[rel] = rows
                tb = el.find(f"{W}tblPr/{W}tblBorders")
                rules = {s: (tb.find(f"{W}{s}").get(f"{W}val"), int(tb.find(f"{W}{s}").get(f"{W}sz") or 0))
                         for s in ("top", "bottom", "left", "right", "insideH", "insideV")}
                head_rules = [[(x.tag.replace(W, ""), x.get(f"{W}val"), int(x.get(f"{W}sz") or 0))
                               for x in list(tc.find(f"{W}tcPr/{W}tcBorders")
                                             if tc.find(f"{W}tcPr/{W}tcBorders") is not None else [])]
                              for tc in el.findall(f"{W}tr")[0].findall(f"{W}tc")]
                body_rules = [tc.find(f"{W}tcPr/{W}tcBorders") for tr in el.findall(f"{W}tr")[1:]
                              for tc in tr.findall(f"{W}tc")]
                ok = (all(rules[s][0] == "single" and rules[s][1] >= 8 for s in ("top", "bottom"))
                      and all(rules[s][0] == "nil" for s in ("left", "right", "insideH", "insideV"))
                      and all(h == [("bottom", "single", h[0][2])] and h[0][2] >= 8 for h in head_rules)
                      and all(b is None for b in body_rules))
                if not ok:
                    bad.append((it["content"], rel, "rules", rules))
                sizes[rel] = sorted({int(s.get(f"{W}val")) for s in el.iter(f"{W}sz")})
            xml = zipfile.ZipFile(io.BytesIO(data)).read("word/document.xml").decode()
            a4 = (16838, 11906) if it["content"] in turned else (11906, 16838)
            if ((page["w"], page["h"]) != a4 or (page["orient"] == "landscape") != (it["content"] in turned)
                    or page["margins"] != [margin]):
                bad.append((it["content"], "page", page))
            if "w:shd" in xml or set(re.findall(r'w:rFonts w:ascii="([^"]+)"', xml)) != {"Arial"}:
                bad.append((it["content"], "shading or font"))
            if {i.date_time for i in zipfile.ZipFile(io.BytesIO(data)).infolist()} != {(1980, 1, 1, 0, 0, 0)}:
                bad.append((it["content"], "member dates"))
            gate[it["number"]] = cells_mode(csvs, OUT, Path(tmp), it["number"].replace(" ", ""))
    check("main_table_files_as_approved_and_ruled_across_only", not bad, bad[:6])
    check("main_table_files_table_type_8pt", all(v == [16] for v in sizes.values()), sizes)
    check("cells_mode_every_number_matches",
          all(g["mismatch"] == g["cannot_bind"] == g["not_in_table"] == 0 and g["match"] > 0 for g in gate.values()),
          gate)
    log.note("cells_mode_on_the_files", gate)
    for k, g in gate.items():
        print(f"    cells mode, {k}: {g['numbers']} numbers, {g['constant']} constants, "
              f"{g['match']} match, {g['mismatch']} mismatch, {g['cannot_bind']} cannot bind, "
              f"{g['not_in_table']} not in the table")

    # ------------------------------------------------------------- residues, numbering
    print("E8.11 residues in what is printed, the files' texts among it")
    texts = printed_texts(OUT) + docx_texts(OUT)
    for name, pat in RESIDUES.items():
        hits = [(w, m.group(0)) for w, t in texts for m in pat.finditer(t)]
        check(f"no_{name.replace(' ', '_').replace('(', '').replace(')', '')}", not hits, hits[:8])
    with tempfile.TemporaryDirectory() as tmp:
        q_texts = q_render(Path(tmp), L) + docx_texts(Path(tmp) / "layout")
    left = sorted({(w.split(" ")[0].split("!")[0], m.group(0), t) for w, t in q_texts
                   for m in NUMBERED.finditer(t)})
    ours = [x[:2] + (x[2][:140],) for x in left
            if OURS.match(x[1]) and not SOURCE_PAPERS.search(x[2]) and not SOURCE_PAPERS_TOO.search(x[2])]
    check("no_number_of_this_paper_left_fixed", not ours, ours[:10])
    q_docx = [t for w, t in q_texts if w.endswith(".docx paragraph 0") or ".docx paragraph" in w]
    check("the_files_name_their_tables_and_others_by_the_configuration",
          any(t.startswith("Table Q1.") for t in q_docx) and any("Supplementary Tables Q8 and Q9" in t for t in q_docx)
          and any("Supplementary Tables Q10–Q12" in t for t in q_docx),
          [t[:60] for t in q_docx][:6])


def s4_label_box(L, cfg, np, old_label: str) -> list[int]:
    """Supplementary Figure S4's x-axis label in the PNG's pixels (x0, y0, x1, y1 from the
    top left): the union of its box as drawn now and as it was, two pixels round."""
    kept, save = [], L._save
    L._save = lambda fig, folder, stem, cfg_: kept.append(fig) or []
    try:
        item = next(i for i in cfg["supplementary"] if i["content"] == "clinvar_coverage_by_position")
        L.figure_clinvar_coverage(item, Path(tempfile.mkdtemp()), cfg)
    finally:
        L._save = save
    fig = kept[0]
    dpi = cfg["artwork"]["dpi"]
    fig.set_dpi(dpi)
    L._grouped_ticks(fig)
    boxes = []
    for text in (None, old_label):
        if text is not None:
            fig.axes[1].set_xlabel(text)
        fig.canvas.draw()
        boxes.append(fig.axes[1].xaxis.label.get_window_extent(fig.canvas.get_renderer()))
    h = fig.get_figheight() * dpi
    return [int(np.floor(min(b.x0 for b in boxes) - 2)), int(np.floor(h - max(b.y1 for b in boxes) - 2)),
            int(np.ceil(max(b.x1 for b in boxes) + 2)), int(np.ceil(h - min(b.y0 for b in boxes) + 2))]


def q_render(tmp: Path, L) -> list:
    """The whole layout drawn again with each item's number in another scheme (Table Q2,
    Supplementary Table Q24, Supplementary Note Q4, ...), into `tmp`: its printed texts."""
    orig_load = L.load_config

    def load(path=L.CONFIG):
        cfg = orig_load(path)
        for it in cfg["main"] + cfg["supplementary"] + cfg.get("supplementary_notes", []):
            it["number"] = re.sub(r"(Table|Figure|Note) S?(\d+)", r"\1 Q\2", it["number"])
        return cfg

    L.load_config = load
    try:
        L.main(out=tmp / "layout")
    finally:
        L.load_config = orig_load
    return printed_texts(tmp / "layout")


if __name__ == "__main__":
    main()
