"""E8.6: what each E8 commit changed in the layout's outputs, cell by cell, row by row and
pixel by pixel, against the commit before it (--base), and the residues left in what the
layout prints. It reads the outputs, and for the numbering check draws the whole layout
again into a temporary folder with every item numbered in another scheme; it writes only
its run record.

  --commit 2  E8.3-E8.4. Every CSV, contents.csv, manifest.csv and every sheet of the
              workbook against commit 1, each change explained: a number of 1,000 or more
              grouped in thousands (evid_layout.grouped), DDX3X's and SBDS's transcripts
              as the MANE summary's RefSeq with the deposit's Ensembl transcript after it,
              the two tables' new note (their sheets one row longer), the MANE summary
              among their sources, a gene's italic marks gone into italic runs, a number
              cell of 1,000 or more shown grouped (its value as before); anything else
              fails. What is printed is a fixed point of the grouping, and the numbers of
              1,000 or more it leaves alone (years, dates, identifiers) are listed. The
              supplementary figures differ from commit 1 only where a tick label of 1,000
              stands; the main figures' line widths before and after (the stage's drawing
              with and without its minimum), Figure 2's shading, crosses and Figure 3's
              letters. --list writes every grouped cell to a TSV.
  --commit 1  E8.1-E8.2. Numbers: the configuration's against E8.1's table of new and old
              numbers. Every CSV the stage writes, contents.csv, manifest.csv and every
              sheet of the workbook, cell by cell against the base, each change explained
              by one of the declared kinds: an item's number (the renumbering, applied to
              every "Supplementary Table/Figure S<n>" the base prints), the "provisional"
              marks, a title's revision code or v5 pointer, a part's replacement of a
              source's text (config/layout.yaml's `replace`), the three repository-only
              texts; anything else fails. Figures pixel for pixel as before. What is
              printed holds no "provisional", revision code, "v5" or "Table 3(b)"; and in
              the layout drawn with every item renumbered (Q1, Q2, ...), no printed text,
              a figure's included, still names an item of this paper by its number: what
              remains are the source papers' numbers, listed.

Run from the repository root:
    .venv/bin/python scripts/stage_E8/check_e8.py --commit 1 [--base d6f7bc5]
or in the reference image:
    ... IMAGE bash scripts/stage_E8/in_image.sh scripts/stage_E8/check_e8.py --commit 1
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

PHASE1 = REPO / "phase1"
OUT = Path("reports/evidence/layout")

# E8.1's table: new number <- content (the instruction's order, S1 first)
NEW_TABLES = ["functional_standards", "label_sources", "analysis_set_by_status",
              "predictor_columns", "training_signal", "column_agreement", "published_cut_points",
              "thresholds_pathogenic", "thresholds_benign", "clinvar_status", "within_gene_arms",
              "arms_at_fitted_thresholds", "external_genes", "new_external_genes",
              "external_column_basis", "external_summary", "practical_metrics",
              "confusion_matrices", "without_bap1", "heldout_intervals", "dilution",
              "label_sensitivity", "family_holdout", "fusion", "clinvar_coverage",
              "clinvar_consistency", "inframe_attribution", "inframe_true_positives",
              "twelve_gene_refit", "prior_sensitivity", "bp4_at_bp7_positions", "vcep_rules",
              "alphagenome_definitions", "by_splice_site_position"]
NEW_FIGURES = ["local_lr_curves", "threshold_dilution", "precision_recall_curves",
               "clinvar_coverage_by_position", "auroc_by_clinvar_status"]
# E8.1's table: the old numbers, in the same order
OLD_TABLES = [1, 17, 2, 3, 14, 4, 5, 6, 7, 10, 11, 15, 12, 19, 20, 21, 16, 31, 29, 26, 8, 24, 30,
              9, 23, 22, 13, 33, 32, 25, 27, 34, 18, 28]
OLD_FIGURES = [1, 2, 5, 4, 3]

RESIDUES = {
    "provisional": re.compile(r"provisional", re.I),
    "revision code": re.compile(r"\b[CE]\d{1,2}[a-z]?(?:\.\d+)?\b|\b[CE]\d{1,2}[a-z]?'s"),
    "v5": re.compile(r"\bv5\b"),
    "Table 3(b)": re.compile(r"Table 3\(b\)"),
}
# a number of this paper, or of a source paper, as printed
NUMBERED = re.compile(r"(?:Supplementary\s+)?(?:Table|Figure|Note|Fig\.)s?\s+S?\d+[a-z]?(?:\([a-z]\))?")


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _git_show(commit: str, path: str) -> bytes | None:
    r = subprocess.run(["git", "-C", str(REPO), "show", f"{commit}:{path}"], capture_output=True)
    return r.stdout if r.returncode == 0 else None


def _git_ls(commit: str, path: str) -> list[str]:
    out = subprocess.run(["git", "-C", str(REPO), "ls-tree", "-r", "--name-only", commit, "--", path],
                         check=True, capture_output=True, text=True).stdout
    return [p for p in out.splitlines() if p]


def csv_cells(data: bytes) -> tuple[list[str], list[str], list[list[str]]]:
    """A CSV as the stage writes it: its comment lines, its header, its rows (text)."""
    text = data.decode("utf-8")
    lines = text.splitlines(keepends=True)
    comments = []
    while lines and lines[0].startswith("#"):
        comments.append(lines.pop(0).rstrip("\r\n"))
    rows = list(csv.reader(io.StringIO("".join(lines))))
    return comments, (rows[0] if rows else []), rows[1:]


def _unescape(s: str) -> str:
    return s.replace("&lt;", "<").replace("&gt;", ">").replace("&quot;", '"').replace("&amp;", "&")


def sheets(data: bytes, styled: bool = False) -> dict[str, dict[tuple[int, str], str]]:
    """The workbook's sheets in order, each {(row, column): text} (an inline string's
    runs joined, a number as written); with `styled`, each number as (value, style) and
    each string as (text, the runs set in italics)."""
    z = zipfile.ZipFile(io.BytesIO(data))
    names = re.findall(r'<sheet name="([^"]+)" sheetId="(\d+)"', z.read("xl/workbook.xml").decode())
    out = {}
    for name, i in names:
        xml = z.read(f"xl/worksheets/sheet{i}.xml").decode()
        cells = {}
        for ref, attrs, body in re.findall(r'<c r="([A-Z]+\d+)"([^>]*)>(.*?)</c>', xml):
            m = re.match(r"([A-Z]+)(\d+)", ref)
            v = re.search(r"<v>(.*?)</v>", body)
            if v:
                s = re.search(r's="(\d+)"', attrs)
                val = (v.group(1), int(s.group(1)) if s else 0) if styled else v.group(1)
            else:
                runs = re.findall(r"<r>(?:<rPr>(.*?)</rPr>)?<t[^>]*>(.*?)</t></r>", body)
                text = ("".join(_unescape(x) for _, x in runs) if runs
                        else _unescape("".join(re.findall(r"<t[^>]*>(.*?)</t>", body))))
                val = (text, [_unescape(x) for p, x in runs if "<i/>" in (p or "")]) if styled else text
            cells[(int(m.group(2)), m.group(1))] = val
        out[_unescape(name)] = cells
    return out


def printed_texts(root: Path) -> list[tuple[str, str]]:
    """Every text the layout prints under `root`: (where, text) for each CSV cell, header
    and comment line (contents.csv and manifest.csv among them) and each workbook cell."""
    out = []
    for p in sorted(root.rglob("*")):
        rel = str(p.relative_to(root))
        if p.suffix == ".csv":
            comments, header, rows = csv_cells(p.read_bytes())
            out += [(f"{rel} comment {i}", c) for i, c in enumerate(comments)]
            out += [(f"{rel} header {j}", c) for j, c in enumerate(header)]
            out += [(f"{rel} row {i} col {j}", c) for i, r in enumerate(rows, 1) for j, c in enumerate(r)]
        elif p.suffix == ".xlsx":
            for name, cells in sheets(p.read_bytes()).items():
                out += [(f"{rel}:{name}!{c}{r}", v) for (r, c), v in cells.items()]
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--commit", type=int, required=True, choices=[1, 2, 3])
    ap.add_argument("--base", default=None, help="the commit before (default: HEAD)")
    ap.add_argument("--list", default=None, help="commit 2: a TSV of every grouped cell")
    args = ap.parse_args()
    base = args.base or subprocess.run(["git", "-C", str(REPO), "rev-parse", "--short", "HEAD"],
                                       check=True, capture_output=True, text=True).stdout.strip()
    if args.list:                               # from the folder it is run in, before the next line
        args.list = str(Path(args.list).resolve())
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
        numbers = L._numbers(cfg)
        old_numbers = {it["content"]: it["number"] for it in old_cfg["main"] + old_cfg["supplementary"]}
        if args.commit == 1:
            commit1(cfg, old_cfg, numbers, old_numbers, base, check, log, np, Image, L)
        if args.commit == 2:
            commit2(cfg, base, check, log, np, Image, L, args.list)
        if failed:
            print(f"{len(failed)} check(s) failed: {failed}")
            raise SystemExit(1)
        print("0 check(s) failed")


# ------------------------------------------------------------------------- commit 1
def commit1(cfg, old_cfg, numbers, old_numbers, base, check, log, np, Image, L) -> None:
    print("E8.1 numbering")
    got_t = [it["content"] for it in cfg["supplementary"] if it["kind"] == "table"]
    got_f = [it["content"] for it in cfg["supplementary"] if it["kind"] == "figure"]
    check("tables_numbered_as_the_instruction", got_t == NEW_TABLES,
          [f"S{k} {c}" for k, c in enumerate(got_t, 1) if k > len(NEW_TABLES) or NEW_TABLES[k - 1] != c])
    check("figures_numbered_as_the_instruction", got_f == NEW_FIGURES)
    check("old_numbers_as_the_instruction",
          [old_numbers[c] for c in NEW_TABLES] == [f"Supplementary Table S{n}" for n in OLD_TABLES]
          and [old_numbers[c] for c in NEW_FIGURES] == [f"Supplementary Figure S{n}" for n in OLD_FIGURES])
    check("main_numbers_unchanged", all(numbers[it["content"]] == old_numbers[it["content"]]
                                        for it in cfg["main"]))
    mapping = {old_numbers[c]: numbers[c] for c in numbers}
    log.note("new_from_old", {numbers[c]: old_numbers[c] for c in numbers if numbers[c] != old_numbers[c]})

    def renumber(text: str) -> str:
        return re.sub(r"Supplementary (Table|Figure) S\d+",
                      lambda m: mapping.get(m.group(0), m.group(0)), text)

    def cleanup(text: str, title: bool = False) -> str:
        text = text.replace("provisional: ", "")
        if title:
            text = re.sub(r" \(C\d{1,2}\)$", "", text)
            text = text.replace(" (v5 Table 2)", "").replace(" (v5 Table 4(b), (c))", "")
        return text

    repo_only = {"ROC and precision-recall curves in one figure (E15); drawn here as":
                 "ROC and precision-recall curves in one figure; drawn here as",
                 "the same curves at slide size (E4.1)": "the same curves at slide size",
                 "E19's draft of ": "an earlier draft of "}

    def repo_text(text: str) -> str:
        for a, b in repo_only.items():
            text = text.replace(a, b)
        return text

    items = {it["content"]: it for it in cfg["main"] + cfg["supplementary"]}
    reps_by_file, reps_by_item = {}, {}
    for it in items.values():
        for part in it.get("parts", []):
            reps = L._replacements(part, numbers)
            if not reps:
                continue
            name = (f"{it['content']}.csv" if len(it["parts"]) == 1
                    else f"{it['content']}.{part['part']}.csv")
            section = "main" if it in cfg["main"] else "supplementary"
            reps_by_file[f"{section}/{name}"] = reps
            reps_by_item.setdefault(it["content"], []).extend(reps)

    def apply(text: str, reps) -> str:
        for a, b, _ in reps:
            text = text.replace(a, b)
        return text

    # ------------------------------------------------------------- every CSV, cell by cell
    print("E8.6 every CSV against the base, cell by cell")
    changed_cells: dict[str, int] = {}
    unexplained: list = []
    base_files = {p.split("layout/", 1)[1] for p in _git_ls(base, f"phase1/{OUT}")}
    here_files = {str(p.relative_to(OUT)) for p in OUT.rglob("*") if p.is_file()}
    check("same_files_as_the_base", base_files == here_files, sorted(base_files ^ here_files))
    for rel in sorted(f for f in base_files & here_files if f.endswith(".csv")
                      and f not in ("contents.csv", "manifest.csv")):
        a = csv_cells(_git_show(base, f"phase1/{OUT}/{rel}"))
        b = csv_cells((OUT / rel).read_bytes())
        reps = reps_by_file.get(rel, [])
        n = 0
        if len(a[0]) != len(b[0]) or len(a[1]) != len(b[1]) or len(a[2]) != len(b[2]):
            unexplained.append((rel, "shape", [len(a[0]), len(a[1]), len(a[2])], [len(b[0]), len(b[1]), len(b[2])]))
            continue
        for x, y in zip(a[0] + a[1] + [c for r in a[2] for c in r], b[0] + b[1] + [c for r in b[2] for c in r]):
            if x != y:
                n += 1
                if apply(x, reps) != y:
                    unexplained.append((rel, x[:120], y[:120]))
        if n:
            changed_cells[rel] = n
    check("csv_changes_explained", not unexplained, unexplained[:10])
    log.note("csv_changed_cells", changed_cells)
    print("    changed cells:", changed_cells)

    # ------------------------------------------------------------- contents.csv
    print("E8.6 contents.csv")
    a = list(csv.DictReader(io.StringIO(_git_show(base, f"phase1/{OUT}/contents.csv").decode())))
    b = list(csv.DictReader(io.StringIO((OUT / "contents.csv").read_text())))
    key = lambda r: r["content"] or r["sources"]                              # noqa: E731
    ka, kb = {key(r): r for r in a}, {key(r): r for r in b}
    check("contents_same_rows", set(ka) == set(kb) and len(a) == len(b))
    want_order = ([it["content"] for it in cfg["main"]] + [it["content"] for it in cfg["supplementary"]]
                  + [L._rel(r["path"]) for r in cfg["repository_only"]])
    check("contents_in_number_order", [key(r) for r in b] == want_order)
    bad, per_field = [], {}
    for k in ka:
        ra, rb = ka[k], kb[k]
        for f in ra:
            x, y = ra[f], rb[f]
            if f == "number":
                want = mapping.get(x, x)
            elif f == "title":
                want = repo_text(cleanup(x, title=True)) if ra["number"] == "repository only" else cleanup(x, title=True)
            elif f in ("description", "notes"):
                want = renumber(cleanup(x))
            else:
                want = x
            if f == "title":
                want = renumber(want)
            if x != y:
                per_field[f] = per_field.get(f, 0) + 1
            if want != y:
                bad.append((k, f, x[:100], y[:100]))
    check("contents_changes_explained", not bad, bad[:10])
    log.note("contents_changed_fields", per_field)
    print("    changed fields:", per_field)

    # ------------------------------------------------------------- manifest.csv
    print("E8.6 manifest.csv")
    a = list(csv.DictReader(io.StringIO(_git_show(base, f"phase1/{OUT}/manifest.csv").decode())))
    b = list(csv.DictReader(io.StringIO((OUT / "manifest.csv").read_text())))
    ka, kb = {r["file"]: r for r in a}, {r["file"]: r for r in b}
    bad, per_field = [], {}
    differs = {f for f in base_files & here_files
               if _git_show(base, f"phase1/{OUT}/{f}") != (OUT / f).read_bytes()}
    for k in ka:
        ra, rb = ka[k], kb.get(k)
        if rb is None:
            bad.append((k, "row gone"))
            continue
        for f in ra:
            x, y = ra[f], rb[f]
            if x != y:
                per_field[f] = per_field.get(f, 0) + 1
            if f == "number":
                ok = y == mapping.get(x, x)
            elif f in ("bytes", "sha256"):
                ok = (x == y) or k.split("layout/", 1)[1] in differs
            elif f in ("sources", "sources_sha256") and k.endswith("supplementary_tables.xlsx"):
                # the workbook's sources are its tables', in the tables' order: the same
                # files and hashes, paired as before, in the new order
                pa = list(zip(ra["sources"].split("; "), ra["sources_sha256"].split("; ")))
                pb = list(zip(rb["sources"].split("; "), rb["sources_sha256"].split("; ")))
                want_order = list(dict.fromkeys(L._rel(s) for i in cfg["supplementary"]
                                                if i["kind"] == "table" for s in L._sources(i)))
                ok = sorted(pa) == sorted(pb) and [s for s, _ in pb] == want_order
            else:
                ok = x == y
            if not ok:
                bad.append((k, f, x[:80], y[:80]))
    written = [r["file"] for r in b]
    check("manifest_changes_explained", not bad and set(ka) == set(kb), bad[:10])
    check("manifest_hashes_true", all(_sha((REPO / r["file"]).read_bytes()) == r["sha256"]
                                      for r in b if not r["file"].endswith("manifest.csv")))
    log.note("manifest_changed_fields", per_field)
    log.note("files_that_differ_from_the_base", sorted(differs))
    print("    changed fields:", per_field)

    # ------------------------------------------------------------- the workbook
    print("E8.6 the workbook, sheet by sheet")
    sa = sheets(_git_show(base, f"phase1/{OUT}/supplementary_tables.xlsx"))
    sb = sheets((OUT / "supplementary_tables.xlsx").read_bytes())
    tables = [it for it in cfg["supplementary"] if it["kind"] == "table"]
    check("sheets_in_number_order", list(sb) == ["README"] + [L._sheet_name(it["number"]) for it in tables])
    bad, per_sheet = [], {}
    for it in tables:
        old_name = L._sheet_name(old_numbers[it["content"]])
        new_name = L._sheet_name(it["number"])
        ca, cb = sa[old_name], sb[new_name]
        if set(ca) != set(cb):
            bad.append((new_name, "cells differ in place", sorted(set(ca) ^ set(cb))[:5]))
            continue
        n = 0
        reps = reps_by_item.get(it["content"], [])
        for ref in ca:
            x, y = ca[ref], cb[ref]
            # the title row and the item's notes are the configuration's (renumbered,
            # cleaned); every other cell is a part's (its source's text, replaced)
            if ref == (1, "A"):
                wants = {renumber(cleanup(x, title=True))}
            else:
                wants = {renumber(cleanup(x)), apply(x, reps)}
            if x != y:
                n += 1
            if y not in wants:
                bad.append((new_name, f"{ref[1]}{ref[0]}", x[:100], y[:100]))
        if n:
            per_sheet[f"{new_name} ({it['content']}, was {old_name})"] = n
    # README: the header lines, then a row per table (by content), then the repository's
    ra, rb = sa["README"], sb["README"]
    rows_a = {r: {c: v for (rr, c), v in ra.items() if rr == r} for r in sorted({r for r, _ in ra})}
    rows_b = {r: {c: v for (rr, c), v in rb.items() if rr == r} for r in sorted({r for r, _ in rb})}
    n_readme = 0
    head_a, head_b = [rows_a.get(r, {}) for r in range(1, 5)], [rows_b.get(r, {}) for r in range(1, 5)]
    want_head = [{"A": f"Supplementary tables ({cfg['supplementary_format']['numbering']}; titles and "
                       "descriptions: config/layout.yaml)"},
                 {"A": head_a[1]["A"].replace("Assembly at submission (provisional):", "Assembly at submission:")},
                 {}, {c: v.replace(" (provisional)", "") for c, v in head_a[3].items()}]
    check("readme_header_as_declared", head_b == want_head, [h for h, w in zip(head_b, want_head) if h != w])
    n_readme += sum(1 for h, w in zip(head_a, head_b) for c in set(h) | set(w) if h.get(c) != w.get(c))
    by_content_a = {r["C"]: r for k, r in rows_a.items() if k >= 5 and "C" in r and r.get("A", "").startswith("Supplementary")}
    by_content_b = {r["C"]: (k, r) for k, r in rows_b.items() if k >= 5 and "C" in r and r.get("A", "").startswith("Supplementary")}
    check("readme_rows_in_number_order", [r[1]["A"] for r in sorted(by_content_b.values())]
          == [it["number"] for it in tables])
    for c, r in by_content_a.items():
        k, s = by_content_b[c]
        want = dict(r)
        want["A"], want["B"] = mapping[r["A"]], L._sheet_name(mapping[r["A"]])
        want["D"] = renumber(cleanup(r["D"], title=True))
        want["F"] = renumber(cleanup(r["F"]))
        n_readme += sum(1 for col in r if r[col] != s.get(col))
        if want != s:
            bad.append(("README", c, {x: (r.get(x), s.get(x)) for x in set(r) | set(s) if want.get(x) != s.get(x)}))
    tail_a = [r for k, r in rows_a.items() if k >= 5 and not r.get("A", "").startswith("Supplementary")]
    tail_b = [r for k, r in rows_b.items() if k >= 5 and not r.get("A", "").startswith("Supplementary")]
    want_tail = [{c: renumber(repo_text(v)) for c, v in r.items()} for r in tail_a]
    check("readme_repository_rows_as_declared", tail_b == want_tail,
          [(x, y) for x, y in zip(tail_b, want_tail) if x != y][:5])
    n_readme += sum(1 for x, y in zip(tail_a, tail_b) for c in x if x[c] != y.get(c))
    per_sheet["README"] = n_readme
    check("workbook_changes_explained", not bad, bad[:10])
    log.note("workbook_changed_cells", per_sheet)
    print("    changed cells:", per_sheet)

    # ------------------------------------------------------------- figures, pixel by pixel
    print("E8.6 figures, pixel for pixel")
    same = {}
    for rel in sorted(f for f in here_files if f.endswith(".png")):
        a = np.asarray(Image.open(io.BytesIO(_git_show(base, f"phase1/{OUT}/{rel}"))).convert("RGBA"))
        b = np.asarray(Image.open(OUT / rel).convert("RGBA"))
        same[rel] = a.shape == b.shape and bool((a == b).all())
    check("figures_pixel_for_pixel_as_before", all(same.values()), [k for k, v in same.items() if not v])

    # ------------------------------------------------------------- residues
    print("E8.6 residues in what is printed")
    texts = printed_texts(OUT)
    for name, pat in RESIDUES.items():
        hits = [(w, m.group(0)) for w, t in texts for m in pat.finditer(t)]
        check(f"no_{name.replace(' ', '_').replace('(', '').replace(')', '')}", not hits, hits[:8])
    # ------------------------------------------------------------- generated numbers only
    print("E8.1 every number of this paper printed through the configuration")
    with tempfile.TemporaryDirectory() as tmp:
        q_texts, fig_texts = q_render(Path(tmp), L)
    left = sorted({(w.split(" ")[0].split("!")[0], m.group(0), t) for w, t in q_texts + fig_texts
                   for m in NUMBERED.finditer(t)})
    ours = [x[:2] + (x[2][:140],) for x in left
            if re.match(r"(Supplementary )?(Table|Figure)s? S?\d", x[1]) and not SOURCE_PAPERS.search(x[2])]
    check("no_number_of_this_paper_left_fixed", not ours, ours[:10])
    log.note("numbers_left_when_renumbered", [[x[0], x[1], x[2][:200]] for x in left])
    print("    left in the renumbered layout (the source papers' numbers):")
    for x in left:
        print(f"      {x[0]}: «{x[1]}»  {x[2][:110]}")


# ------------------------------------------------------------------------- commit 2
def _ungroup(s: str) -> str:
    """A text with the commas between digit groups taken out."""
    return re.sub(r"(?<=\d),(?=\d{3}\b)", "", s)


TX_COLS = {"supplementary/functional_standards.csv": "Transcript",
           "supplementary/new_external_genes.records.csv": "Target transcript accession"}
MANE = "data/evidence/reference_cache/mane/MANE.GRCh38.v1.5.summary.txt.gz"
LEFT_ALONE = re.compile(r"(?<![\w.,/:#+*\-])(\d{4,})(\.\d+)?(?!\w|,\d|\.\d)")


def commit2(cfg, base, check, log, np, Image, L, list_path) -> None:
    mane = L._mane_select(MANE)
    grouped_cells: list = []
    transcript_cells: list = []

    def explained(rel, where, x, y, col=None, gene=None) -> bool:
        if x == y:
            return True
        if col is not None and TX_COLS.get(rel) == col:
            r = mane[gene]
            ok = x == r["Ensembl_nuc"] and y == f"{r['RefSeq_nuc']} ({x})"
            if ok:
                transcript_cells.append((rel, gene, x, y))
            return ok
        if y == L.grouped(x) and _ungroup(y) == _ungroup(x):
            grouped_cells.append((rel, where, x, y))
            return True
        return False

    print("E8.6 every CSV against commit 1, cell by cell")
    base_files = {p.split("layout/", 1)[1] for p in _git_ls(base, f"phase1/{OUT}")}
    here_files = {str(p.relative_to(OUT)) for p in OUT.rglob("*") if p.is_file()}
    check("same_files_as_the_base", base_files == here_files, sorted(base_files ^ here_files))
    bad, per_file = [], {}
    for rel in sorted(f for f in base_files & here_files if f.endswith(".csv")
                      and f not in ("contents.csv", "manifest.csv")):
        a = csv_cells(_git_show(base, f"phase1/{OUT}/{rel}"))
        b = csv_cells((OUT / rel).read_bytes())
        if [len(a[0]), len(a[1]), len(a[2])] != [len(b[0]), len(b[1]), len(b[2])]:
            bad.append((rel, "shape"))
            continue
        n = 0
        gi = a[1].index("Gene") if "Gene" in a[1] else None
        for k, (x, y) in enumerate(zip(a[0], b[0])):
            n += x != y
            if not explained(rel, f"comment {k}", x, y):
                bad.append((rel, f"comment {k}", x[:90], y[:90]))
        for j, (x, y) in enumerate(zip(a[1], b[1])):
            n += x != y
            if not explained(rel, f"header {j}", x, y):
                bad.append((rel, f"header {j}", x[:90], y[:90]))
        for i, (ra, rb) in enumerate(zip(a[2], b[2]), 1):
            for j, (x, y) in enumerate(zip(ra, rb)):
                n += x != y
                if not explained(rel, f"row {i} {a[1][j]}", x, y, a[1][j], ra[gi] if gi is not None else None):
                    bad.append((rel, f"row {i} col {a[1][j]}", x[:90], y[:90]))
        if n:
            per_file[rel] = n
    check("csv_changes_explained", not bad, bad[:10])
    log.note("csv_changed_cells", per_file)
    print("    changed cells:", per_file)
    t4 = csv_cells((OUT / "main/practical_performance.clinvar.csv").read_bytes())
    t4_old = csv_cells(_git_show(base, f"phase1/{OUT}/main/practical_performance.clinvar.csv"))
    five = [(x, y) for ra, rb in zip(t4_old[2], t4[2]) for x, y in zip(ra, rb) if x != y]
    check("table_4b_five_cells", five == [("1410", "1,410"), ("0.87 (1223/1410)", "0.87 (1,223/1,410)"),
                                          ("0.10 (138/1410)", "0.10 (138/1,410)"), ("1410", "1,410"),
                                          ("0.05 (68/1410)", "0.05 (68/1,410)")], five)
    check("transcript_cells", sorted(transcript_cells) == sorted([
        ("supplementary/functional_standards.csv", "DDX3X", "ENST00000644876.2", "NM_001356.5 (ENST00000644876.2)"),
        ("supplementary/new_external_genes.records.csv", "SBDS", "ENST00000246868.7", "NM_016038.4 (ENST00000246868.7)")]),
        transcript_cells)
    same_row = {g: (mane[g]["RefSeq_nuc"], mane[g]["Ensembl_nuc"], mane[g]["MANE_status"]) for g in ("DDX3X", "SBDS")}
    log.note("mane_rows", same_row)
    check("ddx3x_and_sbds_on_one_mane_select_row", same_row == {
        "DDX3X": ("NM_001356.5", "ENST00000644876.2", "MANE Select"),
        "SBDS": ("NM_016038.4", "ENST00000246868.7", "MANE Select")}, same_row)

    print("E8.6 contents.csv")
    a = {r["content"] or r["sources"]: r for r in csv.DictReader(io.StringIO(_git_show(base, f"phase1/{OUT}/contents.csv").decode()))}
    b = {r["content"] or r["sources"]: r for r in csv.DictReader(io.StringIO((OUT / "contents.csv").read_text()))}
    notes = {"functional_standards": "Transcripts are MANE Select (v1.5); for DDX3X, the Ensembl transcript named by the deposit is given in parentheses.",
             "new_external_genes": "Transcripts are MANE Select (v1.5); for SBDS, the Ensembl transcript named by the deposit is given in parentheses."}
    bad, per_field, estimates = [], {}, []
    for k in a:
        for f in a[k]:
            x, y = a[k][f], b[k][f]
            if x == y:
                continue
            per_field[f] = per_field.get(f, 0) + 1
            if f == "notes" and k in notes:
                ok = y == " | ".join(([notes[k]] if True else []) + ([x] if x else []))
            elif f == "sources" and k in notes:
                ok = with_mane(x, y, "phase1/" + MANE)
            elif f == "fits_one_page_estimate":
                ok = y.split(" (")[0] == x.split(" (")[0]
                estimates.append((k, x, y))
            else:
                ok = y == L.grouped(x)
            if not ok:
                bad.append((k, f, x[:90], y[:90]))
    check("contents_changes_explained", not bad and set(a) == set(b), bad[:10])
    check("no_destination_estimate_flips", all(y.split(" (")[0] == x.split(" (")[0] for _, x, y in estimates))
    log.note("contents_changed_fields", per_field)
    log.note("contents_estimates", estimates)
    print("    changed fields:", per_field)

    print("E8.6 manifest.csv")
    a = {r["file"]: r for r in csv.DictReader(io.StringIO(_git_show(base, f"phase1/{OUT}/manifest.csv").decode()))}
    b = {r["file"]: r for r in csv.DictReader(io.StringIO((OUT / "manifest.csv").read_text()))}
    differs = {f for f in base_files & here_files if _git_show(base, f"phase1/{OUT}/{f}") != (OUT / f).read_bytes()}
    mane_sha = hashlib.sha256((PHASE1 / MANE).read_bytes()).hexdigest()
    bad, per_field = [], {}
    for k in a:
        for f in a[k]:
            x, y = a[k][f], b[k][f]
            if x == y:
                continue
            per_field[f] = per_field.get(f, 0) + 1
            rel = k.split("layout/", 1)[1]
            if f in ("bytes", "sha256"):
                ok = rel in differs
            elif f == "sources":
                ok = (b[k]["content"] in notes or rel == "supplementary_tables.xlsx") and \
                    with_mane(x, y, "phase1/" + MANE)
            elif f == "sources_sha256":
                ok = (b[k]["content"] in notes or rel == "supplementary_tables.xlsx") and with_mane(x, y, mane_sha)
            else:
                ok = False
            if not ok:
                bad.append((rel, f, x[:80], y[:80]))
    check("manifest_changes_explained", not bad and set(a) == set(b), bad[:10])
    log.note("manifest_changed_fields", per_field)
    log.note("files_that_differ_from_the_base", sorted(differs))
    print("    changed fields:", per_field, "| files that differ:", len(differs))

    print("E8.6 the workbook, sheet by sheet")
    sa = sheets(_git_show(base, f"phase1/{OUT}/supplementary_tables.xlsx"), styled=True)
    sb = sheets((OUT / "supplementary_tables.xlsx").read_bytes(), styled=True)
    check("same_sheets", list(sa) == list(sb))
    by_sheet = {L._sheet_name(it["number"]): it for it in cfg["supplementary"] if it["kind"] == "table"}
    bad, per_sheet, styled_numbers, italics = [], {}, 0, []
    for name in sa:
        ca, cb = sa[name], dict(sb[name])
        it = by_sheet.get(name)
        if it and it["content"] in notes:
            # the table's note is the sheet's second row; the rows below it move down one
            check(f"{name}_note_row", cb.get((2, "A"), ("",))[0] == notes[it["content"]], cb.get((2, "A")))
            cb = {(r - 1 if r > 2 else r, c): v for (r, c), v in cb.items() if r != 2}
        if set(ca) != set(cb):
            bad.append((name, "cells differ in place", sorted(set(ca) ^ set(cb))[:4]))
            continue
        n = 0
        for ref in ca:
            x, y = ca[ref], cb[ref]
            if x == y:
                continue
            n += 1
            num_x, num_y = re.fullmatch(r"-?\d+(\.\d+)?", str(x[0])) and isinstance(x[1], int), \
                re.fullmatch(r"-?\d+(\.\d+)?", str(y[0])) and isinstance(y[1], int)
            if num_y:
                d = len((y[0].split(".") + [""])[1])
                value_x = x[0] if num_x else _ungroup(x[0])
                ok = (float(value_x) == float(y[0]) and abs(float(y[0])) >= 1000
                      and y[1] == L.GROUPED_STYLE + d and (not num_x or x[1] == L.NUMBER_STYLE + d))
                styled_numbers += ok
            elif x[0] in {mane[g]["Ensembl_nuc"] for g in ("DDX3X", "SBDS")}:
                g = next(g for g in ("DDX3X", "SBDS") if mane[g]["Ensembl_nuc"] == x[0])
                ok = y[0] == f"{mane[g]['RefSeq_nuc']} ({x[0]})" and name == L._sheet_name(
                    next(i["number"] for i in cfg["supplementary"] if i["content"] == {"DDX3X": "functional_standards", "SBDS": "new_external_genes"}[g]))
                if ok:
                    transcript_cells.append(("supplementary_tables.xlsx", g, x[0], y[0]))
            elif name == "README" and ref[1] in ("H", "I") and it is None:
                row_content = cb.get((ref[0], "C"), ("",))[0]
                ok = row_content in notes and with_mane(x[0], y[0], "phase1/" + MANE if ref[1] == "H" else mane_sha)
            else:
                plain = L.ITALIC.sub(r"\1", L.grouped(x[0]))
                ok = y[0] == plain and y[1] == [m.group(1) for m in L.ITALIC.finditer(L.grouped(x[0]))]
                if y[1]:
                    italics.append((name, f"{ref[1]}{ref[0]}", y[1]))
                elif ok and y[0] != x[0]:
                    grouped_cells.append(("supplementary_tables.xlsx", f"{name}!{ref[1]}{ref[0]}", x[0], y[0]))
            if not ok:
                bad.append((name, f"{ref[1]}{ref[0]}", str(x)[:90], str(y)[:90]))
        if n:
            per_sheet[name] = n
    check("workbook_changes_explained", not bad, bad[:10])
    log.note("workbook_changed_cells", per_sheet)
    log.note("workbook_numbers_now_shown_grouped", styled_numbers)
    log.note("workbook_italic_runs", italics)
    print("    changed cells:", per_sheet)
    print("    numbers shown grouped:", styled_numbers, "| cells with italic runs:", len(italics))

    print("E8.6 what is printed, grouped; what the grouping leaves alone")
    texts = printed_texts(OUT)
    # a workbook cell that is a number shows its thousands by its style (checked above)
    unfixed = [(w, t[:100]) for w, t in texts if L.grouped(t) != t and not re.fullmatch(r"[\da-f]{64}(; [\da-f]{64})*", t)
               and not w.startswith("manifest.csv") and not (".xlsx:" in w and re.fullmatch(r"-?\d+(\.\d+)?", t))]
    check("printed_text_is_grouped", not unfixed, unfixed[:8])
    alone = sorted({(w.split(" ")[0].split("!")[0].split(":")[0], m.group(0), t[max(0, m.start() - 30):m.end() + 12])
                    for w, t in texts if not w.startswith("manifest.csv")
                    and not re.fullmatch(r"[\da-f]{64}(; [\da-f]{64})*", t) and not re.fullmatch(r"-?\d+(\.\d+)?", t)
                    for m in LEFT_ALONE.finditer(t)})
    log.note("left_alone", [list(x) for x in alone])
    print(f"    numbers of 1,000 or more left as written: {len(alone)} (years, dates, identifiers; listed in the run record)")
    if list_path:
        Path(list_path).parent.mkdir(parents=True, exist_ok=True)
        with open(list_path, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(["file", "where", "before", "after"])
            w.writerows(grouped_cells)
        print(f"    {len(grouped_cells)} grouped cells listed in {list_path}")
    log.note("grouped_cells", len(grouped_cells))

    print("E8.4 the supplementary figures: only the tick labels of 1,000")
    with tempfile.TemporaryDirectory() as tmp:
        figs = draw(L, Path(tmp), heavier=True)
        before = draw(L, Path(tmp) / "before", heavier=False)
    for rel in sorted(f for f in here_files if f.endswith(".png") and f.startswith("supplementary/")):
        a = np.asarray(Image.open(io.BytesIO(_git_show(base, f"phase1/{OUT}/{rel}"))).convert("RGBA"))
        b = np.asarray(Image.open(OUT / rel).convert("RGBA"))
        stem = Path(rel).stem
        differ = (a != b).any(axis=2) if a.shape == b.shape else None
        allowed = np.zeros(b.shape[:2], dtype=bool)
        boxes = tick_boxes(figs[stem], b.shape[0])
        for x0, y0, x1, y1 in boxes:
            allowed[max(0, y0):y1, max(0, x0):x1] = True
        ok = differ is not None and not (differ & ~allowed).any()
        check(f"{stem}_differs_only_at_its_1000_labels", ok,
              {"pixels_differing": int(differ.sum()) if differ is not None else None, "labels": len(boxes)})
    print("E8.4 the main figures")
    widths = {stem: {"before": line_classes(before[stem]), "after": line_classes(figs[stem])}
              for stem in ("published_cut_point_by_band", "heldout_and_external_lr", "roc_curves")}
    log.note("main_line_widths", widths)
    for stem, w in widths.items():
        print(f"    {stem}:")
        for cls in w["after"]:
            print(f"      {cls:22s} {w['before'].get(cls)} -> {w['after'][cls]}")
        check(f"{stem}_lines_at_least_{cfg['artwork']['main_lines_pt']}_pt",
              all(lo >= cfg["artwork"]["main_lines_pt"] for lo, hi in w["after"].values()))
    from matplotlib.colors import to_hex
    from src import evid_figures as F
    f2 = figs["heldout_and_external_lr"]
    check("figure2_no_row_shading", not [p for ax in f2.axes for p in ax.patches
                                         if to_hex(p.get_facecolor()) == to_hex(F.BAND)])
    key = [h.get_marker() for lg in f2.legends for h in lg.legend_handles if h.get_label().startswith("Fold not counted")]
    check("figure2_folds_not_counted_as_crosses", key == [L.CROSS], key)
    letters = [t.get_text() for ax in figs["roc_curves"].axes for t in ax.texts if t.get_text() in ("a", "b", "c")]
    check("figure3_letters_a_to_c", letters == ["a", "b", "c"], letters)
    for f in figs.values():
        f.canvas.draw()
    labels = {s: sorted({x for x in drawn_texts(f) if re.search(r"\d{1,3}(,\d{3})+", x)}) for s, f in figs.items()}
    log.note("figure_texts_grouped", labels)
    check("figure_numbers_of_1000_grouped", all(not re.search(r"(?<![\d,.])\d{4,}(?![\d,])", x)
                                                for f in figs.values() for x in drawn_texts(f)), labels)


def with_mane(x: str, y: str, entry: str) -> bool:
    """A "; "-joined list with one entry more: the MANE summary's path or hash, where the
    item's parts list it."""
    a, b = x.split("; "), y.split("; ")
    return len(b) == len(a) + 1 and b.count(entry) == 1 and [e for e in b if e != entry] == a


def drawn_ticklabels(ax, axis) -> list:
    """The tick labels an axis draws: visible, set, within the view."""
    lo, hi = sorted(axis.get_view_interval())
    return [t.label1 for t, x in zip(axis.get_major_ticks(), axis.get_majorticklocs())
            if lo - 1e-9 <= x <= hi + 1e-9 and t.label1.get_visible() and t.label1.get_text().strip()]


def drawn_texts(fig) -> list[str]:
    """Every text a figure draws: its axes' drawn tick labels and every other visible text."""
    import matplotlib.axes
    ticks = set()
    out = []
    for ax in fig.findobj(matplotlib.axes.Axes):
        for axis in (ax.xaxis, ax.yaxis):
            for t in axis.get_major_ticks() + axis.get_minor_ticks():
                ticks.add(id(t.label1))
                ticks.add(id(t.label2))
            out += [t.get_text() for t in drawn_ticklabels(ax, axis)]
    out += [t.get_text() for t in fig.findobj(plt_text()) if id(t) not in ticks and t.get_visible()
            and t.get_text().strip()]
    return out


def plt_text():
    import matplotlib.text
    return matplotlib.text.Text


def draw(L, folder: Path, heavier: bool) -> dict:
    """Every figure the stage draws, by its own makers into `folder`, each kept after
    its _save; without the main figures' line minimum if `heavier` is False."""
    import matplotlib
    kept, save, heavy = {}, L._save, L._heavier
    L._save = lambda fig, f, stem, cfg: (kept.__setitem__(stem, fig), save(fig, f, stem, cfg))[1]
    if not heavier:
        L._heavier = lambda fig, pt: None
    try:
        with matplotlib.rc_context():
            cfg = L.load_config()
            for it in cfg["main"] + cfg["supplementary"]:
                if it["kind"] == "figure":
                    L.MAKERS[it["make"]](it, folder / ("main" if it in cfg["main"] else "supplementary"), cfg)
    finally:
        L._save, L._heavier = save, heavy
    return kept


def tick_boxes(fig, height: int) -> list[tuple[int, int, int, int]]:
    """The pixel boxes, in the PNG at 300 dpi (rows from the top), of the tick labels that
    read 1,000 or more grouped, two pixels round."""
    import matplotlib.axes
    screen = fig.dpi
    fig.set_dpi(300)
    try:
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        out = []
        for ax in fig.findobj(matplotlib.axes.Axes):
            for axis in (ax.xaxis, ax.yaxis):
                for t in drawn_ticklabels(ax, axis):
                    if re.search(r"\d,\d{3}", t.get_text()):
                        e = t.get_window_extent(r)
                        out.append((int(e.x0) - 2, int(height - e.y1) - 2, int(e.x1) + 3, int(height - e.y0) + 3))
    finally:
        fig.set_dpi(screen)
    return out


def line_classes(fig) -> dict:
    """A figure's line widths by kind, (smallest, largest): its axes' spines, their major
    ticks, the lines drawn in its panels (a dotted or dashed line is a reference line),
    and its keys' lines."""
    import matplotlib.axes
    from matplotlib.lines import Line2D
    out: dict = {}

    def add(cls, w):
        lo, hi = out.get(cls, (w, w))
        out[cls] = (round(min(lo, w), 3), round(max(hi, w), 3))
    keyed = set()
    for lg in [*fig.legends, *[ax.get_legend() for ax in fig.findobj(matplotlib.axes.Axes) if ax.get_legend()]]:
        for ln in lg.get_lines():
            keyed.add(id(ln))
            if ln.get_linestyle() not in ("None", "", " "):
                add("key lines", ln.get_linewidth())
    for ax in fig.findobj(matplotlib.axes.Axes):
        for s in ax.spines.values():
            if s.get_visible():
                add("spines", s.get_linewidth())
        for axis in (ax.xaxis, ax.yaxis):
            for t in axis.get_major_ticks():
                if t.tick1line.get_visible() and t.tick1line.get_markersize() > 0:
                    add("ticks", t.tick1line.get_markeredgewidth())
        for ln in ax.lines:
            if id(ln) in keyed or ln.get_linestyle() in ("None", "", " "):
                continue
            add("reference lines" if ln.get_linestyle() not in ("-", "solid") else "data, medians, intervals",
                ln.get_linewidth())
    return out


# the source papers' own numbers, which the renumbered layout keeps: the label audit's
# rules and source tables (functional_standards' "Table S1" is this paper's, and keeps its
# number, E8.1), the second BRCA2 dataset's Supplementary Table 5
SOURCE_PAPERS = re.compile(r"Figure 3a legend|Figure 3 legend|Figure 4e legend|sheet 'Table S3'|"
                           r"Supplementary Table 3 'Functional category'|Table S1 cites no paper|"
                           r"Supplementary Table 5|Table 5")


def q_render(tmp: Path, L) -> tuple[list, list]:
    """The whole layout drawn again with each item's number in another scheme (Table Q2,
    Supplementary Table Q24, ...), into `tmp`: its printed texts, and every text its
    figures draw (taken from each figure as the stage saves it)."""
    import matplotlib.text as mtext
    orig_load, orig_save = L.load_config, L._save
    fig_texts: list = []

    def load(path=L.CONFIG):
        cfg = orig_load(path)
        for it in cfg["main"] + cfg["supplementary"]:
            it["number"] = re.sub(r"(Table|Figure) S?(\d+)", r"\1 Q\2", it["number"])
        return cfg

    def save(fig, folder, stem, cfg):
        fig.canvas.draw()
        fig_texts.extend((f"{stem} figure", t.get_text()) for t in fig.findobj(mtext.Text)
                         if t.get_visible() and t.get_text().strip())
        return orig_save(fig, folder, stem, cfg)

    L.load_config, L._save = load, save
    try:
        L.main(out=tmp / "layout")
    finally:
        L.load_config, L._save = orig_load, orig_save
    # a Q number is the configuration's: NUMBERED does not read it
    return printed_texts(tmp / "layout"), fig_texts


if __name__ == "__main__":
    main()
