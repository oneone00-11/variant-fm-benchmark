"""E8.13: the number gate (scripts/check_evidence_manuscript_numbers.py) on the v6 text in
each mode v5 was checked in (scripts/stage_E5/check_bindings_v5.py, which stays as it is),
and the cells mode on the main tables.

- pool mode, the weak check: on the manuscript as it stands (docx) and on the two texts
  build_bindings_v6.py writes (v6_main.md, v6_supplement.md);
- bound mode: each sentence's numbers against the files its source tag names (the keyed
  texts, tagged from the bindings), and the self-test, which plants numbers the named
  files do not hold and passes only if the bound mode reports each;
- cells mode on the two texts with docs/manuscript_bindings_v6.tsv: every measured number
  against its one cell, every constant with its reason;
- cells mode on the manuscript's Tables 1-4 and on the four editable files of E8.11
  (layout/main/<content>.docx): each data cell a paragraph, bound to the same cell of the
  layout's CSV (scripts/stage_E8/check_e8_part2.cells_mode).

Each mode's verdict is recorded whole, and so is a tally of every number written in the
texts by what the gate does with it: read (bound or a constant), or left out (headings,
identifiers, citations, labels, digits within a word or symbol), with the numbers the texts
leave out whole (the reference lists, the notes' list numbers). Exits 1 if the cells mode,
on the texts or on a table, finds a number that does not match its cell or that cannot be
bound.

    .venv/bin/python scripts/stage_E8/check_bindings_v6.py \\
        --manuscript ~/Desktop/calibration/CALInewversion.docx \\
        --texts ~/Desktop/整理文件/E8b_工作文件/v6_texts
"""
from __future__ import annotations

import argparse
import hashlib
import re
import subprocess
import sys
import tempfile
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog  # noqa: E402
from check_e8_part2 import cells_mode, docx_blocks  # noqa: E402

sys.path.insert(0, str(REPO / "scripts"))
import check_evidence_manuscript_numbers as gate  # noqa: E402

TABLE = REPO / "docs" / "manuscript_bindings_v6.tsv"
LAYOUT = REPO / "phase1" / "reports" / "evidence" / "layout"
# each main table's parts, in order, as the layout writes them
PARTS = {"Table 1": ["main/analysis_set_composition.by_gene.csv", "main/analysis_set_composition.by_clinvar_status.csv"],
         "Table 2": ["main/evidence_tiers_and_thresholds.tiers.csv", "main/evidence_tiers_and_thresholds.thresholds.csv"],
         "Table 3": ["main/published_cut_point_by_clinvar_status.csv"],
         "Table 4": ["main/practical_performance.seven_genes.csv", "main/practical_performance.clinvar.csv"]}
FILES = {"Table 1": "analysis_set_composition", "Table 2": "evidence_tiers_and_thresholds",
         "Table 3": "published_cut_point_by_clinvar_status", "Table 4": "practical_performance"}


def manuscript_tables(path: Path) -> dict[str, list[list[list[str]]]]:
    """The manuscript's Tables 1-4 (read only): after each paragraph "Table N.", its parts
    in order, each as its rows of cell texts."""
    import docx
    from docx.oxml.ns import qn
    doc = docx.Document(str(path))
    out: dict[str, list] = {}
    current = None
    for el in doc.element.body:
        if el.tag == qn("w:p"):
            text = "".join(t.text or "" for t in el.iter(qn("w:t"))).strip()
            for name in PARTS:
                if text.startswith(f"{name}. "):
                    current = name
                    out[name] = []
        elif el.tag == qn("w:tbl") and current:
            t = docx.table.Table(el, doc)
            out[current].append([[c.text for c in r.cells] for r in t.rows])
    return out


NUMBER = re.compile(r"\d+(?:[.,]\d+)*")


def excluded_in_text(text: str) -> dict:
    """Every number written in a text's body, by what the gate does with it: read (a
    token: bound or a constant), or left out as it cleans a sentence (gate.clean_sentence,
    in its order: a source tag, an identifier such as a DOI, URN or commit, a citation in
    brackets, a label such as "Supplementary Table S8"), in a heading, or within a word
    or a symbol (BRCA1, PP3, GRCh38, 5th). The tokens counted here must be the gate's."""
    start, end = gate.body_span(text)
    out = {"read": 0, "headings": 0, "source tags": 0, "identifiers": 0, "citations": 0,
           "labels": 0, "within words": 0, "outside the body": len(NUMBER.findall(text[:start] + text[end:]))}
    for _, para in gate._pieces(text[start:end], gate.PARAGRAPH_BREAK):
        if para.lstrip().startswith("#"):
            out["headings"] += len(NUMBER.findall(para))
            continue
        for _, s in gate._pieces(para, gate.SENTENCE_BREAK):
            for key, pat in (("source tags", gate.TAGGED), ("identifiers", gate.IDENTIFIER),
                             ("citations", re.compile(r"\[[^\]]*\]")), ("labels", gate.REF)):
                out[key] += sum(len(NUMBER.findall(m.group(0))) for m in pat.finditer(s))
                s = pat.sub(" ", s)
            read = len(list(gate.TOKEN.finditer(s)))
            out["read"] += read
            out["within words"] += max(0, len(NUMBER.findall(s)) - read)
    if out["read"] != len(gate.read_tokens(text)):
        raise SystemExit("the tally of numbers read is not the gate's")
    return out


def excluded_in_manuscript(path: Path) -> dict:
    """The numbers the texts leave out whole: the reference lists (the main text's and the
    notes'), and the numbers Word gives the notes' numbered items (not in the text)."""
    import docx
    from docx.oxml.ns import qn
    doc = docx.Document(str(path))
    out = {"main references": 0, "notes references": 0, "notes list numbers": 0}
    where = None
    for el in doc.element.body:
        if el.tag != qn("w:p"):
            continue
        p = docx.text.paragraph.Paragraph(el, doc)
        t, heading = p.text.strip(), p.style.name.startswith("Heading")
        if heading:
            where = {"References": "main references", "Supplementary Notes": "notes",
                     "References cited in the Supplementary Notes": "notes references"}.get(t, where)
            if t == "Tables":
                where = None
            continue
        if where in ("main references", "notes references"):
            out[where] += len(NUMBER.findall(t))
        if where == "notes" and el.find(f"{qn('w:pPr')}/{qn('w:numPr')}") is not None:
            out["notes list numbers"] += 1
    return out


def run_gate(*args: str) -> tuple[int, str]:
    r = subprocess.run([sys.executable, str(REPO / "scripts" / "check_evidence_manuscript_numbers.py"), *args],
                       capture_output=True, text=True, cwd=REPO)
    return r.returncode, r.stdout + r.stderr


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manuscript", required=True, type=Path)
    ap.add_argument("--texts", required=True, type=Path)
    args = ap.parse_args()
    ms, folder = args.manuscript.expanduser().resolve(), args.texts.expanduser().resolve()
    failed = []
    with RunLog(__file__) as log:
        log.note("manuscript_sha256", hashlib.sha256(ms.read_bytes()).hexdigest())
        log.note("bindings_sha256", hashlib.sha256(TABLE.read_bytes()).hexdigest())
        for name in ("v6_main.md", "v6_supplement.md", "v6_main_keyed.md", "v6_supplement_keyed.md"):
            log.note(f"{name}_sha256", hashlib.sha256((folder / name).read_bytes()).hexdigest())
        # pool mode
        for label, target in (("pool_manuscript", ms), ("pool_main", folder / "v6_main.md"),
                              ("pool_supplement", folder / "v6_supplement.md")):
            rc, out = run_gate(str(target), *(["--pool"] if target.suffix == ".md" else []))
            log.note(label, {"exit": rc, "output": out})
            print(f"== {label}: exit {rc}\n" + "\n".join(out.strip().splitlines()[-3:]))
        # bound mode and the self-test
        for doc in ("main", "supplement"):
            rc, out = run_gate(str(folder / f"v6_{doc}_keyed.md"))
            log.note(f"bound_{doc}", {"exit": rc, "output": out})
            print(f"== bound mode, {doc}: exit {rc}\n" + "\n".join(
                ln for ln in out.splitlines() if ln.startswith(("bound mode", "  NOT IN", "  UNBOUND"))))
        rc, out = run_gate(str(folder / "v6_main_keyed.md"), "--selftest")
        log.note("selftest", {"exit": rc, "output": out})
        print(f"== self-test: exit {rc}\n" + "\n".join(out.strip().splitlines()[-2:]))
        if rc:
            failed.append("selftest")
        # cells mode on the texts
        bindings = gate.read_bindings(TABLE)
        whole = excluded_in_manuscript(ms)
        log.note("excluded_whole", whole)
        print(f"== left out whole: {whole}")
        for doc in ("main", "supplement"):
            tally = excluded_in_text((folder / f"v6_{doc}.md").read_text())
            log.note(f"numbers_{doc}", tally)
            print(f"== numbers written in the {doc} text: {tally}")
            res = gate.cells_results((folder / f"v6_{doc}.md").read_text(), [b for b in bindings if b.document == doc])
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
        # cells mode on the main tables: the manuscript's, and the editable files'
        ms_tables = manuscript_tables(ms)
        with tempfile.TemporaryDirectory() as tmp:
            for name, parts in PARTS.items():
                got = ms_tables.get(name, [])
                if len(got) != len(parts):
                    failed.append(f"manuscript_{name}")
                    log.note(f"cells_manuscript_{name}", {"parts_found": len(got), "parts": len(parts)})
                    continue
                for label, tables in (
                        ("manuscript", dict(zip(parts, got))),
                        ("file", dict(zip(parts, [b[1] for b in docx_blocks(
                            (LAYOUT / "main" / f"{FILES[name]}.docx").read_bytes())[0] if b[0] == "tbl"])))):
                    g = cells_mode(tables, LAYOUT, Path(tmp), f"{label}_{name.replace(' ', '')}")
                    log.note(f"cells_{label}_{name}", g)
                    print(f"== cells mode, {label} {name}: {g['numbers']} numbers, {g['constant']} constants, "
                          f"{g['match']} match, {g['mismatch']} mismatch, {g['cannot_bind']} cannot bind, "
                          f"{g['not_in_table']} not in the table")
                    if g["mismatch"] or g["cannot_bind"] or g["not_in_table"]:
                        failed.append(f"cells_{label}_{name}")
        log.note("failed", failed)
        print(f"{len(failed)} check(s) failed: {failed}" if failed else "0 check(s) failed")
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
