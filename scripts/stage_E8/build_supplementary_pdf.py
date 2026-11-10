"""E8.12: the supplementary information as one PDF (the author's decision of 2026-11-09),
assembled from the manuscript, read only, and the layout's outputs. The PDF holds the
manuscript's supplementary notes, so it is not committed: this script and its test are.

In order: a front page (the heading, the paper's title, the author, and the contents with
page numbers, a range of numbers never broken at its dash, E8 commit 8); the conventions (config/layout.yaml's supplementary_format); the
supplementary notes word for word, as the manuscript holds them, from its heading
"Supplementary Notes" to the last entry of "References cited in the Supplementary Notes"
(headings, numbered items, run-in headings and italics as they are, the type set in Arial
at 10 pt); the supplementary figures S1-S5, each on a page of its own: the layout's vector
PDF at the size it was drawn (never scaled, which would print its 5 pt type smaller),
centred between side margins of 1.5 cm, with "Supplementary Figure Sx." and its title and
legend below it, or at the top of the next page where they do not fit below (S1, 240 mm
tall); then the tables marked pdf in the configuration, each its title, notes and parts,
the cells the layout's CSVs', ruled across only as the main tables' files are (E8.11),
turned where wide, a wide table's short columns kept whole and its long ones wrapping,
the gene symbols in their headers and cells in italics (evid_layout.italic_genes, E8 commit
8). A4; text pages with margins of 2 cm; page numbers at the foot.

The pages are typeset by LibreOffice from a document built here with python-docx on a copy
of the manuscript, so that the notes keep their own numbering and styles; a figure's page
holds a space as tall as the figure, onto which the figure's PDF is laid with pypdf, moved
and never scaled. A first typesetting, with a mark in each figure's space, finds the
figures' pages and the contents' page numbers; the second is the PDF.

Checks, recorded beside the PDF (<name>.check.json) or where --record says; the run
stops if one fails:
- the notes' text in the PDF is the manuscript's (LibreOffice's text of it, the list numbers
  included), but for white space;
- each of the five figures and the tables marked pdf appears once, in order, after the
  notes, whose headings stand in order; the contents' page numbers are where they are;
- the front page prints its text in order, the joiners aside, and breaks no range across
  two lines, after its dash or before it;
- the figures' pages hold no image, only vector drawing, and each figure is laid at a
  scale of 1 within its page;
- every font is embedded;
- the page count, and the file's size against the supplementary files' limit (150 MB).

    .venv/bin/python scripts/stage_E8/build_supplementary_pdf.py \\
        --manuscript ~/Desktop/calibration/CALInewversion.docx

Needs LibreOffice (soffice) and Arial (macOS: /System/Library/Fonts/Supplemental, given to
LibreOffice through fontconfig). Writes nothing into the repository.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
PHASE1 = REPO / "phase1"

OUT = Path.home() / "Desktop" / "calibration" / "submission" / "Supplementary_Information.pdf"
HEADING = "Supplementary Information"
PAPER = ("Evidence strength of splice-region variant effect predictors against saturation genome "
         "editing data, by distance from the exon and across genes")
AUTHOR = "Ningyi Zhang"
NOTES_FROM = "Supplementary Notes"
NOTES_REFS = "References cited in the Supplementary Notes"
FONT, TEXT_PT, TABLE_PT = "Arial", 10, 8
TEXT_MARGIN_MM, FIGURE_SIDE_MM, FOOTER_MM = 20.0, 15.0, 10.0
A4 = (210.0, 297.0)
LIMIT_MB = 150
# the styles the notes use, at their sizes here (pt)
STYLE_PT = {"Heading 2": 14, "Heading 3": 11}
W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
FONTS_CONF = """<?xml version="1.0"?>
<!DOCTYPE fontconfig SYSTEM "fonts.dtd">
<fontconfig>
  <dir>/System/Library/Fonts/Supplemental</dir>
  <dir>/System/Library/Fonts</dir>
  <dir>/Library/Fonts</dir>
  <cachedir>{cache}</cachedir>
</fontconfig>
"""
MARK = "FIGUREMARK{}"            # in a figure's space, first typesetting only
# a range of numbers, "11–50" or "S29–S34", is kept on one line by a word joiner on each
# side of its dash (U+2060, which prints nothing; E8 commit 8)
JOINER = "\u2060"
RANGE_DASH = re.compile(r"(?<=\d)–(?=[A-Z]?\d)")


# a range broken across two lines of a page's text, after its dash or before it
RANGE_BROKEN = re.compile(r"\d[\u2060 ]*–[\u2060 ]*\n[\u2060 ]*[A-Z]?\d|\d[\u2060 ]*\n[\u2060 ]*–[\u2060 ]*[A-Z]?\d")


def unbroken(text: str) -> str:
    """A text with each range's dash joined to the numbers on either side."""
    return RANGE_DASH.sub(JOINER + "–" + JOINER, text)


def twips(mm: float) -> int:
    return round(mm * 1440 / 25.4)


def qn(tag: str) -> str:
    return f"{{{W}}}{tag}"


def element(xml: str):
    """A WordprocessingML element from its text, the namespace declared on it."""
    from docx.oxml import parse_xml
    return parse_xml(re.sub(r"^<(w:\w+)", rf'<\1 xmlns:w="{W}"', xml, count=1))


def text_of(el) -> str:
    return "".join(t.text or "" for t in el.iter(qn("t")))


def style_of(el) -> str:
    s = el.find(f"{qn('pPr')}/{qn('pStyle')}")
    return s.get(qn("val")) if s is not None else ""


# ------------------------------------------------------------------ the notes, as written
def notes_range(body) -> tuple[int, int]:
    """The notes' span among the body's elements: from the paragraph "Supplementary Notes"
    to the last entry of "References cited in the Supplementary Notes" (the last
    paragraph with text before the next heading or the end)."""
    kids = list(body)
    paras = [(i, text_of(el).strip()) for i, el in enumerate(kids) if el.tag == qn("p")]
    start = [i for i, t in paras if t == NOTES_FROM]
    refs = [i for i, t in paras if t == NOTES_REFS]
    if len(start) != 1 or len(refs) != 1 or refs[0] < start[0]:
        raise SystemExit(f"the manuscript has {len(start)} paragraph(s) '{NOTES_FROM}' and "
                         f"{len(refs)} '{NOTES_REFS}', not one of each in that order")
    end = refs[0]
    for i in range(refs[0] + 1, len(kids)):
        el = kids[i]
        if el.tag != qn("p") or style_of(el).lower().startswith("heading"):
            break
        if text_of(el).strip():
            end = i
    if end == refs[0]:
        raise SystemExit(f"'{NOTES_REFS}' is followed by no entry")
    return start[0], end


def note_headings(body, span: tuple[int, int]) -> list[str]:
    """The notes' own headings, in order: "Supplementary Note 1. Predictor scores", ..."""
    kids = list(body)[span[0]:span[1] + 1]
    return [text_of(el).strip() for el in kids
            if el.tag == qn("p") and re.match(r"Supplementary Note \d+\. ", text_of(el).strip())]


def set_fonts(doc) -> None:
    """Every style, and the document's defaults, in Arial; the notes' paragraph styles at
    10 pt but their headings (STYLE_PT); a run's own font or size, where it sets one, left
    to its style."""
    from docx.shared import Pt
    styles = doc.styles.element
    defaults = styles.find(f"{qn('docDefaults')}/{qn('rPrDefault')}/{qn('rPr')}")
    for st in [defaults] + [s.find(qn("rPr")) for s in styles.findall(qn("style"))]:
        if st is None:
            continue
        for f in st.findall(qn("rFonts")):
            st.remove(f)
        st.insert(0, element(f'<w:rFonts w:ascii="{FONT}" w:hAnsi="{FONT}" w:cs="{FONT}" w:eastAsia="{FONT}"/>'))
    from docx.shared import RGBColor
    for s in doc.styles:
        if s.type == 1:                                         # paragraph styles
            s.font.size = Pt(STYLE_PT.get(s.name, TEXT_PT))
            if s.name in STYLE_PT:                              # the headings: bold, in black
                s.font.bold = True
                s.font.color.rgb = RGBColor(0, 0, 0)
    for r in doc.element.body.iter(qn("rPr")):
        for tag in ("rFonts", "sz", "szCs"):
            for x in r.findall(qn(tag)):
                r.remove(x)


# ------------------------------------------------------------------ what the layout gives
def layout():
    """The layout's configuration and its module, read from phase1/ (the paths in the
    configuration are from there)."""
    sys.path.insert(0, str(PHASE1))
    os.chdir(PHASE1)
    from src import evid_layout as L
    return L, L.load_config()


def figures_and_tables(L, cfg) -> tuple[list[dict], list[dict], list[dict]]:
    figures = [it for it in cfg["supplementary"] if it["kind"] == "figure"]
    tables = [it for it in cfg["supplementary"] if it["kind"] == "table"]
    return figures, [t for t in tables if t["destination"] == "pdf"], [t for t in tables if t["destination"] == "excel"]


def figure_pdf(it: dict) -> Path:
    return PHASE1 / "reports" / "evidence" / "layout" / "supplementary" / f"{it['content']}.pdf"


def figure_size_mm(path: Path) -> tuple[float, float]:
    from pypdf import PdfReader
    box = PdfReader(path).pages[0].mediabox
    return float(box.width) / 72 * 25.4, float(box.height) / 72 * 25.4


def caption(L, it: dict, numbers: dict) -> tuple[str, str]:
    """A figure's caption below it: its number, and its title and legend as printed."""
    return (f"{it['number']}.",
            f" {L.grouped(L._resolve(it['title'], numbers))} {L.grouped(L._resolve(it['legend'], numbers))}")


# ------------------------------------------------------------------ the document
def section(landscape: bool = False, side_mm: float = TEXT_MARGIN_MM) -> str:
    w, h = (A4[1], A4[0]) if landscape else A4
    return (f'<w:sectPr><w:type w:val="nextPage"/><w:pgSz w:w="{twips(w)}" w:h="{twips(h)}"'
            + (' w:orient="landscape"' if landscape else "")
            + f'/><w:pgMar w:top="{twips(TEXT_MARGIN_MM)}" w:right="{twips(side_mm)}" '
              f'w:bottom="{twips(TEXT_MARGIN_MM)}" w:left="{twips(side_mm)}" w:header="{twips(FOOTER_MM)}" '
              f'w:footer="{twips(FOOTER_MM)}" w:gutter="0"/><w:cols w:space="708"/></w:sectPr>')


def ends_section(p, landscape: bool = False, side_mm: float = TEXT_MARGIN_MM) -> None:
    """A paragraph made the last of its section, the next section starting on a new page."""
    ppr = p.find(qn("pPr"))
    if ppr is None:
        ppr = element("<w:pPr/>")
        p.insert(0, ppr)
    ppr.append(element(section(landscape, side_mm)))


def para(L, runs: str, before: float = 0, after: float = 0, keep_next: bool = False,
         keep_lines: bool = False, tabs: str = "") -> str:
    ppr = (("<w:keepNext/>" if keep_next else "") + ("<w:keepLines/>" if keep_lines else "") + tabs
           + f'<w:spacing w:before="{round(before * 20)}" w:after="{round(after * 20)}"/>')
    return f"<w:p><w:pPr>{ppr}</w:pPr>{runs}</w:p>"


def front_texts(L, numbers, headings, figures, pdf_tables, excel_tables) -> list[tuple[str, str]]:
    """The contents' entries, each its label (bold) and the rest, then the sentence that
    names the tables in the Excel file (its label empty)."""
    out = [("Abbreviations and conventions", "")]
    for h in headings:
        num, _, rest = h.partition(". ")
        out.append((f"{num}.", f" {rest}"))
    for it in figures + pdf_tables:
        out.append((f"{it['number']}.", f" {L.grouped(L._resolve(it['title'], numbers))}"))
    named = L._resolve("{" + ", ".join(it["content"] for it in excel_tables) + "}", numbers)
    out.append(("", f"The other {len(excel_tables)} supplementary tables, {named}, are in the accompanying Excel file."))
    return out


def front_page(L, cfg, numbers, headings, figures, pdf_tables, excel_tables, pages) -> list[str]:
    """The front page: the heading, the paper, the author and the contents, each entry's
    page number at the right after a dotted leader (two digits where it is not known yet,
    so that the first typesetting lays the page out as the second), and the sentence that
    names the tables in the Excel file; a range in an entry or the sentence is not broken
    at its dash (unbroken())."""
    run = lambda text, pt=TEXT_PT, bold=False: L._w_runs(unbroken(text), pt, FONT, bold=bold)   # noqa: E731
    width = twips(A4[0] - 2 * TEXT_MARGIN_MM)
    tabs = f'<w:tabs><w:tab w:val="right" w:leader="dot" w:pos="{width}"/></w:tabs>'
    out = [para(L, run(HEADING, 18, bold=True), after=12),
           para(L, run(PAPER, 12, bold=True), after=12),
           para(L, run(AUTHOR), after=24),
           para(L, run("Contents", 11, bold=True), after=6, keep_next=True)]

    def entry(key: str, label: str, rest: str) -> str:
        page = pages.get(key, "00")
        return para(L, run(label, bold=True) + run(rest) + f"<w:r><w:tab/></w:r>" + run(str(page)),
                    after=3, tabs=tabs)
    *entries, (_, sentence) = front_texts(L, numbers, headings, figures, pdf_tables, excel_tables)
    keys = (["conventions"] + [h.partition(". ")[0] for h in headings]
            + [it["number"] for it in figures + pdf_tables])
    out += [entry(key, label, rest) for key, (label, rest) in zip(keys, entries)]
    out.append(para(L, run(sentence), before=6))
    return out


def table_widths(L, tables: list[list[list[str]]]) -> tuple[bool, list[list[float]]]:
    """A table's page and its columns: upright if every part's columns fit the upright
    text width as wide as their text runs (evid_layout.column_widths), else turned. A part
    that fits shares the width in proportion to its columns' natural widths; one that does
    not keeps its narrower columns whole, and its wider ones, which wrap, share what is
    left, each at least as wide as its longest word. The run stops if a part is wider
    than the turned page at its longest words."""
    upright, turned = A4[0] - 2 * TEXT_MARGIN_MM, A4[1] - 2 * TEXT_MARGIN_MM
    sizes = [L.column_widths(t, TABLE_PT) for t in tables]
    landscape = any(sum(n) > upright for n, _ in sizes)
    width = turned if landscape else upright
    out = []
    for natural, least in sizes:
        if sum(least) > width + 1e-9:
            raise SystemExit(f"a table needs {sum(least):.0f} mm at its longest words, more than the "
                             f"{width:.0f} mm of the turned page")
        if sum(natural) <= width:
            out.append([n * width / sum(natural) for n in natural])
            continue
        lo, hi = 0.0, max(natural)
        for _ in range(60):                       # the widest cap at which the columns fit
            cap = (lo + hi) / 2
            if sum(max(m, min(n, cap)) for n, m in zip(natural, least)) <= width:
                lo = cap
            else:
                hi = cap
        w = [max(m, min(n, lo)) for n, m in zip(natural, least)]
        capped = [n > lo for n in natural]
        spare = width - sum(w)
        share = sum(x for x, c in zip(w, capped) if c) or 1.0
        out.append([x + spare * x / share if c else x for x, c in zip(w, capped)])
    return landscape, out


def top_aligned(table_xml: str) -> str:
    """A table's body cells set at the top of their rows, its header at the foot."""
    head, *rows = table_xml.split("<w:tr>")
    return "<w:tr>".join([head, rows[0]] + [r.replace('<w:vAlign w:val="bottom"/>', '<w:vAlign w:val="top"/>')
                                           for r in rows[1:]])


def table_blocks(L, cfg, numbers, it) -> tuple[bool, list[str], int]:
    """A table marked pdf: its title, its notes, and each part's title and table, the cells
    read from the CSVs the layout wrote, their gene symbols in italics
    (evid_layout.italic_genes), ruled across only (evid_layout._w_table); turned where a
    part is wider than the upright page. Returns whether turned, the blocks, and how many
    cells set a gene symbol in italics."""
    folder = PHASE1 / "reports" / "evidence" / "layout" / "supplementary"
    written = [L.written_rows(L._part_path(folder, it, p)) for p in it["parts"]]
    rows = [[[L.italic_genes(c, cfg["gene_symbols"]) for c in r] for r in t] for t in written]
    italic = sum(a != b for t, u in zip(written, rows) for r, q in zip(t, u) for a, b in zip(r, q))
    pdf_cfg = dict(cfg)
    pdf_cfg["docx"] = dict(cfg["docx"], margin_cm=TEXT_MARGIN_MM / 10, table_pt=TABLE_PT)
    landscape, widths = table_widths(L, rows)
    run = lambda text, bold=False: L._w_runs(text, TEXT_PT, FONT, bold=bold)          # noqa: E731
    out = [para(L, run(f"{it['number']}.", bold=True) + run(f" {L.grouped(L._resolve(it['title'], numbers))}"),
                after=6, keep_next=True)]
    out += [para(L, run(n), after=3, keep_next=True) for n in L._notes(it, numbers)]
    for k, (r, w) in enumerate(zip(rows, widths)):
        if len(it["parts"]) > 1:
            out.append(para(L, run(L.part_title(it, k, numbers)), before=6, after=3, keep_next=True))
        out.append(top_aligned(L._w_table(r, w, pdf_cfg)))
    return landscape, out, italic


def assemble(manuscript: Path, dest: Path, L, cfg, pages: dict, marks: bool) -> dict:
    """The document LibreOffice typesets: the manuscript's copy cut to its notes, the
    front page and the conventions before them, the figures' pages and the tables after.
    Returns what was laid out where it is known before typesetting."""
    import docx
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    numbers = L._numbers(cfg)
    doc = docx.Document(str(manuscript))
    body = doc.element.body
    start, end = notes_range(body)
    headings = note_headings(body, (start, end))
    kids = list(body)
    for el in kids[:start] + kids[end + 1:]:
        body.remove(el)                                   # the manuscript's own section too
    set_fonts(doc)
    figures, pdf_tables, excel_tables = figures_and_tables(L, cfg)
    first = list(body)[0]
    blocks = front_page(L, cfg, numbers, headings, figures, pdf_tables, excel_tables, pages)
    blocks += ['<w:p><w:r><w:br w:type="page"/></w:r></w:p>',
               para(L, L._w_runs(L.grouped(L._resolve(cfg["supplementary_format"]["conventions"], numbers)),
                                 TEXT_PT, FONT), after=12)]
    for xml in blocks:
        first.addprevious(element(xml))
    ends_section(list(body)[-1])                          # the notes' last entry ends the text pages
    laid = {"headings": headings, "figures": [], "tables": [],
            "front": front_texts(L, numbers, headings, figures, pdf_tables, excel_tables),
            "front_joiners": sum(x.count(JOINER) for x in blocks[:len(blocks) - 2])}
    for k, it in enumerate(figures):
        w_mm, h_mm = figure_size_mm(figure_pdf(it))
        if w_mm > A4[0] - 2 * FIGURE_SIDE_MM + 1e-6 or h_mm > A4[1] - 2 * TEXT_MARGIN_MM + 1e-6:
            raise SystemExit(f"{it['number']} is {w_mm:.0f} x {h_mm:.0f} mm, larger than its page's "
                             "space: it would have to be scaled")
        mark = MARK.format(k + 1) if marks else " "
        space = (f'<w:p><w:pPr><w:spacing w:before="0" w:after="0" w:line="{round(h_mm / 25.4 * 1440)}" '
                 f'w:lineRule="exact"/></w:pPr><w:r><w:rPr><w:color w:val="FFFFFF"/><w:sz w:val="2"/></w:rPr>'
                 f'<w:t xml:space="preserve">{mark}</w:t></w:r></w:p>')
        label, rest = caption(L, it, numbers)
        cap = element(para(L, L._w_runs(label, TEXT_PT, FONT, bold=True) + L._w_runs(rest, TEXT_PT, FONT),
                           before=8, keep_lines=True))
        body.append(element(space))
        body.append(cap)
        laid["figures"].append({"number": it["number"], "pdf": str(figure_pdf(it).relative_to(REPO)),
                                "width_mm": round(w_mm, 2), "height_mm": round(h_mm, 2)})
        ends_section(cap, side_mm=FIGURE_SIDE_MM)          # the figure's page: side margins of 1.5 cm
    for k, it in enumerate(pdf_tables):
        landscape, xml, italic = table_blocks(L, cfg, numbers, it)
        for x in xml:
            body.append(element(x))
        last = element('<w:p><w:pPr><w:spacing w:before="0" w:after="0"/></w:pPr></w:p>')
        body.append(last)
        laid["tables"].append({"number": it["number"], "turned": landscape, "italic_cells": italic})
        if k < len(pdf_tables) - 1:
            ends_section(last, landscape)
        else:
            body.append(element(section(landscape)))          # the document's last section
    sec = doc.sections[0]
    sec.footer.is_linked_to_previous = False
    p = sec.footer.paragraphs[0]
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p._p.append(element(f'<w:fldSimple w:instr="PAGE"><w:r><w:rPr><w:rFonts w:ascii="{FONT}" w:hAnsi="{FONT}" '
                        f'w:cs="{FONT}"/><w:sz w:val="18"/></w:rPr><w:t>1</w:t></w:r></w:fldSimple>'))
    doc.core_properties.title = HEADING
    doc.save(str(dest))
    return laid


# ------------------------------------------------------------------ typesetting, checking
def typeset(docx_path: Path, work: Path, soffice: str) -> Path:
    env = dict(os.environ)
    conf = work / "fonts.conf"
    conf.write_text(FONTS_CONF.format(cache=work / "fontcache"))
    env["FONTCONFIG_FILE"] = str(conf)
    profile = (work / "profile").as_uri()
    r = subprocess.run([soffice, "--headless", "--norestore", f"-env:UserInstallation={profile}",
                        "--convert-to", "pdf", "--outdir", str(work), str(docx_path)],
                       env=env, capture_output=True, text=True, timeout=600)
    pdf = work / (docx_path.stem + ".pdf")
    if r.returncode or not pdf.is_file():
        raise SystemExit(f"LibreOffice did not typeset {docx_path.name}: {r.stdout} {r.stderr}")
    return pdf


def page_texts(pdf: Path) -> list[str]:
    """Each page's text, the foot of the page (its number) left out."""
    from pypdf import PdfReader
    out = []
    for i, page in enumerate(PdfReader(pdf).pages, start=1):
        w, h = float(page.mediabox.width), float(page.mediabox.height)
        foot = (TEXT_MARGIN_MM - 2) / 25.4 * 72
        r = subprocess.run(["pdftotext", "-raw", "-f", str(i), "-l", str(i), "-x", "0", "-y", "0",
                            "-W", str(int(w)), "-H", str(int(h - foot)), str(pdf), "-"],
                           capture_output=True, text=True, check=True)
        out.append(r.stdout)
    return out


def squeeze(text: str) -> str:
    return " ".join(text.split())


def manuscript_notes_text(manuscript: Path, work: Path, soffice: str) -> str:
    """The notes' text as LibreOffice reads the manuscript (its list numbers included),
    from "Supplementary Notes" to the end of the last reference."""
    env = dict(os.environ)
    env["FONTCONFIG_FILE"] = str(work / "fonts.conf")
    copy_ = work / "manuscript_copy.docx"
    shutil.copyfile(manuscript, copy_)
    subprocess.run([soffice, "--headless", "--norestore", f"-env:UserInstallation={(work / 'profile').as_uri()}",
                    "--convert-to", "txt:Text (encoded):UTF8", "--outdir", str(work), str(copy_)],
                   env=env, capture_output=True, text=True, timeout=600, check=True)
    text = (work / "manuscript_copy.txt").read_text(encoding="utf-8-sig")
    lines = text.splitlines()
    first = [k for k, ln in enumerate(lines) if ln.strip() == NOTES_FROM]
    if len(first) != 1:
        raise SystemExit(f"LibreOffice's text of the manuscript has {len(first)} lines '{NOTES_FROM}'")
    return squeeze("\n".join(lines[first[0]:]))


def images_in(obj, seen=None) -> int:
    """How many images a page's resources hold, its forms' resources included."""
    seen = seen if seen is not None else set()
    n = 0
    res = obj.get("/Resources")
    xo = res.get_object().get("/XObject") if res is not None else None
    if xo is None:
        return 0
    for name in xo.get_object():
        x = xo.get_object()[name].get_object()
        if id(x) in seen:
            continue
        seen.add(id(x))
        if x.get("/Subtype") == "/Image":
            n += 1
        elif x.get("/Subtype") == "/Form":
            n += images_in(x, seen)
    return n


def fonts_of(pdf: Path) -> list[dict]:
    r = subprocess.run(["pdffonts", str(pdf)], capture_output=True, text=True, check=True)
    rows = r.stdout.splitlines()[2:]
    out = []
    for row in rows:
        parts = row.split()
        # name, type (one or two words), encoding, emb, sub, uni, object, id
        emb, sub, uni = parts[-5:-2]
        out.append({"name": parts[0], "embedded": emb == "yes", "subset": sub == "yes"})
    return out


def find_pages(texts: list[str], after: int, wanted: dict[str, str]) -> dict[str, list[int]]:
    """For each key, the pages (from 1) after page `after` whose text holds its string."""
    return {k: [i + 1 for i, t in enumerate(texts) if i + 1 > after and s in squeeze(t)]
            for k, s in wanted.items()}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manuscript", required=True, type=Path)
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--record", type=Path, default=None,
                    help="where the check record goes (default: beside the PDF, <name>.check.json)")
    ap.add_argument("--soffice", default=shutil.which("soffice") or "soffice")
    args = ap.parse_args()
    manuscript, out = args.manuscript.expanduser().resolve(), args.out.expanduser().resolve()
    record_path = (args.record.expanduser().resolve() if args.record else out.with_suffix(".check.json"))
    from pypdf import PdfReader, PdfWriter, Transformation
    L, cfg = layout()
    numbers = L._numbers(cfg)
    figures, pdf_tables, excel_tables = figures_and_tables(L, cfg)
    failed, record = [], {"manuscript": str(manuscript),
                          "manuscript_sha256": hashlib.sha256(manuscript.read_bytes()).hexdigest(),
                          "layout_commit": subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"],
                                                          capture_output=True, text=True).stdout.strip()}

    def check(name, ok, value=None):
        record.setdefault("checks", {})[name] = {"ok": bool(ok), "value": value}
        print(f"  {'ok  ' if ok else 'FAIL'} {name}" + ("" if value is None else f": {str(value)[:300]}"))
        if not ok:
            failed.append(name)

    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        # first typesetting: a mark in each figure's space, the contents' numbers unknown
        laid = assemble(manuscript, work / "first.docx", L, cfg, {}, marks=True)
        texts = page_texts(typeset(work / "first.docx", work, args.soffice))
        conv = squeeze(L._resolve(cfg["supplementary_format"]["conventions"], numbers))[:60]
        conv_page = [i + 1 for i, t in enumerate(texts) if conv in squeeze(t)]
        front = conv_page[0] - 1 if conv_page else 0
        wanted = {"conventions": conv}
        wanted |= {h.partition(". ")[0]: h for h in laid["headings"]}
        wanted |= {f"mark {k + 1}": MARK.format(k + 1) for k in range(len(figures))}
        wanted |= {f"caption {it['number']}": f"{it['number']}. {caption(L, it, numbers)[1].strip()[:40]}"
                   for it in figures}
        wanted |= {it["number"]: f"{it['number']}. {L.grouped(L._resolve(it['title'], numbers))[:40]}"
                   for it in pdf_tables}
        found = find_pages(texts, front, wanted)
        check("first_typesetting_finds_each_part_once", all(len(v) == 1 for v in found.values()),
              {k: v for k, v in found.items() if len(v) != 1})
        if failed:
            raise SystemExit(f"stopped: {failed}")
        pages = {"conventions": found["conventions"][0]}
        pages |= {h.partition(". ")[0]: found[h.partition(". ")[0]][0] for h in laid["headings"]}
        pages |= {it["number"]: found[f"mark {k + 1}"][0] for k, it in enumerate(figures)}
        pages |= {it["number"]: found[it["number"]][0] for it in pdf_tables}
        # the second typesetting, the contents' numbers filled in and no mark
        assemble(manuscript, work / "Supplementary_Information.docx", L, cfg, pages, marks=False)
        second = typeset(work / "Supplementary_Information.docx", work, args.soffice)
        texts2 = page_texts(second)
        found2 = find_pages(texts2, front, {k: v for k, v in wanted.items() if not k.startswith("mark")})
        same = all(found2[k] == found[k] for k in found2)
        check("second_typesetting_lays_out_as_the_first", same and len(texts2) == len(texts),
              {"pages": [len(texts), len(texts2)]})
        # the figures, laid on their spaces, moved and not scaled
        reader = PdfReader(second)
        writer = PdfWriter(clone_from=reader)
        placed = []
        for k, it in enumerate(figures):
            p = pages[it["number"]]
            page = writer.pages[p - 1]
            fig = PdfReader(figure_pdf(it)).pages[0]
            pw, ph = float(page.mediabox.width), float(page.mediabox.height)
            fw, fh = float(fig.mediabox.width), float(fig.mediabox.height)
            x, y = (pw - fw) / 2, ph - TEXT_MARGIN_MM / 25.4 * 72 - fh
            page.merge_transformed_page(fig, Transformation().translate(x, y))
            cap_page = found2[f"caption {it['number']}"][0]
            placed.append({"number": it["number"], "page": p, "caption_page": cap_page, "scale": 1.0,
                           "x_mm": round(x / 72 * 25.4, 2), "top_mm": TEXT_MARGIN_MM,
                           "width_mm": round(fw / 72 * 25.4, 2), "height_mm": round(fh / 72 * 25.4, 2),
                           "within_page": x >= FIGURE_SIDE_MM / 25.4 * 72 - 1e-6 and y >= TEXT_MARGIN_MM / 25.4 * 72 - 1e-6,
                           "images": images_in(page)})
        writer.add_metadata({"/Title": HEADING, "/Author": AUTHOR})
        out.parent.mkdir(parents=True, exist_ok=True)
        with open(out, "wb") as fh_:
            writer.write(fh_)
        final_texts = page_texts(out)
        # the front page: its text, the joiners and white space aside, in order, and no line
        # ending at a range's dash
        bare = lambda t: "".join(t.replace(JOINER, "").split())             # noqa: E731
        front_text = "\n".join(final_texts[:front])
        at, missing = 0, []
        for piece in [HEADING, PAPER, AUTHOR, "Contents"] + [a + b for a, b in laid["front"]]:
            k = bare(front_text).find(bare(L._plain(piece)), at)          # an italic run without its marks
            if k < 0:
                missing.append(piece[:80])
            else:
                at = k + len(bare(L._plain(piece)))
        check("front_page_prints_its_text_in_order_but_for_the_joiners", not missing, missing)
        broken = [m.group(0) for m in RANGE_BROKEN.finditer(front_text)]
        check("front_page_breaks_no_range_at_its_dash", not broken, broken)
        record["front_page_joiners"] = laid["front_joiners"]
        record["front_page_joiners_in_the_pdf_text"] = front_text.count(JOINER)
        record["pages"] = len(final_texts)
        record["bytes"] = out.stat().st_size
        record["starts"] = pages | {"notes": found["conventions"][0], "front": 1}
        record["figures"] = placed
        record["tables"] = laid["tables"]
        # checks on the PDF as written
        check("figures_on_their_own_pages_caption_below_or_next",
              all(f["caption_page"] in (f["page"], f["page"] + 1) for f in placed)
              and len({f["page"] for f in placed}) == len(placed), [(f["number"], f["page"], f["caption_page"]) for f in placed])
        check("figures_laid_at_scale_1_within_their_pages", all(f["scale"] == 1.0 and f["within_page"] for f in placed))
        check("figure_pages_vector_only", all(f["images"] == 0 for f in placed), [f["images"] for f in placed])
        order = [found2[f"caption {it['number']}"][0] for it in figures] + [found2[it["number"]][0] for it in pdf_tables]
        check("five_figures_and_the_pdf_tables_once_in_order", order == sorted(order) and
              all(len(found2[k]) == 1 for k in found2), order)
        heads = [found2[h.partition(". ")[0]][0] for h in laid["headings"]]
        check("notes_headings_in_order_before_the_figures", heads == sorted(heads) and heads[-1] <= order[0], heads)
        fonts = fonts_of(out)
        check("every_font_embedded", fonts and all(f["embedded"] for f in fonts),
              sorted({f["name"] for f in fonts}))
        record["fonts"] = sorted({f["name"] for f in fonts})
        # the notes' text: the PDF's from the notes' heading to the last reference
        want = manuscript_notes_text(manuscript, work, args.soffice)
        notes_pages = range(found["conventions"][0], pages[figures[0]["number"]])
        got = squeeze("\n".join(final_texts[p - 1] for p in notes_pages))
        got = got[got.find(NOTES_FROM):] if NOTES_FROM in got else got
        # white space aside: a line broken after a hyphen ("splice-|region") reads back with
        # a space the text does not have, and nothing else may differ
        bare_got, bare_want = "".join(got.split()), "".join(want.split())
        same = bare_got == bare_want
        diff = None
        if not same:
            k = next((i for i, (a, b) in enumerate(zip(bare_got, bare_want)) if a != b),
                     min(len(bare_got), len(bare_want)))
            diff = {"at": k, "pdf": bare_got[max(0, k - 80):k + 80], "manuscript": bare_want[max(0, k - 80):k + 80],
                    "lengths": [len(bare_got), len(bare_want)]}
        check("notes_text_is_the_manuscripts_but_for_white_space", same, diff)
        record["notes_characters"] = len(want)
        record["notes_characters_without_white_space"] = len(bare_want)
        record["notes_words_pdf_and_manuscript"] = [len(got.split()), len(want.split())]
        check("within_the_size_limit", record["bytes"] < LIMIT_MB * 1024 * 1024,
              f"{record['bytes'] / 1024 / 1024:.2f} MB")
    record["failed"] = failed
    record_path.parent.mkdir(parents=True, exist_ok=True)
    record_path.write_text(json.dumps(record, indent=1, ensure_ascii=False))
    print(f"{out}: {record['pages']} pages, {record['bytes']:,} bytes; {len(failed)} check(s) failed")
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
