"""scripts/stage_E8/build_supplementary_pdf.py (E8.12): the supplementary information as one
PDF, from the manuscript (read only) and the layout. The PDF itself holds the manuscript's
notes and is not committed; these tests build the document it is typeset from, on a
manuscript made here, with the layout's own outputs, and check what goes where. The whole
build, which needs LibreOffice and the manuscript, runs only where E8_MANUSCRIPT names the
manuscript and soffice is installed."""
from __future__ import annotations

import importlib.util
import json
import os
import re
import shutil
import zipfile
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
PHASE1 = REPO / "phase1"
SCRIPT = REPO / "scripts" / "stage_E8" / "build_supplementary_pdf.py"


@pytest.fixture(scope="module")
def B():
    spec = importlib.util.spec_from_file_location("build_supplementary_pdf", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


@pytest.fixture()
def layout(B, monkeypatch):
    monkeypatch.chdir(PHASE1)
    monkeypatch.syspath_prepend(str(PHASE1))
    from src import evid_layout as L
    return L, L.load_config()


def _manuscript(path: Path, extra_after: bool = False) -> Path:
    """A small manuscript as the real one is built: a main text, then the notes under the
    heading "Supplementary Notes", numbered items, a run-in heading and italics, and the
    notes' references."""
    import docx
    d = docx.Document()
    d.add_paragraph("Introduction", style="Heading 2")
    d.add_paragraph("Main text that the PDF leaves out.")
    d.add_paragraph("Supplementary Notes", style="Heading 2")
    d.add_paragraph("Supplementary Note 1. Predictor scores", style="Heading 3")
    p = d.add_paragraph("1. SpliceAI was run on ", style="List Number")
    p.add_run("BRCA1").italic = True
    d.add_paragraph("Supplementary Note 2. External genes", style="Heading 3")
    d.add_paragraph().add_run("At the published cut points.").bold = True
    d.add_paragraph("References cited in the Supplementary Notes", style="Heading 3")
    d.add_paragraph("21. Zhang N. A reference.")
    d.add_paragraph("44. Sahu S. Another reference.")
    d.add_paragraph("")
    if extra_after:
        d.add_paragraph("Appendix", style="Heading 2")
        d.add_paragraph("Not a note.")
    d.save(str(path))
    return path


def test_the_notes_run_from_their_heading_to_the_last_reference(B, tmp_path):
    import docx
    body = docx.Document(str(_manuscript(tmp_path / "m.docx", extra_after=True))).element.body
    start, end = B.notes_range(body)
    texts = [B.text_of(el) for el in list(body)[start:end + 1]]
    assert texts[0] == "Supplementary Notes" and texts[-1] == "44. Sahu S. Another reference."
    assert "Main text that the PDF leaves out." not in texts and "Not a note." not in texts
    assert B.note_headings(body, (start, end)) == ["Supplementary Note 1. Predictor scores",
                                                   "Supplementary Note 2. External genes"]
    d = docx.Document()
    d.add_paragraph("Supplementary Notes")
    with pytest.raises(SystemExit, match="not one of each"):
        B.notes_range(d.element.body)


def test_a_wide_table_keeps_its_short_columns_whole_and_a_too_wide_one_stops(B, layout):
    L, cfg = layout
    upright, turned = 210 - 40, 297 - 40
    short = [["Gene", "N"], ["*RAD51C*", "1,152"]]
    landscape, widths = B.table_widths(L, [short])
    assert not landscape and abs(sum(widths[0]) - upright) < 1e-6
    wide = [["Gene set", "Note"] + [f"Value {k}" for k in range(8)],
            ["Seven genes", "a long text that has to wrap over several lines " * 4] + ["0.3215"] * 8]
    landscape, widths = B.table_widths(L, [wide])
    natural, least = L.column_widths(wide, B.TABLE_PT)
    assert landscape and abs(sum(widths[0]) - turned) < 1e-6
    assert all(abs(w - n) < 1e-6 or w > n for w, n in zip(widths[0][:1] + widths[0][2:], natural[:1] + natural[2:]))
    assert widths[0][1] < natural[1] and all(w >= m - 1e-9 for w, m in zip(widths[0], least))
    with pytest.raises(SystemExit, match="turned page"):
        B.table_widths(L, [[["H"] * 14, ["supercalifragilistic"] * 14]])


def test_a_tables_body_is_set_at_the_top_and_its_header_at_the_foot(B):
    xml = ('<w:tbl><w:tr><w:tc><w:tcPr><w:vAlign w:val="bottom"/></w:tcPr></w:tc></w:tr>'
           '<w:tr><w:tc><w:tcPr><w:vAlign w:val="bottom"/></w:tcPr></w:tc></w:tr></w:tbl>')
    out = B.top_aligned(xml)
    head, body = out.split("<w:tr>")[1], out.split("<w:tr>")[2]
    assert 'w:val="bottom"' in head and 'w:val="top"' in body and 'w:val="bottom"' not in body


def test_the_document_holds_the_parts_in_order_each_figure_its_own_page_at_its_size(B, layout, tmp_path):
    """The document typeset into the PDF: the front page and the conventions, the notes
    copied with their numbering and italics, then a page per figure, its space as tall as
    the figure and its caption (number, title, legend) below, kept whole; then the tables
    marked pdf, ruled across only; A4, margins of 2 cm, the figures' pages 1.5 cm at the
    sides; a page number at the foot; Arial."""
    import docx
    L, cfg = layout
    numbers = L._numbers(cfg)
    out = tmp_path / "assembled.docx"
    laid = B.assemble(_manuscript(tmp_path / "m.docx"), out, L, cfg, {}, marks=False)
    d = docx.Document(str(out))
    body = d.element.body
    paras = [B.text_of(p) for p in body.iter(B.qn("p"))]
    figures, pdf_tables, excel = B.figures_and_tables(L, cfg)
    assert [f["number"] for f in laid["figures"]] == [f"Supplementary Figure S{k}" for k in range(1, 6)]
    assert [t["number"] for t in laid["tables"]] == [t["number"] for t in pdf_tables] == [
        "Supplementary Table S5", "Supplementary Table S6", "Supplementary Table S12",
        "Supplementary Table S25", "Supplementary Table S28"]
    order = [paras.index(x) for x in (B.HEADING, B.PAPER, B.AUTHOR, "Supplementary Notes",
                                      "21. Zhang N. A reference.", "44. Sahu S. Another reference.")]
    assert order == sorted(order) and "Main text that the PDF leaves out." not in paras
    conv = L.grouped(L._resolve(cfg["supplementary_format"]["conventions"], numbers))
    assert paras.index(conv) < paras.index("Supplementary Notes")
    for it, f in zip(figures, laid["figures"]):
        label, rest = B.caption(L, it, numbers)
        # the document's text: an italic run without its marks
        assert L._plain(label + rest) in paras and rest.endswith(L.grouped(L._resolve(it["legend"], numbers)))
        assert f["width_mm"] == 175.0
    heights = [round(int(s.get(B.qn("line"))) / 1440 * 25.4, 1) for s in body.iter(B.qn("spacing"))
               if s.get(B.qn("lineRule")) == "exact"]
    assert heights == [round(f["height_mm"], 1) for f in laid["figures"]] and heights[0] == 240.0
    sects = [s for s in body.iter(B.qn("sectPr"))]
    assert len(sects) == 1 + len(figures) + len(pdf_tables)
    margins = [(int(s.find(B.qn("pgMar")).get(B.qn("left"))), s.find(B.qn("pgSz")).get(B.qn("orient")))
               for s in sects]
    side, text = B.twips(B.FIGURE_SIDE_MM), B.twips(B.TEXT_MARGIN_MM)
    assert [m for m, _ in margins] == [text] + [side] * 5 + [text] * 5
    assert [o == "landscape" for _, o in margins[6:]] == [t["turned"] for t in laid["tables"]]
    for tbl in body.iter(B.qn("tbl")):
        borders = tbl.find(f"{B.qn('tblPr')}/{B.qn('tblBorders')}")
        assert [borders.find(B.qn(s)).get(B.qn("val")) for s in ("top", "bottom", "left", "right", "insideH", "insideV")] \
            == ["single", "single", "nil", "nil", "nil", "nil"]
        assert int(borders.find(B.qn("top")).get(B.qn("sz"))) >= 8
    with zipfile.ZipFile(out) as z:
        footers = [z.read(n).decode() for n in z.namelist() if n.startswith("word/footer")]
        styles = z.read("word/styles.xml").decode()
    assert any('w:instr="PAGE"' in f for f in footers)
    assert set(re.findall(r'w:rFonts w:ascii="([^"]+)"', styles)) == {"Arial"}
    assert "21. Zhang N. A reference." in paras and any(
        r.italic and r.text == "BRCA1" for p in d.paragraphs for r in p.runs)


def test_a_figure_too_large_for_its_page_stops_the_run(B, layout, tmp_path, monkeypatch):
    L, cfg = layout
    monkeypatch.setattr(B, "figure_size_mm", lambda path: (180.5, 100.0))
    with pytest.raises(SystemExit, match="scaled"):
        B.assemble(_manuscript(tmp_path / "m.docx"), tmp_path / "a.docx", L, cfg, {}, marks=False)


@pytest.mark.skipif(not (os.environ.get("E8_MANUSCRIPT") and shutil.which("soffice") and shutil.which("pdftotext")),
                    reason="the whole build needs the manuscript (E8_MANUSCRIPT), LibreOffice and poppler")
def test_the_whole_build_passes_its_checks(B, tmp_path):
    """The script as it is run: every check passes and its record says where each part is."""
    import subprocess
    import sys
    out = tmp_path / "Supplementary_Information.pdf"
    r = subprocess.run([sys.executable, str(SCRIPT), "--manuscript", os.environ["E8_MANUSCRIPT"], "--out", str(out)],
                       capture_output=True, text=True, timeout=1800)
    assert r.returncode == 0, r.stdout[-2000:] + r.stderr[-2000:]
    rec = json.loads(out.with_suffix(".check.json").read_text())
    assert rec["failed"] == [] and all(c["ok"] for c in rec["checks"].values())
    assert all(f["scale"] == 1.0 and f["images"] == 0 for f in rec["figures"])
