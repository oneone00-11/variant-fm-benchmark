"""The layout layer (src/evid_layout.py, E5.3): numbers only in config/layout.yaml, files
named for their content, every table and figure from an output already written, the
printed tables as they stand, the new tables under readable names and the number formats
of the tables v3.0.1 printed, the figures within the artwork guide, and a place for every output added since
release 3.0.1. Figure 1 (c), E6: one axis, the curve Supplementary Figure S1 draws for the
same predictor and band, point for point, in greys, with (a) and (b) where they were."""
from __future__ import annotations

import hashlib
import re
import subprocess
import sys
import zipfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
import yaml

PHASE1 = Path(__file__).resolve().parents[1]
REPO = PHASE1.parent
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tests"))
from _checkout import needs_history  # noqa: E402  one judgement for the tests that read the history (E9)
if str(PHASE1) not in sys.path:
    sys.path.insert(0, str(PHASE1))
from src import evid_layout as L  # noqa: E402

CONFIG = PHASE1 / "config" / "layout.yaml"
OUT = PHASE1 / "reports" / "evidence" / "layout"
V301 = "8e0d609370565fe50933b2c9ec2f62325289e24d"     # release 3.0.1, the layout v5 used
NUMBERED = re.compile(r"(?i)\b(?:figure|table|fig|tab)\s*S?\d|(?:^|[_.])S\d+(?:[_.]|$)")
NOTE_NUMBERED = re.compile(r"(?i)\bnotes?\s+\d")
# a number of another publication that the approved texts name (E8.10): printed as it is
SOURCE_PAPERS = re.compile(r"Supplementary Table 5 of that publication")
# the author's approved texts of 2026-11-09 (E8.10), a copy of the file the instruction
# names, byte for byte
APPROVED = PHASE1 / "tests" / "fixtures" / "approved_texts_20261109.yaml"
APPROVED_SHA256 = "654a1cc95d01ec84f872182ca289b8bd4af44e450e0864c7112f888a7c31b2c8"
MM_PER_PT = 25.4 / 72

need_out = pytest.mark.skipif(not (OUT / "manifest.csv").exists(), reason="the layout stage has not run")


def _cfg() -> dict:
    return yaml.safe_load(CONFIG.read_text())


def _items() -> list[dict]:
    """Every numbered item, in contents.csv's order: the main ones, the supplementary
    ones, and the variant list (E8.5)."""
    c = _cfg()
    return c["main"] + c["supplementary"] + c.get("variants", [])


def _new_parts():
    """(item, part, its source's printed columns) for every part the stage prints itself:
    every column but one a part's `drop` leaves out as a copy of another, or its `omit`
    leaves out."""
    for it in _cfg()["supplementary"]:
        if it["make"] != "table_parts":
            continue
        for part in it["parts"]:
            srcs = part.get("sources") or [part["source"]]
            if isinstance(srcs, dict):
                cols = [part["key"]] + list(L.read_source(PHASE1 / next(iter(srcs.values())))[0].columns)
            else:
                if str(srcs[0]).startswith(L.PRINTED):
                    continue
                cols = list(L.read_source(PHASE1 / srcs[0])[0].columns)
            yield it, part, [c for c in cols if c not in part.get("drop", {})
                             and c not in part.get("omit", [])]


def test_content_names_and_numbers_are_unique_and_names_carry_no_number():
    items = _items()
    contents = [it["content"] for it in items]
    numbers = [it["number"] for it in items]
    assert len(set(contents)) == len(contents)
    assert len(set(numbers)) == len(numbers)
    for it in items:
        assert not NUMBERED.search(it["content"]), it["content"]
        for part in it.get("parts", []):
            assert not NUMBERED.search(part["part"]), part["part"]


def test_every_source_exists_and_every_maker_is_known():
    for it in _items():
        assert it["make"] in L.MAKERS, it["content"]
        for s in L._sources(it):
            assert (PHASE1 / s).is_file(), s
    for r in _cfg().get("repository_only", []):
        assert (PHASE1 / r["path"]).is_file(), r["path"]


def test_every_output_added_since_release_301_has_a_place():
    """Each output an earlier stage added after release 3.0.1 is a source of a numbered
    table or figure, or is listed as kept in the repository only (E5.0 §4.4)."""
    needs_history(V301)
    added = subprocess.run(["git", "-C", str(REPO), "diff", "--name-only", "--diff-filter=A", V301,
                            "HEAD", "--", "phase1/reports/evidence"],
                           capture_output=True, text=True, check=True).stdout.split()
    added = [a[len("phase1/"):] for a in added if not a.startswith("phase1/reports/evidence/layout/")]
    placed = {s for it in _items() for s in L._sources(it)} | {r["path"] for r in _cfg()["repository_only"]}
    assert added
    assert sorted(a for a in added if a not in placed) == []


def test_every_new_column_has_one_readable_name_and_a_known_type():
    """Every column of a table the stage prints has its entry, each entry names a column
    the table has, names are distinct within a part, and none names an item's number,
    which may change while the names stay."""
    cfg = _cfg()
    used = {k: set() for k in cfg["columns"]}
    for it, part, cols in _new_parts():
        names = []
        for c in cols:
            spec = L._spec(it, part, c, cfg)
            assert spec[1] in L.TYPES or spec[1] in cfg["display"]["values"], (it["content"], c, spec)
            assert not NUMBERED.search(spec[0]), (it["content"], c, spec[0])
            names.append(spec[0])
            key = next(k for k in (f"{it['content']}.{part['part']}", it["content"])
                       if c in cfg["columns"].get(k, {}))
            used[key].add(c)
        assert len(set(names)) == len(names), (it["content"], part["part"])
    for key, names in cfg["columns"].items():
        assert key.split(".")[0] in {it["content"] for it in cfg["supplementary"]}, key
        assert set(names) == used[key], (key, sorted(set(names) - used[key]))


@pytest.mark.parametrize("text, kind, row, want", [
    ("570.0", "count", {}, "570"),
    ("0.7947368421052632", "proportion", {}, "0.795"),
    ("0.5", "proportion", {}, "0.500"),
    ("13.363157894736842", "lr", {}, "13.4"),
    ("0.4807692307692308", "lr", {}, "0.481"),
    ("0.006782", "score", {}, "0.006782"),
    ("25.488039593071214", "percent", {}, "25.5"),
    ("0.22699386503067487", "number", {}, "0.227"),
    ("0.05", "constant", {}, "0.05"),
    ("0.2", "threshold", {"threshold_kind": "walker_pp3"}, "0.2"),
    ("0.366", "proportion", {}, "0.366"),
    ("True", "boolean", {}, "Yes"),
    ("False", "boolean", {}, "No"),
    ("spliceai_walker", "tool", {}, "SpliceAI (published basis)"),
    ("fusion_enet", "tool", {}, "Elastic-net combination"),
    ("s3_10", "stratum", {}, "3–10"),
    ("out_of_scope_deep", "stratum", {}, "Beyond 50 (out of scope)"),
    ("very_strong", "tier", {}, "Very strong"),
    ("below supporting", "tier", {}, "Below Supporting"),
    ("not evaluable", "tier", {}, "not evaluable"),
    ("bp4", "side", {}, "BP4"),
    ("walker_pp3", "threshold_kind", {}, "Published cut point (≥ 0.2)"),
    ("logo", "threshold_source", {}, "Held out (leave one gene out)"),
    ("stars_ge_1", "star_level", {}, "At least one star"),
    ("inf", "threshold", {}, "∞"),
    ("", "count", {}, ""),
    ("anything at all", "text", {}, "anything at all"),
])
def test_a_cell_prints_as_the_printed_tables_print_its_type(text, kind, row, want):
    values = _cfg()["display"]["values"]
    assert L._cell(text, kind, row, [], values, "test") == want


def test_a_threshold_keeps_the_digits_that_keep_its_band():
    """evid_supp_tables.fmt_thr, as thresholds_pathogenic prints a threshold: four significant figures, and
    more where four would move a variant across it (here, the side decides which way)."""
    values = _cfg()["display"]["values"]
    from src import evid_supp_tables as S
    x = "0.5108555555343628"
    assert L._cell(x, "threshold", {"tool": "pangolin", "side": "pp3"}, [], values, "t") == \
        S.fmt_thr(float(x), "pangolin", "upper")
    assert L._side({"side": "BP4"}) == "lower" and L._side({"threshold_kind": "bp4_moderate"}) == "lower"
    assert L._side({"side": "pp3", "threshold_kind": "strong"}) == "upper"


@pytest.mark.parametrize("text, kind", [("n/a", "count"), ("1.5", "count"), ("yes", "boolean"),
                                        ("donor", "side"), ("not_a_tool", "tool"),
                                        ("s9_99", "stratum"), ("not_a_kind", "threshold_kind")])
def test_a_cell_its_type_cannot_print_stops_the_run(text, kind):
    with pytest.raises(SystemExit):
        L._cell(text, kind, {}, [], _cfg()["display"]["values"], "test")


def test_supplementary_items_have_a_destination_and_a_short_description():
    cfg = _cfg()
    limit = cfg["supplementary_format"]["description_words"]
    assert limit == 50
    for it in cfg["supplementary"]:
        assert it["destination"] in ("pdf", "excel"), it["content"]
        assert 0 < len(it["description"].split()) <= limit, it["content"]
        if it["kind"] == "figure":
            assert it["destination"] == "pdf", it["content"]
    for f in cfg["supplementary_format"]["files"]:
        assert 0 < len(f["description"].split()) <= limit


def test_notes_name_other_items_by_content():
    cfg = _cfg()
    numbers = L._numbers(cfg)
    texts = [n for it in _items() for n in it.get("notes", [])]
    texts += [n for it in _items() for p in it.get("parts", []) for n in p.get("notes", [])]
    texts += [it["legend"] for it in _items() if "legend" in it]
    texts += [r["what"] for r in cfg["repository_only"]] + [cfg["supplementary_format"]["conventions"]]
    refs = [n.strip() for t in texts for m in L.REFERENCE.findall(t) for n in m.split(",")]
    assert refs and all(r in numbers for r in refs)
    assert {"local_lr_calibration", "external_gene_analyses"} <= set(refs)       # the notes'
    for t in texts:
        assert not NUMBERED.search(SOURCE_PAPERS.sub("", t)), t
        assert not NOTE_NUMBERED.search(t), t
    with pytest.raises(SystemExit):
        L._resolve("see {no_such_item}", numbers)


def test_text_outside_the_guide_is_brought_into_it_or_stops_the_run(tmp_path):
    cfg = _cfg()
    fig, ax = plt.subplots(figsize=(L._width(cfg), 1))
    big = ax.text(0, 0, "a panel letter", fontsize=9)
    L._save(fig, tmp_path, "big", cfg)
    assert big.get_fontsize() == cfg["artwork"]["text_pt"][1]
    fig, ax = plt.subplots(figsize=(L._width(cfg), 1))
    ax.text(0, 0, "too small", fontsize=4)
    with pytest.raises(SystemExit):
        L._save(fig, tmp_path, "small", cfg)
    fig, ax = plt.subplots(figsize=(L._width(cfg) * 1.05, 1))
    with pytest.raises(SystemExit):
        L._save(fig, tmp_path, "wide", cfg)


def test_redrawing_leaves_evid_figures_as_it_was(tmp_path):
    from src import evid_figures as F
    double, save = F.DOUBLE, F._save
    kept = []
    with pytest.raises(RuntimeError):
        with L._drawn_by_figures(1.0, kept):
            assert F.DOUBLE == 1.0
            raise RuntimeError
    assert (F.DOUBLE, F._save) == (double, save)


@need_out
def test_every_figure_is_within_the_artwork_width():
    cfg = _cfg()
    width = cfg["artwork"]["double_column_mm"]
    pdfs = sorted(OUT.rglob("*.pdf"))
    assert len(pdfs) == sum(it["kind"] == "figure" for it in _items())
    for pdf in pdfs:
        box = re.search(rb"/MediaBox \[\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*\]", pdf.read_bytes())
        assert float(box.group(3)) * MM_PER_PT <= width + 0.01, pdf.name
        png = pdf.with_suffix(".png").read_bytes()
        assert int.from_bytes(png[16:20], "big") / cfg["artwork"]["dpi"] * 25.4 <= width + 0.1, pdf.name


@need_out
def test_the_manifest_lists_every_file_written_with_its_sha256():
    m = pd.read_csv(OUT / "manifest.csv", dtype=str, keep_default_na=False)
    on_disk = sorted(p.relative_to(REPO).as_posix() for p in OUT.rglob("*")
                     if p.is_file() and p.name != "manifest.csv")
    assert sorted(m.file) == on_disk
    for r in m.itertuples():
        assert hashlib.sha256((REPO / r.file).read_bytes()).hexdigest() == r.sha256, r.file
    for name in on_disk:
        assert not NUMBERED.search(Path(name).name.split(".")[0]), name


@need_out
def test_the_contents_hold_every_item_and_the_repository_only_outputs():
    c = pd.read_csv(OUT / "contents.csv", dtype=str, keep_default_na=False)
    assert list(c[c.number != "repository only"].content) == [it["content"] for it in _items()]
    assert len(c[c.number == "repository only"]) == len(_cfg()["repository_only"])
    supp = c[(c.section == "supplementary") & (c.kind == "table")]
    assert supp.fits_one_page_estimate.str.match(r"^(yes|no) \((\d{1,3}(,\d{3})*) lines, (\d{1,3}(,\d{3})*) mm wide\)$").all()
    # the destination defaults to the estimate (the author's decision of 2026-10-28)
    assert ((supp.destination == "pdf") == supp.fits_one_page_estimate.str.startswith("yes")).all()
    assert not c.notes.str.contains(r"\{", regex=True).any()


@need_out
def test_the_external_genes_without_a_marker_are_named_for_the_caption():
    """SBDS has no damaging variant at 11-50 bp, so Figure 2 draws no marker for it there
    and the caption says so (E5.3); the stage writes the reason beside the figure."""
    c = pd.read_csv(OUT / "contents.csv", dtype=str, keep_default_na=False)
    notes = c[c.content == "heldout_and_external_lr"].notes.iloc[0]
    assert "no marker: SBDS at 11–50 bp, every row: fewer than 10 labelled on one side (0 damaging" in notes


@need_out
def test_the_printed_tables_are_written_byte_for_byte():
    """But for the folds a rule leaves out of the counts (thresholds_pathogenic, tested
    below), the texts a part prints in place of its source's (its `replace`, E8.1), the
    thousands grouped and the transcripts as RefSeq (E8.3): each cell of the copy is the
    source's with those texts replaced and its numbers grouped (ungrouped, it is the
    source's), or its transcript as the MANE summary gives it; and a copy whose source
    has no such text, number or transcript is the source byte for byte."""
    numbers = L._numbers(_cfg())
    replaced = untouched = 0
    for it in _items():
        for part in it.get("parts", []):
            src = part.get("source")
            if not src or not src.startswith(L.PRINTED) or "not_counted" in part:
                continue
            section = "main" if it in _cfg()["main"] else "supplementary"
            name = f"{it['content']}.csv" if len(it["parts"]) == 1 else f"{it['content']}.{part['part']}.csv"
            text = (PHASE1 / src).read_bytes().decode()
            for old_, new_, _ in L._replacements(part, numbers):
                text, replaced = text.replace(old_, new_), replaced + 1
            got = (OUT / section / name).read_bytes().decode()
            if got == text:
                untouched += 1
                continue
            a, b = L._parse_source(text), L._parse_source(got)
            assert [ungroup(n) for n in b[1]] == [ungroup(n) for n in a[1]], name
            assert [ungroup(c) for c in b[0].columns] == [ungroup(c) for c in a[0].columns], name
            assert a[0].shape == b[0].shape, name
            tx = part.get("transcripts", {}).get("column")
            for col in a[0].columns:
                for x, y in zip(a[0][col], b[0][col]):
                    if col == tx:
                        assert y == x or (y.endswith(f" ({x})") and y.startswith("NM_")), (name, x, y)
                    else:
                        assert ungroup(y) == ungroup(x) and y == L.grouped(x), (name, col, x, y)
    assert replaced == 2            # training_signal's fusion row, external_genes' five headers
    assert untouched >= 5           # Tables 1-3's parts, among others


def ungroup(text: str) -> str:
    """A text with the commas between digit groups taken out."""
    return re.sub(r"(?<=\d),(?=\d{3}\b)", "", text)


@need_out
def test_the_new_tables_keep_every_row_and_column_of_their_sources():
    for it, part, cols in _new_parts():
        name = f"{it['content']}.csv" if len(it["parts"]) == 1 else f"{it['content']}.{part['part']}.csv"
        printed = pd.read_csv(OUT / "supplementary" / name, dtype=str, keep_default_na=False)
        srcs = part.get("sources") or [part["source"]]
        srcs = list(srcs.values()) if isinstance(srcs, dict) else srcs
        rows = sum(len(L.read_source(PHASE1 / s)[0]) for s in srcs)
        assert printed.shape == (rows, len(cols) + len(part.get("marks", []))), name


@need_out
def test_the_practical_table_prints_the_summary_cells():
    """Table 4's first part gives, for its predictors, the held-out values of the
    practical-metrics summary's table (the part before the brackets), the same rule and
    display, and its fold count under a name of its own."""
    md = (REPO / "docs" / "practical-metrics-summary.md").read_text().splitlines()
    rows = {tuple(c.strip() for c in ln.strip("|").split("|")[:2]): [c.strip() for c in ln.strip("|").split("|")]
            for ln in md if ln.startswith("| SpliceAI (published basis) |") or ln.startswith("| Pangolin |")}
    item = next(it for it in _cfg()["main"] if it["make"] == "table_practical_performance")
    t = pd.read_csv(OUT / "main" / f"{item['content']}.{item['parts'][0]['part']}.csv", dtype=str,
                    keep_default_na=False)
    assert t.columns[-1] == "Folds with a fitted threshold (S, M, St)"
    assert len(t) == 6
    for r in t.itertuples(index=False):
        summary = rows[(r[0], r[1])]
        assert list(r)[2:11] == [cell.split(" (")[0] for cell in summary[2:11]], (r[0], r[1])
        assert r[11] == summary[11]


@need_out
def test_the_clinvar_part_leaves_out_the_band_it_cannot_evaluate():
    item = next(it for it in _cfg()["main"] if it["make"] == "table_practical_performance")
    t = pd.read_csv(OUT / "main" / f"{item['content']}.{item['parts'][1]['part']}.csv", dtype=str,
                    keep_default_na=False)
    assert sorted(set(t.Band)) == ["3–10 bp", "3–50 bp"]
    # the approved notes say so, with the band's P/LP count, which the stage checks (E8.11)
    c = pd.read_csv(OUT / "contents.csv", dtype=str, keep_default_na=False)
    assert "The distal band, with 3 P/LP variants, is left out." in \
        c[c.content == item["content"]].notes.iloc[0]


@need_out
def test_the_workbook_holds_every_supplementary_table_with_its_printed_decimals():
    tables = [it for it in _cfg()["supplementary"] if it["kind"] == "table"]
    with zipfile.ZipFile(OUT / "supplementary_tables.xlsx") as z:
        names = z.namelist()
        book = z.read("xl/workbook.xml").decode()
        styles = z.read("xl/styles.xml").decode()
        assert all(info.date_time == (1980, 1, 1, 0, 0, 0) for info in z.infolist())
        sheet = z.read("xl/worksheets/sheet2.xml").decode()
    sheets = re.findall(r'<sheet name="([^"]+)"', book)
    assert sheets == ["README"] + [L._sheet_name(it["number"]) for it in tables]
    assert sum(n.startswith("xl/worksheets/sheet") for n in names) == len(sheets)
    # a cell printed 0.500 is the number 0.5 shown with three decimals
    assert 'formatCode="0.000"' in styles
    cells = re.findall(r'<c r="[A-Z]+\d+" s="(\d+)"><v>([^<]+)</v></c>', sheet)
    assert cells
    for s, v in cells:
        # a number of 1,000 or more shown grouped in thousands (E8.3), with its decimals
        base = L.GROUPED_STYLE if abs(float(v)) >= 1000 else L.NUMBER_STYLE
        assert int(s) - base == len((v.split(".") + [""])[1]), (s, v)
    assert 'formatCode="#,##0"' in styles and 'formatCode="#,##0.000"' in styles


def test_a_share_stored_rounded_prints_at_its_own_decimals():
    """C12 stores its shares rounded to four decimals; printed at three they would be
    rounded twice (0.37746 -> 0.3775 -> 0.378, where inframe_attribution prints
    0.377), so the configuration names four."""
    values = _cfg()["display"]["values"]
    assert L._cell("0.366", "proportion", {}, [4], values, "t") == "0.3660"
    assert L._cell("0.3775", "proportion", {}, [4], values, "t") == "0.3775"
    spec = _cfg()["columns"]["inframe_true_positives"]["in_frame_share"]
    assert spec[1:] == ["proportion", 4]


def test_every_external_gene_has_a_shape_with_an_open_form():
    """A ratio below its boundary is drawn open; a cross has no open form, so TP53 is a
    hexagon (the author's review of 2026-10-28), still in the second colour."""
    assert set(L.GENE_MARKER.values()) <= L.OPEN_SHAPES
    assert len(set(L.GENE_MARKER.values())) == len(L.GENE_MARKER)
    assert L.GENE_MARKER["TP53"] == "h" and L.GENE_COLOUR["TP53"] != L.GENE_COLOUR["DDX3X"]
    assert set(L.MARKER_SIZE) == set(L.GENE_MARKER.values())


def test_a_figure_taller_than_the_limit_stops_the_run(tmp_path):
    cfg = _cfg()
    fig, ax = plt.subplots(figsize=(L._width(cfg), (cfg["artwork"]["max_height_mm"] + 1) / 25.4))
    with pytest.raises(SystemExit):
        L._save(fig, tmp_path, "tall", cfg)


@need_out
def test_every_figure_is_within_the_height_limit():
    """Supplementary Figure S1, 250 mm tall as evid_figures draws it, is drawn at 240 mm."""
    limit = _cfg()["artwork"]["max_height_mm"]
    for pdf in sorted(OUT.rglob("*.pdf")):
        box = re.search(rb"/MediaBox \[\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*\]", pdf.read_bytes())
        assert float(box.group(4)) * MM_PER_PT <= limit + 0.01, pdf.name
    box = re.search(rb"/MediaBox \[\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*\]",
                    (OUT / "supplementary" / "local_lr_curves.pdf").read_bytes())
    assert abs(float(box.group(4)) * MM_PER_PT - limit) < 0.01


def _part_of(content, part):
    it = next(i for i in _cfg()["supplementary"] if i["content"] == content)
    return it, next(p for p in it["parts"] if p["part"] == part)


def test_a_column_left_out_must_repeat_the_one_kept():
    cfg = {"display": {"values": {}}, "columns": {"t": {"a": ["A", "count"], "b": ["B", "count"]}}}
    part = {"part": "p", "drop": {"b": "a"}}
    same = pd.DataFrame({"a": ["1", "2"], "b": ["1", "2"]})
    assert list(L.displayed(same, {"content": "t"}, part, cfg).columns) == ["A"]
    with pytest.raises(SystemExit):
        L.displayed(same.assign(b=["1", "3"]), {"content": "t"}, part, cfg)


@need_out
def test_the_atlas_rows_held_out_on_its_selection_assays_are_marked():
    """family_holdout keeps the combined Atlas score's rows with BRCA1 or RAD51C
    held out and marks each, by evid_common's rule, in (b)'s wording for DDX3X."""
    from src import evid_common as K
    it, part = _part_of("family_holdout", "holdout")
    src = L.read_source(PHASE1 / part["source"])[0]
    printed = pd.read_csv(OUT / "supplementary" / "family_holdout.csv", dtype=str, keep_default_na=False)
    mark = part["marks"][0]
    want = (src.tool == "avi") & src.heldout_gene.isin(K.AVI_SEEN_IN_TRAINING)
    assert want.sum() == 12
    assert (printed[mark["name"]] == mark["text"]).tolist() == want.tolist()
    assert (printed[mark["name"]].isin([mark["text"], ""])).all()
    cols = list(printed.columns)
    assert cols.index(mark["name"]) == cols.index("Note") - 1
    assert mark["text"].startswith("not an independent test: the combined Atlas score's model was selected")


@need_out
@pytest.mark.parametrize("content, part, file, condition", [
    ("heldout_intervals", "intervals", "heldout_intervals.intervals.csv",
     lambda s: s.threshold_from_other_genes == ""),
    ("label_sources", "effect", "label_sources.effect.csv", lambda s: s.side == "bp4"),
])
def test_counts_without_a_definition_print_a_dash(content, part, file, condition):
    """A count the row has no definition for prints "–", not 0: the band counts where no
    threshold was fitted on the other genes (heldout_intervals), the fold counts on BP4
    rows (label_sources).
    The source tables are as they were."""
    it, p = _part_of(content, part)
    src = L.read_source(PHASE1 / p["source"])[0]
    printed = pd.read_csv(OUT / "supplementary" / file, dtype=str, keep_default_na=False)
    hit = condition(src)
    assert hit.any() and (~hit).any()
    for c in p["dashes"][0]["columns"]:
        name = L._spec(it, p, c, _cfg())[0]
        assert (printed.loc[hit, name] == "–").all(), (file, name)
        assert (printed.loc[~hit, name] != "–").all(), (file, name)
        assert (src.loc[hit, c] == "0").all()


# ---------------------------------------------------------------------------
# Figure 1 (c) (E6, the author's decision of 2026-11-04)
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def figure_one(tmp_path_factory):
    """Figure 1 and Supplementary Figure S1 as the stage draws them, into a temporary
    folder, each kept after the stage's own _save."""
    import matplotlib
    out = tmp_path_factory.mktemp("layout")
    kept, save = [], L._save
    mp = pytest.MonkeyPatch()
    mp.chdir(PHASE1)                         # the makers read reports/evidence/...
    mp.setattr(L, "_save", lambda fig, folder, stem, cfg: (kept.append(fig), save(fig, folder, stem, cfg))[1])
    try:
        with matplotlib.rc_context():
            cfg = L.load_config()
            item = next(i for i in cfg["main"] if i["content"] == "published_cut_point_by_band")
            s1 = next(i for i in cfg["supplementary"] if i["content"] == "local_lr_curves")
            L.figure_published_cut_point(item, out, cfg)
            L.redraw_figure(s1, out, cfg)
    finally:
        mp.undo()
    return {"cfg": cfg, "item": item, "fig": kept[0], "s1": kept[1]}


def _curve_drawn(ax):
    """A panel's curve as drawn (score; local likelihood ratio, clipped to the axis) and
    its shade's polygon."""
    from matplotlib.collections import PolyCollection
    line = [ln for ln in ax.lines if len(ln.get_xdata()) > 2]
    shade = [c for c in ax.collections if isinstance(c, PolyCollection)]
    assert len(line) == 1 and len(shade) == 1 and len(shade[0].get_paths()) == 1
    return (np.asarray(line[0].get_xdata(), dtype=float), np.asarray(line[0].get_ydata(), dtype=float),
            shade[0].get_paths()[0].vertices)


def test_figure_one_has_three_axes_and_c_draws_s1s_curve_point_for_point(figure_one):
    """(c) is one axis, so that Figure 1 has three (each axis counts as a figure, and a
    composite of more than three as two), and it draws S1's curve of the same predictor
    and band, every row of its file, with S1's clipping and S1's shade."""
    from src import evid_figures as F
    fig, s1, opts = figure_one["fig"], figure_one["s1"], figure_one["item"]["options"]
    assert len(fig.axes) == 3
    tool, st = opts["c_tool"], opts["c_stratum"]
    panel = s1.axes[F.COLUMNS.index(tool) * len(F.STRATA) + F.STRATA.index(st)]
    c, s = _curve_drawn(fig.axes[2]), _curve_drawn(panel)
    rows = len((PHASE1 / "reports" / "evidence" / f"interval_lr_{tool}_{st}.csv").read_text().splitlines()) - 1
    assert len(c[0]) == rows
    for drawn_c, drawn_s1 in zip(c, s):
        assert np.array_equal(drawn_c, drawn_s1)
    assert fig.axes[2].get_ylim() == panel.get_ylim() == (L.C_FLOOR, L.C_CEIL)


def test_c_marks_each_in_sample_threshold_on_its_tier_at_full_precision(figure_one):
    from src import evid_common as K
    from src import evid_figures as F
    opts = figure_one["item"]["options"]
    t = K.read_back(PHASE1 / "reports" / "evidence" / "evidence_thresholds.csv")
    t = t[(t.tool == opts["c_tool"]) & (t.stratum == opts["c_stratum"]) & (t.status == "ok")
          & t.pp3_threshold_reachable.astype(bool)]
    want = [(float(t[t.tier == name.lower()].pp3_threshold_insample.iloc[0]), v)
            for name, v in F.PATH_CUTS if (t.tier == name.lower()).any()]
    marks = [(float(m.get_xdata()[0]), float(m.get_ydata()[0]))
             for m in figure_one["fig"].axes[2].lines if m.get_marker() == "|"]
    assert want and marks == want


def test_c_is_drawn_in_greys_and_names_its_table_by_the_configurations_number(figure_one):
    """The colours of (a) and (b) mark ClinVar record status; (c) draws every variant, in
    black and grey. Its key names the table of the thresholds by the configuration's
    number, which may change while the content name stays."""
    import matplotlib.colors as mcolors
    import matplotlib.text as mtext
    fig, cfg = figure_one["fig"], figure_one["cfg"]
    ax, key = fig.axes[2], fig.legends[-1]
    colours = [ln.get_color() for ln in ax.lines]
    colours += [tuple(x) for c in ax.collections for x in c.get_facecolor()]
    colours += [t.get_color() for t in ax.findobj(mtext.Text) if t.get_text().strip()]
    colours += [h.get_color() if hasattr(h, "get_color") else h.get_facecolor() for h in key.legend_handles]
    assert len(colours) > 10
    for x in colours:
        r, g, b = mcolors.to_rgb(x)
        assert max(r, g, b) - min(r, g, b) <= 0.03, x
    table = next(i for i in cfg["main"] if i["content"] == "evidence_tiers_and_thresholds")
    letter = "abcdefghij"[[p["part"] for p in table["parts"]].index("thresholds")]
    assert key.get_texts()[-1].get_text() == f"In-sample PP3 threshold ({table['number']}({letter}))"
    renumbered = {"main": [{"content": "t", "number": "Table 9", "parts": [{"part": "x"}, {"part": "y"}]}],
                  "supplementary": []}
    assert L._part_ref(renumbered, "t", "y") == "Table 9(b)"


def test_a_and_b_keep_their_place_from_the_top_and_c_takes_whole_pixels(figure_one):
    """(a) and (b) stay where they were in the figure before (c), measured from the top,
    at the same size, and (c)'s block is a whole number of pixels: the bitmap counts its
    rows from the top, so (a) and (b) fall on the pixels they fell on before."""
    fig, cfg = figure_one["fig"], figure_one["cfg"]
    height = fig.get_size_inches()[1] * 25.4
    fig0, axes0 = plt.subplots(1, 2, figsize=(L._width(cfg), L.AB_MM / 25.4), sharey=True,
                               gridspec_kw={"wspace": 0.08})
    fig0.subplots_adjust(left=0.2, right=0.97, top=0.94, bottom=0.2)
    for ax, ax0 in zip(fig.axes[:2], axes0):
        p, p0 = ax.get_position(), ax0.get_position()
        assert abs(p.x0 - p0.x0) < 1e-12 and abs(p.width - p0.width) < 1e-12
        assert abs((1 - p.y1) * height - (1 - p0.y1) * L.AB_MM) < 1e-9
        assert abs(p.height * height - p0.height * L.AB_MM) < 1e-9
    plt.close(fig0)
    pixels = (height - L.AB_MM) / 25.4 * cfg["artwork"]["dpi"]
    assert abs(pixels - round(pixels)) < 1e-6
    assert height <= cfg["artwork"]["max_height_mm"]


def test_c_stops_the_run_when_the_sources_do_not_list_its_curve(tmp_path):
    """The options name the predictor and band; the curve file they name must be among
    the item's sources, which contents.csv and manifest.csv list."""
    mp = pytest.MonkeyPatch()
    mp.chdir(PHASE1)
    try:
        cfg = L.load_config()
        item = dict(next(i for i in cfg["main"] if i["content"] == "published_cut_point_by_band"))
        item["options"] = dict(item["options"], c_stratum="s3_50")
        with pytest.raises(SystemExit, match="do not list"):
            L.figure_published_cut_point(item, tmp_path, cfg)
    finally:
        mp.undo()
        plt.close("all")


# ---------------------------------------------------------------------------
# E7 (the author's decisions of 2026-11-05): the combined Atlas score's folds that the
# held-out counts leave out, twelve_gene_refit's change columns, S1's caption line, the
# method's name
# ---------------------------------------------------------------------------
AVI = "Atlas combined score"


def _table2a_cells() -> dict:
    """Table 2(a)'s held-out cells of the combined Atlas score, k/n (m), by (band, tier)."""
    t = pd.read_csv(PHASE1 / "reports/evidence/tables/table3_evidence_tiers.csv", dtype=str,
                    keep_default_na=False)
    row = t[t.Predictor.str.startswith(AVI)].iloc[0]
    return {(band, tier): row[f"{band} bp: {tier}, held out"]
            for band in ("3–10", "11–50", "3–50") for tier in ("Moderate", "Strong")}


@need_out
def test_thresholds_pathogenic_counts_the_atlas_rows_without_the_two_folds_as_table_2a_counts_them():
    """thresholds_pathogenic: the combined Atlas score's rows list every fold, BRCA1's
    and RAD51C's marked "†", and count and summarise the others, so that their Moderate
    and Strong counts are Table 2(a)'s k/n (m) cell for cell; every other row, and every
    other cell of these rows, is the printed table's, line for line."""
    from src import evid_common as K
    from src import evid_supp_tables as S
    it, part = _part_of("thresholds_pathogenic", "pp3")
    src = L.read_source(PHASE1 / part["source"])[0]
    printed = pd.read_csv(OUT / "supplementary" / "thresholds_pathogenic.csv", dtype=str,
                          keep_default_na=False)
    assert list(printed.columns) == list(src.columns) and len(printed) == len(src)
    atlas = (src.Tool == AVI).to_numpy()
    assert atlas.sum() == 9
    before = (PHASE1 / part["source"]).read_text().splitlines()
    after = (OUT / "supplementary" / "thresholds_pathogenic.csv").read_text().splitlines()
    assert len(before) == len(after) and before[0] == after[0]
    assert all(b == a for b, a, x in zip(before[1:], after[1:], atlas) if not x)
    folds = K.read_back(PHASE1 / part["not_counted"]["folds"])
    folds = folds[(folds.side == "pp3") & (folds.tool == "avi")
                  & ~folds.heldout_gene.isin(K.AVI_SEEN_IN_TRAINING)]
    band = {S.STRATUM_LABEL[s]: s for s in S.IN_SCOPE}
    cells = _table2a_cells()
    for i in np.flatnonzero(atlas):
        r, s = printed.loc[i], src.loc[i]
        entries = s[L.S6_FOLDS].split("; ")
        assert r[L.S6_FOLDS] == "; ".join(e + ("†" if e.split("=")[0] in K.AVI_SEEN_IN_TRAINING else "")
                                          for e in entries)
        for c in printed.columns:
            if c not in L.S6_OVER_FOLDS and c != L.S6_FOLDS:
                assert r[c] == s[c], c
        f = folds[(folds.stratum == band[r["Stratum (intronic offset, bp)"]]) & (folds.tier == r.Tier.lower())]
        lr = f.heldout_lr.dropna()
        assert r["Held-out genes"] == str(len(f))
        assert r["Folds with an evaluable held-out LR"] == str(len(lr))
        assert r["Held-out genes clearing the cut"] == str(int((lr >= float(r["Tier LR cut"])).sum()))
        assert [r["Held-out LR, median"], r["Held-out LR, minimum"], r["Held-out LR, maximum"]] == \
            [S.fmt_lr(lr.median()), S.fmt_lr(lr.min()), S.fmt_lr(lr.max())]
        if r.Tier != "Supporting":
            n, ev, k = (int(r["Held-out genes"]), int(r["Folds with an evaluable held-out LR"]),
                        int(r["Held-out genes clearing the cut"]))
            assert f"{k}/{n}" + (f" ({n - ev})" if ev < n else "") == \
                cells[(r["Stratum (intronic offset, bp)"], r.Tier)]


def test_thresholds_pathogenic_stops_where_the_printed_row_and_the_folds_disagree_or_no_fold_is_picked(monkeypatch):
    monkeypatch.chdir(PHASE1)
    it, part = _part_of("thresholds_pathogenic", "pp3")
    src = L.read_source(PHASE1 / part["source"])[0]
    out, changed = L.not_counted(src, it, part)
    assert len(changed) == 9 and {k[0] for k in changed} == {AVI}
    bad = src.copy()
    bad.at[bad.index[bad.Tool == AVI][0], "Held-out genes clearing the cut"] = "0"
    with pytest.raises(SystemExit, match="disagree"):
        L.not_counted(bad, it, part)
    monkeypatch.setitem(L.MARK_RULES, "avi_seen_in_training", lambda r: False)
    with pytest.raises(SystemExit, match="picks no fold"):
        L.not_counted(src, it, part)


@need_out
def test_thresholds_benign_is_printed_as_it_stands_with_a_note_on_the_atlas_medians():
    it, part = _part_of("thresholds_benign", "bp4")
    assert (OUT / "supplementary" / "thresholds_benign.csv").read_bytes() == \
        (PHASE1 / part["source"]).read_bytes()
    c = pd.read_csv(OUT / "contents.csv", dtype=str, keep_default_na=False)
    note = c[c.content == "thresholds_benign"].notes.iloc[0]
    assert "Atlas combined score" in note and "*BRCA1* and *RAD51C*" in note and "median" in note


def test_a_column_to_leave_out_must_be_in_the_source():
    cfg = {"display": {"values": {}}, "columns": {"t": {"a": ["A", "count"]}}}
    df = pd.DataFrame({"a": ["1"], "b": ["2"]})
    assert list(L.displayed(df, {"content": "t"}, {"part": "p", "omit": ["b"]}, cfg).columns) == ["A"]
    with pytest.raises(SystemExit, match="no such column"):
        L.displayed(df, {"content": "t"}, {"part": "p", "omit": ["c"]}, cfg)


@need_out
def test_twelve_gene_refit_leaves_out_the_two_change_columns_which_its_source_keeps():
    """The text reports the change of tier only (the author's decision of 2026-10-25): the
    absolute change is the difference of the two thresholds printed beside it, and the
    relative change has no meaning where the threshold changes sign, as it does for
    GPN-MSA and the combined Atlas score at 3-10 bp, BP4 Supporting."""
    from src import evid_common as K
    it, part = _part_of("twelve_gene_refit", "refit")
    assert part["omit"] == ["abs_change", "rel_change"]
    d = K.read_back(PHASE1 / part["source"])
    both = d.threshold_7.notna() & d.threshold_12.notna()
    assert ((d.abs_change - (d.threshold_12 - d.threshold_7)).abs()[both] == 0).all()
    assert d.loc[~both, "abs_change"].isna().all()
    flips = d[both & (np.sign(d.threshold_7) * np.sign(d.threshold_12) < 0)]
    assert sorted(zip(flips.tool, flips.stratum, flips.side, flips.tier)) == \
        [("avi", "s3_10", "bp4", "supporting"), ("gpn_msa", "s3_10", "bp4", "supporting")]
    printed = pd.read_csv(OUT / "supplementary" / "twelve_gene_refit.csv", dtype=str,
                          keep_default_na=False)
    cols = [c for c in L.read_source(PHASE1 / part["source"])[0].columns if c not in part["omit"]]
    assert list(printed.columns) == [L._spec(it, part, c, _cfg())[0] for c in cols]
    assert len(printed) == len(d)


@need_out
def test_the_layout_names_the_method_as_the_text_does():
    """Local likelihood-ratio calibration (the author's decision of 2026-11-04): no
    spelling of the old name, interval calibration, in the configuration, the stage or
    any text the stage writes, the workbook's parts among them (the figures' text is
    checked in their PDFs by scripts/stage_E7/check_e7.py)."""
    old = re.compile(r"(?i)interval[\s_-]*calibration")
    texts = {"config/layout.yaml": CONFIG.read_text(),
             "src/evid_layout.py": (PHASE1 / "src/evid_layout.py").read_text()}
    for p in sorted(OUT.rglob("*.csv")):
        texts[str(p.relative_to(OUT))] = p.read_text()
    for book in ("supplementary_tables.xlsx", "submission/supplementary_tables.xlsx"):
        with zipfile.ZipFile(OUT / book) as z:
            for n in z.namelist():
                texts[f"{book}:{n}"] = z.read(n).decode()
    assert not [k for k, v in texts.items() if old.search(v)]
    c = pd.read_csv(OUT / "contents.csv", dtype=str, keep_default_na=False).set_index("content")
    assert c.loc["thresholds_pathogenic", "title"].startswith(
        "PP3 thresholds by local likelihood-ratio calibration for every predictor column")
    assert c.loc["thresholds_benign", "title"].startswith(
        "BP4 thresholds by local likelihood-ratio calibration for every predictor column")
    head = pd.read_csv(OUT / "supplementary" / "practical_metrics.by_threshold.csv", nrows=0).columns
    assert "Tier by the calibration rule" in head


@pytest.fixture(scope="module")
def figure_two(tmp_path_factory):
    """Figure 2 as the stage draws it, into a temporary folder."""
    import matplotlib
    out = tmp_path_factory.mktemp("layout2")
    kept, save = [], L._save
    mp = pytest.MonkeyPatch()
    mp.chdir(PHASE1)
    mp.setattr(L, "_save", lambda fig, folder, stem, cfg: (kept.append(fig), save(fig, folder, stem, cfg))[1])
    try:
        with matplotlib.rc_context():
            cfg = L.load_config()
            item = dict(next(i for i in cfg["main"] if i["content"] == "heldout_and_external_lr"))
            L.figure_heldout_external(item, out, cfg)
    finally:
        mp.undo()
    return kept[0]


def test_figure_two_draws_the_folds_not_counted_in_grey_and_leaves_them_out_of_the_median(figure_two):
    """(a) and (b): the combined Atlas score's folds with BRCA1 or RAD51C held out are
    drawn where they were, in grey; the median is the other folds'; no other row has a
    grey point; the key names them, in the first column with the other held-out entries."""
    from matplotlib.colors import to_rgba
    from src import evid_common as K
    from src import evid_figures as F
    folds = K.read_back(PHASE1 / "reports/evidence/evidence_thresholds_logo_folds.csv")
    folds = folds[(folds.side == "pp3") & np.isfinite(folds.heldout_lr)]
    grey, blue = to_rgba(F.MUTED), to_rgba(F.BLUE)
    row = F.COLUMNS.index("avi")
    seen_grey = 0
    for k, ax in enumerate(figure_two.axes[:6]):
        tier, st = ("moderate", "strong")[k // 3], F.STRATA[k % 3]
        pts = []
        for coll in ax.collections:
            fc = [tuple(c) for c in coll.get_facecolor()]
            ec = [tuple(c) for c in coll.get_edgecolor()]
            for j, (x, y) in enumerate(coll.get_offsets()):
                colour = fc[j % len(fc)] if fc and fc[j % len(fc)] != (1, 1, 1, 1) else ec[j % len(ec)]
                pts.append((round(float(y)), float(x), colour))
        assert not [p for p in pts if p[0] != row and p[2] == grey]
        assert [p for p in pts if p[2] == blue]
        f = folds[(folds.tool == "avi") & (folds.stratum == st) & (folds.tier == tier)]
        out = f.heldout_gene.isin(K.AVI_SEEN_IN_TRAINING)
        drawn_grey = sorted(x for r_, x, c in pts if r_ == row and c == grey)
        assert drawn_grey == sorted(np.clip(f.heldout_lr[out], 0.3, 1000).tolist())
        seen_grey += len(drawn_grey)
        medians = [ln.get_xdata()[0] for ln in ax.lines
                   if len(ln.get_ydata()) == 2 and abs(np.mean(ln.get_ydata()) - row) < 1e-9]
        counted = np.clip(f.heldout_lr[~out], 0.3, 1000)
        assert medians == ([float(np.median(counted))] if len(counted) else [])
    assert seen_grey == 9          # (a) 2 + 1 + 2, (b) 2 + 0 + 2: n.e. and no threshold draw none
    leg = figure_two.legends[0]
    names = [t.get_text() for t in leg.get_texts()]
    assert names[:4] == ["Held-out gene, clears the boundary", "Held-out gene, below it",
                         "Median of held-out genes", "Fold not counted (†)"]
    import matplotlib
    with matplotlib.rc_context():               # laid out in the stage's font, as drawn
        F._style()
        figure_two.canvas.draw()
        x0 = [t.get_window_extent().x0 for t in leg.get_texts()]
    assert x0[0] == x0[3] < x0[4]


def test_s1_keeps_its_caption_line_two_mm_clear_of_the_panels_at_240_mm(figure_one):
    """Supplementary Figure S1 at 240 mm: the panels' lower edge raised until the caption
    line stands at least 2 mm below the last row's tick labels, by the font boxes of the
    bitmap and of the PDF, no further than the rounding of the edge; every panel as tall
    as the others, the top where it was, the type as it was."""
    import matplotlib
    s1, cfg = figure_one["s1"], figure_one["cfg"]
    assert abs(s1.get_figheight() * 25.4 - 240) < 1e-9 and abs(s1.get_figwidth() * 25.4 - 175) < 1e-9
    with matplotlib.rc_context():               # measured in the stage's font, as drawn
        L.F._style()
        gap = L.caption_gap_mm(s1, cfg["artwork"]["dpi"])
    assert 2 <= gap < 2 + 1e-4 * 240 + 1e-9
    assert s1.subplotpars.top == 0.97 and s1.subplotpars.bottom > 0.03
    heights = {round(ax.get_position().height, 12) for ax in s1.axes}
    assert len(heights) == 1
    assert abs(min(ax.get_position().y0 for ax in s1.axes) - s1.subplotpars.bottom) < 1e-12
    sizes = {t.get_fontsize() for t in s1.findobj(matplotlib.text.Text)
             if t.get_visible() and t.get_text().strip()}
    assert sizes <= {5.5, 6.5, 7.0}


# ---------------------------------------------------------------------------
# E8: the supplementary items numbered in the order the text first cites them, and
# what is printed cleaned of the provisional marks, revision codes and v5 pointers
# (E8.1, E8.2, the author's decisions of 2026-11-07)
# ---------------------------------------------------------------------------
def _printed_config_texts(cfg: dict) -> list[str]:
    """Every text of the configuration the stage prints: titles, descriptions, an item's
    and a part's notes, part titles, a mark's text, what a replacement prints, the
    supplement's rule and numbering, and the repository-only entries."""
    fmt = cfg["supplementary_format"]
    out = [fmt["rule"], fmt["numbering"], fmt["conventions"]] + [f["description"] for f in fmt["files"]]
    out += [str(v) for v in fmt["submission_workbook"].values()]
    for it in cfg["main"] + cfg["supplementary"]:
        out += [it.get("title", ""), it.get("description", ""), it.get("legend", "")] + list(it.get("notes", []))
        for p in it.get("parts", []):
            out += [str(p.get("title", ""))] + list(p.get("notes", []))
            out += [r["new"] for r in p.get("replace", [])] + [m["text"] for m in p.get("marks", [])]
    out += [r["what"] for r in cfg["repository_only"]]
    return [t for t in out if t]


def test_each_list_is_in_its_numbers_order_from_one(tmp_path, monkeypatch):
    """Tables 1-4 and Figures 1-3, and the supplementary tables and figures, each listed in
    its numbers' order from 1 without a gap (E8.1), so that contents.csv and the
    workbook's sheets follow the numbers; a list out of order stops the run."""
    cfg = _cfg()
    for section, prefix, s in (("main", "", ""), ("supplementary", "Supplementary ", "S")):
        for kind in ("figure", "table"):
            have = [it["number"] for it in cfg[section] if it["kind"] == kind]
            assert have == [f"{prefix}{kind.capitalize()} {s}{k}" for k in range(1, len(have) + 1)]
    text = CONFIG.read_text()
    swapped = re.sub(r"number: Supplementary Table S(2|3)\b",
                     lambda m: f"number: Supplementary Table S{5 - int(m.group(1))}", text)
    assert len(re.findall(r"number: Supplementary Table S3\b", swapped)) == 1 and swapped != text
    (tmp_path / "layout.yaml").write_text(swapped)
    monkeypatch.chdir(PHASE1)
    with pytest.raises(SystemExit, match="not numbered in their list's order"):
        L.load_config(tmp_path / "layout.yaml")


def test_what_the_configuration_prints_names_items_only_as_content():
    """No printed text of the configuration writes a number of this paper by hand: every
    item it names is named as {content} and resolves (E8.1). None carries "provisional",
    a revision code (C1-C13, an E number) or a v5 number (E8.2)."""
    cfg = _cfg()
    numbers = L._numbers(cfg)
    texts = _printed_config_texts(cfg)
    assert len(texts) > 100
    for t in texts:
        assert not NUMBERED.search(SOURCE_PAPERS.sub("", t)), t
        assert not NOTE_NUMBERED.search(t), t
        assert not re.search(r"(?i)provisional|\bv5\b", t), t
        assert not re.search(r"\b[CE]\d{1,2}(?:\.\d+)?\b|\b[CE]\d{1,2}'s", t), t
        L._resolve(t, numbers)


def test_a_replacement_prints_its_text_exactly_as_often_as_declared(tmp_path):
    """A part's replacement is made in the headers, the cells and the notes above the
    part, in order; the run stops where the part holds the text more or less often than
    declared. A printed table's copy differs from its source's bytes in the replaced texts
    alone, its line endings and quoting kept."""
    src = tmp_path / "s.csv"
    src.write_bytes(b"# E9 note, see Table 9(b)\r\nA (Table 9(b)),B\r\n\"x, E9's\",Table 9(b)\r\n2,3\r\n")
    reps = L._replacements({"replace": [{"old": "Table 9(b)", "new": "{t}(b)", "count": 3},
                                        {"old": "E9's", "new": "the", "count": 1}]}, {"t": "Table 2"})
    assert reps == [("Table 9(b)", "Table 2(b)", 3), ("E9's", "the", 1)]
    df, notes = L._copy_replaced(src, tmp_path / "out.csv", reps, "t.p")
    assert (tmp_path / "out.csv").read_bytes() == \
        src.read_bytes().replace(b"Table 9(b)", b"Table 2(b)").replace(b"E9's", b"the")
    assert list(df.columns) == ["A (Table 2(b))", "B"]
    assert df.iloc[0].tolist() == ["x, the", "Table 2(b)"] and notes == ["E9 note, see Table 2(b)"]
    with pytest.raises(SystemExit, match="3 times"):
        L._copy_replaced(src, tmp_path / "o.csv", [("Table 9(b)", "x", 2)], "t.p")
    with pytest.raises(SystemExit):                 # across a delimiter: no cell holds it
        L._copy_replaced(src, tmp_path / "o.csv", [(",B", ";B", 1)], "t.p")
    df2, notes2 = L.replaced(pd.DataFrame({"a": ["E9's x"]}), ["E9's"], [("E9's", "the", 2)], "t.p")
    assert df2.a.tolist() == ["the x"] and notes2 == ["the"]
    with pytest.raises(SystemExit, match="1 times"):
        L.replaced(pd.DataFrame({"a": ["E9's x"]}), [], [("E9's", "the", 2)], "t.p")
    unchanged = pd.DataFrame({"a": ["1"]})
    assert L.replaced(unchanged, ["n"], [], "t.p")[0] is unchanged


@need_out
def test_the_replaced_texts_are_gone_from_what_is_printed():
    """Each part's replaced text is printed nowhere in its table, CSV or sheet (its cells'
    text, read as printed); what replaced it is."""
    cfg = _cfg()
    numbers = L._numbers(cfg)
    # each sheet's cells, as printed: a text's runs joined (a gene's italics aside)
    sheet = {name: "\n".join(L._plain(t) for row in rows.values() for t in row.values())
             for name, rows in _book_rows(OUT / "supplementary_tables.xlsx").items()}
    seen = 0
    for it in cfg["supplementary"]:
        for p in it.get("parts", []):
            book = sheet.get(L._sheet_name(it["number"]), "")
            for old, new, count in L._replacements(p, numbers):
                csv_text = (OUT / "supplementary" / L._part_path(Path("."), it, p).name).read_text()
                # the old text may stand inside the new one ("Supplementary Table S1 cites
                # no paper" holds "Table S1 cites no paper"): nowhere else
                assert old not in (csv_text.replace(new, "") if new else csv_text), (it["content"], old)
                assert old not in (book.replace(new, "") if new else book), (it["content"], old)
                if new:
                    assert new in book, (it["content"], new)
                seen += 1
    assert seen >= 30


@need_out
def test_the_readme_states_the_numbering_and_nothing_printed_is_provisional():
    fmt = _cfg()["supplementary_format"]
    with zipfile.ZipFile(OUT / "supplementary_tables.xlsx") as z:
        book = "".join(z.read(n).decode() for n in z.namelist())
    assert f"Supplementary tables ({fmt['numbering']}; titles and descriptions: config/layout.yaml)" in book
    assert "provisional" not in book.lower()
    contents = (OUT / "contents.csv").read_text()
    assert "provisional" not in contents.lower()


# ---------------------------------------------------------------------------
# E8.3: thousands, RefSeq transcripts, italics in the workbook; E8.4: the main figures
# (the author's decisions of 2026-11-07, D2-D4)
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("text, want", [
    ("1410", "1,410"), ("999", "999"), ("-1234.5", "-1,234.5"), ("8853", "8,853"), ("11022", "11,022"),
    ("0.0012704133987426758", "0.0012704133987426758"), ("2014", "2,014"),
    ("0.87 (1223/1410)", "0.87 (1,223/1,410)"), ("0.10 (138/1410)", "0.10 (138/1,410)"),
    ("over 7 genes; 2000 of 2000 draws evaluable", "over 7 genes; 2,000 of 2,000 draws evaluable"),
    ("PP3 LR 1935 of 2000, BP4 LR 0 of 2000", "PP3 LR 1,935 of 2,000, BP4 LR 0 of 2,000"),
    ("209 damaging, 1648 normal", "209 damaging, 1,648 normal"), ("distance 4999; unmasked", "distance 4,999; unmasked"),
    ("no (4475 lines, 1023 mm wide)", "no (4,475 lines, 1,023 mm wide)"),
    # not numbers to group: years, dates and times, identifiers, coordinates, versions
    ("2026-02-24", "2026-02-24"), ("2026-10-19T16:27:59+00:00", "2026-10-19T16:27:59+00:00"),
    ("Waters 2024", "Waters 2024"), ("(Sahu et al. 2025)", "(Sahu et al. 2025)"),
    ("(Google DeepMind, 2026), Methods", "(Google DeepMind, 2026), Methods"),
    ("Superseded on 11 September 2026 by urn:mavedb:00001250-a-3", "Superseded on 11 September 2026 by urn:mavedb:00001250-a-3"),
    ("(PMID 30209399; PMC6181777)", "(PMID 30209399; PMC6181777)"), ("ipynb at 55943112", "ipynb at 55943112"),
    ("for g.32394937, the set's c.9501+4A>C", "for g.32394937, the set's c.9501+4A>C"),
    ("API client 0.10.0", "API client 0.10.0"), ("ENST00000644876.2", "ENST00000644876.2"),
    ("Table 5's 6,551 rows", "Table 5's 6,551 rows"), ("Max Distance 10 000", "Max Distance 10 000"),
])
def test_numbers_of_a_thousand_or_more_are_grouped_but_not_years_dates_or_identifiers(text, want):
    assert L.grouped(text) == want


def test_transcripts_are_read_from_the_mane_summary_not_written_by_hand():
    """A MANE Select RefSeq prints as it stands; the Ensembl transcript a deposit names
    prints after the RefSeq of the same row of the summary; any other stops the run."""
    mane = "data/evidence/reference_cache/mane/MANE.GRCh38.v1.5.summary.txt.gz"
    part = {"transcripts": {"column": "T", "gene": "G", "mane": str(PHASE1 / mane)}}
    df = pd.DataFrame({"G": ["DDX3X", "SBDS", "BRCA1"], "T": ["ENST00000644876.2", "ENST00000246868.7", "NM_007294.4"]})
    assert L.transcripts(df, part, "t")["T"].tolist() == [
        "NM_001356.5 (ENST00000644876.2)", "NM_016038.4 (ENST00000246868.7)", "NM_007294.4"]
    for bad in (["DDX3X", "ENST00000644876.1"], ["BRCA1", "NM_007294.3"], ["NOSUCHGENE", "NM_1.1"]):
        with pytest.raises(SystemExit):
            L.transcripts(pd.DataFrame({"G": [bad[0]], "T": [bad[1]]}), part, "t")
    for content, part_name, column in (("functional_standards", "standards", "Transcript"),
                                       ("new_external_genes", "records", "Target transcript accession")):
        it, p = _part_of(content, part_name)
        assert p["transcripts"]["column"] == column and p["transcripts"]["mane"] == mane
        assert any("Transcripts are MANE Select (version 1.5)" in n for n in it["notes"])


@need_out
def test_the_printed_transcripts_are_refseq_with_the_deposits_ensembl_in_parentheses():
    s1 = pd.read_csv(OUT / "supplementary" / "functional_standards.csv", dtype=str, keep_default_na=False)
    s14 = pd.read_csv(OUT / "supplementary" / "new_external_genes.records.csv", dtype=str, keep_default_na=False)
    got = dict(zip(s1.Gene, s1.Transcript)) | dict(zip(s14.Gene, s14["Target transcript accession"]))
    assert got["DDX3X"] == "NM_001356.5 (ENST00000644876.2)" and got["SBDS"] == "NM_016038.4 (ENST00000246868.7)"
    assert all(v.startswith("NM_") for v in got.values()) and len(got) == 14


@need_out
def test_table_4b_groups_its_counts_in_thousands():
    t = pd.read_csv(OUT / "main" / "practical_performance.clinvar.csv", dtype=str, keep_default_na=False)
    cells = t.values.ravel().tolist()
    for want in ("1,410", "0.87 (1,223/1,410)", "0.10 (138/1,410)", "0.05 (68/1,410)"):
        assert want in cells, want
    assert not any(re.search(r"(?<![\d,.])\d{4,}", c) for c in cells)


def test_a_gene_marked_in_italics_is_an_italic_run_and_a_grouped_number_stays_a_number():
    xml = L._runs_xml("Direct for *BRCA1*, *RAD51C* and *DDX3X*", False)
    assert "*" not in xml and xml.count("<i/>") == 3 and "<t xml:space=\"preserve\">BRCA1</t>" in xml
    assert L._runs_xml("c.*123A>G at 0.2", False) == '<t xml:space="preserve">c.*123A&gt;G at 0.2</t>'
    assert L._cell_xml("A1", "1,410", False) == f'<c r="A1" s="{L.GROUPED_STYLE}"><v>1410</v></c>'
    assert L._cell_xml("A1", "1,040.5", False) == f'<c r="A1" s="{L.GROUPED_STYLE + 1}"><v>1040.5</v></c>'
    assert L._cell_xml("A1", "410", False) == f'<c r="A1" s="{L.NUMBER_STYLE}"><v>410</v></c>'
    assert 't="inlineStr"' in L._cell_xml("A1", "1,410", True)          # a header is text


@need_out
def test_the_workbook_prints_no_italic_mark_and_the_csv_keeps_them():
    with zipfile.ZipFile(OUT / "supplementary_tables.xlsx") as z:
        book = "".join(z.read(n).decode() for n in z.namelist() if n.startswith("xl/worksheets/"))
    assert not L.ITALIC.search(book) and "<i/>" in book
    assert "*BRCA1*" in (OUT / "supplementary" / "training_signal.csv").read_text()


@pytest.fixture(scope="module")
def main_figures(tmp_path_factory):
    """Figures 1-3 and the precision-recall figure as the stage draws them, each kept
    after the stage's own _save, into a temporary folder."""
    import matplotlib
    out = tmp_path_factory.mktemp("main")
    kept, save = {}, L._save
    mp = pytest.MonkeyPatch()
    mp.chdir(PHASE1)
    mp.setattr(L, "_save", lambda fig, folder, stem, cfg: (kept.__setitem__(stem, fig), save(fig, folder, stem, cfg))[1])
    try:
        with matplotlib.rc_context():
            cfg = L.load_config()
            for it in cfg["main"] + cfg["supplementary"]:
                if it["content"] in ("published_cut_point_by_band", "heldout_and_external_lr",
                                     "roc_curves", "precision_recall_curves"):
                    L.MAKERS[it["make"]](it, out, cfg)
    finally:
        mp.undo()
    return kept


def _line_widths(fig) -> list[float]:
    """Every spine's, tick's and drawn line's width in a figure (a marker's edge aside)."""
    import matplotlib.axes
    from matplotlib.lines import Line2D
    out = []
    for ax in fig.findobj(matplotlib.axes.Axes):
        out += [s.get_linewidth() for s in ax.spines.values() if s.get_visible()]
        for axis in (ax.xaxis, ax.yaxis):
            out += [t.tick1line.get_markeredgewidth() for t in axis.get_major_ticks() if t.tick1line.get_visible()]
    out += [ln.get_linewidth() for ln in fig.findobj(Line2D)
            if ln.get_linestyle() not in ("None", "", " ") and ln.get_visible()]
    return out


def test_the_main_figures_lines_are_at_least_one_point_and_the_supplements_as_they_were(main_figures):
    pt = _cfg()["artwork"]["main_lines_pt"]
    assert pt == 1.0
    for stem in ("published_cut_point_by_band", "heldout_and_external_lr", "roc_curves"):
        ws = _line_widths(main_figures[stem])
        assert ws and min(ws) >= pt, (stem, min(ws))
    assert min(_line_widths(main_figures["precision_recall_curves"])) < pt   # drawn as before


def test_figure_2_has_no_row_shading_and_crosses_for_the_folds_not_counted(main_figures):
    from matplotlib.colors import to_hex
    from matplotlib.markers import MarkerStyle
    from src import evid_figures as F
    fig = main_figures["heldout_and_external_lr"]
    band = to_hex(F.BAND)
    assert not [p for ax in fig.axes for p in ax.patches if to_hex(p.get_facecolor()) == band]
    m = MarkerStyle(L.CROSS)
    cross = m.get_path().transformed(m.get_transform()).vertices
    grey = [c for ax in fig.axes for c in ax.collections
            if len(c.get_offsets()) and to_hex(c.get_edgecolor()[0] if not len(c.get_facecolor())
                                              or to_hex(c.get_facecolor()[0]) == "#ffffff"
                                              else c.get_facecolor()[0]) == to_hex(F.MUTED)]
    assert grey and all(np.allclose(c.get_paths()[0].vertices, cross) for c in grey)
    key = [h for lg in fig.legends for h in lg.legend_handles if h.get_label().startswith("Fold not counted")]
    assert len(key) == 1 and key[0].get_marker() == L.CROSS


def test_figure_3_letters_its_panels_and_ticks_of_a_thousand_are_grouped(main_figures):
    fig = main_figures["roc_curves"]
    letters = [t.get_text() for ax in fig.axes for t in ax.texts if t.get_text() in ("a", "b", "c")]
    assert letters == ["a", "b", "c"]
    for stem in ("published_cut_point_by_band", "heldout_and_external_lr"):
        labels = {t.get_text() for t in main_figures[stem].findobj(plt.Text)}
        assert "1,000" in labels and "1000" not in labels, stem
    assert not any(re.fullmatch(r"\d{4,}", t.get_text()) for f in main_figures.values() for t in f.findobj(plt.Text))


# ---------------------------------------------------------------------------
# E8.10: the approved texts, the seven reworded cells, the workbook for submission
# ---------------------------------------------------------------------------
def _approved() -> dict:
    assert hashlib.sha256(APPROVED.read_bytes()).hexdigest() == APPROVED_SHA256
    return yaml.safe_load(APPROVED.read_text())


def _notes_of(entry: dict) -> list[str]:
    n = entry.get("notes") or []
    return [n] if isinstance(n, str) else list(n)


def _unxml(s: str) -> str:
    return s.replace("&lt;", "<").replace("&gt;", ">").replace("&quot;", '"').replace("&amp;", "&")


def _book_rows(path: Path) -> dict[str, dict[int, dict[str, str]]]:
    """A workbook's sheets, each {row: {column: text}}: a number as written, an inline
    string's runs joined with each italic run marked *GENE* again, as the CSV marks it."""
    out = {}
    with zipfile.ZipFile(path) as z:
        for name, i in re.findall(r'<sheet name="([^"]+)" sheetId="(\d+)"', z.read("xl/workbook.xml").decode()):
            rows: dict = {}
            for col, row, body in re.findall(r'<c r="([A-Z]+)(\d+)"[^>]*>(.*?)</c>',
                                             z.read(f"xl/worksheets/sheet{i}.xml").decode()):
                runs = re.findall(r"<r><rPr>(.*?)</rPr><t[^>]*>(.*?)</t></r>", body)
                if runs:
                    text = "".join(f"*{_unxml(t)}*" if "<i/>" in p else _unxml(t) for p, t in runs)
                else:
                    v = re.search(r"<v>(.*?)</v>", body)
                    text = v.group(1) if v else _unxml("".join(re.findall(r"<t[^>]*>(.*?)</t>", body)))
                rows.setdefault(int(row), {})[col] = text
            out[_unxml(name)] = rows
    return out


def _cell_text(cell: str) -> str:
    """A cell's text from its XML (as _sheet_cells gives it), an italic run as *GENE*."""
    runs = re.findall(r"<r><rPr>(.*?)</rPr><t[^>]*>(.*?)</t></r>", cell)
    if runs:
        return "".join(f"*{_unxml(t)}*" if "<i/>" in p else _unxml(t) for p, t in runs)
    v = re.search(r"<v>(.*?)</v>", cell)
    return v.group(1) if v else _unxml("".join(re.findall(r"<t[^>]*>(.*?)</t>", cell)))


def _sheet_cells(path: Path) -> dict[str, list[list[str]]]:
    """A workbook's sheets, each a list of rows, each row its cells' XML with the row
    number taken out of their references: what a cell holds, its style (the number's
    format, bold) and its runs (italics), so that two rows compare cell for cell."""
    out = {}
    with zipfile.ZipFile(path) as z:
        for name, i in re.findall(r'<sheet name="([^"]+)" sheetId="(\d+)"', z.read("xl/workbook.xml").decode()):
            xml = z.read(f"xl/worksheets/sheet{i}.xml").decode()
            out[_unxml(name)] = [re.findall(r'<c r="[A-Z]+\d+"(.*?)</c>', body)
                                 for body in re.findall(r"<row r=\"\d+\"[^>]*>(.*?)</row>", xml)]
    return out


def test_the_configuration_holds_the_approved_texts_word_for_word():
    """E8.10: every supplementary table's title, part titles and notes, every
    supplementary figure's title and legend, and the conventions, as the configuration
    prints them (its {content} references written as numbers), are the author's approved
    texts of 2026-11-09 character for character; the notes are the table's, none a part's."""
    A, cfg = _approved(), _cfg()
    numbers = L._numbers(cfg)

    def printed(t: str) -> str:
        return L.grouped(L._resolve(t, numbers))

    items = {it["content"]: it for it in cfg["supplementary"]}
    assert set(A["supplementary_tables"]) == {c for c, it in items.items() if it["kind"] == "table"}
    assert set(A["supplementary_figures"]) == {c for c, it in items.items() if it["kind"] == "figure"}
    assert printed(cfg["supplementary_format"]["conventions"]) == A["conventions"]
    compared = 1
    for c, a in A["supplementary_figures"].items():
        it = items[c]
        assert (printed(it["title"]), printed(it["legend"])) == (a["title"], a["legend"]), c
        assert not it.get("notes"), c
        compared += 2
    for c, a in A["supplementary_tables"].items():
        it = items[c]
        assert printed(it["title"]) == a["title"], c
        assert [printed(n) for n in it.get("notes", [])] == _notes_of(a), c
        parts = list((a.get("parts") or {}).values())
        if parts:
            assert [printed(p["title"]) for p in it["parts"]] == parts, c
        else:
            assert len(it["parts"]) == 1, c
        assert not any(p.get("notes") for p in it["parts"]), c
        compared += 1 + len(_notes_of(a)) + len(parts)
    assert compared > 100


@need_out
def test_every_printed_title_note_and_legend_is_the_approved_text():
    """E8.10: the same texts where the stage prints them: contents.csv (title, legend,
    notes), the workbook (README's titles; each sheet's title row, notes and part titles)
    and the workbook for submission (its contents sheet; each sheet's title row, notes and
    part titles), an italic run read back as *GENE*."""
    A, cfg = _approved(), _cfg()
    c = pd.read_csv(OUT / "contents.csv", dtype=str, keep_default_na=False).set_index("content")
    book = _book_rows(OUT / "supplementary_tables.xlsx")
    sub = _book_rows(OUT / "submission" / "supplementary_tables.xlsx")
    readme = {r["C"]: r for k, r in book["README"].items() if k >= 5 and "C" in r}
    compared = 0
    for name, a in A["supplementary_figures"].items():
        assert (c.loc[name, "title"], c.loc[name, "legend"], c.loc[name, "notes"]) == (a["title"], a["legend"], "")
        compared += 2
    for name, a in A["supplementary_tables"].items():
        it = next(i for i in cfg["supplementary"] if i["content"] == name)
        notes = _notes_of(a)
        parts = [f"({chr(97 + k)}) {t}" for k, t in enumerate((a.get("parts") or {}).values())]
        assert (c.loc[name, "title"], c.loc[name, "notes"]) == (a["title"], " | ".join(notes)), name
        assert readme[name]["D"] == a["title"], name
        compared += 2 + len(notes)
        sheet = L._sheet_name(it["number"])
        for wb in (book, sub) if it["destination"] == "excel" else (book,):
            rows = wb[sheet]
            assert rows[1] == {"A": f"{it['number']}. {a['title']}"}, name
            assert [rows[k] for k in range(2, 2 + len(notes))] == [{"A": n} for n in notes], name
            assert 2 + len(notes) not in rows, name                     # the empty row after them
            assert [r["A"] for r in rows.values() if r.keys() == {"A"} and r["A"] in parts] == parts, name
            compared += 1 + len(notes) + len(parts)
        if it["destination"] == "pdf":
            assert sheet not in sub, name
    assert sub["Contents"][1] == {"A": "Supplementary Tables"} and sub["Contents"][3] == {"A": A["conventions"]}
    assert compared > 300


@need_out
def test_the_workbook_for_submission_is_the_excel_tables_without_their_sources_comment_lines(tmp_path, monkeypatch):
    """E8.10: a contents sheet, then a sheet for each table marked excel in number order,
    each the repository workbook's sheet cell for cell (each value with its number format
    and italic runs) without the comment lines its sources carry above their headers, and
    nothing else left out; the contents sheet lists the tables by number and title and
    ends by naming the five tables marked pdf, listed from their numbers."""
    cfg = _cfg()
    numbers = L._numbers(cfg)
    tables = [it for it in cfg["supplementary"] if it["kind"] == "table"]
    excel = [it for it in tables if it["destination"] == "excel"]
    pdf = [it for it in tables if it["destination"] == "pdf"]
    assert len(excel) == 29 and [i["content"] for i in pdf] == [
        "training_signal", "column_agreement", "arms_at_fitted_thresholds", "clinvar_coverage",
        "inframe_true_positives"]
    sub, book = _sheet_cells(OUT / "submission" / "supplementary_tables.xlsx"), _sheet_cells(OUT / "supplementary_tables.xlsx")
    assert list(sub) == ["Contents"] + [L._sheet_name(it["number"]) for it in excel]
    monkeypatch.chdir(PHASE1)
    left_out = 0
    for it in excel:
        item = dict(it)
        L.table_parts(item, tmp_path / it["content"], cfg)
        comments = [L.grouped(L._resolve(n, numbers)) for part, _, notes in item["_frames"]
                    for n in notes if n not in part.get("notes", [])]
        a, b = book[L._sheet_name(it["number"])], sub[L._sheet_name(it["number"])]
        gone, j = [], 0
        for row in a:
            if j < len(b) and row == b[j]:
                j += 1
            else:
                gone.append(row)
        assert j == len(b), it["content"]                         # the sheet's rows, in order
        # what was left out: one text cell a row, set as a paragraph, each a comment line of
        # a source
        assert all(len(r) == 1 and f' s="{L.WRAP_STYLE}"' in r[0].split(">")[0] for r in gone), it["content"]
        assert [_cell_text(r[0]) for r in gone] == comments, it["content"]
        left_out += len(gone)
    assert left_out > 0
    rows = _book_rows(OUT / "submission" / "supplementary_tables.xlsx")["Contents"]
    listed = [rows[k] for k in sorted(rows) if k > 4 and "B" in rows[k]]
    assert listed == [{"A": it["number"], "B": L.grouped(L._resolve(it["title"], numbers))} for it in excel]
    assert rows[max(rows)] == {"A": "Supplementary Tables S5, S6, S12, S25 and S28 are in the Supplementary Information PDF."}


def test_several_items_are_listed_as_the_text_lists_them():
    """{a, b, c} prints the items' label once, in the plural, then their numbers in
    order, three or more in a row as a range, the last joined by "and" (E8.10)."""
    T = "Supplementary Table S"
    assert L.listed([f"{T}8", f"{T}9"]) == "Supplementary Tables S8 and S9"
    assert L.listed([f"{T}12", f"{T}10", f"{T}11"]) == "Supplementary Tables S10–S12"
    assert L.listed([f"{T}{k}" for k in (5, 6, 12, 25, 28)]) == "Supplementary Tables S5, S6, S12, S25 and S28"
    assert L.listed(["Supplementary Note 1", "Supplementary Note 2"]) == "Supplementary Notes 1 and 2"
    assert L.listed(["Table Q2", "Table Q3", "Table Q4", "Table Q7"]) == "Tables Q2–Q4 and Q7"
    numbers = {"a": f"{T}8", "b": f"{T}9", "n": "Supplementary Note 4"}
    assert L._resolve("in {a, b}; ({n})", numbers) == "in Supplementary Tables S8 and S9; (Supplementary Note 4)"
    for bad in ([f"{T}8", "Supplementary Figure S2"], ["Table 2", "Table 2"]):
        with pytest.raises(SystemExit):
            L.listed(bad)
    with pytest.raises(SystemExit, match="no item"):
        L._resolve("{a, nothing}", numbers)


def test_the_supplementary_notes_are_numbered_in_order_and_named_by_content(tmp_path, monkeypatch):
    notes = _cfg()["supplementary_notes"]
    assert [n["number"] for n in notes] == [f"Supplementary Note {k}" for k in range(1, 7)]
    assert [n["content"] for n in notes] == ["predictor_scores", "alphagenome_and_atlas", "elastic_net",
                                             "local_lr_calibration", "further_analyses",
                                             "external_gene_analyses"]
    text = CONFIG.read_text()
    assert text.count("number: Supplementary Note 3}") == 1
    (tmp_path / "layout.yaml").write_text(text.replace("number: Supplementary Note 3}", "number: Supplementary Note 9}"))
    monkeypatch.chdir(PHASE1)
    with pytest.raises(SystemExit, match="supplementary notes"):
        L.load_config(tmp_path / "layout.yaml")


def test_a_map_prints_a_number_it_does_not_name_only_when_its_entry_says_so():
    values = {"offset_bin": {"s3_10": "3–10"}}
    asis = [L.NUMBERS_AS_WRITTEN]
    assert [L._cell(t, "offset_bin", {}, asis, values, "w") for t in ("s3_10", "17", ">50")] == ["3–10", "17", ">50"]
    for text, extra in (("17", []), ("s9_9", asis), ("17b", asis)):
        with pytest.raises(SystemExit):
            L._cell(text, "offset_bin", {}, extra, values, "w")


@need_out
def test_the_seven_cells_the_author_reworded_print_as_approved(tmp_path, monkeypatch):
    """E8.10's cell_fixes, made in the layout's maps, replacements and drawing; the
    sources keep their own texts."""
    fixes = _approved()["cell_fixes"]
    numbers = L._numbers(_cfg())
    want = {(f["item"], f["where"].split(",")[0]): L._resolve(f["print"], numbers) for f in fixes}
    sup = OUT / "supplementary"
    audit = (sup / "label_sources.audit.csv").read_text()
    assert audit.count(want[("label_sources", "part audit")]) == 2
    assert audit.replace(want[("label_sources", "part audit")], "").count("Table S1 cites no paper") == 0
    bp = pd.read_csv(sup / "bp4_at_bp7_positions.csv", dtype=str, keep_default_na=False)
    texts = [L._resolve(f["print"], numbers) for f in fixes if f["item"] == "bp4_at_bp7_positions"]
    assert sorted(bp.Condition.value_counts().items()) == sorted((t, 24) for t in texts)
    assert "Offset band (bp)" in pd.read_csv(sup / "clinvar_coverage.csv", nrows=0).columns
    off = pd.read_csv(sup / "new_external_genes.counts_by_offset.csv", dtype=str, keep_default_na=False)
    schemes = next(f for f in fixes if f["where"].endswith("Bin scheme"))["print"].split("; ")
    bands = next(f for f in fixes if "Bin (offset" in f["where"])["print"].split("; ")
    assert set(off["Bin scheme"]) == set(schemes)
    assert set(off.loc[off["Bin scheme"] == schemes[1], "Bin (offset in bp, or stratum)"]) == set(bands)
    assert off.loc[off["Bin scheme"] == schemes[0], "Bin (offset in bp, or stratum)"].str.fullmatch(r"\d+|>50").all()
    src = (PHASE1 / "reports/evidence/suppl_analyses/c6_bp4_bp7.csv").read_text()
    assert "Table 3(b)'s BP4 Supporting threshold" in src                     # the source as it was
    kept, save = [], L._save
    monkeypatch.chdir(PHASE1)
    monkeypatch.setattr(L, "_save", lambda fig, folder, stem, cfg: (kept.append(fig), save(fig, folder, stem, cfg))[1])
    cfg = L.load_config()
    L.figure_clinvar_coverage(next(i for i in cfg["supplementary"] if i["content"] == "clinvar_coverage_by_position"),
                              tmp_path, cfg)
    assert kept[0].axes[1].get_xlabel() == want[("clinvar_coverage_by_position", "figure x-axis label")]


# ---------------------------------------------------------------------------
# E8.11: the main tables' approved texts and their editable files
# ---------------------------------------------------------------------------
W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"


def _docx(path: Path):
    """A main table's file: its body's blocks in order, each ("p", text, bold texts) or
    ("tbl", rows of cell texts, the table element), an italic run as *GENE*; and its
    section's page size, orientation and margins."""
    from xml.etree import ElementTree as ET
    with zipfile.ZipFile(path) as z:
        root = ET.fromstring(z.read("word/document.xml"))

    def runs(p):
        out, bold = [], []
        for r in p.iter(f"{W}r"):
            t = "".join(x.text or "" for x in r.iter(f"{W}t"))
            rpr = r.find(f"{W}rPr")
            it = rpr is not None and rpr.find(f"{W}i") is not None
            out.append(f"*{t}*" if it else t)
            if rpr is not None and rpr.find(f"{W}b") is not None:
                bold.append(t)
        return "".join(out), bold

    blocks = []
    body = root.find(f"{W}body")
    for el in body:
        if el.tag == f"{W}p":
            blocks.append(("p",) + runs(el))
        elif el.tag == f"{W}tbl":
            rows = [[runs(c)[0] for c in tr.findall(f"{W}tc")] for tr in el.findall(f"{W}tr")]
            blocks.append(("tbl", rows, el))
    sect = body.find(f"{W}sectPr")
    size, mar = sect.find(f"{W}pgSz"), sect.find(f"{W}pgMar")
    page = {"w": int(size.get(f"{W}w")), "h": int(size.get(f"{W}h")), "orient": size.get(f"{W}orient"),
            "margins": {k: int(mar.get(f"{W}{k}")) for k in ("top", "right", "bottom", "left")}}
    return blocks, page


def test_the_main_items_hold_the_approved_texts_word_for_word():
    """E8.11: the main tables' titles, part titles (none printed for Table 3, which has
    one part) and notes, and the main figures' titles, as the configuration prints them,
    are the approved texts; the main figures' legends are the manuscript's, not here."""
    A, cfg = _approved(), _cfg()
    numbers = L._numbers(cfg)
    main = {it["content"]: it for it in cfg["main"]}
    printed = lambda t: L.grouped(L._resolve(t, numbers))        # noqa: E731
    for c, a in A["main_tables"].items():
        it = main[c]
        assert it["number"] == a["number"] and printed(it["title"]) == a["title"], c
        assert [printed(n) for n in it["notes"]] == [a["notes"]], c
        titles = [t for t in a["parts"].values()]
        if len(titles) == 1:
            assert titles == [None] and len(it["parts"]) == 1 and "title" not in it["parts"][0], c
        else:
            assert [printed(p["title"]) for p in it["parts"]] == titles, c
    for c, a in A["main_figures"].items():
        assert main[c]["number"] == a["number"] and printed(main[c]["title"]) == a["title"], c
        assert "legend" not in main[c], c


@need_out
def test_each_main_table_is_an_editable_file_as_approved():
    """E8.11: each main table's file holds "Table N." in bold and its title, each part's
    title above its table, its tables cell for cell as the CSVs the stage wrote (italics as
    runs), and the notes in one paragraph below, all the approved texts; it is listed with
    the table's files and in the manifest."""
    A, cfg = _approved(), _cfg()
    c = pd.read_csv(OUT / "contents.csv", dtype=str, keep_default_na=False).set_index("content")
    manifest = pd.read_csv(OUT / "manifest.csv", dtype=str, keep_default_na=False)
    for it in [i for i in cfg["main"] if i["kind"] == "table"]:
        a = A["main_tables"][it["content"]]
        path = OUT / "main" / f"{it['content']}.docx"
        blocks, page = _docx(path)
        kinds = [b[0] for b in blocks]
        parts = [t for t in a["parts"].values() if t]
        want = ["p"] + (["p", "tbl"] * len(parts) if parts else ["tbl"]) + ["p"]
        assert kinds == want, it["content"]
        assert blocks[0][1] == f"{it['number']}. {a['title']}" and blocks[0][2] == [f"{it['number']}."]
        assert [b[1] for b in blocks[1:-1] if b[0] == "p"] == [f"({chr(97 + k)}) {t}" for k, t in enumerate(parts)]
        assert blocks[-1][1] == a["notes"] and not blocks[-1][2]
        tables = [b[1] for b in blocks if b[0] == "tbl"]
        for part, rows in zip(it["parts"], tables):
            assert rows == L.written_rows(OUT / "main" / L._part_path(Path("."), it, part).name), it["content"]
        assert c.loc[it["content"], "files"].endswith(f"phase1/reports/evidence/layout/main/{it['content']}.docx")
        assert c.loc[it["content"], "title"] == a["title"] and c.loc[it["content"], "notes"] == a["notes"]
        assert f"phase1/reports/evidence/layout/main/{it['content']}.docx" in set(manifest.file)


@need_out
def test_the_main_table_files_are_a4_ruled_across_only_and_reproducible():
    """E8.11: A4, margins of 2.5 cm, Tables 2 and 4 turned (wider than the upright page),
    Arial at 8 pt in the table and 9 pt around it; each table ruled at its top, below its
    header and at its foot, every rule at least 1 pt, with no rule down and no shading; the
    file's members dated as the workbook's, so that it is the same bytes on every run."""
    import xml.etree.ElementTree as ET
    cfg = _cfg()
    margin = round(25 * 1440 / 25.4)
    turned = {"evidence_tiers_and_thresholds", "practical_performance"}
    for it in [i for i in cfg["main"] if i["kind"] == "table"]:
        path = OUT / "main" / f"{it['content']}.docx"
        blocks, page = _docx(path)
        a4 = (11906, 16838) if it["content"] not in turned else (16838, 11906)
        assert (page["w"], page["h"]) == a4 and (page["orient"] == "landscape") == (it["content"] in turned)
        assert set(page["margins"].values()) == {margin}
        with zipfile.ZipFile(path) as z:
            assert {i.date_time for i in z.infolist()} == {L.ZIP_TIME}
            xml = z.read("word/document.xml").decode()
        assert "w:shd" not in xml and "w:vMerge" not in xml
        assert set(re.findall(r'w:rFonts w:ascii="([^"]+)"', xml)) == {"Arial"}
        for _, rows, tbl in [b for b in blocks if b[0] == "tbl"]:
            borders = tbl.find(f"{W}tblPr/{W}tblBorders")
            for side in ("top", "bottom"):
                e = borders.find(f"{W}{side}")
                assert e.get(f"{W}val") == "single" and int(e.get(f"{W}sz")) >= 8
            for side in ("left", "right", "insideH", "insideV"):
                assert borders.find(f"{W}{side}").get(f"{W}val") == "nil"
            trs = tbl.findall(f"{W}tr")
            for k, tr in enumerate(trs):
                for tc in tr.findall(f"{W}tc"):
                    b = tc.find(f"{W}tcPr/{W}tcBorders")
                    if k == 0:
                        assert [x.tag for x in b] == [f"{W}bottom"] and int(b[0].get(f"{W}sz")) >= 8
                    else:
                        assert b is None
            sizes = {int(s.get(f"{W}val")) for s in tbl.iter(f"{W}sz")}
            assert sizes == {16}
        outside = ET.fromstring(xml)
        texts = [p for p in outside.find(f"{W}body") if p.tag == f"{W}p"]
        assert {int(s.get(f"{W}val")) for p in texts for s in p.iter(f"{W}sz")} == {18}


def test_a_main_table_wider_than_the_turned_page_stops_the_run():
    cfg = _cfg()
    narrow = [["Gene", "N"], ["*RAD51C*", "1,152"]]
    assert L.page_for([narrow], cfg)[0] is False
    wide = [["Head"] * 12, ["supercalifragilistic"] * 12]
    with pytest.raises(SystemExit, match="turned page"):
        L.page_for([wide], cfg)
    turned = [["Head"] * 9, ["a rather longer cell text"] * 9]
    landscape, widths = L.page_for([turned], cfg)
    assert landscape and abs(sum(widths[0]) - 247) < 1e-6
    # Arial's widths, from Helvetica's metrics: an uppercase gene name is wider than its
    # letters' count suggests, a word of i, t and f narrower
    assert L.text_mm("RAD51C", 8) > L.text_mm("specificity", 8) * 0.8
    assert abs(L.text_mm("0", 10) - 0.556 * 10 * 25.4 / 72) < 1e-9


def test_table_4_stops_if_its_notes_misstate_the_band_it_leaves_out(tmp_path, monkeypatch):
    monkeypatch.chdir(PHASE1)
    cfg = L.load_config()
    item = dict(next(i for i in cfg["main"] if i["content"] == "practical_performance"))
    L.table_practical_performance(dict(item), tmp_path, cfg)                       # as approved
    item["notes"] = [item["notes"][0].replace("with 3 P/LP variants", "with 4 P/LP variants")]
    with pytest.raises(SystemExit, match="do not say"):
        L.table_practical_performance(item, tmp_path, cfg)


# ---------------------------------------------------------------------------
# E8 commit 8: the workbooks' column widths, wrapping and merged paragraphs; the gene
# symbols in the tables' cells in italics (H12); S32's two columns in words (the
# author's decisions of 2026-11-10)
# ---------------------------------------------------------------------------
GENES = ["BAP1", "BARD1", "BRCA1", "BRCA2", "PALB2", "RAD51C", "VHL", "DDX3X", "TP53", "CTCF", "SFPQ",
         "RAD51D", "XRCC2", "SBDS", "ATM"]


def _screen_px(width: float) -> int:
    """A column's width on the screen, in pixels, for Arial at 10 pt (ECMA-376, 18.3.1.13)."""
    return int(((256 * width + int(128 / 7)) / 256) * 7)


def test_the_configuration_lists_the_gene_symbols_and_the_workbooks_layout():
    cfg = _cfg()
    assert cfg["gene_symbols"] == GENES
    w = cfg["workbook"]
    assert w["min_width"] < w["max_width"] < w["max_word_width"] < w["paragraph_width"]
    assert w["header_lines"] >= 2 and w["line_pt"] * 2 <= L.MAX_ROW_PT


@pytest.mark.parametrize("text, want", [
    ("BAP1", "*BAP1*"), ("BAP1;VHL", "*BAP1*;*VHL*"), ("BAP1=10.4; BARD1=n.e.;", "*BAP1*=10.4; *BARD1*=n.e.;"),
    ("RAD51C=11.2†;", "*RAD51C*=11.2†;"),
    ("SNVs on CTCF-dependent HAP1 cell survival", "SNVs on *CTCF*-dependent HAP1 cell survival"),
    ("ENIGMA BRCA1 and BRCA2 Variant Curation Expert Panel",
     "ENIGMA *BRCA1* and *BRCA2* Variant Curation Expert Panel"),
    ("Direct for *BRCA1*, RAD51C and DDX3X", "Direct for *BRCA1*, *RAD51C* and *DDX3X*"),
    ("(BAP1, BARD1)", "(*BAP1*, *BARD1*)"), ("SBDS's deposit", "*SBDS*'s deposit"),
    ("an ATM prime-editing screen", "an *ATM* prime-editing screen"),
    # not a whole symbol; in a URN, a DOI, a web address, a path, a file's name, an identifier
    ("BRCA12 xBRCA1 BRCA1x BRCA1.2", "BRCA12 xBRCA1 BRCA1x BRCA1.2"),
    ("urn:mavedb:BRCA1-a-1", "urn:mavedb:BRCA1-a-1"), ("doi:10.1000/BRCA1", "doi:10.1000/BRCA1"),
    ("https://example.org/BRCA1", "https://example.org/BRCA1"),
    ("data/evidence/BRCA1/scores", "data/evidence/BRCA1/scores"), ("BRCA1.csv", "BRCA1.csv"),
    ("BRCA1_SGE", "BRCA1_SGE"), ("NM_007294.4(BRCA1):c.5074+2T>C", "NM_007294.4(BRCA1):c.5074+2T>C"),
    ("sge_scores BRCA1", "sge_scores *BRCA1*"),
    ("1,410", "1,410"), ("", ""), ("c.*123A>G at 0.2", "c.*123A>G at 0.2"),
])
def test_a_gene_symbol_is_set_in_italics_as_a_whole_word_and_never_in_an_identifier(text, want):
    assert L.italic_genes(text, GENES) == want
    assert L._plain(L.italic_genes(text, GENES)) == L._plain(text)        # the text printed is the same


def test_each_role_gives_its_cells_their_style():
    assert L._cell_xml("A1", "Gene", True, L.HEAD) == \
        f'<c r="A1" s="{L.HEAD_STYLE}" t="inlineStr"><is><t xml:space="preserve">Gene</t></is></c>'
    assert f' s="{L.WRAP_STYLE}" ' in L._cell_xml("A1", "a text", False, L.BODY)
    assert f' s="{L.WRAP_STYLE}" ' in L._cell_xml("A1", "a note", False, L.PARA)
    assert f' s="{L.BOLD_WRAP_STYLE}" ' in L._cell_xml("A1", "A title", True, L.PARA)
    # a number keeps its format whatever its row
    assert L._cell_xml("A1", "1,410", False, L.BODY) == f'<c r="A1" s="{L.GROUPED_STYLE}"><v>1410</v></c>'
    assert L._cell_xml("A1", "0.500", False, L.BODY) == f'<c r="A1" s="{L.NUMBER_STYLE + 3}"><v>0.500</v></c>'


def test_a_sheets_columns_fit_their_text_and_its_paragraphs_are_merged_to_the_paragraph_width():
    w = _cfg()["workbook"]
    path = "reports/evidence/suppl_analyses/c13_vcep_band_lr.csv"
    rows = [(["Supplementary Table S99. " + "A title long enough to wrap. " * 6], True, L.PARA),
            (["a note " * 150], False, L.PARA), ([], False, L.BLANK),
            (["Gene", "A header that wraps at its spaces into lines", "Path"], True, L.HEAD),
            (["*BRCA1*", "1.5", path], False, L.BODY),
            (["BAP1", "a text much longer than the widest a column may be " * 4, "x"], False, L.BODY)]
    widths, heights, merges = L.sheet_layout(rows, w)
    assert len(widths) == 3 and min(widths) >= w["min_width"]
    assert widths[1] <= w["max_width"] + 0.01                             # the long text wraps
    room = _screen_px(widths[2]) - L.PAD_PX
    assert L._lines(path, False, room, w["font_pt"]) == 1                # a path is kept whole
    assert sum(widths) >= w["paragraph_width"] - 0.01
    # the paragraphs, merged from A across the fewest columns as wide as paragraph_width
    end = merges[0].split(":")[1].rstrip("0123456789")
    k = ord(end) - 64
    assert merges == [f"A1:{end}1", f"A2:{end}2"]
    assert sum(widths[:k]) >= w["paragraph_width"] and sum(widths[:k - 1]) < w["paragraph_width"]
    # only the paragraphs are given a height, line_pt a line; the long note's the taller
    assert set(heights) == {1, 2} and heights[2] > heights[1] >= 2 * w["line_pt"]
    assert all(h / w["line_pt"] == int(h / w["line_pt"]) for h in heights.values())


def test_a_narrow_table_is_widened_to_the_paragraph_width_and_a_tall_paragraph_merged_further():
    w = _cfg()["workbook"]
    narrow = [(["Title"], True, L.PARA), ([], False, L.BLANK),
              (["A", "B"], True, L.HEAD), (["1", "a text long enough to wrap " * 3], False, L.BODY)]
    widths, _, merges = L.sheet_layout(narrow, w)
    assert abs(sum(widths) - w["paragraph_width"]) < 0.05 and widths[0] == w["min_width"]   # B, which wraps, takes it
    assert merges == ["A1:B1"]
    # a paragraph taller than a row may be at the paragraph width is merged across more
    wide = [(["word " * 2500], False, L.PARA), ([], False, L.BLANK),
            ([f"Header {j}" for j in range(8)], True, L.HEAD), (["x" * 60] * 8, False, L.BODY)]
    widths, heights, merges = L.sheet_layout(wide, w)
    k = ord(merges[0].split(":")[1].rstrip("0123456789")) - 64
    assert sum(widths[:k - 1]) >= w["paragraph_width"] and heights[1] <= L.MAX_ROW_PT
    # and one that no merge can hold stops the run
    with pytest.raises(SystemExit, match="taller than"):
        L.sheet_layout([(["word " * 20000], False, L.PARA), ([], False, L.BLANK),
                        (["A", "B"], True, L.HEAD), (["1", "2"], False, L.BODY)], w)


def _printed_cell(v: str, genes) -> str | None:
    """A CSV cell as _book_rows reads it back from the workbook: a number as written,
    without its commas; a text with its gene symbols in italics; None if empty."""
    if v == "":
        return None
    if 't="inlineStr"' not in L._cell_xml("A1", v, False, L.BODY):
        return v.replace(",", "")
    return L.italic_genes(v, genes)


@need_out
def test_the_tables_set_their_gene_symbols_in_italics_in_header_and_cells_and_the_csvs_do_not():
    """Every table's sheet, in both workbooks: its header and rows are its CSV's, each
    gene symbol of the configuration set in italics (read back as *GENE*); the CSVs keep
    their text as the stage wrote it."""
    cfg = _cfg()
    genes = cfg["gene_symbols"]
    italic = {}
    for wb in ("supplementary_tables.xlsx", "submission/supplementary_tables.xlsx"):
        book, n = _book_rows(OUT / wb), 0
        for it in (i for i in cfg["supplementary"] if i["kind"] == "table"):
            rows = book.get(L._sheet_name(it["number"]))
            if rows is None:
                assert wb.startswith("submission") and it["destination"] == "pdf", it["content"]
                continue
            after = 0
            for p in it["parts"]:
                header, *body = L.written_rows(OUT / "supplementary" / L._part_path(Path("."), it, p).name)
                want = [{L._col(j): x for j, c in enumerate(r) if (x := _printed_cell(c, genes)) is not None}
                        for r in [header] + body]
                start = next(k for k in sorted(rows) if k > after and rows[k] == want[0])
                assert [rows.get(start + i, {}) for i in range(len(want))] == want, (wb, it["content"], p["part"])
                after = start + len(want) - 1
                n += sum(t.count("*") // 2 for r in want for t in r.values())
        italic[wb] = n
    assert italic["supplementary_tables.xlsx"] > italic["submission/supplementary_tables.xlsx"] > 1000
    csv_text = (OUT / "supplementary" / "vcep_rules.csv").read_text()
    assert "*BRCA1*" not in csv_text and ",BRCA1," in csv_text


@need_out
def test_the_workbooks_size_their_columns_and_merge_and_size_their_paragraphs():
    w = _cfg()["workbook"]
    for wb in ("supplementary_tables.xlsx", "submission/supplementary_tables.xlsx"):
        with zipfile.ZipFile(OUT / wb) as z:
            for name, i in re.findall(r'<sheet name="([^"]+)" sheetId="(\d+)"', z.read("xl/workbook.xml").decode()):
                xml = z.read(f"xl/worksheets/sheet{i}.xml").decode()
                widths = [float(x) for x in re.findall(r'<col min="\d+" max="\d+" width="([\d.]+)" customWidth="1"/>', xml)]
                merges = re.findall(r'<mergeCell ref="A(\d+):([A-Z]+)\1"/>', xml)
                assert widths and min(widths) >= w["min_width"], (wb, name)
                assert sum(widths) >= w["paragraph_width"] - 0.01, (wb, name)
                assert merges and merges[0][0] == "1" and len({c for _, c in merges}) == 1, (wb, name)
                k = sum((ord(ch) - 64) * 26 ** e for e, ch in enumerate(reversed(merges[0][1])))
                assert sum(widths[:k]) >= min(w["paragraph_width"], sum(widths)) - 0.01, (wb, name)
                tall = [r for r, a in re.findall(r'<row r="(\d+)"([^>]*)>', xml) if "customHeight" in a]
                assert set(tall) <= {r for r, _ in merges}, (wb, name)
                # a paragraph's cell alone in its row
                for r, _ in merges:
                    assert re.findall(rf'<c r="([A-Z]+){r}"', xml) == ["A"], (wb, name, r)


@need_out
def test_s32_prints_its_rows_and_the_table_its_rows_were_checked_against_in_words():
    """S32's "Rows" by kind in words, and the BP7 rows' check named by the table's
    number and its condition (a); the source keeps its codes."""
    df = pd.read_csv(OUT / "supplementary" / "vcep_rules.csv", dtype=str, keep_default_na=False)
    assert df["Rows"].value_counts().to_dict() == {"Band likelihood ratio": 24, "BP7": 15, "No specification": 3}
    checked = df["Checked cell for cell against"]
    assert checked.value_counts().to_dict() == {"": 36, "Supplementary Table S31, condition (a)": 6}
    assert L._numbers(_cfg())["bp4_at_bp7_positions"] == "Supplementary Table S31"
    src, _ = L.read_source(PHASE1 / "reports" / "evidence" / "suppl_analyses" / "c13_vcep_band_lr.csv")
    assert src["part"].map(_cfg()["display"]["values"]["vcep_part"]).tolist() == df["Rows"].tolist()
    assert (src["checked_against"] != "").tolist() == (checked != "").tolist()
    assert set(src["checked_against"]) == {"", "c6_bp4_bp7.csv, (a) SpliceAI at or below the published cut, 0.1"}
