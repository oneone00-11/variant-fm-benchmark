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
if str(PHASE1) not in sys.path:
    sys.path.insert(0, str(PHASE1))
from src import evid_layout as L  # noqa: E402

CONFIG = PHASE1 / "config" / "layout.yaml"
OUT = PHASE1 / "reports" / "evidence" / "layout"
V301 = "8e0d609370565fe50933b2c9ec2f62325289e24d"     # release 3.0.1, the layout v5 used
NUMBERED = re.compile(r"(?i)\b(?:figure|table|fig|tab)\s*S?\d|(?:^|[_.])S\d+(?:[_.]|$)")
MM_PER_PT = 25.4 / 72

need_out = pytest.mark.skipif(not (OUT / "manifest.csv").exists(), reason="the layout stage has not run")


def _cfg() -> dict:
    return yaml.safe_load(CONFIG.read_text())


def _items() -> list[dict]:
    c = _cfg()
    return c["main"] + c["supplementary"]


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
    if subprocess.run(["git", "-C", str(REPO), "cat-file", "-e", f"{V301}^{{commit}}"],
                      capture_output=True).returncode != 0:
        pytest.skip("release 3.0.1's commit is not in this clone")
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
    numbers = {it["content"]: it["number"] for it in _items()}
    texts = [n for it in _items() for n in it.get("notes", [])]
    texts += [n for it in _items() for p in it.get("parts", []) for n in p.get("notes", [])]
    texts += [r["what"] for r in cfg["repository_only"]]
    refs = [m for t in texts for m in L.REFERENCE.findall(t)]
    assert refs and all(r in numbers for r in refs)
    for t in texts:
        assert not NUMBERED.search(t), t
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
    assert supp.fits_one_page_estimate.str.match(r"^(yes|no) \(\d+ lines, \d+ mm wide\)$").all()
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
    below), and the texts a part prints in place of its source's (its `replace`, E8): the
    copy is then the source's bytes with those texts replaced, and nothing else."""
    numbers = L._numbers(_cfg())
    replaced = 0
    for it in _items():
        for part in it.get("parts", []):
            src = part.get("source")
            if not src or not src.startswith(L.PRINTED) or "not_counted" in part:
                continue
            section = "main" if it in _cfg()["main"] else "supplementary"
            name = f"{it['content']}.csv" if len(it["parts"]) == 1 else f"{it['content']}.{part['part']}.csv"
            want = (PHASE1 / src).read_bytes()
            for old_, new_, _ in L._replacements(part, numbers):
                want, replaced = want.replace(old_.encode(), new_.encode()), replaced + 1
            assert (OUT / section / name).read_bytes() == want, name
    assert replaced == 2            # training_signal's fusion row, external_genes' five headers


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
    c = pd.read_csv(OUT / "contents.csv", dtype=str, keep_default_na=False)
    assert "(b) leaves out 11–50 bp: 3 P/LP variants of at least one star" in \
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
    assert cells and all(int(s) - L.NUMBER_STYLE == len((v.split(".") + [""])[1]) for s, v in cells)


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
    assert "combined Atlas score" in note and "BRCA1 and RAD51C" in note and "median" in note


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
    with zipfile.ZipFile(OUT / "supplementary_tables.xlsx") as z:
        for n in z.namelist():
            texts[f"supplementary_tables.xlsx:{n}"] = z.read(n).decode()
    assert not [k for k, v in texts.items() if old.search(v)]
    c = pd.read_csv(OUT / "contents.csv", dtype=str, keep_default_na=False).set_index("content")
    assert c.loc["thresholds_pathogenic", "title"].startswith(
        "Local likelihood-ratio calibration thresholds, pathogenic side")
    assert c.loc["thresholds_benign", "title"] == "Local likelihood-ratio calibration thresholds, benign side"
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
    out = [fmt["rule"], fmt["numbering"]] + [f["description"] for f in fmt["files"]]
    for it in cfg["main"] + cfg["supplementary"]:
        out += [it.get("title", ""), it.get("description", "")] + list(it.get("notes", []))
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
        assert not NUMBERED.search(t), t
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
    """Each part's replaced text is printed nowhere, CSV or workbook; what replaced it is."""
    cfg = _cfg()
    numbers = L._numbers(cfg)
    with zipfile.ZipFile(OUT / "supplementary_tables.xlsx") as z:
        book = "".join(z.read(n).decode() for n in z.namelist() if n.startswith("xl/worksheets/"))
    seen = 0
    for it in cfg["supplementary"]:
        for p in it.get("parts", []):
            for old, new, count in L._replacements(p, numbers):
                csv_text = (OUT / "supplementary" / L._part_path(Path("."), it, p).name).read_text()
                assert old not in csv_text, (it["content"], old)
                assert escape_xml(old) not in book, (it["content"], old)
                if new:
                    assert escape_xml(new) in book, (it["content"], new)
                seen += 1
    assert seen >= 30


def escape_xml(s: str) -> str:
    from xml.sax.saxutils import escape
    return escape(s)


@need_out
def test_the_readme_states_the_numbering_and_nothing_printed_is_provisional():
    fmt = _cfg()["supplementary_format"]
    with zipfile.ZipFile(OUT / "supplementary_tables.xlsx") as z:
        book = "".join(z.read(n).decode() for n in z.namelist())
    assert f"Supplementary tables ({fmt['numbering']}; titles and descriptions: config/layout.yaml)" in book
    assert "provisional" not in book.lower()
    contents = (OUT / "contents.csv").read_text()
    assert "provisional" not in contents.lower()
