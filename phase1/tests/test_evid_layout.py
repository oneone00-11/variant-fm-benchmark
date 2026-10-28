"""The layout layer (src/evid_layout.py, E5.3): numbers only in config/layout.yaml, files
named for their content, every table and figure from an output already written, the
printed tables as they stand, the new tables under readable names and S1-S13's number
formats, the figures within the artwork guide, and a place for every output added since
release 3.0.1."""
from __future__ import annotations

import hashlib
import re
import subprocess
import sys
import zipfile
from pathlib import Path

import matplotlib.pyplot as plt
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
    """(item, part, its source's columns) for every part the stage prints itself."""
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
            yield it, part, cols


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
def test_a_cell_prints_as_s1_s13_print_its_type(text, kind, row, want):
    values = _cfg()["display"]["values"]
    assert L._cell(text, kind, row, [], values, "test") == want


def test_a_threshold_keeps_the_digits_that_keep_its_band():
    """evid_supp_tables.fmt_thr, as S6 prints a threshold: four significant figures, and
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
    for it in _items():
        for part in it.get("parts", []):
            src = part.get("source")
            if not src or not src.startswith(L.PRINTED):
                continue
            section = "main" if it in _cfg()["main"] else "supplementary"
            name = f"{it['content']}.csv" if len(it["parts"]) == 1 else f"{it['content']}.{part['part']}.csv"
            assert (OUT / section / name).read_bytes() == (PHASE1 / src).read_bytes(), name


@need_out
def test_the_new_tables_keep_every_row_and_column_of_their_sources():
    for it, part, cols in _new_parts():
        name = f"{it['content']}.csv" if len(it["parts"]) == 1 else f"{it['content']}.{part['part']}.csv"
        printed = pd.read_csv(OUT / "supplementary" / name, dtype=str, keep_default_na=False)
        srcs = part.get("sources") or [part["source"]]
        srcs = list(srcs.values()) if isinstance(srcs, dict) else srcs
        rows = sum(len(L.read_source(PHASE1 / s)[0]) for s in srcs)
        assert printed.shape == (rows, len(cols)), name


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
    rounded twice (0.37746 -> 0.3775 -> 0.378, where Supplementary Table S13 prints
    0.377), so the configuration names four."""
    values = _cfg()["display"]["values"]
    assert L._cell("0.366", "proportion", {}, [4], values, "t") == "0.3660"
    assert L._cell("0.3775", "proportion", {}, [4], values, "t") == "0.3775"
    spec = _cfg()["columns"]["inframe_true_positives"]["in_frame_share"]
    assert spec[1:] == ["proportion", 4]
