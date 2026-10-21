"""Every figure of the evidence pipeline is drawn in the typeface the repository ships.

matplotlib falls back to another font, with no more than a log line, when the one asked
for cannot be found, so the figures could change typeface on a machine that lacks it.
This draws every figure again (the main and supplementary figures, the practical-metrics
figure and slide, the draft PNGs), into a temporary folder, while recording each font
file matplotlib resolves: every one must be a file in phase1/assets/fonts/, and no glyph
may be missing from them. The font files themselves are the release files SOURCE.json
records.
"""
from __future__ import annotations

import hashlib
import json
import sys
import warnings
from pathlib import Path

import pytest

PHASE1 = Path(__file__).resolve().parents[1]
FONTS = PHASE1 / "assets" / "fonts"
EV = PHASE1 / "reports" / "evidence"
FACES = {"LiberationSans-Regular.ttf", "LiberationSans-Bold.ttf",
         "LiberationSans-Italic.ttf", "LiberationSans-BoldItalic.ttf"}

if str(PHASE1) not in sys.path:
    sys.path.insert(0, str(PHASE1))


def test_the_font_files_are_the_release_files_source_json_records():
    src = json.loads((FONTS / "SOURCE.json").read_text())
    assert src["family"] == "Liberation Sans"
    assert set(src["files"]) == FACES | {"LICENSE"}
    for name, f in src["files"].items():
        assert hashlib.sha256((FONTS / name).read_bytes()).hexdigest() == f["sha256"], name
    assert {p.name for p in FONTS.glob("*.ttf")} == FACES
    assert "SIL OPEN FONT LICENSE VERSION 1.1" in " ".join(
        (FONTS / "LICENSE").read_text().upper().split())


@pytest.fixture(scope="module")
def drawn(tmp_path_factory):
    """Every figure drawn again into a temporary folder, with the font files resolved
    on the way and the warnings raised."""
    import matplotlib
    import pandas as pd
    import yaml
    from matplotlib import font_manager
    from src import evid_fig_data as D
    from src import evid_figures as F
    from src import evid_practical_metrics as P
    out = tmp_path_factory.mktemp("figures")
    resolved, opened = [], []
    # Every lookup by name passes through _findfont_cached (mathtext reaches it through
    # a module-level alias of findfont, so findfont itself is not the place), and every
    # font file opened, mathtext's by path included, through _get_font.
    lookup, load = font_manager.FontManager._findfont_cached, font_manager._get_font

    def recording_lookup(self, *args, **kwargs):
        found = lookup(self, *args, **kwargs)
        if isinstance(found, str):
            resolved.append(found)
        return found

    def recording_load(paths, *args, **kwargs):
        opened.extend(str(p) for p in paths)
        return load(paths, *args, **kwargs)

    recording_lookup.cache_clear = lookup.cache_clear      # addfont clears it
    recording_load.cache_clear = load.cache_clear
    mp = pytest.MonkeyPatch()
    mp.setattr(font_manager.FontManager, "_findfont_cached", recording_lookup)
    mp.setattr(font_manager, "_get_font", recording_load)
    mp.chdir(PHASE1)                         # the modules read reports/evidence/...
    mp.setattr(F, "OUT", out)
    mp.setattr(D, "FIG_DIR", out)
    try:
        with matplotlib.rc_context(), warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            F.main()
            D.draft_pngs(D.fig1(), *D.fig2(), D.fig3(), D.fig4())
            cfg = yaml.safe_load(open(P.CONFIG_PATH))
            roc, pr, auc, table = (pd.read_csv(p) for p in
                                   (P.OUT_ROC, P.OUT_PR, P.OUT_AUC, P.OUT_TABLE))
            P.figure(roc, pr, auc, table, cfg)
            P.figure_slides(roc, pr, auc, table, cfg)
    finally:
        mp.undo()
    return {"dir": out, "resolved": resolved, "opened": opened,
            "warnings": [str(w.message) for w in caught]}


def test_every_figure_was_drawn_again(drawn):
    committed = {p.name for d in ("figures", "fig_data") for p in (EV / d).iterdir()
                 if p.suffix in (".png", ".pdf")}
    assert {p.name for p in drawn["dir"].iterdir()} == committed


def test_every_font_file_opened_is_a_repository_file(drawn):
    opened = {Path(p).resolve() for p in drawn["opened"]}
    assert opened, "no font was opened"
    assert {p for p in opened if p.parent != FONTS.resolve()} == set()
    # regular text, bold titles and panel letters, italic gene symbols
    assert {"LiberationSans-Regular.ttf", "LiberationSans-Bold.ttf",
            "LiberationSans-Italic.ttf"} <= {p.name for p in opened}


def test_every_font_looked_up_by_name_is_a_repository_file(drawn):
    """All but one: mathtext looks up Computer Modern's cmex10 by its name, for large
    delimiters, whatever the settings (matplotlib/_mathtext.py, UnicodeFonts); no
    figure has one, and the test above holds that the file is never opened."""
    resolved = {Path(p).resolve() for p in drawn["resolved"]}
    assert resolved, "no font was looked up"
    assert {p.name for p in resolved if p.parent != FONTS.resolve()} <= {"cmex10.ttf"}


def test_no_glyph_is_missing_from_the_font(drawn):
    assert [w for w in drawn["warnings"] if "missing from" in w] == []
