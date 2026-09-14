"""Main Figures 1-3 are generated from tracked tables, not pasted.

Until frozen-matrix-v2 the figures had no generator: their captions were updated
to v2 while the images stayed at v1. phase1/src/phase9_figures.py now draws them
from the report tables, and reproduce_calibration.py runs it. These tests keep
that true and keep the images at the frames the manuscript embeds them in, so a
redraw replaces the picture without distorting it.
"""
from __future__ import annotations

from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
FIG = REPO / "phase1" / "reports" / "phase1" / "figures"
# pixel size of each PNG = the aspect ratio of its frame in the manuscript
GEOMETRY = {"fig1_forest_ranking": (2431, 996),
            "fig2_reliability": (1165, 1135),
            "fig3_calibration_advantage": (2431, 995)}


def test_every_main_figure_is_written_in_three_formats():
    for stem in GEOMETRY:
        for ext in ("png", "pdf", "svg"):
            p = FIG / f"{stem}.{ext}"
            assert p.exists() and p.stat().st_size > 1000, f"{p.name} missing or empty"


def test_figures_are_drawn_at_the_manuscript_frame_geometry():
    from PIL import Image

    for stem, size in GEOMETRY.items():
        assert Image.open(FIG / f"{stem}.png").size == size, stem


def test_the_figure_stage_is_part_of_the_documented_reproduction():
    assert '"src.phase9_figures"' in (REPO / "scripts" / "reproduce_calibration.py").read_text()
