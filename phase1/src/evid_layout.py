"""The layout layer (E5.3): the manuscript's figures and tables, main and supplementary,
under the numbers config/layout.yaml gives them.

Every figure and table here is drawn or written from an output an earlier stage wrote;
nothing is computed that a reader could not find there, and no earlier output is
touched: the v3.0.1 layout under figures/, tables/ and supplement/ stays as it was, so
that release can still be reproduced. One run writes them all, the main figures and
tables, the supplementary figures, every supplementary table as CSV and the workbook,
so that the text and its supplement come from one commit.

Numbers live only in config/layout.yaml, and in what this stage writes from it: the
workbook's README sheet and sheet names, contents.csv and manifest.csv, and every text
that names another item (written {content} in the configuration: a title, a part's
title, a legend, a note, the conventions, or a part's replacement of a source's text;
{a, b} names several, printed as the text lists them, listed()). The supplementary notes
are numbered there too, though the stage prints none of them. Every file is named for
its content, so renumbering figures or tables changes the configuration alone; each
list there is in its number's order (the supplementary items numbered in the order the
text first cites them, E8, 2026-11-08). The supplementary tables' and figures' titles,
part titles, notes and legends, and the conventions, are the author's wording
(approved 2026-11-09, E8.10), printed exactly as the configuration writes them.

Figures follow the configuration's artwork block: no wider than its double-column width
and no taller than its height limit, type no smaller and no larger than its sizes,
bitmaps at its resolution. A larger size set by evid_figures (its titles are 7.5 pt, its
panel letters 9 pt) is brought down to the largest; a smaller one, or a figure too wide
or too tall, stops the run. The three supplementary figures of v3.0.1 (the local
likelihood-ratio curves, the threshold dilution, AUROC by ClinVar record status) are
drawn again by their own functions in evid_figures, from the same outputs, at this width
(the curves also at the height limit, shorter than evid_figures draws them, the panels
raised clear of the caption line at the foot); the files under figures/ are not touched.

Tables. A printed table (tables/, supplement/) is written as it stands, byte for byte,
except where its entry names folds a rule leaves out of the counts: the pathogenic-side
thresholds (thresholds_pathogenic) list the combined Atlas score's folds that Table 2(a)
does not count, marked, and count their rows without them (not_counted()). Every other
table keeps all its columns under the readable names its entry in the configuration
gives, and each cell is printed by its column's type the way the tables v3.0.1 printed
print it (evid_supp_tables' formatters); the source files keep their full precision. A
part's entry can change what is printed, never the source: leave out a column that
repeats another cell for cell, or one the text does not use, print "–" where a count has
no definition, or add a column that marks the rows a rule picks (displayed()); and print
a text of its own in place of a source's, in its headers, cells and notes, each exactly
as often as the entry says (replaced(): a printed table's copy differs from its source
in those texts alone). A table of several parts is written as one CSV per part. A
source's comment lines above its header are kept as the part's notes, in the workbook
above the part.

Outputs, under reports/evidence/layout/:
  main/<content>.pdf, .png                the main figures
  main/<content>[.<part>].csv             the main tables
  main/<content>.docx                     each main table as an editable file (E8.11)
  supplementary/<content>.pdf, .png       the supplementary figures
  supplementary/<content>[.<part>].csv    the supplementary tables
  supplementary_tables.xlsx               the supplementary tables, README sheet first
  submission/supplementary_tables.xlsx    the tables marked excel, for submission: a
                                          contents sheet, then each table without its
                                          sources' comment lines (E8.10)
  contents.csv                            number, content, title, legend, destination,
                                          description, notes, files and sources, and
                                          the outputs kept in the repository only
  manifest.csv                            every file written, with its sha256 and sources

Run (PYTHONPATH=phase1):  python -m src.evid_layout
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import io
import os
import re
import shutil
import zipfile
from contextlib import contextmanager
from pathlib import Path
from xml.sax.saxutils import escape

import matplotlib

matplotlib.use("Agg")
import matplotlib.axes  # noqa: E402
import matplotlib.collections  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.text as mtext  # noqa: E402
import matplotlib.ticker  # noqa: E402
import matplotlib.transforms as mtransforms  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

from . import evid_common as K  # noqa: E402
from . import evid_external as X  # noqa: E402
from . import evid_figures as F  # noqa: E402
from . import evid_practical_metrics as P  # noqa: E402
from . import evid_supp_tables as S  # noqa: E402

CONFIG = Path("config/layout.yaml")
REPORT_DIR = Path("reports/evidence")
OUT = REPORT_DIR / "layout"
REPO_PREFIX = "phase1/"            # paths in the contents and manifest, from the repository root

MAKERS: dict = {}


def maker(fn):
    MAKERS[fn.__name__] = fn
    return fn


def _rel(path) -> str:
    """A path from phase1/ as the repository names it (a file above phase1/ too)."""
    return os.path.normpath(REPO_PREFIX + Path(path).as_posix())


def _sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# ---------------------------------------------------------------------------
# Figures: the artwork block's width and type sizes
# ---------------------------------------------------------------------------
def _width(cfg: dict) -> float:
    """The double-column width, in inches."""
    return cfg["artwork"]["double_column_mm"] * F.MM


def _save(fig, folder: Path, stem: str, cfg: dict) -> list[Path]:
    """PDF and PNG at the artwork block's width, type sizes and resolution, without the
    timestamps and versions that would change every run."""
    art = cfg["artwork"]
    small, large = art["text_pt"]
    if fig.get_figwidth() > _width(cfg) + 1e-9:
        raise SystemExit(f"[layout] {stem}: {fig.get_figwidth() / F.MM:.1f} mm wide, wider "
                         f"than {art['double_column_mm']} mm")
    if fig.get_figheight() > art["max_height_mm"] * F.MM + 1e-9:
        raise SystemExit(f"[layout] {stem}: {fig.get_figheight() / F.MM:.1f} mm tall, taller "
                         f"than {art['max_height_mm']} mm")
    _grouped_ticks(fig)
    fig.canvas.draw()                       # tick labels exist only once drawn
    for t in fig.findobj(mtext.Text):
        if not t.get_visible() or not t.get_text().strip():
            continue
        if t.get_fontsize() > large:
            t.set_fontsize(large)
        elif t.get_fontsize() < small:
            raise SystemExit(f"[layout] {stem}: '{t.get_text()}' is set at "
                             f"{t.get_fontsize()} pt, below {small} pt")
    folder.mkdir(parents=True, exist_ok=True)
    pdf, png = folder / f"{stem}.pdf", folder / f"{stem}.png"
    fig.savefig(pdf, metadata={"CreationDate": None, "ModDate": None, "Producer": None,
                               "Creator": None})
    fig.savefig(png, dpi=art["dpi"], metadata={"Software": None})
    plt.close(fig)
    return [pdf, png]


class _GroupedTicks(matplotlib.ticker.Formatter):
    """An axis's own tick labels, each number of 1,000 or more grouped in thousands
    (1000 -> 1,000), the labels below 1,000 as they were."""

    def __init__(self, inner):
        self.inner = inner

    def __call__(self, x, pos=None):
        return grouped(self.inner(x, pos))

    def format_ticks(self, values):
        return [grouped(s) for s in self.inner.format_ticks(values)]

    def set_locs(self, locs):
        self.inner.set_locs(locs)

    def get_offset(self):
        return self.inner.get_offset()

    def set_axis(self, axis):
        super().set_axis(axis)
        self.inner.set_axis(axis)


def _grouped_ticks(fig) -> None:
    """Every axis of a figure, insets included, with its tick labels grouped in
    thousands (E8.4: the target's style), whatever drew it; a figure without such a
    label is drawn as it was."""
    for ax in fig.findobj(matplotlib.axes.Axes):
        for axis in (ax.xaxis, ax.yaxis):
            for get, put in ((axis.get_major_formatter, axis.set_major_formatter),
                             (axis.get_minor_formatter, axis.set_minor_formatter)):
                f = get()
                if not isinstance(f, _GroupedTicks):
                    put(_GroupedTicks(f))


def _heavier(fig, pt: float) -> None:
    """A main figure's axes, ticks and lines at least `pt` wide (E8.4, the author's
    decision of 2026-11-07, D3): the spines and ticks of every axis, insets included,
    and every drawn line, reference lines, medians, intervals and the key's lines among
    them; a line drawn only as markers keeps its markers' edges, which are not lines."""
    for ax in fig.findobj(matplotlib.axes.Axes):
        for s in ax.spines.values():
            s.set_linewidth(max(s.get_linewidth(), pt))
        for axis in (ax.xaxis, ax.yaxis):
            for tick in axis.get_major_ticks() + axis.get_minor_ticks():
                for line in (tick.tick1line, tick.tick2line):
                    line.set_markeredgewidth(max(line.get_markeredgewidth(), pt))
            axis.set_tick_params(which="both", width=max(pt, 0))
    for line in fig.findobj(Line2D):
        if line.get_linestyle() not in ("None", "", " ") and line.get_linewidth() < pt:
            line.set_linewidth(pt)
    for col in fig.findobj(matplotlib.collections.LineCollection):
        col.set_linewidths([max(w, pt) for w in col.get_linewidths()])


# ---------------------------------------------------------------------------
# Figure: the published cut points, seven genes (a, b); the local likelihood ratio of
# one predictor along its score, in one distance band (c)
# ---------------------------------------------------------------------------
# (a) and (b) keep the block they filled before (c) was added (E6, 2026-11-04): AB_MM
# tall at the top of the figure, every mark in it as far from the top as it was, the
# axes the same size. (c) and its key take about C_MM below, rounded to whole pixels at
# the artwork's resolution: the bitmap's rows are counted from the top, so (a) and (b)
# then fall on the pixel grid they fell on before, pixel for pixel. (c) is one axis, so
# that the figure has three: read strictly, the target's rule on composite figures counts
# each axis as a figure and a composite of more than three as two (the author's decision
# of 2026-11-04).
AB_MM = 92
C_MM = 60
C_AXES_MM = {"left": 16.0, "right": 161.0, "bottom": 15.0, "top": 55.0}   # from the figure's lower left
A_LETTER_X = -0.42                 # (a)'s panel letter, in (a)'s width; (c)'s stands in its column
C_FLOOR, C_CEIL = 0.05, 1000       # Supplementary Figure S1's axis (evid_figures.figure_s1)
C_SHADE = {"color": F.MUTED, "alpha": 0.27, "lw": 0}
C_CUT = {"color": F.INK, "lw": 0.9, "ls": "-."}
C_MARK = {"marker": "|", "ms": 7, "mew": 1.2, "color": F.INK, "ls": ""}


def _part_ref(cfg: dict, content: str, part: str) -> str:
    """An item's part as the text names it, from the configuration: 'Table 2(b)'."""
    it = next(i for i in cfg["main"] + cfg["supplementary"] if i["content"] == content)
    return f"{it['number']}({'abcdefghij'[[p['part'] for p in it['parts']].index(part)]})"


def _score_axis(tool: str) -> str:
    """A predictor's score as an axis label: SpliceAI score (published basis)."""
    name = F.NAME[tool]
    if name.endswith("score"):
        return name
    head, sep, tail = name.partition(" (")
    return f"{head} score{sep}{tail}"


def _local_lr_panel(fig, item: dict, cfg: dict, height: float, ax_a) -> None:
    """(c): the options' predictor and band, its local likelihood ratio along the score,
    as Supplementary Figure S1 draws it for them: the same file read by the same
    function (evid_figures._curve), the same floor and ceiling, a value below the floor
    (zero among them) drawn on it and one above the ceiling at it. In black and grey:
    the colours of (a) and (b) mark ClinVar record status, and (c) draws every variant.
    On each tier's line, a mark at the predictor's in-sample PP3 threshold for that tier,
    where one was fitted (evidence_thresholds.csv, at full precision, as S1 reads it)."""
    tool, st = item["options"]["c_tool"], item["options"]["c_stratum"]
    needs = [F.REPORT_DIR / f"interval_lr_{tool}_{st}.csv", F.REPORT_DIR / "evidence_thresholds.csv"]
    unlisted = [str(p) for p in needs if str(p) not in item["sources"]]
    if unlisted:
        raise SystemExit(f"[layout] {item['content']}: (c) reads {unlisted}, which its sources "
                         f"in config/layout.yaml do not list")
    c = F._curve(tool, st)
    thr = F._thresholds()
    thr = thr[(thr.tool == tool) & (thr.stratum == st)]
    walker = yaml.safe_load(Path(P.CONFIG_PATH).read_text())
    cut = float(walker["thresholds"]["pp3"]["value"])
    width = cfg["artwork"]["double_column_mm"]
    box = C_AXES_MM
    ax = fig.add_axes([box["left"] / width, box["bottom"] / height,
                       (box["right"] - box["left"]) / width, (box["top"] - box["bottom"]) / height])
    ax.set_yscale("log")
    ax.set_ylim(C_FLOOR, C_CEIL)
    ax.set_xlim(0, 1)
    F._tier_lines(ax, F.PATH_CUTS, vertical=False)
    ax.axhline(1.0, color=F.GRID, lw=0.6, zorder=0)
    lo = c["lr_lo_gene_boot"].clip(lower=C_FLOOR).to_numpy()
    hi = c["lr_hi_gene_boot"].clip(upper=C_CEIL).to_numpy()
    ax.fill_between(c.score, lo, hi, zorder=1, **C_SHADE)
    ax.plot(c.score, c.local_lr.clip(lower=C_FLOOR, upper=C_CEIL), color=F.INK, lw=1.0, zorder=2)
    ax.axvline(cut, zorder=1, **C_CUT)
    for name, v in F.PATH_CUTS:
        r = thr[thr.tier == name.lower()]
        if len(r):
            ax.plot(float(r.pp3_threshold_insample.iloc[0]), v, zorder=4, **C_MARK)
    F._logfmt(ax, "y")
    ax.set_xlabel(_score_axis(tool))
    ax.set_ylabel("Local likelihood ratio")
    # the band, as (a) and (b) name their bands
    ax.text(0.012, 0.965, F.STRATUM_LABEL[st], transform=ax.transAxes, fontsize=7,
            fontweight="bold", color=F.INK, ha="left", va="top")
    pa, pc = ax_a.get_position(), ax.get_position()
    F._panel_letter(ax, "c", x=(pa.x0 + A_LETTER_X * pa.width - pc.x0) / pc.width)
    handles = [Line2D([], [], color=F.INK, lw=1.0, label="Local likelihood ratio"),
               Patch(label="5th–95th percentile, gene-clustered bootstrap", **C_SHADE),
               Line2D([], [], label=f"Published PP3 cut point, {cut:g}", **C_CUT),
               Line2D([], [], label="In-sample PP3 threshold "
                      f"({_part_ref(cfg, 'evidence_tiers_and_thresholds', 'thresholds')})", **C_MARK)]
    fig.legend(handles=handles, loc="lower center", ncol=len(handles), bbox_to_anchor=(0.5, 0.0),
               handletextpad=0.4, columnspacing=1.6)


@maker
def figure_published_cut_point(item: dict, folder: Path, cfg: dict) -> list[Path]:
    """v3.0.1's Figure 1 without the external genes, which the figure of the held-out
    and external genes draws instead (a, b), and below it the local likelihood ratio
    along the score (c)."""
    F._style()
    d = K.read_back(item["sources"][0])
    seven = d[(d.source == "seven_genes") & (d.gene == "seven_genes")]
    rows = []
    y = 0.0
    for st in F.STRATA:
        rows.append((y, None, F.STRATUM_LABEL[st], None))
        y += 1
        for arm in F.ARMS:
            r = seven[(seven.stratum == st) & (seven.arm == arm)]
            if len(r):
                rows.append((y, F.ARM_LABEL[arm], None, r.iloc[0]))
                y += 1
        y += 0.4
    ymax = y - 0.4
    dpi = cfg["artwork"]["dpi"]
    c_mm = round(C_MM / 25.4 * dpi) * 25.4 / dpi
    height = AB_MM + c_mm
    fig, axes = plt.subplots(1, 2, figsize=(_width(cfg), height * F.MM), sharey=True,
                             gridspec_kw={"wspace": 0.08})
    specs = [("lr_pp3", "lr_pp3_lo", "lr_pp3_hi", F.PATH_CUTS, (0.8, 60),
              "Likelihood ratio of the PP3 band (score ≥ 0.2)"),
             ("lr_le01", "lr_le01_lo", "lr_le01_hi", F.BEN_CUTS, (0.02, 1.25),
              "Likelihood ratio of the BP4 band (score ≤ 0.1);\nsmaller values are stronger benign evidence")]
    zero = False
    for ax, (col, lo, hi, cuts, xlim, xlabel), letter in zip(axes, specs, "ab"):
        ax.set_xscale("log")
        ax.set_xlim(*xlim)
        F._tier_lines(ax, cuts)
        for yy, lab, head, rec in rows:
            if rec is None or not np.isfinite(rec[col]):
                continue
            colour = F.ARM_COLOUR.get(rec["arm"], F.INK)
            v = rec[col]
            if np.isfinite(rec.get(lo, np.nan)) and np.isfinite(rec.get(hi, np.nan)):
                ax.plot([max(rec[lo], xlim[0]), min(rec[hi], xlim[1])], [yy, yy],
                        color=colour, lw=1.0, solid_capstyle="butt", zorder=2)
            if v <= 0:
                zero = True
                ax.plot(xlim[0], yy, marker="o", ms=4.2, mfc="white", mec=colour, mew=0.9,
                        zorder=3, clip_on=False)
            else:
                ax.plot(v, yy, marker="o", ms=4.2, color=colour, mec="white", mew=0.5,
                        zorder=3)
        ax.set_ylim(ymax - 0.4, -0.8)
        ax.set_xlabel(xlabel)
        F._logfmt(ax)
        ax.tick_params(axis="y", length=0)
        F._panel_letter(ax, letter, x=-0.01 if letter == "b" else A_LETTER_X)
    ax0 = axes[0]
    ax0.set_yticks([r[0] for r in rows if r[1] is not None])
    ax0.set_yticklabels([r[1] for r in rows if r[1] is not None])
    for yy, lab, head, rec in rows:
        if head:
            ax0.text(-0.5, yy, head, transform=ax0.get_yaxis_transform(), fontsize=7,
                     fontweight="bold", color=F.INK, va="center", ha="left")
    handles = [Line2D([], [], marker="o", ls="", color=F.ARM_COLOUR[a], mec="white",
                      label=F.ARM_LABEL[a], ms=4.5) for a in F.ARMS]
    if zero:
        handles += [Line2D([], [], marker="o", ls="", mfc="white", mec=F.INK2,
                           label="Ratio of zero (drawn at the axis floor)", ms=4.5)]
    # the key and the margins (a) and (b) had in a figure AB_MM tall, kept from the top:
    # the key's foot on the block's lower edge
    fig.legend(handles=handles, loc="lower center", ncol=len(handles),
               bbox_to_anchor=(0.6, c_mm / height), handletextpad=0.3, columnspacing=1.2)
    fig.subplots_adjust(left=0.2, right=0.97, top=1 - (1 - 0.94) * AB_MM / height,
                        bottom=1 - (1 - 0.2) * AB_MM / height)
    _local_lr_panel(fig, item, cfg, height, axes[0])
    _heavier(fig, cfg["artwork"]["main_lines_pt"])
    return _save(fig, folder, item["content"], cfg)


# ---------------------------------------------------------------------------
# Figure: held-out genes and external genes
# ---------------------------------------------------------------------------
# genes by marker; TP53's assay differs from the HAP1 deposits (tp53_raw/
# scoreset_meta.json: HCT116 cells, Nutlin-3a), so it alone takes the second colour.
# Every shape has an open form, drawn where a ratio is below its boundary (TP53 had a
# cross until the author's review of 2026-10-28)
EXTERNAL_GENES = ["DDX3X", "CTCF", "SFPQ", "RAD51D", "XRCC2", "SBDS", "TP53"]
GENE_MARKER = {"DDX3X": "o", "CTCF": "s", "SFPQ": "D", "RAD51D": "^", "XRCC2": "v",
               "SBDS": "p", "TP53": "h"}
OPEN_SHAPES = {"o", "s", "D", "^", "v", "p", "h"}
# sizes that give the shapes about the same weight
MARKER_SIZE = {"o": 3.6, "s": 3.3, "D": 3.0, "^": 4.0, "v": 4.0, "p": 4.0, "h": 4.0}
GENE_COLOUR = {g: (F.ORANGE if g == "TP53" else F.BLUE) for g in EXTERNAL_GENES}
TP53_ASSAY = "HCT116 cells, Nutlin-3a"
# the external genes' rows, the author's layout B (2026-10-28): SpliceAI on the published
# basis at the published cut point, judged against the Moderate boundary, and SpliceAI
# and Pangolin at the seven genes' fitted Moderate and Strong thresholds
EXTERNAL_ROWS = [("spliceai_walker", "walker"), ("spliceai_walker", "moderate"),
                 ("spliceai_walker", "strong"), ("pangolin", "moderate"),
                 ("pangolin", "strong")]
TIER_CUT = {"walker": 4.33, "moderate": 4.33, "strong": 18.7}
EXTERNAL_CUTS = [("Moderate", 4.33), ("Strong", 18.7)]


# the held-out folds Table 2(a) does not count, as the pathogenic-side thresholds
# (thresholds_pathogenic) mark them: the combined Atlas score's with BRCA1 or RAD51C held
# out (MARK_RULES); drawn where they were, in grey, and left out of the median (the
# author's decision of 2026-11-05)
FOLD_NOT_COUNTED = "avi_seen_in_training"
CROSS, CROSS_S = "X", 12              # its marker, and its area beside the dots' 9


def _heldout_panel(ax, folds: pd.DataFrame, thr: pd.DataFrame, tier: str, cut: float,
                   st: str, floor: float, ceil: float) -> tuple[bool, bool]:
    """v3.0.1's Figure 3, one panel: each held-out gene's ratio at the fitted threshold,
    and the median of the folds counted. A fold not counted is drawn in grey by the same
    rule, filled where it clears the boundary, open where it does not. Returns whether a
    ratio of zero, and a fold not counted, were drawn."""
    ax.set_xscale("log")
    ax.set_xlim(floor, ceil)
    ax.axvline(cut, color=F.INK2, lw=0.8, ls=(0, (1.5, 1.5)), zorder=0)
    any_zero = any_grey = False
    rule = MARK_RULES[FOLD_NOT_COUNTED]
    for yy, c in enumerate(F.COLUMNS):
        if not len(thr[(thr.tool == c) & (thr.stratum == st) & (thr.tier == tier)]):
            ax.text(floor * 1.15, yy, "not reached in-sample", fontsize=5.5, color=F.MUTED,
                    va="center", zorder=5, bbox=dict(facecolor="white", edgecolor="none", pad=0.4))
            continue
        f = folds[(folds.tool == c) & (folds.stratum == st) & (folds.tier == tier)]
        f = f[np.isfinite(f["heldout_lr"])]
        v = f["heldout_lr"].to_numpy()
        if not len(v):
            ax.text(floor * 1.15, yy, "no evaluable fold", fontsize=5.5, color=F.MUTED,
                    va="center", zorder=5, bbox=dict(facecolor="white", edgecolor="none", pad=0.4))
            continue
        counted = ~np.array([bool(rule(r)) for _, r in f.iterrows()], dtype=bool)
        jitter = (np.arange(len(v)) - (len(v) - 1) / 2) * 0.09
        above, zero = v >= cut, v <= 0
        pos = np.clip(v, floor, ceil)
        draw = counted & ~zero
        ax.scatter(pos[above & draw], yy + jitter[above & draw], s=9, color=F.BLUE,
                   edgecolor="white", linewidth=0.3, zorder=3)
        ax.scatter(pos[~above & draw], yy + jitter[~above & draw], s=9, facecolor="white",
                   edgecolor=F.BLUE, linewidth=0.7, zorder=3)
        if (zero & counted).any():
            any_zero = True
            ax.scatter(np.full((zero & counted).sum(), floor), yy + jitter[zero & counted], s=9,
                       facecolor="white", edgecolor=F.INK2, linewidth=0.7, zorder=3, clip_on=False)
        if (~counted).any():
            # a fold not counted is a grey cross, told from the counted dots by its shape
            # in colour and in grey alike (the author's decision of 2026-11-07, D4):
            # filled where it clears the boundary, open where it does not, as the dots
            any_grey = True
            grey = ~counted & ~zero
            ax.scatter(pos[above & grey], yy + jitter[above & grey], marker=CROSS, s=CROSS_S,
                       color=F.MUTED, edgecolor="white", linewidth=0.3, zorder=3)
            ax.scatter(pos[~above & grey], yy + jitter[~above & grey], marker=CROSS, s=CROSS_S,
                       facecolor="white", edgecolor=F.MUTED, linewidth=0.7, zorder=3)
            if (zero & ~counted).any():
                any_zero = True
                ax.scatter(np.full((zero & ~counted).sum(), floor), yy + jitter[zero & ~counted],
                           marker=CROSS, s=CROSS_S, facecolor="white", edgecolor=F.MUTED,
                           linewidth=0.7, zorder=3, clip_on=False)
        if counted.any():
            ax.plot([np.median(pos[counted])] * 2, [yy - 0.32, yy + 0.32], color=F.INK, lw=1.0,
                    zorder=4)
    F._logfmt(ax)
    ax.tick_params(axis="y", length=0)
    return any_zero, any_grey


def _external_tables() -> dict[str, pd.DataFrame]:
    out = {}
    for g in EXTERNAL_GENES:
        d = K.read_back(REPORT_DIR / f"external_{g.lower()}.csv")
        out[g] = d[d.label_definition == X.PRIMARY_LABEL[g]]
    return out


def _external_value(r: pd.Series, kind: str) -> float:
    col = "walker_lr_pp3" if kind == "walker" else f"pp3_insample_{kind}_lr_here"
    return float(r[col]) if pd.notna(r[col]) else np.nan


def _band_count(r: pd.Series, kind: str) -> int:
    """Labelled variants at or above the row's threshold."""
    p = "walker_pp3" if kind == "walker" else f"pp3_insample_{kind}"
    return int(r[f"{p}_n_pos_band"] + r[f"{p}_n_neg_band"])


def _reachable(r: pd.Series, kind: str) -> bool:
    if kind == "walker":
        return True
    return str(r[f"pp3_insample_{kind}_reachable"]).strip().lower() in ("true", "1", "1.0")


def _cut_lines(ax, cuts: list[tuple[str, float]]) -> None:
    """The boundaries, dotted and named at the top as Figure 1 names its: the lower name
    ends at its line and the upper starts at its own, so that the two do not meet."""
    for k, (name, v) in enumerate(cuts):
        ax.axvline(v, color=F.MUTED, lw=0.6, ls=(0, (1.5, 1.5)), zorder=0)
        ax.text(v, 1.0, name, transform=ax.get_xaxis_transform(), fontsize=6, color=F.INK2,
                ha="right" if k == 0 else "left", va="bottom")


def _external_panel(ax, ext: dict, st: str, rows: list[tuple[str, str]], floor: float,
                    ceil: float) -> list[dict]:
    """One band: a row per (column, threshold); each external gene a marker of its own,
    at a fixed height within the row; filled where the ratio reaches the row's boundary,
    the published cut point's row judged against Moderate. Returns the genes drawn without
    a marker where the row has a threshold, and why."""
    ax.set_xscale("log")
    ax.set_xlim(floor, ceil)
    _cut_lines(ax, EXTERNAL_CUTS)
    slot = {g: (i - (len(EXTERNAL_GENES) - 1) / 2) * 0.1 for i, g in enumerate(EXTERNAL_GENES)}
    missing = []
    for yy, (tool, kind) in enumerate(rows):
        drawn, reachable = 0, False
        for g in EXTERNAL_GENES:
            r = ext[g][(ext[g].stratum == st) & (ext[g].tool == tool)]
            if not len(r):
                continue
            r = r.iloc[0]
            here = _reachable(r, kind)
            reachable |= here
            if not here:
                continue
            why = None
            if r["status"] != "ok":
                why = (f"{r['reason']} ({int(r['n_pos'])} damaging, {int(r['n_neg'])} normal)")
            elif not bool(r["column_basis_comparable"]):
                why = "the column is not the analysis set's variable"
            v = _external_value(r, kind)
            if why is None and not np.isfinite(v) and not np.isinf(v):
                why = (f"{_band_count(r, kind)} labelled variants at or above the threshold, "
                       "too few for a ratio")
            if why:
                missing.append({"gene": g, "stratum": st, "tool": tool, "kind": kind,
                                "status_ok": r["status"] == "ok", "why": why})
                continue
            clears = v >= TIER_CUT[kind]
            x = floor if v <= 0 else min(max(v, floor), ceil)
            # the combined Atlas score's model was selected on the DDX3X assay: not an
            # external test there, drawn in grey as v3.0.1's Figure 4 drew it
            colour = F.MUTED if (g, tool) in F.NOT_EXTERNAL else GENE_COLOUR[g]
            ax.plot(x, yy + slot[g], marker=GENE_MARKER[g], ms=MARKER_SIZE[GENE_MARKER[g]], ls="",
                    mfc=colour if clears and v > 0 else "white", mec=colour if not clears or v <= 0 else "white",
                    mew=0.7 if not clears or v <= 0 else 0.3, zorder=3, clip_on=False)
            drawn += 1
        if not drawn:
            ax.text(floor * 1.15, yy, "not reached in-sample" if not reachable else "no evaluable gene",
                    fontsize=5.5, color=F.MUTED, va="center", zorder=5,
                    bbox=dict(facecolor="white", edgecolor="none", pad=0.4))
    F._logfmt(ax)
    ax.tick_params(axis="y", length=0)
    return missing


def _row_label(tool: str, kind: str) -> str:
    """A row's name: the column and the threshold."""
    if kind == "walker":
        return "SpliceAI, published 0.2"
    name = {"spliceai_walker": "SpliceAI", "pangolin": "Pangolin"}.get(tool, F.NAME[tool])
    return f"{name}, fitted {kind.capitalize()}"


def _missing_notes(missing: list[dict]) -> list[str]:
    """The genes without a marker, for the caption: a gene a band cannot evaluate, once
    for the band; a ratio missing at one threshold, for that row."""
    notes, seen = [], set()
    for m in missing:
        band = F.STRATUM_LABEL[m["stratum"]]
        if not m["status_ok"]:
            key = (m["gene"], m["stratum"])
            if key not in seen:
                seen.add(key)
                notes.append(f"no marker: {m['gene']} at {band}, every row: {m['why']}")
        else:
            notes.append(f"no marker: {m['gene']} at {band}, {_row_label(m['tool'], m['kind'])}: "
                         f"{m['why']}")
    return notes


@maker
def figure_heldout_external(item: dict, folder: Path, cfg: dict) -> list[Path]:
    """v3.0.1's Figure 3 (held-out genes, a and b) with the seven external genes below it
    (c), at the seven genes' in-sample thresholds beside the published cut point:
    SpliceAI on the published basis and Pangolin, a row per threshold."""
    F._style()
    folds = K.read_back(REPORT_DIR / "evidence_thresholds_logo_folds.csv")
    folds = folds[(folds.side == "pp3") & folds.stratum.isin(F.STRATA) & folds.tool.isin(F.COLUMNS)]
    thr = F._thresholds()
    ext = _external_tables()
    heights = [len(F.COLUMNS), len(F.COLUMNS), len(EXTERNAL_ROWS)]
    fig = plt.figure(figsize=(_width(cfg), 190 * F.MM))
    gs = fig.add_gridspec(len(heights), 3, height_ratios=heights, hspace=0.42, wspace=0.08,
                          left=0.24, right=0.98, top=0.965, bottom=0.145)
    floor, ceil = 0.3, 1000
    any_zero = any_grey = False
    letters = iter("abc")
    for i, (tier, cut) in enumerate((("moderate", 4.33), ("strong", 18.7))):
        axes = [fig.add_subplot(gs[i, j]) for j in range(3)]
        for j, (ax, st) in enumerate(zip(axes, F.STRATA)):
            zero, grey = _heldout_panel(ax, folds, thr, tier, cut, st, floor, ceil)
            any_zero, any_grey = any_zero or zero, any_grey or grey
            ax.set_ylim(len(F.COLUMNS) - 0.5, -0.5)
            if i == 0:
                ax.set_title(F.STRATUM_LABEL[st], fontweight="bold")
            if j == 0:
                ax.set_yticks(range(len(F.COLUMNS)))
                ax.set_yticklabels([F.NAME[c] + (" †" if c == "avi" else "") for c in F.COLUMNS])
            else:
                ax.set_yticks([])
        axes[1].text(0.5, -0.1, f"Held-out gene: likelihood ratio at the fitted {tier.capitalize()} "
                     f"threshold (dotted: the {tier.capitalize()} boundary, {cut:g})",
                     transform=axes[1].transAxes, ha="center", va="top", fontsize=6.5)
        F._panel_letter(axes[0], next(letters), x=-0.8)
    axes = [fig.add_subplot(gs[2, j]) for j in range(3)]
    missing = []
    for j, (ax, st) in enumerate(zip(axes, F.STRATA)):
        missing += _external_panel(ax, ext, st, EXTERNAL_ROWS, floor, ceil)
        ax.set_ylim(len(EXTERNAL_ROWS) - 0.5, -0.5)
        if j == 0:
            ax.set_yticks(range(len(EXTERNAL_ROWS)))
            ax.set_yticklabels([_row_label(t, kind) for t, kind in EXTERNAL_ROWS])
        else:
            ax.set_yticks([])
        # above the boundaries' names
        ax.set_title(F.STRATUM_LABEL[st], fontweight="bold", fontsize=7, pad=9)
    axes[1].text(0.5, -0.3, "External gene: likelihood ratio at the published cut point and "
                 "the fitted Moderate and Strong thresholds",
                 transform=axes[1].transAxes, ha="center", va="top", fontsize=6.5)
    F._panel_letter(axes[0], next(letters), x=-0.8)
    handles = [Line2D([], [], marker="o", ls="", color=F.BLUE, mec="white", ms=4,
                      label="Held-out gene, clears the boundary"),
               Line2D([], [], marker="o", ls="", mfc="white", mec=F.BLUE, ms=4,
                      label="Held-out gene, below it"),
               Line2D([], [], color=F.INK, lw=1.0, label="Median of held-out genes")]
    if any_grey:
        # the row's † names why (Table 2(a)'s †); the median leaves these folds out
        handles.append(Line2D([], [], marker=CROSS, ls="", color=F.MUTED, mec="white", ms=4.6,
                              label=f"Fold not counted ({NOT_COUNTED})"))
    if any_zero:
        handles.append(Line2D([], [], marker="o", ls="", mfc="white", mec=F.INK2, ms=4,
                              label="Ratio of zero (at the axis floor)"))
    genes = [Line2D([], [], marker=GENE_MARKER[g], ls="", color=GENE_COLOUR[g], mec="white",
                    ms=MARKER_SIZE[GENE_MARKER[g]] + 0.6,
                    label=F._it(g) + (f" ({TP53_ASSAY})" if g == "TP53" else ""))
             for g in EXTERNAL_GENES]
    genes.append(Line2D([], [], marker="o", ls="", mfc="white", mec=F.INK2, ms=4,
                        label="Open: below the boundary"))
    # the key is filled column by column: the held-out entries fill the first column and
    # the genes the columns after it (four columns for three entries, as before E7)
    fig.legend(handles=handles + genes, loc="lower center",
               ncol=1 + -(-len(genes) // len(handles)),
               bbox_to_anchor=(0.55, 0.0), handletextpad=0.3, columnspacing=1.0)
    item["_notes"] = _missing_notes(missing)
    _heavier(fig, cfg["artwork"]["main_lines_pt"])
    return _save(fig, folder, item["content"], cfg)


# ---------------------------------------------------------------------------
# Figures: ROC curves, precision-recall curves (v3.0.1's practical figure, one row each)
# ---------------------------------------------------------------------------
ROC_LETTER_X = -0.06                  # the ROC panels' letters, in each panel's width


def _curves_figure(item: dict, folder: Path, cfg: dict, row: str) -> list[Path]:
    F._style()
    roc = K.read_back(REPORT_DIR / "practical_roc_curves.csv")
    pr = K.read_back(REPORT_DIR / "practical_pr_curves.csv")
    auc = K.read_back(REPORT_DIR / "practical_curve_auc.csv")
    table = K.read_back(REPORT_DIR / "practical_metrics_by_threshold.csv")
    walker = yaml.safe_load(Path(P.CONFIG_PATH).read_text())
    pp3_cut = float(walker["thresholds"]["pp3"]["value"])
    marks = P._tier_marks(table)
    xmax, ymax = P._inset_range(roc, marks)
    unit = [0, 0.25, 0.5, 0.75, 1]
    fig, axes = plt.subplots(1, 3, figsize=(_width(cfg), (72 if row == "roc" else 92) * F.MM),
                             gridspec_kw={"wspace": 0.28})

    def nudged(ax, tool):
        dx, dy = P.MARKER_NUDGE[tool]
        return mtransforms.offset_copy(ax.transData, fig=fig, x=dx, y=dy, units="points")

    for j, st in enumerate(P.CURVE_STRATA):
        ax = axes[j]
        a = auc[auc.stratum == st].set_index("tool")
        inset = ax.inset_axes(P.INSET_BOUNDS) if row == "roc" else None
        if row == "roc":
            ax.plot([0, 1], [0, 1], color=F.GRID, lw=0.6, ls=(0, (1.5, 1.5)), zorder=0)
        else:
            prevalence = float(a["prevalence"].iloc[0])
            ax.axhline(prevalence, color=F.INK2, lw=0.7, ls=(0, (3, 1.5)), zorder=0)
            ax.text(0.04, prevalence, f"prevalence {prevalence:.2f}", fontsize=5.8,
                    color=F.INK2, ha="left", va="bottom")
        handles = []
        for tool in P.MAIN_TOOLS:
            colour, dash = P.TOOL_COLOUR[tool], P.TOOL_DASH[tool]
            r = roc[(roc.tool == tool) & (roc.stratum == st)]
            p = pr[(pr.tool == tool) & (pr.stratum == st)]
            if row == "roc":
                ax.plot(r.fpr, r.tpr, color=colour, lw=1.0, ls=dash, zorder=2)
                inset.plot(r.fpr, r.tpr, color=colour, lw=0.8, ls=dash, zorder=2)
            else:
                ax.plot(p.recall.iloc[1:], p.precision.iloc[1:], color=colour, lw=1.0, ls=dash,
                        drawstyle="steps-pre", zorder=2)
            for m in marks[(marks.tool == tool) & (marks.stratum == st)].itertuples():
                x, yv, rc, pc = P._point(roc, pr, tool, st, m.threshold)
                shape, size = P.TIER_MARKER[m.threshold_kind]
                style = dict(marker=shape, ms=size, color=colour, mec="white", mew=0.4, ls="",
                             zorder=4, clip_on=False)
                if row == "roc":
                    ax.plot(x, yv, transform=nudged(ax, tool), **style)
                    inset.plot(x, yv, transform=nudged(inset, tool), **style)
                else:
                    ax.plot(rc, pc, transform=nudged(ax, tool), **style)
            if tool == "spliceai_walker":
                x, yv, rc, pc = P._point(roc, pr, tool, st, pp3_cut)
                style = dict(marker="D", ms=3.8, mfc="white", mec=colour, mew=0.9, ls="", zorder=5)
                if row == "roc":
                    ax.plot(x, yv, **style)
                    if x <= xmax:
                        inset.plot(x, yv, **style)
                else:
                    ax.plot(rc, pc, **style)
            # the area under the curve drawn, the seven genes pooled into one curve; not
            # the text's AUROC, which pools per-gene AUROCs (the author's decision of
            # 2026-10-28; the supplementary table of practical metrics gives both)
            key = "auroc_pooled_curve" if row == "roc" else "auprc_pooled_curve"
            handles.append(Line2D([], [], color=colour, lw=1.0, ls=dash,
                                  label=f"{P.SHORT[tool]}  {a.loc[tool, key]:.2f}"))
        if row == "roc":
            ax.legend(handles=handles, loc="lower right", title="Area under the curve",
                      alignment="left", fontsize=5.8, title_fontsize=5.8, handlelength=2.2,
                      handletextpad=0.4, borderaxespad=0.2, labelspacing=0.25)
        else:
            # the curves fill the top of each panel, so the key sits below it
            ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.3), ncol=2,
                      title="Area under the curve", fontsize=5.8, title_fontsize=5.8,
                      handlelength=2.2, handletextpad=0.4, columnspacing=0.8, labelspacing=0.25)
        if row == "roc":
            inset.set_xlim(0, xmax)
            inset.set_ylim(0, ymax)
            xt = [round(v, 10) for v in np.arange(0, xmax + P.INSET_X_STEP / 2, P.INSET_X_STEP)]
            inset.set_xticks(xt)
            inset.set_xticklabels([f"{v:g}" for v in xt], fontsize=5.5)
            yt = [0, round(ymax / 2, 10), ymax]
            inset.set_yticks(yt)
            inset.set_yticklabels([f"{v:g}" for v in yt], fontsize=5.5)
            inset.tick_params(length=1.5, pad=1)
            for spine in inset.spines.values():
                spine.set_visible(True)
                spine.set_color(F.MUTED)
                spine.set_linewidth(0.5)
        ax.set_xlim(-0.01, 1.01)
        ax.set_ylim(-0.01, 1.01)
        ax.set_xticks(unit)
        ax.set_yticks(unit)
        ax.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"])
        ax.set_yticklabels(["0", "0.25", "0.5", "0.75", "1"])
        ax.set_aspect("equal")
        ax.set_title(F.STRATUM_LABEL[st], fontweight="bold")
        if row == "roc":
            F._panel_letter(ax, "abc"[j], x=ROC_LETTER_X)
        ax.set_xlabel("1 − specificity" if row == "roc" else "Sensitivity (recall)")
        if j == 0:
            ax.set_ylabel("Sensitivity" if row == "roc" else "Precision")
    tools = [Line2D([], [], color=P.TOOL_COLOUR[t], lw=1.0, ls=P.TOOL_DASH[t],
                    label=F.NAME[t] + (" †" if t == "avi" else "")) for t in P.MAIN_TOOLS]
    tier_marks = [Line2D([], [], marker=P.TIER_MARKER[t][0], ms=P.TIER_MARKER[t][1], ls="",
                         color=F.INK2, mec="white", mew=0.4,
                         label=f"Fitted {P.TIER_NAME[t]} threshold")
                  for t in ("supporting", "moderate", "strong")]
    tier_marks += [Line2D([], [], marker="D", ms=3.8, ls="", mfc="white", mec=F.INK2, mew=0.9,
                          label=f"Published cut point, {pp3_cut:g}")]
    handles = [h for pair in zip(tools, tier_marks) for h in pair]
    if row == "pr":
        handles.append(Line2D([], [], color=F.INK2, lw=0.7, ls=(0, (3, 1.5)), label="Prevalence"))
    fig.legend(handles=handles, loc="lower center", ncol=5 if row == "pr" else 4,
               bbox_to_anchor=(0.5, 0.0), handletextpad=0.4, columnspacing=1.0, handlelength=2.4)
    fig.subplots_adjust(left=0.07, right=0.99, top=0.92 if row == "roc" else 0.94,
                        bottom=0.25 if row == "roc" else 0.4)
    if row == "roc":
        _heavier(fig, cfg["artwork"]["main_lines_pt"])
    return _save(fig, folder, item["content"], cfg)


@maker
def figure_roc(item: dict, folder: Path, cfg: dict) -> list[Path]:
    return _curves_figure(item, folder, cfg, "roc")


@maker
def figure_precision_recall(item: dict, folder: Path, cfg: dict) -> list[Path]:
    return _curves_figure(item, folder, cfg, "pr")


# ---------------------------------------------------------------------------
# Figure: ClinVar coverage by position (C2), the final form of E19's draft
# ---------------------------------------------------------------------------
CLASS_COLOUR = [("plp", "P/LP", F.BLUE), ("blb", "B/LB", F.ORANGE), ("vus", "VUS", F.AQUA),
                ("conflicting", "Conflicting", "#eda100"), ("other", "Other", F.MUTED)]


@maker
def figure_clinvar_coverage(item: dict, folder: Path, cfg: dict) -> list[Path]:
    F._style()
    d = K.read_back(item["sources"][0])
    d = d[d.gene == "Seven genes"].sort_values("signed_offset")
    fig, axes = plt.subplots(2, 1, figsize=(_width(cfg), 110 * F.MM), sharex=True,
                             gridspec_kw={"hspace": 0.18})
    x = d.signed_offset.to_numpy()
    for ax in axes:
        for b in (-10.5, -2.5, 2.5, 10.5):
            ax.axvline(b, color=F.MUTED, lw=0.6, ls=(0, (1.5, 1.5)), zorder=0)
        ax.axvspan(-2.5, 2.5, color=F.BAND, zorder=-1, lw=0)
        ax.set_ylim(0, 1)
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
    a, b = axes
    a.plot(x, d.share_recorded, color=F.INK, lw=1.0, marker="o", ms=2.2, label="Any ClinVar record")
    a.plot(x, d.share_classified, color=F.BLUE, lw=1.0, ls=(0, (4.0, 1.4)), marker="s", ms=2.0,
           label="Classified P/LP or B/LB")
    a.set_ylabel("Share of SGE variants")
    a.legend(loc="upper right")
    bottom = np.zeros(len(d))
    for key, label, colour in CLASS_COLOUR:
        v = d[f"share_{key}"].to_numpy()
        b.bar(x, v, bottom=bottom, width=0.8, color=colour, label=label, lw=0)
        bottom += v
    b.set_ylabel("Share of SGE variants")
    # distances in bp, as the text gives them (E8.10, the author's decision of 2026-11-09)
    b.set_xlabel("Signed intronic offset (acceptor side < 0 < donor side), bp")
    b.legend(loc="upper right", ncol=5)
    for ax, letter in zip(axes, "ab"):
        F._panel_letter(ax, letter, x=-0.06)
        ax.text(0, 1.0, "1–2", transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                fontsize=5.8, color=F.INK2)
        for lo, hi, lab in ((-10, -3, "3–10"), (3, 10, "3–10"), (-50, -11, "11–50"), (11, 50, "11–50")):
            ax.text((lo + hi) / 2, 1.0, lab, transform=ax.get_xaxis_transform(), ha="center",
                    va="bottom", fontsize=5.8, color=F.INK2)
    fig.subplots_adjust(left=0.08, right=0.99, top=0.95, bottom=0.1)
    return _save(fig, folder, item["content"], cfg)


# ---------------------------------------------------------------------------
# Figures: v3.0.1's supplementary figures, drawn again at this layout's width
# ---------------------------------------------------------------------------
@contextmanager
def _drawn_by_figures(width: float, kept: list):
    """evid_figures' functions draw at its own double-column width and save under
    figures/; while this holds they draw at this layout's width and hand their figure
    to `kept` instead of saving it."""
    double, save = F.DOUBLE, F._save
    F.DOUBLE, F._save = width, (lambda fig, stem: kept.append(fig))
    try:
        yield
    finally:
        F.DOUBLE, F._save = double, save


def caption_gap_mm(fig, dpi: float) -> float:
    """The space between a figure's caption line, the lowest of its own texts, and the
    tick labels of its panels above it, in mm: from the top of the caption's box to the
    foot of the lowest label's, each box spanning its font's full ascent and descent, as
    a PDF reader boxes a word. That is the most such a measure takes off the blank
    between them; matplotlib's own text boxes, and the ink, leave more. Text is laid out
    a little differently in the bitmap, at `dpi`, and in the PDF, by the PDF backend's
    metrics; the lesser of the two is returned."""
    from matplotlib.backends.backend_pdf import PdfFile, RendererPdf
    from matplotlib.font_manager import findfont, get_font
    screen, gaps = fig.dpi, []
    for at, pdf in ((dpi, False), (72, True)):
        fig.set_dpi(at)
        try:
            fig.canvas.draw()                        # places the ticks and their labels
            r = fig.canvas.get_renderer()
            if pdf:
                w, h = fig.get_size_inches()
                r = RendererPdf(PdfFile(io.BytesIO()), 72, h, w)

            def span(t) -> tuple[float, float]:
                """A one-line text's font box, foot and top, in display units:
                matplotlib boxes it from its baseline less the descent of the text or of
                "lp"."""
                prop = t.get_fontproperties()
                font = get_font(findfont(prop))
                em = r.points_to_pixels(prop.get_size_in_points())
                d = max(r.get_text_width_height_descent(s, prop, ismath=False)[2]
                        for s in (t.get_text(), "lp"))
                base = t.get_window_extent(r).y0 + d
                return (base + font.descender / font.units_per_EM * em,
                        base + font.ascender / font.units_per_EM * em)

            labels = []
            for ax in fig.axes:
                lo, hi = sorted(ax.get_xlim())
                labels += [tick.label1 for tick, x in zip(ax.xaxis.get_major_ticks(),
                                                           ax.xaxis.get_majorticklocs())
                           if lo - 1e-9 <= x <= hi + 1e-9 and tick.label1.get_visible()
                           and tick.label1.get_text().strip()]
            cap = min(fig.texts, key=lambda t: t.get_window_extent(r).y0)
            # display y runs upwards: the labels' lowest foot less the caption's top
            gaps.append((min(span(t)[0] for t in labels) - span(cap)[1]) / at * 25.4)
        finally:
            fig.set_dpi(screen)
    return min(gaps)


@maker
def redraw_figure(item: dict, folder: Path, cfg: dict) -> list[Path]:
    """A figure of v3.0.1, by its own function in evid_figures, from the same outputs; at
    the height the options give, if they give one. The figure is set to that height
    before it is saved: its panels are placed in fractions of the figure, so they become
    shorter and its type keeps its size, where scaling the page would shrink the type.
    With a caption gap in the options, the panels' lower edge is raised until the
    caption line at the foot stands that far below their lowest text: every panel and
    the space between them shorter in proportion, the height and the type as they were."""
    F._style()
    kept: list = []
    with _drawn_by_figures(_width(cfg), kept):
        getattr(F, item["options"]["function"])()
    if len(kept) != 1:
        raise SystemExit(f"[layout] {item['content']}: evid_figures.{item['options']['function']} "
                         f"drew {len(kept)} figures, not one")
    fig = kept[0]
    if "height_mm" in item["options"]:
        fig.set_size_inches(_width(cfg), item["options"]["height_mm"] * F.MM)
    want = item["options"].get("caption_gap_mm")
    if want is not None:
        dpi = cfg["artwork"]["dpi"]
        gap = caption_gap_mm(fig, dpi)
        if gap < want:
            raise_by = (want - gap) / (fig.get_figheight() / F.MM)
            # rounded up at the fourth decimal, so that the gap is not short of it by a hair
            fig.subplots_adjust(bottom=np.ceil((fig.subplotpars.bottom + raise_by) * 1e4) / 1e4)
            gap = caption_gap_mm(fig, dpi)
        if gap < want - 1e-9:
            raise SystemExit(f"[layout] {item['content']}: the caption line is {gap:.2f} mm below "
                             f"the panels' lowest text, less than {want} mm")
    return _save(fig, folder, item["content"], cfg)


# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------
# the printed tables, which are written as they stand, byte for byte
PRINTED = (str(REPORT_DIR / "tables") + "/", str(REPORT_DIR / "supplement") + "/")
NUMERIC = {"count", "proportion", "lr", "threshold", "score", "percent", "number", "constant"}
CODED = {"tool", "stratum", "tier", "side", "label", "boolean"}
# the types the code prints; a column can also name one of the configuration's maps of
# codes to names (display: values:)
TYPES = {"text"} | NUMERIC | CODED
# the codes the study's tables print as the tables v3.0.1 printed print them, beyond
# evid_supp_tables'
TOOL_EXTRA = {"fusion_enet": "Elastic-net combination"}
STRATUM_EXTRA = {"out_of_scope_deep": "Beyond 50 (out of scope)"}
TIER_LABEL = {"supporting": "Supporting", "moderate": "Moderate", "strong": "Strong",
              "very_strong": "Very strong", "below supporting": "Below Supporting"}
# the seven genes' label is each deposit's own classification, which the printed
# in-frame attribution table leaves blank
LABEL_BLANK = {"y_assay"}
# the rules a part's marks: name in the configuration. The combined Atlas score's model
# was selected on the BRCA1 and RAD51C assays (evid_common.AVI_SEEN_IN_TRAINING, the
# rule every held-out count follows), so its rows with one of them held out test
# nothing independent of that choice
MARK_RULES = {
    "avi_seen_in_training": lambda r: r.get("tool") == "avi"
                                      and r.get("heldout_gene") in K.AVI_SEEN_IN_TRAINING,
}
DASH = "–"
NOT_COUNTED = "†"          # a fold listed but left out of its row's counts (thresholds_pathogenic, Figure 2)
NUMBERS_AS_WRITTEN = "numbers as written"     # a map's option: see _cell()


def read_source(path: str | Path) -> tuple[pd.DataFrame, list[str]]:
    """A table as text cells, with the comment lines above its header as notes."""
    return _parse_source(Path(path).read_text())


def _parse_source(text: str) -> tuple[pd.DataFrame, list[str]]:
    lines = text.splitlines(keepends=True)
    notes = []
    while lines and lines[0].startswith("#"):
        notes.append(lines.pop(0).lstrip("#").strip())
    df = pd.read_csv(io.StringIO("".join(lines)), dtype=str, keep_default_na=False)
    return df, notes


def _replacements(part: dict, numbers: dict) -> list[tuple[str, str, int]]:
    """A part's `replace` entries: the source's text, the text printed in its place (its
    {content} references written as the items' numbers), and how often the part holds it."""
    return [(r["old"], _resolve(r["new"], numbers), int(r["count"])) for r in part.get("replace", [])]


def replaced(df: pd.DataFrame, notes: list[str], reps: list[tuple[str, str, int]],
             where_: str) -> tuple[pd.DataFrame, list[str]]:
    """A part's table and its source's notes with each replacement made in the headers,
    the cells and the notes, in order, each counted on the text the ones before it left;
    the run stops where a text is not held exactly as often as its entry says."""
    if not reps:
        return df, notes
    cols, data = list(df.columns), [list(r) for r in df.itertuples(index=False)]
    notes = list(notes)
    for old, new, count in reps:
        held = (sum(c.count(old) for c in cols) + sum(v.count(old) for r in data for v in r)
                + sum(n.count(old) for n in notes))
        if held != count:
            raise SystemExit(f"{where_}: '{old}' is to be printed as '{new}' {count} times, and "
                             f"the part holds it {held} times")
        cols = [c.replace(old, new) for c in cols]
        data = [[v.replace(old, new) for v in r] for r in data]
        notes = [n.replace(old, new) for n in notes]
    if len(set(cols)) != len(cols):
        raise SystemExit(f"{where_}: the replacements give two columns one name")
    return pd.DataFrame(data, columns=cols, dtype=str), notes


def _copy_replaced(src: str | Path, dst: Path, reps: list[tuple[str, str, int]],
                   where_: str) -> tuple[pd.DataFrame, list[str]]:
    """A printed table's copy with its part's replacements: the source's bytes with each
    text replaced where it stands, which must read as the source does with the same
    replacements made cell by cell (replaced()), so that nothing else in the copy differs
    from the source. Returns the copy's table and notes."""
    text = Path(src).read_bytes().decode("utf-8")
    out = text
    for old, new, count in reps:
        if out.count(old) != count:
            raise SystemExit(f"{where_}: '{old}' is to be printed as '{new}' {count} times, and "
                             f"the source holds it {out.count(old)} times")
        out = out.replace(old, new)
    want_df, want_notes = replaced(*_parse_source(text), reps, where_)
    got_df, got_notes = _parse_source(out)
    if not (got_df.equals(want_df) and got_notes == want_notes):
        raise SystemExit(f"{where_}: a replacement in the source's text does not read as the same "
                         "replacement in its cells (it spans a delimiter or a quote)")
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_bytes(out.encode("utf-8"))
    return got_df, got_notes


# ---------------------------------------------------------------------------
# How every table finishes what it prints (E8.3): thousands, transcripts
# ---------------------------------------------------------------------------
# A number of 1,000 or more is printed with commas between its thousands, the target's
# style, a count or a value alike: a cell that is a number, a ratio of counts, a number
# in a text. Not a year (four digits from 1900 to 2099 in a text, unless "of" stands
# next to it, as in "2000 of 2000 draws"), a date or a time, a version, a coordinate or
# an identifier: a number of seven digits or more in a text, or one joined to a letter,
# a point, a colon, a slash, a hyphen, an underscore or a quote (PMID 30209399,
# g.32394937, urn:mavedb:00001250-a-3, 2026-10-19, ENST00000644876.2, 0.10.0).
PLAIN_NUMBER = re.compile(r"^(-?)(\d+)(\.\d+)?$")
NUMBER_IN_TEXT = re.compile(r"(?<![\w.,/:#+*\-−–'])(\d{4,6})(\.\d+)?(?![\w'’]|[.,]\d|\s?[-–/:]\s?\d)")
RATIO_IN_TEXT = re.compile(r"(?<![\w.,/:#+*\-−–])(\d+)/(\d+)(?![\w.,/])")
YEARS = (1900, 2099)


def grouped(text: str) -> str:
    """A printed text with every number of 1,000 or more grouped in thousands, by the
    rule above; any other text as it was."""
    m = PLAIN_NUMBER.match(text)
    if m:
        sign, whole, frac = m.groups()
        return text if len(whole) < 4 else f"{sign}{int(whole):,}{frac or ''}"
    text = RATIO_IN_TEXT.sub(lambda r: f"{int(r.group(1)):,}/{int(r.group(2)):,}", text)

    def one(n):
        whole, frac = n.group(1), n.group(2) or ""
        if (len(whole) == 4 and not frac and YEARS[0] <= int(whole) <= YEARS[1]
                and not text[n.end():].startswith(" of ") and not text[:n.start()].endswith("of ")):
            return n.group(0)                                   # a year
        return f"{int(whole):,}{frac}"
    return NUMBER_IN_TEXT.sub(one, text)


def _mane_select(path) -> dict:
    """The MANE summary's MANE Select rows, by gene symbol."""
    with gzip.open(path, "rt") as f:
        header = f.readline().lstrip("#").rstrip("\n").split("\t")
        rows = [dict(zip(header, line.rstrip("\n").split("\t"))) for line in f]
    return {r["symbol"]: r for r in rows if r["MANE_status"] == "MANE Select"}


def transcripts(df: pd.DataFrame, part: dict, where_: str) -> pd.DataFrame:
    """A part's transcript column (its `transcripts` entry) as RefSeq accession.version,
    read from the MANE summary it names, never written by hand: a gene's MANE Select
    RefSeq prints as it stands; its Ensembl transcript, which a deposit may name instead,
    prints after the RefSeq of the same row, in parentheses. Any other transcript, or a
    gene without a MANE Select row, stops the run."""
    spec = part.get("transcripts")
    if not spec:
        return df
    mane = _mane_select(spec["mane"])
    out = []
    for gene, tx in zip(df[spec["gene"]], df[spec["column"]]):
        r = mane.get(gene)
        if r is None:
            raise SystemExit(f"{where_}: {gene} has no MANE Select row in {spec['mane']}")
        if tx == r["RefSeq_nuc"]:
            out.append(tx)
        elif tx == r["Ensembl_nuc"]:
            out.append(f"{r['RefSeq_nuc']} ({tx})")
        else:
            raise SystemExit(f"{where_}: {gene}'s transcript {tx} is neither its MANE Select "
                             f"RefSeq ({r['RefSeq_nuc']}) nor the Ensembl transcript of the same "
                             f"row ({r['Ensembl_nuc']})")
    df = df.copy()
    df[spec["column"]] = out
    return df


def finished(df: pd.DataFrame, notes: list[str], part: dict,
             where_: str) -> tuple[pd.DataFrame, list[str]]:
    """What a part prints, finished: its transcripts as RefSeq, and every number of 1,000
    or more in its headers, cells and notes grouped in thousands."""
    df = transcripts(df, part, where_)
    data = [[grouped(v) for v in r] for r in df.itertuples(index=False)]
    return (pd.DataFrame(data, columns=[grouped(c) for c in df.columns], dtype=str),
            [grouped(n) for n in notes])


def _csv_text(rows: list[list[str]]) -> str:
    out = io.StringIO()
    csv.writer(out, lineterminator="\n", quoting=csv.QUOTE_MINIMAL).writerows(rows)
    return out.getvalue()


def _finish_copy(path: Path, part: dict, where_: str) -> tuple[pd.DataFrame, list[str]]:
    """A printed table's copy, finished cell by cell (finished()), its comment lines with
    their numbers grouped too. The copy must read back as it was written, so that every
    cell finished() leaves alone keeps its bytes; the run stops otherwise. Returns the
    table and notes as printed."""
    text = path.read_bytes().decode("utf-8")
    lines = text.splitlines(keepends=True)
    head = []
    while lines and lines[0].startswith("#"):
        head.append(lines.pop(0))
    body = "".join(lines)
    rows = list(csv.reader(io.StringIO(body, newline="")))
    if _csv_text(rows) != body:
        raise SystemExit(f"{where_}: the printed table does not read back as it is written, so "
                         "its copy cannot be finished cell by cell")
    df = pd.DataFrame(rows[1:], columns=rows[0], dtype=str)
    df, _ = finished(df, [], part, where_)
    head = [grouped(h) for h in head]
    path.write_bytes(("".join(head) + _csv_text([list(df.columns)] + df.values.tolist())).encode("utf-8"))
    return df, [h.lstrip("#").strip() for h in head]


def _spec(item: dict, part: dict, column: str, cfg: dict) -> list:
    """A column's readable name and type: the part's own entry under `columns:` first,
    then the table's."""
    for key in (f"{item['content']}.{part['part']}", item["content"]):
        names = cfg["columns"].get(key, {})
        if column in names:
            return names[column]
    raise SystemExit(f"config/layout.yaml: {item['content']}.{part['part']}: no readable name "
                     f"for the column '{column}'")


def _side(row: dict) -> str:
    """The direction a threshold is used in: BP4's at or below it."""
    side = str(row.get("side", "")).lower()
    kind = str(row.get("threshold_kind", ""))
    return "lower" if side == "bp4" or kind.startswith("bp4") or kind == "walker_bp4" else "upper"


def _cell(text: str, kind: str, row: dict, extra: list, values: dict, where: str) -> str:
    """One cell as the tables v3.0.1 printed print a cell of its type; an empty cell
    stays empty."""
    if text == "" or kind == "text":
        return text
    if kind in NUMERIC:
        try:
            x = float(text)
        except ValueError:
            raise SystemExit(f"{where}: '{text}' is not a number") from None
        if np.isnan(x):
            return ""
        if np.isinf(x):
            return "∞" if x > 0 else "−∞"
        if kind == "count":
            if x != int(x):
                raise SystemExit(f"{where}: '{text}' is not a whole number")
            return S.fmt_int(x)
        if kind == "proportion":
            # three decimals; a source that stores its shares rounded already names its
            # decimals, so that they are not rounded twice
            return S.fmt_dp(x, int(extra[0]) if extra else 3)
        if kind == "lr":
            return S.fmt_lr(x)
        if kind == "threshold":
            # a published cut point is the configuration's constant, printed as written
            if row.get("threshold_source") == "published" or str(row.get("threshold_kind", "")).startswith("walker"):
                return S.fmt_const(x)
            tool = extra[0] if extra else row.get("tool")
            return S.fmt_thr(x, tool if tool in S.tool_labels() else None, _side(row))
        if kind == "score":
            return S.fmt_sig(x, 4)
        if kind == "percent":
            return S.fmt_dp(x, 1)
        if kind == "constant":
            return S.fmt_const(x)
        return S.fmt_sig(x, 3)                                      # number
    if kind == "boolean":
        if text.strip().lower() not in ("true", "false"):
            raise SystemExit(f"{where}: '{text}' is not True or False")
        return "Yes" if text.strip().lower() == "true" else "No"
    if kind == "tool":
        names = {**S.tool_labels(), **TOOL_EXTRA}
        if text not in names:
            raise SystemExit(f"{where}: no printed name for the tool '{text}'")
        return names[text]
    if kind == "stratum":
        # an external gene's band is cut at the depth its deposit reaches, as the printed
        # external-gene tables print it (DDX3X's 11-50 is 11-25)
        gene = row.get("gene") or row.get("gene_set") or row.get("Gene")
        if gene in S.S12_GENES and text in ("s3_10", "s11_50", "s3_50"):
            return _external_band(gene, text)
        names = {**S.STRATUM_LABEL, **STRATUM_EXTRA}
        if text not in names:
            raise SystemExit(f"{where}: no printed name for the stratum '{text}'")
        return names[text]
    if kind == "tier":
        return TIER_LABEL.get(text, text)
    if kind == "side":
        if text.lower() not in ("pp3", "bp4"):
            raise SystemExit(f"{where}: '{text}' is not PP3 or BP4")
        return text.upper()
    if kind == "label":
        if text in LABEL_BLANK:
            return ""
        try:
            return S.label_name(text, row.get("gene") or row.get("gene_set"))
        except KeyError:
            raise SystemExit(f"{where}: no printed name for the label '{text}'") from None
    if kind not in values:
        raise SystemExit(f"{where}: '{kind}' is no type and no map under display: values:")
    if text not in values[kind]:
        # a map's third entry, `numbers as written`, prints a whole number the map does
        # not name, or a bound such as >50, as it stands: an intronic offset beside the
        # bands' codes (E8.10)
        if extra and extra[0] == NUMBERS_AS_WRITTEN and re.fullmatch(r"[<>]?-?\d+", text):
            return text
        raise SystemExit(f"{where}: config/layout.yaml names no '{text}' under display: values: {kind}")
    return values[kind][text]


_BANDS: dict = {}


def _external_band(gene: str, stratum: str) -> str:
    if (gene, stratum) not in _BANDS:
        _BANDS[(gene, stratum)] = S._ext_band(gene, stratum)
    return _BANDS[(gene, stratum)]


def _holds(row: dict, where: dict) -> bool:
    """A row meets a condition: each named column holds the value, or one of the values
    a list gives; an empty string is an empty cell."""
    for c, want in where.items():
        if c not in row:
            raise SystemExit(f"config/layout.yaml: a condition names the column '{c}', which "
                             "the table does not have")
        if row[c].strip() not in (want if isinstance(want, list) else [want]):
            return False
    return True


def displayed(df: pd.DataFrame, item: dict, part: dict, cfg: dict) -> pd.DataFrame:
    """A source table as the supplement prints it: every column, under its readable name,
    each cell by its column's type. The part's options change only what is printed:
    `drop` leaves out a column that repeats another cell for cell (and stops the run if
    it does not); `omit` leaves out the columns it names, which the configuration says
    why the text does not use (and stops the run if the source lacks one); `dashes`
    print "–" in the columns named, in the rows a condition picks; `marks` add a column
    of text on the rows a named rule picks."""
    where_ = f"{item['content']}.{part['part']}"
    for gone, kept in part.get("drop", {}).items():
        if not (df[gone] == df[kept]).all():
            raise SystemExit(f"{where_}: '{gone}' is to be left out as a copy of '{kept}', "
                             f"and {int((df[gone] != df[kept]).sum())} rows differ")
    absent = [c for c in part.get("omit", []) if c not in df.columns]
    if absent:
        raise SystemExit(f"{where_}: the columns {absent} are to be left out, and the source "
                         "has no such column")
    df = df.drop(columns=list(part.get("drop", {})) + list(part.get("omit", [])))
    values = cfg["display"]["values"]
    specs = {c: _spec(item, part, c, cfg) for c in df.columns}
    names = [s[0] for s in specs.values()]
    twice = sorted({n for n in names if names.count(n) > 1})
    if twice:
        raise SystemExit(f"config/layout.yaml: {item['content']}.{part['part']}: names used "
                         f"twice: {twice}")
    for c, s in specs.items():
        if s[1] not in TYPES and s[1] not in values:
            raise SystemExit(f"config/layout.yaml: {item['content']}: the column '{c}' has "
                             f"no known type ({s[1]})")
    rows = df.to_dict("records")
    out = {}
    for c, (name, kind, *extra) in specs.items():
        where = f"{item['content']}.{part['part']}, column '{c}'"
        out[name] = [_cell(r[c], kind, r, extra, values, where) for r in rows]
    for d in part.get("dashes", []):
        hit = [_holds(r, d["where"]) for r in rows]
        if not any(hit):
            raise SystemExit(f"{where_}: no row meets the dashes' condition {d['where']}")
        for c in d["columns"]:
            out[specs[c][0]] = [DASH if h else v for h, v in zip(hit, out[specs[c][0]])]
    for m in part.get("marks", []):
        hit = [MARK_RULES[m["rule"]](r) for r in rows]
        if not any(hit):
            raise SystemExit(f"{where_}: the rule '{m['rule']}' marks no row")
        if m["name"] in names:
            raise SystemExit(f"{where_}: the mark's name '{m['name']}' is a column's already")
        names.insert(names.index(specs[m["before"]][0]), m["name"])
        out[m["name"]] = [m["text"] if h else "" for h in hit]
    return pd.DataFrame(out, columns=names)


# The pathogenic-side thresholds (thresholds_pathogenic) as evid_supp_tables.s6 prints
# them: a row per (tool, band, tier), its counts and summaries over the held-out folds,
# and each fold's ratio in one column
S6_KEY = ("Tool", "Stratum (intronic offset, bp)", "Tier")
S6_CUT = "Tier LR cut"
S6_FOLDS = "Held-out LR per fold"
S6_OVER_FOLDS = {
    "Held-out genes": lambda f, cut: str(len(f)),
    "Folds with a fitted threshold": lambda f, cut: str(int(f.threshold_from_6_genes.notna().sum())),
    "Folds with an evaluable held-out LR": lambda f, cut: str(int(f.heldout_lr.notna().sum())),
    "Held-out genes clearing the cut": lambda f, cut: str(int((f.heldout_lr >= cut).sum())),
    "Held-out LR, median": lambda f, cut: S.fmt_lr(f.heldout_lr.dropna().median()),
    "Held-out LR, minimum": lambda f, cut: S.fmt_lr(f.heldout_lr.dropna().min()),
    "Held-out LR, maximum": lambda f, cut: S.fmt_lr(f.heldout_lr.dropna().max()),
}


def not_counted(df: pd.DataFrame, item: dict, part: dict) -> tuple[pd.DataFrame, list[tuple]]:
    """The pathogenic-side thresholds with the folds a named rule picks listed but not
    counted, by the rule Table 2(a) counts by (evid_tables._folds_excluding): in a row with such a
    fold, the held-out genes, the folds with a fitted threshold, with an evaluable
    held-out ratio and clearing the cut, and the median, minimum and maximum held-out
    ratio are taken over the other folds, from the per-fold file at full precision and
    printed as S6 prints them; each picked fold keeps its entry in the per-fold column,
    followed by "†". Every other row and cell is printed as it stands.

    A row's printed cells are first recomputed over all its folds and must come out as
    printed, so that the two files tell one story; its cut is the tier's boundary, which
    prints in full. The run stops if they do not, if the rule picks no fold, or if it
    picks every fold of a row. Returns the table and the keys of the rows it changed."""
    spec = part["not_counted"]
    where_ = f"{item['content']}.{part['part']}"
    rule = MARK_RULES[spec["rule"]]
    folds = K.read_back(spec["folds"])
    folds = folds[folds.side == spec["side"]]
    picked_all = folds.apply(lambda r: bool(rule(r)), axis=1)
    if not picked_all.any():
        raise SystemExit(f"{where_}: the rule '{spec['rule']}' picks no fold")
    codes = {(S._tool(t), S.STRATUM_LABEL[st], tier.capitalize()): (t, st, tier)
             for t in S.TOOLS for st in S.IN_SCOPE for tier in S.TIERS}
    boundaries = {v for _, v in F.PATH_CUTS}
    out, changed = df.copy(), []
    for i, r in df.iterrows():
        key = tuple(r[c] for c in S6_KEY)
        tool, st, tier = codes[key]
        f = folds[(folds.tool == tool) & (folds.stratum == st) & (folds.tier == tier)]
        picked = picked_all.loc[f.index]
        if not picked.any():
            continue
        if picked.all():
            raise SystemExit(f"{where_}: {key}: the rule '{spec['rule']}' picks every fold")
        cut = float(r[S6_CUT])
        if cut not in boundaries or S.fmt_lr(cut) != r[S6_CUT]:
            raise SystemExit(f"{where_}: {key}: the cut {r[S6_CUT]} is no tier boundary as printed")
        f = f.sort_values("heldout_gene")
        entries = [S._fold_string(f.loc[[j]]) for j in f.index]
        recount = {c: fn(f, cut) for c, fn in S6_OVER_FOLDS.items()} | {S6_FOLDS: "; ".join(entries)}
        if recount != {c: r[c] for c in recount} or recount[S6_FOLDS] != S._fold_string(f):
            raise SystemExit(f"{where_}: {key}: the printed row and {spec['folds']} disagree")
        counted = f[~picked.loc[f.index]]
        for c, fn in S6_OVER_FOLDS.items():
            out.at[i, c] = fn(counted, cut)
        out.at[i, S6_FOLDS] = "; ".join(e + (NOT_COUNTED if picked.loc[j] else "")
                                        for e, j in zip(entries, f.index))
        changed.append(key)
    return out, changed


def _write_csv(df: pd.DataFrame, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, lineterminator="\n")
    return path


def _part_path(folder: Path, item: dict, part: dict) -> Path:
    if len(item["parts"]) == 1:
        return folder / f"{item['content']}.csv"
    return folder / f"{item['content']}.{part['part']}.csv"


@maker
def table_parts(item: dict, folder: Path, cfg: dict) -> list[Path]:
    """Each part from its source, or from several sources of one shape stacked: a list
    of sources as they stand, or a mapping from a value of a new first column ("key",
    its name) to a source, for files that do not name what they are about. A printed
    source is copied; any other is printed as displayed() prints it. A part's
    replacements are made in what is printed, a printed source's copy included, and
    every part is then finished (finished(): thousands grouped, transcripts as RefSeq)."""
    out = []
    item["_frames"] = []
    numbers = _numbers(cfg)
    for part in item["parts"]:
        reps = _replacements(part, numbers)
        where_ = f"{item['content']}.{part['part']}"
        if "sources" in part:
            srcs = part["sources"]
            pairs = srcs.items() if isinstance(srcs, dict) else [(None, x) for x in srcs]
            frames, notes = [], []
            for key, src in pairs:
                df, n = read_source(src)
                if part.get("key"):
                    df.insert(0, part["key"], key)
                frames.append(df)
                notes += [x for x in n if x not in notes]
            if len({tuple(f.columns) for f in frames}) != 1:
                raise SystemExit(f"{item['content']}.{part['part']}: the stacked sources differ in their columns")
            df = pd.concat(frames, ignore_index=True)
            printed = False
        else:
            df, notes = read_source(part["source"])
            printed = str(part["source"]).startswith(PRINTED)
        p = _part_path(folder, item, part)
        if printed and "not_counted" in part:
            df = not_counted(df, item, part)[0]
            df, notes = finished(*replaced(df, notes, reps, where_), part, where_)
            _write_csv(df, p)
        elif printed:
            if reps:
                _copy_replaced(part["source"], p, reps, where_)
            else:
                p.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(part["source"], p)
            df, notes = _finish_copy(p, part, where_)
        else:
            df = displayed(df, item, part, cfg)
            df, notes = finished(*replaced(df, notes, reps, where_), part, where_)
            _write_csv(df, p)
        out.append(p)
        item["_frames"].append((part, df, notes + list(part.get("notes", []))))
    return out


@maker
def variant_list(item: dict, folder: Path, cfg: dict) -> list[Path]:
    """Every analysed variant with its HGVS description, and their counts
    (src/evid_variant_list.py, E8.5): written after the list's six checks, and only if
    the item's sources are the ones the list reads."""
    from . import evid_variant_list as V
    if sorted(item["sources"]) != sorted(V.SOURCES):
        raise SystemExit(f"[layout] {item['content']}: the sources in config/layout.yaml are "
                         f"not the ones the list reads ({sorted(V.SOURCES)})")
    return V.write(folder)


# the distance bands as the text names them (the conventions)
BAND_NAME = {"s3_10": "proximal band", "s11_50": "distal band", "s3_50": "pooled window"}


@maker
def table_practical_performance(item: dict, folder: Path, cfg: dict) -> list[Path]:
    """The practical-performance table: (a) at the study's thresholds in the seven genes,
    held out, as docs/practical-metrics-summary.md prints its held-out values (its
    in-sample values, in brackets there, are in the supplementary table of practical
    metrics); (b) against ClinVar, at least one star, each predictor at the published
    cut points and its in-sample thresholds, in the bands the options name."""
    tools = item["options"]["tools"]
    table = K.read_back(REPORT_DIR / "practical_metrics_by_threshold.csv")
    tiers = ("supporting", "moderate", "strong")
    metrics = {"sensitivity": "sensitivity", "specificity": "specificity",
               "coverage_labelled": "coverage"}
    rows = []
    for tool in tools:
        for st in P.CURVE_STRATA:
            r = {"Predictor": F.NAME[tool], "Band": F.STRATUM_LABEL[st]}
            folds, held_out = [], 0
            for t in tiers:
                ins = P._row(table, tool, st, t, "in_sample")
                lg = P._row(table, tool, st, t, "logo")
                held = lg if P._shown(ins, lg) else None
                for m, name in metrics.items():
                    r[f"{P.TIER_NAME[t]}: {name}"] = P._two(held[m]) if held is not None else "–"
                held_out = max(held_out, P._held_out_folds(lg))
                folds.append(str(int(lg.folds_used)) if lg is not None and np.isfinite(lg.folds_used) else "–")
            # not Table 2(a)'s k/n, which counts the folds with an evaluable ratio
            r["Folds with a fitted threshold (S, M, St)"] = (", ".join(folds)
                                                             + (f" of {held_out}" if held_out else ""))
            rows.append(r)
    sge = pd.DataFrame(rows)
    c1 = K.read_back(REPORT_DIR / "suppl_analyses" / "c1_clinvar_reference_metrics.csv")
    c1 = c1[(c1.star_level == "stars_ge_1") & (c1.analysis == "main") & (c1.measure == "threshold")]
    cols = [("pp3", "walker_pp3", "P/LP at or above: published 0.2"),
            ("pp3", "supporting", "P/LP at or above: Supporting"),
            ("pp3", "moderate", "P/LP at or above: Moderate"),
            ("pp3", "strong", "P/LP at or above: Strong"),
            ("bp4", "walker_bp4", "B/LB at or below: published 0.1"),
            ("bp4", "bp4_supporting", "B/LB at or below: Supporting"),
            ("bp4", "bp4_moderate", "B/LB at or below: Moderate")]
    bands = item["options"]["clinvar_bands"]
    rows = []
    for tool in tools:
        for st in bands:
            sub = c1[(c1.tool == tool) & (c1.stratum == st)]
            if not len(sub):
                continue
            r = {"Predictor": F.NAME[tool], "Band": F.STRATUM_LABEL[st],
                 "P/LP": f"{int(sub.n_plp.iloc[0]):,}", "B/LB": f"{int(sub.n_blb.iloc[0]):,}"}
            for side, kind, name in cols:
                x = sub[(sub.side == side) & (sub.threshold_kind == kind)]
                if not len(x):
                    r[name] = "–"          # no such threshold: see the table's note
                    continue
                x = x.iloc[0]
                if x.status != "ok":
                    r[name] = "not evaluable"
                elif side == "pp3":
                    r[name] = f"{P._two(x.sensitivity)} ({int(x.n_pos_above):,}/{int(x.n_plp):,})"
                else:
                    r[name] = f"{P._two(x.bp4_sensitivity)} ({int(x.bp4_n_neg_at_or_below):,}/{int(x.n_blb):,})"
            rows.append(r)
    clinvar = pd.DataFrame(rows)
    # the bands left out of (b): the approved notes name each, with its P/LP variants of at
    # least one star (E8.11), and the run stops if a count there is not this file's
    said = " ".join(_resolve(n, _numbers(cfg)) for n in item.get("notes", []))
    for st in [s for s in P.CURVE_STRATA if s not in bands]:
        sub = c1[c1.tool.isin(tools) & (c1.stratum == st)]
        n_plp = sorted({int(v) for v in sub.n_plp})
        phrase = (f"The {BAND_NAME[st]}, with {n_plp[0]} P/LP variant{'' if n_plp[0] == 1 else 's'}, "
                  "is left out." if len(n_plp) == 1 else None)
        if phrase is None or phrase not in said:
            raise SystemExit(f"[layout] {item['content']}: (b) leaves out {F.STRATUM_LABEL[st]}, "
                             f"{n_plp} P/LP variants, which its notes do not say as '{phrase}'")
    item["_frames"] = []
    out = []
    for part, df in zip(item["parts"], (sge, clinvar)):
        out.append(_write_csv(df, _part_path(folder, item, part)))
        item["_frames"].append((part, df, []))
    return out


# ---------------------------------------------------------------------------
# One page or more: the default destination of a supplementary table
# ---------------------------------------------------------------------------
LINE_MM = 4.0            # a line of 7 pt type with its spacing
CHAR_MM = 1.3            # an average character of 7 pt type
CELL_CHARS = 30          # a longer cell wraps
PAGE_MM = 247            # an A4 page's text block, the long side; the short side is the width


def one_page(frames: list, width_mm: float) -> tuple[bool, int, int]:
    """An estimate of whether a table fits on one A4 page, upright or turned: its lines
    (each part's title, notes, header and rows) and its widest part (each column as wide
    as its widest cell, wrapped at CELL_CHARS, or the longest word of its header), within
    the text block one way or the other. Returns the verdict, the lines and the width."""
    lines = sum(2 + len(notes) + len(df) for _, df, notes in frames)
    wide = max(sum(min(max([len(str(v)) for v in df[c]] + [max(len(w) for w in str(c).split())]),
                       CELL_CHARS) * CHAR_MM + 2 for c in df.columns) for _, df, _ in frames)
    fits = any(lines * LINE_MM <= tall and wide <= across
               for tall, across in ((PAGE_MM, width_mm), (width_mm, PAGE_MM)))
    return fits, lines, round(wide)


# ---------------------------------------------------------------------------
# The workbook: a minimal Office Open XML writer, so that the file is the same bytes
# on every run (fixed member order and times, no document properties). A cell that
# prints a plain decimal is a number shown with the decimals it prints, and one that
# prints it grouped in thousands (1,410) is the same number shown so; any other is text,
# with a gene's italics, marked *BRCA1* in the text (the CSV keeps the marks), as an
# italic run.
# ---------------------------------------------------------------------------
DECIMAL = re.compile(r"^-?\d+(?:\.(\d+))?$")
GROUPED_DECIMAL = re.compile(r"^-?\d{1,3}(?:,\d{3})+(?:\.(\d+))?$")
ITALIC = re.compile(r"(?<![\w*])\*([A-Z][A-Za-z0-9-]*)\*(?![\w*])")
ZIP_TIME = (1980, 1, 1, 0, 0, 0)
MAX_DECIMALS = 15
NUMBER_STYLE = 2         # cellXfs index of 0 decimals; d decimals at NUMBER_STYLE + d
GROUPED_STYLE = NUMBER_STYLE + MAX_DECIMALS + 1      # the same, grouped in thousands


def _col(i: int) -> str:
    s = ""
    i += 1
    while i:
        i, r = divmod(i - 1, 26)
        s = chr(65 + r) + s
    return s


def _runs_xml(text: str, bold: bool) -> str:
    """An inline string's runs: the text, with each *GENE* as an italic run without its
    marks; a single plain run where there is none."""
    if not ITALIC.search(text):
        return f'<t xml:space="preserve">{escape(text)}</t>'
    font = ("<b/>" if bold else "") + '<sz val="10"/><rFont val="Arial"/>'
    runs, at = [], 0
    for m in ITALIC.finditer(text):
        runs += [(text[at:m.start()], False), (m.group(1), True)]
        at = m.end()
    runs.append((text[at:], False))
    return "".join(f'<r><rPr>{"<i/>" if it else ""}{font}</rPr><t xml:space="preserve">{escape(s)}</t></r>'
                   for s, it in runs if s)


def _cell_xml(ref: str, value, bold: bool) -> str:
    text = "" if value is None else str(value)
    if text == "":
        return ""
    for pattern, style, digits in ((DECIMAL, NUMBER_STYLE, text),
                                   (GROUPED_DECIMAL, GROUPED_STYLE, text.replace(",", ""))):
        m = pattern.match(text)
        if not bold and m and len(digits.lstrip("-").replace(".", "").lstrip("0")) <= 15:
            places = len(m.group(1) or "")
            if places <= MAX_DECIMALS:
                return f'<c r="{ref}" s="{style + places}"><v>{digits}</v></c>'
    style = ' s="1"' if bold else ""
    return f'<c r="{ref}"{style} t="inlineStr"><is>{_runs_xml(text, bold)}</is></c>'


def _sheet_xml(rows: list[tuple[list, bool]]) -> str:
    body = []
    for r, (cells, bold) in enumerate(rows, start=1):
        xs = "".join(_cell_xml(f"{_col(c)}{r}", v, bold) for c, v in enumerate(cells))
        body.append(f'<row r="{r}">{xs}</row>')
    return ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
            '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
            f'<sheetData>{"".join(body)}</sheetData></worksheet>')


def _styles_xml(ns: str) -> str:
    # d decimals at numFmtId 164 + d; the same grouped in thousands after them
    codes = [("0" if d == 0 else "0." + "0" * d) for d in range(MAX_DECIMALS + 1)]
    codes += ["#,##" + c for c in codes]
    formats = "".join(f'<numFmt numFmtId="{164 + k}" formatCode="{c}"/>' for k, c in enumerate(codes))
    numbers = "".join(f'<xf numFmtId="{164 + k}" fontId="0" fillId="0" borderId="0" xfId="0" '
                      'applyNumberFormat="1"/>' for k in range(len(codes)))
    return ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
            f'<styleSheet xmlns="{ns}"><numFmts count="{len(codes)}">{formats}</numFmts>'
            '<fonts count="2"><font><sz val="10"/><name val="Arial"/></font>'
            '<font><b/><sz val="10"/><name val="Arial"/></font></fonts>'
            '<fills count="2"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill></fills>'
            '<borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>'
            '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
            f'<cellXfs count="{GROUPED_STYLE + MAX_DECIMALS + 1}"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>'
            '<xf numFmtId="0" fontId="1" fillId="0" borderId="0" xfId="0" applyFont="1"/>'
            f'{numbers}</cellXfs>'
            '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>'
            "</styleSheet>")


def write_xlsx(path: Path, sheets: list[tuple[str, list[tuple[list, bool]]]]) -> Path:
    ns = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
    rel = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
    files = {
        "[Content_Types].xml": (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
            '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
            '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
            '<Default Extension="xml" ContentType="application/xml"/>'
            '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
            '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>'
            + "".join(f'<Override PartName="/xl/worksheets/sheet{i}.xml" ContentType="application/'
                      'vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
                      for i in range(1, len(sheets) + 1)) + "</Types>"),
        "_rels/.rels": (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            f'<Relationship Id="rId1" Type="{rel}/officeDocument" Target="xl/workbook.xml"/>'
            "</Relationships>"),
        "xl/workbook.xml": (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
            f'<workbook xmlns="{ns}" xmlns:r="{rel}"><sheets>'
            + "".join(f'<sheet name="{escape(name)}" sheetId="{i}" r:id="rId{i}"/>'
                      for i, (name, _) in enumerate(sheets, start=1)) + "</sheets></workbook>"),
        "xl/_rels/workbook.xml.rels": (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            + "".join(f'<Relationship Id="rId{i}" Type="{rel}/worksheet" Target="worksheets/sheet{i}.xml"/>'
                      for i in range(1, len(sheets) + 1))
            + f'<Relationship Id="rId{len(sheets) + 1}" Type="{rel}/styles" Target="styles.xml"/>'
            "</Relationships>"),
        "xl/styles.xml": _styles_xml(ns),
    }
    for i, (_, rows) in enumerate(sheets, start=1):
        files[f"xl/worksheets/sheet{i}.xml"] = _sheet_xml(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(path, "w") as z:
        for name, text in files.items():
            info = zipfile.ZipInfo(name, date_time=ZIP_TIME)
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            z.writestr(info, text.encode("utf-8"), compresslevel=6)
    return path


def _sheet_name(number: str) -> str:
    return number.replace("Supplementary Table ", "")


def part_title(item: dict, k: int, numbers: dict) -> str:
    """A part's title as printed above it, its letter first: "(a) By gene"."""
    part = item["parts"][k]
    return f"({chr(97 + k)}) {grouped(_resolve(str(part.get('title', part['part'])), numbers))}"


def _sheet_rows(it: dict, numbers: dict, comments: bool = True) -> list[tuple[list, bool]]:
    """A table's sheet: its number and title, its notes, then each part after an empty
    row, its title (for a table of several parts), the comment lines its source carries
    above its header (unless `comments` is false) and the part's own notes, its header
    and its rows."""
    rows: list[tuple[list, bool]] = [([f"{it['number']}. {grouped(_resolve(it['title'], numbers))}"], True)]
    for n in _notes(it, numbers):
        rows.append(([n], False))
    for k, (part, df, notes) in enumerate(it["_frames"]):
        rows.append(([], False))
        if len(it["parts"]) > 1:
            rows.append(([part_title(it, k, numbers)], True))
        for n in notes if comments else part.get("notes", []):
            rows.append(([grouped(_resolve(n, numbers))], False))
        rows.append((list(df.columns), True))
        rows += [(list(r), False) for r in df.itertuples(index=False)]
    return rows


def workbook(items: list[dict], cfg: dict, numbers: dict, path: Path) -> Path:
    fmt = cfg["supplementary_format"]
    readme: list[tuple[list, bool]] = [
        ([f"Supplementary tables ({fmt['numbering']}; titles and descriptions: config/layout.yaml)"], True),
        ([_resolve(fmt["rule"], numbers)], False),
        ([], False),
        (["Table", "Sheet", "Content", "Title", "Destination",
          "Description", "Parts", "Source files", "Source sha256"], True)]
    sheets = []
    for it in items:
        srcs = _sources(it)
        title = grouped(_resolve(it["title"], numbers))
        readme.append(([it["number"], _sheet_name(it["number"]), it["content"], title,
                        it["destination"], grouped(_resolve(it["description"], numbers)),
                        "; ".join(p["part"] for p in it["parts"]) if len(it["parts"]) > 1 else "",
                        "; ".join(_rel(s) for s in srcs), "; ".join(_sha256(s) for s in srcs)], False))
        sheets.append((_sheet_name(it["number"]), _sheet_rows(it, numbers)))
    readme.append(([], False))
    readme.append((["Kept in the repository only (not numbered)"], True))
    readme.append((["Path", "What it holds"], True))
    readme += [([_rel(r["path"]), _resolve(r["what"], numbers)], False)
               for r in cfg.get("repository_only", [])]
    return write_xlsx(path, [("README", readme)] + sheets)


def submission_workbook(items: list[dict], cfg: dict, numbers: dict, path: Path) -> Path:
    """The workbook for submission (E8.10): a contents sheet, then a sheet for each table
    marked excel, as the repository's workbook prints it (its cells, numbers and italics
    the same) without the comment lines a source carries above its header. The contents
    sheet opens with its heading and the conventions, lists the tables the file holds by
    number and title, and ends with the tables marked pdf, listed from their numbers."""
    fmt = cfg["supplementary_format"]
    sub = fmt["submission_workbook"]
    excel = [it for it in items if it["destination"] == "excel"]
    pdf = [it for it in items if it["destination"] == "pdf"]
    contents: list[tuple[list, bool]] = [([sub["heading"]], True), ([], False),
                                         ([grouped(_resolve(fmt["conventions"], numbers))], False),
                                         ([], False)]
    contents += [([it["number"], grouped(_resolve(it["title"], numbers))], False) for it in excel]
    if pdf:
        named = _resolve("{" + ", ".join(it["content"] for it in pdf) + "}", numbers)
        contents += [([], False), ([f"{named} {sub['pdf_tables']}"], False)]
    sheets = [(sub["sheet"], contents)]
    sheets += [(_sheet_name(it["number"]), _sheet_rows(it, numbers, comments=False)) for it in excel]
    return write_xlsx(path, sheets)


# ---------------------------------------------------------------------------
# The main tables as editable files (E8.11): a minimal WordprocessingML writer, so that
# a file is the same bytes on every run, as the workbook is (fixed member order and
# times, no document properties). A table's cells are its CSV's, read back from the file
# the stage wrote; a gene's italics, marked *BRCA1* there, are italic runs.
# ---------------------------------------------------------------------------
W_NS = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
A4_MM = (210.0, 297.0)
TWIPS_PER_MM = 1440 / 25.4
CELL_MARGIN_MM = 1.9            # a cell's left and right margins (Word's default, 0.19 cm)
CELL_SLACK_MM = 0.4             # and a little more, so that a word measured to fit does
MM_PER_PT = 25.4 / 72
NO_GLYPH = 600                  # a glyph Helvetica's metrics lack (≥, ≤, ρ): its width, in 1/1000 em


def _twips(mm: float) -> int:
    return round(mm * TWIPS_PER_MM)


_AFM: dict = {}


def text_mm(text: str, pt: float, bold: bool = False) -> float:
    """A text's width in Arial at `pt`, in mm, from the metrics of Helvetica, which
    Arial's match glyph for glyph (matplotlib's own copy, the same wherever the stage
    runs); italics as wide as upright."""
    from fontTools.agl import UV2AGL
    from matplotlib import _afm, get_data_path
    if bold not in _AFM:
        with open(Path(get_data_path()) / "fonts" / "afm" / ("phvb8a.afm" if bold else "phvr8a.afm"), "rb") as f:
            _AFM[bold] = _afm.AFM(f)
    units = 0.0
    for ch in text:
        try:
            units += _AFM[bold].get_width_from_char_name(UV2AGL[ord(ch)])
        except KeyError:
            units += NO_GLYPH
    return units / 1000 * pt * MM_PER_PT


def _w_runs(text: str, pt: float, font: str, bold: bool = False) -> str:
    """A paragraph's runs: the text, each *GENE* an italic run without its marks."""
    def props(italic: bool) -> str:
        return (f'<w:rPr><w:rFonts w:ascii="{font}" w:hAnsi="{font}" w:cs="{font}"/>'
                + ("<w:b/>" if bold else "") + ("<w:i/>" if italic else "")
                + f'<w:sz w:val="{round(pt * 2)}"/><w:szCs w:val="{round(pt * 2)}"/></w:rPr>')
    runs, at = [], 0
    for m in ITALIC.finditer(text):
        runs += [(text[at:m.start()], False), (m.group(1), True)]
        at = m.end()
    runs.append((text[at:], False))
    return "".join(f'<w:r>{props(it)}<w:t xml:space="preserve">{escape(s)}</w:t></w:r>'
                   for s, it in runs if s)


def _w_p(runs: str, before_pt: float = 0, after_pt: float = 0, keep_next: bool = False) -> str:
    ppr = (("<w:keepNext/>" if keep_next else "")
           + f'<w:spacing w:before="{round(before_pt * 20)}" w:after="{round(after_pt * 20)}"/>')
    return f"<w:p><w:pPr>{ppr}</w:pPr>{runs}</w:p>"


def _plain(text: str) -> str:
    return ITALIC.sub(lambda m: m.group(1), text)


def column_widths(rows: list[list[str]], pt: float) -> tuple[list[float], list[float]]:
    """Each column's width in mm, its margins included: as wide as its widest cell, its
    bold header wrapping at its spaces (natural); and as wide as its longest word (least)."""
    pad = 2 * CELL_MARGIN_MM + CELL_SLACK_MM
    head, body = rows[0], rows[1:]
    natural, least = [], []
    for j in range(len(head)):
        cells = [_plain(r[j]) for r in body]
        head_word = max([text_mm(w, pt, bold=True) for w in _plain(head[j]).split()] or [0])
        body_word = max([text_mm(w, pt) for c in cells for w in c.split()] or [0])
        natural.append(max([text_mm(c, pt) for c in cells] + [head_word]) + pad)
        least.append(max(head_word, body_word) + pad)
    return natural, least


def page_for(tables: list[list[list[str]]], cfg: dict) -> tuple[bool, list[list[float]]]:
    """A main table's page and its columns: upright if every part fits the upright
    text width at its natural widths, else turned; each part's columns then share the
    width in proportion to their natural widths, none below its least. The run stops if
    a part is wider than the turned page at its least widths."""
    d = cfg["docx"]
    upright, turned = (A4_MM[0] - 2 * d["margin_cm"] * 10), (A4_MM[1] - 2 * d["margin_cm"] * 10)
    sizes = [column_widths(t, d["table_pt"]) for t in tables]
    landscape = any(sum(n) > upright for n, _ in sizes)
    width = turned if landscape else upright
    out = []
    for natural, least in sizes:
        if sum(least) > width + 1e-9:
            raise SystemExit(f"[layout] a main table needs {sum(least):.0f} mm, more than the "
                             f"{width:.0f} mm of the turned page")
        w = [n * width / sum(natural) for n in natural]
        fixed = [False] * len(w)
        while True:                                   # a column below its least is set to it
            low = [not f and x < m - 1e-9 for x, m, f in zip(w, least, fixed)]
            free = sum(n for n, f, lo in zip(natural, fixed, low) if not (f or lo))
            if not any(low) or not free:
                break
            fixed = [f or lo for f, lo in zip(fixed, low)]
            rest = width - sum(m for m, f in zip(least, fixed) if f)
            w = [m if f else n * rest / free for n, m, f in zip(natural, least, fixed)]
        out.append(w)
    return landscape, out


def _w_table(rows: list[list[str]], widths_mm: list[float], cfg: dict) -> str:
    d = cfg["docx"]
    rule = f'w:val="single" w:sz="{round(d["rule_pt"] * 8)}" w:space="0" w:color="000000"'
    tw = [_twips(w) for w in widths_mm]
    out = [f'<w:tbl><w:tblPr><w:tblW w:w="{sum(tw)}" w:type="dxa"/>'
           f'<w:tblBorders><w:top {rule}/><w:left w:val="nil"/><w:bottom {rule}/><w:right w:val="nil"/>'
           '<w:insideH w:val="nil"/><w:insideV w:val="nil"/></w:tblBorders>'
           '<w:tblLayout w:type="fixed"/>'
           f'<w:tblCellMar><w:top w:w="0" w:type="dxa"/><w:left w:w="{_twips(CELL_MARGIN_MM)}" w:type="dxa"/>'
           f'<w:bottom w:w="0" w:type="dxa"/><w:right w:w="{_twips(CELL_MARGIN_MM)}" w:type="dxa"/></w:tblCellMar>'
           '<w:tblLook w:val="0000" w:firstRow="0" w:lastRow="0" w:firstColumn="0" w:lastColumn="0" '
           'w:noHBand="1" w:noVBand="1"/></w:tblPr>',
           "<w:tblGrid>" + "".join(f'<w:gridCol w:w="{w}"/>' for w in tw) + "</w:tblGrid>"]
    for i, row in enumerate(rows):
        head = i == 0
        cells = []
        for text, w in zip(row, tw):
            borders = f"<w:tcBorders><w:bottom {rule}/></w:tcBorders>" if head else ""
            cells.append(f'<w:tc><w:tcPr><w:tcW w:w="{w}" w:type="dxa"/>{borders}<w:vAlign w:val="bottom"/></w:tcPr>'
                         + _w_p(_w_runs(text, d["table_pt"], d["font"], bold=head), before_pt=1, after_pt=1)
                         + "</w:tc>")
        trpr = "<w:trPr><w:cantSplit/><w:tblHeader/></w:trPr>" if head else "<w:trPr><w:cantSplit/></w:trPr>"
        out.append(f"<w:tr>{trpr}{''.join(cells)}</w:tr>")
    return "".join(out) + "</w:tbl>"


def written_rows(path: Path) -> list[list[str]]:
    """A table as the stage wrote it: its header and rows, the comment lines above them
    left out."""
    lines = [l for l in path.read_text().splitlines(keepends=True) if not l.startswith("#")]
    return list(csv.reader(io.StringIO("".join(lines))))


def main_table_docx(item: dict, cfg: dict, numbers: dict, folder: Path) -> Path:
    """A main table as an editable file: A4, "Table N." in bold and its title, each
    part's title above its table, the table ruled across at its top, below its header and
    at its foot, and the notes in one paragraph below; the cells read back from the CSVs
    the stage wrote for the table."""
    d = cfg["docx"]
    tables = [written_rows(_part_path(folder, item, p)) for p in item["parts"]]
    landscape, widths = page_for(tables, cfg)
    body = [_w_p(_w_runs(f"{item['number']}.", d["text_pt"], d["font"], bold=True)
                 + _w_runs(f" {grouped(_resolve(item['title'], numbers))}", d["text_pt"], d["font"]),
                 after_pt=6, keep_next=True)]
    for k, (rows, w) in enumerate(zip(tables, widths)):
        if len(item["parts"]) > 1:
            body.append(_w_p(_w_runs(part_title(item, k, numbers), d["text_pt"], d["font"]),
                             before_pt=6 if k else 0, after_pt=3, keep_next=True))
        body.append(_w_table(rows, w, cfg))
    body.append(_w_p(_w_runs(" ".join(_notes(item, numbers)), d["text_pt"], d["font"]), before_pt=6))
    page_w, page_h = (A4_MM[1], A4_MM[0]) if landscape else A4_MM
    m = _twips(d["margin_cm"] * 10)
    body.append(f'<w:sectPr><w:pgSz w:w="{_twips(page_w)}" w:h="{_twips(page_h)}"'
                + (' w:orient="landscape"' if landscape else "")
                + f'/><w:pgMar w:top="{m}" w:right="{m}" w:bottom="{m}" w:left="{m}" w:header="709" '
                  'w:footer="709" w:gutter="0"/></w:sectPr>')
    rel = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
    files = {
        "[Content_Types].xml": (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
            '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
            '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
            '<Default Extension="xml" ContentType="application/xml"/>'
            '<Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>'
            '<Override PartName="/word/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.styles+xml"/>'
            "</Types>"),
        "_rels/.rels": (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            f'<Relationship Id="rId1" Type="{rel}/officeDocument" Target="word/document.xml"/></Relationships>'),
        "word/_rels/document.xml.rels": (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            f'<Relationship Id="rId1" Type="{rel}/styles" Target="styles.xml"/></Relationships>'),
        "word/styles.xml": (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
            f'<w:styles xmlns:w="{W_NS}"><w:docDefaults><w:rPrDefault><w:rPr>'
            f'<w:rFonts w:ascii="{d["font"]}" w:hAnsi="{d["font"]}" w:cs="{d["font"]}" w:eastAsia="{d["font"]}"/>'
            f'<w:sz w:val="{round(d["text_pt"] * 2)}"/><w:szCs w:val="{round(d["text_pt"] * 2)}"/>'
            '<w:lang w:val="en-GB"/></w:rPr></w:rPrDefault><w:pPrDefault><w:pPr>'
            '<w:spacing w:after="0" w:line="240" w:lineRule="auto"/></w:pPr></w:pPrDefault></w:docDefaults>'
            '<w:style w:type="paragraph" w:default="1" w:styleId="Normal"><w:name w:val="Normal"/></w:style>'
            "</w:styles>"),
        "word/document.xml": (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
            f'<w:document xmlns:w="{W_NS}"><w:body>{"".join(body)}</w:body></w:document>'),
    }
    path = folder / f"{item['content']}.docx"
    path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(path, "w") as z:
        for name, text in files.items():
            info = zipfile.ZipInfo(name, date_time=ZIP_TIME)
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            z.writestr(info, text.encode("utf-8"), compresslevel=6)
    return path


# ---------------------------------------------------------------------------
# The run
# ---------------------------------------------------------------------------
REFERENCE = re.compile(r"\{([a-z0-9_]+(?:,\s*[a-z0-9_]+)*)\}")
# a number as its label, its letter and its rank: "Supplementary Table", "S", 12
NUMBER_PARTS = re.compile(r"^(.*\S)\s+([A-Z]?)(\d+)$")


SECTIONS = ("main", "supplementary", "variants")   # each a folder of layout/


def _numbers(cfg: dict) -> dict:
    """Each item's number by its content name, the supplementary notes' among them."""
    items = [it for s in SECTIONS for it in cfg.get(s, [])] + cfg.get("supplementary_notes", [])
    return {it["content"]: it["number"] for it in items}


def listed(numbers_: list[str]) -> str:
    """Several items' numbers as the text lists them (E8.10): their label once, in the
    plural, then their numbers in order, three or more in a row as a range and the last
    joined by "and": Supplementary Tables S8 and S9, Supplementary Tables S10–S12,
    Supplementary Tables S5, S6, S12, S25 and S28. Items of different kinds, or one item
    twice, stop the run."""
    parts = [NUMBER_PARTS.match(n) for n in numbers_]
    if not all(parts) or len({(m.group(1), m.group(2)) for m in parts}) != 1:
        raise SystemExit(f"config/layout.yaml: {numbers_} cannot be listed under one label")
    ranks = sorted(int(m.group(3)) for m in parts)
    if len(set(ranks)) != len(ranks):
        raise SystemExit(f"config/layout.yaml: {numbers_} lists an item twice")
    label, letter = parts[0].group(1), parts[0].group(2)
    runs: list[list[int]] = []
    for k in ranks:
        if runs and k == runs[-1][-1] + 1:
            runs[-1].append(k)
        else:
            runs.append([k])
    words = []
    for r in runs:
        words += [f"{letter}{r[0]}–{letter}{r[-1]}"] if len(r) >= 3 else [f"{letter}{k}" for k in r]
    return f"{label}s " + (words[0] if len(words) == 1 else ", ".join(words[:-1]) + " and " + words[-1])


def _resolve(text: str, numbers: dict) -> str:
    """A text's {content} references, as the items' numbers; {a, b, c} lists several
    (listed())."""
    def number(m):
        names = [n.strip() for n in m.group(1).split(",")]
        unknown = [n for n in names if n not in numbers]
        if unknown:
            raise SystemExit(f"config/layout.yaml: a text names {unknown}, which is no item")
        return numbers[names[0]] if len(names) == 1 else listed([numbers[n] for n in names])
    return REFERENCE.sub(number, text)


def _notes(item: dict, numbers: dict) -> list[str]:
    return [grouped(_resolve(n, numbers)) for n in list(item.get("notes", [])) + item.get("_notes", [])]


def _part_notes(item: dict, numbers: dict) -> list[str]:
    """The configuration's notes on a part, named by the part (the workbook writes them
    above it)."""
    return [f"{p['part']}: {grouped(_resolve(n, numbers))}" for p in item.get("parts", [])
            for n in p.get("notes", [])]


def _sources(item: dict) -> list[str]:
    out = list(item.get("sources", []))
    for part in item.get("parts", []):
        if "source" in part:
            out.append(part["source"])
        srcs = part.get("sources", [])
        out += list(srcs.values()) if isinstance(srcs, dict) else list(srcs)
        if "not_counted" in part:                    # the per-fold file it counts from
            out.append(part["not_counted"]["folds"])
        if "transcripts" in part:                    # the MANE summary it reads them from
            out.append(part["transcripts"]["mane"])
    return list(dict.fromkeys(out))


def load_config(path: Path = CONFIG) -> dict:
    cfg = yaml.safe_load(path.read_text())
    items = [it for s in SECTIONS for it in cfg.get(s, [])]
    notes_ = cfg.get("supplementary_notes", [])
    for key in ("content", "number"):
        seen = [it[key] for it in items + notes_]
        twice = sorted({c for c in seen if seen.count(c) > 1})
        if twice:
            raise SystemExit(f"{path}: {key}s used twice: {twice}")
    # the supplementary notes, numbered from 1 in their list's order as the tables are
    if [n["number"] for n in notes_] != [f"Supplementary Note {k}" for k in range(1, len(notes_) + 1)]:
        raise SystemExit(f"{path}: the supplementary notes are listed as "
                         f"{[n['number'] for n in notes_]}, not numbered in their list's order")
    missing = [s for it in items for s in _sources(it) if not Path(s).is_file()]
    missing += [r["path"] for r in cfg.get("repository_only", []) if not Path(r["path"]).is_file()]
    if missing:
        raise SystemExit(f"{path}: sources not found: {missing}")
    # each list in its number's order, numbered from 1 without a gap: the place of an item
    # in contents.csv and of a table's sheet in the workbook follow its number
    for section, prefix in (("main", ""), ("supplementary", "Supplementary ")):
        for kind in ("figure", "table"):
            have = [it["number"] for it in cfg[section] if it["kind"] == kind]
            s = "S" if section == "supplementary" else ""
            want = [f"{prefix}{kind.capitalize()} {s}{k}" for k in range(1, len(have) + 1)]
            if have != want:
                raise SystemExit(f"{path}: the {section} {kind}s are listed as {have}, not "
                                 f"numbered in their list's order")
    numbers = _numbers(cfg)
    for it in items:
        for part in it.get("parts", []):
            for r in part.get("replace", []):
                if (set(r) != {"old", "new", "count"} or not isinstance(r["old"], str) or not r["old"]
                        or not isinstance(r["new"], str) or not isinstance(r["count"], int)
                        or r["count"] < 1):
                    raise SystemExit(f"{path}: {it['content']}.{part['part']}: a replacement is "
                                     "{old: text, new: text, count: a whole number from 1}")
                _resolve(r["new"], numbers)
            for text in [part.get("title") or ""] + list(part.get("notes", [])):
                _resolve(str(text), numbers)
        for key in ("title", "description", "legend"):
            if key in it:
                _resolve(it[key], numbers)
        for n in it.get("notes", []):
            _resolve(n, numbers)
    fmt = cfg["supplementary_format"]
    if not fmt.get("numbering"):
        raise SystemExit(f"{path}: supplementary_format names no numbering")
    for key in ("conventions", "rule"):
        _resolve(fmt[key], numbers)
    limit = cfg["supplementary_format"]["description_words"]
    for it in items:
        if it["make"] not in MAKERS:
            raise SystemExit(f"{path}: {it['content']} names no known maker ({it['make']})")
    for it in cfg["supplementary"]:
        if it.get("destination") not in ("pdf", "excel"):
            raise SystemExit(f"{path}: {it['content']}: destination is pdf or excel")
        if not it.get("description") or len(it["description"].split()) > limit:
            raise SystemExit(f"{path}: {it['content']}: a description of at most {limit} words")
    return cfg


def main(config: Path = CONFIG, out: Path = OUT) -> None:
    cfg = load_config(config)
    numbers = _numbers(cfg)
    if out.exists():
        shutil.rmtree(out)
    written: list[tuple[Path, dict]] = []
    for section in SECTIONS:
        folder = out / section
        for it in cfg.get(section, []):
            for p in MAKERS[it["make"]](it, folder, cfg):
                written.append((p, it))
            if section == "main" and it["kind"] == "table":      # its editable file (E8.11)
                written.append((main_table_docx(it, cfg, numbers, folder), it))
            print(f"[layout] {it['number']}: {it['content']}")
    supp_tables = [it for it in cfg["supplementary"] if it["kind"] == "table"]
    book = workbook(supp_tables, cfg, numbers, out / "supplementary_tables.xlsx")
    submission = submission_workbook(supp_tables, cfg, numbers, out / "submission" / "supplementary_tables.xlsx")
    width = cfg["artwork"]["double_column_mm"]
    contents = []
    for section in SECTIONS:
        for it in cfg.get(section, []):
            files = [p for p, i in written if i is it]
            estimate = ""
            if it["kind"] == "table" and section == "supplementary":
                fits, lines, wide = one_page(it["_frames"], width)
                estimate = f"{'yes' if fits else 'no'} ({lines:,} lines, {wide:,} mm wide)"
            contents.append({"number": it["number"], "section": section, "kind": it["kind"],
                             "content": it["content"], "title": grouped(_resolve(it["title"], numbers)),
                             "legend": grouped(_resolve(it.get("legend", ""), numbers)),
                             "destination": it.get("destination", ""),
                             "fits_one_page_estimate": estimate,
                             "description": grouped(_resolve(it.get("description", ""), numbers)),
                             "notes": " | ".join(_notes(it, numbers) + _part_notes(it, numbers)),
                             "files": "; ".join(_rel(p) for p in files),
                             "sources": "; ".join(_rel(s) for s in _sources(it)),
                             "make": it["make"] + (f" ({it['options']['function']})"
                                                   if it["make"] == "redraw_figure" else "")})
    for r in cfg.get("repository_only", []):
        contents.append({"number": "repository only", "section": "", "kind": "", "content": "",
                         "title": _resolve(r["what"], numbers), "legend": "", "destination": "",
                         "fits_one_page_estimate": "", "description": "", "notes": "",
                         "files": "", "sources": _rel(r["path"]), "make": ""})
    _write_csv(pd.DataFrame(contents), out / "contents.csv")
    manifest = []
    for p, it in written + [(book, None), (submission, None), (out / "contents.csv", None)]:
        srcs = (_sources(it) if it else [x for i in supp_tables for x in _sources(i)]
                if p in (book, submission) else [])
        manifest.append({"file": _rel(p), "bytes": p.stat().st_size, "sha256": _sha256(p),
                         "number": it["number"] if it else "", "content": it["content"] if it else "",
                         "sources": "; ".join(_rel(s) for s in dict.fromkeys(srcs)),
                         "sources_sha256": "; ".join(_sha256(s) for s in dict.fromkeys(srcs))})
    _write_csv(pd.DataFrame(manifest), out / "manifest.csv")
    print(f"[layout] wrote {len(written) + 4} files to {out}")


if __name__ == "__main__":
    main()
