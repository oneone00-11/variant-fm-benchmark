"""The layout layer (E5.3): the manuscript's figures and tables, main and supplementary,
under the numbers config/layout.yaml gives them.

Every figure and table here is drawn or written from an output an earlier stage wrote;
nothing is computed that a reader could not find there, and no earlier output is
touched: the v3.0.1 layout under figures/, tables/ and supplement/ stays as it was, so
that release can still be reproduced. One run writes them all, the main figures and
tables, the supplementary figures, every supplementary table as CSV and the workbook,
so that the text and its supplement come from one commit.

Numbers live only in config/layout.yaml, and in what this stage writes from it: the
workbook's README sheet and sheet names, contents.csv and manifest.csv, and the notes
that name another item (written {content} in the configuration). Every file is named
for its content, so reordering figures or tables once the text is final changes the
configuration alone. Titles, descriptions and notes in the configuration are
placeholders until then.

Figures follow the configuration's artwork block: no wider than its double-column width,
type no smaller and no larger than its sizes, bitmaps at its resolution. A larger size
set by evid_figures (its titles are 7.5 pt, its panel letters 9 pt) is brought down to
the largest; a smaller one stops the run. Supplementary Figures S1-S3 are v3.0.1's,
drawn again by their own functions in evid_figures, from the same outputs, at this
width; the files under figures/ are not touched.

Tables. A printed table (tables/, supplement/) is written as it stands, byte for byte.
Every other table keeps all its columns, under the readable names its entry in the
configuration gives, and each cell is printed by its column's type the way
Supplementary Tables S1-S13 print it (evid_supp_tables' formatters); the source files
keep their full precision. A table of several parts is written as one CSV per part. A
source's comment lines above its header are kept as the part's notes, in the workbook
above the part.

Outputs, under reports/evidence/layout/:
  main/<content>.pdf, .png                the main figures
  main/<content>[.<part>].csv             the main tables
  supplementary/<content>.pdf, .png       the supplementary figures
  supplementary/<content>[.<part>].csv    the supplementary tables
  supplementary_tables.xlsx               the supplementary tables, README sheet first
  contents.csv                            number, content, title, destination,
                                          description, notes, files and sources, and
                                          the outputs kept in the repository only
  manifest.csv                            every file written, with its sha256 and sources

Run (PYTHONPATH=phase1):  python -m src.evid_layout
"""
from __future__ import annotations

import hashlib
import io
import re
import shutil
import zipfile
from contextlib import contextmanager
from pathlib import Path
from xml.sax.saxutils import escape

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.text as mtext  # noqa: E402
import matplotlib.transforms as mtransforms  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

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
    return REPO_PREFIX + Path(path).as_posix()


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


# ---------------------------------------------------------------------------
# Figure: the published cut points, seven genes
# ---------------------------------------------------------------------------
@maker
def figure_published_cut_point(item: dict, folder: Path, cfg: dict) -> list[Path]:
    """v3.0.1's Figure 1 without the external genes, which the figure of the held-out
    and external genes draws instead."""
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
    fig, axes = plt.subplots(1, 2, figsize=(_width(cfg), 92 * F.MM), sharey=True,
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
        F._panel_letter(ax, letter, x=-0.01 if letter == "b" else -0.42)
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
    fig.legend(handles=handles, loc="lower center", ncol=len(handles),
               bbox_to_anchor=(0.6, 0.0), handletextpad=0.3, columnspacing=1.2)
    fig.subplots_adjust(left=0.2, right=0.97, top=0.94, bottom=0.2)
    return _save(fig, folder, item["content"], cfg)


# ---------------------------------------------------------------------------
# Figure: held-out genes and external genes
# ---------------------------------------------------------------------------
# genes by marker; TP53's assay differs from the HAP1 deposits (tp53_raw/
# scoreset_meta.json: HCT116 cells, Nutlin-3a), so it alone takes the second colour
EXTERNAL_GENES = ["DDX3X", "CTCF", "SFPQ", "RAD51D", "XRCC2", "SBDS", "TP53"]
GENE_MARKER = {"DDX3X": "o", "CTCF": "s", "SFPQ": "D", "RAD51D": "^", "XRCC2": "v",
               "SBDS": "p", "TP53": "X"}
# sizes that give the shapes about the same weight
MARKER_SIZE = {"o": 3.6, "s": 3.3, "D": 3.0, "^": 4.0, "v": 4.0, "p": 4.0, "X": 4.0}
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


def _heldout_panel(ax, folds: pd.DataFrame, thr: pd.DataFrame, tier: str, cut: float,
                   st: str, floor: float, ceil: float) -> bool:
    """v3.0.1's Figure 3, one panel: each held-out gene's ratio at the fitted threshold."""
    ax.set_xscale("log")
    ax.set_xlim(floor, ceil)
    ax.axvline(cut, color=F.INK2, lw=0.8, ls=(0, (1.5, 1.5)), zorder=0)
    any_zero = False
    for yy, c in enumerate(F.COLUMNS):
        if yy % 2 == 0:
            ax.axhspan(yy - 0.5, yy + 0.5, color=F.BAND, zorder=-1, lw=0)
        if not len(thr[(thr.tool == c) & (thr.stratum == st) & (thr.tier == tier)]):
            ax.text(floor * 1.15, yy, "not reached in-sample", fontsize=5.5, color=F.MUTED,
                    va="center", zorder=5, bbox=dict(facecolor="white", edgecolor="none", pad=0.4))
            continue
        v = folds[(folds.tool == c) & (folds.stratum == st) & (folds.tier == tier)]["heldout_lr"]
        v = v[np.isfinite(v)].to_numpy()
        if not len(v):
            ax.text(floor * 1.15, yy, "no evaluable fold", fontsize=5.5, color=F.MUTED,
                    va="center", zorder=5, bbox=dict(facecolor="white", edgecolor="none", pad=0.4))
            continue
        jitter = (np.arange(len(v)) - (len(v) - 1) / 2) * 0.09
        above, zero = v >= cut, v <= 0
        pos = np.clip(v, floor, ceil)
        ax.scatter(pos[above & ~zero], yy + jitter[above & ~zero], s=9, color=F.BLUE,
                   edgecolor="white", linewidth=0.3, zorder=3)
        ax.scatter(pos[~above & ~zero], yy + jitter[~above & ~zero], s=9, facecolor="white",
                   edgecolor=F.BLUE, linewidth=0.7, zorder=3)
        if zero.any():
            any_zero = True
            ax.scatter(np.full(zero.sum(), floor), yy + jitter[zero], s=9, facecolor="white",
                       edgecolor=F.INK2, linewidth=0.7, zorder=3, clip_on=False)
        ax.plot([np.median(pos)] * 2, [yy - 0.32, yy + 0.32], color=F.INK, lw=1.0, zorder=4)
    F._logfmt(ax)
    ax.tick_params(axis="y", length=0)
    return any_zero


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
        if yy % 2 == 0:
            ax.axhspan(yy - 0.5, yy + 0.5, color=F.BAND, zorder=-1, lw=0)
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
    any_zero = False
    letters = iter("abc")
    for i, (tier, cut) in enumerate((("moderate", 4.33), ("strong", 18.7))):
        axes = [fig.add_subplot(gs[i, j]) for j in range(3)]
        for j, (ax, st) in enumerate(zip(axes, F.STRATA)):
            any_zero |= _heldout_panel(ax, folds, thr, tier, cut, st, floor, ceil)
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
    if any_zero:
        handles.append(Line2D([], [], marker="o", ls="", mfc="white", mec=F.INK2, ms=4,
                              label="Ratio of zero (at the axis floor)"))
    genes = [Line2D([], [], marker=GENE_MARKER[g], ls="", color=GENE_COLOUR[g], mec="white",
                    ms=MARKER_SIZE[GENE_MARKER[g]] + 0.6,
                    label=F._it(g) + (f" ({TP53_ASSAY})" if g == "TP53" else ""))
             for g in EXTERNAL_GENES]
    genes.append(Line2D([], [], marker="o", ls="", mfc="white", mec=F.INK2, ms=4,
                        label="Open: below the boundary"))
    fig.legend(handles=handles + genes, loc="lower center", ncol=4,
               bbox_to_anchor=(0.55, 0.0), handletextpad=0.3, columnspacing=1.0)
    item["_notes"] = _missing_notes(missing)
    return _save(fig, folder, item["content"], cfg)


# ---------------------------------------------------------------------------
# Figures: ROC curves, precision-recall curves (v3.0.1's practical figure, one row each)
# ---------------------------------------------------------------------------
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
    b.set_xlabel("Signed intronic offset (acceptor side < 0 < donor side), nucleotides")
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


@maker
def redraw_figure(item: dict, folder: Path, cfg: dict) -> list[Path]:
    """A figure of v3.0.1, by its own function in evid_figures, from the same outputs."""
    F._style()
    kept: list = []
    with _drawn_by_figures(_width(cfg), kept):
        getattr(F, item["options"]["function"])()
    if len(kept) != 1:
        raise SystemExit(f"[layout] {item['content']}: evid_figures.{item['options']['function']} "
                         f"drew {len(kept)} figures, not one")
    return _save(kept[0], folder, item["content"], cfg)


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
# the codes the study's tables print as S1-S13 print them, beyond evid_supp_tables'
TOOL_EXTRA = {"fusion_enet": "Elastic-net combination"}
STRATUM_EXTRA = {"out_of_scope_deep": "Beyond 50 (out of scope)"}
TIER_LABEL = {"supporting": "Supporting", "moderate": "Moderate", "strong": "Strong",
              "very_strong": "Very strong", "below supporting": "Below Supporting"}
# the seven genes' label is each deposit's own classification, which S13 leaves blank
LABEL_BLANK = {"y_assay"}


def read_source(path: str | Path) -> tuple[pd.DataFrame, list[str]]:
    """A table as text cells, with the comment lines above its header as notes."""
    lines = Path(path).read_text().splitlines(keepends=True)
    notes = []
    while lines and lines[0].startswith("#"):
        notes.append(lines.pop(0).lstrip("#").strip())
    df = pd.read_csv(io.StringIO("".join(lines)), dtype=str, keep_default_na=False)
    return df, notes


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
    """One cell as S1-S13 print a cell of its type; an empty cell stays empty."""
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
        # an external gene's band is cut at the depth its deposit reaches, as S12 and S13
        # print it (DDX3X's 11-50 is 11-25)
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
        raise SystemExit(f"{where}: config/layout.yaml names no '{text}' under display: values: {kind}")
    return values[kind][text]


_BANDS: dict = {}


def _external_band(gene: str, stratum: str) -> str:
    if (gene, stratum) not in _BANDS:
        _BANDS[(gene, stratum)] = S._ext_band(gene, stratum)
    return _BANDS[(gene, stratum)]


def displayed(df: pd.DataFrame, item: dict, part: dict, cfg: dict) -> pd.DataFrame:
    """A source table as the supplement prints it: every column, under its readable name,
    each cell by its column's type."""
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
    return pd.DataFrame(out, columns=names)


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
    source is copied; any other is printed as displayed() prints it."""
    out = []
    item["_frames"] = []
    for part in item["parts"]:
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
        if printed:
            p.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(part["source"], p)
        else:
            df = displayed(df, item, part, cfg)
            _write_csv(df, p)
        out.append(p)
        item["_frames"].append((part, df, notes + list(part.get("notes", []))))
    return out


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
                 "P/LP": str(int(sub.n_plp.iloc[0])), "B/LB": str(int(sub.n_blb.iloc[0]))}
            for side, kind, name in cols:
                x = sub[(sub.side == side) & (sub.threshold_kind == kind)]
                if not len(x):
                    r[name] = "–"          # no such threshold: see the table's note
                    continue
                x = x.iloc[0]
                if x.status != "ok":
                    r[name] = "not evaluable"
                elif side == "pp3":
                    r[name] = f"{P._two(x.sensitivity)} ({int(x.n_pos_above)}/{int(x.n_plp)})"
                else:
                    r[name] = f"{P._two(x.bp4_sensitivity)} ({int(x.bp4_n_neg_at_or_below)}/{int(x.n_blb)})"
            rows.append(r)
    clinvar = pd.DataFrame(rows)
    # the bands left out of (b), from the same file: how many P/LP variants, and why
    notes = []
    for st in [s for s in P.CURVE_STRATA if s not in bands]:
        sub = c1[c1.tool.isin(tools) & (c1.stratum == st)]
        n_plp = sorted({int(v) for v in sub.n_plp})
        why = sorted({str(v) for v in sub.note.dropna()} | {str(v) for v in sub.status if v != "ok"})
        notes.append(f"(b) leaves out {F.STRATUM_LABEL[st]}: "
                     f"{' or '.join(map(str, n_plp))} P/LP variants of at least one star "
                     f"({'; '.join(why)})")
    item["_notes"] = notes
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
# prints a plain decimal is a number shown with the decimals it prints; any other is text.
# ---------------------------------------------------------------------------
DECIMAL = re.compile(r"^-?\d+(?:\.(\d+))?$")
ZIP_TIME = (1980, 1, 1, 0, 0, 0)
MAX_DECIMALS = 15
NUMBER_STYLE = 2         # cellXfs index of 0 decimals; d decimals at NUMBER_STYLE + d


def _col(i: int) -> str:
    s = ""
    i += 1
    while i:
        i, r = divmod(i - 1, 26)
        s = chr(65 + r) + s
    return s


def _cell_xml(ref: str, value, bold: bool) -> str:
    text = "" if value is None else str(value)
    if text == "":
        return ""
    m = DECIMAL.match(text)
    if not bold and m and len(text.lstrip("-").replace(".", "").lstrip("0")) <= 15:
        places = len(m.group(1) or "")
        if places <= MAX_DECIMALS:
            return f'<c r="{ref}" s="{NUMBER_STYLE + places}"><v>{text}</v></c>'
    style = ' s="1"' if bold else ""
    return f'<c r="{ref}"{style} t="inlineStr"><is><t xml:space="preserve">{escape(text)}</t></is></c>'


def _sheet_xml(rows: list[tuple[list, bool]]) -> str:
    body = []
    for r, (cells, bold) in enumerate(rows, start=1):
        xs = "".join(_cell_xml(f"{_col(c)}{r}", v, bold) for c, v in enumerate(cells))
        body.append(f'<row r="{r}">{xs}</row>')
    return ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
            '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
            f'<sheetData>{"".join(body)}</sheetData></worksheet>')


def _styles_xml(ns: str) -> str:
    formats = "".join(f'<numFmt numFmtId="{164 + d}" formatCode="{"0" if d == 0 else "0." + "0" * d}"/>'
                      for d in range(MAX_DECIMALS + 1))
    numbers = "".join(f'<xf numFmtId="{164 + d}" fontId="0" fillId="0" borderId="0" xfId="0" '
                      'applyNumberFormat="1"/>' for d in range(MAX_DECIMALS + 1))
    return ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
            f'<styleSheet xmlns="{ns}"><numFmts count="{MAX_DECIMALS + 1}">{formats}</numFmts>'
            '<fonts count="2"><font><sz val="10"/><name val="Arial"/></font>'
            '<font><b/><sz val="10"/><name val="Arial"/></font></fonts>'
            '<fills count="2"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill></fills>'
            '<borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>'
            '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
            f'<cellXfs count="{NUMBER_STYLE + MAX_DECIMALS + 1}"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>'
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


def workbook(items: list[dict], cfg: dict, numbers: dict, path: Path) -> Path:
    fmt = cfg["supplementary_format"]
    readme: list[tuple[list, bool]] = [
        (["Supplementary tables (provisional numbering, titles and descriptions: config/layout.yaml)"], True),
        ([_resolve(fmt["rule"], numbers)], False),
        ([], False),
        (["Table", "Sheet", "Content", "Title (provisional)", "Destination",
          "Description (provisional)", "Parts", "Source files", "Source sha256"], True)]
    sheets = []
    for it in items:
        srcs = _sources(it)
        readme.append(([it["number"], _sheet_name(it["number"]), it["content"], it["title"],
                        it["destination"], it["description"],
                        "; ".join(p["part"] for p in it["parts"]) if len(it["parts"]) > 1 else "",
                        "; ".join(_rel(s) for s in srcs), "; ".join(_sha256(s) for s in srcs)], False))
        rows: list[tuple[list, bool]] = [([f"{it['number']}. {it['title']}"], True)]
        for n in _notes(it, numbers):
            rows.append(([n], False))
        for k, (part, df, notes) in enumerate(it["_frames"]):
            rows.append(([], False))
            if len(it["parts"]) > 1:
                rows.append(([f"({chr(97 + k)}) {part.get('title', part['part'])}"], True))
            for n in notes:
                rows.append(([_resolve(n, numbers)], False))
            rows.append((list(df.columns), True))
            rows += [(list(r), False) for r in df.itertuples(index=False)]
        sheets.append((_sheet_name(it["number"]), rows))
    readme.append(([], False))
    readme.append((["Kept in the repository only (not numbered)"], True))
    readme.append((["Path", "What it holds"], True))
    readme += [([_rel(r["path"]), _resolve(r["what"], numbers)], False)
               for r in cfg.get("repository_only", [])]
    return write_xlsx(path, [("README", readme)] + sheets)


# ---------------------------------------------------------------------------
# The run
# ---------------------------------------------------------------------------
REFERENCE = re.compile(r"\{([a-z0-9_]+)\}")


def _resolve(text: str, numbers: dict) -> str:
    """A note's {content} references, as the items' numbers."""
    def number(m):
        if m.group(1) not in numbers:
            raise SystemExit(f"config/layout.yaml: a note names '{m.group(1)}', which is no item")
        return numbers[m.group(1)]
    return REFERENCE.sub(number, text)


def _notes(item: dict, numbers: dict) -> list[str]:
    return [_resolve(n, numbers) for n in list(item.get("notes", [])) + item.get("_notes", [])]


def _part_notes(item: dict, numbers: dict) -> list[str]:
    """The configuration's notes on a part, named by the part (the workbook writes them
    above it)."""
    return [f"{p['part']}: {_resolve(n, numbers)}" for p in item.get("parts", [])
            for n in p.get("notes", [])]


def _sources(item: dict) -> list[str]:
    out = list(item.get("sources", []))
    for part in item.get("parts", []):
        if "source" in part:
            out.append(part["source"])
        srcs = part.get("sources", [])
        out += list(srcs.values()) if isinstance(srcs, dict) else list(srcs)
    return list(dict.fromkeys(out))


def load_config(path: Path = CONFIG) -> dict:
    cfg = yaml.safe_load(path.read_text())
    items = cfg["main"] + cfg["supplementary"]
    for key in ("content", "number"):
        seen = [it[key] for it in items]
        twice = sorted({c for c in seen if seen.count(c) > 1})
        if twice:
            raise SystemExit(f"{path}: {key}s used twice: {twice}")
    missing = [s for it in items for s in _sources(it) if not Path(s).is_file()]
    missing += [r["path"] for r in cfg.get("repository_only", []) if not Path(r["path"]).is_file()]
    if missing:
        raise SystemExit(f"{path}: sources not found: {missing}")
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
    numbers = {it["content"]: it["number"] for it in cfg["main"] + cfg["supplementary"]}
    if out.exists():
        shutil.rmtree(out)
    written: list[tuple[Path, dict]] = []
    for section in ("main", "supplementary"):
        folder = out / section
        for it in cfg[section]:
            for p in MAKERS[it["make"]](it, folder, cfg):
                written.append((p, it))
            print(f"[layout] {it['number']}: {it['content']}")
    supp_tables = [it for it in cfg["supplementary"] if it["kind"] == "table"]
    book = workbook(supp_tables, cfg, numbers, out / "supplementary_tables.xlsx")
    width = cfg["artwork"]["double_column_mm"]
    contents = []
    for section in ("main", "supplementary"):
        for it in cfg[section]:
            files = [p for p, i in written if i is it]
            estimate = ""
            if it["kind"] == "table" and section == "supplementary":
                fits, lines, wide = one_page(it["_frames"], width)
                estimate = f"{'yes' if fits else 'no'} ({lines} lines, {wide} mm wide)"
            contents.append({"number": it["number"], "section": section, "kind": it["kind"],
                             "content": it["content"], "title": it["title"],
                             "destination": it.get("destination", ""),
                             "fits_one_page_estimate": estimate,
                             "description": it.get("description", ""),
                             "notes": " | ".join(_notes(it, numbers) + _part_notes(it, numbers)),
                             "files": "; ".join(_rel(p) for p in files),
                             "sources": "; ".join(_rel(s) for s in _sources(it)),
                             "make": it["make"] + (f" ({it['options']['function']})"
                                                   if it["make"] == "redraw_figure" else "")})
    for r in cfg.get("repository_only", []):
        contents.append({"number": "repository only", "section": "", "kind": "", "content": "",
                         "title": _resolve(r["what"], numbers), "destination": "",
                         "fits_one_page_estimate": "", "description": "", "notes": "",
                         "files": "", "sources": _rel(r["path"]), "make": ""})
    _write_csv(pd.DataFrame(contents), out / "contents.csv")
    manifest = []
    for p, it in written + [(book, None), (out / "contents.csv", None)]:
        srcs = _sources(it) if it else [x for i in supp_tables for x in _sources(i)] if p == book else []
        manifest.append({"file": _rel(p), "bytes": p.stat().st_size, "sha256": _sha256(p),
                         "number": it["number"] if it else "", "content": it["content"] if it else "",
                         "sources": "; ".join(_rel(s) for s in dict.fromkeys(srcs)),
                         "sources_sha256": "; ".join(_sha256(s) for s in dict.fromkeys(srcs))})
    _write_csv(pd.DataFrame(manifest), out / "manifest.csv")
    print(f"[layout] wrote {len(written) + 3} files to {out}")


if __name__ == "__main__":
    main()
