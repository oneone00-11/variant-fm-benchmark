"""Final figures for the evidence-strength study, main and supplementary.

Every panel is drawn from a table the pipeline has already written; nothing is
computed here that a reader could not find in reports/evidence/. The draft PNGs
that evid_fig_data writes are layout sketches and are not used.

Design choices that change what a reader sees, recorded so they are not mistaken
for accidents:

  * Likelihood ratios are on a logarithmic axis, because the tier boundaries are
    multiplicative (2.08, 4.33, 18.7) and a linear axis would compress the whole
    Supporting-to-Moderate range into a sliver.
  * Tier boundaries are drawn as dotted reference lines and named once, at the
    top of each panel; they are the reading frame of every figure.
  * Only the in-scope strata appear. The canonical dinucleotides are outside the
    PP3/BP4 recommendation and are reported in the supplementary tables instead.
  * Colour carries one job per figure: ClinVar record status in Figure 1, label
    definition in Figure 4. The three colours are the first three slots of a
    palette validated for colour-vision deficiency on all pairs, and every
    coloured mark has a second channel (position or marker shape), so no reading
    depends on hue alone.
  * A likelihood ratio of zero cannot sit on a log axis. Such points are drawn
    at the axis floor as open markers and the legend says so.

Output is deterministic: PDF and PNG metadata that would carry a timestamp or a
software version are cleared, so a rerun on the same machine writes identical
files. Font rendering depends on the fonts installed, so byte identity across
machines is not claimed; the content is.

Run (PYTHONPATH=phase1):  python -m src.evid_figures
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

REPORT_DIR = Path("reports/evidence")
FIG_DATA = REPORT_DIR / "fig_data"
OUT = REPORT_DIR / "figures"

MM = 1 / 25.4
DOUBLE = 183 * MM          # double-column width

INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a8984"
GRID = "#d9d8d4"
BAND = "#f3f2ef"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"

PATH_CUTS = [("Supporting", 2.08), ("Moderate", 4.33), ("Strong", 18.7)]
BEN_CUTS = [("Supporting", 1 / 2.08), ("Moderate", 1 / 4.33), ("Strong", 1 / 18.7)]

STRATA = ["s3_10", "s11_50", "s3_50"]
STRATUM_LABEL = {"s3_10": "3–10 bp", "s11_50": "11–50 bp", "s3_50": "3–50 bp"}
ARMS = ["all", "classified", "recorded_unclassified", "unrecorded"]
ARM_LABEL = {"all": "All variants", "classified": "Classified in ClinVar",
             "recorded_unclassified": "Recorded, unclassified",
             "unrecorded": "Not in ClinVar"}
ARM_COLOUR = {"all": INK, "classified": BLUE, "recorded_unclassified": ORANGE,
              "unrecorded": AQUA}

# display order and names for score columns, shared with the tables
COLUMNS = ["spliceai_walker", "spliceai", "pangolin", "alphagenome", "avi",
           "avi_splice_sites", "avi_splice_site_usage", "avi_splice_junctions",
           "cadd", "phylop", "phastcons", "gpn_msa", "nt"]
NAME = {
    "spliceai_walker": "SpliceAI (published basis)",
    "spliceai": "SpliceAI (distance 50)",
    "pangolin": "Pangolin",
    "alphagenome": "AlphaGenome splice score",
    "avi": "Atlas combined score",
    "avi_splice_sites": "Atlas splice sites",
    "avi_splice_site_usage": "Atlas splice-site usage",
    "avi_splice_junctions": "Atlas splice junctions",
    "cadd": "CADD",
    "phylop": "phyloP",
    "phastcons": "phastCons",
    "gpn_msa": "GPN-MSA",
    "nt": "Nucleotide Transformer",
}
SPLICE_AWARE = ["spliceai_walker", "pangolin", "alphagenome", "avi",
                "avi_splice_sites", "avi_splice_site_usage", "avi_splice_junctions"]
LABEL_DEF = {"y_control_anchored": "Control-anchored labels",
             "y_fdr": "FDR labels"}
# The external genes' two label definitions are told apart by fill, not colour:
# filled for control-anchored, open for FDR. Colour already carries ClinVar record
# status in Figure 1, and one colour must not mean two things in one figure.
LABEL_FILLED = {"y_control_anchored": True, "y_fdr": False}


def _style() -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7,
        "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
        "axes.edgecolor": INK2, "axes.linewidth": 0.6, "axes.labelcolor": INK,
        "xtick.color": INK2, "ytick.color": INK2, "xtick.major.width": 0.6,
        "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "pdf.fonttype": 42, "svg.hashsalt": "evidence",
        "legend.frameon": False,
    })


def _save(fig, stem: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{stem}.pdf", metadata={"CreationDate": None, "ModDate": None,
                                                "Producer": None, "Creator": None})
    fig.savefig(OUT / f"{stem}.png", dpi=300, metadata={"Software": None})
    plt.close(fig)
    print(f"[figures] wrote {stem}.pdf / .png")


def _panel_letter(ax, letter: str, x: float = -0.02) -> None:
    ax.text(x, 1.02, letter, transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="bottom", ha="right", color=INK)


def _tier_lines(ax, cuts, ymax_frac: float = 1.0, label: bool = True,
                vertical: bool = True) -> None:
    for name, v in cuts:
        if vertical:
            ax.axvline(v, color=MUTED, lw=0.6, ls=(0, (1.5, 1.5)), zorder=0)
            if label:
                ax.text(v, 1.0, name, transform=ax.get_xaxis_transform(), fontsize=6,
                        color=INK2, ha="left", va="bottom", rotation=0)
        else:
            ax.axhline(v, color=MUTED, lw=0.6, ls=(0, (1.5, 1.5)), zorder=0)
            if label:
                ax.text(1.0, v, " " + name, transform=ax.get_yaxis_transform(),
                        fontsize=6, color=INK2, ha="left", va="center")


def _logfmt(ax, axis: str = "x") -> None:
    fmt = matplotlib.ticker.FuncFormatter(
        lambda v, _: (f"{v:g}" if v >= 0.01 else f"{v:.0e}"))
    getattr(ax, f"{axis}axis").set_major_formatter(fmt)
    getattr(ax, f"{axis}axis").set_minor_formatter(matplotlib.ticker.NullFormatter())
    getattr(ax, f"{axis}axis").set_minor_locator(matplotlib.ticker.NullLocator())


# ---------------------------------------------------------------------------
# Figure 1 -- the published cut points by band and ClinVar record status
# ---------------------------------------------------------------------------
def figure1() -> None:
    d = pd.read_csv(FIG_DATA / "fig1_walker_cutpoints.csv")
    seven = d[(d.source == "seven_genes") & (d.gene == "seven_genes")]
    ext = d[(d.source == "external")
            & ((d.gene == "DDX3X") | ((d.gene == "TP53")
                                      & (d.label_definition == "y_control_anchored")
                                      & (d.stratum == "s3_10")))]

    rows = []           # (y, label, group header or None, record)
    y = 0.0
    for st in STRATA:
        rows.append((y, None, f"Seven genes, {STRATUM_LABEL[st]}", None)); y += 1
        for arm in ARMS:
            r = seven[(seven.stratum == st) & (seven.arm == arm)]
            if len(r):
                rows.append((y, ARM_LABEL[arm], None, r.iloc[0])); y += 1
        y += 0.4
    rows.append((y, None, "DDX3X, all variants", None)); y += 1
    for st in STRATA:
        for ld in ("y_control_anchored", "y_fdr"):
            r = ext[(ext.gene == "DDX3X") & (ext.stratum == st) & (ext.label_definition == ld)]
            if len(r):
                short = "control-anchored" if ld == "y_control_anchored" else "FDR"
                rows.append((y, f"{STRATUM_LABEL[st]}, {short}", None, r.iloc[0])); y += 1
    y += 0.4
    rows.append((y, None, "TP53, all variants", None)); y += 1
    r = ext[ext.gene == "TP53"]
    if len(r):
        rows.append((y, "3–10 bp, control-anchored", None, r.iloc[0])); y += 1
    ymax = y

    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE, 118 * MM), sharey=True,
                             gridspec_kw={"wspace": 0.08})
    specs = [("lr_pp3", "lr_pp3_lo", "lr_pp3_hi", PATH_CUTS, (0.8, 60),
              "Likelihood ratio of the PP3 band (score ≥ 0.2)"),
             ("lr_le01", "lr_le01_lo", "lr_le01_hi", BEN_CUTS, (0.02, 1.25),
              "Likelihood ratio of the BP4 band (score ≤ 0.1);\nsmaller values are stronger benign evidence")]
    for ax, (col, lo, hi, cuts, xlim, xlabel), letter in zip(axes, specs, "ab"):
        ax.set_xscale("log")
        ax.set_xlim(*xlim)
        _tier_lines(ax, cuts)
        for yy, lab, head, rec in rows:
            if rec is None:
                continue
            arm = rec["arm"] if rec["source"] == "seven_genes" else "all"
            colour = ARM_COLOUR.get(arm, INK) if rec["source"] == "seven_genes" else INK2
            filled = (rec["source"] == "seven_genes"
                      or LABEL_FILLED.get(rec["label_definition"], True))
            v = rec[col]
            if not np.isfinite(v):
                continue
            if np.isfinite(rec.get(lo, np.nan)) and np.isfinite(rec.get(hi, np.nan)):
                ax.plot([max(rec[lo], xlim[0]), min(rec[hi], xlim[1])], [yy, yy],
                        color=colour, lw=1.0, solid_capstyle="butt", zorder=2)
            if v <= 0:
                ax.plot(xlim[0], yy, marker="o", ms=4.2, mfc="white", mec=colour,
                        mew=0.9, zorder=3, clip_on=False)
            else:
                marker = "o" if rec["source"] == "seven_genes" else "D"
                ax.plot(v, yy, marker=marker, ms=4.2 if marker == "o" else 3.8,
                        color=colour, mfc=colour if filled else "white",
                        mec="white" if filled else colour, mew=0.5 if filled else 0.9,
                        zorder=3)
        ax.set_ylim(ymax - 0.4, -0.8)
        ax.set_xlabel(xlabel)
        _logfmt(ax)
        ax.tick_params(axis="y", length=0)
        _panel_letter(ax, letter, x=-0.01 if letter == "b" else -0.58)
    ax0 = axes[0]
    ax0.set_yticks([r[0] for r in rows if r[1] is not None])
    ax0.set_yticklabels([r[1] for r in rows if r[1] is not None])
    for yy, lab, head, rec in rows:
        if head:
            ax0.text(-0.56, yy, head, transform=ax0.get_yaxis_transform(),
                     fontsize=7, fontweight="bold", color=INK, va="center", ha="left")
    handles = [Line2D([], [], marker="o", ls="", color=ARM_COLOUR[a], mec="white",
                      label=ARM_LABEL[a], ms=4.5) for a in ARMS]
    handles += [Line2D([], [], marker="D", ls="", color=INK2,
                       mfc=INK2 if LABEL_FILLED[k] else "white",
                       mec="white" if LABEL_FILLED[k] else INK2,
                       label=f"External gene, {v.replace('Control', 'control')}", ms=4)
                for k, v in LABEL_DEF.items()]
    if any(rec is not None and (rec["lr_pp3"] == 0 or rec["lr_le01"] == 0)
           for _, _, _, rec in rows):
        handles += [Line2D([], [], marker="o", ls="", mfc="white", mec=INK2,
                           label="Ratio of zero (drawn at the axis floor)", ms=4.5)]
    fig.legend(handles=handles, loc="lower center", ncol=3, bbox_to_anchor=(0.6, -0.005),
               handletextpad=0.3, columnspacing=1.2)
    fig.subplots_adjust(left=0.27, right=0.97, top=0.95, bottom=0.19)
    _save(fig, "figure1")


# ---------------------------------------------------------------------------
# Figure 2 / Supplementary Figure S1 -- local likelihood ratio against score
# ---------------------------------------------------------------------------
def _curve(tool: str, stratum: str) -> pd.DataFrame | None:
    p = REPORT_DIR / f"interval_lr_{tool}_{stratum}.csv"
    if not p.exists():
        return None
    c = pd.read_csv(p).sort_values("score")
    return c[np.isfinite(c.local_lr)]


def _thresholds() -> pd.DataFrame:
    e = pd.read_csv(REPORT_DIR / "evidence_thresholds.csv")
    return e[(e.status == "ok") & e.pp3_threshold_reachable.astype(bool)]


def _curve_panel(ax, tool: str, stratum: str, thr: pd.DataFrame, floor: float,
                 ceil: float, walker_mark: bool) -> None:
    c = _curve(tool, stratum)
    ax.set_yscale("log")
    ax.set_ylim(floor, ceil)
    _tier_lines(ax, PATH_CUTS, vertical=False, label=False)
    ax.axhline(1.0, color=GRID, lw=0.6, zorder=0)
    if c is None or not len(c):
        ax.text(0.5, 0.5, "not scored", transform=ax.transAxes, ha="center",
                color=MUTED)
        return
    lo = c["lr_lo_gene_boot"].clip(lower=floor).to_numpy()
    hi = c["lr_hi_gene_boot"].clip(upper=ceil).to_numpy()
    ax.fill_between(c.score, lo, hi, color=BLUE, alpha=0.18, lw=0, zorder=1)
    ax.plot(c.score, c.local_lr.clip(lower=floor, upper=ceil), color=BLUE, lw=1.0,
            zorder=2)
    t = thr[(thr.tool == tool) & (thr.stratum == stratum)]
    for tier, ls in (("moderate", (0, (3, 1.5))), ("strong", "-")):
        r = t[t.tier == tier]
        if len(r):
            x = float(r.pp3_threshold_insample.iloc[0])
            ax.axvline(x, color=INK2, lw=0.7, ls=ls, zorder=1)
    if walker_mark:
        ax.axvline(0.2, color=ORANGE, lw=0.9, zorder=1)
    _logfmt(ax, "y")


def figure2() -> None:
    tools = ["spliceai_walker", "pangolin", "alphagenome", "avi"]
    thr = _thresholds()
    fig, axes = plt.subplots(len(tools), 3, figsize=(DOUBLE, 150 * MM),
                             gridspec_kw={"hspace": 0.55, "wspace": 0.18})
    for i, tool in enumerate(tools):
        for j, st in enumerate(STRATA):
            ax = axes[i, j]
            _curve_panel(ax, tool, st, thr, 0.05, 1000, tool == "spliceai_walker")
            if i == 0:
                ax.set_title(STRATUM_LABEL[st], fontweight="bold", pad=8)
            if j == 0:
                ax.set_ylabel(f"{NAME[tool]}\nlocal likelihood ratio")
            else:
                ax.set_yticklabels([])
            if j == 2:
                for name, v in PATH_CUTS:
                    ax.text(1.01, v, name, transform=ax.get_yaxis_transform(),
                            fontsize=5.8, color=INK2,
                            va={"Supporting": "top", "Moderate": "bottom"}.get(name, "center"))
    handles = [Line2D([], [], color=BLUE, lw=1.0, label="Local likelihood ratio"),
               matplotlib.patches.Patch(color=BLUE, alpha=0.18, lw=0,
                                        label="95% interval, gene-clustered bootstrap"),
               Line2D([], [], color=INK2, lw=0.7, ls=(0, (3, 1.5)),
                      label="Fitted Moderate threshold"),
               Line2D([], [], color=INK2, lw=0.7, label="Fitted Strong threshold"),
               Line2D([], [], color=ORANGE, lw=0.9, label="Published cut point (0.2)")]
    fig.legend(handles=handles, loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.0))
    fig.text(0.5, 0.075, "Score, on each predictor's own scale", ha="center", fontsize=7)
    fig.subplots_adjust(left=0.1, right=0.91, top=0.95, bottom=0.13)
    _save(fig, "figure2")


def figure_s1() -> None:
    thr = _thresholds()
    fig, axes = plt.subplots(len(COLUMNS), 3, figsize=(DOUBLE, 250 * MM),
                             gridspec_kw={"hspace": 0.75, "wspace": 0.18})
    for i, tool in enumerate(COLUMNS):
        for j, st in enumerate(STRATA):
            ax = axes[i, j]
            _curve_panel(ax, tool, st, thr, 0.05, 1000, tool in ("spliceai", "spliceai_walker"))
            if i == 0:
                ax.set_title(STRATUM_LABEL[st], fontweight="bold", pad=6)
            if j == 0:
                ax.set_ylabel(NAME[tool], rotation=0, ha="right", va="center", fontsize=6.5)
            else:
                ax.set_yticklabels([])
            ax.tick_params(labelsize=5.5)
    fig.text(0.5, 0.005, "Score (each column on its own scale); y axis: local likelihood "
             "ratio, log scale, with Supporting, Moderate and Strong boundaries dotted",
             ha="center", fontsize=6.5, color=INK2)
    fig.subplots_adjust(left=0.2, right=0.98, top=0.97, bottom=0.03)
    _save(fig, "figureS1")


# ---------------------------------------------------------------------------
# Figure 3 -- held-out likelihood ratios, one point per held-out gene
# ---------------------------------------------------------------------------
def figure3() -> None:
    f = pd.read_csv(REPORT_DIR / "evidence_thresholds_logo_folds.csv")
    f = f[(f.side == "pp3") & f.stratum.isin(STRATA) & f.tool.isin(COLUMNS)]
    thr = _thresholds()
    tiers = [("moderate", "Moderate", 4.33), ("strong", "Strong", 18.7)]
    floor, ceil = 0.3, 1000
    fig, axes = plt.subplots(2, 3, figsize=(DOUBLE, 150 * MM), sharey=True,
                             gridspec_kw={"hspace": 0.3, "wspace": 0.08})
    ypos = {c: i for i, c in enumerate(COLUMNS)}
    any_zero = False
    for i, (tier, tname, cut) in enumerate(tiers):
        for j, st in enumerate(STRATA):
            ax = axes[i, j]
            ax.set_xscale("log")
            ax.set_xlim(floor, ceil)
            ax.axvline(cut, color=INK2, lw=0.8, ls=(0, (1.5, 1.5)), zorder=0)
            for c in COLUMNS:
                yy = ypos[c]
                reached = len(thr[(thr.tool == c) & (thr.stratum == st) & (thr.tier == tier)])
                if not reached:
                    ax.text(floor * 1.15, yy, "not reached in-sample", fontsize=5.5,
                            color=MUTED, va="center")
                    continue
                v = f[(f.tool == c) & (f.stratum == st) & (f.tier == tier)]["heldout_lr"]
                v = v[np.isfinite(v)].to_numpy()
                if not len(v):
                    ax.text(floor * 1.15, yy, "no evaluable fold", fontsize=5.5,
                            color=MUTED, va="center")
                    continue
                jitter = (np.arange(len(v)) - (len(v) - 1) / 2) * 0.09
                above = v >= cut
                pos = np.clip(v, floor, ceil)
                zero = v <= 0
                ax.scatter(pos[above & ~zero], yy + jitter[above & ~zero], s=9,
                           color=BLUE, edgecolor="white", linewidth=0.3, zorder=3)
                ax.scatter(pos[~above & ~zero], yy + jitter[~above & ~zero], s=9,
                           facecolor="white", edgecolor=BLUE, linewidth=0.7, zorder=3)
                if zero.any():
                    any_zero = True
                    ax.scatter(np.full(zero.sum(), floor), yy + jitter[zero], s=9,
                               facecolor="white", edgecolor=INK2, linewidth=0.7,
                               zorder=3, clip_on=False)
                ax.plot([np.median(pos)] * 2, [yy - 0.32, yy + 0.32], color=INK,
                        lw=1.0, zorder=4)
            if i == 0:
                ax.set_title(STRATUM_LABEL[st], fontweight="bold")
            _logfmt(ax)
            ax.tick_params(axis="y", length=0)
            for k in range(len(COLUMNS)):
                if k % 2 == 0:
                    ax.axhspan(k - 0.5, k + 0.5, color=BAND, zorder=-1, lw=0)
    axes[0, 0].set_yticks(range(len(COLUMNS)))
    axes[0, 0].set_yticklabels([NAME[c] for c in COLUMNS])
    axes[0, 0].set_ylim(len(COLUMNS) - 0.5, -0.5)
    _panel_letter(axes[0, 0], "a", x=-0.62)
    _panel_letter(axes[1, 0], "b", x=-0.62)
    for i, (tier, tname, cut) in enumerate(tiers):
        mid = axes[i, 1]
        mid.text(0.5, -0.13, f"Likelihood ratio in the held-out gene at the fitted {tname} "
                 f"threshold (dotted line: {tname} boundary, {cut:g})",
                 transform=mid.transAxes, ha="center", va="top", fontsize=7)
    handles = [Line2D([], [], marker="o", ls="", color=BLUE, mec="white", ms=4,
                      label="Held-out gene, ratio clears the cut"),
               Line2D([], [], marker="o", ls="", mfc="white", mec=BLUE, ms=4,
                      label="Held-out gene, ratio below the cut"),
               Line2D([], [], color=INK, lw=1.0, label="Median across held-out genes")]
    if any_zero:
        handles.insert(2, Line2D([], [], marker="o", ls="", mfc="white", mec=INK2, ms=4,
                                 label="Ratio of zero (drawn at the axis floor)"))
    fig.legend(handles=handles, loc="lower center", ncol=4, bbox_to_anchor=(0.55, 0.0),
               handletextpad=0.3, columnspacing=1.0)
    fig.subplots_adjust(left=0.19, right=0.98, top=0.96, bottom=0.1)
    _save(fig, "figure3")


# ---------------------------------------------------------------------------
# Figure 4 -- thresholds carried to genes outside the seven
# ---------------------------------------------------------------------------
def figure4() -> None:
    d = pd.read_csv(FIG_DATA / "fig4_external.csv")
    basis = {g: pd.read_csv(REPORT_DIR / f"external_{g.lower()}_column_basis.csv")
             for g in ("DDX3X", "TP53")}
    marks = {"walker": ("o", "Published cut point (0.2)"),
             "moderate": ("s", "Fitted Moderate threshold"),
             "strong": ("^", "Fitted Strong threshold")}
    cols = ["spliceai_walker", "spliceai", "pangolin", "alphagenome", "avi",
            "avi_splice_sites", "avi_splice_site_usage", "avi_splice_junctions"]
    panels = [("DDX3X", st, ["y_control_anchored", "y_fdr"]) for st in STRATA]
    panels.append(("TP53", "s3_10", ["y_control_anchored"]))
    floor, ceil = 0.8, 300
    fig, axes = plt.subplots(2, 2, figsize=(DOUBLE, 150 * MM), sharey=True,
                             gridspec_kw={"wspace": 0.06, "hspace": 0.42})
    ypos = {c: i for i, c in enumerate(cols)}
    for ax, (gene, st, lds), letter in zip(axes.ravel(), panels, "abcd"):
        sub = d[(d.gene == gene) & (d.stratum == st) & d.label_definition.isin(lds)]
        ax.set_xscale("log")
        ax.set_xlim(floor, ceil)
        for name, v in PATH_CUTS:
            ax.axvline(v, color=MUTED, lw=0.6, ls=(0, (1.5, 1.5)), zorder=0)
            ax.text(v * 1.06, -0.45, name, fontsize=5.8, color=INK2, rotation=90,
                    ha="left", va="top")
        for k in range(len(cols)):
            if k % 2 == 0:
                ax.axhspan(k - 0.5, k + 0.5, color=BAND, zorder=-1, lw=0)
        b = basis[gene].set_index("tool")
        for c in cols:
            r = sub[sub.tool == c]
            if not len(r) or not bool(r.column_basis_comparable.iloc[0]):
                if c in b.index and not bool(b.loc[c, "comparable"]):
                    note = "excluded by the column-basis check"
                elif gene == "TP53" and c.startswith("avi"):
                    note = "not evaluated in this gene"
                else:
                    note = "not scored in this gene"
                ax.text(ceil / 1.1, ypos[c], note, fontsize=5.6, color=MUTED,
                        ha="right", va="center")
                continue
            for li, ld in enumerate(lds):
                rr = r[r.label_definition == ld]
                if not len(rr):
                    continue
                dy = (li - (len(lds) - 1) / 2) * 0.28
                filled = LABEL_FILLED[ld]
                style = dict(ls="", ms=4.4, color=BLUE, mfc=BLUE if filled else "white",
                             mec="white" if filled else BLUE, mew=0.4 if filled else 0.9,
                             zorder=3)
                w = rr.iloc[0]["lr_at_walker_cut"]
                if np.isfinite(w):
                    ax.plot(w, ypos[c] + dy, marker=marks["walker"][0], **style)
                for tier in ("moderate", "strong"):
                    t = rr[rr.tier == tier]
                    if len(t) and np.isfinite(t.iloc[0]["lr_at_logo_threshold"]):
                        ax.plot(t.iloc[0]["lr_at_logo_threshold"], ypos[c] + dy,
                                marker=marks[tier][0], **style)
        ax.set_title(f"{gene}, {STRATUM_LABEL[st]}", fontweight="bold", loc="left")
        ax.set_xlabel("Likelihood ratio in the external gene")
        _logfmt(ax)
        ax.tick_params(axis="y", length=0)
        _panel_letter(ax, letter, x=-0.02 if letter in "bd" else -0.45)
    axes[0, 0].set_yticks(range(len(cols)))
    axes[0, 0].set_yticklabels([NAME[c] for c in cols])
    axes[0, 0].set_ylim(len(cols) - 0.5, -0.5)
    axes[1, 0].set_yticklabels([NAME[c] for c in cols])
    handles = [Line2D([], [], marker=m, ls="", color=INK2, mec="white", ms=4.4, label=lab)
               for m, lab in marks.values()]
    handles += [matplotlib.patches.Patch(facecolor=BLUE, edgecolor=BLUE,
                                         label="Filled: control-anchored labels"),
                matplotlib.patches.Patch(facecolor="white", edgecolor=BLUE, lw=0.9,
                                         label="Open: FDR labels")]
    fig.legend(handles=handles, loc="lower center", ncol=3, bbox_to_anchor=(0.58, 0.0),
               handlelength=1.0)
    fig.subplots_adjust(left=0.19, right=0.98, top=0.95, bottom=0.14)
    _save(fig, "figure4")


# ---------------------------------------------------------------------------
# Supplementary Figure S3 -- AUROC by ClinVar record status
# ---------------------------------------------------------------------------
def figure_s3() -> None:
    a = pd.read_csv(REPORT_DIR / "territory_metrics_arms_noBRCA1.csv")
    a = a[a.status == "ok"]
    pools = [("s3_50", "all_genes", "3–50 bp, seven genes"),
             ("s3_50", "no_BRCA1", "3–50 bp, without BRCA1"),
             ("all_1_50", "all_genes", "1–50 bp including ±1,2 (out of scope)")]
    fig, axes = plt.subplots(1, 3, figsize=(DOUBLE, 95 * MM), sharey=True,
                             gridspec_kw={"wspace": 0.14})
    arms = ["classified", "recorded_unclassified", "unrecorded"]
    ypos = {c: i for i, c in enumerate(COLUMNS)}
    for ax, (st, gs, title) in zip(axes, pools):
        sub = a[(a.stratum == st) & (a.gene_set == gs)]
        for k in range(len(COLUMNS)):
            if k % 2 == 0:
                ax.axhspan(k - 0.5, k + 0.5, color=BAND, zorder=-1, lw=0)
        for ai, arm in enumerate(arms):
            dy = (ai - 1) * 0.25
            s2 = sub[sub.clinvar_arm == arm]
            for c in COLUMNS:
                r = s2[s2.tool == c]
                if not len(r):
                    continue
                r = r.iloc[0]
                ax.plot([r.auroc_lo, r.auroc_hi], [ypos[c] + dy] * 2, color=ARM_COLOUR[arm],
                        lw=0.8, zorder=2)
                ax.plot(r.auroc, ypos[c] + dy, "o", ms=3.2, color=ARM_COLOUR[arm],
                        mec="white", mew=0.3, zorder=3)
        ax.set_title(title, fontweight="bold")
        ax.set_xlim(0.5, 1.0)
        ax.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
        ax.set_xlabel("AUROC, pooled across genes (95% CI)")
        ax.tick_params(axis="y", length=0)
    axes[0].set_yticks(range(len(COLUMNS)))
    axes[0].set_yticklabels([NAME[c] for c in COLUMNS])
    axes[0].set_ylim(len(COLUMNS) - 0.5, -0.5)
    handles = [Line2D([], [], marker="o", ls="-", lw=0.8, color=ARM_COLOUR[a2],
                      mec="white", ms=4, label=ARM_LABEL[a2]) for a2 in arms]
    fig.legend(handles=handles, loc="lower center", ncol=3, bbox_to_anchor=(0.55, 0.0))
    fig.subplots_adjust(left=0.19, right=0.98, top=0.93, bottom=0.17)
    _save(fig, "figureS3")


# ---------------------------------------------------------------------------
# Supplementary Figure S2 -- how many training genes a Strong threshold needs
# ---------------------------------------------------------------------------
# Four lines per panel, so four categorical slots in the palette's fixed order;
# the palette validates its first four slots on the adjacent pairlist that line
# charts use, and every line also carries a distinct marker and a direct label.
DILUTION_TOOLS = [("spliceai_walker", BLUE, "o"), ("pangolin", ORANGE, "s"),
                  ("alphagenome", AQUA, "^"), ("avi", "#eda100", "D")]


def figure_s2() -> None:
    p = REPORT_DIR / "threshold_dilution_summary.csv"
    if not p.exists():
        print("[figures] threshold_dilution_summary.csv missing; figureS2 not drawn")
        return
    d = pd.read_csv(p)
    d = d[d.tier == "strong"]
    metrics = [("frac_subsets_reached",
                "Training subsets on which\nthe Strong tier is reached"),
               ("frac_pairs_cleared",
                "Held-out genes clearing the\nStrong cut, where evaluable")]
    fig, axes = plt.subplots(2, 3, figsize=(DOUBLE, 105 * MM), sharey=True,
                             gridspec_kw={"hspace": 0.35, "wspace": 0.08})
    for i, (col, ylab) in enumerate(metrics):
        for j, st in enumerate(STRATA):
            ax = axes[i, j]
            for tool, colour, marker in DILUTION_TOOLS:
                t = d[(d.tool == tool) & (d.stratum == st)].sort_values("k")
                t = t[np.isfinite(t[col])]
                if not len(t):
                    continue
                ax.plot(t.k, t[col], color=colour, marker=marker, ms=3.5, lw=1.2,
                        mec="white", mew=0.4)
            ax.set_ylim(-0.04, 1.04)
            ax.set_xticks(range(2, 7))
            ax.set_xlim(1.7, 6.3)
            if i == 0:
                ax.set_title(STRATUM_LABEL[st], fontweight="bold")
            if i == 1:
                ax.set_xlabel("Training genes")
            if j == 0:
                ax.set_ylabel(ylab)
            ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0))
    handles = [Line2D([], [], color=c, marker=m, ms=4, lw=1.2, mec="white",
                      label=NAME[t]) for t, c, m in DILUTION_TOOLS]
    fig.legend(handles=handles, loc="lower center", ncol=4, bbox_to_anchor=(0.55, 0.0))
    fig.subplots_adjust(left=0.13, right=0.98, top=0.93, bottom=0.2)
    _save(fig, "figureS2")


def main() -> None:
    _style()
    figure1()
    figure2()
    figure3()
    figure4()
    figure_s1()
    figure_s2()
    figure_s3()


if __name__ == "__main__":
    main()
