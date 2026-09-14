"""Main figures 1-3, drawn from the tracked report tables and nothing else.

Until frozen-matrix-v2 the three main figures were the only delivered product with
no generator in the repository: they were drawn once, pasted into the manuscript,
and could not be re-derived when the analysis set changed -- so their captions
could be updated while the images stayed at the old numbers. This module closes
that gap. Every value plotted is read from a CSV under REPORT_DIR; no statistic is
recomputed here, so a figure cannot disagree with the table it illustrates.

    Figure 1  phase2_per_gene_rho.csv + phase2_leaderboard.csv
              pooled per-gene Spearman rho with 95% CI for the elastic-net fusion
              and the three top single tools, the seven per-gene values as faint
              dots, the best single tool as a dashed line, I^2 beside each row.
    Figure 2  phase3_reliability_fusion.csv + phase3_reliability_best_single.csv
              calibrated-probability bins (mean predicted vs observed fraction
              damaging) in the primary condition, marker area by bin occupancy.
    Figure 3  phase3_H3_headline.csv (isotonic rows)
              dBrier (best single - fusion) and dYield (fusion - best) with
              gene-clustered bootstrap intervals across the four conditions; a
              marker is filled when the headline's own test (fusion_better_on)
              says the interval excludes zero.

Output (REPORT_DIR/figures/, three formats each):
    fig1_forest_ranking, fig2_reliability, fig3_calibration_advantage

Reruns are byte-identical: the backend is Agg, the font is matplotlib's bundled
DejaVu Sans, and the date metadata every format would otherwise stamp is pinned.

Usage (PYTHONPATH=phase1):  python -m src.phase9_figures
"""
from __future__ import annotations

import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import config as C

DPI = 350
# Okabe-Ito, the palette the published figures used
BLUE, ORANGE, GREEN, GOLD, GREY = "#0072B2", "#D55E00", "#009E73", "#E69F00", "#7F7F7F"
PRETTY = {"M1_enet": "Fusion (elastic net)", "single:pangolin": "Pangolin",
          "single:spliceai": "SpliceAI", "single:alphagenome": "AlphaGenome"}
FIG1_ROWS = [("M1_enet", BLUE, "D"), ("single:pangolin", ORANGE, "o"),
             ("single:spliceai", GREEN, "o"), ("single:alphagenome", GOLD, "o")]
COND = [("y_assay/BRCA1_included", "Functional · +BRCA1"),
        ("y_assay/BRCA1_excluded", "Functional · −BRCA1"),
        ("y_clinvar/BRCA1_included", "ClinVar · +BRCA1"),
        ("y_clinvar/BRCA1_excluded", "ClinVar · −BRCA1")]

RC = {"font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"], "font.size": 9,
      "axes.linewidth": 0.8, "axes.edgecolor": "#444444", "xtick.color": "#444444",
      "ytick.color": "#444444", "xtick.labelsize": 9, "ytick.labelsize": 9,
      "axes.labelcolor": "#222222", "svg.fonttype": "path", "svg.hashsalt": "variant-fm-benchmark",
      "pdf.compression": 6, "path.simplify": False}


def _half_up(x, nd=0):
    from decimal import Decimal, ROUND_HALF_UP
    q = Decimal(1).scaleb(-nd)
    return Decimal(repr(float(x))).quantize(q, rounding=ROUND_HALF_UP)


def _ci(s):
    """"[lo,hi]" -> (lo, hi); the tables print intervals as one string."""
    if not isinstance(s, str):
        return np.nan, np.nan
    m = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", s.replace("−", "-"))
    return (float(m[0]), float(m[1])) if len(m) >= 2 else (np.nan, np.nan)


def _save(fig, stem):
    d = C.REPORT_DIR / "figures"
    d.mkdir(parents=True, exist_ok=True)
    fig.savefig(d / f"{stem}.png", dpi=DPI, metadata={"Software": "phase9_figures"})
    fig.savefig(d / f"{stem}.pdf", metadata={"CreationDate": None})
    fig.savefig(d / f"{stem}.svg", metadata={"Date": None})
    plt.close(fig)
    print(f"[phase9] {d / stem}.{{png,pdf,svg}}")


# ---------------------------------------------------------------------------
def figure1(rep):
    pg = pd.read_csv(rep / "phase2_per_gene_rho.csv")
    # the leaderboard stores rho, lo, hi and I^2 at full precision; each is rounded once here
    pool = pd.read_csv(rep / "phase2_leaderboard.csv").set_index("model").rename(columns={"pooled_rho": "rho"})
    singles = pool[pool.index.str.startswith("single:")].rho.dropna().sort_values(ascending=False)
    top = [o for o, _, _ in FIG1_ROWS[1:]]
    if list(singles.index[:3]) != top:
        sys.exit(f"[phase9] the three top single tools are now {list(singles.index[:3])}, "
                 f"not {top}; update FIG1_ROWS and the Figure 1 caption together")
    dots = {o: pg[pg.object == o].rho.to_numpy() for o, _, _ in FIG1_ROWS}
    lo_all = min(v.min() for v in dots.values())
    hi_all = max(v.max() for v in dots.values())
    x0, x1 = np.floor((lo_all - 0.012) * 50) / 50, np.ceil((hi_all + 0.004) * 50) / 50
    label_x = x1 + 0.004
    fig, ax = plt.subplots(figsize=(2431 / DPI, 996 / DPI))
    best = singles.iloc[0]
    ax.axvline(best, color=GREY, ls="--", lw=1.2, zorder=1)
    for i, (obj, colour, marker) in enumerate(FIG1_ROWS):
        y = len(FIG1_ROWS) - 1 - i
        row = pool.loc[obj]
        lo, hi = row.lo, row.hi
        ax.scatter(dots[obj], np.full(dots[obj].shape, y), s=22, color=colour,
                   alpha=0.30, linewidths=0, zorder=2)
        ax.errorbar(row.rho, y, xerr=[[row.rho - lo], [hi - row.rho]],
                    fmt=marker, color=colour, markersize=9 if marker == "D" else 8,
                    elinewidth=2.4, capsize=4, capthick=2.0, zorder=3)
        head = obj == "M1_enet"
        ax.text(label_x, y, f"{_half_up(row.rho, 3)}  (I²={_half_up(row.I2)}%)", va="center", ha="left",
                fontsize=9, fontweight="bold" if head else "normal",
                color="#222222" if not head else BLUE)
    ax.set_yticks(range(len(FIG1_ROWS)))
    ax.set_yticklabels([PRETTY[o] for o, _, _ in FIG1_ROWS][::-1])
    for t, (obj, _, _) in zip(ax.get_yticklabels()[::-1], FIG1_ROWS):
        if obj == "M1_enet":
            t.set_fontweight("bold")
            t.set_color(BLUE)
    ax.set_ylim(-0.6, len(FIG1_ROWS) - 0.4)
    ax.set_xlim(x0, x1 + 0.10)
    ax.set_xticks(np.arange(x0, x1 + 0.001, 0.02))
    ax.set_xlabel("Per-gene Spearman ρ vs functional pathogenicity  (pooled, 95% CI)")
    ax.grid(axis="x", color="#DDDDDD", lw=0.7)
    ax.set_axisbelow(True)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.tick_params(axis="y", length=0)
    fig.tight_layout()
    _save(fig, "fig1_forest_ranking")


def figure2(rep):
    fus = pd.read_csv(rep / "phase3_reliability_fusion.csv")
    best = pd.read_csv(rep / "phase3_reliability_best_single.csv")
    fig, ax = plt.subplots(figsize=(1165 / DPI, 1135 / DPI))
    ax.plot([0, 1], [0, 1], ls="--", color="#AAAAAA", lw=1.6, label="perfect calibration", zorder=1)
    area = lambda n: 18 + 260 * (n / max(fus.n.max(), best.n.max()))
    ax.scatter(best.mean_pred, best.frac_pos, s=area(best.n), color=ORANGE, marker="o",
               label="Pangolin (calibrated)", zorder=2)
    ax.scatter(fus.mean_pred, fus.frac_pos, s=area(fus.n), color=BLUE, marker="D",
               label="Fusion (calibrated)", zorder=3)
    ax.set_xlabel("Mean predicted probability (calibrated)")
    ax.set_ylabel("Observed fraction damaging")
    ax.set_xlim(-0.06, 1.06)
    ax.set_ylim(-0.06, 1.06)
    ax.set_xticks(np.arange(0, 1.01, 0.2))
    ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.grid(color="#EEEEEE", lw=0.7)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], ls="--", color="#AAAAAA", lw=1.6, label="perfect calibration"),
               Line2D([], [], ls="", marker="o", ms=6, color=ORANGE, label="Pangolin (calibrated)"),
               Line2D([], [], ls="", marker="D", ms=6, color=BLUE, label="Fusion (calibrated)")]
    ax.legend(handles=handles, loc="upper left", frameon=False, fontsize=8.5,
              handletextpad=0.4, borderpad=0.2)
    fig.tight_layout()
    _save(fig, "fig2_reliability")


def figure3(rep):
    h = pd.read_csv(rep / "phase3_H3_headline.csv")
    h = h[h.calib == "isotonic"].set_index("set")
    # (column, interval column, headline flag, colour, x label, axis starts at zero)
    # Both point estimates are the headline's own observed statistics, printed at the
    # four decimals they are stored and quoted at; a marker is filled when the
    # headline's test (fusion_better_on) says the interval excludes zero.
    panels = [("dBrier", "dBrier_ci_geneclust", "Brier", BLUE, "ΔBrier  (best single − fusion)", True),
              ("dYield", "dYield_ci", "yield", GOLD, "ΔHigh-confidence fraction  (fusion − best)", False)]
    fig, axes = plt.subplots(1, 2, figsize=(2431 / DPI, 995 / DPI), layout="constrained")
    for ax, (col, cicol, flag, colour, xlabel, from_zero) in zip(axes, panels):
        ci = {cond: _ci(h.loc[cond][cicol]) for cond, _ in COND}
        lo_all = 0.0 if from_zero else min(lo for lo, _ in ci.values())
        hi_all = max(hi for _, hi in ci.values())
        span = hi_all - lo_all
        for i, (cond, _) in enumerate(COND):
            y = len(COND) - 1 - i
            v = float(h.loc[cond, col])
            lo, hi = ci[cond]
            excl = flag in str(h.loc[cond, "fusion_better_on"]).split(",")
            ax.errorbar(v, y, xerr=[[v - lo], [hi - v]], fmt="D", color=colour,
                        markerfacecolor=colour if excl else "white",
                        markeredgecolor=colour, markeredgewidth=1.6, markersize=9,
                        elinewidth=2.4, capsize=0, alpha=1.0 if excl else 0.55, zorder=3)
            ax.text(hi + 0.04 * span, y, f"{v:+.4f}".replace("-", "−"), va="center", ha="left",
                    fontsize=9, fontweight="bold" if excl else "normal", color="#222222")
        if not from_zero:
            ax.axvline(0, color=GREY, lw=1.0, zorder=1)
        ax.set_yticks(range(len(COND)))
        # the row labels belong to the left panel only: repeating them on the right
        # steals the width its x label needs
        ax.set_yticklabels([lab for _, lab in COND][::-1] if from_zero else [])
        ax.set_ylim(-0.6, len(COND) - 0.4)
        ax.set_xlim(0.0 if from_zero else lo_all - 0.08 * span, hi_all + 0.34 * span)
        ax.set_xlabel(xlabel, fontsize=8.5)
        ax.grid(axis="x", color="#DDDDDD", lw=0.7)
        ax.set_axisbelow(True)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.tick_params(axis="y", length=0)
    fig.get_layout_engine().set(w_pad=0.06, h_pad=0.02)
    _save(fig, "fig3_calibration_advantage")


def run():
    rep = C.REPORT_DIR
    need = ["phase2_per_gene_rho.csv", "phase2_leaderboard.csv", "phase3_reliability_fusion.csv",
            "phase3_reliability_best_single.csv", "phase3_H3_headline.csv"]
    missing = [n for n in need if not (rep / n).exists()]
    if missing:
        sys.exit(f"[phase9] missing inputs {missing}; run the earlier stages first")
    with plt.rc_context(RC):
        figure1(rep)
        figure2(rep)
        figure3(rep)


if __name__ == "__main__":
    run()
