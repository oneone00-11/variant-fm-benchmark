"""E2.9 -- one tidy table per figure, plus a rough PNG to look at the data.

The tables are the deliverable. The PNGs are matplotlib defaults and exist only so
the shape of each figure can be checked before anyone spends time on styling; they
are not the figures.

Every value is read from a pipeline output. Nothing is recomputed here, so a figure
and the table behind it cannot drift apart.

Run (PYTHONPATH=phase1):  python -m src.evid_fig_data
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                        # noqa: E402
import pandas as pd                       # noqa: E402
import yaml                               # noqa: E402

from . import evid_common as K            # noqa: E402

REPORT_DIR = Path("reports/evidence")
FIG_DIR = REPORT_DIR / "fig_data"
IN_SCOPE = ["s3_10", "s11_50", "s3_50"]


# ---------------------------------------------------------------- fig 1
def fig1() -> pd.DataFrame:
    """Walker's fixed cut points: LR+ and LR(<=0.1) by stratum and by ClinVar arm,
    with the no-BRCA1 control and the two external genes as further rows."""
    rows = []
    w = pd.read_csv(REPORT_DIR / "walker_thresholds.csv")
    col = "spliceai_walker" if "spliceai_walker" in set(w.score_column) else "spliceai"
    w = w[(w.score_column == col) & (w.scope == "pooled") & (w.status == "ok")]
    for r in w.itertuples():
        if r.stratum not in IN_SCOPE:
            continue
        rows.append({"source": "seven_genes", "gene": "seven_genes",
                     "label_definition": "assay_official", "stratum": r.stratum,
                     "arm": r.clinvar_arm, "n_pos": r.n_pos, "n_neg": r.n_neg,
                     "lr_pp3": r.lr_pp3, "lr_pp3_lo": r.lr_pp3_lo,
                     "lr_pp3_hi": r.lr_pp3_hi, "tier_pp3": r.tier_pp3,
                     "lr_le01": r.lr_bp4, "lr_le01_lo": getattr(r, "lr_bp4_lo", np.nan),
                     "lr_le01_hi": getattr(r, "lr_bp4_hi", np.nan),
                     "tier_le01": r.tier_bp4})
    arms = REPORT_DIR / "territory_metrics_arms_noBRCA1.csv"
    if arms.exists():
        a = pd.read_csv(arms)
        a = a[(a.gene_set == "no_BRCA1") & (a.tool == col) & (a.stratum.isin(IN_SCOPE))
              & a.walker_lr_pp3.notna()]
        for r in a.itertuples():
            rows.append({"source": "seven_genes", "gene": "seven_genes_noBRCA1",
                         "label_definition": "assay_official", "stratum": r.stratum,
                         "arm": f"{r.clinvar_arm}_noBRCA1",
                         "n_pos": np.nan, "n_neg": np.nan,
                         "lr_pp3": r.walker_lr_pp3, "lr_pp3_lo": np.nan,
                         "lr_pp3_hi": np.nan, "tier_pp3": r.walker_tier_pp3,
                         "lr_le01": r.walker_lr_bp4, "lr_le01_lo": np.nan,
                         "lr_le01_hi": np.nan, "tier_le01": ""})
    for gene in ("ddx3x", "tp53"):
        f = REPORT_DIR / f"external_{gene}.csv"
        if not f.exists():
            continue
        e = pd.read_csv(f)
        e = e[(e.status == "ok") & (e.tool == col) & (e.stratum.isin(IN_SCOPE))
              & e.walker_lr_pp3.notna()]
        for r in e.itertuples():
            rows.append({"source": "external", "gene": gene.upper(),
                         "label_definition": r.label_definition, "stratum": r.stratum,
                         "arm": "all", "n_pos": r.n_pos, "n_neg": r.n_neg,
                         # one gene: a variant-level bootstrap within the gene
                         # (evid_external._band_detail), not the gene-clustered
                         # interval of the seven-gene rows
                         "lr_pp3": r.walker_lr_pp3,
                         "lr_pp3_lo": getattr(r, "walker_lr_pp3_lo", np.nan),
                         "lr_pp3_hi": getattr(r, "walker_lr_pp3_hi", np.nan),
                         "tier_pp3": r.walker_tier_pp3,
                         "lr_le01": r.walker_lr_bp4,
                         "lr_le01_lo": getattr(r, "walker_lr_bp4_lo", np.nan),
                         "lr_le01_hi": getattr(r, "walker_lr_bp4_hi", np.nan),
                         "tier_le01": r.walker_tier_bp4})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- fig 2
def fig2() -> tuple[pd.DataFrame, pd.DataFrame]:
    """The score-to-evidence curves, and the tier lines as constants beside them."""
    frames = []
    for f in sorted(REPORT_DIR.glob("interval_lr_*.csv")):
        c = pd.read_csv(f)
        if c.empty or c["stratum"].iloc[0] not in IN_SCOPE:
            continue
        keep = ["tool", "stratum", "score", "local_lr", "lr_lo_gene_boot",
                "lr_hi_gene_boot", "n_damaging_in_window", "n_normal_in_window"]
        frames.append(c[[k for k in keep if k in c.columns]])
    curves = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    curves = curves.rename(columns={"local_lr": "local_lr",
                                    "lr_lo_gene_boot": "lo", "lr_hi_gene_boot": "hi"})
    cfg = yaml.safe_load(Path("config/walker2023.yaml").read_text())
    path_cuts, ben_cuts = K.acmg_bands(cfg)
    lines = pd.DataFrame(
        [{"direction": "pathogenic", "tier": k, "lr": v} for k, v in path_cuts.items()]
        + [{"direction": "benign", "tier": k, "lr": v} for k, v in ben_cuts.items()]
        + [{"direction": "walker_cutpoint", "tier": "pp3",
            "lr": cfg["thresholds"]["pp3"]["value"]},
           {"direction": "walker_cutpoint", "tier": "bp4",
            "lr": cfg["thresholds"]["bp4"]["value"]}])
    return curves, lines


# ---------------------------------------------------------------- fig 3
def fig3() -> pd.DataFrame:
    """Evidence by territory and ClinVar arm: the tier a held-out gene supports,
    beside the threshold-free metrics on the same cells."""
    tm = pd.read_csv(REPORT_DIR / "territory_metrics.csv")
    tm = tm[(tm.status == "ok") & tm.stratum.isin(IN_SCOPE)]
    tier = pd.read_csv(REPORT_DIR / "tier_logo_table.csv")
    tier = tier[tier.status == "ok"]
    order = {"supporting": 1, "moderate": 2, "strong": 3}
    tier["rank"] = tier.tier.map(order)
    # the strongest tier the stratum supports in-sample, and the strongest whose
    # threshold clears its own cut on a majority of held-out genes
    ins = (tier[tier.in_sample_tier_reached.fillna(False)]
           .sort_values("rank").groupby(["tool", "stratum"]).tail(1)
           [["tool", "stratum", "tier"]].rename(columns={"tier": "max_tier_in_sample"}))
    logo = tier.copy()
    logo["majority"] = logo.folds_heldout_lr_above_cut >= (logo.n_folds / 2.0)
    lg = (logo[logo.majority].sort_values("rank").groupby(["tool", "stratum"]).tail(1)
          [["tool", "stratum", "tier"]].rename(columns={"tier": "max_tier_logo"}))
    out = (tm[["tool", "stratum", "clinvar_arm", "auroc", "auroc_lo", "auroc_hi",
               "auprc_median", "n_pos", "n_neg", "k_genes"]]
           .rename(columns={"clinvar_arm": "arm", "auprc_median": "prauc"})
           .merge(ins, on=["tool", "stratum"], how="left")
           .merge(lg, on=["tool", "stratum"], how="left"))
    out["max_tier_in_sample"] = out["max_tier_in_sample"].fillna("none")
    out["max_tier_logo"] = out["max_tier_logo"].fillna("none")
    return out


# ---------------------------------------------------------------- fig 4
def fig4() -> pd.DataFrame:
    """External genes: the fixed cut point and the fitted threshold, side by side."""
    rows = []
    for gene in ("ddx3x", "tp53"):
        f = REPORT_DIR / f"external_{gene}.csv"
        if not f.exists():
            continue
        e = pd.read_csv(f)
        e = e[(e.status == "ok") & e.stratum.isin(IN_SCOPE)]
        for r in e.itertuples():
            for tier in ("supporting", "moderate", "strong"):
                rows.append({
                    "gene": gene.upper(), "label_definition": r.label_definition,
                    "stratum": r.stratum, "tool": r.tool, "tier": tier,
                    "n_pos": r.n_pos, "n_neg": r.n_neg,
                    "column_basis_comparable": getattr(r, "column_basis_comparable", True),
                    "lr_at_walker_cut": getattr(r, "walker_lr_pp3", np.nan),
                    "tier_at_walker_cut": getattr(r, "walker_tier_pp3", ""),
                    "logo_threshold": getattr(r, f"e3_{tier}_threshold", np.nan),
                    "lr_at_logo_threshold": getattr(r, f"e3_{tier}_lr_here", np.nan),
                    "tier_at_logo_threshold": getattr(r, f"e3_{tier}_tier_here", ""),
                    "n_pos_above_logo_threshold": getattr(r, f"e3_{tier}_n_pos_band", np.nan),
                    "n_neg_above_logo_threshold": getattr(r, f"e3_{tier}_n_neg_band", np.nan),
                    "lr_at_logo_threshold_lo": getattr(r, f"e3_{tier}_lr_lo", np.nan),
                    "lr_at_logo_threshold_hi": getattr(r, f"e3_{tier}_lr_hi", np.nan),
                    "tier_at_logo_threshold_lower_bound":
                        getattr(r, f"e3_{tier}_tier_at_bound", ""),
                    "lr_at_walker_cut_lo": getattr(r, "walker_lr_pp3_lo", np.nan),
                    "lr_at_walker_cut_hi": getattr(r, "walker_lr_pp3_hi", np.nan),
                })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- drafts
def draft_pngs(f1, f2c, f2l, f3, f4) -> None:
    if len(f1):
        d = f1[(f1.source == "seven_genes") & (f1.gene == "seven_genes")]
        fig, ax = plt.subplots(figsize=(9, 4))
        for i, (arm, g) in enumerate(d.groupby("arm")):
            g = g.set_index("stratum").reindex(IN_SCOPE)
            ax.errorbar(np.arange(len(IN_SCOPE)) + i * 0.12, g.lr_pp3,
                        yerr=[g.lr_pp3 - g.lr_pp3_lo, g.lr_pp3_hi - g.lr_pp3],
                        fmt="o", capsize=3, label=arm)
        for cut, name in [(2.08, "Supporting"), (4.33, "Moderate"), (18.7, "Strong")]:
            ax.axhline(cut, ls=":", lw=0.8, color="grey")
            ax.text(len(IN_SCOPE) - 0.4, cut, name, fontsize=7, va="bottom")
        ax.set_xticks(range(len(IN_SCOPE))); ax.set_xticklabels(IN_SCOPE)
        ax.set_yscale("log"); ax.set_ylabel("LR+ at SpliceAI >= 0.2"); ax.legend(fontsize=7)
        ax.set_title("fig1 draft: Walker cut point by stratum and ClinVar arm")
        fig.tight_layout(); fig.savefig(FIG_DIR / "fig1_draft.png", dpi=110); plt.close(fig)

    if len(f2c):
        tools = ["spliceai_walker", "pangolin", "alphagenome", "avi"]
        tools = [t for t in tools if t in set(f2c.tool)]
        fig, axes = plt.subplots(1, len(IN_SCOPE), figsize=(13, 3.6), sharey=True)
        for ax, st in zip(np.atleast_1d(axes), IN_SCOPE):
            for tool in tools:
                c = f2c[(f2c.tool == tool) & (f2c.stratum == st)].sort_values("score")
                if len(c):
                    ax.plot(c.score, c.local_lr, lw=1, label=tool)
            for cut in (2.08, 4.33, 18.7):
                ax.axhline(cut, ls=":", lw=0.8, color="grey")
            ax.set_yscale("log"); ax.set_title(st, fontsize=9); ax.set_xlabel("score")
        np.atleast_1d(axes)[0].set_ylabel("local LR"); np.atleast_1d(axes)[0].legend(fontsize=6)
        fig.suptitle("fig2 draft: local likelihood ratio against score", fontsize=10)
        fig.tight_layout(); fig.savefig(FIG_DIR / "fig2_draft.png", dpi=110); plt.close(fig)

    if len(f3):
        d = f3[f3.arm == "all"]
        fig, ax = plt.subplots(figsize=(9, 4))
        for st in IN_SCOPE:
            g = d[d.stratum == st].sort_values("auroc")
            ax.plot(g.auroc, g.tool, "o", label=st)
        ax.set_xlabel("pooled AUROC"); ax.legend(fontsize=7)
        ax.set_title("fig3 draft: threshold-free performance by territory")
        fig.tight_layout(); fig.savefig(FIG_DIR / "fig3_draft.png", dpi=110); plt.close(fig)

    if len(f4):
        d = f4[(f4.tier == "moderate") & f4.lr_at_walker_cut.notna()]
        fig, ax = plt.subplots(figsize=(9, 4))
        for i, (key, g) in enumerate(d.groupby(["gene", "label_definition"])):
            ax.plot(g.lr_at_walker_cut, g.tool + " | " + g.stratum, "o",
                    label=" / ".join(key))
        for cut, name in [(2.08, "Sup"), (4.33, "Mod"), (18.7, "Str")]:
            ax.axvline(cut, ls=":", lw=0.8, color="grey")
        ax.set_xscale("log"); ax.set_xlabel("LR+ at the fixed cut point")
        ax.legend(fontsize=6); ax.set_title("fig4 draft: external genes")
        fig.tight_layout(); fig.savefig(FIG_DIR / "fig4_draft.png", dpi=110); plt.close(fig)


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    f1 = fig1()
    f2c, f2l = fig2()
    f3 = fig3()
    f4 = fig4()
    written = []
    for name, d in (("fig1_walker_cutpoints.csv", f1),
                    ("fig2_interval_lr.csv", f2c),
                    ("fig2_tier_lines.csv", f2l),
                    ("fig3_evidence_by_territory.csv", f3),
                    ("fig4_external.csv", f4)):
        d.to_csv(FIG_DIR / name, index=False)
        written.append((name, len(d)))
    draft_pngs(f1, f2c, f2l, f3, f4)
    print(f"[E2.9] wrote {FIG_DIR}/")
    for name, n in written:
        print(f"  {name:34s} {n:>7,} rows")
    for p in sorted(FIG_DIR.glob("*.png")):
        print(f"  {p.name:34s} {p.stat().st_size // 1024:>7} KB")


if __name__ == "__main__":
    main()
