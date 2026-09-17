"""E8 -- what moved between the published analysis set and the evidence-strength one.

A list, not an argument. For every quantity the published study reports that still
has a counterpart on the new set, the two values sit side by side. Quantities with
no counterpart are listed as such rather than silently dropped: the evidence-strength design
drops the Brier/ECE/confidence-band axis entirely, so nothing in it is comparable.

Run (PYTHONPATH=phase1):  python -m src.evid_delta
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from . import config as C
from . import evid_common as K
from .phase2_model import dl_pool, per_gene_rho
from .phase5_likelihood_ratios import lr_at_specificity, lr_interp_at_specificity

OLD_DIR = Path("reports/phase1")
OUT_MD = Path("../docs/evidence-delta.md")
OUT_CSV = Path("reports/evidence/evid_delta.csv")

# published column name -> new panel name
OLD_TO_NEW = {
    "spliceai": "spliceai", "pangolin": "pangolin", "alphagenome": "alphagenome",
    "cadd": "cadd", "phylop": "phylop", "phastcons": "phastcons",
    "gpn_msa": "gpn_msa", "nt": "nt",
}


def _md_table(df: pd.DataFrame) -> str:
    """A markdown table without adding a dependency for one call."""
    def fmt(v):
        if isinstance(v, float):
            return "n/a" if not np.isfinite(v) else f"{v:.4f}"
        return str(v)
    head = "| " + " | ".join(df.columns) + " |"
    rule = "|" + "|".join("---" for _ in df.columns) + "|"
    body = ["| " + " | ".join(fmt(v) for v in row) + " |"
            for row in df.itertuples(index=False)]
    return "\n".join([head, rule, *body])


def _old_leaderboard() -> pd.DataFrame:
    p = OLD_DIR / "phase2_leaderboard.csv"
    return pd.read_csv(p) if p.exists() else pd.DataFrame()


def _old_lr() -> pd.DataFrame:
    p = OLD_DIR / "phase5_likelihood_ratios.csv"
    return pd.read_csv(p) if p.exists() else pd.DataFrame()


def new_rho(df: pd.DataFrame, tool: str) -> tuple[float, int]:
    sub = df.dropna(subset=[tool, "func_pathogenicity"])
    pg = per_gene_rho(sub[tool].to_numpy(float),
                      sub["func_pathogenicity"].to_numpy(float),
                      sub["gene"].astype(str).to_numpy())
    if len(pg) < 2:
        return np.nan, len(sub)
    return float(dl_pool(pg)["rho"]), len(sub)


def new_lr(df: pd.DataFrame, tool: str) -> tuple[float, float, int]:
    """(scanned LR+, interpolated LR+, n). The likelihood ratio needs a binary
    label, so it runs on the labelled subset -- a different and smaller sample than
    the rank correlation, which needs only the continuous score. The count is
    returned with the estimate so the two cannot drift apart in the table."""
    sub = df.dropna(subset=[tool, "y_assay"])
    y = sub["y_assay"].to_numpy(float)
    s = sub[tool].to_numpy(float)
    return lr_at_specificity(y, s)[0], lr_interp_at_specificity(y, s)[0], len(sub)


def build() -> None:
    df = K.load_set()
    old_lb = _old_leaderboard()
    old_lr = _old_lr()

    rows = []
    for old_name, new_name in OLD_TO_NEW.items():
        if new_name not in df.columns:
            continue
        rho_new, n_new = new_rho(df, new_name)
        scan, interp, n_lr = new_lr(df, new_name)
        old_rho = np.nan
        if len(old_lb):
            m = old_lb[old_lb["model"].astype(str) == f"single:{old_name}"]
            if len(m):
                old_rho = float(m["pooled_rho"].iloc[0])
        old_lrplus = old_lrscan = np.nan
        if len(old_lr):
            # the published primary condition: functional labels, BRCA1 included,
            # raw scores (the calibrated rows are a monotone transform of the same
            # order, so the raw row is the comparable one)
            m = old_lr[(old_lr["object"].astype(str) == f"single:{old_name}")
                       & (old_lr["condition"] == "y_assay/BRCA1_included")
                       & (old_lr["calibration"] == "raw")]
            if len(m):
                old_lrplus = float(m["lr_plus_interp"].iloc[0])
                old_lrscan = float(m["lr_plus"].iloc[0])
        rows.append({
            "quantity": "pooled Spearman rho", "object": new_name,
            "published_set_1781": old_rho, "evid_set_8853": rho_new,
            "n_new": n_new, "comparable": "same estimator, different variant set",
        })
        rows.append({
            "quantity": "LR+ at 95% specificity (interpolated)", "object": new_name,
            "published_set_1781": old_lrplus, "evid_set_8853": interp,
            "n_new": n_lr, "comparable": "same estimator, different variant set",
        })
        rows.append({
            "quantity": "LR+ at 95% specificity (scanned)", "object": new_name,
            "published_set_1781": old_lrscan, "evid_set_8853": scan,
            "n_new": n_lr, "comparable": "same estimator, different variant set",
        })

    out = pd.DataFrame(rows)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_CSV, index=False)

    lines = [
        "# What moved between the published analysis set and the evidence-strength one",
        "",
        "A list, with no interpretation. Generated by `phase1/src/evid_delta.py`.",
        "",
        "## The two sets",
        "",
        "| | published | evidence-strength |",
        "|---|---|---|",
        "| window | intron-side, \\|offset\\| <= 8 | intron-side, 1 <= \\|offset\\| <= 50 |",
        "| strata | core (<=2) / region (3-8) | pm12 / s3_10 / s11_50 |",
        "| SNVs | 1,781 | 8,853 |",
        "| ClinVar record required | yes | no |",
        "| complete predictor row required | yes | no |",
        "| panel | 10 predictors | 8 (AlphaMissense and gnomAD AF removed) |",
        "| ClinVar arms available | 2 | 3 |",
        "",
        "Of the published 1,781, 1,774 carry over; 7 are reclassified coding/UTR by the",
        "atlas region classifier (`reports/evidence/old_set_disposition.csv`).",
        "",
        "## Quantities with a counterpart",
        "",
        _md_table(out),
        "",
        "## Quantities with no counterpart",
        "",
        "The evidence-strength design drops the probability-calibration axis, so these published",
        "quantities have nothing to compare against and are not recomputed:",
        "",
        "- Brier score and Delta-Brier, and the Murphy decomposition into reliability",
        "  and resolution (`phase3_H3_headline.csv`, `phase8_lr_decomposition.csv`)",
        "- expected calibration error (`phase3_calibration_summary.csv`)",
        "- clinical yield at the 0.90/0.10 confidence band (`phase5b_evidence_yield.csv`)",
        "- the out-of-gene isotonic and Platt calibrators themselves",
        "- Delta-rho of fusion against best single tool, and the H2 evolution-axis",
        "  ablation (`phase2_H1_stratified.csv`, `phase2_H2_ablation.csv`): the panel",
        "  changed, so the fusion is not the same object",
        "",
    ]
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(lines))
    print(f"[E8] wrote {OUT_MD} and {OUT_CSV}")
    print(out.to_string(index=False))


if __name__ == "__main__":
    build()
