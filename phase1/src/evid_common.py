"""Shared pieces for the evidence-strength stages (E2-E7).

Kept separate so the panel, the orientation convention and the fusion protocol are
defined once. Nothing here fits on data it is then evaluated on.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from . import config as C

SET_PATH = Path("data/evidence/analysis_set_v1.parquet")
REPORT_DIR = Path("reports/evidence")

# The E1 panel. AlphaMissense (no splice coverage) and the gnomAD allele-frequency
# columns (BA1/BS1/PM2 evidence, not PP3) are out by design.
PANEL = ["spliceai", "pangolin", "alphagenome", "cadd", "phylop", "phastcons",
         "gpn_msa", "nt"]
# scored separately; present only once their stage has run
OPTIONAL = ["spliceai_walker", "avi", "avi_splice"]

STRATA = ["pm12", "s3_10", "s11_50", "all_1_50"]
ARMS = ["all", "classified", "recorded_unclassified", "unrecorded"]

FUSION = "fusion_enet"

MIN_POS = MIN_NEG = 10


def load_set(path: Path | str = SET_PATH) -> pd.DataFrame:
    return pd.read_parquet(path)


def panel_of(df: pd.DataFrame) -> list[str]:
    """Panel columns actually present, in a fixed order."""
    return [c for c in PANEL + OPTIONAL if c in df.columns]


def stratum_frame(df: pd.DataFrame, stratum: str) -> pd.DataFrame:
    return df if stratum == "all_1_50" else df[df["stratum"] == stratum]


def arm_frame(df: pd.DataFrame, arm: str) -> pd.DataFrame:
    return df if arm == "all" else df[df["clinvar_arm"] == arm]


def rank_within_gene(df: pd.DataFrame, feats: list[str]) -> pd.DataFrame:
    """Within-gene rank to the unit interval. NaN preserved; self-contained per
    gene, so it is LOGO-safe.

    No feature is flipped here. config.REVERSED_FEATURES exists because the
    published frozen matrix stored GPN-MSA and allele frequency on an
    anti-correlated convention; the atlas columns this set is built from are
    already oriented, and E1 measures and records the sign of every panel column
    (reports/evidence/feature_orientation.csv) rather than assuming either.
    """
    out = df[feats].copy()
    for _, idx in df.groupby("gene", observed=True).groups.items():
        for f in feats:
            v = df.loc[idx, f].to_numpy(dtype=float)
            m = ~np.isnan(v)
            r = np.full_like(v, np.nan)
            if m.sum() > 1:
                r[m] = (rankdata(v[m]) - 0.5) / m.sum()
            elif m.sum() == 1:
                r[m] = 0.5
            out.loc[idx, f] = r
    return out


def logo_fusion(df: pd.DataFrame, feats: list[str]) -> np.ndarray:
    """Out-of-fold elastic-net fusion, leave-one-gene-out, on the within-gene rank
    of functional pathogenicity -- the published study's protocol, on the new panel.

    Returned as an out-of-fold score aligned to df's rows. Genes are never a
    feature and every in-fold statistic is fit on training genes only.
    """
    from .phase2_model import logo_oof, rank_target_within_gene
    xrn = rank_within_gene(df, feats)
    y = rank_target_within_gene(df, "func_pathogenicity")
    genes = df["gene"].astype(str).reset_index(drop=True)
    xrn = xrn.reset_index(drop=True)
    keep = y.reset_index(drop=True).notna().to_numpy()
    oof = np.full(len(df), np.nan)
    sub = logo_oof(xrn.loc[keep].reset_index(drop=True),
                   y.reset_index(drop=True).loc[keep].reset_index(drop=True),
                   genes.loc[keep].reset_index(drop=True), feats, "enet")
    oof[keep] = sub
    return oof


def band_lr(y: np.ndarray, s: np.ndarray, thr: float, side: str = "upper") -> float:
    """P(band | damaging) / P(band | normal) for the band a threshold cuts off.

    side='upper' is the PP3 band (s >= thr), side='lower' the BP4 band (s <= thr).
    A zero denominator is bounded away from zero rather than returned as infinity,
    the convention phase5 already uses.
    """
    pos, neg = s[y == 1], s[y == 0]
    if len(pos) < MIN_POS or len(neg) < MIN_NEG or not np.isfinite(thr):
        return np.nan
    if side == "upper":
        p, n = float((pos >= thr).mean()), float((neg >= thr).mean())
    else:
        p, n = float((pos <= thr).mean()), float((neg <= thr).mean())
    n = max(n, 1.0 / (len(neg) + 1))
    p = max(p, 1.0 / (len(pos) + 1)) if p == 0 else p
    return p / n


def acmg_bands(cfg: dict) -> tuple[dict, dict]:
    return (cfg["evidence_bands"]["pathogenic_lr_thresholds"],
            cfg["evidence_bands"]["benign_lr_thresholds"])


def tier_of(lr: float, bands: dict, direction: str) -> str:
    if lr is None or not np.isfinite(lr):
        return "not evaluable"
    order = [("Very strong", bands["very_strong"]), ("Strong", bands["strong"]),
             ("Moderate", bands["moderate"]), ("Supporting", bands["supporting"])]
    for name, t in order:
        if (lr >= t) if direction == "pathogenic" else (lr <= t):
            return name
    return "below supporting"
