"""Shared pieces for the evidence-strength stages (E2-E7).

Kept separate so the panel, the orientation convention and the fusion protocol are
defined once. Nothing here fits on data it is then evaluated on.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from . import config as C

# The companion atlas repository, read by several stages. EVID_ATLAS_REPO overrides
# it; the default is a checkout beside this one, so that a fresh clone of both
# repositories into one directory works without configuration.
ATLAS_REPO = Path(os.environ.get(
    "EVID_ATLAS_REPO",
    str(Path(__file__).resolve().parents[2].parent / "functional-standard-atlas")))

SET_PATH = Path("data/evidence/analysis_set_v1.parquet")
REPORT_DIR = Path("reports/evidence")

# The E1 panel. AlphaMissense (no splice coverage) and the gnomAD allele-frequency
# columns (BA1/BS1/PM2 evidence, not PP3) are out by design.
PANEL = ["spliceai", "pangolin", "alphagenome", "cadd", "phylop", "phastcons",
         "gpn_msa", "nt"]
# The fusion's features: the panel without the AlphaGenome splice score. The
# AlphaGenome terms of service forbid using its outputs to train another model
# (LICENSE-DATA), and the elastic net is a model trained on its inputs, with its
# fitted coefficients in the repository. An earlier version fused all eight.
FUSION_FEATURES = [t for t in PANEL if t != "alphagenome"]
# The combined Atlas score's checkpoint was selected on the BRCA1 and RAD51C assays
# used here (predictor_training_provenance.csv, row avi), so its held-out results
# are counted over the other genes only, wherever held-out genes are counted.
AVI_SEEN_IN_TRAINING = ("BRCA1", "RAD51C")
# scored separately; present only once their stage has run
# Carried beside the panel: scored separately, evaluated like any other tool, but
# not part of the fusion's feature set (the fusion is the published study's panel,
# refit on the new strata, and adding a column would make it a different object).
OPTIONAL = ["spliceai_walker", "avi", "avi_splice_sites",
            "avi_splice_site_usage", "avi_splice_junctions"]

# pm12 is OUT of the ClinGen SVI recommendation's scope: those variants go through
# the PVS1 decision tree, not PP3/BP4. So there are two pools, not one. `s3_50` is
# the in-scope pool and is the one a Walker-threshold result should be read off;
# `all_1_50` is the descriptive total and carries half its positives from pm12,
# which inflates sensitivity and moves the BP4 band.
STRATA = ["pm12", "s3_10", "s11_50", "s3_50", "all_1_50"]
IN_SCOPE_STRATA = {"s3_10", "s11_50", "s3_50"}
ARMS = ["all", "classified", "recorded_unclassified", "unrecorded"]

FUSION = "fusion_enet"

MIN_POS = MIN_NEG = 10
# a band needs occupants before it gets a ratio -- see band_lr
MIN_BAND = 10


def load_set(path: Path | str = SET_PATH) -> pd.DataFrame:
    return pd.read_parquet(path)


def panel_of(df: pd.DataFrame) -> list[str]:
    """Panel columns actually present, in a fixed order."""
    return [c for c in PANEL + OPTIONAL if c in df.columns]


def stratum_frame(df: pd.DataFrame, stratum: str) -> pd.DataFrame:
    if stratum == "all_1_50":
        return df
    if stratum == "s3_50":
        return df[df["stratum"] != "pm12"]
    return df[df["stratum"] == stratum]


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

    One consequence to carry downstream: this is a SINGLE out-of-fold column, so
    gene h's score came from a model trained on every gene except h -- including
    gene g. Anything that later holds g out and fits on the other six is therefore
    not free of g. That is fine for per-gene rank correlation, which compares each
    gene's score with its own labels, and it is NOT fine for fitting a threshold on
    six genes and reporting it on the seventh; E3 flags the affected rows.
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

    A zero DENOMINATOR is bounded away from zero, so a perfectly separating score
    gives a large finite ratio rather than infinity -- the convention phase5 uses.
    A zero NUMERATOR is not bounded, and a band with fewer than MIN_BAND occupants
    is not evaluated at all.
    Flooring the numerator would turn "no damaging variant is in this band" into a
    positive likelihood ratio computed from the class sizes alone, with no
    observation from the band in it; on the PP3 side that manufactures evidence for
    pathogenicity out of a band that contains none.
    """
    pos, neg = s[y == 1], s[y == 0]
    if len(pos) < MIN_POS or len(neg) < MIN_NEG or not np.isfinite(thr):
        return np.nan
    in_band = (s >= thr) if side == "upper" else (s <= thr)
    if int(in_band.sum()) < MIN_BAND:
        return np.nan            # too few occupants to estimate a ratio from
    if side == "upper":
        p, n = float((pos >= thr).mean()), float((neg >= thr).mean())
    else:
        p, n = float((pos <= thr).mean()), float((neg <= thr).mean())
    if p == 0.0:
        # Zero damaging variants in the band. On the PP3 side that is simply no
        # evidence. On the BP4 side the point estimate is also zero, but zero clears
        # every benign cut, so a TIER must not be read off it -- callers assigning a
        # benign tier use band_lr_benign_bound instead.
        return 0.0
    return p / max(n, 1.0 / (len(neg) + 1))


def band_lr_benign_bound(y: np.ndarray, s: np.ndarray, thr: float) -> float:
    """The value a BP4 tier should be read off: the point estimate, except where the
    band holds no damaging variant, where it is the rule-of-three upper bound."""
    pos, neg = s[y == 1], s[y == 0]
    if len(pos) < MIN_POS or len(neg) < MIN_NEG or not np.isfinite(thr):
        return np.nan
    if int((s <= thr).sum()) < MIN_BAND:
        return np.nan
    p = float((pos <= thr).mean())
    if p > 0:
        return p / max(float((neg <= thr).mean()), 1.0 / (len(neg) + 1))
    return rule_of_three_lr(len(pos), float((neg <= thr).mean()), len(neg))


def rule_of_three_lr(n_pos_total: int, frac_neg_in_band: float,
                     n_neg_total: int) -> float:
    """Upper bound on a band likelihood ratio whose damaging count is zero.

    Observing zero of n has a one-sided 95% upper bound of 3/n on the rate (Hanley
    and Lippman-Hand 1983). Dividing by the band's normal fraction gives the
    matching bound on the ratio. Without it a band holding no damaging variant has a
    ratio of exactly zero, and zero clears every benign cut -- so an empty-of-
    positives band would be published as Very strong benign evidence, which is the
    mirror of the flooring defect the PP3 side had.
    """
    fn = max(frac_neg_in_band, 1.0 / (n_neg_total + 1))
    return (3.0 / n_pos_total) / fn


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
