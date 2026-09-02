"""Clinical yield redefined by evidence tier rather than by an arbitrary cut-off.

The manuscript defines an "actionable" call as a calibrated confidence of at
least 0.90 or at most 0.10. Those two numbers have no source: they appear in
phase3_calibration.py as `HI, LO = 0.90, 0.10` with the comment "calibrated-
confidence thresholds", carry no citation in the text, and correspond to nothing
in ACMG/AMP. They are a plausible-looking high-confidence band, chosen.

This computes the same quantity against a defined scale instead: the fraction of
variants a predictor places at ACMG Moderate evidence or above, in either
direction. The thresholds are the score values at which the likelihood ratio
reaches the Moderate band of the Tavtigian point system at a prior of 0.10
(LR+ >= 4.33 towards pathogenic, LR- <= 1/4.33 towards benign), and they are
fitted on the training genes and applied to the held-out gene, exactly as the
calibrator is, so the yield is not read off the same variants that set it.

Both definitions are reported side by side. The old one is not deleted, because
the point is to show what changes when the cut-off is given a meaning.

Usage (PYTHONPATH=phase1):  python -m src.phase5b_evidence_yield
"""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd

from . import config as C
from .phase2_model import (usable_features, rank_within_gene,
                           rank_target_within_gene, logo_oof)
from .phase3_calibration import logo_calibrate, HI, LO
from .phase5_likelihood_ratios import PRIOR, MIN_POS, MIN_NEG

MODERATE_PATH = 4.33          # Tavtigian points, prior 0.10
MODERATE_BEN = 1 / 4.33


def implied_lr(prob, prevalence):
    """Turn a calibrated probability into the likelihood ratio it implies.

    A calibrated probability is only a probability under the prevalence of the
    set it was calibrated on. Here that is 0.49 by construction, nothing like a
    clinical prior. Dividing the posterior odds by the prior odds of the
    calibration set removes it:

        LR(v) = [p / (1 - p)] / [pi / (1 - pi)]

    what remains is a property of the score's discrimination, which is the
    quantity ACMG/AMP evidence is expressed in and the quantity that transfers to
    a population with a different prevalence. This is why the 0.49 argues for
    reporting likelihood ratios rather than against reporting anything.
    """
    p = np.clip(prob, 1e-6, 1 - 1e-6)
    prior_odds = prevalence / (1.0 - prevalence)
    return (p / (1.0 - p)) / prior_odds


def evidence_yield(prob, y, genes):
    """Fraction of variants the score places at Moderate evidence or above.

    The prevalence divided out is the training genes' own positive fraction, so
    it is out-of-gene in the same sense the calibrator is.
    """
    called = np.zeros(len(prob), dtype=bool)
    seen = np.zeros(len(prob), dtype=bool)
    lrs = np.full(len(prob), np.nan)
    for g in pd.unique(genes):
        tr, te = (genes != g), (genes == g)
        ytr = y[tr]
        ytr = ytr[~np.isnan(ytr)]
        if len(ytr) < MIN_POS + MIN_NEG:
            continue
        prev = float((ytr == 1).mean())
        if not 0 < prev < 1:
            continue
        m = te & ~np.isnan(prob)
        lr = implied_lr(prob[m], prev)
        lrs[m] = lr
        seen |= m
        idx = np.where(m)[0]
        called[idx] = (lr >= MODERATE_PATH) | (lr <= MODERATE_BEN)
    n = int(seen.sum())
    return (float(called[seen].mean()) if n else np.nan), int(called.sum()), n, lrs


def run():
    p = C.OUTPUT_DIR / f"frozen_matrix_{C.FROZEN_VERSION}.parquet"
    if not p.exists():
        sys.exit(f"[phase5b] frozen matrix not found at {p}")
    df = pd.read_parquet(p)
    df = df[df["is_splice"] & df["func_pathogenicity"].notna()].reset_index(drop=True)
    feats = usable_features(df)
    genes = df["gene"].astype(str)
    Xrn = rank_within_gene(df, feats)
    ytr = (rank_target_within_gene(df)
           if getattr(C, "TRAIN_TARGET", "rank_within_gene") == "rank_within_gene"
           else df["func_pathogenicity"])
    from .phase2_model import per_gene_rho, dl_pool
    yraw = df["func_pathogenicity"].to_numpy(); gv = genes.to_numpy()
    single_pool = {f: dl_pool(per_gene_rho(Xrn[f].to_numpy(), yraw, gv))["rho"] for f in feats}
    best = max(single_pool, key=lambda k: single_pool[k] if np.isfinite(single_pool[k]) else -9)

    scores = {"fusion_M1": logo_oof(Xrn, ytr, genes, feats, "enet"),
              "mean_M0b": logo_oof(Xrn, ytr, genes, feats, "mean")}
    for f in feats:
        scores[f"single:{f}"] = Xrn[f].to_numpy()
    for k in scores:
        s = scores[k]; lo, hi = np.nanmin(s), np.nanmax(s)
        scores[k] = (s - lo) / (hi - lo) if hi > lo else s

    rows = []
    for label_col in ["y_assay", "y_clinvar"]:
        if label_col not in df.columns or df[label_col].notna().sum() < 50:
            continue
        for brca in ["included", "excluded"]:
            sub = df if brca == "included" else df[df["gene"] != "BRCA1"].reset_index(drop=True)
            idx = (df.index if brca == "included" else df[df["gene"] != "BRCA1"].index).to_numpy()
            y = sub[label_col].astype(float).to_numpy()
            g = sub["gene"].astype(str).to_numpy()
            for key, full in scores.items():
                s = full[idx]
                cal = logo_calibrate(s, y, g, method="isotonic")
                # phase3 computes yield on labelled variants only; the two
                # definitions must be compared on exactly the same set
                m = ~np.isnan(cal) & ~np.isnan(y)
                old = float(((cal[m] >= HI) | (cal[m] <= LO)).mean()) if m.sum() else np.nan
                calm = np.where(m, cal, np.nan)
                new, n_called, n_seen, _ = evidence_yield(calm, y, g)
                rows.append({"condition": f"{label_col}/BRCA1_{brca}", "object": key,
                             "n": int(m.sum()),
                             "yield_confidence_0.90_0.10": old,
                             "yield_acmg_moderate_or_above": new,
                             "n_called_moderate": n_called, "n_evaluable": n_seen,
                             "delta": (new - old) if np.isfinite(new) and np.isfinite(old) else np.nan})
    out = pd.DataFrame(rows)
    C.REPORT_DIR.mkdir(parents=True, exist_ok=True)
    out.to_csv(C.REPORT_DIR / "phase5b_evidence_yield.csv", index=False)
    print(f"[phase5b] best single = single:{best}; prior {PRIOR}; "
          f"Moderate band LR+ >= {MODERATE_PATH} / LR- <= {MODERATE_BEN:.3f}")
    for c, gdf in out.groupby("condition"):
        print(f"\n=== {c} ===")
        print(gdf.sort_values("yield_acmg_moderate_or_above", ascending=False)[
            ["object", "n", "yield_confidence_0.90_0.10",
             "yield_acmg_moderate_or_above", "delta"]].round(4).to_string(index=False))
    print(f"\n[phase5b] wrote phase5b_evidence_yield.csv to {C.REPORT_DIR}")


if __name__ == "__main__":
    run()
