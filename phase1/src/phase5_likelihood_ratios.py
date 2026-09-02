"""Likelihood ratios by evidence tier, before and after calibration.

Reviewer 2's third point: the analysis set's positive fraction is 0.49, so a
calibrated probability from it is not a clinical probability. ACMG/AMP evidence
is expressed as a likelihood ratio, which is a property of the score's
discrimination and does not depend on the prevalence of the set it was measured
on. This stage reports LR+ and LR- for every object H3 evaluates, in the same
four conditions, with the same gene-clustered bootstrap.

Two questions this is built to answer honestly rather than favourably:

  1. Does the isotonic calibration change the likelihood ratios at all?
     It is a monotone map, so at a threshold defined by a specificity it cannot
     change which variants are called -- except where it creates ties, which can
     make a specificity unreachable. Any difference found is a tie artefact, not
     an improvement in evidence.

  2. Does the fusion's Brier advantage survive translation into LR?
     If it does not, that is the result. The operating point is fixed in advance
     at 95% specificity, matching the companion atlas paper, and is not moved.

Usage (PYTHONPATH=phase1):  python -m src.phase5_likelihood_ratios
"""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd

from . import config as C
from .phase2_model import (usable_features, rank_within_gene,
                           rank_target_within_gene, logo_oof)
from .phase3_calibration import logo_calibrate, CAL_METHODS

RNG = np.random.default_rng(C.RANDOM_SEED)
N_BOOT = 2000
TARGET_SPEC = 0.95          # fixed in advance; matches the companion atlas paper
TARGET_SENS = 0.95          # rule-out point, where LR- is the meaningful figure
MIN_POS = MIN_NEG = 10

# Tavtigian point system at a prior of 0.10, as used by the companion paper.
# Pathogenic tiers are LR+ thresholds; benign tiers are their reciprocals, which
# is what the point system implies for BP4-direction evidence.
PRIOR = 0.10
ACMG_TIERS = [(350.0, "Very strong"), (18.7, "Strong"),
              (4.33, "Moderate"), (2.08, "Supporting")]
BENIGN_TIERS = [(1 / 350.0, "Very strong"), (1 / 18.7, "Strong"),
                (1 / 4.33, "Moderate"), (1 / 2.08, "Supporting")]


def acmg_tier(lr: float) -> str:
    """Ported from atlas.clinical_evidence so the two papers band identically."""
    if lr is None or not np.isfinite(lr):
        # no threshold reaches the target specificity -- the score is too coarsely
        # quantised to operate at that point, which is not the same as being weak
        return "not evaluable"
    for thresh, name in ACMG_TIERS:
        if lr >= thresh:
            return name
    return "below supporting"


def benign_tier(lrm: float) -> str:
    if lrm is None or not np.isfinite(lrm):
        return "not evaluable"
    for thresh, name in BENIGN_TIERS:
        if lrm <= thresh:
            return name
    return "below supporting"


def _rates(y, s, thr):
    pos, neg = s[y == 1], s[y == 0]
    tpr = float((pos >= thr).mean())
    fpr = float((neg >= thr).mean())
    return tpr, fpr, len(pos), len(neg)


def lr_at_specificity(y, s, spec=TARGET_SPEC):
    """(LR+, LR-, threshold) at the most sensitive threshold meeting `spec`.

    Thresholds are scanned over observed score values rather than taken as a
    quantile of the negatives: several of these scores are heavily quantised, and
    a quantile threshold can land on a large tie group, driving the false-positive
    rate to 1 and collapsing the ratio even for a well-separating score.
    """
    if (y == 1).sum() < MIN_POS or (y == 0).sum() < MIN_NEG:
        return np.nan, np.nan, np.nan, np.nan, np.nan
    target_fpr = 1.0 - spec
    neg = s[y == 0]
    for thr in np.unique(s):
        if float((neg >= thr).mean()) <= target_fpr:
            tpr, fpr, _, n_neg = _rates(y, s, thr)
            # a zero observed FPR is bounded away from zero rather than reported
            # as an infinite likelihood ratio
            fpr = max(fpr, 1.0 / (n_neg + 1))
            spec_obs = 1.0 - fpr
            lr_minus = (1.0 - tpr) / spec_obs if spec_obs > 0 else np.nan
            return tpr / fpr, lr_minus, float(thr), tpr, spec_obs
    # quantisation prevents reaching this specificity
    return np.nan, np.nan, np.nan, np.nan, np.nan


def lr_at_sensitivity(y, s, sens=TARGET_SENS):
    """LR- at a rule-out point. At 95% specificity LR- is close to 1 by
    construction, so it carries almost no rule-out information there; the
    meaningful figure is at a threshold chosen for sensitivity instead."""
    if (y == 1).sum() < MIN_POS or (y == 0).sum() < MIN_NEG:
        return np.nan, np.nan
    pos = s[y == 1]
    for thr in np.unique(s)[::-1]:
        if float((pos >= thr).mean()) >= sens:
            tpr, fpr, _, n_neg = _rates(y, s, thr)
            spec_obs = 1.0 - fpr
            if spec_obs <= 0:
                return np.nan, float(thr)
            tpr_c = min(tpr, 1.0 - 1.0 / (len(pos) + 1))   # bound away from LR- = 0
            return (1.0 - tpr_c) / spec_obs, float(thr)
    return np.nan, np.nan


def _boot_gene(fn, genes):
    ug = pd.unique(genes)
    if len(ug) == 0:
        # AlphaMissense scores essentially no splice variant, so its evaluable
        # set can be empty; that is reported as unevaluable, not resampled
        return []
    out = []
    for _ in range(N_BOOT):
        idx = np.concatenate([np.where(genes == g)[0]
                              for g in RNG.choice(ug, len(ug), replace=True)])
        out.append(fn(idx))
    return out


def _ci(vals):
    v = np.asarray(vals, dtype=float)
    if v.size == 0:
        return np.nan, np.nan, 0
    v = v[np.isfinite(v)]
    if v.size < N_BOOT * 0.5:      # more than half the resamples were unevaluable
        return np.nan, np.nan, int(v.size)
    lo, hi = np.nanpercentile(v, [2.5, 97.5])
    return float(lo), float(hi), int(v.size)


def evaluate(sub, scores, label_col, tag, best_key, method="isotonic"):
    y = sub[label_col].astype(float).to_numpy()
    genes = sub["gene"].astype(str).to_numpy()
    keep = ~np.isnan(y)
    rows, cal_store = [], {}
    for key, raw in scores.items():
        cal = logo_calibrate(raw, y, genes, method=method)
        cal_store[key] = cal
        for stage, s in (("raw", raw), (f"calibrated_{method}", cal)):
            m = keep & ~np.isnan(s)
            yy, ss, gg = y[m].astype(int), s[m], genes[m]
            lrp, lrm, thr, tpr, spec_obs = lr_at_specificity(yy, ss)
            lrm_ro, thr_ro = lr_at_sensitivity(yy, ss)
            lo, hi, n_ok = _ci(_boot_gene(
                lambda ix: lr_at_specificity(yy[ix], ss[ix])[0], gg))
            rows.append({
                "condition": tag, "calibration": stage, "object": key,
                "n": int(m.sum()), "n_pos": int((yy == 1).sum()),
                "n_neg": int((yy == 0).sum()),
                "threshold_at_spec95": thr,
                "sensitivity_at_that_point": tpr,
                "specificity_achieved": spec_obs,
                "lr_plus": lrp, "lr_plus_lo": lo, "lr_plus_hi": hi,
                "n_boot_evaluable": n_ok,
                "lr_minus_at_spec95": lrm,
                "lr_minus_at_sens95": lrm_ro, "threshold_at_sens95": thr_ro,
                "acmg_tier": acmg_tier(lrp),
                "benign_tier_at_sens95": benign_tier(lrm_ro),
            })
    return pd.DataFrame(rows), cal_store


def headline(sub, scores, cal_store, label_col, tag, best_key, method):
    """AA-3: does the fusion's calibration advantage survive as evidence?

    The paired difference in LR+ between the fusion and the same best single tool
    H3 uses, on the same variants, resampled by gene.
    """
    y = sub[label_col].astype(float).to_numpy()
    genes = sub["gene"].astype(str).to_numpy()
    out = []
    for stage, src in (("raw", scores), (f"calibrated_{method}", cal_store)):
        f, s = src["fusion_M1"], src[best_key]
        m = ~np.isnan(y) & ~np.isnan(f) & ~np.isnan(s)
        yy, ff, ss, gg = y[m].astype(int), f[m], s[m], genes[m]

        def d(ix):
            a = lr_at_specificity(yy[ix], ff[ix])[0]
            b = lr_at_specificity(yy[ix], ss[ix])[0]
            return a - b

        obs = d(np.arange(len(yy)))
        lo, hi, n_ok = _ci(_boot_gene(d, gg))
        lrf = lr_at_specificity(yy, ff)[0]
        lrs = lr_at_specificity(yy, ss)[0]
        out.append({"condition": tag, "calibration": stage,
                    "best_single": best_key, "n": int(m.sum()),
                    "lr_plus_fusion": lrf, "lr_plus_best_single": lrs,
                    "delta_lr_plus": obs, "lo": lo, "hi": hi,
                    "n_boot_evaluable": n_ok,
                    "fusion_better": bool(np.isfinite(lo) and lo > 0),
                    "tier_fusion": acmg_tier(lrf),
                    "tier_best_single": acmg_tier(lrs),
                    "crosses_a_tier_boundary": acmg_tier(lrf) != acmg_tier(lrs)})
    return pd.DataFrame(out)


def summarise_calibration_effect():
    """What isotonic calibration does to the likelihood ratios, as a table.

    Reads the ratios already computed rather than refitting: this is pure
    post-processing of phase5_likelihood_ratios.csv. It exists because the
    Results paragraph quotes a median change and a count of band changes, and
    a number quoted from a summary nothing stores has no source -- the failure
    this pipeline has had before.
    """
    src = C.REPORT_DIR / "phase5_likelihood_ratios.csv"
    if not src.exists():
        sys.exit(f"[phase5] {src} not found; run the main stage first")
    d = pd.read_csv(src)
    rows = []
    for method, g in d.groupby("cal_method"):
        cal = f"calibrated_{method}"
        lr = g.pivot_table(index=["condition", "object"], columns="calibration",
                           values="lr_plus").dropna()
        tier = g.pivot_table(index=["condition", "object"], columns="calibration",
                            values="acmg_tier", aggfunc="first").dropna()
        delta = (lr[cal] - lr["raw"]).abs()
        rows.append({"cal_method": method,
                     "cells_compared": int(len(lr)),
                     "cells_identical": int((delta < 0.001).sum()),
                     "median_abs_change_lr_plus": float(delta.median()),
                     "max_abs_change_lr_plus": float(delta.max()),
                     "cells_with_a_tier": int(len(tier)),
                     "cells_changing_tier": int((tier["raw"] != tier[cal]).sum())})
    out = pd.DataFrame(rows)
    out.to_csv(C.REPORT_DIR / "phase5_calibration_effect.csv", index=False)
    print("=== effect of calibration on the likelihood ratios ===")
    print(out.round(3).to_string(index=False))
    return out


def run():
    p = C.OUTPUT_DIR / f"frozen_matrix_{C.FROZEN_VERSION}.parquet"
    if not p.exists():
        sys.exit(f"[phase5] frozen matrix not found at {p}")
    df = pd.read_parquet(p)
    df = df[df["is_splice"] & df["func_pathogenicity"].notna()].reset_index(drop=True)

    # identical to phase3: same features, same target, same fusion, same best single
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
    best_key = f"single:{best}"
    for k in scores:
        s = scores[k]; lo, hi = np.nanmin(s), np.nanmax(s)
        scores[k] = (s - lo) / (hi - lo) if hi > lo else s

    C.REPORT_DIR.mkdir(parents=True, exist_ok=True)
    tables, heads = [], []
    for method in CAL_METHODS:
        for label_col in ["y_assay", "y_clinvar"]:
            if label_col not in df.columns or df[label_col].notna().sum() < 50:
                continue
            for brca in ["included", "excluded"]:
                sub = df if brca == "included" else df[df["gene"] != "BRCA1"].reset_index(drop=True)
                idx = df.index if brca == "included" else df[df["gene"] != "BRCA1"].index
                sc = {k: v[idx.to_numpy()] for k, v in scores.items()}
                tag = f"{label_col}/BRCA1_{brca}"
                tab, cal_store = evaluate(sub, sc, label_col, tag, best_key, method=method)
                tab["cal_method"] = method
                tables.append(tab)
                h = headline(sub, sc, cal_store, label_col, tag, best_key, method)
                h["cal_method"] = method
                heads.append(h)

    full = pd.concat(tables, ignore_index=True)
    head = pd.concat(heads, ignore_index=True)
    full.to_csv(C.REPORT_DIR / "phase5_likelihood_ratios.csv", index=False)
    head.to_csv(C.REPORT_DIR / "phase5_lr_headline.csv", index=False)

    print(f"[phase5] best single = {best_key}; prior {PRIOR}; "
          f"operating point = {TARGET_SPEC:.0%} specificity; {N_BOOT} gene-clustered resamples")
    print(f"[phase5] wrote phase5_likelihood_ratios.csv ({len(full)} rows) and "
          f"phase5_lr_headline.csv ({len(head)} rows) to {C.REPORT_DIR}")
    iso = head[head.cal_method == "isotonic"]
    print("\n=== AA-3: does the fusion's advantage survive as evidence? ===")
    print(iso[["condition", "calibration", "lr_plus_fusion", "lr_plus_best_single",
               "delta_lr_plus", "lo", "hi", "fusion_better",
               "crosses_a_tier_boundary"]].round(3).to_string(index=False))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "effect":
        summarise_calibration_effect()
    else:
        run()
        summarise_calibration_effect()
