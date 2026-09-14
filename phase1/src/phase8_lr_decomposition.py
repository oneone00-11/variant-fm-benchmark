"""(g) LR+ from frozen-matrix-v1 to v2: score change versus threshold placement.

Replacing the SpliceAI, Pangolin and Nucleotide Transformer columns moved two
evidence-tier assignments at the 95%-specificity operating point (Pangolin, isotonic-
calibrated: 16.4 -> 19.9 in the primary condition, 21.1 -> 17.4 without BRCA1) while
the interpolated values barely moved. This stage makes that decomposition a tracked
product for every object x condition x calibration stage:

    delta_obs        = LR+ at the scanned threshold, v2 - v1
    delta_interp     = LR+ interpolated to exactly 95% specificity, v2 - v1  (score change)
    delta_placement  = delta_obs - delta_interp                              (threshold placement)

and, within each version, placement = observed - interpolated, which is what a cell
such as AlphaGenome's 24.3 (BRCA1 excluded, calibrated) is made of: at exactly 95%
specificity the ratio is bounded by 1 / 0.05 = 20 by construction, so anything above
20 is the achieved specificity running past the nominal one.

Both matrices are rebuilt through the identical phase-5 code path (features, LOGO
fusion, min-max scaling, isotonic calibration, seed), so the v1 side reproduces the
published phase5_likelihood_ratios.csv exactly.

Usage (PYTHONPATH=phase1):  python -m src.phase8_lr_decomposition
Writes reports/phase1/phase8_lr_placement_decomposition.csv
"""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd

from . import config as C
from .phase3_calibration import logo_calibrate
from .phase5_likelihood_ratios import (lr_at_specificity, lr_interp_at_specificity,
                                       acmg_tier, TARGET_SPEC)
from .phase7_label_contrast import build_scores

V1 = C.OUTPUT_DIR / "frozen_matrix_v1.parquet"
V2 = C.OUTPUT_DIR / f"frozen_matrix_{C.FROZEN_VERSION}.parquet"


def _splice(p):
    df = pd.read_parquet(p)
    return df[df["is_splice"] & df["func_pathogenicity"].notna()].reset_index(drop=True)


def _tie(y, s, thr):
    """Tie-group sizes at the chosen threshold and at the next lower value (the
    candidate the scan rejected), with the specificity each would achieve."""
    if not np.isfinite(thr):
        return {}
    u = np.unique(s)
    i = int(np.searchsorted(u, thr))
    neg, pos = s[y == 0], s[y == 1]
    out = {}
    for name, j in (("chosen", i), ("below", i - 1)):
        if 0 <= j < len(u):
            v = u[j]
            out[f"tie_{name}_all"] = int((s == v).sum())
            out[f"tie_{name}_neg"] = int((neg == v).sum())
            out[f"tie_{name}_pos"] = int((pos == v).sum())
            out[f"spec_{name}"] = float(1 - (neg >= v).mean())
        else:
            out[f"tie_{name}_all"] = out[f"tie_{name}_neg"] = out[f"tie_{name}_pos"] = np.nan
            out[f"spec_{name}"] = np.nan
    return out


def cells(df, scores):
    rows = {}
    for label_col in ("y_assay", "y_clinvar"):
        for brca in ("included", "excluded"):
            keep_g = np.ones(len(df), bool) if brca == "included" else (df["gene"] != "BRCA1").to_numpy()
            sub = df[keep_g].reset_index(drop=True)
            sc = {k: v[keep_g] for k, v in scores.items()}
            tag = f"{label_col}/BRCA1_{brca}"
            y_all = sub[label_col].astype(float).to_numpy()
            g_all = sub["gene"].astype(str).to_numpy()
            keep = ~np.isnan(y_all)
            for key, raw in sc.items():
                cal = logo_calibrate(raw, y_all, g_all, method="isotonic")
                for stage, s in (("raw", raw), ("calibrated_isotonic", cal)):
                    m = keep & ~np.isnan(s)
                    yy, ss = y_all[m].astype(int), s[m]
                    lrp, _, thr, tpr, spec_obs = lr_at_specificity(yy, ss)
                    lri, tpr_i = lr_interp_at_specificity(yy, ss)
                    rows[(tag, stage, key)] = {
                        "n": int(m.sum()), "n_pos": int((yy == 1).sum()), "n_neg": int((yy == 0).sum()),
                        "threshold": thr, "lr_plus_obs": lrp, "spec_achieved": spec_obs,
                        "sens_obs": tpr, "lr_plus_interp": lri, "sens_interp": tpr_i,
                        "tier_obs": acmg_tier(lrp), "tier_interp": acmg_tier(lri),
                        **_tie(yy, ss, thr),
                    }
    return rows


def run():
    if C.FROZEN_VERSION == "v1":
        sys.exit("[phase8g] FROZEN_VERSION is v1; this stage compares v1 with the current "
                 "version and needs FROZEN_VERSION=v2 (the default)")
    for p in (V1, V2):
        if not p.exists():
            sys.exit(f"[phase8g] frozen matrix not found at {p}")
    a = cells(*(lambda d: (d, build_scores(d)))(_splice(V1)))
    b = cells(*(lambda d: (d, build_scores(d)))(_splice(V2)))
    out = []
    for key in a:
        ra, rb = a[key], b[key]
        tag, stage, obj = key
        d_obs = rb["lr_plus_obs"] - ra["lr_plus_obs"]
        d_int = rb["lr_plus_interp"] - ra["lr_plus_interp"]
        out.append({
            "condition": tag, "calibration": stage, "object": obj,
            "n": rb["n"], "n_pos": rb["n_pos"], "n_neg": rb["n_neg"],
            "lr_plus_obs_v1": ra["lr_plus_obs"], "spec_achieved_v1": ra["spec_achieved"],
            "lr_plus_interp_v1": ra["lr_plus_interp"],
            "tier_obs_v1": ra["tier_obs"], "tier_interp_v1": ra["tier_interp"],
            "lr_plus_obs_v2": rb["lr_plus_obs"], "spec_achieved_v2": rb["spec_achieved"],
            "lr_plus_interp_v2": rb["lr_plus_interp"],
            "tier_obs_v2": rb["tier_obs"], "tier_interp_v2": rb["tier_interp"],
            "delta_obs": d_obs, "delta_interp_score_change": d_int,
            "delta_placement": d_obs - d_int,
            "placement_within_v1": ra["lr_plus_obs"] - ra["lr_plus_interp"],
            "placement_within_v2": rb["lr_plus_obs"] - rb["lr_plus_interp"],
            "threshold_v1": ra["threshold"], "threshold_v2": rb["threshold"],
            **{f"{k}_v1": v for k, v in ra.items() if k.startswith(("tie_", "spec_chosen", "spec_below"))},
            **{f"{k}_v2": v for k, v in rb.items() if k.startswith(("tie_", "spec_chosen", "spec_below"))},
        })
    out = pd.DataFrame(out)
    C.REPORT_DIR.mkdir(parents=True, exist_ok=True)
    out.to_csv(C.REPORT_DIR / "phase8_lr_placement_decomposition.csv", index=False)
    ev = out.dropna(subset=["lr_plus_obs_v1", "lr_plus_obs_v2"])
    print(f"[phase8g] {len(out)} cells ({len(ev)} evaluable in both versions); "
          f"target specificity {TARGET_SPEC:.0%}")
    print(f"  |delta_obs| max {ev.delta_obs.abs().max():.3f}, median {ev.delta_obs.abs().median():.3f}; "
          f"|delta_interp| max {ev.delta_interp_score_change.abs().max():.3f}, "
          f"median {ev.delta_interp_score_change.abs().median():.3f}")
    print(f"  tier changes v1->v2 at the observed threshold: "
          f"{int((ev.tier_obs_v1 != ev.tier_obs_v2).sum())}; at the interpolated point: "
          f"{int((ev.tier_interp_v1 != ev.tier_interp_v2).sum())}")
    big = ev.reindex(ev.placement_within_v2.abs().sort_values(ascending=False).index).head(6)
    print("  largest within-v2 placement effects (observed - interpolated):")
    print(big[["condition", "calibration", "object", "lr_plus_obs_v2", "spec_achieved_v2",
               "lr_plus_interp_v2", "placement_within_v2"]].round(3).to_string(index=False))


if __name__ == "__main__":
    run()
