"""E19 C5 -- the held-out likelihood ratios with intervals, Table 3's fold counts by
the point estimate and by the lower bound, and the Strong thresholds across the folds.

  1. intervals (E19.7 item 1). For every column x band x tier x held-out gene of E3's
     leave-one-gene-out folds (evidence_thresholds_logo_folds.csv: the PP3 side, the
     threshold fitted on the other genes), the held-out ratio at that fold's threshold
     with a 95% interval from the resampling E7 already uses on the external genes:
     evid_external._band_detail, the held-out gene's variants resampled within each
     class (2,000 draws, class sizes fixed), seeded from the cell's identity as E7 seeds
     its rows (the gene, the label, the band, the column, then the tier). Its point
     estimate is E3's held-out ratio (evid_common.band_lr on the same variants); the
     stage stops if one differs from the per-fold file. The fusion's rows carry E3's
     leakage flag.
  2. counts (item 2). Table 3's k/n (m), once by the point estimate (a fold counts when
     its held-out ratio is at or above the tier's edge, Table 3's rule), which must
     equal the printed Table 3 cell for cell or the stage stops, and once by the lower
     bound (a fold counts when the interval's lower bound is at or above the edge; m
     then counts the folds with no lower bound). Table 3's cells only: its thirteen
     columns, the three bands, Moderate and Strong, a dash where the tier is not reached
     in-sample, and the combined Atlas score without the folds of BRCA1 and RAD51C.
     Table 3 itself is not changed.
  3. jackknife (item 3). From the per-fold file only, nothing refitted: each column's
     and band's Strong thresholds over the folds (how many, the smallest, the median,
     the largest), and the largest shift from the in-sample Strong threshold with the
     fold it comes from.

The variant-level and gene-level bootstrap thresholds stand side by side in
Supplementary Tables S6 and S7 already; no file here repeats them.

Outputs, in reports/evidence/suppl_analyses/:
  c5_heldout_lr_ci.csv                per column x band x tier x held-out gene
  c5_fold_counts_point_vs_bound.csv   per Table 3 cell: both counts, and Table 3's own
  c5_jackknife_strong.csv             per column x band

Run (PYTHONPATH=phase1):  python -m src.evid_suppl_c5
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from . import evid_common as K
from . import evid_external as X
from . import evid_interval_lr as E3
from . import evid_tables as T

OUT_DIR = Path("reports/evidence/suppl_analyses")
OUT_CI = OUT_DIR / "c5_heldout_lr_ci.csv"
OUT_COUNTS = OUT_DIR / "c5_fold_counts_point_vs_bound.csv"
OUT_JACK = OUT_DIR / "c5_jackknife_strong.csv"
FOLDS = E3.REPORT_DIR / "evidence_thresholds_logo_folds.csv"
FROZEN = E3.REPORT_DIR / "evidence_thresholds.csv"
TIER_LOGO = E3.REPORT_DIR / "tier_logo_table.csv"
TABLE3 = T.OUT / "table3_evidence_tiers.csv"
BANDS = T.STRATA
TABLE3_TIERS = [("moderate", "Moderate, held out"), ("strong", "Strong, held out")]
LEAKAGE = ("fusion out-of-fold scores in the training folds were produced by models "
           "trained on the held-out gene; held-out LR is optimistic")


def _same(a: float, b: float) -> bool:
    return (np.isnan(a) and np.isnan(b)) or a == b


def intervals(df: pd.DataFrame, folds: pd.DataFrame, path_cuts: dict) -> pd.DataFrame:
    rows = []
    for r in folds.itertuples():
        sub = K.stratum_frame(df, r.stratum).dropna(subset=["y_assay", r.tool])
        te = sub[sub["gene"].astype(str) == r.heldout_gene]
        y = te["y_assay"].to_numpy(dtype=float)
        v = te[r.tool].to_numpy(dtype=float)
        thr = float(r.threshold_from_6_genes)
        d = X._band_detail(y, v, thr, "upper", path_cuts, "pathogenic",
                           (r.heldout_gene, "y_assay", r.stratum, r.tool, f"heldout_{r.tier}"))
        if not _same(d["lr"], float(r.heldout_lr)):
            raise SystemExit(f"[C5] {r.tool} {r.stratum} {r.tier} held out {r.heldout_gene}: the "
                             f"held-out ratio {d['lr']!r} is not the per-fold file's "
                             f"{r.heldout_lr!r}: stopped")
        cut = path_cuts[r.tier]
        rows.append({
            "tool": r.tool, "stratum": r.stratum, "tier": r.tier, "cut_lr": cut,
            "heldout_gene": r.heldout_gene, "threshold_from_other_genes": thr,
            "heldout_n_pos": int((y == 1).sum()), "heldout_n_neg": int((y == 0).sum()),
            "n_pos_band": d["n_pos_band"], "n_neg_band": d["n_neg_band"],
            "heldout_lr": d["lr"], "lr_lo": d["lr_lo"], "lr_hi": d["lr_hi"],
            "tier_at_point": d["tier"], "tier_at_bound": d["tier_at_bound"],
            "clears_at_point": bool(np.isfinite(d["lr"]) and d["lr"] >= cut),
            "clears_at_bound": bool(np.isfinite(d["lr_lo"]) and d["lr_lo"] >= cut),
            "logo_threshold_leakage": LEAKAGE if r.tool == K.FUSION else ""})
    return pd.DataFrame(rows)


def _cell(n: int, ev: int, k: int) -> str:
    return f"{k}/{n}" + (f" ({n - ev})" if ev < n else "")


def counts(ci: pd.DataFrame) -> pd.DataFrame:
    tl = K.read_back(TIER_LOGO)
    t3 = K.read_back(TABLE3, dtype=str, keep_default_na=False).set_index("Predictor")
    rows = []
    for c in T.COLUMNS:
        pred = T.NAME[c] + (" †" if c == "avi" else "")
        excl = K.AVI_SEEN_IN_TRAINING if c == "avi" else ()
        for st in BANDS:
            for tier, head in TABLE3_TIERS:
                printed = t3.loc[pred, f"{T.STRATUM_LABEL[st]}: {head}"]
                r = tl[(tl.tool == c) & (tl.stratum == st) & (tl.tier == tier)]
                base = {"tool": c, "predictor": pred, "stratum": st, "tier": tier,
                        "excluded_folds": ";".join(excl), "table3_cell": printed}
                if not len(r) or r.status.iloc[0] != "ok" or not bool(r.in_sample_tier_reached.iloc[0]):
                    point = bound = "–"
                    rows.append(base | {"in_sample_tier_reached": False, "point_cell": point,
                                        "bound_cell": bound})
                else:
                    f = ci[(ci.tool == c) & (ci.stratum == st) & (ci.tier == tier)
                           & ~ci.heldout_gene.isin(excl)]
                    n = len(f)
                    p_ev, p_k = int(f.heldout_lr.notna().sum()), int(f.clears_at_point.sum())
                    b_ev, b_k = int(f.lr_lo.notna().sum()), int(f.clears_at_bound.sum())
                    point, bound = _cell(n, p_ev, p_k), _cell(n, b_ev, b_k)
                    rows.append(base | {"in_sample_tier_reached": True, "n_folds": n,
                                        "point_evaluable": p_ev, "point_k": p_k,
                                        "point_m": n - p_ev, "point_cell": point,
                                        "bound_evaluable": b_ev, "bound_k": b_k,
                                        "bound_m": n - b_ev, "bound_cell": bound})
                if point != printed:
                    raise SystemExit(f"[C5] {pred} {st} {tier}: the count by the point estimate, "
                                     f"{point!r}, is not Table 3's {printed!r}: stopped (E19.7)")
    cols = ["tool", "predictor", "stratum", "tier", "in_sample_tier_reached", "excluded_folds",
            "n_folds", "point_evaluable", "point_k", "point_m", "point_cell", "table3_cell",
            "bound_evaluable", "bound_k", "bound_m", "bound_cell"]
    return pd.DataFrame(rows)[cols]


def jackknife(folds: pd.DataFrame) -> pd.DataFrame:
    ev = K.read_back(FROZEN)
    rows = []
    for (tool, st), f in folds[folds.tier == "strong"].groupby(["tool", "stratum"], sort=False):
        e = ev[(ev.tool == tool) & (ev.stratum == st) & (ev.tier == "strong")].iloc[0]
        ins = float(e.pp3_threshold_insample)
        taus = f.set_index("heldout_gene")["threshold_from_6_genes"].astype(float)
        fin = taus[np.isfinite(taus)]
        shift = (fin - ins).abs() if np.isfinite(ins) else pd.Series(dtype=float)
        rows.append({"tool": tool, "stratum": st,
                     "in_sample_reached": bool(e.pp3_threshold_reachable),
                     "in_sample_threshold": ins, "n_folds": len(taus),
                     "n_folds_reaching": len(fin),
                     "threshold_min": float(fin.min()) if len(fin) else np.nan,
                     "threshold_median": float(fin.median()) if len(fin) else np.nan,
                     "threshold_max": float(fin.max()) if len(fin) else np.nan,
                     "max_abs_shift": float(shift.max()) if len(shift) else np.nan,
                     "fold_of_max_shift": str(shift.idxmax()) if len(shift) else ""})
    return pd.DataFrame(rows)


def main() -> None:
    cfg = yaml.safe_load(E3.CONFIG_PATH.read_text())
    path_cuts, _ = K.acmg_bands(cfg)
    df = K.load_set()
    df[K.FUSION] = K.logo_fusion(df, [t for t in K.FUSION_FEATURES if t in df.columns])
    folds = K.read_back(FOLDS)
    folds = folds[(folds.side == "pp3") & folds.stratum.isin(BANDS)].reset_index(drop=True)

    ci = intervals(df, folds, path_cuts)
    print(f"[C5] {len(ci)} folds: every held-out ratio is the per-fold file's")
    cn = counts(ci)
    print("[C5] the counts by the point estimate are Table 3's, cell for cell")
    jk = jackknife(folds)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ci.to_csv(OUT_CI, index=False)
    cn.to_csv(OUT_COUNTS, index=False)
    jk.to_csv(OUT_JACK, index=False)
    diff = cn[cn.point_cell != cn.bound_cell]
    print(f"[C5] {len(diff)} of {len(cn)} Table 3 cells count differently at the lower bound"
          + "".join(f"\n      {r.predictor} {r.stratum} {r.tier}: {r.point_cell} -> {r.bound_cell}"
                    for r in diff.itertuples()))


if __name__ == "__main__":
    main()
