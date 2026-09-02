"""Cross-fold weight stability, and how few training genes the calibration needs.

Two questions left open by review.

AB-1, Reviewer 1 methodological point 3: are the fusion's weights stable across
the leave-one-gene-out folds? The pipeline reported only the cross-fold mean
|coefficient|, which cannot answer it -- a feature whose weight swings between
+0.4 and -0.4 has the same mean magnitude as one that sits at 0.4 in every fold.
This reports the weights fold by fold and their dispersion.

AB-2, the second half of Reviewer 1's main point 3: does the isotonic
calibration's advantage depend on the number of training genes? The LOGO
protocol is repeated with the training set thinned to 5, 4 and 3 genes, drawn
without the held-out gene, and the Brier score tracked as it shrinks.

Usage (PYTHONPATH=phase1):  python -m src.phase6_sensitivity
"""
from __future__ import annotations

import itertools
import sys

import numpy as np
import pandas as pd

from . import config as C
from .phase2_model import (usable_features, rank_within_gene,
                           rank_target_within_gene, logo_oof, _en_design)
from .phase3_calibration import logo_calibrate, brier

RNG = np.random.default_rng(C.RANDOM_SEED)
MAX_SUBSETS = 20      # per (held-out gene, size): all of them when fewer exist


def weight_stability(Xrn, ytr, genes, feats):
    """Per-fold elastic-net weights and their spread across the seven folds."""
    _, per_fold = logo_oof(Xrn, ytr, genes, feats, "enet", collect_coefs=True)
    per_fold = per_fold.sort_index()
    summary = pd.DataFrame({
        "mean": per_fold.mean(axis=1),
        "sd": per_fold.std(axis=1, ddof=1),
        "min": per_fold.min(axis=1),
        "max": per_fold.max(axis=1),
        "range": per_fold.max(axis=1) - per_fold.min(axis=1),
        "mean_abs": per_fold.abs().mean(axis=1),
        "n_folds_nonzero": (per_fold != 0).sum(axis=1),
        "n_folds_sign_flip": per_fold.apply(
            lambda r: int((r > 0).any() and (r < 0).any()), axis=1),
    })
    # a weight whose spread exceeds its own magnitude is not a stable weight
    summary["range_over_mean_abs"] = np.where(
        summary["mean_abs"] > 0, summary["range"] / summary["mean_abs"], np.nan)
    return per_fold, summary.sort_values("mean_abs", ascending=False)


def _fit_once(Xrn, ytr, feats, tr_mask, te_mask):
    """One elastic-net fit, predicting on both the training and the held-out
    genes. The calibrator needs the in-training predictions and the evaluation
    needs the held-out ones; fitting twice would only cost time."""
    from sklearn.linear_model import ElasticNetCV
    med = Xrn.loc[tr_mask, feats].median().fillna(0.5)
    Xtr = _en_design(Xrn, tr_mask, feats, med)
    en = ElasticNetCV(l1_ratio=[.1, .5, .9], n_alphas=50, cv=5,
                      random_state=C.RANDOM_SEED, max_iter=5000)
    en.fit(Xtr.to_numpy(), ytr[tr_mask].to_numpy())
    return (en.predict(Xtr.to_numpy()),
            en.predict(_en_design(Xrn, te_mask, feats, med).to_numpy()))


def training_gene_sensitivity(df, Xrn, ytr, genes, feats, label_col="y_assay"):
    """Thin the training set to k genes and re-run the whole LOGO chain.

    Both the fusion and its calibrator are refitted on the k training genes only,
    so the number reported is what a k-gene study would actually have obtained,
    not a k-gene model wearing a seven-gene calibrator.
    """
    y = df[label_col].astype(float).to_numpy()
    gv = genes.to_numpy()
    ug = list(pd.unique(gv))
    rows = []
    for k in range(len(ug) - 1, 2, -1):          # 6, 5, 4, 3 training genes
        for held in ug:
            pool = [g for g in ug if g != held]
            combos = list(itertools.combinations(pool, k))
            if len(combos) > MAX_SUBSETS:
                pick = RNG.choice(len(combos), MAX_SUBSETS, replace=False)
                combos = [combos[i] for i in sorted(pick)]
            for combo in combos:
                tr = np.isin(gv, combo)
                te = gv == held
                if tr.sum() < 20 or te.sum() == 0:
                    continue
                # calibrate on the same k training genes, applied to the held-out gene
                s_tr, pred = _fit_once(Xrn, ytr, feats, tr, te)
                m = ~np.isnan(y[tr])
                if m.sum() < 20:
                    continue
                from sklearn.isotonic import IsotonicRegression
                iso = IsotonicRegression(out_of_bounds="clip", y_min=0.0, y_max=1.0)
                iso.fit(s_tr[m], y[tr][m])
                prob = iso.predict(pred)
                mt = ~np.isnan(y[te])
                if mt.sum() < 10:
                    continue
                rows.append({"n_training_genes": k, "held_out_gene": held,
                             "training_genes": "+".join(combo),
                             "n_test": int(mt.sum()),
                             "brier": brier(prob[mt], y[te][mt])})
    return pd.DataFrame(rows)


def run():
    p = C.OUTPUT_DIR / f"frozen_matrix_{C.FROZEN_VERSION}.parquet"
    if not p.exists():
        sys.exit(f"[phase6] frozen matrix not found at {p}")
    df = pd.read_parquet(p)
    df = df[df["is_splice"] & df["func_pathogenicity"].notna()].reset_index(drop=True)
    feats = usable_features(df)
    genes = df["gene"].astype(str)
    Xrn = rank_within_gene(df, feats)
    ytr = (rank_target_within_gene(df)
           if getattr(C, "TRAIN_TARGET", "rank_within_gene") == "rank_within_gene"
           else df["func_pathogenicity"])

    C.REPORT_DIR.mkdir(parents=True, exist_ok=True)
    per_fold, summary = weight_stability(Xrn, ytr, genes, feats)
    per_fold.to_csv(C.REPORT_DIR / "phase6_enet_weights_per_fold.csv")
    summary.to_csv(C.REPORT_DIR / "phase6_enet_weight_stability.csv")
    print("=== AB-1: elastic-net weights, fold by fold (held-out gene as column) ===")
    print(per_fold.round(4).to_string())
    print("\n=== dispersion across the seven folds ===")
    print(summary.round(4).to_string())

    sens = training_gene_sensitivity(df, Xrn, ytr, genes, feats)
    sens.to_csv(C.REPORT_DIR / "phase6_training_gene_sensitivity.csv", index=False)
    agg = sens.groupby("n_training_genes")["brier"].agg(["count", "mean", "std", "min", "max"])
    print("\n=== AB-2: Brier of the calibrated fusion by number of training genes ===")
    print(agg.round(4).to_string())
    print(f"\n[phase6] wrote phase6_enet_weights_per_fold.csv, "
          f"phase6_enet_weight_stability.csv, phase6_training_gene_sensitivity.csv "
          f"to {C.REPORT_DIR}")


if __name__ == "__main__":
    run()
