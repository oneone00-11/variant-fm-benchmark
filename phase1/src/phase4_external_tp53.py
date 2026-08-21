"""
Phase 4 -- external validation on TP53 (a fully held-out gene).

The single honest question: does the calibration advantage (H3, established on the
seven training genes) REPLICATE on a gene the model has never seen?

Design (leakage-proof by construction):
  1. Fit the 8-feature fusion (M1 elastic net) and the isotonic calibrator on the
     SEVEN genes ONLY -> FREEZE them. NT and Pangolin are excluded (their scoring
     could not be byte-identically reproduced for a new gene; see manuscript).
  2. Re-run the control-anchored DIRECTIONALITY GATE on TP53 before any comparison
     (nonsense vs synonymous). If TP53 is mis-oriented, ABORT -- a flipped sign
     would silently invert every downstream number. (We already know TP53 needs a
     flip; the gate confirms the loaded matrix is correctly flipped.)
  3. APPLY the frozen fusion + frozen calibrator to TP53's 192 splice variants.
  4. Compare calibrated Brier: fusion vs the best single tool, on TP53 only, with a
     VARIANT-LEVEL bootstrap (the only valid scheme for a single gene -- there is no
     within-gene cluster structure here, unlike the 7-gene analysis).

Success = directional replication (fusion Brier lower on TP53, consistent with the
7-gene result), NOT a second independent p<0.05. n=192, single gene: limited power.

Run:  python -m src.phase4_external_tp53
"""
from __future__ import annotations
import sys
import numpy as np
import pandas as pd
from sklearn.isotonic import IsotonicRegression
from . import config as C
from .phase2_model import usable_features, rank_within_gene, rank_target_within_gene, logo_oof
from .phase3_calibration import ece, brier, yield_metrics, reliability

RNG = np.random.default_rng(C.RANDOM_SEED)
N_BOOT = 2000
TP53_PATH = C.OUTPUT_DIR.parent / "external" / "tp53_splice_scored.parquet"

# 8-feature set: full predictor/evo set minus the two irreproducible tools
EXCLUDE = {"nt", "pangolin"}


def features_8(df):
    return [f for f in usable_features(df) if f not in EXCLUDE]


# ---------------------------------------------------------------------------
# directionality gate (control-anchored) -- must pass before any comparison
# ---------------------------------------------------------------------------
def directionality_ok(df, tag):
    r = df["region"].astype(str).str.lower() if "region" in df.columns else pd.Series("", index=df.index)
    # TP53 splice-only table may not carry nonsense/synonymous; fall back to consequence
    cons = df["consequence"].astype(str).str.lower() if "consequence" in df.columns else r
    lof = df.loc[cons.str.contains("nonsense|stop_gained|frameshift", na=False), "func_pathogenicity"].dropna()
    neu = df.loc[cons.str.contains("synonymous|silent", na=False), "func_pathogenicity"].dropna()
    if len(lof) < 5 or len(neu) < 5:
        print(f"[phase4] {tag}: too few controls in splice table to gate "
              f"(nonsense={len(lof)}, synonymous={len(neu)}) -- "
              f"relying on the phase-1 directionality check already run on full TP53.")
        return True
    from sklearn.metrics import roc_auc_score
    y = np.r_[np.ones(len(lof)), np.zeros(len(neu))]
    auroc = roc_auc_score(y, np.r_[lof.values, neu.values])
    print(f"[phase4] {tag}: control AUROC = {auroc:.4f} "
          f"(nonsense median {lof.median():.3f} vs synonymous {neu.median():.3f})")
    if auroc < 0.5:
        sys.exit(f"[phase4] ABORT: TP53 appears mis-oriented (AUROC {auroc:.3f} < 0.5). "
                 f"Add TP53 to FLIP_GENES / fix the parquet before comparing.")
    return True


# ---------------------------------------------------------------------------
# fit-and-freeze on the 7 genes, apply to TP53
# ---------------------------------------------------------------------------
def fit_frozen_fusion(train, feats):
    """Train M1 elastic net on ALL 7 genes with within-gene-rank target; return a
    frozen predictor callable score(X_rank_df) -> raw fusion score."""
    from sklearn.linear_model import ElasticNetCV
    Xr = rank_within_gene(train, feats)
    y = rank_target_within_gene(train)
    med = Xr[feats].median()
    def design(Xrank):
        base = Xrank[feats].fillna(med)
        isna = Xrank[feats].isna().astype(float); isna.columns = [f"{c}_isna" for c in feats]
        return pd.concat([base, isna], axis=1).to_numpy()
    en = ElasticNetCV(l1_ratio=[.1, .5, .9], n_alphas=50, cv=5,
                      random_state=C.RANDOM_SEED, max_iter=5000)
    en.fit(design(Xr), y)
    return en, design, med


def frozen_calibrator(scores, labels):
    """Fit isotonic on the 7-gene (score,label) pairs -> frozen calibrator."""
    m = ~np.isnan(scores) & ~np.isnan(labels)
    iso = IsotonicRegression(out_of_bounds="clip", y_min=0.0, y_max=1.0)
    iso.fit(scores[m], labels[m])
    return iso


def boot_brier_diff_variant(prob_f, prob_s, label):
    """Variant-level bootstrap of (Brier_single - Brier_fusion) on TP53.
    Variant-level is the ONLY valid scheme for one gene (no cluster structure)."""
    m = ~np.isnan(prob_f) & ~np.isnan(prob_s) & ~np.isnan(label)
    pf, ps, y = prob_f[m], prob_s[m], label[m]
    obs = brier(ps, y) - brier(pf, y)
    ds = []
    n = len(y)
    for _ in range(N_BOOT):
        idx = RNG.integers(0, n, n)
        ds.append(brier(ps[idx], y[idx]) - brier(pf[idx], y[idx]))
    lo, hi = np.nanpercentile(ds, [2.5, 97.5])
    return {"obs": obs, "lo": lo, "hi": hi, "fusion_better": bool(lo > 0)}


# ---------------------------------------------------------------------------
# orchestration
# ---------------------------------------------------------------------------
def run():
    frozen = pd.read_parquet(C.OUTPUT_DIR / f"frozen_matrix_{C.FROZEN_VERSION}.parquet")
    train = frozen[frozen["is_splice"] & frozen["func_pathogenicity"].notna()].reset_index(drop=True)
    if not TP53_PATH.exists():
        sys.exit(f"[phase4] TP53 scored matrix not found at {TP53_PATH}")
    tp53 = pd.read_parquet(TP53_PATH)
    tp53 = tp53[tp53["is_splice"] & tp53["func_pathogenicity"].notna()].reset_index(drop=True)

    feats = features_8(train)
    feats = [f for f in feats if f in tp53.columns]
    print(f"[phase4] 8-feature set ({len(feats)}): {feats}")
    print(f"[phase4] train (7 genes) splice n={len(train)} | TP53 splice n={len(tp53)}")

    # gate: confirm TP53 orientation (already flipped in the parquet)
    directionality_ok(tp53, "TP53")

    # best single tool by 7-gene pooled per-gene Spearman (same criterion as H1/H3)
    from .phase2_model import per_gene_rho, dl_pool
    Xr_tr = rank_within_gene(train, feats)
    yraw = train["func_pathogenicity"].to_numpy(); gv = train["gene"].astype(str).to_numpy()
    best = max(feats, key=lambda f: dl_pool(per_gene_rho(Xr_tr[f].to_numpy(), yraw, gv))["rho"] or -9)
    print(f"[phase4] best single tool (7-gene) = {best}")

    # ---- fit & freeze fusion + calibrator on 7 genes; get 7-gene OOF for calibrator fit
    en, design, med = fit_frozen_fusion(train, feats)
    oof7 = logo_oof(Xr_tr, rank_target_within_gene(train),
                    train["gene"].astype(str), feats, "enet")
    y7_assay = train["y_assay"].astype(float).to_numpy()
    # min-max the fusion score to [0,1] using 7-gene range, then freeze that mapping
    lo7, hi7 = np.nanmin(oof7), np.nanmax(oof7)
    def mm(s): return (s - lo7) / (hi7 - lo7) if hi7 > lo7 else s
    cal_fusion = frozen_calibrator(mm(oof7), y7_assay)
    # best single: freeze its calibrator on 7-gene rank scores
    s7_best = Xr_tr[best].to_numpy()
    lo7b, hi7b = np.nanmin(s7_best), np.nanmax(s7_best)
    def mmb(s): return (s - lo7b) / (hi7b - lo7b) if hi7b > lo7b else s
    cal_best = frozen_calibrator(mmb(s7_best), y7_assay)

    # ---- APPLY frozen models to TP53 (never fit on TP53) ----
    Xr_te = rank_within_gene(tp53, feats)
    fusion_raw_tp53 = en.predict(design(Xr_te))
    prob_fusion = cal_fusion.predict(np.clip(mm(fusion_raw_tp53), 0, 1))
    prob_best = cal_best.predict(np.clip(mmb(Xr_te[best].to_numpy()), 0, 1))
    rfs = tp53["func_pathogenicity"].to_numpy()   # = +func_score = RFS (syn -1 .. nonsense +1)

    def _eval(prob_f, prob_s, y, best_name):
        res = pd.DataFrame([
            {"model": "fusion_8feat", "ECE": round(ece(prob_f, y), 4),
             "Brier": round(brier(prob_f, y), 4), **yield_metrics(prob_f, y)},
            {"model": f"best_single({best_name})", "ECE": round(ece(prob_s, y), 4),
             "Brier": round(brier(prob_s, y), 4), **yield_metrics(prob_s, y)},
        ])
        return res, boot_brier_diff_variant(prob_f, prob_s, y)

    def _verdict(d):
        return ("YES — fusion Brier lower on held-out TP53 (CI excludes 0)"
                if d["fusion_better"] and d["obs"] > 0 else
                "DIRECTIONAL — sign consistent with 7 genes" if d["obs"] > 0 else
                "NOT replicated — fusion not better on TP53 (report honestly)")

    C.REPORT_DIR.mkdir(parents=True, exist_ok=True)

    # ---- MAIN analysis: label (1) control-anchored RFS>0, all 192 (approved 2026-07-10) ----
    # RFS anchored synonymous=-1 / nonsense=+1 (Funk et al. 2025); RFS>0 => damaging. No
    # official per-variant classification was released. Label used ONLY for the on-TP53
    # Brier/ECE -- never for training or calibration.
    y_main = (rfs > 0).astype(float)
    res, d = _eval(prob_fusion, prob_best, y_main, best)
    res.to_csv(C.REPORT_DIR / "phase4_tp53_external.csv", index=False)
    reliability(prob_fusion, y_main).to_csv(C.REPORT_DIR / "phase4_tp53_reliability.csv", index=False)
    print(f"\n=== TP53 EXTERNAL VALIDATION -- MAIN (label 1: control-anchored RFS>0, "
          f"n={len(y_main)}, damaging={int(y_main.sum())}) ===")
    print(res.to_string(index=False))
    print(f"\ndelta Brier (best_single - fusion) = {d['obs']:.4f}   "
          f"variant-level 95% CI [{d['lo']:.4f}, {d['hi']:.4f}]")
    print("REPLICATION:", _verdict(d))

    # ---- LABEL-DEFINITION ROBUSTNESS: does the delta-Brier survive all 3 label rules? ----
    y_med = (rfs > np.median(rfs)).astype(float)                # median split (top 50% by RFS)
    _, d_med = _eval(prob_fusion, prob_best, y_med, best)
    mask05 = np.abs(rfs) >= 0.5                                 # NA mid-band (|RFS|>=0.5)
    _, d_nab = _eval(prob_fusion[mask05], prob_best[mask05], (rfs[mask05] > 0).astype(float), best)
    print("\n=== LABEL-DEFINITION ROBUSTNESS (delta-Brier = best_single - fusion) ===")
    print(f"  {'label rule':24} {'n':>4} {'dmg':>4} {'dBrier':>8} {'variant 95% CI':>20} {'fusion_better':>13}")
    for name, dd, nn, nd in [("control-anchored RFS>0", d, len(y_main), int(y_main.sum())),
                             ("median split (top 50%)", d_med, len(y_med), int(y_med.sum())),
                             ("NA mid-band |RFS|>=0.5", d_nab, int(mask05.sum()), int((rfs[mask05] > 0).sum()))]:
        print(f"  {name:24} {nn:>4} {nd:>4} {dd['obs']:>8.4f} "
              f"[{dd['lo']:+.4f}, {dd['hi']:+.4f}] {str(dd['fusion_better']):>13}")

    # ---- SENSITIVITY: label (3) NA mid-band -- drop |RFS|<band (pLOF-like ambiguous),
    #      keep clear WT-like (RFS<=-band) vs LOF (RFS>=+band). Proves the delta-Brier
    #      conclusion does NOT depend on how mid-band variants are handled.
    print("\n=== SENSITIVITY (label 3: NA mid-band) -- delta-Brier stability vs mid-band handling ===")
    print(f"  {'band|RFS|>=':11} {'n':>4} {'dmg':>4} {'fusionBrier':>11} {'bestBrier':>10} "
          f"{'dBrier':>8} {'variant 95% CI':>20} {'fusion_better':>13}")
    primary = None
    for band in (0.25, 0.5, 0.75):
        mask = np.abs(rfs) >= band
        yb = (rfs[mask] > 0).astype(float)
        res_b, d_b = _eval(prob_fusion[mask], prob_best[mask], yb, best)
        fb = float(res_b.loc[res_b.model == "fusion_8feat", "Brier"].iloc[0])
        sb = float(res_b.loc[res_b.model.str.startswith("best"), "Brier"].iloc[0])
        print(f"  {band:<11} {int(mask.sum()):>4} {int(yb.sum()):>4} {fb:>11.4f} {sb:>10.4f} "
              f"{d_b['obs']:>8.4f} [{d_b['lo']:+.4f}, {d_b['hi']:+.4f}] {str(d_b['fusion_better']):>13}")
        if band == 0.5:
            primary = (res_b, d_b, mask, yb)
    res_b, d_b, mask, yb = primary   # band 0.5 is the reported appendix sensitivity
    res_b.to_csv(C.REPORT_DIR / "phase4_tp53_external_naband.csv", index=False)
    reliability(prob_fusion[mask], yb).to_csv(C.REPORT_DIR / "phase4_tp53_reliability_naband.csv", index=False)
    print(f"\n  appendix (band 0.5, n={int(mask.sum())}): REPLICATION: {_verdict(d_b)}")

    print("\n[phase4] main=label(1); appendix sensitivity=label(3). n small, single gene: "
          "read as directional replication, not a new significance test.")


if __name__ == "__main__":
    run()
