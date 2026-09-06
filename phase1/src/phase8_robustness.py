"""Phase 8 -- reviewer-requested robustness analyses (revision round 1).

Six analyses, each probing a different threat to the two headline claims
(H1: the fusion ranks better than the best single tool; H3: it is better
calibrated). Everything reuses the earlier stages' own machinery -- the same
frozen matrix, the same LOGO out-of-fold predictions, the same out-of-gene
calibration, the same LR+ threshold scan -- so a number here is directly
comparable to the number it is a robustness check on.

  (a) Exact sign-flip randomisation test. The gene-cluster bootstrap is
      asymptotic in the number of CLUSTERS, and there are only 7. Enumerating
      all 2^7 = 128 sign assignments of the per-gene effect gives the exact
      two-sided p under the sharp null; the wild cluster bootstrap (Rademacher
      weights, 2000 draws) is the same idea without enumeration. Done for
      delta-rho and for delta-Brier in all 8 calibrator x label x +/-BRCA1
      conditions (the isotonic / y_assay / +BRCA1 primary condition included).
  (b) Leave-two-genes-out: delta-rho over all 21 gene pairs, to show the H1
      gain is not carried by any single pair of genes.
  (c) Hartung-Knapp CIs for the whole Table-1 leaderboard. The DL normal
      interval is anti-conservative at k=7; dl_pool(method="hk") reports the
      calibrated t(k-1) interval alongside.
  (d) Murphy decomposition of the Brier score (reliability / resolution /
      uncertainty), fusion vs best single, under 3 binnings (10 equal-width,
      10 equal-frequency, 15 equal-width) x 4 conditions (isotonic x
      y_assay/y_clinvar x +/-BRCA1), with gene-clustered bootstrap CIs on the
      component deltas -- so the H3 Brier gain is attributed to refinement
      (RES) rather than relabelled to recalibration (REL) by assertion.
  (e) Stratified selection test (phase7's decisive test for ClinVar-record
      selection), re-run within +/-BRCA1, within the two offset strata
      (splice-core |off|<=2 vs splice-region 3-8), and with the assay-only arm
      split by record type (VUS vs Conflicting/Other).
  (f) LR+ at four fixed operating points (specificity 0.90 / 0.95 / 0.975 /
      0.99) for every evaluable object (9 single tools + M1 + M0b), on the
      primary label (y_assay), with and without BRCA1. The threshold scan over
      observed score values is unchanged from phase5; only the target
      specificity varies.

Bootstrap streams: each analysis that resamples draws from its OWN generator
seeded with config.RANDOM_SEED, in a fixed documented order, so the output is
deterministic and independent of any other stage's stream.

Run (PYTHONPATH=phase1):  python -m src.phase8_robustness
"""
from __future__ import annotations

import itertools
import sys

import numpy as np
import pandas as pd

from . import config as C
from .phase2_model import (usable_features, rank_within_gene,
                           rank_target_within_gene, logo_oof,
                           per_gene_rho, dl_pool)
from .phase3_calibration import logo_calibrate, brier, CAL_METHODS
from .phase5_likelihood_ratios import lr_at_specificity, acmg_tier

N_BOOT = 2000
# (f) fixed operating points; 0.95 is the phase5 primary, kept for continuity
SPEC_POINTS = (0.90, 0.95, 0.975, 0.99)
# (d) binning schemes: (name, n_bins, strategy)
MURPHY_BINNINGS = [("width10", 10, "width"),
                   ("freq10", 10, "frequency"),
                   ("width15", 15, "width")]
EPS = 1e-15   # tolerance for |T_perm| >= |T_obs| comparisons (float fuzz)


# ---------------------------------------------------------------------------
# shared per-gene effect decompositions
# ---------------------------------------------------------------------------
def _dl_pool_z(d, v):
    """DL pooling of per-gene Fisher-z DIFFERENCES. tau^2 is recomputed on
    every call, which is what makes this valid inside sign-flip assignments
    and wild-weight resamples (the heterogeneity changes with the signs)."""
    d = np.asarray(d, dtype=float)
    w = 1.0 / v
    zfe = (w * d).sum() / w.sum()
    Q = (w * (d - zfe) ** 2).sum()
    k = len(d)
    Cc = w.sum() - (w ** 2).sum() / w.sum()
    tau2 = max(0.0, (Q - (k - 1)) / Cc) if Cc > 0 else 0.0
    ws = 1.0 / (v + tau2)
    return float((ws * d).sum() / ws.sum())


def _n_weighted_T(d, n):
    """Size-weighted mean of per-gene effects -- the sign-flip statistic for
    delta-Brier, matching the pooled Brier difference on the full set."""
    return float(np.sum(d) / np.sum(n))


def sign_flip_exact(d, v, T_fn):
    """Exact two-sided randomisation p over all 2^k sign assignments of the
    per-gene effects, under the sharp null of symmetrically distributed
    per-gene effects. k clusters -> 2^k assignments; the smallest achievable
    p is 2/2^k (the all-plus and all-minus assignments)."""
    d = np.asarray(d, dtype=float)
    k = len(d)
    t_obs = T_fn(d, v)
    n_assign = 2 ** k
    t_perm = np.empty(n_assign)
    for i, signs in enumerate(itertools.product((1.0, -1.0), repeat=k)):
        t_perm[i] = T_fn(np.array(signs) * d, v)
    return {"T_obs": t_obs, "k": k, "n_assignments": n_assign,
            "p_exact": float((np.abs(t_perm) >= abs(t_obs) - EPS).mean())}


def wild_cluster_p(d, v, T_fn, rng, n_boot=N_BOOT):
    """Wild cluster bootstrap two-sided p (Rademacher +/-1 weights on the
    per-gene effects) -- the non-enumerated counterpart of the sign-flip
    test, reported as its asymptotic cross-check."""
    d = np.asarray(d, dtype=float)
    t_obs = T_fn(d, v)
    cnt = 0
    for _ in range(n_boot):
        wgt = rng.choice((1.0, -1.0), size=len(d))
        if abs(T_fn(wgt * d, v)) >= abs(t_obs) - EPS:
            cnt += 1
    return float(cnt / n_boot)


def _rho_diff_inputs(fusion_pg, single_pg):
    """Per-gene Fisher-z delta-rho and its variance, on the genes both
    per-gene tables share (sorted, so assignment order is stable)."""
    f = fusion_pg.set_index("gene")
    s = single_pg.set_index("gene")
    common = sorted(set(f.index) & set(s.index))
    zf = np.arctanh(f.loc[common, "rho"].clip(-0.999, 0.999).to_numpy())
    zs = np.arctanh(s.loc[common, "rho"].clip(-0.999, 0.999).to_numpy())
    nf = f.loc[common, "n"].to_numpy()
    ns = s.loc[common, "n"].to_numpy()
    return common, zf - zs, 1.0 / (nf - 3) + 1.0 / (ns - 3)


def _brier_diff_inputs(pf, ps, label, genes):
    """Per-gene (Brier_single - Brier_fusion) * n_g, on variants where both
    calibrated probabilities and the label exist. Positive = fusion better."""
    m = ~np.isnan(pf) & ~np.isnan(ps) & ~np.isnan(label)
    ug = sorted(pd.unique(genes[m]))
    d, n = [], []
    for g in ug:
        gm = m & (genes == g)
        d.append((brier(ps[gm], label[gm]) - brier(pf[gm], label[gm])) * gm.sum())
        n.append(int(gm.sum()))
    return ug, np.asarray(d), np.asarray(n)


# ---------------------------------------------------------------------------
# (a) sign-flip exact tests + wild cluster bootstrap
# ---------------------------------------------------------------------------
def sign_flip_rho(fusion_pg, single_pg, wild_rng):
    """Delta-rho rows: all genes, and BRCA1 excluded (6 clusters -> 2^6)."""
    rows = []
    for brca in ("included", "excluded"):
        fpg, spg = fusion_pg, single_pg
        if brca == "excluded":
            fpg = fpg[fpg.gene != "BRCA1"]
            spg = spg[spg.gene != "BRCA1"]
        common, d, v = _rho_diff_inputs(fpg, spg)
        r = sign_flip_exact(d, v, _dl_pool_z)
        rows.append({"analysis": "delta_rho", "condition": f"BRCA1_{brca}",
                     "statistic": "DL-pooled Fisher-z delta", **r,
                     "p_wild_rademacher": wild_cluster_p(d, v, _dl_pool_z, wild_rng),
                     "delta_rho_obs": dl_pool(fpg)["rho"] - dl_pool(spg)["rho"]})
    return rows


def sign_flip_brier(cal, wild_rng):
    """Delta-Brier rows for every calibrated condition. `cal` maps
    (method, label_col, brca) -> (p_fusion, p_single, label, genes)."""
    rows = []
    for (method, label_col, brca), (pf, ps, label, genes) in cal.items():
        ug, d, n = _brier_diff_inputs(pf, ps, label, genes)
        T_fn = _n_weighted_T
        r = sign_flip_exact(d, n, T_fn)
        rows.append({"analysis": "delta_brier",
                     "condition": f"{method}/{label_col}/BRCA1_{brca}",
                     "statistic": "n-weighted per-gene Brier delta", **r,
                     "p_wild_rademacher": wild_cluster_p(d, n, T_fn, wild_rng),
                     "delta_brier_obs": float(T_fn(d, n))})
    return rows


# ---------------------------------------------------------------------------
# (b) leave-two-genes-out
# ---------------------------------------------------------------------------
def leave_two_genes_out(fusion_pg, single_pg):
    common = sorted(set(fusion_pg.gene) & set(single_pg.gene))
    rows = []
    for drop in itertools.combinations(common, 2):
        keep = [g for g in common if g not in drop]
        ff = fusion_pg[fusion_pg.gene.isin(keep)]
        ss = single_pg[single_pg.gene.isin(keep)]
        rows.append({"dropped": "+".join(drop), "k_remaining": len(keep),
                     "fusion_rho": dl_pool(ff)["rho"],
                     "single_rho": dl_pool(ss)["rho"],
                     "delta_rho": dl_pool(ff)["rho"] - dl_pool(ss)["rho"]})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# (c) Hartung-Knapp leaderboard
# ---------------------------------------------------------------------------
def leaderboard_hk(board):
    """DL and HK intervals side by side for every Table-1 row. `board` is a
    list of (model_name, per_gene_rho_df); unpoolable rows keep NaN intervals."""
    rows = []
    for name, pg in board:
        dl = dl_pool(pg, method="dl")
        hk = dl_pool(pg, method="hk")
        rows.append({"model": name, "pooled_rho": dl["rho"],
                     "dl_lo": dl["lo"], "dl_hi": dl["hi"],
                     "hk_lo": hk["lo"], "hk_hi": hk["hi"],
                     "I2": dl["I2"], "k": dl["k"]})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# (d) Murphy decomposition of the Brier score
# ---------------------------------------------------------------------------
def murphy_parts(prob, label, n_bins=10, binning="width"):
    """Brier = REL - RES + UNC (Murphy 1973), binned on the probability axis.
    binning="width": equal-width bins on [0,1]; "frequency": equal-frequency
    bins (quantile edges of the predictions, de-duplicated)."""
    m = ~np.isnan(prob) & ~np.isnan(label)
    p, y = prob[m], label[m]
    ybar = y.mean()
    unc = ybar * (1.0 - ybar)
    if binning == "width":
        edges = np.linspace(0, 1, n_bins + 1)
    elif binning == "frequency":
        edges = np.unique(np.quantile(p, np.linspace(0, 1, n_bins + 1)))
    else:
        raise ValueError(f"unknown binning {binning!r}")
    nb = len(edges) - 1
    idx = np.clip(np.digitize(p, edges) - 1, 0, nb - 1)
    rel = res = 0.0
    for b in range(nb):
        mb = idx == b
        if mb.sum():
            w = mb.sum() / len(p)
            rel += w * (p[mb].mean() - y[mb].mean()) ** 2
            res += w * (y[mb].mean() - ybar) ** 2
    return rel, res, unc


def cluster_boot_diff(prob_f, prob_s, label, genes, fn, lower_is_better, rng,
                      n_boot=N_BOOT):
    """The phase3 gene-clustered percentile bootstrap with an explicit
    generator, so phase8's stream is independent of other stages'. Same
    resampling scheme as phase3_calibration.boot_diff: whole genes drawn with
    replacement, the metric recomputed on the concatenated resample."""
    ug = pd.unique(genes)

    def diff(idx):
        a, b = fn(prob_s[idx], label[idx]), fn(prob_f[idx], label[idx])
        return (a - b) if lower_is_better else (b - a)  # positive => fusion better

    obs = diff(np.arange(len(label)))
    ds = [diff(np.concatenate([np.where(genes == g)[0]
                               for g in rng.choice(ug, len(ug), replace=True)]))
          for _ in range(n_boot)]
    lo, hi = np.nanpercentile(ds, [2.5, 97.5])
    return {"obs": float(obs), "lo": float(lo), "hi": float(hi)}


def murphy_table(cal, boot_rng):
    """(d) over the 4 isotonic conditions x 3 binnings. `cal` is the same
    cache the sign-flip test filled; only isotonic entries are used.

    Bootstrap draw order is FIXED and documented: conditions iterate
    (y_assay/included, y_assay/excluded, y_clinvar/included,
    y_clinvar/excluded), binnings iterate MURPHY_BINNINGS, and within a cell
    the draws are Brier, REL, RES. The primary cell is therefore the first
    consumer of the stream.
    """
    rows = []
    for label_col in ("y_assay", "y_clinvar"):
        for brca in ("included", "excluded"):
            pf, ps, label, genes = cal[("isotonic", label_col, brca)]
            cond = f"isotonic/{label_col}/BRCA1_{brca}"
            for bname, nb, strat in MURPHY_BINNINGS:
                rel_f, res_f, unc = murphy_parts(pf, label, nb, strat)
                rel_s, res_s, _ = murphy_parts(ps, label, nb, strat)
                d_bri = cluster_boot_diff(pf, ps, label, genes, brier, True, boot_rng)
                d_rel = cluster_boot_diff(pf, ps, label, genes,
                                          lambda p, y: murphy_parts(p, y, nb, strat)[0],
                                          True, boot_rng)
                d_res = cluster_boot_diff(pf, ps, label, genes,
                                          lambda p, y: murphy_parts(p, y, nb, strat)[1],
                                          False, boot_rng)
                rows.append({
                    "condition": cond, "binning": bname,
                    "brier_fusion": brier(pf, label), "brier_single": brier(ps, label),
                    "rel_fusion": rel_f, "rel_single": rel_s,
                    "res_fusion": res_f, "res_single": res_s, "unc": unc,
                    "d_brier": d_bri["obs"], "d_brier_lo": d_bri["lo"], "d_brier_hi": d_bri["hi"],
                    # REL: lower is better -> single - fusion; >0 favours fusion
                    "d_rel": d_rel["obs"], "d_rel_lo": d_rel["lo"], "d_rel_hi": d_rel["hi"],
                    # RES: higher is better -> fusion - single; >0 = fusion sharper
                    "d_res": d_res["obs"], "d_res_lo": d_res["lo"], "d_res_hi": d_res["hi"],
                })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# (e) stratified selection test
# ---------------------------------------------------------------------------
def stratified_selection_test(df, scores):
    """Phase 7's selection test, re-run inside strata. Both arms are scored
    against the SAME functional labels; only set membership varies. Strata:
    overall, -BRCA1, the two offset bins, and the assay-only arm split by the
    variant's ClinVar record type (VUS vs Conflicting/Other)."""
    y = df["y_assay"].astype(float).to_numpy()
    rec = (df["y_assay"].notna() & df["y_clinvar"].notna()).to_numpy()
    only = (df["y_assay"].notna() & df["y_clinvar"].isna()).to_numpy()
    clv = df["clinvar"].astype(str)
    strata = [("overall", rec, only)]
    nb = (df["gene"] != "BRCA1").to_numpy()
    strata.append(("minus_BRCA1", rec & nb, only & nb))
    for b in ("core_like", "region"):
        mb = (df["splice_bin"] == b).to_numpy()
        strata.append((f"offset_{b}", rec & mb, only & mb))
    for cls, msk in (("VUS", clv.eq("VUS").to_numpy()),
                     ("Conflicting/Other", clv.isin(["Conflicting", "Other"]).to_numpy())):
        strata.append((f"assay_only_{cls}", rec, only & msk))

    rows = []
    for tag, mrec, monly in strata:
        for subset, mask in (("ClinVar-recorded", mrec), ("assay-only", monly)):
            for key, s in scores.items():
                name = key.replace("single:", "")
                m = mask & ~np.isnan(y) & ~np.isnan(s)
                lrp, _, thr, tpr, spec = lr_at_specificity(y[m].astype(int), s[m])
                rows.append({"stratum": tag, "subset": subset, "object": name,
                             "n": int(m.sum()), "n_pos": int((y[m] == 1).sum()),
                             "n_neg": int((y[m] == 0).sum()),
                             "lr_plus": lrp, "sensitivity": tpr,
                             "specificity": spec, "threshold": thr,
                             "tier": acmg_tier(lrp)})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# (f) LR+ at multiple fixed operating points
# ---------------------------------------------------------------------------
def lr_plus_operating_points(df, scores, label_col="y_assay"):
    """LR+ / LR- at SPEC_POINTS for every object, same observed-value
    threshold scan as phase5 (lr_at_specificity), only the target varies."""
    y = df[label_col].astype(float).to_numpy()
    rows = []
    for brca in ("included", "excluded"):
        kg = (np.ones(len(df), bool) if brca == "included"
              else (df["gene"] != "BRCA1").to_numpy())
        for key, s in scores.items():
            name = key.replace("single:", "")
            m = kg & ~np.isnan(y) & ~np.isnan(s)
            yy, ss = y[m].astype(int), s[m]
            for spec in SPEC_POINTS:
                lrp, lrm, thr, tpr, spec_obs = lr_at_specificity(yy, ss, spec=spec)
                rows.append({"condition": f"{label_col}/BRCA1_{brca}",
                             "object": name, "target_spec": spec,
                             "threshold": thr, "sensitivity": tpr,
                             "specificity_achieved": spec_obs,
                             "n": int(m.sum()), "n_pos": int((yy == 1).sum()),
                             "n_neg": int((yy == 0).sum()),
                             "lr_plus": lrp, "lr_minus": lrm,
                             "acmg_tier": acmg_tier(lrp)})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# orchestration
# ---------------------------------------------------------------------------
def run():
    p = C.OUTPUT_DIR / f"frozen_matrix_{C.FROZEN_VERSION}.parquet"
    if not p.exists():
        sys.exit(f"[phase8] frozen matrix not found at {p}")
    df = pd.read_parquet(p)
    df = df[df["is_splice"] & df["func_pathogenicity"].notna()].reset_index(drop=True)

    feats = usable_features(df)
    genes = df["gene"].astype(str)
    Xrn = rank_within_gene(df, feats)
    ytr = (rank_target_within_gene(df)
           if getattr(C, "TRAIN_TARGET", "rank_within_gene") == "rank_within_gene"
           else df["func_pathogenicity"])
    y = df["func_pathogenicity"].astype(float).to_numpy()
    gv = genes.to_numpy()

    # same fusion and same best-single criterion as phases 2/3/5/7
    oof_enet = logo_oof(Xrn, ytr, genes, feats, "enet")
    oof_gbt = logo_oof(Xrn, ytr, genes, feats, "gbt")
    oof_mean = logo_oof(Xrn, ytr, genes, feats, "mean")
    single_pg = {f: per_gene_rho(Xrn[f].to_numpy(), y, gv) for f in feats}
    single_pool = {f: dl_pool(single_pg[f])["rho"] for f in feats}
    best = max(single_pool, key=lambda k: single_pool[k] if np.isfinite(single_pool[k]) else -9)
    fus_pg = {"M1_enet": per_gene_rho(oof_enet, y, gv),
              "M2_gbt": per_gene_rho(oof_gbt, y, gv),
              "M0b_mean": per_gene_rho(oof_mean, y, gv)}

    # min-maxed scores on the full splice set, exactly as phase3/5/7 build them
    scores = {"fusion_M1": oof_enet, "mean_M0b": oof_mean}
    for f in feats:
        scores[f"single:{f}"] = Xrn[f].to_numpy()
    for k in scores:
        s = scores[k]
        lo, hi = np.nanmin(s), np.nanmax(s)
        scores[k] = (s - lo) / (hi - lo) if hi > lo else s
    best_key = f"single:{best}"

    C.REPORT_DIR.mkdir(parents=True, exist_ok=True)

    # ---- calibrated (fusion, best single) per condition: one cache shared by
    # (a) delta-Brier sign-flip (all 8 conditions) and (d) Murphy (isotonic 4).
    cal = {}
    for method in CAL_METHODS:
        for label_col in ("y_assay", "y_clinvar"):
            for brca in ("included", "excluded"):
                idx = (np.arange(len(df)) if brca == "included"
                       else np.where((df["gene"] != "BRCA1").to_numpy())[0])
                lab = df[label_col].astype(float).to_numpy()[idx]
                gsub = gv[idx]
                pf = logo_calibrate(scores["fusion_M1"][idx], lab, gsub, method=method)
                ps = logo_calibrate(scores[best_key][idx], lab, gsub, method=method)
                cal[(method, label_col, brca)] = (pf, ps, lab, gsub)

    # ---- (a) sign-flip exact + wild cluster bootstrap ----
    wild_rng = np.random.default_rng(C.RANDOM_SEED)
    sf_rows = sign_flip_rho(fus_pg["M1_enet"], single_pg[best], wild_rng)
    sf_rows += sign_flip_brier(cal, wild_rng)
    sf_df = pd.DataFrame(sf_rows)
    sf_df.to_csv(C.REPORT_DIR / "phase8_sign_flip_exact.csv", index=False)
    print("=== (a) exact sign-flip randomisation test (2^k assignments) + wild cluster bootstrap ===")
    print(sf_df.round(4).to_string(index=False))

    # ---- (b) leave-two-genes-out ----
    l2go = leave_two_genes_out(fus_pg["M1_enet"], single_pg[best])
    l2go.to_csv(C.REPORT_DIR / "phase8_leave_two_genes_out.csv", index=False)
    imin, imax = l2go.delta_rho.idxmin(), l2go.delta_rho.idxmax()
    print(f"\n=== (b) leave-two-genes-out delta-rho ({len(l2go)} pairs) ===")
    print(f"  all positive: {bool((l2go.delta_rho > 0).all())} | "
          f"min {l2go.delta_rho.min():+.4f} (drop {l2go.loc[imin, 'dropped']}) | "
          f"median {l2go.delta_rho.median():+.4f} | "
          f"max {l2go.delta_rho.max():+.4f} (drop {l2go.loc[imax, 'dropped']})")

    # ---- (c) Hartung-Knapp leaderboard (new CSV; phase2 output untouched) ----
    board = [(f"single:{f}", single_pg[f]) for f in feats] + list(fus_pg.items())
    hk_df = leaderboard_hk(board)
    hk_df = hk_df.sort_values("pooled_rho", ascending=False,
                              na_position="last").reset_index(drop=True)
    hk_df.to_csv(C.REPORT_DIR / "phase8_leaderboard_hk.csv", index=False)
    print("\n=== (c) Table-1 leaderboard: DL vs Hartung-Knapp 95% CI ===")
    print(hk_df.round(3).to_string(index=False))

    # ---- (d) Murphy decomposition (first consumer of its own stream) ----
    boot_rng = np.random.default_rng(C.RANDOM_SEED)
    mur = murphy_table(cal, boot_rng)
    mur.to_csv(C.REPORT_DIR / "phase8_murphy_decomposition.csv", index=False)
    print("\n=== (d) Murphy decomposition of the Brier score ===")
    print(mur.round(4).to_string(index=False))

    # ---- (e) stratified selection test ----
    sel = stratified_selection_test(df, scores)
    sel.to_csv(C.REPORT_DIR / "phase8_selection_test_stratified.csv", index=False)
    print("\n=== (e) selection test, stratified (LR+ at 95% specificity, assay labels) ===")
    piv = sel.pivot_table(index=["stratum", "object"], columns="subset", values="lr_plus")
    print(piv.round(2).to_string())

    # ---- (f) LR+ at four operating points ----
    lrpt = lr_plus_operating_points(df, scores)
    lrpt.to_csv(C.REPORT_DIR / "phase8_lr_plus_operating_points.csv", index=False)
    print(f"\n=== (f) LR+ at specificity {SPEC_POINTS} (y_assay, +/- BRCA1) ===")
    piv2 = lrpt.pivot_table(index=["condition", "object"], columns="target_spec",
                            values="lr_plus")
    print(piv2.round(2).to_string())

    print(f"\n[phase8] wrote phase8_sign_flip_exact.csv, "
          f"phase8_leave_two_genes_out.csv, phase8_leaderboard_hk.csv, "
          f"phase8_murphy_decomposition.csv, phase8_selection_test_stratified.csv, "
          f"phase8_lr_plus_operating_points.csv to {C.REPORT_DIR}")


if __name__ == "__main__":
    run()
