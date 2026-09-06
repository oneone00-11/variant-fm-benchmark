"""
Phase 2 -- fusion vs single tool (H1), stratified by splice_bin, plus the
evolution-axis ablation (H2). Pinned to frozen-matrix-v1.

Leakage-safe under leave-one-gene-out (LOGO): within-gene rank-normalisation is
self-contained per gene; imputation medians fit on training genes only; feature
orientation is a fixed monotonic flip (no target peeking inside folds).

Models: M0b equal-weight mean | M1 elastic-net (primary) | M2 HGBR (robustness).
Metric: per-gene Spearman rho, DerSimonian-Laird pooled (Fisher-z), I^2.
H1: delta-rho vs best single, gene-clustered bootstrap 95% CI -- reported OVERALL
    and stratified by splice_bin (core_like vs region).
H2: elastic-net ablation of the evolutionary axis, two versions
    (drop conservation only; drop conservation+alignment), delta-rho + coefs.

Run:  python -m src.phase2_model
"""
from __future__ import annotations
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata, t as tdist
from sklearn.linear_model import ElasticNetCV
from sklearn.ensemble import HistGradientBoostingRegressor
from . import config as C

RNG = np.random.default_rng(C.RANDOM_SEED)
N_BOOT = 2000


# ---------------------------------------------------------------------------
# feature assembly + orientation
# ---------------------------------------------------------------------------
def usable_features(df):
    return [f for f in dict.fromkeys(C.PREDICTOR_FEATURES + C.EVO_FEATURES)
            if C.COLUMNS.get(f) and f in df.columns]

def rank_within_gene(df, feats):
    """Within-gene rank -> [0,1]; then flip REVERSED_FEATURES to 1-r so larger =
    more damaging for all. NaN preserved. Self-contained per gene -> LOGO-safe."""
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
    for f in feats:
        if f in getattr(C, "REVERSED_FEATURES", []):
            out[f] = 1.0 - out[f]
    return out


# ---------------------------------------------------------------------------
# LOGO out-of-fold predictions for a given feature set + model
# ---------------------------------------------------------------------------
def _en_design(Xrn, mask, feats, med):
    base = Xrn.loc[mask, feats].fillna(med)
    isna = Xrn.loc[mask, feats].isna().astype(float)
    isna.columns = [f"{c}_isna" for c in feats]
    return pd.concat([base, isna], axis=1)

def rank_target_within_gene(df, col="func_pathogenicity"):
    """Within-gene rank of the target -> [0,1]. Self-contained per gene (LOGO-safe);
    training on this aligns the objective with the per-gene rank evaluation."""
    out = pd.Series(np.nan, index=df.index, dtype=float)
    for _, idx in df.groupby("gene", observed=True).groups.items():
        v = df.loc[idx, col].to_numpy(dtype=float); m = ~np.isnan(v)
        r = np.full_like(v, np.nan)
        if m.sum() > 1: r[m] = (rankdata(v[m]) - 0.5) / m.sum()
        elif m.sum() == 1: r[m] = 0.5
        out.loc[idx] = r
    return out


def logo_oof(Xrn, y, genes, feats, model, collect_coefs=False):
    oof = np.full(len(y), np.nan)
    coefs = []
    for gtest in genes.unique():
        tr = (genes != gtest).to_numpy(); te = (genes == gtest).to_numpy()
        if model == "mean":
            oof[te] = np.nanmean(Xrn.loc[te, feats].to_numpy(), axis=1)
        elif model == "enet":
            med = Xrn.loc[tr, feats].median().fillna(0.5)
            Xtr = _en_design(Xrn, tr, feats, med)
            en = ElasticNetCV(l1_ratio=[.1, .5, .9], n_alphas=50, cv=5,
                              random_state=C.RANDOM_SEED, max_iter=5000)
            en.fit(Xtr.to_numpy(), y[tr].to_numpy())
            oof[te] = en.predict(_en_design(Xrn, te, feats, med).to_numpy())
            if collect_coefs:
                # keyed by held-out gene: the per-fold weights are what a
                # stability claim needs, and the cross-fold mean hides exactly
                # the dispersion the question is about
                coefs.append(pd.Series(en.coef_, index=Xtr.columns, name=gtest))
        elif model == "gbt":
            gb = HistGradientBoostingRegressor(max_depth=3, max_iter=300,
                    learning_rate=0.05, l2_regularization=1.0,
                    random_state=C.RANDOM_SEED, early_stopping=True)
            gb.fit(Xrn.loc[tr, feats].to_numpy(), y[tr].to_numpy())
            oof[te] = gb.predict(Xrn.loc[te, feats].to_numpy())
    if collect_coefs:
        # DataFrame: feature x held-out gene. Callers wanting the old summary
        # take .abs().mean(axis=1).sort_values(ascending=False) themselves.
        return oof, pd.concat(coefs, axis=1)
    return oof


# ---------------------------------------------------------------------------
# per-gene Spearman + DerSimonian-Laird pooling
# ---------------------------------------------------------------------------
def per_gene_rho(score, y, genes, min_n=10):
    rows = []
    for g in pd.unique(genes):
        m = (genes == g) & ~np.isnan(score) & ~np.isnan(y)
        if m.sum() >= min_n:
            rows.append({"gene": g, "rho": spearmanr(score[m], y[m]).statistic, "n": int(m.sum())})
    return pd.DataFrame(rows, columns=["gene", "rho", "n"])

def dl_pool(pg, method="dl"):
    """DerSimonian-Laird random-effects pooling of per-gene Spearman rhos
    (Fisher-z). method="dl" gives the historical normal-theory 95% CI;
    method="hk" gives the Hartung-Knapp interval -- the same DL tau^2 and
    weights, but the pooled variance is the HK quadratic form and the critical
    value is t(k-1) instead of 1.96. At k=7 genes the normal CI is
    anti-conservative; the HK interval is the calibrated one and is reported
    alongside DL in phase8's leaderboard (phase8_leaderboard_hk.csv)."""
    if method not in ("dl", "hk"):
        raise ValueError(f"unknown pooling method {method!r}; expected 'dl' or 'hk'")
    if pg is None or "rho" not in pg.columns or len(pg) == 0:
        return {"rho": np.nan, "lo": np.nan, "hi": np.nan, "I2": np.nan, "k": 0}
    d = pg.dropna(subset=["rho"]); d = d[d["n"] > 3]
    if len(d) < 2:
        return {"rho": d["rho"].mean() if len(d) else np.nan,
                "lo": np.nan, "hi": np.nan, "I2": np.nan, "k": len(d)}
    z = np.arctanh(d["rho"].clip(-0.999, 0.999).to_numpy())
    v = 1.0 / (d["n"].to_numpy() - 3); w = 1.0 / v
    z_fe = np.sum(w * z) / np.sum(w)
    Q = np.sum(w * (z - z_fe) ** 2); k = len(d)
    Cc = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0.0, (Q - (k - 1)) / Cc) if Cc > 0 else 0.0
    ws = 1.0 / (v + tau2)
    z_re = np.sum(ws * z) / np.sum(ws)
    I2 = max(0.0, (Q - (k - 1)) / Q) * 100 if Q > 0 else 0.0
    if method == "hk":
        q_hk = np.sum(ws * (z - z_re) ** 2) / (k - 1)
        se = np.sqrt(q_hk / np.sum(ws))
        crit = float(tdist.ppf(0.975, k - 1))
    else:
        se = np.sqrt(1.0 / np.sum(ws))
        crit = 1.96
    return {"rho": np.tanh(z_re), "lo": np.tanh(z_re - crit * se),
            "hi": np.tanh(z_re + crit * se), "I2": I2, "k": k}

def boot_delta(fusion_pg, single_pg):
    obs = dl_pool(fusion_pg)["rho"] - dl_pool(single_pg)["rho"]
    genes = np.intersect1d(fusion_pg["gene"], single_pg["gene"])
    if len(genes) < 2 or not np.isfinite(obs):
        return {"delta": obs, "lo": np.nan, "hi": np.nan, "excludes_zero": False}
    deltas = []
    for _ in range(N_BOOT):
        samp = RNG.choice(genes, size=len(genes), replace=True)
        fb = pd.concat([fusion_pg[fusion_pg.gene == g] for g in samp])
        sb = pd.concat([single_pg[single_pg.gene == g] for g in samp])
        deltas.append(dl_pool(fb)["rho"] - dl_pool(sb)["rho"])
    lo, hi = np.nanpercentile(deltas, [2.5, 97.5])
    return {"delta": obs, "lo": lo, "hi": hi, "excludes_zero": bool(lo > 0 or hi < 0)}

def h1_on(mask, oof_fusion, single_score, y, genes, label):
    fpg = per_gene_rho(oof_fusion[mask], y[mask], genes[mask])
    spg = per_gene_rho(single_score[mask], y[mask], genes[mask])
    b = boot_delta(fpg, spg)
    return {"stratum": label, "n": int(mask.sum()),
            "fusion_rho": round(dl_pool(fpg)["rho"], 4),
            "single_rho": round(dl_pool(spg)["rho"], 4),
            "delta": round(b["delta"], 4) if np.isfinite(b["delta"]) else np.nan,
            "ci95": f"[{b['lo']:.4f}, {b['hi']:.4f}]" if np.isfinite(b["lo"]) else "n/a",
            "H1": "SUPPORTED" if b["excludes_zero"] and b["delta"] > 0 else "ns"}


# ---------------------------------------------------------------------------
# orchestration
# ---------------------------------------------------------------------------
def run():
    p = C.OUTPUT_DIR / f"frozen_matrix_{C.FROZEN_VERSION}.parquet"
    if not p.exists():
        sys.exit(f"[phase2] frozen matrix not found at {p}")
    df = pd.read_parquet(p)
    df = df[df["is_splice"] & df["func_pathogenicity"].notna()].reset_index(drop=True)

    feats = usable_features(df)
    y = df["func_pathogenicity"].astype(float)                 # raw -> EVALUATION only
    if getattr(C, "TRAIN_TARGET", "rank_within_gene") == "rank_within_gene":
        ytrain = rank_target_within_gene(df)                   # aligned -> TRAINING
    else:
        ytrain = y
    genes = df["gene"].astype(str)
    sbin = df["splice_bin"].astype(str)
    Xrn = rank_within_gene(df, feats)

    # ---- single-tool pooled rho (oriented) -> best single ----
    single_pg = {f: per_gene_rho(Xrn[f].to_numpy(), y.to_numpy(), genes.to_numpy()) for f in feats}
    single_pool = {f: dl_pool(single_pg[f])["rho"] for f in feats}
    best = max(single_pool, key=lambda k: single_pool[k] if np.isfinite(single_pool[k]) else -9)

    # ---- fusion models (full feature set) ----
    oof_enet, coefs = logo_oof(Xrn, ytrain, genes, feats, "enet", collect_coefs=True)
    oof_gbt = logo_oof(Xrn, ytrain, genes, feats, "gbt")
    oof_mean = logo_oof(Xrn, ytrain, genes, feats, "mean")
    fus_pg = {"M1_enet": per_gene_rho(oof_enet, y.to_numpy(), genes.to_numpy()),
              "M2_gbt":  per_gene_rho(oof_gbt,  y.to_numpy(), genes.to_numpy()),
              "M0b_mean":per_gene_rho(oof_mean, y.to_numpy(), genes.to_numpy())}

    # ---- leaderboard ----
    board = [("single:" + f, dl_pool(single_pg[f])) for f in feats]
    board += [(m, dl_pool(fus_pg[m])) for m in fus_pg]
    board_df = pd.DataFrame([{"model": n, "pooled_rho": round(d["rho"], 4),
        "ci95": f"[{d['lo']:.3f},{d['hi']:.3f}]" if np.isfinite(d['lo']) else "n/a",
        "I2": round(d["I2"], 1), "k": d["k"]} for n, d in board]
        ).sort_values("pooled_rho", ascending=False, na_position="last").reset_index(drop=True)

    # ---- stratified H1 (M1 primary vs best single) ----
    ally = y.to_numpy(); allg = genes.to_numpy()
    strat = [h1_on(np.ones(len(df), bool), oof_enet, Xrn[best].to_numpy(), ally, allg, "overall")]
    for b in ["core_like", "region"]:
        m = (sbin == b).to_numpy()
        if m.sum() > 0:
            strat.append(h1_on(m, oof_enet, Xrn[best].to_numpy(), ally, allg, b))
    strat_df = pd.DataFrame(strat)

    # ---- H2 evolution-axis ablation (M1 primary; M2 robustness) ----
    cons = [f for f in C.EVO_CONSERVATION if f in feats]
    align = [f for f in C.EVO_ALIGNMENT if f in feats]
    sets = {"full": feats,
            "drop_conservation": [f for f in feats if f not in cons],
            "drop_all_evo": [f for f in feats if f not in cons + align]}
    h2 = []
    for model in ("enet", "gbt"):
        oof_by = {name: logo_oof(Xrn, ytrain, genes, fs, model) for name, fs in sets.items()}
        pg_by = {name: per_gene_rho(oof_by[name], ally, allg) for name in sets}
        full_pool = dl_pool(pg_by["full"])["rho"]
        for name in ("drop_conservation", "drop_all_evo"):
            b = boot_delta(pg_by["full"], pg_by[name])
            h2.append({"model": "M1_enet" if model == "enet" else "M2_gbt",
                       "ablation": name, "full_rho": round(full_pool, 4),
                       "ablated_rho": round(dl_pool(pg_by[name])["rho"], 4),
                       "delta_full_minus_ablated": round(b["delta"], 4),
                       "ci95": f"[{b['lo']:.4f}, {b['hi']:.4f}]" if np.isfinite(b['lo']) else "n/a",
                       "evo_contributes": "YES" if b["excludes_zero"] and b["delta"] > 0 else "ns"})
    h2_df = pd.DataFrame(h2)

    # evo coefficients in the full elastic-net
    evo_all = cons + align
    per_fold = coefs                      # feature x held-out gene
    coefs = coefs.abs().mean(axis=1).sort_values(ascending=False)
    coef_evo = coefs[[c for c in coefs.index if c in evo_all]]
    coef_rank = {c: int(np.where(coefs.index == c)[0][0]) + 1 for c in coefs.index if c in evo_all}

    C.REPORT_DIR.mkdir(parents=True, exist_ok=True)
    board_df.to_csv(C.REPORT_DIR / "phase2_leaderboard.csv", index=False)
    strat_df.to_csv(C.REPORT_DIR / "phase2_H1_stratified.csv", index=False)
    h2_df.to_csv(C.REPORT_DIR / "phase2_H2_ablation.csv", index=False)

    print(f"[phase2] splice set n={len(df)} | genes={genes.nunique()} | features={len(feats)} | best single={best}\n")
    print("=== LEADERBOARD ==="); print(board_df.to_string(index=False))
    print("\n=== H1 stratified (M1 elastic net vs best single) ==="); print(strat_df.to_string(index=False))
    print("\n=== H2 evolution-axis ablation ==="); print(h2_df.to_string(index=False))
    print("\n[H2] elastic-net |coef| of evo features (full model):")
    for c in coef_evo.index:
        print(f"   {c:12s} |coef|={coef_evo[c]:.4f}   rank {coef_rank[c]}/{len(coefs)}")


if __name__ == "__main__":
    run()
