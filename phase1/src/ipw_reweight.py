"""
Inverse-probability weighting for the sampling-frame selection.

WHY
---
The analysis set is an inner join of the functional (SGE/MAVE) standard with
ClinVar *records* (`scripts/51_clean_and_merge.py`), so a variant enters only if
ClinVar has heard of it. The Introduction criticises ClinVar benchmarks for
ascertainment bias; the sampling frame carries a variant of that same bias, and
it is not independent of effect size: within the proximal splice window the
analysis set is enriched for functionally extreme variants
(`ascertainment_bias_check`: mean |z| 0.726 vs 0.582 on the 5-gene master).

This module asks the corresponding counterfactual: if the analysis set
had the *full functional set's* effect-size distribution, would any conclusion
move? No model is re-scored and no model is re-trained -- the frozen
out-of-fold predictions and the out-of-gene calibrators are exactly those of
`phase2_model` / `phase3_calibration`; only the EVALUATION is reweighted.

METHOD
------
target population : full functional SNV set, proximal splice window |offset| <= 8
                    (7 genes; `functional_scores_master.tsv` +
                    `candidate_functional_master.tsv`), n = 3,291
observed sample   : the analysis set's splice variants inside that frame,
                    n = 1,768 of 1,781 (13 have no intron offset -- 7 exon-edge
                    and 6 audited UTR-intron -- so they cannot be placed in the
                    frame; they are dropped and the drop is reported)

  1. effect size   z = functional score standardised WITHIN GENE over the target
                   frame (deposit scales differ ~100x), oriented as pathogenicity;
                   the selection variable is |z|.
  2. bins          quantile cutpoints of |z| taken on the TARGET frame, never on
                   the observed sample, so p_full(k) = 1/K by construction.
  3. weights       w_k = p_full(k) / p_int(k), renormalised to mean(w) = 1 over
                   the observed sample.
  4. metrics       weighted Spearman rho (rank-transform, then weighted Pearson
                   on the ranks), per gene, DerSimonian-Laird pooled with the
                   Fisher-z variance taken at the per-gene Kish ESS rather than n;
                   weighted Brier  sum w (p-y)^2 / sum w; weighted ECE, same bins.
  5. intervals     gene-level cluster bootstrap (1,000 resamples of the 7 genes),
                   matching the main analysis's cluster structure. Weights are
                   held fixed at their point estimate (design-based IPW).

Run:  python -m src.ipw_reweight
"""
from __future__ import annotations
import sys
import zlib
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr, ks_2samp, mannwhitneyu
from . import config as C
from .phase2_model import (usable_features, rank_within_gene,
                           rank_target_within_gene, logo_oof, per_gene_rho, dl_pool)
from .phase3_calibration import logo_calibrate, ece, brier, yield_metrics, HI, LO

N_BOOT = 1000                 # B2 spec; the main analysis uses 2,000
SPLICE_WINDOW = C.SPLICE_REGION_MAX          # |intron offset| <= 8
BIN_SCHEMES = {"quintile": 5, "decile": 10, "vigintile": 20}
PRIMARY_SCHEME = "decile"
TRUNC = (0.1, 10.0)           # sensitivity truncation for exploding weights
WEIGHT_ALERT = 10.0

FULL_SETS = ["data/processed/functional_scores_master.tsv",       # 5 original genes
             "data/processed/candidate_functional_master.tsv"]    # VHL, BAP1


# ---------------------------------------------------------------------------
# weighted estimators
# ---------------------------------------------------------------------------
def weighted_pearson(x, y, w):
    w = np.asarray(w, dtype=float)
    sw = w.sum()
    if sw <= 0:
        return np.nan
    mx, my = (w * x).sum() / sw, (w * y).sum() / sw
    cov = (w * (x - mx) * (y - my)).sum() / sw
    vx = (w * (x - mx) ** 2).sum() / sw
    vy = (w * (y - my) ** 2).sum() / sw
    return cov / np.sqrt(vx * vy) if vx > 0 and vy > 0 else np.nan


def weighted_spearman(x, y, w):
    """Rank-transform (ties -> average ranks, as scipy's spearmanr does), then a
    weighted Pearson on the ranks. With w == 1 this reduces to spearmanr."""
    x, y, w = np.asarray(x, float), np.asarray(y, float), np.asarray(w, float)
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(w)
    if m.sum() < 3:
        return np.nan
    return weighted_pearson(rankdata(x[m]), rankdata(y[m]), w[m])


def weighted_brier(prob, label, w):
    m = np.isfinite(prob) & np.isfinite(label) & np.isfinite(w)
    if m.sum() == 0:
        return np.nan
    return float((w[m] * (prob[m] - label[m]) ** 2).sum() / w[m].sum())


def weighted_ece(prob, label, w, n_bins=10):
    m = np.isfinite(prob) & np.isfinite(label) & np.isfinite(w)
    p, y, ww = prob[m], label[m], w[m]
    if len(p) == 0:
        return np.nan
    bins = np.linspace(0, 1, n_bins + 1)
    idx = np.clip(np.digitize(p, bins) - 1, 0, n_bins - 1)
    tot, e = ww.sum(), 0.0
    for b in range(n_bins):
        mb = idx == b
        if mb.sum() and ww[mb].sum() > 0:
            wb = ww[mb].sum()
            e += (wb / tot) * abs((ww[mb] * p[mb]).sum() / wb - (ww[mb] * y[mb]).sum() / wb)
    return float(e)


def weighted_yield(prob, label, w):
    m = np.isfinite(prob) & np.isfinite(label) & np.isfinite(w)
    p, y, ww = prob[m], label[m], w[m]
    if len(p) == 0:
        return {"actionable_frac": np.nan, "actionable_acc": np.nan}
    path, ben = p >= HI, p <= LO
    act = path | ben
    frac = ww[act].sum() / ww.sum() if ww.sum() > 0 else np.nan
    ok = (path & (y == 1)) | (ben & (y == 0))
    acc = ww[ok].sum() / ww[act].sum() if ww[act].sum() > 0 else np.nan
    return {"actionable_frac": float(frac), "actionable_acc": float(acc)}


def ess(w):
    """Kish effective sample size."""
    w = np.asarray(w, float)
    return float(w.sum() ** 2 / (w ** 2).sum()) if (w ** 2).sum() > 0 else np.nan


# ---------------------------------------------------------------------------
# weighted per-gene rho + DL pooling (variance at the per-gene ESS)
# ---------------------------------------------------------------------------
def per_gene_rho_w(score, y, genes, w, min_n=10):
    rows = []
    for g in pd.unique(genes):
        m = (genes == g) & np.isfinite(score) & np.isfinite(y) & np.isfinite(w)
        if m.sum() >= min_n:
            rows.append({"gene": g,
                         "rho": weighted_spearman(score[m], y[m], w[m]),
                         "n": ess(w[m])})     # DL variance uses the ESS, not n
    return pd.DataFrame(rows, columns=["gene", "rho", "n"])


def boot_delta_pg(pg_a, pg_b, rng, n_boot=N_BOOT):
    """Gene-cluster bootstrap on the pooled difference, mirroring
    phase2_model.boot_delta (resamples the per-gene rho rows)."""
    obs = dl_pool(pg_a)["rho"] - dl_pool(pg_b)["rho"]
    genes = np.intersect1d(pg_a["gene"], pg_b["gene"])
    if len(genes) < 2 or not np.isfinite(obs):
        return {"delta": obs, "lo": np.nan, "hi": np.nan, "excludes_zero": False}
    d = []
    for _ in range(n_boot):
        samp = rng.choice(genes, size=len(genes), replace=True)
        a = pd.concat([pg_a[pg_a.gene == g] for g in samp])
        b = pd.concat([pg_b[pg_b.gene == g] for g in samp])
        d.append(dl_pool(a)["rho"] - dl_pool(b)["rho"])
    lo, hi = np.nanpercentile(d, [2.5, 97.5])
    return {"delta": obs, "lo": lo, "hi": hi, "excludes_zero": bool(lo > 0 or hi < 0)}


def boot_diff_w(prob_f, prob_s, label, genes, w, fn, rng, lower_is_better=True,
                n_boot=N_BOOT):
    """Gene-clustered bootstrap on a weighted metric difference (positive =>
    fusion better), mirroring phase3_calibration.boot_diff."""
    ug = pd.unique(genes)

    def diff(idx):
        a, b = fn(prob_s[idx], label[idx], w[idx]), fn(prob_f[idx], label[idx], w[idx])
        return (a - b) if lower_is_better else (b - a)

    obs = diff(np.arange(len(label)))
    ds = [diff(np.concatenate([np.where(genes == g)[0]
                               for g in rng.choice(ug, len(ug), replace=True)]))
          for _ in range(n_boot)]
    lo, hi = np.nanpercentile(ds, [2.5, 97.5])
    return {"obs": obs, "lo": lo, "hi": hi, "excludes_zero": bool(lo > 0 or hi < 0)}


# ---------------------------------------------------------------------------
# target frame + weights
# ---------------------------------------------------------------------------
def _key(chrom, pos, ref, alt):
    return (pd.Series(chrom).astype(str).str.replace(r"\.0$", "", regex=True) + ":"
            + pd.Series(pos).astype(float).astype("int64").astype(str) + ":"
            + pd.Series(ref).astype(str) + ":" + pd.Series(alt).astype(str))


def load_target_frame(repo_root: "Path"):
    """Full functional SNV set restricted to the proximal splice window."""
    parts = []
    for rel in FULL_SETS:
        p = repo_root / rel
        if not p.exists():
            sys.exit(f"[ipw] full functional set not found: {p}")
        d = pd.read_csv(p, sep="\t", low_memory=False)
        parts.append(d)
    full = pd.concat(parts, ignore_index=True)
    full = full[(full["is_snv"] == True) & full["functional_score"].notna()
                & full["pos"].notna()].copy()
    full["key"] = _key(full.chrom, full.pos, full.ref, full.alt).to_numpy()
    full = full.drop_duplicates("key")
    frame = full[full["intron_offset"].abs() <= SPLICE_WINDOW].copy()
    # functional pathogenicity = -score, standardised within gene over the FRAME
    g = frame.groupby("gene")["functional_score"]
    frame["z"] = -(frame["functional_score"] - g.transform("mean")) / g.transform("std")
    frame["abs_z"] = frame["z"].abs()
    return frame


def make_weights(abs_z_frame, abs_z_obs, n_bins):
    """Quantile bins fixed on the TARGET frame; w_k = p_full(k)/p_int(k),
    renormalised to mean(w) = 1 on the observed sample."""
    qs = np.linspace(0, 1, n_bins + 1)[1:-1]
    cuts = np.unique(np.quantile(abs_z_frame, qs))
    edges = np.concatenate([[-np.inf], cuts, [np.inf]])
    k_frame = np.digitize(abs_z_frame, edges) - 1
    k_obs = np.digitize(abs_z_obs, edges) - 1
    K = len(edges) - 1
    p_full = np.array([(k_frame == k).mean() for k in range(K)])
    p_int = np.array([(k_obs == k).mean() for k in range(K)])
    with np.errstate(divide="ignore", invalid="ignore"):
        w_k = np.where(p_int > 0, p_full / p_int, np.nan)
    w = w_k[k_obs]
    w = w / np.nanmean(w)                     # mean(w) == 1 on the observed sample
    per_bin = pd.DataFrame({"bin": np.arange(K), "p_full": p_full, "p_int": p_int,
                            "n_frame": [int((k_frame == k).sum()) for k in range(K)],
                            "n_obs": [int((k_obs == k).sum()) for k in range(K)],
                            "w_raw": w_k})
    per_bin["w_norm"] = per_bin["w_raw"] / np.nanmean(w_k[k_obs])
    return w, per_bin, edges


def frame_balance(frame, observed_keys):
    """How enriched is the observed sample within the target frame? The published
    ascertainment check measured this over WHOLE genes; the analysis stratum is
    the proximal splice window, and the answer there is what the weights see."""
    f = frame.assign(observed=frame["key"].isin(observed_keys))
    rows = []
    for g, d in list(f.groupby("gene")) + [("POOLED", f)]:
        a, r = d[d.observed], d[~d.observed]
        lo, hi = d["z"].quantile([.10, .90])
        row = {"gene": g, "n_frame": len(d), "n_observed": len(a),
               "n_remainder": len(r),
               "pct_covered": round(100 * len(a) / len(d), 1),
               "mean_abs_z_frame": d["abs_z"].mean(),
               "mean_abs_z_observed": a["abs_z"].mean(),
               "mean_abs_z_remainder": r["abs_z"].mean() if len(r) else np.nan,
               "extreme_decile_frac_observed": ((a.z < lo) | (a.z > hi)).mean(),
               "extreme_decile_frac_remainder": (((r.z < lo) | (r.z > hi)).mean()
                                                 if len(r) else np.nan)}
        if len(r) >= 30:
            row["ks_p"] = f"{ks_2samp(a.z, r.z).pvalue:.3e}"
            row["mwu_abs_z_p"] = f"{mannwhitneyu(a.abs_z, r.abs_z).pvalue:.3e}"
        else:
            row["ks_p"] = row["mwu_abs_z_p"] = "n/a (fully covered)"
        rows.append(row)
    return pd.DataFrame(rows)


def weight_diagnostics(w, tag):
    return {"scheme": tag, "n": int(len(w)), "w_min": float(np.min(w)),
            "w_median": float(np.median(w)), "w_max": float(np.max(w)),
            "ESS": ess(w), "ESS_over_n": ess(w) / len(w),
            "n_w_gt_alert": int((w > WEIGHT_ALERT).sum())}


# ---------------------------------------------------------------------------
# guardrails -- these gate the analysis; see also tests/test_ipw_reweight.py
# ---------------------------------------------------------------------------
def guardrails(rng=None):
    rng = rng or np.random.default_rng(0)
    n = 400
    x = rng.normal(size=n)
    y = 0.6 * x + rng.normal(size=n)
    one = np.ones(n)
    # (a) weighted Spearman degenerates to the unweighted one when w == 1
    assert abs(weighted_spearman(x, y, one) - spearmanr(x, y).statistic) < 1e-10, \
        "weighted Spearman does not reduce to spearmanr at w == 1"
    # ... including with ties, where average ranks matter
    xt, yt = np.round(x, 1), np.round(y, 1)
    assert abs(weighted_spearman(xt, yt, one) - spearmanr(xt, yt).statistic) < 1e-10, \
        "weighted Spearman disagrees with spearmanr under ties"
    # (b) weighted Brier / ECE degenerate too
    p = rng.uniform(size=n); lab = (rng.uniform(size=n) < p).astype(float)
    assert abs(weighted_brier(p, lab, one) - brier(p, lab)) < 1e-12
    assert abs(weighted_ece(p, lab, one) - ece(p, lab)) < 1e-12
    assert abs(weighted_yield(p, lab, one)["actionable_frac"]
               - (yield_metrics(p, lab)["actionable_frac"] or 0.0)) < 1e-4
    # (c) a weight of 2 is exactly a duplicated observation
    ww = np.ones(n); ww[:50] = 2.0
    assert abs(weighted_brier(np.concatenate([p, p[:50]]),
                              np.concatenate([lab, lab[:50]]), np.ones(n + 50))
               - weighted_brier(p, lab, ww)) < 1e-12, \
        "weighted Brier is not equivalent to replication"
    # (d) mean(w) == 1 after normalisation, on real-shaped inputs
    frame = rng.normal(size=5000) ** 2
    obs = rng.normal(size=1500) ** 2 + 0.3
    for k in (5, 10, 20):
        w, per_bin, _ = make_weights(frame, obs, k)
        assert abs(np.nanmean(w) - 1.0) < 1e-9, f"mean(w) != 1 at {k} bins"
        assert np.isfinite(w).all(), f"non-finite weight at {k} bins"
    return True


# ---------------------------------------------------------------------------
# orchestration
# ---------------------------------------------------------------------------
def run():
    from pathlib import Path
    guardrails()
    here = Path(__file__).resolve()
    repo_root = here.parents[2]

    p = C.OUTPUT_DIR / f"frozen_matrix_{C.FROZEN_VERSION}.parquet"
    if not p.exists():
        sys.exit(f"[ipw] frozen matrix not found at {p}")
    df = pd.read_parquet(p)
    df = df[df["is_splice"] & df["func_pathogenicity"].notna()].reset_index(drop=True)
    df["key"] = _key(df.chrom, df.pos, df.ref, df.alt).to_numpy()

    frame = load_target_frame(repo_root)

    # --- observed sample = analysis set inside the target frame ------------
    fz = frame.set_index("key")["abs_z"]
    df["abs_z"] = df["key"].map(fz)
    in_frame = df["abs_z"].notna().to_numpy()
    n_drop = int((~in_frame).sum())
    print(f"[ipw] target frame (full functional, |offset| <= {SPLICE_WINDOW}): "
          f"{len(frame):,} variants, {frame.gene.nunique()} genes")
    print(f"[ipw] analysis splice set: {len(df):,}; inside the frame: "
          f"{int(in_frame.sum()):,}; dropped (no intron offset): {n_drop}")
    print("[ipw] per-gene coverage of the frame by the analysis set:")
    cov = (frame.assign(observed=frame.key.isin(set(df.key)))
                .groupby("gene")["observed"].agg(["sum", "size", "mean"]))
    print(cov.to_string())

    # --- frozen predictions: identical to phase2/phase3, no re-training ----
    feats = usable_features(df)
    genes_s = df["gene"].astype(str)
    Xrn = rank_within_gene(df, feats)
    ytr = rank_target_within_gene(df)
    oof_enet = logo_oof(Xrn, ytr, genes_s, feats, "enet")
    oof_mean = logo_oof(Xrn, ytr, genes_s, feats, "mean")
    yraw = df["func_pathogenicity"].to_numpy()
    gv = genes_s.to_numpy()
    single_pool = {f: dl_pool(per_gene_rho(Xrn[f].to_numpy(), yraw, gv))["rho"]
                   for f in feats}
    best = max(single_pool, key=lambda k: single_pool[k]
               if np.isfinite(single_pool[k]) else -9)
    print(f"[ipw] best single tool (unweighted, pooled per-gene Spearman) = {best}")

    scores = {"fusion_M1": oof_enet, "mean_M0b": oof_mean}
    for f in feats:
        scores[f"single:{f}"] = Xrn[f].to_numpy()
    scores_raw = {k: v.copy() for k, v in scores.items()}      # for rho
    for k in scores:                                           # min-max for calibration
        s = scores[k]
        lo, hi = np.nanmin(s), np.nanmax(s)
        scores[k] = (s - lo) / (hi - lo) if hi > lo else s
    best_key = f"single:{best}"

    # --- weights ----------------------------------------------------------
    obs_abs_z = df.loc[in_frame, "abs_z"].to_numpy()
    weights, per_bin, diag_rows = {}, {}, []
    for tag, k in BIN_SCHEMES.items():
        w, pb, _ = make_weights(frame["abs_z"].to_numpy(), obs_abs_z, k)
        weights[tag], per_bin[tag] = w, pb.assign(scheme=tag)
        diag_rows.append(weight_diagnostics(w, tag))
    w_pri = weights[PRIMARY_SCHEME]
    w_trunc = np.clip(w_pri, *TRUNC)
    w_trunc = w_trunc / w_trunc.mean()
    diag_rows.append(weight_diagnostics(w_trunc, f"{PRIMARY_SCHEME}_truncated"))
    diag_rows.append(weight_diagnostics(np.ones(len(w_pri)), "unweighted"))
    diag = pd.DataFrame(diag_rows)
    wsets = {"unweighted": np.ones(len(w_pri)), **weights,
             f"{PRIMARY_SCHEME}_truncated": w_trunc}

    C.REPORT_DIR.mkdir(parents=True, exist_ok=True)
    frame_balance(frame, set(df["key"])).to_csv(
        C.REPORT_DIR / "ipw_frame_balance.csv", index=False)
    diag.to_csv(C.REPORT_DIR / "ipw_weight_diagnostics.csv", index=False)
    pd.concat(per_bin.values(), ignore_index=True).to_csv(
        C.REPORT_DIR / "ipw_weight_bins.csv", index=False)
    print("\n=== selection inside the target frame (proximal splice window) ===")
    print(frame_balance(frame, set(df["key"])).round(4).to_string(index=False))
    print("\n=== weight diagnostics ===")
    print(diag.to_string(index=False))

    # everything below is evaluated on the in-frame subset
    sub = df[in_frame].reset_index(drop=True)
    gsub = sub["gene"].astype(str).to_numpy()
    ysub = sub["func_pathogenicity"].to_numpy()
    idx = np.where(in_frame)[0]
    raw_sub = {k: v[idx] for k, v in scores_raw.items()}
    cal_in = {k: v[idx] for k, v in scores.items()}

    # --- Q1/Q2: ranking ----------------------------------------------------
    # Each scheme gets its own deterministic stream, so two schemes that produce
    # identical weights produce identical intervals (a truncation that changes
    # nothing must LOOK like it changed nothing).
    def stream(w, what):
        # Seeded from the WEIGHT VECTOR, not the scheme name, so two schemes that
        # produce identical weights produce identical intervals -- a truncation
        # that changed nothing has to look like it changed nothing. crc32 rather
        # than hash(): str hashing is salted per process and would not reproduce.
        wsig = zlib.crc32(np.round(np.asarray(w, float), 12).tobytes())
        return np.random.default_rng([C.RANDOM_SEED, wsig,
                                      zlib.crc32(what.encode())])

    rank_rows = []
    for scheme, w in wsets.items():
        pg = {k: per_gene_rho_w(v, ysub, gsub, w) for k, v in raw_sub.items()}
        pooled = {k: dl_pool(v)["rho"] for k, v in pg.items()}
        order = sorted([k for k in pooled if np.isfinite(pooled[k])],
                       key=lambda k: -pooled[k])
        for r, k in enumerate(order, 1):
            rank_rows.append({"scheme": scheme, "model": k, "rho": pooled[k], "rank": r})
        if scheme in ("unweighted", PRIMARY_SCHEME, f"{PRIMARY_SCHEME}_truncated",
                      "quintile", "vigintile"):
            b = boot_delta_pg(pg["fusion_M1"], pg[best_key], stream(w, "rho"))
            rank_rows.append({"scheme": scheme, "model": f"DELTA fusion-{best}",
                              "rho": b["delta"], "rank": np.nan,
                              "ci_lo": b["lo"], "ci_hi": b["hi"],
                              "excludes_zero": b["excludes_zero"]})
    rank_df = pd.DataFrame(rank_rows)
    rank_df.to_csv(C.REPORT_DIR / "ipw_ranking.csv", index=False)

    # --- Q3: calibration ---------------------------------------------------
    cal_rows, head_rows = [], []
    for label_col in ["y_assay", "y_clinvar"]:
        for brca in ["included", "excluded"]:
            keep = np.ones(len(sub), bool) if brca == "included" else (gsub != "BRCA1")
            lab = sub[label_col].astype(float).to_numpy()
            g_k, lab_k = gsub[keep], lab[keep]
            cal = {k: logo_calibrate(v[keep], lab_k, g_k, method="isotonic")
                   for k, v in cal_in.items()}
            for scheme, w in wsets.items():
                wk = w[keep] / w[keep].mean()
                for k, pr in cal.items():
                    cal_rows.append({"set": f"{label_col}/BRCA1_{brca}", "scheme": scheme,
                                     "model": k,
                                     "Brier": weighted_brier(pr, lab_k, wk),
                                     "ECE": weighted_ece(pr, lab_k, wk),
                                     **weighted_yield(pr, lab_k, wk)})
                tag = f"{label_col}/{brca}"
                d_b = boot_diff_w(cal["fusion_M1"], cal[best_key], lab_k, g_k, wk,
                                  weighted_brier, stream(wk, tag + "/brier"),
                                  lower_is_better=True)
                d_e = boot_diff_w(cal["fusion_M1"], cal[best_key], lab_k, g_k, wk,
                                  weighted_ece, stream(wk, tag + "/ece"),
                                  lower_is_better=True)
                head_rows.append({"set": f"{label_col}/BRCA1_{brca}", "scheme": scheme,
                                  "dBrier": d_b["obs"], "dBrier_lo": d_b["lo"],
                                  "dBrier_hi": d_b["hi"],
                                  "dBrier_excludes_zero": d_b["excludes_zero"],
                                  "dECE": d_e["obs"], "dECE_lo": d_e["lo"],
                                  "dECE_hi": d_e["hi"]})
    cal_df = pd.DataFrame(cal_rows)
    head_df = pd.DataFrame(head_rows)
    cal_df.to_csv(C.REPORT_DIR / "ipw_calibration.csv", index=False)
    head_df.to_csv(C.REPORT_DIR / "ipw_headline.csv", index=False)

    # --- B5 deliverable: one row per predictor -----------------------------
    unw = rank_df[(rank_df.scheme == "unweighted") & rank_df["rank"].notna()]
    wgt = rank_df[(rank_df.scheme == PRIMARY_SCHEME) & rank_df["rank"].notna()]
    cu = cal_df[(cal_df.set == "y_assay/BRCA1_included") & (cal_df.scheme == "unweighted")]
    cw = cal_df[(cal_df.set == "y_assay/BRCA1_included") & (cal_df.scheme == PRIMARY_SCHEME)]
    tab = (unw.set_index("model")[["rho", "rank"]].rename(columns={"rho": "rho_unw", "rank": "rank_unw"})
             .join(wgt.set_index("model")[["rho", "rank"]].rename(columns={"rho": "rho_wgt", "rank": "rank_wgt"}))
             .join(cu.set_index("model")[["Brier"]].rename(columns={"Brier": "Brier_unw"}))
             .join(cw.set_index("model")[["Brier"]].rename(columns={"Brier": "Brier_wgt"})))
    tab["d_rho"] = tab.rho_wgt - tab.rho_unw
    tab["d_Brier"] = tab.Brier_wgt - tab.Brier_unw
    tab["rank_change"] = tab.rank_unw - tab.rank_wgt
    tab = tab.sort_values("rank_unw").reset_index().rename(columns={"index": "model"})
    tab = tab[["model", "rho_unw", "rho_wgt", "d_rho", "Brier_unw", "Brier_wgt",
               "d_Brier", "rank_unw", "rank_wgt", "rank_change"]]
    tab.to_csv(C.REPORT_DIR / "ipw_predictor_table.csv", index=False)

    pd.set_option("display.width", 200, "display.max_columns", 40)
    print("\n=== B5 predictor table (rho on the splice set; Brier on "
          "y_assay/BRCA1_included, isotonic) ===")
    print(tab.round(4).to_string(index=False))
    print("\n=== Q1 ranking under each weighting scheme ===")
    print(rank_df[rank_df["rank"].notna()].pivot(index="model", columns="scheme",
          values="rank").to_string())
    print("\n=== Q2 ranking saturation: delta-rho fusion - best single ===")
    print(rank_df[rank_df["rank"].isna()].round(4).to_string(index=False))
    print("\n=== Q3 calibration headline: delta-Brier (best single - fusion) ===")
    print(head_df.round(4).to_string(index=False))
    print(f"\n[ipw] wrote ipw_weight_diagnostics.csv, ipw_weight_bins.csv, "
          f"ipw_frame_balance.csv, ipw_ranking.csv, ipw_calibration.csv, ipw_headline.csv, "
          f"ipw_predictor_table.csv to {C.REPORT_DIR}\n")


if __name__ == "__main__":
    run()
