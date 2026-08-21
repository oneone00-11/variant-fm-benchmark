"""
Phase 3 -- calibration (H3). The make-or-break test, orthogonal to H1/H2.

Question: after mapping scores to pathogenicity PROBABILITIES, is the fusion
better *calibrated* than the best single tool -- i.e. when it says 80% it is
really ~80%? This can be true even when ranking (rho) is a dead heat.

FAIRNESS (the crux): every object -- fusion M1, best single, equal-weight mean --
is put through the IDENTICAL out-of-gene isotonic calibration. We compare
calibrated-vs-calibrated, never calibrated-vs-raw. Calibration is fit under LOGO
(fit on the other genes, applied to the held-out gene), matching H1/H2 geometry.

Label: y_assay (official functional labels) primary; y_clinvar sensitivity.
Metrics: ECE, Brier, and clinical yield (actionable fraction at a controlled
confidence). Headline: delta-ECE (best single - fusion), gene-clustered bootstrap.
Everything reported with and without BRCA1.

Run:  python -m src.phase3_calibration
"""
from __future__ import annotations
import sys
import numpy as np
import pandas as pd
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LogisticRegression
from scipy.stats import norm
from . import config as C
from .phase2_model import (usable_features, rank_within_gene,
                           rank_target_within_gene, logo_oof)

RNG = np.random.default_rng(C.RANDOM_SEED)
N_BOOT = 2000
HI, LO = 0.90, 0.10   # calibrated-confidence thresholds for "actionable" calls
CAL_METHODS = ("isotonic", "platt")   # robustness: results must hold under both


# ---------------------------------------------------------------------------
# calibration metrics
# ---------------------------------------------------------------------------
def ece(prob, label, n_bins=10):
    m = ~np.isnan(prob) & ~np.isnan(label)
    p, y = prob[m], label[m]
    if len(p) == 0:
        return np.nan
    bins = np.linspace(0, 1, n_bins + 1)
    idx = np.clip(np.digitize(p, bins) - 1, 0, n_bins - 1)
    e = 0.0
    for b in range(n_bins):
        mb = idx == b
        if mb.sum():
            e += (mb.sum() / len(p)) * abs(p[mb].mean() - y[mb].mean())
    return e

def brier(prob, label):
    m = ~np.isnan(prob) & ~np.isnan(label)
    return np.mean((prob[m] - label[m]) ** 2) if m.sum() else np.nan

def reliability(prob, label, n_bins=10):
    m = ~np.isnan(prob) & ~np.isnan(label)
    p, y = prob[m], label[m]
    bins = np.linspace(0, 1, n_bins + 1)
    idx = np.clip(np.digitize(p, bins) - 1, 0, n_bins - 1)
    rows = []
    for b in range(n_bins):
        mb = idx == b
        if mb.sum():
            rows.append({"bin_mid": (bins[b] + bins[b+1]) / 2,
                         "mean_pred": p[mb].mean(), "frac_pos": y[mb].mean(),
                         "n": int(mb.sum())})
    return pd.DataFrame(rows)

def yield_metrics(prob, label):
    m = ~np.isnan(prob) & ~np.isnan(label)
    p, y = prob[m], label[m]
    path, ben = p >= HI, p <= LO
    act = path | ben
    correct = np.sum((path & (y == 1)) | (ben & (y == 0)))
    return {"actionable_frac": round(float(act.mean()), 4) if len(p) else np.nan,
            "actionable_acc": round(float(correct / act.sum()), 4) if act.sum() else np.nan,
            "n_actionable": int(act.sum())}


# ---------------------------------------------------------------------------
# identical out-of-gene isotonic calibration for any score
# ---------------------------------------------------------------------------
def logo_calibrate(score, label, genes, method="isotonic"):
    """Identical out-of-gene calibration for any score. method='isotonic'
    (non-parametric) or 'platt' (logistic / Platt scaling, lower variance).
    Both are run so H3 cannot be an artefact of one calibrator."""
    cal = np.full(len(score), np.nan)
    for g in pd.unique(genes):
        tr = (genes != g) & ~np.isnan(score) & ~np.isnan(label)
        te = (genes == g) & ~np.isnan(score)
        if tr.sum() < 20 or te.sum() == 0:
            continue
        if method == "isotonic":
            cal_fn = IsotonicRegression(out_of_bounds="clip", y_min=0.0, y_max=1.0)
            cal_fn.fit(score[tr], label[tr])
            cal[te] = cal_fn.predict(score[te])
        elif method == "platt":
            lr = LogisticRegression(C=1e6, solver="lbfgs", max_iter=1000)
            lr.fit(score[tr].reshape(-1, 1), label[tr].astype(int))
            cal[te] = lr.predict_proba(score[te].reshape(-1, 1))[:, 1]
    return cal


# ---------------------------------------------------------------------------
# gene-clustered bootstrap on a metric difference
# ---------------------------------------------------------------------------
def _diff_fn(prob_f, prob_s, label, fn, lower_is_better):
    def diff(idx):
        a, b = fn(prob_s[idx], label[idx]), fn(prob_f[idx], label[idx])
        return (a - b) if lower_is_better else (b - a)   # positive => fusion better
    return diff


def boot_diff(prob_f, prob_s, label, genes, fn, lower_is_better=True):
    """PRIMARY inference: gene-clustered percentile bootstrap (respects the
    non-independence of variants within a gene). This is the headline CI."""
    ug = pd.unique(genes)
    diff = _diff_fn(prob_f, prob_s, label, fn, lower_is_better)
    obs = diff(np.arange(len(label)))
    ds = [diff(np.concatenate([np.where(genes == g)[0]
                               for g in RNG.choice(ug, len(ug), replace=True)]))
          for _ in range(N_BOOT)]
    lo, hi = np.nanpercentile(ds, [2.5, 97.5])
    return {"obs": obs, "lo": lo, "hi": hi, "fusion_better": bool(lo > 0)}


def boot_diff_variant(prob_f, prob_s, label, fn, lower_is_better=True):
    """SENSITIVITY ONLY: variant-level bootstrap. Assumes variants are
    independent, which they are NOT (they cluster within genes), so this CI is
    ANTI-CONSERVATIVE (too narrow). Report only as a lower bound / lenient check,
    never as the primary result."""
    n = len(label)
    diff = _diff_fn(prob_f, prob_s, label, fn, lower_is_better)
    obs = diff(np.arange(n))
    ds = [diff(RNG.integers(0, n, n)) for _ in range(N_BOOT)]
    lo, hi = np.nanpercentile(ds, [2.5, 97.5])
    return {"obs": obs, "lo": lo, "hi": hi, "fusion_better": bool(lo > 0)}


def boot_diff_bca(prob_f, prob_s, label, genes, fn, lower_is_better=True):
    """SENSITIVITY: gene-clustered bootstrap with BCa (bias-corrected and
    accelerated) interval; acceleration from a leave-one-gene-out jackknife."""
    ug = pd.unique(genes)
    diff = _diff_fn(prob_f, prob_s, label, fn, lower_is_better)
    obs = diff(np.arange(len(label)))
    boots = np.array([diff(np.concatenate([np.where(genes == g)[0]
                      for g in RNG.choice(ug, len(ug), replace=True)]))
                      for _ in range(N_BOOT)])
    boots = boots[np.isfinite(boots)]
    if len(boots) < 100:
        return {"obs": obs, "lo": np.nan, "hi": np.nan, "fusion_better": False}
    # bias correction
    z0 = norm.ppf(np.clip(np.mean(boots < obs), 1e-6, 1 - 1e-6))
    # acceleration via leave-one-gene-out jackknife
    jack = np.array([diff(np.where(genes != g)[0]) for g in ug])
    jack = jack[np.isfinite(jack)]
    jbar = jack.mean()
    num = np.sum((jbar - jack) ** 3)
    den = 6.0 * (np.sum((jbar - jack) ** 2) ** 1.5)
    a = num / den if den != 0 else 0.0
    def adj(alpha):
        z = norm.ppf(alpha)
        return norm.cdf(z0 + (z0 + z) / (1 - a * (z0 + z)))
    lo = np.nanpercentile(boots, 100 * adj(0.025))
    hi = np.nanpercentile(boots, 100 * adj(0.975))
    return {"obs": obs, "lo": lo, "hi": hi, "fusion_better": bool(lo > 0)}


# ---------------------------------------------------------------------------
# one full evaluation for a given label column and gene set
# ---------------------------------------------------------------------------
def evaluate(df, scores: dict, label_col: str, tag: str, best_key: str, method="isotonic"):
    genes = df["gene"].astype(str).to_numpy()
    label = df[label_col].astype(float).to_numpy()

    cal = {name: logo_calibrate(s, label, genes, method=method) for name, s in scores.items()}
    summary = []
    for name, p in cal.items():
        summary.append({"set": tag, "calib": method, "model": name,
                        "ECE": round(ece(p, label), 4), "Brier": round(brier(p, label), 4),
                        **yield_metrics(p, label)})
    summ_df = pd.DataFrame(summary)

    fus, best = cal["fusion_M1"], cal[best_key]
    d_ece = boot_diff(fus, best, label, genes, ece, lower_is_better=True)
    d_bri = boot_diff(fus, best, label, genes, brier, lower_is_better=True)
    d_yld = boot_diff(fus, best, label, genes,
                      lambda p, y: yield_metrics(p, y)["actionable_frac"] or 0.0,
                      lower_is_better=False)
    # (three) -- second/third CI schemes on the Brier headline
    d_bri_var = boot_diff_variant(fus, best, label, brier, lower_is_better=True)
    d_bri_bca = boot_diff_bca(fus, best, label, genes, brier, lower_is_better=True)
    supported = [name for name, d in
                 [("Brier", d_bri), ("ECE", d_ece), ("yield", d_yld)] if d["fusion_better"]]
    head = {"set": tag, "calib": method,
            "dECE": round(d_ece["obs"], 4), "dECE_ci": f"[{d_ece['lo']:.4f},{d_ece['hi']:.4f}]",
            "dBrier": round(d_bri["obs"], 4),
            "dBrier_ci_geneclust": f"[{d_bri['lo']:.4f},{d_bri['hi']:.4f}]",
            "dBrier_ci_BCa": f"[{d_bri_bca['lo']:.4f},{d_bri_bca['hi']:.4f}]",
            "dBrier_ci_variant(sens)": f"[{d_bri_var['lo']:.4f},{d_bri_var['hi']:.4f}]",
            "dYield": round(d_yld["obs"], 4), "dYield_ci": f"[{d_yld['lo']:.4f},{d_yld['hi']:.4f}]",
            "fusion_better_on": ",".join(supported) if supported else "none"}

    # (three) -- per-tool Brier delta CI: does fusion beat EVERY single tool?
    per_tool = []
    for name, p in cal.items():
        if not name.startswith("single:"):
            continue
        d = boot_diff(fus, p, label, genes, brier, lower_is_better=True)
        per_tool.append({"set": tag, "calib": method, "tool": name.replace("single:", ""),
                         "dBrier_vs_fusion": round(d["obs"], 4),
                         "ci95": f"[{d['lo']:.4f},{d['hi']:.4f}]",
                         "fusion_lower": d["fusion_better"]})
    per_tool_df = pd.DataFrame(per_tool)
    return summ_df, head, cal, per_tool_df


# ---------------------------------------------------------------------------
# orchestration
# ---------------------------------------------------------------------------
def run():
    p = C.OUTPUT_DIR / f"frozen_matrix_{C.FROZEN_VERSION}.parquet"
    if not p.exists():
        sys.exit(f"[phase3] frozen matrix not found at {p}")
    df = pd.read_parquet(p)
    df = df[df["is_splice"] & df["func_pathogenicity"].notna()].reset_index(drop=True)

    feats = usable_features(df)
    genes = df["gene"].astype(str)
    Xrn = rank_within_gene(df, feats)
    ytr = rank_target_within_gene(df) if getattr(C, "TRAIN_TARGET", "rank_within_gene") == "rank_within_gene" \
          else df["func_pathogenicity"]

    # best single by the SAME criterion as Phase 2 (DL-pooled per-gene Spearman),
    # so the H1 and H3 comparators are identical (avoids picking a weaker tool).
    from .phase2_model import per_gene_rho, dl_pool
    yraw = df["func_pathogenicity"].to_numpy(); gv = genes.to_numpy()
    single_pool = {f: dl_pool(per_gene_rho(Xrn[f].to_numpy(), yraw, gv))["rho"] for f in feats}
    best = max(single_pool, key=lambda k: single_pool[k] if np.isfinite(single_pool[k]) else -9)

    oof_enet = logo_oof(Xrn, ytr, genes, feats, "enet")
    oof_mean = logo_oof(Xrn, ytr, genes, feats, "mean")
    # calibrate EVERY single tool too, so no one can claim the comparator was cherry-picked
    scores = {"fusion_M1": oof_enet, "mean_M0b": oof_mean}
    for f in feats:
        scores[f"single:{f}"] = Xrn[f].to_numpy()
    best_key = f"single:{best}"
    # min-max each score to [0,1] (monotonic, harmless) for a common isotonic support
    for k in scores:
        s = scores[k]; lo, hi = np.nanmin(s), np.nanmax(s)
        scores[k] = (s - lo) / (hi - lo) if hi > lo else s

    C.REPORT_DIR.mkdir(parents=True, exist_ok=True)
    all_summary, all_head, all_pertool = [], [], []
    for method in CAL_METHODS:                       # (one) isotonic AND platt
        for label_col in ["y_assay", "y_clinvar"]:
            if label_col not in df.columns or df[label_col].notna().sum() < 50:
                continue
            for brca in ["included", "excluded"]:
                sub = df if brca == "included" else df[df["gene"] != "BRCA1"].reset_index(drop=True)
                idx = df.index if brca == "included" else df[df["gene"] != "BRCA1"].index
                sc = {k: v[idx.to_numpy()] for k, v in scores.items()}
                tag = f"{label_col}/BRCA1_{brca}"
                summ, head, cal, ptool = evaluate(sub, sc, label_col, tag, best_key, method=method)
                all_summary.append(summ); all_head.append(head); all_pertool.append(ptool)
                if method == "isotonic" and label_col == "y_assay" and brca == "included":
                    reliability(cal["fusion_M1"], sub[label_col].astype(float).to_numpy()
                                ).to_csv(C.REPORT_DIR / "phase3_reliability_fusion.csv", index=False)

    summ_df = pd.concat(all_summary, ignore_index=True)
    head_df = pd.DataFrame(all_head)
    pertool_df = pd.concat(all_pertool, ignore_index=True)
    summ_df.to_csv(C.REPORT_DIR / "phase3_calibration_summary.csv", index=False)
    head_df.to_csv(C.REPORT_DIR / "phase3_H3_headline.csv", index=False)
    pertool_df.to_csv(C.REPORT_DIR / "phase3_pertool_brier_ci.csv", index=False)

    print(f"[phase3] best single (pooled per-gene Spearman) = {best} | "
          f"labelled splice (y_assay) = {int(df['y_assay'].notna().sum())}/{len(df)}\n")
    print("=== per-model calibration (ALL objects identically calibrated) ===")
    print(summ_df.to_string(index=False))
    print(f"\n=== H3 headline: fusion vs best single ({best}) -- 3 CI schemes, 2 calibrators ===")
    print(head_df.to_string(index=False))
    print("\n=== per-tool: does fusion's Brier beat EVERY single tool? (gene-clustered CI) ===")
    print(pertool_df.to_string(index=False))

    # robustness one-liners
    bri_ok = head_df.apply(lambda r: r["dBrier"] > 0 and "[-" not in r["dBrier_ci_geneclust"].split(",")[0], axis=1)
    print("\n[phase3] ROBUSTNESS SUMMARY")
    print(f"  Brier favours fusion with gene-clustered CI excluding 0 in "
          f"{int(bri_ok.sum())}/{len(head_df)} (calibrator x label x BRCA1) conditions.")
    print("  Primary CI = gene-clustered; BCa = sensitivity; variant-level = anti-conservative lower bound only.")
    print("  per-tool table shows whether fusion's Brier is significantly below EACH single tool.")


if __name__ == "__main__":
    run()
