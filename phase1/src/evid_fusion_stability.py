"""Fusion stability -- are the fused model's weights stable across the leave-one-
gene-out folds, and does the fusion reach a tier that no single column reaches?

Are the elastic-net fusion's weights stable? The fusion is
refitted inside every leave-one-gene-out fold, so there is one set of weights per
held-out gene, and E3 reports only the out-of-fold score those fits produce. This
stage refits exactly the same fusion -- the same columns (evid_common.FUSION_FEATURES,
the panel without the AlphaGenome splice score, in the order E3 passes them), the same rows, the same within-gene rank target and
the same call into phase2_model.logo_oof that evid_common.logo_fusion makes -- with
the per-fold coefficients kept, and writes them fold by fold beside their spread.
The cross-fold mean alone cannot answer the question: a weight that swings between
+0.4 and -0.4 has the same mean magnitude as one that sits at 0.4 in every fold,
which is the point phase6_sensitivity made for the published fusion.

The refit is checked rather than assumed to be the same object. Its out-of-fold
score must equal logo_fusion's bit for bit, and it must reproduce the row counts E3
wrote for the fusion and the held-out likelihood ratio of every fusion row in
evidence_thresholds_logo_folds.csv. If either check fails the stage stops. A
coefficient table from a different fit would be worse than no table, and so would a
Table S9 that joins a fresh fit to stale E3 output.

What stability across these folds does and does not mean. Any two folds share all
but one of their training genes, so the fits are jackknife replicates, not
independent ones. Agreement says that removing any single gene does not move the
weights much; it does not say that the weights would come back on a new set of
genes. The features are within-gene ranks on the unit interval and the target is a
within-gene rank, so a coefficient is on the same scale in every fold and comparing
folds is meaningful. The model is fitted on the whole analysis set, canonical-site
stratum included, exactly as E3's is: there is one set of weights per fold, not one
per stratum. Each fold's weights are fitted without the held-out gene, so the
coefficients themselves carry no leakage from it.

"Non-zero" means exactly non-zero. The elastic net zeroes a weight by soft-
thresholding, so a zero is a decision of the fit and not a rounding artefact. A
term counts as changing sign only when one fold gives it a strictly positive weight
and another a strictly negative one. A drop to zero is counted separately, because
it is a different kind of instability. Negative zeros returned by the solver are
written as zero. The standard deviation is taken across folds with ddof=1, as in
phase6. Terms are ranked by absolute weight within each fold, and the best and worst
rank across folds are reported, because a reader asking about stability usually
means the ordering as well as the values. Ties take the worse rank, so a zeroed
weight ranks last. The fit draws no random numbers (cyclic coordinate descent and an
unshuffled inner cross-validation), so there is nothing to seed beyond the
random_state phase2_model already fixes.

The design also carries one missingness indicator per panel column
(phase2_model._en_design). Where a panel column has no missing value in the rows
used, its indicator is a constant zero and the elastic net cannot weight it. The
coefficient file keeps every term of the fitted design, with the count of missing
values behind each one. The print block drops only indicator terms whose weight is
zero in every fold.

Table S9 has three blocks, told apart by its `block` column:

  a_insample_tier  per in-scope stratum, the highest PP3 tier the fusion reaches
                   in-sample (E3's pp3_threshold_reachable) against the highest
                   tier any single column reaches, and which columns reach it. Two
                   comparison sets are given. The fusion's own inputs are the fair
                   test of whether combining them adds a tier; every single column
                   E3 evaluated also includes those carried beside the panel.
  b_logo_fold      per in-scope stratum and Supporting/Moderate/Strong: the fusion
                   threshold fitted on the stratum's other genes and the likelihood
                   ratio of the band above it in the held-out gene. The held-out
                   counts explain every missing ratio. Values are kept at E3's
                   precision, so each row can be matched to its source by
                   equality. E3 names its column threshold_from_6_genes, but a
                   stratum in which a gene has no variant fits on one gene fewer,
                   so the number of training genes is counted and written beside
                   it.
  c_coefficients   the coefficient summary, rounded for print; the full-precision
                   values are in fusion_fold_coefficients.csv.

The leakage block b carries, stated here as it is in evid_interval_lr's docstring.
logo_fusion produces ONE out-of-fold column. Gene g's score comes from a model
trained on the other genes, and gene h's score comes from a model whose training
set INCLUDED g. When E3 holds g out and fits a threshold on the other genes'
scores, g's labels have therefore already shaped the scores that threshold is
fitted on. The threshold is not free of the held-out gene, and the fusion's held-out
likelihood ratios are optimistic. What block b shows is that the fusion's fitted
thresholds can seldom be tested in the held-out gene: across the in-scope strata
and the three tiers, a held-out ratio exists in fewer than half of the folds, and
where one exists it nearly always clears the tier's cut. The contamination can only
make those ratios look better than they are.

Stability is measured by leaving one gene out at a time, so any two folds share all
but one training gene. "Stable" therefore means that removing a single gene barely
moves the weights, not that the weights would reappear on new genes. The inner
cross-validation that picks the penalty splits rows in their stored order, so the
selected penalty, and with it which small weights reach exactly zero, depends on the
row order of the analysis set, which is fixed by its manifest. A clean estimate
would need a nested design, with the fusion refitted inside each outer fold before
any threshold is fitted, and that design is not run here. Every block b row carries
E3's `logo_threshold_leakage` flag.

Block a has a related caveat. The fusion's in-sample tier is read off a curve that
pools out-of-fold scores from differently fitted models, each on its own scale.

Run (PYTHONPATH=phase1):  python -m src.evid_fusion_stability
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from . import evid_common as K
from .evid_interval_lr import TIERS

REPORT_DIR = K.REPORT_DIR
SUPP_DIR = REPORT_DIR / "supplement"
COEF_OUT = REPORT_DIR / "fusion_fold_coefficients.csv"
TABLE_OUT = SUPP_DIR / "tableS9_fusion.csv"
THRESHOLDS = "evidence_thresholds.csv"
LOGO_FOLDS = "evidence_thresholds_logo_folds.csv"

STRATA = [s for s in K.STRATA if s in K.IN_SCOPE_STRATA]
FOLD_TIERS = ["supporting", "moderate", "strong"]
PRINT_DECIMALS = 3
FOLD_PREFIX = "coef_heldout_"
NO_TIER = "none"

COLS_A = ["stratum", "comparison_set", "n_single_columns_compared",
          "fusion_insample_tier", "best_single_insample_tier",
          "single_columns_reaching_best_tier", "fusion_vs_best_single"]
COLS_B = ["stratum", "tier", "tier_lr_cut", "heldout_gene", "n_training_genes",
          "threshold_from_training_genes", "heldout_lr", "heldout_lr_at_or_above_cut",
          "heldout_n", "heldout_n_damaging", "heldout_n_normal", "heldout_n_in_band",
          "why_no_heldout_lr", "logo_threshold_leakage"]
COLS_C = ["term", "term_kind", "n_folds", "mean", "sd", "min", "max", "mean_abs",
          "rank_mean_abs", "n_folds_nonzero", "n_folds_nonzero_same_sign_as_mean",
          "changes_sign_across_folds", "rank_in_fold_best", "rank_in_fold_worst"]
INT_COLS = ["n_single_columns_compared", "n_training_genes", "heldout_n",
            "heldout_n_damaging", "heldout_n_normal", "heldout_n_in_band",
            "n_folds", "rank_mean_abs", "n_folds_nonzero",
            "n_folds_nonzero_same_sign_as_mean",
            "rank_in_fold_best", "rank_in_fold_worst", "n_missing_in_rows_used",
            "n_folds_positive", "n_folds_negative", "n_folds_zero"]
BOOL_COLS = ["heldout_lr_at_or_above_cut", "changes_sign_across_folds"]


def fusion_features(df: pd.DataFrame) -> list[str]:
    """The fusion's inputs, chosen exactly as evid_interval_lr.main chooses them."""
    return [t for t in K.FUSION_FEATURES if t in df.columns]


def read_e3(report_dir: Path, name: str) -> pd.DataFrame:
    # round_trip so a stored threshold equals the score it was read off, bit for bit
    return pd.read_csv(Path(report_dir) / name, float_precision="round_trip")


# One set of writer options for the files and for csv_text, so the text a test
# compares is the text main() writes; the line ending is fixed rather than left
# to the platform.
CSV_OPTS = {"index": False, "lineterminator": "\n"}


def csv_text(d: pd.DataFrame) -> str:
    """The exact bytes written, so a test can compare a rerun with the file."""
    return d.to_csv(**CSV_OPTS)


# ---------------------------------------------------------------------------
# the refit
# ---------------------------------------------------------------------------
def rows_used(df: pd.DataFrame) -> np.ndarray:
    """The rows logo_fusion fits on: those with a within-gene rank target."""
    from .phase2_model import rank_target_within_gene
    y = rank_target_within_gene(df, "func_pathogenicity")
    return y.reset_index(drop=True).notna().to_numpy()


def refit_with_coefficients(df: pd.DataFrame,
                            feats: list[str]) -> tuple[np.ndarray, pd.DataFrame]:
    """evid_common.logo_fusion, line for line, with collect_coefs=True.

    The body is repeated rather than imported because logo_fusion does not pass
    collect_coefs through. That repetition is why build() compares the two scores
    bit for bit before anything is written.
    """
    from .phase2_model import logo_oof, rank_target_within_gene
    xrn = K.rank_within_gene(df, feats)
    y = rank_target_within_gene(df, "func_pathogenicity")
    genes = df["gene"].astype(str).reset_index(drop=True)
    xrn = xrn.reset_index(drop=True)
    keep = y.reset_index(drop=True).notna().to_numpy()
    oof = np.full(len(df), np.nan)
    sub, coefs = logo_oof(xrn.loc[keep].reset_index(drop=True),
                          y.reset_index(drop=True).loc[keep].reset_index(drop=True),
                          genes.loc[keep].reset_index(drop=True), feats, "enet",
                          collect_coefs=True)
    oof[keep] = sub
    return oof, coefs


def coefficient_table(coefs: pd.DataFrame, feats: list[str],
                      n_missing: pd.Series) -> pd.DataFrame:
    """Term x held-out gene, then the cross-fold summary, in design order."""
    per_fold = coefs.astype(float) + 0.0          # -0.0 -> 0.0; values unchanged
    folds = [str(g) for g in per_fold.columns]
    v = per_fold.to_numpy()
    mean = v.mean(axis=1)
    pos, neg, zero = (v > 0).sum(axis=1), (v < 0).sum(axis=1), (v == 0).sum(axis=1)
    same = ((np.sign(v) == np.sign(mean)[:, None]) & (v != 0)).sum(axis=1)
    # Rank of |coef| within each fold, 1 = largest. Ties take the WORSE rank: in
    # practice the only ties are exact zeros, and a zeroed weight should rank last
    # rather than just behind the last non-zero one, where it would look placed.
    kind, column = [], []
    for term in per_fold.index:
        base = term[:-len("_isna")] if term.endswith("_isna") else term
        is_ind = term.endswith("_isna") and base in feats
        kind.append("missingness indicator" if is_ind else "within-gene rank")
        column.append(base if is_ind else term)
    # ranks are taken among the within-gene rank terms only; the indicators are
    # given no rank
    is_rank_term = pd.Series([k == "within-gene rank" for k in kind], index=per_fold.index)
    fold_rank = (per_fold[is_rank_term].abs()
                 .rank(axis=0, ascending=False, method="max")
                 .reindex(per_fold.index))

    out = pd.DataFrame({
        "term": list(per_fold.index),
        "panel_column": column,
        "term_kind": kind,
        "n_missing_in_rows_used": [int(n_missing[c]) for c in column],
    })
    for j, g in enumerate(folds):
        out[f"{FOLD_PREFIX}{g}"] = v[:, j]
    out["n_folds"] = v.shape[1]
    out["mean"] = mean + 0.0
    out["sd"] = v.std(axis=1, ddof=1)
    out["min"] = v.min(axis=1)
    out["max"] = v.max(axis=1)
    out["mean_abs"] = np.abs(v).mean(axis=1)
    out["rank_mean_abs"] = (out["mean_abs"].where(is_rank_term.to_numpy())
                            .rank(ascending=False, method="max").astype("Int64"))
    out["n_folds_positive"] = pos
    out["n_folds_negative"] = neg
    out["n_folds_zero"] = zero
    out["n_folds_nonzero"] = pos + neg
    out["n_folds_nonzero_same_sign_as_mean"] = same
    out["changes_sign_across_folds"] = (pos > 0) & (neg > 0)
    out["rank_in_fold_best"] = fold_rank.min(axis=1).astype("Int64").to_numpy()
    out["rank_in_fold_worst"] = fold_rank.max(axis=1).astype("Int64").to_numpy()
    return out


# ---------------------------------------------------------------------------
# agreement with E3
# ---------------------------------------------------------------------------
def _scored(df: pd.DataFrame, stratum: str) -> pd.DataFrame:
    """The rows E3 evaluates the fusion on in a stratum (run_tool_stratum)."""
    return K.stratum_frame(df, stratum).dropna(subset=["y_assay", K.FUSION])


def check_e3_is_current(df: pd.DataFrame, ev: pd.DataFrame,
                        folds: pd.DataFrame) -> None:
    """Stop if E3's fusion outputs were not produced by this fit on this set."""
    bad = []
    fe = ev[ev["tool"] == K.FUSION]
    if not len(fe):
        bad.append(f"{THRESHOLDS} has no {K.FUSION} rows")
    for stratum in K.STRATA:
        sub = _scored(df, stratum)
        y = sub["y_assay"].to_numpy(dtype=float)
        want = {"n": len(sub), "n_pos": int((y == 1).sum()), "n_neg": int((y == 0).sum())}
        for r in fe[fe["stratum"] == stratum].itertuples():
            got = {k: int(getattr(r, k)) for k in want}
            if got != want:
                bad.append(f"{THRESHOLDS} {stratum}/{r.tier}: {got} != {want}")
    ff = folds[folds["tool"] == K.FUSION]
    if not len(ff):
        bad.append(f"{LOGO_FOLDS} has no {K.FUSION} rows")
    for r in ff.itertuples():
        sub = _scored(df, r.stratum)
        te = (sub["gene"].astype(str) == r.heldout_gene).to_numpy()
        lr = K.band_lr(sub["y_assay"].to_numpy(dtype=float)[te],
                       sub[K.FUSION].to_numpy(dtype=float)[te],
                       r.threshold_from_6_genes, "upper")
        stored = r.heldout_lr
        if not ((np.isnan(lr) and np.isnan(stored)) or lr == stored):
            bad.append(f"{LOGO_FOLDS} {r.stratum}/{r.tier}/{r.heldout_gene}: "
                       f"stored {stored!r}, refit gives {lr!r}")
    if bad:
        raise SystemExit(
            "[fusion-stability] E3's fusion outputs do not match the current fit "
            f"({len(bad)} disagreements); rerun python -m src.evid_interval_lr.\n  "
            + "\n  ".join(bad[:20]))


# ---------------------------------------------------------------------------
# Table S9
# ---------------------------------------------------------------------------
def _reached(ev: pd.DataFrame) -> pd.Series:
    """tool x stratum -> highest PP3 tier reached in-sample, as a rank (-1 = none)."""
    ok = ev[ev["status"] == "ok"].copy()
    ok["reached"] = ok["pp3_threshold_reachable"].map(lambda v: str(v) == "True")
    ok["rank"] = ok["tier"].map({t: i for i, t in enumerate(TIERS)})
    hit = ok[ok["reached"]].groupby(["tool", "stratum"])["rank"].max()
    idx = pd.MultiIndex.from_product([sorted(ev["tool"].unique()),
                                      sorted(ev["stratum"].unique())],
                                     names=["tool", "stratum"])
    return hit.reindex(idx, fill_value=-1).astype(int)


def _tier_name(rank: int) -> str:
    return TIERS[rank] if rank >= 0 else NO_TIER


def block_a(ev: pd.DataFrame) -> pd.DataFrame:
    reached = _reached(ev)
    tools = list(ev["tool"].unique())
    singles = ([t for t in K.PANEL + K.OPTIONAL if t in tools]
               + sorted(t for t in tools
                        if t not in K.PANEL + K.OPTIONAL and t != K.FUSION))
    sets = {"fusion inputs": [t for t in singles if t in K.FUSION_FEATURES],
            "all single columns in E3": singles}
    rows = []
    for stratum in STRATA:
        f = int(reached[(K.FUSION, stratum)])
        for name, cols in sets.items():
            ranks = {c: int(reached[(c, stratum)]) for c in cols}
            best = max(ranks.values())
            rows.append({
                "stratum": stratum, "comparison_set": name,
                "n_single_columns_compared": len(cols),
                "fusion_insample_tier": _tier_name(f),
                "best_single_insample_tier": _tier_name(best),
                "single_columns_reaching_best_tier": (
                    "; ".join(c for c in cols if ranks[c] == best) if best >= 0 else ""),
                "fusion_vs_best_single": ("fusion higher" if f > best else
                                          "same tier" if f == best else "fusion lower"),
            })
    return pd.DataFrame(rows, columns=COLS_A)


def block_b(df: pd.DataFrame, ev: pd.DataFrame, folds: pd.DataFrame) -> pd.DataFrame:
    fe = ev[ev["tool"] == K.FUSION].set_index(["stratum", "tier"])
    f = folds[(folds["tool"] == K.FUSION) & (folds["side"] == "pp3")
              & folds["stratum"].isin(STRATA) & folds["tier"].isin(FOLD_TIERS)].copy()
    f["_s"] = f["stratum"].map({s: i for i, s in enumerate(STRATA)})
    f["_t"] = f["tier"].map({t: i for i, t in enumerate(FOLD_TIERS)})
    f = f.sort_values(["_s", "_t"], kind="mergesort")
    rows = []
    for r in f.itertuples():
        sub = _scored(df, r.stratum)
        g = sub["gene"].astype(str)
        te = (g == r.heldout_gene).to_numpy()
        y = sub["y_assay"].to_numpy(dtype=float)[te]
        s = sub[K.FUSION].to_numpy(dtype=float)[te]
        tau, lr = r.threshold_from_6_genes, r.heldout_lr
        n_pos, n_neg = int((y == 1).sum()), int((y == 0).sum())
        n_band = int((s >= tau).sum()) if np.isfinite(tau) else pd.NA
        # the reasons mirror evid_common.band_lr's own refusals, in order
        if not np.isfinite(tau):
            why = "tier not reached on the training genes"
        elif n_pos < K.MIN_POS or n_neg < K.MIN_NEG:
            why = (f"held-out gene has fewer than {K.MIN_POS} damaging or "
                   f"{K.MIN_NEG} normal variants in this stratum")
        elif n_band < K.MIN_BAND:
            why = (f"fewer than {K.MIN_BAND} held-out variants score at or "
                   f"above the threshold")
        else:
            why = ""
        e = fe.loc[(r.stratum, r.tier)]
        rows.append({
            "stratum": r.stratum, "tier": r.tier, "tier_lr_cut": e["path_cut_lr"],
            "heldout_gene": r.heldout_gene,
            "n_training_genes": int(g[~te].nunique()),
            "threshold_from_training_genes": tau,
            "heldout_lr": lr,
            "heldout_lr_at_or_above_cut": (bool(lr >= e["path_cut_lr"])
                                           if np.isfinite(lr) else pd.NA),
            "heldout_n": int(te.sum()), "heldout_n_damaging": n_pos,
            "heldout_n_normal": n_neg, "heldout_n_in_band": n_band,
            "why_no_heldout_lr": why,
            "logo_threshold_leakage": e["logo_threshold_leakage"],
        })
    return pd.DataFrame(rows, columns=COLS_B)


def block_c(coef: pd.DataFrame) -> pd.DataFrame:
    shown = coef[~((coef["term_kind"] == "missingness indicator")
                   & (coef["n_folds_nonzero"] == 0))]
    shown = shown.sort_values("rank_mean_abs", kind="mergesort")
    out = shown[COLS_C].copy()
    for c in ("mean", "sd", "min", "max", "mean_abs"):
        # + 0.0 so a small negative mean does not print as -0.0
        out[c] = np.round(out[c].to_numpy(dtype=float), PRINT_DECIMALS) + 0.0
    return out


def _typed(d: pd.DataFrame) -> pd.DataFrame:
    """Integer and boolean columns keep their type through the block union, so a
    count prints as 7 and not 7.0."""
    for c in d.columns:
        if c in INT_COLS:
            d[c] = d[c].astype("Int64")
        elif c in BOOL_COLS:
            d[c] = d[c].astype("boolean")
    return d


def table_s9(a: pd.DataFrame, b: pd.DataFrame, c: pd.DataFrame) -> pd.DataFrame:
    cols = ["block"] + list(dict.fromkeys(COLS_A + COLS_B + COLS_C))
    parts = [p.assign(block=name) for name, p in
             (("a_insample_tier", a), ("b_logo_fold", b), ("c_coefficients", c))]
    return _typed(pd.concat(parts, ignore_index=True).reindex(columns=cols))


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------
def build(df: pd.DataFrame, report_dir: Path | str = REPORT_DIR
          ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(fusion_fold_coefficients, tableS9_fusion) from the set and E3's outputs."""
    feats = fusion_features(df)
    oof, coefs = refit_with_coefficients(df, feats)
    ref = K.logo_fusion(df, feats)
    if not np.array_equal(oof, ref, equal_nan=True):
        raise SystemExit("[fusion-stability] the refit's out-of-fold score differs "
                         "from evid_common.logo_fusion; the two code paths have "
                         "drifted apart")
    df = df.copy()
    df[K.FUSION] = oof

    ev = read_e3(report_dir, THRESHOLDS)
    folds = read_e3(report_dir, LOGO_FOLDS)
    check_e3_is_current(df, ev, folds)

    n_missing = df.loc[rows_used(df), feats].isna().sum()
    coef = _typed(coefficient_table(coefs, feats, n_missing))
    return coef, table_s9(block_a(ev), block_b(df, ev, folds), block_c(coef))


# Printed headers for Table S9; build() keeps the internal names the tests use.
S9_HEADERS = {
    "block": "Block", "stratum": "Stratum", "comparison_set": "Comparison set",
    "n_single_columns_compared": "Single columns compared",
    "fusion_insample_tier": "Fusion, in-sample tier",
    "best_single_insample_tier": "Best single column, in-sample tier",
    "single_columns_reaching_best_tier": "Single columns reaching that tier",
    "fusion_vs_best_single": "Fusion against the best single column",
    "tier": "Tier", "tier_lr_cut": "Tier boundary (LR)", "heldout_gene": "Held-out gene",
    "n_training_genes": "Training genes",
    "threshold_from_training_genes": "Threshold from the training genes",
    "heldout_lr": "Held-out LR", "heldout_lr_at_or_above_cut": "Held-out LR clears the boundary",
    "heldout_n": "Held-out variants", "heldout_n_damaging": "Held-out damaging",
    "heldout_n_normal": "Held-out normal", "heldout_n_in_band": "Held-out variants above the threshold",
    "why_no_heldout_lr": "Why no held-out LR",
    "logo_threshold_leakage": "Leakage into the threshold",
    "term": "Term", "term_kind": "Term kind", "n_folds": "Folds",
    "mean": "Mean coefficient", "sd": "SD", "min": "Minimum", "max": "Maximum",
    "mean_abs": "Mean absolute coefficient", "rank_mean_abs": "Rank by mean absolute value",
    "n_folds_nonzero": "Folds with a non-zero coefficient",
    "n_folds_nonzero_same_sign_as_mean": "Non-zero folds with the sign of the mean",
    "changes_sign_across_folds": "Changes sign across folds",
    "rank_in_fold_best": "Best rank in a fold", "rank_in_fold_worst": "Worst rank in a fold",
}


def printed_s9(s9: pd.DataFrame) -> pd.DataFrame:
    """Table S9 as printed: readable headers, and Yes/No for the two flags."""
    out = s9.copy()
    for c in BOOL_COLS:
        out[c] = out[c].map({True: "Yes", False: "No"}).astype("string")
    return out.rename(columns=S9_HEADERS)


def main() -> None:
    df = K.load_set()
    coef, s9 = build(df)
    missing = [c for c in s9.columns if c not in S9_HEADERS]
    if missing:
        raise SystemExit(f"[fusion-stability] no printed header for {missing}")
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    SUPP_DIR.mkdir(parents=True, exist_ok=True)
    # literal names in the to_csv calls, the form the reproduction-coverage check
    # looks for when it works back from an output to the module that writes it
    coef.to_csv(REPORT_DIR / "fusion_fold_coefficients.csv", **CSV_OPTS)
    printed_s9(s9).to_csv(SUPP_DIR / "tableS9_fusion.csv", **CSV_OPTS)
    n_folds = int(coef["n_folds"].iloc[0])
    print(f"[fusion-stability] refit matches logo_fusion and E3; {n_folds} folds, "
          f"{len(coef)} design terms")
    print(f"[fusion-stability] wrote {COEF_OUT} and {TABLE_OUT} ({len(s9)} rows)\n")

    c = s9[s9["block"] == "c_coefficients"]
    print("--- per-fold weights, by mean |coefficient| ---")
    print(c[["term", "mean", "sd", "min", "max", "mean_abs", "n_folds_nonzero",
             "n_folds_nonzero_same_sign_as_mean", "changes_sign_across_folds",
             "rank_in_fold_best", "rank_in_fold_worst"]].to_string(index=False))
    flips = coef.loc[coef["changes_sign_across_folds"].astype(bool), "term"].tolist()
    print(f"\nterms changing sign across folds: {', '.join(flips) if flips else 'none'}")

    a = s9[s9["block"] == "a_insample_tier"]
    print("\n--- highest PP3 tier reached in-sample: fusion against single columns ---")
    print(a[["stratum", "comparison_set", "fusion_insample_tier",
             "best_single_insample_tier", "fusion_vs_best_single",
             "single_columns_reaching_best_tier"]].to_string(index=False))

    b = s9[s9["block"] == "b_logo_fold"]
    print("\n--- fusion threshold from the other genes, ratio in the held-out one ---")
    print(b.groupby(["stratum", "tier"], sort=False).agg(
        folds=("heldout_gene", "size"),
        with_threshold=("threshold_from_training_genes", "count"),
        with_heldout_lr=("heldout_lr", "count"),
        at_or_above_cut=("heldout_lr_at_or_above_cut", "sum"),
        lr_min=("heldout_lr", "min"), lr_max=("heldout_lr", "max"),
    ).round(3).to_string())


if __name__ == "__main__":
    main()
