"""Guardrails on the fusion-stability stage (src/evid_fusion_stability.py).

The stage refits the fusion to recover its per-fold weights, so the first thing to
pin is that the refit is the fusion E3 used and not a near relative. That is shown
three ways: bit for bit against evid_common.logo_fusion, against the score grid E3
wrote, and against E3's stored held-out ratios. The guard that enforces the last of
these must actually fire, so a perturbed copy is fed to it. The weight summary is
checked against an independent code path (phase6_sensitivity.weight_stability).
Table S9 is checked against the E3 tables it is read from. Finally, the outputs
must be identical on a rerun, identical to the files on disk, and independent of
the BLAS thread count.
"""
from __future__ import annotations

import hashlib
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

PHASE1 = Path(__file__).resolve().parents[1]
if str(PHASE1) not in sys.path:
    sys.path.insert(0, str(PHASE1))

from src import evid_common as K  # noqa: E402
from src import evid_fusion_stability as M  # noqa: E402
from src.evid_interval_lr import TIERS  # noqa: E402

REPORTS = PHASE1 / "reports/evidence"
SET = PHASE1 / "data/evidence/analysis_set_v1.parquet"
E3_FILES = [REPORTS / M.THRESHOLDS, REPORTS / M.LOGO_FOLDS]

pytestmark = pytest.mark.skipif(
    not SET.exists() or not all(p.exists() for p in E3_FILES),
    reason="analysis set or E3 outputs not built")


@pytest.fixture(scope="module")
def df():
    return pd.read_parquet(SET)


@pytest.fixture(scope="module")
def refit(df):
    feats = M.fusion_features(df)
    oof, coefs = M.refit_with_coefficients(df, feats)
    return feats, oof, coefs


@pytest.fixture(scope="module")
def built(df):
    return M.build(df, REPORTS)


@pytest.fixture(scope="module")
def e3():
    return M.read_e3(REPORTS, M.THRESHOLDS), M.read_e3(REPORTS, M.LOGO_FOLDS)


def _block(s9, name):
    return s9[s9["block"] == name].reset_index(drop=True)


# ---------------------------------------------------------------- the refit
def test_the_refit_is_logo_fusion_bit_for_bit(df, refit):
    feats, oof, coefs = refit
    assert feats == [t for t in K.PANEL if t in df.columns]
    ref = K.logo_fusion(df, feats)
    assert np.array_equal(oof, ref, equal_nan=True)
    assert np.isfinite(oof).any()
    # one rank term and one missingness indicator per feature, one column per gene
    genes = list(pd.unique(df["gene"].astype(str)))
    assert list(coefs.columns) == genes
    assert list(coefs.index) == feats + [f"{f}_isna" for f in feats]


def test_the_refit_reproduces_the_fusion_score_grid_e3_wrote(df, refit):
    """E3's curve files list every distinct fusion score in a stratum. The refit
    must land on exactly those values, not on values close to them."""
    _, oof, _ = refit
    d = df.assign(**{K.FUSION: oof})
    seen = 0
    for stratum in K.STRATA:
        path = REPORTS / f"interval_lr_{K.FUSION}_{stratum}.csv"
        if not path.exists():
            continue
        grid = pd.read_csv(path, float_precision="round_trip")["score"].to_numpy()
        sub = K.stratum_frame(d, stratum).dropna(subset=["y_assay", K.FUSION])
        assert np.array_equal(np.unique(sub[K.FUSION].to_numpy()), grid), stratum
        seen += 1
    if not seen:
        pytest.skip("E3 curve files not present in this checkout")


def test_the_staleness_guard_fires_on_a_changed_heldout_ratio(df, refit, e3):
    _, oof, _ = refit
    ev, folds = e3
    d = df.assign(**{K.FUSION: oof})
    M.check_e3_is_current(d, ev, folds)            # the real tables pass
    bad = folds.copy()
    i = bad.index[(bad["tool"] == K.FUSION) & bad["heldout_lr"].notna()][0]
    bad.loc[i, "heldout_lr"] = bad.loc[i, "heldout_lr"] * (1 + 1e-9)
    with pytest.raises(SystemExit):
        M.check_e3_is_current(d, ev, bad)
    worse = ev.copy()
    worse.loc[worse["tool"] == K.FUSION, "n_pos"] += 1
    with pytest.raises(SystemExit):
        M.check_e3_is_current(d, worse, folds)


# ------------------------------------------------------ the coefficient file
def test_the_weight_summary_agrees_with_phase6s_independent_summary(df, refit, built):
    """phase6_sensitivity.weight_stability summarises logo_oof's weights by its own
    code. Fed the same design, it must give the same summary numbers."""
    from src.phase2_model import rank_target_within_gene
    from src.phase6_sensitivity import weight_stability
    feats = refit[0]
    _, p6 = weight_stability(K.rank_within_gene(df, feats),
                             rank_target_within_gene(df, "func_pathogenicity"),
                             df["gene"].astype(str), feats)
    coef = built[0].set_index("term")
    p6 = p6.loc[coef.index]
    assert np.allclose(coef["mean"], p6["mean"], rtol=0, atol=1e-15)
    assert np.allclose(coef["sd"], p6["sd"], rtol=0, atol=1e-15)
    assert (coef["min"].to_numpy() == p6["min"].to_numpy()).all()
    assert (coef["max"].to_numpy() == p6["max"].to_numpy()).all()
    assert (coef["n_folds_nonzero"].to_numpy() == p6["n_folds_nonzero"].to_numpy()).all()
    assert (coef["changes_sign_across_folds"].astype(int).to_numpy()
            == p6["n_folds_sign_flip"].to_numpy()).all()


def test_the_fold_counts_add_up(df, built):
    coef = built[0]
    fold_cols = [c for c in coef.columns if c.startswith(M.FOLD_PREFIX)]
    assert len(fold_cols) == df["gene"].nunique()
    assert (coef["n_folds"] == len(fold_cols)).all()
    assert (coef["n_folds_positive"] + coef["n_folds_negative"]
            + coef["n_folds_zero"] == coef["n_folds"]).all()
    assert (coef["n_folds_nonzero_same_sign_as_mean"] <= coef["n_folds_nonzero"]).all()
    v = coef[fold_cols].to_numpy()
    assert not np.signbit(v[v == 0]).any(), "a negative zero reached the file"
    assert ((coef["rank_in_fold_best"] <= coef["rank_in_fold_worst"]).all())
    # a missingness indicator for a column with no missing value is a constant
    # zero in every fold, and the elastic net cannot give it weight
    ind = coef[(coef["term_kind"] == "missingness indicator")
               & (coef["n_missing_in_rows_used"] == 0)]
    assert (ind[fold_cols].to_numpy() == 0).all()


# ------------------------------------------------------------------ Table S9
def test_block_a_is_the_highest_tier_e3_reached(built, e3):
    ev, _ = e3
    a = _block(built[1], "a_insample_tier")

    def best(tool, stratum):
        r = ev[(ev.tool == tool) & (ev.stratum == stratum) & (ev.status == "ok")]
        got = [t for t in TIERS
               if (r.loc[r.tier == t, "pp3_threshold_reachable"].astype(str) == "True").any()]
        return got[-1] if got else M.NO_TIER

    singles = [t for t in ev.tool.unique() if t != K.FUSION]
    assert sorted(a["stratum"].unique()) == sorted(K.IN_SCOPE_STRATA)
    for r in a.itertuples():
        pool = [t for t in singles if t in K.PANEL] if r.comparison_set == "fusion inputs" \
            else singles
        assert r.n_single_columns_compared == len(pool)
        assert r.fusion_insample_tier == best(K.FUSION, r.stratum)
        tiers = {t: best(t, r.stratum) for t in pool}
        order = [M.NO_TIER] + TIERS
        top = max(tiers.values(), key=order.index)
        assert r.best_single_insample_tier == top
        assert set(r.single_columns_reaching_best_tier.split("; ")) == {
            t for t, v in tiers.items() if v == top}


def test_block_b_is_e3s_fold_table_for_the_fusion(built, e3, df):
    ev, folds = e3
    b = _block(built[1], "b_logo_fold")
    src = folds[(folds.tool == K.FUSION) & folds.stratum.isin(K.IN_SCOPE_STRATA)
                & folds.tier.isin(["supporting", "moderate", "strong"])]
    key = ["stratum", "tier", "heldout_gene"]
    m = b.merge(src, on=key, how="outer", indicator=True, suffixes=("", "_e3"))
    assert (m["_merge"] == "both").all() and len(b) == len(src)
    # the block's threshold column is E3's threshold_from_6_genes, renamed because
    # the distal band trains on five genes
    for mine, theirs in (("threshold_from_training_genes", "threshold_from_6_genes"),
                         ("heldout_lr", "heldout_lr_e3")):
        x, y = m[mine].to_numpy(dtype=float), m[theirs].to_numpy(dtype=float)
        assert np.array_equal(x, y, equal_nan=True), mine
    # every missing ratio is explained, and no present ratio is
    has_lr = b["heldout_lr"].notna()
    assert (b.loc[~has_lr, "why_no_heldout_lr"].fillna("") != "").all()
    assert (b.loc[has_lr, "why_no_heldout_lr"].fillna("") == "").all()
    assert (b.loc[has_lr, "heldout_n_in_band"] >= K.MIN_BAND).all()
    assert (b["heldout_n_damaging"] + b["heldout_n_normal"] == b["heldout_n"]).all()
    assert b["logo_threshold_leakage"].fillna("").str.len().gt(0).all()
    for s, g in b.groupby("stratum"):
        n_genes = K.stratum_frame(df, s).dropna(subset=["y_assay"])["gene"].nunique()
        assert (g["n_training_genes"] == n_genes - 1).all(), s


def test_block_c_is_the_coefficient_file_rounded_for_print(built):
    coef, s9 = built
    c = _block(s9, "c_coefficients").set_index("term")
    full = coef.set_index("term")
    dropped = set(full.index) - set(c.index)
    assert all(full.loc[t, "term_kind"] == "missingness indicator"
               and full.loc[t, "n_folds_nonzero"] == 0 for t in dropped)
    assert set(full.index[full["term_kind"] == "within-gene rank"]) <= set(c.index)
    for col in ("mean", "sd", "min", "max", "mean_abs"):
        want = np.round(full.loc[c.index, col].to_numpy(dtype=float), M.PRINT_DECIMALS)
        assert np.array_equal(c[col].to_numpy(dtype=float), want), col
    assert list(c["rank_mean_abs"]) == sorted(c["rank_mean_abs"])


# ------------------------------------------------------------- reproducibility
def test_a_rebuild_is_byte_identical_and_matches_the_files_on_disk(df, built):
    again = M.build(df, REPORTS)
    for first, second in zip(built, again):
        assert M.csv_text(first) == M.csv_text(second)
    for d, path in ((built[0], REPORTS / M.COEF_OUT.name),
                    (built[1], REPORTS / "supplement" / M.TABLE_OUT.name)):
        if path.exists():
            assert path.read_text() == M.csv_text(d), f"{path.name} is stale"


def test_the_fit_does_not_depend_on_the_blas_thread_count(refit):
    """Summation order in a multithreaded BLAS can move the last bit of a
    coefficient. The fit is run single-threaded in a fresh process and must hash
    the same as the default-threaded fit in this one."""
    feats, oof, coefs = refit
    here = hashlib.sha256(oof.tobytes() + coefs.to_numpy().tobytes()).hexdigest()
    code = (
        "import hashlib; from src import evid_common as K, evid_fusion_stability as M; "
        f"df = K.load_set(); o, c = M.refit_with_coefficients(df, {feats!r}); "
        "print(hashlib.sha256(o.tobytes() + c.to_numpy().tobytes()).hexdigest())")
    env = dict(os.environ, PYTHONPATH=str(PHASE1), OMP_NUM_THREADS="1",
               OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
               VECLIB_MAXIMUM_THREADS="1")
    out = subprocess.run([sys.executable, "-W", "ignore", "-c", code], cwd=PHASE1,
                         env=env, capture_output=True, text=True, check=True)
    assert out.stdout.strip() == here
