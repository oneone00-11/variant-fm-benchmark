"""
Guardrails for the IPW reweighting.

The degeneracy tests are the point of this file: every weighted estimator must
reduce EXACTLY to the unweighted one when all weights are 1. If it does not, a
"weighting changed nothing" result would be indistinguishable from a bug in the
weighting itself.

Run:  cd phase1 && python -m pytest tests/test_ipw_reweight.py -q
  or: cd phase1 && python -m tests.test_ipw_reweight
"""
from __future__ import annotations
import numpy as np
from scipy.stats import spearmanr

from src.ipw_reweight import (weighted_pearson, weighted_spearman, weighted_brier,
                              weighted_ece, weighted_yield, ess, make_weights,
                              guardrails)
from src.phase3_calibration import ece, brier, yield_metrics

RNG = np.random.default_rng(20260830)
N = 500


def _sample():
    x = RNG.normal(size=N)
    y = 0.6 * x + RNG.normal(size=N)
    p = RNG.uniform(size=N)
    lab = (RNG.uniform(size=N) < p).astype(float)
    return x, y, p, lab


# --- degeneracy: w == 1 must reproduce the unweighted estimator exactly -----
def test_weighted_spearman_degenerates_to_spearmanr():
    x, y, _, _ = _sample()
    assert abs(weighted_spearman(x, y, np.ones(N)) - spearmanr(x, y).statistic) < 1e-12


def test_weighted_spearman_degenerates_under_ties():
    """Average ranks, not ordinal ranks -- ties are the usual place this breaks."""
    x, y, _, _ = _sample()
    xt, yt = np.round(x, 1), np.round(y, 1)
    assert abs(weighted_spearman(xt, yt, np.ones(N)) - spearmanr(xt, yt).statistic) < 1e-12


def test_weighted_pearson_degenerates_to_corrcoef():
    x, y, _, _ = _sample()
    assert abs(weighted_pearson(x, y, np.ones(N)) - np.corrcoef(x, y)[0, 1]) < 1e-12


def test_weighted_brier_degenerates():
    _, _, p, lab = _sample()
    assert abs(weighted_brier(p, lab, np.ones(N)) - brier(p, lab)) < 1e-12


def test_weighted_ece_degenerates():
    _, _, p, lab = _sample()
    assert abs(weighted_ece(p, lab, np.ones(N)) - ece(p, lab)) < 1e-12


def test_weighted_yield_degenerates():
    _, _, p, lab = _sample()
    got = weighted_yield(p, lab, np.ones(N))
    ref = yield_metrics(p, lab)
    assert abs(got["actionable_frac"] - ref["actionable_frac"]) < 1e-4
    assert abs(got["actionable_acc"] - ref["actionable_acc"]) < 1e-4


# --- weights behave like replication ---------------------------------------
def test_integer_weight_equals_replication():
    _, _, p, lab = _sample()
    w = np.ones(N)
    w[:60] = 2.0
    dup_p, dup_lab = np.concatenate([p, p[:60]]), np.concatenate([lab, lab[:60]])
    assert abs(weighted_brier(p, lab, w)
               - weighted_brier(dup_p, dup_lab, np.ones(N + 60))) < 1e-12
    assert abs(weighted_ece(p, lab, w)
               - weighted_ece(dup_p, dup_lab, np.ones(N + 60))) < 1e-12


# --- normalisation and ESS --------------------------------------------------
def test_weights_normalise_to_mean_one():
    frame = RNG.normal(size=8000) ** 2
    obs = RNG.normal(size=2000) ** 2 + 0.4          # deliberately shifted
    for k in (5, 10, 20):
        w, per_bin, _ = make_weights(frame, obs, k)
        assert abs(np.mean(w) - 1.0) < 1e-9
        assert np.isfinite(w).all()
        assert len(per_bin) == k


def test_cutpoints_come_from_the_frame_not_the_sample():
    """p_full must be ~1/K by construction; if cutpoints were taken on the
    observed sample it would be p_int that came out flat instead."""
    frame = RNG.normal(size=8000) ** 2
    obs = RNG.normal(size=2000) ** 2 + 0.4
    _, per_bin, _ = make_weights(frame, obs, 10)
    assert np.allclose(per_bin["p_full"], 0.1, atol=0.002)
    assert not np.allclose(per_bin["p_int"], 0.1, atol=0.002)


def test_ess_bounds():
    assert abs(ess(np.ones(100)) - 100) < 1e-9         # equal weights -> n
    w = np.concatenate([np.full(99, 0.01), [99.01]])
    assert ess(w) < 5                                   # one dominating weight


def test_module_guardrails_pass():
    """The same checks the analysis script gates itself on."""
    assert guardrails() is True


if __name__ == "__main__":
    import sys
    fails = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"PASS {name}")
            except AssertionError as e:
                fails += 1
                print(f"FAIL {name}: {e}")
    sys.exit(1 if fails else 0)
