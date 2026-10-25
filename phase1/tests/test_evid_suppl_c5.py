"""E19 C5: held-out ratios with E7's within-gene intervals, Table 3's counts by the point
estimate and by the lower bound, and the Strong thresholds across the folds.

Every fold of E3's per-fold file (PP3, the three bands) has a row whose point estimate
is the file's; the point count is the printed Table 3 cell; the bound count is read off
the same rows; and the jackknife is the per-fold file's Strong thresholds.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

PHASE1 = Path(__file__).resolve().parents[1]
REPO = PHASE1.parent
OUT = PHASE1 / "reports/evidence/suppl_analyses"
EV = PHASE1 / "reports/evidence"

if str(PHASE1) not in sys.path:
    sys.path.insert(0, str(PHASE1))

from src import evid_common as K  # noqa: E402
from src import evid_suppl_c5 as C5  # noqa: E402

need = pytest.mark.skipif(not (OUT / "c5_heldout_lr_ci.csv").exists(),
                          reason="run the pipeline first")


@need
def test_every_fold_has_its_ratio_and_an_interval_around_it():
    ci = K.read_back(OUT / "c5_heldout_lr_ci.csv", keep_default_na=False, na_values=[""])
    folds = K.read_back(EV / "evidence_thresholds_logo_folds.csv")
    folds = folds[(folds.side == "pp3") & folds.stratum.isin(C5.BANDS)]
    key = ["tool", "stratum", "tier", "heldout_gene"]
    m = folds.merge(ci, on=key, how="outer", indicator=True, validate="one_to_one")
    assert (m._merge == "both").all() and len(m) == len(folds)
    same = (m.heldout_lr_x == m.heldout_lr_y) | (m.heldout_lr_x.isna() & m.heldout_lr_y.isna())
    assert same.all()
    thr = (m.threshold_from_6_genes == m.threshold_from_other_genes) | (
        m.threshold_from_6_genes.isna() & m.threshold_from_other_genes.isna())
    assert thr.all()
    fin = ci[ci.lr_lo.notna()]
    assert (fin.lr_lo <= fin.heldout_lr).all() and (fin.heldout_lr <= fin.lr_hi).all()
    assert ((ci.tool == K.FUSION) == ci.logo_threshold_leakage.notna()).all()
    assert (ci.clears_at_point == (ci.heldout_lr >= ci.cut_lr)).all()
    assert (ci.clears_at_bound == (ci.lr_lo >= ci.cut_lr)).all()


@need
def test_the_point_count_is_table_3s_and_the_bound_count_is_read_off_the_folds():
    cn = K.read_back(OUT / "c5_fold_counts_point_vs_bound.csv")
    cn["excluded_folds"] = cn.excluded_folds.fillna("")
    ci = K.read_back(OUT / "c5_heldout_lr_ci.csv", keep_default_na=False, na_values=[""])
    assert len(cn) == 13 * 3 * 2
    assert (cn.point_cell == cn.table3_cell).all()
    for r in cn[cn.in_sample_tier_reached].itertuples():
        excl = set(r.excluded_folds.split(";")) if r.excluded_folds else set()
        f = ci[(ci.tool == r.tool) & (ci.stratum == r.stratum) & (ci.tier == r.tier)
               & ~ci.heldout_gene.isin(excl)]
        assert r.n_folds == len(f)
        assert (r.bound_k, r.bound_evaluable) == (int(f.clears_at_bound.sum()),
                                                  int(f.lr_lo.notna().sum()))
        assert r.bound_k <= r.point_k
    assert set(cn[cn.tool == "avi"].excluded_folds) == {";".join(K.AVI_SEEN_IN_TRAINING)}
    assert (cn[~cn.in_sample_tier_reached].bound_cell == "–").all()


@need
def test_the_jackknife_reads_the_per_fold_strong_thresholds():
    jk = K.read_back(OUT / "c5_jackknife_strong.csv", keep_default_na=False, na_values=[""])
    folds = K.read_back(EV / "evidence_thresholds_logo_folds.csv")
    ev = K.read_back(EV / "evidence_thresholds.csv")
    for r in jk.itertuples():
        f = folds[(folds.tool == r.tool) & (folds.stratum == r.stratum) & (folds.tier == "strong")
                  & (folds.side == "pp3")]
        t = f.threshold_from_6_genes.dropna()
        assert (r.n_folds, r.n_folds_reaching) == (len(f), len(t))
        e = ev[(ev.tool == r.tool) & (ev.stratum == r.stratum) & (ev.tier == "strong")].iloc[0]
        ins = e.pp3_threshold_insample
        assert (np.isnan(ins) and np.isnan(r.in_sample_threshold)) or ins == r.in_sample_threshold
        if len(t):
            assert (r.threshold_min, r.threshold_max) == (t.min(), t.max())
            assert r.threshold_median == t.median()
        if len(t) and np.isfinite(ins):
            assert r.max_abs_shift == (t - ins).abs().max()


def test_the_stage_runs_after_c4():
    spec = importlib.util.spec_from_file_location(
        "reproduce_evidence", REPO / "scripts" / "reproduce_evidence.py")
    E = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(E)
    mods = [m for m, _, _ in E.STAGES]
    assert mods.index("src.evid_suppl_c5") == mods.index("src.evid_suppl_c4") + 1
