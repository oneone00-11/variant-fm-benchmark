"""The phase-8 robustness numbers the revision reports are the numbers in the
pipeline output.

Same guardrail as test_reported_numbers.py: each pinned value below is read
back from the committed `phase1/reports/phase1/phase8_*.csv`, so a change in
the robustness pipeline that moves a reported number fails here. Tolerances
are half the last printed digit (see test_reported_numbers.approx for why not
round()).
"""
from __future__ import annotations

from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
REPORTS = REPO / "phase1" / "reports" / "phase1"


def _csv(name):
    p = REPORTS / name
    if not p.exists():
        pytest.skip(f"{name} not built in this checkout")
    import pandas as pd
    return pd.read_csv(p)


def approx(value, printed):
    dp = len(str(printed).split(".")[1]) if "." in str(printed) else 0
    return abs(float(value) - float(printed)) <= 0.5 * 10 ** (-dp) + 1e-12


# ---------------------------------------------------------------------------
# (a) exact sign-flip tests
# ---------------------------------------------------------------------------
def test_sign_flip_delta_rho_exact_p():
    """2^7 = 128 assignments enumerated; the H1 gain's exact two-sided p."""
    d = _csv("phase8_sign_flip_exact.csv")
    r = d[(d.analysis == "delta_rho") & (d.condition == "BRCA1_included")].iloc[0]
    assert r.k == 7 and r.n_assignments == 128
    assert approx(r.p_exact, 0.1250)
    assert r.p_exact == pytest.approx(16 / 128)  # 16 assignments beat |T_obs|
    r6 = d[(d.analysis == "delta_rho") & (d.condition == "BRCA1_excluded")].iloc[0]
    assert r6.k == 6 and r6.n_assignments == 64
    assert approx(r6.p_exact, 0.2500)


def test_sign_flip_delta_brier_primary_condition():
    """The calibration headline under the exact test: p = 6/128 = 0.0469."""
    d = _csv("phase8_sign_flip_exact.csv")
    r = d[(d.analysis == "delta_brier") &
          (d.condition == "isotonic/y_assay/BRCA1_included")].iloc[0]
    assert r.n_assignments == 128
    assert r.p_exact == pytest.approx(6 / 128)
    assert approx(r.p_exact, 0.0469)
    assert approx(r.T_obs, 0.0085)


def test_sign_flip_delta_brier_all_eight_conditions():
    """The full calibrator x label x +/-BRCA1 grid, pinned exactly."""
    d = _csv("phase8_sign_flip_exact.csv").set_index("condition")
    expected = {
        "isotonic/y_assay/BRCA1_included": (0.0469, 0.0085),
        "isotonic/y_assay/BRCA1_excluded": (0.0625, 0.0116),
        "isotonic/y_clinvar/BRCA1_included": (0.1406, 0.0119),
        "isotonic/y_clinvar/BRCA1_excluded": (0.0938, 0.0192),
        "platt/y_assay/BRCA1_included": (0.0312, 0.0086),
        "platt/y_assay/BRCA1_excluded": (0.0625, 0.0108),
        "platt/y_clinvar/BRCA1_included": (0.1094, 0.0120),
        "platt/y_clinvar/BRCA1_excluded": (0.0938, 0.0169),
    }
    b = d[d.analysis == "delta_brier"]
    assert len(b) == 8
    for cond, (p, obs) in expected.items():
        assert approx(b.loc[cond, "p_exact"], p), cond
        assert approx(b.loc[cond, "delta_brier_obs"], obs), cond
    # 7-cluster conditions cannot go below 2/128; 6-cluster not below 2/64
    assert (b[b.k == 7].p_exact >= 2 / 128).all()
    assert (b[b.k == 6].p_exact >= 2 / 64).all()


def test_wild_cluster_bootstrap_corroborates_exact_p():
    """Rademacher wild bootstrap (2000 draws) sits next to the exact p."""
    d = _csv("phase8_sign_flip_exact.csv")
    r = d[(d.analysis == "delta_brier") &
          (d.condition == "isotonic/y_assay/BRCA1_included")].iloc[0]
    assert approx(r.p_wild_rademacher, 0.0560)
    r = d[(d.analysis == "delta_rho") & (d.condition == "BRCA1_included")].iloc[0]
    assert approx(r.p_wild_rademacher, 0.1310)


# ---------------------------------------------------------------------------
# (b) leave-two-genes-out
# ---------------------------------------------------------------------------
def test_leave_two_genes_out_all_positive():
    d = _csv("phase8_leave_two_genes_out.csv")
    assert len(d) == 21
    assert (d.delta_rho > 0).all(), "the H1 gain must survive every gene pair"
    imin, imax = d.delta_rho.idxmin(), d.delta_rho.idxmax()
    assert d.loc[imin, "dropped"] == "BARD1+BRCA1"
    assert d.loc[imax, "dropped"] == "BAP1+BRCA2"
    assert approx(d.delta_rho.min(), 0.0061)
    assert approx(d.delta_rho.median(), 0.0113)
    assert approx(d.delta_rho.max(), 0.0187)


# ---------------------------------------------------------------------------
# (c) Hartung-Knapp leaderboard
# ---------------------------------------------------------------------------
def test_hartung_knapp_widens_the_headline_rows():
    d = _csv("phase8_leaderboard_hk.csv").set_index("model")
    f = d.loc["M1_enet"]
    assert approx(f.pooled_rho, 0.773)
    assert approx(f.dl_lo, 0.740) and approx(f.dl_hi, 0.802)
    assert approx(f.hk_lo, 0.729) and approx(f.hk_hi, 0.810)
    nt = d.loc["single:nt"]
    assert approx(nt.dl_lo, 0.374) and approx(nt.dl_hi, 0.598)
    assert approx(nt.hk_lo, 0.367) and approx(nt.hk_hi, 0.603)
    # HK widens the heterogeneous rows (I2 > 0)...
    het = d.dropna(subset=["hk_lo"])
    het = het[het.I2 > 0]
    assert ((het.hk_hi - het.hk_lo) >= (het.dl_hi - het.dl_lo) - 1e-12).all()
    # ...but on a tau^2 = 0 row the t-critical value does not compensate for a
    # small quadratic form, so HK can land NARROWER than DL; pinned so this is
    # seen rather than discovered (gpn_msa, I2 = 0)
    g = d.loc["single:gpn_msa"]
    assert g.I2 == 0.0
    assert approx(g.hk_lo, 0.641) and approx(g.hk_hi, 0.689)
    assert (g.hk_hi - g.hk_lo) < (g.dl_hi - g.dl_lo)


# ---------------------------------------------------------------------------
# (d) Murphy decomposition
# ---------------------------------------------------------------------------
def _murphy_primary():
    d = _csv("phase8_murphy_decomposition.csv")
    return d[(d.condition == "isotonic/y_assay/BRCA1_included")
             & (d.binning == "width10")].iloc[0]


def test_murphy_primary_components():
    r = _murphy_primary()
    assert approx(r.brier_fusion, 0.0631) and approx(r.brier_single, 0.0717)
    assert approx(r.rel_fusion, 0.0047) and approx(r.rel_single, 0.0048)
    assert approx(r.res_fusion, 0.1911) and approx(r.res_single, 0.1822)
    assert approx(r.unc, 0.2499)
    # the binned decomposition closes up to the within-bin variance residual
    # (bin-mean predictions are not the predictions); it must be small
    assert abs(r.brier_fusion - (r.rel_fusion - r.res_fusion + r.unc)) < 1e-3
    assert abs(r.brier_single - (r.rel_single - r.res_single + r.unc)) < 1e-3


def test_murphy_primary_deltas_and_intervals():
    """The Brier gain is a RESOLUTION gain: dRES excludes 0, dREL does not."""
    r = _murphy_primary()
    assert approx(r.d_brier, 0.0085)
    assert approx(r.d_brier_lo, 0.0031) and approx(r.d_brier_hi, 0.0166)
    assert approx(r.d_rel, 0.0001)
    assert approx(r.d_rel_lo, -0.0008) and approx(r.d_rel_hi, 0.0067)
    assert approx(r.d_res, 0.0090)
    assert approx(r.d_res_lo, 0.0006) and approx(r.d_res_hi, 0.0128)
    assert r.d_res_lo > 0 and r.d_brier_lo > 0


def test_murphy_grid_is_complete_and_brier_is_binning_invariant():
    d = _csv("phase8_murphy_decomposition.csv")
    assert len(d) == 12  # 4 conditions x 3 binnings
    assert d.condition.nunique() == 4 and d.binning.nunique() == 3
    # Brier itself does not depend on the binning, only the REL/RES split does
    for cond, g in d.groupby("condition"):
        assert g.d_brier.nunique() == 1
        assert g.brier_fusion.nunique() == 1 and g.brier_single.nunique() == 1


# ---------------------------------------------------------------------------
# (e) stratified selection test
# ---------------------------------------------------------------------------
def _sel(stratum, obj, subset):
    d = _csv("phase8_selection_test_stratified.csv")
    r = d[(d.stratum == stratum) & (d.object == obj) & (d.subset == subset)]
    assert len(r) == 1
    return r.iloc[0]


def test_selection_test_baseline_and_minus_brca1():
    assert approx(_sel("overall", "fusion_M1", "ClinVar-recorded").lr_plus, 20.3)
    assert approx(_sel("overall", "fusion_M1", "assay-only").lr_plus, 14.1)
    assert approx(_sel("overall", "cadd", "ClinVar-recorded").lr_plus, 18.7)
    assert approx(_sel("overall", "cadd", "assay-only").lr_plus, 6.6)
    assert approx(_sel("minus_BRCA1", "fusion_M1", "ClinVar-recorded").lr_plus, 20.2)
    assert approx(_sel("minus_BRCA1", "fusion_M1", "assay-only").lr_plus, 13.7)
    assert approx(_sel("minus_BRCA1", "cadd", "ClinVar-recorded").lr_plus, 19.1)
    assert approx(_sel("minus_BRCA1", "cadd", "assay-only").lr_plus, 9.9)


def test_selection_test_offset_strata():
    # core recorded arm is degenerate (too few negatives for a 95%-spec point);
    # the value is pinned so a change in that degeneracy is visible, and the
    # small negative count is asserted so the number is never read at face value
    core = _sel("offset_core_like", "fusion_M1", "ClinVar-recorded")
    assert approx(core.lr_plus, 1.7) and core.n_neg < 40
    assert approx(_sel("offset_core_like", "fusion_M1", "assay-only").lr_plus, 4.2)
    assert approx(_sel("offset_core_like", "cadd", "ClinVar-recorded").lr_plus, 3.9)
    assert approx(_sel("offset_core_like", "cadd", "assay-only").lr_plus, 5.0)
    assert approx(_sel("offset_region", "fusion_M1", "ClinVar-recorded").lr_plus, 16.7)
    assert approx(_sel("offset_region", "fusion_M1", "assay-only").lr_plus, 13.3)
    assert approx(_sel("offset_region", "cadd", "ClinVar-recorded").lr_plus, 13.8)
    assert approx(_sel("offset_region", "cadd", "assay-only").lr_plus, 1.9)


def test_selection_test_assay_only_record_type_split():
    assert approx(_sel("assay_only_VUS", "fusion_M1", "assay-only").lr_plus, 12.1)
    assert approx(_sel("assay_only_VUS", "cadd", "assay-only").lr_plus, 4.4)
    assert approx(_sel("assay_only_Conflicting/Other", "fusion_M1", "assay-only").lr_plus, 17.6)
    assert approx(_sel("assay_only_Conflicting/Other", "cadd", "assay-only").lr_plus, 9.3)


# ---------------------------------------------------------------------------
# (f) LR+ at four operating points
# ---------------------------------------------------------------------------
def _lr(obj, spec, brca="included"):
    d = _csv("phase8_lr_plus_operating_points.csv")
    r = d[(d.condition == f"y_assay/BRCA1_{brca}") & (d.object == obj)
          & (d.target_spec == spec)]
    assert len(r) == 1
    return r.iloc[0]


def test_lr_plus_grid_covers_every_object_and_point():
    d = _csv("phase8_lr_plus_operating_points.csv")
    assert d.target_spec.nunique() == 4
    assert set(d.target_spec.unique()) == {0.9, 0.95, 0.975, 0.99}
    assert d.condition.nunique() == 2
    # 10 single tools + M1 + M0b; alphamissense is present but unevaluable
    # (it scores essentially no splice variant)
    assert d.object.nunique() == 12


def test_lr_plus_fusion_four_points():
    for spec, printed in ((0.90, 9.40), (0.95, 18.08), (0.975, 31.75), (0.99, 67.48)):
        assert approx(_lr("fusion_M1", spec).lr_plus, printed), spec
    for spec, printed in ((0.90, 9.31), (0.95, 17.95), (0.975, 32.44), (0.99, 74.21)):
        assert approx(_lr("fusion_M1", spec, brca="excluded").lr_plus, printed), spec


def test_lr_plus_key_singles_and_degeneracy():
    assert approx(_lr("pangolin", 0.95).lr_plus, 17.06)
    assert approx(_lr("spliceai", 0.95).lr_plus, 16.92)
    assert approx(_lr("cadd", 0.99).lr_plus, 57.87)
    assert approx(_lr("nt", 0.95).lr_plus, 6.46)
    assert approx(_lr("mean_M0b", 0.95).lr_plus, 16.57)
    # quantised scores cannot reach the strictest points: NaN, not a fake number
    import math
    assert math.isnan(_lr("phastcons", 0.99).lr_plus)
    assert math.isnan(_lr("alphamissense", 0.95).lr_plus)
