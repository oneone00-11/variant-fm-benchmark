"""The phase-8 robustness numbers the manuscript reports are the numbers in the
pipeline output.

Same guardrail as test_reported_numbers.py: each pinned value is read back from
the committed `phase1/reports/phase1/phase8_*.csv`, so a change in the
robustness pipeline that moves a reported number fails here. The expected
values come from `tests/fixtures/expected_values.json`, generated from the
outputs by `scripts/build_expected_values.py` at printed precision; tolerances
are half the last printed digit (see test_reported_numbers.approx for why not
round()).
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
REPORTS = REPO / "phase1" / "reports" / "phase1"
FIXTURE = REPO / "tests" / "fixtures" / "expected_values.json"


def _csv(name):
    p = REPORTS / name
    if not p.exists():
        pytest.skip(f"{name} not built in this checkout")
    import pandas as pd
    return pd.read_csv(p)


def E(name):
    return json.loads(FIXTURE.read_text())["values"][name]["value"]


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
    assert approx(r.p_exact, E("sf.rho.incl.p"))
    num = E("sf.rho.incl.num")
    assert r.p_exact == pytest.approx(num["numerator"] / num["denominator"])
    r6 = d[(d.analysis == "delta_rho") & (d.condition == "BRCA1_excluded")].iloc[0]
    assert r6.k == 6 and r6.n_assignments == 64
    assert approx(r6.p_exact, E("sf.rho.excl.p"))


def test_sign_flip_delta_brier_primary_condition():
    """The calibration headline under the exact test: p as a multiple of 2/128."""
    d = _csv("phase8_sign_flip_exact.csv")
    r = d[(d.analysis == "delta_brier") &
          (d.condition == "isotonic/y_assay/BRCA1_included")].iloc[0]
    assert r.n_assignments == 128
    num = E("sf.brier.primary.num")
    assert r.p_exact == pytest.approx(num["numerator"] / num["denominator"])
    assert num["numerator"] % 2 == 0, "a two-sided exact p over 128 assignments is an even numerator"
    assert approx(r.p_exact, E("sf.brier.primary.p"))
    assert approx(r.T_obs, E("sf.brier.primary.T_obs"))
    assert r.p_exact < 0.05


def test_sign_flip_delta_brier_all_eight_conditions():
    """The full calibrator x label x +/-BRCA1 grid, pinned exactly."""
    d = _csv("phase8_sign_flip_exact.csv").set_index("condition")
    b = d[d.analysis == "delta_brier"]
    assert len(b) == 8
    for cond in ("isotonic/y_assay/BRCA1_included", "isotonic/y_assay/BRCA1_excluded",
                 "isotonic/y_clinvar/BRCA1_included", "isotonic/y_clinvar/BRCA1_excluded",
                 "platt/y_assay/BRCA1_included", "platt/y_assay/BRCA1_excluded",
                 "platt/y_clinvar/BRCA1_included", "platt/y_clinvar/BRCA1_excluded"):
        assert approx(b.loc[cond, "p_exact"], E(f"sf.brier.{cond}.p")), cond
        assert approx(b.loc[cond, "delta_brier_obs"], E(f"sf.brier.{cond}.obs")), cond
    # 7-cluster conditions cannot go below 2/128; 6-cluster not below 2/64
    assert (b[b.k == 7].p_exact >= 2 / 128).all()
    assert (b[b.k == 6].p_exact >= 2 / 64).all()


def test_wild_cluster_bootstrap_corroborates_exact_p():
    """Rademacher wild bootstrap (2000 draws) sits next to the exact p."""
    d = _csv("phase8_sign_flip_exact.csv")
    r = d[(d.analysis == "delta_brier") &
          (d.condition == "isotonic/y_assay/BRCA1_included")].iloc[0]
    assert approx(r.p_wild_rademacher, E("sf.brier.primary.wild"))
    r = d[(d.analysis == "delta_rho") & (d.condition == "BRCA1_included")].iloc[0]
    assert approx(r.p_wild_rademacher, E("sf.rho.incl.wild"))


# ---------------------------------------------------------------------------
# (b) leave-two-genes-out
# ---------------------------------------------------------------------------
def test_leave_two_genes_out_all_positive():
    d = _csv("phase8_leave_two_genes_out.csv")
    assert len(d) == 21
    assert (d.delta_rho > 0).all(), "the H1 gain must survive every gene pair"
    imin, imax = d.delta_rho.idxmin(), d.delta_rho.idxmax()
    assert d.loc[imin, "dropped"] == E("l2go.min.dropped")
    assert d.loc[imax, "dropped"] == E("l2go.max.dropped")
    assert approx(d.delta_rho.min(), E("l2go.min"))
    assert approx(d.delta_rho.median(), E("l2go.median"))
    assert approx(d.delta_rho.max(), E("l2go.max"))


# ---------------------------------------------------------------------------
# (c) Hartung-Knapp leaderboard
# ---------------------------------------------------------------------------
def test_hartung_knapp_widens_the_headline_rows():
    d = _csv("phase8_leaderboard_hk.csv").set_index("model")
    f = d.loc["M1_enet"]
    for c in ("pooled_rho", "dl_lo", "dl_hi", "hk_lo", "hk_hi"):
        assert approx(f[c], E(f"hk.M1_enet.{c}")), c
    nt = d.loc["single:nt"]
    for c in ("dl_lo", "dl_hi", "hk_lo", "hk_hi"):
        assert approx(nt[c], E(f"hk.nt.{c}")), c
    # HK widens the heterogeneous rows (I2 > 0)...
    het = d.dropna(subset=["hk_lo"])
    het = het[het.I2 > 0]
    assert ((het.hk_hi - het.hk_lo) >= (het.dl_hi - het.dl_lo) - 1e-12).all()
    # ...but on a tau^2 = 0 row the t-critical value does not compensate for a
    # small quadratic form, so HK can land NARROWER than DL; pinned so this is
    # seen rather than discovered (gpn_msa, I2 = 0)
    g = d.loc["single:gpn_msa"]
    assert g.I2 == 0.0
    assert approx(g.hk_lo, E("hk.gpn_msa.hk_lo")) and approx(g.hk_hi, E("hk.gpn_msa.hk_hi"))
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
    for c in ("brier_fusion", "brier_single", "rel_fusion", "rel_single", "res_fusion", "res_single", "unc"):
        assert approx(r[c], E(f"murphy.{c}")), c
    # the binned decomposition closes up to the within-bin variance residual
    # (bin-mean predictions are not the predictions); it must be small
    assert abs(r.brier_fusion - (r.rel_fusion - r.res_fusion + r.unc)) < 1e-3
    assert abs(r.brier_single - (r.rel_single - r.res_single + r.unc)) < 1e-3


def test_murphy_primary_deltas_and_intervals():
    """The Brier gain is a RESOLUTION gain: dRES excludes 0, dREL does not."""
    r = _murphy_primary()
    for c in ("d_brier", "d_brier_lo", "d_brier_hi", "d_rel", "d_rel_lo", "d_rel_hi", "d_res", "d_res_lo", "d_res_hi"):
        assert approx(r[c], E(f"murphy.{c}")), c
    assert r.d_res_lo > 0 and r.d_brier_lo > 0
    assert r.d_rel_lo < 0 < r.d_rel_hi


def test_murphy_grid_is_complete_and_brier_is_binning_invariant():
    d = _csv("phase8_murphy_decomposition.csv")
    assert len(d) == 12  # 4 conditions x 3 binnings
    assert d.condition.nunique() == 4 and d.binning.nunique() == 3
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
    for stratum in ("overall", "minus_BRCA1"):
        for obj in ("fusion_M1", "cadd"):
            for subset in ("ClinVar-recorded", "assay-only"):
                assert approx(_sel(stratum, obj, subset).lr_plus, E(f"sel8.{stratum}.{obj}.{subset}")), (stratum, obj, subset)
    assert _sel("overall", "fusion_M1", "ClinVar-recorded").lr_plus > _sel("overall", "fusion_M1", "assay-only").lr_plus


def test_selection_test_offset_strata():
    # core recorded arm is degenerate (too few negatives for a 95%-spec point);
    # the value is pinned so a change in that degeneracy is visible, and the
    # small negative count is asserted so the number is never read at face value
    core = _sel("offset_core_like", "fusion_M1", "ClinVar-recorded")
    assert approx(core.lr_plus, E("sel8.offset_core_like.fusion_M1.ClinVar-recorded")) and core.n_neg < 40
    for stratum in ("offset_core_like", "offset_region"):
        for obj in ("fusion_M1", "cadd"):
            for subset in ("ClinVar-recorded", "assay-only"):
                assert approx(_sel(stratum, obj, subset).lr_plus, E(f"sel8.{stratum}.{obj}.{subset}")), (stratum, obj, subset)


def test_selection_test_assay_only_record_type_split():
    for stratum in ("assay_only_VUS", "assay_only_Conflicting/Other"):
        for obj in ("fusion_M1", "cadd"):
            assert approx(_sel(stratum, obj, "assay-only").lr_plus, E(f"sel8.{stratum}.{obj}.assay-only")), (stratum, obj)


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
    assert (d.tier_basis == "interpolated to the nominal specificity").all()


def test_lr_plus_fusion_four_points():
    for brca in ("included", "excluded"):
        for spec in (0.9, 0.95, 0.975, 0.99):
            r = _lr("fusion_M1", spec, brca=brca)
            assert approx(r.lr_plus, E(f"op.fusion.{brca}.{spec}")), (brca, spec)
            assert approx(r.lr_plus_interp, E(f"op.fusion.{brca}.{spec}.interp")), (brca, spec)
            # the interpolated ratio is bounded by 1 / (1 - spec) by construction
            assert r.lr_plus_interp <= 1.0 / (1.0 - spec) + 1e-9


def test_lr_plus_key_singles_and_degeneracy():
    for obj, spec in (("pangolin", 0.95), ("spliceai", 0.95), ("cadd", 0.99), ("nt", 0.95), ("mean_M0b", 0.95)):
        r = _lr(obj, spec)
        assert approx(r.lr_plus, E(f"op.{obj}.included.{spec}")), (obj, spec)
        assert approx(r.lr_plus_interp, E(f"op.{obj}.included.{spec}.interp")), (obj, spec)
    # quantised scores cannot reach the strictest points at an observed value:
    # NaN, not a fake number -- while the interpolated value stays defined
    import math
    assert math.isnan(_lr("phastcons", 0.99).lr_plus)
    assert math.isfinite(_lr("phastcons", 0.99).lr_plus_interp)
    assert math.isnan(_lr("alphamissense", 0.95).lr_plus)


# ---------------------------------------------------------------------------
# (h) per-gene deltas behind the "six of seven genes" sentences (3.1, 3.3)
# ---------------------------------------------------------------------------
def test_per_gene_deltas_name_the_dissenting_genes():
    d = _csv("phase8_per_gene_deltas.csv")
    assert len(d) == 7 and set(d.best_single) == {E("pergene.best_single")}
    assert int(d.rho_favours_fusion.sum()) == E("pergene.rho_favours_n")
    assert int(d.brier_favours_fusion.sum()) == E("pergene.brier_favours_n")
    assert list(d.gene[~d.rho_favours_fusion]) == [E("pergene.rho_dissent")]
    assert list(d.gene[~d.brier_favours_fusion]) == [E("pergene.brier_dissent")]
    # the two axes do not fail on the same gene
    assert E("pergene.rho_dissent") != E("pergene.brier_dissent")
    # the columns are what their names say
    assert ((d.delta_rho > 0) == d.rho_favours_fusion).all()
    assert ((d.delta_brier > 0) == d.brier_favours_fusion).all()
