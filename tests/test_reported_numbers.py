"""The numbers in the manuscript are the numbers in the pipeline output.

Reviewer 2's complaint about this repository was that the analysis behind the
paper could not be found in it. The analysis is here; what was missing was
anything asserting that the reported values and the pipeline agree. Each value
below is quoted in the manuscript, and each is read from
`phase1/reports/phase1/`, so a change in the pipeline that moves a headline
number fails here rather than silently diverging from the text.

Tolerances are half the last printed digit of the value as it appears in the
manuscript.
"""
from __future__ import annotations

import re
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


def _ci(s):
    return [float(x) for x in re.findall(r"-?\d+\.?\d*", str(s))]


def approx(value, printed):
    """Half the last printed digit, as the manuscript prints it. Python's round()
    is half-to-even and 0.0785 is stored just below the tie, so `round(x, 3)`
    would fail a value the manuscript correctly prints as 0.079."""
    dp = len(str(printed).split(".")[1]) if "." in str(printed) else 0
    return abs(float(value) - float(printed)) <= 0.5 * 10 ** (-dp) + 1e-12


def test_h1_ranking_headline():
    """Results 3.1 and the Abstract: rho 0.773 vs 0.761, delta +0.012 [+0.002, +0.021]."""
    d = _csv("phase2_H1_stratified.csv")
    o = d[d.stratum == "overall"].iloc[0]
    assert approx(o.fusion_rho, 0.773)
    assert approx(o.single_rho, 0.761)
    assert approx(o.delta, 0.012)
    lo, hi = _ci(o.ci95)
    assert approx(lo, 0.002) and approx(hi, 0.021)


def test_h1_stratified_gain_is_larger_in_the_core():
    """Results 3.1 states this explicitly; it was wrong in an earlier draft."""
    d = _csv("phase2_H1_stratified.csv").set_index("stratum")
    assert approx(d.loc["core_like", "delta"], 0.0372)
    assert approx(d.loc["region", "delta"], 0.0217)
    assert d.loc["core_like", "delta"] > d.loc["region", "delta"]
    for s in ("core_like", "region"):
        assert _ci(d.loc[s, "ci95"])[0] < 0, f"{s} interval no longer spans zero"


def test_h3_calibration_headline():
    """Results 3.3: delta-Brier by condition, isotonic."""
    d = _csv("phase3_H3_headline.csv")
    d = d[d.calib == "isotonic"].set_index("set")
    for cond, val in (("y_assay/BRCA1_included", 0.0085), ("y_assay/BRCA1_excluded", 0.0116),
                      ("y_clinvar/BRCA1_included", 0.0119), ("y_clinvar/BRCA1_excluded", 0.0192)):
        assert approx(d.loc[cond, "dBrier"], val), cond
        assert _ci(d.loc[cond, "dBrier_ci_geneclust"])[0] > 0, f"{cond} no longer excludes zero"


def test_clinical_yield_numbers_and_direction():
    """Results 3.3, the sentence a reviewer flagged: 75.1% at 96.8% versus
    75.2% at 96.3%. The fusion's actionable fraction is LOWER here."""
    d = _csv("phase3_calibration_summary.csv")
    d = d[(d.set == "y_assay/BRCA1_included") & (d.calib == "isotonic")].set_index("model")
    fus, pan = d.loc["fusion_M1"], d.loc["single:pangolin"]
    assert approx(fus.actionable_frac * 100, 75.1)
    assert approx(fus.actionable_acc * 100, 96.8)
    assert approx(pan.actionable_frac * 100, 75.2)
    assert approx(pan.actionable_acc * 100, 96.3)
    assert fus.actionable_frac < pan.actionable_frac, (
        "the fusion's actionable fraction is no longer the lower of the two; the "
        "manuscript sentence describing them as indistinguishable needs rechecking")


def test_yield_is_significant_in_three_of_four_conditions():
    """The Abstract, Introduction, Results and Discussion all say three of four."""
    d = _csv("phase3_H3_headline.csv")
    d = d[d.calib == "isotonic"]
    sig = sum(1 for _, r in d.iterrows() if _ci(r.dYield_ci)[0] > 0)
    assert sig == 3, f"yield is significant in {sig} of four conditions, not three"


def test_tp53_external_validation():
    """Results 3.4: Brier 0.079 versus 0.122."""
    d = _csv("phase4_tp53_external.csv").set_index("model")
    fus = d.loc["fusion_8feat", "Brier"]
    best = [i for i in d.index if i.startswith("best_single")][0]
    assert approx(fus, 0.079)
    assert approx(d.loc[best, "Brier"], 0.122)
    assert fus < d.loc[best, "Brier"]


def test_sampling_frame_reweighting():
    """Methods 2.6 'Sampling frame': weights 0.72 to 1.24, Kish ESS 97% of nominal."""
    d = _csv("ipw_weight_diagnostics.csv").set_index("scheme")
    dec = d.loc["decile"]
    assert approx(dec.w_min, 0.72)
    assert approx(dec.w_max, 1.24)
    assert approx(dec.ESS_over_n, 0.97)
    assert int(dec.n_w_gt_alert) == 0


def test_reweighting_does_not_reorder_the_predictors():
    """Methods 2.6: 'the predictor ordering is unchanged for all eleven models'."""
    d = _csv("ipw_predictor_table.csv")
    assert len(d) == 11
    assert (d.rank_unw == d.rank_wgt).all()


def test_sampling_frame_enrichment_is_larger_across_whole_genes():
    """The Methods sentence contrasts recruitment across whole genes with
    recruitment inside the splice window. It reported 0.713 versus 0.578 for the
    whole-gene comparison, which no output held: ipw_frame_balance.csv is
    window-restricted, so that comparison had never been computed. Both figures
    are pinned here, and so is the direction the sentence turns on."""
    whole = _csv("ipw_frame_balance_wholegene.csv").set_index("gene")
    window = _csv("ipw_frame_balance.csv").set_index("gene")
    w, n = whole.loc["POOLED"], window.loc["POOLED"]
    assert approx(w.mean_abs_z_observed, 1.26) and approx(w.mean_abs_z_remainder, 0.62)
    assert approx(n.mean_abs_z_observed, 0.92) and approx(n.mean_abs_z_remainder, 0.80)
    # "Within the splice window ... the enrichment is much smaller."
    assert (w.mean_abs_z_observed - w.mean_abs_z_remainder) > \
           (n.mean_abs_z_observed - n.mean_abs_z_remainder)
    # The whole-gene frame is the full functional SNV set, not the analysis set.
    assert int(w.n_observed) == 1781


# --- evidence strength (H4) --------------------------------------------------

def _lr():
    return _csv("phase5_likelihood_ratios.csv")


def test_nothing_reaches_strong_on_the_functional_standard():
    """The central claim of 3.5. If a future change lets an object cross 18.7
    against the assay labels, the Results paragraph is wrong, not the test."""
    d = _lr()
    fun = d[(d.cal_method == "isotonic") & (d.condition.str.startswith("y_assay"))]
    assert len(fun) > 0
    assert not (fun.acmg_tier.isin(["Strong", "Very strong"]) &
                (fun.calibration == "raw")).any()


def test_clinvar_labels_reach_strong_where_the_assay_does_not():
    """The contrast the paragraph turns on: same predictors, same variants,
    same operating point, different labels."""
    d = _lr()
    raw = d[(d.cal_method == "isotonic") & (d.calibration == "raw")]
    cv = raw[raw.condition.str.startswith("y_clinvar")]
    n_strong = cv.groupby("condition").acmg_tier.apply(lambda s: (s == "Strong").sum())
    assert n_strong.min() >= 7 and n_strong.max() <= 8, n_strong.to_dict()


def test_the_fusion_lr_advantage_does_not_resolve():
    """Pinned because it is the result that kept the framing honest: the Brier
    advantage does not survive translation into evidence."""
    h = _csv("phase5_lr_headline.csv")
    iso = h[h.cal_method == "isotonic"]
    assert len(iso) == 8
    assert not iso.fusion_better.any(), iso[iso.fusion_better].to_dict("records")
    row = iso[(iso.condition == "y_assay/BRCA1_included") & (iso.calibration == "raw")].iloc[0]
    assert approx(row.lr_plus_fusion, 18.1) and approx(row.lr_plus_best_single, 17.1)
    assert row.lo < 0 < row.hi or row.lo < 0


def test_the_old_yield_definition_still_reproduces_the_published_number():
    """75.1% is what the manuscript has always printed for the fusion. The
    evidence-based definition is reported beside it, not instead of it, and this
    pins the bridge between the two."""
    y = _csv("phase5b_evidence_yield.csv")
    r = y[(y.condition == "y_assay/BRCA1_included") & (y.object == "fusion_M1")].iloc[0]
    assert approx(r["yield_confidence_0.90_0.10"], 0.751)
    assert approx(r["yield_acmg_moderate_or_above"], 0.868)


def test_elastic_net_weights_do_not_flip_sign_across_folds():
    """The stability claim. A sign flip in any feature would contradict it."""
    if not (REPORTS / "phase6_enet_weight_stability.csv").exists():
        pytest.skip("phase6 not built in this checkout")
    import pandas as pd
    w = pd.read_csv(REPORTS / "phase6_enet_weight_stability.csv", index_col=0)
    assert int(w.n_folds_sign_flip.sum()) == 0
    for f in ("alphagenome", "pangolin", "spliceai"):
        assert int(w.loc[f, "n_folds_nonzero"]) == 7


def test_thinning_the_training_genes_degrades_brier_only_slightly():
    g = _csv("phase6_training_gene_summary.csv").set_index("n_training_genes")
    assert approx(g.loc[6, "mean"], 0.0589) and approx(g.loc[3, "mean"], 0.0608)
    # monotone in the direction claimed, and worse in every held-out gene
    assert g.loc[3, "mean"] > g.loc[6, "mean"]
    assert int(g.loc[3, "n_genes_worse_than_6"]) == 7
