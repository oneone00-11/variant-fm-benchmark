"""The numbers in the manuscript are the numbers in the pipeline output.

The analysis is in this repository; what was missing was anything asserting
that the reported values and the pipeline agree. Each value below is quoted in
the manuscript and is read from `phase1/reports/phase1/`, so a change in the
pipeline that moves a headline number fails here rather than silently diverging
from the text.

The expected values are not typed into this file. They live in
`tests/fixtures/expected_values.json`, generated from the pipeline outputs by
`scripts/build_expected_values.py` at the manuscript's printed precision; each
entry records the output file, row and column it came from. Tolerances are half
the last printed digit of the value as the manuscript prints it, as before.
"""
from __future__ import annotations

import json
import re
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


def _ci(s):
    return [float(x) for x in re.findall(r"-?\d+\.?\d*", str(s))]


def E(name):
    """The expected value as the manuscript prints it (a string at printed precision,
    an int, a str, or a small structure), from the generated fixture."""
    return json.loads(FIXTURE.read_text())["values"][name]["value"]


def approx(value, printed):
    """Half the last printed digit, as the manuscript prints it. Python's round()
    is half-to-even and 0.0785 is stored just below the tie, so `round(x, 3)`
    would fail a value the manuscript correctly prints as 0.079."""
    dp = len(str(printed).split(".")[1]) if "." in str(printed) else 0
    return abs(float(value) - float(printed)) <= 0.5 * 10 ** (-dp) + 1e-12


def test_the_fixture_is_generated_from_the_current_outputs():
    """Regenerating the fixture must be a no-op: every expected value is the
    current pipeline output at printed precision. (A stale fixture is what a
    moved number looks like from the other side.)"""
    import subprocess, sys
    r = subprocess.run([sys.executable, str(REPO / "scripts" / "build_expected_values.py"), "--check"],
                       capture_output=True, text=True, cwd=REPO)
    assert r.returncode == 0, r.stdout + r.stderr


def test_h1_ranking_headline():
    """Results 3.1 and the Abstract: rho fusion vs best single, delta and its CI."""
    d = _csv("phase2_H1_stratified.csv")
    o = d[d.stratum == "overall"].iloc[0]
    assert approx(o.fusion_rho, E("h1.overall.fusion_rho"))
    assert approx(o.single_rho, E("h1.overall.single_rho"))
    assert approx(o.delta, E("h1.overall.delta"))
    lo, hi = _ci(o.ci95)
    assert approx(lo, E("h1.overall.ci_lo")) and approx(hi, E("h1.overall.ci_hi"))
    assert lo > 0, "the H1 interval no longer excludes zero"


def test_h1_stratified_gain_is_larger_in_the_core():
    """Results 3.1 states this explicitly; it was wrong in an earlier draft."""
    d = _csv("phase2_H1_stratified.csv").set_index("stratum")
    assert approx(d.loc["core_like", "delta"], E("h1.core_like.delta"))
    assert approx(d.loc["region", "delta"], E("h1.region.delta"))
    assert d.loc["core_like", "delta"] > d.loc["region", "delta"]
    for s in ("core_like", "region"):
        assert _ci(d.loc[s, "ci95"])[0] < 0, f"{s} interval no longer spans zero"


def test_h3_calibration_headline():
    """Results 3.3: delta-Brier by condition, isotonic."""
    d = _csv("phase3_H3_headline.csv")
    d = d[d.calib == "isotonic"].set_index("set")
    for cond in ("y_assay/BRCA1_included", "y_assay/BRCA1_excluded",
                 "y_clinvar/BRCA1_included", "y_clinvar/BRCA1_excluded"):
        assert approx(d.loc[cond, "dBrier"], E(f"h3.dBrier.{cond}")), cond
        assert _ci(d.loc[cond, "dBrier_ci_geneclust"])[0] > 0, f"{cond} no longer excludes zero"


def test_clinical_yield_numbers_and_direction():
    """Results 3.3: the fusion's and Pangolin's high-confidence fractions and
    accuracies in the primary condition, and which of the two is higher."""
    d = _csv("phase3_calibration_summary.csv")
    d = d[(d.set == "y_assay/BRCA1_included") & (d.calib == "isotonic")].set_index("model")
    fus, pan = d.loc["fusion_M1"], d.loc["single:pangolin"]
    assert approx(fus.actionable_frac * 100, E("h3.yield.fusion.frac_pct"))
    assert approx(fus.actionable_acc * 100, E("h3.yield.fusion.acc_pct"))
    assert approx(pan.actionable_frac * 100, E("h3.yield.pangolin.frac_pct"))
    assert approx(pan.actionable_acc * 100, E("h3.yield.pangolin.acc_pct"))
    direction = E("h3.yield.direction")
    got = "fusion_higher" if fus.actionable_frac > pan.actionable_frac else (
        "pangolin_higher" if pan.actionable_frac > fus.actionable_frac else "equal")
    assert got == direction, (
        f"the fusion/Pangolin high-confidence ordering is {got}, the manuscript describes "
        f"{direction}; the Results sentence needs rechecking")


def test_yield_is_significant_in_three_of_four_conditions():
    """The Abstract, Introduction, Results and Discussion all say three of four."""
    d = _csv("phase3_H3_headline.csv")
    d = d[d.calib == "isotonic"]
    sig = sum(1 for _, r in d.iterrows() if _ci(r.dYield_ci)[0] > 0)
    assert sig == E("h3.yield.n_significant_isotonic") == 3, f"yield is significant in {sig} conditions"


def test_tp53_external_validation():
    """Results 3.4: Brier of the frozen fusion against the best single tool on TP53,
    and the names of the two rows (the fusion carries all ten predictors)."""
    d = _csv("phase4_tp53_external.csv").set_index("model")
    fus_row = [i for i in d.index if i.startswith("fusion")][0]
    best = [i for i in d.index if i.startswith("best_single")][0]
    assert fus_row == E("tp53.fusion.row") and best == E("tp53.best.row")
    assert approx(d.loc[fus_row, "Brier"], E("tp53.fusion.brier"))
    assert approx(d.loc[best, "Brier"], E("tp53.best.brier"))
    assert d.loc[fus_row, "Brier"] < d.loc[best, "Brier"]


def test_sampling_frame_reweighting():
    """Methods 2.4 'Sampling frame': weight range and Kish ESS."""
    d = _csv("ipw_weight_diagnostics.csv").set_index("scheme")
    dec = d.loc["decile"]
    assert approx(dec.w_min, E("ipw.decile.w_min"))
    assert approx(dec.w_max, E("ipw.decile.w_max"))
    assert approx(dec.ESS_over_n, E("ipw.decile.ESS_over_n"))
    assert int(dec.n_w_gt_alert) == 0


def test_reweighting_does_not_reorder_the_predictors():
    """Methods 2.4: 'the predictor ordering is unchanged for all eleven models'."""
    d = _csv("ipw_predictor_table.csv")
    assert len(d) == 11
    assert (d.rank_unw == d.rank_wgt).all()


def test_sampling_frame_enrichment_is_larger_across_whole_genes():
    """The Methods sentence contrasts recruitment across whole genes with
    recruitment inside the splice window, and the direction it turns on."""
    whole = _csv("ipw_frame_balance_wholegene.csv").set_index("gene")
    window = _csv("ipw_frame_balance.csv").set_index("gene")
    w, n = whole.loc["POOLED"], window.loc["POOLED"]
    assert approx(w.mean_abs_z_observed, E("ipw.whole.observed")) and approx(w.mean_abs_z_remainder, E("ipw.whole.remainder"))
    assert approx(n.mean_abs_z_observed, E("ipw.window.observed")) and approx(n.mean_abs_z_remainder, E("ipw.window.remainder"))
    assert (w.mean_abs_z_observed - w.mean_abs_z_remainder) > \
           (n.mean_abs_z_observed - n.mean_abs_z_remainder)
    assert int(w.n_observed) == 1781


# --- evidence strength (H4) --------------------------------------------------

def _lr():
    return _csv("phase5_likelihood_ratios.csv")


def test_nothing_reaches_strong_on_the_functional_standard():
    """The central claim of 3.3 on the tier basis the manuscript uses -- the
    likelihood ratio interpolated to exactly 95% specificity: no object, raw or
    calibrated, reaches Strong against the assay labels. The observed-threshold
    crossings (a tie group landing past the 5% line) are counted beside it."""
    d = _lr()
    fun = d[(d.cal_method == "isotonic") & (d.condition.str.startswith("y_assay"))]
    assert len(fun) > 0
    assert (fun.tier_basis == "interpolated to the nominal specificity").all()
    n_interp = int(fun.acmg_tier.isin(["Strong", "Very strong"]).sum())
    assert n_interp == E("h4.assay_strong_interp_any_stage") == 0
    assert int(fun.acmg_tier_observed.isin(["Strong", "Very strong"]).sum()) == E("h4.assay_strong_observed_any_stage")


def test_clinvar_labels_reach_strong_where_the_assay_does_not():
    """The contrast the paragraph turns on: same predictors, same variants,
    same operating point, different labels."""
    d = _lr()
    raw = d[(d.cal_method == "isotonic") & (d.calibration == "raw")]
    cv = raw[raw.condition.str.startswith("y_clinvar")]
    for cond, g in cv.groupby("condition"):
        assert int((g.acmg_tier == "Strong").sum()) == E(f"h4.clinvar_strong_raw_interp.{cond}"), cond
        assert int((g.acmg_tier_observed == "Strong").sum()) == E(f"h4.clinvar_strong_raw_observed.{cond}"), cond
        assert int((g.acmg_tier == "Strong").sum()) >= 5, cond


def test_the_fusion_lr_advantage_and_the_tiers():
    """The paired difference in LR+ between the fusion and the best single tool,
    on both operating-point conventions, and the fact that carries the framing:
    on the functional standard no comparison moves a tier."""
    h = _csv("phase5_lr_headline.csv")
    iso = h[h.cal_method == "isotonic"]
    assert len(iso) == 8
    rows_obs = sorted(f"{r.condition}|{r.calibration}" for r in iso[iso.fusion_better].itertuples())
    rows_int = sorted(f"{r.condition}|{r.calibration}" for r in iso[iso.fusion_better_interp].itertuples())
    assert rows_obs == E("h4.observed_resolves")
    assert rows_int == E("h4.interp_resolves")
    assert int(iso.crosses_a_tier_boundary.sum()) == E("h4.interp_tier_crossings")
    assert int(iso.crosses_a_tier_boundary_observed.sum()) == E("h4.observed_tier_crossings")
    row = iso[(iso.condition == "y_assay/BRCA1_included") & (iso.calibration == "raw")].iloc[0]
    assert approx(row.lr_plus_fusion, E("h4.raw.fusion.lr_plus")) and approx(row.lr_plus_best_single, E("h4.raw.best.lr_plus"))
    assert approx(row.lr_plus_fusion_interp, E("h4.raw.fusion.lr_plus_interp"))
    assert approx(row.lr_plus_best_single_interp, E("h4.raw.best.lr_plus_interp"))
    assert row.tier_fusion == row.tier_best_single == "Moderate"


def test_calibration_effect_on_the_ratios():
    """Results 3.3: the median change in LR+ on calibrating and the band changes,
    on the interpolated basis."""
    c = _csv("phase5_calibration_effect.csv").set_index("cal_method").loc["isotonic"]
    assert approx(c.median_abs_change_lr_plus_interp, E("h4.cal_effect.median_abs_change_interp"))
    assert int(c.cells_changing_tier) == E("h4.cal_effect.cells_changing_tier")
    assert int(c.cells_with_a_tier) == E("h4.cal_effect.cells_with_a_tier")


def test_the_old_yield_definition_still_reproduces_the_published_number():
    """The high-confidence fraction the manuscript prints for the fusion, and the
    evidence-based definition reported beside it."""
    y = _csv("phase5b_evidence_yield.csv")
    r = y[(y.condition == "y_assay/BRCA1_included") & (y.object == "fusion_M1")].iloc[0]
    assert approx(r["yield_confidence_0.90_0.10"], E("h5b.fusion.yield_conf"))
    assert approx(r["yield_acmg_moderate_or_above"], E("h5b.fusion.yield_moderate"))


def test_elastic_net_weights_do_not_flip_sign_across_folds():
    """The stability claim. A sign flip in any feature would contradict it."""
    if not (REPORTS / "phase6_enet_weight_stability.csv").exists():
        pytest.skip("phase6 not built in this checkout")
    import pandas as pd
    w = pd.read_csv(REPORTS / "phase6_enet_weight_stability.csv", index_col=0)
    assert int(w.n_folds_sign_flip.sum()) == 0
    for f in ("alphagenome", "pangolin", "spliceai"):
        assert int(w.loc[f, "n_folds_nonzero"]) == 7
        assert approx(w.loc[f, "mean"], E(f"weights.{f}.mean")) and approx(w.loc[f, "sd"], E(f"weights.{f}.sd"))


def test_thinning_the_training_genes_degrades_brier_only_slightly():
    g = _csv("phase6_training_gene_summary.csv").set_index("n_training_genes")
    assert approx(g.loc[6, "mean"], E("thin.mean6")) and approx(g.loc[3, "mean"], E("thin.mean3"))
    assert g.loc[3, "mean"] > g.loc[6, "mean"]
    assert int(g.loc[3, "n_genes_worse_than_6"]) == E("thin.n_genes_worse")


# --- the label-standard contrast (3.3) ---------------------------------------

def test_the_two_label_sets_were_never_the_same_variants():
    b = _csv("phase7_base.csv").iloc[0]
    assert int(b.assay_labelled) == E("labels.assay_labelled") and int(b.clinvar_labelled) == E("labels.clinvar_labelled")
    assert int(b.both) == E("labels.both")
    assert int(b.assay_only) == E("labels.assay_only")


def test_the_label_sets_agree_where_they_overlap():
    """Rules out label divergence as the explanation."""
    a = _csv("phase7_label_agreement.csv").iloc[0]
    assert int(a.n) == E("labels.both")
    assert approx(a.agreement, E("labels.agreement")) and approx(a.kappa, E("labels.kappa"))
    assert approx(a.positive_rate_assay, E("labels.pos_assay")) and approx(a.positive_rate_clinvar, E("labels.pos_clinvar"))


def test_the_gap_almost_closes_on_the_shared_variants():
    c = _csv("phase7_label_contrast.csv")
    assert (c.n_fun.max() == E("labels.both"))
    strong_fun = int((c.tier_fun == "Strong").sum())
    strong_cv = int((c.tier_cv == "Strong").sum())
    assert (strong_fun, strong_cv) == (E("labels.strong_fun"), E("labels.strong_cv")), (strong_fun, strong_cv)


def test_recorded_variants_are_easier_under_the_same_functional_labels():
    """The decisive test for selection: no ClinVar label is involved on either
    side, only whether the variant carries a record."""
    s = _csv("phase7_selection_test.csv")
    piv = s.pivot_table(index="object", columns="subset", values="lr_plus")
    both = piv.dropna()
    assert len(both) >= 10
    assert (both["ClinVar-recorded"] > both["assay-only"]).all(), both.to_dict()
    sens = s.pivot_table(index="object", columns="subset", values="sensitivity")
    assert approx(sens.loc["fusion_M1", "ClinVar-recorded"], E("sel.fusion.recorded.sens"))
    assert approx(sens.loc["fusion_M1", "assay-only"], E("sel.fusion.assay_only.sens"))


def test_no_predictor_in_this_panel_carries_clinical_training_labels():
    """Why circularity cannot be the mechanism here: the group is empty."""
    g = _csv("phase7_delta_by_training_group.csv").set_index("group")
    assert int(g.loc["clinical labels", "n_objects"]) == 0


def test_the_disclosed_test_count_is_the_measured_one():
    """The disclosure states a number that lives in no reports/ file, so it is
    compared against the measurement recorded when it was taken. It went stale
    once already: the manuscript said thirty-three while a bare clone ran
    forty-seven, and nothing noticed because the checker had no source to
    compare against."""
    import re

    facts = REPO / "phase1" / "config" / "pipeline_facts.json"
    if not facts.exists():
        pytest.skip("pipeline facts not recorded in this checkout")
    f = json.loads(facts.read_text())
    want = int(f.get("guardrail_tests_passing_from_archive", f.get("guardrail_tests_passing_from_bare_clone")))

    manuscript = _manuscript_path()
    if manuscript is None:
        pytest.skip("manuscript not available in this checkout")
    from docx import Document

    units = ["", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine"]
    tens = {20: "twenty", 30: "thirty", 40: "forty", 50: "fifty", 60: "sixty",
            70: "seventy", 80: "eighty", 90: "ninety"}
    words = {(t if not u else f"{t}-{u}"): n + i
             for n, t in tens.items() for i, u in enumerate(units)}
    text = " ".join(p.text for p in Document(str(manuscript)).paragraphs)
    m = re.search(r"([A-Za-z-]+) automated tests", text)
    assert m, "the disclosure no longer states a test count"
    said = words.get(m.group(1).lower())
    assert said is not None, f"unrecognised count word {m.group(1)!r}; add it here"
    assert said == want, (f"the disclosure says {m.group(1)} ({said}) but the recorded "
                          f"measurement is {want}")


def _manuscript_path():
    import os

    env = os.environ.get("VARIANT_FM_MANUSCRIPT")
    if env:
        p = Path(env).expanduser()
        return p if p.exists() else None
    desk = Path.home() / "Desktop"
    hits = sorted([*desk.glob("calibration_draft*.docx"),
                   *(desk / "calibration").glob("calibration_draft*.docx"),
                   *desk.glob("draft_reframed_*.docx")],
                  key=lambda q: q.stat().st_mtime, reverse=True)
    return hits[0] if hits else None


def test_tp53_label_definitions_split_as_reported():
    """Table 4 / 3.4: the control-anchored and mid-band intervals exclude zero,
    the median-split interval spans it (a v2 fact; under frozen-matrix-v1 all
    three excluded zero)."""
    d = _csv("phase4_tp53_label_definitions.csv").set_index("label_rule")
    assert int(d.fusion_better.sum()) == E("tp53.rules.fusion_better_n")
    for rule, tag in (("control-anchored RFS>0", "control"), ("median split (top 50%)", "median"),
                      ("NA mid-band |RFS|>=0.5", "midband")):
        assert approx(d.loc[rule, "dBrier"], E(f"tp53.rules.{tag}.dBrier"))
        assert approx(d.loc[rule, "ci_lo"], E(f"tp53.rules.{tag}.ci_lo"))
        assert approx(d.loc[rule, "ci_hi"], E(f"tp53.rules.{tag}.ci_hi"))
        assert bool(d.loc[rule, "fusion_better"]) == (float(d.loc[rule, "ci_lo"]) > 0)


def test_leaderboard_is_the_pooling_of_its_per_gene_points():
    """Figure 1 draws phase2_per_gene_rho.csv against phase2_leaderboard.csv and Table 1
    prints the leaderboard: pooling the per-gene points must give the leaderboard's
    unrounded rho, interval and I^2 exactly, and the interval string must be those
    ends formatted once."""
    import sys
    import pandas as pd
    sys.path.insert(0, str(REPO / "phase1"))
    from src.phase2_model import dl_pool

    lb = _csv("phase2_leaderboard.csv").set_index("model")
    pg = _csv("phase2_per_gene_rho.csv")
    for m, r in lb.iterrows():
        g = pg[pg.object == m][["gene", "rho", "n"]]
        if pd.isna(r.pooled_rho):
            assert pd.isna(r.ci95) or r.ci95 == "n/a"   # pandas reads "n/a" as missing
            continue
        assert len(g) == int(r.k)
        d = dl_pool(g)
        for col, key in (("pooled_rho", "rho"), ("lo", "lo"), ("hi", "hi"), ("I2", "I2")):
            assert abs(d[key] - float(r[col])) < 1e-12, (m, col)
        assert f"[{r.lo:.3f},{r.hi:.3f}]" == r.ci95


def test_reliability_bins_cover_the_labelled_primary_set_for_both_curves():
    """Figure 2 plots both curves; each must bin every assay-labelled variant once."""
    cb = _csv("class_balance.csv")
    n = int(cb[(cb.label == "y_assay") & (cb.brca1 == "included")].n_labelled.iloc[0])
    for name in ("phase3_reliability_fusion.csv", "phase3_reliability_best_single.csv"):
        d = _csv(name)
        assert int(d.n.sum()) == n, name
        assert d.mean_pred.between(0, 1).all() and d.frac_pos.between(0, 1).all(), name


def test_table_1_rounds_the_unrounded_pooled_values_once():
    """Table 1 printed the fusion's I^2 as 60 for a true 59.49: the leaderboard then
    stored I^2 to one decimal (59.5) and the table rounded that again. The number check
    could not see it -- 60 lies within half a unit of 59.5. The leaderboard now stores
    its estimates unrounded and every Table 1 cell is checked against them, rounded
    half-up once."""
    from decimal import Decimal, ROUND_HALF_UP
    import pandas as pd

    manuscript = _manuscript_path()
    if manuscript is None:
        pytest.skip("manuscript not available in this checkout")
    from docx import Document

    def once(x, nd):
        return str(Decimal(repr(float(x))).quantize(Decimal(1).scaleb(-nd), rounding=ROUND_HALF_UP))

    names = {"Fusion — elastic net (M1)": "M1_enet", "Pangolin": "single:pangolin",
             "Fusion — gradient-boosted trees (M2)": "M2_gbt", "SpliceAI": "single:spliceai",
             "Fusion — equal-weight mean (M0b)": "M0b_mean", "AlphaGenome": "single:alphagenome",
             "CADD": "single:cadd", "GPN-MSA": "single:gpn_msa", "phyloP": "single:phylop",
             "phastCons": "single:phastcons", "Nucleotide Transformer": "single:nt",
             "gnomAD AF": "single:gnomad_af", "AlphaMissense": "single:alphamissense"}
    pool = _csv("phase2_leaderboard.csv").set_index("model").rename(columns={"pooled_rho": "rho"})
    tables = [t for t in Document(str(manuscript)).tables
              if t.rows[0].cells[0].text.strip().startswith("Predictor")
              and "Pooled" in t.rows[0].cells[2].text]
    assert len(tables) == 1, "Table 1 not found"
    seen, bad = set(), []
    for row in tables[0].rows[1:]:
        c = [x.text.strip() for x in row.cells]
        m = names[c[0]]
        seen.add(m)
        p = pool.loc[m]
        if pd.isna(p.rho):
            if c[2] != "—":
                bad.append(f"{c[0]}: {c[2]!r} for a non-evaluable row")
            continue
        want = (once(p.rho, 3), f"{once(p.lo, 3)} – {once(p.hi, 3)}", once(p.I2, 0))
        if tuple(c[2:5]) != want:
            bad.append(f"{c[0]}: printed {tuple(c[2:5])}, should be {want}")
    assert seen == set(names.values())
    assert not bad, "Table 1 cells that do not round the pooled estimates once:\n  " + "\n  ".join(bad)


def _round_once(x, nd):
    from decimal import Decimal, ROUND_HALF_UP
    return str(Decimal(repr(float(x))).quantize(Decimal(1).scaleb(-nd), rounding=ROUND_HALF_UP))


def test_table_3_rounds_the_calibration_summary_once():
    """Table 3 printed four cells from a four-decimal summary rounded again: Pangolin's
    ECE (0.032495 -> 0.0325 -> 0.033), SpliceAI's high-confidence share (1180/1675 =
    70.448% -> 0.7045 -> 70.5), and Nucleotide Transformer's ECE and share. The summary
    now stores its values unrounded; every cell is checked against it, rounded half-up
    once."""
    import pandas as pd

    manuscript = _manuscript_path()
    if manuscript is None:
        pytest.skip("manuscript not available in this checkout")
    from docx import Document

    names = {"Fusion — elastic net (M1)": "fusion_M1", "Pangolin": "single:pangolin",
             "SpliceAI": "single:spliceai", "Fusion — equal-weight mean (M0b)": "mean_M0b",
             "AlphaGenome": "single:alphagenome", "CADD": "single:cadd", "GPN-MSA": "single:gpn_msa",
             "phyloP": "single:phylop", "phastCons": "single:phastcons",
             "Nucleotide Transformer": "single:nt", "gnomAD AF": "single:gnomad_af",
             "AlphaMissense": "single:alphamissense"}
    pr = _csv("phase3_calibration_summary.csv")
    pr = pr[(pr.set == "y_assay/BRCA1_included") & (pr.calib == "isotonic")].set_index("model")
    tables = [t for t in Document(str(manuscript)).tables if t.rows[0].cells[2].text.strip() == "ECE"]
    assert len(tables) == 1, "Table 3 not found"
    seen, bad = set(), []
    for row in tables[0].rows[1:]:
        c = [x.text.strip() for x in row.cells]
        m = names[c[0]]
        seen.add(m)
        p = pr.loc[m]
        if pd.isna(p.ECE):
            if c[2] != "—":
                bad.append(f"{c[0]}: {c[2]!r} for a non-evaluable row")
            continue
        want = (_round_once(p.ECE, 3), _round_once(p.Brier, 3),
                _round_once(100 * p.actionable_frac, 1), _round_once(100 * p.actionable_acc, 1))
        if tuple(c[2:6]) != want:
            bad.append(f"{c[0]}: printed {tuple(c[2:6])}, should be {want}")
    assert seen == set(names.values())
    assert not bad, "Table 3 cells that do not round the summary once:\n  " + "\n  ".join(bad)


# Every column the manuscript prints at fewer decimals than the pipeline computes.
# Storing a rounded copy and rounding it again for the text misprinted six cells in
# the v2 revision (Table 1's I^2, three Table 3 cells, the Figure 2 caption, a 3.4
# share); these columns are now written unrounded, and this keeps them so. Interval
# strings and the bootstrap statistic behind dYield are formatted or rounded once by
# design and are not listed.
FULL_PRECISION = {
    "phase2_leaderboard.csv": ["pooled_rho", "I2", "lo", "hi"],
    "phase2_H1_stratified.csv": ["fusion_rho", "single_rho", "delta"],
    "phase2_H2_ablation.csv": ["full_rho", "ablated_rho", "delta_full_minus_ablated"],
    "phase3_calibration_summary.csv": ["ECE", "Brier", "actionable_frac", "actionable_acc"],
    "phase3_H3_headline.csv": ["dECE", "dBrier"],
    "phase3_pertool_brier_ci.csv": ["dBrier_vs_fusion"],
    "phase4_tp53_external.csv": ["ECE", "Brier", "actionable_frac", "actionable_acc"],
    "phase4_tp53_external_naband.csv": ["ECE", "Brier", "actionable_frac", "actionable_acc"],
    "class_balance.csv": ["pos_frac"],
    "coverage_by_region.csv": None,     # every region column
}


def test_reporting_tables_store_point_estimates_unrounded():
    import pandas as pd

    for name, cols in FULL_PRECISION.items():
        raw = pd.read_csv(REPORTS / name, dtype=str)
        for c in (cols or [x for x in raw.columns if x != "predictor"]):
            vals = [v for v in raw[c].dropna() if v not in ("", "n/a")]
            decimals = [len(v.split(".")[1]) if "." in v else 0 for v in vals]
            assert max(decimals) > 6, f"{name}:{c} is stored rounded ({max(decimals)} decimals)"


def test_tp53_orientation_gate_is_persisted_and_passes():
    """3.4 quotes the post-orientation control AUROC on the full TP53 SGE. It used to be
    printed once and copied; the manuscript and its supplement then disagreed on the
    third decimal (0.998 against 0.997). It is now written on every run."""
    g = _csv("phase4_tp53_orientation.csv").set_index("orientation")
    run = [i for i in g.index if i.startswith("as run")][0]
    dft = [i for i in g.index if i.startswith("default")][0]
    assert g.loc[run, "n_nonsense"] > 1000 and g.loc[run, "n_synonymous"] > 100
    assert abs(g.loc[run, "control_auroc"] + g.loc[dft, "control_auroc"] - 1) < 1e-12
    assert g.loc[run, "control_auroc"] > 0.99
    assert approx(g.loc[run, "control_auroc"], E("tp53.orientation.auroc"))


def test_every_measured_number_in_the_manuscript_is_bound():
    """Methods 2.6 says the tests fail if a reported result has no output behind it. This runs
    the number check on the manuscript and requires what the check requires: every measured
    number bound to a cell of an output its paragraph, table or caption declares, or to a
    derivation computed from one; no derivation taken on trust; no parameter rule whose code
    has moved; no pin left untaken; no cardinality mismatch. The coverage counts themselves
    are no longer printed in the manuscript."""
    import sys

    manuscript = _manuscript_path()
    if manuscript is None:
        pytest.skip("manuscript not available in this checkout")
    sys.path.insert(0, str(REPO / "phase1"))
    import src.check_manuscript_numbers as cmn

    r = cmn.classify(manuscript)
    unbound = [(t["uid"], t["token"]) for k in ("pool", "whitelisted", "no_source", "scoped_mismatch") for t in r[k]]
    assert not unbound, f"measured numbers with no bound output: {unbound}"
    assert r["coverage"]["bound"] == r["coverage"]["measured"]
    declared = [(t["uid"], t["token"]) for t in r["bound"] if "(declared, not computed)" in t["cells"][0]]
    assert not declared, f"derivations bound on trust rather than computed: {declared}"
    assert not (r["rule_problems"] or r["pin_problems"] or r["cardinality_mismatch"]), \
        (r["rule_problems"], r["pin_problems"], r["cardinality_mismatch"])
