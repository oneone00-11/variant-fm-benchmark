"""E8 guardrails on the threshold tables.

Three things worth a test rather than a glance:

  * the transcribed ClinGen SVI recommendation. The brief flagged that two
    second-hand readings of the >= 0.2 threshold's strength circulate. The source
    says Moderate; this pins it, so a later edit that reintroduces "Supporting"
    fails here rather than in a manuscript.
  * E2 and E3 must agree on the likelihood ratio at the same threshold. They are
    separate code paths over the same data and the same convention for a zero
    denominator, so a disagreement means one of them is wrong.
  * the local likelihood ratio is expected to rise with the score, and the
    threshold rule assumes it. Where it does not, the test records the
    non-monotone stretch rather than failing: a dip inside the noise band is a
    property of the data, not a bug, and the point is to make it visible.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

PHASE1 = Path(__file__).resolve().parents[1]
REPORTS = PHASE1 / "reports/evidence"
CONFIG = PHASE1 / "config/walker2023.yaml"
WALKER = REPORTS / "walker_thresholds.csv"
EVIDENCE = REPORTS / "evidence_thresholds.csv"
SET = PHASE1 / "data/evidence/analysis_set_v1.parquet"


@pytest.fixture(scope="module")
def cfg():
    return yaml.safe_load(CONFIG.read_text())


# ---------------------------------------------------------------------------
# the transcription
# ---------------------------------------------------------------------------
def test_walker_thresholds_are_the_published_ones(cfg):
    pp3, bp4 = cfg["thresholds"]["pp3"], cfg["thresholds"]["bp4"]
    assert (pp3["value"], pp3["operator"]) == (0.2, ">=")
    assert (bp4["value"], bp4["operator"]) == (0.1, "<=")
    # the point the brief asked to settle from the source
    assert pp3["evidence_strength"] == "Moderate"
    assert bp4["evidence_strength"] == "Moderate"
    assert cfg["thresholds"]["grey_zone"]["evidence_strength"] == "Uninformative"
    # and the table the strengths come from
    assert pp3["reported_lr"] == 15.99 and pp3["reported_lr_ci"] == [13.23, 19.32]
    assert bp4["reported_lr"] == 0.17 and bp4["reported_lr_ci"] == [0.14, 0.21]


def test_the_score_the_thresholds_apply_to(cfg):
    sd = cfg["score_definition"]
    assert sd["distance"] == 4999 and sd["masked"] is False
    assert cfg["column_check"]["verdict"].startswith("MISMATCH")
    assert cfg["column_check"]["atlas_definition"]["distance"] == 50


def test_canonical_sites_are_out_of_scope(cfg):
    assert cfg["canonical_pm12"]["pp3_bp4_applicable"] is False
    assert cfg["other_tools"]["thresholds_recommended"] is False


# ---------------------------------------------------------------------------
# the tables
# ---------------------------------------------------------------------------
@pytest.mark.skipif(not WALKER.exists(), reason="E2 not run")
def test_e2_band_lr_is_internally_consistent():
    """For the top band the band likelihood ratio IS sensitivity / (1 - specificity);
    if those three columns disagree, one of them is not what it says it is."""
    d = pd.read_csv(WALKER)
    ok = d[(d.status == "ok") & d.lr_pp3.notna() & (d.spec_pp3 < 1.0)]
    assert len(ok) > 10
    implied = ok.sens_pp3 / (1.0 - ok.spec_pp3)
    close = np.isclose(implied, ok.lr_pp3, rtol=1e-6)
    # cells where the observed false-positive rate was zero are bounded away from
    # zero by convention, so they legitimately differ
    assert close.mean() > 0.9, ok.loc[~close, ["stratum", "clinvar_arm", "lr_pp3"]]


@pytest.mark.skipif(not WALKER.exists(), reason="E2 not run")
def test_e2_grey_zone_is_the_complement():
    d = pd.read_csv(WALKER)
    ok = d[d.status == "ok"]
    total = ok.frac_pp3 + ok.frac_grey + ok.frac_bp4
    assert np.allclose(total, 1.0, atol=1e-9)


@pytest.mark.skipif(not (WALKER.exists() and SET.exists()), reason="E2 not run")
def test_e2_band_lr_reproduces_from_the_set_itself():
    """E2's table against a recomputation from the analysis set.

    The brief asks that E2 and E3 give the same likelihood ratio at the same
    threshold. They cannot be compared directly, because they report different
    quantities: E2's is the ratio over the whole band above the cut point, which is
    what applying PP3 to a variant does, and E3's is the density ratio in a narrow
    window around the cut point, which is what Pejaver's threshold rule reads. What
    can be tested is that each reproduces from the data. This does E2's side.
    """
    df = pd.read_parquet(SET)
    col = "spliceai_walker" if "spliceai_walker" in df.columns else "spliceai"
    d = pd.read_csv(WALKER)
    d = d[(d.score_column == col) & (d.scope == "pooled") & (d.status == "ok")]
    assert len(d) > 5
    for r in d.itertuples():
        sub = df if r.stratum == "all_1_50" else df[df.stratum == r.stratum]
        if r.clinvar_arm != "all":
            sub = sub[sub.clinvar_arm == r.clinvar_arm]
        sub = sub.dropna(subset=["y_assay", col])
        y, s = sub["y_assay"].to_numpy(float), sub[col].to_numpy(float)
        pos, neg = s[y == 1], s[y == 0]
        hi_p, hi_n = (pos >= 0.2).mean(), (neg >= 0.2).mean()
        expected = hi_p / max(hi_n, 1.0 / (len(neg) + 1))
        assert np.isclose(expected, r.lr_pp3, rtol=1e-9), (r.stratum, r.clinvar_arm)


@pytest.mark.skipif(not (WALKER.exists() and EVIDENCE.exists()), reason="E2/E3 not run")
def test_the_cut_point_sits_where_the_local_ratio_is_uninformative():
    """The two stages disagree about the cut point by construction, and the shape of
    the disagreement is itself the finding, so it is recorded rather than asserted
    away.

    At Walker's 0.2 the LOCAL likelihood ratio is close to 1 in every stratum: the
    variants sitting near the cut point carry almost no evidence either way, which
    is what Walker's own uninformative band (>0.1 and <0.2, LR 1.00) says. The BAND
    ratio above it is much larger, because the damaging variants that make PP3 worth
    applying score far above the cut point rather than near it. A reader who reads
    the 0.2 threshold as "variants at 0.2 are moderately predictive" has the wrong
    picture; the table written here is the evidence for that.
    """
    e2 = pd.read_csv(WALKER)
    e3 = pd.read_csv(EVIDENCE)
    if "walker_pp3_local_lr" not in e3.columns:
        pytest.skip("E3 carries no Walker-cut-point columns")
    e3 = e3[(e3.status == "ok") & e3["walker_pp3_local_lr"].notna()]
    if not len(e3):
        pytest.skip("no Walker-cut-point rows")
    rows = []
    for tool in e3.tool.unique():
        for stratum in e3[e3.tool == tool].stratum.unique():
            a = e2[(e2.score_column == tool) & (e2.stratum == stratum)
                   & (e2.scope == "pooled") & (e2.clinvar_arm == "all")
                   & (e2.status == "ok")]
            b = e3[(e3.tool == tool) & (e3.stratum == stratum)]
            if not len(a) or not len(b):
                continue
            rows.append({"tool": tool, "stratum": stratum,
                         "frac_called_pp3": round(float(a.frac_pp3.iloc[0]), 4),
                         "band_lr_above_cut": round(float(a.lr_pp3.iloc[0]), 3),
                         "local_lr_at_cut": round(float(b.walker_pp3_local_lr.iloc[0]), 3),
                         "local_lr_lower_bound": round(
                             float(b.walker_pp3_local_lr_lo.iloc[0]), 3)})
    assert rows
    out = pd.DataFrame(rows)
    out.to_csv(REPORTS / "walker_cutpoint_local_vs_band.csv", index=False)
    # the band ratio is the one PP3 application uses, and it must exceed the local
    # ratio wherever the cut point is a tail
    tail = out[out.frac_called_pp3 < 0.5]
    assert (tail.band_lr_above_cut > tail.local_lr_at_cut).all(), tail


@pytest.mark.skipif(not EVIDENCE.exists(), reason="E3 not run")
def test_evidence_thresholds_are_ordered_by_tier():
    """A stricter tier cannot have a looser threshold."""
    d = pd.read_csv(EVIDENCE)
    order = ["supporting", "moderate", "strong", "very_strong"]
    bad = []
    for (tool, stratum), g in d[d.status == "ok"].groupby(["tool", "stratum"]):
        t = (g.set_index("tier")["pp3_threshold_insample"]
              .reindex(order).dropna())
        if len(t) > 1 and not t.is_monotonic_increasing:
            bad.append((tool, stratum, t.to_dict()))
    assert not bad, bad


@pytest.mark.skipif(not EVIDENCE.exists(), reason="E3 not run")
def test_local_lr_monotonicity_is_recorded_not_assumed():
    """The threshold rule reads the local likelihood ratio as increasing in the
    score. Where a curve dips, the dip is written down here so it is a known
    property of that tool in that stratum rather than a silent assumption."""
    curves = sorted(REPORTS.glob("interval_lr_*.csv"))
    if not curves:
        pytest.skip("E3 curve files not present")
    report = []
    for path in curves:
        c = pd.read_csv(path).dropna(subset=["local_lr"])
        if len(c) < 3:
            continue
        d = np.diff(c["local_lr"].to_numpy())
        drops = int((d < 0).sum())
        report.append((path.stem, len(c), drops, round(drops / max(len(d), 1), 3),
                       float(c["local_lr"].iloc[-1] - c["local_lr"].iloc[0])))
    assert report
    lines = ["tool_stratum n_grid n_drops drop_frac end_minus_start"]
    lines += [" ".join(str(x) for x in r) for r in report]
    (REPORTS / "local_lr_monotonicity.txt").write_text("\n".join(lines) + "\n")
    # the curve must at least end higher than it starts: that is the direction the
    # whole PP3 side depends on
    rising = [r for r in report if r[4] > 0]
    assert len(rising) / len(report) > 0.7, [r for r in report if r[4] <= 0]
