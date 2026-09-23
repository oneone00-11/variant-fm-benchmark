"""E2.10 guardrails for the round that added AVI, the tier table and the arm controls.

Each test pins something a reader would otherwise have to take on trust: that the
AVI column is present and carries its provenance, that the tier table keeps the
in-sample claim and the held-out estimate apart, that the no-BRCA1 control really
is the same analysis minus one gene, and that the figure tables are derived rather
than retyped.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

PHASE1 = Path(__file__).resolve().parents[1]
DATA = PHASE1 / "data/evidence"
REPORTS = PHASE1 / "reports/evidence"
ATLAS_REPO = Path(os.environ.get(
    "EVID_ATLAS_REPO",
    str(Path(__file__).resolve().parents[2].parent / "functional-standard-atlas")))
if str(PHASE1) not in sys.path:
    sys.path.insert(0, str(PHASE1))
SET = DATA / "analysis_set_v1.parquet"

AVI_COLS = ["avi", "avi_splice_sites", "avi_splice_site_usage", "avi_splice_junctions"]


# --------------------------------------------------------------------- AVI
@pytest.mark.skipif(not SET.exists(), reason="analysis set not built")
def test_the_avi_columns_are_present_and_complete():
    df = pd.read_parquet(SET)
    for c in AVI_COLS:
        assert c in df.columns, f"{c} missing from the analysis set"
        assert df[c].notna().all(), f"{c} has {int(df[c].isna().sum())} gaps"


@pytest.mark.skipif(not (DATA / "avi.parquet.provenance.json").exists(),
                    reason="AVI not scored")
def test_the_avi_provenance_records_what_a_reader_needs():
    p = json.loads((DATA / "avi.parquet.provenance.json").read_text())
    for field in ("client_version", "scorers_requested", "accessed_utc",
                  "n_variants", "n_scored", "licence", "licence_source",
                  "output_sha256", "aggregation", "failed_batches"):
        assert field in p, f"provenance has no {field}"
    assert p["client_version"].startswith("0.9"), p["client_version"]
    # the licence is more restrictive than the repo's default and must say so
    assert "NON-COMMERCIAL" in p["licence"].upper()
    assert "train" in p["licence"].lower()
    # every requested scorer produced a complete column
    assert all(v == p["n_variants"] for v in p["n_scored"].values()), p["n_scored"]
    assert p["failed_batches"] == [], p["failed_batches"]


@pytest.mark.skipif(not list(DATA.glob("avi_scorers_*.json")), reason="no scorer list")
def test_the_scorer_list_was_recorded_verbatim():
    """The AVI scorer's key is served at runtime, so the list is the record of
    which scorer was chosen and what else was on offer."""
    path = sorted(DATA.glob("avi_scorers_*.json"))[-1]
    rec = json.loads(path.read_text())
    assert rec["n_scorers"] == len(rec["scorers"])
    assert rec["n_scorers"] > 10, "suspiciously short scorer list"
    assert "AVI_SCORE" in rec["scorers"]


@pytest.mark.skipif(not (PHASE1.parent / "LICENSE-DATA").exists(), reason="no LICENSE-DATA")
def test_the_avi_licence_is_declared_per_column():
    text = (PHASE1.parent / "LICENSE-DATA").read_text()
    for c in AVI_COLS:
        assert c in text, f"{c} is not named in LICENSE-DATA"
    assert "NON-COMMERCIAL USE ONLY" in text.upper()


# ------------------------------------------------------------- tier table
@pytest.mark.skipif(not (REPORTS / "tier_logo_table.csv").exists(), reason="E2.2 not run")
def test_the_tier_table_keeps_the_two_quantities_apart():
    d = pd.read_csv(REPORTS / "tier_logo_table.csv")
    ok = d[d.status == "ok"]
    assert {"in_sample_tier_reached", "heldout_lr_median", "folds_reached",
            "folds_heldout_lr_above_cut"} <= set(d.columns)
    # the fusion is excluded by design -- its threshold does not carry across genes
    assert "fusion_enet" not in set(d.tool)
    # an unreached tier must not carry an in-sample threshold
    unreached = ok[~ok.in_sample_tier_reached.fillna(False)]
    assert unreached["in_sample_threshold"].isna().all()
    # the two quantities must be able to disagree, and on this data they do:
    # a tier reached on the stratum whose held-out folds never clear its own cut
    disagree = ok[(ok.in_sample_tier_reached.fillna(False))
                  & (ok.folds_heldout_lr_above_cut == 0)]
    assert len(disagree), ("no cell where in-sample and held-out disagree; either "
                           "the data changed or the two columns are the same number")


@pytest.mark.skipif(not (REPORTS / "tier_logo_table.csv").exists(), reason="E2.2 not run")
def test_the_tier_table_is_derived_from_the_threshold_tables():
    """It is a view, not a retyping: every row must match evidence_thresholds.csv."""
    d = pd.read_csv(REPORTS / "tier_logo_table.csv")
    ev = pd.read_csv(REPORTS / "evidence_thresholds.csv")
    ok = d[d.status == "ok"]
    m = ok.merge(ev, on=["tool", "stratum", "tier"], suffixes=("", "_ev"))
    assert len(m) == len(ok)
    assert (m["in_sample_tier_reached"].astype(bool)
            == m["pp3_threshold_reachable"].astype(bool)).all()
    got = m["in_sample_threshold"].to_numpy(dtype=float)
    want = m["pp3_threshold_insample"].to_numpy(dtype=float)
    assert np.allclose(got, want, equal_nan=True)


# --------------------------------------------------------------- the arms
@pytest.mark.skipif(not (REPORTS / "territory_metrics_arms_noBRCA1.csv").exists(),
                    reason="E2.3 not run")
def test_the_noBRCA1_control_is_the_same_analysis_minus_one_gene():
    d = pd.read_csv(REPORTS / "territory_metrics_arms_noBRCA1.csv")
    assert set(d.gene_set) == {"all_genes", "no_BRCA1"}
    a = d[d.gene_set == "all_genes"]
    b = d[d.gene_set == "no_BRCA1"]
    # the same grid of cells is attempted in both
    key = ["stratum", "clinvar_arm", "tool"]
    assert set(map(tuple, a[key].values)) == set(map(tuple, b[key].values))
    # dropping a gene cannot add variants to a cell
    m = a.merge(b, on=key, suffixes=("_all", "_no"))
    assert (m["n_no"] <= m["n_all"]).all()
    # BRCA1 has no unrecorded variants, so that arm must be unchanged by the drop
    un = m[(m.clinvar_arm == "unrecorded")]
    assert (un["n_no"] == un["n_all"]).all(), \
        "dropping BRCA1 changed the unrecorded arm, which should hold no BRCA1"


@pytest.mark.skipif(not (REPORTS / "arms_within_gene.csv").exists(), reason="E2.3 not run")
def test_the_within_gene_table_has_genes_carrying_all_three_arms():
    w = pd.read_csv(REPORTS / "arms_within_gene.csv")
    ok = w[w.status == "ok"]
    per_gene = ok[ok.stratum == "s3_50"].groupby("gene")["clinvar_arm"].nunique()
    assert (per_gene == 3).any(), \
        "no gene contributes all three arms, so the confound cannot be removed within gene"


@pytest.mark.skipif(not (REPORTS / "clinvar_arm_no_assertion.csv").exists(),
                    reason="E2.4 not run")
def test_the_no_assertion_records_are_visible_and_inside_the_middle_arm():
    d = pd.read_csv(REPORTS / "clinvar_arm_no_assertion.csv")
    strict = d[d.clinvar_arm_strict == "recorded_no_assertion"]
    assert len(strict), "no_assertion records are not broken out"
    # they belong to the middle arm, not to classified or unrecorded
    assert set(strict.clinvar_arm) == {"recorded_unclassified"}
    assert strict["n"].sum() > 100


# ------------------------------------------------------------ figure data
@pytest.mark.skipif(not (REPORTS / "fig_data").exists(), reason="E2.9 not run")
def test_every_figure_table_exists_and_is_non_empty():
    fd = REPORTS / "fig_data"
    for name in ("fig1_walker_cutpoints.csv", "fig2_interval_lr.csv",
                 "fig2_tier_lines.csv", "fig3_evidence_by_territory.csv",
                 "fig4_external.csv"):
        p = fd / name
        assert p.exists(), f"{name} missing"
        assert len(pd.read_csv(p)), f"{name} is empty"


@pytest.mark.skipif(not (REPORTS / "fig_data/fig1_walker_cutpoints.csv").exists(),
                    reason="E2.9 not run")
def test_the_figure_tables_only_carry_in_scope_strata():
    """A figure must not invite a threshold reading off a pool that contains the
    canonical sites the recommendation excludes."""
    fd = REPORTS / "fig_data"
    for name in ("fig1_walker_cutpoints.csv", "fig2_interval_lr.csv",
                 "fig3_evidence_by_territory.csv", "fig4_external.csv"):
        d = pd.read_csv(fd / name)
        if "stratum" in d.columns:
            assert set(d.stratum) <= {"s3_10", "s11_50", "s3_50"}, \
                f"{name} carries an out-of-scope stratum: {set(d.stratum)}"


# ---------------------------------------------------------------------------
# E9 -- the predictor training-signal table
# ---------------------------------------------------------------------------
def test_every_evaluated_column_has_a_training_provenance_row():
    """The table answers the question of hidden dependencies, so a
    column that is evaluated anywhere and missing here is the failure mode."""
    t = pd.read_csv(REPORTS / "predictor_training_provenance.csv")
    have = set(t.score_column)
    evaluated = set(pd.read_csv(REPORTS / "evidence_thresholds.csv").tool)
    missing = evaluated - have
    assert not missing, f"evaluated but no training provenance row: {sorted(missing)}"


def test_the_training_rows_that_come_from_the_atlas_are_not_restated():
    """Nine columns are carried from the companion atlas's curated table. If this
    module ever restates one instead of importing it, the two papers can disagree
    about what a model was trained on."""
    import importlib.util
    import sys as _sys

    t = pd.read_csv(REPORTS / "predictor_training_provenance.csv")
    spec = importlib.util.spec_from_file_location(
        "_atlas_pr", ATLAS_REPO / "src/atlas/predictor_resources.py")
    mod = importlib.util.module_from_spec(spec)
    _sys.modules.setdefault("_atlas_pr", mod)
    spec.loader.exec_module(mod)

    from src.evid_training_provenance import FROM_ATLAS
    checked = 0
    for col, name in FROM_ATLAS.items():
        row = t[t.score_column == col]
        assert len(row) == 1, col
        assert row.training_data.iloc[0] == mod.TRAINING[name][0], col
        assert row.row_origin.iloc[0].startswith("companion atlas"), col
        checked += 1
    assert checked == len(FROM_ATLAS)


def test_no_panel_column_is_trained_on_clinical_classifications():
    """The claim the manuscript makes from this table. It is a statement about
    the panel as assembled, so it is asserted against the file, not in prose."""
    t = pd.read_csv(REPORTS / "predictor_training_provenance.csv")
    offenders = t.loc[t.contains_clinical_labels != "No", "score_column"].tolist()
    assert not offenders, f"trained on clinical classifications: {offenders}"
