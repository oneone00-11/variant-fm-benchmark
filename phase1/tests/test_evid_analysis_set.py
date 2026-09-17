"""E8 guardrails on the evidence-strength analysis set.

These check the things a reader would have to take on trust otherwise: that the
set on disk is the set the manifest describes, that every variant of the published
study is accounted for, and that no panel column is stored on the wrong sign.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest

PHASE1 = Path(__file__).resolve().parents[1]
SET = PHASE1 / "data/evidence/analysis_set_v1.parquet"
MANIFEST = PHASE1 / "data/evidence/analysis_set_v1.manifest.json"
REPORTS = PHASE1 / "reports/evidence"

pytestmark = pytest.mark.skipif(
    not SET.exists(), reason="evidence-strength analysis set not built (scripts/reproduce_evidence.py)")


@pytest.fixture(scope="module")
def manifest():
    return json.loads(MANIFEST.read_text())


@pytest.fixture(scope="module")
def df():
    return pd.read_parquet(SET)


def test_set_matches_its_manifest(df, manifest):
    assert len(df) == manifest["n_rows"]
    h = hashlib.sha256()
    with open(SET, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    assert h.hexdigest() == manifest["sha256"]


def test_selection_rule_holds(df):
    assert set(df["stratum"]) == {"pm12", "s3_10", "s11_50"}
    assert (df["ref"].str.len() == 1).all() and (df["alt"].str.len() == 1).all()
    off = df["intron_offset_abs"]
    assert off.between(1, 50).all(), "a variant outside 1 <= |offset| <= 50 got in"
    by_stratum = df.groupby("stratum")["intron_offset_abs"].agg(["min", "max"])
    assert by_stratum.loc["pm12", "max"] <= 2
    assert by_stratum.loc["s3_10", "min"] >= 3 and by_stratum.loc["s3_10", "max"] <= 10
    assert by_stratum.loc["s11_50", "min"] >= 11


def test_clinvar_arms_are_the_three_declared(df):
    assert set(df["clinvar_arm"]) == {"classified", "recorded_unclassified", "unrecorded"}
    # the arm must follow from the ClinVar class, not from set membership
    cls = df.loc[df["clinvar_arm"] == "classified", "clnsig_class"]
    assert set(cls) <= {"P/LP", "B/LB"}
    un = df.loc[df["clinvar_arm"] == "unrecorded", "clnsig_class"]
    assert un.isna().all()


def test_every_published_splice_variant_is_accounted_for(manifest):
    disp = pd.read_csv(REPORTS / "old_set_disposition.csv")
    assert len(disp) == 1781, "the published study's splice set is 1,781 variants"
    assert manifest["checks"]["old_splice_set_in_new"] \
        + manifest["checks"]["old_splice_set_explained"] == 1781
    # nothing may be unexplained
    assert not (disp["disposition"] == "other").any()
    assert disp.loc[disp["disposition"] != "in_new_set", "reason"].str.len().gt(0).all()


def test_hgvs_join_key_agrees_with_the_published_matrix(manifest):
    """The label join runs on the atlas's c. notation. It is only valid because the
    two products agree on that notation wherever they overlap."""
    assert manifest["checks"]["hgvs_agreement_with_frozen_v2_on_shared"] == 1.0
    assert manifest["checks"]["n_shared_with_frozen_v2"] > 3000


def test_no_panel_column_is_stored_on_the_wrong_sign():
    orient = pd.read_csv(REPORTS / "feature_orientation.csv")
    measured = orient.dropna(subset=["spearman_vs_pathogenicity"])
    wrong = measured[~measured["oriented_larger_is_damaging"].astype(bool)]
    assert wrong.empty, f"anti-correlated panel columns: {list(wrong['feature'])}"


def test_orientation_gate_passes_in_every_gene():
    gate = pd.read_csv(REPORTS / "directionality_check.csv")
    assert len(gate) == 7
    assert (gate["verdict"] == "OK").all(), gate[["gene", "verdict", "control_auroc"]]


def test_set_counts_agree_with_the_set(df):
    counts = pd.read_csv(REPORTS / "set_counts.csv")
    assert counts["n"].sum() == len(df)
    got = counts.groupby("stratum")["n"].sum().sort_index()
    want = df.groupby("stratum").size().sort_index()
    assert got.equals(want.astype(got.dtype))
