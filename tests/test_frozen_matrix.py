"""The frozen analysis set is what the manuscript says it is.

Every number in the paper is computed from `phase1/data/frozen/frozen_matrix_v1.parquet`.
`scripts/reproduce_calibration.py` verifies its content hash before running, but
that check only fires when someone runs the whole pipeline. These tests make the
same guarantees available in seconds, so a clone can confirm the analysis set
without a full rebuild.

The hash is over a canonical CSV serialisation with rows sorted by variant_id --
parquet itself is not byte-reproducible across writer versions.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
FROZEN = REPO / "phase1" / "data" / "frozen" / "frozen_matrix_v1.parquet"

# Published content hash of frozen-matrix-v1 (git tag frozen-matrix-v1,
# commit 2c78945a41; build commit 10d94bbe01), as cited in Supplementary Note S9.
EXPECTED_SHA256 = "2a0e249b44906f11ed91ce4746aba7389d8fff70d791a17710a3e26ec66d3199"
EXPECTED_ROWS, EXPECTED_GENES, EXPECTED_SPLICE = 21410, 7, 1781
GENES = {"BRCA1", "BRCA2", "BARD1", "PALB2", "RAD51C", "VHL", "BAP1"}


@pytest.fixture(scope="module")
def frozen():
    if not FROZEN.exists():
        pytest.skip("frozen matrix not built in this checkout")
    import pandas as pd
    return pd.read_parquet(FROZEN)


def test_content_hash_matches_the_published_value(frozen):
    canonical = frozen.sort_values("variant_id").reset_index(drop=True).to_csv(index=False)
    got = hashlib.sha256(canonical.encode()).hexdigest()
    assert got == EXPECTED_SHA256, (
        "the frozen matrix is not the one the manuscript used:\n"
        f"  expected {EXPECTED_SHA256}\n  got      {got}")


def test_shape_is_what_the_manuscript_reports(frozen):
    assert len(frozen) == EXPECTED_ROWS
    assert frozen["gene"].nunique() == EXPECTED_GENES
    assert set(frozen["gene"]) == GENES
    assert int(frozen["is_splice"].sum()) == EXPECTED_SPLICE


def test_splice_strata_sum_to_the_splice_set(frozen):
    """583 + 1,191 + 7 = 1,781, the decomposition Methods 2.1 states."""
    s = frozen[frozen["is_splice"]]["splice_subclass"].value_counts()
    assert s.get("splice_core", 0) == 583
    assert s.get("splice_region", 0) == 1191
    assert s.get("splice_exon_edge", 0) == 7
    assert s.sum() == EXPECTED_SPLICE


def test_class_balance_is_what_the_manuscript_reports(frozen):
    """Methods 2.1: 1,675 of 1,781 splice variants labelled, positive fraction 0.49.
    Reviewer 2 read the near-balance off this number, so it is pinned."""
    s = frozen[frozen["is_splice"]]
    lab = s["y_assay"].notna()
    assert int(lab.sum()) == 1675
    assert round(float(s.loc[lab, "y_assay"].mean()), 2) == 0.49


def test_the_deep_intronic_tier_is_empty_by_construction(frozen):
    """The manuscript's generalisation limit rests on this."""
    s = frozen[frozen["is_splice"]]
    assert (s["intron_offset"].abs() > 8).sum() == 0
