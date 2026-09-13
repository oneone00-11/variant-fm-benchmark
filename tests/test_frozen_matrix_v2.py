"""frozen-matrix-v2 is frozen-matrix-v1 with exactly three columns replaced.

The build (phase1/src/phase1_build_frozen_matrix_v2.py) swaps the SpliceAI, Pangolin
and Nucleotide Transformer columns for the pinned full-precision scorings and must
leave everything else untouched: same variants, same order, every other column
value-identical. These tests pin that contract, the manifest's content hash, and the
hash the reproduce scripts verify against, so that a rebuilt v2 that drifts from the
one the v1-vs-v2 comparison was made on cannot pass silently.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
FROZEN = REPO / "phase1" / "data" / "frozen"
V1, V2 = FROZEN / "frozen_matrix_v1.parquet", FROZEN / "frozen_matrix_v2.parquet"
MANIFEST = FROZEN / "manifest_v2.json"
REPLACED = {"spliceai": "spliceai_ds_fullprec", "pangolin": "pangolin_fullprec", "nt": "nucleotide_transformer"}


@pytest.fixture(scope="module")
def pair():
    if not (V1.exists() and V2.exists()):
        pytest.skip("frozen-matrix-v2 not built in this checkout")
    return pd.read_parquet(V1), pd.read_parquet(V2)


def _sha(df: pd.DataFrame) -> str:
    canon = df.sort_values("variant_id").reset_index(drop=True).to_csv(index=False)
    return hashlib.sha256(canon.encode()).hexdigest()


def test_same_variants_in_the_same_order(pair):
    v1, v2 = pair
    assert len(v1) == len(v2) == 21410
    assert (v1["variant_id"].to_numpy() == v2["variant_id"].to_numpy()).all()
    assert list(v1.columns) == list(v2.columns)


def test_every_untouched_column_is_value_identical(pair):
    v1, v2 = pair
    touched = set(REPLACED) | {f"{c}_isna" for c in REPLACED}
    for col in v1.columns:
        if col in touched:
            continue
        assert v1[col].equals(v2[col]), f"{col} differs between v1 and v2"


def test_replaced_columns_come_from_the_tracked_inputs(pair):
    _, v2 = pair
    atlas = pd.read_csv(REPO / "phase1" / "data" / "rescore" / "atlas_columns_v2.tsv", sep="\t").set_index("variant_id")
    for canon, col in REPLACED.items():
        only16 = pd.read_csv(REPO / "data" / "rescore" / f"{canon}_frozen_only16.tsv", sep="\t").set_index("variant_id")
        src = pd.concat([atlas[col], only16[col]])
        expected = src.loc[v2["variant_id"]].to_numpy(dtype=float)
        got = v2[canon].to_numpy(dtype=float)
        both = ~(np.isnan(expected) | np.isnan(got))
        assert (np.isnan(expected) == np.isnan(got)).all(), f"{canon}: missingness differs from the inputs"
        assert np.allclose(expected[both], got[both], rtol=0, atol=1e-12), f"{canon}: values differ from the inputs"
        assert (v2[f"{canon}_isna"].to_numpy() == np.isnan(got).astype(int)).all()


def test_full_precision_re_rounds_to_the_v1_cli_values(pair):
    """The v1 SpliceAI/Pangolin values are the CLI's 2-d.p. output (Pangolin's stored as float32);
    the v2 values are the same computation without the rounding, so no v1 value can sit more than
    half a unit in the last printed digit from its v2 value. Exact re-rounding is not the invariant:
    one SpliceAI variant (17-58696697-C-T) lies on a rounding boundary (v2 = 0.2149994, v1 printed
    0.22 by the WSL run that produced v1), a cross-platform float32 difference of ~1e-6."""
    v1, v2 = pair
    for canon in ("spliceai", "pangolin"):
        a, b = v1[canon].to_numpy(dtype=float), v2[canon].to_numpy(dtype=float)
        ok = ~(np.isnan(a) | np.isnan(b))
        assert np.nanmax(np.abs(b[ok] - a[ok])) <= 0.005 + 2e-6, canon


def test_manifest_hash_matches_the_matrix_and_the_reproduce_pin(pair):
    _, v2 = pair
    m = json.loads(MANIFEST.read_text())
    assert m["version"] == "v2" and m["derived_from"]["version"] == "v1"
    assert _sha(v2) == m["sha256"]
    sys.path.insert(0, str(REPO / "scripts"))
    from reproduce_calibration import PINS  # noqa: E402
    assert PINS["v2"][0] == m["sha256"], "scripts/reproduce_calibration.py pins a different v2 hash"
    assert (m["n_rows"], m["n_genes"], m["n_splice"]) == (21410, 7, 1781)


def test_v1_is_untouched():
    """v2 is an addition; the v1 pin must still hold."""
    if not V1.exists():
        pytest.skip("v1 not built")
    v1 = pd.read_parquet(V1)
    assert _sha(v1) == "2a0e249b44906f11ed91ce4746aba7389d8fff70d791a17710a3e26ec66d3199"
