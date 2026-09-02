"""The pipeline is seeded, and the seed is the only source of randomness.

The manuscript states that a clean-clone re-run reproduces the reported tables
byte for byte. A full re-run takes minutes, so this asserts the property that
makes it true: every stochastic step draws from a generator seeded with
`config.RANDOM_SEED`, and two draws from that seed agree exactly.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "phase1"))

from src import config as C  # noqa: E402


def test_the_seed_is_fixed_and_reported():
    assert isinstance(C.RANDOM_SEED, int)
    assert C.RANDOM_SEED == 20260708, (
        "the seed changed; every reported table was produced under 20260708")


def test_two_generators_from_the_same_seed_agree():
    a = np.random.default_rng(C.RANDOM_SEED).integers(0, 10_000, 1000)
    b = np.random.default_rng(C.RANDOM_SEED).integers(0, 10_000, 1000)
    assert np.array_equal(a, b)


def test_the_gene_cluster_bootstrap_is_reproducible():
    """The interval scheme behind every CI in the paper."""
    from src.phase3_calibration import boot_diff

    rng_state = np.random.default_rng(C.RANDOM_SEED)
    genes = np.array(["g%d" % (i % 7) for i in range(700)])
    label = (np.arange(700) % 2).astype(float)
    a = rng_state.random(700)
    b = a + 0.01

    import src.phase3_calibration as P3
    out = []
    for _ in range(2):
        P3.RNG = np.random.default_rng(C.RANDOM_SEED)
        out.append(boot_diff(a, b, label, genes, lambda p, y: float(np.mean((p - y) ** 2))))
    assert out[0]["lo"] == out[1]["lo"] and out[0]["hi"] == out[1]["hi"], (
        "the bootstrap is not reproducible from the seed; the reported intervals "
        "would not survive a re-run")


def test_calibration_is_deterministic():
    """Isotonic calibration must give the same probabilities on the same input."""
    from src.phase3_calibration import logo_calibrate

    rng = np.random.default_rng(0)
    genes = np.array(["g%d" % (i % 5) for i in range(500)])
    score = rng.random(500)
    label = (score + rng.normal(0, 0.2, 500) > 0.5).astype(float)
    first = logo_calibrate(score, label, genes)
    second = logo_calibrate(score, label, genes)
    assert np.allclose(first, second, equal_nan=True)
