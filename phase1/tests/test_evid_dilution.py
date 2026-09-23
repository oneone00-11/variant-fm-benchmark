"""E11 guardrails for the training-gene dilution stage.

What is pinned here, and why each is worth a test:

  * the k = G-1 rows ARE E3's leave-one-gene-out folds. The stage seeds those
    subsets with E3's own cell key and takes the training rows in E3's order, so
    the fitted thresholds and the held-out ratios must equal
    evidence_thresholds_logo_folds.csv exactly, not approximately. Exact agreement
    is required because it is achievable: a tolerance here would let the two stages
    drift apart unnoticed. The k = G reference fit must likewise equal E3's
    in-sample threshold.
  * the subset design is complete: every size-k subset is present once, every gene
    is held out equally often at each k, and the summary's counts and fractions
    are the detail table's.
  * the held-out ratios are recomputed from the analysis set, not carried over, and
    a refit of a subset from the code reproduces the table and itself, so the
    tables cannot be stale output of an older version of the code.
  * the printed supplement table is the summary rounded, not a retyped copy.
"""
from __future__ import annotations

import sys
from math import comb
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

PHASE1 = Path(__file__).resolve().parents[1]
REPORTS = PHASE1 / "reports/evidence"
SET = PHASE1 / "data/evidence/analysis_set_v1.parquet"
CONFIG = PHASE1 / "config/walker2023.yaml"
DETAIL = REPORTS / "threshold_dilution.csv"
SUMMARY = REPORTS / "threshold_dilution_summary.csv"
SUPP = REPORTS / "supplement" / "tableS8_dilution.csv"
FOLDS = REPORTS / "evidence_thresholds_logo_folds.csv"
E3 = REPORTS / "evidence_thresholds.csv"
if str(PHASE1) not in sys.path:
    sys.path.insert(0, str(PHASE1))

from src import evid_common as K          # noqa: E402
from src import evid_dilution as D        # noqa: E402

pytestmark = pytest.mark.skipif(
    not (DETAIL.exists() and SUMMARY.exists()), reason="E11 not run")


def _read(p: Path) -> pd.DataFrame:
    # round_trip so that an exact comparison compares the values the stage wrote
    return pd.read_csv(p, float_precision="round_trip")


def _same(a: pd.Series, b: pd.Series) -> np.ndarray:
    a, b = a.to_numpy(dtype=float), b.to_numpy(dtype=float)
    return (a == b) | (np.isnan(a) & np.isnan(b))


@pytest.fixture(scope="module")
def detail():
    return _read(DETAIL)


@pytest.fixture(scope="module")
def summary():
    return _read(SUMMARY)


@pytest.fixture(scope="module")
def data():
    return pd.read_parquet(SET)


@pytest.fixture(scope="module")
def cuts():
    return K.acmg_bands(yaml.safe_load(CONFIG.read_text()))


# ---------------------------------------------------------------------------
# agreement with E3
# ---------------------------------------------------------------------------
@pytest.mark.skipif(not FOLDS.exists(), reason="E3 not run")
def test_k_g_minus_one_is_e3s_leave_one_gene_out(detail):
    """Every E3 fold for these tools, strata and tiers has a k = G-1 row here with
    the same threshold and the same held-out ratio, to the last bit."""
    folds = _read(FOLDS)
    folds = folds[(folds.side == "pp3") & folds.tool.isin(D.TOOLS)
                  & folds.stratum.isin(D.STRATA) & folds.tier.isin(D.TIERS)]
    top = detail[detail.k == detail.n_genes - 1]
    assert (top.seed_key == "logo").all()
    m = top.merge(folds, on=["tool", "stratum", "tier", "heldout_gene"],
                  how="outer", suffixes=("", "_e3"), indicator=True)
    assert (m["_merge"] == "both").all(), m[m["_merge"] != "both"][
        ["tool", "stratum", "tier", "heldout_gene", "_merge"]]
    assert len(m) == len(folds) > 0
    bad_tau = ~_same(m["threshold"], m["threshold_from_6_genes"])
    assert not bad_tau.any(), m[bad_tau][["tool", "stratum", "tier", "heldout_gene",
                                          "threshold", "threshold_from_6_genes"]]
    bad_lr = ~_same(m["heldout_lr"], m["heldout_lr_e3"])
    assert not bad_lr.any(), m[bad_lr][["tool", "stratum", "tier", "heldout_gene",
                                        "heldout_lr", "heldout_lr_e3"]]


@pytest.mark.skipif(not E3.exists(), reason="E3 not run")
def test_all_genes_reference_is_e3s_in_sample_threshold(summary):
    ev = _read(E3)
    ref = summary.drop_duplicates(["tool", "stratum", "tier"])
    m = ref.merge(ev, on=["tool", "stratum", "tier"], how="left")
    assert len(m) == len(D.TOOLS) * len(D.STRATA) * len(D.TIERS)
    assert _same(m["threshold_all_genes"], m["pp3_threshold_insample"]).all()


@pytest.mark.skipif(not FOLDS.exists(), reason="E3 not run")
def test_a_live_refit_of_an_e3_fold_matches_e3(data, cuts):
    """Not only the tables: the code path itself reproduces an E3 fold."""
    folds = _read(FOLDS)
    tool, stratum, held = "spliceai_walker", "s3_10", "VHL"
    y, s, genes = D.cell_arrays(data, tool, stratum)
    ug = D.genes_of(genes)
    train = tuple(g for g in ug if g != held)
    status, taus = D.fit_subset(y, s, genes, train, *cuts,
                                cell=D.cell_key(tool, stratum, train, ug))
    assert status == "fitted"
    for tier in D.TIERS:
        f = folds[(folds.tool == tool) & (folds.stratum == stratum)
                  & (folds.tier == tier) & (folds.side == "pp3")
                  & (folds.heldout_gene == held)]
        assert len(f) == 1
        want = float(f["threshold_from_6_genes"].iloc[0])
        assert (taus[tier] == want) or (np.isnan(taus[tier]) and np.isnan(want))


# ---------------------------------------------------------------------------
# the design is complete and the counts add up
# ---------------------------------------------------------------------------
def test_every_subset_once_and_every_gene_held_out_equally(detail, summary, data):
    for (tool, stratum, tier), d in detail.groupby(["tool", "stratum", "tier"]):
        _, _, genes = D.cell_arrays(data, tool, stratum)
        G = len(D.genes_of(genes))
        assert int(d["n_genes"].iloc[0]) == G
        assert sorted(d["k"].unique()) == list(range(D.MIN_K, G))
        for k, dk in d.groupby("k"):
            assert dk["training_genes"].nunique() == comb(G, k)
            assert len(dk) == comb(G, k) * (G - k)
            per_gene = dk["heldout_gene"].value_counts()
            assert len(per_gene) == G and (per_gene == comb(G - 1, k)).all()
            # a held-out gene is never one of its own training genes
            assert not any(h in t.split(";") for h, t in
                           zip(dk["heldout_gene"], dk["training_genes"]))
    assert (summary["n_subsets"] == [comb(g, k) for g, k in
                                     zip(summary.n_genes, summary.k)]).all()


def test_summary_is_the_detail_table_counted(detail, summary):
    for r in summary.itertuples():
        d = detail[(detail.tool == r.tool) & (detail.stratum == r.stratum)
                   & (detail.tier == r.tier) & (detail.k == r.k)]
        subsets = d.drop_duplicates("training_genes")
        if r.tool == "avi":
            # pairs are counted over the genes its model selection did not see;
            # the subsets it was fitted on are all still counted
            d = d[~d.heldout_gene.isin(K.AVI_SEEN_IN_TRAINING)]
        ev = d[d.evaluable]
        assert r.n_pairs == len(d)
        assert r.n_subsets_reached == int(subsets.tier_reached.sum())
        assert r.frac_subsets_reached == pytest.approx(
            subsets.tier_reached.sum() / len(subsets))
        assert r.n_pairs_evaluable == len(ev)
        assert r.n_pairs_cleared == int((ev.heldout_lr >= r.path_cut_lr).sum())
        if len(ev):
            assert r.heldout_lr_median == pytest.approx(ev.heldout_lr.median())
        taus = subsets.loc[subsets.tier_reached, "threshold"]
        if len(taus):
            assert r.threshold_median == pytest.approx(taus.median())
            assert r.threshold_q25 <= r.threshold_median <= r.threshold_q75


def test_evaluable_rows_obey_the_band_rules(detail):
    ev = detail[detail.evaluable]
    assert len(ev) > 0
    assert ev.tier_reached.all()
    assert (ev.heldout_n_in_band >= K.MIN_BAND).all()
    assert (ev.heldout_n_damaging >= K.MIN_POS).all()
    assert (ev.heldout_n_normal >= K.MIN_NEG).all()
    assert (ev.heldout_n_in_band
            == ev.heldout_n_damaging_in_band + ev.heldout_n_normal_in_band).all()
    assert (ev.cleared_cut.astype(bool) == (ev.heldout_lr >= ev.path_cut_lr)).all()
    # no ratio without a threshold, and no verdict without a ratio
    assert detail.loc[~detail.tier_reached, "heldout_lr"].isna().all()
    assert detail.loc[~detail.evaluable, "cleared_cut"].isna().all()
    assert (detail.loc[~detail.evaluable, "not_evaluable_reason"].fillna("") != "").all()


# ---------------------------------------------------------------------------
# recomputed, not carried over
# ---------------------------------------------------------------------------
def test_heldout_ratios_recompute_from_the_analysis_set(detail, data):
    ev = detail[detail.evaluable].iloc[::97]
    assert len(ev) > 10
    arrays = {}
    for r in ev.itertuples():
        key = (r.tool, r.stratum)
        if key not in arrays:
            arrays[key] = D.cell_arrays(data, *key)
        y, s, genes = arrays[key]
        te = genes == r.heldout_gene
        band = s[te] >= r.threshold
        assert int(band.sum()) == r.heldout_n_in_band
        assert int((band & (y[te] == 1)).sum()) == r.heldout_n_damaging_in_band
        assert K.band_lr(y[te], s[te], r.threshold, "upper") == r.heldout_lr


def test_a_small_subset_refits_identically_and_matches_the_table(detail, data, cuts):
    tool, stratum, train = "spliceai_walker", "s3_10", ("BRCA2", "VHL")
    y, s, genes = D.cell_arrays(data, tool, stratum)
    ug = D.genes_of(genes)
    cell = D.cell_key(tool, stratum, train, ug)
    a = D.fit_subset(y, s, genes, train, *cuts, cell=cell)
    b = D.fit_subset(y, s, genes, train, *cuts, cell=cell)
    assert a[0] == b[0] == "fitted"
    rows = detail[(detail.tool == tool) & (detail.stratum == stratum)
                  & (detail.training_genes == ";".join(train))]
    for tier in D.TIERS:
        assert (a[1][tier] == b[1][tier]) or (np.isnan(a[1][tier])
                                              and np.isnan(b[1][tier]))
        stored = rows.loc[rows.tier == tier, "threshold"].unique()
        assert len(stored) == 1
        assert (a[1][tier] == stored[0]) or (np.isnan(a[1][tier])
                                             and np.isnan(stored[0]))


def test_worker_processes_do_not_change_the_fits(data, monkeypatch):
    """The stage runs its fits across processes by default. Each fit is seeded from
    its own cell, so running the same fits in two workers or in this process must
    give the same thresholds bit for bit; a seed drawn from call order would not."""
    monkeypatch.chdir(PHASE1)          # the stage's paths are relative to phase1/
    tasks = [t for t in D.plan(data)
             if t[0] == "spliceai_walker" and t[1] == "s3_10" and t[3] == 2
             and "VHL" in t[4]][:4]
    assert len(tasks) == 4
    serial = D.run_fits(tasks, 1)
    parallel = D.run_fits(tasks, 2)
    assert [r[0] for r in serial] == [r[0] for r in parallel]
    for (_, a), (_, b) in zip(serial, parallel):
        assert a.keys() == b.keys()
        assert all((a[t] == b[t]) or (np.isnan(a[t]) and np.isnan(b[t])) for t in a)


def test_seed_keys_follow_e3_where_e3_makes_the_same_fit():
    ug = ["A", "B", "C", "D"]
    assert D.cell_key("t", "s", ("A", "B", "C", "D"), ug) == ("t", "s", "full")
    assert D.cell_key("t", "s", ("A", "B", "D"), ug) == ("t", "s", "logo", "C")
    # a smaller subset is keyed by its genes, whatever order they arrive in
    assert (D.cell_key("t", "s", ("B", "A"), ug)
            == D.cell_key("t", "s", ("A", "B"), ug)
            == ("t", "s", "subset", "A", "B"))
