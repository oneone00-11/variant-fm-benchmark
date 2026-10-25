"""E19 C3: the in-sample thresholds refitted under other labels, and read for the
subsets without BARD1 and PALB2.

The label variants change exactly the variants they name and no other; the join of a
label file is attach_labels' own (the current files give the set's labels back); the
refit files carry the frozen fit as evidence_thresholds.csv holds it and every change is
read off the two; and the subsets are E11's fits, the one-gene-out ones E3's own folds.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

PHASE1 = Path(__file__).resolve().parents[1]
REPO = PHASE1.parent
OUT = PHASE1 / "reports/evidence/suppl_analyses"
EV = PHASE1 / "reports/evidence"

if str(PHASE1) not in sys.path:
    sys.path.insert(0, str(PHASE1))

from src import config as C  # noqa: E402
from src import evid_common as K  # noqa: E402
from src import evid_suppl_c3 as C3  # noqa: E402

need = pytest.mark.skipif(not (OUT / "c3_intermediate.csv").exists(),
                          reason="run the pipeline first")
KEYS = ["tool", "stratum", "side", "tier"]


@pytest.fixture(scope="module")
def df() -> pd.DataFrame:
    with pytest.MonkeyPatch.context() as mp:
        mp.chdir(PHASE1)                     # the stage's paths are from phase1/
        return K.load_set()


def _eq(a: pd.Series, b: pd.Series) -> bool:
    a, b = a.astype("Float64"), b.astype("Float64")
    return bool(((a == b) | (a.isna() & b.isna())).fillna(False).all())


def test_a_label_file_joins_as_attach_labels_joins_it(df):
    for gene, spec in C.ASSAY_LABELS.items():
        got = C3.labels_from_file(df, gene, PHASE1 / C.ASSAY_LABEL_DIR / spec["file"])
        assert _eq(got, df.loc[df.gene == gene, "y_assay"]), gene


def test_the_intermediate_variants_are_the_337_and_nothing_else_moves(df):
    for value in (1.0, 0.0):
        y = C3.intermediate_labels(df, value)
        moved = ~((y == df.y_assay) | (y.isna() & df.y_assay.isna())).fillna(False)
        assert moved.sum() == 337
        assert (df.loc[moved, "y_assay_source"] == C3.INDETERMINATE).all()
        assert (df.loc[moved, "stratum"] != "pm12").all()
        assert (y[moved] == value).all()


def test_only_palb2_and_bard1_change_under_the_superseded_labels(df, monkeypatch):
    monkeypatch.chdir(PHASE1)
    y = C3.superseded_labels(df)
    other = ~df.gene.isin(C3.SUPERSEDED_GENES)
    assert _eq(y[other], df.loc[other, "y_assay"])
    assert not _eq(y[~other], df.loc[~other, "y_assay"])


@need
def test_the_refit_files_carry_the_frozen_fit():
    e = K.read_back(EV / "evidence_thresholds.csv")
    e = e[e.stratum.isin(C3.BANDS)].set_index(["tool", "stratum", "tier"])
    for name, n_analyses in (("c3_intermediate.csv", 2), ("c3_superseded_labels.csv", 1)):
        t = K.read_back(OUT / name)
        assert len(t) == n_analyses * 14 * 3 * 2 * 4
        assert not t.duplicated(["analysis"] + KEYS).any()
        for r in t.itertuples():
            f = e.loc[(r.tool, r.stratum, r.tier)]
            thr = f[f"{r.side}_threshold_insample"]
            assert (np.isnan(thr) and np.isnan(r.frozen_threshold)) or thr == r.frozen_threshold
            assert r.frozen_reached == f[f"{r.side}_threshold_reachable"]
            assert (r.frozen_n_pos, r.frozen_n_neg) == (f.n_pos, f.n_neg)
            assert r.cut_lr == f["path_cut_lr" if r.side == "pp3" else "ben_cut_lr"]


@need
def test_every_change_is_read_off_the_two_fits():
    for name in ("c3_intermediate.csv", "c3_superseded_labels.csv"):
        t = K.read_back(OUT / name)
        assert (t.reached == np.isfinite(t.threshold)).all()
        want = np.where(t.reached & ~t.frozen_reached, "gained",
                        np.where(~t.reached & t.frozen_reached, "lost", "unchanged"))
        assert (t.tier_change == want).all()
        same = (t.threshold == t.frozen_threshold) | (t.threshold.isna() & t.frozen_threshold.isna())
        assert (t.threshold_changed == ~same).all()


@need
def test_the_intermediate_refits_add_the_337(df):
    t = K.read_back(OUT / "c3_intermediate.csv")
    ind = df[(df.y_assay_source == C3.INDETERMINATE) & (df.stratum != "pm12")]
    for r in t.drop_duplicates(["analysis", "stratum"]).itertuples():
        want = len(K.stratum_frame(ind, r.stratum))
        assert (r.n_added, r.n_removed, r.n_switched) == (want, 0, 0)
    s = t[t.stratum == "s3_50"].drop_duplicates("analysis")
    assert (s.n_added == 337).all()
    # every column's added variants land on the side the analysis names
    dam = t[t.analysis == "intermediate as damaging"]
    assert (dam.n_pos >= dam.frozen_n_pos).all() and (dam.n_neg == dam.frozen_n_neg).all()
    nor = t[t.analysis == "intermediate as normal"]
    assert (nor.n_neg >= nor.frozen_n_neg).all() and (nor.n_pos == nor.frozen_n_pos).all()


@need
def test_the_subsets_are_e11s_fits_and_the_one_gene_out_ones_e3s_folds():
    wo = K.read_back(OUT / "c3_without_bard1_palb2.csv")
    e = K.read_back(EV / "evidence_thresholds.csv").set_index(["tool", "stratum", "tier"])
    folds = K.read_back(EV / "evidence_thresholds_logo_folds.csv")
    folds = folds[folds.side == "pp3"].set_index(["tool", "stratum", "tier", "heldout_gene"])
    assert len(wo) == 24 * 4
    for r in wo.itertuples():
        ref = e.loc[(r.tool, r.stratum, r.tier), "pp3_threshold_insample"]
        assert (np.isnan(ref) and np.isnan(r.threshold_every_gene)) or ref == r.threshold_every_gene
        genes = set(r.training_genes.split(";"))
        if r.subset == C3.ALL_GENES:
            assert r.n_genes == (6 if r.stratum == "s11_50" else 7)
            continue
        drop = dict(C3.SUBSETS)[r.subset]
        assert not genes & drop and r.n_genes == len(genes)
        if len(drop) == 1:
            f = folds.loc[(r.tool, r.stratum, r.tier, next(iter(drop))), "threshold_from_6_genes"]
            assert (np.isnan(f) and np.isnan(r.threshold)) or f == r.threshold
        assert r.tier_reached == bool(np.isfinite(r.threshold))


def test_the_stage_runs_after_c2():
    spec = importlib.util.spec_from_file_location(
        "reproduce_evidence", REPO / "scripts" / "reproduce_evidence.py")
    E = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(E)
    mods = [m for m, _, _ in E.STAGES]
    assert mods.index("src.evid_suppl_c3") == mods.index("src.evid_suppl_c2") + 1
