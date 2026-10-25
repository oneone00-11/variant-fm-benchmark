"""E19.10: E3's in-sample thresholds refitted on the seven genes and the five new ones.

The twelve-gene frame is the analysis set's seven genes and the five genes' scored rows
under their primary labels, nothing else; the seven-gene columns are the frozen fit;
every change is read off the two; NT and the fusion are not refitted and say why.
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

from src import evid_common as K  # noqa: E402
from src import evid_external as X  # noqa: E402
from src import evid_suppl_twelve_gene as TG  # noqa: E402

need = pytest.mark.skipif(not (OUT / "twelve_gene_refit.csv").exists(),
                          reason="run the pipeline first")


def test_the_frame_is_the_seven_genes_and_the_five_under_their_primary_labels(monkeypatch):
    monkeypatch.chdir(PHASE1)
    df = K.load_set()
    tools = [t for t in K.panel_of(df) if t not in TG.NOT_REFITTED]
    assert len(tools) == 12 and "nt" not in tools
    d = TG.twelve_gene_frame(df, tools)
    assert list(pd.unique(d.gene)) == list(pd.unique(df.gene)) + X.NEW_GENES
    seven = d.iloc[:len(df)]
    assert (seven.gene.to_numpy() == df.gene.to_numpy()).all()
    same = (seven.y_assay.to_numpy() == df.y_assay.to_numpy()) | (
        seven.y_assay.isna().to_numpy() & df.y_assay.isna().to_numpy())
    assert same.all()
    s = pd.read_parquet(PHASE1 / X.NEWGENES_SCORED)
    for g in X.NEW_GENES:
        a, b = d[d.gene == g], s[s.gene == g]
        assert len(a) == len(b) and set(a.stratum) <= {"s3_10", "s11_50"}
        assert ((a.y_assay.to_numpy() == b.label.to_numpy()) |
                (a.y_assay.isna().to_numpy() & b.label.isna().to_numpy())).all()


@need
def test_the_seven_gene_columns_are_the_frozen_fit():
    t = K.read_back(OUT / "twelve_gene_refit.csv")
    e = K.read_back(EV / "evidence_thresholds.csv")
    e = e[e.stratum.isin(["s3_10", "s11_50", "s3_50"])].set_index(["tool", "stratum", "tier"])
    assert len(t) == 14 * 3 * 2 * 4
    for r in t.itertuples():
        f = e.loc[(r.tool, r.stratum, r.tier)]
        thr = f[f"{r.side}_threshold_insample"]
        assert (np.isnan(thr) and np.isnan(r.threshold_7)) or thr == r.threshold_7
        assert r.reached_7 == f[f"{r.side}_threshold_reachable"]
        assert (r.n_pos_7, r.n_neg_7) == (f.n_pos, f.n_neg)


@need
def test_every_change_is_read_off_the_two_fits():
    t = K.read_back(OUT / "twelve_gene_refit.csv")
    fit = t[t.status_12 != "not refitted"]
    assert len(fit) == 12 * 3 * 2 * 4
    assert (fit.reached_12 == np.isfinite(fit.threshold_12)).all()
    d = fit.threshold_12 - fit.threshold_7
    assert ((fit.abs_change == d) | (fit.abs_change.isna() & d.isna())).all()
    r = d / fit.threshold_7.abs()
    assert ((fit.rel_change == r) | (fit.rel_change.isna() & r.isna())).all()
    want = np.where(fit.reached_12 & ~fit.reached_7.astype(bool), "gained",
                    np.where(~fit.reached_12 & fit.reached_7.astype(bool), "lost", "unchanged"))
    assert (fit.tier_change == want).all()
    assert (fit.n_genes_12 == fit.n_genes_7 + 5).all()
    rest = t[t.status_12 == "not refitted"]
    assert set(rest.tool) == set(TG.NOT_REFITTED) and rest.reason.notna().all()


def test_the_stage_runs_last():
    spec = importlib.util.spec_from_file_location(
        "reproduce_evidence", REPO / "scripts" / "reproduce_evidence.py")
    E = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(E)
    mods = [m for m, _, _ in E.STAGES]
    assert mods.index("src.evid_suppl_twelve_gene") > mods.index("src.evid_suppl_c5")
