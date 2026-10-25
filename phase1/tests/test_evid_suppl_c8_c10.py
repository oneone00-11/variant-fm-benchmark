"""E19 C8, C9, C10: read from existing outputs.

The distal band's composition is the analysis set's; every subset row is E11's row and
the one-gene-out ones E3's folds; the external genes' two groups add up to the seven
external genes' published cells; and every MCC and F1 is the standard formula on E15's
counts.
"""
from __future__ import annotations

import importlib.util
import math
import re
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
from src import evid_suppl_c8_c10 as M  # noqa: E402

need = pytest.mark.skipif(not (OUT / "c8_without_bap1.csv").exists(),
                          reason="run the pipeline first")


def _dil() -> pd.DataFrame:
    return K.read_back(EV / "threshold_dilution.csv")


def _row(dil, tool, st, tier, train, held):
    hit = dil[(dil.tool == tool) & (dil.stratum == st) & (dil.tier == tier)
              & (dil.training_genes.map(lambda t: set(t.split(";"))) == set(train))
              & (dil.heldout_gene == held)]
    assert len(hit) == 1
    return hit.iloc[0]


def _same(a, b) -> bool:
    return (pd.isna(a) and pd.isna(b)) or a == b


@need
def test_the_distal_composition_is_the_sets():
    t = K.read_back(OUT / "c8_without_bap1.csv")
    s = pd.read_parquet(PHASE1 / "data/evidence/analysis_set_v1.parquet")
    d = s[(s.stratum == "s11_50") & s.y_assay.notna()]
    c = t[t.part == "composition"].set_index("gene")
    assert set(c.index) == set(d.gene)
    for g, r in c.iterrows():
        x = d[d.gene == g]
        assert (r.n_labelled, r.n_damaging) == (len(x), int((x.y_assay == 1).sum()))
    assert c.share_of_band_damaging.sum() == pytest.approx(1.0)


@need
def test_without_bap1_is_e11s_subset_and_e3s_fold():
    t = K.read_back(OUT / "c8_without_bap1.csv")
    dil, folds = _dil(), K.read_back(EV / "evidence_thresholds_logo_folds.csv")
    w = t[(t.part == "thresholds") & (t.subset == "without BAP1")]
    assert len(w) == 4 * 2 * 2
    for r in w.itertuples():
        train = r.training_genes.split(";")
        assert "BAP1" not in train and r.heldout_gene == "BAP1"
        d = _row(dil, r.tool, r.stratum, r.tier, train, "BAP1")
        assert _same(d.threshold, r.threshold) and _same(d.heldout_lr, r.heldout_lr)
        f = folds[(folds.tool == r.tool) & (folds.stratum == r.stratum) & (folds.tier == r.tier)
                  & (folds.side == "pp3") & (folds.heldout_gene == "BAP1")].iloc[0]
        assert _same(f.threshold_from_6_genes, r.threshold) and _same(f.heldout_lr, r.heldout_lr)


@need
def test_the_family_splits_are_e11s_rows():
    t = K.read_back(OUT / "c9_family_holdout.csv")
    a = t[t.part == "(a) family and function held out"]
    dil = _dil()
    for r in a.itertuples():
        train = set(r.training_genes.split(";"))
        if r.split.startswith("homologous"):
            assert train == set(M.HR) - ({"BRCA2"} if r.stratum == "s11_50" else set())
            assert r.heldout_gene in M.FUNCTION
        else:
            assert train == set(M.FUNCTION) and r.heldout_gene in M.HR
        d = _row(dil, r.tool, r.stratum, r.tier, train, r.heldout_gene)
        assert _same(d.threshold, r.threshold) and _same(d.heldout_lr, r.heldout_lr)
    # five held out per stratum and cell where BRCA2 has variants, four where it has none
    n = a[a.split.str.startswith("VHL")].groupby(["tool", "stratum", "tier"]).size()
    assert all(v == (4 if k[1] == "s11_50" else 5) for k, v in n.items())


@need
def test_the_two_external_groups_add_up_to_the_seven():
    t = K.read_back(OUT / "c9_family_holdout.csv", keep_default_na=False)
    b = t[t.part == "(b) external genes"]
    cell = re.compile(r"^(\d+)/(\d+)(?: \((\d+)\))?$")
    key = ["tool", "stratum_label", "side", "tier"]
    groups = {s: g.set_index(key)["genes_reaching_the_tier"] for s, g in b.groupby("split")}
    seven = groups.pop("Seven external genes")
    assert len(groups) == 2
    (_, x), (_, y) = groups.items()
    for k, v in seven.items():
        m, mx, my = cell.match(v), cell.match(x[k]), cell.match(y[k])
        if m and mx and my:
            got = [int(a or 0) + int(b_ or 0) for a, b_ in zip(mx.groups(), my.groups())]
            assert got == [int(z or 0) for z in m.groups()]
        else:
            assert x[k] == y[k] == v or not m


@need
def test_mcc_and_f1_are_the_standard_formulas_on_e15s_counts():
    head = [x for x in (OUT / "c10_confusion_mcc_f1.csv").read_text().splitlines()
            if x.startswith("#")]
    assert any("MCC = (TP x TN - FP x FN)" in x for x in head)
    assert any("F1 = 2 TP / (2 TP + FP + FN)" in x for x in head)
    t = K.read_back(OUT / "c10_confusion_mcc_f1.csv", comment="#")
    p = K.read_back(EV / "practical_metrics_by_threshold.csv")
    p = p[p.threshold_source.isin(["in_sample", "logo"])]
    assert len(t) == len(p)
    for r, q in zip(t.itertuples(), p.itertuples()):
        assert (r.tool, r.stratum, r.threshold_kind, r.threshold_source) == (
            q.tool, q.stratum, q.threshold_kind, q.threshold_source)
        want = ((q.n_pos_above, q.n_neg_above, q.n_pos_below, q.n_neg_below) if q.side == "pp3"
                else (q.bp4_n_neg_at_or_below, q.bp4_n_pos_at_or_below, q.bp4_n_neg_above,
                      q.bp4_n_pos_above))
        assert all(_same(a, b) for a, b in zip((r.tp, r.fp, r.fn, r.tn), want))
        if not any(pd.isna(w) for w in want):
            tp, fp, fn, tn = want
            den = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
            if den:
                assert r.mcc == pytest.approx((tp * tn - fp * fn) / math.sqrt(den), rel=1e-12)
            assert r.f1 == pytest.approx(2 * tp / (2 * tp + fp + fn), rel=1e-12)


def test_the_stage_runs_after_c5():
    spec = importlib.util.spec_from_file_location(
        "reproduce_evidence", REPO / "scripts" / "reproduce_evidence.py")
    E = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(E)
    mods = [m for m, _, _ in E.STAGES]
    assert mods.index("src.evid_suppl_c8_c10") > mods.index("src.evid_suppl_c5")
