"""E20 C13: the VCEP rules applied to this study's data.

Every evaluable row's counts, ratio and tier at the point estimate follow from the analysis
set, the rule's cut point and column, by E2's band ratio; the BP7 rows are c6_bp4_bp7.csv's
own; every PP3 and BP4 rule of the rules file has its three bands; and the header names the
sources, the column, what is not evaluated and the labels.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

PHASE1 = Path(__file__).resolve().parents[1]
REPO = PHASE1.parent
OUT = PHASE1 / "reports/evidence/suppl_analyses/c13_vcep_band_lr.csv"

if str(PHASE1) not in sys.path:
    sys.path.insert(0, str(PHASE1))

from src import evid_common as K  # noqa: E402

need = pytest.mark.skipif(not OUT.exists(), reason="run the pipeline first")
BANDS = {"s3_10": ["s3_10"], "s11_50": ["s11_50"], "s3_50": ["s3_10", "s11_50"]}


def _out() -> pd.DataFrame:
    return K.read_back(OUT, comment="#")


@need
def test_every_rule_has_its_three_bands():
    rules = json.loads((PHASE1 / "data/evidence/vcep/vcep_rules.json").read_text())["rules"]
    o = _out()
    b = o[o["part"] == "band_lr"]
    for r in rules:
        if r["criterion"] in ("PP3", "BP4"):
            rows = b[(b["gene"] == r["gene"]) & (b["criterion"] == r["criterion"])]
            assert sorted(rows["stratum"]) == sorted(BANDS)


@need
def test_the_evaluable_rows_follow_from_the_data():
    df = pd.read_parquet(PHASE1 / "data/evidence/analysis_set_v1.parquet")
    cfg = yaml.safe_load((PHASE1 / "config/walker2023.yaml").read_text())
    path_bands, ben_bands = K.acmg_bands(cfg)
    o = _out()
    ok = o[(o["part"] == "band_lr") & (o["status"] == "ok")]
    assert len(ok)
    for r in ok.itertuples():
        sub = df[(df["gene"] == r.gene) & df["stratum"].isin(BANDS[r.stratum]) & df["y_assay"].notna()
                 & df[r.score_column].notna()]
        y, s = sub["y_assay"].to_numpy(float), sub[r.score_column].to_numpy(float)
        cut = float(r.rule.split()[-1])
        side = "upper" if r.criterion == "PP3" else "lower"
        assert (int((y == 1).sum()), int((y == 0).sum())) == (r.n_damaging, r.n_normal)
        lr = K.band_lr(y, s, cut, side)
        assert lr == r.lr
        if side == "upper":
            assert r.tier_at_point == K.tier_of(lr, path_bands, "pathogenic")
        else:
            assert r.tier_at_point == K.tier_of(K.band_lr_benign_bound(y, s, cut), ben_bands, "benign")


@need
def test_the_bp7_rows_are_c6s():
    c6 = K.read_back(PHASE1 / "reports/evidence/suppl_analyses/c6_bp4_bp7.csv", comment="#")
    c6 = c6[c6["condition"].str.startswith("(a)")].set_index(["gene", "side"])
    o = _out()
    for r in o[o["part"] == "bp7_from_c6"].itertuples():
        c = c6.loc[(r.gene, r.side)]
        assert (r.n_at_positions, r.n_meeting, r.n_damaging, r.n_normal) == (
            int(c["n_at_positions"]), int(c["n_meeting"]), int(c["n_damaging"]), int(c["n_normal"]))


@need
def test_the_header_names_sources_column_and_what_is_not_evaluated():
    head = "\n".join(x for x in OUT.read_text().splitlines() if x.startswith("#"))
    for gn in ("GN092", "GN097", "GN077", "GN078"):
        assert gn in head
    assert "spliceai_walker" in head and "VarSeak" in head
    assert "SGE labels" in head and "circularity" in head


def test_the_stage_runs_after_c12():
    spec = importlib.util.spec_from_file_location("reproduce_evidence", REPO / "scripts/reproduce_evidence.py")
    E = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(E)
    mods = [m for m, _, _ in E.STAGES]
    assert mods.index("src.evid_suppl_c13") > mods.index("src.evid_suppl_c12")
