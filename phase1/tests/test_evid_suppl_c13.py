"""E20 C13: the VCEP rules applied to this study's data.

Every evaluable row's counts, ratio and tier at the point estimate follow from the analysis
set, the rule's cut point and column, by E2's band ratio; a band without labelled variants
says so; BP7 follows each specification's own definition (E20b.3): ENIGMA's rows are
c6_bp4_bp7.csv's condition (a), PALB2's and VHL's recount from the analysis set; every PP3
and BP4 rule of the rules file has its three bands; and the header names the sources, the
column, ENIGMA's setting, PhyloP's source, what is not evaluated and the labels.
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
def test_a_band_without_labelled_variants_says_so():
    b = _out().query("part == 'band_lr'")
    empty = b[(b["n_damaging"] == 0) & (b["n_normal"] == 0)]
    assert len(empty) == 2 and set(empty["gene"]) == {"BRCA2"} and set(empty["stratum"]) == {"s11_50"}
    assert (empty["reason"] == "no labelled variants at this distance in this gene's assay").all()
    rest = b[(b["status"] == "not evaluable") & ~b.index.isin(empty.index)]
    assert not rest["reason"].str.startswith("no labelled").any()


def _bp7() -> pd.DataFrame:
    return _out().query("part == 'bp7'")


def _labelled() -> pd.DataFrame:
    from src import evid_suppl_c6_c7 as C6C7
    return C6C7.labelled(pd.read_parquet(PHASE1 / "data/evidence/analysis_set_v1.parquet"))


def _row_counts(r) -> tuple:
    return (r.n_at_positions, r.n_damaging_at_positions, r.n_normal_at_positions, r.n_meeting,
            r.n_damaging, r.n_normal)


def _want(at: pd.DataFrame, meet: pd.DataFrame) -> tuple:
    return (len(at), int((at["y_assay"] == 1).sum()), int((at["y_assay"] == 0).sum()), len(meet),
            int((meet["y_assay"] == 1).sum()), int((meet["y_assay"] == 0).sum()))


@need
def test_every_bp7_rule_has_its_rows():
    b = _bp7()
    assert b.groupby("gene").size().to_dict() == {"BRCA1": 3, "BRCA2": 3, "PALB2": 6, "VHL": 3}


@need
def test_enigmas_bp7_rows_are_c6s_condition_a():
    c6 = K.read_back(PHASE1 / "reports/evidence/suppl_analyses/c6_bp4_bp7.csv", comment="#")
    c6 = c6[c6["condition"].str.startswith("(a)")].set_index(["gene", "side"])
    rows = _bp7()[_bp7()["gene"].isin(["BRCA1", "BRCA2"])]
    assert len(rows) == 6
    for r in rows.itertuples():
        c = c6.loc[(r.gene, r.side)]
        assert _row_counts(r) == tuple(int(c[f]) for f in (
            "n_at_positions", "n_damaging_at_positions", "n_normal_at_positions", "n_meeting",
            "n_damaging", "n_normal"))
        assert np.array_equal([r.share_damaging], [c["share_damaging"]], equal_nan=True)


@need
def test_palb2s_bp7_rows_follow_from_the_data():
    p = _labelled().query("gene == 'PALB2'")
    o = p["signed_offset"]
    sides = {"donor, +8 and beyond": o >= 8, "acceptor, -22 and beyond": o <= -22}
    sides["both sides"] = sides["donor, +8 and beyond"] | sides["acceptor, -22 and beyond"]
    rows = _bp7().query("gene == 'PALB2'")
    assert sorted(rows["bp7_condition"].str[:4].tolist()) == ["(i) "] * 3 + ["(ii)"] * 3
    for r in rows.itertuples():
        at = p[sides[r.side]]
        meet = at if r.bp7_condition.startswith("(i) ") else at[at["spliceai_walker"] <= 0.1]
        assert _row_counts(r) == _want(at, meet)


@need
def test_vhls_bp7_rows_follow_from_the_data():
    v = _labelled().query("gene == 'VHL'")
    off = v["signed_offset"].abs()
    rows = _bp7().query("gene == 'VHL'")
    assert sorted(rows["stratum"]) == sorted(BANDS)
    for r in rows.itertuples():
        lo, hi = {"s3_10": (3, 10), "s11_50": (11, 50), "s3_50": (3, 50)}[r.stratum]
        band = off.between(lo, hi)
        at = v[band & v["phylop"].notna()]
        meet = at[(at["spliceai_walker"] <= 0.1) & (at["phylop"] <= 0.2)]
        assert _row_counts(r) == _want(at, meet)
        assert r.n_phylop_missing == int((band & v["phylop"].isna()).sum())


@need
def test_the_header_names_sources_column_and_what_is_not_evaluated():
    head = "\n".join(x for x in OUT.read_text().splitlines() if x.startswith("#"))
    for gn in ("GN092", "GN097", "GN077", "GN078"):
        assert gn in head
    assert "spliceai_walker" in head and "VarSeak" in head
    assert "SGE labels" in head and "circularity" in head
    assert "Appendix J" in head and "Max Distance 10 000" in head and "'unmasked'" in head
    assert "UCSC phyloP100way, hg38" in head and "names no version" in head
    assert "+8 and beyond" in head and "-22 and beyond" in head and "condition (a)" in head


def test_the_stage_runs_after_c12():
    spec = importlib.util.spec_from_file_location("reproduce_evidence", REPO / "scripts/reproduce_evidence.py")
    E = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(E)
    mods = [m for m, _, _ in E.STAGES]
    assert mods.index("src.evid_suppl_c13") > mods.index("src.evid_suppl_c12")
