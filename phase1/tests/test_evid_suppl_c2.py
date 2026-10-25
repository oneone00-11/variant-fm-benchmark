"""E19 C2: ClinVar's coverage of the analysis set by signed offset and by band.

Every one of the 8,853 variants is counted once, at its signed offset, whose magnitude
is the set's intron_offset_abs; the counts and shares are recomputed from the set; the
band summary adds up to the positions; and the stage runs after every existing stage.
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
SET = PHASE1 / "data/evidence/analysis_set_v1.parquet"

if str(PHASE1) not in sys.path:
    sys.path.insert(0, str(PHASE1))

from src import evid_common as K  # noqa: E402
from src import evid_suppl_c2 as C2  # noqa: E402

need = pytest.mark.skipif(not (OUT / "c2_clinvar_coverage_by_position.csv").exists(),
                          reason="run the pipeline first")


def _set() -> pd.DataFrame:
    d = pd.read_parquet(SET)
    d["signed_offset"] = d["hgvs_c"].map(C2.signed_offset)
    return d


def test_the_signed_offset_keeps_the_sets_magnitude():
    d = _set()
    assert len(d) == 8853
    assert (d.signed_offset.abs() == d.intron_offset_abs).all()
    assert C2.signed_offset("c.123-8T>G") == -8 and C2.signed_offset("c.9501+5G>C") == 5
    assert C2.signed_offset("c.-19+2T>A") == 2 and C2.signed_offset("c.*7-3C>T") == -3


@need
def test_every_count_is_the_sets():
    d = _set()
    pos = K.read_back(OUT / "c2_clinvar_coverage_by_position.csv")
    assert len(pos) == 8 * 100
    allg = pos[pos.gene == C2.ALL]
    assert allg.n_variants.sum() == len(d)
    for r in pos.itertuples():
        g = d if r.gene == C2.ALL else d[d.gene == r.gene]
        c = g[g.signed_offset == r.signed_offset]
        assert r.n_variants == len(c)
        assert r.n_recorded == int((c.clinvar_arm != "unrecorded").sum())
        assert r.n_classified == int((c.clinvar_arm == "classified").sum())
        for cls, short in C2.CLASSES:
            assert getattr(r, f"n_{short}") == int((c.clnsig_class == cls).sum())
        if len(c):
            assert r.share_recorded == r.n_recorded / len(c)
        else:
            assert np.isnan(r.share_recorded)
    # the five classes partition the recorded variants
    k = sum(allg[f"n_{s}"] for _, s in C2.CLASSES)
    assert (k == allg.n_recorded).all()


@need
def test_the_band_summary_adds_up_to_the_positions():
    pos = K.read_back(OUT / "c2_clinvar_coverage_by_position.csv")
    st = K.read_back(OUT / "c2_clinvar_coverage_by_stratum.csv")
    band_of = lambda o: "pm12" if abs(o) <= 2 else "s3_10" if abs(o) <= 10 else "s11_50"  # noqa: E731
    pos["stratum"] = pos.signed_offset.map(band_of)
    agg = pos.groupby(["gene", "stratum"])[["n_variants", "n_recorded", "n_classified"]].sum()
    for r in st.itertuples():
        a = agg.loc[(r.gene, r.stratum)]
        assert r.n_variants == a.n_variants
        assert r.n_unrecorded + r.n_recorded_unclassified + r.n_classified == r.n_variants
        assert r.n_classified == a.n_classified
        assert r.n_variants - r.n_unrecorded == a.n_recorded


def test_the_stage_runs_after_every_existing_stage():
    spec = importlib.util.spec_from_file_location(
        "reproduce_evidence", REPO / "scripts" / "reproduce_evidence.py")
    E = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(E)
    mods = [m for m, _, _ in E.STAGES]
    i = mods.index("src.evid_suppl_c2")
    assert i > mods.index("src.evid_alphagenome_definitions")
