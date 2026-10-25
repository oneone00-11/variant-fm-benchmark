"""E19 C6 and C7: BP4 at the BP7 positions, and SpliceAI's published cut points by position.

Every count is the analysis set's; condition (b) uses evidence_thresholds.csv's fitted
BP4 Supporting thresholds at full precision; the header names exactly the variants the
printed Table 3(b) values would leave out; and C7's ratios are E2's _band_stats on the
same variants.
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
from src import evid_suppl_c2 as C2  # noqa: E402
from src import evid_suppl_c6_c7 as M  # noqa: E402
from src import evid_walker_thresholds as W  # noqa: E402

need = pytest.mark.skipif(not (OUT / "c6_bp4_bp7.csv").exists(), reason="run the pipeline first")


def _set() -> pd.DataFrame:
    d = pd.read_parquet(PHASE1 / "data/evidence/analysis_set_v1.parquet")
    d = d.dropna(subset=["y_assay", M.SCORE]).copy()
    d["signed_offset"] = d["hgvs_c"].map(C2.signed_offset)
    return d


def _fitted() -> dict[str, float]:
    e = K.read_back(EV / "evidence_thresholds.csv")
    return {st: float(e[(e.tool == M.SCORE) & (e.stratum == st) & (e.tier == "supporting")]
                      .bp4_threshold_insample.iloc[0]) for st in ("s3_10", "s11_50")}


@need
def test_the_header_gives_the_fitted_thresholds_and_the_printed_values_difference():
    head = [x for x in (OUT / "c6_bp4_bp7.csv").read_text().splitlines() if x.startswith("#")]
    f = _fitted()
    assert any(f"3-10 nt {f['s3_10']!r}, 11-50 nt {f['s11_50']!r}" in x for x in head)
    d = _set()
    bp7 = d[(d.signed_offset >= 7) | (d.signed_offset <= -21)]
    band = np.where(bp7.signed_offset.abs() <= 10, "s3_10", "s11_50")
    fitted = pd.Series(band, index=bp7.index).map(f)
    printed = pd.Series(band, index=bp7.index).map({"s3_10": 0.0260, "s11_50": 0.00127})
    left = bp7[(bp7[M.SCORE] <= fitted) & (bp7[M.SCORE] > printed)]
    assert sorted(zip(left.gene, left.hgvs_c)) == [("BAP1", "c.255+9G>C"),
                                                   ("BARD1", "c.1810+43A>G")]
    assert (left.y_assay == 0).all()
    assert any("BAP1 c.255+9G>C (3-10 nt); BARD1 c.1810+43A>G (11-50 nt)" in x
               and "The damaging counts are the same" in x for x in head)


@need
def test_every_c6_count_is_the_sets():
    t = K.read_back(OUT / "c6_bp4_bp7.csv", comment="#")
    d, f = _set(), _fitted()
    assert len(t) == 2 * 3 * 8
    sides = dict(M.C6_SIDES)
    for r in t.itertuples():
        g = d if r.gene == M.ALL else d[d.gene == r.gene]
        g = g[sides[r.side](g.signed_offset)]
        if r.condition.startswith("(a)"):
            cut = 0.1
        else:
            cut = pd.Series(np.where(g.signed_offset.abs() <= 10, "s3_10", "s11_50"),
                            index=g.index).map(f)
        m = g[g[M.SCORE] <= cut]
        assert (r.n_at_positions, r.n_meeting) == (len(g), len(m))
        assert (r.n_damaging, r.n_normal) == (int((m.y_assay == 1).sum()), int((m.y_assay == 0).sum()))
    both = t[t.side == "both sides"].set_index(["condition", "gene"])
    parts = t[t.side != "both sides"].groupby(["condition", "gene"])[["n_meeting", "n_damaging"]].sum()
    assert (both[["n_meeting", "n_damaging"]].sort_index() == parts.sort_index()).all().all()


@need
def test_c7_is_e2s_band_statistics_by_position_group():
    t = K.read_back(OUT / "c7_by_position.csv")
    d = _set()
    assert list(t.group) == [g for _, g, _, _ in M.C7_GROUPS]
    for r, (_, _, lo, hi) in zip(t.itertuples(), M.C7_GROUPS):
        g = d[(d.signed_offset >= lo) & (d.signed_offset <= hi)]
        y, s = g.y_assay.to_numpy(float), g[M.SCORE].to_numpy(float)
        assert (r.n_damaging, r.n_normal) == (int((y == 1).sum()), int((y == 0).sum()))
        st = W._band_stats(y, s, 0.2, 0.1)
        if not st:
            assert r.status == "counts only"
            continue
        assert r.status == "ok"
        for k in ("lr_pp3", "lr_bp4"):
            v = getattr(r, k)
            assert (np.isnan(st[k]) and np.isnan(v)) or st[k] == pytest.approx(v, rel=1e-12)
            lo_, hi_ = getattr(r, f"{k}_lo"), getattr(r, f"{k}_hi")
            if np.isfinite(lo_):
                assert lo_ <= hi_


def test_the_stage_runs_after_the_twelve_gene_refit():
    spec = importlib.util.spec_from_file_location(
        "reproduce_evidence", REPO / "scripts" / "reproduce_evidence.py")
    E = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(E)
    mods = [m for m, _, _ in E.STAGES]
    assert mods.index("src.evid_suppl_c6_c7") > mods.index("src.evid_suppl_twelve_gene")
