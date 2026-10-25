"""E19 C4: E3's in-sample thresholds at the priors 0.05, 0.10 and 0.185.

Every edge keeps its posterior across the priors (the fixed-posterior rule written in
the files' header); the fit at 0.10 is the frozen one, cell for cell; and every change
and every Strong flag is read off the two fits.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest
import yaml

PHASE1 = Path(__file__).resolve().parents[1]
REPO = PHASE1.parent
OUT = PHASE1 / "reports/evidence/suppl_analyses"
EV = PHASE1 / "reports/evidence"

if str(PHASE1) not in sys.path:
    sys.path.insert(0, str(PHASE1))

from src import evid_common as K  # noqa: E402
from src import evid_interval_lr as E3  # noqa: E402
from src import evid_suppl_c4 as C4  # noqa: E402

need = pytest.mark.skipif(not (OUT / "c4_prior_sensitivity.csv").exists(),
                          reason="run the pipeline first")


def _cuts() -> dict:
    cfg = yaml.safe_load((PHASE1 / E3.CONFIG_PATH).read_text())
    p, b = K.acmg_bands(cfg)
    return {"path": p, "ben": b}


def test_every_edge_keeps_its_posterior():
    cuts = _cuts()
    for p in C4.PRIORS:
        e = C4.edges_at(p, cuts)
        for key in ("path", "ben"):
            for t in E3.TIERS:
                assert e[key][t] * C4._odds(p) == pytest.approx(
                    cuts[key][t] * C4._odds(C4.BASE_PRIOR), rel=1e-12)
    assert C4.edges_at(C4.BASE_PRIOR, cuts) == cuts


@need
def test_the_files_carry_the_formula_and_the_edges():
    cuts = _cuts()
    for name in ("c4_band_edges.csv", "c4_prior_sensitivity.csv"):
        head = [x for x in (OUT / name).read_text().splitlines() if x.startswith("#")]
        assert any(C4.FORMULA in x for x in head)
        for p in C4.PRIORS:
            e = C4.edges_at(p, cuts)
            assert any(x.endswith(", ".join(repr(e["path"][t]) for t in E3.TIERS)) for x in head)
    be = K.read_back(OUT / "c4_band_edges.csv", comment="#")
    assert len(be) == 3 * 2 * 4
    for (side, tier), g in be.groupby(["side", "tier"]):
        assert g.posterior.nunique() == 1 or np.ptp(g.posterior) < 1e-15


@need
def test_the_fit_at_0_10_is_the_frozen_one():
    t = K.read_back(OUT / "c4_prior_sensitivity.csv", comment="#")
    e = K.read_back(EV / "evidence_thresholds.csv").set_index(["tool", "stratum", "tier"])
    assert len(t) == 3 * 2 * len(e)
    base = t[t.prior == C4.BASE_PRIOR]
    for r in base.itertuples():
        f = e.loc[(r.tool, r.stratum, r.tier)]
        thr = f[f"{r.side}_threshold_insample"]
        assert (np.isnan(thr) and np.isnan(r.threshold)) or thr == r.threshold
        assert r.reached == f[f"{r.side}_threshold_reachable"]
        assert (r.n_pos, r.n_neg) == (f.n_pos, f.n_neg)
    assert (base.tier_change == "unchanged").all() and not base.threshold_changed.any()


@need
def test_every_change_and_strong_flag_is_read_off_the_fits():
    t = K.read_back(OUT / "c4_prior_sensitivity.csv", comment="#")
    assert (t.reached == np.isfinite(t.threshold)).all()
    want = np.where(t.reached & ~t.reached_at_0_10, "gained",
                    np.where(~t.reached & t.reached_at_0_10, "lost", "unchanged"))
    assert (t.tier_change == want).all()
    key = ["prior", "tool", "stratum", "side"]
    strong = t[t.tier == "strong"].set_index(key)
    for k, g in t.groupby(key):
        s = strong.loc[k]
        assert (g.strong_reached == s.reached).all()
        assert (g.strong_reached_at_0_10 == s.reached_at_0_10).all()
        if s.reached_at_0_10:
            assert (g.strong_still_reached == s.reached).all()
        else:
            assert g.strong_still_reached.isna().all()


def test_the_stage_runs_after_c3():
    spec = importlib.util.spec_from_file_location(
        "reproduce_evidence", REPO / "scripts" / "reproduce_evidence.py")
    E = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(E)
    mods = [m for m, _, _ in E.STAGES]
    assert mods.index("src.evid_suppl_c4") == mods.index("src.evid_suppl_c3") + 1
