"""E19 C3, BRCA2 against a second dataset: the concordance and the two refits.

The 3 x 3 table counts the 201 variants the two sets share at 3-50 nt; the refits change
exactly the labels their definitions name; every change is read off the frozen fit and
the refit; the 11-50 nt band cannot change; and both files say in their headers which
second set, which range and which rows of Table 5 the author's decisions concern.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

PHASE1 = Path(__file__).resolve().parents[1]
REPO = PHASE1.parent
OUT = PHASE1 / "reports/evidence/suppl_analyses"
EV = PHASE1 / "reports/evidence"

if str(PHASE1) not in sys.path:
    sys.path.insert(0, str(PHASE1))

from src import evid_common as K  # noqa: E402
from src import evid_suppl_c3_brca2 as B  # noqa: E402

need = pytest.mark.skipif(not (OUT / "c3_brca2_concordance.csv").exists(),
                          reason="run the pipeline first")


def _read(name: str):
    return K.read_back(OUT / name, comment="#")


@need
def test_the_concordance_counts_the_shared_variants():
    c = _read("c3_brca2_concordance.csv")
    assert len(c) == 9 and set(c.huang_class) == set(B.ORDER) == set(c.second_set_class)
    n = int(c.n.sum())
    assert n == int(c.n_shared_3_50.iloc[0]) == 201
    agree = int(c[c.huang_class == c.second_set_class].n.sum())
    assert agree == int(c.n_agree.iloc[0])
    assert np.allclose(c.share_of_shared, c.n / n)


@need
def test_each_refit_changes_exactly_the_labels_it_names():
    c = _read("c3_brca2_concordance.csv").set_index(["huang_class", "second_set_class"]).n
    labelled = ("abnormal", "normal")
    differ_labelled = sum(int(c[(h, s)]) for h in labelled for s in B.ORDER if h != s)
    switched = int(c[("abnormal", "normal")] + c[("normal", "abnormal")])
    to_inter = sum(int(c[(h, "intermediate")]) for h in labelled)
    from_inter = sum(int(c[("intermediate", s)]) for s in labelled)
    r = _read("c3_brca2_second_set.csv")
    s = r[r.stratum == "s3_50"].drop_duplicates("analysis").set_index("analysis")
    assert tuple(s.loc[B.AGREE, ["n_added", "n_removed", "n_switched"]]) == (0, differ_labelled, 0)
    assert tuple(s.loc[B.SECOND, ["n_added", "n_removed", "n_switched"]]) == (from_inter, to_inter,
                                                                            switched)


@need
def test_every_change_is_read_off_the_two_fits_and_the_distal_band_cannot_change():
    r = _read("c3_brca2_second_set.csv")
    e = K.read_back(EV / "evidence_thresholds.csv").set_index(["tool", "stratum", "tier"])
    assert len(r) == 2 * 14 * 3 * 2 * 4
    for x in r.itertuples():
        f = e.loc[(x.tool, x.stratum, x.tier)]
        thr = f[f"{x.side}_threshold_insample"]
        assert (np.isnan(thr) and np.isnan(x.frozen_threshold)) or thr == x.frozen_threshold
    want = np.where(r.reached & ~r.frozen_reached, "gained",
                    np.where(~r.reached & r.frozen_reached, "lost", "unchanged"))
    assert (r.tier_change == want).all()
    distal = r[r.stratum == "s11_50"]
    assert (distal.tier_change == "unchanged").all() and not distal.threshold_changed.any()
    assert (distal[["n_added", "n_removed", "n_switched"]] == 0).all().all()


@need
def test_both_headers_name_the_set_its_range_and_the_rows_of_table_5():
    for name in ("c3_brca2_concordance.csv", "c3_brca2_second_set.csv"):
        head = "\n".join(x for x in (OUT / name).read_text().splitlines() if x.startswith("#"))
        assert "urn:mavedb:00001242-a-1" in head
        assert "the donor side +3 ... +9 and the acceptor side -3 ... -7" in head
        assert "201 by Table 5" in head and "198, by MaveDB's score table, is unchanged" in head
        assert "c.9505A>C/G/T thus stands twice" in head and "lacks all six rows" in head
        assert "c.8610A>C at g.32371078A>C (GRCh38 G)" in head
        assert head.count("c. miswritten") == 3


def test_the_stage_runs_after_c1_low():
    spec = importlib.util.spec_from_file_location(
        "reproduce_evidence", REPO / "scripts" / "reproduce_evidence.py")
    E = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(E)
    mods = [m for m, _, _ in E.STAGES]
    assert mods.index("src.evid_suppl_c3_brca2") > mods.index("src.evid_suppl_c1_low")
