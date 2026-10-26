"""E20 C12: the true positives' in-frame attribution beside the false positives'.

The false-positive rows are E6's and E2.7's summaries, cell for cell; the true-positive
rows count every variant of the subsets on their bands, label by label, and none of the
anchors; the shares and their Wilson intervals follow from the counts; and the header says
what E20.2 asks it to say.
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
OUT = PHASE1 / "reports/evidence/suppl_analyses/c12_inframe_tp_fp.csv"
EV = PHASE1 / "data/evidence"
REP = PHASE1 / "reports/evidence"

if str(PHASE1) not in sys.path:
    sys.path.insert(0, str(PHASE1))

from src import evid_common as K  # noqa: E402
from src import evid_practical_metrics as P  # noqa: E402

need = pytest.mark.skipif(not OUT.exists(), reason="run the pipeline first")
FRAMES = ["in_frame", "out_of_frame", "undetermined", "n"]


def _out() -> pd.DataFrame:
    return K.read_back(OUT, comment="#")


@need
def test_the_false_positives_are_e6s_and_e27s_summaries():
    o = _out()
    fp = o[o["group"] == "false_positive"]
    seven = K.read_back(REP / "inframe_attribution_summary.csv")
    for r in fp[fp["gene_set"] == "seven genes"].itertuples():
        s = seven[(seven["tool"] == r.tool) & (seven["stratum"] == r.stratum)].iloc[0]
        assert [int(s[f]) for f in FRAMES] == [getattr(r, f) for f in FRAMES]
    for gene in ("ddx3x", "tp53"):
        summ = K.read_back(REP / f"inframe_attribution_{gene}_summary.csv")
        rows = fp[fp["gene_set"] == gene.upper()]
        assert set(rows["tool"]) == {"spliceai"}
        for r in rows.itertuples():
            s = summ[(summ["label_definition"] == r.label_definition) & (summ["tool"] == r.tool)
                     & (summ["stratum"] == r.stratum)]
            want = [int(s[f].iloc[0]) for f in FRAMES] if len(s) else [0, 0, 0, 0]
            assert want == [getattr(r, f) for f in FRAMES]


@need
def test_the_true_positives_are_the_subsets_without_the_anchors():
    o = _out()
    tp = o[o["group"] == "true_positive"].set_index(["gene_set", "label_definition", "tool", "stratum"])
    for tool, name in (("spliceai", "inframe_subset_tp.parquet"), ("pangolin", "inframe_subset_pangolin_tp.parquet")):
        sub = pd.read_parquet(EV / name)
        for band, n in sub.groupby("stratum").size().items():
            assert tp.loc[("seven genes", "y_assay", tool, band), "n"] == n
    for gene, name in (("DDX3X", "ddx3x_inframe_subset_tp.parquet"), ("TP53", "tp53_inframe_subset_12nt_tp.parquet")):
        sub = pd.read_parquet(EV / "external" / name)
        labels = [c for c in sub.columns if c.startswith("y_")]
        for lab in labels:
            for band in ("s3_10", "s11_50"):
                n = int(((sub[lab] == 1) & (sub["stratum"] == band)).sum())
                assert tp.loc[(gene, lab, "spliceai", band), "n"] == n
    assert (tp[["in_frame", "out_of_frame", "undetermined"]].sum(axis=1) == tp["n"]).all()


@need
def test_the_shares_follow_from_the_counts():
    o = _out()
    for r in o.itertuples():
        if r.n == 0:
            assert np.isnan(r.in_frame_share) and np.isnan(r.in_frame_share_lo)
            continue
        assert r.in_frame_share == round(r.in_frame / r.n, 4)
        lo, hi = P._wilson(r.in_frame, r.n)
        assert (r.in_frame_share_lo, r.in_frame_share_hi) == (round(lo, 4), round(hi, 4))
        if r.in_frame + r.out_of_frame:
            assert r.resolved_in_frame_share == round(r.in_frame / (r.in_frame + r.out_of_frame), 4)


@need
def test_the_header_says_what_e20_2_asks():
    head = "\n".join(x for x in OUT.read_text().splitlines() if x.startswith("#"))
    assert "not measured on transcripts" in head and "nonsense-mediated decay" in head
    assert "urn:mavedb:00001277-a-1" in head and "checked in E18" in head
    assert "3-10 and 11-50 nt" in head and "DDX3X 11-25 nt" in head
    assert "drawn again under the current labels in E20b (2026-10-26)" in head


def test_the_stage_runs_last():
    spec = importlib.util.spec_from_file_location("reproduce_evidence", REPO / "scripts/reproduce_evidence.py")
    E = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(E)
    mods = [m for m, _, _ in E.STAGES]
    assert mods.index("src.evid_suppl_c12") > mods.index("src.evid_suppl_c3_brca2")
