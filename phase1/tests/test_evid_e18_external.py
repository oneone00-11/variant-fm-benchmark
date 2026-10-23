"""E18.3: the external test reads the analysis set's variables.

DDX3X and TP53 are tested on E18.2's columns (docs/external-column-identity.md): E7
reads each of them value for value; every column it tests is the analysis set's
variable, and its column check says so for twelve columns of each gene; TP53's
declared definition mismatch is withdrawn, its client-0.6.1 AlphaGenome kept beside the
official definition and not tested, its GPN-MSA turned to the analysis set's direction;
NT, which has no column of the same variable on either gene, is named in the column
check with its reason and tested on neither; and LICENSE-DATA covers every score file
stage E18 committed.
"""
from __future__ import annotations

import fnmatch
import importlib.util
import json
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

PHASE1 = Path(__file__).resolve().parents[1]
REPO = PHASE1.parent
EXT = PHASE1 / "data/evidence/external"
REPORTS = PHASE1 / "reports/evidence"

if str(PHASE1) not in sys.path:
    sys.path.insert(0, str(PHASE1))

from src import evid_common as K  # noqa: E402
from src import evid_external as X  # noqa: E402


def _build_columns():
    stage = REPO / "scripts/stage_E18"
    if str(stage) not in sys.path:
        sys.path.insert(0, str(stage))
    spec = importlib.util.spec_from_file_location("e18_build_columns", stage / "build_columns.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture
def in_phase1(monkeypatch):
    monkeypatch.chdir(PHASE1)


@pytest.mark.parametrize("gene", ["DDX3X", "TP53"])
def test_the_column_check_finds_twelve_comparable_columns_and_names_nt(gene):
    b = pd.read_csv(REPORTS / f"external_{gene.lower()}_column_basis.csv")
    assert list(b.tool) == K.PANEL + K.OPTIONAL
    tested = b[b.tool != "nt"]
    assert tested.comparable.all() and len(tested) == 12
    nt = b[b.tool == "nt"].iloc[0]
    assert not nt.comparable and nt.reason.startswith("no column of the same variable")
    ext = pd.read_csv(REPORTS / f"external_{gene.lower()}.csv")
    assert set(ext.tool) == set(tested.tool)


def test_no_definition_mismatch_is_declared():
    assert X.DECLARED_MISMATCH == {}


def _alphagenome(gene: str) -> pd.Series:
    ag = pd.read_csv(PHASE1 / "data/evidence/alphagenome_official/external_scores.tsv", sep="\t",
                     usecols=["gene", "variant_id", "alphagenome_splicing"],
                     float_precision="round_trip")
    return ag[ag.gene == gene].set_index("variant_id")["alphagenome_splicing"]


def _same(a, b) -> bool:
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    return bool(((a == b) | (np.isnan(a) & np.isnan(b))).all())


@pytest.mark.parametrize("gene", ["DDX3X", "TP53"])
def test_e7_reads_e18s_columns_value_for_value(gene, in_phase1):
    df, _ = X._load_external(gene)
    df = df.set_index("variant_id")
    files = X.DDX3X_SCORE_FILES if gene == "DDX3X" else X.TP53_SCORE_FILES
    for fname, cols in files.items():
        if not (EXT / fname).exists():
            assert cols == ["nt"] and "nt" not in df.columns
            continue
        d = pd.read_parquet(EXT / fname).set_index("variant_id")
        for c in cols:
            assert _same(df.loc[d.index, c], d[c]), (fname, c)
    assert _same(df["alphagenome"], _alphagenome(gene).loc[df.index])
    if gene == "TP53":
        t = pd.read_parquet(X.TP53_TABLE).set_index("variant_id").loc[df.index]
        assert _same(df["gpn_msa"], -t["gpn_msa"])
        assert _same(df["alphagenome_v061"], t["alphagenome"])
        for c in ("cadd", "phylop", "phastcons"):
            assert _same(df[c], t[c]), c


def test_licence_data_names_every_score_file_stage_e18_committed():
    B = _build_columns()
    text = (REPO / "LICENSE-DATA").read_text()
    flat = " ".join(text.split())
    sec1 = flat.split(" 1. The frozen analysis set")[1].split(" 2. AlphaGenome outputs")[0]
    sec2 = flat.split(" 2. AlphaGenome outputs")[1].split(" 3. What a reuser")[0]
    patterns = re.findall(r"`([^`]+\.(?:parquet|tsv))`", sec1)
    ag = {"newgenes_avi.parquet", "tp53_avi_rescored.parquet", "e18_anchors_avi.parquet",
          "newgenes_scored.parquet", "e18_anchors.parquet"}
    files = {B.file_name(r, g) for r, spec in B.RUNS.items() for g in spec["genes"]
             if (EXT / B.file_name(r, g)).exists()}
    files |= {p.name for p in EXT.glob("e18_anchors_*.parquet")} | {"newgenes_scored.parquet"}
    for f in sorted(files):
        if f in ag:
            assert f"`{f}`" in sec2 or f in sec2, f
        else:
            assert any(fnmatch.fnmatch(f, p.replace("<gene>", "*").replace("<run>", "*"))
                       for p in patterns), f
    assert "nt_gpu_run/anchor_check.tsv" in sec1
    for f in ("external_scores.tsv", "drift_2b_start_scores.tsv", "drift_2b_end_scores.tsv"):
        assert f in sec2, f
