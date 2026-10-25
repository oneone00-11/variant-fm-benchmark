"""E19 C1: the predictors against ClinVar's germline classifications, as a concordance check.

The star conversion is the saved page's table for the germline classification of
aggregate records, and grades every review status of the set; every count is the
analysis set's; a pooled AUROC is missing exactly where a gene's AUROC is 0 or 1; no
likelihood ratio or tier appears; and the cross-tabulation counts are the set's.
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
PAGE = PHASE1 / "data/evidence/e19/clinvar_review_status"

if str(PHASE1) not in sys.path:
    sys.path.insert(0, str(PHASE1))

from src import evid_common as K  # noqa: E402
from src import evid_suppl_c1 as C1  # noqa: E402

need = pytest.mark.skipif(not (OUT / "c1_clinvar_reference_metrics.csv").exists(),
                          reason="run the pipeline first")


def _set() -> pd.DataFrame:
    d = pd.read_parquet(PHASE1 / "data/evidence/analysis_set_v1.parquet")
    d["clinvar_stars"] = C1.stars_of(d, C1.star_table(PAGE / "review_status.html"))
    d["y_clinvar"] = d["clnsig_class"].map({"P/LP": 1.0, "B/LB": 0.0})
    return d


def test_the_saved_page_is_the_documentation_as_its_manifest_records():
    import hashlib
    import json
    man = json.loads((PAGE / "manifest.json").read_text())
    data = (PAGE / "review_status.html").read_bytes()
    assert hashlib.sha256(data).hexdigest() == man["sha256"] and len(data) == man["bytes"]
    assert man["source_url"] == "https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/"
    assert "saved by the author with a browser" in man["how"]
    assert b"recaptcha" not in data.lower()


def test_the_star_table_is_the_germline_aggregate_one_and_grades_every_status():
    t = C1.star_table(PAGE / "review_status.html")
    assert (t.table == " / ".join(C1.TABLE_HEADINGS)).all()
    lut = dict(zip(t.review_status, t.stars))
    assert lut["criteria provided, multiple submitters, no conflicts"] == 2
    assert lut["criteria provided, conflicting classifications"] == 1
    assert lut["reviewed by expert panel"] == 3 and lut["practice guideline"] == 4
    assert t.stars.isin([0, 1, 2, 3, 4]).all() and not t.review_status.duplicated().any()
    d = _set()
    graded = d.clnrevstat.fillna("").str.len() > 0
    assert d.loc[graded, "clinvar_stars"].notna().all()
    assert d.loc[~graded, "clinvar_stars"].isna().all()


@need
def test_the_conversion_file_counts_the_sets_statuses():
    c = K.read_back(OUT / "c1_clinvar_star_conversion.csv")
    d = pd.read_parquet(PHASE1 / "data/evidence/analysis_set_v1.parquet")
    assert int(c.n_variants_in_set.sum()) == int((d.clnrevstat.fillna("").str.len() > 0).sum())


@need
def test_the_reference_counts_are_the_sets_and_no_ratio_appears():
    m = K.read_back(OUT / "c1_clinvar_reference_metrics.csv")
    assert not [c for c in m.columns
                if c.endswith("lr") or "_lr" in c or "tier" in c or "posterior" in c]
    d = _set()
    for level, k, _ in C1.STAR_LEVELS:
        ref = d[d.y_clinvar.notna() & (d.clinvar_stars >= k)]
        for st in C1.BANDS:
            r = m[(m.star_level == level) & (m.stratum == st) & (m.measure == "auroc")
                  & (m.tool == "cadd")].iloc[0]
            f = K.stratum_frame(ref, st).dropna(subset=["cadd"])
            assert (r.n_plp, r.n_blb) == (int((f.y_clinvar == 1).sum()), int((f.y_clinvar == 0).sum()))
            ok = r.n_plp >= K.MIN_POS and r.n_blb >= K.MIN_NEG
            assert (r.status != "not evaluable") == ok


@need
def test_a_pooled_auroc_is_missing_exactly_where_a_gene_is_0_or_1():
    m = K.read_back(OUT / "c1_clinvar_reference_metrics.csv", keep_default_na=False,
                    na_values=[""])
    a = m[(m.measure == "auroc") & m.status.isin(["ok", "undefined"])]
    for r in a.itertuples():
        per = [float(x.split()[1]) for x in r.auroc_per_gene.split("; ")]
        degenerate = any(v in (0.0, 1.0) for v in per)
        assert (r.status == "undefined") == degenerate
        assert np.isfinite(r.auroc) == (not degenerate)
        if not degenerate:
            assert r.auroc_lo <= r.auroc <= r.auroc_hi


@need
def test_the_crosstab_counts_are_the_sets():
    x = K.read_back(OUT / "c1_clinvar_sge_crosstab.csv")
    d = _set()
    r = x[(x.star_level == "stars_ge_1") & (x.tool == "spliceai_walker") & (x.stratum == "s3_10")
          & (x.threshold_kind == "walker_pp3")]
    for row in r.itertuples():
        if row.clinvar_group == "unrecorded":
            g = d[d.clinvar_arm == "unrecorded"]
        else:
            g = d[(d.clnsig_class == row.clinvar_group) & (d.clinvar_stars >= 1)]
        g = K.stratum_frame(g, "s3_10").dropna(subset=["y_assay", "spliceai_walker"])
        above = g.spliceai_walker >= row.threshold
        assert row.n_pos_above == int(((g.y_assay == 1) & above).sum())
        assert row.n_neg_below == int(((g.y_assay == 0) & ~above).sum())


def test_the_stage_runs_after_c6_c7():
    spec = importlib.util.spec_from_file_location(
        "reproduce_evidence", REPO / "scripts" / "reproduce_evidence.py")
    E = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(E)
    mods = [m for m, _, _ in E.STAGES]
    assert mods.index("src.evid_suppl_c1") > mods.index("src.evid_suppl_c6_c7")
