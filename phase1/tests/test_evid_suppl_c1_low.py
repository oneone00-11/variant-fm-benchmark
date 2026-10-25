"""E19 C1, low circularity: assertions last evaluated before the gene's SGE data were public.

Each gene's first public date is the earliest of its dated sources, a preprint's being
the earlier of Crossref's posted and created; the preprint searches are recorded for
every gene and their matches follow the rule; the archive is compared with the VCF as
the docstring says; every variant's reason is the first that applies; and the subset's
composition is its variants'.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from datetime import date
from pathlib import Path

import pandas as pd
import pytest

PHASE1 = Path(__file__).resolve().parents[1]
REPO = PHASE1.parent
OUT = PHASE1 / "reports/evidence/suppl_analyses"
DATES = PHASE1 / "data/evidence/e19/first_public_dates"

if str(PHASE1) not in sys.path:
    sys.path.insert(0, str(PHASE1))

from src import evid_common as K  # noqa: E402
from src import evid_suppl_c1_low as L  # noqa: E402

need = pytest.mark.skipif(not (OUT / "c1_low_circularity_variants.csv").exists(),
                          reason="run the pipeline first")
SEVEN = ["BAP1", "BARD1", "BRCA1", "BRCA2", "PALB2", "RAD51C", "VHL"]


def test_every_saved_record_is_the_one_its_manifest_names():
    for name in ("manifest.json", "preprint_search.json"):
        man = json.loads((DATES / name).read_text())
        for f in man["files"]:
            if f.get("file"):
                assert hashlib.sha256((REPO / f["file"]).read_bytes()).hexdigest() == f["sha256"]


@need
def test_each_first_public_date_is_the_earliest_of_its_sources():
    d = K.read_back(OUT / "c1_first_public_dates.csv", keep_default_na=False, na_values=[""])
    assert set(d.gene) == set(SEVEN)
    for gene, g in d.groupby("gene"):
        first = g[g.source == "first_public_date"]
        assert len(first) == 1
        dated = g[(g.source != "first_public_date") & g.date.notna()]
        assert first.date.iloc[0] == dated.date.min()
    b = d[(d.gene == "BRCA2") & (d.source == "first_public_date")].iloc[0]
    assert (b.identifier, b.date) == ("10.1101/2023.12.14.571597", "2023-12-15")


@need
def test_a_preprint_date_is_the_earlier_of_posted_and_created():
    d = K.read_back(OUT / "c1_first_public_dates.csv", keep_default_na=False, na_values=[""])
    for r in d[d.source == "preprint"].itertuples():
        m = json.loads((PHASE1 / L._repo_path(r.record)).read_text())["message"]
        got = [date(*m[k]["date-parts"][0][:3]) for k in ("posted", "created") if m.get(k)]
        assert r.date == min(got).isoformat()


@need
def test_every_gene_has_its_search_record_and_the_matches_follow_the_rule():
    s = K.read_back(OUT / "c1_preprint_search.csv", keep_default_na=False)
    summary = s[s.preprint == ""].set_index("gene")
    assert set(summary.index) == set(SEVEN)
    hits = s[s.preprint != ""]
    rule = (hits.title_names_gene.astype(str) == "True") & \
           (hits.names_saturation_genome_editing.astype(str) == "True") & \
           (hits.authors_met.astype(str) == "True")
    assert ((hits.is_the_studys.astype(str) == "True") == rule).all()
    for gene, r in summary.iterrows():
        found = sorted(set(hits[(hits.gene == gene) & rule].doi))
        assert r.title == ("found: " + "; ".join(found) if found else "not found")
    assert summary.loc["BRCA2", "title"].startswith("found: 10.1101/2023.12.14.571597")
    for g in ("BAP1", "RAD51C", "PALB2"):
        assert summary.loc[g, "title"] == "not found"


@need
def test_every_variants_reason_is_the_first_that_applies():
    v = K.read_back(OUT / "c1_low_circularity_variants.csv", keep_default_na=False)
    assert ((v.why_not == "") == (v.in_subset.astype(str) == "True")).all()
    assert set(v.why_not) <= {"", "not P/LP or B/LB", "fewer than one star", "archive and VCF differ",
                              "no LastEvaluated date", "evaluated on or after the first public date"}
    agree = v.agree.astype(str) == "True"
    assert agree.mean() >= L.MIN_AGREEMENT
    sub = v[v.in_subset.astype(str) == "True"]
    assert (sub.clnsig_class.isin(["P/LP", "B/LB"])).all() and (sub.agree.astype(str) == "True").all()
    assert (pd.to_datetime(sub.last_evaluated) < pd.to_datetime(sub.first_public_date)).all()


@need
def test_the_composition_is_the_subsets_variants():
    v = K.read_back(OUT / "c1_low_circularity_variants.csv", keep_default_na=False)
    c = K.read_back(OUT / "c1_low_circularity_composition.csv")
    sub = v[v.in_subset.astype(str) == "True"]
    for r in c[c.gene != "Seven genes"].itertuples():
        g = sub[sub.gene == r.gene]
        g = K.stratum_frame(g, r.stratum)
        assert (r.subset_n_plp, r.subset_n_blb) == (int((g.clnsig_class == "P/LP").sum()),
                                                    int((g.clnsig_class == "B/LB").sum()))


def test_the_stage_runs_after_c1():
    spec = importlib.util.spec_from_file_location(
        "reproduce_evidence", REPO / "scripts" / "reproduce_evidence.py")
    E = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(E)
    mods = [m for m, _, _ in E.STAGES]
    assert mods.index("src.evid_suppl_c1_low") == mods.index("src.evid_suppl_c1") + 1
