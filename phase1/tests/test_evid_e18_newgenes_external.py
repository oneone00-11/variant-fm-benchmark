"""E18.6: the five genes of E16 and E16b enter the external test.

Their labels are the ones their <gene>_variants.parquet files carry, value for value,
four from the one investigator-provided calibration their records carry and SBDS's
from its class column, each named for its source; every column E7 tests on them is the
analysis set's variable, NT named with its reason; only the three in-scope strata are
evaluated, at least ten on each side; S12 carries their rows; the normal-only strata of
all seven external genes count the normal variants at or above each frozen PP3 threshold
and Walker's cut point; the external summary's k/n (m) is what E7's own tables give;
none of the five assays is recorded in a predictor's training or model selection; and
the entry point runs the stage after DDX3X's and before anything that reads it.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

PHASE1 = Path(__file__).resolve().parents[1]
REPO = PHASE1.parent
EXT = PHASE1 / "data/evidence/external"
REPORTS = PHASE1 / "reports/evidence"
SUPP = REPORTS / "supplement"

if str(PHASE1) not in sys.path:
    sys.path.insert(0, str(PHASE1))

from src import evid_common as K  # noqa: E402
from src import evid_external as X  # noqa: E402
from src import evid_supp_tables as S  # noqa: E402

IN_SCOPE = ["s3_10", "s11_50", "s3_50"]
RANGE_GENES = ["CTCF", "SFPQ", "RAD51D", "XRCC2"]
need = pytest.mark.skipif(not (REPORTS / "external_summary.csv").exists(),
                          reason="run the pipeline first")


@pytest.fixture
def in_phase1(monkeypatch):
    monkeypatch.chdir(PHASE1)


def _sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


@pytest.mark.parametrize("gene", X.NEW_GENES)
def test_the_labels_are_the_ones_the_variants_file_carries(gene, in_phase1):
    df, labels = X._load_external(gene)
    assert labels == [X.NEW_GENE_LABEL[gene]]
    v = pd.read_parquet(EXT / f"{gene.lower()}_variants.parquet")
    v = v[(v["offset"].abs() >= 3) & (v["offset"].abs() <= 50)].reset_index(drop=True)
    assert list(df.variant_id) == list(v.variant_id)
    assert np.array_equal(df[labels[0]].to_numpy(float), v["label"].to_numpy(float),
                          equal_nan=True)
    man = json.loads((EXT / "newgenes_scored.manifest.json").read_text())
    path = f"phase1/data/evidence/external/{gene.lower()}_variants.parquet"
    assert man["inputs"][path] == _sha256(REPO / path)
    assert set(df.urn) == {X.PRIMARY_LABEL_SOURCE[gene][1]}


def test_each_label_is_named_for_the_source_its_record_names():
    """The four range genes: the one calibration their records carry is investigator-
    provided, titled "Investigator-provided functional classes", and E16 selected it.
    SBDS: its record carries none; its labels are the score table's class column."""
    for gene in RANGE_GENES:
        meta = json.loads((EXT / f"{gene.lower()}_raw/scoreset_meta.json").read_text())
        cals = meta["scoreCalibrations"]
        assert [(c["title"], c["investigatorProvided"]) for c in cals] == [
            ("Investigator-provided functional classes", True)]
        sel = json.loads((EXT / f"{gene.lower()}_raw/score_ranges.json").read_text())
        assert sel["selection"] == "investigatorProvided is true"
        assert [c["title"] for c in sel["investigator_provided"]] == [
            "Investigator-provided functional classes"]
        assert X.PRIMARY_LABEL_SOURCE[gene] == ("ranges", meta["urn"])
        assert X.primary_label_name(gene) == f"Investigator-provided range ({meta['urn']})"
        assert X.NEW_GENE_LABEL[gene] == "y_investigator_range"
    meta = json.loads((EXT / "sbds_raw/scoreset_meta.json").read_text())
    assert not meta.get("scoreCalibrations")
    assert X.PRIMARY_LABEL_SOURCE["SBDS"] == ("classification", meta["urn"])
    assert X.primary_label_name("SBDS") == "Deposit classification (urn:mavedb:00001284-a-1)"
    assert X.NEW_GENE_LABEL["SBDS"] == "y_deposit"


@need
@pytest.mark.parametrize("gene", X.NEW_GENES)
def test_the_column_check_finds_twelve_comparable_columns_and_names_nt(gene):
    b = pd.read_csv(REPORTS / f"external_{gene.lower()}_column_basis.csv")
    assert list(b.tool) == K.PANEL + K.OPTIONAL
    tested = b[b.tool != "nt"]
    assert tested.comparable.all() and len(tested) == 12
    nt = b[b.tool == "nt"].iloc[0]
    assert not nt.comparable and nt.reason.startswith("no column of the same variable")


@need
@pytest.mark.parametrize("gene", X.NEW_GENES)
def test_e7_evaluates_the_three_in_scope_strata_at_ten_a_side(gene, in_phase1):
    e = pd.read_csv(REPORTS / f"external_{gene.lower()}.csv")
    assert set(e.stratum) == set(IN_SCOPE) and set(e.label_definition) == {X.NEW_GENE_LABEL[gene]}
    df, (lab,) = X._load_external(gene)
    for r in e.itertuples():
        sub = df if r.stratum == "s3_50" else df[df.stratum == r.stratum]
        s = sub.dropna(subset=[lab, r.tool])
        assert (r.n_pos, r.n_neg) == (int((s[lab] == 1).sum()), int((s[lab] == 0).sum()))
        ok = r.n_pos >= K.MIN_POS and r.n_neg >= K.MIN_NEG
        assert (r.status == "ok") == ok


@need
def test_s12_carries_the_new_genes_under_their_label_names():
    s12 = pd.read_csv(SUPP / "tableS12_external_genes.csv", dtype=str, keep_default_na=False)
    assert list(dict.fromkeys(s12.Gene)) == S.S12_GENES
    for gene in X.NEW_GENES:
        rows = s12[s12.Gene == gene]
        assert set(rows["Label definition"]) == {X.primary_label_name(gene)}
        assert set(rows["Stratum (intronic offset, bp)"]) == {S._ext_band(gene, st)
                                                              for st in IN_SCOPE}
        assert S.tool_labels()["nt"] not in set(rows.Tool)


@need
def test_the_normal_only_strata_are_counted_from_the_data(in_phase1):
    out = K.read_back(PHASE1 / "reports/evidence/external_normal_only.csv")
    assert X.NORMAL_ONLY_GENES == ["DDX3X", "TP53"] + X.NEW_GENES
    ev = K.read_back(REPORTS / "evidence_thresholds.csv")
    pp3 = 0.2
    want, normal_only = set(), set()
    for gene in X.NORMAL_ONLY_GENES:
        df, _ = X._load_external(gene)
        lab = X.PRIMARY_LABEL[gene]
        basis = pd.read_csv(REPORTS / f"external_{gene.lower()}_column_basis.csv")
        tools = [t for t in K.PANEL + K.OPTIONAL
                 if t in set(basis.loc[basis.comparable.astype(bool), "tool"])]
        for st in IN_SCOPE:
            sub = (df[df.stratum != "pm12"] if st == "s3_50" else df[df.stratum == st])
            for tool in tools:
                s = sub.dropna(subset=[lab, tool])
                y, v = s[lab].to_numpy(float), s[tool].to_numpy(float)
                if not ((y == 1).sum() < K.MIN_POS and (y == 0).sum() >= K.MIN_NEG):
                    continue
                normal_only.add((gene, st))
                rows = out[(out.gene == gene) & (out.stratum == st) & (out.tool == tool)]
                kinds = []
                if tool in ("spliceai", "spliceai_walker"):
                    kinds.append(("walker_pp3", pp3))
                e = ev[(ev.tool == tool) & (ev.stratum == st) & (ev.status == "ok")]
                for tier in ("supporting", "moderate", "strong"):
                    t = e[e.tier == tier]
                    if len(t) and bool(t.pp3_threshold_reachable.iloc[0]):
                        kinds.append((tier, float(t.pp3_threshold_insample.iloc[0])))
                # a column with no frozen PP3 threshold in the stratum, and not SpliceAI,
                # has nothing to count
                if kinds:
                    want.add((gene, st, tool))
                assert list(rows.threshold_kind) == [k for k, _ in kinds]
                for r, (_, thr) in zip(rows.itertuples(), kinds):
                    assert r.threshold == thr and r.n_neg == int((y == 0).sum())
                    assert r.n_neg_at_or_above == int((v[y == 0] >= thr).sum())
    assert set(zip(out.gene, out.stratum, out.tool)) == want
    # SBDS 11-50 nt holds no damaging variant, XRCC2 fewer than ten in either band, and
    # TP53's 11-12 nt none under its primary label
    assert {("SBDS", "s11_50"), ("XRCC2", "s3_10"), ("XRCC2", "s11_50"),
            ("TP53", "s11_50")} == normal_only == {(g, st) for g, st, _ in want}


RANK = {"below supporting": 0, "Supporting": 1, "Moderate": 2, "Strong": 3, "Very strong": 4}


@need
def test_the_external_summary_is_what_e7s_tables_give():
    """Recomputed from E7's own tables, not from S12: k/n (m) per gene set, stratum,
    column, side and tier, n counting every gene counted as in Table 3; DDX3X x the
    combined Atlas score not counted."""
    summ = pd.read_csv(REPORTS / "external_summary.csv", dtype=str, keep_default_na=False)
    assert list(summ.columns) == ["Gene set", "Stratum (intronic offset, bp)", "Tool", "Side",
                                  "Tier", "External genes reaching the tier, k/n (m)", "Note"]
    sets = {"Seven external genes": S.S12_GENES, "DDX3X and TP53": ["DDX3X", "TP53"]}
    assert list(dict.fromkeys(summ["Gene set"])) == list(sets)
    assert (summ.groupby("Gene set").size() == 3 * 13 * 5).all()
    ev = K.read_back(REPORTS / "evidence_thresholds.csv")
    ext = {g: K.read_back(REPORTS / f"external_{g.lower()}.csv") for g in S.S12_GENES}
    labels = {v: k for k, v in S.tool_labels().items()}
    strata = {v: k for k, v in S.STRATUM_LABEL.items()}
    for r in summ.itertuples(index=False):
        genes = sets[r[0]]
        st, tool = strata[r[1]], labels[r[2]]
        side, tier = r[3], r[4].lower()
        stat, note = r[5], r[6]
        if tool == "nt":
            assert stat == "No column of the same variable"
            continue
        e = ev[(ev.tool == tool) & (ev.stratum == st) & (ev.tier == tier) & (ev.status == "ok")]
        reach = "pp3_threshold_reachable" if side == "PP3" else "bp4_threshold_reachable"
        if not len(e) or not bool(e[reach].iloc[0]):
            assert stat == "–"
            continue
        col = (f"pp3_insample_{tier}_tier_here" if side == "PP3"
               else f"bp4_insample_{tier}_tier_here")
        k = n = m = 0
        for g in genes:
            if (g, tool) == ("DDX3X", "avi"):
                assert "DDX3X not counted" in note
                continue
            n += 1
            x = ext[g]
            x = x[(x.label_definition == X.PRIMARY_LABEL[g]) & (x.stratum == st)
                  & (x.tool == tool)]
            if not len(x) or x.status.iloc[0] != "ok" or x[col].iloc[0] not in RANK:
                m += 1
                continue
            k += RANK[x[col].iloc[0]] >= RANK[tier.capitalize()]
        assert stat == f"{k}/{n}" + (f" ({m})" if m else ""), (r, k, n, m)


@need
def test_no_new_gene_is_in_a_predictors_training_or_model_selection():
    hits = S.training_overlap(S.S12_GENES)
    assert hits == {("DDX3X", "avi")}
    assert S.not_external() == {("DDX3X", "avi")}


def test_the_entry_point_runs_the_new_genes_after_ddx3x_and_before_their_readers():
    spec = importlib.util.spec_from_file_location(
        "reproduce_evidence", REPO / "scripts" / "reproduce_evidence.py")
    E = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(E)
    rows = [(m, a) for m, _, a in E.STAGES]
    i = rows.index(("src.evid_external", ["--apply-new"]))
    assert rows[i - 1] == ("src.evid_external", ["--apply", "DDX3X"])
    mods = [m for m, _ in rows]
    assert i < mods.index("src.evid_supp_tables")
