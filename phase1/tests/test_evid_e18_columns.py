"""E18.2: the external genes' columns are the analysis set's variables.

The external score files of E18.2 are written by scripts/stage_E18/build_columns.py
from the scorer runs' caches, which are not committed. From the committed files alone:
each run's anchors equal the analysis set's stored values bit for bit, in every column
of the run the analysis set has (E18.0 rule 4(a)); the anchors are the seeded draw
scripts/stage_E18/anchors.py makes; each gene's file holds exactly that gene's E18
variants, in e18_score_input.parquet's order, and its record names its sha256, its
counts and its scoring times, every one after the tag; the five further genes' table and
the evidence table rebuild from the committed files byte for byte (build_columns.py
--offline); every cell of the evidence table, 7 genes x 13 columns, has evidence (a) or
(b), or is NT's, in a commit of its own; and run_avi.py's guard stops the Atlas run on
an authentication error and on a batch whose every attempt fails, with the key cut out
of what it records.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import re
import sys
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

PHASE1 = Path(__file__).resolve().parents[1]
REPO = PHASE1.parent
EXT = PHASE1 / "data/evidence/external"
SET = PHASE1 / "data/evidence/analysis_set_v1.parquet"
STAGE = REPO / "scripts/stage_E18"
# the tagger time of the tag frozen-thresholds-v1, from its tag object
# (`git cat-file -p frozen-thresholds-v1`: 1792686865 +0800)
TAG_UTC = "2026-10-22T16:34:25+00:00"

if str(PHASE1) not in sys.path:
    sys.path.insert(0, str(PHASE1))


def _load(name: str):
    if str(STAGE) not in sys.path:
        sys.path.insert(0, str(STAGE))
    spec = importlib.util.spec_from_file_location(f"e18_{name}", STAGE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


B = _load("build_columns")


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _record(path: Path) -> dict:
    return json.loads(Path(str(path) + ".provenance.json").read_text())


def _present(run: str) -> bool:
    return (EXT / f"e18_anchors_{run}.parquet").exists()


def _files(run: str) -> list[tuple[str, Path]]:
    """(gene or "anchors", path) for every file build_columns writes for the run."""
    spec = B.RUNS[run]
    out = [("anchors", EXT / f"e18_anchors_{run}.parquet")]
    out += [(g, EXT / B.file_name(run, g)) for g in spec["genes"] if not (run == "avi" and g in B.NEW)]
    if run == "avi":
        out.append(("NEW", EXT / "newgenes_avi.parquet"))
    return out


RUNS = [r for r in B.RUNS if _present(r) or not B.RUNS[r].get("optional")]


@pytest.mark.parametrize("run", RUNS)
def test_every_file_of_the_run_is_there_with_its_records(run):
    for _, path in _files(run):
        for p in (path, Path(str(path) + ".provenance.json"), Path(str(path) + ".environment.json")):
            assert p.is_file(), p


@pytest.mark.parametrize("run", RUNS)
def test_the_runs_anchors_equal_the_analysis_set_bit_for_bit(run):
    got = pd.read_parquet(EXT / f"e18_anchors_{run}.parquet")
    anchors = pd.read_parquet(EXT / "e18_anchors.parquet")
    assert list(got["variant_id"]) == list(anchors["variant_id"])
    s = pd.read_parquet(SET).set_index("variant_id").loc[got["variant_id"]]
    cols = [c for c in B.RUNS[run]["columns"] if c in s.columns]
    assert set(B.RUNS[run]["columns"]) & set(B.COLUMNS) <= set(cols)
    rec = _record(EXT / f"e18_anchors_{run}.parquet")
    for c in cols:
        g, w = got[c].to_numpy(dtype=float), s[c].to_numpy(dtype=float)
        assert ((g == w) | (np.isnan(g) & np.isnan(w))).all(), c
        assert rec["same_variable"]["n_identical"][c] == len(got) == 39


def test_the_anchors_are_the_seeded_draw():
    A = _load("anchors")
    drawn = A.draw(pd.read_parquet(SET))
    anchors = pd.read_parquet(EXT / "e18_anchors.parquet")
    assert list(drawn["variant_id"]) == list(anchors["variant_id"])
    assert len(anchors) == 39 and anchors.groupby(["gene", "stratum"]).size().max() == 3


@pytest.mark.parametrize("run", RUNS)
def test_each_file_holds_exactly_its_genes_variants_in_order(run):
    inp = pd.read_parquet(EXT / "e18_score_input.parquet")
    for g, path in _files(run):
        got = pd.read_parquet(path)
        if g == "anchors":
            want = inp[inp["set"] == "anchor"]
        elif g == "NEW":
            want = inp[(inp["set"] == "external") & inp["gene"].isin(B.NEW)]
        else:
            want = inp[(inp["set"] == "external") & (inp["gene"] == g)]
        assert list(got["variant_id"]) == list(want["variant_id"]), path.name
        assert list(got[["gene", "chrom", "pos", "ref", "alt"]].itertuples(index=False)) == \
            list(want[["gene", "chrom", "pos", "ref", "alt"]].itertuples(index=False))


@pytest.mark.parametrize("run", RUNS)
def test_each_record_matches_its_file_and_every_score_postdates_the_tag(run):
    tag = datetime.fromisoformat(TAG_UTC)
    for _, path in _files(run):
        df, rec = pd.read_parquet(path), _record(path)
        assert rec["sha256"] == _sha(path) and rec["n"] == len(df) and rec["run"] == run
        assert rec["n_scored"] == {c: int(df[c].notna().sum()) for c in B.RUNS[run]["columns"]}
        first, last = rec["scored_utc_range"]
        if "scored_utc" in df.columns:
            assert [first, last] == [min(df["scored_utc"]), max(df["scored_utc"])]
            assert all(datetime.fromisoformat(t) > tag for t in df["scored_utc"])
        assert tag < datetime.fromisoformat(first) <= datetime.fromisoformat(last)
        env = json.loads(Path(str(path) + ".environment.json").read_text())
        assert env["file"] == rec["file"] and env["pip_freeze"]


def test_the_assembled_tables_rebuild_byte_for_byte():
    """build_columns.py --offline: from the committed files alone."""
    df, manifest = B.assemble()
    for path, data in B.write_assembled(df, manifest, {}).items():
        assert path.read_bytes() == data, path


def test_the_five_genes_table_is_their_variants_with_every_column():
    df = pd.read_parquet(B.SCORED)
    man = json.loads(B.MANIFEST.read_text())
    assert man["sha256"] == _sha(B.SCORED) and man["n"] == len(df) == 4522
    assert man["n_by_gene"] == {"CTCF": 1370, "SFPQ": 1019, "RAD51D": 1296, "XRCC2": 390, "SBDS": 447}
    present = [c for c in B.COLUMNS if c not in man["not_scored"]]
    assert set(man["not_scored"]) <= {"nt"} and set(present) <= set(df.columns)
    assert datetime.fromisoformat(man["scored_utc_range"][0]) > datetime.fromisoformat(TAG_UTC)


def test_every_cell_of_the_evidence_table_has_evidence():
    rows = [r for r in B.DOC.read_text().splitlines() if r.startswith("| ") and not r.startswith("| gene")]
    cells = [[c.strip() for c in re.split(r"(?<!\\)\|", r.strip()[1:-1])] for r in rows]
    assert [(c[0], c[1]) for c in cells] == [(g, col) for g in B.GENES for col in B.COLUMNS]
    for c in cells:
        ev = c[6]
        assert ev.startswith("(a)") or ev.startswith("(b)") or (
            c[1] == "nt" and ev.startswith("not part of this commit")), c
        if "re-scored" in ev:
            assert "anchors 39/39 identical" in ev, c


class _Fake:
    def __init__(self, errors):
        self.errors, self.calls = list(errors), 0

    def scorer_metadata(self):
        return {"AVI_SCORE": None}

    def query_variants(self, *a, **k):
        self.calls += 1
        if self.errors:
            e = self.errors.pop(0)
            if e is not None:
                raise e
        return "scores"


def _run_like_evid_score_avi(guard, attempts):
    """The shape of evid_score_avi.run()'s retry loop: `except Exception`, then retry."""
    for attempt in range(attempts):
        try:
            return guard.query_variants([])
        except Exception:                                       # noqa: BLE001
            if attempt == attempts - 1:
                return "the batch failed and run() went on"


@pytest.fixture
def run_avi(monkeypatch):
    R = _load("run_avi")
    monkeypatch.setattr(R, "BACKOFF", 0.0)
    return R


def test_the_guard_stops_on_an_authentication_error_with_the_key_cut_out(run_avi):
    A = SimpleNamespace(read_key=lambda: "s3cr3t", MAX_RETRIES=6)
    g = run_avi.Guarded(_Fake([RuntimeError("UNAUTHENTICATED: bad key s3cr3t")]), A)
    with pytest.raises(run_avi.Stop) as e:
        _run_like_evid_score_avi(g, 6)
    assert "authentication" in str(e.value) and "s3cr3t" not in str(e.value) and "<key>" in str(e.value)
    assert g._c.calls == 1


def test_the_guard_stops_when_a_batch_fails_every_attempt(run_avi):
    A = SimpleNamespace(read_key=lambda: None, MAX_RETRIES=6)
    g = run_avi.Guarded(_Fake([RuntimeError("UNAVAILABLE")] * 6), A)
    with pytest.raises(run_avi.Stop):
        _run_like_evid_score_avi(g, 6)
    assert g._c.calls == 6 and len(g.errors) == 6


def test_the_guard_lets_a_batch_that_recovers_go_on(run_avi):
    A = SimpleNamespace(read_key=lambda: None, MAX_RETRIES=6)
    g = run_avi.Guarded(_Fake([RuntimeError("UNAVAILABLE")] * 5 + [None]), A)
    assert _run_like_evid_score_avi(g, 6) == "scores" and g.in_a_row == 0
