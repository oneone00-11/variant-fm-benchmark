"""E17: AlphaGenome's splicing score by its official definition, and the records it rests on.

The definition is quoted from AlphaGenome's own notebook at a pinned commit, item by
item with the lines it rests on, and the notebook's words around it beside it; the
scoring environment is pinned with hashes, and the model a request with no version gets
is quoted from the client; the composite is the notebook's combination of its three
components; the scores cover the seven training genes only, so that no external gene is
scored before the thresholds are frozen (E17.0 rule 10); the pilot is the seeded draw it
says it is, and its variants scored again (ten of them twice after the pilot, once after
the full run) agree to the last bit; no record holds the API key; the scoring stops
when it should; and no AlphaGenome column enters the fusion model (the AlphaGenome
terms; E17.0 rule 5). Since the switch (E17.5): the analysis set reads the new scores
bit for bit, oriented as they come; the old composite stays beside them, outside every
model; the comparison of the two definitions holds only the study's own quantities;
and Supplementary Table S3 and LICENSE-DATA describe the new definition.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import re
import sys
import threading
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

PHASE1 = Path(__file__).resolve().parents[1]
REPO = PHASE1.parent
DIR = PHASE1 / "data/evidence/alphagenome_official"
SOURCE = DIR / "NOTEBOOK_SOURCE.json"
CONTEXT = DIR / "NOTEBOOK_CONTEXT.json"
ENV = DIR / "ENVIRONMENT.json"
REQ = REPO / "scripts/stage_E17/requirements.txt"
SET = PHASE1 / "data/evidence/analysis_set_v1.parquet"
COMPONENTS = ["SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS"]
COMPOSITE = "alphagenome_splicing"
FULL, DRIFT = "all_scores.tsv", "drift_scores.tsv"
SCORE_FILES = {"examples_scores.tsv", "pilot_scores.tsv", "repeat_scores.tsv", FULL, DRIFT}

if str(PHASE1) not in sys.path:
    sys.path.insert(0, str(PHASE1))


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _scores(name: str, as_text: bool = False) -> pd.DataFrame:
    """The scores as written: as text, or as numbers read back exactly (pandas' default
    float parser is not round-trip exact)."""
    kw = ({"dtype": str, "keep_default_na": False} if as_text
          else {"float_precision": "round_trip"})
    return pd.read_csv(DIR / name, sep="\t", **kw)


def _prov(name: str) -> dict:
    return json.loads((DIR / f"{name}.provenance.json").read_text())


def _present() -> list[str]:
    return sorted(p.name for p in DIR.glob("*_scores.tsv") if p.name in SCORE_FILES)


def _need(name: str) -> None:
    if not (DIR / name).exists():
        pytest.skip(f"{name} is written by a later step")


def _scorer():
    spec = importlib.util.spec_from_file_location(
        "score_official", REPO / "scripts/stage_E17/score_official.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# --------------------------------------------------------------------------
# the definition and the environment
# --------------------------------------------------------------------------
def test_every_item_of_the_definition_rests_on_a_quoted_line():
    src = json.loads(SOURCE.read_text())
    text = {c["cell"]: c["source"] for c in src["quoted_cells"]}
    for d in src["definition"]:
        assert d["quoted"], d["item"]
        for q in d["quoted"]:
            assert q["holds"] in text[q["cell"]].split("\n")[q["line"] - 1], (d["item"], q)
    assert {d["item"] for d in src["definition"]} >= {
        "client call", "scorers", "input window", "organism", "tissue or ontology",
        "score field", "aggregation", "combination", "missing values", "strand"}
    assert src["client"]["SEQUENCE_LENGTH_1MB"] == 2 ** 20
    assert [s["is_signed"] for s in src["client"]["scorers"]] == [False] * 3


def test_the_notebooks_words_are_quoted_from_the_same_file():
    """Gate 1 review R6b: markdown cells 1, 10 and 16 of the very notebook the code is
    quoted from (same commit, path and sha256), verbatim."""
    src, ctx = json.loads(SOURCE.read_text()), json.loads(CONTEXT.read_text())
    for k in ("repository", "commit", "path", "raw_url", "sha256", "bytes"):
        assert ctx[k] == src[k], k
    assert ctx["code_record"]["sha256"] == _sha(SOURCE)
    cells = {c["cell"]: c for c in ctx["quoted_cells"]}
    assert sorted(cells) == [1, 10, 16]
    for c in cells.values():
        assert c["type"] == "markdown"
        assert c["sha256"] == hashlib.sha256(c["source"].encode()).hexdigest()
    assert "the recommended method" in cells[1]["source"]
    assert "We use a 1MB sequence length" in cells[10]["source"]
    assert "We are working on a recommended variant score cutoff" in cells[16]["source"]


def test_the_scoring_environment_is_pinned_with_hashes():
    lines = [l for l in REQ.read_text().splitlines() if l.strip() and not l.startswith("#")]
    for line in lines:
        assert re.fullmatch(r"[A-Za-z0-9_.\-]+==\S+ --hash=sha256:[0-9a-f]{64}", line), line
    env = json.loads(ENV.read_text())
    pins = {line.split()[0] for line in lines}
    assert pins == set(env["pip_freeze"])
    assert f"alphagenome=={env['client']['version']}" in pins
    assert env["client"]["version"] == json.loads(SOURCE.read_text())["client"]["version"]
    assert env["requirements_sha256"] == _sha(REQ)


def test_the_default_model_is_quoted_from_the_client():
    """Gate 1 review R6a: the request names no model version, as in the notebook, and
    the records say which one that means, in the client's words. ENVIRONMENT.json has
    held the quote since the full run."""
    _need(FULL)
    env = json.loads(ENV.read_text())
    assert env["model"]["model_version_requested"] is None
    d = env["model"]["default_when_not_requested"]
    assert d["file"] == "alphagenome/models/v1/dna_model.py"
    assert "ALL_FOLDS = enum.auto()" in d["text"] and "Default model version" in d["text"]
    assert d["file_sha256_in_wheel_record"] is True and len(d["file_sha256"]) == 64
    assert d["wheel"] == f"alphagenome=={env['client']['version']}"
    assert f"{d['wheel']} --hash=sha256:{d['wheel_sha256']}" in REQ.read_text().splitlines()
    for name in (FULL, DRIFT):
        if (DIR / name).exists():
            assert _prov(name)["model_version_default"] == d, name


# --------------------------------------------------------------------------
# the scores
# --------------------------------------------------------------------------
@pytest.mark.parametrize("name", sorted(SCORE_FILES))
def test_the_composite_is_the_notebooks_combination(name):
    """The client's raw scores are float32, so the notebook's cell 13 adds and divides
    in float32: the components are float32 values, and the composite is
    (sites + usage) + junctions / 5 in float32, to the last bit."""
    import numpy as np
    _need(name)
    df = _scores(name)
    f32 = {c: df[c].to_numpy(dtype=np.float32) for c in COMPONENTS}
    for c in COMPONENTS:
        assert (df[c].to_numpy() == f32[c].astype(np.float64)).all(), c
    expect = (f32["SPLICE_SITES"] + f32["SPLICE_SITE_USAGE"]) + f32["SPLICE_JUNCTIONS"] / np.float32(5.0)
    assert (df[COMPOSITE].to_numpy() == expect.astype(np.float64)).all()
    for c in COMPONENTS:                 # unsigned scorers: no raw score below zero
        assert (df[f"min_raw_{c}"].dropna() >= 0).all(), c


def test_the_new_scores_cover_the_seven_training_genes_only():
    """E17.0 rule 10: nothing for DDX3X, TP53 or the five further genes before the
    thresholds are frozen. The notebook's own three examples are the one exception, and
    they are exactly those."""
    aset = pd.read_parquet(SET, columns=["variant_id", "gene"])
    from src import config as C
    assert set(aset.gene) == set(C.ASSAY_LABELS)
    for p in DIR.iterdir():
        assert p.name in SCORE_FILES | {f"{n}.provenance.json" for n in SCORE_FILES} | {
            "NOTEBOOK_SOURCE.json", "NOTEBOOK_CONTEXT.json", "ENVIRONMENT.json",
            "ATLAS_COLUMNS_CHECK.json"}, p.name
    for name in _present():
        df = _scores(name, as_text=True)
        if name == "examples_scores.tsv":
            cell15 = next(c["source"] for c in json.loads(SOURCE.read_text())["quoted_cells"]
                          if c["cell"] == 15)
            notebook = {f"{c}:{p}:{r}>{a}" for c, p, r, a in re.findall(
                r"chromosome='(\w+)',\s+position=(\d+),\s+reference_bases='(\w+)',\s+"
                r"alternate_bases='(\w+)'", cell15)}
            assert set(df.variant_id) == notebook and len(notebook) == 3
            continue
        assert set(df.variant_id) <= set(aset.variant_id), name
        assert set(df.gene) <= set(C.ASSAY_LABELS), name


def test_the_full_run_scores_every_row_of_the_analysis_set():
    _need(FULL)
    aset = pd.read_parquet(SET, columns=["variant_id", "gene", "stratum"])
    df = _scores(FULL, as_text=True)
    assert list(df.variant_id) == list(aset.variant_id)
    assert list(df.gene) == list(aset.gene) and list(df.stratum) == list(aset.stratum)
    prov = _prov(FULL)
    assert prov["n_variants"] == len(aset) and prov["genes"] == sorted(set(aset.gene))


def test_the_pilot_is_the_seeded_draw():
    S = _scorer()
    aset = pd.read_parquet(SET, columns=["variant_id", "gene", "stratum"])
    pilot, repeat, share = S.pilot_sample(aset)
    assert len(pilot) == S.PILOT_N and len(repeat) == S.REPEAT_N
    assert sorted(_scores("pilot_scores.tsv").variant_id) == pilot
    for name in ("repeat_scores.tsv", DRIFT):
        if (DIR / name).exists():
            assert sorted(_scores(name).variant_id) == repeat, name
    prov = json.loads((DIR / "pilot_scores.tsv.provenance.json").read_text())
    assert prov["sample"]["share"] == share


def test_the_ten_variants_scored_twice_agree_to_the_last_bit():
    _need("repeat_scores.tsv")
    pilot = _scores("pilot_scores.tsv", as_text=True).set_index("variant_id")
    repeat = _scores("repeat_scores.tsv", as_text=True).set_index("variant_id")
    cols = COMPONENTS + [COMPOSITE]
    assert (pilot.loc[repeat.index, cols] == repeat[cols]).all().all()


def test_the_pilot_scored_again_in_the_full_run_agrees_to_the_last_bit():
    """Gate 1 review R4: the full run has its own cache, so it scored the pilot's 120
    afresh; they equal the pilot's scores as written."""
    _need(FULL)
    pilot = _scores("pilot_scores.tsv", as_text=True).set_index("variant_id")
    full = _scores(FULL, as_text=True).set_index("variant_id")
    cols = COMPONENTS + [COMPOSITE]
    assert (full.loc[pilot.index, cols] == pilot[cols]).all().all()
    check = _prov(FULL)["pilot_cross_check"]
    assert check["identical"] == check["variants"] == len(pilot) and check["differing"] == []
    assert check["against_sha256"] == _sha(DIR / "pilot_scores.tsv")


def test_the_drift_check_scored_the_ten_afresh_after_the_full_run():
    """Gate 1 review R1: the ten were each a fresh call (none from a cache), made after
    the full run ended, and equal the pilot's scores as written."""
    _need(DRIFT)
    S = _scorer()
    prov, full = _prov(DRIFT), _prov(FULL)
    assert prov["n_scored_this_run"] == S.REPEAT_N and prov["n_from_cache"] == 0
    drift = _scores(DRIFT, as_text=True)
    end = datetime.fromisoformat(full["end_utc"])
    assert all(datetime.fromisoformat(t) > end for t in drift.scored_utc)
    assert prov["after_full_run"]["end_utc"] == full["end_utc"]
    assert prov["after_full_run"]["sha256"] == _sha(DIR / f"{FULL}.provenance.json")
    pilot = _scores("pilot_scores.tsv", as_text=True).set_index("variant_id")
    cols = COMPONENTS + [COMPOSITE]
    assert (pilot.loc[drift.variant_id, cols].to_numpy() == drift[cols].to_numpy()).all()
    assert prov["drift_check"]["identical"] == S.REPEAT_N


@pytest.mark.parametrize("name", sorted(SCORE_FILES))
def test_each_record_matches_its_scores_and_holds_no_key(name):
    prov_path = DIR / f"{name}.provenance.json"
    if not prov_path.exists():
        pytest.skip(f"{name} is written by a later step")
    prov = json.loads(prov_path.read_text())
    assert prov["output_sha256"] == _sha(DIR / name)
    assert prov["definition_sha256"] == _sha(SOURCE)
    assert prov["api_key_present"] is True
    assert prov["api_key_source"] in ("ALPHAGENOME_API_KEY", "~/.config/alphagenome/key")
    assert prov["failures"] == []
    text = prov_path.read_text() + ENV.read_text()
    assert not re.search(r"AIza[0-9A-Za-z_\-]{30,}", text)


# --------------------------------------------------------------------------
# the scoring stops when it should (gate 1 review R3): a fake client that only
# raises, and no call to the service
# --------------------------------------------------------------------------
def _hold_until_stopped(workers: int, error: BaseException, at: int = 5):
    """A fake client: calls before `at` return at once; call `at` waits until every
    worker holds a call, then raises `error`; the calls then in flight return once the
    run has stopped. Counts the calls, and those begun after the stop."""
    stop, cond = threading.Event(), threading.Condition()
    seen = {"calls": 0, "active": 0, "after_stop": 0, "all_held": None}

    def score(r):
        with cond:
            seen["calls"] += 1
            n = seen["calls"]
            seen["after_stop"] += stop.is_set()
            seen["active"] += 1
            cond.notify_all()
        try:
            if n < at:
                return {"score": 1.0}
            if n == at:
                with cond:
                    seen["all_held"] = cond.wait_for(lambda: seen["active"] == workers, 30)
                raise error
            stop.wait(30)
            return {"score": 1.0}
        finally:
            with cond:
                seen["active"] -= 1
    return score, stop, seen


@pytest.mark.parametrize("workers", [1, 2, 4])
def test_an_authentication_error_stops_the_run_after_the_calls_in_flight(workers):
    """The calls in flight end, nothing queued starts, and no call begins after the
    stop."""
    S = _scorer()
    score, stop, seen = _hold_until_stopped(
        workers, RuntimeError("UNAUTHENTICATED: the request does not have valid credentials"))
    rows = [SimpleNamespace(variant_id=f"v{i}") for i in range(200)]
    kept = []
    with pytest.raises(S.Stopped, match="authentication error") as e:
        S.run(rows, score, lambda r, out: kept.append(r.variant_id), workers=workers,
              backoff=0.0, stop=stop)
    assert seen["all_held"] is True
    assert seen["calls"] == 5 + workers - 1          # the fatal call and those in flight
    assert seen["after_stop"] == 0
    assert len(kept) == seen["calls"] - 1 and e.value.state["failures"] == []


@pytest.mark.parametrize("workers", [1, 2, 4])
def test_ctrl_c_stops_the_run_by_the_same_path(workers):
    """Ctrl+C, here a KeyboardInterrupt that reaches the loop through a future as Ctrl+C
    reaches the main thread. Until the loop sets the flag, the worker that raised it may
    begin one more call; so at most one call per worker is in flight at the stop, and
    nothing queued starts."""
    S = _scorer()
    score, stop, seen = _hold_until_stopped(workers, KeyboardInterrupt())
    rows = [SimpleNamespace(variant_id=f"v{i}") for i in range(200)]
    kept = []
    with pytest.raises(S.Stopped, match=re.escape("Ctrl+C")):
        S.run(rows, score, lambda r, out: kept.append(r.variant_id), workers=workers,
              backoff=0.0, stop=stop)
    assert seen["all_held"] is True
    assert 5 + workers - 1 <= seen["calls"] <= 5 + workers
    assert seen["after_stop"] <= 1
    assert len(kept) == seen["calls"] - 1


@pytest.mark.parametrize("workers", [1, 2, 4])
def test_three_variants_in_a_row_failing_stop_the_run(workers):
    """Errors that do not clear (here RESOURCE_EXHAUSTED on every call): the run stops
    once three variants in a row have used all their attempts, with at most one more
    variant per other worker begun, and at most one call each begun after the stop."""
    S = _scorer()
    stop, lock = threading.Event(), threading.Lock()
    seen = {"calls": 0, "after_stop": 0}

    def score(r):
        with lock:
            seen["calls"] += 1
            seen["after_stop"] += stop.is_set()
        raise RuntimeError("RESOURCE_EXHAUSTED: quota exceeded")

    rows = [SimpleNamespace(variant_id=f"v{i}") for i in range(200)]
    with pytest.raises(S.Stopped, match="3 variants in a row") as e:
        S.run(rows, score, lambda r, out: None, workers=workers, attempts=3, backoff=0.0,
              stop=stop)
    failed = e.value.state["failures"]
    assert 3 <= len(failed) <= 3 + workers - 1
    assert seen["calls"] <= (3 + workers - 1) * 3
    assert seen["after_stop"] <= workers - 1
    if workers == 1:
        assert [f["variant_id"] for f in failed] == ["v0", "v1", "v2"] and seen["calls"] == 9


def test_failures_not_in_a_row_do_not_stop_the_run():
    """A variant scored between failures starts the count again; the failures are
    returned, and the caller stops after the run."""
    S = _scorer()
    rows = [SimpleNamespace(variant_id=f"v{i}") for i in range(12)]
    bad = {"v0", "v1", "v3", "v4", "v6", "v7"}

    def score(r):
        if r.variant_id in bad:
            raise RuntimeError("UNAVAILABLE: the service is unavailable")
        return {"score": 1.0}

    kept = []
    state = S.run(rows, score, lambda r, out: kept.append(r.variant_id), workers=1,
                  attempts=2, backoff=0.0)
    assert [f["variant_id"] for f in state["failures"]] == sorted(bad)
    assert kept == [f"v{i}" for i in range(12) if f"v{i}" not in bad]
    assert state["retries"] == len(bad) and state["calls"] == 2 * len(bad) + len(kept)


# --------------------------------------------------------------------------
# the pipeline reads the new scores (E17.5)
# --------------------------------------------------------------------------
OLD = "alphagenome_atlas_quantile"
COMPARISON = PHASE1 / "reports/evidence/alphagenome_definition_comparison.csv"
SUPP = PHASE1 / "reports/evidence/supplement"


def test_the_analysis_set_reads_the_new_scores_bit_for_bit():
    """N1: the analysis set's `alphagenome` is all_scores.tsv's composite, read back
    exactly."""
    import numpy as np
    aset = pd.read_parquet(SET, columns=["variant_id", "alphagenome"])
    full = _scores(FULL).set_index("variant_id")[COMPOSITE]
    a = aset["alphagenome"].to_numpy(dtype=float)
    b = full.loc[aset["variant_id"]].to_numpy(dtype=float)
    assert np.array_equal(a.view(np.int64), b.view(np.int64))


def test_the_old_composite_is_kept_beside_it_and_outside_every_model():
    from src import evid_common as K
    aset = pd.read_parquet(SET)
    cols = list(aset.columns)
    assert cols.index(OLD) == cols.index("alphagenome") + 1
    assert OLD not in K.PANEL + K.OPTIONAL + K.FUSION_FEATURES


def test_the_old_composite_is_the_atlas_column():
    import numpy as np
    from src import evid_build_set as B
    if not B.ATLAS_MATRIX.exists():
        pytest.skip("the companion atlas release is not here")
    atlas = (pd.read_parquet(B.ATLAS_MATRIX, columns=["variant_id", "alphagenome"])
             .drop_duplicates("variant_id").set_index("variant_id")["alphagenome"])
    aset = pd.read_parquet(SET, columns=["variant_id", OLD])
    assert np.array_equal(aset[OLD].to_numpy(dtype=float),
                          atlas.loc[aset["variant_id"]].to_numpy(dtype=float), equal_nan=True)


def test_the_new_column_is_larger_is_damaging():
    """E17.5: measured, not assumed; no sign is flipped."""
    o = pd.read_csv(PHASE1 / "reports/evidence/feature_orientation.csv").set_index("feature")
    assert o.loc["alphagenome", "spearman_vs_pathogenicity"] > 0
    assert bool(o.loc["alphagenome", "oriented_larger_is_damaging"])


def test_the_comparison_table_holds_only_the_studys_own_quantities():
    """Three strata x three tiers x two sides. Held out: Table 3's k/n on the PP3 side;
    none on the BP4 side, which carries E3's median held-out ratio and says why."""
    from src import evid_practical_metrics as P
    t = pd.read_csv(COMPARISON, dtype=str, keep_default_na=False)
    assert len(t) == 18 and not t.duplicated(["stratum", "tier", "side"]).any()
    assert set(t.stratum) == {"s3_10", "s11_50", "s3_50"} and set(t.side) == {"PP3", "BP4"}
    pp3, bp4 = t[t.side == "PP3"], t[t.side == "BP4"]
    for n in ("alphagenome", OLD):
        assert (bp4[f"{n}_heldout_k_of_n"] == "").all() and (pp3[f"{n}_heldout_lr_median"] == "").all()
        assert (pp3[f"{n}_heldout_k_of_n"] != "").all()
    assert (bp4["note"] == P.BP4_NO_LOGO).all() and (pp3["note"] == "").all()


def test_the_comparison_tables_new_definition_is_the_studys():
    """Float for float what E3, E2.2 and E4 wrote for `alphagenome` (the stage checks
    E15's band ratios as well)."""
    from src import evid_common as K
    t = K.read_back(COMPARISON, keep_default_na=False)
    ev = K.read_back(PHASE1 / "reports/evidence/evidence_thresholds.csv")
    tm = K.read_back(PHASE1 / "reports/evidence/territory_metrics.csv")
    for r in t.itertuples():
        side = r.side.lower()
        e = ev[(ev.tool == "alphagenome") & (ev.stratum == r.stratum) & (ev.tier == r.tier)].iloc[0]
        got = getattr(r, "alphagenome_in_sample_threshold")
        want = e[f"{side}_threshold_insample"]
        assert (got == "" and pd.isna(want)) or float(got) == want, (r.stratum, r.tier, side)
        a = tm[(tm.tool == "alphagenome") & (tm.stratum == r.stratum) & (tm.clinvar_arm == "all")]
        assert float(r.alphagenome_pooled_auroc) == a["auroc"].iloc[0]


def test_s3_describes_the_new_definition():
    """The new column's row comes from E17's records, not from the atlas's resource
    row; the old composite has no row."""
    s3 = pd.read_csv(SUPP / "tableS3_predictor_columns.csv", dtype=str, keep_default_na=False)
    r = s3[s3["Score column"] == "alphagenome"].iloc[0]
    text = " ".join(r.to_numpy())
    for old in ("merged-quantile", "16-kb", "0.7.0"):
        assert old not in text, old
    src = json.loads(SOURCE.read_text())
    assert f"API client {src['client']['version']}" in r["Version, resource and scoring configuration"]
    assert src["commit"][:8] in r["Version, resource and scoring configuration"]
    assert OLD not in set(s3["Score column"])


def test_licence_data_section_2_names_every_alphagenome_output():
    """Gate 1 review R5 and the author's decision of 2026-10-22: one section for every
    AlphaGenome value, with the anchors the supplementary tables read it by."""
    text = (REPO / "LICENSE-DATA").read_text()
    head = next(l for l in text.splitlines() if l.startswith("2. "))
    flat = " ".join(text.split())
    sec2 = flat.split(" 2. AlphaGenome outputs")[1].split(" 3. What a reuser")[0]
    assert "AVI" in head and "NON-COMMERCIAL" in sec2
    for col in ["avi", "avi_splice_sites", "avi_splice_site_usage", "avi_splice_junctions",
                "alphagenome", OLD, "alphagenome_v061", *COMPONENTS, COMPOSITE]:
        assert f"`{col}`" in sec2, col
    for s in ["phase1/data/evidence/alphagenome_official/", "docs/alphagenome-official-pilot.md",
              "client 0.10.0", "dna_client.score_variant", "client 0.7.0", "client 0.6.1",
              "alphagenome_definition_comparison.csv"]:
        assert s in sec2, s


# --------------------------------------------------------------------------
# the AlphaGenome terms: no AlphaGenome column trains a model
# --------------------------------------------------------------------------
def test_no_alphagenome_column_enters_the_fusion_model():
    from src import evid_common as K
    assert not [f for f in K.FUSION_FEATURES if f.startswith("alphagenome")]
