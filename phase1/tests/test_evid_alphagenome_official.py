"""E17: AlphaGenome's splicing score by its official definition, and the records it rests on.

The definition is quoted from AlphaGenome's own notebook at a pinned commit, item by
item with the lines it rests on; the scoring environment is pinned with hashes; the
composite is the notebook's combination of its three components; the scores cover the
seven training genes only, so that no external gene is scored before the thresholds are
frozen (E17.0 rule 10); the pilot is the seeded draw it says it is, and its ten variants
scored twice agree to the last bit; no record holds the API key; and no AlphaGenome
column enters the fusion model (the AlphaGenome terms; E17.0 rule 5).
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import re
import sys
from pathlib import Path

import pandas as pd
import pytest

PHASE1 = Path(__file__).resolve().parents[1]
REPO = PHASE1.parent
DIR = PHASE1 / "data/evidence/alphagenome_official"
SOURCE = DIR / "NOTEBOOK_SOURCE.json"
ENV = DIR / "ENVIRONMENT.json"
REQ = REPO / "scripts/stage_E17/requirements.txt"
SET = PHASE1 / "data/evidence/analysis_set_v1.parquet"
COMPONENTS = ["SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS"]
COMPOSITE = "alphagenome_splicing"
SCORE_FILES = {"examples_scores.tsv", "pilot_scores.tsv", "repeat_scores.tsv", "scores.tsv"}

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


def _present() -> list[str]:
    return sorted(p.name for p in DIR.glob("*_scores.tsv") if p.name in SCORE_FILES)


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


# --------------------------------------------------------------------------
# the scores
# --------------------------------------------------------------------------
@pytest.mark.parametrize("name", sorted(SCORE_FILES))
def test_the_composite_is_the_notebooks_combination(name):
    """The client's raw scores are float32, so the notebook's cell 13 adds and divides
    in float32: the components are float32 values, and the composite is
    (sites + usage) + junctions / 5 in float32, to the last bit."""
    import numpy as np
    if not (DIR / name).exists():
        pytest.skip(f"{name} is written by a later step")
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
            "NOTEBOOK_SOURCE.json", "ENVIRONMENT.json", "ATLAS_COLUMNS_CHECK.json"}, p.name
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


def test_the_pilot_is_the_seeded_draw():
    S = _scorer()
    aset = pd.read_parquet(SET, columns=["variant_id", "gene", "stratum"])
    pilot, repeat, share = S.pilot_sample(aset)
    assert len(pilot) == S.PILOT_N and len(repeat) == S.REPEAT_N
    assert sorted(_scores("pilot_scores.tsv").variant_id) == pilot
    if (DIR / "repeat_scores.tsv").exists():
        assert sorted(_scores("repeat_scores.tsv").variant_id) == repeat
    prov = json.loads((DIR / "pilot_scores.tsv.provenance.json").read_text())
    assert prov["sample"]["share"] == share


def test_the_ten_variants_scored_twice_agree_to_the_last_bit():
    if not (DIR / "repeat_scores.tsv").exists():
        pytest.skip("the determinism check has not run yet")
    pilot = _scores("pilot_scores.tsv", as_text=True).set_index("variant_id")
    repeat = _scores("repeat_scores.tsv", as_text=True).set_index("variant_id")
    cols = COMPONENTS + [COMPOSITE]
    assert (pilot.loc[repeat.index, cols] == repeat[cols]).all().all()


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
# the AlphaGenome terms: no AlphaGenome column trains a model
# --------------------------------------------------------------------------
def test_no_alphagenome_column_enters_the_fusion_model():
    from src import evid_common as K
    assert not [f for f in K.FUSION_FEATURES if f.startswith("alphagenome")]
