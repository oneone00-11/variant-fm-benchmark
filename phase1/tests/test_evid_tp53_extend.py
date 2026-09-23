"""E7 guardrails on the TP53 table extended to 1 <= |offset| <= 12.

These check what a reader would otherwise have to take on trust: that the 192 rows
scored for the published study are carried unchanged, that the 96 added rows are
exactly the deposit's offsets 9 to 12, that every panel column is complete or its
gap is documented, that the added rows were scored by the same scorers (the
re-scored frozen rows reproduce), and that the table on disk is the one its
manifest and its score files describe.
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

PHASE1 = Path(__file__).resolve().parents[1]
REPO = PHASE1.parent
EXT = PHASE1 / "data/evidence/external"
TABLE = EXT / "tp53_splice_scored_12nt.parquet"
MANIFEST = EXT / "tp53_splice_scored_12nt.manifest.json"
FROZEN = PHASE1 / "data/external/tp53_splice_scored_v2.parquet"

if str(PHASE1) not in sys.path:
    sys.path.insert(0, str(PHASE1))

PANEL = ["spliceai", "spliceai_walker", "pangolin", "cadd", "phylop", "phastcons", "nt",
         "avi", "avi_splice_sites", "avi_splice_site_usage", "avi_splice_junctions"]

pytestmark = pytest.mark.skipif(
    not TABLE.exists(),
    reason="extended TP53 table not built (cd phase1 && python -m src.evid_tp53_extend)")


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


@pytest.fixture(scope="module")
def df():
    return pd.read_parquet(TABLE)


@pytest.fixture(scope="module")
def frozen():
    return pd.read_parquet(FROZEN)


@pytest.fixture(scope="module")
def manifest():
    return json.loads(MANIFEST.read_text())


def test_288_rows_sorted_and_unique(df):
    assert len(df) == 288
    assert df["variant_id"].is_unique
    assert df["variant_id"].is_monotonic_increasing


def test_24_rows_at_every_offset_1_to_12(df):
    counts = df["intron_offset"].abs().astype(int).value_counts().sort_index()
    assert counts.to_dict() == {k: 24 for k in range(1, 13)}


def test_frozen_rows_are_identical_to_v2(df, frozen):
    assert len(frozen) == 192
    part = df[df["variant_id"].isin(frozen["variant_id"])][list(frozen.columns)]
    pd.testing.assert_frame_equal(part.reset_index(drop=True),
                                  frozen.reset_index(drop=True), check_exact=True)


def test_new_rows_are_exactly_offsets_9_to_12(df, frozen):
    new = df[~df["variant_id"].isin(frozen["variant_id"])]
    assert len(new) == 96
    assert set(new["intron_offset"].abs().astype(int)) == {9, 10, 11, 12}
    assert (frozen["intron_offset"].abs() <= 8).all()


def test_panel_columns_complete_or_gap_documented(df, manifest):
    not_scored = manifest.get("not_scored", {})
    for c in PANEL:
        if c in not_scored:
            assert len(not_scored[c]) > 40, f"{c}: the reason for the gap is not recorded"
            continue
        assert c in df.columns, f"{c} missing from the table"
        assert df[c].notna().all(), f"{c}: {int(df[c].isna().sum())} unscored rows"
    for c in manifest["required_complete"]:
        assert df[c].notna().all(), f"{c}: declared complete but has gaps"


def test_the_out_of_panel_gap_is_documented(df, manifest):
    """AlphaMissense is empty for every intron-side SNV; the manifest must say why."""
    assert df["alphamissense"].isna().all()
    assert "missense" in manifest["not_scored"]["alphamissense"]


def test_new_rows_follow_the_frozen_derivation(df, frozen):
    new = df[~df["variant_id"].isin(frozen["variant_id"])]
    # TP53 is flipped (config.FLIP_GENES): pathogenicity is the raw score
    assert np.array_equal(new["func_pathogenicity"], new["func_score"])
    side = np.where(new["intron_offset"] > 0, "donor", "acceptor")
    assert (new["splice_side"].to_numpy() == side).all()
    assert new["is_splice"].all() and (new["region"] == "splice").all()
    assert (new["splice_subclass"] == "intronic").all()      # |offset| > SPLICE_REGION_MAX
    assert (new["offset_source"] == "upstream").all()
    for c in [c for c in df.columns if c.endswith("_isna")]:
        base = c[:-len("_isna")]
        assert (df[c].to_numpy() == df[base].isna().astype("int8").to_numpy()).all(), c
    ids = new["variant_id"].str.split("-", expand=True)
    assert (ids[0] == new["chrom"]).all()
    assert (ids[1].astype(int) == new["pos"]).all()
    assert (new["hgvs_nt"].str.startswith("NM_000546.6:c.")).all()


def test_table_matches_its_manifest(df, manifest):
    from src.phase1_build_frozen_matrix_v2 import canonical_sha256
    assert manifest["parquet_sha256"] == _sha(TABLE)
    assert manifest["sha256"] == canonical_sha256(df)
    assert manifest["n_rows"] == len(df) == manifest["n_frozen_rows"] + manifest["n_new_rows"]
    assert manifest["rows_per_abs_offset"] == {str(k): 24 for k in range(1, 13)}
    for c, v in manifest["columns"].items():
        assert v["n_non_null"] == int(df[c].notna().sum()), c


def test_added_rows_come_from_the_same_scorers(manifest):
    """Every scorer re-scored frozen rows beside the new ones; all must reproduce."""
    checks = manifest["reproduction_checks"]
    assert {c["scorer"] for c in checks} >= {"spliceai", "spliceai_walker", "pangolin", "nt",
                                             "cadd", "conservation", "gpn_msa", "avi",
                                             "alphagenome"}
    for c in checks:
        assert c["n_checked"] >= 24, c
        assert c["n_identical"] == c["n_checked"], c


def test_score_files_are_the_ones_the_manifest_records(manifest):
    for name, spec in manifest["score_files"].items():
        path = REPO / spec["file"]
        assert path.exists(), path
        assert manifest["inputs_sha256"][spec["file"]] == _sha(path), name
        prov = json.loads(Path(str(path) + ".provenance.json").read_text())
        recorded = prov.get("sha256") or prov.get("output_sha256")
        assert recorded == _sha(path), f"{name}: provenance hash is stale"
    for rel in ("data/external/tp53_mavedb_scores.tsv",
                "phase1/data/external/tp53_splice_scored_v2.parquet"):
        assert manifest["inputs_sha256"][rel] == _sha(REPO / rel), rel


def test_rebuild_from_cached_scores_is_byte_identical(tmp_path, monkeypatch):
    """No network and no model: the table is a function of the cached score files."""
    monkeypatch.chdir(PHASE1)
    from src import evid_tp53_extend as T
    table, _, _ = T.assemble()
    out = tmp_path / "rebuilt.parquet"
    table.to_parquet(out, index=False)
    assert _sha(out) == _sha(TABLE)


def test_inframe_subset_extends_the_frozen_one(manifest):
    inf = manifest["inframe"]
    sub = pd.read_parquet(REPO / inf["subset"])
    ev = pd.read_parquet(REPO / inf["events"])
    old = pd.read_parquet(EXT / "tp53_inframe_subset.parquet")
    assert set(old["variant_id"]) <= set(sub["variant_id"])
    assert set(ev["variant_id"]) == set(sub["variant_id"])
    assert (sub["spliceai_walker"] >= 0.2).all() and (sub["y_tp53_control_anchored"] == 0).all()
    assert inf["subset_sha256"] == _sha(REPO / inf["subset"])
    assert inf["events_sha256"] == _sha(REPO / inf["events"])


def test_no_credential_in_what_this_stage_wrote():
    """These outputs are not tracked yet, so tests/test_no_secrets.py cannot see them."""
    pat = re.compile(r"AIza[0-9A-Za-z_\-]{30,}")
    hits = [p.name for p in EXT.glob("tp53_*12nt*.json") if pat.search(p.read_text())]
    assert not hits, f"credential-shaped string in {hits}"
