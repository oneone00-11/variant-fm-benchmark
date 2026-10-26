"""E20.1 (C12): the true positives' event records.

The subsets are E6's and E2.7's definition with "called damaging" for "called normal", on
the 3-10 and 11-50 nt bands; the score inputs add anchor 2's ten variants per tool, drawn
from E6's event files at equal steps; every event file is the one its provenance names,
with the environment beside it; and both anchors hold: every row's largest score stands in
the established relation to the stored column, and every anchor gives E6's record back in
every field.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

PHASE1 = Path(__file__).resolve().parents[1]
REPO = PHASE1.parent
EV = PHASE1 / "data/evidence"
EXT = EV / "external"

BANDS = ["s3_10", "s11_50"]
SEVEN = {"spliceai": ("spliceai_walker", EV / "inframe_subset_tp.parquet",
                      EV / "inframe_score_input_tp.parquet", EV / "inframe_events_tp.parquet",
                      EV / "inframe_events.parquet"),
         "pangolin": ("pangolin", EV / "inframe_subset_pangolin_tp.parquet",
                      EV / "inframe_score_input_pangolin_tp.parquet",
                      EV / "inframe_events_pangolin_tp.parquet", EV / "inframe_events_pangolin.parquet")}
EXTERNAL = {"ddx3x": (EXT / "ddx3x_inframe_subset_tp.parquet", EXT / "ddx3x_inframe_events_tp.parquet",
                      EXT / "ddx3x_scored.parquet"),
            "tp53": (EXT / "tp53_inframe_subset_12nt_tp.parquet",
                     EXT / "tp53_inframe_events_12nt_tp.parquet", EXT / "tp53_splice_scored_12nt.parquet")}
EVENT_FILES = [v[3] for v in SEVEN.values()] + [v[1] for v in EXTERNAL.values()]


def _set():
    return pd.read_parquet(EV / "analysis_set_v1.parquet")


def test_the_seven_genes_subsets_follow_the_definition():
    df = _set()
    for tool, (col, subset, *_rest) in SEVEN.items():
        want = df[(df[col] >= 0.2) & (df["y_assay"] == 1) & df["stratum"].isin(BANDS)]
        got = pd.read_parquet(subset)
        assert list(got["variant_id"]) == list(want["variant_id"])
        assert set(got["stratum"]) <= set(BANDS) and (got["score_column"] == col).all()
    assert len(pd.read_parquet(SEVEN["spliceai"][1])) == 553
    assert len(pd.read_parquet(SEVEN["pangolin"][1])) == 539


def test_the_external_subsets_are_the_union_of_damaging_labels():
    d = pd.read_parquet(EXT / "ddx3x_splice.parquet").merge(
        pd.read_parquet(EXT / "ddx3x_scored.parquet"), on="variant_id", how="left")
    labels = ["y_deposit", "y_control_anchored_gated", "y_control_anchored", "y_fdr"]
    damaging = np.zeros(len(d), dtype=bool)
    for c in labels:
        damaging |= (d[c] == 1).to_numpy()
    want = d[(d["spliceai_walker"] >= 0.2).to_numpy() & damaging & d["stratum"].isin(BANDS).to_numpy()]
    got = pd.read_parquet(EXTERNAL["ddx3x"][0])
    assert list(got["variant_id"]) == list(want["variant_id"]) and len(got) == 66
    t = pd.read_parquet(EXT / "tp53_splice_scored_12nt.parquet")
    t = t[t["is_splice"]]
    off = t["intron_offset"].abs()
    band = np.where(off <= 2, "pm12", np.where(off <= 10, "s3_10", "s11_50"))
    keep = (t["spliceai_walker"] >= 0.2).to_numpy() & (t["func_pathogenicity"].to_numpy(float) > 0) \
        & np.isin(band, BANDS)
    got = pd.read_parquet(EXTERNAL["tp53"][0])
    assert list(got["variant_id"]) == list(t.loc[keep, "variant_id"]) and len(got) == 15


def _added_in_e20b(tool: str) -> set:
    """The variants E20b added to E6's subsets (its subset run record): E20.1 drew the anchors
    from E6's event files as they stood before, 609 and 524 variants."""
    rec = sorted((REPO / "scripts/stage_E20b/run_log").glob("write_e6_subsets_*.json"))
    return set(json.loads(rec[-1].read_text())["notes"]["checks"][tool]["added"]) if rec else set()


def test_the_score_inputs_add_the_anchors_at_equal_steps():
    for tool, (_col, subset, inp, _new, old) in SEVEN.items():
        s, i = pd.read_parquet(subset), pd.read_parquet(inp)
        tp = i[i["role"] == "true_positive"]
        assert list(tp["variant_id"]) == list(s["variant_id"])
        ev = pd.read_parquet(old)
        ev = ev[~ev["variant_id"].isin(_added_in_e20b(tool))]
        ev = ev.sort_values("variant_id", kind="mergesort").reset_index(drop=True)
        n = len(ev)
        assert n == {"spliceai": 609, "pangolin": 524}[tool]
        want = [ev.loc[k * n // 10, "variant_id"] for k in range(10)]
        assert list(i.loc[i["role"] == "anchor", "variant_id"]) == want


def test_every_event_file_is_the_one_its_records_name():
    for p in EVENT_FILES:
        prov = json.loads(Path(str(p) + ".provenance.json").read_text())
        env = json.loads(Path(str(p) + ".environment.json").read_text())
        assert prov["output_sha256"] == hashlib.sha256(p.read_bytes()).hexdigest()
        assert prov["n_scored"] == prov["n_variants"] == len(pd.read_parquet(p))
        assert prov["distance"] == (50 if "pangolin" in p.name else 4999)
        assert env["file"] == str(p.relative_to(REPO)) and env["pip_freeze"] and env["command"]
        assert env["scorer_sha256"] == hashlib.sha256((REPO / env["scorer"]).read_bytes()).hexdigest()


def _cli_round(x: float) -> float:
    return float(round(np.float32(x), 2))


def test_anchor_1_every_row_stands_in_the_relation_to_the_column():
    df = _set().set_index("variant_id")
    for tool, (col, _s, _i, new, old) in SEVEN.items():
        for p in (old, new):
            ev = pd.read_parquet(p)
            stored = df.loc[ev["variant_id"], col].to_numpy(dtype=float)
            if tool == "spliceai":
                assert (ev["spliceai_event_max"].to_numpy() == stored).all()
            else:
                want = np.array([_cli_round(x) for x in stored])
                assert (ev["pangolin_event_max"].to_numpy() == want).all()
    for gene, (_s, new, scored) in EXTERNAL.items():
        table = pd.read_parquet(scored).set_index("variant_id")
        ev = pd.read_parquet(new)
        assert (ev["spliceai_event_max"].to_numpy()
                == table.loc[ev["variant_id"], "spliceai_walker"].to_numpy(dtype=float)).all()


def test_anchor_2_the_anchors_give_e6s_records_back():
    for tool, (_col, _s, inp, new, old) in SEVEN.items():
        anc = pd.read_parquet(inp).query("role == 'anchor'")["variant_id"]
        a = pd.read_parquet(old).set_index("variant_id").loc[anc]
        b = pd.read_parquet(new).set_index("variant_id").loc[anc]
        for c in a.columns:
            if a[c].dtype.kind == "f":
                assert np.array_equal(a[c].to_numpy().view(np.int64), b[c].to_numpy().view(np.int64)), c
            else:
                assert (a[c].astype(str) == b[c].astype(str)).all(), c


def test_the_pangolin_work_folder_is_complete():
    work = EV / "pangolin_events_work_tp/chunks"
    chunks = sorted(work.glob("chunk_*.csv"))
    chunks = [c for c in chunks if not c.name.endswith(".out.csv")]
    assert len(chunks) == 5
    for c in chunks:
        assert c.with_suffix(".done").read_text() == hashlib.md5(c.read_bytes()).hexdigest()
        assert c.with_suffix(".out.csv").is_file()
