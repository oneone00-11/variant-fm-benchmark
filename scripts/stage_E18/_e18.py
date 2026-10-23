"""Shared by the E18.2 scorer runs: where things are, the trial, the anchor check.

Every scorer run reads phase1/data/evidence/external/e18_score_input.parquet (the 39
anchors, then the external variants; scripts/stage_E18/score_input.py) and caches what
it scores under phase1/data/e18_cache/ (gitignored). A trial (--trial) scores the
anchors and the first TRIAL_N external variants into a cache of its own, to measure
throughput before the full run (E18.2: "先试吞吐"); the full run starts from an empty
cache, so its anchors are scored in the same run as its external variants.

An anchor is identical when the value scored here equals the analysis set's stored
value exactly, as float64, NaN equal to NaN (E18.0 rule 4(a)); the stored values are
e18_anchors.parquet's, copied from analysis_set_v1.parquet by scripts/stage_E18/anchors.py.
Every column of the analysis set reached it as the scorer's float64 output (the atlas's
score files and matrix are parquet; spliceai_walker.parquet and avi.parquet too), so
the comparison is direct; the probes of E18 section 4 compared the same way.

The atlas's scorers are the release's (v2.5.0, unpacked by
scripts/reproduce_evidence.py --fetch-inputs under phase1/data/evidence/companion_atlas/);
a run that loads a scorer from a checkout first checks it is byte for byte the
release's file.

Standard library, numpy and pandas only: each scorer's environment has them.
"""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from _runlog import REPO, rel

EXT = REPO / "phase1/data/evidence/external"
INPUT = EXT / "e18_score_input.parquet"
ANCHORS = EXT / "e18_anchors.parquet"
CACHE = REPO / "phase1/data/e18_cache"
RELEASE = REPO / "phase1/data/evidence/companion_atlas/atlas"
TRIAL_N = 50
LIMIT_HOURS = 24.0          # E18.2: a projected full run longer than this stops the stage


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def release_checked(path: Path, release_rel: str) -> str:
    """The sha256 of `path`, which must equal the release's `release_rel`."""
    got, want = sha(path), sha(RELEASE / release_rel)
    if got != want:
        raise SystemExit(f"{rel(path)} is not the atlas release's {release_rel}")
    return got


def load_input(trial: bool, keep=None) -> pd.DataFrame:
    """The anchors and the external variants (`keep(df)` selects among the latter);
    in a trial, the anchors and the first TRIAL_N of those."""
    df = pd.read_parquet(INPUT)
    anchors, ext = df[df["set"] == "anchor"], df[df["set"] == "external"]
    if keep is not None:
        ext = ext[keep(ext)]
    if trial:
        ext = ext.head(TRIAL_N)
    return pd.concat([anchors, ext], ignore_index=True)


def cache_path(name: str, trial: bool, suffix: str) -> Path:
    return CACHE / f"{name}{'_trial' if trial else ''}{suffix}"


def anchor_check(column: str, values: dict) -> dict:
    """values: variant_id -> scored value, for whichever anchors were scored."""
    stored = pd.read_parquet(ANCHORS).set_index("variant_id")[column]
    ids = [i for i in stored.index if i in values]
    got = np.array([values[i] for i in ids], dtype=float)
    want = stored.loc[ids].to_numpy(dtype=float)
    same = (got == want) | (np.isnan(got) & np.isnan(want))
    diff = np.abs(got - want)
    return {"column": column, "n_anchors": int(len(stored)), "n_checked": len(ids),
            "n_identical": int(same.sum()),
            "differing": [{"variant_id": i, "stored": repr(float(w)), "scored": repr(float(g))}
                          for i, g, w, s in zip(ids, got, want, same) if not s],
            "max_abs_diff": float(np.nanmax(diff)) if np.isfinite(diff).any() else 0.0}


def throughput(seconds: float, n_scored: int, n_full: int, processes: int = 1) -> dict:
    per = seconds / max(n_scored, 1)
    hours = per * n_full / processes / 3600
    return {"seconds": round(seconds, 1), "n_scored": n_scored,
            "seconds_per_variant": round(per, 3), "n_full_run": n_full,
            "processes": processes, "projected_full_run_hours": round(hours, 2),
            "within_limit": hours <= LIMIT_HOURS, "limit_hours": LIMIT_HOURS}


def report(log, check: dict | None, tp: dict | None) -> None:
    if check is not None:
        log.note("anchors", check)
        print(f"anchors {check['column']}: {check['n_identical']}/{check['n_checked']} identical"
              + (f" (of {check['n_anchors']})" if check["n_checked"] != check["n_anchors"] else "")
              + (f"; max |diff| {check['max_abs_diff']:.3g}" if check["differing"] else ""),
              flush=True)
    if tp is not None:
        log.note("throughput", tp)
        print(f"throughput: {tp['seconds_per_variant']} s per variant; the full run of "
              f"{tp['n_full_run']} on {tp['processes']} process(es) projects to "
              f"{tp['projected_full_run_hours']} h"
              + ("" if tp["within_limit"] else f" -- over {LIMIT_HOURS:.0f} h: stop and ask"),
              flush=True)


def write_json(path: Path, obj: dict) -> None:
    path.write_text(json.dumps(obj, indent=1, default=str) + "\n")


class Clock:
    def __init__(self):
        self.t0 = time.time()

    def seconds(self) -> float:
        return time.time() - self.t0
