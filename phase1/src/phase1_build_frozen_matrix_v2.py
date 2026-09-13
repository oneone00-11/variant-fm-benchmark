"""Build frozen-matrix-v2: frozen-matrix-v1 with three predictor columns replaced.

    spliceai  <- SpliceAI 1.3.1 at full float precision (the CLI rounds to 2 d.p. only
                 when printing; the atlas patched the format call and validated that
                 re-rounding reproduces the stock output)
    pangolin  <- Pangolin pinned to git 5cf94b8, full float precision, same recipe
    nt        <- Nucleotide Transformer scored under the atlas's pinned environment
                 (6,000-bp window, variant aligned inside one 6-mer)

Everything else -- variant set, row order, labels, offsets, the other seven predictors
and all annotation columns -- is copied from v1 unchanged, and that is asserted here and
in tests/test_frozen_matrix_v2.py.

Inputs (all tracked):
    data/frozen/frozen_matrix_v1.parquet            (verified against the v1 pin first)
    data/rescore/atlas_columns_v2.tsv               the atlas's columns for the 21,394 shared SNVs
    <repo>/data/rescore/{spliceai,pangolin,nt}_frozen_only16.tsv
                                                    the 16 RAD51C SNVs the atlas lacks, scored with
                                                    the same pinned scorers (scripts/phase4b_score_pinned.py)
    <repo>/data/tp53/{spliceai,pangolin,nt}_tp53.tsv
                                                    the 192 TP53 splice SNVs, same scorers

Outputs:
    data/frozen/frozen_matrix_v2.parquet + manifest_v2.json   (v1 is never touched)
    data/frozen/frozen_matrix_v2_column_report.tsv            per-column before/after audit
    data/external/tp53_splice_scored_v2.parquet               TP53 with all ten predictors

Run from phase1/:  python -m src.phase1_build_frozen_matrix_v2
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from . import config as C

PHASE1 = Path(__file__).resolve().parents[1]
REPO = PHASE1.parent
V1 = PHASE1 / "data" / "frozen" / "frozen_matrix_v1.parquet"
V1_SHA = "2a0e249b44906f11ed91ce4746aba7389d8fff70d791a17710a3e26ec66d3199"
ATLAS_COLS = PHASE1 / "data" / "rescore" / "atlas_columns_v2.tsv"
ONLY16 = {m: REPO / "data" / "rescore" / f"{m}_frozen_only16.tsv" for m in ("spliceai", "pangolin", "nt")}
TP53_V1 = PHASE1 / "data" / "external" / "tp53_splice_scored.parquet"
TP53_IN = {m: REPO / "data" / "tp53" / f"{m}_tp53.tsv" for m in ("spliceai", "pangolin", "nt")}
OUT = PHASE1 / "data" / "frozen" / "frozen_matrix_v2.parquet"
MANIFEST = PHASE1 / "data" / "frozen" / "manifest_v2.json"
REPORT = PHASE1 / "data" / "frozen" / "frozen_matrix_v2_column_report.tsv"
TP53_OUT = PHASE1 / "data" / "external" / "tp53_splice_scored_v2.parquet"

REPLACED = {"spliceai": "spliceai_ds_fullprec", "pangolin": "pangolin_fullprec", "nt": "nucleotide_transformer"}


def canonical_sha256(df: pd.DataFrame) -> str:
    """Same recipe as v1: rows sorted by variant_id, CSV serialisation."""
    canon = df.sort_values("variant_id").reset_index(drop=True).to_csv(index=False)
    return hashlib.sha256(canon.encode()).hexdigest()


def file_sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def max_spearman_given_ties(x: np.ndarray) -> float:
    """The atlas's arithmetic ceiling (atlas.robustness), copied verbatim so this build
    has no import dependency on the companion repository."""
    n = len(x)
    if n < 3:
        return float("nan")
    rx = stats.rankdata(x, method="average")
    if rx.std() == 0:
        return 0.0
    order = np.argsort(x, kind="stable")
    ideal = np.empty(n, dtype=float)
    ideal[order] = np.arange(1, n + 1, dtype=float)
    return float(np.corrcoef(rx, ideal)[0, 1])


def replacement_frame() -> pd.DataFrame:
    a = pd.read_csv(ATLAS_COLS, sep="\t")
    parts = [a[["variant_id", *REPLACED.values()]]]
    only = None
    for m, col in REPLACED.items():
        t = pd.read_csv(ONLY16[m], sep="\t")[["variant_id", col]]
        only = t if only is None else only.merge(t, on="variant_id", how="outer")
    parts.append(only)
    rep = pd.concat(parts, ignore_index=True)
    assert rep["variant_id"].is_unique, "duplicate variant ids across replacement inputs"
    return rep.set_index("variant_id")


def build() -> dict:
    v1 = pd.read_parquet(V1)
    got = canonical_sha256(v1)
    if got != V1_SHA:
        sys.exit(f"[v2] frozen-matrix-v1 does not match its pin ({got} != {V1_SHA}); refusing to build v2 on it")
    rep = replacement_frame()
    missing = set(v1["variant_id"]) - set(rep.index)
    extra = set(rep.index) - set(v1["variant_id"])
    assert not missing and not extra, f"replacement rows do not cover v1 exactly: missing {len(missing)}, extra {len(extra)}"

    v2 = v1.copy()
    report_rows = []
    for canon, col in REPLACED.items():
        new = rep.loc[v2["variant_id"], col].to_numpy(dtype=float)
        old = v2[canon].to_numpy(dtype=float)
        ok = ~(np.isnan(old) | np.isnan(new))
        sp = v2["is_splice"].to_numpy()
        report_rows.append({
            "column": canon, "n_scored_v1": int((~np.isnan(old)).sum()), "n_scored_v2": int((~np.isnan(new)).sum()),
            "spearman_v1_v2": float(stats.spearmanr(old[ok], new[ok]).statistic),
            "distinct_v1": int(pd.Series(old[ok]).nunique()), "distinct_v2": int(pd.Series(new[ok]).nunique()),
            "splice_distinct_v1": int(pd.Series(old[sp & ok]).nunique()), "splice_distinct_v2": int(pd.Series(new[sp & ok]).nunique()),
            "splice_tie_ceiling_v1": max_spearman_given_ties(old[sp & ok]),
            "splice_tie_ceiling_v2": max_spearman_given_ties(new[sp & ok]),
            "max_abs_diff_after_rounding_to_2dp": (float(np.nanmax(np.abs(np.round(new[ok], 2) - old[ok])))
                                                   if canon != "nt" else float("nan")),
        })
        v2[canon] = new
        v2[f"{canon}_isna"] = np.isnan(new).astype("int8")

    # invariants: same rows, same order, everything else identical
    assert (v2["variant_id"].to_numpy() == v1["variant_id"].to_numpy()).all()
    touched = set(REPLACED) | {f"{c}_isna" for c in REPLACED}
    for col in v1.columns:
        if col in touched:
            continue
        a, b = v1[col], v2[col]
        assert a.equals(b), f"column {col} changed during the v2 build"

    C.OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    v2.to_parquet(OUT, index=False)
    sha = canonical_sha256(v2)
    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip()
    except Exception:
        commit = "unknown"
    provenance = {
        "atlas_columns_v2.tsv": file_sha256(ATLAS_COLS),
        "atlas_columns_v2.provenance.json": json.loads((ATLAS_COLS.with_suffix(".provenance.json")).read_text()),
        **{f"{m}_frozen_only16.tsv": file_sha256(p) for m, p in ONLY16.items()},
        **{f"{m}_frozen_only16.provenance": json.loads(p.with_name(p.name.replace(".tsv", ".provenance.json")).read_text())
           for m, p in ONLY16.items() if p.with_name(p.name.replace(".tsv", ".provenance.json")).exists()},
    }
    built_at = datetime.now(timezone.utc).isoformat()
    if MANIFEST.exists():
        prior = json.loads(MANIFEST.read_text())
        if prior.get("sha256") == sha:        # unchanged matrix: keep the original build record
            built_at, commit = prior.get("built_at_utc", built_at), prior.get("git_commit", commit)
    manifest = {
        "version": "v2", "built_at_utc": built_at, "git_commit": commit,
        "derived_from": {"version": "v1", "sha256": V1_SHA},
        "replaced_columns": {
            "spliceai": "SpliceAI 1.3.1, -A grch38 -D 50, max of four deltas, full float precision",
            "pangolin": "Pangolin git 5cf94b8db938c658391b4305cd7ce33297d44ff7, -d 50, max(gain,|loss|), full float precision",
            "nt": "InstaDeepAI/nucleotide-transformer-v2-500m-multi-species, masked 6-mer LLR, 6,000-bp window",
        },
        "n_rows": int(len(v2)), "n_genes": int(v2["gene"].nunique()), "n_splice": int(v2["is_splice"].sum()),
        "sha256": sha, "output_path": str(OUT.relative_to(PHASE1)), "provenance": provenance,
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
    pd.DataFrame(report_rows).to_csv(REPORT, sep="\t", index=False, float_format="%.6g")
    return manifest


def build_tp53() -> dict:
    t1 = pd.read_parquet(TP53_V1)
    t2 = t1.copy()
    for canon, col in REPLACED.items():
        s = pd.read_csv(TP53_IN[canon], sep="\t").set_index("variant_id")[col]
        assert set(s.index) == set(t2["variant_id"]), f"TP53 {canon} rows do not match"
        new = s.loc[t2["variant_id"]].to_numpy(dtype=float)
        t2[canon] = new
        t2[f"{canon}_isna"] = np.isnan(new).astype("int8")
    assert (t2["variant_id"].to_numpy() == t1["variant_id"].to_numpy()).all()
    t2.to_parquet(TP53_OUT, index=False)
    info = {"n_rows": int(len(t2)), "sha256": canonical_sha256(t2),
            "inputs": {m: file_sha256(p) for m, p in TP53_IN.items()},
            "scored": {m: int(t2[m].notna().sum()) for m in REPLACED}}
    (TP53_OUT.with_suffix(".manifest.json")).write_text(json.dumps(info, indent=2) + "\n")
    return info


def main() -> None:
    m = build()
    print(json.dumps({k: v for k, v in m.items() if k != "provenance"}, indent=2))
    print(pd.read_csv(REPORT, sep="\t").to_string(index=False))
    t = build_tp53()
    print("TP53 v2:", json.dumps(t))


if __name__ == "__main__":
    main()
