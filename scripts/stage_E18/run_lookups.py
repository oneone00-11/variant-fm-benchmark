"""E18.2: CADD, phyloP, phastCons and GPN-MSA for the external genes, by the analysis
set's own scorers.

The analysis set's columns are the companion atlas's, from its models/cadd/score.py
(CADD GRCh38-v1.7, the precomputed whole-genome SNV file by remote tabix),
models/conservation/score.py (the UCSC track API, phyloP100way and phastCons100way) and
models/gpn_msa/score.py (the precomputed GPN-MSA hg38 scores of
songlab/gpn-msa-hg38-scores, flipped so that larger is more damaging, the analysis set's
direction). Each runs here unchanged, in the atlas's pinned environment (requests
2.32.5, pandas 3.0.5, pyarrow 21.0.0, pysam 0.24.0), on the 39 anchors and the external
variants that need it, in one run: DDX3X's 1,857 and the five further genes' 4,522.
TP53's four columns come from the published study's scorers, which
docs/column-provenance.md section 2 records as value for value the atlas's (evidence b),
so TP53 is not scored again here.

A real run of these scripts writes its run log beside itself, so each runs from a copy
in phase1/data/e18_cache/atlas_scorers/<scorer>/ (gitignored), byte for byte the
release's file, with GPN-MSA's pinned tabix index beside it, and with a cache of its
own, empty at the start; the atlas checkout is not written to. Before scoring, the
remote sources' versions are read and recorded: GPN-MSA's dataset revision must be
cf1718a9 (the revision the analysis set's column records) and its file the size the
scorer pins, or the run stops; CADD's file headers and the UCSC tracks' data times are
recorded. Outputs: phase1/data/e18_cache/lookup_<column>.parquet, the scripts' own run
logs beside the copies, lookups.sources.json and an environment record. --trial scores
the anchors and the first 50 of these external variants, with caches of its own.

Run from the repository root, in the atlas's pinned environment:
    <atlas checkout>/.venv/bin/python scripts/stage_E18/run_lookups.py
"""
from __future__ import annotations

import argparse
import importlib.util
import platform
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import RunLog, rel, utc_now  # noqa: E402
import _e18 as E  # noqa: E402

COPIES = E.CACHE / "atlas_scorers"
# column -> (scorer folder, extra arguments, the scorer's output column)
RUNS = {"cadd": ("cadd", [], "cadd"),
        "phylop": ("conservation", ["--track", "phyloP100way"], "phylop100way"),
        "phastcons": ("conservation", ["--track", "phastCons100way"], "phastcons100way"),
        "gpn_msa": ("gpn_msa", [], "gpn_msa")}
GENES = {"DDX3X", "CTCF", "SFPQ", "RAD51D", "XRCC2", "SBDS"}
GPN_REVISION = "cf1718a9"       # docs/column-provenance.md section 1, the analysis set's


def _load(path: Path):
    spec = importlib.util.spec_from_file_location(path.parent.name + "_score", path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def sources(df) -> dict:
    """The remote sources' versions, as served now."""
    import requests
    out = {"read_utc": utc_now()}
    gpn = _load(COPIES / "gpn_msa" / "score.py")
    h = requests.head(gpn.RESOLVE_URL, allow_redirects=False, timeout=60).headers
    out["gpn_msa"] = {"url": gpn.RESOLVE_URL, "x_repo_commit": h.get("X-Repo-Commit"),
                      "x_linked_size": h.get("X-Linked-Size"),
                      "x_linked_etag": h.get("X-Linked-Etag")}
    if not str(h.get("X-Repo-Commit", "")).startswith(GPN_REVISION) or \
            int(h.get("X-Linked-Size") or -1) != gpn.TOTAL_SIZE:
        raise SystemExit(f"GPN-MSA is not revision {GPN_REVISION} at the pinned size: {out['gpn_msa']}")
    cadd = _load(COPIES / "cadd" / "score.py")
    out["cadd"] = {"api_version": cadd.API_VERSION, "files": {}}
    for url in (cadd.MIRRORS["US"], cadd.MIRRORS["US"] + ".tbi"):
        h = requests.head(url, allow_redirects=True, timeout=60).headers
        out["cadd"]["files"][url] = {k: h.get(k) for k in ("Content-Length", "Last-Modified", "ETag")}
    cons = _load(COPIES / "conservation" / "score.py")
    a = df[df["set"] == "anchor"].iloc[0]
    out["ucsc"] = {}
    for track in cons.TRACKS:
        r = requests.get(cons.API, params={"genome": "hg38", "track": track,
                                           "chrom": f"chr{a.chrom}", "start": int(a.pos) - 1,
                                           "end": int(a.pos)}, timeout=120)
        r.raise_for_status()
        j = r.json()
        out["ucsc"][track] = {k: j.get(k) for k in ("dataTime", "dataTimeStamp", "trackType",
                                                    "bigDataUrl", "genome")}
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--trial", action="store_true",
                    help=f"the anchors and the first {E.TRIAL_N} external variants only")
    args = ap.parse_args()
    import pandas as pd
    tag = "_trial" if args.trial else ""
    with RunLog(__file__) as log:
        log.input(E.INPUT)
        df = E.load_input(args.trial, keep=lambda x: x["gene"].isin(GENES))
        E.CACHE.mkdir(parents=True, exist_ok=True)
        inp = E.cache_path("lookup_input", args.trial, ".parquet")
        df[["variant_id", "chrom", "pos", "ref", "alt"]].to_parquet(inp, index=False)
        log.output(inp)
        shas = {}
        for folder in sorted({f for f, _, _ in RUNS.values()}):
            dst = COPIES / folder
            dst.mkdir(parents=True, exist_ok=True)
            for name in ["score.py"] + (["gpn_scores.tsv.bgz.tbi"] if folder == "gpn_msa" else []):
                src = E.RELEASE / "models" / folder / name
                shutil.copyfile(src, dst / name)
                shas[f"models/{folder}/{name}"] = E.release_checked(dst / name, f"models/{folder}/{name}")
                log.input(src)
            cache = dst / f"cache{tag}"
            if cache.exists():
                shutil.rmtree(cache)
        src_versions = sources(df)
        srcf = E.cache_path("lookups", args.trial, ".sources.json")
        E.write_json(srcf, src_versions)
        log.output(srcf)
        log.note("sources", src_versions)
        start, clock, checks = utc_now(), E.Clock(), {}
        for col, (folder, extra, out_col) in RUNS.items():
            out = E.cache_path(f"lookup_{col}", args.trial, ".parquet")
            cmd = [sys.executable, str(COPIES / folder / "score.py"), "--input", str(inp),
                   "--output", str(out), "--cache-dir", str(COPIES / folder / f"cache{tag}"), *extra]
            print("$", " ".join(rel(c) if c.startswith("/") else c for c in cmd), flush=True)
            subprocess.run(cmd, check=True)
            log.output(out)
            got = pd.read_parquet(out)
            checks[col] = E.anchor_check(col, dict(zip(got["variant_id"], got[out_col].astype(float))))
            E.report(log, checks[col], None)
            log.note(f"n_scored_{col}", int(len(got)))
        for folder in sorted({f for f, _, _ in RUNS.values()}):
            for rl in sorted((COPIES / folder).glob("run_log*.json")):
                log.output(rl)
        seconds = clock.seconds()
        log.note("anchors", checks)
        if args.trial:
            E.report(log, None, E.throughput(
                seconds, len(df), len(E.load_input(False, keep=lambda x: x["gene"].isin(GENES)))))
        freeze = subprocess.run([sys.executable, "-m", "pip", "freeze"], check=True,
                                capture_output=True, text=True).stdout.split()
        env = {"scorers": "the atlas's models/{cadd,conservation,gpn_msa}/score.py, run "
                          "unchanged from byte-identical copies of the release's files",
               "scorer_sha256": shas, "interpreter": sys.executable,
               "python": platform.python_version(), "platform": platform.platform(),
               "pip_freeze": freeze, "run_started_utc": start, "run_ended_utc": utc_now(),
               "seconds": round(seconds, 1), "n_variants": len(df)}
        envf = E.cache_path("lookups", args.trial, ".environment.json")
        E.write_json(envf, env)
        log.output(envf)
        print(f"done: {len(df)} variants, {', '.join(RUNS)}")


if __name__ == "__main__":
    main()
