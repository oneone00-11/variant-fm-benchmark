"""E18.2: the external genes' committed score files, from the scorer runs' caches.

For each scorer run of E18.2 (run_spliceai.py --mode d50 and --mode walker,
run_pangolin.py, run_lookups.py, run_avi.py, and NT's run on the GPU) this

  1. checks the run scored exactly the 39 anchors and the variants it was asked for,
     once each, and that every anchor equals the analysis set's stored value bit for
     bit in every column of the run (E18.0 rule 4(a)): the value as the scorer wrote it
     against analysis_set_v1.parquet's, float64 against float64, NaN equal to NaN; the
     analysis set's columns reached it as their scorers' float64 output (the atlas's
     score files and matrix are parquet, as are spliceai_walker.parquet and avi.parquet),
     so that is the storage path. A difference, a missing anchor or a variant scored
     twice stops the build;
  2. writes, in phase1/data/evidence/external/, one parquet per gene with that gene's
     variants (as many rows as the run scored for it, in e18_score_input.parquet's
     order) and one with the run's anchors, e18_anchors_<run>.parquet, each with a
     .provenance.json and an .environment.json beside it, as the tp53_12nt_* files have;
     the five further genes' Atlas columns go to one file, newgenes_avi.parquet. Names
     say what a file holds beside the existing ones (E18.0 rule 7): <gene>_<run> for a
     column the gene did not have, ddx3x_pangolin_fullprec for the full-precision
     Pangolin, <gene>_<run>_rescored for a column scored again by the analysis set's
     scorer;
  3. assembles the five further genes' 4,522 variants (3 <= |offset| <= 50) with every
     column, newgenes_scored.parquet, and newgenes_scored.manifest.json; and
  4. generates docs/external-column-identity.md, the same-variable evidence table:
     7 external genes x 13 columns, each with its scorer, version, precision,
     direction and evidence, (a) anchors or (b) a committed value-for-value record.

--offline does steps 3 and 4 from the committed files of step 2 alone, reading no cache
and no network, and fails unless every file it writes is byte for byte the committed
one. Nothing written carries the time of the build, so a rebuild is byte-identical.
NT was scored on an RTX A6000 in the original run's recorded environment (the author's
decision of 2026-10-23), and its anchors decide alone: all 39 identical, and its files
are written like any other run's. Otherwise, by the rule the author fixed before the run,
NT does not enter the external test, with no tolerance, and no NT file per gene is
written; the run's evidence goes to phase1/data/evidence/external/nt_gpu_run/: the pod's
environment record, its three logs and both manifests as served, the 39 anchors' stored
and re-scored values with their difference (anchor_check.tsv), and record.json, which
says what was committed and names the two result tables that were not by their sha256
(the 6,667 external values are not the analysis set's variable and stay in the
gitignored cache; the author's decision of 2026-10-23). Every other run stops on an
anchor that differs.

Run in the reference image (修改单 §0 rule 9; scripts/stage_E18/in_image.sh):
    docker run ... IMAGE bash scripts/stage_E18/in_image.sh scripts/stage_E18/build_columns.py [--offline]
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel  # noqa: E402

EXT = REPO / "phase1/data/evidence/external"
CACHE = REPO / "phase1/data/e18_cache"
RUN_LOGS = REPO / "scripts/stage_E18/run_log"
INPUT = EXT / "e18_score_input.parquet"
SET = REPO / "phase1/data/evidence/analysis_set_v1.parquet"
AVI_RUN = REPO / "phase1/data/evidence/avi.parquet"
AG_EXTERNAL = REPO / "phase1/data/evidence/alphagenome_official/external_scores.tsv"
AG_DRIFT = {"start": REPO / "phase1/data/evidence/alphagenome_official/drift_2b_start_scores.tsv",
            "end": REPO / "phase1/data/evidence/alphagenome_official/drift_2b_end_scores.tsv"}
COLUMN_PROVENANCE = REPO / "docs/column-provenance.md"
DOC = REPO / "docs/external-column-identity.md"
SCORED = EXT / "newgenes_scored.parquet"
MANIFEST = EXT / "newgenes_scored.manifest.json"
# NT's evidence when its anchors do not hold (the author's decisions of 2026-10-23)
NT_DIR = EXT / "nt_gpu_run"
NT_SERVED = ["out/environment.json", "out/MANIFEST.sha256", "out_trial/MANIFEST.sha256",
             "run.log", "run_trial.log", "run_trial_failed_hf_transfer.log"]
NT_RULE = ("by the rule the author fixed before the run (2026-10-23), NT enters the external "
           "test only if all 39 anchors are identical bit for bit; otherwise it does not, "
           "with no tolerance and no re-scoring of the seven genes, and every external "
           "gene's NT cell reads: no column of the same variable")

NEW = ["CTCF", "SFPQ", "RAD51D", "XRCC2", "SBDS"]
GENES = ["DDX3X", "TP53"] + NEW
LOOKUP_GENES = ["DDX3X"] + NEW
AVI_GENES = ["TP53"] + NEW
KEY = ["variant_id", "gene", "chrom", "pos", "ref", "alt"]
WALKER = ["spliceai_walker", "spliceai_walker_ds_ag", "spliceai_walker_ds_al",
          "spliceai_walker_ds_dg", "spliceai_walker_ds_dl"]
AVI = ["avi", "avi_splice_sites", "avi_splice_site_usage", "avi_splice_junctions"]
# the thirteen predictor columns, in scripts/stage_E18/anchors.py's order
COLUMNS = ["spliceai", "spliceai_walker", "pangolin", "alphagenome", "nt", "cadd",
           "phylop", "phastcons", "gpn_msa"] + AVI

# run -> what it is. "cache": the files run_*.py wrote; "logs": the run log name and the
# argument that marks the run (trial runs excluded).
RUNS = {
    "spliceai": dict(
        columns=["spliceai"], genes=GENES, cache=["spliceai_d50.tsv"],
        environments=["spliceai_d50.environment.json"], logs=("run_spliceai", "d50"),
        tool="SpliceAI 1.3.1 (tensorflow 2.19.1, CPU)",
        scorer="the atlas's models/spliceai/score_fullprec.py",
        configuration={"annotation": "grch38", "distance": 50, "mask": 0,
                       "score_definition": "max(DS_AG, DS_AL, DS_DG, DS_DL)"},
        precision="full float precision (the CLI's {:.2f} replaced by {:.17g})",
        direction="larger = more damaging, as the analysis set's column"),
    "spliceai_walker": dict(
        columns=WALKER, genes=GENES, cache=["spliceai_walker.tsv"],
        environments=["spliceai_walker.environment.json"], logs=("run_spliceai", "walker"),
        tool="SpliceAI 1.3.1 (tensorflow 2.19.1, CPU)",
        scorer="phase1/src/evid_score_spliceai_walker.py",
        configuration={"annotation": "grch38", "distance": 4999, "mask": 0,
                       "score_definition": "max(DS_AG, DS_AL, DS_DG, DS_DL); the four "
                                           "components at the gene of the maximum"},
        precision="full float precision",
        direction="larger = more damaging, as the analysis set's column"),
    "pangolin": dict(
        columns=["pangolin"], genes=GENES, cache=["pangolin_w*.tsv"],
        environments=["pangolin_w*.environment.json"], logs=("run_pangolin", None),
        tool="Pangolin at git 5cf94b8 (torch 2.13.0, CPU, one thread)",
        scorer="the atlas's models/pangolin/score_fullprec.py",
        configuration={"distance": 50, "mask": "True (Pangolin's default)",
                       "score_definition": "max(splice gain, |splice loss|) over the "
                                           "transcript records"},
        precision="full float precision (the two round(..., 2) removed)",
        direction="larger = more damaging, as the analysis set's column"),
    "cadd": dict(
        columns=["cadd"], genes=LOOKUP_GENES, cache=["lookup_cadd.parquet"],
        environments=["lookups.environment.json"], logs=("run_lookups", None),
        tool="CADD GRCh38-v1.7, the precomputed whole-genome SNV file (remote tabix)",
        scorer="the atlas's models/cadd/score.py",
        configuration={"source": "precomputed", "mirror": "US"},
        precision="as served (PHRED as printed in the file)",
        direction="larger = more damaging, as the analysis set's column"),
    "conservation": dict(
        columns=["phylop", "phastcons"], genes=LOOKUP_GENES,
        cache=["lookup_phylop.parquet", "lookup_phastcons.parquet"],
        environments=["lookups.environment.json"], logs=("run_lookups", None),
        tool="UCSC phyloP100way and phastCons100way, hg38 (track API)",
        scorer="the atlas's models/conservation/score.py",
        configuration={"tracks": ["phyloP100way", "phastCons100way"]},
        precision="as served",
        direction="larger = more conserved = more damaging, as the analysis set's columns"),
    "gpn_msa": dict(
        columns=["gpn_msa"], genes=LOOKUP_GENES, cache=["lookup_gpn_msa.parquet"],
        environments=["lookups.environment.json"], logs=("run_lookups", None),
        tool="GPN-MSA hg38 precomputed scores, songlab/gpn-msa-hg38-scores at cf1718a9",
        scorer="the atlas's models/gpn_msa/score.py",
        configuration={"revision": "cf1718a926db35d078bcb08b9dcebc803494bf37"},
        precision="as served",
        direction="gpn_msa = -raw: larger = more damaging, as the analysis set's column"),
    "avi": dict(
        columns=AVI, genes=AVI_GENES, cache=["avi_e18.parquet"],
        environments=["avi_e18.environment.json"], logs=("run_avi", None),
        tool="AlphaGenome Atlas, client 0.9.0, scorers AVI_SCORE, SPLICE_SITES, "
             "SPLICE_SITE_USAGE, SPLICE_JUNCTIONS",
        scorer="phase1/src/evid_score_avi.py::run",
        configuration={"aggregation": "the maximum over tracks and over rows"},
        precision="as returned (float32), stored as float64",
        direction="larger = more damaging, as the analysis set's columns"),
    "nt": dict(
        columns=["nt"], genes=GENES, cache=["nt/out/results.tsv"],
        environments=["nt/out/environment.json"], logs=None,
        tool="InstaDeepAI/nucleotide-transformer-v2-500m-multi-species at 06615c16, "
             "NVIDIA RTX A6000, torch 2.11.0+cu128, transformers 4.46.3",
        scorer="the atlas's models/nt/score.py (masked_llr)",
        configuration={"window_nt": 6000, "dtype": "float32", "batch": 1},
        precision="float32 model output as float64",
        direction="larger = more damaging, as the analysis set's column", optional=True),
}
COLUMN_RUN = {c: r for r, spec in RUNS.items() for c in spec["columns"] if c in COLUMNS}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def file_name(run: str, gene: str) -> str:
    g = gene.lower()
    if gene in NEW:
        return "newgenes_avi.parquet" if run == "avi" else f"{g}_{run}.parquet"
    if gene == "DDX3X":
        # not ddx3x_spliceai_d50 or ddx3x_nt: evid_external.DDX3X_SCORE_FILES reads those
        # names if they exist, and this stage changes what E7 reads only in commit 3
        return {"spliceai": "ddx3x_spliceai.parquet", "nt": "ddx3x_nt_a6000.parquet",
                "pangolin": "ddx3x_pangolin_fullprec.parquet"}.get(run, f"ddx3x_{run}_rescored.parquet")
    return f"tp53_{run}_rescored.parquet"


def _floats(s: pd.Series) -> np.ndarray:
    """Python's own float() of each repr: exact, unlike pandas's default parser."""
    return np.array([float(x) for x in s], dtype=float)


def _tsv(path: Path, columns: list[str]) -> pd.DataFrame:
    t = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    names = ["value"] + [c for c in t.columns if c.startswith("ds_")]
    out = pd.DataFrame({"variant_id": t["variant_id"], "scored_utc": t["scored_utc"]})
    for name, col in zip(names, columns):
        out[col] = _floats(t[name])
    return out


def read_run(run: str) -> tuple[pd.DataFrame, list[Path], dict]:
    """(variant_id, the run's columns, scored_utc) for every variant the run scored; the
    cache files read; the run's time range when it records no time per variant."""
    spec = RUNS[run]
    paths = sorted(p for pat in spec["cache"] for p in CACHE.glob(pat))
    span = {}
    if run in ("spliceai", "spliceai_walker", "pangolin"):
        df = pd.concat([_tsv(p, spec["columns"]) for p in paths], ignore_index=True)
    elif run in ("cadd", "conservation", "gpn_msa"):
        src = {"cadd": "cadd", "phylop": "phylop100way", "phastcons": "phastcons100way",
               "gpn_msa": "gpn_msa"}
        df = None
        for col, p in zip(spec["columns"], paths if run != "conservation" else
                          [CACHE / "lookup_phylop.parquet", CACHE / "lookup_phastcons.parquet"]):
            d = pd.read_parquet(p)[["variant_id", src[col]]].rename(columns={src[col]: col})
            df = d if df is None else df.merge(d, on="variant_id", how="outer")
        env = json.loads((CACHE / "lookups.environment.json").read_text())
        span = {"scored_utc_range": [env["run_started_utc"], env["run_ended_utc"]]}
    elif run == "avi":
        d = pd.read_parquet(paths[0])
        df = d[["variant_id"] + AVI].copy()
        seg = json.loads((CACHE / "avi_e18.segments.json").read_text())
        span = {"scored_utc_range": [seg[0]["start_utc"], seg[-1]["end_utc"]],
                "segments": seg}
    else:                                                       # nt
        t = pd.read_csv(paths[0], sep="\t", dtype=str, keep_default_na=False)
        df = pd.DataFrame({"variant_id": t["variant_id"],
                           "nt": [float(x) if x != "" else np.nan for x in t["nt"]]})
        env = json.loads((CACHE / "nt/out/environment.json").read_text())
        span = {"scored_utc_range": [env["start_utc"], env["end_utc"]]}
    return df, paths, span


def anchors_check(run: str, df: pd.DataFrame, anchors: pd.DataFrame, stored: pd.DataFrame) -> dict:
    cols = [c for c in RUNS[run]["columns"] if c in stored.columns]
    got = df.set_index("variant_id").reindex(anchors["variant_id"])
    want = stored.set_index("variant_id").loc[anchors["variant_id"], cols]
    out = {"n_anchors": int(len(anchors)), "columns": {}}
    for c in cols:
        g, w = got[c].to_numpy(dtype=float), want[c].to_numpy(dtype=float)
        same = (g == w) | (np.isnan(g) & np.isnan(w))
        out["columns"][c] = {"n_identical": int(same.sum()),
                             "differing": [v for v, s in zip(anchors["variant_id"], same) if not s]}
    return out


def run_logs(run: str) -> list[str]:
    spec = RUNS[run]
    if spec["logs"] is None:
        return []
    name, mark = spec["logs"]
    out = []
    for p in sorted(RUN_LOGS.glob(f"{name}_*.json")):
        argv = json.loads(p.read_text())["argv"]
        if "--trial" in argv or (mark is not None and mark not in argv):
            continue
        out.append(rel(p))
    return out


def environment(run: str) -> dict:
    paths = sorted(p for pat in RUNS[run]["environments"] for p in CACHE.glob(pat))
    envs = [json.loads(p.read_text()) for p in paths]
    keep = ("interpreter", "python", "platform", "pip_freeze", "scorer_sha256",
            "torch_threads", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
            "workers", "client_version", "atlas_module", "api_key_present",
            "api_key_source", "gpu_name", "gpu_capability", "torch", "cuda", "cudnn",
            "transformers", "nvidia_smi", "nvidia_smi_query", "pip_freeze_all",
            "checkpoint", "model_revision", "model_class", "dtype", "device")
    first = {k: envs[0][k] for k in keep if k in envs[0]}
    for e in envs[1:]:
        for k in ("interpreter", "python", "pip_freeze", "scorer_sha256"):
            if e.get(k) != envs[0].get(k):
                raise SystemExit(f"{run}: the environments of its processes differ in {k}")
    return first | {"processes": len(envs)}


def _write_parquet(df: pd.DataFrame, path: Path) -> None:
    df.to_parquet(path, index=False)


def _write_json(obj, path: Path) -> None:
    path.write_text(json.dumps(obj, indent=1, ensure_ascii=False) + "\n")


def _manifest(path: Path) -> dict[str, str]:
    """name -> sha256, from a MANIFEST.sha256 the pod wrote."""
    return {n: h for h, n in (ln.split("  ", 1) for ln in path.read_text().splitlines() if ln)}


def nt_evidence(df: pd.DataFrame, anchors: pd.DataFrame, stored: pd.DataFrame, log) -> dict:
    """NT's anchors do not hold: its run's evidence, and no NT file per gene."""
    src = CACHE / "nt"
    for sub in ("out", "out_trial"):
        for name, h in _manifest(src / sub / "MANIFEST.sha256").items():
            if sha(src / sub / name) != h:
                raise SystemExit(f"[nt] {sub}/{name} is not the file the pod's manifest names")
    committed = {}
    for r in NT_SERVED:
        dst = NT_DIR / r
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_bytes((src / r).read_bytes())
        committed[r] = sha(dst)
        log.output(dst)
    got = df.set_index("variant_id")["nt"]
    st = stored.set_index("variant_id")
    gene = anchors.set_index("variant_id")["gene"]
    rows = []
    for v in anchors["variant_id"]:
        g, w = float(got[v]), float(st.loc[v, "nt"])
        rows.append({"variant_id": v, "gene": gene[v], "stratum": st.loc[v, "stratum"],
                     "stored": repr(w), "rescored": repr(g), "difference": repr(g - w),
                     "identical": g == w})
    table = pd.DataFrame(rows)
    out = NT_DIR / "anchor_check.tsv"
    table.to_csv(out, sep="\t", index=False)
    log.output(out)
    diff = np.array([abs(float(x)) for x in table["difference"]])
    rel_diff = diff / np.abs(np.array([float(x) for x in table["stored"]]))

    def results(sub):
        t = pd.read_csv(src / sub / "results.tsv", sep="\t", dtype=str, keep_default_na=False)
        return t.set_index("variant_id")
    full, trial = results("out"), results("out_trial")
    both = full.loc[trial.index, "nt"] == trial["nt"]          # repr against repr: bit for bit
    in_trial = trial["set"] == "anchor"
    env = json.loads((src / "out/environment.json").read_text())
    query = [ln for ln in env["nvidia_smi_query"].splitlines() if ln and not ln.startswith("name")]
    cuda_max = re.search(r"CUDA Version:\s*([0-9.]+)", env["nvidia_smi"])
    k = int(table["identical"].sum())
    rec = {
        "column": "nt",
        "rule": NT_RULE,
        "result": (f"{k}/{len(table)} anchors identical bit for bit, largest difference "
                   f"{diff.max():.3g}: NT does not enter the external test"),
        "anchors": {"file": rel(out), "n": len(table), "n_identical": k,
                    "largest_difference": repr(float(diff.max())),
                    "largest_relative_difference": repr(float(rel_diff.max()))},
        "trial_and_full_run_agree": {
            "anchors": f"{int(both[in_trial].sum())}/{int(in_trial.sum())}",
            "external": f"{int(both[~in_trial].sum())}/{int((~in_trial).sum())}"},
        "run": {k_: env[k_] for k_ in ("gpu_name", "gpu_capability", "torch", "cuda", "cudnn",
                                        "transformers", "checkpoint", "model_revision",
                                        "model_class", "dtype", "device", "start_utc", "end_utc",
                                        "n_variants", "n_scored", "package_manifest_sha256")}
               | {"nvidia_smi_query": query, "driver_cuda_version": cuda_max and cuda_max.group(1)},
        "committed_as_served": committed,
        "not_committed": {
            "files": {f"{sub}/results.tsv": _manifest(src / sub / "MANIFEST.sha256")["results.tsv"]
                      for sub in ("out", "out_trial")},
            "kept_in": "phase1/data/e18_cache/nt/ (gitignored)",
            "why": ("the 6,667 external values are not the analysis set's variable and are not "
                    "used (the author's decision of 2026-10-23)")},
        "package_run_logs": [rel(q) for q in sorted(RUN_LOGS.glob("nt_gpu_package_*.json"))],
        "scored_utc_range": [env["start_utc"], env["end_utc"]],
    }
    _write_json(rec, NT_DIR / "record.json")
    log.output(NT_DIR / "record.json")
    print(f"[nt] {rec['result']}; evidence in {rel(NT_DIR)}")
    return {"record": rel(NT_DIR / "record.json")}


def _nt_cell() -> str | None:
    """The NT cell's text when NT's anchors did not hold, else None."""
    rec = NT_DIR / "record.json"
    if not rec.exists():
        return None
    r = json.loads(rec.read_text())["anchors"]
    return (f"no column of the same variable: {r['n_identical']}/{r['n']} anchors identical bit "
            f"for bit, largest difference {float(r['largest_difference']):.3g}; by the rule "
            "fixed before the run (the author's decision of 2026-10-23) NT does not enter the "
            "external test")


def build_run(run: str, inp: pd.DataFrame, stored: pd.DataFrame, log) -> dict | None:
    spec = RUNS[run]
    paths = [p for pat in spec["cache"] for p in CACHE.glob(pat)]
    if not paths:
        if spec.get("optional"):
            print(f"[{run}] no results in the cache yet: left out")
            return None
        raise SystemExit(f"[{run}] no cache: run its scorer first")
    df, paths, span = read_run(run)
    for p in paths:
        log.input(p)
    if df["variant_id"].duplicated().any():
        raise SystemExit(f"[{run}] a variant was scored twice")
    anchors = inp[inp["set"] == "anchor"]
    expect = set(anchors["variant_id"]) | set(inp.loc[(inp["set"] == "external")
                                                      & inp["gene"].isin(spec["genes"]), "variant_id"])
    got = set(df["variant_id"])
    if run in ("cadd", "conservation", "gpn_msa"):       # the atlas's scorers drop the unscored
        if not got <= expect or not set(anchors["variant_id"]) <= got:
            raise SystemExit(f"[{run}] scored variants it was not asked for, or missed an anchor")
    elif got != expect:
        raise SystemExit(f"[{run}] scored {len(got)} variants, not the {len(expect)} it was asked for")
    check = anchors_check(run, df, anchors, stored)
    bad = {c: v for c, v in check["columns"].items() if v["n_identical"] != check["n_anchors"]}
    if bad and run == "nt":
        return {"nt_evidence": nt_evidence(df, anchors, stored, log), "anchors": check}
    if bad:
        raise SystemExit(f"[{run}] anchors not identical (E18.0 rule 4: stop): {bad}")
    env = environment(run)
    logs = run_logs(run)
    cache = {rel(p): sha(p) for p in paths}
    extra = {}
    if run in ("cadd", "conservation", "gpn_msa"):
        s = json.loads((CACHE / "lookups.sources.json").read_text())
        key = {"cadd": "cadd", "conservation": "ucsc", "gpn_msa": "gpn_msa"}[run]
        extra["source"] = {"read_utc": s["read_utc"], key: s[key]}
    cols = spec["columns"]
    per_row = "scored_utc" in df.columns
    files = {}
    groups = [("anchors", anchors)] + [
        (g, inp[(inp["set"] == "external") & (inp["gene"] == g)]) for g in spec["genes"]]
    if run == "avi":
        groups = [("anchors", anchors), ("TP53", groups[1][1]),
                  ("NEW", inp[(inp["set"] == "external") & inp["gene"].isin(NEW)])]
    for g, rows in groups:
        name = (f"e18_anchors_{run}.parquet" if g == "anchors" else
                "newgenes_avi.parquet" if g == "NEW" else file_name(run, g))
        out = EXT / name
        part = rows.merge(df, on="variant_id", how="left")
        keep = (["variant_id", "gene", "chrom", "pos", "ref", "alt"] + cols
                + (["scored_utc"] if per_row else []))
        if g == "anchors":
            st = stored.set_index("variant_id")
            part = part.assign(stratum=st.loc[part["variant_id"], "stratum"].to_numpy())
            keep = keep[:2] + ["stratum"] + keep[2:]
        part = part[keep].reset_index(drop=True)
        _write_parquet(part, out)
        times = (sorted(part["scored_utc"].dropna()) if per_row else span["scored_utc_range"])
        prov = {
            "file": rel(out),
            "genes": sorted(set(part["gene"])), "set": "anchors" if g == "anchors" else "external",
            "columns": cols, "run": run, "tool": spec["tool"], "scorer": spec["scorer"],
            "configuration": spec["configuration"], "precision": spec["precision"],
            "direction": spec["direction"],
            "same_variable": {
                "rule": "E18.0 rule 4(a): the run's anchors equal the analysis set's stored "
                        "values bit for bit",
                "anchors": rel(EXT / f"e18_anchors_{run}.parquet"),
                "n_anchors": check["n_anchors"],
                "n_identical": {c: v["n_identical"] for c, v in check["columns"].items()}},
            **extra,
            "environment": rel(Path(str(out) + ".environment.json")),
            "run_logs": logs, "cache": cache,
            "n": int(len(part)),
            "n_scored": {c: int(part[c].notna().sum()) for c in cols},
            "scored_utc_range": [times[0], times[-1]],
            "sha256": sha(out),
        }
        if "segments" in span:
            prov["segments"] = span["segments"]
        _write_json(prov, Path(str(out) + ".provenance.json"))
        _write_json({"file": rel(out)} | env, Path(str(out) + ".environment.json"))
        for p in (out, Path(str(out) + ".provenance.json"), Path(str(out) + ".environment.json")):
            log.output(p)
        files[g] = rel(out)
        print(f"[{run}] {rel(out)}: {len(part)} rows; "
              + ", ".join(f"{c} {prov['n_scored'][c]}" for c in cols))
    return {"files": files, "anchors": check}


# ---------------------------------------------------------------------------
# assembly (steps 3 and 4): from the committed files only
# ---------------------------------------------------------------------------
def _alphagenome() -> pd.DataFrame:
    """The official definition's composite as the analysis set reads it
    (evid_build_set.attach_alphagenome: the round-trip parser)."""
    ag = pd.read_csv(AG_EXTERNAL, sep="\t", usecols=["variant_id", "alphagenome_splicing"],
                     float_precision="round_trip")
    return ag.rename(columns={"alphagenome_splicing": "alphagenome"})


def _committed(run: str, gene: str) -> Path:
    return EXT / ("newgenes_avi.parquet" if run == "avi" and gene in NEW else file_name(run, gene))


def assemble() -> tuple[pd.DataFrame, dict]:
    parts, inputs = [], {}
    for g in NEW:
        src = EXT / f"{g.lower()}_variants.parquet"
        v = pd.read_parquet(src)
        inputs[rel(src)] = sha(src)
        parts.append(v[v["offset"].abs().between(3, 50)])
    df = pd.concat(parts, ignore_index=True)
    ag = _alphagenome()
    inputs[rel(AG_EXTERNAL)] = sha(AG_EXTERNAL)
    columns, ranges, absent = {}, [], {}
    for col in COLUMNS:
        if col == "alphagenome":
            df = df.merge(ag, on="variant_id", how="left")
            prov = json.loads(Path(str(AG_EXTERNAL) + ".provenance.json").read_text())
            ranges.append(prov["scored_utc_range"])
            columns[col] = {"source": rel(AG_EXTERNAL)}
            continue
        run = COLUMN_RUN[col]
        files = sorted({_committed(run, g) for g in NEW})
        if not all(f.exists() for f in files):
            if RUNS[run].get("optional"):
                absent[col] = (f"{_nt_cell()} ({rel(NT_DIR / 'anchor_check.tsv')})" if _nt_cell()
                               else "not part of this commit: NT is scored on the original GPU "
                                    "and has a commit of its own (the author's decision of "
                                    "2026-10-23)")
                continue
            raise SystemExit(f"{col}: {[rel(f) for f in files if not f.exists()]} missing")
        take = [col] + ([c for c in WALKER if c != col] if col == "spliceai_walker" else [])
        vals = pd.concat([pd.read_parquet(f, columns=["variant_id"] + take) for f in files],
                         ignore_index=True)
        df = df.merge(vals, on="variant_id", how="left")
        for f in files:
            inputs[rel(f)] = sha(f)
            ranges.append(json.loads(Path(str(f) + ".provenance.json").read_text())["scored_utc_range"])
        columns[col] = {"source": [rel(f) for f in files]}
    for col, info in columns.items():
        info["n_scored_by_gene"] = {g: int(df.loc[df["gene"] == g, col].notna().sum()) for g in NEW}
    manifest = {
        "file": rel(SCORED),
        "n": int(len(df)), "n_by_gene": {g: int((df["gene"] == g).sum()) for g in NEW},
        "rows": "each gene's <gene>_variants.parquet rows at 3 <= |offset| <= 50 (E18.1), "
                "in that file's order, with every one of its columns unchanged: the labels "
                "are E16's and, for SBDS, E16b's",
        "columns": columns,
        "not_scored": absent,
        "same_variable": rel(DOC),
        "scored_utc_range": [min(r[0] for r in ranges), max(r[1] for r in ranges)],
        "inputs": inputs,
    }
    return df, manifest


def _provenance_section2() -> dict:
    """docs/column-provenance.md section 2, the committed value-for-value record."""
    text = COLUMN_PROVENANCE.read_text()
    sec = text.split("## 2.", 1)[1].split("## 3.", 1)[0]
    rows = [r for r in sec.splitlines() if r.startswith("| ") and not r.startswith("| column")
            and not r.startswith("|---")]
    head = [h.strip() for h in next(r for r in sec.splitlines() if r.startswith("| column")).strip("|").split("|")]
    out = {}
    for r in rows:
        cells = [c.strip() for c in r.strip("|").split("|")]
        out[cells[0]] = dict(zip(head, cells))
    return out


def _avi_same_run() -> str:
    a = pd.read_parquet(AVI_RUN)
    s = pd.read_parquet(SET, columns=["variant_id"] + AVI).set_index("variant_id")
    x = a[a["set"] == "analysis"].set_index("variant_id").loc[s.index, AVI]
    same = int(((x.to_numpy() == s.to_numpy()) | (np.isnan(x.to_numpy()) & np.isnan(s.to_numpy())))
               .all(axis=1).sum())
    n = {k: int(v) for k, v in a["set"].value_counts().items()}
    return (f"(a) the analysis set's own run: avi.parquet scored the analysis set's "
            f"{n['analysis']:,}, DDX3X's {n['ddx3x']:,} and TP53's first {n['tp53']} in one run; "
            f"the analysis set's values are that run's ({same:,}/{len(s):,} identical)")


def identity_table(built: dict) -> str:
    sec2 = _provenance_section2()
    drift = {}
    for k, p in AG_DRIFT.items():
        d = json.loads(Path(str(p) + ".provenance.json").read_text())["drift_check"]
        drift[k] = f"{d['identical']}/{d['variants']}"
    avi_run = _avi_same_run()
    rows = []
    for g in GENES:
        for col in COLUMNS:
            if col == "alphagenome":
                rows.append((g, col, "scripts/stage_E17/score_official.py (the official "
                             "definition, notebook 55943112)", "client 0.10.0, default model",
                             "float32, read back round-trip", "as the analysis set's column",
                             f"(a) drift checks, the pilot's 10 variants before and after the "
                             f"external run: {drift['start']} and {drift['end']} identical",
                             rel(AG_EXTERNAL)))
                continue
            run = COLUMN_RUN[col]
            spec = RUNS[run]
            if g in spec["genes"]:
                f = _committed(run, g)
                if not f.exists() and run == "nt" and _nt_cell():
                    rows.append((g, col, spec["scorer"], spec["tool"], spec["precision"],
                                 spec["direction"], _nt_cell(), rel(NT_DIR / "anchor_check.tsv")))
                    continue
                if not f.exists():
                    rows.append((g, col, spec["scorer"], spec["tool"], spec["precision"],
                                 spec["direction"], "not part of this commit: scored on the "
                                 "original GPU, in a commit of its own (the author's decision "
                                 "of 2026-10-23)", "—"))
                    continue
                prov = json.loads(Path(str(f) + ".provenance.json").read_text())
                k = prov["same_variable"]["n_identical"][col]
                n = prov["same_variable"]["n_anchors"]
                rows.append((g, col, spec["scorer"], spec["tool"], spec["precision"],
                             spec["direction"], f"(a) E18 re-scored, anchors {k}/{n} identical "
                             f"({prov['same_variable']['anchors']})", rel(f)))
            elif g == "DDX3X" and col in AVI:
                rows.append((g, col, RUNS["avi"]["scorer"], RUNS["avi"]["tool"],
                             RUNS["avi"]["precision"], RUNS["avi"]["direction"], avi_run,
                             "phase1/data/evidence/avi.parquet"))
            elif g == "TP53" and col in ("cadd", "phylop", "phastcons", "gpn_msa"):
                r = sec2[col]
                frac = (f"identical_up_to_sign_frac {r['identical_up_to_sign_frac']}" if col == "gpn_msa"
                        else f"identical_frac {r['identical_frac']}")
                how = ("the published study's scorer, the same tabix source at the pinned "
                       "revision cf1718a9; raw sign, turned to the analysis set's direction "
                       "when the external test reads it (E18.0 rule 5)" if col == "gpn_msa" else
                       "the published study's scorer (scripts/phase4a_score_tp53.py calls "
                       + ("75_cadd.fetch)" if col == "cadd" else "73_conservation.fetch_track)"))
                rows.append((g, col, how, "as the analysis set's column", "as served",
                             "raw sign (negative = more deleterious)" if col == "gpn_msa"
                             else "as the analysis set's column",
                             f"(b) docs/column-provenance.md section 2: {frac} over "
                             f"{int(r['n_shared_scored']):,} shared variants",
                             "phase1/data/evidence/external/tp53_splice_scored_12nt.parquet"))
            else:
                raise SystemExit(f"{g} x {col}: no evidence and not scored")
    head = ("| gene | column | scorer | version | precision | direction | evidence and result | file |\n"
            "|---|---|---|---|---|---|---|---|\n")
    body = "".join("| " + " | ".join(str(c).replace("|", "\\|") for c in r) + " |\n" for r in rows)
    n_a = sum(1 for r in rows if r[6].startswith("(a)"))
    n_b = sum(1 for r in rows if r[6].startswith("(b)"))
    n_p = sum(1 for r in rows if r[6].startswith("not part of this commit"))
    n_none = sum(1 for r in rows if r[6].startswith("no column of the same variable"))
    return (
        "# External genes: the same variable as the analysis set, column by column\n\n"
        "Generated by `scripts/stage_E18/build_columns.py` (E18.2); do not edit by hand.\n\n"
        "Every predictor column of the seven external genes (DDX3X, TP53 and, at "
        "3 <= |offset| <= 50, CTCF, SFPQ, RAD51D, XRCC2 and SBDS) must come from the scorer "
        "that produced the analysis set's column: the same version, precision and direction "
        "(E18.0 rule 4). The evidence is one of two kinds: (a) anchors, analysis-set variants "
        "re-scored in the same run, equal to the analysis set's stored values bit for bit "
        "(three per training gene and stratum s3_10, s11_50, drawn with config.RANDOM_SEED by "
        "`scripts/stage_E18/anchors.py`, 39 in all; for AlphaGenome, the pilot's ten variants "
        "of the drift checks); or (b) a committed value-for-value record. No tolerance.\n\n"
        f"{len(rows)} cells: {n_a} with evidence (a), {n_b} with evidence (b)"
        + (f", {n_p} not part of this commit" if n_p else "")
        + (f", {n_none} with no column of the same variable (NT)" if n_none else "")
        + ".\n\n" + head + body)


def write_assembled(df: pd.DataFrame, manifest: dict, built: dict) -> dict[Path, bytes]:
    buf = io.BytesIO()
    df.to_parquet(buf, index=False)
    data = buf.getvalue()
    manifest = dict(manifest, sha256=hashlib.sha256(data).hexdigest())
    return {SCORED: data,
            MANIFEST: (json.dumps(manifest, indent=1, ensure_ascii=False) + "\n").encode(),
            DOC: identity_table(built).encode()}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--offline", action="store_true",
                    help="steps 3 and 4 from the committed files only; fail unless byte-identical")
    ap.add_argument("--runs", default=",".join(RUNS),
                    help="the runs to build from the cache (default: all)")
    args = ap.parse_args()
    with RunLog(__file__) as log:
        built = {}
        if not args.offline:
            inp = pd.read_parquet(INPUT)
            log.input(INPUT)
            # every column a run writes that the analysis set has: the thirteen, and the
            # walker's four components beside its maximum
            have = set(pd.read_parquet(SET).columns)
            stored = pd.read_parquet(SET, columns=["variant_id", "stratum"] + [
                c for r in RUNS.values() for c in r["columns"] if c in have])
            log.input(SET)
            for run in args.runs.split(","):
                built[run] = build_run(run, inp, stored, log)
            log.note("anchors", {r: b["anchors"] for r, b in built.items() if b})
        df, manifest = assemble()
        for p in manifest["inputs"]:
            log.input(REPO / p)
        out = write_assembled(df, manifest, built)
        if args.offline:
            same = {rel(p): p.exists() and p.read_bytes() == data for p, data in out.items()}
            log.note("byte_identical", same)
            print("\n".join(f"  {'identical' if s else 'DIFFERS'}  {p}" for p, s in same.items()))
            if not all(same.values()):
                raise SystemExit("--offline: not byte for byte the committed files")
            print("--offline: every assembled file is byte for byte the committed one")
            return
        for p, data in out.items():
            p.write_bytes(data)
            log.output(p)
        print(f"{rel(SCORED)}: {manifest['n']} rows; {rel(DOC)} written")


if __name__ == "__main__":
    main()
