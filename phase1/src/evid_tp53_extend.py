"""E7 support -- TP53 extended from |offset| <= 8 to |offset| <= 12.

The TP53 table the external-gene stages read (data/external/tp53_splice_scored_v2
.parquet) holds 192 intron-side SNVs, 24 at each |offset| 1 to 8. The window was a
constant of the scoring script (scripts/phase4a_score_tp53.py, SPLICE_REGION_MAX
= 8), not a property of the data: the deposit (MaveDB urn:mavedb:00001213-a-1, a
copy at data/external/tp53_mavedb_scores.tsv) carries intron-side SNVs to |offset|
12, 24 at each, and the 96 at 9 to 12 were never scored. Those 96 are what carries
TP53 across the 10/11 bp boundary into the 11-50 bp stratum, which the frozen 192
cannot populate.

This stage builds the 288-row table. It changes nothing that exists:

  * the 192 frozen rows are copied from the v2 table and asserted identical to it,
    value for value;
  * the 96 new rows are parsed and mapped by the scoring script's own functions,
    imported by path with the window widened at run time, and every derived field
    (orientation, offset side, region, missingness flags) comes from the same
    transforms. The 192 are rebuilt the same way and must come out identical to the
    frozen rows, which is the check that the rule applied to the 96 is the rule
    that built the 192;
  * every score column of the new rows comes from the scorer that produced that
    column for the 192, called through its own entry point (scripts/phase4b_score_
    pinned.py for SpliceAI at distance 50, Pangolin and the Nucleotide Transformer;
    scripts/phase4a_score_tp53.py for CADD, phyloP, phastCons, GPN-MSA and gnomAD;
    the AlphaGenome client 0.6.1 scoring path for `alphagenome`; this project's
    Walker-basis scorer for `spliceai_walker`; the Atlas lookup for the four Atlas
    columns). Each model runs in its own pinned environment, as the 192 were.

A score is only as good as the claim that it is the same variable, so each scorer
also re-scores frozen rows beside the new ones -- 24 of the 192 (three at each
offset) for the per-variant scorers, all 192 where one query covers them anyway --
and the build stops unless every re-scored value equals the frozen one exactly.
Exactly means as stored: the frozen SpliceAI, Pangolin and NT values reached the v2
table through a TSV written by pandas and read back with its default float parser,
which moves the last bits of a double, so the new rows' values take the same path
(TSV_STORED below). Their float32 model outputs agree before that step as well.

AlphaMissense stays empty for the new rows, as it is for the 192: it scores
missense variants only, so an intronic SNV has no entry, and the file the scoring
script reads is not on this machine. The manifest says so rather than leaving a
silent gap.

For the in-frame attribution (evid_inframe, E2.7) the false-positive subset is
recomputed over the 288 by the same rule (Walker-basis score at or above the PP3
cut point, control-anchored label normal) and SpliceAI's event positions are scored
for it at -D 4999. Both go to new files; the 192-row subset and event table are
left as they are, and their rows must be reproduced exactly inside the new ones.

Each scorer writes a parquet, a provenance json and an environment json (the
interpreter and its full `pip freeze`) beside it. A scorer whose output already
exists is skipped, so a rerun needs no network, no GPU and no model environment,
and produces a byte-identical table.

Outputs (data/evidence/external/):
    tp53_splice_scored_12nt.parquet        the 288-row table, rows sorted by variant_id
    tp53_splice_scored_12nt.manifest.json  counts, coverage, definitions, checks, hashes
    tp53_12nt_*.parquet (+ .provenance.json, .environment.json) per-scorer outputs
    tp53_avi_12nt.parquet                  the four Atlas columns for all 288
    tp53_inframe_subset_12nt.parquet       E2.7 false-positive subset over the 288
    tp53_inframe_events_12nt.parquet       SpliceAI events for that subset

Model environments (override with the environment variable in brackets):
    SpliceAI     <atlas>/models/spliceai/.venv            [EVID_SPLICEAI_PYTHON]
    Pangolin     <atlas>/models/pangolin/.venv            [EVID_PANGOLIN_PYTHON]
    NT           <atlas>/../envs/nt                        [EVID_NT_PYTHON]
    AlphaGenome  an env with alphagenome==0.6.1 (ag-env)   [EVID_AG061_PYTHON]
    Atlas        an env with alphagenome>=0.9.0            [EVID_AG_ATLAS_PYTHON]
The AlphaGenome key is read at run time from $ALPHAGENOME_API_KEY or
~/.config/alphagenome/key (evid_score_avi.read_key) and is never written anywhere.

Run (from phase1/):  python -m src.evid_tp53_extend
                     python -m src.evid_tp53_extend --offline   # cached scores only
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import io
import json
import os
import platform
import subprocess
import sys
import time
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

PHASE1 = Path(__file__).resolve().parents[1]
REPO = PHASE1.parent
SCRIPTS = REPO / "scripts"
EXT_DIR = PHASE1 / "data" / "evidence" / "external"

# frozen inputs, read only
DEPOSIT = REPO / "data" / "external" / "tp53_mavedb_scores.tsv"
FROZEN_TABLE = PHASE1 / "data" / "external" / "tp53_splice_scored_v2.parquet"
FROZEN_MANIFEST = PHASE1 / "data" / "external" / "tp53_splice_scored_v2.manifest.json"
FROZEN_WALKER = EXT_DIR / "tp53_spliceai_walker.parquet"
FROZEN_AVI = PHASE1 / "data" / "evidence" / "avi.parquet"
FROZEN_SUBSET = EXT_DIR / "tp53_inframe_subset.parquet"
FROZEN_EVENTS = EXT_DIR / "tp53_inframe_events.parquet"
WALKER_CFG = PHASE1 / "config" / "walker2023.yaml"
SCORE_TP53 = SCRIPTS / "phase4a_score_tp53.py"
SCORE_PINNED = SCRIPTS / "phase4b_score_pinned.py"
AG_V061 = SCRIPTS / "phase4a_ag_drift_check.py"
WALKER_SCORER = PHASE1 / "src" / "evid_score_spliceai_walker.py"
EVENT_SCORER = PHASE1 / "src" / "evid_score_spliceai_events.py"

# outputs
OUT = EXT_DIR / "tp53_splice_scored_12nt.parquet"
MANIFEST = EXT_DIR / "tp53_splice_scored_12nt.manifest.json"
ALL_VARIANTS = EXT_DIR / "tp53_12nt_variants.parquet"
SCORE_INPUT = EXT_DIR / "tp53_12nt_score_input.parquet"
AVI_OUT = EXT_DIR / "tp53_avi_12nt.parquet"
SUBSET_OUT = EXT_DIR / "tp53_inframe_subset_12nt.parquet"
EVENTS_OUT = EXT_DIR / "tp53_inframe_events_12nt.parquet"
WORK = EXT_DIR / "tp53_12nt_work"

OFFSET_MAX = 12
PER_OFFSET = 24
CHECK_PICKS = (0, 8, 16)       # per |offset| 1..8, by variant_id: three sites, three alts
BASE = ["variant_id", "gene", "chrom", "pos", "ref", "alt"]
WALKER_COLS = ["spliceai_walker", "spliceai_walker_ds_ag", "spliceai_walker_ds_al",
               "spliceai_walker_ds_dg", "spliceai_walker_ds_dl"]
AVI_COLS = ["avi", "avi_splice_sites", "avi_splice_site_usage", "avi_splice_junctions"]

# The scorers. `input` is the variant list the scorer is run on: "check" = the 96 new
# rows plus 24 frozen rows re-scored as a check, "all" = all 288 (for scorers whose
# cost does not depend on the number of variants). `columns` maps the scorer's output
# column to the raw predictor column of scripts/phase4a_score_tp53.py (and so, via
# config.COLUMNS, to the table column).
SCORERS = {
    "spliceai": dict(file="tp53_12nt_spliceai_d50.parquet", env="spliceai", input="check",
                     columns={"spliceai_ds_fullprec": "spliceai_ds"}),
    "pangolin": dict(file="tp53_12nt_pangolin.parquet", env="pangolin", input="check",
                     columns={"pangolin_fullprec": "pangolin_score"}),
    "nt": dict(file="tp53_12nt_nt.parquet", env="nt", input="check",
               columns={"nucleotide_transformer": "nucleotide_transformer"}),
    "alphagenome": dict(file="tp53_12nt_alphagenome_v061.parquet", env="alphagenome_v061",
                        input="check", columns={"alphagenome_splice": "alphagenome_splice"}),
    "cadd": dict(file="tp53_12nt_cadd.parquet", env="analysis", input="check",
                 columns={"cadd_phred": "cadd_phred"}),
    "conservation": dict(file="tp53_12nt_conservation.parquet", env="analysis", input="all",
                         columns={"phylop100way": "phylop100way",
                                  "phastcons100way": "phastcons100way"}),
    "gpn_msa": dict(file="tp53_12nt_gpn_msa.parquet", env="analysis", input="all",
                    columns={"gpn_msa_score": "gpn_msa_score"}),
    "gnomad": dict(file="tp53_12nt_gnomad.parquet", env="analysis", input="all",
                   columns={"gnomad_af_global": "gnomad_af_global"}),
    "spliceai_walker": dict(file="tp53_12nt_spliceai_walker.parquet", env="spliceai",
                            input="check", columns={c: c for c in WALKER_COLS}),
    "avi": dict(file=AVI_OUT.name, env="atlas", input="all",
                columns={c: c for c in AVI_COLS}),
}

# What each column of the table is, for the manifest. The frozen rows' definitions are
# the ones recorded in the v2 build (phase1_build_frozen_matrix_v2.py) and the scoring
# script's BENCHMARK_SCORERS; the new rows are scored by the same entry points.
DEFINITIONS = {
    "spliceai": ("SpliceAI 1.3.1, -A grch38, -D 50, mask 0 (raw), max(DS_AG, DS_AL, DS_DG, "
                 "DS_DL) over the annotated genes, full float precision; "
                 "scripts/phase4b_score_pinned.py::score_spliceai -> "
                 "<atlas>/models/spliceai/score_fullprec.py"),
    "pangolin": ("Pangolin git 5cf94b8, -d 50, default masking, no score cutoff, "
                 "max(gain, |loss|) over transcript records, full float precision; "
                 "scripts/phase4b_score_pinned.py::score_pangolin -> "
                 "<atlas>/models/pangolin/score_fullprec.py"),
    "alphagenome": ("AlphaGenome dna_client, client 0.6.1, 1 Mb interval, the recommended "
                    "SPLICE_* variant scorers, max |raw_score| printed to 6 d.p.; "
                    "scripts/phase4a_ag_drift_check.py::build_scorer/score_one, the "
                    "verbatim copy of scripts/91_score_alphagenome.py::score_one that "
                    "scored the frozen rows (phase4a_tp53_external.py). A different "
                    "variable from the analysis set's client 0.7.0 column "
                    "(evid_external.DECLARED_MISMATCH)"),
    "gpn_msa": ("GPN-MSA precomputed hg38 scores (songlab/gpn-msa-hg38-scores at pinned "
                "revision cf1718a9), raw sign convention (config.REVERSED_FEATURES); "
                "scripts/phase4a_score_tp53.py::_score_gpn"),
    "nt": ("InstaDeepAI/nucleotide-transformer-v2-500m-multi-species, masked 6-mer LLR, "
           "6,000-bp window, CPU; scripts/phase4b_score_pinned.py::score_nt -> "
           "<atlas>/models/nt/score.py"),
    "cadd": ("CADD GRCh38-v1.7 PHRED via the CADD REST API; "
             "scripts/phase4a_score_tp53.py::_score_cadd -> scripts/75_cadd.py::fetch"),
    "alphamissense": "not scored for the new rows; see not_scored",
    "phylop": ("UCSC phyloP100way (hg38), getData/track API; "
               "scripts/phase4a_score_tp53.py::_score_conservation -> "
               "scripts/73_conservation.py::fetch_track"),
    "phastcons": ("UCSC phastCons100way (hg38), getData/track API; "
                  "scripts/phase4a_score_tp53.py::_score_conservation -> "
                  "scripts/73_conservation.py::fetch_track"),
    "gnomad_af": ("gnomAD r4 (gnomad_r4) global AF, max over genome and exome, GraphQL API; "
                  "scripts/phase4a_score_tp53.py::_score_gnomad -> "
                  "scripts/72_gnomad.py::fetch_gene('TP53')"),
    "spliceai_walker": ("SpliceAI 1.3.1, -A grch38, -D 4999, mask 0, max of the four deltas, "
                        "full float precision; the four components at the winning gene "
                        "beside it; src/evid_score_spliceai_walker.py"),
    "avi": ("AlphaGenome Atlas, client 0.9.0, scorer AVI_SCORE, max over tracks and rows; "
            "src/evid_score_avi.py::run"),
    "avi_splice_sites": "AlphaGenome Atlas, client 0.9.0, scorer SPLICE_SITES; as avi",
    "avi_splice_site_usage": "AlphaGenome Atlas, client 0.9.0, scorer SPLICE_SITE_USAGE; as avi",
    "avi_splice_junctions": "AlphaGenome Atlas, client 0.9.0, scorer SPLICE_JUNCTIONS; as avi",
}

# The frozen spliceai / pangolin / nt values reached the v2 table through a TSV:
# phase4b_score_pinned.py wrote <model>_tp53.tsv with DataFrame.to_csv and
# phase1_build_frozen_matrix_v2.build_tp53 read it back with pd.read_csv, whose default
# float parser is not round-trip exact -- it moves the last bits of a full-precision
# double (0.00040796271059662104 is stored as 0.0004079627105966). The model output,
# float32, is untouched. The new rows' values take the same path, so each column has
# one storage convention throughout; the score files keep the unrounded output.
TSV_STORED = {"spliceai", "pangolin", "nt"}

NOT_SCORED = {
    "alphamissense": (
        "AlphaMissense scores missense variants only, so an intron-side SNV has no entry "
        "by construction (the 192 frozen rows are 0/192 as well), and the file the "
        "scoring script reads (data/raw/scores/AlphaMissense_hg38.tsv.gz, "
        "phase4a_score_tp53._score_alphamissense) is not on this machine. Out of the "
        "E1 panel (evid_common.PANEL)."),
}

# Panel columns that must be complete on every row (the test reads this list).
REQUIRED_COMPLETE = ["spliceai", "spliceai_walker", "pangolin", "alphagenome", "gpn_msa",
                     "nt", "cadd", "phylop", "phastcons", *AVI_COLS]


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------
def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def _rel(path: Path) -> str:
    try:
        return str(Path(path).resolve().relative_to(REPO))
    except ValueError:
        return str(path)


def _load_script(path: Path, name: str):
    """Import a script by path. Numbered and phase scripts are not importable as
    modules; these three are side-effect free at import."""
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def _scoring_script():
    return _load_script(SCORE_TP53, "tp53_scoring_script")


@contextmanager
def _window(mod, max_offset: int):
    """The scoring script's own parser, with its window constant widened for the call."""
    old = mod.SPLICE_REGION_MAX
    mod.SPLICE_REGION_MAX = max_offset
    try:
        yield mod
    finally:
        mod.SPLICE_REGION_MAX = old


@contextmanager
def _cache_root(mod, root: Path):
    """_score_cadd keeps its request cache at ROOT/data/external/, which is frozen
    input. Point ROOT elsewhere for the call so the cache lands in this stage's work
    directory and nothing under data/external/ is written."""
    old = mod.ROOT
    (root / "data" / "external").mkdir(parents=True, exist_ok=True)
    mod.ROOT = root
    try:
        yield mod
    finally:
        mod.ROOT = old


@contextmanager
def _no_stray_dirs():
    """The numbered scorer scripts create data/raw/scores/ when imported. Remove it
    again if it did not exist and nothing was written into it."""
    d = REPO / "data" / "raw" / "scores"
    existed = d.exists()
    try:
        yield
    finally:
        if not existed and d.exists() and not any(d.iterdir()):
            d.rmdir()


def _atlas_repo() -> Path:
    from .evid_common import ATLAS_REPO
    return Path(ATLAS_REPO)


def _python(env: str, atlas: Path) -> str:
    """The interpreter of a scorer's pinned environment: the environment variable,
    or the environment beside the atlas checkout. Only a fresh scoring run needs
    these; the stage itself runs offline from the cached score files."""
    options = {
        "spliceai": ("EVID_SPLICEAI_PYTHON",
                     [atlas / "models/spliceai/.venv/bin/python"]),
        "pangolin": ("EVID_PANGOLIN_PYTHON",
                     [atlas / "models/pangolin/.venv/bin/python"]),
        "nt": ("EVID_NT_PYTHON", [atlas.parent / "envs/nt/bin/python"]),
        "alphagenome_v061": ("EVID_AG061_PYTHON",
                             [REPO / "ag-env/bin/python"]),
        "atlas": ("EVID_AG_ATLAS_PYTHON", [atlas.parent / "envs/alphagenome-atlas/bin/python"]),
    }
    if env == "analysis":
        return sys.executable
    var, candidates = options[env]
    if os.environ.get(var):
        return os.environ[var]
    for c in candidates:
        if Path(c).exists():
            return str(c)
    raise SystemExit(f"[tp53-12nt] no {env} environment found; set {var} "
                     f"(looked in {', '.join(str(c) for c in candidates)})")


ENV_SUFFIX = ".environment.json"
# the packages worth naming in the manifest for each environment; the sidecar keeps
# the full `pip freeze --all`
KEY_PACKAGES = {
    "spliceai": ["spliceai", "tensorflow", "keras", "numpy", "pandas", "pyfaidx", "pysam"],
    "pangolin": ["pangolin", "torch", "gffutils", "pyfaidx", "numpy", "pandas"],
    "nt": ["torch", "transformers", "pyfaidx", "numpy", "pandas"],
    "alphagenome_v061": ["alphagenome", "numpy", "pandas"],
    "atlas": ["alphagenome", "numpy", "pandas"],
    "analysis": ["pandas", "pyarrow", "numpy", "pysam", "requests"],
}


def record_environment(path: Path, python: str) -> None:
    """Which interpreter scored `path`, and every package version it had."""
    ver = subprocess.run([python, "-c", "import platform; print(platform.python_version())"],
                         capture_output=True, text=True, check=True).stdout.strip()
    freeze = subprocess.run([python, "-m", "pip", "freeze", "--all"],
                            capture_output=True, text=True).stdout.splitlines()
    rec = {"file": _rel(path), "interpreter": python, "python_version": ver,
           "pip_freeze": sorted(x for x in freeze if x.strip())}
    Path(str(path) + ENV_SUFFIX).write_text(json.dumps(rec, indent=2) + "\n")


def _environment_summary(path: Path, env: str) -> dict:
    side = Path(str(path) + ENV_SUFFIX)
    if not side.exists():
        return {"recorded": False}
    rec = json.loads(side.read_text())
    want = {p.lower() for p in KEY_PACKAGES[env]}
    key = [x for x in rec["pip_freeze"]
           if x.split("==")[0].split(" @ ")[0].strip().lower() in want]
    return {"interpreter": rec["interpreter"], "python_version": rec["python_version"],
            "key_packages": key, "full_record": _rel(side), "full_record_sha256": sha256_file(side)}


def _write_scored(out: pd.DataFrame, path: Path, *, tool: str, version, configuration,
                  scorer: str, input_path: Path, t0: float, record: dict | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    out.reset_index(drop=True).to_parquet(path, index=False)
    score_cols = [c for c in out.columns if c not in BASE]
    prov = {
        "tool": tool, "version": version, "configuration": configuration,
        "scorer": scorer,
        "n": int(len(out)),
        "n_scored": {c: int(out[c].notna().sum()) for c in score_cols},
        "input": _rel(input_path), "input_sha256": sha256_file(input_path),
        "sha256": sha256_file(path),
        "scored_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "runtime_seconds": round(time.time() - t0, 1),
        "python": platform.python_version(), "interpreter": sys.executable,
    }
    if record:
        prov["scorer_record"] = record
    Path(str(path) + ".provenance.json").write_text(
        json.dumps(prov, indent=2, default=str) + "\n")
    print(f"[tp53-12nt] wrote {_rel(path)} ({prov['n_scored']})", flush=True)


# ---------------------------------------------------------------------------
# the variant set
# ---------------------------------------------------------------------------
def deposit_rows(mod=None) -> pd.DataFrame:
    """Every intron-side SNV with a score at 1 <= |offset| <= 12, parsed and mapped by
    the scoring script's own load_tp53_splice() (HGVS(cDNA) for the offset, the
    deposit's HGVS(genomic) g.-coordinate on GRCh38 chr17 for the variant)."""
    mod = mod or _scoring_script()
    with _window(mod, OFFSET_MAX):
        raw = mod.load_tp53_splice()
    return raw.sort_values("variant_id").reset_index(drop=True)


def variant_lists(raw: pd.DataFrame, frozen: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    allv = raw[BASE + ["intron_offset"]].copy()
    allv["pos"] = allv["pos"].astype(int)
    allv["in_frozen_table"] = allv["variant_id"].isin(frozen["variant_id"])
    allv = allv.sort_values("variant_id").reset_index(drop=True)
    new = allv[~allv["in_frozen_table"]].assign(role="new")
    fz = allv[allv["in_frozen_table"]].copy()
    fz["abs_offset"] = fz["intron_offset"].abs().astype(int)
    picks = [g.sort_values("variant_id").iloc[list(CHECK_PICKS)]
             for _, g in fz.groupby("abs_offset")]
    check = pd.concat(picks).drop(columns="abs_offset").assign(role="reproduction_check")
    inp = pd.concat([new, check]).sort_values("variant_id").reset_index(drop=True)
    return allv, inp[BASE + ["intron_offset", "role"]]


def write_variant_lists() -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = deposit_rows()
    frozen = pd.read_parquet(FROZEN_TABLE)
    allv, inp = variant_lists(raw, frozen)
    EXT_DIR.mkdir(parents=True, exist_ok=True)
    allv.to_parquet(ALL_VARIANTS, index=False)
    inp.to_parquet(SCORE_INPUT, index=False)
    return allv, inp


def _check_reference_alleles(allv: pd.DataFrame, atlas: Path) -> dict:
    """The deposit's genomic REF base against GRCh38 at every one of the 288."""
    import pysam
    fasta = atlas / "data" / "refs" / "grch38_subset.fa"
    if not fasta.exists():
        return {"checked": False, "reason": f"{fasta} not found"}
    fa = pysam.FastaFile(str(fasta))
    bad = [r.variant_id for r in allv.itertuples()
           if fa.fetch(str(r.chrom), int(r.pos) - 1, int(r.pos)).upper() != r.ref]
    fa.close()
    if bad:
        raise SystemExit(f"[tp53-12nt] REF does not match GRCh38 at {bad[:5]} ...")
    return {"checked": True, "n": int(len(allv)), "mismatches": 0,
            "reference": "Ensembl r112 GRCh38 chr17 (<atlas>/data/refs/grch38_subset.fa)"}


# ---------------------------------------------------------------------------
# workers: each runs inside its scorer's environment
# ---------------------------------------------------------------------------
def _input(which: str) -> tuple[pd.DataFrame, Path]:
    path = SCORE_INPUT if which == "check" else ALL_VARIANTS
    df = pd.read_parquet(path)
    df["chrom"] = df["chrom"].astype(str)
    df["pos"] = df["pos"].astype(int)
    return df, path


def _pinned(name: str, atlas: Path):
    df, path = _input(SCORERS[name]["input"])
    p4b = _load_script(SCORE_PINNED, "tp53_pinned_scorers")
    prov: dict = {}
    t0 = time.time()
    frame = df[BASE].copy()
    if name == "spliceai":
        vals = p4b.score_spliceai(frame, atlas, prov)
        tool, version = "SpliceAI", prov.get("tool_version")
        entry = "scripts/phase4b_score_pinned.py::score_spliceai"
    elif name == "pangolin":
        vals = p4b.score_pangolin(frame, atlas, prov)
        tool, version = "Pangolin", f"{prov.get('tool_version')} (git {prov.get('git_commit')})"
        entry = "scripts/phase4b_score_pinned.py::score_pangolin"
    else:
        vals = p4b.score_nt(frame, atlas, prov, "cpu")
        tool, version = prov.get("tool"), prov.get("checkpoint_revision_sha")
        entry = "scripts/phase4b_score_pinned.py::score_nt"
    col = next(iter(SCORERS[name]["columns"]))
    out = df[BASE].copy()
    out[col] = np.asarray(vals, dtype=float)
    _write_scored(out, EXT_DIR / SCORERS[name]["file"], tool=tool, version=version,
                  configuration=prov.get("parameters"), scorer=entry, input_path=path,
                  t0=t0, record=prov)


def _alphagenome_v061(atlas: Path) -> None:
    import alphagenome
    ver = getattr(alphagenome, "__version__", "?")
    if ver != "0.6.1":
        raise SystemExit(f"[tp53-12nt] the frozen column is the client 0.6.1 definition; "
                         f"this environment has {ver}")
    from .evid_score_avi import read_key
    key = read_key()
    if not key:
        raise SystemExit("[tp53-12nt] no AlphaGenome key: set ALPHAGENOME_API_KEY or "
                         "write ~/.config/alphagenome/key")
    # build_scorer reads the key from this process's environment; it is not passed
    # on to any child process and not written anywhere
    os.environ["ALPHAGENOME_API_KEY"] = key
    del key
    ag = _load_script(AG_V061, "tp53_alphagenome_v061")
    df, path = _input("check")
    t0 = time.time()
    model, splice, dna_client, variant_scorers = ag.build_scorer()
    vals, failed = [], []
    for r in df.itertuples():
        s = ag.score_one((str(r.chrom), int(r.pos), r.ref, r.alt), model, splice,
                         dna_client, variant_scorers)
        if s is None:
            failed.append(r.variant_id)
        vals.append(float(s) if s not in (None, "") else np.nan)
    out = df[BASE].copy()
    out["alphagenome_splice"] = vals
    _write_scored(out, EXT_DIR / SCORERS["alphagenome"]["file"], tool="AlphaGenome",
                  version=f"alphagenome client {ver}",
                  configuration={"interval": "SEQUENCE_LENGTH_1MB",
                                 "variant_scorers": [str(s) for s in splice],
                                 "aggregate": "max |raw_score| over tidy_scores rows, "
                                              "formatted {:.6f}"},
                  scorer="scripts/phase4a_ag_drift_check.py::build_scorer + score_one",
                  input_path=path, t0=t0, record={"api_failures": failed})


def _analysis(name: str) -> None:
    """CADD, conservation, GPN-MSA and gnomAD: the scoring script's own functions."""
    mod = _scoring_script()
    df, path = _input(SCORERS[name]["input"])
    t0 = time.time()
    with _no_stray_dirs():
        if name == "cadd":
            with _cache_root(mod, WORK / "cadd"):
                cols = mod._score_cadd(df)
            meta = dict(tool="CADD", version="GRCh38-v1.7 (REST API)",
                        configuration={"api": "https://cadd.gs.washington.edu/api/v1.0/"
                                              "GRCh38-v1.7", "value": "PHRED, max over "
                                              "returned annotations"},
                        scorer="scripts/phase4a_score_tp53.py::_score_cadd")
        elif name == "conservation":
            cols = mod._score_conservation(df)
            meta = dict(tool="UCSC phyloP100way / phastCons100way", version="hg38",
                        configuration={"api": "https://api.genome.ucsc.edu/getData/track",
                                       "span": [int(df.pos.min()) - 1, int(df.pos.max())]},
                        scorer="scripts/phase4a_score_tp53.py::_score_conservation")
        elif name == "gpn_msa":
            cols = mod._score_gpn(df)
            meta = dict(tool="GPN-MSA", version=f"songlab/gpn-msa-hg38-scores@{mod.GPN_PIN}",
                        configuration={"index": _rel(Path(mod.GPN_TBI)),
                                       "sign": "raw (config.REVERSED_FEATURES)"},
                        scorer="scripts/phase4a_score_tp53.py::_score_gpn")
        elif name == "gnomad":
            cols = mod._score_gnomad(df)
            meta = dict(tool="gnomAD", version="gnomad_r4 (GraphQL API)",
                        configuration={"value": "max(genome af, exome af)",
                                       "query": "gene TP53, GRCh38"},
                        scorer="scripts/phase4a_score_tp53.py::_score_gnomad")
        else:
            raise SystemExit(f"unknown analysis scorer {name}")
    out = df[BASE].copy()
    for c, vals in cols.items():
        out[c] = pd.to_numeric(pd.Series(vals, index=out.index, dtype="object"),
                               errors="coerce").astype(float)
    _write_scored(out, EXT_DIR / SCORERS[name]["file"], input_path=path, t0=t0, **meta)


def _atlas_lookup() -> None:
    """The Atlas columns through evid_score_avi's own run(), pointed at the 288."""
    from . import evid_score_avi as A
    A.TARGETS = dict(A.TARGETS, tp53=ALL_VARIANTS)
    rc = A.run(None, ["tp53"], True, 250, AVI_OUT)
    if rc:
        raise SystemExit(rc)


def worker(name: str, atlas: Path) -> None:
    if name in ("spliceai", "pangolin", "nt"):
        _pinned(name, atlas)
    elif name == "alphagenome":
        _alphagenome_v061(atlas)
    elif name == "avi":
        _atlas_lookup()
    else:
        _analysis(name)


# ---------------------------------------------------------------------------
# orchestration: run what is missing, each in its own environment
# ---------------------------------------------------------------------------
def _child_env(extra: dict | None = None) -> dict:
    env = dict(os.environ)
    env.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")
    env.update(extra or {})
    return env


def _run(cmd: list[str], extra_env: dict | None = None) -> None:
    shown = " ".join(str(c) for c in cmd)
    print(f"[tp53-12nt] $ {shown}", flush=True)
    r = subprocess.run([str(c) for c in cmd], cwd=PHASE1, env=_child_env(extra_env))
    if r.returncode != 0:
        raise SystemExit(f"[tp53-12nt] failed ({r.returncode}): {shown}")


def ensure_scored(name: str, atlas: Path, offline: bool) -> Path:
    spec = SCORERS[name]
    path = EXT_DIR / spec["file"]
    env = spec["env"]
    if path.exists() and Path(str(path) + ".provenance.json").exists():
        if not offline and not Path(str(path) + ENV_SUFFIX).exists():
            record_environment(path, _python(env, atlas))
        return path
    if offline:
        raise SystemExit(f"[tp53-12nt] {_rel(path)} is missing and --offline was given")
    py = _python(env, atlas)
    if name == "spliceai_walker":
        _run([py, WALKER_SCORER, "--variants", SCORE_INPUT, "--out", path, "--run"])
    elif env == "analysis":
        worker(name, atlas)
    else:
        extra = {"HF_HUB_OFFLINE": "1", "TRANSFORMERS_OFFLINE": "1"} if env == "nt" else None
        _run([py, "-m", "src.evid_tp53_extend", "--worker", name, "--atlas", atlas], extra)
    if not path.exists():
        raise SystemExit(f"[tp53-12nt] {name} produced no {_rel(path)}")
    record_environment(path, py)
    return path


# ---------------------------------------------------------------------------
# assembly
# ---------------------------------------------------------------------------
def _as_stored(d: pd.DataFrame) -> pd.DataFrame:
    """The TSV write and default-parser read the frozen values went through."""
    buf = io.StringIO()
    d.to_csv(buf, sep="\t", index=False)
    buf.seek(0)
    return pd.read_csv(buf, sep="\t")


def _read_scores(name: str, expected: pd.DataFrame) -> pd.DataFrame:
    path = EXT_DIR / SCORERS[name]["file"]
    d = pd.read_parquet(path)
    got, want = set(d["variant_id"]), set(expected["variant_id"])
    if got != want:
        raise SystemExit(f"[tp53-12nt] {_rel(path)} covers a different variant list "
                         f"({len(got)} vs {len(want)}); delete it and rerun")
    if name in TSV_STORED:
        d = _as_stored(d)
    return d.set_index("variant_id")


def _identical(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    return (a == b) | (np.isnan(a) & np.isnan(b))


def _check(name: str, col: str, ids: list[str], got: pd.Series, want: pd.Series) -> dict:
    g = got.loc[ids].to_numpy(dtype=float)
    w = want.loc[ids].to_numpy(dtype=float)
    same = _identical(g, w)
    diff = np.abs(g - w)
    return {"scorer": name, "column": col, "n_checked": int(len(ids)),
            "n_identical": int(same.sum()),
            "max_abs_diff": float(np.nanmax(diff)) if np.isfinite(diff).any() else 0.0,
            "differing": [i for i, s in zip(ids, same) if not s]}


def assemble() -> tuple[pd.DataFrame, dict, pd.DataFrame]:
    """The 288-row table from cached score files. No network, no model."""
    from . import config as C
    from .phase1_build_frozen_matrix_v2 import canonical_sha256

    mod = _scoring_script()
    frozen = pd.read_parquet(FROZEN_TABLE)
    raw = deposit_rows(mod)
    allv, inp = variant_lists(raw, frozen)
    fz_ids = frozen["variant_id"].tolist()
    new_ids = allv.loc[~allv["in_frozen_table"], "variant_id"].tolist()
    check_ids = inp.loc[inp["role"] == "reproduction_check", "variant_id"].tolist()
    if len(raw) != OFFSET_MAX * PER_OFFSET or len(new_ids) != len(raw) - len(frozen):
        raise SystemExit(f"[tp53-12nt] deposit gives {len(raw)} rows, {len(new_ids)} new")
    if not set(fz_ids) <= set(raw["variant_id"]):
        raise SystemExit("[tp53-12nt] a frozen row is not in the deposit parse")

    # --- score files: new-row values and the reproduction checks
    canon_of = {src: canon for canon, src in C.COLUMNS.items() if src is not None}
    fz = frozen.set_index("variant_id")
    walker_fz = pd.read_parquet(FROZEN_WALKER).set_index("variant_id")
    avi_fz = pd.read_parquet(FROZEN_AVI)
    avi_fz = avi_fz[avi_fz["set"] == "tp53"].set_index("variant_id")
    new_raw = pd.DataFrame(index=pd.Index(new_ids, name="variant_id"))
    checks, extra = [], {}
    for name, spec in SCORERS.items():
        expected = inp if spec["input"] == "check" else allv
        d = _read_scores(name, expected)
        ids = check_ids if spec["input"] == "check" else fz_ids
        for col, target in spec["columns"].items():
            if name == "spliceai_walker":
                ref = walker_fz[col]
            elif name == "avi":
                ref = avi_fz[col]
            else:
                ref = fz[canon_of[target]]
            checks.append(_check(name, col, ids, d[col], ref))
            if name in ("spliceai_walker", "avi"):
                extra[col] = d[col]
            else:
                new_raw[target] = d.loc[new_ids, col].to_numpy(dtype=float)
    failed = [c for c in checks if c["n_identical"] != c["n_checked"]]
    if failed:
        lines = [f"  {c['scorer']}/{c['column']}: {c['n_identical']}/{c['n_checked']} "
                 f"identical, max |diff| {c['max_abs_diff']:.3g}" for c in failed]
        raise SystemExit("[tp53-12nt] re-scored frozen rows do not reproduce the frozen "
                         "values -- not the same scorer:\n" + "\n".join(lines))

    # --- the frozen rows rebuilt through the scoring script's transforms must come
    # out identical to the frozen table; that is the check on the rule used for the 96
    filled = raw.copy()
    for src in mod.PREDICTOR_COLS:
        canon = canon_of[src]
        vals = pd.Series(np.nan, index=filled["variant_id"], dtype=float)
        vals.loc[fz_ids] = fz.loc[fz_ids, canon].to_numpy(dtype=float)
        if src in new_raw.columns:
            vals.loc[new_ids] = new_raw[src].to_numpy(dtype=float)
        filled[src] = vals.to_numpy()
    # through parquet, as the frozen rows were: missing text becomes None, not pd.NA
    rebuilt = pd.read_parquet(io.BytesIO(mod.to_frozen_schema(filled).to_parquet(index=False)))
    reb_fz = rebuilt[rebuilt["variant_id"].isin(fz_ids)].reset_index(drop=True)
    pd.testing.assert_frame_equal(reb_fz, frozen.reset_index(drop=True), check_exact=True)
    new_rows = rebuilt[~rebuilt["variant_id"].isin(fz_ids)]

    table = pd.concat([frozen, new_rows], ignore_index=True)
    table = table.sort_values("variant_id").reset_index(drop=True)
    walker = pd.concat([walker_fz.loc[fz_ids, WALKER_COLS],
                        pd.DataFrame({c: extra[c].loc[new_ids] for c in WALKER_COLS})])
    avi = pd.DataFrame({c: extra[c] for c in AVI_COLS})
    for c in WALKER_COLS:
        table[c] = walker.loc[table["variant_id"], c].to_numpy(dtype=float)
    for c in AVI_COLS:
        table[c] = avi.loc[table["variant_id"], c].to_numpy(dtype=float)

    # --- invariants
    part = table[table["variant_id"].isin(fz_ids)][list(frozen.columns)].reset_index(drop=True)
    pd.testing.assert_frame_equal(part, frozen.reset_index(drop=True), check_exact=True)
    counts = table["intron_offset"].abs().astype(int).value_counts().sort_index()
    if list(counts.index) != list(range(1, OFFSET_MAX + 1)) or (counts != PER_OFFSET).any():
        raise SystemExit(f"[tp53-12nt] rows per |offset| are not {PER_OFFSET}: {counts}")
    for c in REQUIRED_COMPLETE:
        if table[c].isna().any():
            raise SystemExit(f"[tp53-12nt] {c}: {int(table[c].isna().sum())} rows unscored")

    new_mask = ~table["variant_id"].isin(fz_ids)
    columns = {}
    for c in table.columns:
        if c in DEFINITIONS or c in WALKER_COLS:
            columns[c] = {
                "n_non_null": int(table[c].notna().sum()),
                "n_non_null_frozen_rows": int(table.loc[~new_mask, c].notna().sum()),
                "n_non_null_new_rows": int(table.loc[new_mask, c].notna().sum()),
                "definition": DEFINITIONS.get(c, DEFINITIONS["spliceai_walker"]),
            }
    stratum = table["intron_offset"].abs().map(
        lambda o: "pm12" if o <= 2 else ("s3_10" if o <= 10 else "s11_50"))
    derived = {
        "variant_id / chrom / pos / ref / alt": (
            "phase4a_score_tp53.load_tp53_splice: the deposit's HGVS(genomic) "
            "g.-substitution on GRCh38 chr17, plus strand"),
        "hgvs_nt": "the deposit's HGVS(cDNA), NM_000546.6",
        "intron_offset": "signed offset parsed from HGVS(cDNA) (c.N+k / c.N-k)",
        "func_score": "the deposit's `score` column (relative fitness score)",
        "func_pathogenicity": (
            f"+func_score (TP53 is in config.FLIP_GENES = {C.FLIP_GENES}); "
            "phase1_build_frozen_matrix.orient_functional"),
        "offset_source / splice_side": (
            "phase1_build_frozen_matrix.derive_intron_offset: 'upstream'; donor if the "
            "offset is positive, acceptor if negative"),
        "region / is_splice": (
            "phase1_build_frozen_matrix.assign_region from region_label 'splice' "
            "(phase4a_score_tp53 sets it for every row): True for all 288"),
        "splice_subclass / splice_bin / analysis_region": (
            f"assign_region: |offset| <= {C.SPLICE_CORE_MAX} splice_core, "
            f"{C.SPLICE_CORE_MAX + 1}..{C.SPLICE_REGION_MAX} splice_region, "
            f"> {C.SPLICE_REGION_MAX} intronic (so the 96 new rows are 'intronic', "
            "splice_bin 'region', analysis_region 'intronic')"),
        "y_clinvar / y_assay / y_assay_source": "missing, as for the frozen rows",
        "*_isna": "phase1_build_frozen_matrix.add_missingness_flags",
    }
    inputs = {_rel(p): sha256_file(p) for p in [
        DEPOSIT, FROZEN_TABLE, FROZEN_MANIFEST, FROZEN_WALKER, FROZEN_AVI, FROZEN_SUBSET,
        FROZEN_EVENTS, WALKER_CFG, SCORE_TP53, SCORE_PINNED, AG_V061, WALKER_SCORER,
        EVENT_SCORER, PHASE1 / "src/evid_score_avi.py", PHASE1 / "src/config.py",
        PHASE1 / "src/phase1_build_frozen_matrix.py",
        PHASE1 / "data/frozen/frozen_matrix_v1.parquet", Path(__file__),
        *[EXT_DIR / s["file"] for s in SCORERS.values()],
        *[Path(str(EXT_DIR / s["file"]) + ".provenance.json") for s in SCORERS.values()],
    ]}
    manifest = {
        "product": ("TP53 intron-side SNVs at 1 <= |offset| <= 12 from MaveDB "
                    "urn:mavedb:00001213-a-1 on NM_000546.6, GRCh38: the 192 frozen rows "
                    "of tp53_splice_scored_v2.parquet unchanged plus the 96 at "
                    "|offset| 9-12, scored by the same scorers"),
        "n_rows": int(len(table)), "n_frozen_rows": int((~new_mask).sum()),
        "n_new_rows": int(new_mask.sum()),
        "rows_per_abs_offset": {str(k): int(v) for k, v in counts.items()},
        "rows_per_stratum": {k: int(v) for k, v in stratum.value_counts().sort_index().items()},
        "frozen_rows": {"source": _rel(FROZEN_TABLE),
                        "identical_to_source": True,
                        "rebuilt_by_the_same_transforms_and_identical": True},
        "reference_alleles": None,       # filled by build() when the fasta is present
        "columns": columns,
        "not_scored": NOT_SCORED,
        "required_complete": REQUIRED_COMPLETE,
        "derived_fields": derived,
        "reproduction_checks": checks,
        "score_files": {name: {"file": _rel(EXT_DIR / s["file"]),
                               "run_on": ("the 96 new rows and 24 frozen rows re-scored "
                                          "as a check" if s["input"] == "check" else
                                          "all 288 rows"),
                               "columns": list(s["columns"]),
                               "stored_as": ("TSV write + pandas.read_csv default parse, "
                                             "as the frozen rows (see TSV_STORED)"
                                             if name in TSV_STORED else "as scored"),
                               "environment": _environment_summary(EXT_DIR / s["file"],
                                                                   s["env"])}
                        for name, s in SCORERS.items()},
        "inframe": None,                 # filled by build() once the events exist
        "inputs_sha256": inputs,
        "sha256": canonical_sha256(table),
        "sha256_recipe": ("rows sorted by variant_id, CSV serialisation "
                          "(phase1_build_frozen_matrix_v2.canonical_sha256), the same "
                          "content hash the v2 manifest carries"),
    }
    return table, manifest, allv


# ---------------------------------------------------------------------------
# E2.7: the false-positive subset and its events
# ---------------------------------------------------------------------------
def inframe_subset(table: pd.DataFrame) -> pd.DataFrame:
    """evid_inframe.external_subset's TP53 rule, applied to the 288.

    external_subset merges the 192-row Walker file onto whatever table it is given,
    so it cannot be pointed at a table that already carries the column; the rule is
    restated here and must reproduce the existing 192-row subset exactly."""
    import yaml
    pp3 = yaml.safe_load(WALKER_CFG.read_text())["thresholds"]["pp3"]["value"]
    df = table[table["is_splice"]].copy()
    df["stratum"] = df["intron_offset"].abs().map(
        lambda o: "pm12" if o <= 2 else ("s3_10" if o <= 10 else "s11_50"))
    df["y_tp53_control_anchored"] = (
        df["func_pathogenicity"].to_numpy(dtype=float) > 0).astype(float)
    col = "spliceai_walker"
    sub = df[(df[col] >= pp3) & (df["y_tp53_control_anchored"] == 0)].copy()
    sub["score_column"] = col
    keep = ["variant_id", "gene", "chrom", "pos", "ref", "alt", "stratum", col,
            "score_column", "y_tp53_control_anchored"]
    sub = sub[keep].reset_index(drop=True)
    old = pd.read_parquet(FROZEN_SUBSET)
    mine = sub[sub["variant_id"].isin(pd.read_parquet(FROZEN_TABLE)["variant_id"])]
    pd.testing.assert_frame_equal(mine.reset_index(drop=True),
                                  old[keep].reset_index(drop=True), check_exact=True)
    return sub


def ensure_events(sub: pd.DataFrame, atlas: Path, offline: bool) -> dict:
    SUBSET_OUT.parent.mkdir(parents=True, exist_ok=True)
    sub.to_parquet(SUBSET_OUT, index=False)
    if EVENTS_OUT.exists():
        ev = pd.read_parquet(EVENTS_OUT)
        if set(ev["variant_id"]) != set(sub["variant_id"]):
            raise SystemExit(f"[tp53-12nt] {_rel(EVENTS_OUT)} was scored for a different "
                             "subset; delete it and rerun")
    elif offline:
        raise SystemExit(f"[tp53-12nt] {_rel(EVENTS_OUT)} is missing and --offline was given")
    else:
        py = _python("spliceai", atlas)
        _run([py, EVENT_SCORER, "--variants", SUBSET_OUT,
              "--distance", "4999", "--out", EVENTS_OUT, "--run"])
        record_environment(EVENTS_OUT, py)
        ev = pd.read_parquet(EVENTS_OUT)
    old = pd.read_parquet(FROZEN_EVENTS).set_index("variant_id")
    evi = ev.set_index("variant_id")
    num = [c for c in old.columns if c.startswith("spliceai_") and c != "spliceai_event_symbol"]
    ids = list(old.index)
    same = all(_identical(evi.loc[ids, c], old.loc[ids, c]).all() for c in num)
    same = same and (evi.loc[ids, "spliceai_event_symbol"] == old.loc[ids, "spliceai_event_symbol"]).all()
    if not same:
        raise SystemExit("[tp53-12nt] the event table does not reproduce the 192-row one")
    # the event maximum is the Walker-basis maximum: same scorer, same distance
    if not _identical(evi.loc[sub["variant_id"], "spliceai_event_max"],
                      sub.set_index("variant_id")["spliceai_walker"]).all():
        raise SystemExit("[tp53-12nt] event maxima disagree with the Walker column")
    return {"subset": _rel(SUBSET_OUT), "subset_sha256": sha256_file(SUBSET_OUT),
            "events": _rel(EVENTS_OUT), "events_sha256": sha256_file(EVENTS_OUT),
            "rule": ("spliceai_walker >= the Walker PP3 cut point (config/walker2023.yaml) "
                     "and func_pathogenicity <= 0 (control-anchored label normal), "
                     "as evid_inframe.external_subset"),
            "n_subset": int(len(sub)),
            "by_stratum": {k: int(v) for k, v in sub["stratum"].value_counts().sort_index().items()},
            "frozen_rows_reproduced": int(len(ids)),
            "new_rows": int((~sub["variant_id"].isin(ids)).sum()),
            "events_environment": _environment_summary(EVENTS_OUT, "spliceai"),
            "pangolin_events": ("not scored: there is no Pangolin event table for the 192 "
                                "frozen rows either (evid_inframe's tp53_pangolin.parquet "
                                "does not exist), so TP53 attribution stays SpliceAI-only")}


# ---------------------------------------------------------------------------
def build(offline: bool = False) -> dict:
    atlas = _atlas_repo()
    write_variant_lists()
    for name in SCORERS:
        ensure_scored(name, atlas, offline)
    table, manifest, allv = assemble()
    manifest["inframe"] = ensure_events(inframe_subset(table), atlas, offline)
    manifest["reference_alleles"] = _check_reference_alleles(allv, atlas)
    table.to_parquet(OUT, index=False)
    manifest["parquet_sha256"] = sha256_file(OUT)
    MANIFEST.write_text(json.dumps(manifest, indent=2, default=str) + "\n")
    return manifest


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--offline", action="store_true",
                    help="use cached score files only; fail if one is missing")
    ap.add_argument("--worker", choices=sorted(SCORERS), help=argparse.SUPPRESS)
    ap.add_argument("--atlas", default=None, help=argparse.SUPPRESS)
    args = ap.parse_args()
    if args.worker:
        worker(args.worker, Path(args.atlas) if args.atlas else _atlas_repo())
        return 0
    m = build(args.offline)
    print(f"\n[tp53-12nt] {m['n_rows']} rows ({m['n_frozen_rows']} frozen + "
          f"{m['n_new_rows']} new); per |offset| {m['rows_per_abs_offset']}")
    print(f"     strata {m['rows_per_stratum']}")
    for c, v in m["columns"].items():
        print(f"     {c:24s} {v['n_non_null']:>4}/{m['n_rows']}  "
              f"(new rows {v['n_non_null_new_rows']}/{m['n_new_rows']})")
    ok = sum(c["n_identical"] == c["n_checked"] for c in m["reproduction_checks"])
    print(f"     reproduction checks: {ok}/{len(m['reproduction_checks'])} columns identical")
    print(f"     E2.7 subset {m['inframe']['n_subset']} ({m['inframe']['by_stratum']})")
    print(f"     wrote {_rel(OUT)}  sha256 {m['parquet_sha256'][:12]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
