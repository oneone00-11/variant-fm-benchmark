#!/usr/bin/env python3
"""
evidence-strength reframe -- one command for every analysis stage of the reworked study.

What this reproduces, and what it does not. The analysis stages below run in order
and are deterministic under the seed in `phase1/src/config.py`; several of them,
including the build of the analysis set, also read the companion atlas, release
v2.5.0: EVID_ATLAS_REPO, or a checkout named functional-standard-atlas beside this
one, or the release archive that --fetch-inputs downloads from Zenodo
(10.5281/zenodo.22751081, md5 checked) and unpacks under phase1/data/evidence/
(gitignored). The atlas's results/ files ship in that archive and not in its git
repository, so a bare clone of the atlas is not enough; the entry point checks for
every atlas file the stages read before it starts. Five kinds of input are NOT
rebuilt here, because they are third-party downloads or model re-scores that need
their own environments; each is pinned by sha256 in a manifest or provenance file:

  * the ClinVar GRCh38 VCF of 15 June 2026 (192 MB; `phase1/data/evidence/clinvar/`,
    gitignored; --fetch-inputs downloads it from NCBI's archive and checks its md5
    against NCBI's file and its sha256 against the analysis-set manifest);
  * the `spliceai_walker` column, re-scored at Walker's -D 4999 in the pinned
    SpliceAI 1.3.1 environment (`src/evid_score_spliceai_walker.py`), and the
    SpliceAI event records for the in-frame attribution
    (`src/evid_score_spliceai_events.py`);
  * the DDX3X deposit and its labels (`src/evid_external.py --prepare ddx3x`, which
    fetches the MaveDB score sets and the deposit's classification) and its panel
    scores;
  * the TP53 scores for offsets 9-12 (`src/evid_tp53_extend.py`, whose model scores
    are cached; the stage below reruns it offline from that cache);
  * the AlphaGenome Atlas columns (`src/evid_score_avi.py`), which need an API key
    outside the tree and a separate venv carrying alphagenome>=0.9.0 -- the pinned
    .venv keeps 0.7.0 because that is the provenance of the alphagenome column.

Each is tracked with a checksum or a provenance record, except four DDX3X score
lookups (CADD, GPN-MSA, phyloP, phastCons). The ClinVar VCF and the SpliceAI event
records print the command that produces them when they are missing; the module
named above documents the others.

Stages
------
  1  evid_build_set          rebuild the 1 <= |offset| <= 50 analysis set (E1)
  2  evid_walker_thresholds  the ClinGen fixed cut points, by stratum and arm (E2)
  3  evid_interval_lr        score-to-evidence intervals, Pejaver's procedure (E3)
  4  evid_territory_metrics  AUROC / PR-AUC by territory and ClinVar arm (E4)
  5  evid_inframe            what the false positives are predicting (E6)
  6  evid_inframe DDX3X      the same attribution on DDX3X (E2.7)
  7  evid_tp53_extend        TP53 to |offset| 12, offline from cached scores (E7)
  8  evid_inframe TP53       the same attribution on TP53 (E2.7)
  9  evid_external TP53      fixed and fitted thresholds on TP53 (E7)
 10  evid_external DDX3X     merge the scored columns (E7)
 11  evid_external DDX3X     fixed and fitted thresholds on DDX3X (E7)
 12  evid_diagnostics        cut-point, monotonicity, concordance and depth tables (E8)
 13  evid_tier_logo          in-sample tier vs held-out ratio, per fold (E2.2)
 14  evid_arms               ClinVar arms without BRCA1, and within gene (E2.3/E2.4)
 15  evid_fig_data           the figure tables and draft PNGs (E2.9)
 16  evid_training_provenance where each predictor's training signal comes from (E9)
 17  evid_fusion_stability   the fusion's coefficients in every fold (E10)
 18  evid_dilution           thresholds refitted on every subset of training genes (E11)
 19  evid_tables             the four main tables as printed (E12)
 20  evid_supp_tables        the supplementary tables as printed (E13)
 21  evid_figures            main and supplementary figures, PDF and PNG (E14)
 22  evid_practical_metrics  sensitivity, specificity and ROC at the evidence thresholds (E15)
 23  evid_delta              every published quantity with a counterpart (E8)

Outputs land in `phase1/reports/evidence/`; stages 22 and 23 also write a page under
`docs/`.

    bash reproduce.sh                                       # from a fresh clone: everything
    bash reproduce.sh --check                               # five minutes, no download
    python scripts/reproduce_evidence.py --fetch-inputs --verify
    python scripts/reproduce_evidence.py
    python scripts/reproduce_evidence.py --from 3      # resume at a stage
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import http.client
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.request
import zipfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PHASE1 = REPO / "phase1"

# The two public inputs this repository does not carry, and where --fetch-inputs
# gets them. NCBI's archive keeps every release; the weekly directory is the
# fallback. The sha256 is the one phase1/data/evidence/analysis_set_v1.manifest.json
# records (a test keeps the two equal).
CLINVAR_VCF = PHASE1 / "data" / "evidence" / "clinvar" / "clinvar_20260615.vcf.gz"
# Tried in order: (url, url of the source's own md5 file, or None). A mirror without
# an md5 file is checked by the sha256 alone, which is the check that matters.
_NCBI = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/"
CLINVAR_SOURCES = [
    (_NCBI + "archive_2.0/2026/clinvar_20260615.vcf.gz",
     _NCBI + "archive_2.0/2026/clinvar_20260615.vcf.gz.md5"),
    (_NCBI + "weekly/clinvar_20260615.vcf.gz",
     _NCBI + "weekly/clinvar_20260615.vcf.gz.md5"),
    # NCBI's file, unchanged, in this release's Zenodo record (10.5281/zenodo.23051009)
    ("https://zenodo.org/api/records/23051009/files/clinvar_20260615.vcf.gz/content", None),
]
CLINVAR_URLS = [url for url, _ in CLINVAR_SOURCES]
CLINVAR_SHA256 = "10d86b892aae1f035e1950844e13fb039dad50be1087a0d1445c60d29191a342"
ATLAS_ARCHIVE = {
    "version": "2.5.0",
    "doi": "10.5281/zenodo.22751081",
    "url": ("https://zenodo.org/api/records/22751081/files/"
            "functional-standard-atlas-v2.5.0.zip/content"),
    "md5": "0858fe136346428b764e8b1c3c35b223",      # Zenodo's checksum of the file
}
# What --verify compares: every file the stages write, as it stood before the run.
# Fetched inputs are left out; manifests may differ only in their wall-clock fields,
# the set evid_supp_tables also leaves out when it hashes them.
OUTPUT_ROOTS = [PHASE1 / "reports" / "evidence", PHASE1 / "data" / "evidence",
                REPO / "docs" / "evidence-delta.md",
                REPO / "docs" / "practical-metrics-summary.md"]
NOT_OUTPUTS = [PHASE1 / "data" / "evidence" / "clinvar",
               PHASE1 / "data" / "evidence" / "companion_atlas"]
WALL_CLOCK_KEYS = {"built_utc", "scored_utc", "retrieved_utc", "run_at"}
# Tracked files under the output folders that the stages only read; every other
# tracked file there must be rewritten by a full run.
RUN_INPUTS = REPO / "scripts" / "reproduce_inputs.txt"
# The sha256 of every atlas file the stages read, as the release has them.
ATLAS_PINS = REPO / "scripts" / f"atlas_release_v{ATLAS_ARCHIVE['version']}.sha256"
# In a tree without .git (a release archive), where --verify keeps the tree as first
# unpacked, so that a rerun after an interrupted run is still compared with it.
BASELINE = REPO / ".reproduce_baseline"

# One rule for finding the atlas, shared with the stages and the tests.
sys.path.insert(0, str(PHASE1))
from src.evid_common import (ATLAS_BESIDE, ATLAS_FETCHED, ATLAS_NEEDS,  # noqa: E402
                             ATLAS_RELEASE_COMMIT, resolve_atlas)

STAGES = [
    ("src.evid_build_set",          "E1  rebuild the analysis set", []),
    ("src.evid_walker_thresholds",  "E2  Walker fixed cut points", []),
    ("src.evid_interval_lr",        "E3  score-to-evidence intervals", []),
    ("src.evid_territory_metrics",  "E4  territory and ClinVar-arm metrics", []),
    ("src.evid_inframe",            "E6  in-frame attribution", ["--attribute"]),
    ("src.evid_inframe",            "E2.7 in-frame attribution: DDX3X", ["--external-attribute", "ddx3x"]),
    ("src.evid_tp53_extend",        "E7  TP53 extended to |offset| 12 (offline)", ["--offline"]),
    ("src.evid_inframe",            "E2.7 in-frame attribution: TP53", ["--external-attribute", "tp53"]),
    ("src.evid_external",           "E7  external gene: TP53", ["--apply", "TP53"]),
    ("src.evid_external",           "E7  external gene: DDX3X, merge scored columns",
     ["--merge-scores", "ddx3x"]),
    ("src.evid_external",           "E7  external gene: DDX3X", ["--apply", "DDX3X"]),
    ("src.evid_diagnostics",        "E8  cut-point and monotonicity diagnostics", []),
    ("src.evid_tier_logo",          "E2.2 in-sample tier against held-out ratio", []),
    ("src.evid_arms",               "E2.3/E2.4 ClinVar arms without the gene confound", []),
    ("src.evid_fig_data",           "E2.9 one tidy table per figure, plus draft PNGs", []),
    ("src.evid_training_provenance", "E9  predictor training signals and overlaps", []),
    ("src.evid_fusion_stability",   "E10 fusion coefficients across folds", []),
    ("src.evid_dilution",           "E11 thresholds refitted on subsets of training genes", []),
    ("src.evid_tables",             "E12 the four main tables, as printed", []),
    ("src.evid_supp_tables",        "E13 the supplementary tables, as printed", []),
    ("src.evid_figures",            "E14 main and supplementary figures", []),
    ("src.evid_practical_metrics",  "E15 practical metrics", []),
    ("src.evid_delta",              "E8  old/new quantity list -> docs/evidence-delta.md", []),
]

# Stages that need an input this script does not produce, and what produces it.
NEEDS = {
    "src.evid_tier_logo": [
        ("phase1/reports/evidence/evidence_thresholds_logo_folds.csv",
         "python -m src.evid_interval_lr (stage 3)"),
    ],
    "src.evid_practical_metrics": [
        ("phase1/reports/evidence/evidence_thresholds_logo_folds.csv",
         "python -m src.evid_interval_lr (stage 3)"),
    ],
    "src.evid_build_set": [
        ("phase1/data/evidence/clinvar/clinvar_20260615.vcf.gz",
         "python scripts/reproduce_evidence.py --fetch-inputs, or curl -O "
         + CLINVAR_URLS[0] + "   (into phase1/data/evidence/clinvar/)"),
    ],
    "src.evid_inframe": [
        ("phase1/data/evidence/inframe_events.parquet",
         "python -m src.evid_inframe --subset, then the SpliceAI-environment command "
         "it prints"),
    ],
}


def _digest(path: Path, algo: str) -> str:
    h = hashlib.new(algo)
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


# --------------------------------------------------------------------------
# --fetch-inputs: the two public inputs this repository does not carry
# --------------------------------------------------------------------------
def _permanent(e: urllib.error.HTTPError) -> bool:
    """A client error that another attempt will not cure (a 408 or 429 may)."""
    return 400 <= e.code < 500 and e.code not in (408, 429)


def _expected_size(r, offset: int) -> int | None:
    """The full size of the file the response belongs to, when the server says."""
    cr = r.headers.get("Content-Range", "")               # bytes a-b/total
    if "/" in cr and cr.rsplit("/", 1)[1].isdigit():
        return int(cr.rsplit("/", 1)[1])
    cl = r.headers.get("Content-Length", "")
    return offset + int(cl) if cl.isdigit() else None


def _download(url: str, dest: Path, attempts: int = 30) -> None:
    """Download `url` to `dest` through a .part file that survives a stall or a
    dropped connection: each retry asks for the remaining bytes only (NCBI and Zenodo
    both honour HTTP ranges). A minute without data counts as a stall; a client error
    such as 404 ends at once. The caller checks the checksum, so a resumed file is
    never trusted on its own."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    for attempt in range(1, attempts + 1):
        have = part.stat().st_size if part.exists() else 0
        headers = {"Range": f"bytes={have}-"} if have else {}
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=headers),
                                        timeout=60) as r:
                resume = bool(have) and r.status == 206   # else it sends it all
                total = _expected_size(r, have if resume else 0)
                with open(part, "ab" if resume else "wb") as f:
                    shutil.copyfileobj(r, f, 1 << 20)
            # a connection closed early ends the read without an error, so the size
            # the server announced is the only sign the body is incomplete
            if total is None or part.stat().st_size >= total:
                part.replace(dest)
                return
            err = f"connection closed before {total} bytes"
        except urllib.error.HTTPError as e:
            if e.code == 416 and have:                 # nothing left to send
                part.replace(dest)
                return
            if _permanent(e):
                raise
            err = e
        except (OSError, http.client.HTTPException) as e:   # stalls, resets, a
            err = e                                           # connection cut short
        if attempt == attempts:
            break
        got = part.stat().st_size if part.exists() else 0
        print(f"  interrupted at {got / 2**20:.0f} MB ({err}); resuming, attempt "
              f"{attempt + 1} of {attempts}", flush=True)
        time.sleep(min(30, 2 * attempt))
    raise OSError(f"gave up on {url} after {attempts} attempts")


def _fetch_md5(url: str, attempts: int = 3) -> str:
    """A source's own md5 file, read with a few retries; anything that is not an md5
    (an error page, say) counts as the source failing."""
    for attempt in range(1, attempts + 1):
        try:
            text = urllib.request.urlopen(url, timeout=60).read().decode("ascii", "replace")
            m = re.match(r"\s*([0-9a-fA-F]{32})\b", text)
            if m:
                return m.group(1).lower()
            err = f"no md5 in {url}"
        except urllib.error.HTTPError as e:
            if _permanent(e):
                raise
            err = e
        except (OSError, http.client.HTTPException) as e:
            err = e
        if attempt < attempts:
            time.sleep(2 * attempt)
    raise OSError(f"could not read an md5 from {url} ({err})")


def fetch_clinvar() -> None:
    """The ClinVar release, from the first source that serves the exact file: its own
    md5, where the source has one, and the sha256 the analysis-set manifest records.
    A source whose file fails either is set aside for the next one. An existing file is
    only checked."""
    rel = CLINVAR_VCF.relative_to(REPO)
    if CLINVAR_VCF.exists():
        if _digest(CLINVAR_VCF, "sha256") != CLINVAR_SHA256:
            sys.exit(f"{rel} is not the release the analysis set was built from (sha256 "
                     "differs); remove it and rerun")
        print(f"ClinVar release present, sha256 checked: {rel}")
        return
    part = CLINVAR_VCF.with_name(CLINVAR_VCF.name + ".part")
    for url, md5_url in CLINVAR_SOURCES:
        try:
            md5 = _fetch_md5(md5_url) if md5_url else None
            print(f"downloading {url} (192 MB)", flush=True)
            _download(url, CLINVAR_VCF)
        except (OSError, http.client.HTTPException) as e:
            print(f"  not available there ({e})")
            part.unlink(missing_ok=True)          # never resume one source's bytes
            continue                              # from another
        if md5 is not None and _digest(CLINVAR_VCF, "md5") != md5:
            print("  the file does not match that source's own md5; trying the next source")
        elif _digest(CLINVAR_VCF, "sha256") != CLINVAR_SHA256:
            print("  that source's file is not the release the analysis set was built "
                  "from (sha256 differs); trying the next source")
        else:
            print(("  md5 matches the source's, " if md5 else "  ")
                  + "sha256 matches the analysis-set manifest")
            return
        CLINVAR_VCF.unlink(missing_ok=True)
    sys.exit("no source served the ClinVar release the analysis set was built from")


def fetch_atlas() -> None:
    """The atlas release archive, from Zenodo, checked against Zenodo's md5 and
    unpacked to phase1/data/evidence/companion_atlas/atlas/."""
    a = ATLAS_ARCHIVE
    root = ATLAS_FETCHED.parent
    zpath = root / f"functional-standard-atlas-v{a['version']}.zip"
    if not (zpath.exists() and _digest(zpath, "md5") == a["md5"]):
        print(f"downloading the atlas release v{a['version']} (doi {a['doi']}, 50 MB)",
              flush=True)
        try:
            _download(a["url"], zpath)
        except (OSError, http.client.HTTPException) as e:
            sys.exit(f"could not download the atlas release archive ({e}). Download "
                     f"https://doi.org/{a['doi']} by hand and set EVID_ATLAS_REPO to "
                     "its unpacked atlas/ folder.")
    if _digest(zpath, "md5") != a["md5"]:
        zpath.unlink()
        sys.exit("the downloaded atlas archive failed its md5 check and was removed")
    with zipfile.ZipFile(zpath) as z:
        for name in z.namelist():            # nothing may land outside root
            if not (root / name).resolve().is_relative_to(root.resolve()):
                sys.exit(f"unexpected path in the atlas archive: {name}")
        z.extractall(root)
    print(f"  md5 matches Zenodo's; unpacked to {ATLAS_FETCHED.relative_to(REPO)}")


# --------------------------------------------------------------------------
# which atlas the stages read
# --------------------------------------------------------------------------
def _pins() -> dict[str, str]:
    pins = {}
    for line in ATLAS_PINS.read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            sha, rel = line.split(maxsplit=1)
            pins[rel.strip()] = sha
    return pins


def _git_head(path: Path) -> str | None:
    """HEAD of `path` if `path` is itself the top of a git checkout; an archive
    unpacked inside this repository is not."""
    try:
        run = lambda *a: subprocess.run(["git", "-C", str(path), *a], check=True,  # noqa: E731
                                        capture_output=True, text=True).stdout.strip()
        if Path(run("rev-parse", "--show-toplevel")).resolve() != Path(path).resolve():
            return None
        return run("rev-parse", "HEAD")
    except (OSError, subprocess.CalledProcessError):
        return None


def atlas_release_mismatch(path: Path) -> list[str]:
    """What keeps `path` from being the atlas release the outputs were made from: files
    that differ from the release's pinned sha256, and, for a git checkout, a HEAD other
    than the release commit, which the analysis-set manifest records."""
    path = Path(path)
    bad = [f"{rel} {'missing' if not (path / rel).is_file() else 'differs'}"
           for rel, sha in _pins().items()
           if not (path / rel).is_file() or _digest(path / rel, "sha256") != sha]
    head, want = _git_head(path), ATLAS_RELEASE_COMMIT[ATLAS_ARCHIVE["version"]]
    if head is not None and head != want:
        bad.append(f"a checkout at {head[:12]}, not the release commit {want[:12]}")
    return bad


def atlas_version(atlas: Path) -> str:
    for line in (atlas / "CITATION.cff").read_text().splitlines():
        if line.startswith("version:"):
            return line.split(":", 1)[1].strip().strip('"').strip("'")
    return "unknown"


def choose_atlas(release_only: bool, fetch: bool) -> Path:
    """The atlas every stage will read. With `release_only` (--verify), only a copy
    that is the release file for file will do: EVID_ATLAS_REPO if it names one, else
    the checkout beside this repository if it is one, else the fetched archive, which
    --fetch-inputs downloads when needed. Otherwise the shared rule of
    evid_common.resolve_atlas, which only asks that no file be missing."""
    env = os.environ.get("EVID_ATLAS_REPO")
    ver = ATLAS_ARCHIVE["version"]
    if not release_only:
        atlas, missing = resolve_atlas()
        if missing and fetch and not env:
            fetch_atlas()
            atlas, missing = resolve_atlas()
        if missing:
            hint = ("point EVID_ATLAS_REPO at an unpacked copy of the release archive, or "
                    "unset it and rerun with --fetch-inputs" if env else
                    "rerun with --fetch-inputs to download the release archive")
            sys.exit(f"The companion atlas at {atlas} lacks {', '.join(missing)}.\n"
                     "Its results/ files ship in the release archive "
                     f"(doi {ATLAS_ARCHIVE['doi']}), not in its git repository; {hint}.")
        return atlas
    if env:
        atlas = Path(env).expanduser().resolve()
        bad = atlas_release_mismatch(atlas)
        if bad:
            sys.exit(f"EVID_ATLAS_REPO={env} is not the atlas release {ver} file for "
                     f"file ({'; '.join(bad[:3])}{'; ...' if len(bad) > 3 else ''}). "
                     "Unset it, and --fetch-inputs downloads the release archive.")
        return atlas
    for atlas in (ATLAS_BESIDE, ATLAS_FETCHED):
        if atlas.is_dir() and not atlas_release_mismatch(atlas):
            return atlas
    if not fetch:
        sys.exit(f"No copy of the atlas release {ver} is here; rerun with --fetch-inputs.")
    fetch_atlas()
    bad = atlas_release_mismatch(ATLAS_FETCHED)
    if bad:
        sys.exit(f"the unpacked atlas archive does not match the release's pinned "
                 f"checksums: {'; '.join(bad[:3])}")
    return ATLAS_FETCHED


# --------------------------------------------------------------------------
# --verify: every output compared with the committed one
# --------------------------------------------------------------------------
# How --verify reads a file that is not byte-identical. Tables are compared cell by
# cell, floats to a relative tolerance and integers and strings exactly, so that a
# report from another platform says whether a difference is in the numbers; figures
# are named as such. Neither refinement changes the verdict: REPRODUCED means every
# output is byte-identical, a manifest's build time aside, and was rewritten by the run.
TABLE_DELIMITERS = {".csv": ",", ".tsv": "\t"}
COPY_SUFFIXES = {".csv", ".tsv", ".json", ".parquet", ".txt"}
FIGURE_SUFFIXES = {".png", ".pdf", ".svg", ".tif", ".tiff"}
RTOL = 1e-9
REPORT_NAME = "reproduce_report.txt"      # at the repository root, gitignored
OK_STATUSES = ("identical", "identical (build time aside)")
# The second tier. On another machine the maths libraries (Apple's Accelerate, which
# numpy calls on macOS, and the OpenBLAS scipy bundles, which scikit-learn's coordinate
# descent calls) can move a float's last digits: on a GitHub M1 runner every printed
# table and figure matched byte for byte, while the tables that carry the elastic-net
# combination's scores agreed only to 1e-9.
# REPRODUCED NUMERICALLY allows exactly that: files outside the printed tables and
# figures may differ within RTOL, and a manifest may differ only in the hashes of such
# files. A printed table or figure, or any other difference, still fails.
NUMERIC_STATUSES = OK_STATUSES + ("differs (numerically equal)",
                                  "differs (hashes of numerically equal files)")
PRINTED = [PHASE1 / "reports" / "evidence" / "tables",
           PHASE1 / "reports" / "evidence" / "figures"]
_WALL_CLOCK = re.compile(r'("(?:' + "|".join(map(re.escape, sorted(WALL_CLOCK_KEYS)))
                         + r')"\s*:\s*)"[^"]*"')


def _in_roots(p: Path) -> bool:
    return (any(p == r or p.is_relative_to(r) for r in OUTPUT_ROOTS)
            and not any(p.is_relative_to(x) for x in NOT_OUTPUTS))


def _outputs() -> list[Path]:
    files = []
    for root in OUTPUT_ROOTS:
        for p in ([root] if root.is_file() else sorted(root.rglob("*"))):
            if (p.is_file() and not p.name.startswith(".") and "__pycache__" not in p.parts
                    and _in_roots(p)):
                files.append(p)
    return files


def _is_git_tree() -> bool:
    return (REPO / ".git").exists() and _git_head(REPO) is not None


def _git_files(*args: str) -> list[str]:
    roots = [str(r.relative_to(REPO)) for r in OUTPUT_ROOTS]
    out = subprocess.run(["git", "-C", str(REPO), "ls-files", "-z", *args, "--", *roots],
                         check=True, capture_output=True, text=True).stdout
    return [f for f in out.split("\0") if f and _in_roots(REPO / f)]


def _head_blobs(paths: list[str]):
    """(path, bytes) for each path as HEAD has it, through one `git cat-file --batch`."""
    request = "".join(f"HEAD:{p}\n" for p in paths).encode()
    out = subprocess.run(["git", "-C", str(REPO), "cat-file", "--batch"], input=request,
                         check=True, capture_output=True).stdout
    pos = 0
    for p in paths:
        nl = out.index(b"\n", pos)
        header = out[pos:nl].split()
        if header[-1] == b"missing":
            raise SystemExit(f"--verify: HEAD has no {p}")
        size = int(header[2])
        yield p, out[nl + 1: nl + 1 + size]
        pos = nl + 1 + size + 1


def baseline(keep: Path) -> dict:
    """The committed outputs the run is compared against, with copies of the tables
    and manifests. In a clone, the files as HEAD has them, whatever the working tree
    holds, so an interrupted earlier run cannot become the baseline. In a release
    archive, which has no .git, the tree as first unpacked: recorded in
    .reproduce_baseline/ on the first --verify run and reused by every later one."""
    files = {}
    if _is_git_tree():
        for rel, blob in _head_blobs(_git_files()):
            entry = {"sha256": hashlib.sha256(blob).hexdigest()}
            if Path(rel).suffix in COPY_SUFFIXES:
                entry["copy"] = keep / rel
                entry["copy"].parent.mkdir(parents=True, exist_ok=True)
                entry["copy"].write_bytes(blob)
            files[rel] = entry
        return {"files": files, "git": True, "start": time.time()}
    record = BASELINE / "baseline.json"
    if not record.exists():
        print(f"recording the tree as unpacked in {BASELINE.name}/ (no .git here)")
        for p in _outputs():
            rel = str(p.relative_to(REPO))
            files[rel] = {"sha256": _digest(p, "sha256")}
            if p.suffix in COPY_SUFFIXES:
                (BASELINE / "copies" / rel).parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(p, BASELINE / "copies" / rel)
        record.write_text(json.dumps(files, indent=0, sort_keys=True))
    files = json.loads(record.read_text())
    for rel, entry in files.items():
        if Path(rel).suffix in COPY_SUFFIXES:
            entry["copy"] = BASELINE / "copies" / rel
    return {"files": files, "git": False, "start": time.time()}


def _is_na(x) -> bool:
    return x is None or type(x).__name__ in ("NAType", "NaTType")


def _cell(x):
    """A table cell as a number when it reads as one."""
    if isinstance(x, str):
        for cast in (int, float):
            try:
                return cast(x)
            except ValueError:
                pass
    return x


def _equal(a, b) -> bool:
    if _is_na(a) or _is_na(b):
        return _is_na(a) and _is_na(b)
    a, b = _cell(a), _cell(b)
    if isinstance(a, bool) or isinstance(b, bool) or not (
            isinstance(a, (int, float)) and isinstance(b, (int, float))):
        return bool(a == b)
    if isinstance(a, int) and isinstance(b, int):
        return a == b
    if math.isnan(a) or math.isnan(b):
        return math.isnan(a) and math.isnan(b)
    return math.isclose(a, b, rel_tol=RTOL, abs_tol=0.0)


def _json_equal(a, b) -> bool:
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(_json_equal(a[k], b[k]) for k in a)
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(_json_equal(x, y) for x, y in zip(a, b))
    return _equal(a, b)


def _without_wall_clock(path: Path):
    def strip(x):
        if isinstance(x, dict):
            return {k: strip(v) for k, v in x.items() if k not in WALL_CLOCK_KEYS}
        if isinstance(x, list):
            return [strip(v) for v in x]
        return x
    return strip(json.loads(path.read_text()))


def _numerically_equal(before: Path, after: Path) -> bool:
    try:
        if after.suffix in TABLE_DELIMITERS:
            d = TABLE_DELIMITERS[after.suffix]
            with open(before, newline="") as f, open(after, newline="") as g:
                x, y = list(csv.reader(f, delimiter=d)), list(csv.reader(g, delimiter=d))
            return len(x) == len(y) and all(
                len(r) == len(t) and all(_equal(u, v) for u, v in zip(r, t))
                for r, t in zip(x, y))
        if after.suffix == ".json":
            return _json_equal(_without_wall_clock(before), _without_wall_clock(after))
        if after.suffix == ".txt":                # a text report, word by word
            x, y = before.read_text().split(), after.read_text().split()
            return len(x) == len(y) and all(_equal(u, v) for u, v in zip(x, y))
        if after.suffix == ".parquet":
            import pandas as pd
            x, y = pd.read_parquet(before), pd.read_parquet(after)
            if list(x.columns) != list(y.columns) or x.shape != y.shape:
                return False
            return all(_equal(u, v) for c in x.columns
                       for u, v in zip(x[c].tolist(), y[c].tolist()))
    except Exception:            # a file that cannot be read as a table is "differs"
        return False
    return False


def _status(p: Path, b: dict) -> str:
    if not p.exists():
        return "missing"
    if _digest(p, "sha256") == b["sha256"]:
        return "identical"
    copy = b.get("copy")
    # only the values of the wall-clock fields may differ, byte for byte otherwise
    if (p.suffix == ".json" and copy is not None
            and _WALL_CLOCK.sub(r'\1""', p.read_text())
            == _WALL_CLOCK.sub(r'\1""', copy.read_text())):
        return "identical (build time aside)"
    if p.suffix in FIGURE_SUFFIXES:
        return "differs (figure)"
    if copy is not None and _numerically_equal(copy, p):
        return "differs (numerically equal)"
    return "differs"


def _printed(p: Path) -> bool:
    """A table or figure as printed: the main tables, the figures and the
    supplementary tables (not their manifest)."""
    supp = PHASE1 / "reports" / "evidence" / "supplement"
    return (any(p.is_relative_to(r) for r in PRINTED)
            or (p.is_relative_to(supp) and p.name.startswith("table")))


def _hash_only_differences(rows: list, base: dict) -> list:
    """Re-grade a differing text file whose only differences are the hashes of files
    that are themselves numerically equal (a manifest recording them)."""
    swaps = {}
    for status, rel in rows:
        if status == "differs (numerically equal)":
            p, b = REPO / rel, base["files"][rel]
            swaps[_digest(p, "sha256")] = b["sha256"]
            if p.suffix == ".json" and b.get("copy") is not None:
                swaps[_source_sha256(p)] = _source_sha256(b["copy"])
    if not swaps:
        return rows
    out = []
    for status, rel in rows:
        b = base["files"].get(rel, {})
        if status == "differs" and b.get("copy") is not None and (REPO / rel).exists():
            try:
                text = (REPO / rel).read_text()
                for new, old in swaps.items():
                    text = text.replace(new, old)
                if text == b["copy"].read_text():
                    status = "differs (hashes of numerically equal files)"
            except (OSError, UnicodeDecodeError):
                pass
        out.append((status, rel))
    return out


def _run_inputs() -> set[str]:
    return {line.strip() for line in RUN_INPUTS.read_text().splitlines()
            if line.strip() and not line.startswith("#")}


def verify(base: dict) -> bool:
    """Every committed output compared with what this run wrote; one line per file in
    reproduce_report.txt, and a summary here. A committed output the run left
    untouched fails as "not rewritten", unless scripts/reproduce_inputs.txt lists it
    as an input the stages only read."""
    inputs, start = _run_inputs(), base["start"]
    rows = []
    for rel, b in sorted(base["files"].items()):
        p = REPO / rel
        status = _status(p, b)
        if status in OK_STATUSES and rel not in inputs and p.stat().st_mtime < start:
            status = "not rewritten"
        rows.append((status, rel))
    if base["git"]:
        new = _git_files("--others", "--exclude-standard")
    else:
        new = [str(p.relative_to(REPO)) for p in _outputs()
               if str(p.relative_to(REPO)) not in base["files"]]
    rows += [("new", r) for r in sorted(new)]
    rows = _hash_only_differences(rows, base)
    counts = {}
    for status, _ in rows:
        counts[status] = counts.get(status, 0) + 1
    byte_ok = all(status in OK_STATUSES for status, _ in rows)
    numeric_ok = all(status in NUMERIC_STATUSES and (status in OK_STATUSES
                                                    or not _printed(REPO / rel))
                     for status, rel in rows)
    ok = byte_ok or numeric_ok
    n_inputs = sum(1 for _, rel in rows if rel in inputs)
    n_last = sum(1 for status, _ in rows if status not in OK_STATUSES)
    verdict = ("REPRODUCED: every output matches the committed one byte for byte."
               if byte_ok else
               f"REPRODUCED NUMERICALLY: every printed table and figure matches byte for "
               f"byte; {n_last} other file(s) differ only in the last digits of their "
               f"numbers (relative {RTOL:g})" if numeric_ok else
               "NOT REPRODUCED: every file not marked identical")
    with open(REPO / REPORT_NAME, "w") as f:
        f.write(f"# reproduce_evidence.py --verify, against "
                f"{'HEAD' if base['git'] else 'the tree as unpacked'}\n# {verdict}\n")
        f.write("# " + "; ".join(f"{k}: {v}" for k, v in sorted(counts.items())) + "\n")
        f.write(f"# {n_inputs} of these are inputs the stages only read "
                f"({RUN_INPUTS.relative_to(REPO)})\n")
        for status, rel in rows:
            f.write(f"{status}\t{rel}\n")
    print(f"\n{'=' * 74}\nVerification against the committed outputs "
          f"({'HEAD' if base['git'] else 'the tree as unpacked'})\n{'=' * 74}")
    print(f"  {len(base['files'])} committed files under the output folders, {n_inputs} "
          "of them inputs the stages only read")
    for status in ("identical", "identical (build time aside)", "differs (numerically equal)",
                   "differs (hashes of numerically equal files)", "differs (figure)",
                   "differs", "not rewritten", "missing", "new"):
        if counts.get(status):
            print(f"  {status + ':':30s} {counts[status]}")
            if status != "identical":
                for st, rel in rows:
                    if st == status:
                        print(f"      {rel}")
    print(f"  one line per file: {REPORT_NAME}")
    print("\n" + (verdict if ok else verdict + " (listed above)"))
    return ok


# --------------------------------------------------------------------------
# --check: the quick path, with no download and no stage run
# --------------------------------------------------------------------------
def _canonical_sha256(path: Path) -> str:
    """The row-order-free content hash the frozen and TP53 manifests record."""
    import pandas as pd
    from src.phase1_build_frozen_matrix_v2 import canonical_sha256
    return canonical_sha256(pd.read_parquet(path))


def _source_sha256(path: Path) -> str:
    """The hash evid_supp_tables records for a source: a JSON without its wall-clock
    fields, keys sorted; any other file by its bytes."""
    if path.suffix != ".json":
        return _digest(path, "sha256")
    blob = json.dumps(_without_wall_clock(path), sort_keys=True).encode()
    return hashlib.sha256(blob).hexdigest()


def recorded_checksums() -> tuple[list[tuple[str, Path, str, str]], list[str]]:
    """(record, file, expected, how) for every checksum a tracked manifest or
    provenance record carries about a file in this repository, and about the ClinVar
    release and the atlas when they are here; and the names of those two when they are
    not. `how` is "file" (sha256 of the bytes), "content" (the row-order-free hash of a
    table) or "source" (evid_supp_tables' hash of a source). Hashes of files that only
    existed where the model scores were made (reference FASTAs, scorer modules,
    interpreters) are left out: nothing here can check them."""
    ev, out, absent = PHASE1 / "data" / "evidence", [], []
    atlas, atlas_missing = resolve_atlas()
    if atlas_missing:
        absent.append("the atlas")

    def add(record: Path, target: Path, sha: str, how: str = "file") -> None:
        out.append((str(record.relative_to(REPO)), target, sha, how))

    # provenance sidecars, in phase1/data and data/: <output>.provenance.json beside
    # the output it describes (a .tsv where the output name has no suffix)
    if _is_git_tree():                 # what the checkout carries, not stray local files
        out_ = subprocess.run(["git", "-C", str(REPO), "ls-files", "-z", "--",
                               "phase1/data/*.provenance.json", "data/*.provenance.json"],
                              check=True, capture_output=True, text=True).stdout
        records = [REPO / f for f in out_.split("\0") if f]
    else:
        records = [*(PHASE1 / "data").rglob("*.provenance.json"),
                   *(REPO / "data").rglob("*.provenance.json")]
    for rec in sorted(records):
        if any(rec.is_relative_to(x) for x in NOT_OUTPUTS):
            continue                                       # fetched inputs
        d = json.loads(rec.read_text())
        target = rec.with_name(rec.name[: -len(".provenance.json")])
        if not target.suffix:
            target = target.with_name(target.name + ".tsv")
        for key in ("output_sha256", "sha256"):
            if isinstance(d.get(key), str):
                add(rec, target, d[key])
        source = d.get("input") or re.sub(r"^.*/variant-fm-benchmark/", "",
                                          str(d.get("input_file", "")))
        if isinstance(d.get("input_sha256"), str) and source and (REPO / source).is_file():
            add(rec, REPO / source, d["input_sha256"])

    # the analysis set and the sources it was built from
    rec = ev / "analysis_set_v1.manifest.json"
    d = json.loads(rec.read_text())
    add(rec, ev / "analysis_set_v1.parquet", d["sha256"])
    for key, src in d["sources"].items():
        if key == "assay_labels":
            for spec in src.values():
                add(rec, PHASE1 / "data" / "assay_labels" / spec["file"], spec["sha256"])
        elif src.get("path", "").startswith("atlas:"):
            if not atlas_missing:
                add(rec, atlas / src["path"][len("atlas:"):], src["sha256"])
        elif key == "clinvar_vcf":                 # CLINVAR_VCF is this path (a test)
            if CLINVAR_VCF.exists():
                add(rec, CLINVAR_VCF, src["sha256"])
            else:
                absent.append("the ClinVar release")
        elif "path" in src and "sha256" in src:
            add(rec, PHASE1 / src["path"], src["sha256"])

    # the external genes
    rec = ev / "external" / "ddx3x_splice.manifest.json"
    add(rec, ev / "external" / "ddx3x_splice.parquet", json.loads(rec.read_text())["sha256"])
    rec = ev / "external" / "tp53_splice_scored_12nt.manifest.json"
    d = json.loads(rec.read_text())
    table = ev / "external" / "tp53_splice_scored_12nt.parquet"
    add(rec, table, d["parquet_sha256"])
    add(rec, table, d["sha256"], "content")
    add(rec, REPO / d["reference_alleles"]["file"], d["reference_alleles"]["sha256"])
    for part in ("subset", "events"):
        add(rec, REPO / d["inframe"][part], d["inframe"][part + "_sha256"])
    envs = [s["environment"] for s in d["score_files"].values()]
    envs.append(d["inframe"]["events_environment"])
    for env in envs:
        add(rec, REPO / env["full_record"], env["full_record_sha256"])
    for rel, sha in d["inputs_sha256"].items():
        add(rec, REPO / rel, sha)
    rec = PHASE1 / "data" / "external" / "tp53_splice_scored_v2.manifest.json"
    d = json.loads(rec.read_text())
    add(rec, rec.with_name("tp53_splice_scored_v2.parquet"), d["sha256"], "content")
    for model, sha in d["inputs"].items():
        add(rec, REPO / "data" / "tp53" / f"{model}_tp53.tsv", sha)

    # the frozen matrices the set joins to, and the rescored columns of the second
    frozen = PHASE1 / "data" / "frozen"
    for n in (1, 2):
        rec = frozen / f"manifest_v{n}.json"
        d = json.loads(rec.read_text())
        add(rec, PHASE1 / d["output_path"], d["sha256"], "content")
    add(rec, frozen / "frozen_matrix_v1.parquet", d["derived_from"]["sha256"], "content")
    for name, sha in d["provenance"].items():
        if isinstance(sha, str) and name.endswith(".tsv"):
            here = [x / name for x in (PHASE1 / "data" / "rescore", REPO / "data" / "rescore")
                    if (x / name).is_file()]
            add(rec, here[0] if here else REPO / "data" / "rescore" / name, sha)

    # the reference cache, and the supplementary tables with their sources
    rec = ev / "reference_cache" / "MANIFEST.json"
    for rel, meta in json.loads(rec.read_text())["files"].items():
        add(rec, ev / "reference_cache" / rel, meta["sha256"])
    rec = PHASE1 / "reports" / "evidence" / "supplement" / "supplement_tables_manifest.csv"
    with open(rec, newline="") as f:
        for row in csv.DictReader(f):
            add(rec, rec.with_name(row["File"]), row["File sha256"])
            files = [x for x in row["Source files"].split("; ") if x]
            shas = [x for x in row["Source sha256"].split("; ") if x]
            if len(files) != len(shas):                # S9: "written by <stage>"
                continue
            for src, sha in zip(files, shas):
                if sha == "n/a (code)":
                    continue
                if src.startswith("atlas:"):
                    if not atlas_missing:
                        add(rec, atlas / src[len("atlas:"):], sha, "source")
                    continue
                add(rec, REPO / src, sha, "source")
    return out, absent


def check() -> bool:
    """Every recorded checksum, then the test suite. No download, no stage run."""
    print(f"{'=' * 74}\nChecksums recorded in the manifests and provenance records\n"
          f"{'=' * 74}")
    rows, absent = recorded_checksums()
    bad = []
    for record, target, sha, how in rows:
        if not target.exists():
            bad.append(f"missing  {target}  (recorded in {record})")
            continue
        got = (_canonical_sha256(target) if how == "content" else
               _source_sha256(target) if how == "source" else _digest(target, "sha256"))
        if got != sha:
            bad.append(f"differs  {target}  (recorded in {record})")
    inside = {t for _, t, _, _ in rows if t.is_relative_to(REPO)}
    outside = {t for _, t, _, _ in rows} - inside
    print(f"  {len(rows)} checksums over {len(inside)} files in this repository"
          + (f" and {len(outside)} outside it" if outside else "")
          + f": {len(rows) - len(bad)} match, {len(bad)} do not")
    for b in bad:
        print(f"      {b}")
    for name in absent:
        print(f"  {name} is not here, so its checksums wait for the full run")
    print(f"\n{'=' * 74}\nThe test suite (pytest -q)\n{'=' * 74}", flush=True)
    env = dict(os.environ)
    atlas, missing = resolve_atlas()
    if not missing:                        # the tests read the atlas the checksums did
        env["EVID_ATLAS_REPO"] = str(atlas)
    tests = subprocess.run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider"],
                           cwd=REPO, env=env)
    ok = not bad and tests.returncode == 0
    print("\nCHECKED: every recorded checksum matches and every test passes." if ok else
          "\nCHECK FAILED: see above.")
    return ok


def run_stage(module: str, label: str, n: int, total: int, args: list) -> bool:
    """Run one stage. Returns False if it was skipped for a missing input.

    A missing input skips its own stage and lets the rest run. Exiting the whole
    script would mean an un-scored model column silently prevents the stages after
    it from regenerating at all, which is how outputs drift from the code that is
    supposed to produce them.
    """
    print(f"\n{'=' * 74}\n[{n}/{total}] {module}  --  {label}\n{'=' * 74}", flush=True)
    missing = [(rel, how) for rel, how in NEEDS.get(module, [])
               if not (REPO / rel).exists()]
    if missing:
        for rel, how in missing:
            print(f"  SKIPPED -- missing input {rel}\n    produce it with:  {how}")
        return False
    r = subprocess.run([sys.executable, "-m", module, *args], cwd=PHASE1)
    if r.returncode != 0:
        sys.exit(f"\nFAILED at stage {n}/{total} ({module}); exit code {r.returncode}")
    return True


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--from", dest="start", type=int, default=1)
    ap.add_argument("--fetch-inputs", action="store_true",
                    help="download the ClinVar release and the atlas release archive "
                         "when they are missing, and check their checksums")
    ap.add_argument("--check", action="store_true",
                    help="quick path, about five minutes, no data download: check every "
                         "checksum the manifests and provenance records carry, then run "
                         "the test suite; no stage runs")
    ap.add_argument("--verify", action="store_true",
                    help="compare every output with the committed one: REPRODUCED if all "
                         "match byte for byte, REPRODUCED NUMERICALLY if every printed "
                         "table and figure does and other files differ only in their last "
                         "digits; exit non-zero otherwise. Uses only an atlas that is the "
                         "release file for file")
    args = ap.parse_args()
    if args.check:
        sys.exit(0 if check() else 1)
    if args.verify and args.start != 1:
        sys.exit("--verify needs the full run (no --from)")
    if args.fetch_inputs:
        fetch_clinvar()
    atlas = choose_atlas(release_only=args.verify, fetch=args.fetch_inputs)
    os.environ["EVID_ATLAS_REPO"] = str(atlas)    # every stage reads this one atlas
    version = atlas_version(atlas)
    print(f"companion atlas: {atlas} (release {version}"
          + (", matches the release file for file)" if args.verify else ")"))
    if version != ATLAS_ARCHIVE["version"]:
        print(f"  WARNING: the analysis set was built from release "
              f"{ATLAS_ARCHIVE['version']}; outputs will differ")
    with tempfile.TemporaryDirectory(prefix="reproduce_before_") as keep:
        base = baseline(Path(keep)) if args.verify else None
        total = len(STAGES)
        skipped = []
        for i, (mod, label, extra) in enumerate(STAGES, 1):
            if i < args.start:
                print(f"[{i}/{total}] {mod} -- skipped (--from)")
                continue
            if not run_stage(mod, label, i, total, extra):
                skipped.append(f"{i}. {mod} -- {label}")
        print("\nTables: phase1/reports/evidence/")
        if skipped:
            print("\nStages skipped for a missing input -- their outputs are whatever "
                  "the last successful run left:")
            for s in skipped:
                print(f"  {s}")
            sys.exit(1)
        print("All stages ran.")
        if base is not None and not verify(base):
            sys.exit(1)
    print()


if __name__ == "__main__":
    main()
