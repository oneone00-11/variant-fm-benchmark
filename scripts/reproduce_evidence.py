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
 22  evid_delta              every published quantity with a counterpart (E8)

Outputs land in `phase1/reports/evidence/`.

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
                REPO / "docs" / "evidence-delta.md"]
NOT_OUTPUTS = [PHASE1 / "data" / "evidence" / "clinvar",
               PHASE1 / "data" / "evidence" / "companion_atlas"]
WALL_CLOCK_KEYS = {"built_utc", "scored_utc", "retrieved_utc", "run_at"}
ATLAS_BESIDE = REPO.parent / "functional-standard-atlas"
ATLAS_FETCHED = PHASE1 / "data" / "evidence" / "companion_atlas" / "atlas"
# Every atlas file the stages read. The three under results/ ship in the release
# archive and are not tracked in the atlas's git repository.
ATLAS_NEEDS = [
    "CITATION.cff",
    "src/atlas/__init__.py", "src/atlas/clinical_evidence.py", "src/atlas/evaluate.py",
    "src/atlas/manifest.py", "src/atlas/mapping.py", "src/atlas/predictor_resources.py",
    "results/score_matrix_atlas_v2.parquet", "results/alphagenome_v061_scores.parquet",
    "results/table1_atlas_composition.tsv",
]

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
    ("src.evid_delta",              "E8  old/new quantity list -> docs/evidence-delta.md", []),
]

# Stages that need an input this script does not produce, and what produces it.
NEEDS = {
    "src.evid_tier_logo": [
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
    both honour HTTP ranges). A minute without data counts as a stall. The caller
    checks the checksum, so a resumed file is never trusted on its own."""
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
            err = e
        except (OSError, http.client.HTTPException) as e:   # stalls, resets, a
            err = e                                           # connection cut short
        got = part.stat().st_size if part.exists() else 0
        print(f"  interrupted at {got / 2**20:.0f} MB ({err}); resuming, attempt "
              f"{attempt + 1} of {attempts}", flush=True)
        time.sleep(min(30, 2 * attempt))
    raise OSError(f"gave up on {url} after {attempts} attempts")


def fetch_clinvar() -> None:
    """The ClinVar release, from NCBI, checked against NCBI's md5 and the sha256 the
    analysis-set manifest records. An existing file is only checked."""
    if CLINVAR_VCF.exists():
        if _digest(CLINVAR_VCF, "sha256") != CLINVAR_SHA256:
            sys.exit(f"{CLINVAR_VCF.relative_to(REPO)} is not the release the "
                     "analysis set was built from (sha256 differs); remove it and rerun")
        print(f"ClinVar release present, sha256 checked: {CLINVAR_VCF.relative_to(REPO)}")
        return
    for url, md5_url in CLINVAR_SOURCES:
        try:
            md5 = (urllib.request.urlopen(md5_url, timeout=60).read().decode().split()[0]
                   if md5_url else None)
            print(f"downloading {url} (192 MB)", flush=True)
            _download(url, CLINVAR_VCF)
            break
        except (OSError, http.client.HTTPException) as e:
            print(f"  not available there ({e})")
    else:
        sys.exit("could not download the ClinVar release from any source")
    if ((md5 is not None and _digest(CLINVAR_VCF, "md5") != md5)
            or _digest(CLINVAR_VCF, "sha256") != CLINVAR_SHA256):
        CLINVAR_VCF.unlink()     # a rerun downloads it afresh
        sys.exit("the downloaded ClinVar file failed its checksums and was removed")
    print(("  md5 matches the source's, " if md5 else "  ")
          + "sha256 matches the analysis-set manifest")


def fetch_atlas() -> None:
    """The atlas release archive, from Zenodo, checked against Zenodo's md5 and
    unpacked to phase1/data/evidence/companion_atlas/atlas/."""
    a = ATLAS_ARCHIVE
    root = ATLAS_FETCHED.parent
    zpath = root / f"functional-standard-atlas-v{a['version']}.zip"
    if not (zpath.exists() and _digest(zpath, "md5") == a["md5"]):
        print(f"downloading the atlas release v{a['version']} (doi {a['doi']}, 50 MB)",
              flush=True)
        _download(a["url"], zpath)
    if _digest(zpath, "md5") != a["md5"]:
        zpath.unlink()
        sys.exit("the downloaded atlas archive failed its md5 check and was removed")
    with zipfile.ZipFile(zpath) as z:
        for name in z.namelist():            # nothing may land outside root
            if not (root / name).resolve().is_relative_to(root.resolve()):
                sys.exit(f"unexpected path in the atlas archive: {name}")
        z.extractall(root)
    print(f"  md5 matches Zenodo's; unpacked to {ATLAS_FETCHED.relative_to(REPO)}")


def _missing(atlas: Path) -> list[str]:
    return [f for f in ATLAS_NEEDS if not (atlas / f).is_file()]


def resolve_atlas() -> tuple[Path, list[str]]:
    """The atlas every stage will read, and the files it lacks: EVID_ATLAS_REPO if
    set; otherwise the checkout beside this repository if it is complete, else the
    fetched release archive."""
    if os.environ.get("EVID_ATLAS_REPO"):
        path = Path(os.environ["EVID_ATLAS_REPO"])
        return path, _missing(path)
    for path in (ATLAS_BESIDE, ATLAS_FETCHED):
        if path.is_dir() and not _missing(path):
            return path, []
    path = ATLAS_BESIDE if ATLAS_BESIDE.is_dir() else ATLAS_FETCHED
    return path, _missing(path)


def atlas_version(atlas: Path) -> str:
    for line in (atlas / "CITATION.cff").read_text().splitlines():
        if line.startswith("version:"):
            return line.split(":", 1)[1].strip().strip('"').strip("'")
    return "unknown"


def _outputs() -> list[Path]:
    files = []
    for root in OUTPUT_ROOTS:
        for p in ([root] if root.is_file() else sorted(root.rglob("*"))):
            if (p.is_file() and not p.name.startswith(".") and "__pycache__" not in p.parts
                    and not any(p.is_relative_to(x) for x in NOT_OUTPUTS)):
                files.append(p)
    return files


def _without_wall_clock(path: Path):
    def strip(x):
        if isinstance(x, dict):
            return {k: strip(v) for k, v in x.items() if k not in WALL_CLOCK_KEYS}
        if isinstance(x, list):
            return [strip(v) for v in x]
        return x
    return strip(json.loads(path.read_text()))


# How --verify reads a file that is not byte-identical. Tables are compared cell by
# cell, floats to a relative tolerance and integers and strings exactly, so that a
# report from another platform says whether a difference is in the numbers; figures
# are named as such. Neither refinement changes the verdict: REPRODUCED means every
# file is byte-identical, the build time in a manifest aside.
TABLE_DELIMITERS = {".csv": ",", ".tsv": "\t"}
FIGURE_SUFFIXES = {".png", ".pdf", ".svg", ".tif", ".tiff"}
RTOL = 1e-9
REPORT_NAME = "reproduce_report.txt"      # at the repository root, gitignored


def snapshot() -> dict:
    """The outputs as the checkout carries them, before any stage runs. Tables and
    manifests are also copied aside, so a changed one can be compared cell by cell."""
    keep = Path(tempfile.mkdtemp(prefix="reproduce_before_"))
    snap = {}
    for p in _outputs():
        entry = {"sha256": _digest(p, "sha256"), "mtime": p.stat().st_mtime}
        if p.suffix in TABLE_DELIMITERS or p.suffix in (".json", ".parquet"):
            entry["copy"] = keep / p.relative_to(REPO)
            entry["copy"].parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(p, entry["copy"])
        snap[p] = entry
    return snap


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
    a, b = _cell(a), _cell(b)
    if isinstance(a, bool) or isinstance(b, bool) or not (
            isinstance(a, (int, float)) and isinstance(b, (int, float))):
        return a == b
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
        if after.suffix == ".parquet":
            import pandas as pd
            x, y = pd.read_parquet(before), pd.read_parquet(after)
            if list(x.columns) != list(y.columns) or x.shape != y.shape:
                return False
            return all(_equal(u, v) for c in x.columns
                       for u, v in zip(x[c].tolist(), y[c].tolist()))
    except (OSError, ValueError, csv.Error):
        return False
    return False


def _status(p: Path, b: dict) -> str:
    if not p.exists():
        return "missing"
    if _digest(p, "sha256") == b["sha256"]:
        return "identical"
    if p.suffix == ".json" and _without_wall_clock(p) == _without_wall_clock(b["copy"]):
        return "identical (build time aside)"
    if p.suffix in FIGURE_SUFFIXES:
        return "differs (figure)"
    if "copy" in b and _numerically_equal(b["copy"], p):
        return "differs (numerically equal)"
    return "differs"


def verify(before: dict) -> bool:
    """Every output the checkout carried, compared with what this run wrote; one line
    per file in reproduce_report.txt, and a summary here."""
    rows = [(_status(p, b), p.relative_to(REPO), p.exists() and p.stat().st_mtime == b["mtime"])
            for p, b in before.items()]
    new = [p.relative_to(REPO) for p in _outputs() if p not in before]
    rows += [("new", r, False) for r in new]
    counts = {}
    for status, _, _ in rows:
        counts[status] = counts.get(status, 0) + 1
    ok = all(st in ("identical", "identical (build time aside)") for st, _, _ in rows)
    stale = sum(1 for _, _, untouched in rows if untouched)
    verdict = ("REPRODUCED: every output matches the committed one." if ok else
               "NOT REPRODUCED: every file not marked identical")
    report = REPO / REPORT_NAME
    with open(report, "w") as f:
        f.write(f"# reproduce_evidence.py --verify\n# {verdict}\n")
        f.write("# " + "; ".join(f"{k}: {v}" for k, v in sorted(counts.items())) + "\n")
        f.write(f"# {stale} files under the output folders are inputs the run only reads\n")
        for status, rel, _ in sorted(rows, key=lambda r: str(r[1])):
            f.write(f"{status}\t{rel}\n")
    print(f"\n{'=' * 74}\nVerification against the outputs this checkout carried\n"
          f"{'=' * 74}")
    print(f"  {len(before)} files under the output folders: {len(before) - stale} "
          f"rewritten by this run, {stale} inputs it only reads")
    for status in ("identical", "identical (build time aside)", "differs (numerically equal)",
                   "differs (figure)", "differs", "missing", "new"):
        if counts.get(status):
            print(f"  {status + ':':30s} {counts[status]}")
            if status != "identical":
                for st, rel, _ in rows:
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
    if str(PHASE1) not in sys.path:
        sys.path.insert(0, str(PHASE1))
    import pandas as pd
    from src.phase1_build_frozen_matrix_v2 import canonical_sha256
    return canonical_sha256(pd.read_parquet(path))


def recorded_checksums() -> list[tuple[str, Path, str, str]]:
    """(record, file, expected sha256, "file" or "content") for every checksum a
    tracked manifest or provenance record carries about a file in this repository,
    and about the ClinVar release and the atlas files when they are present. A
    "content" hash is taken over the table's rows, not the file's bytes. Hashes of
    files that live only where the model scores were made (reference FASTAs, scorer
    modules) are left out: nothing here can check them."""
    ev, out = PHASE1 / "data" / "evidence", []

    def add(record: Path, target: Path, sha: str, how: str = "file") -> None:
        out.append((str(record.relative_to(REPO)), target, sha, how))

    # provenance sidecars: <output>.provenance.json beside the file it describes
    for rec in sorted((PHASE1 / "data").rglob("*.provenance.json")):
        d = json.loads(rec.read_text())
        stem = rec.name[: -len(".provenance.json")]
        target = rec.with_name(stem)
        if not target.exists():                    # atlas_columns_v2 -> .tsv
            target = next(iter(sorted(p for p in rec.parent.glob(stem + ".*")
                                      if p != rec)), target)
        for key in ("output_sha256", "sha256"):
            if isinstance(d.get(key), str):
                add(rec, target, d[key])
        if isinstance(d.get("input_sha256"), str) and isinstance(d.get("input"), str):
            add(rec, REPO / d["input"], d["input_sha256"])

    # the analysis set and the sources it was built from
    rec = ev / "analysis_set_v1.manifest.json"
    d = json.loads(rec.read_text())
    add(rec, ev / "analysis_set_v1.parquet", d["sha256"])
    for key, src in d["sources"].items():
        if key == "assay_labels":
            for spec in src.values():
                add(rec, PHASE1 / "data" / "assay_labels" / spec["file"], spec["sha256"])
        elif src.get("path", "").startswith("atlas:"):
            atlas, missing = resolve_atlas()
            if not missing:
                add(rec, atlas / src["path"][len("atlas:"):], src["sha256"])
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
    rec = PHASE1 / "data" / "external" / "tp53_splice_scored_v2.manifest.json"
    add(rec, rec.with_name("tp53_splice_scored_v2.parquet"),
        json.loads(rec.read_text())["sha256"], "content")

    # the frozen matrices the set joins to
    for n in (1, 2):
        rec = PHASE1 / "data" / "frozen" / f"manifest_v{n}.json"
        d = json.loads(rec.read_text())
        add(rec, PHASE1 / d["output_path"], d["sha256"], "content")

    # the reference cache and the supplementary tables
    rec = ev / "reference_cache" / "MANIFEST.json"
    for rel, meta in json.loads(rec.read_text())["files"].items():
        add(rec, ev / "reference_cache" / rel, meta["sha256"])
    rec = PHASE1 / "reports" / "evidence" / "supplement" / "supplement_tables_manifest.csv"
    with open(rec, newline="") as f:
        for row in csv.DictReader(f):
            add(rec, rec.with_name(row["File"]), row["File sha256"])
    return out


def check() -> bool:
    """Every recorded checksum, then the test suite. No download, no stage run."""
    print(f"{'=' * 74}\nChecksums recorded in the manifests and provenance records\n"
          f"{'=' * 74}")
    rows = recorded_checksums()
    bad = []
    for record, target, sha, how in rows:
        if not target.exists():
            bad.append(f"missing  {target.relative_to(REPO)}  (recorded in {record})")
            continue
        got = _canonical_sha256(target) if how == "content" else _digest(target, "sha256")
        if got != sha:
            bad.append(f"differs  {target.relative_to(REPO)}  (recorded in {record})")
    n_files = len({t for _, t, _, _ in rows})
    print(f"  {len(rows)} checksums over {n_files} files: "
          f"{len(rows) - len(bad)} match, {len(bad)} do not")
    for b in bad:
        print(f"      {b}")
    if not CLINVAR_VCF.exists():
        print("  the ClinVar release and the atlas are not here, so their checksums wait "
              "for the full run")
    print(f"\n{'=' * 74}\nThe test suite (pytest -q)\n{'=' * 74}", flush=True)
    tests = subprocess.run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider"],
                           cwd=REPO)
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
                    help="quick path, about five minutes, no download: check every "
                         "checksum the manifests and provenance records carry, then run "
                         "the test suite; no stage runs")
    ap.add_argument("--verify", action="store_true",
                    help="compare every output with the one the checkout carried before "
                         "the run, and exit non-zero unless all match")
    args = ap.parse_args()
    if args.check:
        sys.exit(0 if check() else 1)
    if args.verify and args.start != 1:
        sys.exit("--verify needs the full run (no --from)")
    if args.fetch_inputs:
        fetch_clinvar()
    atlas, missing = resolve_atlas()
    if missing and args.fetch_inputs and not os.environ.get("EVID_ATLAS_REPO"):
        fetch_atlas()
        atlas, missing = resolve_atlas()
    if missing:
        sys.exit(f"The companion atlas at {atlas} lacks {', '.join(missing)}.\n"
                 "Its results/ files ship in the release archive, not in its git "
                 "repository. Rerun with --fetch-inputs to download the archive "
                 f"(doi {ATLAS_ARCHIVE['doi']}), or set EVID_ATLAS_REPO to an unpacked "
                 "copy of it.")
    os.environ["EVID_ATLAS_REPO"] = str(atlas)    # every stage reads this one atlas
    version = atlas_version(atlas)
    print(f"companion atlas: {atlas} (release {version})")
    if version != ATLAS_ARCHIVE["version"]:
        print(f"  WARNING: the analysis set was built from release "
              f"{ATLAS_ARCHIVE['version']}; outputs will differ")
    before = snapshot() if args.verify else None
    total = len(STAGES)
    skipped = []
    for i, (mod, label, extra) in enumerate(STAGES, 1):
        if i < args.start:
            print(f"[{i}/{total}] {mod} -- skipped (--from)")
            continue
        if not run_stage(mod, label, i, total, extra):
            skipped.append(f"{i}. {mod} -- {label}")
    print(f"\nTables: phase1/reports/evidence/")
    if skipped:
        print("\nStages skipped for a missing input -- their outputs are whatever "
              "the last successful run left:")
        for s in skipped:
            print(f"  {s}")
        sys.exit(1)
    print("All stages ran.")
    if before is not None and not verify(before):
        sys.exit(1)
    print()


if __name__ == "__main__":
    main()
