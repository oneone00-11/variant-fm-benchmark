"""The seven training-gene label files built from their recorded sources (stage E16d):
what scripts/stage_E16d/build_labels.py must reproduce, and what relabelling PALB2 and
BARD1 may change.

The generator writes into phase1/data/assay_labels/_generated/. Five of its outputs
must be the files in use, byte for byte. PALB2's and BARD1's must be the investigator-
provided ranges of the cited a-2, variant by variant, and must differ from the files
that entered the repository on 2026-07-08 on exactly the variants the label source
audit (stage E16c) listed: 359 and 13. These tests read the July files from
superseded/ once the switch has moved them there, and from their own place before.
"""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import io
import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

PHASE1 = Path(__file__).resolve().parents[1]
REPO = PHASE1.parent
LABELS = PHASE1 / "data/assay_labels"
GENERATED = LABELS / "_generated"
SUPERSEDED = LABELS / "superseded"
PROVENANCE = LABELS / "PROVENANCE.json"
AUDIT = PHASE1 / "data/evidence/label_audit"

if str(PHASE1) not in sys.path:
    sys.path.insert(0, str(PHASE1))
from src import config as C              # noqa: E402
from src import evid_newgenes as N       # noqa: E402
from src.evid_build_set import _norm_hgvs  # noqa: E402

UNCHANGED = ["BRCA1", "BRCA2", "RAD51C", "VHL", "BAP1"]
RELABELLED = ["PALB2", "BARD1"]
# The label source audit as committed (E16d commit 1): label_audit_diff.csv there lists
# the variants on which R1 differs from the July files; label_audit_summary.csv counts
# them (359 and 13).
AUDIT_COMMIT = "76ed7c22d28a0a38503038d248563715e3fcac32"
AUDIT_COUNTS = {"PALB2": 359, "BARD1": 13}
JULY_COMMIT = "10d94bb"                  # where all seven label files entered (2026-07-08)
SWITCH_SUBJECT = "Use the archived classification for PALB2 and BARD1"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _module(name: str):
    spec = importlib.util.spec_from_file_location(name, REPO / "scripts" / "stage_E16d" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


B = _module("build_labels")


def _classes(text: str, name: str, gene: str) -> dict:
    spec = C.ASSAY_LABELS[gene]
    rows = csv.DictReader(io.StringIO(text), delimiter="\t" if name.endswith(".tsv") else ",")
    return {_norm_hgvs(r[spec["key"]]): r[spec["class_col"]] for r in rows}


def _july(gene: str) -> Path:
    name = C.ASSAY_LABELS[gene]["file"]
    return SUPERSEDED / name if (SUPERSEDED / name).exists() else LABELS / name


def _git(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, text=True)


# --------------------------------------------------------------------------
# the record and the generator
# --------------------------------------------------------------------------
@pytest.mark.parametrize("gene", list(C.ASSAY_LABELS))
def test_the_provenance_record_matches_the_generated_files_and_their_sources(gene):
    doc = json.loads(PROVENANCE.read_text())
    rec = doc["genes"][gene]
    name = C.ASSAY_LABELS[gene]["file"]
    gen = REPO / rec["generated"]["path"]
    assert gen == GENERATED / name and _sha(gen) == rec["generated"]["sha256"]
    assert rec["generated"]["rows"] == len(gen.read_text().splitlines()) - 1
    assert rec["label_file"] == f"phase1/data/assay_labels/{name}"
    assert rec["in_use"] == ((LABELS / name).read_bytes() == gen.read_bytes())
    for f in rec["source"]["files"]:
        if f["committed"] or (REPO / f["path"]).exists():
            assert _sha(REPO / f["path"]) == f["sha256"], f["path"]
    assert rec["source"]["type"] in ("MaveDB score set version", "paper supplementary table")
    assert rec["source"]["rule"].startswith(B.RULES[gene]["rule"])
    assert doc["generator_sha256"] == _sha(REPO / "scripts/stage_E16d/build_labels.py")


@pytest.mark.parametrize("gene", list(C.ASSAY_LABELS))
def test_the_generator_rebuilds_its_committed_output(gene):
    """Skipped, with the reason, for a gene whose source is not here (the BRCA2 table
    sits in the label audit's cache, which --fetch-inputs fills)."""
    try:
        blob, _ = B.build(gene)
    except B.Skip as e:
        print(f"skipped {gene}: {e}")
        pytest.skip(f"skipped {gene}: {e}")
    assert blob == (GENERATED / C.ASSAY_LABELS[gene]["file"]).read_bytes()


@pytest.mark.parametrize("gene", UNCHANGED)
def test_the_five_unchanged_genes_are_rebuilt_byte_for_byte(gene):
    name = C.ASSAY_LABELS[gene]["file"]
    assert (GENERATED / name).read_bytes() == (LABELS / name).read_bytes()


def test_every_wording_table_maps_into_the_configuration():
    for gene, words in B.WORDING.items():
        assert set(words.values()) <= set(C.ASSAY_LABELS[gene]["map"]), gene
        got = set(_classes((GENERATED / C.ASSAY_LABELS[gene]["file"]).read_text(),
                           C.ASSAY_LABELS[gene]["file"], gene).values())
        assert got <= set(words.values()), gene


# --------------------------------------------------------------------------
# PALB2 and BARD1
# --------------------------------------------------------------------------
def _r1(gene: str) -> dict:
    """The cited a-2's investigator-provided ranges applied here, independently of the
    generator, in the generator's wording."""
    urn = B.RULES[gene]["version"]
    man = json.loads((AUDIT / f"{gene.lower()}_raw/manifest.json").read_text())
    vid = urn.split(":")[-1]
    path = {Path(e["path"]).name: REPO / e["path"] for e in man["files"]}
    meta = json.loads(path[f"scoreset_meta_{vid}.json"].read_text())
    cal = [c for c in meta["scoreCalibrations"] if c.get("investigatorProvided") is True]
    assert len(cal) == 1
    classes = cal[0]["functionalClassifications"]
    out = {}
    with open(path[f"urn_mavedb_{vid}.csv"], newline="") as f:
        for r in csv.DictReader(f):
            hit = [c["label"] for c in classes if N._matches(float(r["score"]), c)]
            assert len(hit) == 1, (r["hgvs_nt"], r["score"], hit)   # none, or a shared point
            out[_norm_hgvs(r["hgvs_nt"])] = B.WORDING[gene][hit[0]]
    return out


@pytest.mark.parametrize("gene", RELABELLED)
def test_palb2_and_bard1_are_the_cited_ranges_variant_by_variant(gene):
    name = C.ASSAY_LABELS[gene]["file"]
    assert _classes((GENERATED / name).read_text(), name, gene) == _r1(gene)


@pytest.mark.parametrize("gene", RELABELLED)
def test_palb2_and_bard1_differ_from_the_july_files_exactly_where_the_audit_said(gene):
    name = C.ASSAY_LABELS[gene]["file"]
    gen = _classes((GENERATED / name).read_text(), name, gene)
    july = _classes(_july(gene).read_text(), name, gene)
    assert set(gen) == set(july)
    differ = {k for k in gen if gen[k] != july[k]}
    assert len(differ) == AUDIT_COUNTS[gene]
    shown = _git("show", f"{AUDIT_COMMIT}:phase1/reports/evidence/label_audit_diff.csv")
    if shown.returncode:
        pytest.skip(f"the audit's commit {AUDIT_COMMIT[:7]} is not in this clone; the count is checked")
    listed = {_norm_hgvs(r["hgvs_nt"]) for r in csv.DictReader(io.StringIO(shown.stdout))
              if r["gene"] == gene and r["rule"] == "R1" and r["version"] == B.RULES[gene]["version"]}
    assert differ == listed


def test_no_score_of_a_ranged_source_is_unassigned_or_on_a_shared_end_point():
    """R1 for BRCA1, BARD1 and PALB2: every score of the cited table lies in exactly one
    range (_r1 asserts it row by row)."""
    for gene in ("BRCA1", "BARD1", "PALB2"):
        rows = len((GENERATED / C.ASSAY_LABELS[gene]["file"]).read_text().splitlines()) - 1
        assert len(_r1(gene)) == rows, gene


# --------------------------------------------------------------------------
# the switch: data only, and the July files kept
# --------------------------------------------------------------------------
def test_the_switch_commit_changes_no_code():
    """E16d.0 rule 2, on the commit that replaces the two label files, once it exists
    and its parent is in the clone."""
    if not (REPO / ".git").exists():
        pytest.skip("a release archive has no .git")
    log = _git("log", "--format=%H\t%s").stdout.splitlines()
    hits = [h for h, s in (l.split("\t", 1) for l in log) if s.startswith(SWITCH_SUBJECT)]
    if not hits:
        pytest.skip("the switch commit has not been made yet")
    if _git("rev-parse", "--verify", "--quiet", f"{hits[0]}^").returncode:
        pytest.skip("the switch commit's parent is not in this clone (shallow)")
    check = _module("check_switch_commit")
    assert not [p for p in check.changed(hits[0]) if check.is_code(p)]


def test_the_superseded_files_are_the_july_files():
    if not SUPERSEDED.exists():
        pytest.skip("the switch has not been made yet")
    assert (SUPERSEDED / "README.md").is_file()
    for gene in RELABELLED:
        name = C.ASSAY_LABELS[gene]["file"]
        shown = subprocess.run(["git", "-C", str(REPO), "show",
                                f"{JULY_COMMIT}:phase1/data/assay_labels/{name}"], capture_output=True)
        if shown.returncode:
            pytest.skip(f"commit {JULY_COMMIT} is not in this clone")
        assert (SUPERSEDED / name).read_bytes() == shown.stdout, name


def test_only_the_comparison_report_is_exempt_from_being_rewritten():
    """--verify expects every tracked output to be rewritten by a full run; the one file
    the E16d comparison writes once is the only exception, and is not an input."""
    spec = importlib.util.spec_from_file_location("reproduce_evidence",
                                                  REPO / "scripts" / "reproduce_evidence.py")
    E = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(E)
    assert E.COMPARISON_REPORTS == ["phase1/reports/evidence/label_change_effect.csv"]
    assert set(E.COMPARISON_REPORTS).isdisjoint(E._run_inputs())
