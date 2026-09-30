"""Guardrails on reproducing the evidence-strength outputs from public sources.

A clean-room run (a fresh clone, the atlas release archive from Zenodo, a Python
environment built from requirements-evidence.lock.txt, ClinVar from NCBI) found three
committed files that depended on the machine that wrote them: the analysis-set
manifest read the atlas commit from `git` and recorded "unknown" for the archive, the
TP53 manifest checked its REF bases against a 706 MB local FASTA that no archive
carries, and the supplementary-table manifest hashes both. These tests hold the fixes,
and the one command built on them: what it fetches is what the manifests record, it
finds the atlas by one documented rule and accepts only the release when it verifies,
its downloads survive the network, and its two verdicts (REPRODUCED, CHECKED) cannot
be reached by a run that should fail them.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd
import pytest

PHASE1 = Path(__file__).resolve().parents[1]
REPO = PHASE1.parent
DATA = PHASE1 / "data/evidence"
SET_MANIFEST = DATA / "analysis_set_v1.manifest.json"
TP53_MANIFEST = DATA / "external/tp53_splice_scored_12nt.manifest.json"
TP53_VARIANTS = DATA / "external/tp53_12nt_variants.parquet"
REF_CACHE = DATA / "reference_cache"

if str(PHASE1) not in sys.path:
    sys.path.insert(0, str(PHASE1))
from src import evid_build_set as B  # noqa: E402
from src import evid_common as K     # noqa: E402


def _entry_point():
    spec = importlib.util.spec_from_file_location(
        "reproduce_evidence", REPO / "scripts" / "reproduce_evidence.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


E = _entry_point()


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _no_local_path(obj) -> bool:
    text = json.dumps(obj)
    return not any(p in text for p in ("/Users/", "/home/", "/private/", "/tmp/"))


# --------------------------------------------------------------------------
# the inputs the entry point fetches are the ones the analysis set was built from
# --------------------------------------------------------------------------
def test_the_fetched_clinvar_release_is_the_one_the_manifest_records():
    sources = json.loads(SET_MANIFEST.read_text())["sources"]
    assert E.CLINVAR_SHA256 == sources["clinvar_vcf"]["sha256"]
    assert E.CLINVAR_VCF == PHASE1 / sources["clinvar_vcf"]["path"]
    assert all("/" + E.CLINVAR_VCF.name in u for u in E.CLINVAR_URLS)
    assert E.CLINVAR_URLS[0].startswith("https://ftp.ncbi.nlm.nih.gov/")   # NCBI first
    assert "zenodo.org" in E.CLINVAR_URLS[-1]                               # the mirror last


def test_the_fetched_atlas_release_is_the_one_the_manifest_records():
    sources = json.loads(SET_MANIFEST.read_text())["sources"]
    version = E.ATLAS_ARCHIVE["version"]
    assert sources["atlas_matrix"]["atlas_git_head"] == K.ATLAS_RELEASE_COMMIT[version]
    pins = E._pins()
    for key in ("atlas_matrix", "atlas_alphagenome_v061"):
        rel = sources[key]["path"][len("atlas:"):]
        assert rel in K.ATLAS_NEEDS and pins[rel] == sources[key]["sha256"]


def test_every_atlas_file_the_stages_need_is_pinned_to_the_release():
    for p in (B.ATLAS_MATRIX, B.ATLAS_V061):
        assert str(p.relative_to(B.ATLAS_REPO)) in K.ATLAS_NEEDS
    assert set(K.ATLAS_NEEDS) <= set(E._pins())


def test_neither_manifest_carries_a_local_path():
    assert _no_local_path(json.loads(SET_MANIFEST.read_text())["sources"])
    assert _no_local_path(json.loads(TP53_MANIFEST.read_text())["reference_alleles"])


# --------------------------------------------------------------------------
# where the atlas is found
# --------------------------------------------------------------------------
def _atlas(root: Path, complete: bool = True) -> Path:
    for f in K.ATLAS_NEEDS if complete else K.ATLAS_NEEDS[:3]:
        (root / f).parent.mkdir(parents=True, exist_ok=True)
        (root / f).write_text('version: "2.5.0"\n' if f == "CITATION.cff" else "")
    return root


def test_the_atlas_is_found_in_the_documented_order(tmp_path, monkeypatch):
    beside, fetched = tmp_path / "beside", tmp_path / "fetched"
    monkeypatch.delenv("EVID_ATLAS_REPO", raising=False)
    find = lambda: K.resolve_atlas(beside, fetched)       # noqa: E731

    path, missing = find()                               # nothing anywhere
    assert missing == K.ATLAS_NEEDS

    _atlas(beside, complete=False)                       # a bare clone beside
    path, missing = find()
    assert path == beside and "results/score_matrix_atlas_v2.parquet" in missing

    _atlas(fetched)                                      # the fetched archive wins
    assert find() == (fetched, [])

    _atlas(beside)                                       # a complete checkout beside
    assert find() == (beside, [])

    other = _atlas(tmp_path / "other", complete=False)   # an explicit path is never
    monkeypatch.setenv("EVID_ATLAS_REPO", str(other))    # silently replaced
    path, missing = find()
    assert path == other and missing

    monkeypatch.chdir(tmp_path)                          # and a relative one is made
    monkeypatch.setenv("EVID_ATLAS_REPO", "beside")      # absolute, so every stage
    assert find() == (beside.resolve(), [])              # reads the same place


def test_verify_accepts_only_the_release_file_for_file(tmp_path, monkeypatch):
    copy = _atlas(tmp_path / "not-the-release")
    monkeypatch.setenv("EVID_ATLAS_REPO", str(copy))
    assert E.atlas_release_mismatch(copy)
    with pytest.raises(SystemExit, match="not the atlas release"):
        E.choose_atlas(release_only=True, fetch=True)
    assert E.choose_atlas(release_only=False, fetch=False) == copy.resolve()


def test_the_atlas_here_is_the_release_when_there_is_one():
    atlas, missing = K.resolve_atlas()
    if missing:
        pytest.skip("no atlas here")
    bad = E.atlas_release_mismatch(atlas)
    assert not bad or all("checkout at" in b for b in bad), bad


# --------------------------------------------------------------------------
# the atlas commit is the same from a checkout and from the release archive
# --------------------------------------------------------------------------
_GIT_ENV = {**os.environ, "GIT_CONFIG_GLOBAL": os.devnull, "GIT_CONFIG_NOSYSTEM": "1",
            "GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@t",
            "GIT_COMMITTER_NAME": "t", "GIT_COMMITTER_EMAIL": "t@t"}


def _git(*args, cwd):
    """git without the user's configuration: no signing, no hooks."""
    return subprocess.run(["git", *args], cwd=cwd, check=True, capture_output=True,
                          text=True, env=_GIT_ENV).stdout.strip()


def test_an_archive_records_the_commit_its_release_was_cut_from(tmp_path):
    outer = tmp_path / "vfm"
    outer.mkdir()
    _git("init", "-q", cwd=outer)                       # the archive is unpacked inside
    (outer / "x").write_text("x")                       # another checkout, whose HEAD
    _git("add", "x", cwd=outer)                         # must not be taken for the
    _git("commit", "-q", "-m", "x", cwd=outer)          # atlas's
    archive = _atlas(outer / "phase1/data/evidence/companion_atlas/atlas")
    assert B.git_head(archive) == "unknown"
    assert B.atlas_commit(archive) == K.ATLAS_RELEASE_COMMIT["2.5.0"]


def test_a_checkout_records_its_own_head(tmp_path):
    atlas = _atlas(tmp_path / "functional-standard-atlas")
    _git("init", "-q", cwd=atlas)
    _git("add", "CITATION.cff", cwd=atlas)
    _git("commit", "-q", "-m", "c", cwd=atlas)
    head = _git("rev-parse", "HEAD", cwd=atlas)
    assert B.atlas_commit(atlas) == head != K.ATLAS_RELEASE_COMMIT["2.5.0"]


# --------------------------------------------------------------------------
# TP53's REF bases are checked against a sequence the repository carries
# --------------------------------------------------------------------------
def test_the_tp53_reference_check_ran_on_every_variant():
    ref = json.loads(TP53_MANIFEST.read_text())["reference_alleles"]
    assert ref["checked"] is True and ref["n"] == 288 and ref["mismatches"] == 0


def test_the_tp53_reference_is_tracked_registered_and_agrees_with_every_ref():
    ref = json.loads(TP53_MANIFEST.read_text())["reference_alleles"]
    path = REPO / ref["file"]
    rel = str(path.relative_to(REF_CACHE))
    registered = json.loads((REF_CACHE / "MANIFEST.json").read_text())["files"][rel]
    assert _sha(path) == registered["sha256"] == ref["sha256"]
    if (REPO / ".git").exists():                        # a release archive has no .git
        tracked = subprocess.run(["git", "ls-files", "--error-unmatch", str(path)],
                                 cwd=REPO, capture_output=True)
        assert tracked.returncode == 0, f"{rel} is not tracked"

    _, _, chrom, start, end = path.stem.split("_")      # seq_GRCh38_chr17_<s>_<e>
    start, end, seq = int(start), int(end), path.read_text().strip().upper()
    assert len(seq) == end - start + 1
    v = pd.read_parquet(TP53_VARIANTS)
    assert len(v) == 288 and (v["chrom"].astype(str) == chrom.removeprefix("chr")).all()
    got = [seq[p - start] for p in v["pos"].astype(int)]
    assert got == list(v["ref"])


def test_atlas_sources_are_named_the_same_wherever_the_atlas_sits(monkeypatch):
    from src import evid_supp_tables as S
    inside = PHASE1 / "data/evidence/companion_atlas/atlas"      # --fetch-inputs
    monkeypatch.setattr(S, "ATLAS_REPO", inside)
    assert S._display(inside / "src/atlas/predictor_resources.py") == \
        "atlas:src/atlas/predictor_resources.py"
    assert S._display(PHASE1 / "src/evid_common.py") == "phase1/src/evid_common.py"


# --------------------------------------------------------------------------
# --verify: the verdict the one command ends with
# --------------------------------------------------------------------------
def _outputs_tree(tmp_path, monkeypatch, git: bool = False):
    """A miniature repository: one output folder, one input the stages only read."""
    repo = tmp_path / "repo"
    out = repo / "reports"
    out.mkdir(parents=True)
    (out / "t.csv").write_text("a,b\n1,2\n")
    (out / "m.manifest.json").write_text('{\n  "n": 1,\n  "built_utc": "2026-01-01"\n}\n')
    (out / "input.parquet").write_bytes(b"read, never written")
    inputs = repo / "inputs.txt"
    inputs.write_text("# inputs\nreports/input.parquet\n")
    monkeypatch.setattr(E, "REPO", repo)
    monkeypatch.setattr(E, "OUTPUT_ROOTS", [out])
    monkeypatch.setattr(E, "NOT_OUTPUTS", [])
    monkeypatch.setattr(E, "RUN_INPUTS", inputs)
    monkeypatch.setattr(E, "BASELINE", repo / ".reproduce_baseline")
    if git:
        _git("init", "-q", cwd=repo)
        _git("add", "-A", cwd=repo)
        _git("commit", "-q", "-m", "outputs", cwd=repo)
    return repo, out


def _rerun(out: Path, csv="a,b\n1,2\n",
           manifest='{\n  "n": 1,\n  "built_utc": "2026-09-29"\n}\n', skip=()):
    """What a run writes: every output rewritten (mtime after the baseline's start)."""
    time.sleep(0.01)
    if "t.csv" not in skip:
        (out / "t.csv").write_text(csv)
    (out / "m.manifest.json").write_text(manifest)


def _report(repo: Path) -> str:
    return (repo / "reproduce_report.txt").read_text()


def _verdict(repo: Path) -> str:
    return _report(repo).splitlines()[1]


def test_verify_passes_only_when_every_output_matches(tmp_path, monkeypatch, capsys):
    repo, out = _outputs_tree(tmp_path, monkeypatch)
    base = E.baseline(tmp_path / "keep")
    _rerun(out)                                             # only the build time moved
    assert E.verify(base)
    assert "REPRODUCED: every output matches" in capsys.readouterr().out
    assert "identical (build time aside)\treports/m.manifest.json" in _report(repo)

    base = E.baseline(tmp_path / "keep")                    # same number, other bytes
    _rerun(out, csv="a,b\n1.0000000000001,2\n")             # in an intermediate table:
    assert E.verify(base)                                   # the second tier, named
    assert "differs (numerically equal)\treports/t.csv" in _report(repo)
    assert _verdict(repo).startswith("# REPRODUCED NUMERICALLY")

    base = E.baseline(tmp_path / "keep")                    # a real difference fails
    _rerun(out, csv="a,b\n1.001,2\n")                       # both tiers
    assert not E.verify(base)
    assert _verdict(repo).startswith("# NOT REPRODUCED")


def test_verify_holds_manifests_to_their_bytes_bar_the_build_time(tmp_path, monkeypatch):
    repo, out = _outputs_tree(tmp_path, monkeypatch)
    base = E.baseline(tmp_path / "keep")
    _rerun(out, manifest='{\n  "n": 1.0,\n  "built_utc": "2026-09-29"\n}\n')   # 1 -> 1.0
    assert E.verify(base)                                   # not byte for byte
    assert "differs (numerically equal)\treports/m.manifest.json" in _report(repo)
    assert _verdict(repo).startswith("# REPRODUCED NUMERICALLY")


def test_a_printed_table_must_match_byte_for_byte(tmp_path, monkeypatch):
    repo, out = _outputs_tree(tmp_path, monkeypatch)
    printed = out / "tables"
    printed.mkdir()
    (printed / "table1.csv").write_text("x\n0.5\n")
    monkeypatch.setattr(E, "PRINTED", [printed])
    base = E.baseline(tmp_path / "keep")
    _rerun(out)
    (printed / "table1.csv").write_text("x\n0.50000000000001\n")   # numerically equal
    assert not E.verify(base)
    assert _verdict(repo).startswith("# NOT REPRODUCED")


def test_a_manifest_may_differ_only_in_the_hash_of_a_numerically_equal_file(
        tmp_path, monkeypatch):
    repo, out = _outputs_tree(tmp_path, monkeypatch)
    (out / "sums.csv").write_text(f"file,sha256\nt.csv,{_sha(out / 't.csv')}\n")
    base = E.baseline(tmp_path / "keep")
    _rerun(out, csv="a,b\n1.0000000000001,2\n")
    (out / "sums.csv").write_text(f"file,sha256\nt.csv,{_sha(out / 't.csv')}\n")
    assert E.verify(base)
    assert ("differs (hashes of numerically equal files)\treports/sums.csv"
            in _report(repo))
    base = E.baseline(tmp_path / "keep")                    # any other change fails
    _rerun(out, csv="a,b\n1.0000000000001,2\n")
    (out / "sums.csv").write_text(f"file,sha256\nt.csv,{_sha(out / 't.csv')}\nx,y\n")
    assert not E.verify(base)


def test_verify_fails_an_output_the_run_did_not_rewrite(tmp_path, monkeypatch):
    repo, out = _outputs_tree(tmp_path, monkeypatch)
    base = E.baseline(tmp_path / "keep")
    _rerun(out, skip=("t.csv",))                            # a stale committed output
    assert not E.verify(base)
    report = _report(repo)
    assert "not rewritten\treports/t.csv" in report
    assert "identical\treports/input.parquet" in report      # an input may stay untouched


def test_verify_names_a_changed_figure_and_a_new_file(tmp_path, monkeypatch):
    repo, out = _outputs_tree(tmp_path, monkeypatch)
    (out / "f.png").write_bytes(b"\x89PNG one")
    base = E.baseline(tmp_path / "keep")
    _rerun(out)
    (out / "f.png").write_bytes(b"\x89PNG two")
    (out / "extra.csv").write_text("x")
    assert not E.verify(base)
    report = _report(repo)
    assert "differs (figure)\treports/f.png" in report and "new\treports/extra.csv" in report


def test_an_archive_keeps_its_first_baseline_across_runs(tmp_path, monkeypatch):
    """No .git: an interrupted run must not become the next run's baseline."""
    repo, out = _outputs_tree(tmp_path, monkeypatch)
    E.baseline(tmp_path / "keep")                            # first run: records the tree
    (out / "t.csv").write_text("a,b\n9,9\n")                 # ... and is interrupted
    base = E.baseline(tmp_path / "keep2")                    # the rerun compares with the
    _rerun(out, csv="a,b\n9,9\n")                            # tree as unpacked
    assert not E.verify(base)
    assert "differs\treports/t.csv" in _report(repo)


def test_a_clone_is_compared_with_head_not_the_working_tree(tmp_path, monkeypatch):
    repo, out = _outputs_tree(tmp_path, monkeypatch, git=True)
    (out / "t.csv").write_text("a,b\n9,9\n")                 # left by an interrupted run
    base = E.baseline(tmp_path / "keep")
    assert base["git"]
    _rerun(out, csv="a,b\n9,9\n")
    assert not E.verify(base)
    assert "differs\treports/t.csv" in _report(repo)
    base = E.baseline(tmp_path / "keep2")
    _rerun(out)
    assert E.verify(base)


def test_verify_reads_a_parquet_with_missing_values(tmp_path):
    a, b = tmp_path / "a.parquet", tmp_path / "b.parquet"
    frame = pd.DataFrame({"x": pd.array([1.0, None], dtype="Float64"), "y": ["p", "q"]})
    frame.to_parquet(a)
    frame.to_parquet(b, compression="gzip")                  # other bytes, same cells
    assert E._numerically_equal(a, b)
    frame.assign(y=["p", "r"]).to_parquet(b)
    assert not E._numerically_equal(a, b)


def test_the_one_command_fetches_and_verifies():
    script = (REPO / "reproduce.sh").read_text()
    assert "requirements-evidence.lock.txt" in script
    assert "reproduce_evidence.py --fetch-inputs --verify" in script
    assert "reproduce_evidence.py --check" in script


def test_every_listed_input_is_a_committed_file():
    if not (REPO / ".git").exists():
        pytest.skip("needs the list of tracked files")
    inputs = E._run_inputs()
    tracked = set(subprocess.run(
        ["git", "ls-files", "--", "phase1/reports/evidence", "phase1/data/evidence"],
        cwd=REPO, capture_output=True, text=True, check=True).stdout.split())
    assert inputs <= tracked, sorted(inputs - tracked)[:5]


# --------------------------------------------------------------------------
# downloads that survive the network
# --------------------------------------------------------------------------
def _serve(files: dict, cut_first: set = frozenset()):
    """A local server; paths in `cut_first` announce the whole body and send half on
    their first request, as a dropped connection does. Records each request's Range."""
    import http.server
    import threading
    seen = []

    class Handler(http.server.BaseHTTPRequestHandler):
        def log_message(self, *a):
            pass

        def do_GET(self):
            body, rng = files.get(self.path), self.headers.get("Range")
            seen.append((self.path, rng))
            if body is None:
                self.send_error(404)
                return
            if rng:
                start = int(rng.split("=")[1].rstrip("-"))
                self.send_response(206)
                self.send_header("Content-Range", f"bytes {start}-{len(body) - 1}/{len(body)}")
                self.send_header("Content-Length", str(len(body) - start))
                self.end_headers()
                self.wfile.write(body[start:])
                return
            self.send_response(200)
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            if self.path in cut_first and sum(p == self.path for p, _ in seen) == 1:
                self.wfile.write(body[: len(body) // 2])
                self.close_connection = True
                return
            self.wfile.write(body)

    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    return server, f"http://127.0.0.1:{server.server_address[1]}", seen


def test_a_download_cut_short_resumes_from_where_it_stopped(tmp_path, monkeypatch):
    payload = bytes(range(256)) * 4096                   # 1 MiB
    server, base, seen = _serve({"/f": payload}, cut_first={"/f"})
    monkeypatch.setattr(E.time, "sleep", lambda s: None)
    try:
        E._download(base + "/f", tmp_path / "file.bin")
    finally:
        server.shutdown()
    assert (tmp_path / "file.bin").read_bytes() == payload
    assert seen == [("/f", None), ("/f", f"bytes={len(payload) // 2}-")]
    assert not (tmp_path / "file.bin.part").exists()


def test_a_missing_file_is_not_retried(tmp_path, monkeypatch):
    server, base, seen = _serve({})
    monkeypatch.setattr(E.time, "sleep", lambda s: None)
    try:
        with pytest.raises(OSError):
            E._download(base + "/gone", tmp_path / "file.bin")
    finally:
        server.shutdown()
    assert len(seen) == 1


def test_clinvar_moves_on_from_a_source_that_serves_the_wrong_file(tmp_path, monkeypatch):
    right, wrong = b"the release" * 1000, b"a re-issued file" * 1000
    md5 = hashlib.md5(wrong).hexdigest().encode()
    server, base, seen = _serve({"/ncbi/c.vcf.gz": wrong, "/ncbi/c.vcf.gz.md5": md5,
                                 "/mirror/c.vcf.gz": right})
    monkeypatch.setattr(E, "CLINVAR_VCF", tmp_path / "c.vcf.gz")
    monkeypatch.setattr(E, "REPO", tmp_path)
    monkeypatch.setattr(E, "CLINVAR_SHA256", hashlib.sha256(right).hexdigest())
    monkeypatch.setattr(E, "CLINVAR_SOURCES", [
        (base + "/gone/c.vcf.gz", base + "/gone/c.vcf.gz.md5"),     # 404
        (base + "/ncbi/c.vcf.gz", base + "/ncbi/c.vcf.gz.md5"),     # its own md5 fits,
        (base + "/mirror/c.vcf.gz", None)])                          # the sha256 not
    monkeypatch.setattr(E.time, "sleep", lambda s: None)
    try:
        E.fetch_clinvar()
        assert (tmp_path / "c.vcf.gz").read_bytes() == right
        (tmp_path / "c.vcf.gz").unlink()
        monkeypatch.setattr(E, "CLINVAR_SHA256", "0" * 64)          # no source serves it
        with pytest.raises(SystemExit):
            E.fetch_clinvar()
        assert not (tmp_path / "c.vcf.gz").exists()
    finally:
        server.shutdown()


# --------------------------------------------------------------------------
# --check: the quick path
# --------------------------------------------------------------------------
def test_every_recorded_checksum_names_a_file_of_its_own():
    rows, _ = E.recorded_checksums()
    for record, target, sha, how in rows:
        assert target != REPO / record, f"{record} points at itself"
        assert len(sha) == 64 and how in ("file", "content", "source")
    inside = {t for _, t, _, _ in rows if t.is_relative_to(REPO)}
    assert len(inside) >= 140


def test_the_quick_check_waits_for_an_absent_clinvar(tmp_path, monkeypatch):
    monkeypatch.setattr(E, "CLINVAR_VCF", tmp_path / "not-downloaded.vcf.gz")
    rows, absent = E.recorded_checksums()
    assert "the ClinVar release" in absent
    assert not any(t == tmp_path / "not-downloaded.vcf.gz" for _, t, _, _ in rows)
