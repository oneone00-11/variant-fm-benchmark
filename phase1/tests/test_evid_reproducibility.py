"""Guardrails on reproducing the evidence-strength outputs from public sources.

A clean-room run (a fresh clone, the atlas release archive from Zenodo, a Python
environment built from requirements-evidence.lock.txt, ClinVar from NCBI) found three
committed files that depended on the machine that wrote them: the analysis-set
manifest read the atlas commit from `git` and recorded "unknown" for the archive, the
TP53 manifest checked its REF bases against a 706 MB local FASTA that no archive
carries, and the supplementary-table manifest hashes both. These tests hold the fixes:
the entry point's inputs agree with what the manifests record, it finds the atlas the
way it documents, and neither manifest carries a local path.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import subprocess
import sys
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
    assert all(u.endswith("/" + E.CLINVAR_VCF.name) for u in E.CLINVAR_URLS)


def test_the_fetched_atlas_release_is_the_one_the_manifest_records():
    sources = json.loads(SET_MANIFEST.read_text())["sources"]
    version = E.ATLAS_ARCHIVE["version"]
    assert sources["atlas_matrix"]["atlas_git_head"] == B.ATLAS_RELEASE_COMMIT[version]
    for key in ("atlas_matrix", "atlas_alphagenome_v061"):
        path = sources[key]["path"]
        assert path.startswith("atlas:") and path[len("atlas:"):] in E.ATLAS_NEEDS


def test_every_atlas_file_the_set_build_reads_is_checked_before_the_run():
    for p in (B.ATLAS_MATRIX, B.ATLAS_V061):
        assert str(p.relative_to(B.ATLAS_REPO)) in E.ATLAS_NEEDS


def test_neither_manifest_carries_a_local_path():
    assert _no_local_path(json.loads(SET_MANIFEST.read_text())["sources"])
    assert _no_local_path(json.loads(TP53_MANIFEST.read_text())["reference_alleles"])


# --------------------------------------------------------------------------
# where the entry point looks for the atlas
# --------------------------------------------------------------------------
def _atlas(root: Path, complete: bool = True) -> Path:
    for f in E.ATLAS_NEEDS if complete else E.ATLAS_NEEDS[:3]:
        (root / f).parent.mkdir(parents=True, exist_ok=True)
        (root / f).write_text('version: "2.5.0"\n' if f == "CITATION.cff" else "")
    return root


def test_the_atlas_is_found_in_the_documented_order(tmp_path, monkeypatch):
    beside, fetched = tmp_path / "beside", tmp_path / "fetched"
    monkeypatch.setattr(E, "ATLAS_BESIDE", beside)
    monkeypatch.setattr(E, "ATLAS_FETCHED", fetched)
    monkeypatch.delenv("EVID_ATLAS_REPO", raising=False)

    path, missing = E.resolve_atlas()                   # nothing anywhere
    assert missing == E.ATLAS_NEEDS

    _atlas(beside, complete=False)                      # a bare clone beside
    path, missing = E.resolve_atlas()
    assert path == beside and "results/score_matrix_atlas_v2.parquet" in missing

    _atlas(fetched)                                     # the fetched archive wins
    assert E.resolve_atlas() == (fetched, [])

    _atlas(beside)                                      # a complete checkout beside
    assert E.resolve_atlas() == (beside, [])

    other = _atlas(tmp_path / "other", complete=False)  # an explicit path is never
    monkeypatch.setenv("EVID_ATLAS_REPO", str(other))   # silently replaced
    path, missing = E.resolve_atlas()
    assert path == other and missing


# --------------------------------------------------------------------------
# the atlas commit is the same from a checkout and from the release archive
# --------------------------------------------------------------------------
def _git(*args, cwd):
    subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", *args], cwd=cwd,
                   check=True, capture_output=True)


def test_an_archive_records_the_commit_its_release_was_cut_from(tmp_path):
    outer = tmp_path / "vfm"
    outer.mkdir()
    _git("init", "-q", cwd=outer)                       # the archive is unpacked inside
    (outer / "x").write_text("x")                       # another checkout, whose HEAD
    _git("add", "x", cwd=outer)                         # must not be taken for the
    _git("commit", "-q", "-m", "x", cwd=outer)          # atlas's
    archive = _atlas(outer / "phase1/data/evidence/companion_atlas/atlas")
    assert B.git_head(archive) == "unknown"
    assert B.atlas_commit(archive) == B.ATLAS_RELEASE_COMMIT["2.5.0"]


def test_a_checkout_records_its_own_head(tmp_path):
    atlas = _atlas(tmp_path / "functional-standard-atlas")
    _git("init", "-q", cwd=atlas)
    _git("add", "CITATION.cff", cwd=atlas)
    _git("commit", "-q", "-m", "c", cwd=atlas)
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=atlas, check=True,
                          capture_output=True, text=True).stdout.strip()
    assert B.atlas_commit(atlas) == head != B.ATLAS_RELEASE_COMMIT["2.5.0"]


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


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))


# --------------------------------------------------------------------------
# --verify: the verdict the one command ends with
# --------------------------------------------------------------------------
def test_verify_passes_only_when_every_output_matches(tmp_path, monkeypatch, capsys):
    out = tmp_path / "reports"
    out.mkdir()
    (out / "t.csv").write_text("a,b\n1,2\n")
    (out / "m.manifest.json").write_text('{"n": 1, "built_utc": "2026-01-01"}')
    (out / "input.parquet").write_text("read, never written")
    monkeypatch.setattr(E, "REPO", tmp_path)
    monkeypatch.setattr(E, "OUTPUT_ROOTS", [out])
    monkeypatch.setattr(E, "NOT_OUTPUTS", [])

    def rewrite(csv="a,b\n1,2\n", manifest='{"n": 1, "built_utc": "2026-09-29"}'):
        before = E.snapshot()
        for f in ("t.csv", "m.manifest.json"):
            (out / f).touch()
        (out / "t.csv").write_text(csv)
        (out / "m.manifest.json").write_text(manifest)
        return E.verify(before)

    assert rewrite()                                   # only the build time moved
    assert "REPRODUCED" in capsys.readouterr().out
    report = (tmp_path / "reproduce_report.txt").read_text()
    assert "identical (build time aside)\treports/m.manifest.json" in report
    assert not rewrite(csv="a,b\n1.0000000000001,2\n")  # same number, other bytes:
    report = (tmp_path / "reproduce_report.txt").read_text()   # named, still a fail
    assert "differs (numerically equal)\treports/t.csv" in report
    assert not rewrite(csv="a,b\n1,3\n")               # a table changed
    assert not rewrite(manifest='{"n": 2, "built_utc": "2026-09-29"}')
    (out / "extra.csv").write_text("x")                # a file no commit carries
    before = E.snapshot()
    (out / "extra.csv").unlink()
    assert not E.verify(before)                        # an output went missing


def test_the_one_command_fetches_and_verifies():
    script = (REPO / "reproduce.sh").read_text()
    assert "requirements-evidence.lock.txt" in script
    assert "reproduce_evidence.py --fetch-inputs --verify" in script


def test_atlas_sources_are_named_the_same_wherever_the_atlas_sits(monkeypatch):
    from src import evid_supp_tables as S
    inside = PHASE1 / "data/evidence/companion_atlas/atlas"      # --fetch-inputs
    monkeypatch.setattr(S, "ATLAS_REPO", inside)
    assert S._display(inside / "src/atlas/predictor_resources.py") == \
        "atlas:src/atlas/predictor_resources.py"
    assert S._display(PHASE1 / "src/evid_common.py") == "phase1/src/evid_common.py"


# --------------------------------------------------------------------------
# a download cut short resumes where it stopped
# --------------------------------------------------------------------------
def test_a_download_cut_short_resumes_from_where_it_stopped(tmp_path, monkeypatch):
    import http.server
    import threading

    payload = bytes(range(256)) * 4096                   # 1 MiB
    requests = []

    class Handler(http.server.BaseHTTPRequestHandler):
        def log_message(self, *a):
            pass

        def do_GET(self):
            rng = self.headers.get("Range")
            requests.append(rng)
            if rng is None:                                # the first try: announce
                self.send_response(200)                    # it all, send half, hang up
                self.send_header("Content-Length", str(len(payload)))
                self.end_headers()
                self.wfile.write(payload[: len(payload) // 2])
                self.wfile.flush()
                self.close_connection = True
                return
            start = int(rng.split("=")[1].rstrip("-"))
            self.send_response(206)
            self.send_header("Content-Range", f"bytes {start}-{len(payload) - 1}/{len(payload)}")
            self.send_header("Content-Length", str(len(payload) - start))
            self.end_headers()
            self.wfile.write(payload[start:])

    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    monkeypatch.setattr(E.time, "sleep", lambda s: None)
    try:
        dest = tmp_path / "file.bin"
        E._download(f"http://127.0.0.1:{server.server_address[1]}/f", dest)
    finally:
        server.shutdown()
    assert dest.read_bytes() == payload
    assert requests[0] is None and requests[1] == f"bytes={len(payload) // 2}-"
    assert not (tmp_path / "file.bin.part").exists()


def test_verify_names_a_changed_figure_as_one(tmp_path, monkeypatch):
    out = tmp_path / "figures"
    out.mkdir()
    (out / "f.png").write_bytes(b"\x89PNG one")
    monkeypatch.setattr(E, "REPO", tmp_path)
    monkeypatch.setattr(E, "OUTPUT_ROOTS", [out])
    monkeypatch.setattr(E, "NOT_OUTPUTS", [])
    before = E.snapshot()
    (out / "f.png").write_bytes(b"\x89PNG two")
    assert not E.verify(before)
    assert "differs (figure)\tfigures/f.png" in (tmp_path / "reproduce_report.txt").read_text()


# --------------------------------------------------------------------------
# a ClinVar source without an md5 file (a mirror) is checked by sha256 alone
# --------------------------------------------------------------------------
def _serve(files: dict):
    import http.server
    import threading

    class Handler(http.server.BaseHTTPRequestHandler):
        def log_message(self, *a):
            pass

        def do_GET(self):
            body = files.get(self.path)
            if body is None:
                self.send_error(404)
                return
            self.send_response(200)
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    return server, f"http://127.0.0.1:{server.server_address[1]}"


def test_clinvar_falls_back_to_a_mirror_and_checks_its_sha256(tmp_path, monkeypatch):
    body = b"a stand-in for the release" * 1000
    server, base = _serve({"/mirror/clinvar.vcf.gz": body})
    monkeypatch.setattr(E, "CLINVAR_VCF", tmp_path / "clinvar.vcf.gz")
    monkeypatch.setattr(E, "REPO", tmp_path)
    monkeypatch.setattr(E, "CLINVAR_SOURCES", [
        (base + "/ncbi/clinvar.vcf.gz", base + "/ncbi/clinvar.vcf.gz.md5"),   # 404
        (base + "/mirror/clinvar.vcf.gz", None)])
    monkeypatch.setattr(E.time, "sleep", lambda s: None)
    try:
        monkeypatch.setattr(E, "CLINVAR_SHA256", hashlib.sha256(body).hexdigest())
        E.fetch_clinvar()
        assert (tmp_path / "clinvar.vcf.gz").read_bytes() == body
        (tmp_path / "clinvar.vcf.gz").unlink()
        monkeypatch.setattr(E, "CLINVAR_SHA256", "0" * 64)          # a wrong file
        with pytest.raises(SystemExit):
            E.fetch_clinvar()
        assert not (tmp_path / "clinvar.vcf.gz").exists()
    finally:
        server.shutdown()
