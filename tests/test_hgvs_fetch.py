"""E9.4: the script that asked VariantValidator about every listed variant
(scripts/stage_E9/fetch_hgvs_check.py) and the one that takes the excerpt from the stored
replies (extract_hgvs_check.py), with a stand-in for the service: no test here asks it.

What is held: how the descriptions are cut into requests and written into the address;
that an entry is found by the description submitted, never by its key; that a run stores
each reply once, can be started again without asking twice, tries again after a 429, a
5xx or a reply that lacks an entry, and stops on a status it does not retry, on three
rounds of failures, or when the versions change; that the excerpt copies and does not
judge; and that the archive of the raw replies depends on the files alone.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import sys
import tarfile
import urllib.parse
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts" / "stage_E9"))
import extract_hgvs_check as X  # noqa: E402
import fetch_hgvs_check as F  # noqa: E402

META = {"variantvalidator_hgvs_version": "4.0.0", "variantvalidator_version": "4.0.1", "vvdb_version": "vvdb_1",
        "vvseqrepo_db": "VV_SR_1", "vvta_version": "vvta_1"}
COLUMNS = ["Gene", "RefSeq transcript (MANE Select)", "HGVS c.", "Ensembl transcript (MANE Select)"]


def _rows(n_bap1: int = 12, n_sbds: int = 3) -> list[dict]:
    rows = [dict(zip(COLUMNS, ["BAP1", "NM_004656.4", f"c.37+{k}G>A", "ENST00000460680.6"])) for k in range(3, 3 + n_bap1)]
    rows += [dict(zip(COLUMNS, ["DDX3X", "NM_001356.5", "c.45+3G>A", "ENST00000644876.2"]))]
    rows += [dict(zip(COLUMNS, ["SBDS", "NM_016038.4", f"c.129-{k}C>T", "ENST00000246868.7"])) for k in range(3, 3 + n_sbds)]
    return rows


def _write_list(path: Path, rows: list[dict]) -> Path:
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    return path


class Service:
    """A stand-in for F.get: answers /hello/, the transcripts' lists and the descriptions
    of a request; `script` holds answers to give first, in order, to the requests for
    descriptions: (status, body, headers)."""

    def __init__(self, script=(), meta=None, rewrite=None):
        self.script, self.meta, self.rewrite = list(script), dict(meta or META), rewrite or {}
        self.asked: list[str] = []

    def __call__(self, url: str):
        self.asked.append(url)
        if url.endswith("/hello/"):
            return 200, {}, json.dumps({"metadata": self.meta, "status": "hello_world"}).encode(), 0.01
        if "/tools/gene2transcripts_v2/" in url:
            return 200, {}, b"[]", 0.01
        if self.script:
            status, body, headers = self.script.pop(0)
            return status, headers, body, 0.01
        return 200, {}, self.reply(url), 0.01

    def reply(self, url: str, leave_out: int = 0) -> bytes:
        joined = urllib.parse.unquote(url.split(f"/{F.BUILD}/")[1].split("/")[0])
        out = {}
        for d in joined.split("|")[leave_out:]:
            shown = self.rewrite.get(d, d)
            out[shown] = {"submitted_variant": d, "hgvs_transcript_variant": shown, "validation_warnings": [],
                          "primary_assembly_loci": {"grch38": {"hgvs_genomic_description": "NC_1:g.1A>C",
                                                               "vcf": {"chr": "1", "pos": "1", "ref": "A", "alt": "C"}}},
                          "annotations": {"mane_select": True}}
        return json.dumps(out | {"flag": "gene_variant", "metadata": self.meta | {"vvseqrepo_db": self.meta["vvseqrepo_db"] + "/master"}}).encode()

    def requests(self) -> list[str]:
        return [u for u in self.asked if f"/{F.BUILD}/" in u and "/tools/" not in u]


@pytest.fixture
def quiet(monkeypatch):
    monkeypatch.setattr(F.time, "sleep", lambda s: None)


def _args(tmp_path: Path, rows=None) -> argparse.Namespace:
    return argparse.Namespace(list=str(_write_list(tmp_path / "list.csv", rows or _rows())), raw=str(tmp_path / "raw"), base=F.BASE)


def _journal(tmp_path: Path) -> list[dict]:
    return [json.loads(ln) for ln in (tmp_path / "raw" / "journal.jsonl").read_text().splitlines()]


def test_every_row_is_submitted_on_its_refseq_transcript_and_sbds_again_on_its_ensembl_transcript():
    items = F.descriptions(_rows())
    assert [it["description"] for it in items[:2]] == ["NM_004656.4:c.37+3G>A", "NM_004656.4:c.37+4G>A"]
    assert [(it["interface"], it["gene"]) for it in items[16:]] == [(F.ENSEMBL, "SBDS")] * 3
    assert items[16]["description"] == "ENST00000246868.7:c.129-3C>T" and items[16]["row"] == 14
    assert not any(it["interface"] == F.ENSEMBL and it["gene"] == "DDX3X" for it in items)
    assert [it["row"] for it in items[:16]] == list(range(1, 17))


def test_a_request_holds_at_most_ten_descriptions_of_one_transcript_in_order():
    reqs = F.requests_of(F.descriptions(_rows()))
    assert [len(r) for r in reqs] == [10, 2, 1, 3, 3]
    assert all(len({(it["interface"], it["transcript"]) for it in r}) == 1 for r in reqs)
    assert [it["row"] for r in reqs for it in r] == list(range(1, 17)) + [14, 15, 16]


def test_the_address_carries_the_descriptions_encoded_and_names_the_one_transcript():
    req = F.requests_of(F.descriptions(_rows()))[1]
    url = F.url_of(req)
    assert url == (f"{F.BASE}/VariantValidator/variantvalidator/GRCh38/NM_004656.4%3Ac.37%2B13G%3EA%7CNM_004656.4%3Ac.37%2B14G%3EA"
                   "/NM_004656.4?content-type=application%2Fjson")
    assert F.url_of(F.requests_of(F.descriptions(_rows()))[-1]).startswith(f"{F.BASE}/VariantValidator/variantvalidator_ensembl/GRCh38/ENST")
    assert F.reply_path(Path("raw"), req) == Path("raw/replies/variantvalidator/NM_004656.4/00011_00012.json")


def test_an_entry_is_found_by_the_description_submitted_not_by_its_key():
    """Where the reference base does not match, the tool rewrites the description and keys
    the entry by the rewritten one: the entry is still the submitted description's."""
    service = Service(rewrite={"NM_004656.4:c.37+13G>A": "NM_004656.4:c.37+13="})
    req = F.requests_of(F.descriptions(_rows()))[1]
    body = service.reply(F.url_of(req))
    assert "NM_004656.4:c.37+13G>A" not in json.loads(body)
    found = F.entries(body)
    assert found["NM_004656.4:c.37+13G>A"]["hgvs_transcript_variant"] == "NM_004656.4:c.37+13="
    assert F.whole(body, req) == ""
    assert F.whole(service.reply(F.url_of(req), leave_out=1), req).startswith("no entry for 1 of 2 descriptions")
    assert F.whole(b"<html>busy</html>", req).startswith("not the JSON expected")
    assert F.versions(body)["vvseqrepo_db"] == "VV_SR_1/master"


def test_a_run_stores_each_reply_once_and_a_second_run_asks_nothing(tmp_path, monkeypatch, quiet, capsys):
    service = Service()
    monkeypatch.setattr(F, "get", service)
    args = _args(tmp_path)
    F.run(args)
    raw = tmp_path / "raw"
    stored = sorted(p.relative_to(raw / "replies").as_posix() for p in (raw / "replies").rglob("*.json"))
    assert stored == ["variantvalidator/NM_001356.5/00013_00013.json", "variantvalidator/NM_004656.4/00001_00010.json",
                      "variantvalidator/NM_004656.4/00011_00012.json", "variantvalidator/NM_016038.4/00014_00016.json",
                      "variantvalidator_ensembl/ENST00000246868.7/00014_00016.json"]
    assert len(service.requests()) == 5 and len(list(raw.glob("hello_*.json"))) >= 1
    assert sorted(p.name for p in raw.glob("transcripts_*")) == ["transcripts_DDX3X_ensembl.json", "transcripts_SBDS_ensembl.json"]
    for url in service.asked:                    # nothing but the descriptions: no key, no name, no address of a person
        assert "@" not in url and "key" not in url.lower() and "email" not in url.lower()
    before = len(service.asked)
    F.run(args)
    assert len(service.asked) == before          # every reply is stored: not even /hello/ is asked
    assert "5 already stored" in capsys.readouterr().out
    assert [r["outcome"] for r in _journal(tmp_path) if r["what"] == "request"] == ["stored"] * 5


def test_a_run_tries_again_after_a_429_a_5xx_and_a_reply_that_lacks_an_entry(tmp_path, monkeypatch, quiet):
    service = Service()
    first = F.url_of(F.requests_of(F.descriptions(_rows()))[0])
    service.script = [(429, b'{"message":"Rate limit hit"}', {"Retry-After": "7"}), (500, b"An internal server error occurred.", {}),
                      (200, service.reply(first, leave_out=1), {}), (0, b"URLError: timed out", {})]
    monkeypatch.setattr(F, "get", service)
    F.run(_args(tmp_path))
    outcomes = [r["outcome"] for r in _journal(tmp_path) if r["what"] == "request"]
    assert outcomes[:5] == ["failed (status 429); again after 7 s", "failed (status 500); again after 10 s",
                            "failed (no entry for 1 of 10 descriptions, the first NM_004656.4:c.37+3G>A); again after 20 s",
                            "failed (status 0); again after 40 s", "stored"]
    assert outcomes[5:] == ["stored"] * 4
    assert len(F.entries((tmp_path / "raw/replies/variantvalidator/NM_004656.4/00001_00010.json").read_bytes())) == 10


def test_a_run_stops_on_a_status_it_does_not_retry(tmp_path, monkeypatch, quiet):
    service = Service(script=[(404, b"not found", {})])
    monkeypatch.setattr(F, "get", service)
    with pytest.raises(SystemExit) as stop:
        F.run(_args(tmp_path))
    assert stop.value.code == 4 and not list((tmp_path / "raw").glob("replies/**/*.json"))


def test_a_run_pauses_after_five_failures_in_a_row_and_stops_at_the_third_round(tmp_path, monkeypatch):
    slept: list[float] = []
    monkeypatch.setattr(F.time, "sleep", slept.append)
    service = Service(script=[(503, b"busy", {})] * 15)
    monkeypatch.setattr(F, "get", service)
    with pytest.raises(SystemExit) as stop:
        F.run(_args(tmp_path))
    assert stop.value.code == 3 and len(service.requests()) == 15
    assert slept.count(F.ROUND_PAUSE) == 2
    assert [r["what"] for r in _journal(tmp_path)][-1] == "run stops"


def test_a_run_stops_when_the_versions_change_on_the_way(tmp_path, monkeypatch, quiet):
    service = Service()
    monkeypatch.setattr(F, "get", service)
    args = _args(tmp_path)
    F.run(args)
    (tmp_path / "raw/replies/variantvalidator/NM_016038.4/00014_00016.json").unlink()
    later = Service(meta=META | {"vvta_version": "vvta_2"})
    monkeypatch.setattr(F, "get", later)
    with pytest.raises(SystemExit) as stop:
        F.run(args)
    assert stop.value.code == 5 and not later.requests()        # /hello/ no longer answers as the first time


def test_the_second_asking_compares_each_entry_with_the_stored_one(tmp_path, monkeypatch, quiet, capsys):
    monkeypatch.setattr(F, "RECHECK_ROWS", 4)
    service = Service()
    monkeypatch.setattr(F, "get", service)
    args = _args(tmp_path)
    F.run(args)
    F.recheck(args)
    result = json.loads((tmp_path / "raw/recheck/comparison.json").read_text())
    assert len(result["rows"]) == 4 and all(e["same_entry"] for e in result["entries"])
    assert "the same as stored" in capsys.readouterr().out
    # an entry that comes back otherwise is reported with the fields that differ
    for p in (tmp_path / "raw/recheck").glob("variantvalidator*.json"):
        p.unlink()
    changed = Service(rewrite={it["description"]: "changed" for it in F.descriptions(_rows())})
    monkeypatch.setattr(F, "get", changed)
    F.recheck(args)
    result = json.loads((tmp_path / "raw/recheck/comparison.json").read_text())
    assert not any(e["same_entry"] for e in result["entries"])
    assert result["entries"][0]["differs_in"] == ["hgvs_transcript_variant"]


def test_the_excerpt_copies_the_entrys_fields_and_judges_nothing():
    item = {"description": "NM_1.1:c.1+3G>A", "interface": F.REFSEQ}
    entry = {"submitted_variant": "NM_1.1:c.1+3G>A", "hgvs_transcript_variant": "NM_1.1:c.1+3=",
             "validation_warnings": ["a warning, with a comma", "é: another"],
             "primary_assembly_loci": {"grch38": {"hgvs_genomic_description": "NC_1.1:g.5=",
                                                  "vcf": {"chr": "1", "pos": "5", "ref": "G", "alt": "G"}}},
             "annotations": {"mane_select": False}}
    row = X.excerpt_row(item, entry, "gene_variant")
    assert list(row) == X.COLUMNS
    assert row == {"submitted_description": "NM_1.1:c.1+3G>A", "interface": F.REFSEQ, "transcript_description": "NM_1.1:c.1+3=",
                   "genomic_description_grch38": "NC_1.1:g.5=", "vcf_chr": "1", "vcf_pos": "5", "vcf_ref": "G", "vcf_alt": "G",
                   "mane_select": "false", "flag": "gene_variant", "warnings": '["a warning, with a comma", "é: another"]'}
    # what the tool did not return is left empty, never filled in
    failed = X.excerpt_row(item, {"submitted_variant": "NM_1.1:c.1+3G>A", "hgvs_transcript_variant": "",
                                  "validation_warnings": ["The transcript is not in our database"],
                                  "primary_assembly_loci": {}, "annotations": {}}, "warning")
    assert [failed[k] for k in X.COLUMNS[2:10]] == ["", "", "", "", "", "", "", "warning"]


def test_the_excerpt_is_read_from_the_stored_replies_in_the_order_submitted(tmp_path, monkeypatch, quiet):
    service = Service(rewrite={"NM_016038.4:c.129-4C>T": "NM_016038.4:c.129-4="})
    monkeypatch.setattr(F, "get", service)
    args = _args(tmp_path)
    F.run(args)
    rows = F.read_list(Path(args.list))
    excerpt, versions = X.read_replies(tmp_path / "raw", rows)
    assert [e["submitted_description"] for e in excerpt] == [it["description"] for it in F.descriptions(rows)]
    assert excerpt[14]["transcript_description"] == "NM_016038.4:c.129-4=" and versions["vvta_version"] == "vvta_1"
    (tmp_path / "raw/replies/variantvalidator/NM_001356.5/00013_00013.json").unlink()
    with pytest.raises(SystemExit, match="no reply is stored for rows 13-13"):
        X.read_replies(tmp_path / "raw", rows)


def test_the_archive_of_the_raw_replies_depends_on_the_files_alone(tmp_path):
    raw = tmp_path / "raw"
    (raw / "replies" / "b").mkdir(parents=True)
    (raw / "replies" / "b" / "2.json").write_text("two")
    (raw / "replies" / "a.json").write_text("one")
    (raw / "journal.jsonl").write_text("{}\n")
    one = X.write_archive(raw, tmp_path / "x" / "one.tar.gz")
    (raw / "replies" / "a.json").touch()                         # another time, the same bytes
    two = X.write_archive(raw, tmp_path / "y" / "two.tar.gz")
    assert (one["sha256"], one["bytes"], one["files"]) == (two["sha256"], two["bytes"], 3)
    data = (tmp_path / "x" / "one.tar.gz").read_bytes()
    assert hashlib.sha256(data).hexdigest() == one["sha256"] and data[4:8] == b"\x00\x00\x00\x00"       # no time in gzip's header
    with tarfile.open(fileobj=io.BytesIO(gzip.decompress(data))) as tar:
        members = tar.getmembers()
    assert [m.name for m in members] == ["hgvs_raw/journal.jsonl", "hgvs_raw/replies/a.json", "hgvs_raw/replies/b/2.json"]
    assert {(m.mtime, m.uid, m.gid, m.mode) for m in members} == {(0, 0, 0, 0o644)}
    assert hashlib.sha256(gzip.decompress(data)).hexdigest() == one["tar_sha256"]
