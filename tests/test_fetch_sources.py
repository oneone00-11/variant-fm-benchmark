"""E9.1: the label audit's uncommitted files are fetched from sources tried in order, and a
reply that is not the file is a failed attempt. No test here touches the network: the
entry point's _download is replaced by a stand-in that serves prepared bytes, or fails,
for each address, and records the order in which the addresses were asked."""
from __future__ import annotations

import hashlib
import io
import json
import sys
import urllib.error
import zipfile
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import reproduce_evidence as E  # noqa: E402

THE_FILE = b"the file the manifest records\n" * 40
SHA = hashlib.sha256(THE_FILE).hexdigest()
MEMBER = "table.xlsx"


def zipped(members: dict[str, bytes]) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        for name, data in members.items():
            z.writestr(name, data)
    return buf.getvalue()


class Served:
    """Stands in for _download: `replies` maps an address to what it gives, attempt after
    attempt (the last repeats): bytes, written to the destination as a download would, or
    an exception, raised."""

    def __init__(self, replies: dict[str, list]):
        self.replies, self.asked = replies, []

    def __call__(self, url: str, dest: Path, attempts: int = 30) -> None:
        assert attempts == 1, "one download an attempt: the attempts are counted by the caller"
        n = self.asked.count(url)
        self.asked.append(url)
        reply = self.replies[url][min(n, len(self.replies[url]) - 1)]
        if isinstance(reply, Exception):
            raise reply
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(reply)


@pytest.fixture
def tree(tmp_path, monkeypatch):
    """A repository of its own, and no waiting between attempts."""
    monkeypatch.setattr(E, "REPO", tmp_path)
    waits = []
    return tmp_path, {"path": "cache/gene/table.xlsx", "sha256": SHA,
                      "source_url": "https://first.example/zip", "zip_member": MEMBER}, waits


def source(url: str, zip_member: str | None = None) -> dict:
    return {"url": url, "zip_member": zip_member, "what": f"the stand-in at {url}"}


def left_in(folder: Path) -> list[str]:
    return sorted(p.name for p in folder.iterdir()) if folder.is_dir() else []


def test_the_sources_are_tried_in_the_order_written_and_the_first_that_serves_the_file_ends_it(tree, monkeypatch):
    root, entry, waits = tree
    served = Served({"https://a.example/f": [b"another file"], "https://b.example/f": [THE_FILE],
                     "https://c.example/f": [THE_FILE]})
    monkeypatch.setattr(E, "_download", served)
    got = E.fetch_in_order(entry, [source("https://a.example/f"), source("https://b.example/f"),
                                   source("https://c.example/f")], attempts=3, pause=waits.append)
    assert got == "https://b.example/f"
    assert served.asked == ["https://a.example/f"] * 3 + ["https://b.example/f"]     # c is never asked
    assert waits == [2, 4]                                # between a's attempts, not after its last
    assert (root / entry["path"]).read_bytes() == THE_FILE
    assert left_in(root / "cache/gene") == ["table.xlsx"]


def test_a_file_counts_only_with_the_manifests_sha256_whichever_source_served_it(tree, monkeypatch):
    root, entry, waits = tree
    almost = THE_FILE[:-1] + b"?"
    for good_at in (0, 1, 2):
        urls = [f"https://s{i}.example/f" for i in range(3)]
        served = Served({u: [THE_FILE if i == good_at else almost] for i, u in enumerate(urls)})
        monkeypatch.setattr(E, "_download", served)
        (root / entry["path"]).unlink(missing_ok=True)
        assert E.fetch_in_order(entry, [source(u) for u in urls], attempts=1, pause=waits.append) == urls[good_at]
        assert hashlib.sha256((root / entry["path"]).read_bytes()).hexdigest() == SHA
        assert served.asked == urls[:good_at + 1]


def test_one_attempt_that_gets_other_bytes_leaves_nothing_behind_and_says_why(tree, monkeypatch):
    root, entry, _ = tree
    monkeypatch.setattr(E, "_download", Served({"https://a.example/f": [b"an error page"]}))
    why, pointless = E._take(source("https://a.example/f"), root / entry["path"], SHA)
    assert "sha256 differs" in why and "13 bytes" in why and pointless is False
    assert left_in(root / "cache/gene") == []


@pytest.mark.parametrize("reply, said", [
    (b"<html>Service temporarily unavailable</html>", "the reply is not a zip"),
    (zipped({"another.xlsx": THE_FILE}), f"the zip served has no member {MEMBER}"),
    (zipped({MEMBER: b"a different table"}), "sha256 differs"),
])
def test_a_reply_that_is_not_the_zip_of_the_file_is_a_failed_attempt_and_the_next_source_is_tried(tree, monkeypatch, reply, said):
    root, entry, waits = tree
    served = Served({"https://first.example/zip": [reply], "https://b.example/f": [THE_FILE]})
    monkeypatch.setattr(E, "_download", served)
    why, pointless = E._take(source("https://first.example/zip", MEMBER), root / entry["path"], SHA)
    assert said in why and pointless is False and left_in(root / "cache/gene") == []
    served.asked.clear()
    got = E.fetch_in_order(entry, [source("https://first.example/zip", MEMBER), source("https://b.example/f")],
                           attempts=2, pause=waits.append)
    assert got == "https://b.example/f" and served.asked == ["https://first.example/zip"] * 2 + ["https://b.example/f"]
    assert (root / entry["path"]).read_bytes() == THE_FILE and left_in(root / "cache/gene") == ["table.xlsx"]


def test_a_zip_that_holds_the_file_gives_the_member_and_a_later_attempt_can_succeed(tree, monkeypatch):
    root, entry, waits = tree
    served = Served({"https://first.example/zip": [b"not a zip this time", zipped({"readme.txt": b"x", MEMBER: THE_FILE})]})
    monkeypatch.setattr(E, "_download", served)
    assert E.fetch_in_order(entry, [source("https://first.example/zip", MEMBER)], attempts=5, pause=waits.append) == "https://first.example/zip"
    assert served.asked == ["https://first.example/zip"] * 2 and waits == [2]
    assert (root / entry["path"]).read_bytes() == THE_FILE and left_in(root / "cache/gene") == ["table.xlsx"]


def test_a_connection_that_fails_is_tried_again_and_a_client_error_is_not(tree, monkeypatch):
    root, entry, waits = tree
    gone = urllib.error.HTTPError("https://gone.example/f", 404, "Not Found", None, None)
    served = Served({"https://gone.example/f": [gone],
                     "https://down.example/f": [OSError("gave up on https://down.example/f after 1 attempt (timed out)")],
                     "https://b.example/f": [THE_FILE]})
    monkeypatch.setattr(E, "_download", served)
    got = E.fetch_in_order(entry, [source("https://gone.example/f"), source("https://down.example/f"),
                                   source("https://b.example/f")], attempts=4, pause=waits.append)
    assert got == "https://b.example/f"
    assert served.asked == ["https://gone.example/f"] + ["https://down.example/f"] * 4 + ["https://b.example/f"]
    assert waits == [2, 4, 6]


def test_when_no_source_serves_the_file_the_run_stops_and_says_what_each_one_gave(tree, monkeypatch):
    root, entry, waits = tree
    gone = urllib.error.HTTPError("https://gone.example/f", 404, "Not Found", None, None)
    served = Served({"https://first.example/zip": [b"<html>busy</html>"], "https://gone.example/f": [gone],
                     "https://down.example/f": [OSError("gave up on https://down.example/f after 1 attempt (timed out)")],
                     "https://other.example/f": [b"other bytes"]})
    monkeypatch.setattr(E, "_download", served)
    sources = [source("https://first.example/zip", MEMBER), source("https://gone.example/f"),
               source("https://down.example/f"), source("https://other.example/f")]
    with pytest.raises(SystemExit) as stop:
        E.fetch_in_order(entry, sources, attempts=2, pause=waits.append)
    said = str(stop.value)
    assert said.startswith("cache/gene/table.xlsx: no source served the file its manifest records.")
    for s, why in zip(sources, ("the reply is not a zip", "HTTP 404 Not Found", "timed out", "sha256 differs")):
        line = next(ln for ln in said.splitlines() if s["url"] in ln)
        assert why in line and s["what"] in line
    assert [said.index(s["url"]) for s in sources] == sorted(said.index(s["url"]) for s in sources)   # in the order tried
    assert served.asked == ["https://first.example/zip"] * 2 + ["https://gone.example/f"] + ["https://down.example/f"] * 2 + ["https://other.example/f"] * 2
    assert left_in(root / "cache/gene") == []             # no file, no half-download, no zip


def test_a_half_download_of_one_source_is_never_resumed_from_the_next(tree, monkeypatch):
    root, entry, waits = tree

    def cut_short(url, dest, attempts=30):
        if "a.example" in url:
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.with_name(dest.name + ".part").write_bytes(THE_FILE[:100])     # as _download leaves it
            raise OSError(f"gave up on {url} after 1 attempt (connection closed before {len(THE_FILE)} bytes)")
        assert not dest.with_name(dest.name + ".part").exists(), "the next source would resume another's bytes"
        dest.write_bytes(THE_FILE)

    monkeypatch.setattr(E, "_download", cut_short)
    assert E.fetch_in_order(entry, [source("https://a.example/f"), source("https://b.example/f")], attempts=2, pause=waits.append) == "https://b.example/f"
    assert left_in(root / "cache/gene") == ["table.xlsx"]


# ------------------------------------------------------------------ the table itself
def cached_entries() -> list[dict]:
    out = []
    for man in sorted(E.LABEL_AUDIT.glob("*_raw/manifest.json")) + sorted(E.SECOND_SET.glob("*_raw/label_source_manifest.json")):
        out += [e for e in json.loads(man.read_text())["files"] if e.get("cached")]
    return out


def test_the_table_lists_every_cached_file_and_nothing_else_with_its_manifests_source_once():
    table, entries = json.loads(E.FETCH_SOURCES.read_text()), cached_entries()
    assert sorted(table["files"]) == sorted(e["path"] for e in entries) and len(entries) == 5
    for e in entries:
        listed = table["files"][e["path"]]
        assert sum(1 for s in listed if s.get("manifest")) == 1, e["path"]
        for s in listed:
            assert s["what"] and (s.get("manifest") is True) != ("url" in s), s
            assert s.get("manifest") or s["url"].startswith("https://")
            assert s["checked"] == "" or len(s["checked"]) == 10        # a day, once the source has been checked
        sources = E.label_audit_sources(e, table)
        assert len(sources) == len(listed) and len({s["url"] for s in sources}) == len(sources)
        own = next(s for s, t in zip(sources, listed) if t.get("manifest"))
        assert own["url"] == e["source_url"] and own["zip_member"] == e.get("zip_member")
    assert set(table["release_archive"]["urls"]) <= {e["path"] for e in entries}
    assert all(u.startswith("https://") for u in table["release_archive"]["urls"].values())


def test_the_release_archive_is_tried_first_for_the_files_it_holds_once_its_addresses_are_filled_in():
    table, entries = json.loads(E.FETCH_SOURCES.read_text()), cached_entries()
    e, other = entries[0], entries[1]
    before = E.label_audit_sources(e, table)
    table["release_archive"]["urls"] = {e["path"]: "https://archive.example/records/1/files/the-file/content"}
    after = E.label_audit_sources(e, table)
    assert after[0] == {"url": "https://archive.example/records/1/files/the-file/content", "zip_member": None,
                        "what": "this release's archive of the inputs"}
    assert after[1:] == before
    assert E.label_audit_sources(other, table) == E.label_audit_sources(other, json.loads(E.FETCH_SOURCES.read_text()))


def test_a_cached_file_the_table_does_not_list_stops_the_run():
    table = json.loads(E.FETCH_SOURCES.read_text())
    with pytest.raises(SystemExit) as stop:
        E.label_audit_sources({"path": "phase1/data/evidence/label_audit_cache/new/file.json", "source_url": "https://x.example/f"}, table)
    assert "lists no source for phase1/data/evidence/label_audit_cache/new/file.json" in str(stop.value)


def test_the_fetch_asks_the_tables_sources_only_for_a_file_that_is_missing_or_differs(tmp_path, monkeypatch):
    """fetch_label_audit_cache over a manifest of its own: a file present with the manifest's
    sha256 is only checked; one that differs is removed and fetched from the table's sources."""
    monkeypatch.setattr(E, "REPO", tmp_path)
    audit = tmp_path / "audit"
    (audit / "g_raw").mkdir(parents=True)
    files = [{"path": "cache/g/here.json", "sha256": SHA, "source_url": "https://m.example/here", "cached": True},
             {"path": "cache/g/stale.json", "sha256": SHA, "source_url": "https://m.example/stale", "cached": True},
             {"path": "committed.json", "sha256": "0" * 64}]
    (audit / "g_raw" / "manifest.json").write_text(json.dumps({"files": files}))
    (tmp_path / "cache/g").mkdir(parents=True)
    (tmp_path / "cache/g/here.json").write_bytes(THE_FILE)
    (tmp_path / "cache/g/stale.json").write_bytes(b"an older file")
    table = {"release_archive": {"urls": {}}, "files": {
        "cache/g/here.json": [{"manifest": True, "what": "its own"}],
        "cache/g/stale.json": [{"url": "https://new.example/stale", "what": "a copy"}, {"manifest": True, "what": "its own"}]}}
    (tmp_path / "sources.json").write_text(json.dumps(table))
    monkeypatch.setattr(E, "LABEL_AUDIT", audit)
    monkeypatch.setattr(E, "SECOND_SET", tmp_path / "none")
    monkeypatch.setattr(E, "FETCH_SOURCES", tmp_path / "sources.json")
    from src import evid_label_sources as LS
    monkeypatch.setattr(LS, "extract_differences", lambda man, repo: ([], []))
    served = Served({"https://new.example/stale": [THE_FILE]})
    monkeypatch.setattr(E, "_download", served)
    E.fetch_label_audit_cache(attempts=2)
    assert served.asked == ["https://new.example/stale"]
    assert (tmp_path / "cache/g/stale.json").read_bytes() == THE_FILE
