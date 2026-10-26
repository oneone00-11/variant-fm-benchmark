"""E20b.2 (commit 2): the ENIGMA BRCA1/BRCA2 specification files that hold Figure 1A and
Appendix J, as the Criteria Specification Registry serves them. No predictor score is read.

Where. The registry's records (E20.3) name "Specifications Figure 1A and Appendix J"; the
API's JSON carries no file, but the registry's page for each specification lists the files
attached to its rule set (GN092 for BRCA1, GN097 for BRCA2), each version 1.2:
Specifications, Appendix, Supplementary tables, Tables 4 and 9. ClinGen's page for the
panel (the records' affiliation, 50087) is read for the links it carries, recorded by
sha256 and not saved. The registry is the first source E20b.2 names, so the published
paper's supplement (the fallback) is not needed. Fetched here: both registry pages, and
from each specification the two files Figure 1A and Appendix J are in
(Specifications_V1.2, Appendix_V1.2); the three tables are listed, not fetched.

Licence. ClinGen's terms of use (raw/clingen_terms_of_use.html) put the curated content it
publishes under CC0; the panel's specification as the journal published it (Parsons et al.
2024, doi:10.1016/j.ajhg.2024.07.013, Data S2; its version not read here) is under the
publisher's licences, which Crossref's record of the paper lists (saved in raw/, Crossref's
metadata being open): none is an open licence. The files are therefore treated as of
unclear licence (E20b.2): kept in the gitignored phase1/data/e20b_cache/, recorded by URL
and sha256, not committed. The registry's pages are committed. PMC's page for the paper is
not requested: it served this work a challenge page on 2026-10-26, which was not answered.

A registry response that is not HTTP 200 or that is a challenge page stops the run; no
challenge is answered. ClinGen's panel page and Crossref's record are only read for what
they say: either one not served is recorded as such and the run goes on. Every file is recorded
in raw/manifest.json: the committed ones in "files" (E20.3's form), the cached ones in
"enigma_specification" with their cache paths.

Run from the repository root, with network access:
    .venv/bin/python scripts/stage_E20b/fetch_enigma_specification.py
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
import urllib.error
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel, utc_now  # noqa: E402

RAW = REPO / "phase1/data/evidence/vcep/raw"
MANIFEST = RAW / "manifest.json"
CACHE = REPO / "phase1/data/e20b_cache/enigma_specification"
CSPEC = "https://cspec.genome.network/cspec"
SPECS = {"GN092": "BRCA1", "GN097": "BRCA2"}
WANTED = ["Specifications_V1.2", "Appendix_V1.2"]
AFFILIATION = "https://clinicalgenome.org/affiliation/50087/"
CROSSREF = "https://api.crossref.org/works/10.1016/j.ajhg.2024.07.013"
PAPER = {"doi": "10.1016/j.ajhg.2024.07.013", "pmid": "39142283", "pmcid": "PMC11393667",
         "supplement_item": "Data S2. Approved specifications developed by the ClinGen ENIGMA BRCA1 "
                            "and BRCA2 Variant Curation Expert Panel (mmc7.zip)"}
AGENT = "variant-fm-benchmark E20b record fetch (urllib)"
CHALLENGE = re.compile(rb"(?i)just a moment|captcha|proof.of.work|POW_CHALLENGE|Preparing to download")


class NotServed(Exception):
    """Not HTTP 200, or a challenge page (never answered)."""


def get(url: str) -> tuple[bytes, dict]:
    req = urllib.request.Request(url, headers={"User-Agent": AGENT})
    try:
        with urllib.request.urlopen(req, timeout=180) as r:
            body = r.read()
            info = {"url": url, "final_url": r.geturl(), "http_status": r.status,
                    "content_type": r.headers.get("Content-Type"),
                    "content_disposition": r.headers.get("Content-Disposition"), "fetched_utc": utc_now()}
    except urllib.error.HTTPError as e:
        raise NotServed(f"HTTP {e.code}") from None
    if info["http_status"] != 200:
        raise NotServed(f"HTTP {info['http_status']}")
    if CHALLENGE.search(body[:20000]):
        raise NotServed("a challenge page")
    return body, info


def required(url: str) -> tuple[bytes, dict]:
    try:
        return get(url)
    except NotServed as e:
        raise SystemExit(f"{url}: {e}: stopped") from None


def optional(url: str) -> tuple[bytes | None, dict]:
    try:
        return get(url)
    except NotServed as e:
        return None, {"url": url, "not_served": str(e), "fetched_utc": utc_now()}


ITEM = re.compile(r'<span class="file-label">\s*([^<:]+?):\s*</span>\s*'
                  r'<span class="description">([^<]*)</span>.*?'
                  r'href="https://cspec\.genome\.network/cspec/File/id/([0-9a-f-]+)/data"', re.S)


def attachments(page: str) -> list[dict]:
    """(name, description, file id) of each file on a registry page's files panel: each item
    is a file-label span, a description span, then its download link."""
    out = [{"name": n.strip(), "description": d.strip(), "file_id": i} for n, d, i in ITEM.findall(page)]
    if len(out) != len(re.findall(r"/cspec/File/id/[0-9a-f-]+/data", page)):
        raise SystemExit("a file on the page has no label and description: stopped")
    return out


def main() -> None:
    man = json.loads(MANIFEST.read_text())
    files, section = [], {"what": "the ENIGMA specification files that hold Figure 1A and Appendix J "
                                  "(E20b.2)", "written_by": rel(Path(__file__)), "specifications": {}}
    with RunLog(__file__) as log:
        CACHE.mkdir(parents=True, exist_ok=True)
        for gn, gene in SPECS.items():
            url = f"{CSPEC}/ui/svi/doc/{gn}"
            body, info = required(url)
            path = RAW / f"{gn}_page.html"
            path.write_bytes(body)
            log.output(path)
            files.append(info | {"file": rel(path), "bytes": len(body),
                                 "sha256": hashlib.sha256(body).hexdigest()})
            listed = attachments(body.decode("utf-8", errors="replace"))
            names = [a["name"] for a in listed]
            for w in WANTED:
                if names.count(w) != 1:
                    raise SystemExit(f"{gn}: {names.count(w)} files named {w} on the page: stopped")
            entries = []
            for a in listed:
                entry = dict(a) | {"url": f"{CSPEC}/File/id/{a['file_id']}/data"}
                if a["name"] in WANTED:
                    data, finfo = required(entry["url"])
                    fname = re.search(r'filename="?([^";]+)"?', finfo["content_disposition"] or "")
                    fname = fname.group(1) if fname else f"{a['name']}.bin"
                    dest = CACHE / gn / fname
                    dest.parent.mkdir(parents=True, exist_ok=True)
                    dest.write_bytes(data)
                    log.output(dest)
                    entry |= finfo | {"cached_at": rel(dest), "committed": False, "bytes": len(data),
                                      "sha256": hashlib.sha256(data).hexdigest(), "filename": fname}
                else:
                    entry |= {"fetched": False, "why": "not needed: Figure 1A and Appendix J are in "
                                                       "the Specifications and Appendix files"}
                entries.append(entry)
            record = json.loads((RAW / f"{gn}_record.json").read_text())["data"]["entContent"]
            version = record["version"]
            for e in entries:
                if e.get("cached_at") and not (e["name"].endswith(f"_V{version}")
                                               and f"_V{version}_" in e["filename"]):
                    raise SystemExit(f"{gn}: {e['name']} ({e['filename']}) is not version {version}, "
                                     "the registry record's: stopped (E20b.6)")
            notes = record.get("releaseNotes") or ""
            section["specifications"][gn] = {
                "gene": gene, "page": rel(path), "record_version": version,
                "record_approved_on": record.get("approvedOn"),
                "release_notes_on_spliceai": [x.strip() for x in notes.splitlines() if "SpliceAI" in x],
                "attachments": entries}
        a92 = {e["name"]: e.get("sha256") for e in section["specifications"]["GN092"]["attachments"]}
        a97 = {e["name"]: e.get("sha256") for e in section["specifications"]["GN097"]["attachments"]}
        section["same_files_for_brca1_and_brca2"] = {w: a92[w] == a97[w] for w in WANTED}
        body, info = optional(CROSSREF)
        journal = PAPER | {"supplement_read": False,
                           "pmc_page": "not requested: it served this work a challenge page on 2026-10-26, "
                                       "which was not answered"}
        if body:
            path = RAW / f"crossref_{PAPER['doi'].replace('/', '_')}.json"
            path.write_bytes(body)
            log.output(path)
            files.append(info | {"file": rel(path), "bytes": len(body),
                                 "sha256": hashlib.sha256(body).hexdigest()})
            journal |= {"crossref_record": rel(path),
                        "licences": [{"url": x.get("URL"), "content_version": x.get("content-version"),
                                      "from": (x.get("start") or {}).get("date-time")}
                                     for x in json.loads(body)["message"].get("license", [])]}
        else:
            journal |= {"crossref_record": info}
        section["licence"] = {
            "clingen_terms": "raw/clingen_terms_of_use.html: curated content published by ClinGen under CC0 1.0",
            "the_specification_as_the_journal_published_it": journal,
            "decision": ("unclear: kept in the gitignored cache, recorded by URL and sha256, not committed "
                         "(E20b.2, as Sahu 2025's PDF in E19)")}
        body, info = optional(AFFILIATION)
        links = sorted(set(re.findall(r'href="([^"]+)"', body.decode("utf-8", errors="replace")))) if body else []
        section["searched"] = {
            "registry_records": "raw/GN092_record.json, raw/GN097_record.json (E20.3): no file in the JSON",
            "registry_pages": "the files panel of each specification's page: found",
            "clingen_panel_page": info | ({
                "bytes": len(body), "sha256": hashlib.sha256(body).hexdigest(), "page_saved": False,
                "links_to_the_registry": [x for x in links if "cspec.genome.network" in x],
                "links_to_a_file": [x for x in links if re.search(r"(?i)\.(docx?|pdf|xlsx?|zip)(\?|$)", x)],
                "links_to_the_paper": [x for x in links if "10.1016/j.ajhg.2024.07.013" in x
                                       or "39142283" in x]} if body else {}),
            "published_paper": "not needed (the registry holds the files)"}
        man["files"] = [f for f in man["files"] if f["file"] not in {x["file"] for x in files}] + files
        man["enigma_specification"] = section
        MANIFEST.write_text(json.dumps(man, indent=1, ensure_ascii=False) + "\n")
        log.output(MANIFEST)
        print(json.dumps({gn: [(e["name"], e.get("filename"), e.get("bytes")) for e in s["attachments"]]
                          for gn, s in section["specifications"].items()}, indent=1))
        print("same files:", section["same_files_for_brca1_and_brca2"], "| journal licences:", section[
            "licence"]["the_specification_as_the_journal_published_it"].get("licences"))


if __name__ == "__main__":
    main()
