"""E19 C1 (E19.4; the author's review of 2026-10-25, item 1): each training gene's SGE
study searched for a preprint in Europe PMC's preprint index, not left to Crossref's
relation fields or a reference list. Every response is saved as served, with
Crossref's record of every preprint the searches return.

The searches (Europe PMC REST search, SRC:PPR, every result page, core results):
  * for a gene whose cited score set records a paper (first_public_dates/
    manifest.json, MaveDB's publication identifiers): TITLE:"<gene>" AND AUTH:"<the
    paper's first author>", and TITLE:"<gene>" AND AUTH:"<its last author>", the
    authors' family names as Crossref's record of the paper gives them (and without
    their accents too, where they have any);
  * for a gene whose score set records no paper (BARD1 and PALB2): TITLE:"<gene>" AND
    AUTH:"<family name>" for every contributor a version of the score set names;
  * for every gene: TITLE:"<gene>" AND (TITLE:"saturation genome editing" OR
    ABSTRACT:"saturation genome editing").
Saved under phase1/data/evidence/e19/first_public_dates/: europepmc/ (every page of
every response) and crossref/ (each returned preprint DOI's record; a DOI already there
from fetch_first_public_dates.py is not fetched again), all listed in
preprint_search.json with URL, query, HTTP status, size, sha256 and time.

Nothing is decided here. Which preprint, if any, is a study's is decided offline by the
C1 stage (evid_suppl_c1_low.preprint_matches), by the rule written there.

Run from the repository root (the network is used here only):
    .venv/bin/python scripts/stage_E19/fetch_preprints.py
"""
from __future__ import annotations

import json
import sys
import time
import unicodedata
import urllib.parse
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel, sha256_of, utc_now  # noqa: E402

EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
CROSSREF = "https://api.crossref.org/works/"
DATES = REPO / "phase1/data/evidence/e19/first_public_dates"
OUT_EPMC = DATES / "europepmc"
OUT_CROSSREF = DATES / "crossref"
SGE = '(TITLE:"saturation genome editing" OR ABSTRACT:"saturation genome editing")'
PAGE = 1000

S = requests.Session()
S.trust_env = False


def _safe(name: str) -> str:
    return urllib.parse.quote(name, safe="").replace("%", "_")


def _get(url: str, params: dict | None = None) -> requests.Response:
    for attempt in range(5):
        try:
            r = S.get(url, params=params, timeout=120)
            if r.status_code in (429, 500, 502, 503, 504) and attempt < 4:
                time.sleep(2 ** attempt)
                continue
            return r
        except requests.RequestException:
            if attempt == 4:
                raise
            time.sleep(2 ** attempt)
    raise AssertionError("unreachable")


def _fold(s: str) -> str:
    return "".join(c for c in unicodedata.normalize("NFKD", s) if not unicodedata.combining(c))


def queries(man: dict) -> dict[str, list[dict]]:
    out = {}
    saved = {f["url"]: f.get("file") for f in man["files"]}
    for gene, g in man["genes"].items():
        qs = []
        papers = [d for d in g["dois"] if not d.startswith("10.1101/")]
        if papers:
            for doi in papers:
                path = saved.get(man["crossref_api"] + urllib.parse.quote(doi, safe="/"))
                msg = json.loads((REPO / path).read_text())["message"]
                people = [a["family"] for a in msg.get("author") or [] if a.get("family")]
                for role, name in (("first author", people[0]), ("last author", people[-1])):
                    for n in dict.fromkeys([name, _fold(name)]):
                        qs.append({"why": f"{role} of {doi}", "query":
                                   f'SRC:PPR AND TITLE:"{gene}" AND AUTH:"{n}"'})
        else:
            names = {}
            for urn in g["versions"]:
                meta = json.loads((REPO / saved[f"{man['mavedb_api']}/score-sets/{urn}"]).read_text())
                for c in meta.get("contributors") or []:
                    if c.get("familyName"):
                        names.setdefault(c["familyName"], urn)
            for name, urn in names.items():
                qs.append({"why": f"contributor named on {urn}", "query":
                           f'SRC:PPR AND TITLE:"{gene}" AND AUTH:"{name}"'})
        qs.append({"why": "saturation genome editing", "query":
                   f'SRC:PPR AND TITLE:"{gene}" AND {SGE}'})
        out[gene] = qs
    return out


def main() -> None:
    man = json.loads((DATES / "manifest.json").read_text())
    files, record = [], {}
    with RunLog(__file__) as log:
        log.input(DATES / "manifest.json")
        have = {f["url"] for f in man["files"] if f.get("file")}
        dois = set()
        for gene, qs in queries(man).items():
            record[gene] = []
            for i, q in enumerate(qs, 1):
                cursor, page, hits = "*", 1, []
                while True:
                    params = {"query": q["query"], "format": "json", "resultType": "core",
                              "pageSize": PAGE, "cursorMark": cursor}
                    r = _get(EPMC, params)
                    entry = {"url": r.url, "query": q["query"], "why": q["why"],
                             "http_status": r.status_code, "fetched_utc": utc_now()}
                    if r.status_code != 200:
                        raise SystemExit(f"{q['query']}: HTTP {r.status_code}: stopped (E19.0 rule 4)")
                    path = OUT_EPMC / f"{gene}_{i:02d}_page{page}.json"
                    path.parent.mkdir(parents=True, exist_ok=True)
                    path.write_bytes(r.content)
                    log.output(path)
                    files.append(entry | {"file": rel(path), "bytes": len(r.content),
                                          "sha256": sha256_of(path)})
                    d = r.json()
                    res = d.get("resultList", {}).get("result", [])
                    hits += res
                    nxt = d.get("nextCursorMark")
                    if not res or not nxt or nxt == cursor or len(hits) >= int(d.get("hitCount", 0)):
                        break
                    cursor, page = nxt, page + 1
                record[gene].append({"query": q["query"], "why": q["why"],
                                     "hit_count": len(hits),
                                     "ids": [h.get("id") for h in hits]})
                dois |= {h["doi"].lower() for h in hits if h.get("doi")}
                print(f"{gene}: {q['query']} -> {len(hits)}")
        for doi in sorted(dois):
            url = CROSSREF + urllib.parse.quote(doi, safe="/")
            if url in have:
                continue
            r = _get(url)
            entry = {"url": url, "http_status": r.status_code, "fetched_utc": utc_now()}
            if r.status_code == 404:
                files.append(entry | {"file": None})
                continue
            if r.status_code != 200:
                raise SystemExit(f"{url}: HTTP {r.status_code}: stopped (E19.0 rule 4)")
            path = OUT_CROSSREF / f"{_safe(doi)}.json"
            path.write_bytes(r.content)
            log.output(path)
            files.append(entry | {"file": rel(path), "bytes": len(r.content),
                                  "sha256": sha256_of(path)})
        out = {"what": ("Europe PMC preprint searches for each training gene's SGE study, and "
                        "Crossref's record of every preprint they return, as served"),
               "written_by": rel(Path(__file__)), "europepmc_api": EPMC,
               "crossref_api": CROSSREF, "searches": record, "files": files}
        (DATES / "preprint_search.json").write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n")
        log.output(DATES / "preprint_search.json")
        log.note("searches", {g: [x["hit_count"] for x in v] for g, v in record.items()})


if __name__ == "__main__":
    main()
