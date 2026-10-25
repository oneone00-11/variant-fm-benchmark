"""E19 C3, BRCA2 (E19.5 item 3, as the author specified it on 2026-10-25): the records a
second BRCA2 dataset is chosen from, saved as served.

A score set is a BRCA2 candidate when its title, a target's name or a target's
identifiers point to BRCA2 (the author's rule, 2026-10-25, item 1):
  * title       the symbol BRCA2 as a word (case folded);
  * target name the symbol BRCA2 as a word; or the whole name, without its version,
                one of BRCA2's identifiers below, or one of HGNC's previous or alias
                symbols for it;
  * identifiers the target's accession (and the gene its accession record names), an
                external identifier (Ensembl, RefSeq, HGNC) or MaveDB's own HGNC
                mapping of the target (mappedHgncName), without version, one of
                BRCA2's identifiers.
BRCA2's identifiers come from three records saved here: HGNC's record of the symbol
(HGNC ID, Ensembl gene, RefSeq accession, MANE Select pair), Ensembl's record of that
gene with its transcripts, and NCBI's list of the gene's RefSeq transcripts. The rule
is applied to every published score set MaveDB serves, not to a search's answer: the
first search (by gene name) missed the score sets whose target is named by its
transcript, ENST00000380152.8.

How every score set is read. MaveDB's search serves at most 100 score sets a request
and orders them by their experiment's title (src/mavedb/lib/score_sets.py, release
v2026.3.2, the release the API reports), a key many score sets share, so pages taken by
offset can repeat a score set and skip another (the first listing did: 00001270-a-1 and
-a-9 came on two pages). The listing is therefore read in parts that each fit in one
request: the text search, which matches a score set's URN among its fields, for each
URN prefix urn:mavedb:0000NNN (ten experiment sets), split down to one experiment set,
then to one experiment, where a part holds more than 100. A part keeps only the score
sets whose own URN carries its prefix (the text search also matches a URN written in
another record's fields), so the parts do not overlap; their union must be exactly
the number of published score sets MaveDB counts, before and after, or the script
stops. The search leaves out a score set a published version has superseded, as the
listing of current versions should.

Saved under phase1/data/evidence/second_set/brca2_raw/, every file listed in
manifest.json with its URL, request body, HTTP status, size, sha256 and time of fetch:
  identifiers/  the three records above
  search/       every part of the listing as served, gzip-compressed without loss (the
                sha256 recorded is that of the bytes as served, and of the .gz beside
                it); the count of published score sets before and after; the
                listing's filter options (every target name and accession MaveDB
                holds, with counts); and three targeted searches (text BRCA2; targets
                BRCA2 and ENST00000380152.8; target accessions ENST00000380152.8 and
                NM_000059.4), kept as cross-checks of the listing
  listing.csv   one row per published score set: the fields the rule reads and which
                of them point to BRCA2
  candidates/   each candidate's score set metadata, its experiment's metadata and its
                score table, as served
  papers/       Crossref's record of every paper MaveDB records on a candidate or its
                experiment, and of Huang 2025, the analysis set's BRCA2 source:
                condition (b) is judged on their author lists

Nothing is decided here; scripts/stage_E19/screen_brca2_second_set.py applies the four
conditions to these files, offline, in the reference image.

Run from the repository root (the network is used here only):
    .venv/bin/python scripts/stage_E19/fetch_brca2_second_set.py
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
import sys
import time
import urllib.parse
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel, sha256_of, utc_now  # noqa: E402

MAVEDB = "https://api.mavedb.org/api/v1"
HGNC = "https://rest.genenames.org/fetch/symbol/{symbol}"
ENSEMBL = "https://rest.ensembl.org/lookup/id/{gene}?expand=1"
NCBI = "https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id/{entrez}/product_report"
CROSSREF = "https://api.crossref.org/works/"
OUT = REPO / "phase1/data/evidence/second_set/brca2_raw"
GENE = "BRCA2"
HUANG_DOI = "10.1038/s41586-024-08388-8"      # Supplementary Table S1's BRCA2 source
PAGE = 100                                     # the most MaveDB's search serves a page
MAX_BYTES = 25 * 1024 * 1024                   # E19.0 rule 4
TARGETED = [("text", {"text": "BRCA2"}),
            ("targets", {"targets": ["BRCA2", "ENST00000380152.8"]}),
            ("target_accessions", {"targetAccessions": ["ENST00000380152.8", "NM_000059.4"]})]
LISTING_COLUMNS = ["urn", "title", "published_date", "target_names", "target_accessions",
                   "accession_genes", "mapped_hgnc_names", "external_identifiers",
                   "title_hit", "name_hit", "identifier_hit", "points_to_brca2"]

S = requests.Session()
S.trust_env = False


def _request(method: str, url: str, *, body=None, headers=None) -> requests.Response:
    for attempt in range(5):
        try:
            r = S.request(method, url, json=body, headers=headers, timeout=600)
            if r.status_code in (429, 500, 502, 503, 504) and attempt < 4:
                time.sleep(2 ** attempt)
                continue
            return r
        except requests.RequestException:
            if attempt == 4:
                raise
            time.sleep(2 ** attempt)
    raise AssertionError("unreachable")


def _save(r: requests.Response, path: Path, files: list, log, *, body=None,
          compress: bool = False) -> None:
    """The response as served, at path (gzip-compressed without loss if asked); a
    response that is not 200 stops the script: no record may be missing."""
    entry = {"url": r.url, "method": r.request.method, "http_status": r.status_code,
             "fetched_utc": utc_now(), "bytes": len(r.content),
             "sha256": hashlib.sha256(r.content).hexdigest()}
    if body is not None:
        entry["request_body"] = body
    if r.status_code != 200:
        raise SystemExit(f"{r.request.method} {r.url}: HTTP {r.status_code}: stopped "
                         "(E19.0 rule 4: a record that cannot be fetched stops the work)")
    if len(r.content) > MAX_BYTES:
        raise SystemExit(f"{r.url}: {len(r.content):,} bytes, above 25 MB: stopped")
    path.parent.mkdir(parents=True, exist_ok=True)
    if compress:
        path = path.with_name(path.name + ".gz")
        path.write_bytes(gzip.compress(r.content, mtime=0))
        entry["stored_as"] = "gzip of the bytes as served"
        entry["sha256_gz"] = sha256_of(path)
    else:
        path.write_bytes(r.content)
    entry["file"] = rel(path)
    files.append(entry)
    log.output(path)


def _bare(x: str | None) -> str:
    return re.sub(r"\.\d+$", "", (x or "").strip()).upper()


def identifiers(files: list, log) -> dict:
    r = _request("GET", HGNC.format(symbol=GENE), headers={"Accept": "application/json"})
    _save(r, OUT / "identifiers" / "hgnc_BRCA2.json", files, log)
    doc = r.json()["response"]["docs"]
    if len(doc) != 1 or doc[0]["symbol"] != GENE:
        raise SystemExit("HGNC's record for BRCA2 is not one record of that symbol")
    doc = doc[0]
    r = _request("GET", ENSEMBL.format(gene=doc["ensembl_gene_id"]),
                 headers={"Content-Type": "application/json"})
    _save(r, OUT / "identifiers" / f"ensembl_{doc['ensembl_gene_id']}.json", files, log)
    ens = r.json()
    r = _request("GET", NCBI.format(entrez=doc["entrez_id"]))
    _save(r, OUT / "identifiers" / f"ncbi_gene_{doc['entrez_id']}_product_report.json",
          files, log)
    ncbi = r.json()["reports"][0]["product"]
    ids = ({doc["hgnc_id"], doc["symbol"], doc["ensembl_gene_id"]}
           | set(doc.get("refseq_accession") or []) | set(doc.get("mane_select") or [])
           | {t["id"] for t in ens.get("Transcript") or []}
           | {t["accession_version"] for t in ncbi.get("transcripts") or []})
    return {"identifiers": sorted({_bare(i) for i in ids}),
            "other_symbols": sorted({s.upper() for s in (doc.get("prev_symbol") or [])
                                     + (doc.get("alias_symbol") or [])}),
            "mane_select": doc.get("mane_select") or []}


WORD = re.compile(rf"(?<![A-Za-z0-9]){GENE}(?![A-Za-z0-9])", re.I)


def listing_row(ss: dict, ids: dict) -> dict:
    ident = set(ids["identifiers"])
    names, accs, acc_genes, mapped, ext = [], [], [], [], []
    for t in ss.get("targetGenes") or []:
        names.append(t.get("name") or "")
        acc = t.get("targetAccession") or {}
        if acc.get("accession"):
            accs.append(acc["accession"])
        if acc.get("gene"):
            acc_genes.append(acc["gene"])
        if t.get("mappedHgncName"):
            mapped.append(t["mappedHgncName"])
        for e in t.get("externalIdentifiers") or []:
            i = e.get("identifier") or {}
            ext.append(f"{i.get('dbName')}:{i.get('identifier')}")
    title_hit = bool(WORD.search(ss.get("title") or ""))
    name_hit = any(WORD.search(n) or _bare(n) in ident or n.strip().upper() in ids["other_symbols"]
                   for n in names)
    id_hit = (any(_bare(a) in ident for a in accs + acc_genes + mapped)
              or any(_bare(x.split(":", 1)[1]) in ident for x in ext))
    return {"urn": ss["urn"], "title": ss.get("title") or "",
            "published_date": ss.get("publishedDate") or "",
            "target_names": ";".join(names), "target_accessions": ";".join(accs),
            "accession_genes": ";".join(acc_genes), "mapped_hgnc_names": ";".join(mapped),
            "external_identifiers": ";".join(ext), "title_hit": title_hit,
            "name_hit": name_hit, "identifier_hit": id_hit,
            "points_to_brca2": title_hit or name_hit or id_hit}


URN = re.compile(r"^urn:mavedb:(\d{8})-([a-z0-9]+)-\d+$")
MAX_PREFIX = 9999999                           # every 8-digit experiment set number


def _count(files: list, log, name: str) -> int:
    body = {"limit": 1}
    r = _request("POST", f"{MAVEDB}/score-sets/search", body=body)
    _save(r, OUT / "search" / f"count_{name}.json", files, log, body=body)
    return int(r.json()["numScoreSets"])


def _part(prefix: str, files: list, log, rows: list, seen: set, ids: dict,
          splits: list) -> None:
    """The score sets whose URN starts with prefix, in one request, or split further."""
    body = {"text": prefix, "limit": PAGE}
    r = _request("POST", f"{MAVEDB}/score-sets/search", body=body)
    if r.status_code != 200:
        raise SystemExit(f"search {body}: HTTP {r.status_code}: stopped")
    page = r.json()
    if page["numScoreSets"] > PAGE:
        splits.append({"prefix": prefix, "num_score_sets": page["numScoreSets"]})
        m = re.fullmatch(r"urn:mavedb:(\d+)", prefix)
        if m and len(m.group(1)) < 8:
            for d in range(10):
                _part(f"{prefix}{d}", files, log, rows, seen, ids, splits)
            return
        if m:                                  # one experiment set: by experiment
            es = _request("GET", f"{MAVEDB}/experiment-sets/{prefix}")
            _save(es, OUT / "search" / f"experiment_set_{prefix.replace(':', '_')}.json",
                  files, log)
            for e in es.json().get("experiments") or []:
                _part(f"{e['urn']}-", files, log, rows, seen, ids, splits)
            return
        raise SystemExit(f"{prefix}: more than {PAGE} score sets in one experiment: stopped")
    safe = prefix.replace(":", "_")
    _save(r, OUT / "search" / f"listing_{safe}.json", files, log, body=body, compress=True)
    for ss in page["scoreSets"]:
        if not ss["urn"].startswith(prefix):
            continue
        if not URN.match(ss["urn"]):
            raise SystemExit(f"{ss['urn']}: not a URN this listing can place: stopped")
        if ss["urn"] in seen:
            raise SystemExit(f"{ss['urn']} in two parts of the listing: stopped")
        seen.add(ss["urn"])
        rows.append(listing_row(ss, ids))


def listing(files: list, ids: dict, log) -> tuple[list[dict], dict]:
    before = _count(files, log, "before")
    rows, seen, splits, n = [], set(), [], 0
    while len(seen) < before and n <= MAX_PREFIX // 10:
        _part(f"urn:mavedb:{n:07d}", files, log, rows, seen, ids, splits)
        n += 1
    after = _count(files, log, "after")
    if not (len(seen) == before == after):
        raise SystemExit(f"the listing gave {len(seen)} score sets; MaveDB counts {before} "
                         f"before and {after} after: stopped")
    rows.sort(key=lambda x: x["urn"])
    return rows, {"published_score_sets_before": before, "published_score_sets_after": after,
                  "listed": len(seen), "parts_read": n, "parts_split": splits}


def main() -> None:
    files: list[dict] = []
    with RunLog(__file__) as log:
        r = _request("GET", f"{MAVEDB}/api/version")
        r.raise_for_status()
        version = r.json()["version"]
        log.note("api_version", version)
        ids = identifiers(files, log)
        rows, how = listing(files, ids, log)
        log.note("listing", how)
        with open(OUT / "listing.csv", "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=LISTING_COLUMNS)
            w.writeheader()
            w.writerows(rows)
        log.output(OUT / "listing.csv")
        body = {}
        r = _request("POST", f"{MAVEDB}/score-sets/search/filter-options", body=body)
        _save(r, OUT / "search" / "filter_options.json", files, log, body=body)
        for name, body in TARGETED:
            r = _request("POST", f"{MAVEDB}/score-sets/search", body=body | {"limit": PAGE})
            _save(r, OUT / "search" / f"search_{name}.json", files, log,
                  body=body | {"limit": PAGE})
        cands = [x["urn"] for x in rows if x["points_to_brca2"]]
        dois = {HUANG_DOI}
        for urn in cands:
            safe = urn.replace(":", "_")
            r = _request("GET", f"{MAVEDB}/score-sets/{urn}")
            _save(r, OUT / "candidates" / f"{safe}.json", files, log)
            meta = r.json()
            exp = (meta.get("experiment") or {}).get("urn")
            r = _request("GET", f"{MAVEDB}/experiments/{exp}")
            _save(r, OUT / "candidates" / f"{exp.replace(':', '_')}.json", files, log)
            exp_meta = r.json()
            r = _request("GET", f"{MAVEDB}/score-sets/{urn}/scores")
            _save(r, OUT / "candidates" / f"{safe}_scores.csv", files, log)
            for rec in (meta, exp_meta):
                for field in ("primaryPublicationIdentifiers", "secondaryPublicationIdentifiers"):
                    for p in rec.get(field) or []:
                        if p.get("doi"):
                            dois.add(p["doi"].lower())
        for doi in sorted(dois):
            r = _request("GET", CROSSREF + urllib.parse.quote(doi, safe="/"))
            _save(r, OUT / "papers" / f"crossref_{doi.replace('/', '_')}.json", files, log)
        manifest = {"what": ("the records a second BRCA2 dataset is chosen from (E19.5 item 3; "
                             "the author's specification of 2026-10-25), as served"),
                    "written_by": rel(Path(__file__)), "mavedb_api": MAVEDB,
                    "mavedb_api_version": version, "rule": __doc__.split("Saved under")[0].strip(),
                    "brca2_identifiers": ids, "listing": how,
                    "candidates": cands, "papers": sorted(dois), "files": files}
        (OUT / "manifest.json").write_text(json.dumps(manifest, indent=1, ensure_ascii=False) + "\n")
        log.output(OUT / "manifest.json")
        log.note("candidates", cands)
        print(f"{len(rows):,} published score sets listed; {len(cands)} point to BRCA2:")
        for x in rows:
            if x["points_to_brca2"]:
                hits = [k for k in ("title_hit", "name_hit", "identifier_hit") if x[k]]
                print(f"  {x['urn']}  {x['title'][:70]!r}  ({', '.join(hits)})")


if __name__ == "__main__":
    main()
