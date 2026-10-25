"""E19 C1 (E19.4): the records each training gene's first public date is read from.

The low-circularity subset keeps an assertion only if ClinVar last evaluated it before
the gene's SGE results were first public. For each of the seven training genes that
date is the earliest of three, each from a record saved here as served (E19.4):

  * the earliest public version of the score set the analysis cites (Supplementary
    Table S1's MaveDB accession), found by walking MaveDB's supersededScoreSet chain
    back to its first version: each version's metadata and its experiment's;
  * the journal article's earliest online date, from Crossref, for every DOI MaveDB
    records on those score sets and experiments (their publication identifiers'
    DOIs and their own doiIdentifiers);
  * a preprint's posting date, from Crossref, for a preprint DOI MaveDB records or a
    "has-preprint" relation in a Crossref record lists.

This script only fetches and saves, into phase1/data/evidence/e19/first_public_dates/
(mavedb/, crossref/ and a manifest of every file's URL, HTTP status, size, sha256 and
time of fetch); the dates are read from the saved records by the C1 stage, offline.
A gene for which no record gives any date stops here, as E19.12 says; nothing is
inferred or filled in from memory.

Run from the repository root (the network is used here only):
    .venv/bin/python scripts/stage_E19/fetch_first_public_dates.py
"""
from __future__ import annotations

import json
import sys
import time
import urllib.parse
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel, sha256_of, utc_now  # noqa: E402

MAVEDB = "https://api.mavedb.org/api/v1"
CROSSREF = "https://api.crossref.org/works/"
OUT = REPO / "phase1/data/evidence/e19/first_public_dates"
S1 = REPO / "phase1/reports/evidence/supplement/tableS1_functional_standards.csv"
SEVEN = ["BAP1", "BARD1", "BRCA1", "BRCA2", "PALB2", "RAD51C", "VHL"]
PREPRINT_PREFIXES = ("10.1101/",)      # bioRxiv and medRxiv, Cold Spring Harbor's prefix

S = requests.Session()
S.trust_env = False


def _safe(name: str) -> str:
    return urllib.parse.quote(name, safe="").replace("%", "_")


def _get(url: str, path: Path, files: list, log) -> dict | None:
    """GET url, save the body as served at path; 404 is recorded, not raised."""
    for attempt in range(5):
        try:
            r = S.get(url, timeout=120)
            break
        except requests.RequestException as e:
            if attempt == 4:
                raise
            time.sleep(2 ** attempt)
    entry = {"url": url, "http_status": r.status_code, "fetched_utc": utc_now()}
    if r.status_code == 404:
        files.append(entry | {"file": None})
        return None
    r.raise_for_status()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(r.content)
    log.output(path)
    files.append(entry | {"file": rel(path), "bytes": len(r.content), "sha256": sha256_of(path)})
    return r.json()


def _dois(record: dict) -> list[tuple[str, str]]:
    """(doi, where MaveDB records it) for a score set's or an experiment's record."""
    out = []
    for field in ("primaryPublicationIdentifiers", "secondaryPublicationIdentifiers"):
        for p in record.get(field) or []:
            doi = p.get("doi")
            if not doi and str(p.get("dbName", "")).lower() in ("biorxiv", "medrxiv"):
                doi = p.get("identifier")
            if doi:
                out.append((doi.lower(), f"{field} ({p.get('dbName')} {p.get('identifier')})"))
    for d in record.get("doiIdentifiers") or []:
        if d.get("identifier"):
            out.append((d["identifier"].lower(), "doiIdentifiers"))
    return out


def main() -> None:
    s1 = pd.read_csv(S1, dtype=str, keep_default_na=False).set_index("Gene")
    files: list[dict] = []
    with RunLog(__file__) as log:
        log.input(S1)
        found = {}
        for gene in SEVEN:
            urn = s1.loc[gene, "MaveDB accession"].strip()
            chain, dois, urn_now, missing = [], [], urn, []
            while urn_now:
                ss = _get(f"{MAVEDB}/score-sets/{urn_now}",
                          OUT / "mavedb" / f"{_safe(urn_now)}.json", files, log)
                if ss is None:
                    if urn_now == urn:
                        raise SystemExit(f"{gene}: MaveDB has no score set {urn}: stopped")
                    # an earlier version the chain names but MaveDB does not serve: the
                    # chain's earliest version is then unknown, which the stage reports
                    missing.append(urn_now)
                    break
                chain.append(urn_now)
                dois += [(d, f"{urn_now} {w}") for d, w in _dois(ss)]
                exp = (ss.get("experiment") or {}).get("urn")
                if exp:
                    er = _get(f"{MAVEDB}/experiments/{exp}",
                              OUT / "mavedb" / f"{_safe(exp)}.json", files, log)
                    if er is not None:
                        dois += [(d, f"{exp} {w}") for d, w in _dois(er)]
                urn_now = (ss.get("supersededScoreSet") or {}).get("urn")
            seen, queue = set(), [d for d, _ in dois]
            while queue:
                doi = queue.pop(0)
                if doi in seen:
                    continue
                seen.add(doi)
                cr = _get(CROSSREF + urllib.parse.quote(doi, safe="/"),
                          OUT / "crossref" / f"{_safe(doi)}.json", files, log)
                if cr is None:
                    continue
                rel_ = (cr.get("message") or {}).get("relation") or {}
                for kind in ("has-preprint",):
                    for x in rel_.get(kind) or []:
                        if x.get("id-type") == "doi" and x.get("id"):
                            queue.append(x["id"].lower())
            found[gene] = {"cited_urn": urn, "versions": chain,
                           "versions_not_served": missing, "dois": sorted(seen),
                           "doi_sources": dois}
            print(f"{gene}: {len(chain)} version(s) {chain}"
                  + (f", not served: {missing}" if missing else "") + f"; DOIs {sorted(seen)}")
        OUT.mkdir(parents=True, exist_ok=True)
        manifest = {"what": "the records E19.4's first public dates are read from, as served",
                    "mavedb_api": MAVEDB, "crossref_api": CROSSREF, "genes": found,
                    "files": files}
        (OUT / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")
        log.output(OUT / "manifest.json")
        log.note("genes", {g: {k: v[k] for k in ("versions", "dois")} for g, v in found.items()})


if __name__ == "__main__":
    main()
