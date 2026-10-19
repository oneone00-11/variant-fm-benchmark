"""E16: fetch four HAP1 saturation genome editing deposits from MaveDB, as served,
with the evidence for which score set was taken.

The network is used here and in the other fetch_* scripts only; the stage
phase1/src/evid_newgenes.py reads what they write and nothing else.

For each gene this writes, under phase1/data/evidence/external/<gene>_raw/:
  scoreset_meta.json          the score set's metadata, as served
  scoreset_meta_<id>.json     the metadata of each earlier version it supersedes
  experiment_meta.json        the metadata of its experiment
  urn_mavedb_<id>.csv         the score table, byte for byte as served
  score_ranges.json           the investigator-provided calibrations, extracted from
                              scoreset_meta.json without change
  variant_details.ndjson      the mapped-variant records (VRS pre/postMapped), byte for
                              byte as served by /score-sets/{urn}/variant-details
  manifest.json               every file's path, size, sha256, source URL, time of
                              fetch and API version, with the licence, URN and row count

and under scripts/stage_E16/run_log/: each search request with its response, the
candidate list and the reason for the choice, and this script's run record.

Which score set: CTCF and SFPQ are taken by the URN the author verified on 2026-10-19.
RAD51D and XRCC2 are taken by one rule over the search results: the score sets whose
only target is the gene, whose short description names saturation genome editing and
HAP1, and which no other score set supersedes. Exactly one must remain; otherwise the
script stops and prints the candidates.

Run from the repository root:
    .venv/bin/python scripts/stage_E16/fetch_mavedb.py
"""
from __future__ import annotations

import csv
import io
import json
import sys
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import LOG_DIR, REPO, RunLog, rel, sha256_of, utc_now  # noqa: E402

BASE = "https://api.mavedb.org/api/v1"
EXT = REPO / "phase1" / "data" / "evidence" / "external"
GENES = ["CTCF", "SFPQ", "RAD51D", "XRCC2"]
# verified by the author on 2026-10-19; the other two are found by the rule below
GIVEN_URN = {"CTCF": "urn:mavedb:00001262-a-3", "SFPQ": "urn:mavedb:00001265-a-2"}
# a file above this size is not committed as served (E16.2); the script stops instead
MAX_BYTES = 25 * 1024 * 1024

S = requests.Session()
S.trust_env = False


def _get(url: str, **kw) -> requests.Response:
    r = S.get(url, timeout=600, **kw)
    return r


def _id(urn: str) -> str:
    return urn.split(":")[-1]


def _save(path: Path, blob: bytes, log: RunLog) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(blob)
    log.output(path)
    if len(blob) > MAX_BYTES:
        raise SystemExit(f"{rel(path)} is {len(blob):,} bytes, above the {MAX_BYTES:,}-byte "
                         "limit for committing a file as served; stopped for the author to "
                         "decide how to store it (E16.2)")


def api_version() -> str:
    r = _get(f"{BASE}/api/version")
    r.raise_for_status()
    return r.json()["version"]


def search(gene: str, stamp: str, log: RunLog) -> list[dict]:
    """Both search forms, each request and response kept as served."""
    found = {}
    for mode, body in (("text", {"text": gene}), ("targets", {"targets": [gene]})):
        r = S.post(f"{BASE}/score-sets/search", json=body, timeout=300)
        stem = LOG_DIR / "searches" / f"search_{gene}_{mode}_{stamp}"
        _save(stem.with_suffix(".request.json"),
              (json.dumps({"method": "POST", "url": f"{BASE}/score-sets/search",
                           "body": body, "status": r.status_code,
                           "fetched_utc": utc_now()}, indent=1) + "\n").encode(), log)
        _save(stem.with_suffix(".response.json"), r.content, log)
        r.raise_for_status()
        for item in r.json().get("scoreSets", []):
            found[item["urn"]] = item
    return list(found.values())


def _summary(item: dict) -> dict:
    sup = item.get("supersedingScoreSet")
    return {
        "urn": item["urn"],
        "title": item.get("title"),
        "short_description": item.get("shortDescription"),
        "targets": [t.get("name") for t in item.get("targetGenes", [])],
        "num_variants": item.get("numVariants"),
        "published_date": item.get("publishedDate"),
        "superseding_score_set": sup.get("urn") if isinstance(sup, dict) else sup,
    }


def choose(gene: str, candidates: list[dict]) -> tuple[str, dict]:
    summ = [_summary(c) for c in candidates]
    if gene in GIVEN_URN:
        urn = GIVEN_URN[gene]
        return urn, {"how": "URN given by the author (verified 2026-10-19)",
                     "given_urn_in_search_results": any(s["urn"] == urn for s in summ),
                     "candidates": summ}
    def keep(s: dict) -> tuple[bool, list[str]]:
        why = []
        if s["targets"] != [gene]:
            why.append(f"targets are {s['targets']}, not only {gene}")
        text = (s["short_description"] or "").lower()
        if "saturation genome editing" not in text:
            why.append("short description does not name saturation genome editing")
        if "hap1" not in text:
            why.append("short description does not name HAP1")
        if s["superseding_score_set"]:
            why.append(f"superseded by {s['superseding_score_set']}")
        return not why, why
    verdicts = [(s, *keep(s)) for s in summ]
    kept = [s for s, ok, _ in verdicts if ok]
    record = {"how": ("rule: only target is the gene; short description names saturation "
                      "genome editing and HAP1; not superseded"),
              "candidates": [s | {"kept": ok, "excluded_because": why}
                             for s, ok, why in verdicts]}
    if len(kept) != 1:
        print(json.dumps(record, indent=1))
        raise SystemExit(f"{gene}: {len(kept)} score sets meet the rule; stopped (E16.10)")
    return kept[0]["urn"], record


def fetch_gene(gene: str, version: str, stamp: str, log: RunLog) -> dict:
    out = EXT / f"{gene.lower()}_raw"
    files = []

    def take(name: str, url: str, extra: dict | None = None, **kw) -> bytes:
        r = _get(url, **kw)
        r.raise_for_status()
        path = out / name
        _save(path, r.content, log)
        entry = {"path": rel(path), "bytes": len(r.content), "sha256": sha256_of(path),
                 "source_url": url, "fetched_utc": utc_now(), "api_version": version,
                 "http_status": r.status_code}
        if extra:
            entry.update(extra(r) if callable(extra) else extra)
        files.append(entry)
        return r.content

    candidates = search(gene, stamp, log)
    urn, choice = choose(gene, candidates)
    choice_path = LOG_DIR / "searches" / f"choice_{gene}_{stamp}.json"
    _save(choice_path, (json.dumps({"gene": gene, "chosen": urn} | choice, indent=1,
                                   ensure_ascii=False) + "\n").encode(), log)

    meta = json.loads(take("scoreset_meta.json", f"{BASE}/score-sets/{urn}"))
    if meta.get("supersedingScoreSet"):
        raise SystemExit(f"{gene}: {urn} is superseded by "
                         f"{meta['supersedingScoreSet']}; stopped (E16.10)")
    chain, prev = [], meta.get("supersededScoreSet")
    while prev:                                   # every earlier version, newest first
        purn = prev["urn"]
        pmeta = json.loads(take(f"scoreset_meta_{_id(purn)}.json",
                                f"{BASE}/score-sets/{purn}"))
        chain.append(purn)
        prev = pmeta.get("supersededScoreSet")
    exp_urn = meta["experiment"]["urn"]
    take("experiment_meta.json", f"{BASE}/experiments/{exp_urn}")

    scores = take(f"urn_mavedb_{_id(urn)}.csv", f"{BASE}/score-sets/{urn}/scores")
    n_rows = sum(1 for _ in csv.DictReader(io.StringIO(scores.decode("utf-8"))))

    cals = meta.get("scoreCalibrations") or []
    legacy = meta.get("scoreRanges")
    ranges = {
        "extracted_from": rel(out / "scoreset_meta.json"),
        "key": "scoreCalibrations",
        "selection": "investigatorProvided is true",
        "all_calibration_titles": [c.get("title") for c in cals],
        "legacy_scoreRanges_key_present": legacy is not None,
        "investigator_provided": [c for c in cals if c.get("investigatorProvided") is True],
    }
    path = out / "score_ranges.json"
    _save(path, (json.dumps(ranges, indent=1, ensure_ascii=False) + "\n").encode(), log)
    files.append({"path": rel(path), "bytes": path.stat().st_size, "sha256": sha256_of(path),
                  "source_url": f"{BASE}/score-sets/{urn}",
                  "derived_from": rel(out / "scoreset_meta.json"),
                  "note": "extract of scoreCalibrations, entries unchanged"})

    take("variant_details.ndjson", f"{BASE}/score-sets/{urn}/variant-details",
         extra=lambda r: {"response_headers": {k: r.headers.get(k) for k in (
             "X-As-Of", "X-Total-Count", "X-Stream-Type", "X-Processing-Started")}})

    lic = meta.get("license") or {}
    manifest = {
        "gene": gene, "urn": urn, "experiment_urn": exp_urn,
        "superseded_versions": chain,
        "license": {k: lic.get(k) for k in ("shortName", "longName", "version", "link")},
        "n_rows": n_rows,
        "api": {"base": BASE, "version": version},
        "selection": {"how": choice["how"], "record": rel(choice_path)},
        "files": files,
    }
    mpath = out / "manifest.json"
    mpath.write_text(json.dumps(manifest, indent=1, ensure_ascii=False) + "\n")
    log.output(mpath)
    print(f"{gene}: {urn} ({n_rows:,} rows, licence {lic.get('shortName')}); "
          f"{len(chain)} earlier version(s): {chain or 'none'}")
    return {"gene": gene, "urn": urn, "n_rows": n_rows, "superseded_versions": chain}


def main() -> None:
    with RunLog(__file__) as log:
        stamp = utc_now().replace(":", "").replace("-", "").replace("+0000", "Z")
        version = api_version()
        log.note("api_version", version)
        # the route the earlier probes used for postMapped records, kept as evidence
        # of why this script reads /variant-details instead
        old = _get(f"{BASE}/score-sets/{GIVEN_URN['CTCF']}/mapped-variants")
        evidence = LOG_DIR / "searches" / f"mapped_variants_route_{stamp}.json"
        _save(evidence, (json.dumps({"url": old.url, "status": old.status_code,
                                     "body": old.text, "fetched_utc": utc_now()},
                                    indent=1) + "\n").encode(), log)
        log.note("mapped_variants_route_status", old.status_code)
        log.note("genes", [fetch_gene(g, version, stamp, log) for g in GENES])


if __name__ == "__main__":
    main()
