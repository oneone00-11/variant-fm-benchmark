"""E16b: fetch the IGVF portal copy of the SBDS saturation genome editing scores, as
served, for a cross-check only (E16b.5): the SBDS labels come from MaveDB, never from
this file (decided by the author on 2026-10-20, before any predictor result).

Writes, under phase1/data/evidence/external/igvf_sbds/:
  IGVFFI9790OCYB.tsv.gz    the score file, byte for byte as served (anonymous download)
  IGVFFI9790OCYB.json      the portal record of that file
  IGVFDS8094DYWC.json      the portal record of its analysis set
  manifest.json            each file's URL, time of fetch, size and sha256, the IGVF
                           accessions and the release dates the records give

A download that fails is recorded in the manifest and the script goes on: the
cross-check is then listed as pending (E16b.2), and nothing else depends on it.

Run from the repository root:
    .venv/bin/python scripts/stage_E16/fetch_igvf.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel, sha256_of, utc_now  # noqa: E402

API = "https://api.data.igvf.org"
OUT = REPO / "phase1" / "data" / "evidence" / "external" / "igvf_sbds"
FILE_ACC, SET_ACC = "IGVFFI9790OCYB", "IGVFDS8094DYWC"   # from the author's SGE survey
S = requests.Session()
S.trust_env = False


def take(url: str, name: str, log: RunLog, files: list, failures: list) -> bytes | None:
    try:
        r = S.get(url, timeout=300, headers={"Accept": "application/json"} if name.endswith(".json") else {})
        r.raise_for_status()
    except requests.RequestException as e:
        failures.append({"url": url, "error": str(e), "fetched_utc": utc_now()})
        print(f"  not fetched: {url} ({e})")
        return None
    path = OUT / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(r.content)
    log.output(path)
    files.append({"path": rel(path), "bytes": len(r.content), "sha256": sha256_of(path),
                  "source_url": url, "fetched_utc": utc_now(), "http_status": r.status_code})
    return r.content


def main() -> None:
    with RunLog(__file__) as log:
        files, failures = [], []
        frec = take(f"{API}/tabular-files/{FILE_ACC}/?format=json", f"{FILE_ACC}.json", log, files, failures)
        srec = take(f"{API}/analysis-sets/{SET_ACC}/?format=json", f"{SET_ACC}.json", log, files, failures)
        take(f"{API}/tabular-files/{FILE_ACC}/@@download/{FILE_ACC}.tsv.gz", f"{FILE_ACC}.tsv.gz",
             log, files, failures)
        f = json.loads(frec) if frec else {}
        s = json.loads(srec) if srec else {}
        manifest = {
            "use": "cross-check only (E16b.5); never used to build labels (author, 2026-10-20)",
            "file_accession": FILE_ACC, "analysis_set_accession": SET_ACC,
            "file_status": f.get("status", "not_provided"),
            "file_release_timestamp": f.get("release_timestamp", "not_provided"),
            "file_md5sum": f.get("md5sum", "not_provided"),
            "analysis_set_status": s.get("status", "not_provided"),
            "analysis_set_release_timestamp": s.get("release_timestamp", "not_provided"),
            "files": files, "failures": failures,
        }
        OUT.mkdir(parents=True, exist_ok=True)
        (OUT / "manifest.json").write_text(json.dumps(manifest, indent=1, ensure_ascii=False) + "\n")
        log.output(OUT / "manifest.json")
        log.note("summary", {k: manifest[k] for k in manifest if k not in ("files",)})
        print(json.dumps({k: manifest[k] for k in manifest if k != "files"}, indent=1))


if __name__ == "__main__":
    main()
