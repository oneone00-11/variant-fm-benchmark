"""E18.4: the TP53 score set whose investigator-provided ranges classify E7's TP53 rows,
from MaveDB, as served.

urn:mavedb:00001277-a-1, named by the instruction (E18.4): its metadata, which carries
the investigator-provided score ranges, and its score table, fetched as
scripts/stage_E16/fetch_mavedb.py fetches a score set (the same API, the same requests,
each file saved byte for byte as served and refused at 25 MB or more) into
phase1/data/evidence/external/tp53_raw/:

  scoreset_meta.json                the score set's metadata, as served
  urn_mavedb_00001277-a-1.csv       its score table, as served
  manifest.json                     each file's path, size, sha256, source URL, time of
                                    fetch and API version, with the licence, URN and
                                    row count

Nothing else is fetched (the author's decision of 2026-10-23): not the mapped-variant
records, whose file is 40.4 MB, above the limit, so the cross-check E16 makes with them
is not made; scripts/stage_E18/check_tp53_ranges.py requires instead that every TP53 row
E7 reads maps to its existing coordinates.

Run from the repository root:
    .venv/bin/python scripts/stage_E18/fetch_tp53_ranges.py
"""
from __future__ import annotations

import csv
import io
import json
import sys
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel, sha256_of, utc_now  # noqa: E402

BASE = "https://api.mavedb.org/api/v1"
URN = "urn:mavedb:00001277-a-1"
OUT = REPO / "phase1/data/evidence/external/tp53_raw"
MAX_BYTES = 25 * 1024 * 1024          # E16.2, unchanged

S = requests.Session()
S.trust_env = False


def main() -> None:
    with RunLog(__file__) as log:
        r = S.get(f"{BASE}/api/version", timeout=600)
        r.raise_for_status()
        version = r.json()["version"]
        OUT.mkdir(parents=True, exist_ok=True)
        files = []

        def take(name: str, url: str) -> bytes:
            resp = S.get(url, timeout=600)
            resp.raise_for_status()
            if len(resp.content) >= MAX_BYTES:
                raise SystemExit(f"{name} is {len(resp.content):,} bytes, not under the "
                                 f"{MAX_BYTES:,}-byte limit (E16.2); not saved")
            path = OUT / name
            path.write_bytes(resp.content)
            log.output(path)
            files.append({"path": rel(path), "bytes": len(resp.content),
                          "sha256": sha256_of(path), "source_url": url,
                          "fetched_utc": utc_now(), "api_version": version,
                          "http_status": resp.status_code})
            return resp.content

        meta = json.loads(take("scoreset_meta.json", f"{BASE}/score-sets/{URN}"))
        if meta.get("supersedingScoreSet"):
            raise SystemExit(f"{URN} is superseded by {meta['supersedingScoreSet']}; stopped")
        table = take(f"urn_mavedb_{URN.split(':')[-1]}.csv", f"{BASE}/score-sets/{URN}/scores")
        n_rows = sum(1 for _ in csv.DictReader(io.StringIO(table.decode("utf-8"))))
        lic = meta.get("license") or {}
        manifest = {
            "gene": "TP53", "urn": URN,
            "title": meta.get("title"), "published": meta.get("publishedDate"),
            "license": {k: lic.get(k) for k in ("shortName", "longName", "version", "link")},
            "n_rows": n_rows,
            "api": {"base": BASE, "version": version},
            "selection": "the URN named by the instruction (E18.4)",
            "not_fetched": ("the mapped-variant records (/variant-details, 40.4 MB, above the "
                            "25 MB limit): the author's decision of 2026-10-23"),
            "files": files,
        }
        (OUT / "manifest.json").write_text(json.dumps(manifest, indent=1, ensure_ascii=False) + "\n")
        log.output(OUT / "manifest.json")
        log.note("n_rows", n_rows)
        print(f"TP53: {URN} ({n_rows:,} rows, licence {lic.get('shortName')}); "
              + ", ".join(f"{Path(f['path']).name} {f['bytes']:,} bytes" for f in files))


if __name__ == "__main__":
    main()
