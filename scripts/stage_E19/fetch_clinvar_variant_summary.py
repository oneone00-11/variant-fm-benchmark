"""E19 C1 (E19.4): ClinVar's monthly variant_summary archive for June 2026, as served.

The low-circularity subset needs, for each ClinVar assertion, the date its aggregate
germline classification was last evaluated. The VCF the analysis set was built from
(2026-06-15) carries no such date; the tab-delimited monthly archive does, in its
LastEvaluated column (E19.4: variant_summary_2026-06.txt.gz in tab_delimited/archive/,
joined on VariationID).

The archive is larger than 25 MB, so it is not committed and CI does not fetch it
(E19.0 rule 4): it is downloaded once into the gitignored phase1/data/e19_cache/, through
scripts/reproduce_evidence.py's own resumable download, and checked against the md5 NCBI
publishes beside it when there is one. Its URL, size, sha256, md5, Last-Modified and
time of fetch go into phase1/data/evidence/e19/clinvar_variant_summary/manifest.json,
which is committed. scripts/stage_E19/extract_clinvar_variant_summary.py then writes,
in the reference image, the small table of the analysis set's rows that is committed in
its place.

Run from the repository root (the network is used here only):
    .venv/bin/python scripts/stage_E19/fetch_clinvar_variant_summary.py
"""
from __future__ import annotations

import gzip
import hashlib
import json
import sys
import urllib.error
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel, sha256_of, utc_now  # noqa: E402

sys.path.insert(0, str(REPO / "scripts"))
import reproduce_evidence as R  # noqa: E402

NAME = "variant_summary_2026-06.txt.gz"
URL = f"https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/{NAME}"
CACHE = REPO / "phase1/data/e19_cache" / NAME
OUT = REPO / "phase1/data/evidence/e19/clinvar_variant_summary"


def _head(url: str) -> dict:
    req = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(req, timeout=120) as r:
        return {"http_status": r.status, "content_length": r.headers.get("Content-Length"),
                "last_modified": r.headers.get("Last-Modified")}


def _md5_published(url: str) -> str | None:
    try:
        with urllib.request.urlopen(url + ".md5", timeout=120) as r:
            return r.read().decode("ascii", errors="replace").split()[0]
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return None
        raise


def main() -> None:
    with RunLog(__file__) as log:
        head = _head(URL)
        log.note("head", head)
        previous = (json.loads((OUT / "manifest.json").read_text())
                    if (OUT / "manifest.json").exists() else {})
        if not CACHE.exists():
            R._download(URL, CACHE)
            fetched = utc_now()
        else:
            fetched = previous.get("fetched_utc", "before this run; not recorded")
        size = CACHE.stat().st_size
        if head["content_length"] and int(head["content_length"]) != size:
            raise SystemExit(f"{rel(CACHE)} is {size:,} bytes; the server announces "
                             f"{head['content_length']}")
        md5 = hashlib.md5()
        with open(CACHE, "rb") as f:
            for block in iter(lambda: f.read(1 << 20), b""):
                md5.update(block)
        published = _md5_published(URL)
        if published is not None and published != md5.hexdigest():
            raise SystemExit(f"{rel(CACHE)}: md5 {md5.hexdigest()} is not NCBI's {published}")
        with gzip.open(CACHE, "rt", encoding="utf-8", errors="replace") as f:
            header = f.readline().rstrip("\n").split("\t")
        if "LastEvaluated" not in header:
            raise SystemExit(f"{NAME} has no LastEvaluated column: stopped (E19.12)")
        OUT.mkdir(parents=True, exist_ok=True)
        manifest = {
            "source_url": URL, "file": NAME, "cached_at": rel(CACHE), "committed": False,
            "why_not_committed": "larger than 25 MB (E19.0 rule 4)",
            "bytes": size, "sha256": sha256_of(CACHE), "md5": md5.hexdigest(),
            "md5_published": published if published is not None else "not_provided",
            "last_modified": head["last_modified"], "fetched_utc": fetched,
            "checked_utc": utc_now(),
            "columns": header,
        }
        (OUT / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")
        log.output(OUT / "manifest.json")
        log.note("manifest", {k: v for k, v in manifest.items() if k != "columns"})
        print(f"{NAME}: {size:,} bytes, sha256 {manifest['sha256']}, md5 "
              f"{'= NCBI' if published else 'not published by NCBI'}")


if __name__ == "__main__":
    main()
