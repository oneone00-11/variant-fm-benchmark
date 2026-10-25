"""E19 C1: ClinVar's own documentation of review status and its gold stars, as served.

C1 grades the analysis set's ClinVar assertions by stars (main analysis at least one,
sensitivity analysis at least two), converted from CLNREVSTAT by the table ClinVar's
documentation gives (E19.3: "换算表照 ClinVar 官方文档并入库"). ClinVar's FTP README
and XML schemas point to that page for the stars and carry no table of their own, so
the page is the source. This saves it byte for byte, with a manifest of its URL, size,
sha256 and time, into phase1/data/evidence/e19/clinvar_review_status/; the stage
evid_suppl_c1 reads the saved page offline and writes the conversion table it finds
there. Nothing here is typed from memory.

Two ways in, and a page is kept only if it is the documentation: its title must be
"Review status in ClinVar", it must hold the table headed "Number of gold stars" and
name every review status that occurs in the analysis set's CLNREVSTAT, written as
ClinVar's pages write them (underscores as spaces), and it must carry no reCAPTCHA;
anything else, such as a "Checking your browser" challenge, is refused and not saved.

  * fetched by this script (the default). The first attempt, on 2026-10-25, was served
    NCBI's reCAPTCHA challenge page instead of the documentation; its run record is
    kept. The challenge is not worked around.
  * --from-file PATH: the page as the author saved it from a browser ("Webpage, HTML
    only"), by the author's decision of 2026-10-25 (saved with Chrome on 2026-10-26);
    the manifest says so ("saved by the author with a browser") and records the name it
    was saved under, the file's own modification time (UTC and local), the URL macOS
    records it came from (the kMDItemWhereFroms attribute) and the browser its
    quarantine attribute names, beside the time it was recorded here.

Run from the repository root:
    .venv/bin/python scripts/stage_E19/fetch_clinvar_review_status.py [--from-file PATH]
"""
from __future__ import annotations

import argparse
import json
import plistlib
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel, sha256_of, utc_now  # noqa: E402

URL = "https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/"
OUT = REPO / "phase1/data/evidence/e19/clinvar_review_status"
MAX_BYTES = 25 * 1024 * 1024          # E16.2 and E19.0 rule 4: nothing larger is committed
SET = REPO / "phase1/data/evidence/analysis_set_v1.parquet"


def review_statuses() -> list[str]:
    """The analysis set's CLNREVSTAT values, as ClinVar's pages write them."""
    import pandas as pd
    v = pd.read_parquet(SET, columns=["clnrevstat"])["clnrevstat"]
    return sorted({x.replace("_", " ") for x in v.dropna() if x})


def _is_documentation(data: bytes) -> bool:
    raw = data.decode("utf-8", errors="replace")
    title = re.search(r"<title>(.*?)</title>", raw, re.S | re.I)
    text = " ".join(raw.lower().split())
    return (title is not None and " ".join(title.group(1).split()) == "Review status in ClinVar"
            and "number of gold stars" in text and "recaptcha" not in text
            and all(p.lower() in text for p in review_statuses()))


def _xattr(path: Path, name: str) -> bytes | None:
    r = subprocess.run(["xattr", "-px", name, str(path)], capture_output=True, text=True)
    return bytes.fromhex("".join(r.stdout.split())) if r.returncode == 0 else None


def _saved_from(path: Path) -> dict:
    """What macOS recorded about the saved file: where it came from, which browser."""
    where = _xattr(path, "com.apple.metadata:kMDItemWhereFroms")
    quarantine = _xattr(path, "com.apple.quarantine")
    return {"where_from": plistlib.loads(where) if where else "not recorded",
            "browser": (quarantine.decode().split(";")[2] if quarantine and
                        quarantine.decode().count(";") >= 2 else "not recorded")}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--from-file", type=Path, default=None,
                    help="the page as the author saved it from a browser")
    args = ap.parse_args()
    with RunLog(__file__) as log:
        if args.from_file is None:
            s = requests.Session()
            s.trust_env = False
            r = s.get(URL, timeout=600)
            r.raise_for_status()
            data = r.content
            how = {"how": "fetched by this script", "final_url": r.url,
                   "http_status": r.status_code, "content_type": r.headers.get("Content-Type"),
                   "fetched_utc": utc_now()}
        else:
            src = args.from_file.expanduser()
            log.input(src)
            data = src.read_bytes()
            mtime = datetime.fromtimestamp(src.stat().st_mtime, timezone.utc)
            how = {"how": ("saved by the author with a browser (作者用浏览器保存), as "
                           "'Webpage, HTML only', the scripted fetch having been served a "
                           "reCAPTCHA challenge (the author's decision of 2026-10-25)"),
                   "saved_as": src.name,
                   "file_saved_utc": mtime.isoformat(timespec="seconds"),
                   "file_saved_local": mtime.astimezone().isoformat(timespec="seconds"),
                   **_saved_from(src), "recorded_utc": utc_now()}
        log.note("how", how)
        if not _is_documentation(data):
            raise SystemExit(f"{len(data):,} bytes that are not the review-status documentation "
                             "(a phrase of the star table is missing): not saved")
        if len(data) >= MAX_BYTES:
            raise SystemExit(f"{len(data):,} bytes, not under the 25 MB limit: not saved")
        OUT.mkdir(parents=True, exist_ok=True)
        page = OUT / "review_status.html"
        page.write_bytes(data)
        log.output(page)
        manifest = {"source_url": URL, **how, "file": rel(page), "bytes": len(data),
                    "sha256": sha256_of(page)}
        (OUT / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")
        log.output(OUT / "manifest.json")
        log.note("manifest", manifest)
        print(f"saved {rel(page)}: {len(data):,} bytes, sha256 {manifest['sha256']}")


if __name__ == "__main__":
    main()
