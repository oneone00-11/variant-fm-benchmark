"""E19 C3, BRCA2: the second set's label source under rule R6, fetched as served.

R4 does not apply to Sahu et al. 2025 (paper_rule.json: its classes are cut on a
mixture model's probability, which the score table does not carry), so the source is
the per-variant table published with the paper, as the author allowed (2026-10-25, item
4): Supplementary Table 5, "Function scores for 6,551 SNVs across exons 15 to 26",
taken as Huang 2025's Supplementary Table 3 was taken. Its link is read from the
article page's list of supplementary information by its title, "Supplementary Table 5";
the file is fetched from that link, as served, into the gitignored cache
(phase1/data/e19_cache/brca2_second_set/), and it must be the table the entry describes:
one sheet of 6,551 rows with the columns c.nom, Function.score and Classification. The
article page is not saved (it is the publisher's); its URL, sha256 and time are
recorded.

Written: phase1/data/evidence/second_set/brca2_raw/label_source_manifest.json, in the
label audit's manifest form (scripts/stage_E16c/fetch_sources.py), the cached file
marked "cached", so scripts/reproduce_evidence.py --fetch-inputs fetches it again and
checks its sha256. The compact extract committed in its place is written offline in
the reference image by scripts/stage_E19/extract_brca2_supplementary.py.

Run from the repository root (the network is used here only):
    .venv/bin/python scripts/stage_E19/fetch_brca2_supplementary.py
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel, sha256_of, utc_now  # noqa: E402

sys.path.insert(0, str(REPO / "phase1"))
from src import evid_label_sources as LS  # noqa: E402

ARTICLE = "https://www.nature.com/articles/s41586-024-08349-1"
ENTRY = "Supplementary Table 5"
RAW = REPO / "phase1/data/evidence/second_set/brca2_raw"
CACHE = REPO / "phase1/data/e19_cache/brca2_second_set"
MANIFEST = RAW / "label_source_manifest.json"
SHEET, HEADER_ROW, N_ROWS = "Sheet 1", 1, 6551
NEEDED = ("c.nom", "Function.score", "Classification")


def main() -> None:
    s = requests.Session()
    s.trust_env = False
    with RunLog(__file__) as log:
        r = s.get(ARTICLE, timeout=300)
        r.raise_for_status()
        page = {"url": ARTICLE, "final_url": r.url, "http_status": r.status_code,
                "bytes": len(r.content), "sha256": hashlib.sha256(r.content).hexdigest(),
                "fetched_utc": utc_now(), "saved": False,
                "why_not_saved": "the publisher's page; only its record is kept"}
        # each supplementary item: a link whose text is its title, then its description
        items = re.findall(r'<a [^>]*href="(https://[^"]+_ESM\.[a-z]+)"[^>]*>\s*([^<(]+?)\s*\(download'
                           r'.*?c-article-supplementary__description[^>]*>(.*?)</div>', r.text, re.S)
        hits = [(u, re.sub(r"<[^>]+>", "", d).strip()) for u, title, d in items if title == ENTRY]
        if len(hits) != 1:
            raise SystemExit(f"the article page names {len(hits)} files '{ENTRY}' "
                             f"({len(items)} supplementary items in all): stopped")
        url, description = hits[0]
        page["supplementary_items"] = len(items)
        page["entry"] = {"title": ENTRY, "link": url, "description": description}
        f = s.get(url, timeout=300)
        f.raise_for_status()
        CACHE.mkdir(parents=True, exist_ok=True)
        path = CACHE / url.rsplit("/", 1)[1]
        path.write_bytes(f.content)
        log.output(path)
        sheet = LS.xlsx_sheet(path, SHEET)
        header = sheet[HEADER_ROW]
        if not all(c in header.values() for c in NEEDED) or len(sheet) - 1 != N_ROWS:
            raise SystemExit(f"{path.name}: not the table '{ENTRY}' describes "
                             f"({len(sheet) - 1} rows; header {sorted(header.values())}): stopped")
        entry = {"path": rel(path), "bytes": len(f.content), "sha256": sha256_of(path),
                 "source_url": url, "fetched_utc": utc_now(), "http_status": f.status_code,
                 "cached": True, "article_entry": ENTRY, "sheet": SHEET,
                 "header_row": header}
        manifest = {"gene": "BRCA2", "score_set": "urn:mavedb:00001242-a-1",
                    "what": ("the second BRCA2 set's label source under rule R6 (the author's "
                             "specification of 2026-10-25, item 4), as served"),
                    "written_by": rel(Path(__file__)), "article_page": page,
                    "files": [entry]}
        MANIFEST.write_text(json.dumps(manifest, indent=1, ensure_ascii=False) + "\n")
        log.output(MANIFEST)
        print(f"{ENTRY}: {url} -> {rel(path)} ({len(f.content):,} bytes, sha256 {entry['sha256']})")


if __name__ == "__main__":
    main()
