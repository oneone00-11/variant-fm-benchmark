"""E19 C3, BRCA2: check every quote of the second set's paper_rule.json against the paper.

The paper, Sahu et al. 2025 (doi:10.1038/s41586-024-08349-1), is not open access: the
author obtained the publisher's PDF through the NUS library on 2026-10-26, and it is
kept in the gitignored cache (phase1/data/e19_cache/brca2_second_set/), never
committed (the author's review of 2026-10-25, item 3). This script checks that the
cached file is the one paper_rule.json records (its sha256 and size), extracts its text
with pypdf (the version the reference image carries), and checks that every quote
occurs in it, with Unicode compatibility forms, line-end hyphenation and white space
folded on both sides, nothing else. It writes paper_quote_check.json beside
paper_rule.json: the PDF's sha256 and size, pypdf's version, the sha256 of the
paper_rule.json it checked, and per quote whether it was found. The PDF's text is not
saved.

    docker run ... IMAGE bash scripts/stage_E19/in_image.sh scripts/stage_E19/check_brca2_paper_quotes.py
"""
from __future__ import annotations

import json
import re
import sys
import unicodedata
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel, sha256_of, utc_now  # noqa: E402

RAW = REPO / "phase1/data/evidence/second_set/brca2_raw"
RULE = RAW / "paper_rule.json"
OUT = RAW / "paper_quote_check.json"
FOLDING = ("Unicode NFKC; a hyphen at a line end joined to a following lower-case letter; all "
           "white space removed; the same on both sides")


def _fold(s: str) -> str:
    s = unicodedata.normalize("NFKC", s)
    s = re.sub(r"-\n(?=[a-z])", "", s)
    return re.sub(r"\s+", "", s)


def main() -> None:
    import pypdf
    with RunLog(__file__) as log:
        log.input(RULE)
        rule = json.loads(RULE.read_text())
        pdf = rule["source"]["pdf"]
        path = REPO / pdf["cache_path"]
        if not path.is_file():
            raise SystemExit(f"{pdf['cache_path']} is not here: the author's copy is needed")
        if sha256_of(path) != pdf["sha256"] or path.stat().st_size != pdf["bytes"]:
            raise SystemExit(f"{pdf['cache_path']} is not the file paper_rule.json records: stopped")
        reader = pypdf.PdfReader(str(path))
        text = _fold("\n".join(p.extract_text() or "" for p in reader.pages))
        quotes = [{"where": q["where"], "found": _fold(q["text"]) in text,
                   "first_words": " ".join(q["text"].split()[:8])} for q in rule["quotes"]]
        record = {"checked_utc": utc_now(), "written_by": rel(Path(__file__)),
                  "paper_rule": rel(RULE), "paper_rule_sha256": sha256_of(RULE),
                  "pdf": {"cache_path": pdf["cache_path"], "sha256": pdf["sha256"],
                          "bytes": pdf["bytes"], "pages": len(reader.pages)},
                  "text_extraction": f"pypdf {pypdf.__version__}", "folding": FOLDING,
                  "quotes": quotes}
        OUT.write_text(json.dumps(record, indent=1, ensure_ascii=False) + "\n")
        log.output(OUT)
        n = sum(q["found"] for q in quotes)
        print(f"{n}/{len(quotes)} quotes found in {pdf['cache_path']}")
        if n != len(quotes):
            raise SystemExit("a quote was not found in the paper: fix paper_rule.json")


if __name__ == "__main__":
    main()
