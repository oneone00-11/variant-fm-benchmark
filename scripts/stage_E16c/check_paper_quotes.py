"""E16c: check, against each cited paper as served, the rule label_audit/<gene>_raw/
paper_rule.json quotes from it (candidate rule R4 of the label source audit).

paper_rule.json is written by hand: the passages of the paper that state how each
variant was classed, copied verbatim with their place, and, where the passages give a
rule a program can apply, that rule in machine form. This script makes the quotes
checkable. For each gene with a paper_rule.json it
  - resolves the paper's PMID, PMCID and DOI through NCBI's ID converter and stops if
    they disagree with the file;
  - fetches the full text (Europe PMC's fullTextXML, else NCBI's PMC efetch) and
    checks that every quote occurs in it, whitespace and Unicode compatibility forms
    folded, nothing else;
  - writes label_audit/<gene>_raw/paper_quote_check.json: the URLs, fetch times, the
    sha256 and size of each document as served, the sha256 of the paper_rule.json it
    checked, and, per quote, whether it was found.
The papers themselves are not saved (their licences differ); the sha256 identifies
the document that was read. The audit stage refuses a paper_rule.json that this check
has not passed.

Run from the repository root:
    .venv/bin/python scripts/stage_E16c/check_paper_quotes.py
"""
from __future__ import annotations

import hashlib
import html
import json
import re
import sys
import unicodedata
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel, sha256_of, utc_now  # noqa: E402

AUDIT = REPO / "phase1" / "data" / "evidence" / "label_audit"
IDCONV = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids={id}&format=json"
EUROPEPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id={num}"
S = requests.Session()
S.trust_env = False


# Elements that end a run of text in JATS (a paragraph, a title, a table cell...). Any
# other tag is inline markup (italic, bold, sub, sup, a cross-reference) and is dropped
# without a space, so "<italic>P</italic><sub>nf</sub>" reads "Pnf" on both sides.
BLOCK = re.compile(r"</?(?:p|sec|title|caption|label|td|th|tr|table|table-wrap|fig|list|"
                   r"list-item|abstract|body|back|front|ref|ref-list|article-title|"
                   r"disp-quote|boxed-text|fn|fn-group|def-list|def|term|statement)\b[^>]*>",
                   re.I)


def _text(markup: str) -> str:
    """Text with block elements as spaces, inline markup dropped, entities resolved,
    and Unicode compatibility forms (thin and no-break spaces among them) and
    whitespace folded. A quote, typed with <i>, <sub> and <sup>, gets the same."""
    t = re.sub(r"<[^>]+>", "", BLOCK.sub(" ", markup))
    return " ".join(unicodedata.normalize("NFKC", html.unescape(t)).split())


def _quote_text(quote: str) -> str:
    """A quote's text: its own <i>, <sub> and <sup> are markup; any other '<' is the
    character itself (the documents write it &lt;)."""
    return _text(re.sub(r"<(?!/?(?:i|sub|sup)>)", "&lt;", quote))


def ids(pmid: str) -> dict:
    r = S.get(IDCONV.format(id=pmid), timeout=120)
    r.raise_for_status()
    rec = r.json()["records"][0]
    return {"url": r.url, "fetched_utc": utc_now(), "pmid": rec.get("pmid"),
            "pmcid": rec.get("pmcid"), "doi": rec.get("doi")}


def full_text(pmcid: str) -> dict:
    tried = []
    for url in (EUROPEPMC.format(pmcid=pmcid), EFETCH.format(num=pmcid.removeprefix("PMC"))):
        r = S.get(url, timeout=300)
        tried.append({"url": url, "http_status": r.status_code, "bytes": len(r.content),
                      "fetched_utc": utc_now()})
        if r.status_code == 200 and "<body" in r.text:
            return {"url": url, "fetched_utc": tried[-1]["fetched_utc"], "bytes": len(r.content),
                    "sha256": hashlib.sha256(r.content).hexdigest(), "text": _text(r.text),
                    "tried": tried}
    return {"url": None, "tried": tried, "text": None}


def check(gene_dir: Path, log: RunLog) -> dict:
    path = gene_dir / "paper_rule.json"
    log.input(path)
    pr = json.loads(path.read_text())
    src = pr["source"]
    got = ids(src["pmid"])
    for k in ("pmcid", "doi"):
        if (got.get(k) or "").lower() != (src.get(k) or "").lower():
            raise SystemExit(f"{rel(path)}: {k} {src.get(k)!r}, NCBI's ID converter "
                             f"{got.get(k)!r}")
    doc = full_text(src["pmcid"])
    quotes = []
    for q in pr.get("quotes", []):
        found = doc["text"] is not None and _quote_text(q["text"]) in doc["text"]
        quotes.append({"where": q["where"], "found": found,
                       "first_words": " ".join(q["text"].split()[:8])})
    record = {"gene": pr["gene"], "checked_utc": utc_now(),
              "paper_rule": rel(path), "paper_rule_sha256": sha256_of(path),
              "id_converter": got,
              "document": {k: doc.get(k) for k in ("url", "fetched_utc", "bytes", "sha256")},
              "documents_tried": doc["tried"],
              "folding": ("block elements read as spaces, inline markup dropped, entities "
                          "resolved, Unicode NFKC, whitespace collapsed; the same on both sides"),
              "quotes": quotes}
    out = gene_dir / "paper_quote_check.json"
    out.write_text(json.dumps(record, indent=1, ensure_ascii=False) + "\n")
    log.output(out)
    print(f"{pr['gene']}: {src['pmid']} / {src['pmcid']}: "
          f"{sum(q['found'] for q in quotes)}/{len(quotes)} quotes found in "
          f"{doc['url'] or 'no full text served'}")
    return {"gene": pr["gene"], "found": sum(q["found"] for q in quotes), "of": len(quotes)}


def main() -> None:
    with RunLog(__file__) as log:
        done = [check(d, log) for d in sorted(AUDIT.glob("*_raw"))
                if (d / "paper_rule.json").is_file()]
        log.note("papers", done)
        if any(d["found"] != d["of"] for d in done):
            raise SystemExit("a quote was not found in its paper; fix paper_rule.json")


if __name__ == "__main__":
    main()
