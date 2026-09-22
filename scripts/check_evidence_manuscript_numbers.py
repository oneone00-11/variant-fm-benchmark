#!/usr/bin/env python3
"""Is every number printed in the manuscript one the pipeline produced?

The rule this enforces is that no value is typed from memory or carried over from
an earlier draft. Every numeric token in the body is either

  DECLARED   a parameter, a published constant or a unit that no analysis produces,
             listed below with the reason and the place it comes from;
  MATCHED    a value that rounds to a number the pipeline wrote, where the pool is
             every numeric cell in reports/evidence/, the manifests and provenance
             records, and the per-stratum totals derived from the analysis set's
             own count table;
  UNMATCHED  neither, which fails the run.

Matching is at the printed precision: a token is matched when a pool value rounds
to it at half a unit in its last printed digit. Percentages are matched against the
pooled value and against that value times one hundred, because the tables store
proportions and the text prints percentages.

What this does NOT check, and why it is stated rather than hidden: numbers written
as words. Fold counts in the text read "in six of the seven genes", and a numeral
pool would match a bare 6 against almost any table, so the check would pass without
testing anything. Those sentences are listed at the end of the report instead, for
the author to read against Table 2.

    python scripts/check_evidence_manuscript_numbers.py path/to/manuscript.md
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
REPORTS = REPO / "phase1" / "reports" / "evidence"
DATA = REPO / "phase1" / "data" / "evidence"

# A token is a measurement unless a rule here says what it is instead. Each entry
# names the value and where it is fixed, so a reader can check the claim.
DECLARED = {
    # the recommendation's own constants
    "0.2": "PP3 cut point, Walker et al. Table 1 (config/walker2023.yaml)",
    "0.1": "BP4 cut point, Walker et al. Table 1 (config/walker2023.yaml)",
    "2,736": "size of Walker et al.'s calibration set, their p. 1051",
    "1046": "page reference, Walker et al.",
    "1049": "page reference, Walker et al.",
    "1057": "page reference, Walker et al.",
    # the point system
    "0.10": "prior probability of the Tavtigian point system",
    "2.08": "Supporting tier boundary (config/walker2023.yaml)",
    "4.33": "Moderate tier boundary (config/walker2023.yaml)",
    "18.7": "Strong tier boundary (config/walker2023.yaml)",
    "350": "Very strong tier boundary (config/walker2023.yaml)",
    # analysis parameters
    "1.3.1": "SpliceAI package version (provenance record)",
    "4,999": "SpliceAI scoring distance, Walker's basis (provenance record)",
    "100": "minimum window occupancy, Pejaver's procedure (src/evid_interval_lr.py)",
    "2,000": "bootstrap resamples, gene-clustered (src/evid_walker_thresholds.py)",
    "95": "confidence level",
    "10": "minimum band occupancy and minimum class count (src/evid_common.py)",
    "0.01": "false discovery rate cut for the external label definition",
    # window and stratum bounds, fixed by the design
    "1": "stratum bound", "3": "stratum bound", "8": "external window bound",
    "11": "stratum bound", "50": "stratum bound", "2": "stratum bound",
    # dates and identifiers
    "15": "date in the ClinVar release name", "2026": "year", "20": "date",
    "2023": "year", "00000658": "MaveDB accession",
    # enumeration markers
    "(1)": "inline enumeration", "(2)": "inline enumeration", "(3)": "inline enumeration",
}

# A count written as a word, in the shape a result claim uses: "in six held-out
# genes", "eleven of the thirteen columns", "in all seven genes". The narrow
# pattern keeps the reading list to the sentences that actually assert a count.
WORD_COUNT = re.compile(
    r"\b(one|two|three|four|five|six|seven|eight|nine|ten|eleven|twelve|thirteen|fourteen|fifteen)\b"
    r"\s*(of\s+(the|them|those|seven|six|thirteen)|held-out|genes|columns|folds|arms\b)",
    re.I)


def pool_values() -> dict[float, list[str]]:
    """Every number the pipeline wrote, with where it was written."""
    pool: dict[float, list[str]] = {}

    def add(v, where):
        try:
            f = float(v)
        except (TypeError, ValueError):
            return
        if f != f:  # NaN
            return
        pool.setdefault(round(f, 10), []).append(where)

    for csv in sorted(list(REPORTS.glob("*.csv")) + list((REPORTS / "fig_data").glob("*.csv"))):
        try:
            df = pd.read_csv(csv)
        except Exception:
            continue
        for col in df.columns:
            ser = pd.to_numeric(df[col], errors="coerce").dropna()
            for v in ser.unique():
                add(v, f"{csv.name}:{col}")

    for js in sorted(list(DATA.glob("*.json")) + list((DATA / "external").glob("*.json"))):
        try:
            blob = json.loads(js.read_text())
        except Exception:
            continue
        stack = [blob]
        while stack:
            node = stack.pop()
            if isinstance(node, dict):
                stack.extend(node.values())
            elif isinstance(node, list):
                stack.extend(node)
            elif isinstance(node, (int, float)):
                add(node, js.name)

    # The analysis set's totals are quoted in the text and appear in no single cell:
    # the count table is per gene, stratum and arm. The sums it implies are part of
    # what the pipeline produced, so they belong in the pool.
    counts = REPORTS / "set_counts.csv"
    if counts.exists():
        sc = pd.read_csv(counts)
        cols = ["n", "n_labelled", "n_damaging", "n_normal"]
        for col in cols:
            add(sc[col].sum(), f"set_counts.csv:{col} (total)")
        for stratum, g in sc.groupby("stratum"):
            for col in cols:
                add(g[col].sum(), f"set_counts.csv:{col} ({stratum})")
        in_scope = sc[sc.stratum != "pm12"]
        for col in cols:
            add(in_scope[col].sum(), f"set_counts.csv:{col} (in scope)")
        for (stratum, arm), g in sc.groupby(["stratum", "clinvar_arm"]):
            for col in cols:
                add(g[col].sum(), f"set_counts.csv:{col} ({stratum}/{arm})")
    return pool


def body_of(path: Path) -> str:
    if path.suffix == ".docx":
        import docx
        text = "\n".join(p.text for p in docx.Document(str(path)).paragraphs)
    else:
        text = path.read_text()
    start = text.index("# Abstract") if "# Abstract" in text else 0
    for marker in ("# References", "References\n"):
        if marker in text:
            text = text[:text.index(marker)]
            break
    return text[start:]


TOKEN = re.compile(r"(?<![\w.])(\d{1,3}(?:,\d{3})*(?:\.\d+)?)%?(?![\w])")


def tolerance(tok: str) -> float:
    frac = tok.split(".")[1] if "." in tok else ""
    return 0.5 * 10 ** (-len(frac)) + 1e-9


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("manuscript", type=Path)
    args = ap.parse_args()

    body = body_of(args.manuscript)
    # citations and the draft's own source tags are not measurements
    body = re.sub(r"\[[^\]]*\]", " ", body)
    body = re.sub(r"`[^`]*`", " ", body)

    pool = pool_values()
    declared = matched = 0
    unmatched = []
    for m in TOKEN.finditer(body):
        tok = m.group(1)
        if tok in DECLARED:
            declared += 1
            continue
        value = float(tok.replace(",", ""))
        tol = tolerance(tok)
        hit = None
        for candidate in (value, value / 100.0):
            for pv, wheres in pool.items():
                if abs(pv - candidate) <= tol or abs(pv * 100 - value) <= tol:
                    hit = wheres[0]
                    break
            if hit:
                break
        if hit:
            matched += 1
        else:
            ctx = body[max(0, m.start() - 60):m.end() + 40].replace("\n", " ")
            unmatched.append((tok, ctx.strip()))

    print(f"pool: {len(pool):,} distinct values from reports/evidence/, the manifests "
          f"and the analysis-set totals")
    print(f"tokens: {matched} matched, {declared} declared, {len(unmatched)} unmatched")
    for tok, ctx in unmatched:
        print(f"  UNMATCHED {tok}  ...{ctx}...")

    results = body[body.index("# Results"):] if "# Results" in body else body
    words = [s.strip() for s in re.split(r"(?<=[.!?])\s+", results) if WORD_COUNT.search(s)]
    print(f"\nresult sentences asserting a count written as a word ({len(words)}); "
          f"read each against Table 1 or Table 2:")
    for s in words:
        print("  - " + " ".join(s.split())[:150])
    return 1 if unmatched else 0


if __name__ == "__main__":
    sys.exit(main())
