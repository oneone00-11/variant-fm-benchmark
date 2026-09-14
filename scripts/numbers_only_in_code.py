"""Measured numbers in the manuscript that no output carries, and copies of them in code.

A value check can be satisfied by coincidence. In the frozen-matrix-v2 revision the TP53
orientation AUROC was printed as 0.998, a value that existed only in a config comment,
and the check matched it to BRCA2's 0.9975. The number check now binds each measured
number to a cell of an output its paragraph or table declares, or to a derivation it
computes from those outputs. This script audits the result from the other side:

  1. every measured number that is not bound that way. Such a number has no persisted
     source and must be written to an output or removed from the text; the run fails
     while any remain.
  2. every distinctive measured number (three or more significant digits, or 100 or
     more) also printed in a tracked file outside the pipeline outputs -- source,
     configuration, tests, README, hand-written docs -- beside how the check binds it.
     A copy next to a bound number restates an output; it does not stand in for one.

The page also records what the binding round found and how each case was resolved.

    python scripts/numbers_only_in_code.py <manuscript.docx> [--out docs/numbers-only-in-code.md]
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "phase1"))
import src.check_manuscript_numbers as cmn  # noqa: E402

OUTPUT_PREFIXES = ("phase1/reports/", "phase1/data/", "data/", "results/")
GENERATED = ("docs/v1-v2-delta.md", "docs/manuscript-replacements", "docs/manuscript-number-diff",
             "docs/weak-token-spotcheck", "docs/manuscript-binding.md", "docs/numbers-only-in-code.md",
             "docs/supplement-verify.md", "tests/fixtures/")
SOURCE_SUFFIXES = (".py", ".json", ".md", ".yml", ".yaml", ".txt", ".cfg", ".toml", ".sh")
MEASURED = ("bound", "pool", "whitelisted", "no_source", "scoped_mismatch")

# (as printed now, what it is, where it lived, how it was resolved)
RESOLVED = [
    ("0.997", "TP53 control AUROC after orientation; printed 0.998 until the rounding round",
     "a config comment; the value check matched 0.998 to BRCA2's 0.9975",
     "phase4_external_tp53.py writes it to phase4_tp53_orientation.csv (0.996925) on every run and the text prints "
     "0.997; a cell now binds only if it rounds to the print, so 0.9975 no longer does, and the "
     "paragraph's claim pins the number to the orientation cell"),
    ("21,409", "SpliceAI values that re-round to the v1 print (Methods 2.5)",
     "an assertion in tests/test_frozen_matrix_v2.py, and the whitelist",
     "phase1_build_frozen_matrix_v2.py writes it to frozen_matrix_v2_column_report.tsv (n_rerounded_to_v1_print)"),
    ("0.2149994", "full-precision score of the one re-rounding exception", "the same test, and the whitelist",
     "written to frozen_matrix_v2_column_report.tsv (exception_v2_value)"),
    ("0.9997", "rank correlation of the re-scored Nucleotide Transformer column with v1",
     "frozen_matrix_v2_column_report.tsv, outside the outputs the check read; whitelisted",
     "the paragraph declares the column report; bound to its nt row"),
    ("0.672", "agreement of AlphaGenome client v0.6.1 with the companion's v0.7.0",
     "docs/column-provenance.md, the page a one-off audit script writes",
     "scripts/column_provenance_audit.py also writes data/rescore/column_provenance_shared.tsv (alphagenome, spearman)"),
    ("72", "fusion-versus-single-tool comparisons whose interval excludes zero",
     "no cell; the value check matched it to unrelated 0.72 values read as a percentage",
     "a percentage reading now needs a printed %; the count is computed from phase3_pertool_brier_ci.csv"),
    ("30", "disagreements between the two label sets", "the claims file, as typed arithmetic (18 + 12)",
     "computed from phase7_label_confusion.csv (its off-diagonal cells)"),
    ("1.44", "smallest LR+ ratio, ClinVar-recorded over assay-only",
     "the claims file, as typed arithmetic (20.333 / 14.141)", "computed from phase7_selection_test.csv (min_ratio)"),
    ("4.20", "largest such ratio", "the claims file, as typed arithmetic (12.518 / 2.979)",
     "computed from phase7_selection_test.csv (max_ratio)"),
    ("46,392", "SNVs in the companion atlas's matrix", "the whitelist",
     "a property of the companion study's data, not this pipeline's: set aside as cited, with its reference"),
]


def _distinctive(tok: str) -> bool:
    digits = re.sub(r"\D", "", tok).lstrip("0")
    return len(digits) >= 3 or abs(cmn._norm(tok)) >= 100


def _status(t: dict, key: str) -> str:
    if key == "bound":
        c = t["cells"][0]
        if c.startswith("derived:"):
            return "declared derivation" if "(declared, not computed)" in c else "computed derivation"
        return "single cell" if t["n_cells"] == 1 else f"one of {t['n_cells']} cells"
    if key == "non_measurement":
        return f"not a measurement ({t['kind']})"
    return {"pool": "pool match only", "whitelisted": "whitelisted", "no_source": "NO SOURCE",
            "scoped_mismatch": "SCOPED MISMATCH"}[key]


def _persisted(status: str) -> bool:
    return status in ("single cell", "computed derivation") or status.startswith("one of ")


def _md(s: str) -> str:
    return s.replace("|", "\\|")


def _sources() -> list[str]:
    out = subprocess.run(["git", "ls-files", "-z"], cwd=REPO, capture_output=True, text=True, check=True).stdout
    return [f for f in out.split("\0") if f.endswith(SOURCE_SUFFIXES)
            and not f.startswith(OUTPUT_PREFIXES) and not f.startswith(GENERATED)]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("docx")
    ap.add_argument("--out", default="docs/numbers-only-in-code.md")
    a = ap.parse_args()
    r = cmn.classify(Path(a.docx))
    measured = [(t, _status(t, k)) for k in MEASURED for t in r[k]]
    unpersisted = [(t, st) for t, st in measured if not _persisted(st)]
    now: dict[str, list[str]] = {}
    for k in (*MEASURED, "non_measurement"):
        for t in r[k]:
            now.setdefault(t["token"], []).append(f"{t['uid']}: {_status(t, k)}")

    lines_of = {f: (REPO / f).read_text(errors="ignore").splitlines() for f in _sources() if (REPO / f).is_file()}
    distinctive = [(t, st) for t, st in measured if _distinctive(t["token"])]
    rows = []
    for t, st in distinctive:
        tok = t["token"].replace("−", "-").lstrip("+")
        forms = sorted({tok, tok.replace(",", ""), tok.lstrip("-"), tok.replace(",", "").lstrip("-")},
                       key=len, reverse=True)
        pat = re.compile(r"(?<![\d.])(?:%s)(?![\d])" % "|".join(map(re.escape, forms)))
        copies = [f"{f}:{i}" for f, lines in lines_of.items() for i, line in enumerate(lines, 1) if pat.search(line)]
        if copies:
            rows.append((t, st, copies))

    out = ["# Measured numbers that no output carries, and copies of them in code", "",
           f"Manuscript `{Path(a.docx).name}`, {r['coverage']['measured']} measured numbers. Written by "
           "`python scripts/numbers_only_in_code.py <docx>`. The number check's own record of every token is "
           "`docs/manuscript-binding.md`.", "",
           f"## 1. Measured numbers with no persisted source: {len(unpersisted)}", "",
           "A measured number is persisted when the number check binds it to a cell of an output its paragraph or "
           "table declares, or to a derivation the check computes from those outputs. A pool match, a whitelisted "
           "value, or a derivation declared in the claims file without being computed has no output that carries it.",
           ""]
    if unpersisted:
        out += ["| block | number | status | context |", "|---|---|---|---|"]
        out += [f"| {t['uid']} | {t['token']} | {st} | …{_md(t['context'])}… |" for t, st in unpersisted]
    else:
        out.append("None.")
    out += ["", "## 2. Found and resolved in the binding round", "",
            "| printed | what it is | where it lived | resolution | binding now |", "|---|---|---|---|---|"]
    out += [f"| {tok} | {what} | {lived} | {res} | {'; '.join(now.get(tok, ['not printed']))} |"
            for tok, what, lived, res in RESOLVED]
    out += ["", f"## 3. Distinctive measured numbers with a copy outside the outputs: {len(rows)} of {len(distinctive)}",
            "", f"The {len(distinctive)} measured numbers with three or more significant digits, or of 100 or more, were "
            f"searched for in {len(lines_of)} tracked source, configuration, test, README and hand-written documentation "
            "files; pipeline outputs and generated pages are excluded. A shorter literal such as 7 or 0.5 occurs in "
            "hundreds of unrelated lines, so a copy of it says nothing about where a number came from; section 1 "
            "covers those numbers too.", "",
            "| block | number | binding | copies (first three) | context |", "|---|---|---|---|---|"]
    for t, st, copies in rows:
        more = f" (+{len(copies) - 3})" if len(copies) > 3 else ""
        out.append(f"| {t['uid']} | {t['token']} | {st} | {', '.join(copies[:3])}{more} | …{_md(t['context'][:80])}… |")
    dest = Path(a.out) if Path(a.out).is_absolute() else REPO / a.out
    dest.write_text("\n".join(out) + "\n")
    print(f"{len(unpersisted)} measured numbers with no persisted source; {len(rows)} of {len(distinctive)} distinctive "
          f"measured numbers have a copy outside the outputs -> {a.out}")
    return 1 if unpersisted else 0


if __name__ == "__main__":
    sys.exit(main())
