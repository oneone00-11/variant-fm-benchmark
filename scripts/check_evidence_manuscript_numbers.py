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
to it at half a unit in its last printed digit. A number printed as a percentage
("%" or "per cent") is also matched against a stored proportion times one hundred,
because the tables store proportions and the text prints percentages; a number not
printed as a percentage is matched only as it stands.

A token is declared by what it is in its sentence, not by its digits alone: "95" is
the confidence level in "95% CI" and a measurement in "at least 95% of the cases",
and "1" is a stratum bound in "1-2 nt" and an enumeration in "(1)". The context rules
come first, then the plain list.

What this does NOT check, and why it is stated rather than hidden: numbers written
as words. Fold counts in the text read "in six of the seven genes", and a numeral
pool would match a bare 6 against almost any table, so the check would pass without
testing anything. Those sentences are listed at the end of the report instead, for
the author to read against Table 2.

Three modes, from weakest to strongest:

    python scripts/check_evidence_manuscript_numbers.py path/to/manuscript.docx
        pool mode: every number against every value the pipeline wrote
    python scripts/check_evidence_manuscript_numbers.py path/to/manuscript.md
        bound mode: every number against the file its own sentence names
    python scripts/check_evidence_manuscript_numbers.py path/to/manuscript.md \\
            --cells docs/manuscript_bindings_v5.tsv --document main
        cells mode: every number against the one cell a binding table names

`--selftest` plants numbers no named file holds and passes only if the bound mode
reports every one of them; `--declared` lists each declared token with its reason.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from dataclasses import dataclass, field
from decimal import ROUND_HALF_EVEN, ROUND_HALF_UP, Decimal, InvalidOperation
from pathlib import Path

import pandas as pd
import yaml

REPO = Path(__file__).resolve().parents[1]
REPORTS = REPO / "phase1" / "reports" / "evidence"
DATA = REPO / "phase1" / "data" / "evidence"

# A token is a measurement unless a rule here says what it is instead. Each entry
# names the value and where it is fixed, so a reader can check the claim. The
# context rules below (CONTEXT_RULES) are tried first.
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
    "10": "minimum band occupancy and minimum class count (src/evid_common.py)",
    "0.01": "false discovery rate cut for the external label definition",
    # window and stratum bounds, fixed by the design
    "1": "stratum bound", "3": "stratum bound", "8": "external window bound",
    "11": "stratum bound", "50": "stratum bound", "2": "stratum bound",
    # dates and identifiers
    "15": "date in the ClinVar release name", "2026": "year", "20": "date",
    "2023": "year", "00000658": "MaveDB accession",
}

MONTHS = ("January|February|March|April|May|June|July|August|September|October|"
          "November|December")
# "95% CI", "95% confidence interval", "one-sided 95% lower bound", "the one-sided 95%
# bound", "Confidence intervals (CI) are 95%": the confidence level, written right
# beside the interval or bound it qualifies.
CONFIDENCE_AFTER = re.compile(r"%\s*(?:CI\b|confidence\b|(?:lower|upper)\s+bound\b|bound\b)", re.I)
CONFIDENCE_BEFORE = re.compile(r"\b(?:intervals?|CI\)?)\s+(?:are|is)\s+$", re.I)


@dataclass(frozen=True)
class Context:
    """A token as its sentence prints it: the cleaned sentence, the token's span in it
    (digits only) and whether the number is printed as a percentage."""
    clean: str
    start: int
    end: int
    percent: bool

    @property
    def before(self) -> str:
        return self.clean[:self.start]

    @property
    def after(self) -> str:
        return self.clean[self.end:]


def _enumeration(tok: str, c: Context) -> str | None:
    if tok.isdigit() and int(tok) < 10 and c.before.endswith("(") and c.after.startswith(")"):
        return "inline enumeration, (1) (2) (3)"
    return None


def _date(tok: str, c: Context) -> str | None:
    m = re.match(rf"\s+({MONTHS})\b(?:\s+(\d{{4}}))?", c.after)
    if tok.isdigit() and 1 <= int(tok) <= 31 and m:
        return f"day of the date {tok} {m.group(1)}" + (f" {m.group(2)}" if m.group(2) else "")
    return None


def _confidence_level(tok: str, c: Context) -> str | None:
    if tok != "95" or not c.percent:
        return None
    if CONFIDENCE_AFTER.match(c.after) or CONFIDENCE_BEFORE.search(c.before):
        return "confidence level"
    return None


def _band_bound(tok: str, c: Context) -> str | None:
    if tok == "10" and re.search(r"\b3\s*[–-]\s*$", c.before):
        return "stratum bound, the upper end of 3-10 nt"
    return None


def _scoring_distance(tok: str, c: Context) -> str | None:
    if tok == "50" and re.search(r"\bdistance(?:\s+of)?\s*$", c.before, re.I):
        return "SpliceAI scoring distance of 50 nt (the atlas column, the fusion input)"
    return None


def _allele_frequency(tok: str, c: Context) -> str | None:
    if tok == "0.1" and c.percent and re.search(r"allele frequency", c.clean, re.I):
        return "gnomAD filtering allele frequency of 0.1% (the Atlas model's labels)"
    return None


def _scale_anchor(tok: str, c: Context) -> str | None:
    if tok == "1" and c.before[-1:] in ("−", "+") and re.search(r"\banchored\b", c.clean):
        return "anchor of the score scale, -1 synonymous and +1 nonsense controls"
    return None


CONTEXT_RULES = [_enumeration, _date, _confidence_level, _band_bound, _scoring_distance,
                 _allele_frequency, _scale_anchor]
# declared only through a context rule: printed elsewhere, the number is a measurement
CONTEXT_ONLY = {"95"}


def declared_reason(tok: str, c: Context) -> str | None:
    """Why a token is not a measurement, or None when it is one."""
    for rule in CONTEXT_RULES:
        reason = rule(tok, c)
        if reason:
            return reason
    if tok in CONTEXT_ONLY:
        return None
    return DECLARED.get(tok)


# A count written as a word, in the shape a result claim uses: "in six held-out
# genes", "eleven of the thirteen columns", "in all seven genes". The narrow
# pattern keeps the reading list to the sentences that actually assert a count.
WORD_COUNT = re.compile(
    r"\b(one|two|three|four|five|six|seven|eight|nine|ten|eleven|twelve|thirteen|fourteen|fifteen)\b"
    r"\s*(of\s+(the|them|those|seven|six|thirteen)|held-out|genes|columns|folds|arms\b)",
    re.I)

TOKEN = re.compile(r"(?<![\w.])(\d{1,3}(?:,\d{3})*(?:\.\d+)?)%?(?![\w])(?!\.\d)")
PER_CENT = re.compile(r"\s*per\s?cent\b", re.I)
# identifiers are not measurements: DOIs, URLs, accessions, commits
IDENTIFIER = re.compile(r"(?:https?://\S+|doi:\S+|10\.\d{4,}/\S+|urn:\S+|\bzenodo\.\d+|"
                        r"\bNM_\d+(?:\.\d+)?|\b[0-9a-f]{7,40}\b)", re.I)


def is_percent(printed: str, after: str) -> bool:
    """Printed as a percentage: with a percent sign, or followed by "per cent"."""
    return printed.endswith("%") or bool(PER_CENT.match(after))


def tolerance(tok: str) -> float:
    frac = tok.split(".")[1] if "." in tok else ""
    return 0.5 * 10 ** (-len(frac)) + 1e-9


def _match_value(tok: str, pairs, percent: bool):
    """The first (value, where) whose value rounds to the token; a proportion times one
    hundred only when the token is printed as a percentage."""
    value = float(tok.replace(",", ""))
    tol = tolerance(tok)
    for pv, where in pairs:
        if abs(pv - value) <= tol or (percent and abs(pv * 100 - value) <= tol):
            return where
    return None


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

    for csv_path in sorted(list(REPORTS.glob("*.csv")) + list((REPORTS / "fig_data").glob("*.csv"))):
        try:
            df = pd.read_csv(csv_path)
        except Exception:
            continue
        for col in df.columns:
            ser = pd.to_numeric(df[col], errors="coerce").dropna()
            for v in ser.unique():
                add(v, f"{csv_path.name}:{col}")

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


def pool_check(path: Path) -> int:
    body = body_of(path)
    # citations and the draft's own source tags are not measurements
    body = re.sub(r"\[[^\]]*\]", " ", body)
    body = re.sub(r"`[^`]*`", " ", body)

    pool = pool_values()
    pairs = [(pv, wheres[0]) for pv, wheres in pool.items()]
    declared = matched = 0
    unmatched = []
    # sentence by sentence, so that a context rule reads the token's own sentence
    for _, para in _pieces(body, PARAGRAPH_BREAK):
        for _, sent in _pieces(para, SENTENCE_BREAK):
            for m in TOKEN.finditer(sent):
                tok = m.group(1)
                ctx = Context(sent, m.start(1), m.end(1), is_percent(m.group(0), sent[m.end():]))
                if declared_reason(tok, ctx):
                    declared += 1
                    continue
                if _match_value(tok, pairs, ctx.percent):
                    matched += 1
                else:
                    near = sent[max(0, m.start() - 60):m.end() + 40].replace("\n", " ")
                    unmatched.append((tok, near.strip()))

    print(f"pool: {len(pool):,} distinct values from reports/evidence/, the manifests "
          f"and the analysis-set totals")
    print(f"tokens: {matched} matched, {declared} declared, {len(unmatched)} unmatched")
    for tok, near in unmatched:
        print(f"  UNMATCHED {tok}  ...{near}...")

    results = body[body.index("# Results"):] if "# Results" in body else body
    words = [s.strip() for s in re.split(r"(?<=[.!?])\s+", results) if WORD_COUNT.search(s)]
    print(f"\nresult sentences asserting a count written as a word ({len(words)}); "
          f"read each against Table 1 or Table 2:")
    for s in words:
        print("  - " + " ".join(s.split())[:150])
    return 1 if unmatched else 0


# ---------------------------------------------------------------------------
# reading the Markdown source: one record per numeric token
# ---------------------------------------------------------------------------
# The Markdown source carries a tag naming the output each number was read from.
# The bound and cells modes read the same tokens, in the same order: the text from
# "# Abstract" (or the top, when there is no such heading, as in the supplementary
# notes) to "# References" (or the end), paragraph by paragraph and sentence by
# sentence, with the source tags, identifiers, citations and table, figure and note
# labels taken out before numbers are read.
SEARCH_DIRS = [REPORTS, REPORTS / "fig_data", REPORTS / "tables", REPORTS / "supplement",
               DATA, DATA / "external", REPORTS.parents[1] / "config"]
REF = re.compile(r"(?:Supplementary\s+)?(?:Tables?|Figures?|Notes?)\s+S?\d+[a-d]?"
                 r"(?:\s*(?:,|and|to|–)\s*S?\d+[a-d]?)*")
TAGGED = re.compile(r"\[([^\]]*?\.(?:csv|json|yaml)(?:[^\]]*))\]")
PARAGRAPH_BREAK = re.compile(r"\n\s*\n")
SENTENCE_BREAK = re.compile(r"(?<=[.!?])\s+(?=[A-Z*(])")


@dataclass
class Token:
    line: int                 # the line of the source file the sentence starts on
    printed: str              # as printed, with its percent sign
    tok: str                  # the digits
    context: Context          # the cleaned sentence and the token's place in it
    sentence: str             # the sentence as written, whitespace collapsed
    tags: list[str]           # the files the sentence names, or the paragraph's last named
    occurrence: int = 1       # its rank among the tokens printed alike on its line
    reason: str | None = None  # why it is declared, or None for a measurement

    @property
    def excerpt_space(self) -> str:
        return " ".join(self.context.clean.split())


def body_span(text: str) -> tuple[int, int]:
    start = text.index("# Abstract") if "# Abstract" in text else 0
    end = text.index("# References") if "# References" in text else len(text)
    return start, end


def _pieces(text: str, sep: re.Pattern):
    """(offset, piece) for every piece of `text` between matches of `sep`, the pieces
    re.split would return."""
    pos = 0
    for m in sep.finditer(text):
        yield pos, text[pos:m.start()]
        pos = m.end()
    yield pos, text[pos:]


def _tags(sentence: str) -> list[str]:
    tags = []
    for m in TAGGED.finditer(sentence):
        tags += [t.strip() for t in re.split(r"[;,]", m.group(1))
                 if re.search(r"\.(csv|json|yaml)$", t.strip())]
    return tags


def clean_sentence(sentence: str) -> str:
    clean = TAGGED.sub(" ", sentence)
    clean = IDENTIFIER.sub(" ", clean)              # DOIs, URLs, accessions
    clean = re.sub(r"\[[^\]]*\]", " ", clean)      # citations
    return REF.sub(" ", clean)                      # table, figure, note labels


def read_tokens(text: str) -> list[Token]:
    start, end = body_span(text)
    body = text[start:end]
    out: list[Token] = []
    seen: dict[tuple[int, str], int] = {}
    for p_off, para in _pieces(body, PARAGRAPH_BREAK):
        if para.lstrip().startswith("#"):
            continue
        last_tags: list[str] = []
        for s_off, sent in _pieces(para, SENTENCE_BREAK):
            line = text.count("\n", 0, start + p_off + s_off) + 1
            tags = _tags(sent)
            if tags:
                last_tags = tags
            clean = clean_sentence(sent)
            for m in TOKEN.finditer(clean):
                tok = m.group(1)
                ctx = Context(clean, m.start(1), m.end(1), is_percent(m.group(0), clean[m.end():]))
                key = (line, m.group(0))
                seen[key] = seen.get(key, 0) + 1
                out.append(Token(line=line, printed=m.group(0), tok=tok, context=ctx,
                                 sentence=" ".join(sent.split()), tags=list(last_tags),
                                 occurrence=seen[key], reason=declared_reason(tok, ctx)))
    return out


# ---------------------------------------------------------------------------
# bound mode: each number against the file its own sentence names
# ---------------------------------------------------------------------------
# The pool mode above matches a token against every value the pipeline wrote, which
# is more than half a million values: at two or three significant figures almost
# any number finds a match, so a pass there shows only that a number is not
# impossible. Bound mode holds every number to the file named in its own sentence,
# or, where a sentence carries no tag, in the nearest earlier sentence of the same
# paragraph. A number whose paragraph names no file is reported as unbound. It still
# shows only that the named file holds the value somewhere: the cells mode below
# holds each number to one cell.


def _file_values(name: str) -> dict[float, str]:
    for d in SEARCH_DIRS:
        f = d / name
        if f.exists():
            break
    else:
        return {}
    out: dict[float, str] = {}
    if f.suffix == ".csv":
        df = pd.read_csv(f)
        for col in df.columns:
            for v in pd.to_numeric(df[col], errors="coerce").dropna().unique():
                out.setdefault(round(float(v), 10), f"{name}:{col}")
        if name == "set_counts.csv":
            cols = ["n", "n_labelled", "n_damaging", "n_normal"]
            groups = [("total", df)] + [(str(k), g) for k, g in df.groupby("stratum")] + \
                     [("in scope", df[df.stratum != "pm12"])] + \
                     [(f"{a}/{b}", g) for (a, b), g in df.groupby(["stratum", "clinvar_arm"])]
            for label, g in groups:
                for c in cols:
                    out.setdefault(round(float(g[c].sum()), 10), f"{name}:{c} ({label})")
    else:
        blob = (yaml.safe_load(f.read_text()) if f.suffix == ".yaml"
                else json.loads(f.read_text()))
        stack = [blob]
        while stack:
            node = stack.pop()
            if isinstance(node, dict):
                stack.extend(node.values())
            elif isinstance(node, list):
                stack.extend(node)
            elif isinstance(node, (int, float)):
                out.setdefault(round(float(node), 10), name)
    return out


def _matches(tok: str, values: dict[float, str], percent: bool = False) -> str | None:
    return _match_value(tok, values.items(), percent)


@dataclass
class BoundResult:
    bound: list[Token] = field(default_factory=list)
    declared: list[Token] = field(default_factory=list)
    fails: list[Token] = field(default_factory=list)
    unbound: list[Token] = field(default_factory=list)


def bound_results(text: str) -> BoundResult:
    res = BoundResult()
    cache: dict[str, dict[float, str]] = {}
    for t in read_tokens(text):
        if t.reason:
            res.declared.append(t)
        elif not t.tags:
            res.unbound.append(t)
        else:
            vals: dict[float, str] = {}
            for name in t.tags:
                if name not in cache:
                    cache[name] = _file_values(name)
                vals.update(cache[name])
            (res.bound if _matches(t.tok, vals, t.context.percent) else res.fails).append(t)
    return res


def bound_check(md_path: Path, show_declared: bool = False) -> int:
    res = bound_results(md_path.read_text())
    print(f"bound mode on {md_path.name}: {len(res.bound)} numbers matched in the file their "
          f"sentence names, {len(res.declared)} declared constants, {len(res.fails)} not found "
          f"in their named file, {len(res.unbound)} in a paragraph that names no file")
    for t in res.fails:
        print(f"  NOT IN NAMED FILE {t.tok} [{', '.join(t.tags)}]  ...{t.sentence[:140]}...")
    for t in res.unbound:
        print(f"  UNBOUND {t.tok}  ...{t.sentence[:140]}...")
    if show_declared:
        for t in res.declared:
            print(f"  DECLARED {t.printed} line {t.line}: {t.reason}")
    return 1 if (res.fails or res.unbound) else 0


# The self-test plants numbers, one sentence each, right after the first sentence of
# the body that names a .csv file, and passes only if the bound mode reports every
# planted number. Before judging, it confirms that each planted number is neither a
# declared constant nor a value the named files hold, so that a number can only be
# caught or missed. Whether the text's own numbers pass has no bearing on the verdict.
PLANTED = ["0.3712", "12.34", "987", "0.9137", "45.6%", "7,777"]
PLANT_WORD = "Selftest-planted"


def plant(text: str, planted: list[str]) -> tuple[str, list[str], int] | None:
    """The text with the numbers planted, the files the planting sentence names and
    its line; None when no sentence of the body names a .csv file."""
    start, end = body_span(text)
    body = text[start:end]
    for p_off, para in _pieces(body, PARAGRAPH_BREAK):
        if para.lstrip().startswith("#"):
            continue
        for s_off, sent in _pieces(para, SENTENCE_BREAK):
            named = _tags(sent)
            if any(x.endswith(".csv") for x in named):
                at = start + p_off + s_off + len(sent)
                line = text.count("\n", 0, start + p_off + s_off) + 1
                sown = "".join(f" {PLANT_WORD} {p}." for p in planted)
                return text[:at] + sown + text[at:], named, line
    return None


def selftest(md_path: Path, planted: list[str] = PLANTED) -> int:
    sown = plant(md_path.read_text(), planted)
    if sown is None:
        print("SELFTEST cannot run: no sentence of the body names a .csv file to plant beside")
        return 2
    text, named, line = sown
    mine = [t for t in read_tokens(text) if PLANT_WORD in t.sentence and t.printed in planted]
    declared = [t.printed for t in mine if t.reason]
    vals: dict[float, str] = {}
    for name in named:
        vals.update(_file_values(name))
    held = [t.printed for t in mine if not t.reason and _matches(t.tok, vals, t.context.percent)]
    if sorted(t.printed for t in mine) != sorted(planted) or declared or held:
        print(f"SELFTEST cannot run: of the planted {', '.join(planted)}, "
              f"read {', '.join(t.printed for t in mine) or 'none'}; declared constants "
              f"{', '.join(declared) or 'none'}; held by {', '.join(named)}: "
              f"{', '.join(held) or 'none'}. Plant other numbers.")
        return 2
    print(f"bound mode with {len(planted)} planted numbers beside the sentence on line {line}, "
          f"which names {', '.join(named)}; checked first: no planted number is a declared "
          f"constant or a value those files hold")
    res = bound_results(text)
    reported = {t.printed for t in res.fails + res.unbound if PLANT_WORD in t.sentence}
    caught = [p for p in planted if p in reported]
    for p in planted:
        print(f"  {'caught' if p in caught else 'MISSED'} {p}")
    verdict = "PASSED" if len(caught) == len(planted) else "FAILED"
    print(f"SELFTEST caught {len(caught)} / {len(planted)}: {verdict}")
    return 0 if verdict == "PASSED" else 1


# ---------------------------------------------------------------------------
# cells mode: each number against the one cell a binding table names
# ---------------------------------------------------------------------------
# A binding table (docs/manuscript_bindings_<version>.tsv) has one row per number of
# the text, in these columns:
#
#   document    which text the row belongs to: main or supplement
#   line        the line of the source file the number's sentence starts on
#   occurrence  1 for the first number printed this way on that line, 2 for the next
#   printed     the number as printed, with its percent sign
#   excerpt     a piece of its sentence as this script reads it (source tags and
#               citations out, spaces collapsed), so that a moved line is caught
#   kind        constant (not compared; the reason says why) or measured
#   reason      for a constant, what it is and where it is fixed; for a measured
#               value, an optional note
#   file        the output that holds the value, from the repository root
#   rows        which rows: conditions joined by "; ", each "column=value",
#               "column!=value", "column>=number" (also <=, >, <),
#               "column in a|b|c", "column not in a|b" or "column startswith text";
#               in a JSON or YAML file they test the key a "*" in the path matched
#               ("key in a|b")
#   column      the column of a table, or the path of keys in a JSON or YAML file
#               (dots between keys, [n] for an item of a list, ["key"] for a key with
#               a dot or a colon, * for every key at that level)
#   derive      value (exactly one cell), sum, min, max, count (rows), count_nonnull,
#               all (every cell must round to the number) or text (the cell is text,
#               which must contain the number)
#   scale       1, or 100 for a stored proportion printed as a percentage
#
# A measured number is compared with its derived value at the printed precision: the
# value, times the scale, rounded half up or half to even to as many decimals as the
# number prints must equal it. Three lists come out, and any of them failing the run:
# numbers whose bound value differs (MISMATCH), rows that cannot be resolved because a
# file, row or column is not there, the excerpt is not the sentence's, or no number
# stands at the row's place (CANNOT BIND), and numbers of the text with no row (NOT IN
# TABLE).
BINDING_COLUMNS = ["document", "line", "occurrence", "printed", "excerpt", "kind", "reason",
                   "file", "rows", "column", "derive", "scale"]
DERIVES = {"value", "sum", "min", "max", "count", "count_nonnull", "all", "text"}
CONDITION = re.compile(r"^(?P<col>.+?)\s*(?P<op>!=|>=|<=|=|>|<)\s*(?P<val>.*)$")
WORD_CONDITION = re.compile(r"^(?P<col>.+?)\s+(?P<op>not in|in|startswith)\s+(?P<val>.+)$")


class BindError(Exception):
    """A binding row that cannot be resolved; the message says which part is missing."""


@dataclass
class Binding:
    document: str
    line: int
    occurrence: int
    printed: str
    excerpt: str
    kind: str
    reason: str
    file: str
    rows: str
    column: str
    derive: str
    scale: int

    @property
    def key(self) -> tuple[int, str, int]:
        return self.line, self.printed, self.occurrence

    def where(self) -> str:
        parts = [self.file, self.rows, self.column]
        tail = f" ({self.derive}" + (f", x{self.scale})" if self.scale != 1 else ")")
        return " | ".join(p for p in parts if p) + (tail if self.derive != "value" or self.scale != 1 else "")


def read_bindings(path: Path) -> list[Binding]:
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, quoting=csv.QUOTE_MINIMAL)
    missing = [c for c in BINDING_COLUMNS if c not in df.columns]
    if missing:
        raise SystemExit(f"{path}: the binding table has no column {', '.join(missing)}")
    out = []
    for _, r in df.iterrows():
        out.append(Binding(document=r["document"], line=int(r["line"]), occurrence=int(r["occurrence"]),
                           printed=r["printed"], excerpt=r["excerpt"], kind=r["kind"],
                           reason=r["reason"], file=r["file"], rows=r["rows"], column=r["column"],
                           derive=r["derive"] or "value", scale=int(r["scale"] or 1)))
    keys = [(b.document, *b.key) for b in out]
    dup = sorted({k for k in keys if keys.count(k) > 1})
    if dup:
        raise SystemExit(f"{path}: more than one row for {dup}")
    return out


def write_bindings(path: Path, rows: list[dict]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=BINDING_COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in BINDING_COLUMNS})


def _decimal(x) -> Decimal | None:
    if isinstance(x, bool) or x is None:
        return None
    if isinstance(x, Decimal):
        return x
    try:
        d = Decimal(str(x).strip().replace(",", ""))
    except (InvalidOperation, ValueError):
        return None
    return d if d.is_finite() else None


def _conditions(spec: str) -> list[tuple[str, str, str]]:
    out = []
    for part in [p.strip() for p in spec.split(";") if p.strip()]:
        m = WORD_CONDITION.match(part) or CONDITION.match(part)
        if not m:
            raise BindError(f"cannot read the row condition '{part}'")
        out.append((m["col"].strip(), m["op"], m["val"].strip()))
    return out


def _test(cell: str, op: str, val: str) -> bool:
    if op == "=":
        return cell == val
    if op == "!=":
        return cell != val
    if op == "in":
        return cell in val.split("|")
    if op == "not in":
        return cell not in val.split("|")
    if op == "startswith":
        return cell.startswith(val)
    a, b = _decimal(cell), _decimal(val)
    if a is None or b is None:
        return False
    return {">=": a >= b, "<=": a <= b, ">": a > b, "<": a < b}[op]


PATH_PART = re.compile(r'\["((?:[^"\\]|\\.)*)"\]|\[(\d+)\]|([^.\[\]]+)')


def _path(spec: str) -> list:
    parts, pos = [], 0
    spec = spec.strip()
    while pos < len(spec):
        if spec[pos] == ".":
            pos += 1
            continue
        m = PATH_PART.match(spec, pos)
        if not m:
            raise BindError(f"cannot read the path '{spec}'")
        parts.append(m.group(1) if m.group(1) is not None else
                     int(m.group(2)) if m.group(2) is not None else m.group(3))
        pos = m.end()
    return parts


def _walk(node, parts: list, conds: list, key=None):
    """(key the '*' matched, value) for every node the path reaches."""
    if not parts:
        yield key, node
        return
    head, rest = parts[0], parts[1:]
    if head == "*":
        if not isinstance(node, dict):
            raise BindError("a '*' in the path meets something that is not a mapping")
        for k, v in node.items():
            if all(_test(str(k), op, val) for col, op, val in conds):
                yield from _walk(v, rest, conds, k)
        return
    if isinstance(head, int):
        if not isinstance(node, list) or head >= len(node):
            raise BindError(f"no item [{head}] in the path")
        yield from _walk(node[head], rest, conds, key)
        return
    if not isinstance(node, dict) or head not in node:
        raise BindError(f"no key '{head}' in the path")
    yield from _walk(node[head], rest, conds, key)


_LOADED: dict[tuple, object] = {}


def _load(path: Path):
    """A table as text cells, or a JSON or YAML record with exact decimals; read once
    per file version."""
    st = path.stat()
    key = (str(path), st.st_mtime_ns, st.st_size)
    if key not in _LOADED:
        if path.suffix == ".csv":
            _LOADED[key] = pd.read_csv(path, dtype=str, keep_default_na=False)
        elif path.suffix == ".json":
            _LOADED[key] = json.loads(path.read_text(), parse_float=Decimal, parse_int=Decimal)
        else:
            _LOADED[key] = yaml.safe_load(path.read_text())
    return _LOADED[key]


def resolve(b: Binding, root: Path = REPO) -> list:
    """The cells a measured binding reads: Decimals, or for derive=text one string."""
    path = root / b.file
    if not b.file or not path.is_file():
        raise BindError(f"no file {b.file or '(none named)'}")
    conds = _conditions(b.rows)
    if path.suffix == ".csv":
        df = _load(path)
        for col, _, _ in conds:
            if col not in df.columns:
                raise BindError(f"no column '{col}' in {b.file}")
        keep = pd.Series(True, index=df.index)
        for col, op, val in conds:
            keep &= df[col].map(lambda cell, op=op, val=val: _test(cell, op, val))
        sel = df[keep]
        if b.derive == "count":
            return [Decimal(len(sel))]
        if b.column not in df.columns:
            raise BindError(f"no column '{b.column}' in {b.file}")
        cells = list(sel[b.column])
    elif path.suffix in (".json", ".yaml", ".yml"):
        blob = _load(path)
        if conds and any(col != "key" for col, _, _ in conds):
            raise BindError("in a JSON or YAML file the rows can only test the key")
        found = list(_walk(blob, _path(b.column), conds))
        if b.derive == "count":
            return [Decimal(len(found))]
        cells = [v for _, v in found]
    else:
        raise BindError(f"{b.file} is not a table, a JSON or a YAML file")
    if b.derive == "text":
        if len(cells) != 1 or not isinstance(cells[0], str):
            raise BindError(f"derive=text needs one text cell, the rows select {len(cells)}")
        return cells
    if b.derive == "count_nonnull":
        return [Decimal(sum(1 for c in cells if str(c).strip() != ""))]
    if b.derive == "value" and len(cells) != 1:
        raise BindError(f"the rows select {len(cells)} cells, not one")
    if not cells:
        raise BindError("the rows select no cell")
    nums = [_decimal(c) for c in cells]
    if any(n is None for n in nums):
        raise BindError(f"a selected cell is not a number: {next(c for c, n in zip(cells, nums) if n is None)!r}")
    if b.derive in ("value", "all"):
        return nums
    return [{"sum": sum(nums, Decimal(0)), "min": min(nums), "max": max(nums)}[b.derive]]


def _places(printed: str) -> int:
    tok = printed.rstrip("%").replace(",", "")
    return len(tok.split(".")[1]) if "." in tok else 0


def rounds_to(value: Decimal, printed: str, scale: int = 1, negative: bool = False) -> bool:
    """Does the value, times the scale, round to the printed number at its decimals?"""
    target = Decimal(printed.rstrip("%").replace(",", ""))
    if negative:
        target = -target
    q = Decimal(1).scaleb(-_places(printed))
    v = value * scale
    return target in (v.quantize(q, rounding=ROUND_HALF_UP), v.quantize(q, rounding=ROUND_HALF_EVEN))


def fmt_like(value: Decimal, printed: str, scale: int = 1) -> str:
    v = (value * scale).quantize(Decimal(1).scaleb(-_places(printed)), rounding=ROUND_HALF_UP)
    s = f"{v:,f}" if "," in printed else f"{v:f}"
    return s + ("%" if printed.endswith("%") else "")


def _number_in_text(printed: str, text: str) -> bool:
    plain = printed.rstrip("%").replace(",", "")
    return any(m.group(0).replace(",", "") == plain
               for m in re.finditer(r"(?<![\w.])\d[\d,]*(?:\.\d+)?(?![\w])", text))


@dataclass
class CellResult:
    token: Token | None
    binding: Binding | None
    status: str               # constant, match, mismatch, cannot_bind, not_in_table
    now: str = ""             # the bound value at the printed precision
    values: list = field(default_factory=list)
    why: str = ""


def evaluate(b: Binding, t: Token | None, root: Path = REPO) -> CellResult:
    if b.kind == "constant":
        if not b.reason:
            return CellResult(t, b, "cannot_bind", why="a constant row gives no reason")
        return CellResult(t, b, "constant", now=b.printed)
    if b.kind != "measured":
        return CellResult(t, b, "cannot_bind", why=f"kind '{b.kind}' is neither constant nor measured")
    if b.derive not in DERIVES:
        return CellResult(t, b, "cannot_bind", why=f"derive '{b.derive}' is not one of {sorted(DERIVES)}")
    if b.scale not in (1, 100):
        return CellResult(t, b, "cannot_bind", why=f"scale {b.scale} is neither 1 nor 100")
    percent = t.context.percent if t is not None else b.printed.endswith("%")
    if b.scale == 100 and not percent:
        return CellResult(t, b, "cannot_bind", why="scale 100 on a number not printed as a percentage")
    try:
        cells = resolve(b, root)
    except BindError as e:
        return CellResult(t, b, "cannot_bind", why=str(e))
    negative = t is not None and t.context.before.endswith("−")
    if b.derive == "text":
        ok = _number_in_text(b.printed, cells[0])
        first = next((ln.strip() for ln in cells[0].splitlines() if ln.strip()), "")
        return CellResult(t, b, "match" if ok else "mismatch", now=f"(text) {first[:100]}", values=cells)
    if b.derive == "all":
        ok = all(rounds_to(v, b.printed, b.scale, negative) for v in cells)
        lo, hi = fmt_like(min(cells), b.printed, b.scale), fmt_like(max(cells), b.printed, b.scale)
        now = lo if lo == hi else f"{lo}–{hi}"
        return CellResult(t, b, "match" if ok else "mismatch", now=now, values=cells)
    v = cells[0]
    ok = rounds_to(v, b.printed, b.scale, negative)
    return CellResult(t, b, "match" if ok else "mismatch", now=fmt_like(v, b.printed, b.scale),
                      values=cells)


def cells_results(text: str, bindings: list[Binding], root: Path = REPO) -> list[CellResult]:
    tokens = read_tokens(text)
    by_key = {b.key: b for b in bindings}
    used: set = set()
    out = []
    for t in tokens:
        b = by_key.get((t.line, t.printed, t.occurrence))
        if b is None:
            out.append(CellResult(t, None, "not_in_table"))
            continue
        used.add(b.key)
        if b.excerpt and b.excerpt not in t.excerpt_space:
            out.append(CellResult(t, b, "cannot_bind",
                                  why=f"the excerpt is not in the sentence on line {t.line}"))
            continue
        out.append(evaluate(b, t, root))
    for b in bindings:
        if b.key not in used:
            out.append(CellResult(None, b, "cannot_bind",
                                  why=f"no number {b.printed} (occurrence {b.occurrence}) on line {b.line}"))
    return out


def cells_check(md_path: Path, table: Path, document: str | None = None,
                root: Path | None = None) -> int:
    bindings = read_bindings(table)
    docs = sorted({b.document for b in bindings})
    if document is None:
        if len(docs) != 1:
            raise SystemExit(f"{table.name} binds {', '.join(docs)}: name one with --document")
        document = docs[0]
    bindings = [b for b in bindings if b.document == document]
    res = cells_results(md_path.read_text(), bindings, root or REPO)
    n = {k: sum(1 for r in res if r.status == k)
         for k in ("constant", "match", "mismatch", "cannot_bind", "not_in_table")}
    in_text = sum(1 for r in res if r.token is not None)
    measured = in_text - n["constant"] - n["not_in_table"]
    print(f"cells mode on {md_path.name} ({document}; bindings {table.name}): {in_text} numbers in "
          f"the text; {n['constant']} constants; of {measured} measured, {n['match']} match their bound "
          f"cell, {n['mismatch']} do not, {n['cannot_bind']} cannot be bound; "
          f"{n['not_in_table']} not in the table")
    for r in res:
        b, t = r.binding, r.token
        if r.status == "mismatch":
            print(f"  MISMATCH line {b.line}: {b.printed}, now {r.now}  [{b.where()}]  ...{b.excerpt}...")
    for r in res:
        if r.status == "cannot_bind":
            b = r.binding
            print(f"  CANNOT BIND line {b.line}: {b.printed}  ({r.why})  [{b.where()}]")
    for r in res:
        if r.status == "not_in_table":
            t = r.token
            print(f"  NOT IN TABLE line {t.line}: {t.printed} (occurrence {t.occurrence})  ...{t.sentence[:120]}...")
    return 1 if (n["mismatch"] or n["cannot_bind"] or n["not_in_table"]) else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("manuscript", type=Path)
    ap.add_argument("--pool", action="store_true",
                    help="the weak check: match against every value the pipeline wrote")
    ap.add_argument("--selftest", action="store_true",
                    help="plant numbers no named file holds and confirm the bound mode reports each")
    ap.add_argument("--declared", action="store_true",
                    help="bound mode: also list each declared token with its reason")
    ap.add_argument("--cells", type=Path, metavar="TABLE",
                    help="cells mode: hold each number to the cell this binding table names")
    ap.add_argument("--document", help="cells mode: the binding table's rows for this text "
                                       "(main or supplement)")
    args = ap.parse_args()

    if args.manuscript.suffix == ".md" and not args.pool:
        if args.selftest:
            return selftest(args.manuscript)
        if args.cells:
            return cells_check(args.manuscript, args.cells, args.document)
        rc = bound_check(args.manuscript, show_declared=args.declared)
        _word_counts(args.manuscript.read_text())
        return rc
    print("NOTE: pool mode. It matches each number against every value the pipeline "
          "wrote, so a pass shows only that a number is not impossible. Run the "
          "source-tagged Markdown for the bound check.")
    return pool_check(args.manuscript)


def _word_counts(text: str) -> None:
    start, end = body_span(text)
    body = re.sub(r"\[[^\]]*\]", " ", text[start:end])
    results = body[body.index("# Results"):] if "# Results" in body else body
    words = [s.strip() for s in re.split(r"(?<=[.!?])\s+", results) if WORD_COUNT.search(s)]
    print(f"\nresult sentences asserting a count written as a word ({len(words)}); "
          f"read each against the tables:")
    for s in words:
        print("  - " + " ".join(s.split())[:150])


if __name__ == "__main__":
    sys.exit(main())
