"""Is every number in the manuscript bound to the output it comes from?

Ported from the companion atlas repository (outputs live in `phase1/reports/phase1/`;
configuration is JSON because PyYAML is not pinned here). It began as a value check:
every numeric token matched against every value the pipeline emits, at half a unit in
the last printed digit. A value check can be satisfied by coincidence. The TP53
orientation AUROC was printed as 0.998, a value no output held, and it matched another
gene's 0.9975.

Every numeric token in the .docx (body paragraphs, table cells and captions) is now one of

  not a measurement  set aside by a named rule in config/manuscript_token_kinds.json:
                     headings and cross-references, labels, identifiers, dates, and the
                     parameters and constants, each rule naming the code that sets the value
  BOUND              its block declares outputs in config/analysis_claims.json, and a named
                     cell of one of them rounds to the printed number, or a derivation the
                     check computes from them does
  POOL               its block declares nothing; the value occurs somewhere in the pipeline
  whitelisted        listed in config/manuscript_number_whitelist.json
  SCOPED MISMATCH    declared outputs, none of which holds it           (fails the run)
  NO SOURCE          nothing declared, and nothing in the pipeline holds it  (fails the run)

and the run states how many measured tokens are bound, to a single cell or to one of
several cells that round to the same printed number.

    python -m src.check_manuscript_numbers path/to/manuscript.docx [--binding-report docs/manuscript-binding.md]
"""
from __future__ import annotations

import argparse
import json
import math
import re
import sys
from decimal import ROUND_HALF_UP, Decimal, InvalidOperation
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]   # <repo>/phase1/src/ -> <repo>
RESULTS = REPO / "phase1" / "reports" / "phase1"
# The analysis set itself is a pipeline output: its row count, its strata sizes
# and the held-out gene's size are quoted in the manuscript and appear in no CSV.
DATA_DIRS = (REPO / "phase1" / "data" / "frozen", REPO / "data" / "external")
WHITELIST = REPO / "phase1" / "config" / "manuscript_number_whitelist.json"
CLAIMS = REPO / "phase1" / "config" / "analysis_claims.json"
FACTS = REPO / "phase1" / "config" / "pipeline_facts.json"
KINDS = REPO / "phase1" / "config" / "manuscript_token_kinds.json"

# Tokens are read as printed, so the tolerance can follow the printed precision.
TOKEN = re.compile(r"(?<![\w.])([-−+]?\d{1,3}(?:,\d{3})*(?:\.\d+)?)(?![\w])")
MAX_ABS = 1e7          # beyond this a token is an accession or an identifier
REF_HEADING = re.compile(r"^\s*references\s*$", re.I)
# Identifiers are not measurements: DOIs, URLs and MaveDB URNs carry digits that
# no pipeline output should be expected to explain.
IDENTIFIER = re.compile(r"(?:https?://\S+|doi:\S+|10\.\d{4,}/\S+|urn:\S+|"
                        r"\bzenodo\.\d+|\bNM_\d+(?:\.\d+)?)", re.I)
# Numbered citations point at the reference list, not at pipeline output.
CITATION = re.compile(r"\[\s*\d{1,2}(?:\s*[,–—-]\s*\d{1,2})*\s*\]")


def _norm(tok: str) -> float:
    return float(tok.replace(",", "").replace("−", "-").replace("+", ""))


def _tolerance(tok: str) -> float:
    frac = tok.split(".")[1] if "." in tok else ""
    return 0.5 * 10 ** (-len(frac)) + 1e-12


def pipeline_values(results: Path = RESULTS) -> np.ndarray:
    """Every numeric value appearing under results/, plus the row count of each
    table -- counts like "2,452 of 2,803 score sets" are genuine pipeline
    quantities that appear nowhere inside the files themselves.

    Generated prose (`Supplemental_Note.md`) is deliberately excluded. It is
    produced by `atlas.supplement`, and any number hardcoded in that generator
    would otherwise validate itself through it: the check would confirm that a
    transcribed literal matches the same transcribed literal.
    """
    vals: set[float] = set()
    for p in sorted(results.rglob("*")):
        if not p.is_file() or p.suffix.lower() not in {".tsv", ".csv", ".json"}:
            continue
        try:
            text = p.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        for m in re.finditer(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", text):
            try:
                v = float(m.group(0))
            except ValueError:
                continue
            if abs(v) <= MAX_ABS:
                vals.add(v)
        if p.suffix.lower() in {".tsv", ".csv"}:
            n = sum(1 for _ in text.splitlines() if _.strip())
            vals.add(float(max(n - 1, 0)))          # data rows, header excluded
            vals.add(float(n))
    # Parquet outputs are binary, so only their shape is readable as a value --
    # and shape is often exactly what the manuscript quotes ("10,888 calls in
    # total" is the row count of the definition sweep).
    for p in [q for base in DATA_DIRS if base.is_dir()
              for q in sorted(base.rglob("*.parquet"))] + sorted(results.rglob("*.parquet")):
        try:
            import pyarrow.parquet as pq
            md = pq.ParquetFile(p).metadata
        except Exception:
            continue
        vals.add(float(md.num_rows))
        vals.add(float(md.num_columns))
        # Category sizes within the analysis set -- the splice strata the
        # manuscript decomposes (583 + 1,191 + 7 = 1,781) live only here.
        try:
            import pandas as pd
            df = pd.read_parquet(p)
            for col in df.columns:
                if df[col].dtype == object or str(df[col].dtype) == "bool":
                    if df[col].nunique(dropna=True) <= 40:
                        for n in df[col].value_counts().values:
                            vals.add(float(n))
                        for keep in (True,):
                            if str(df[col].dtype) == "bool":
                                vals.add(float(df[col].sum()))
        except Exception:
            pass
    # Counts the manuscript cites that live in no results/ file -- the guardrail
    # suite size above all. Without them a stale count matches an unrelated
    # value by coincidence and is reported as verified.
    if FACTS.exists():
        for v in json.loads(FACTS.read_text()).values():
            if isinstance(v, (int, float)) and not isinstance(v, bool):
                vals.add(float(v))
    return np.array(sorted(vals))


def manuscript_tokens(docx_path: Path) -> list[dict]:
    """Tokens from body paragraphs, table cells and captions; the reference
    list is excluded because its years and page ranges are not measurements."""
    import docx

    d = docx.Document(str(docx_path))
    units: list[tuple[str, str, str]] = []
    in_refs = False
    for i, p in enumerate(d.paragraphs):
        if REF_HEADING.match(p.text.strip()):
            in_refs = True
        if in_refs or not p.text.strip():
            continue
        units.append((f"P{i}", p.text, p.style.name if p.style is not None else ""))
    for ti, tb in enumerate(d.tables):
        for ri, row in enumerate(tb.rows):
            for ci, cell in enumerate(row.cells):
                if cell.text.strip():
                    units.append((f"T{ti + 1}r{ri}c{ci}", cell.text, "table header" if ri == 0 else "table"))

    out: list[dict] = []
    seen: set[tuple[str, str]] = set()
    for uid, raw, style in units:
        # blank identifiers and citations with spaces of the same length, so a token's
        # offsets are offsets into the printed text the kind rules match against
        blank = lambda m: " " * len(m.group(0))
        text = CITATION.sub(blank, IDENTIFIER.sub(blank, raw))
        for m in TOKEN.finditer(text):
            tok = m.group(1)
            try:
                val = _norm(tok)
            except ValueError:
                continue
            if abs(val) > MAX_ABS:
                continue
            key = (uid, tok)
            if key in seen:
                continue
            seen.add(key)
            lo = max(0, m.start() - 55)
            out.append({"uid": uid, "token": tok, "value": val, "start": m.start(1), "end": m.end(1),
                        "style": style,
                        "context": re.sub(r"\s+", " ", text[lo:m.end() + 35]).strip()})
    return out


def load_whitelist(path: Path = WHITELIST) -> tuple[set[str], list[re.Pattern], dict]:
    spec = json.loads(path.read_text()) if path.exists() else {}
    values: set[str] = set()
    patterns: list[re.Pattern] = []
    for key, group in spec.items():
        if key.startswith("_") or not isinstance(group, dict):
            continue   # `_note` and other prose keys are documentation
        for v in group.get("values", []) or []:
            values.add(str(v))
        if group.get("pattern"):
            patterns.append(re.compile(group["pattern"]))
    return values, patterns, spec


def heading_numbers(docx_path: Path) -> set[str]:
    """Section numbers the document actually has ("2.4", "3", ...).

    The whitelist exempts section numbers by pattern -- a digit or two, a point, a
    digit -- and a pattern cannot tell the cross-reference "(2.4)" from the
    measurement "8.8". It let a stale Table 3 value, 8.8%, through the v1 -> v2
    replacement pass as a section number. The pattern is therefore honoured only for
    numbers that open a heading paragraph.
    """
    import docx

    out: set[str] = set()
    for p in docx.Document(str(docx_path)).paragraphs:
        if p.style is not None and p.style.name.lower().startswith("heading"):
            m = re.match(r"\s*(\d{1,2}(?:\.\d)?)\.?\s", p.text)
            if m:
                out.add(m.group(1))
    return out


# A token printed to few decimals sits in a dense part of the pool, so a match
# carries no information: "0.12" is within half a last-digit of hundreds of
# unrelated pipeline values. Such matches are reported separately as `weak`,
# because treating them as verification is exactly the mistake that lets a
# sourceless number pass.
WEAK_MATCH_MIN = 25          # distinct pool values inside the tolerance window


NUMBER_WORDS = {
    "one": 1, "two": 2, "three": 3, "four": 4, "five": 5, "six": 6, "seven": 7,
    "eight": 8, "nine": 9, "ten": 10, "eleven": 11, "twelve": 12,
    "thirteen": 13, "fourteen": 14, "fifteen": 15, "sixteen": 16,
    "seventeen": 17, "eighteen": 18, "nineteen": 19, "twenty": 20,
}
# A claim may assert how many things an output covers. That is a statement about
# the output's SHAPE, and no value-membership check can test it: "the nineteen
# predictors span two axes" cites a matrix of nine, and every number printed in
# that sentence still corresponds to some cell. Cardinality is the check that
# catches it.
CARDINALITY = re.compile(
    r"\b(" + "|".join(NUMBER_WORDS) + r"|\d{1,3})\s+"
    r"(?:[a-z-]+\s+){0,3}?"                    # "nine broad-scope and splice-aware predictors"
    r"(predictors?|models?|tools?|genes?|strata|stratum|columns?|rows?)\b", re.I)


def paragraph_cardinalities(docx_path: Path) -> dict[str, list[tuple[int, str, str]]]:
    """paragraph id -> [(count, noun, phrase)] for every cardinality asserted."""
    import docx

    d = docx.Document(str(docx_path))
    out: dict[str, list[tuple[int, str, str]]] = {}
    in_refs = False
    for i, para in enumerate(d.paragraphs):
        if REF_HEADING.match(para.text.strip()):
            in_refs = True
        if in_refs or not para.text.strip():
            continue
        for m in CARDINALITY.finditer(para.text):
            tok = m.group(1).lower()
            n = NUMBER_WORDS.get(tok, None)
            if n is None:
                try:
                    n = int(tok)
                except ValueError:
                    continue
            out.setdefault(f"P{i}", []).append((n, m.group(2).lower(), m.group(0)))
    return out


def _blocks(docx_path: Path) -> list[tuple[str, str]]:
    """(id, text) for every paragraph and table cell, in document order."""
    import docx

    d = docx.Document(str(docx_path))
    out = [(f"P{i}", para.text) for i, para in enumerate(d.paragraphs)]
    for ti, t in enumerate(d.tables):
        for ri, row in enumerate(t.rows):
            for ci, cell in enumerate(row.cells):
                out.append((f"T{ti+1}r{ri+1}c{ci+1}", cell.text))
    return out


def resolve_anchors(docx_path: Path, path: Path = CLAIMS) -> dict[int, str]:
    """Row index in analysis_claims.json -> the block id its anchor lands in.

    Claims bind by a unique substring of the sentence that makes them, not by
    paragraph number. Rewriting the abstract shifted every later paragraph by
    two and left all six `para` values pointing at the wrong text with nothing
    to notice. A row whose anchor hits zero blocks, or more than one, is
    omitted here and reported by the tests rather than guessed at.
    """
    spec = json.loads(path.read_text()) if path.exists() else {}
    blocks = _blocks(docx_path)
    resolved: dict[int, str] = {}
    for idx, c in enumerate(spec.get("claims", [])):
        anchor = c.get("anchor")
        if not anchor:
            continue
        total = sum(t.count(anchor) for _, t in blocks)
        if total == 1:
            resolved[idx] = next(uid for uid, t in blocks if anchor in t)
    return resolved


def claim_cardinalities(docx_path: Path, path: Path = CLAIMS) -> dict[str, list[dict]]:
    """paragraph id -> the cardinality assertions recorded for it."""
    if not path.exists():
        return {}
    spec = json.loads(path.read_text()) if path.exists() else {}
    where = resolve_anchors(docx_path, path)
    out: dict[str, list[dict]] = {}
    for idx, c in enumerate(spec.get("claims", [])):
        if idx not in where:
            continue
        if c.get("cardinality"):
            out.setdefault(where[idx], []).append(
                {"count": int(c["cardinality"]["count"]),
                 "noun": c["cardinality"]["noun"],
                 "output": c.get("output")})
    return out


def claim_derivations(docx_path: Path, path: Path = CLAIMS) -> dict[str, list[dict]]:
    """block id -> values the paragraph states that no output stores directly.

    A summary of stored values -- a sum, a median, a count of cells -- has no
    literal source, and widening the scope cannot give it one. Each such value
    is declared individually, with the arithmetic that produces it, so it can be
    read and checked. This is deliberately not a paragraph-level exemption: a
    row here licenses exactly one number and says where it comes from.
    """
    spec = json.loads(path.read_text()) if path.exists() else {}
    where = resolve_anchors(docx_path, path)
    out: dict[str, list[dict]] = {}
    for idx, c in enumerate(spec.get("claims", [])):
        if c.get("derived") and idx in where:
            out.setdefault(where[idx], []).extend(c["derived"])
    return out


def claim_pins(docx_path: Path, path: Path = CLAIMS) -> dict[str, dict[str, str]]:
    """block id -> {printed token: the one cell it must bind to}.

    A paragraph's scope can hold several cells that round to the same printed number:
    TP53's post-orientation control AUROC, 0.997, shares its paragraph's scope with the
    per-variant SpliceAI scores, twelve of which also print as 0.997. A pin names the cell
    the number is, so the token binds to that cell or is reported as a mismatch."""
    spec = json.loads(path.read_text()) if path.exists() else {}
    where = resolve_anchors(docx_path, path)
    out: dict[str, dict[str, str]] = {}
    for idx, c in enumerate(spec.get("claims", [])):
        if c.get("pin") and c.get("status") == "verified" and idx in where:
            out.setdefault(where[idx], {}).update(c["pin"])
    return out


def claim_scopes(docx_path: Path, path: Path = CLAIMS) -> dict[str, list[str]]:
    """paragraph id -> the outputs config/analysis_claims.yaml binds it to.

    A number inside a paragraph that declares its own source should be checked
    against THAT source, not against every value the pipeline has ever emitted.
    Checking against the whole pool is what let the complementarity range pass:
    0.12 matches 299 unrelated values, so it never reached `no source`.
    """
    if not path.exists():
        return {}
    spec = json.loads(path.read_text()) if path.exists() else {}
    where = resolve_anchors(docx_path, path)
    scopes: dict[str, list[str]] = {}
    for idx, c in enumerate(spec.get("claims", [])):
        if idx not in where:
            continue
        out = c.get("output")
        if out and c.get("status") == "verified":
            # a claim may name several outputs; all of them scope the block
            outs = out if isinstance(out, list) else [out]
            scopes.setdefault(where[idx], []).extend(outs)
            # a table caption's claim may also scope the table that follows it
            if c.get("table") and where[idx].startswith("P"):
                for cell_uid in _table_after(docx_path, int(where[idx][1:])):
                    scopes.setdefault(cell_uid, []).extend(outs)
    return scopes


def _table_after(docx_path: Path, para_index: int) -> list[str]:
    """Token ids of every cell of the first table after body paragraph `para_index`."""
    import docx
    from docx.oxml.ns import qn
    from docx.table import Table

    d = docx.Document(str(docx_path))
    pi = ti = 0
    seen = False
    for child in d.element.body.iterchildren():
        if child.tag == qn("w:p"):
            seen = seen or pi == para_index
            pi += 1
        elif child.tag == qn("w:tbl"):
            if seen:
                tb = Table(child, d)
                return [f"T{ti + 1}r{ri}c{ci}" for ri, row in enumerate(tb.rows) for ci in range(len(row.cells))]
            ti += 1
    return []


def scoped_values(outputs: list[str]) -> np.ndarray:
    """Every numeric value in the specific outputs a claim names."""
    vals: set[float] = set()
    for rel in outputs:
        p = REPO / rel
        if not p.exists():
            continue
        if p.suffix.lower() == ".parquet":
            try:
                import pyarrow.parquet as pq
                md = pq.ParquetFile(p).metadata
                vals.add(float(md.num_rows)); vals.add(float(md.num_columns))
                tbl = pq.read_table(p)
                for col in tbl.column_names:
                    for v in tbl.column(col).to_pylist()[:200000]:
                        if isinstance(v, (int, float)) and not isinstance(v, bool) \
                                and abs(v) <= MAX_ABS:
                            vals.add(float(v))
            except Exception:
                pass
            continue
        try:
            text = p.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        for m in re.finditer(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", text):
            try:
                v = float(m.group(0))
            except ValueError:
                continue
            if abs(v) <= MAX_ABS:
                vals.add(v)
        if p.suffix.lower() in {".tsv", ".csv"}:
            n = sum(1 for _ in text.splitlines() if _.strip())
            vals.add(float(max(n - 1, 0))); vals.add(float(n))
    return np.array(sorted(vals))


def _hits(pool: np.ndarray, v: float, tol: float) -> int:
    """Matches for a token, trying it as printed and as a percentage.

    Both signs are tried. The token's magnitude is what reaches here, so a
    pipeline value of -0.7303 would never match a printed -0.73 if only the
    positive form were searched; the global fallback used to hide that, and
    removing the fallback exposed it.
    """
    if not pool.size:
        return 0
    n = 0
    for cand in (v, -v):
        n += int(np.sum(np.abs(pool - cand) <= tol))
        n += int(np.sum(np.abs(pool - cand / 100.0) <= tol / 100))
    return n


def token_kinds(docx_path: Path, tokens: list[dict], path: Path = KINDS,
                headings: set[str] | None = None) -> tuple[dict[int, tuple[str, str]], list[dict]]:
    """Which tokens are not measurements, and the rule that says so.

    Returns {token index: (kind, rule name)} and the rules whose declared source no
    longer defines the value (a parameter rule names the line of code that sets it)."""
    import docx

    spec = json.loads(path.read_text()) if path.exists() else {}
    d = docx.Document(str(docx_path))
    texts = {f"P{i}": p.text for i, p in enumerate(d.paragraphs)}
    for ti, tb in enumerate(d.tables):
        for ri, row in enumerate(tb.rows):
            for ci, cell in enumerate(row.cells):
                texts[f"T{ti + 1}r{ri}c{ci}"] = cell.text
    heads = heading_numbers(docx_path) if headings is None else headings
    rules = [(r, re.compile(r["pattern"])) for r in spec.get("rules", [])]
    problems = []
    for r, _ in rules:
        if r.get("source"):
            src = REPO / r["source"]
            if not src.exists() or not re.search(r["defined_by"], src.read_text(errors="ignore")):
                problems.append({"rule": r["name"], "source": r["source"], "defined_by": r["defined_by"]})
    out: dict[int, tuple[str, str]] = {}
    for i, t in enumerate(tokens):
        text = texts.get(t["uid"], "")
        if t.get("style", "").lower().startswith("heading"):
            out[i] = ("heading number", "numbered heading")
            continue
        bare = t["token"].lstrip("+-−")
        if re.fullmatch(r"\d{1,2}\.\d", bare) and bare in heads:
            before = text[max(0, t["start"] - 9):t["start"]]
            after = text[t["end"]:t["end"] + 1]
            if re.search(r"(?:\(|\bin |\bof |, |; |\bsee )$", before) or after in (")", ";"):
                out[i] = ("section reference", f"cross-reference to section {bare}")
                continue
        for r, rx in rules:
            spans = []
            for m in rx.finditer(text):
                groups = [m.span(g) for g in range(1, rx.groups + 1) if m.group(g) is not None]
                spans.extend(groups or [m.span(0)])
            if any(lo <= t["start"] and t["end"] <= hi for lo, hi in spans):
                out[i] = (r["kind"], r["name"])
                break
    return out, problems


_NUM = re.compile(r"[-+−]?\d*\.?\d+(?:[eE][-+]?\d+)?")


def cell_index(outputs: list[str]) -> list[tuple[str, float]]:
    """(cell id, value) for every number in the declared outputs.

    CSV/TSV: every cell, and each number inside a text cell (an interval string gives
    two), plus the row count. JSON: every numeric leaf by path. Parquet: row and
    column counts, category counts of text and boolean columns, and every numeric cell
    of a table small enough to enumerate."""
    import pandas as pd

    cells: list[tuple[str, float]] = []
    for rel in outputs:
        p = REPO / rel
        if not p.exists():
            continue
        name, suf = p.name, p.suffix.lower()
        if suf in (".csv", ".tsv"):
            df = pd.read_csv(p, sep="\t" if suf == ".tsv" else ",", dtype=str, keep_default_na=False)
            cells.append((f"{name}: row count", float(len(df))))
            # a row is named by up to three of its text-valued columns (condition,
            # calibration, object ...); numeric-text columns such as intervals are not names
            textcols = [c for c in df.columns
                        if df[c].map(lambda x: bool(re.search(r"[A-DF-Za-df-z_]", str(x)))).mean() > 0.5][:3]
            for i in range(len(df)):
                key = "/".join(str(df.iloc[i][c]) for c in textcols) or str(i)
                for c in df.columns:
                    val = str(df.iloc[i][c]).replace("−", "-")
                    # a number is read out of a text cell only when the cell is numeric
                    # text (an interval such as "[0.0016, 0.0239]"); digits inside a label
                    # like y_assay/BRCA1_included or fusion_M1 are not values
                    if re.search(r"[A-DF-Za-df-z_]", val):
                        continue
                    nums = _NUM.findall(val)
                    for k, n in enumerate(nums):
                        try:
                            x = float(n.replace("−", "-"))
                        except ValueError:
                            continue
                        if abs(x) <= MAX_ABS:
                            tag = "" if len(nums) == 1 else f"#{k + 1}"
                            cells.append((f"{name}[{key}].{c}{tag}", x))
        elif suf == ".json":
            def walk(o, path):
                if isinstance(o, bool):
                    return
                if isinstance(o, (int, float)):
                    cells.append((f"{name}:{path}", float(o)))
                elif isinstance(o, dict):
                    for k, v in o.items():
                        walk(v, f"{path}.{k}" if path else k)
                elif isinstance(o, list):
                    for k, v in enumerate(o):
                        walk(v, f"{path}[{k}]")
            walk(json.loads(p.read_text()), "")
        elif suf == ".parquet":
            df = pd.read_parquet(p)
            cells.append((f"{name}: row count", float(len(df))))
            cells.append((f"{name}: column count", float(df.shape[1])))
            for c in df.columns:
                dt = str(df[c].dtype)
                if dt in ("object", "bool", "category", "string"):
                    vc = df[c].value_counts(dropna=True)
                    if len(vc) <= 40:
                        for cat, n in vc.items():
                            cells.append((f"{name}: count of {c} = {cat}", float(n)))
                elif len(df) <= 5000:
                    for i, x in enumerate(df[c].to_numpy()):
                        try:
                            xf = float(x)
                        except (TypeError, ValueError):
                            continue
                        if np.isfinite(xf):
                            cells.append((f"{name}[{i}].{c}", xf))
    return cells


def _rounds_to(x: float, tok: str, as_percent: bool = False) -> bool:
    """Does x print as the token? x (times 100 for a percentage) is rounded to the token's
    decimals, half-up on its shortest decimal form or as Python formats the float; at a
    decimal tie the two differ and either print is correct. A symmetric window of half a
    unit is not this test: it let BRCA2's control AUROC 0.9975, which prints as 0.998, bind
    a printed 0.997 -- the coincidence that once let 0.998 pass for TP53's 0.9969."""
    if not math.isfinite(x):
        return False
    bare = tok.replace(",", "").lstrip("+-−")
    dp = len(bare.split(".")[1]) if "." in bare else 0
    xd, xf = abs(Decimal(repr(float(x)))), abs(float(x))
    if as_percent:
        xd, xf = xd * 100, xf * 100
    try:
        return (xd.quantize(Decimal(1).scaleb(-dp), ROUND_HALF_UP) == Decimal(bare)
                or Decimal(f"{xf:.{dp}f}") == Decimal(bare))
    except InvalidOperation:
        return False


def _bind(index: list[tuple[str, float]], tok: str, as_percent: bool = True) -> list[str]:
    """Cells whose value rounds to the token as printed, or -- for a token printed as a
    percentage -- to the token over 100; each named with the value it holds."""
    v, tol = abs(_norm(tok)), _tolerance(tok)
    out = []
    for cid, x in index:
        # the window is a fast pre-filter; _rounds_to decides
        if abs(abs(x) - v) <= tol and _rounds_to(x, tok):
            out.append(f"{cid} = {x:.6g}")
        elif as_percent and abs(abs(x) - v / 100.0) <= tol / 100.0 and _rounds_to(x, tok, True):
            out.append(f"{cid} = {x:.6g} (as %)")
    return out


def _printed_as_percent(t: dict, texts: dict[str, str]) -> bool:
    """A token is read as a percentage only if it is printed as one: followed by %,
    near "percentage point(s)", or in a table column whose header carries %."""
    text = texts.get(t["uid"], "")
    after = text[t["end"]:t["end"] + 25]
    before = text[max(0, t["start"] - 45):t["start"]]
    if re.match(r"\s?%", after) or "percentage point" in after or "percentage point" in before:
        return True
    m = re.match(r"(T\d+)r\d+c(\d+)$", t["uid"])
    return bool(m and "%" in texts.get(f"{m.group(1)}r0c{m.group(2)}", ""))


def _load_table(rel: str):
    import pandas as pd
    p = REPO / rel
    if p.suffix.lower() == ".parquet":
        return pd.read_parquet(p)
    return pd.read_csv(p, sep="\t" if p.suffix.lower() == ".tsv" else ",")


def _eval_spec(spec: dict) -> float:
    """One number computed from a declared output after its `where` filters: count (rows);
    sum / min / max / mean / median (of `column`); rank (1-based, descending by `column`, of
    the row whose `key_column` is `key`); off_diagonal (the cells of a confusion table whose
    row label differs from their column label); min_ratio / max_ratio (the smallest or
    largest ratio of two columns after pivoting)."""
    df = _load_table(spec["from"])
    for col, want in (spec.get("where") or {}).items():
        df = df[df[col].astype(str).str.lower() == str(want).lower()]
    agg = spec["agg"]
    if agg == "count":
        return float(len(df))
    if agg == "rank":
        keycol = spec.get("key_column", df.columns[0])
        order = df.sort_values(spec["column"], ascending=False)[keycol].astype(str).tolist()
        return float(order.index(spec["key"]) + 1)
    if agg == "off_diagonal":
        rows = df.set_index(df.columns[0])
        return float(sum(rows.at[r, c] for r in rows.index for c in rows.columns if str(r) != str(c)))
    if agg in ("min_ratio", "max_ratio"):
        # one row per `index`, one column per value of `columns`, holding `values`
        pv = df.pivot_table(index=spec["index"], columns=spec["columns"], values=spec["values"], aggfunc="first")
        ratio = (pv[spec["numerator"]] / pv[spec["denominator"]]).replace([np.inf, -np.inf], np.nan).dropna()
        return float(ratio.min() if agg == "min_ratio" else ratio.max())
    col = df[spec["column"]].astype(float)
    if agg in ("sum", "min", "max", "mean", "median"):
        return float(getattr(col, agg)())
    raise ValueError(f"unknown aggregate {agg!r}")


def _derived_value(dv: dict) -> tuple[float, str]:
    """A derivation the check evaluates against the outputs it names. An entry without
    `agg` is only declared -- its value is taken on trust -- and is reported as such."""
    agg = dv.get("agg")
    if agg is None:
        return float(dv["value"]), f"{dv['how']} (declared, not computed)"
    if agg in ("ratio_percent", "complement_percent"):
        r = 100.0 * _eval_spec(dv["numerator"]) / _eval_spec(dv["denominator"])
        return (r if agg == "ratio_percent" else 100.0 - r), dv.get("how", agg)
    return _eval_spec(dv), dv.get("how", agg)


def classify(docx_path: Path, results: Path = RESULTS,
             whitelist: Path = WHITELIST, claims: Path = CLAIMS, kinds: Path = KINDS) -> dict:
    """Every numeric token is one of:

      not a measurement   a heading or cross-reference number, a label, an identifier, a
                          date, or a parameter or constant named by a rule in
                          config/manuscript_token_kinds.json
      BOUND               its paragraph (or table) declares outputs, and a named cell of
                          one of them rounds to the printed number, or a derivation does;
                          a pin in the claim narrows the candidates to one named cell
      POOL                no declared outputs; the value occurs somewhere in the pipeline
      whitelisted         listed in config/manuscript_number_whitelist.json
      SCOPED MISMATCH     declared outputs, and none of their cells holds it  (fails)
      NO SOURCE           undeclared, and nothing in the pipeline holds it    (fails)
    """
    pool = pipeline_values(results)
    wl_values, wl_patterns, wl_spec = load_whitelist(whitelist)
    sect_src = (wl_spec.get("section_numbers") or {}).get("pattern")
    headings = heading_numbers(docx_path)
    scopes = claim_scopes(docx_path, claims)
    derivations = claim_derivations(docx_path, claims)
    pins = claim_pins(docx_path, claims)
    cell_cache: dict[str, list[tuple[str, float]]] = {}

    # Cardinality: does the count a paragraph asserts match what the claims
    # manifest records for the output it cites? No value check can ask this.
    asserted = paragraph_cardinalities(docx_path)
    recorded = claim_cardinalities(docx_path, claims)
    cardinality_mismatch = []
    for para, entries in recorded.items():
        for e in entries:
            for n, noun, phrase in asserted.get(para, []):
                if noun.rstrip("s") != e["noun"].rstrip("s"):
                    continue
                if n != e["count"]:
                    cardinality_mismatch.append(
                        {"uid": para, "phrase": phrase, "asserted": n,
                         "recorded": e["count"], "noun": noun, "output": e["output"]})

    tokens = manuscript_tokens(docx_path)
    kind_of, rule_problems = token_kinds(docx_path, tokens, kinds, headings)
    import docx as _docx
    _d = _docx.Document(str(docx_path))
    texts = {f"P{i}": p.text for i, p in enumerate(_d.paragraphs)}
    for ti, tb in enumerate(_d.tables):
        for ri, row in enumerate(tb.rows):
            for ci, cell in enumerate(row.cells):
                texts[f"T{ti + 1}r{ri}c{ci}"] = cell.text
    matched, weak, no_source, whitelisted, scoped_mismatch = [], [], [], [], []
    bound, pool_only, non_measurement = [], [], []
    for i, t in enumerate(tokens):
        tol = _tolerance(t["token"])
        v = abs(t["value"])
        bare = t["token"].lstrip("+-−")
        white = (bare in wl_values
                 or any(p.match(bare) for p in wl_patterns if p.pattern != sect_src)
                 or (sect_src is not None and re.match(sect_src, bare) is not None
                     and bare in headings))
        if i in kind_of:
            t["kind"], t["rule"] = kind_of[i]
            non_measurement.append(t)
            continue
        outputs = scopes.get(t["uid"])
        if outputs:
            t["scope"] = outputs
            # a declared derivation names the value's origin exactly, so it is tried first
            for dv in derivations.get(t["uid"], []):
                val, how = _derived_value(dv)
                if _rounds_to(val, t["token"]):
                    t["derived_from"] = how
                    t["cells"], t["n_cells"] = [f"derived: {how} = {val:.6g}"], 1
                    break
            if not t.get("cells"):
                key = "|".join(sorted(outputs))
                if key not in cell_cache:
                    cell_cache[key] = cell_index(outputs)
                hits = _bind(cell_cache[key], t["token"], _printed_as_percent(t, texts))
                pin = pins.get(t["uid"], {}).get(t["token"])
                if pin is not None:
                    t["pin"] = pin
                    hits = [h for h in hits if h.split(" = ")[0] == pin]
                if hits:
                    t["cells"], t["n_cells"] = hits[:6], len(hits)
            if t.get("cells"):
                t["n_scoped_matches"] = t["n_pool_matches"] = t["n_cells"]
                bound.append(t)
                matched.append(t)
                continue
            if white:
                whitelisted.append(t)
                continue
            # Absent from the output its own paragraph declares. This is a fault: a
            # paragraph that declares a source asserts that its numbers come from it.
            # A number that belongs to a different output means the scope is too
            # narrow; widen the declaration in config/analysis_claims.json.
            t["n_scoped_matches"] = 0
            t["n_pool_matches"] = _hits(pool, v, tol)
            scoped_mismatch.append(t)
            continue
        n_hits = _hits(pool, v, tol)
        t["n_pool_matches"] = n_hits
        if n_hits:
            (weak if n_hits >= WEAK_MATCH_MIN else matched).append(t)
            pool_only.append(t)
        elif white:
            whitelisted.append(t)
        else:
            no_source.append(t)
    # a pin no scoped token took names a cell for a number the block no longer prints
    applied = {(t["uid"], t["token"]) for t in tokens if t.get("pin")}
    pin_problems = [{"uid": uid, "token": tok, "cell": cell} for uid, pinned in pins.items()
                    for tok, cell in pinned.items() if (uid, tok) not in applied]
    measured = len(bound) + len(pool_only) + len(whitelisted) + len(no_source) + len(scoped_mismatch)
    single = sum(1 for t in bound if t["n_cells"] == 1)
    by_kind: dict[str, int] = {}
    for t in non_measurement:
        by_kind[t["kind"]] = by_kind.get(t["kind"], 0) + 1
    return {"matched": matched, "weak": weak, "whitelisted": whitelisted,
            "no_source": no_source, "scoped_mismatch": scoped_mismatch,
            "cardinality_mismatch": cardinality_mismatch, "stale_counts": [],
            "bound": bound, "pool": pool_only, "non_measurement": non_measurement,
            "non_measurement_by_kind": by_kind, "rule_problems": rule_problems, "pin_problems": pin_problems,
            "coverage": {"tokens": len(tokens), "measured": measured, "bound": len(bound),
                         "single_cell": single, "several_cells": len(bound) - single},
            "n_pool": int(pool.size), "n_scoped_paragraphs": len(scopes)}


def binding_report(r: dict, docx_name: str) -> str:
    """Markdown record of every token's classification, for docs/."""
    cov = r["coverage"]
    _md = lambda x: x.replace("|", "\\|")   # a | in the context would split the table row
    lines = [f"# Number binding: `{docx_name}`", "",
             f"{cov['measured']} of the {cov['tokens']} numeric tokens are measurements. {cov['bound']} of them are "
             "bound to a named cell of an output their paragraph or table declares (or to a derivation the check "
             f"computes from such outputs): {cov['single_cell']} to a single cell, and {cov['several_cells']} to a cell "
             "that the declared outputs share with at least one other cell rounding to the same printed number, so the "
             f"check confirms the value but cannot say which of those cells the text means. {cov['measured'] - cov['bound']} "
             f"are unbound. The other {len(r['non_measurement'])} tokens are not measurements and are listed with the "
             "rule that sets each aside. Written by `python -m src.check_manuscript_numbers <docx> --binding-report "
             "<path>`; a token is one distinct number per paragraph or table cell.", "",
             "## Measured numbers not bound to a cell", ""]
    unbound = ([("pool only", t) for t in r["pool"]] + [("whitelisted", t) for t in r["whitelisted"]]
               + [("NO SOURCE", t) for t in r["no_source"]] + [("SCOPED MISMATCH", t) for t in r["scoped_mismatch"]])
    if unbound:
        lines += ["| block | token | status | context |", "|---|---|---|---|"]
        lines += [f"| {t['uid']} | {t['token']} | {st} | …{_md(t['context'])}… |" for st, t in unbound]
    else:
        lines.append("None.")
    lines += ["", "## Bound numbers", "",
              "Single-cell bindings first, then those with several candidate cells.", "",
              "| block | token | cells holding the value (first six) | context |", "|---|---|---|---|"]
    for t in sorted(r["bound"], key=lambda t: t["n_cells"] > 1):
        more = f" (+{t['n_cells'] - len(t['cells'])} more)" if t["n_cells"] > len(t["cells"]) else ""
        lines.append(f"| {t['uid']} | {t['token']} | {'; '.join(t['cells'])}{more}{' (pinned)' if t.get('pin') else ''} | …{_md(t['context'])}… |")
    lines += ["", "## Not measurements", "", "| block | token | kind | rule | context |", "|---|---|---|---|---|"]
    lines += [f"| {t['uid']} | {t['token']} | {t['kind']} | {t['rule']} | …{_md(t['context'])}… |" for t in r["non_measurement"]]
    return "\n".join(lines) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("docx")
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--json", dest="as_json")
    ap.add_argument("--binding-report", dest="binding_report")
    a = ap.parse_args()
    r = classify(Path(a.docx))
    cov = r["coverage"]
    print(f"[numbers] {Path(a.docx).name}: {cov['tokens']} numeric tokens against "
          f"{r['n_pool']:,} pipeline values, {r['n_scoped_paragraphs']} blocks scoped to declared outputs")
    print(f"  measured                   : {cov['measured']}")
    print(f"    bound to a named cell    : {cov['bound']}  (one cell {cov['single_cell']}, "
          f"several candidate cells {cov['several_cells']})")
    print(f"    pool match only          : {len(r['pool'])}")
    print(f"    whitelisted              : {len(r['whitelisted'])}")
    print(f"    NO SOURCE                : {len(r['no_source'])}")
    print(f"    SCOPED MISMATCH          : {len(r['scoped_mismatch'])}")
    kinds = ", ".join(f"{k} {n}" for k, n in sorted(r["non_measurement_by_kind"].items()))
    print(f"  not measurements           : {len(r['non_measurement'])}  ({kinds})")
    print(f"  CARDINALITY MISMATCH       : {len(r['cardinality_mismatch'])}")
    print(f"  RULE SOURCE MISSING        : {len(r['rule_problems'])}")
    print(f"  PIN NOT APPLIED            : {len(r['pin_problems'])}")
    print(f"{cov['bound']} of {cov['measured']} measured tokens are bound to a named cell of a declared output: "
          f"{cov['single_cell']} to a single cell, {cov['several_cells']} to one of several cells that round to the "
          f"printed number; {cov['measured'] - cov['bound']} unbound.")
    for c in r["cardinality_mismatch"]:
        print(f"    CARDINALITY {c['uid']}: text says \"{c['phrase']}\" but "
              f"{c['output']} carries {c['recorded']} {c['noun']}")
    for q in r["rule_problems"]:
        print(f"    RULE SOURCE MISSING  {q['rule']}: {q['source']} no longer matches {q['defined_by']!r}")
    for q in r["pin_problems"]:
        print(f"    PIN NOT APPLIED  {q['uid']} {q['token']}: no scoped token of that block takes {q['cell']}")
    for t in r["scoped_mismatch"]:
        src = f"the pinned cell {t['pin']}" if t.get("pin") else ", ".join(t["scope"])
        print(f"    SCOPED MISMATCH {t['uid']:>8s}  {t['token']:>10s}  absent from {src}\n        …{t['context']}…")
    for label, key in (("pool only", "pool"), ("whitelisted", "whitelisted"), ("NO SOURCE", "no_source")):
        for t in r[key]:
            print(f"    unbound: {label:11s} {t['uid']:>8s}  {t['token']:>10s}   …{t['context'][:90]}…")
    if a.verbose:
        for t in r["non_measurement"]:
            print(f"    not measured  {t['uid']:>8s}  {t['token']:>10s}  {t['kind']} ({t['rule']})")
    if a.as_json:
        Path(a.as_json).write_text(json.dumps(r, indent=1, default=str))
    if a.binding_report:
        Path(a.binding_report).write_text(binding_report(r, Path(a.docx).name))
    sys.exit(1 if (r["no_source"] or r["cardinality_mismatch"]
                   or r["scoped_mismatch"] or r["rule_problems"] or r["pin_problems"]) else 0)


if __name__ == "__main__":
    main()
