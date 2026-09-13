"""Read-only prose-density lint for the manuscript (.docx or markdown).

Reports, without proposing rewrites:
  1. every sentence longer than 40 words, with section and word count (descending);
  2. every sentence carrying two or more qualifying constructions -- dash asides,
     semicolons, parentheses, and concessive connectors (although / though / while /
     whereas / rather than / not ... but) -- with the count per sentence;
  3. the abstract on its own: sentences per paragraph, mean sentence length, and
     independent claims per sentence (split heuristically on semicolons, dashes and
     coordinated independent clauses);
  4. a per-section summary table (sentences, median length, share over 40 words).

Usage:
    python scripts/prose_density.py path/to/manuscript.docx [--out docs/prose-density.md]
"""
from __future__ import annotations

import argparse
import re
import statistics
from pathlib import Path

HEAD_RE = re.compile(r"^(\d+(?:\.\d+)*\.?\s+\S.*|Abstract|Keywords:?|Data availability|Ethics approval|"
                     r"Funding|Conflict of interest|References|Introduction|Methods|Results|Discussion|"
                     r"Conclusions?|Background|Declarations|Abbreviations|Additional files)$")
CONCESSIVE = re.compile(r"\b(although|though|while|whereas|rather than|nonetheless|nevertheless|albeit|"
                        r"even so|even though)\b", re.I)
NOT_BUT = re.compile(r"\bnot\b[^.;]{0,80}?\bbut\b", re.I)
DASH_ASIDE = re.compile(r"\s[—–]\s")
CLAUSE_SPLIT = re.compile(r";|\s[—–]\s|:\s|,\s(?:and|but|so|while|whereas)\s")


def load_paragraphs(path: Path):
    """[(section, text)] in reading order, references excluded."""
    if path.suffix.lower() == ".docx":
        import docx
        d = docx.Document(path)
        items = []
        for p in d.paragraphs:
            items.append((p.style.name, p.text))
    else:
        items = [("", line) for line in path.read_text(encoding="utf8").splitlines()]
    out, section = [], "front matter"
    for style, text in items:
        t = text.strip()
        if not t:
            continue
        is_heading = style.lower().startswith("heading") or HEAD_RE.match(t) or t.lstrip("#").strip() != t
        if is_heading:
            section = t.lstrip("# ").strip()
            if section.lower().startswith("references"):
                break
            continue
        if section.lower().startswith("references"):
            break
        out.append((section, t))
    return out


def sentences(text: str):
    text = re.sub(r"\s+", " ", text)
    # protect common abbreviations and decimals from the splitter
    prot = re.sub(r"\b(et al|Fig|vs|e\.g|i\.e|cf|approx|no)\.", lambda m: m.group(0).replace(".", "§"), text)
    prot = re.sub(r"(\d)\.(\d)", r"\1¶\2", prot)
    parts = re.split(r"(?<=[.!?])\s+(?=[A-Z(\"'*\[])", prot)
    return [p.replace("§", ".").replace("¶", ".").strip() for p in parts if p.strip()]


def words(s: str) -> int:
    return len(re.findall(r"[A-Za-z0-9ρΔ][\w’'\-.%±/]*", s))


def qualifiers(s: str) -> dict:
    return {
        "dash_asides": len(DASH_ASIDE.findall(s)),
        "semicolons": s.count(";"),
        "parentheses": s.count("("),
        "concessive": len(CONCESSIVE.findall(s)) + len(NOT_BUT.findall(s)),
    }


def claims(s: str) -> int:
    return len([c for c in CLAUSE_SPLIT.split(s) if words(c) >= 4])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("manuscript")
    ap.add_argument("--out", default="docs/prose-density.md")
    ap.add_argument("--long", type=int, default=40)
    a = ap.parse_args()
    paras = load_paragraphs(Path(a.manuscript))

    rows = []            # (section, sentence, n_words, qualifiers)
    for section, text in paras:
        for s in sentences(text):
            rows.append((section, s, words(s), qualifiers(s)))

    # 4. per-section summary
    by_sec: dict[str, list[int]] = {}
    order: list[str] = []
    for sec, _, n, _ in rows:
        if sec not in by_sec:
            by_sec[sec] = []
            order.append(sec)
        by_sec[sec].append(n)
    summary = ["| section | sentences | median words | share > %d words |" % a.long, "|---|---|---|---|"]
    for sec in order:
        ns = by_sec[sec]
        summary.append(f"| {sec} | {len(ns)} | {statistics.median(ns):.0f} | "
                       f"{100 * sum(n > a.long for n in ns) / len(ns):.0f}% |")
    allns = [n for _, _, n, _ in rows]
    summary.append(f"| **all** | {len(allns)} | {statistics.median(allns):.0f} | "
                   f"{100 * sum(n > a.long for n in allns) / len(allns):.0f}% |")

    # 1. long sentences
    long_rows = sorted([r for r in rows if r[2] > a.long], key=lambda r: -r[2])
    long_md = ["| words | section | sentence |", "|---|---|---|"]
    for sec, s, n, _ in long_rows:
        long_md.append(f"| {n} | {sec} | {s.replace('|', '\\|')} |")

    # 2. qualifier-dense sentences
    q_rows = []
    for sec, s, n, q in rows:
        total = sum(q.values())
        if total >= 2:
            q_rows.append((total, sec, s, n, q))
    q_rows.sort(key=lambda r: (-r[0], -r[3]))
    q_md = ["| qualifiers | dash | ; | ( | concessive | words | section | sentence |", "|---|---|---|---|---|---|---|---|"]
    for total, sec, s, n, q in q_rows:
        q_md.append(f"| {total} | {q['dash_asides']} | {q['semicolons']} | {q['parentheses']} | "
                    f"{q['concessive']} | {n} | {sec} | {s.replace('|', '\\|')} |")

    # 3. abstract
    abs_paras = [t for sec, t in paras if sec.lower().startswith("abstract")]
    abs_md = ["| paragraph | sentences | mean words | claims per sentence |", "|---|---|---|---|"]
    abs_detail = ["| paragraph | sentence | words | claims |", "|---|---|---|---|"]
    for i, t in enumerate(abs_paras, 1):
        ss = sentences(t)
        if not ss:
            continue
        cl = [claims(s) for s in ss]
        abs_md.append(f"| {i} | {len(ss)} | {statistics.mean(words(s) for s in ss):.1f} | "
                      f"{', '.join(map(str, cl))} |")
        for s, c in zip(ss, cl):
            abs_detail.append(f"| {i} | {s.replace('|', '\\|')} | {words(s)} | {c} |")

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join([
        f"# Prose density report — `{Path(a.manuscript).name}`", "",
        f"Read-only statistics. Sentences split on terminal punctuation; words are whitespace tokens "
        f"starting with a letter, digit or Greek symbol. Long = more than {a.long} words. "
        f"Qualifier = dash aside, semicolon, parenthesis, or concessive connector "
        f"(although / though / while / whereas / rather than / not…but). Claims per sentence split on "
        f"semicolons, dash asides, colons and coordinated independent clauses.", "",
        "## Summary by section", "", *summary, "",
        f"## 1. Sentences longer than {a.long} words ({len(long_rows)})", "", *long_md, "",
        f"## 2. Sentences with two or more qualifying constructions ({len(q_rows)})", "", *q_md, "",
        "## 3. Abstract", "", *abs_md, "", *abs_detail, "",
    ]), encoding="utf8")
    print(f"{len(rows)} sentences in {len(order)} sections; {len(long_rows)} over {a.long} words; "
          f"{len(q_rows)} with >=2 qualifiers -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
