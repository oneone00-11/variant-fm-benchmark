"""E5 (E19.0, E20.0, E17.0 rule 3): tests/test_no_secrets.py's patterns, unchanged, over the
files this stage writes, by scripts/stage_E18/scan_secrets.py's method: this stage's
scripts and run logs (scripts/stage_E5), the documents and tests it adds or changes, the
numbers gate, and the layout stage's code and outputs once they exist. Prints every match
(the pattern's name, never the matched text) and the count of files scanned; exits 1 on
a match. Extra folders or files (a preview made outside the repository, say) can be named
on the command line.

    .venv/bin/python scripts/stage_E5/scan_secrets.py [extra folder or file ...]
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
E5 = ["scripts/stage_E5", "docs", "tests", "scripts/check_evidence_manuscript_numbers.py",
      "scripts/reproduce_evidence.py", "README.md", "phase1/config", "phase1/src",
      "phase1/reports/evidence/layout"]
SKIP_SUFFIXES = {".db", ".fa", ".fai", ".gtf", ".bgz", ".png", ".pdf", ".xlsx", ".pyc"}


def _files(paths: list[Path]) -> list[Path]:
    out = []
    for p in paths:
        if p.is_dir():
            out += [q for q in sorted(p.rglob("*")) if q.is_file()]
        elif p.is_file():
            out.append(p)
    return [p for p in dict.fromkeys(out)
            if p.suffix not in SKIP_SUFFIXES and "__pycache__" not in p.parts]


def main() -> int:
    spec = importlib.util.spec_from_file_location("no_secrets", REPO / "tests/test_no_secrets.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    paths = [REPO / f for f in E5] + [Path(a).expanduser() for a in sys.argv[1:]]
    files = _files(paths)
    hits = []
    for p in files:
        text = p.read_bytes().decode("utf-8", errors="replace")
        hits += [(p, label) for label, pat in m.PATTERNS.items() if pat.search(text)]
    for p, label in hits:
        print(f"MATCH {label}: {p}")
    shown = [str(f.relative_to(REPO)) if f.is_relative_to(REPO) else str(f) for f in paths if f.exists()]
    print(f"scanned {len(files)} files in {', '.join(shown)}: {len(hits)} match(es)")
    return 1 if hits else 0


if __name__ == "__main__":
    raise SystemExit(main())
