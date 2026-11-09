"""E8 (E19.0, E20.0, E17.0 rule 3): tests/test_no_secrets.py's patterns, unchanged, over the
files this stage writes, by scripts/stage_E5/scan_secrets.py's method (a copy of
scripts/stage_E7/scan_secrets.py): this stage's scripts and run logs (scripts/stage_E8),
the change log, the layout stage's code, configuration and tests, the variant list's
stage and test, and the layout's outputs, the workbook's members unzipped and read as
text. Prints every match (the pattern's name and the file, never the matched text) and
the count of files scanned; exits 1 on a match. Extra folders or files (the report, made
outside the repository) can be named on the command line.

    .venv/bin/python scripts/stage_E8/scan_secrets.py [extra folder or file ...]
"""
from __future__ import annotations

import importlib.util
import sys
import zipfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
E8 = ["scripts/stage_E8", "docs/evidence-change-log.md", "phase1/config/layout.yaml",
      "phase1/src/evid_layout.py", "phase1/tests/test_evid_layout.py",
      "phase1/src/evid_variant_list.py", "phase1/tests/test_evid_variant_list.py",
      "phase1/reports/evidence/layout"]
SKIP_SUFFIXES = {".db", ".fa", ".fai", ".gtf", ".bgz", ".png", ".pdf", ".pyc"}


def _files(paths: list[Path]) -> list[Path]:
    out = []
    for p in paths:
        if p.is_dir():
            out += [q for q in sorted(p.rglob("*")) if q.is_file()]
        elif p.is_file():
            out.append(p)
    return [p for p in dict.fromkeys(out)
            if p.suffix not in SKIP_SUFFIXES and "__pycache__" not in p.parts]


def _texts(p: Path) -> list[tuple[str, str]]:
    """A file's text, or each member's of a workbook: (name, text)."""
    if p.suffix == ".xlsx":
        with zipfile.ZipFile(p) as z:
            return [(f"{p}:{n}", z.read(n).decode("utf-8", errors="replace")) for n in z.namelist()]
    return [(str(p), p.read_bytes().decode("utf-8", errors="replace"))]


def main() -> int:
    spec = importlib.util.spec_from_file_location("no_secrets", REPO / "tests/test_no_secrets.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    paths = [REPO / f for f in E8] + [Path(a).expanduser() for a in sys.argv[1:]]
    files = _files(paths)
    hits, scanned = [], 0
    for p in files:
        for name, text in _texts(p):
            scanned += 1
            hits += [(name, label) for label, pat in m.PATTERNS.items() if pat.search(text)]
    for name, label in hits:
        print(f"MATCH {label}: {name}")
    shown = [str(f.relative_to(REPO)) if f.is_relative_to(REPO) else str(f) for f in paths if f.exists()]
    print(f"scanned {scanned} texts in {len(files)} files in {', '.join(shown)}: {len(hits)} match(es)")
    return 1 if hits else 0


if __name__ == "__main__":
    raise SystemExit(main())
