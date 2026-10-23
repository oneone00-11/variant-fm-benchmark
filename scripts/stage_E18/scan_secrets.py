"""E18.0 rule 6 (E17.0 rule 3): tests/test_no_secrets.py's patterns, unchanged, over the
files this stage writes, tracked or not: the run logs, the score files and their
records, the caches and the NT package. Prints every file scanned and every match (the
pattern's name, never the matched text); exits 1 on a match.

    .venv/bin/python scripts/stage_E18/scan_secrets.py [extra folder ...]
"""
from __future__ import annotations

import gzip
import importlib.util
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
FOLDERS = ["scripts/stage_E18", "phase1/data/evidence/external", "phase1/data/evidence/alphagenome_official",
           "phase1/data/e18_cache", "phase1/data/alphagenome_official_cache", "docs"]
# Pangolin's gene databases, reference sequences, and GPN-MSA's sparse local copy of the
# 37 GB scores file (holes but for the ranges fetched): downloaded data, never written by
# anything that reads the key
SKIP_SUFFIXES = {".db", ".fa", ".fai", ".gtf", ".bgz"}


def main() -> int:
    spec = importlib.util.spec_from_file_location("no_secrets", REPO / "tests/test_no_secrets.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    folders = [REPO / f for f in FOLDERS] + [Path(a).expanduser() for a in sys.argv[1:]]
    n, hits = 0, []
    for folder in folders:
        if not folder.exists():
            continue
        for p in sorted(folder.rglob("*")):
            if not p.is_file() or p.suffix in SKIP_SUFFIXES or "__pycache__" in p.parts:
                continue
            data = p.read_bytes()
            if p.suffix == ".gz":
                try:
                    data = gzip.decompress(data)
                except OSError:
                    pass
            text = data.decode("utf-8", errors="replace")
            n += 1
            hits += [(p, label) for label, pat in m.PATTERNS.items() if pat.search(text)]
    for p, label in hits:
        print(f"MATCH {label}: {p}")
    print(f"scanned {n} files in {', '.join(str(f) for f in folders if f.exists())}: "
          f"{len(hits)} match(es)")
    return 1 if hits else 0


if __name__ == "__main__":
    raise SystemExit(main())
