"""E20b (E19.0, E20.0, E17.0 rule 3): tests/test_no_secrets.py's patterns, unchanged, over the
files this stage writes, by scripts/stage_E18/scan_secrets.py's method: E18's folders
and E19's, E20's and this stage's scripts and run logs (scripts/stage_E20, scripts/stage_E20b), the VCEP records
(phase1/data/evidence/vcep), the true-positive scorer's work folder, and the
true-positive subsets, score inputs and event records written beside E6's in
phase1/data/evidence (the files named *_tp*; the external genes' are in
phase1/data/evidence/external, one of E18's folders). Files in the gitignored caches are
downloaded data, read by no code that reads a key, and are skipped as E18's scanner
skips them. Prints every match (the pattern's name, never the matched text) and the
count of files scanned; exits 1 on a match.

    .venv/bin/python scripts/stage_E20b/scan_secrets.py [extra folder ...]
"""
from __future__ import annotations

import gzip
import importlib.util
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
E18 = ["scripts/stage_E18", "phase1/data/evidence/external", "phase1/data/evidence/alphagenome_official",
       "phase1/data/e18_cache", "phase1/data/alphagenome_official_cache", "docs"]
E19 = ["scripts/stage_E19", "phase1/reports/evidence/suppl_analyses", "phase1/data/evidence/e19",
       "phase1/data/evidence/second_set"]
E20 = ["scripts/stage_E20", "scripts/stage_E20b", "phase1/data/evidence/vcep",
       "phase1/data/evidence/pangolin_events_work_tp", "phase1/data/evidence/pangolin_events_work"]
E20_FILES = ["phase1/data/evidence/*_tp.parquet", "phase1/data/evidence/*_tp.parquet.*.json",
             "phase1/data/evidence/inframe_*.parquet", "phase1/data/evidence/inframe_*.parquet.*.json"]
SKIP_SUFFIXES = {".db", ".fa", ".fai", ".gtf", ".bgz"}


def _files(folders: list[Path]) -> list[Path]:
    out = []
    for folder in folders:
        if folder.is_dir():
            out += [p for p in sorted(folder.rglob("*")) if p.is_file()]
    for pattern in E20_FILES:
        out += sorted(REPO.glob(pattern))
    return [p for p in dict.fromkeys(out)
            if p.suffix not in SKIP_SUFFIXES and "__pycache__" not in p.parts]


def main() -> int:
    spec = importlib.util.spec_from_file_location("no_secrets", REPO / "tests/test_no_secrets.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    folders = [REPO / f for f in E18 + E19 + E20] + [Path(a).expanduser() for a in sys.argv[1:]]
    files = _files(folders)
    hits = []
    for p in files:
        data = p.read_bytes()
        if p.suffix == ".gz":
            try:
                data = gzip.decompress(data)
            except OSError:
                pass
        text = data.decode("utf-8", errors="replace")
        hits += [(p, label) for label, pat in m.PATTERNS.items() if pat.search(text)]
    for p, label in hits:
        print(f"MATCH {label}: {p}")
    shown = [str(f.relative_to(REPO)) if f.is_relative_to(REPO) else str(f)
             for f in folders if f.exists()]
    print(f"scanned {len(files)} files in {', '.join(shown)} and {', '.join(E20_FILES)}: "
          f"{len(hits)} match(es)")
    return 1 if hits else 0


if __name__ == "__main__":
    raise SystemExit(main())
