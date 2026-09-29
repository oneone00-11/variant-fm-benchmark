"""The files this checkout carries: `git ls-files` in a clone, and a walk of the tree in
a release archive, which has no .git. The walk leaves out what a run adds to a tree
(.venv, the fetched inputs, caches, the --verify report and baseline)."""
from __future__ import annotations

import fnmatch
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
_SKIP_DIRS = {".git", ".venv", "ag-env", "__pycache__", ".pytest_cache", ".reproduce_baseline"}
_SKIP_PATHS = ("phase1/data/evidence/companion_atlas/", "phase1/data/evidence/clinvar/")
_SKIP_FILES = {"reproduce_report.txt"}


def is_clone() -> bool:
    return (REPO / ".git").exists()


def checkout_files(pattern: str | None = None) -> list[Path]:
    """Every file the checkout carries, or those matching a git pathspec: a glob
    ("*.py", which matches at any depth) or a directory prefix."""
    if is_clone():
        args = ["git", "-C", str(REPO), "ls-files", "-z"] + ([pattern] if pattern else [])
        out = subprocess.run(args, capture_output=True, text=True, check=True).stdout
        return [REPO / p for p in out.split("\0") if p]
    files = []
    for p in REPO.rglob("*"):
        rel = p.relative_to(REPO).as_posix()
        if (not p.is_file() or set(p.relative_to(REPO).parts) & _SKIP_DIRS
                or rel.startswith(_SKIP_PATHS) or p.name in _SKIP_FILES):
            continue
        if pattern is None or (fnmatch.fnmatch(rel, pattern) if any(c in pattern for c in "*?[")
                               else rel == pattern or rel.startswith(pattern.rstrip("/") + "/")):
            files.append(p)
    return sorted(files)
