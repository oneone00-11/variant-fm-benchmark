"""The files this checkout carries: `git ls-files` in a clone, and a walk of the tree in
a release archive, which has no .git. The walk leaves out what a run adds to a tree
(.venv, the fetched inputs, caches, the --verify report and baseline)."""
from __future__ import annotations

import fnmatch
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
_SKIP_DIRS = {".git", ".venv", "ag-env", "__pycache__", ".pytest_cache", ".reproduce_baseline"}
_SKIP_PATHS = ("phase1/data/evidence/companion_atlas/", "phase1/data/evidence/clinvar/",
               "phase1/data/evidence/refs/", "phase1/data/evidence/label_audit_cache/")
_SKIP_FILES = {"reproduce_report.txt"}


def is_clone() -> bool:
    return (REPO / ".git").exists()


def missing_history(*commits: str, repo: Path = REPO) -> str | None:
    """Why a test that reads this repository's history cannot run here, or None when it
    can. One judgement for every such test (E9): a release archive has no .git, and git is
    not asked at all (an archive unpacked inside another clone would otherwise read that
    clone's history); a shallow clone holds the last commits only, and there git takes
    the oldest commit it holds for the one that added every file; and a full clone must
    hold each commit the test names. A git that cannot be run, or that answers anything
    but "false" to --is-shallow-repository, counts as no history: never an error."""
    if not (repo / ".git").exists():
        return "a release archive has no .git"
    try:
        shallow = subprocess.run(["git", "-C", str(repo), "rev-parse", "--is-shallow-repository"],
                                 capture_output=True, text=True)
        if shallow.returncode or shallow.stdout.strip() != "false":
            return "a shallow clone holds the last commits only"
        for c in commits:
            if subprocess.run(["git", "-C", str(repo), "cat-file", "-e", f"{c}^{{commit}}"],
                              capture_output=True).returncode:
                return f"commit {c[:7]} is not in this clone"
    except OSError as e:
        return f"git cannot be run here ({e})"
    return None


def needs_history(*commits: str) -> None:
    """Skip the calling test, giving missing_history()'s reason, where it gives one."""
    import pytest
    why = missing_history(*commits)
    if why:
        pytest.skip(why)


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
