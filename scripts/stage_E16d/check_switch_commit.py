"""E16d: check that a commit changes no code, configuration or workflow file (E16d.0
rule 2), so that every output it changes can only come from its data.

Code, configuration and workflow mean: a .py, .sh, .yml, .yaml, .toml, .cfg, .ini or
.lock file; config.py; .gitattributes or .gitignore; any .txt under scripts/ (the
entry point's input lists and pins); and anything under .github/. Prints the commit's
git diff --stat and exits non-zero if such a file is in it.

    .venv/bin/python scripts/stage_E16d/check_switch_commit.py <commit>
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
CODE = re.compile(r"(\.py|\.sh|\.ya?ml|\.toml|\.cfg|\.ini|\.lock)$")


def is_code(path: str) -> bool:
    name = Path(path).name
    return bool(CODE.search(path) or name in ("config.py", ".gitattributes", ".gitignore")
                or path.startswith(".github/")
                or (path.startswith("scripts/") and path.endswith(".txt")))


def changed(commit: str) -> list[str]:
    out = subprocess.run(["git", "-C", str(REPO), "diff", "--name-only", f"{commit}^", commit],
                         capture_output=True, text=True, check=True).stdout
    return [p for p in out.splitlines() if p]


def main() -> None:
    commit = sys.argv[1]
    stat = subprocess.run(["git", "-C", str(REPO), "diff", "--stat", f"{commit}^", commit],
                          capture_output=True, text=True, check=True).stdout
    print(stat)
    bad = [p for p in changed(commit) if is_code(p)]
    if bad:
        raise SystemExit(f"{commit} changes code, configuration or workflow files: {bad}")
    print(f"{commit}: no code, configuration or workflow file changed")


if __name__ == "__main__":
    main()
