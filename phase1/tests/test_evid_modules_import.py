"""Every stage module must parse and import.

This exists because a history rewrite that renamed paths and text across the branch
also rewrote an identifier inside one module, turning EVID_DIR into a name with a
hyphen in it. The file stopped parsing, and the whole suite still passed, because no
test imported that module -- the tests read its OUTPUTS, which were already on disk.

A suite that only reads outputs cannot tell a working pipeline from a broken one
whose last run happened to succeed. This closes that gap: parse every tracked Python
file, and import every stage module.
"""
from __future__ import annotations

import ast
import importlib
import re
import subprocess
import sys
from pathlib import Path

import pytest

PHASE1 = Path(__file__).resolve().parents[1]
REPO = PHASE1.parent


def tracked(pattern: str) -> list[Path]:
    out = subprocess.run(["git", "-C", str(REPO), "ls-files", pattern],
                         capture_output=True, text=True, check=True).stdout
    return [REPO / p for p in out.split("\n") if p.strip()]


@pytest.mark.parametrize("path", tracked("*.py"), ids=lambda p: str(p.name))
def test_every_tracked_python_file_parses(path):
    ast.parse(path.read_text(encoding="utf-8", errors="ignore"), filename=str(path))


STAGE_MODULES = sorted(p.stem for p in (PHASE1 / "src").glob("evid_*.py"))


@pytest.mark.parametrize("name", STAGE_MODULES)
def test_every_stage_module_imports(name):
    """Importing runs the module-level code -- the paths, the constants, the
    configuration reads -- which is where the rewrite damage was."""
    sys.path.insert(0, str(PHASE1))
    try:
        mod = importlib.import_module(f"src.{name}")
    except ImportError as e:
        # a scorer that needs a model environment the analysis venv does not carry
        if any(k in str(e) for k in ("spliceai", "pangolin", "alphagenome", "gffutils")):
            pytest.skip(f"{name} needs a model environment: {e}")
        raise
    assert mod is not None


def test_no_identifier_carries_a_hyphen_from_the_rename():
    """The rewrite's substitution was textual, so it could land inside a name.
    A hyphen cannot appear in a Python identifier or in a shell variable name,
    which makes it the signature of exactly that damage."""
    bad = []
    pat = re.compile(r"[A-Za-z_][A-Za-z0-9_]*-[A-Za-z0-9_]*(?:_DIR|_REPO|_PATH|_FILE)\b")
    for path in tracked("*.py") + tracked("*.md") + tracked("*.sh"):
        if path.name == Path(__file__).name or not path.is_file():
            continue
        for i, line in enumerate(path.read_text(encoding="utf-8", errors="ignore")
                                 .splitlines(), 1):
            if pat.search(line):
                bad.append(f"{path.relative_to(REPO)}:{i}: {line.strip()[:90]}")
    assert not bad, "identifier damaged by a textual rename:\n" + "\n".join(bad)
