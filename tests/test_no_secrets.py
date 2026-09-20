"""No credential may enter the repository, and therefore its history.

This runs over `git ls-files`, not over the working tree: what matters is what is
tracked, because a key that reaches a commit stays reachable after the file is
deleted. It is a gate, not a scan of the disk.

The AlphaGenome key this project needs lives outside the tree -- in
$ALPHAGENOME_API_KEY, or in ~/.config/alphagenome/key, which is where
src/evid_score_avi.py looks.
"""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]

# Patterns for credentials that would be recoverable from a clone. Each is anchored
# on a provider's own prefix, so an ordinary word cannot trip it.
PATTERNS = {
    "Google API key (AlphaGenome, Ensembl, …)": re.compile(r"AIza[0-9A-Za-z_\-]{30,}"),
    "AWS access key id": re.compile(r"AKIA[0-9A-Z]{16}"),
    "GitHub token": re.compile(r"gh[pousr]_[0-9A-Za-z]{30,}"),
    "OpenAI key": re.compile(r"sk-[A-Za-z0-9]{40,}"),
    "private key block": re.compile(r"-----BEGIN (RSA |EC |OPENSSH |PGP )?PRIVATE KEY"),
}

# Binary and large data files are skipped: they cannot carry a pasted key and
# reading them all makes the gate slow enough that people stop running it.
SKIP_SUFFIXES = {".parquet", ".png", ".pdf", ".docx", ".xlsx", ".gz", ".zip",
                 ".fa", ".fai", ".db", ".jpg", ".jpeg", ".ico", ".woff", ".woff2"}
MAX_BYTES = 4_000_000


def tracked_files() -> list[Path]:
    out = subprocess.run(["git", "-C", str(REPO), "ls-files", "-z"],
                         capture_output=True, text=True, check=True).stdout
    return [REPO / p for p in out.split("\0") if p]


def test_no_credentials_in_tracked_files():
    hits = []
    for path in tracked_files():
        if path.suffix.lower() in SKIP_SUFFIXES or not path.is_file():
            continue
        # this test's own patterns are literals describing what to look for
        if path.resolve() == Path(__file__).resolve():
            continue
        if path.stat().st_size > MAX_BYTES:
            continue
        try:
            text = path.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        for label, pattern in PATTERNS.items():
            for m in pattern.finditer(text):
                line = text[:m.start()].count("\n") + 1
                hits.append(f"{path.relative_to(REPO)}:{line}  {label}")
    assert not hits, "credential-shaped strings in tracked files:\n" + "\n".join(hits)


def test_secret_bearing_paths_are_ignored():
    """The places a key is conventionally left must not be committable."""
    ignored = subprocess.run(
        ["git", "-C", str(REPO), "check-ignore", "-v",
         ".env", "secrets.env", "foo.key", "phase1/bar.key"],
        capture_output=True, text=True)
    covered = {line.split("\t")[-1] for line in ignored.stdout.splitlines()}
    for probe in (".env", "secrets.env", "foo.key", "phase1/bar.key"):
        assert probe in covered, f"{probe} is not gitignored"


def test_the_key_file_this_project_uses_is_outside_the_repo():
    from_home = Path.home() / ".config" / "alphagenome" / "key"
    assert REPO not in from_home.parents, "the key path resolves inside the repo"
    if from_home.exists():
        assert not any(
            from_home.resolve() == p.resolve() for p in tracked_files()), \
            "the key file is tracked"


@pytest.mark.parametrize("sample", [
    "AIzaSyA" + "b" * 32,
    "AKIA" + "A" * 16,
    "ghp_" + "x" * 36,
])
def test_the_gate_would_catch_a_real_key(sample):
    """A gate that cannot fail is not a gate."""
    assert any(p.search(sample) for p in PATTERNS.values()), sample
