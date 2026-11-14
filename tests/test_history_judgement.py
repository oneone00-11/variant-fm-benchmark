"""E9: the one judgement by which the tests that read the repository's history run or
skip (tests/_checkout.missing_history): run in a full clone, skip in a shallow clone, skip
in a release archive without .git and never an error there. Each case is made here in a
small repository of its own, so that the result does not depend on how this checkout was
made; and every test that reads the history is held to that judgement."""
from __future__ import annotations

import ast
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "tests"))
import _checkout as K  # noqa: E402

# the tests that read the history, each by the one judgement (E9; the author's decision
# of 2026-11-14): its file and its function
HISTORY_TESTS = [
    ("tests/test_container.py", "test_the_rebaselined_pngs_differ_in_their_bytes_only"),
    ("tests/test_evidence_manuscript_numbers.py", "test_the_v5_table_reproduces_e5_0_on_the_outputs_it_read"),
    ("tests/test_evidence_manuscript_numbers.py", "test_the_v5_texts_read_number_by_number"),
    ("phase1/tests/test_evid_layout.py", "test_every_output_added_since_release_301_has_a_place"),
    ("phase1/tests/test_evid_newgenes.py", "test_without_the_argument_the_stage_writes_e16s_files_as_they_were"),
    ("phase1/tests/test_label_build.py", "test_palb2_and_bard1_differ_from_the_july_files_exactly_where_the_audit_said"),
    ("phase1/tests/test_label_build.py", "test_the_switch_commit_changes_no_code"),
    ("phase1/tests/test_label_build.py", "test_the_superseded_files_are_the_july_files"),
    ("phase1/tests/test_evid_label_audit.py", "test_the_label_files_entered_the_repository_on_the_date_the_audit_uses"),
]


def _git(cwd: Path, *args: str) -> str:
    return subprocess.run(["git", "-c", "user.name=test", "-c", "user.email=test@example.invalid",
                           "-c", "init.defaultBranch=main", "-C", str(cwd), *args],
                          capture_output=True, text=True, check=True).stdout.strip()


@pytest.fixture(scope="module")
def full(tmp_path_factory) -> tuple[Path, str, str]:
    """A repository of two commits: its folder, and the first commit and the second."""
    root = tmp_path_factory.mktemp("full")
    _git(root, "init", "-q")
    for k in (1, 2):
        (root / f"file{k}.txt").write_text(f"{k}\n")
        _git(root, "add", f"file{k}.txt")
        _git(root, "commit", "-q", "-m", f"commit {k}")
    first, second = _git(root, "rev-list", "--reverse", "HEAD").split()
    return root, first, second


def test_a_full_clone_has_its_history(full):
    root, first, second = full
    assert K.missing_history(repo=root) is None
    assert K.missing_history(first, second, repo=root) is None
    assert K.missing_history("0" * 40, repo=root) == "commit 0000000 is not in this clone"


def test_a_shallow_clone_has_not(full, tmp_path):
    root, first, second = full
    shallow = tmp_path / "shallow"
    subprocess.run(["git", "clone", "-q", "--depth", "1", f"file://{root}", str(shallow)],
                   capture_output=True, check=True)
    assert (shallow / ".git").exists() and (shallow / "file1.txt").exists()
    assert K.missing_history(repo=shallow) == "a shallow clone holds the last commits only"
    assert K.missing_history(second, repo=shallow) == "a shallow clone holds the last commits only"
    # what the judgement guards against: there git takes the one commit it holds for the
    # one that added every file
    added = _git(shallow, "log", "--diff-filter=A", "--format=%H", "--", "file1.txt").split()
    assert added == [second]


def test_a_release_archive_has_no_git_and_git_is_not_asked(tmp_path, monkeypatch):
    (tmp_path / "file1.txt").write_text("1\n")

    def no_git(*a, **k):
        raise AssertionError("git was asked about a tree without .git")
    monkeypatch.setattr(K.subprocess, "run", no_git)
    assert K.missing_history(repo=tmp_path) == "a release archive has no .git"
    assert K.missing_history("abc1234", repo=tmp_path) == "a release archive has no .git"


def test_a_git_that_cannot_be_run_is_no_history_and_no_error(full, monkeypatch):
    root, first, _ = full

    def missing(*a, **k):
        raise FileNotFoundError("git")
    monkeypatch.setattr(K.subprocess, "run", missing)
    assert K.missing_history(first, repo=root).startswith("git cannot be run here")


def test_needs_history_skips_with_the_reason(monkeypatch, tmp_path):
    monkeypatch.setitem(K.missing_history.__kwdefaults__, "repo", tmp_path)
    with pytest.raises(pytest.skip.Exception, match="a release archive has no .git"):
        K.needs_history()


@pytest.mark.parametrize("path, name", HISTORY_TESTS)
def test_every_test_that_reads_the_history_uses_the_one_judgement(path, name):
    """Its skips for want of history go through needs_history, and only through it."""
    tree = ast.parse((REPO / path).read_text())
    fn = next(n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == name)
    body = fn.body[1:] if ast.get_docstring(fn) is not None else fn.body
    src = "\n".join(ast.unparse(s) for s in body)
    assert "needs_history(" in src, (path, name)
    for old in ("is not in this clone", "has no .git", "cat-file", "rev-parse", "shallow"):
        assert old not in src, (path, name, old)
