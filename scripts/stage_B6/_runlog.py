"""Run records for the B6 one-off scripts, in the format of scripts/stage_E16/_runlog.py
(a copy, so that the run records of the earlier stages keep the sha256 they carry).

Every script in this folder writes run_log/<script>_<UTC timestamp>.json when it
finishes, whether it succeeded or not: its argv, the commit it ran at, the sha256 of
the script and of this module as they ran, the interpreter, the platform and the
versions of the packages that touch the numbers, each input and output with its size
and sha256, start and end times, the exception if there was one, and any notes the
script adds; an input is hashed when it is read, an output when the script ends. Paths
are recorded relative to the repository, or to the home folder (~/) when they lie
outside it. Some of these scripts run inside the reference image (container/), where
the repository is mounted and git may be absent: the commit is then "unavailable".
"""
from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
import traceback
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
LOG_DIR = Path(__file__).resolve().parent / "run_log"
PACKAGES = ["numpy", "scipy", "pandas", "matplotlib", "threadpoolctl"]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def rel(path: Path | str) -> str:
    p = Path(path).resolve()
    for base, prefix in ((REPO, ""), (Path.home().resolve(), "~/")):
        try:
            return prefix + str(p.relative_to(base))
        except ValueError:
            pass
    return str(p)


def _version(pkg: str) -> str:
    try:
        return metadata.version(pkg)
    except metadata.PackageNotFoundError:
        return "not installed"


def _git(*args: str) -> str:
    try:
        return subprocess.run(["git", "-C", str(REPO), *args], check=True,
                              capture_output=True, text=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def _file_entry(path: Path) -> dict:
    path = Path(path)
    if not path.is_file():
        return {"path": rel(path), "exists": False}
    return {"path": rel(path), "bytes": path.stat().st_size, "sha256": sha256_of(path)}


class RunLog:
    """with RunLog(__file__) as log: ... log.input(p); log.output(p); log.note(k, v)"""

    def __init__(self, script: str):
        self.script = Path(script).stem
        self.script_path = rel(Path(script).resolve())
        self.code = {rel(p): sha256_of(p) for p in (Path(script).resolve(),
                                                     Path(__file__).resolve())}
        self.inputs: list[dict] = []
        self.outputs: list[Path] = []
        self.notes: dict = {}
        self.start = utc_now()
        self.path: Path | None = None

    def input(self, path: Path | str) -> None:
        """Hashed now, as read: a script may rewrite a file it reads (a manifest)."""
        self.inputs.append(_file_entry(Path(path)))

    def output(self, path: Path | str) -> None:
        self.outputs.append(Path(path))

    def note(self, key: str, value) -> None:
        self.notes[key] = value

    def __enter__(self) -> "RunLog":
        return self

    def __exit__(self, exc_type, exc, tb) -> bool:
        end = utc_now()
        record = {
            "script": self.script_path,
            "argv": sys.argv,
            "git_commit": _git("rev-parse", "HEAD"),
            "git_worktree_clean": _git("status", "--porcelain") == "",
            "code_sha256": self.code,
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "packages": {p: _version(p) for p in PACKAGES},
            "start_utc": self.start,
            "end_utc": end,
            "inputs": self.inputs,
            "outputs": [_file_entry(p) for p in self.outputs],
            "exception": ("".join(traceback.format_exception(exc_type, exc, tb))
                          if exc_type is not None else None),
            "notes": self.notes,
        }
        LOG_DIR.mkdir(parents=True, exist_ok=True)
        stamp = end.replace(":", "").replace("-", "").replace("+0000", "Z")
        self.path = LOG_DIR / f"{self.script}_{stamp}.json"
        self.path.write_text(json.dumps(record, indent=1, ensure_ascii=False) + "\n")
        print(f"[run log] {rel(self.path)}")
        return False
