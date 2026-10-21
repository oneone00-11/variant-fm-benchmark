"""B6: the reference image's record of itself, written from inside it to the inside
section of container/ENVIRONMENT.json, and, once the fonts are in the repository
(phase1/assets/fonts/), to its fonts section.

  * inside: the interpreter (python -VV), pip, the installed packages (pip freeze),
    compared name==version with both lock files (requirements-evidence.lock.txt and
    container/requirements.linux-amd64.txt), pip check, every Debian package and its
    version (dpkg-query), the C library and the OS release;
  * fonts: each font file the figures load, with its sha256, as the files are here.

Run inside the image:  bash container/run.sh environment
"""
from __future__ import annotations

import os
import platform
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel, sha256_of  # noqa: E402
import _environment  # noqa: E402

LOCKS = [REPO / "requirements-evidence.lock.txt",
         REPO / "container" / "requirements.linux-amd64.txt"]
FONTS = REPO / "phase1" / "assets" / "fonts"


def _run(*cmd: str) -> str:
    return subprocess.run(cmd, check=True, capture_output=True, text=True).stdout


def lock_set(path: Path) -> set[str]:
    return {line.split()[0] for line in path.read_text().splitlines()
            if line.strip() and not line.startswith("#")}


def os_release() -> dict:
    out = {}
    for line in Path("/etc/os-release").read_text().splitlines():
        key, _, value = line.partition("=")
        if key in ("PRETTY_NAME", "ID", "VERSION_ID", "VERSION_CODENAME"):
            out[key] = value.strip('"')
    return out


def main() -> None:
    with RunLog(__file__) as log:
        freeze = sorted(_run(sys.executable, "-m", "pip", "freeze").split())
        locks = {}
        for path in LOCKS:
            log.input(path)
            want = lock_set(path)
            locks[rel(path)] = {"same_set": set(freeze) == want,
                                "only_installed": sorted(set(freeze) - want),
                                "only_in_lock": sorted(want - set(freeze))}
        dpkg = {}
        for line in _run("dpkg-query", "-W", "-f", "${Package}\t${Version}\n").splitlines():
            pkg, version = line.split("\t")
            dpkg[pkg] = version
        inside = {
            "python": _run(sys.executable, "-VV").strip(),
            "pip": _run(sys.executable, "-m", "pip", "--version").strip(),
            "pip_check": _run(sys.executable, "-m", "pip", "check").strip(),
            "pip_freeze": freeze,
            "pip_freeze_against_locks": locks,
            "dpkg_packages": dpkg,
            "libc": " ".join(platform.libc_ver()),
            "os_release": os_release(),
            "machine": os.uname().machine,
        }
        _environment.update("inside", inside, __file__, log)
        log.note("locks", locks)
        print(f"{inside['python']}\n{inside['pip']}\n{len(freeze)} packages; "
              + "; ".join(f"{k}: {'same set' if v['same_set'] else 'DIFFERS'}"
                          for k, v in locks.items())
              + f"\n{len(dpkg)} Debian packages; git {dpkg.get('git', 'absent')}")
        if FONTS.is_dir():
            files = sorted(p for p in FONTS.iterdir() if p.suffix.lower() == ".ttf")
            _environment.update("fonts", {"folder": rel(FONTS), "files": {
                p.name: sha256_of(p) for p in files}}, __file__, log)
            for p in files:
                log.input(p)
        if not all(v["same_set"] for v in locks.values()):
            raise SystemExit("the installed packages are not the locks' (see the record)")


if __name__ == "__main__":
    main()
