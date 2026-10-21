"""B6.1: the exact versions the reference image installs, read in the pinned base image.

Runs inside the base image at the digest container/ENVIRONMENT.json records
(base_image.reference), as root, with this repository mounted:

    docker run --rm --platform linux/amd64 -v "$PWD:/repo" -w /repo <reference> \\
        python scripts/stage_B6/pin_versions.py

  * apt: the candidate version of each Debian package container/Dockerfile installs
    (APT, with what it is for), from `apt-get update` and `apt-cache policy`, written
    to the apt_pins section of container/ENVIRONMENT.json; the Dockerfile's pins are
    taken from there, and a test holds the two equal.
  * pip: every package of requirements-evidence.lock.txt, at its version, downloaded
    as the wheel pip picks for this interpreter and platform (--only-binary=:all:
    --no-deps); each wheel's sha256 is computed here from the file, not copied from an
    index, and written with --hash to container/requirements.linux-amd64.txt, in the
    order of the lock. A package with no wheel for this platform stops the script:
    no compiler goes into the image (B6.1).
"""
from __future__ import annotations

import hashlib
import re
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel  # noqa: E402
import _environment  # noqa: E402

LOCK = REPO / "requirements-evidence.lock.txt"
OUT = REPO / "container" / "requirements.linux-amd64.txt"
APT = {"git": "scripts/reproduce_evidence.py --verify reads the committed outputs "
              "with git (rev-parse, ls-files, cat-file)"}


def _run(*cmd: str) -> str:
    return subprocess.run(cmd, check=True, capture_output=True, text=True).stdout


def lock_lines() -> list[str]:
    return [line.strip() for line in LOCK.read_text().splitlines()
            if line.strip() and not line.startswith("#")]


def apt_candidates() -> dict:
    _run("apt-get", "update")
    out = {}
    for pkg, use in APT.items():
        policy = _run("apt-cache", "policy", pkg)
        m = re.search(r"Candidate:\s*(\S+)", policy)
        if not m or m.group(1) == "(none)":
            raise SystemExit(f"apt has no candidate for {pkg}:\n{policy}")
        out[pkg] = {"version": m.group(1), "use": use, "apt_cache_policy": policy.strip()}
    return out


def wheel(spec: str, log) -> dict:
    with tempfile.TemporaryDirectory() as d:
        r = subprocess.run([sys.executable, "-m", "pip", "download", "--no-deps",
                            "--only-binary=:all:", "--no-cache-dir",
                            "--disable-pip-version-check", "--dest", d, spec],
                           capture_output=True, text=True)
        if r.returncode:
            raise SystemExit(f"{spec}: no wheel for this platform (B6.1 stop)\n"
                             f"{r.stdout[-2000:]}{r.stderr[-2000:]}")
        files = list(Path(d).iterdir())
        if len(files) != 1 or files[0].suffix != ".whl":
            raise SystemExit(f"{spec}: pip download gave {[f.name for f in files]}")
        f = files[0]
        return {"spec": spec, "file": f.name, "bytes": f.stat().st_size,
                "sha256": hashlib.sha256(f.read_bytes()).hexdigest()}


def main() -> None:
    with RunLog(__file__) as log:
        log.input(LOCK)
        env = _environment.read().get("base_image", {})
        log.note("python", _run(sys.executable, "-VV").strip())
        log.note("pip", _run(sys.executable, "-m", "pip", "--version").strip())
        log.note("uname_m", _run("uname", "-m").strip())
        apt = apt_candidates()
        wheels = [wheel(spec, log) for spec in lock_lines()]
        for w in wheels:
            print(f"  {w['file']}  {w['sha256']}")
        header = [
            "# The packages of requirements-evidence.lock.txt, each at the same version, for the",
            "# reference image (container/Dockerfile, linux/amd64). Written by",
            "# scripts/stage_B6/pin_versions.py in the pinned base image",
            f"# {env.get('reference', '(base_image.reference in container/ENVIRONMENT.json)')},",
            "# which downloaded each wheel and hashed the file; do not edit by hand.",
            "# Installed with pip install --require-hashes --no-deps --only-binary=:all:.",
        ]
        OUT.write_text("\n".join(header + [f"{w['spec']} --hash=sha256:{w['sha256']}"
                                           for w in wheels]) + "\n")
        log.output(OUT)
        _environment.update("apt_pins", {
            "base_image": env.get("reference"),
            "packages": {k: {"version": v["version"], "use": v["use"]} for k, v in apt.items()},
        }, __file__, log)
        log.note("apt", apt)
        log.note("wheels", wheels)
        print(f"wrote {rel(OUT)} ({len(wheels)} packages); apt: "
              + ", ".join(f"{k}={v['version']}" for k, v in apt.items()))


if __name__ == "__main__":
    main()
