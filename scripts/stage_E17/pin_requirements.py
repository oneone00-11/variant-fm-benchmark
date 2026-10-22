"""E17.2: the AlphaGenome scoring environment, pinned with hashes.

The official notebook installs the client unpinned (`pip install alphagenome`), so the
client is the newest release on PyPI on the day, read from PyPI's JSON API with its
upload time. This resolves that release's dependencies in a scratch environment of the
interpreter given, then downloads, for every package of the resolved set, the file pip
picks on this machine (a wheel where there is one) and hashes it, and writes
scripts/stage_E17/requirements.txt: `name==version --hash=sha256:...` for each, nothing
typed by hand. .venv-ag-official/ is built from it with `pip install --require-hashes
--no-deps`; it only scores, it never runs the pipeline.

Run from the repository root:
    python3.12 scripts/stage_E17/pin_requirements.py --python python3.12
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import tempfile
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel, utc_now  # noqa: E402

OUT = REPO / "scripts" / "stage_E17" / "requirements.txt"
PYPI = "https://pypi.org/pypi/alphagenome/json"


def _run(*cmd: str) -> str:
    return subprocess.run(cmd, check=True, capture_output=True, text=True).stdout


def latest() -> dict:
    with urllib.request.urlopen(PYPI, timeout=60) as r:
        d = json.loads(r.read())
    v = d["info"]["version"]
    return {"version": v, "queried_utc": utc_now(),
            "files": [{"filename": f["filename"], "upload_time": f["upload_time_iso_8601"],
                       "sha256": f["digests"]["sha256"]} for f in d["releases"][v]]}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--python", default="python3.12", help="the interpreter of the venv")
    args = ap.parse_args()
    with RunLog(__file__) as log:
        rel_ = latest()
        version = rel_["version"]
        log.note("pypi", rel_)
        with tempfile.TemporaryDirectory() as d:
            venv = Path(d) / "resolve"
            _run(args.python, "-m", "venv", str(venv))
            py = str(venv / "bin" / "python")
            log.note("interpreter", _run(py, "-VV").strip())
            log.note("pip", _run(py, "-m", "pip", "--version").strip())
            _run(py, "-m", "pip", "install", "--quiet", "--disable-pip-version-check",
                 f"alphagenome=={version}")
            freeze = sorted(_run(py, "-m", "pip", "freeze").split())
            files = []
            for spec in freeze:
                dest = Path(d) / "files" / spec
                dest.mkdir(parents=True)
                r = subprocess.run([py, "-m", "pip", "download", "--no-deps", "--quiet",
                                    "--disable-pip-version-check", "--prefer-binary",
                                    "--dest", str(dest), spec],
                                   capture_output=True, text=True)
                if r.returncode:
                    raise SystemExit(f"{spec}: pip download failed\n{r.stderr[-1500:]}")
                got = list(dest.iterdir())
                if len(got) != 1:
                    raise SystemExit(f"{spec}: downloaded {[g.name for g in got]}")
                files.append({"spec": spec, "file": got[0].name,
                              "kind": "wheel" if got[0].suffix == ".whl" else "sdist",
                              "sha256": hashlib.sha256(got[0].read_bytes()).hexdigest()})
        header = [
            "# The AlphaGenome scoring environment of stage E17 (not the pipeline's): the",
            f"# client alphagenome=={version}, the newest release on PyPI on "
            f"{rel_['queried_utc'][:10]} (the official",
            "# notebook installs it unpinned), and the packages it resolved to, each with the",
            "# sha256 of the file pip picks here. Written by scripts/stage_E17/pin_requirements.py;",
            "# do not edit by hand. Installed with pip install --require-hashes --no-deps.",
        ]
        OUT.write_text("\n".join(header + [f"{f['spec']} --hash=sha256:{f['sha256']}"
                                           for f in files]) + "\n")
        log.output(OUT)
        log.note("files", files)
        print(f"alphagenome=={version}: {len(files)} packages "
              f"({sum(f['kind'] == 'sdist' for f in files)} as sdist) -> {rel(OUT)}")
        for f in files:
            print(f"  {f['file']}")


if __name__ == "__main__":
    main()
