"""B6.2: the figures' typeface, Liberation Sans, from its official release, into
phase1/assets/fonts/ with the licence it is released under.

Asks GitHub's API for the latest release of liberationfonts/liberation-fonts, the
project's own repository, and takes the TrueType archive its release notes link to
(liberation-fonts-ttf-<version>.tar.gz). From that archive, unchanged: the four
Liberation Sans faces (regular, bold, italic, bold italic) and the licence file, the
SIL Open Font License 1.1 in full. phase1/assets/fonts/SOURCE.json records the release,
the archive's URL, size and sha256, each file's place in the archive and sha256, and
the version each font names itself. No system font is read or copied.

Run from the repository root:
    .venv/bin/python scripts/stage_B6/fetch_fonts.py
"""
from __future__ import annotations

import hashlib
import io
import json
import re
import sys
import tarfile
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel, sha256_of, utc_now  # noqa: E402

PROJECT = "liberationfonts/liberation-fonts"
RELEASE_API = f"https://api.github.com/repos/{PROJECT}/releases/latest"
OUT = REPO / "phase1" / "assets" / "fonts"
FACES = ["LiberationSans-Regular.ttf", "LiberationSans-Bold.ttf",
         "LiberationSans-Italic.ttf", "LiberationSans-BoldItalic.ttf"]
LICENCE = "LICENSE"


def _get(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "variant-fm-benchmark B6",
                                               "Accept": "application/vnd.github+json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.read()


def font_version(path: Path) -> str:
    from fontTools.ttLib import TTFont
    names = TTFont(path)["name"]
    return str(names.getDebugName(5))            # name ID 5: the version string


def main() -> None:
    with RunLog(__file__) as log:
        release = json.loads(_get(RELEASE_API))
        version = release["tag_name"]
        links = re.findall(r"\((https://github\.com/" + re.escape(PROJECT)
                           + r"/files/\d+/liberation-fonts-ttf-" + re.escape(version)
                           + r"\.tar\.gz)\)", release.get("body") or "")
        if len(links) != 1:
            raise SystemExit(f"release {version}: {len(links)} TrueType archive links in "
                             "its notes")
        url = links[0]
        blob = _get(url)
        archive_sha = hashlib.sha256(blob).hexdigest()
        print(f"{PROJECT} {version}: {url}\n  {len(blob):,} bytes, sha256 {archive_sha}")
        OUT.mkdir(parents=True, exist_ok=True)
        files = {}
        with tarfile.open(fileobj=io.BytesIO(blob), mode="r:gz") as t:
            members = {Path(m.name).name: m for m in t.getmembers() if m.isfile()}
            for name in FACES + [LICENCE]:
                if name not in members:
                    raise SystemExit(f"{name} is not in the archive")
                data = t.extractfile(members[name]).read()
                (OUT / name).write_bytes(data)
                log.output(OUT / name)
                files[name] = {"member": members[name].name, "bytes": len(data),
                               "sha256": sha256_of(OUT / name)}
        for name in FACES:
            files[name]["version_string"] = font_version(OUT / name)
        licence = (OUT / LICENCE).read_text(encoding="utf-8", errors="replace")
        if "SIL OPEN FONT LICENSE VERSION 1.1" not in " ".join(licence.upper().split()):
            raise SystemExit("the archive's licence file is not the SIL OFL 1.1")
        source = {
            "family": "Liberation Sans",
            "project": f"https://github.com/{PROJECT}",
            "release": release["html_url"], "version": version,
            "published": release.get("published_at"),
            "archive": {"url": url, "bytes": len(blob), "sha256": archive_sha},
            "licence": {"file": LICENCE, "name": "SIL Open Font License 1.1",
                        "note": "the archive's licence file, unchanged"},
            "files": files,
            "fetched_utc": utc_now(),
            "written_by": rel(Path(__file__).resolve()),
        }
        (OUT / "SOURCE.json").write_text(json.dumps(source, indent=1, ensure_ascii=False)
                                         + "\n")
        log.output(OUT / "SOURCE.json")
        log.note("source", source)
        for name, f in files.items():
            print(f"  {name:32s} {f['sha256']}  {f.get('version_string', '')}")


if __name__ == "__main__":
    main()
