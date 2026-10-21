"""B6.1: the base image of the reference environment, by digest, as the registry serves
it now.

Queries Docker Hub's registry for python:3.12.13-slim-bookworm, the official image of
the Python the committed outputs were made with (3.12.13), on Debian bookworm. Each
document is checked against the digest it is served under: the image index, the
linux/amd64 manifest it lists (the only platform the reference image is built for) and
that manifest's config, whose environment names the Python version. Writes the
base_image section of container/ENVIRONMENT.json. container/Dockerfile's FROM line
names the linux/amd64 manifest digest; a test holds the two equal.

Run from the repository root:
    python3 scripts/stage_B6/base_image.py
"""
from __future__ import annotations

import hashlib
import json
import sys
import urllib.error
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import RunLog  # noqa: E402
import _environment  # noqa: E402

REPOSITORY = "library/python"
TAG = "3.12.13-slim-bookworm"
PLATFORM = {"os": "linux", "architecture": "amd64"}
REGISTRY = "https://registry-1.docker.io/v2"
TOKEN_URL = ("https://auth.docker.io/token?service=registry.docker.io"
             f"&scope=repository:{REPOSITORY}:pull")
INDEX_TYPES = ["application/vnd.oci.image.index.v1+json",
               "application/vnd.docker.distribution.manifest.list.v2+json"]
MANIFEST_TYPES = ["application/vnd.oci.image.manifest.v1+json",
                  "application/vnd.docker.distribution.manifest.v2+json"]


class _NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, *args, **kwargs):
        return None


def _get(url: str, token: str, accept: list[str] | None = None) -> tuple[bytes, dict]:
    headers = {"Authorization": f"Bearer {token}"}
    if accept:
        headers["Accept"] = ", ".join(accept)
    opener = urllib.request.build_opener(_NoRedirect)
    try:
        with opener.open(urllib.request.Request(url, headers=headers), timeout=60) as r:
            return r.read(), dict(r.headers)
    except urllib.error.HTTPError as e:
        if e.code not in (301, 302, 303, 307, 308):
            raise
        # a blob is served from a CDN, which takes no registry token
        with urllib.request.urlopen(e.headers["Location"], timeout=60) as r:
            return r.read(), dict(r.headers)


def _checked(body: bytes, digest: str, what: str) -> None:
    got = "sha256:" + hashlib.sha256(body).hexdigest()
    if got != digest:
        raise SystemExit(f"{what}: served as {digest}, but its bytes hash to {got}")


def main() -> None:
    with RunLog(__file__) as log:
        token = json.loads(urllib.request.urlopen(TOKEN_URL, timeout=60).read())["token"]
        body, headers = _get(f"{REGISTRY}/{REPOSITORY}/manifests/{TAG}", token, INDEX_TYPES)
        index_digest = headers.get("Docker-Content-Digest") or headers.get("docker-content-digest")
        _checked(body, index_digest, f"the index of python:{TAG}")
        index = json.loads(body)
        if index.get("mediaType") not in INDEX_TYPES:
            raise SystemExit(f"python:{TAG} is not served as an image index "
                             f"({index.get('mediaType')})")
        hits = [m for m in index["manifests"]
                if {k: m.get("platform", {}).get(k) for k in PLATFORM} == PLATFORM
                and not m["platform"].get("variant")]
        if len(hits) != 1:
            raise SystemExit(f"python:{TAG} lists {len(hits)} linux/amd64 manifests")
        amd64 = hits[0]["digest"]
        body, headers = _get(f"{REGISTRY}/{REPOSITORY}/manifests/{amd64}", token,
                             MANIFEST_TYPES)
        _checked(body, amd64, "the linux/amd64 manifest")
        manifest = json.loads(body)
        config_digest = manifest["config"]["digest"]
        body, _ = _get(f"{REGISTRY}/{REPOSITORY}/blobs/{config_digest}", token)
        _checked(body, config_digest, "the linux/amd64 config")
        config = json.loads(body)
        env = dict(e.split("=", 1) for e in config.get("config", {}).get("Env", []))
        if env.get("PYTHON_VERSION") != "3.12.13":
            raise SystemExit(f"python:{TAG} carries Python {env.get('PYTHON_VERSION')}")
        data = {
            "registry": "registry-1.docker.io",
            "repository": REPOSITORY,
            "tag": TAG,
            "reference": f"python:{TAG}@{amd64}",
            "index_digest": index_digest,
            "index_media_type": index["mediaType"],
            "platforms_in_index": sorted(
                f"{m['platform']['os']}/{m['platform']['architecture']}"
                + (f"/{m['platform']['variant']}" if m["platform"].get("variant") else "")
                for m in index["manifests"] if m.get("platform", {}).get("os") != "unknown"),
            "linux_amd64_manifest_digest": amd64,
            "linux_amd64_manifest_media_type": manifest.get("mediaType"),
            "linux_amd64_config_digest": config_digest,
            "config_created": config.get("created"),
            "config_architecture": config.get("architecture"),
            "config_os": config.get("os"),
            "config_env": {k: env[k] for k in sorted(env) if k.startswith(("PYTHON", "LANG", "GPG"))},
        }
        _environment.update("base_image", data, __file__, log)
        log.note("base_image", data)
        print(f"python:{TAG}\n  index          {index_digest}\n  linux/amd64    {amd64}\n"
              f"  config         {config_digest}\n  PYTHON_VERSION {env['PYTHON_VERSION']}")


if __name__ == "__main__":
    main()
