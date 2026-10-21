"""B6.1: the pushed reference image, as the build workflow reported it and as it loads
here, written to the image section of container/ENVIRONMENT.json.

Reads the artifacts of one run of .github/workflows/container.yml (task build):
reference-image-record (image.json, manifest.json) and reference-image (image.tar, what
`docker save` wrote on the runner). Checks, before anything is recorded:

  * the registry manifest the workflow kept hashes to the manifest digest it reported,
    and names the config digest it reported (the image ID);
  * the archive's sha256 is the one the workflow reported, and the config it carries
    hashes to that config digest;
  * with --load: `docker load` of the archive, then `docker image inspect` of the image
    it names (architecture, labels, ID; under containerd's image store the ID is the
    digest of the archive's own manifest, which is recorded beside the config digest);
  * who can read the package: an anonymous token's request for the manifest, from here
    (the workflow made the same request from the runner).

Run on the host, with the Docker context that runs the image (DOCKER_CONTEXT):
    python scripts/stage_B6/record_image.py --record DIR --archive DIR/image.tar --load
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import tarfile
import urllib.error
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import RunLog, rel, sha256_of  # noqa: E402
import _environment  # noqa: E402

ACCEPT = ("application/vnd.oci.image.manifest.v1+json, "
          "application/vnd.docker.distribution.manifest.v2+json")


def anonymous_status(name: str, tag: str) -> dict:
    """The HTTP status of an anonymous request for the manifest: 200 means anyone may
    pull the image; a private package answers 401, 403 or 404."""
    token = ""
    try:
        with urllib.request.urlopen(
                f"https://ghcr.io/token?scope=repository:{name}:pull", timeout=60) as r:
            token = json.loads(r.read()).get("token", "")
        token_status = 200
    except urllib.error.HTTPError as e:
        token_status = e.code
    req = urllib.request.Request(f"https://ghcr.io/v2/{name}/manifests/{tag}",
                                 headers={"Authorization": f"Bearer {token}",
                                          "Accept": ACCEPT}, method="HEAD")
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            status = r.status
    except urllib.error.HTTPError as e:
        status = e.code
    return {"token_http_status": token_status, "manifest_http_status": status,
            "anyone_can_pull": status == 200}


def archive_config(path: Path) -> tuple[str, str, list]:
    """The config digest the archive's manifest.json names, checked against the
    config's bytes, and the tags it carries."""
    with tarfile.open(path) as t:
        manifest = json.load(t.extractfile("manifest.json"))
        if len(manifest) != 1:
            raise SystemExit(f"{path}: {len(manifest)} images in the archive")
        config_path = manifest[0]["Config"]
        blob = t.extractfile(config_path).read()
        index = json.load(t.extractfile("index.json")) if "index.json" in t.getnames() else {}
    digest = "sha256:" + hashlib.sha256(blob).hexdigest()
    if not config_path.endswith(digest.split(":")[1]):
        raise SystemExit(f"{path}: the config {config_path} hashes to {digest}")
    own = [m.get("digest") for m in index.get("manifests", [])]
    return digest, manifest[0].get("RepoTags", []), own


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--record", type=Path, required=True,
                    help="the reference-image-record artifact's folder")
    ap.add_argument("--archive", type=Path, required=True, help="image.tar")
    ap.add_argument("--load", action="store_true", help="docker load it and inspect it")
    args = ap.parse_args()
    with RunLog(__file__) as log:
        image = json.loads((args.record / "image.json").read_text())
        log.input(args.record / "image.json")
        raw = (args.record / "manifest.json").read_bytes()
        log.input(args.record / "manifest.json")
        if "sha256:" + hashlib.sha256(raw).hexdigest() != image["manifest_digest"]:
            raise SystemExit("manifest.json does not hash to the reported manifest digest")
        if json.loads(raw)["config"]["digest"] != image["config_digest"]:
            raise SystemExit("the manifest names another config than the reported one")
        if image.get("manifest_digest_push_reported") not in (None, image["manifest_digest"]):
            raise SystemExit("docker push reported another manifest digest")
        archive_sha = sha256_of(args.archive)
        log.input(args.archive)
        if archive_sha != image["docker_save_sha256"]:
            raise SystemExit("the archive's sha256 is not the one the workflow reported")
        config_digest, tags, own_manifests = archive_config(args.archive)
        if config_digest != image["config_digest"]:
            raise SystemExit(f"the archive carries config {config_digest}, "
                             f"not {image['config_digest']}")
        data = {
            "name": image["name"], "tag": image["tag"],
            "reference": f"{image['name']}@{image['manifest_digest']}",
            "manifest_digest": image["manifest_digest"],
            "manifest_media_type": image["manifest_media_type"],
            "config_digest": image["config_digest"],
            "image_id_on_runner": image["image_id_on_runner"],
            "architecture": image["architecture"], "os": image["os"],
            "labels": image["labels"], "size_bytes": image["size_bytes"],
            "built_at_commit": image["built_at_commit"], "build_run": image["run_url"],
            "docker_save_sha256": image["docker_save_sha256"],
            "docker_save_bytes": image["docker_save_bytes"],
            "docker_save_config_digest": config_digest,
            "docker_save_manifests": own_manifests, "docker_save_tags": tags,
            "anonymous_manifest_request_from_runner":
                image["anonymous_manifest_request_http_status"],
            "anonymous_manifest_request_from_here":
                anonymous_status(image["name"].removeprefix("ghcr.io/"), image["tag"]),
        }
        if args.load:
            out = subprocess.run(["docker", "load", "-i", str(args.archive)], check=True,
                                 capture_output=True, text=True).stdout.strip()
            ref = f"{image['name']}:{image['tag']}"
            inspect = json.loads(subprocess.run(["docker", "image", "inspect", ref],
                                                check=True, capture_output=True,
                                                text=True).stdout)[0]
            data["loaded_here"] = {
                "docker_load": out, "image": ref, "inspect_id": inspect["Id"],
                "id_is": ("the config digest" if inspect["Id"] == config_digest else
                          "the archive's own manifest digest" if inspect["Id"] in own_manifests
                          else "neither the config nor the archive's manifest digest"),
                "architecture": inspect.get("Architecture"), "os": inspect.get("Os"),
                "labels": inspect.get("Config", {}).get("Labels"),
                "docker_server": subprocess.run(
                    ["docker", "version", "--format", "{{.Server.Version}} "
                     "{{.Server.Os}}/{{.Server.Arch}}"], check=True, capture_output=True,
                    text=True).stdout.strip(),
            }
            if data["loaded_here"]["id_is"].startswith("neither"):
                raise SystemExit(f"the loaded image's ID {inspect['Id']} matches nothing "
                                 "recorded")
            if (inspect.get("Architecture"), inspect.get("Os")) != ("amd64", "linux"):
                raise SystemExit("the loaded image is not linux/amd64")
        _environment.update("image", data, __file__, log)
        log.note("image", data)
        print(json.dumps(data, indent=1))
        print(f"recorded in {rel(_environment.PATH)}")


if __name__ == "__main__":
    main()
