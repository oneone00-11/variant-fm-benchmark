"""Download the ClinVar GRCh38 VCF (weekly rolling release).

Streams clinvar.vcf.gz to data/raw/ and records the release date (##fileDate)
from the gzip header for reproducibility. Bypasses the system proxy
(trust_env=False) because the local Clash proxy is unreliable.

Run:  python scripts/01_download_clinvar.py
"""
import gzip
import sys

import requests
from tqdm import tqdm

from config import CLINVAR_VCF_URL, CLINVAR_VCF_GZ, DATA_RAW


def download(url: str, dest, chunk=1 << 20):
    sess = requests.Session()
    sess.trust_env = False  # ignore broken registry/env proxy
    with sess.get(url, stream=True, timeout=60) as r:
        r.raise_for_status()
        total = int(r.headers.get("Content-Length", 0))
        tmp = dest.with_suffix(dest.suffix + ".part")
        with open(tmp, "wb") as f, tqdm(
            total=total, unit="B", unit_scale=True, desc=dest.name
        ) as bar:
            for block in r.iter_content(chunk_size=chunk):
                f.write(block)
                bar.update(len(block))
        tmp.replace(dest)
    return dest


def read_filedate(vcf_gz):
    """Pull ##fileDate and ##source from the VCF header."""
    meta = {}
    with gzip.open(vcf_gz, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.startswith("##"):
                break
            if line.startswith("##fileDate"):
                meta["fileDate"] = line.strip().split("=", 1)[1]
            elif line.startswith("##source"):
                meta["source"] = line.strip().split("=", 1)[1]
            elif line.lower().startswith("##reference"):
                meta["reference"] = line.strip().split("=", 1)[1]
    return meta


def main():
    print(f"Downloading ClinVar VCF -> {CLINVAR_VCF_GZ}")
    download(CLINVAR_VCF_URL, CLINVAR_VCF_GZ)
    # also grab the tabix index (handy later)
    try:
        download(CLINVAR_VCF_URL + ".tbi", CLINVAR_VCF_GZ.with_suffix(".gz.tbi"))
    except Exception as e:  # non-fatal
        print(f"  (tbi index download skipped: {e})")

    meta = read_filedate(CLINVAR_VCF_GZ)
    print("ClinVar release metadata:")
    for k, v in meta.items():
        print(f"  {k}: {v}")
    (DATA_RAW / "clinvar_release.txt").write_text(
        "\n".join(f"{k}={v}" for k, v in meta.items()) + "\n", encoding="utf-8"
    )
    size_mb = CLINVAR_VCF_GZ.stat().st_size / 1e6
    print(f"Done. {CLINVAR_VCF_GZ.name} = {size_mb:.1f} MB")


if __name__ == "__main__":
    sys.exit(main())
