# The reference environment of the evidence-strength study (stage B6): Python 3.12.13
# on Debian bookworm, linux/amd64, with exactly the packages of
# requirements-evidence.lock.txt, each checked against the sha256 of its wheel. It holds
# no data and no code of this repository: the repository is mounted when the image
# runs, and container/run.sh runs the pipeline in it.
#
# The pins below are read from container/ENVIRONMENT.json, which
# scripts/stage_B6/base_image.py (the registry) and scripts/stage_B6/pin_versions.py
# (the base image itself) wrote; tests/test_container.py holds this file to it. Built
# for linux/amd64 only, by .github/workflows/container.yml.
FROM python:3.12.13-slim-bookworm@sha256:6e13e65c55e33adf203d77ee371cf8bf5d81bd4902ef07565721f46bf44917af

SHELL ["/bin/bash", "-euo", "pipefail", "-c"]

RUN apt-get update \
 && apt-get install -y --no-install-recommends \
    # git: scripts/reproduce_evidence.py --verify reads the committed outputs with git
    git=1:2.39.5-0+deb12u3 \
 && rm -rf /var/lib/apt/lists/*

# Wheels only, each checked against its hash, nothing resolved: the lock is the whole
# environment. Then the same check as reproduce.sh: what is installed is the lock.
COPY requirements.linux-amd64.txt /opt/reference/requirements.linux-amd64.txt
RUN python -m pip install --no-cache-dir --disable-pip-version-check --root-user-action=ignore \
        --require-hashes --no-deps --only-binary=:all: \
        -r /opt/reference/requirements.linux-amd64.txt \
 && python -m pip check \
 && diff <(grep -v '^#' /opt/reference/requirements.linux-amd64.txt | grep -v '^$' \
               | cut -d' ' -f1 | sort) \
         <(python -m pip freeze | sort)
