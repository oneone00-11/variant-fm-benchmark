#!/usr/bin/env bash
# One command from a fresh clone to every output of the evidence-strength study,
# checked against the committed outputs:
#
#     bash reproduce.sh            # everything: about 25 minutes of computation
#     bash reproduce.sh --check    # about five minutes, no data download: every
#                                  # recorded checksum, then the test suite
#
# It builds .venv with Python 3.12 from requirements-evidence.lock.txt (once), downloads
# and checks the two public inputs the repository does not carry (the ClinVar release
# from NCBI, the companion atlas release archive from Zenodo) and the four Ensembl
# chromosomes the E16 genes' reference-base check reads, runs all stages of
# scripts/reproduce_evidence.py, and compares every output with the committed one. It
# exits non-zero unless all of them match. First use downloads about 530 MB: about
# 120 MB of Python packages for .venv, and 410 MB of data (ClinVar 192 MB, the atlas
# 50 MB, the chromosomes 165 MB), which can take an hour when NCBI is slow. Needs
# python3.12 and network access.
set -euo pipefail
cd "$(dirname "$0")"

PY=""
for c in python3.12 python3; do
    if command -v "$c" >/dev/null 2>&1 &&
       "$c" -c 'import sys; sys.exit(sys.version_info[:2] != (3, 12))'; then
        PY="$c"; break
    fi
done
if [ -z "$PY" ]; then
    echo "Python 3.12 is required (the environment is pinned to it); install it and rerun." >&2
    exit 1
fi

if [ ! -x .venv/bin/python ]; then
    echo "building .venv from requirements-evidence.lock.txt with $("$PY" --version)"
    "$PY" -m venv .venv
    .venv/bin/python -m pip install --quiet --disable-pip-version-check \
        -r requirements-evidence.lock.txt
fi

# the run is only comparable in the environment the outputs were made in
if ! diff -q <(grep -v '^#' requirements-evidence.lock.txt | grep -v '^$' | sort) \
             <(.venv/bin/python -m pip freeze | sort) >/dev/null; then
    echo ".venv differs from requirements-evidence.lock.txt; remove .venv and rerun." >&2
    exit 1
fi

if [ "${1:-}" = "--check" ]; then     # the quick path: no download, no stage run
    exec .venv/bin/python scripts/reproduce_evidence.py --check
fi

# Nothing on the machine may reach the outputs: single-threaded BLAS, so floating-point
# sums run in one order on any number of cores (the outputs are identical this way and
# with the default threads), and no user matplotlib configuration. The atlas is checked
# by the entry point: --verify reads only a copy that is the release file for file.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
       VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
MPLCONFIGDIR="$(mktemp -d)"
export MPLCONFIGDIR
trap 'rm -rf "$MPLCONFIGDIR"' EXIT
unset MATPLOTLIBRC

# The run's verdict is the script's exit status; the note printed after it changes
# nothing.
status=0
.venv/bin/python scripts/reproduce_evidence.py --fetch-inputs --verify "$@" || status=$?
cat <<'EOF'

This was the native path, with this machine's own Python. Since stage B6 the reference
environment is a Docker image, in which every committed output was written:
    IMAGE=ghcr.io/oneone00-11/variant-fm-benchmark-reference@sha256:0dc03436cd94cf46effe25160efa3015ab9808852179d58908c72e16e0537776
    docker run --rm --platform linux/amd64 --user "$(id -u):$(id -g)" \
        -v "$PWD:/repo" -w /repo "$IMAGE" bash container/run.sh fetch      # online
    docker run --rm --platform linux/amd64 --network none --user "$(id -u):$(id -g)" \
        -v "$PWD:/repo" -w /repo "$IMAGE" bash container/run.sh verify     # no network
There the run ends REPRODUCED. Run natively, expect three differences from it (README,
the reproduction section):
  - PNG figures whose compressed bytes differ while every pixel is the same;
  - the TP53 manifest's exact count of re-scored values identical to the archived ones
    (24 of 24 on arm64, 23 of 24 in the image);
  - numbers that differ in their last digits.
To compare the PNG figures pixel by pixel with the committed ones:
    .venv/bin/python scripts/stage_B6/compare_png_pixels.py
EOF
exit "$status"
