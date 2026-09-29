#!/usr/bin/env bash
# One command from a fresh clone to every output of the evidence-strength study,
# checked against the outputs this checkout carries:
#
#     bash reproduce.sh
#
# It builds .venv with Python 3.12 from requirements-evidence.lock.txt (once),
# downloads and checks the two public inputs the repository does not carry (the
# ClinVar release from NCBI, the companion atlas release archive from Zenodo), runs
# all stages of scripts/reproduce_evidence.py, and compares every output with the
# committed one. It exits non-zero unless all of them match. About 25 minutes and
# 250 MB of downloads on first use; needs python3.12 and network access.
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

# Nothing on the machine may reach the outputs: single-threaded BLAS, so floating-point
# sums run in one order on any number of cores (the outputs are identical this way and
# with the default threads), and no user matplotlib configuration.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
       VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
MPLCONFIGDIR="$(mktemp -d)"
export MPLCONFIGDIR
unset MATPLOTLIBRC EVID_ATLAS_REPO

exec .venv/bin/python scripts/reproduce_evidence.py --fetch-inputs --verify "$@"
