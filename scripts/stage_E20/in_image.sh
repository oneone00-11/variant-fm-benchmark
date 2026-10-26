#!/usr/bin/env bash
# E20: one of this folder's scripts in the reference image, so that every file it
# commits is written by the image's interpreter and packages (修改单 §0 rule 9). A copy
# of scripts/stage_E19/in_image.sh for this folder; the settings are container/run.sh's,
# the same lines:
#
#   docker run --rm --platform linux/amd64 --network none --user "$(id -u):$(id -g)" \
#       -v "$PWD:/repo" -w /repo IMAGE bash scripts/stage_E20/in_image.sh \
#       scripts/stage_E20/pipeline_from.py 39
set -euo pipefail
cd "$(dirname "$0")/../.."

if [ "$(id -u)" = 0 ]; then
    echo "run as the host user (docker run --user \"\$(id -u):\$(id -g)\")" >&2
    exit 2
fi
case "${1:-}" in
    scripts/stage_E20/*.py) ;;
    *) echo "usage: bash scripts/stage_E20/in_image.sh scripts/stage_E20/<script>.py [arguments]" >&2
       exit 2 ;;
esac

HOME="$(mktemp -d)"
MPLCONFIGDIR="$(mktemp -d)"
export HOME MPLCONFIGDIR
trap 'rm -rf "$HOME" "$MPLCONFIGDIR"' EXIT
git config --global --add safe.directory "$PWD"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
       VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
unset MATPLOTLIBRC EVID_ATLAS_REPO
export PYTHONHASHSEED=0 TZ=UTC LC_ALL=C.UTF-8
unset OPENBLAS_CORETYPE NPY_DISABLE_CPU_FEATURES NPY_ENABLE_CPU_FEATURES
export OPENBLAS_CORETYPE="Nehalem"
export NPY_DISABLE_CPU_FEATURES="AVX F16C FMA3 AVX2 AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"

LOCK=container/requirements.linux-amd64.txt
if ! diff -q <(grep -v '^#' "$LOCK" | grep -v '^$' | cut -d' ' -f1 | sort) \
             <(python -m pip freeze | sort) >/dev/null; then
    echo "the interpreter's packages differ from $LOCK: run this inside the reference image" >&2
    exit 1
fi

exec python "$@"
