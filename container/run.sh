#!/usr/bin/env bash
# The evidence-strength pipeline inside the reference image (container/Dockerfile). It
# runs in the image, with this repository mounted as the working directory, as the
# host user, in two steps with the same image:
#
#   docker run --rm --platform linux/amd64 --user "$(id -u):$(id -g)" \
#       -v "$PWD:/repo" -w /repo IMAGE bash container/run.sh fetch      # online
#   docker run --rm --platform linux/amd64 --network none --user "$(id -u):$(id -g)" \
#       -v "$PWD:/repo" -w /repo IMAGE bash container/run.sh verify     # offline
#
#   fetch        download the inputs the repository does not carry into its gitignored
#                folders and check them (reproduce_evidence.py --fetch-only); no stage
#   verify       every stage, then every output compared with the committed one
#                (reproduce_evidence.py --verify, whose verdict decides the exit code)
#   probe NAME   the CPU, the BLAS kernels and numpy's dispatch under these settings,
#                to container/cpu_probe/NAME.json (scripts/stage_B6/cpu_probe.py)
#   environment  the image's interpreter and packages, to container/ENVIRONMENT.json
#                (scripts/stage_B6/record_environment.py)
set -euo pipefail
cd "$(dirname "$0")/.."

if [ "$(id -u)" = 0 ]; then
    echo "run as the host user (docker run --user \"\$(id -u):\$(id -g)\"), so that the" \
         "files written into the repository are the user's" >&2
    exit 2
fi

# The host user has no home folder in the image: a temporary one, for git's
# configuration, and a temporary matplotlib cache, both removed on exit.
HOME="$(mktemp -d)"
MPLCONFIGDIR="$(mktemp -d)"
export HOME MPLCONFIGDIR
trap 'rm -rf "$HOME" "$MPLCONFIGDIR"' EXIT
git config --global --add safe.directory "$PWD"     # mounted, so owned by another uid

# reproduce.sh's settings: single-threaded BLAS, so that floating-point sums run in one
# order on any number of cores, and no user matplotlib configuration. Then the hash
# seed, the time zone and the locale, fixed rather than inherited.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
       VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
unset MATPLOTLIBRC EVID_ATLAS_REPO
export PYTHONHASHSEED=0 TZ=UTC LC_ALL=C.UTF-8

# B6.4: the OpenBLAS kernel type and the numpy dispatch every CPU is held to, so that
# two machines run the same instructions. Empty: each library chooses for the CPU it
# finds. The values, and the probe records they rest on (container/cpu_probe/), are
# set before the outputs are rebaselined in the image.
OPENBLAS_CORETYPE_PIN=""
NPY_DISABLE_CPU_FEATURES_PIN=""
unset OPENBLAS_CORETYPE NPY_DISABLE_CPU_FEATURES NPY_ENABLE_CPU_FEATURES
if [ -n "$OPENBLAS_CORETYPE_PIN" ]; then
    export OPENBLAS_CORETYPE="$OPENBLAS_CORETYPE_PIN"
fi
if [ -n "$NPY_DISABLE_CPU_FEATURES_PIN" ]; then
    export NPY_DISABLE_CPU_FEATURES="$NPY_DISABLE_CPU_FEATURES_PIN"
fi

# the run is only comparable in the image's own environment: what is installed must be
# the lock, as the image's build checked
LOCK=container/requirements.linux-amd64.txt
if ! diff -q <(grep -v '^#' "$LOCK" | grep -v '^$' | cut -d' ' -f1 | sort) \
             <(python -m pip freeze | sort) >/dev/null; then
    echo "the interpreter's packages differ from $LOCK: run this inside the reference image" >&2
    exit 1
fi

case "${1:-}" in
    fetch)       python scripts/reproduce_evidence.py --fetch-only ;;
    verify)      python scripts/reproduce_evidence.py --verify ;;
    probe)       python scripts/stage_B6/cpu_probe.py "${2:?probe needs a name}" ;;
    environment) python scripts/stage_B6/record_environment.py ;;
    *)           echo "usage: bash container/run.sh fetch | verify | probe NAME | environment" >&2
                 exit 2 ;;
esac
