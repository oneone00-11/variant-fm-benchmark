"""B6.4: what the CPU offers and what the maths libraries choose on it, inside the
reference image, under container/run.sh's settings.

The image fixes the software but not the CPU: the OpenBLAS copies numpy and scipy
carry pick a kernel for the CPU they find, and numpy dispatches its own loops by
instruction set. Two machines give the same outputs byte for byte only if both run the
same code paths, so this records, for one machine, to container/cpu_probe/<name>.json:

  * /proc/cpuinfo: the vendor, the model name and the flags of the first processor;
  * threadpoolctl.threadpool_info() in full, once numpy, scipy and scikit-learn are
    imported: each OpenBLAS copy reports the kernel architecture it chose;
  * numpy.show_runtime() as printed, and numpy's own lists of its SIMD baseline, its
    dispatch targets and the features it found and will use here;
  * pyarrow.runtime_info() (the SIMD level Arrow's kernels use) and Pillow's codec
    versions, which also choose code by CPU;
  * the settings in force: OPENBLAS_CORETYPE, NPY_DISABLE_CPU_FEATURES and the thread
    counts, and the packages installed (pip freeze).

Run inside the image:  bash container/run.sh probe <name>
"""
from __future__ import annotations

import contextlib
import io
import json
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel, utc_now  # noqa: E402

OUT_DIR = REPO / "container" / "cpu_probe"
SETTINGS = ["OPENBLAS_CORETYPE", "NPY_DISABLE_CPU_FEATURES", "NPY_ENABLE_CPU_FEATURES",
            "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "PYTHONHASHSEED", "TZ",
            "LC_ALL", "ARROW_USER_SIMD_LEVEL"]


def cpuinfo() -> dict:
    first: dict = {}
    for line in Path("/proc/cpuinfo").read_text().splitlines():
        if not line.strip():
            if first:
                break
            continue
        key, _, value = line.partition(":")
        first[key.strip()] = value.strip()
    return {"vendor_id": first.get("vendor_id"), "model_name": first.get("model name"),
            "cpu_family": first.get("cpu family"), "model": first.get("model"),
            "flags": sorted(first.get("flags", "").split()),
            "processors": sum(1 for line in Path("/proc/cpuinfo").read_text().splitlines()
                              if line.startswith("processor"))}


def main() -> None:
    name = sys.argv[1]
    with RunLog(__file__) as log:
        import numpy as np
        import scipy
        import scipy.linalg  # noqa: F401  (loads scipy's OpenBLAS)
        import sklearn
        import sklearn.linear_model  # noqa: F401
        import threadpoolctl
        from numpy._core import _multiarray_umath as mu
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            np.show_runtime()
        try:
            import pyarrow
            arrow = {k: str(v) for k, v in pyarrow.runtime_info()._asdict().items()}
        except Exception as e:                     # recorded, not fatal
            arrow = {"error": repr(e)}
        from PIL import features
        record = {
            "name": name,
            "probed_utc": utc_now(),
            "cpuinfo": cpuinfo(),
            "uname": {k: getattr(os.uname(), k)            # not the container's hostname
                      for k in ("sysname", "release", "version", "machine")},
            "settings": {k: os.environ.get(k) for k in SETTINGS},
            "threadpool_info": threadpoolctl.threadpool_info(),
            "numpy_show_runtime": buf.getvalue(),
            "numpy_cpu": {
                "baseline": list(mu.__cpu_baseline__),
                "dispatch": list(mu.__cpu_dispatch__),
                "found": sorted(k for k, v in mu.__cpu_features__.items() if v),
                "dispatch_in_use": [d for d in mu.__cpu_dispatch__
                                    if mu.__cpu_features__.get(d)],
            },
            "pyarrow_runtime_info": arrow,
            "pillow": {f: features.version(f) for f in ("zlib", "jpg", "freetype2")},
            "versions": {"numpy": np.__version__, "scipy": scipy.__version__,
                         "scikit-learn": sklearn.__version__,
                         "threadpoolctl": threadpoolctl.__version__},
            "pip_freeze": subprocess.run([sys.executable, "-m", "pip", "freeze"],
                                         capture_output=True, text=True,
                                         check=True).stdout.split(),
        }
        OUT_DIR.mkdir(parents=True, exist_ok=True)
        out = OUT_DIR / f"{name}.json"
        out.write_text(json.dumps(record, indent=1, ensure_ascii=False) + "\n")
        log.output(out)
        log.note("cpu", record["cpuinfo"]["model_name"])
        log.note("openblas", [(t.get("prefix"), t.get("architecture"))
                              for t in record["threadpool_info"]
                              if t.get("internal_api") == "openblas"])
        log.note("numpy_dispatch_in_use", record["numpy_cpu"]["dispatch_in_use"])
        print(f"{rel(out)}: {record['cpuinfo']['model_name']}")
        for t in record["threadpool_info"]:
            print(f"  {t.get('internal_api')} {t.get('prefix')} {t.get('version')} "
                  f"architecture={t.get('architecture')} threads={t.get('num_threads')}")
        print(f"  numpy dispatch in use: {' '.join(record['numpy_cpu']['dispatch_in_use'])}")
        print(f"  pyarrow: {arrow}")


if __name__ == "__main__":
    main()
