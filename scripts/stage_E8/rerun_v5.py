"""E8.13: v5's checks run again, before and after the change to the number gate in E8
commit 7 (a number followed by an en dash or a hyphen and a number with "%" is printed as
a percentage, as the 55 of "55–67%"). Each mode must give the same output after the
change as before it.

- scripts/stage_E5/check_bindings_v5.py, run as it stands (the cells mode on the v5 main
  text and supplementary notes, number by number against E5.0, and the bindings read at
  8e0d609); its run record is written to this folder's run_log/, as
  scripts/stage_E8/check_scope.py runs E18's script;
- the gate's other modes on the same texts, as E5.0 ran them: the bound mode on the main
  text and on the notes (a copy of the notes between "# Abstract" and "# References", the
  span the gate reads, E5.0 §1), the self-test on the main text, and the pool mode on both.

Each output is recorded whole. --compare names two records of this script, one made
before the change and one after: every mode's output must be the same, line for line (a
line naming the run record's own file aside); exit 1 otherwise.

    .venv/bin/python scripts/stage_E8/rerun_v5.py \\
        --main ~/Desktop/整理文件/calibration_manuscript_v5_round2_keyed.md \\
        --supplement ~/Desktop/整理文件/calibration_supplement_notes_v5_keyed.md
    .venv/bin/python scripts/stage_E8/rerun_v5.py --compare run_log/rerun_v5_<before>.json run_log/rerun_v5_<after>.json
"""
from __future__ import annotations

import argparse
import contextlib
import io
import json
import runpy
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from _runlog import REPO, RunLog  # noqa: E402  this folder's run record, found first by E5's script

GATE = REPO / "scripts" / "check_evidence_manuscript_numbers.py"
V5 = REPO / "scripts" / "stage_E5" / "check_bindings_v5.py"


def gate(*args: str) -> dict:
    r = subprocess.run([sys.executable, str(GATE), *args], capture_output=True, text=True, cwd=REPO)
    return {"exit": r.returncode, "output": (r.stdout + r.stderr).splitlines()}


def run(main_md: Path, supp_md: Path) -> None:
    with RunLog(__file__) as log:
        log.input(GATE)
        log.input(main_md)
        log.input(supp_md)
        # E5's script, unchanged, in this process: its RunLog is this folder's
        out = io.StringIO()
        argv = sys.argv
        sys.argv = [str(V5), "--main", str(main_md), "--supplement", str(supp_md)]
        code = 0
        try:
            with contextlib.redirect_stdout(out):
                runpy.run_path(str(V5), run_name="__main__")
        except SystemExit as e:
            code = e.code if isinstance(e.code, int) else (0 if e.code is None else 1)
        finally:
            sys.argv = argv
        print(out.getvalue())
        modes = {"check_bindings_v5": {"exit": code, "output": out.getvalue().splitlines()}}
        with tempfile.TemporaryDirectory() as tmp:
            wrapped = Path(tmp) / supp_md.name
            wrapped.write_text("# Abstract\n\n" + supp_md.read_text() + "\n\n# References\n")
            modes["bound_main"] = gate(str(main_md))
            modes["bound_supplement"] = gate(str(wrapped))
            modes["selftest_main"] = gate(str(main_md), "--selftest")
            modes["pool_main"] = gate(str(main_md), "--pool")
            modes["pool_supplement"] = gate(str(wrapped), "--pool")
        for k, v in modes.items():
            print(f"== {k}: exit {v['exit']}")
            print("\n".join(v["output"][-4:]))
        log.note("modes", modes)


def compare(a: Path, b: Path) -> int:
    ma, mb = (json.loads(Path(p).read_text())["notes"]["modes"] for p in (a, b))
    differ = 0
    for k in sorted(set(ma) | set(mb)):
        xa = [ln for ln in ma.get(k, {}).get("output", []) if not ln.startswith("[run log]")]
        xb = [ln for ln in mb.get(k, {}).get("output", []) if not ln.startswith("[run log]")]
        same = ma.get(k, {}).get("exit") == mb.get(k, {}).get("exit") and xa == xb
        differ += not same
        print(f"{'same     ' if same else 'DIFFERENT'} {k}: exit {ma.get(k, {}).get('exit')} / "
              f"{mb.get(k, {}).get('exit')}, {len(xa)} / {len(xb)} lines")
        if not same:
            for x, y in zip(xa + [""] * len(xb), xb + [""] * len(xa)):
                if x != y:
                    print(f"   before: {x[:160]}\n   after:  {y[:160]}")
                    break
    print(f"{differ} mode(s) differ")
    return 1 if differ else 0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--main", type=Path)
    ap.add_argument("--supplement", type=Path)
    ap.add_argument("--compare", nargs=2, type=Path)
    args = ap.parse_args()
    if args.compare:
        raise SystemExit(compare(*args.compare))
    run(args.main.expanduser().resolve(), args.supplement.expanduser().resolve())


if __name__ == "__main__":
    main()
