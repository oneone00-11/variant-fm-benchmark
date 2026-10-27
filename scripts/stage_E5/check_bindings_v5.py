"""E5.1 (commit 1): the gate's cells mode on the v5 texts, number by number against E5.0.

Runs scripts/check_evidence_manuscript_numbers.py --cells docs/manuscript_bindings_v5.tsv on
the v5 main text and on its supplementary notes, as a reader would, and records its output
whole. Then compares its verdict on every number with the verdict E5.0 reached (E5.0 §2;
tests/fixtures/manuscript_bindings_v5_expected.tsv, from build_bindings_v5.py): the status
must be the same, and where the number does not match, the value now must be E5.0's at the
printed precision. One number more or less fails it.

It also reads every measured binding at 8e0d609, the release v5 cites, from that commit's
files (git show, into a temporary folder): there each one must round to its printed
number, which shows the binding is the cell the sentence was written from. Two cannot be
read there and are expected not to be: the median of the folds' thresholds, whose column
E18 renamed (e3_ to pp3_logo_median_, cells unchanged), and the AlphaGenome scoring record,
which E17 added.

Run from the repository root (no network):
    .venv/bin/python scripts/stage_E5/check_bindings_v5.py \\
        --main ~/Desktop/整理文件/calibration_manuscript_v5_round2_keyed.md \\
        --supplement ~/Desktop/整理文件/calibration_supplement_notes_v5_keyed.md
"""
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
import tempfile
from decimal import Decimal
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog  # noqa: E402

sys.path.insert(0, str(REPO / "scripts"))
import check_evidence_manuscript_numbers as gate  # noqa: E402

TABLE = REPO / "docs" / "manuscript_bindings_v5.tsv"
EXPECTED = REPO / "tests" / "fixtures" / "manuscript_bindings_v5_expected.tsv"
V301 = "8e0d609"
NOT_AT_V301 = {("main", 101, "51.8"): "column e3_strong_lr_here, renamed pp3_logo_median_strong_lr_here by E18",
               ("supplement", 18, "16,384"): "the AlphaGenome scoring record, added by E17"}


def numbers(text: str) -> list[Decimal]:
    """The numbers of a value as E5.0 wrote it: one, or the two ends of a range."""
    out = []
    for part in text.replace("%", "").replace(",", "").split("–"):
        try:
            out.append(Decimal(part.strip()))
        except Exception:  # noqa: BLE001
            return []
    return out


def same_value(r: gate.CellResult, e: dict) -> bool:
    if r.now == e["expected_now"]:
        return True
    want = numbers(e["expected_now"])
    if not want:                      # E5.0 wrote words (the scoring interval): status only
        return True
    got = [min(r.values), max(r.values)] if len(want) == 2 else r.values[:1]
    if len(got) != len(want):
        return False
    scale = Decimal(r.binding.scale)
    q = [Decimal(1).scaleb(-(len(str(w).split(".")[1]) if "." in str(w) else 0)) for w in want]
    return all((g * scale).quantize(qq) == w for g, w, qq in zip(got, want, q))


def at_commit(commit: str, bindings: list[gate.Binding], log: RunLog) -> list[gate.CellResult]:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        for f in sorted({b.file for b in bindings if b.kind == "measured"}):
            run = subprocess.run(["git", "-C", str(REPO), "show", f"{commit}:{f}"], capture_output=True)
            if run.returncode == 0:
                (root / f).parent.mkdir(parents=True, exist_ok=True)
                (root / f).write_bytes(run.stdout)
        return [gate.evaluate(b, None, root) for b in bindings if b.kind == "measured"]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--main", type=Path, required=True)
    ap.add_argument("--supplement", type=Path, required=True)
    args = ap.parse_args()
    texts = {"main": args.main.expanduser(), "supplement": args.supplement.expanduser()}
    with RunLog(__file__) as log:
        log.input(REPO / "scripts/check_evidence_manuscript_numbers.py")
        for p in (TABLE, EXPECTED, *texts.values()):
            log.input(p)
        bindings = gate.read_bindings(TABLE)
        with open(EXPECTED, newline="") as fh:
            expected = {(r["document"], int(r["line"]), r["printed"], int(r["occurrence"])): r
                        for r in csv.DictReader(fh, delimiter="\t")}
        failures, counts, cli = [], {}, {}
        for doc, md in texts.items():
            run = subprocess.run([sys.executable, str(REPO / "scripts/check_evidence_manuscript_numbers.py"),
                                  str(md), "--cells", str(TABLE), "--document", doc],
                                 cwd=REPO, capture_output=True, text=True)
            cli[doc] = {"returncode": run.returncode, "output": run.stdout.splitlines()}
            print(run.stdout)
            mine = [b for b in bindings if b.document == doc]
            results = gate.cells_results(md.read_text(), mine)
            for r in results:
                t, b = r.token, r.binding
                key = (doc, b.line if b else t.line, b.printed if b else t.printed,
                       b.occurrence if b else t.occurrence)
                e = expected.pop(key, None)
                counts[f"{doc} / {r.status}"] = counts.get(f"{doc} / {r.status}", 0) + 1
                if e is None:
                    failures.append(f"{key}: {r.status}, not a number of E5.0's table")
                elif r.status != e["expected_status"]:
                    failures.append(f"{key}: {r.status}, E5.0 {e['e5_0_category']} -> expected "
                                    f"{e['expected_status']} ({r.why})")
                elif r.status == "mismatch" and not same_value(r, e):
                    failures.append(f"{key}: now {r.now}, E5.0 now {e['expected_now']}")
        failures += [f"{k}: in E5.0's table, not read" for k in expected]
        v301 = at_commit(V301, bindings, log)
        v301_bad = []
        for r in v301:
            b = r.binding
            key = (b.document, b.line, b.printed)
            if r.status == "match":
                continue
            if r.status == "cannot_bind" and key in NOT_AT_V301:
                continue
            v301_bad.append(f"{key}: {r.status} {r.now or r.why}")
        print(f"at {V301}: {sum(r.status == 'match' for r in v301)} of {len(v301)} measured bindings round "
              f"to their printed number; {len(NOT_AT_V301)} cannot be read there as expected "
              f"({'; '.join(NOT_AT_V301.values())}); other: {len(v301_bad)}")
        for x in v301_bad:
            print("  " + x)
        for k, v in sorted(counts.items()):
            print(f"{v:4d}  {k}")
        print(f"against E5.0, number by number: {len(failures)} difference(s)")
        for x in failures:
            print("  " + x)
        log.note("cli", cli)
        log.note("status_counts", counts)
        log.note("differences_from_e5_0", failures)
        log.note("v301", {"commit": V301, "measured": len(v301),
                          "match": sum(r.status == "match" for r in v301),
                          "expected_not_readable": {" ".join(map(str, k)): v for k, v in NOT_AT_V301.items()},
                          "other": v301_bad})
        if failures or v301_bad:
            raise SystemExit("the cells mode does not reproduce E5.0: see the run record")
        print("the cells mode reproduces E5.0 §2 number by number")
    return 0


if __name__ == "__main__":
    sys.exit(main())
