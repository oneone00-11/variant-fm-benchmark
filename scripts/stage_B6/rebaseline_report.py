"""B6.3: what running every stage in the reference image changed against the committed
outputs, from the comparator's own report, written to docs/rebaseline-report.md.

Reads reproduce_report.txt, the per-file report `scripts/reproduce_evidence.py --verify`
writes, from a run in the image against the committed outputs (HEAD), and the files
that run left in the working tree. The comparator alone grades the files; this lists:

  * its verdict line and its count of files per status, quoted as it wrote them;
  * every file it did not call identical, with its status;
  * for each file it called "differs (numerically equal)": the columns (keys, for a
    JSON file) whose values moved, how many cells moved, and the largest absolute and
    relative difference seen between HEAD's file and the working tree's, read the way
    the comparator reads them;
  * for each file it called "differs (hashes of numerically equal files)": the files
    whose hashes changed in it;
  * the image the run used (container/ENVIRONMENT.json: the reference by manifest
    digest, the image ID) and the machine it ran on (the CPU probe record given on the
    command line: model, flags, the OpenBLAS kernels, numpy's dispatch).

Run from the repository root, after the run in the image and before the commit:
    .venv/bin/python scripts/stage_B6/rebaseline_report.py --probe container/cpu_probe/NAME.json
"""
from __future__ import annotations

import argparse
import csv
import difflib
import hashlib
import importlib.util
import io
import json
import math
import re
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel  # noqa: E402
import _environment  # noqa: E402

REPORT = REPO / "reproduce_report.txt"
OUT = REPO / "docs" / "rebaseline-report.md"
NUMERIC = "differs (numerically equal)"
HASHES = "differs (hashes of numerically equal files)"
FIGURE = "differs (figure)"


def comparator():
    spec = importlib.util.spec_from_file_location(
        "reproduce_evidence", REPO / "scripts" / "reproduce_evidence.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def head(path: str) -> bytes:
    return subprocess.run(["git", "-C", str(REPO), "show", f"HEAD:{path}"], check=True,
                          capture_output=True).stdout


class Moves:
    """Per column or key: cells that moved, largest absolute and relative difference."""

    def __init__(self):
        self.cols: dict[str, list] = {}

    def add(self, col: str, a, b, E) -> None:
        a, b = E._cell(a), E._cell(b)
        if not (isinstance(a, (int, float)) and isinstance(b, (int, float))) or a == b:
            return
        if isinstance(a, bool) or isinstance(b, bool) or math.isnan(a) or math.isnan(b):
            return
        d = abs(a - b)
        r = d / max(abs(a), abs(b))
        c = self.cols.setdefault(col, [0, 0.0, 0.0])
        c[0] += 1
        c[1], c[2] = max(c[1], d), max(c[2], r)


def table_moves(rel_path: str, E) -> Moves:
    m, new_path = Moves(), REPO / rel_path
    suffix = new_path.suffix
    if suffix in E.TABLE_DELIMITERS:
        d = E.TABLE_DELIMITERS[suffix]
        x = list(csv.reader(io.StringIO(head(rel_path).decode()), delimiter=d))
        y = list(csv.reader(io.StringIO(new_path.read_text()), delimiter=d))
        names = x[0] if x else []
        for r, t in zip(x[1:], y[1:]):
            for i, (u, v) in enumerate(zip(r, t)):
                m.add(names[i] if i < len(names) else f"column {i + 1}", u, v, E)
    elif suffix == ".json":
        def walk(a, b, key):
            if isinstance(a, dict) and isinstance(b, dict):
                for k in a.keys() & b.keys():
                    if k not in E.WALL_CLOCK_KEYS:
                        walk(a[k], b[k], f"{key}.{k}" if key else k)
            elif isinstance(a, list) and isinstance(b, list):
                for u, v in zip(a, b):
                    walk(u, v, f"{key}[]")
            else:
                m.add(key or "(value)", a, b, E)
        walk(json.loads(head(rel_path)), json.loads(new_path.read_text()), "")
    elif suffix == ".parquet":
        import pandas as pd
        x = pd.read_parquet(io.BytesIO(head(rel_path)))
        y = pd.read_parquet(new_path)
        for col in x.columns:
            for u, v in zip(x[col].tolist(), y[col].tolist()):
                m.add(col, u, v, E)
    elif suffix == ".txt":
        for u, v in zip(head(rel_path).decode().split(), new_path.read_text().split()):
            m.add("(words of the text)", u, v, E)
    return m


def _hashes_of(path: str, E) -> set[str]:
    """Every sha256 that can stand for this file in a manifest, as HEAD has it and as
    the working tree has it: of its bytes; of a JSON's content without its wall-clock
    fields (evid_supp_tables' hash of a source); of a table's rows in any order (the
    content hash the frozen and TP53 manifests record)."""
    out = set()
    with tempfile.TemporaryDirectory() as d:
        for label, blob in (("head", head(path)), ("now", (REPO / path).read_bytes())):
            tmp = Path(d) / f"{label}{Path(path).suffix}"
            tmp.write_bytes(blob)
            out.add(hashlib.sha256(blob).hexdigest())
            if tmp.suffix == ".json":
                out.add(E._source_sha256(tmp))
            elif tmp.suffix == ".parquet":
                out.add(E._canonical_sha256(tmp))
    return out


def hashes_changed(rel_path: str, numeric: list[str], E) -> list[str]:
    """The files whose hashes the manifest records and which changed: every hash, as
    HEAD's manifest has it, that the new manifest no longer carries, by file."""
    old_text, new_text = head(rel_path).decode(), (REPO / rel_path).read_text()
    by_hash = {h: f for f in numeric for h in _hashes_of(f, E)}
    gone = set(re.findall(r"\b[0-9a-f]{64}\b", old_text)) - set(
        re.findall(r"\b[0-9a-f]{64}\b", new_text))
    return sorted({by_hash.get(h, f"(a hash of no numerically equal file: {h[:12]}...)")
                   for h in gone})


def fmt(x: float) -> str:
    return f"{x:.3g}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--probe", type=Path, required=True,
                    help="the CPU probe record of the machine the run used")
    args = ap.parse_args()
    with RunLog(__file__) as log:
        E = comparator()
        log.input(REPORT)
        lines = REPORT.read_text().splitlines()
        header = [l[2:] for l in lines if l.startswith("# ")]
        rows = [l.split("\t", 1) for l in lines if l and not l.startswith("#")]
        verdict, counts = header[1], header[2]
        numeric = [p for s, p in rows if s == NUMERIC]
        env = _environment.read()
        probe = json.loads(args.probe.read_text())
        log.input(args.probe)
        base = subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"], check=True,
                              capture_output=True, text=True).stdout.strip()
        md = ["# Rebaseline in the reference image", "",
              "Written by `scripts/stage_B6/rebaseline_report.py` from the per-file report of "
              "`scripts/reproduce_evidence.py --verify` (`reproduce_report.txt`), after every "
              "stage ran in the reference image (`container/run.sh`: inputs fetched online, "
              "stages and `--verify` with no network), against the committed outputs of "
              f"commit {base[:7]}; do not edit by hand. The comparator alone grades the files; "
              "the differences below are computed here and listed, not judged.", "",
              "## The comparator's verdict", "", "Quoted as it wrote them:", "",
              f"> {verdict}", ">", f"> {counts}", "",
              "## Files not identical", "",
              "| status | file |", "|---|---|"]
        md += [f"| {s} | `{p}` |" for s, p in rows
               if s not in ("identical", "identical (build time aside)")]
        md += ["", "## How far the numerically equal files moved", "",
               "Per file, the columns (for JSON, the keys; `[]` stands for any list "
               "position) whose numbers moved, the cells that moved, and the largest "
               "absolute and relative difference between the committed file and the "
               "image's (relative: the difference over the larger magnitude).", "",
               "| file | column | cells moved | largest absolute | largest relative |",
               "|---|---|---|---|---|"]
        worst = 0.0
        for p in numeric:
            moves = table_moves(p, E).cols
            if not moves:
                md.append(f"| `{p}` | (no number moved; the text did) | | | |")
            for col, (n, d, r) in sorted(moves.items()):
                md.append(f"| `{p}` | `{col}` | {n} | {fmt(d)} | {fmt(r)} |")
                worst = max(worst, r)
        md += ["", f"Largest relative difference over all of them: {fmt(worst)}.", ""]
        plain = [p for s, p in rows if s == "differs"]
        if plain:
            md += ["## Files that differ in more than their numbers' last digits", "",
                   "The comparator called these \"differs\". Their changed lines, HEAD's "
                   "first (`-`), then the image's (`+`); a sha256 that is the hash of a "
                   "file above is named after it.", ""]
            names = {h: f for f in numeric for h in _hashes_of(f, E)}
            names.update({h: f for f in plain for h in _hashes_of(f, E)})
            def named(text: str) -> str:
                return re.sub(r"\b[0-9a-f]{64}\b",
                              lambda m: (f"<sha256 of {names[m.group(0)]}>" if m.group(0)
                                         in names else m.group(0)[:12] + "..."), text)
            for p in plain:
                if Path(p).suffix in E.TABLE_DELIMITERS:      # cell by cell
                    d = E.TABLE_DELIMITERS[Path(p).suffix]
                    x = list(csv.reader(io.StringIO(head(p).decode()), delimiter=d))
                    y = list(csv.reader(io.StringIO((REPO / p).read_text()), delimiter=d))
                    diff = [f"row {r[0]}, {x[0][i]}: {named(u)} -> {named(v)}"
                            for r, t in zip(x[1:], y[1:]) for i, (u, v) in enumerate(zip(r, t))
                            if u != v]
                    if len(x) != len(y):
                        diff.append(f"rows: {len(x) - 1} -> {len(y) - 1}")
                    what = "changed cells"
                else:
                    old = head(p).decode().splitlines()
                    new = (REPO / p).read_text().splitlines()
                    diff = [named(l) for l in difflib.unified_diff(old, new, lineterm="", n=0)
                            if not l.startswith(("---", "+++", "@@"))]
                    what = "changed lines"
                md += [f"`{p}`, {len(diff)} {what}:", "", "```"]
                md += [line if len(line) <= 400 else line[:400] + " ..." for line in diff[:40]]
                if len(diff) > 40:
                    md.append(f"... {len(diff) - 40} more")
                md += ["```", ""]
        hashed = [p for s, p in rows if s == HASHES]
        if hashed:
            md += ["## Manifests that changed only in hashes", "",
                   "| manifest | the files whose hashes changed |", "|---|---|"]
            for p in hashed:
                md.append(f"| `{p}` | " + "<br>".join(f"`{f}`" for f in
                                                       hashes_changed(p, numeric, E)) + " |")
            md.append("")
        img, cpu = env.get("image", {}), probe["cpuinfo"]
        blas = [f"{t.get('prefix')} {t.get('version')}: {t.get('architecture')}"
                for t in probe["threadpool_info"] if t.get("internal_api") == "openblas"]
        md += ["## The image and the machine", "",
               f"- Image: `{img.get('reference')}` (image ID, the config digest: "
               f"`{img.get('config_digest')}`), built at commit "
               f"{str(img.get('built_at_commit'))[:7]} by {img.get('build_run')} "
               "(`container/ENVIRONMENT.json`, section image).",
               f"- CPU: {cpu['model_name']} ({cpu['vendor_id']}), {cpu['processors']} "
               f"processors; flags: {' '.join(cpu['flags'])} (`{rel(args.probe)}`).",
               f"- OpenBLAS kernels: {'; '.join(blas)}.",
               f"- numpy dispatch in use: {' '.join(probe['numpy_cpu']['dispatch_in_use'])}.",
               "- Settings in force: " + ", ".join(
                   f"{k}={v}" for k, v in probe["settings"].items() if v is not None) + ".",
               ""]
        OUT.write_text("\n".join(md))
        log.output(OUT)
        log.note("verdict", verdict)
        log.note("counts", counts)
        print(f"wrote {rel(OUT)}: {verdict}")


if __name__ == "__main__":
    main()
