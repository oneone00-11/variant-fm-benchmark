"""B6: the PNG figures compared pixel by pixel with the committed ones.

scripts/reproduce_evidence.py --verify compares files byte for byte and calls a PNG that
differs "differs (figure)". A PNG's bytes also depend on how its image data was
compressed, and Pillow compresses it into other bytes on macOS than in the reference
image's Linux. So this decodes each PNG and says which of three it is:

  bytes identical                  the committed file, byte for byte
  bytes differ, pixels identical   the same image (size and every pixel, as RGBA) in
                                   other bytes
  pixels differ                    another image: how many pixels differ, and the box
                                   they lie in

It reads every PNG the repository tracks under the comparator's output folders, as the
old side has it (default: HEAD) and as the new side has it (default: the working tree,
what the last run wrote), and prints a Markdown table, which --out also writes to a
file. It changes no file the comparator reads and leaves the comparator's verdict as it
is: a report, not a verdict.

Run from the repository root:
    .venv/bin/python scripts/stage_B6/compare_png_pixels.py              # after a run
    .venv/bin/python scripts/stage_B6/compare_png_pixels.py --old A --new B --out FILE
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import io
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel  # noqa: E402

SAME_BYTES = "bytes identical"
SAME_PIXELS = "bytes differ, pixels identical"
OTHER_PIXELS = "pixels differ"


def _git(*args: str) -> bytes:
    return subprocess.run(["git", "-C", str(REPO), *args], check=True,
                          capture_output=True).stdout


def _roots() -> tuple[list[str], object]:
    spec = importlib.util.spec_from_file_location(
        "reproduce_evidence", REPO / "scripts" / "reproduce_evidence.py")
    E = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(E)
    return [str(r.relative_to(REPO)) for r in E.OUTPUT_ROOTS], E


def pngs(commit: str) -> list[str]:
    """The PNGs a commit tracks under the comparator's output folders."""
    roots, E = _roots()
    names = _git("ls-tree", "-r", "--name-only", commit, "--", *roots).decode().split("\n")
    return sorted(n for n in names if n.endswith(".png") and E._in_roots(REPO / n))


def read(side: str | None, path: str) -> bytes | None:
    if side is None:
        p = REPO / path
        return p.read_bytes() if p.is_file() else None
    try:
        return _git("show", f"{side}:{path}")
    except subprocess.CalledProcessError:
        return None


def classify(old: bytes, new: bytes) -> dict:
    """One PNG's two versions: their status, sizes and modes, and where pixels differ."""
    import numpy as np
    from PIL import Image
    out = {"bytes_old": len(old), "bytes_new": len(new)}
    if hashlib.sha256(old).digest() == hashlib.sha256(new).digest():
        return out | {"status": SAME_BYTES}
    a, b = Image.open(io.BytesIO(old)), Image.open(io.BytesIO(new))
    out |= {"size_old": a.size, "size_new": b.size, "mode_old": a.mode, "mode_new": b.mode}
    if a.size != b.size:
        return out | {"status": OTHER_PIXELS, "pixels_differing": None, "box": None}
    x, y = np.asarray(a.convert("RGBA")), np.asarray(b.convert("RGBA"))
    differ = (x != y).any(axis=2)
    if not differ.any():
        return out | {"status": SAME_PIXELS}
    rows, cols = np.nonzero(differ)
    return out | {"status": OTHER_PIXELS, "pixels_differing": int(differ.sum()),
                  "box": (int(cols.min()), int(rows.min()), int(cols.max()) + 1,
                          int(rows.max()) + 1)}


def table(results: dict, old_label: str, new_label: str) -> str:
    counts = {s: sum(r["status"] == s for r in results.values())
              for s in (SAME_BYTES, SAME_PIXELS, OTHER_PIXELS)}
    lines = [f"PNG figures, pixel by pixel: old = {old_label}, new = {new_label}.", "",
             "; ".join(f"{s}: {n}" for s, n in counts.items()), "",
             "| file | status | bytes old -> new | image | where pixels differ |",
             "|---|---|---|---|---|"]
    for path, r in results.items():
        if r["status"] == "absent":
            lines.append(f"| `{path}` | absent on the {r['side']} side | | | |")
            continue
        image = (f"{r['size_old'][0]}x{r['size_old'][1]} {r['mode_old']}"
                 + ("" if (r["size_new"], r["mode_new"]) == (r["size_old"], r["mode_old"])
                    else f" -> {r['size_new'][0]}x{r['size_new'][1]} {r['mode_new']}")
                 if "size_old" in r else "")
        where = ("" if r["status"] != OTHER_PIXELS else
                 "sizes differ" if r["pixels_differing"] is None else
                 f"{r['pixels_differing']:,} pixels, within {r['box']}")
        lines.append(f"| `{path}` | {r['status']} | {r['bytes_old']:,} -> {r['bytes_new']:,} "
                     f"| {image} | {where} |")
    return "\n".join(lines) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--old", default="HEAD", help="the commit to compare against (HEAD)")
    ap.add_argument("--new", default=None, help="the commit to compare (the working tree)")
    ap.add_argument("--out", type=Path, default=None, help="also write the table here")
    ap.add_argument("--note", default="", help="a paragraph to put above the table")
    args = ap.parse_args()
    with RunLog(__file__) as log:
        old_label = f"commit {_git('rev-parse', '--short', args.old).decode().strip()}"
        new_label = ("the working tree" if args.new is None else
                     f"commit {_git('rev-parse', '--short', args.new).decode().strip()}")
        paths = sorted(set(pngs(args.old)) | (set() if args.new is None else set(pngs(args.new))))
        results = {}
        for path in paths:
            old, new = read(args.old, path), read(args.new, path)
            if old is None or new is None:
                results[path] = {"status": "absent", "side": "old" if old is None else "new"}
            else:
                results[path] = classify(old, new)
        text = table(results, old_label, new_label)
        if args.out:
            head = ["# PNG figures, pixel by pixel", "",
                    "Written by `scripts/stage_B6/compare_png_pixels.py`; do not edit by hand.",
                    ""]
            if args.note:
                head += [args.note, ""]
            args.out.write_text("\n".join(head) + "\n" + text)
            log.output(args.out)
            print(f"wrote {rel(args.out)}", file=sys.stderr)
        log.note("results", {p: r["status"] for p, r in results.items()})
        print(text)


if __name__ == "__main__":
    main()
