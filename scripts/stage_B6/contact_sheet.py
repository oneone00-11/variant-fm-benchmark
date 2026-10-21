"""B6.2, B6.3: the old and the new version of every figure that changed, side by side
(old on the left, new on the right), for the author to check by eye that nothing in
the layout moved. A review aid, not an output: it is written outside the repository.

Takes each PNG the pipeline writes (phase1/reports/evidence/figures/ and the draft
PNGs in fig_data/) as the old commit has it and as the new one (or the working tree)
has it, and for each one that differs writes <out>/<stem>.png: a caption (the file,
both sizes, the share of pixels that differ and the box they lie in), then the two
images. All pairs also go, scaled to one width, into <out>/contact_sheet.pdf, one per
page. The PDF figures are drawn from the same figures as their PNGs and are listed in
<out>/README.txt with whether their bytes changed. Captions are set in the repository's
own Liberation Sans (phase1/assets/fonts/).

Run from the repository root:
    .venv/bin/python scripts/stage_B6/contact_sheet.py --old COMMIT [--new COMMIT] --out DIR
"""
from __future__ import annotations

import argparse
import io
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel  # noqa: E402

FIGURES = ["phase1/reports/evidence/figures", "phase1/reports/evidence/fig_data"]
FONT = REPO / "phase1" / "assets" / "fonts" / "LiberationSans-Regular.ttf"
GAP, CAPTION, PAGE_WIDTH = 40, 120, 3000


def _git(*args: str) -> bytes:
    return subprocess.run(["git", "-C", str(REPO), *args], check=True,
                          capture_output=True).stdout


def _files(commit: str | None, suffix: str) -> dict[str, bytes]:
    if commit is None:
        return {rel(p): p.read_bytes() for d in FIGURES for p in sorted((REPO / d).glob(f"*{suffix}"))}
    names = _git("ls-tree", "-r", "--name-only", commit, "--", *FIGURES).decode().split()
    return {n: _git("show", f"{commit}:{n}") for n in names if n.endswith(suffix)}


def pair(name: str, old: bytes, new: bytes, labels: tuple[str, str]):
    from PIL import Image, ImageChops, ImageDraw, ImageFont
    a = Image.open(io.BytesIO(old)).convert("RGB")
    b = Image.open(io.BytesIO(new)).convert("RGB")
    if a.size == b.size:
        diff = ImageChops.difference(a, b)
        box = diff.getbbox()
        r, g, bl = diff.split()                     # a pixel differs in any channel
        changed = a.size[0] * a.size[1] - ImageChops.lighter(
            ImageChops.lighter(r, g), bl).histogram()[0]
        share = f"{100 * changed / (a.size[0] * a.size[1]):.2f}% of pixels differ"
        where = f", within {box}" if box else ""
    else:
        share, where = "the sizes differ", ""
    w, h = a.size[0] + GAP + b.size[0], CAPTION + max(a.size[1], b.size[1])
    sheet = Image.new("RGB", (w, h), "white")
    sheet.paste(a, (0, CAPTION))
    sheet.paste(b, (a.size[0] + GAP, CAPTION))
    draw = ImageDraw.Draw(sheet)
    font = ImageFont.truetype(str(FONT), 34)
    draw.text((10, 10), f"{name}    {a.size[0]}x{a.size[1]} -> {b.size[0]}x{b.size[1]};"
                        f" {share}{where}", fill="black", font=font)
    draw.text((10, 60), f"left: {labels[0]}", fill="black", font=font)
    draw.text((a.size[0] + GAP + 10, 60), f"right: {labels[1]}", fill="black", font=font)
    return sheet, share + where


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--old", required=True, help="the commit whose figures are on the left")
    ap.add_argument("--new", default=None, help="the commit on the right (default: the "
                                                "working tree)")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    with RunLog(__file__) as log:
        labels = (f"commit {args.old}", f"commit {args.new}" if args.new else "working tree")
        old, new = _files(args.old, ".png"), _files(args.new, ".png")
        args.out.mkdir(parents=True, exist_ok=True)
        pages, lines = [], []
        for name in sorted(set(old) | set(new)):
            if name not in old or name not in new:
                lines.append(f"{name}: only in {'the new' if name in new else 'the old'}")
                continue
            if old[name] == new[name]:
                lines.append(f"{name}: identical")
                continue
            sheet, what = pair(name, old[name], new[name], labels)
            path = args.out / f"{Path(name).stem}.png"
            sheet.save(path)
            log.output(path)
            scale = PAGE_WIDTH / sheet.size[0]
            pages.append(sheet.resize((PAGE_WIDTH, round(sheet.size[1] * scale))))
            lines.append(f"{name}: {what} -> {path.name}")
        pdfs_old, pdfs_new = _files(args.old, ".pdf"), _files(args.new, ".pdf")
        for name in sorted(set(pdfs_old) | set(pdfs_new)):
            same = pdfs_old.get(name) == pdfs_new.get(name)
            lines.append(f"{name}: {'identical' if same else 'differs'} (drawn from the "
                         "same figure as its PNG)")
        if pages:
            pages[0].save(args.out / "contact_sheet.pdf", save_all=True,
                          append_images=pages[1:], resolution=150)
            log.output(args.out / "contact_sheet.pdf")
        readme = args.out / "README.txt"
        readme.write_text(f"Old (left): {labels[0]}\nNew (right): {labels[1]}\n\n"
                          + "\n".join(lines) + "\n")
        log.output(readme)
        log.note("pairs", lines)
        print("\n".join(lines))


if __name__ == "__main__":
    main()
