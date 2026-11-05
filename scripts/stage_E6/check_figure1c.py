"""E6.6, item 3: Figure 1 (c) against Supplementary Figure S1, point by point, and the
figure's measurements; with (a) and (b) against the previous commit's figure.

Both figures are drawn by the layout stage's own makers, into a temporary folder (the
stage's outputs are read, never written), and every array is read back from what was
drawn:
  * (c)'s curve (score; local likelihood ratio, clipped to the axis) and its shade's
    lower and upper edges, against S1's panel of the same predictor and band, element
    for element (no tolerance), and against the file read again by the function both
    use (evid_figures._curve), clipped as S1 clips it;
  * the file's rows, against the points drawn;
  * each threshold mark against evidence_thresholds.csv read at full precision
    (evid_common.read_back), and its height against its tier's boundary; the cut point's
    line against config/walker2023.yaml;
  * the axes: three in Figure 1; (a) and (b) as far from the top, and as large, as in
    the figure before (c), AB_MM tall; (c)'s block a whole number of pixels;
  * the type: every text drawn between the artwork's sizes, inside the figure, and none
    on another;
  * the colours of (c) and of its key: greys only;
  * the files: the stage's Figure 1 equals this drawing byte for byte, and the PNG's
    first rows, as many as the previous commit's PNG has, equal that PNG pixel for pixel.

Run from the repository root; it writes only its run record:
    .venv/bin/python scripts/stage_E6/check_figure1c.py [--base 6ce66eb]
or in the reference image:
    ... IMAGE bash scripts/stage_E6/in_image.sh scripts/stage_E6/check_figure1c.py
"""
from __future__ import annotations

import argparse
import hashlib
import io
import itertools
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog  # noqa: E402

PHASE1 = REPO / "phase1"
CONTENT, S1 = "published_cut_point_by_band", "local_lr_curves"
STEM = Path("reports/evidence/layout/main") / CONTENT
GREY = 0.03          # a grey's channels differ by at most this (the palette's greys: 0.024)


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default="6ce66eb", help="the commit before (c) was added")
    args = ap.parse_args()
    os.chdir(PHASE1)
    sys.path.insert(0, str(PHASE1))
    import matplotlib.colors as mcolors
    import matplotlib.pyplot as plt
    import matplotlib.text as mtext
    import numpy as np
    import yaml
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.collections import PolyCollection
    from PIL import Image
    from src import evid_common as K
    from src import evid_figures as F
    from src import evid_layout as L

    failed: list[str] = []

    def check(name: str, ok: bool, value=None) -> None:
        log.note(name, {"ok": bool(ok), "value": value} if value is not None else bool(ok))
        print(f"  {'ok  ' if ok else 'FAIL'} {name}" + (f": {value}" if value is not None else ""))
        if not ok:
            failed.append(name)

    with RunLog(__file__) as log:
        cfg = L.load_config()
        item = next(i for i in cfg["main"] if i["content"] == CONTENT)
        s1 = next(i for i in cfg["supplementary"] if i["content"] == S1)
        tool, st = item["options"]["c_tool"], item["options"]["c_stratum"]
        curve_file = F.REPORT_DIR / f"interval_lr_{tool}_{st}.csv"
        for p in [L.CONFIG, Path(L.P.CONFIG_PATH), *map(Path, L._sources(item)),
                  STEM.with_suffix(".pdf"), STEM.with_suffix(".png")]:
            log.input(PHASE1 / p)
        log.note("base", args.base)
        log.note("predictor_and_band", {"c_tool": tool, "c_stratum": st})

        # both figures, as the stage draws them, kept after the stage's own _save
        kept, save = [], L._save
        L._save = lambda fig, folder, stem, c: (kept.append(fig), save(fig, folder, stem, c))[1]
        tmp = Path(tempfile.mkdtemp())
        try:
            L.figure_published_cut_point(item, tmp / "main", cfg)
            L.redraw_figure(s1, tmp / "supplementary", cfg)
        finally:
            L._save = save
        fig, sfig = kept
        art = cfg["artwork"]
        dpi = art["dpi"]

        print("Figure 1")
        width_mm, height_mm = (v / F.MM for v in fig.get_size_inches())
        check("axes_in_figure_1", len(fig.axes) == 3, len(fig.axes))
        check("figure_mm", width_mm <= art["double_column_mm"] + 1e-9
              and height_mm <= art["max_height_mm"] + 1e-9,
              {"width": round(width_mm, 4), "height": round(height_mm, 4)})
        c_px = (height_mm - L.AB_MM) / 25.4 * dpi
        check("c_block_whole_pixels", abs(c_px - round(c_px)) < 1e-6,
              {"mm": round(height_mm - L.AB_MM, 4), "pixels": round(c_px, 6)})
        # (a) and (b) where the figure before (c) had them: its own numbers, 175 x AB_MM
        fig0, axes0 = plt.subplots(1, 2, figsize=(L._width(cfg), L.AB_MM * F.MM), sharey=True,
                                   gridspec_kw={"wspace": 0.08})
        fig0.subplots_adjust(left=0.2, right=0.97, top=0.94, bottom=0.2)
        geometry = {}
        for name, ax, ax0 in zip("ab", fig.axes[:2], axes0):
            p, p0 = ax.get_position(), ax0.get_position()
            now = [p.x0 * width_mm, p.width * width_mm, (1 - p.y1) * height_mm, p.height * height_mm]
            was = [p0.x0 * width_mm, p0.width * width_mm, (1 - p0.y1) * L.AB_MM, p0.height * L.AB_MM]
            geometry[name] = {"left_width_top_height_mm": [round(v, 4) for v in now],
                              "largest_difference_mm": max(abs(a - b) for a, b in zip(now, was))}
        plt.close(fig0)
        check("ab_axes_as_before", all(g["largest_difference_mm"] < 1e-9 for g in geometry.values()),
              geometry)
        pc = fig.axes[2].get_position()
        check("c_axes_mm", True, {"left": round(pc.x0 * width_mm, 4), "width": round(pc.width * width_mm, 4),
                                  "bottom": round(pc.y0 * height_mm, 4),
                                  "height": round(pc.height * height_mm, 4)})

        # the type, on a fresh canvas at the artwork's resolution
        canvas = FigureCanvasAgg(fig)
        fig.set_dpi(dpi)
        canvas.draw()
        r = canvas.get_renderer()

        def drawn(t) -> bool:
            if not t.get_visible() or not t.get_text().strip():
                return False
            for ax in fig.axes:          # a tick label whose tick lies outside the view is not drawn
                for axis, lim in ((ax.xaxis, ax.get_xlim()), (ax.yaxis, ax.get_ylim())):
                    for tk in axis.get_major_ticks():
                        if t in (tk.label1, tk.label2):
                            lo, hi = sorted(lim)
                            return tk.get_loc() is not None and lo - 1e-12 <= tk.get_loc() <= hi + 1e-12
            return True

        texts = [(t.get_text(), t.get_fontsize(), t.get_window_extent(r))
                 for t in fig.findobj(mtext.Text) if drawn(t)]
        sizes = sorted({s for _, s, _ in texts})
        small, large = art["text_pt"]
        check("type_pt", small <= sizes[0] and sizes[-1] <= large, sizes)
        w, h = fig.bbox.width, fig.bbox.height
        outside = [s for s, _, b in texts if b.x0 < -0.5 or b.y0 < -0.5 or b.x1 > w + 0.5 or b.y1 > h + 0.5]
        check("texts_inside_the_figure", not outside, {"texts": len(texts), "outside": outside})
        touching = [(a[0], b[0]) for a, b in itertools.combinations(texts, 2) if a[2].overlaps(b[2])]
        check("texts_apart", not touching, touching)
        cut_text = [s for s, _, _ in texts if s.startswith("In-sample")]
        check("threshold_key_names_its_table", cut_text ==
              [f"In-sample PP3 threshold ({L._part_ref(cfg, 'evidence_tiers_and_thresholds', 'thresholds')})"],
              cut_text)

        print("(c) against Supplementary Figure S1")
        c_ax = fig.axes[2]
        # S1 draws a row per predictor and a column per band; it names the bands above
        # its first row and the predictors beside its first column
        i, j = F.COLUMNS.index(tool), F.STRATA.index(st)
        s_ax = sfig.axes[i * len(F.STRATA) + j]
        named = {"row": i, "column": j, "title": s_ax.get_title(), "ylabel": s_ax.get_ylabel()}
        check("s1_panel_is_the_options", (i != 0 or named["title"] == F.STRATUM_LABEL[st])
              and (j != 0 or named["ylabel"] == F.NAME[tool]), named)

        def curve_of(ax):
            lines = [ln for ln in ax.lines if len(ln.get_xdata()) > 2]
            polys = [c for c in ax.collections if isinstance(c, PolyCollection)]
            assert len(lines) == 1 and len(polys) == 1, (len(lines), len(polys))
            x = np.asarray(lines[0].get_xdata(), dtype=float)
            y = np.asarray(lines[0].get_ydata(), dtype=float)
            paths = polys[0].get_paths()
            assert len(paths) == 1, len(paths)
            v = paths[0].vertices
            n = len(x)
            return {"score": x, "local_lr": y, "shade_x": v[1:n + 1, 0], "shade_lower": v[1:n + 1, 1],
                    "shade_upper_x": v[n + 2:2 * n + 2, 0][::-1],
                    "shade_upper": v[n + 2:2 * n + 2, 1][::-1], "vertices": v}

        cc, ss = curve_of(c_ax), curve_of(s_ax)
        c = F._curve(tool, st)
        expect = {"score": c.score.to_numpy(), "local_lr": c.local_lr.clip(lower=L.C_FLOOR, upper=L.C_CEIL).to_numpy(),
                  "shade_x": c.score.to_numpy(), "shade_upper_x": c.score.to_numpy(),
                  "shade_lower": c.lr_lo_gene_boot.clip(lower=L.C_FLOOR).to_numpy(),
                  "shade_upper": c.lr_hi_gene_boot.clip(upper=L.C_CEIL).to_numpy()}
        for key in ("score", "local_lr", "shade_lower", "shade_upper"):
            check(f"{key}_equals_s1", np.array_equal(cc[key], ss[key]),
                  {"points": len(cc[key]), "differing": int((cc[key] != ss[key]).sum())
                   if cc[key].shape == ss[key].shape else "shapes differ"})
            check(f"{key}_equals_the_file", np.array_equal(cc[key], expect[key]))
        check("shade_x_is_the_score", np.array_equal(cc["shade_x"], expect["score"])
              and np.array_equal(cc["shade_upper_x"], expect["score"]))
        check("shade_polygon_equals_s1", np.array_equal(cc["vertices"], ss["vertices"]),
              {"vertices": len(cc["vertices"])})
        lines = (PHASE1 / curve_file).read_text().splitlines()
        check("rows_drawn", len(lines) - 1 == len(cc["score"]) == len(c),
              {"file_rows": len(lines) - 1, "finite_local_lr": len(c), "points_drawn": len(cc["score"])})
        floor, ceil = L.C_FLOOR, L.C_CEIL
        log.note("clipped_to_the_axis", {
            "floor_ceiling": [floor, ceil],
            "local_lr_below_floor": int((c.local_lr < floor).sum()),
            "local_lr_zero": int((c.local_lr == 0).sum()),
            "local_lr_above_ceiling": int((c.local_lr > ceil).sum()),
            "lower_edge_below_floor": int((c.lr_lo_gene_boot < floor).sum()),
            "upper_edge_above_ceiling": int((c.lr_hi_gene_boot > ceil).sum())})
        exact = K.read_back(curve_file)
        digits = {}
        for col in ("score", "local_lr", "lr_lo_gene_boot", "lr_hi_gene_boot"):
            a, b = c[col].to_numpy(), exact[col].to_numpy()
            rel_ = np.abs(a - b) / np.where(b == 0, 1, np.abs(b))
            digits[col] = {"values_differing": int((a != b).sum()), "largest_relative": float(rel_.max())}
        log.note("default_parser_against_read_back", digits)
        print(f"  info  the curve as S1 reads it (pandas' default parser) against read_back: {digits}")

        print("thresholds and cut point")
        thr = K.read_back(F.REPORT_DIR / "evidence_thresholds.csv")
        thr = thr[(thr.tool == tool) & (thr.stratum == st) & (thr.status == "ok")
                  & thr.pp3_threshold_reachable.astype(bool)]
        want = [(float(thr[thr.tier == n.lower()].pp3_threshold_insample.iloc[0]), v)
                for n, v in F.PATH_CUTS if (thr.tier == n.lower()).any()]
        marks = [(float(m.get_xdata()[0]), float(m.get_ydata()[0])) for m in c_ax.lines if m.get_marker() == "|"]
        check("threshold_marks", marks == want, [{"threshold": repr(x), "on_boundary": y} for x, y in marks])
        walker = yaml.safe_load(Path(L.P.CONFIG_PATH).read_text())
        cut = float(walker["thresholds"]["pp3"]["value"])
        vlines = [ln for ln in c_ax.lines if ln.get_linestyle() == "-."]
        check("cut_point_line", len(vlines) == 1 and list(vlines[0].get_xdata()) == [cut, cut], cut)

        print("colours of (c) and its key")
        colours = set()
        for ln in c_ax.lines:
            colours |= {ln.get_color(), ln.get_markeredgecolor(), ln.get_markerfacecolor()}
        for col in c_ax.collections:
            colours |= {tuple(x) for x in col.get_facecolor()} | {tuple(x) for x in col.get_edgecolor()}
        colours |= {t.get_color() for t in c_ax.findobj(mtext.Text) if t.get_text().strip()}
        colours |= {s.get_edgecolor() for s in c_ax.spines.values() if s.get_visible()}
        colours |= {tl.get_color() for tl in c_ax.xaxis.get_ticklines() + c_ax.yaxis.get_ticklines()}
        key = fig.legends[-1]
        for hd in key.legend_handles:
            if hasattr(hd, "get_color"):
                colours |= {hd.get_color(), hd.get_markeredgecolor()}
            else:
                colours |= {hd.get_facecolor(), hd.get_edgecolor()}
        colours |= {t.get_color() for t in key.get_texts()}
        rgba = sorted({mcolors.to_hex(mcolors.to_rgba(x), keep_alpha=True) for x in colours
                       if not (isinstance(x, str) and x == "none") and mcolors.to_rgba(x)[3] > 0})
        greys = [x for x in rgba if max(mcolors.to_rgb(x)) - min(mcolors.to_rgb(x)) <= GREY]
        check("c_in_greys_only", greys == rgba, rgba)
        check("c_without_the_record_status_colours",
              not {F.BLUE, F.ORANGE, F.AQUA} & {x[:7] for x in rgba})

        print("files")
        for suffix in (".pdf", ".png"):
            stage = (PHASE1 / STEM.with_suffix(suffix)).read_bytes()
            here = (tmp / "main" / f"{CONTENT}{suffix}").read_bytes()
            check(f"stage_{suffix[1:]}_is_this_drawing", stage == here, _sha(stage))
        box = re.search(rb"/MediaBox \[\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*\]",
                        (PHASE1 / STEM.with_suffix(".pdf")).read_bytes())
        check("pdf_mm", True, [round(float(box.group(3)) * 25.4 / 72, 3), round(float(box.group(4)) * 25.4 / 72, 3)])
        old_png = subprocess.run(["git", "-C", str(REPO), "show", f"{args.base}:phase1/{STEM}.png"],
                                 check=True, capture_output=True).stdout
        log.note("base_png_sha256", _sha(old_png))
        old = np.asarray(Image.open(io.BytesIO(old_png)))
        new = np.asarray(Image.open(PHASE1 / STEM.with_suffix(".png")))
        rows = old.shape[0]
        same = new.shape[1:] == old.shape[1:] and np.array_equal(new[:rows], old)
        check("ab_pixels_as_before", same, {"base_png": list(old.shape), "png": list(new.shape),
                                             "rows_compared": rows})
        print(f"\n{len(failed)} check(s) failed" + (f": {failed}" if failed else ""))
        log.note("failed", failed)
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
