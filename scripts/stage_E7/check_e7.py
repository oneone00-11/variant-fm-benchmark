"""E7.8, items 2-4: what E7's commit changed in the outputs, cell by cell, row by row and
pixel by pixel, against the commit before it (--base), and the measurements the report
gives. It reads the outputs and draws Figure 2 and Supplementary Figure S1 again by the
layout stage's own makers into a temporary folder; it writes only its run record.

  E7.1  S32: the two change columns gone, every other cell as before; the source keeps
        them, and they are the difference of the printed thresholds and that difference
        over the seven-gene threshold's magnitude.
  E7.2  S6: the combined Atlas score's nine rows counted without the folds with BRCA1 or
        RAD51C held out, their Moderate and Strong counts Table 2(a)'s k/n (m) cell for
        cell, the two folds marked "†"; every other row byte for byte the printed table;
        the rows' counts and medians before and after. S7 byte for byte, with its note.
        Figure 2: the two folds in grey in the combined Atlas score's rows, the medians
        over the other folds, the key's new entry in its first column; against the
        previous PNG, pixels differ only in those rows and in the key.
  E7.3  S1: 175 x 240 mm; the space between the caption line and the last row's tick
        labels, by the stage's measure (font boxes, bitmap and PDF), by matplotlib's
        text boxes and by the ink of the PNG, before and after; the panels' heights.
  E7.4  no spelling of "interval calibration" in layout.yaml, evid_layout.py, any text
        the stage writes (CSV, the workbook's parts) or any text its figures draw.
  E7.5  predictor_training_provenance.csv: one cell, the fusion row's training data;
        supplement_tables_manifest.csv: one entry of S3's source hashes, that file's;
        the superseded labels' README: the one sentence.
  and   contents.csv, manifest.csv and the workbook: every cell that changed, listed.

Run from the repository root:
    .venv/bin/python scripts/stage_E7/check_e7.py [--base 9b8eae4]
or in the reference image:
    ... IMAGE bash scripts/stage_E7/in_image.sh scripts/stage_E7/check_e7.py
"""
from __future__ import annotations

import argparse
import difflib
import hashlib
import io
import os
import re
import subprocess
import sys
import tempfile
import zipfile
from contextlib import ExitStack
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog  # noqa: E402

PHASE1 = REPO / "phase1"
OUT = Path("reports/evidence/layout")
OLD_NAME = re.compile(r"(?i)interval[\s_-]*calibration")
AVI_LABEL = "Atlas combined score"


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _git_show(commit: str, path: str) -> bytes:
    return subprocess.run(["git", "-C", str(REPO), "show", f"{commit}:{path}"], check=True,
                          capture_output=True).stdout


def _sheets(data: bytes) -> dict[str, dict[tuple[int, str], str]]:
    """The workbook's sheets by name, each {(row, column): text}, as the stage writes them
    (inline strings and plain numbers)."""
    z = zipfile.ZipFile(io.BytesIO(data))
    names = re.findall(r'<sheet name="([^"]+)" sheetId="(\d+)"', z.read("xl/workbook.xml").decode())
    out = {}
    for name, i in names:
        xml = z.read(f"xl/worksheets/sheet{i}.xml").decode()
        cells = {}
        for ref, body in re.findall(r'<c r="([A-Z]+\d+)"[^>]*>(.*?)</c>', xml):
            m = re.match(r"([A-Z]+)(\d+)", ref)
            text = re.search(r"<t[^>]*>(.*?)</t>|<v>(.*?)</v>", body)
            val = text.group(1) if text.group(1) is not None else text.group(2)
            cells[(int(m.group(2)), m.group(1))] = (val.replace("&amp;", "&").replace("&lt;", "<")
                                                    .replace("&gt;", ">").replace("&quot;", '"'))
        out[name] = cells
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default="9b8eae4", help="the commit before E7")
    args = ap.parse_args()
    os.chdir(PHASE1)
    sys.path.insert(0, str(PHASE1))
    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    from matplotlib.colors import to_rgba
    from PIL import Image
    from src import evid_common as K
    from src import evid_figures as F
    from src import evid_layout as L
    from src import evid_supp_tables as S

    failed: list[str] = []

    def check(name: str, ok: bool, value=None) -> None:
        log.note(name, {"ok": bool(ok), "value": value} if value is not None else bool(ok))
        shown = "" if value is None else f": {value}"
        print(f"  {'ok  ' if ok else 'FAIL'} {name}{shown[:300]}")
        if not ok:
            failed.append(name)

    def base_table(rel: str) -> pd.DataFrame:
        return pd.read_csv(io.BytesIO(_git_show(args.base, f"phase1/{rel}")), dtype=str,
                           keep_default_na=False)

    def table(rel: str) -> pd.DataFrame:
        return pd.read_csv(rel, dtype=str, keep_default_na=False)

    with RunLog(__file__) as log:
        log.note("base", args.base)
        cfg = L.load_config()
        items = {i["content"]: i for i in cfg["main"] + cfg["supplementary"]}

        # ---------------------------------------------------------------- E7.1 S32
        print("E7.1 Supplementary Table S32")
        it = items["twelve_gene_refit"]
        part = it["parts"][0]
        rel = OUT / "supplementary/twelve_gene_refit.csv"
        old, new = base_table(str(rel)), table(rel)
        gone = [c for c in old.columns if c not in new.columns]
        kept = [c for c in old.columns if c not in gone]
        check("s32_columns_left_out", gone == ["Threshold change, twelve-gene minus seven-gene",
                                               "Threshold change relative to the seven-gene threshold's magnitude"], gone)
        check("s32_other_cells_as_before", list(new.columns) == kept and new.equals(old[kept]),
              {"rows": len(new), "columns": [len(old.columns), len(new.columns)]})
        d = K.read_back(part["source"])
        both = d.threshold_7.notna() & d.threshold_12.notna()
        nonzero = both & (d.threshold_7 != 0)
        flips = d[both & (np.sign(d.threshold_7) * np.sign(d.threshold_12) < 0)]
        zero = d[both & (d.threshold_7 == 0)]
        check("s32_source_keeps_them_derived", bool(((d.abs_change - (d.threshold_12 - d.threshold_7)).abs()[both] == 0).all()
              and ((d.rel_change - d.abs_change / d.threshold_7.abs()).abs()[nonzero] == 0).all()
              and d.loc[~both, "abs_change"].isna().all() and d.loc[~nonzero, "rel_change"].isna().all()),
              {"sign_changes": [f"{t} {s} {sd} {tr}: {a:.6g} -> {b:.6g}" for t, s, sd, tr, a, b in
                                zip(flips.tool, flips.stratum, flips.side, flips.tier, flips.threshold_7, flips.threshold_12)],
               "seven_gene_threshold_zero_relative_change_empty": [f"{t} {s} {sd} {tr}" for t, s, sd, tr in
                                                                   zip(zero.tool, zero.stratum, zero.side, zero.tier)]})
        check("s32_source_unchanged", (PHASE1 / part["source"]).read_bytes() == _git_show(args.base, f"phase1/{part['source']}"))

        # ---------------------------------------------------------------- E7.2 S6, S7
        print("E7.2 Supplementary Tables S6 and S7")
        part = items["thresholds_pathogenic"]["parts"][0]
        rel = OUT / "supplementary/thresholds_pathogenic.csv"
        src_bytes = (PHASE1 / part["source"]).read_bytes()
        check("s6_source_unchanged", src_bytes == _git_show(args.base, f"phase1/{part['source']}"))
        old, new = base_table(str(rel)), table(rel)
        check("s6_before_was_the_source", _git_show(args.base, f"phase1/{rel}") == src_bytes)
        before_lines = src_bytes.decode().splitlines()
        after_lines = (PHASE1 / rel).read_text().splitlines()
        atlas = (old.Tool == AVI_LABEL).to_numpy()
        same_others = (len(before_lines) == len(after_lines) and before_lines[0] == after_lines[0]
                       and all(b == a for b, a, x in zip(before_lines[1:], after_lines[1:], atlas) if not x))
        check("s6_other_rows_byte_for_byte", same_others, {"rows": len(old), "atlas_rows": int(atlas.sum())})
        recount = set(L.S6_OVER_FOLDS) | {L.S6_FOLDS}
        moved = sorted({c for i in np.flatnonzero(atlas) for c in old.columns if old.at[i, c] != new.at[i, c]},
                       key=list(old.columns).index)
        check("s6_atlas_rows_change_only_the_counts_summaries_and_per_fold", set(moved) <= recount, moved)
        marks_ok = all(new.at[i, L.S6_FOLDS] == "; ".join(
            e + ("†" if e.split("=")[0] in K.AVI_SEEN_IN_TRAINING else "") for e in old.at[i, L.S6_FOLDS].split("; "))
            for i in np.flatnonzero(atlas))
        check("s6_two_folds_marked", marks_ok)
        t3 = pd.read_csv("reports/evidence/tables/table3_evidence_tiers.csv", dtype=str, keep_default_na=False)
        t3row = t3[t3.Predictor.str.startswith(AVI_LABEL)].iloc[0]
        rows, agree = [], True
        for i in np.flatnonzero(atlas):
            o, n = old.loc[i], new.loc[i]
            band, tier = o["Stratum (intronic offset, bp)"], o.Tier

            def cell(r):
                nn, ev, k = (int(r["Held-out genes"]), int(r["Folds with an evaluable held-out LR"]),
                             int(r["Held-out genes clearing the cut"]))
                return f"{k}/{nn}" + (f" ({nn - ev})" if ev < nn else "")
            t2 = t3row[f"{band} bp: {tier}, held out"] if tier != "Supporting" else ""
            if tier != "Supporting":
                agree &= cell(n) == t2
            rows.append({"band": band, "tier": tier,
                         "before": {c: o[c] for c in L.S6_OVER_FOLDS} | {"k/n (m)": cell(o)},
                         "after": {c: n[c] for c in L.S6_OVER_FOLDS} | {"k/n (m)": cell(n)},
                         "table_2a": t2, "per_fold_after": n[L.S6_FOLDS]})
        log.note("s6_atlas_rows", rows)
        for r in rows:
            print(f"      {r['band']:6s} {r['tier']:10s} before {r['before']['k/n (m)']:9s} "
                  f"median {r['before']['Held-out LR, median'] or '–':5s} | after {r['after']['k/n (m)']:9s} "
                  f"median {r['after']['Held-out LR, median'] or '–':5s} | Table 2(a) {r['table_2a'] or '(no column)'}")
        check("s6_moderate_strong_counts_equal_table_2a", agree)
        part7 = items["thresholds_benign"]["parts"][0]
        check("s7_byte_for_byte", (PHASE1 / OUT / "supplementary/thresholds_benign.csv").read_bytes()
              == (PHASE1 / part7["source"]).read_bytes() == _git_show(args.base, f"phase1/{OUT}/supplementary/thresholds_benign.csv"))

        # ---------------------------------------------------------------- E7.4 S16's column
        print("E7.4 Supplementary Table S16's column")
        rel = OUT / "supplementary/practical_metrics.by_threshold.csv"
        old, new = base_table(str(rel)), table(rel)
        renamed = [(a, b) for a, b in zip(old.columns, new.columns) if a != b]
        check("s16_one_column_renamed", renamed == [("Tier by interval calibration", "Tier by the calibration rule")], renamed)
        check("s16_cells_as_before", old.shape == new.shape and (old.values == new.values).all(), list(new.shape))

        # ---------------------------------------------------------------- contents, manifest, workbook
        print("contents.csv, manifest.csv, the workbook")
        a = _git_show(args.base, f"phase1/{OUT}/contents.csv").decode().splitlines()
        b = (PHASE1 / OUT / "contents.csv").read_text().splitlines()
        changed = [ln[1:] for ln in difflib.unified_diff(a, b, lineterm="", n=0)
                   if ln[:1] in "+-" and not ln.startswith(("---", "+++"))]
        heads = sorted({ln.split(",", 4)[3] for ln in changed})
        check("contents_lines_changed", heads == ["thresholds_benign", "thresholds_pathogenic", "twelve_gene_refit"], heads)
        old, new = base_table(str(OUT / "contents.csv")), table(OUT / "contents.csv")
        cells = {}
        for i in range(len(old)):
            diff = [c for c in old.columns if old.at[i, c] != new.at[i, c]]
            if diff:
                cells[old.at[i, "number"]] = {c: [old.at[i, c], new.at[i, c]] for c in diff}
        log.note("contents_cells", cells)
        for k, v in cells.items():
            print(f"      {k}: {sorted(v)}")
        old, new = base_table(str(OUT / "manifest.csv")), table(OUT / "manifest.csv")
        check("manifest_same_files", list(old.file) == list(new.file), len(new))
        mrows = {}
        for i in range(len(old)):
            diff = [c for c in old.columns if old.at[i, c] != new.at[i, c]]
            if diff:
                mrows[old.at[i, "file"]] = diff
        log.note("manifest_rows", mrows)
        for k, v in mrows.items():
            print(f"      {k.split('/layout/')[1]}: {v}")
        ws_old = _sheets(_git_show(args.base, f"phase1/{OUT}/supplementary_tables.xlsx"))
        ws_new = _sheets((PHASE1 / OUT / "supplementary_tables.xlsx").read_bytes())
        check("workbook_same_sheets", list(ws_old) == list(ws_new), len(ws_new))

        def as_rows(cells: dict) -> list[tuple]:
            """A sheet as its rows, each the texts of its cells in column order."""
            n_rows = max(r for r, _ in cells) if cells else 0
            out_ = [[] for _ in range(n_rows)]
            for (r, c), v in sorted(cells.items(), key=lambda kv: (kv[0][0], len(kv[0][1]), kv[0][1])):
                out_[r - 1].append((c, v))
            return [tuple(v for _, v in row) for row in out_]

        sheets = {}
        for name in ws_new:
            o, n = as_rows(ws_old[name]), as_rows(ws_new[name])
            if o != n:
                ops = []
                for tag, i1, i2, j1, j2 in difflib.SequenceMatcher(None, o, n, autojunk=False).get_opcodes():
                    if tag != "equal":
                        ops.append({"op": tag, "rows_before": [i1 + 1, i2], "rows_after": [j1 + 1, j2],
                                    "first_before": list(o[i1])[:3] if i2 > i1 else None,
                                    "first_after": list(n[j1])[:3] if j2 > j1 else None})
                sheets[name] = {"rows_before": len(o), "rows_after": len(n), "changes": ops}
        log.note("workbook_rows", sheets)
        for k, v in sheets.items():
            print(f"      {k}: rows {v['rows_before']} -> {v['rows_after']}; " + "; ".join(
                f"{c['op']} {c['rows_before'][0]}-{c['rows_before'][1]} -> {c['rows_after'][0]}-{c['rows_after'][1]}"
                for c in v["changes"]))
        check("workbook_sheets_changed", sorted(sheets) == ["README", "S16", "S32", "S6", "S7"], sorted(sheets))
        # each changed table's sheet holds its CSV, header and rows, text for text
        same = {}
        for name, csv_ in (("S6", "thresholds_pathogenic.csv"), ("S7", "thresholds_benign.csv"),
                           ("S32", "twelve_gene_refit.csv")):
            t = table(OUT / "supplementary" / csv_)
            rows_ = as_rows(ws_new[name])
            want = [tuple(t.columns)] + [tuple(v for v in r if v != "") for r in t.itertuples(index=False)]
            start = rows_.index(tuple(t.columns))
            same[name] = rows_[start:start + len(want)] == want
        t = table(OUT / "supplementary/practical_metrics.by_threshold.csv")
        same["S16"] = tuple(t.columns) in as_rows(ws_new["S16"])
        check("workbook_tables_are_the_csvs", all(same.values()), same)

        # ---------------------------------------------------------------- E7.2 Figure 2
        # every figure is drawn and measured under the style the stage draws in: its font
        # is found through rcParams, and another font measures text differently
        print("E7.2 Figure 2")
        style = ExitStack()
        style.enter_context(matplotlib.rc_context())
        F._style()
        tmp = Path(tempfile.mkdtemp())
        kept, save = [], L._save
        L._save = lambda fig, folder, stem, c: (kept.append(fig), save(fig, folder, stem, c))[1]
        try:
            with matplotlib.rc_context():
                L.figure_heldout_external(dict(items["heldout_and_external_lr"]), tmp / "main", cfg)
                fig2 = kept[-1]
                L.redraw_figure(items["local_lr_curves"], tmp / "supplementary", cfg)
                s1 = kept[-1]
                s1_before_item = dict(items["local_lr_curves"], options={
                    k: v for k, v in items["local_lr_curves"]["options"].items() if k != "caption_gap_mm"})
                L.redraw_figure(s1_before_item, tmp / "before", cfg)
                s1_before = kept[-1]
        finally:
            L._save = save
        for stem, folder in (("heldout_and_external_lr", "main"), ("local_lr_curves", "supplementary")):
            for suffix in (".pdf", ".png"):
                stage = (PHASE1 / OUT / folder / f"{stem}{suffix}").read_bytes()
                here = (tmp / folder / f"{stem}{suffix}").read_bytes()
                check(f"{stem}{suffix}_is_this_drawing", stage == here, _sha(stage))
        folds = K.read_back("reports/evidence/evidence_thresholds_logo_folds.csv")
        folds = folds[(folds.side == "pp3") & np.isfinite(folds.heldout_lr)]
        grey, blue = to_rgba(F.MUTED), to_rgba(F.BLUE)
        row = F.COLUMNS.index("avi")
        drawn = []
        ok_grey = ok_median = True
        for k, ax in enumerate(fig2.axes[:6]):
            tier, st = ("moderate", "strong")[k // 3], F.STRATA[k % 3]
            pts = []
            for coll in ax.collections:
                fc, ec = [tuple(c) for c in coll.get_facecolor()], [tuple(c) for c in coll.get_edgecolor()]
                for j, (x, y) in enumerate(coll.get_offsets()):
                    face = fc[j % len(fc)] if fc else None
                    pts.append((round(float(y)), float(x), face if face and face != (1, 1, 1, 1) else ec[j % len(ec)]))
            f = folds[(folds.tool == "avi") & (folds.stratum == st) & (folds.tier == tier)]
            out = f.heldout_gene.isin(K.AVI_SEEN_IN_TRAINING)
            g = sorted(x for r_, x, c in pts if r_ == row and c == grey)
            ok_grey &= g == sorted(np.clip(f.heldout_lr[out], 0.3, 1000).tolist())
            ok_grey &= not [p for p in pts if p[0] != row and p[2] == grey]
            med = [ln.get_xdata()[0] for ln in ax.lines
                   if len(ln.get_ydata()) == 2 and abs(np.mean(ln.get_ydata()) - row) < 1e-9]
            counted = np.clip(f.heldout_lr[~out], 0.3, 1000)
            ok_median &= med == ([float(np.median(counted))] if len(counted) else [])
            allf = np.clip(f.heldout_lr, 0.3, 1000)
            drawn.append({"panel": f"{'ab'[k // 3]} {F.STRATUM_LABEL[st]}", "grey": [
                f"{gname}={v:.3g}" for gname, v in zip(f.heldout_gene[out], f.heldout_lr[out])],
                "median_before": round(float(np.median(allf)), 4) if len(allf) else None,
                "median_after": round(float(med[0]), 4) if med else None})
        log.note("figure2_atlas_rows", drawn)
        for dd in drawn:
            print(f"      {dd['panel']:12s} grey {dd['grey']}  median {dd['median_before']} -> {dd['median_after']}")
        check("figure2_two_folds_in_grey_only_in_the_atlas_rows", ok_grey)
        check("figure2_medians_over_the_counted_folds", ok_median)
        leg = fig2.legends[0]
        fig2.set_dpi(300)
        fig2.canvas.draw()
        texts = leg.get_texts()
        x0 = [t.get_window_extent().x0 for t in texts]
        check("figure2_key_entry_in_the_first_column",
              [t.get_text() for t in texts][:4] == ["Held-out gene, clears the boundary", "Held-out gene, below it",
                                                     "Median of held-out genes", "Fold not counted (†)"]
              and x0[0] == x0[3] < x0[4], [t.get_text() for t in texts])
        # pixels: where the PNG may differ from the previous one
        old_png = np.asarray(Image.open(io.BytesIO(_git_show(args.base, f"phase1/{OUT}/main/heldout_and_external_lr.png"))).convert("RGB"))
        new_png = np.asarray(Image.open(PHASE1 / OUT / "main/heldout_and_external_lr.png").convert("RGB"))
        H = new_png.shape[0]
        allowed = np.zeros(new_png.shape[:2], dtype=bool)
        for ax in fig2.axes[:6]:
            bb = ax.get_window_extent()
            lo_y, hi_y = sorted(ax.transData.transform([[1, row - 0.5], [1, row + 0.5]])[:, 1])
            allowed[int(H - hi_y) - 1:int(H - lo_y) + 2, int(bb.x0) - 1:int(bb.x1) + 2] = True
        # the key: every row below the lowest text of (c), its axis label
        r_ = fig2.canvas.get_renderer()
        key_from = int(H - min(ax.get_tightbbox(r_).y0 for ax in fig2.axes[6:9]))
        allowed[key_from:, :] = True
        differ = (old_png != new_png).any(axis=2) if old_png.shape == new_png.shape else None
        outside = int((differ & ~allowed).sum()) if differ is not None else -1
        check("figure2_pixels_differ_only_in_the_atlas_rows_and_the_key",
              differ is not None and outside == 0,
              {"png": list(new_png.shape), "differing": int(differ.sum()) if differ is not None else None,
               "outside": outside, "key_from_row": key_from})
        log.note("figure2_png_shapes", {"before": list(old_png.shape), "after": list(new_png.shape)})

        # ---------------------------------------------------------------- E7.3 S1
        print("E7.3 Supplementary Figure S1")
        box = re.search(rb"/MediaBox \[\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*\]",
                        (PHASE1 / OUT / "supplementary/local_lr_curves.pdf").read_bytes())
        mm = [round(float(box.group(3)) * 25.4 / 72, 3), round(float(box.group(4)) * 25.4 / 72, 3)]
        check("s1_pdf_mm", mm == [175.0, 240.0], mm)
        dpi = cfg["artwork"]["dpi"]

        def boxes_gap(fig):
            fig.set_dpi(dpi)
            fig.canvas.draw()
            r = fig.canvas.get_renderer()
            cap = min(fig.texts, key=lambda t: t.get_window_extent(r).y0)
            g = (min(ax.get_tightbbox(r).y0 for ax in fig.axes) - cap.get_window_extent(r).y1) / dpi * 25.4
            fig.set_dpi(100)
            return g

        def ink_gap(png: np.ndarray) -> float:
            rows_ = np.where((png.min(axis=2) < 200).any(axis=1))[0]
            breaks = np.where(np.diff(rows_) > 1)[0]
            return (rows_[breaks[-1] + 1] - rows_[breaks[-1]] - 1) / dpi * 25.4

        old_s1 = np.asarray(Image.open(io.BytesIO(_git_show(args.base, f"phase1/{OUT}/supplementary/local_lr_curves.png"))).convert("RGB"))
        new_s1 = np.asarray(Image.open(PHASE1 / OUT / "supplementary/local_lr_curves.png").convert("RGB"))
        gaps = {"stage_measure_font_boxes_bitmap_and_pdf": [round(L.caption_gap_mm(s1_before, dpi), 3),
                                                             round(L.caption_gap_mm(s1, dpi), 3)],
                "matplotlib_text_boxes_at_300_dpi": [round(boxes_gap(s1_before), 3), round(boxes_gap(s1), 3)],
                "ink_in_the_png": [round(ink_gap(old_s1), 3), round(ink_gap(new_s1), 3)]}
        check("s1_caption_gap_at_least_2_mm", min(v[1] for v in gaps.values()) >= 2 - 1e-9, gaps)
        heights = [sorted({round(ax.get_position().height * 240, 4) for ax in f.axes}) for f in (s1_before, s1)]
        check("s1_panels_one_height_shorter_in_proportion",
              len(heights[0]) == len(heights[1]) == 1 and s1.subplotpars.top == s1_before.subplotpars.top,
              {"panel_mm": [heights[0][0], heights[1][0]], "bottom": [s1_before.subplotpars.bottom, s1.subplotpars.bottom],
               "top": s1.subplotpars.top, "png_px": [list(old_s1.shape[:2]), list(new_s1.shape[:2])]})
        sizes = [sorted({t.get_fontsize() for t in f.findobj(matplotlib.text.Text)
                         if t.get_visible() and t.get_text().strip()}) for f in (s1_before, s1)]
        check("s1_type_as_before", sizes[0] == sizes[1], sizes[1])

        # ---------------------------------------------------------------- E7.4 the old name
        print("E7.4 the old name")
        hits, scanned = [], 0
        for p in [Path("config/layout.yaml"), Path("src/evid_layout.py")] + sorted(OUT.rglob("*.csv")):
            scanned += 1
            if OLD_NAME.search(p.read_text()):
                hits.append(str(p))
        with zipfile.ZipFile(OUT / "supplementary_tables.xlsx") as z:
            for n in z.namelist():
                scanned += 1
                if OLD_NAME.search(z.read(n).decode()):
                    hits.append(f"supplementary_tables.xlsx:{n}")
        fig_texts = 0
        figs = []
        kept.clear()
        L._save = lambda fig, folder, stem, c: (kept.append(fig), save(fig, folder, stem, c))[1]
        try:
            with matplotlib.rc_context():
                for it in cfg["main"] + cfg["supplementary"]:
                    if it["kind"] == "figure":
                        L.MAKERS[it["make"]](dict(it), tmp / "all", cfg)
                        figs.append((it["content"], kept[-1]))
        finally:
            L._save = save
        for content, fig in figs:
            for t in fig.findobj(matplotlib.text.Text):
                fig_texts += 1
                if OLD_NAME.search(t.get_text()):
                    hits.append(f"{content}: '{t.get_text()}'")
        check("no_interval_calibration", not hits,
              {"files_and_workbook_parts": scanned, "figures": len(figs), "figure_texts": fig_texts, "hits": hits})
        plt.close("all")
        style.close()

        # ---------------------------------------------------------------- E7.5
        print("E7.5 the two texts")
        rel = "reports/evidence/predictor_training_provenance.csv"
        old, new = base_table(rel), table(rel)
        diff = [(old.at[i, "score_column"], c) for i in range(len(old)) for c in old.columns if old.at[i, c] != new.at[i, c]]
        check("provenance_one_cell", diff == [("fusion_enet", "training_data")] and old.shape == new.shape, diff)
        i = int(np.flatnonzero(old.score_column == "fusion_enet")[0])
        log.note("provenance_cell", {"before": old.at[i, "training_data"], "after": new.at[i, "training_data"]})
        rel = "reports/evidence/supplement/supplement_tables_manifest.csv"
        old, new = base_table(rel), table(rel)
        diff = [(old.at[i, "Table"], c) for i in range(len(old)) for c in old.columns if old.at[i, c] != new.at[i, c]]
        j = int(np.flatnonzero(old.Table == "S3")[0])
        files_ = old.at[j, "Source files"].split("; ")
        pa, pb = old.at[j, "Source sha256"].split("; "), new.at[j, "Source sha256"].split("; ")
        entries = [(f_, x, y) for f_, x, y in zip(files_, pa, pb) if x != y]
        prov_sha = _sha((PHASE1 / "reports/evidence/predictor_training_provenance.csv").read_bytes())
        old_prov_sha = _sha(_git_show(args.base, "phase1/reports/evidence/predictor_training_provenance.csv"))
        check("supplement_manifest_one_entry", diff == [("S3", "Source sha256")] and entries ==
              [("phase1/reports/evidence/predictor_training_provenance.csv", old_prov_sha, prov_sha)], entries)
        readme = "data/assay_labels/superseded/README.md"
        a = _git_show(args.base, f"phase1/{readme}").decode().splitlines()
        b = (PHASE1 / readme).read_text().splitlines()
        lines = [ln for ln in difflib.unified_diff(a, b, lineterm="", n=0) if not ln.startswith(("---", "+++"))]
        log.note("readme_diff", lines)
        check("readme_names_the_stage", "No stage reads" not in "\n".join(b) and "Stage 30" in "\n".join(b), lines)

        print(f"\n{len(failed)} check(s) failed" + (f": {failed}" if failed else ""))
        log.note("failed", failed)
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
