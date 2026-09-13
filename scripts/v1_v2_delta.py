"""frozen-matrix-v1 versus v2: value-by-value comparison of the pipeline outputs.

Three views, written to one report:

  A. the conclusion-bearing quantities (H1 leaderboard and delta-rho, sign-flip p,
     leave-two-genes-out, no-offset drop; H2 ablation; H3 delta-Brier, Murphy terms,
     high-confidence yield; H4 likelihood ratios, evidence tiers and operating points;
     TP53; fusion weights; tie structure of the splice subset);
  B. every number printed in the manuscript: each token is traced to the v1 output
     cell(s) it matches (the number checker's own tokeniser and tolerance) and the
     same cell in the v2 run is read back;
  C. every numeric cell in every result table, v1 against v2, sorted by |delta|.

    python scripts/v1_v2_delta.py --v2 phase1/reports/phase1_v2 \
        --manuscript ~/Desktop/calibration_draft.docx --out docs/v1-v2-delta.md
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
V1 = REPO / "phase1" / "reports" / "phase1"


def load_checker():
    spec = importlib.util.spec_from_file_location("cmn", REPO / "phase1" / "src" / "check_manuscript_numbers.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def read(d: Path, name: str) -> pd.DataFrame | None:
    p = d / name
    return pd.read_csv(p) if p.exists() else None


def ci(s) -> tuple[float, float] | None:
    m = re.findall(r"[-+]?\d*\.?\d+(?:e[-+]?\d+)?", str(s))
    return (float(m[0]), float(m[1])) if len(m) >= 2 else None


def excl0(c) -> str:
    return "n/a" if c is None else ("excludes 0" if (c[0] > 0 or c[1] < 0) else "spans 0")


def fmt(v) -> str:
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return "—"
    if isinstance(v, float):
        return f"{v:.4f}" if abs(v) < 100 else f"{v:.2f}"
    return str(v)


class Delta:
    def __init__(self):
        self.rows: list[dict] = []

    def add(self, name, where, v1, v2, note1="", note2=""):
        try:
            d = float(v2) - float(v1)
        except (TypeError, ValueError):
            d = float("nan")
        flag = ""
        if note1 and note2 and note1 != note2:
            flag = f"{note1} → {note2}"
        self.rows.append({"quantity": name, "paper": where, "v1": v1, "v2": v2, "abs_change": abs(d) if not np.isnan(d) else float("nan"),
                          "change": d, "flag": flag or (note1 if note1 == note2 else "")})

    def table(self, sort=True) -> str:
        rows = sorted(self.rows, key=lambda r: -(r["abs_change"] if not np.isnan(r["abs_change"]) else -1)) if sort else self.rows
        out = ["| quantity | paper location | v1 | v2 | Δ (v2−v1) | interval / tier |", "|---|---|---|---|---|---|"]
        for r in rows:
            out.append(f"| {r['quantity']} | {r['paper']} | {fmt(r['v1'])} | {fmt(r['v2'])} | "
                       f"{'' if np.isnan(r['change']) else f'{r['change']:+.4f}'} | {r['flag']} |")
        return "\n".join(out)


def curated(v1d: Path, v2d: Path) -> tuple[str, list[str]]:
    D = Delta()
    highlights: list[str] = []
    # ---- H1
    a, b = read(v1d, "phase2_leaderboard.csv"), read(v2d, "phase2_leaderboard.csv")
    if a is not None and b is not None:
        b = b.set_index("model")
        for r in a.itertuples():
            if r.model in b.index:
                D.add(f"pooled ρ — {r.model}", "Table 1 / Abstract / 3.1", r.pooled_rho, b.loc[r.model, "pooled_rho"],
                      f"CI {r.ci95}", f"CI {b.loc[r.model, 'ci95']}")
        order1 = list(a["model"]); order2 = list(b.sort_values("pooled_rho", ascending=False).index)
        highlights.append(f"Leaderboard order v1: {' > '.join(order1)}")
        highlights.append(f"Leaderboard order v2: {' > '.join(order2)}")
    a, b = read(v1d, "phase2_H1_stratified.csv"), read(v2d, "phase2_H1_stratified.csv")
    if a is not None and b is not None:
        b = b.set_index("stratum")
        for r in a.itertuples():
            if r.stratum in b.index:
                D.add(f"Δρ fusion − best single — {r.stratum}", "Abstract / 3.1 (H1)", r.delta, b.loc[r.stratum, "delta"],
                      f"{excl0(ci(r.ci95))} {r.ci95} [{r.H1}]", f"{excl0(ci(b.loc[r.stratum, 'ci95']))} {b.loc[r.stratum, 'ci95']} [{b.loc[r.stratum, 'H1']}]")
                D.add(f"best single ρ — {r.stratum}", "3.1", r.single_rho, b.loc[r.stratum, "single_rho"])
                D.add(f"fusion ρ — {r.stratum}", "3.1", r.fusion_rho, b.loc[r.stratum, "fusion_rho"])
    a, b = read(v1d, "phase8_sign_flip_exact.csv"), read(v2d, "phase8_sign_flip_exact.csv")
    if a is not None and b is not None:
        for r in a.itertuples():
            m = b[(b.analysis == r.analysis) & (b.condition == r.condition)]
            if len(m):
                D.add(f"sign-flip exact p — {r.analysis} {r.condition}", "3.1 / 3.3 (Table S16)", r.p_exact, float(m.p_exact.iloc[0]),
                      "p<0.05" if r.p_exact < 0.05 else "ns", "p<0.05" if float(m.p_exact.iloc[0]) < 0.05 else "ns")
    a, b = read(v1d, "phase8_leave_two_genes_out.csv"), read(v2d, "phase8_leave_two_genes_out.csv")
    if a is not None and b is not None:
        for stat, f in (("median", np.median), ("min", np.min), ("max", np.max)):
            D.add(f"leave-two-genes-out Δρ — {stat}", "3.1", float(f(a.delta_rho)), float(f(b.delta_rho)),
                  "all positive" if (a.delta_rho > 0).all() else "some ≤ 0", "all positive" if (b.delta_rho > 0).all() else "some ≤ 0")
    a, b = read(v1d, "phase2_no_offset_drop.csv"), read(v2d, "phase2_no_offset_drop.csv")
    if a is not None and b is not None:
        b = b.set_index("quantity")
        for r in a.itertuples():
            if r.quantity in b.index:
                try:
                    D.add(f"no-offset drop — {r.quantity}", "3.1", float(r.value), float(b.loc[r.quantity, "value"]))
                except (TypeError, ValueError):
                    pass
    # ---- H2
    a, b = read(v1d, "phase2_H2_ablation.csv"), read(v2d, "phase2_H2_ablation.csv")
    if a is not None and b is not None:
        for r in a.itertuples():
            m = b[(b.model == r.model) & (b.ablation == r.ablation)]
            if len(m):
                D.add(f"H2 ablation Δρ — {r.model} {r.ablation}", "Table 2 / 3.2", r.delta_full_minus_ablated,
                      float(m.delta_full_minus_ablated.iloc[0]), f"{excl0(ci(r.ci95))}", f"{excl0(ci(m.ci95.iloc[0]))}")
    # ---- H3
    a, b = read(v1d, "phase3_H3_headline.csv"), read(v2d, "phase3_H3_headline.csv")
    if a is not None and b is not None:
        for r in a.itertuples():
            m = b[(b.set == r.set) & (b.calib == r.calib)]
            if len(m):
                m = m.iloc[0]
                D.add(f"ΔBrier best single − fusion — {r.set} {r.calib}", "3.3 (H3) / Fig 3", r.dBrier, m.dBrier,
                      f"{excl0(ci(r.dBrier_ci_geneclust))} {r.dBrier_ci_geneclust}", f"{excl0(ci(m.dBrier_ci_geneclust))} {m.dBrier_ci_geneclust}")
                D.add(f"ΔECE — {r.set} {r.calib}", "3.3 / S3b", r.dECE, m.dECE, excl0(ci(r.dECE_ci)), excl0(ci(m.dECE_ci)))
                D.add(f"Δ high-confidence yield — {r.set} {r.calib}", "3.3 / Fig 3", r.dYield, m.dYield, excl0(ci(r.dYield_ci)), excl0(ci(m.dYield_ci)))
    a, b = read(v1d, "phase3_calibration_summary.csv"), read(v2d, "phase3_calibration_summary.csv")
    if a is not None and b is not None:
        sel = (a.set == "y_assay/BRCA1_included") & (a.calib == "isotonic")
        for r in a[sel].itertuples():
            m = b[(b.set == r.set) & (b.calib == r.calib) & (b.model == r.model)]
            if len(m):
                m = m.iloc[0]
                for col, where in (("ECE", "Table 3"), ("Brier", "Table 3"), ("actionable_frac", "Table 3 (high-confidence %)"), ("actionable_acc", "Table 3 (accuracy %)")):
                    D.add(f"{col} — {r.model} (primary condition)", where, getattr(r, col), getattr(m, col))
    a, b = read(v1d, "phase8_murphy_decomposition.csv"), read(v2d, "phase8_murphy_decomposition.csv")
    if a is not None and b is not None:
        r = a[(a.condition == "isotonic/y_assay/BRCA1_included") & (a.binning == "width10")]
        m = b[(b.condition == "isotonic/y_assay/BRCA1_included") & (b.binning == "width10")]
        if len(r) and len(m):
            r, m = r.iloc[0], m.iloc[0]
            for col in ("brier_fusion", "brier_single", "rel_fusion", "rel_single", "res_fusion", "res_single", "unc"):
                D.add(f"Murphy — {col}", "3.3 (Murphy decomposition) / S15", r[col], m[col])
            D.add("Murphy — ΔREL", "3.3", r.d_rel, m.d_rel, excl0((r.d_rel_lo, r.d_rel_hi)), excl0((m.d_rel_lo, m.d_rel_hi)))
            D.add("Murphy — ΔRES", "3.3", r.d_res, m.d_res, excl0((r.d_res_lo, r.d_res_hi)), excl0((m.d_res_lo, m.d_res_hi)))
    # ---- H4
    a, b = read(v1d, "phase5_likelihood_ratios.csv"), read(v2d, "phase5_likelihood_ratios.csv")
    if a is not None and b is not None:
        for cond in ("y_assay/BRCA1_included", "y_assay/BRCA1_excluded"):
            for cal in ("raw", "calibrated_isotonic"):
                sub = a[(a.condition == cond) & (a.calibration == cal)]
                for r in sub.itertuples():
                    m = b[(b.condition == cond) & (b.calibration == cal) & (b.object == r.object)]
                    if len(m):
                        m = m.iloc[0]
                        D.add(f"LR+ at 95% spec — {r.object} ({cond}, {cal})", "3.3 (H4) / S11", r.lr_plus, m.lr_plus,
                              f"{r.acmg_tier}", f"{m.acmg_tier}")
    a, b = read(v1d, "phase5_lr_headline.csv"), read(v2d, "phase5_lr_headline.csv")
    if a is not None and b is not None:
        for r in a.itertuples():
            m = b[(b.condition == r.condition) & (b.calibration == r.calibration)]
            if len(m):
                m = m.iloc[0]
                D.add(f"ΔLR+ fusion − best single — {r.condition} {r.calibration}", "3.3 (H4)", r.delta_lr_plus, m.delta_lr_plus,
                      f"{excl0((r.lo, r.hi))}; {r.tier_fusion} vs {r.tier_best_single}", f"{excl0((m.lo, m.hi))}; {m.tier_fusion} vs {m.tier_best_single}")
    a, b = read(v1d, "phase8_lr_plus_operating_points.csv"), read(v2d, "phase8_lr_plus_operating_points.csv")
    if a is not None and b is not None:
        sub = a[a.condition == "y_assay/BRCA1_included"]
        for r in sub.itertuples():
            m = b[(b.condition == r.condition) & (b.object == r.object) & (np.isclose(b.target_spec, r.target_spec))]
            if len(m):
                m = m.iloc[0]
                D.add(f"LR+ at {r.target_spec:.3f} spec — {r.object}", "3.3 (operating points) / S19", r.lr_plus, m.lr_plus, r.acmg_tier, m.acmg_tier)
    # ---- TP53
    a, b = read(v1d, "phase4_tp53_external.csv"), read(v2d, "phase4_tp53_external.csv")
    if a is not None and b is not None:
        key = lambda s: "fusion" if str(s).startswith("fusion") else "best_single"
        b2 = {key(r.model): r for r in b.itertuples()}
        for r in a.itertuples():
            m = b2.get(key(r.model))
            if m is not None:
                for col in ("ECE", "Brier", "actionable_frac", "actionable_acc"):
                    D.add(f"TP53 {col} — {r.model} → {m.model}", "3.4 / Table 4", getattr(r, col), getattr(m, col))
    a, b = read(v1d, "phase4_tp53_label_definitions.csv"), read(v2d, "phase4_tp53_label_definitions.csv")
    if a is not None and b is not None and len(a) == len(b):
        for i, r in a.iterrows():
            for col in a.columns:
                if col in b.columns and isinstance(r[col], (int, float)) and not isinstance(r[col], bool):
                    D.add(f"TP53 label definitions row {i} — {col}", "Table 4 / S8", r[col], b.iloc[i][col])
    # ---- weights
    a, b = read(v1d, "phase6_enet_weight_stability.csv"), read(v2d, "phase6_enet_weight_stability.csv")
    if a is not None and b is not None:
        a, b = a.set_index(a.columns[0]), b.set_index(b.columns[0])
        for f in ("alphagenome", "pangolin", "spliceai", "gpn_msa", "phastcons", "phylop", "nt", "cadd"):
            if f in a.index and f in b.index:
                D.add(f"elastic-net weight mean — {f}", "4 (weight stability) / S12", a.loc[f, "mean"], b.loc[f, "mean"])
    return D.table(), highlights


def every_manuscript_number(v1d: Path, v2d: Path, docx: Path, cmn) -> tuple[str, dict]:
    """Trace each manuscript token to v1 cells and read the v2 value of those cells."""
    def cells(d: Path) -> dict[str, pd.DataFrame]:
        out = {}
        for p in sorted(d.glob("*.csv")):
            try:
                out[p.name] = pd.read_csv(p)
            except Exception:
                pass
        return out
    c1, c2 = cells(v1d), cells(v2d)
    index = []  # (file, row, col, value)
    for name, df in c1.items():
        for col in df.columns:
            s = pd.to_numeric(df[col], errors="coerce")
            for i, v in s.items():
                if pd.notna(v):
                    index.append((name, i, col, float(v)))
    vals = np.array([x[3] for x in index])
    toks = cmn.manuscript_tokens(docx)
    rows = []
    stats = {"tokens": len(toks), "traced": 0, "changed": 0, "ambiguous": 0, "unmatched": 0}
    for t in toks:
        tol = cmn._tolerance(t["token"])
        hit = np.where(np.abs(vals - t["value"]) <= tol + 1e-12)[0]
        # also try percentage
        hitp = np.where(np.abs(vals * 100 - t["value"]) <= tol + 1e-12)[0]
        cand = [(index[i], False) for i in hit] + [(index[i], True) for i in hitp]
        distinct = len(set(round(index[i][3], 10) for i in hit) | set(round(index[i][3] * 100, 10) for i in hitp))
        if not cand:
            stats["unmatched"] += 1
            continue
        if distinct >= cmn.WEAK_MATCH_MIN:
            stats["ambiguous"] += 1
            continue
        stats["traced"] += 1
        best = None
        for (name, i, col, v1v), pct in cand:
            df2 = c2.get(name)
            if df2 is None or col not in df2.columns or i >= len(df2):
                continue
            v2v = pd.to_numeric(df2[col], errors="coerce").iloc[i]
            if pd.isna(v2v):
                continue
            v2v = float(v2v) * (100 if pct else 1)
            v1s = v1v * (100 if pct else 1)
            d = v2v - v1s
            if best is None or abs(d) > abs(best[-1]):
                best = (name, i, col, v1s, v2v, d)
        if best is None:
            continue
        name, i, col, v1s, v2v, d = best
        moved = abs(d) > tol
        if moved:
            stats["changed"] += 1
        rows.append({"uid": t["uid"], "token": t["token"], "context": t["context"], "cell": f"{name}[{i}].{col}",
                     "v1": v1s, "v2": v2v, "delta": d, "moves_printed_value": moved})
    rows.sort(key=lambda r: -abs(r["delta"]))
    out = ["| block | printed | v1 cell | v1 | v2 | Δ | printed value moves? | context |", "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        out.append(f"| {r['uid']} | {r['token']} | {r['cell']} | {fmt(r['v1'])} | {fmt(r['v2'])} | {r['delta']:+.4g} | "
                   f"{'YES' if r['moves_printed_value'] else 'no'} | {r['context'].replace('|', '/')[:110]} |")
    return "\n".join(out), stats


def all_cells(v1d: Path, v2d: Path, top: int = 150) -> tuple[str, dict]:
    rows = []
    n_cells = n_changed = 0
    missing = []
    for p in sorted(v1d.glob("*.csv")):
        q = v2d / p.name
        if not q.exists():
            missing.append(p.name)
            continue
        a, b = pd.read_csv(p), pd.read_csv(q)
        for col in a.columns:
            if col not in b.columns:
                continue
            s1, s2 = pd.to_numeric(a[col], errors="coerce"), pd.to_numeric(b[col], errors="coerce")
            n = min(len(s1), len(s2))
            for i in range(n):
                x, y = s1.iloc[i], s2.iloc[i]
                if pd.isna(x) or pd.isna(y):
                    continue
                n_cells += 1
                d = float(y) - float(x)
                if d != 0:
                    n_changed += 1
                    rows.append((abs(d), p.name, i, col, float(x), float(y), d))
    rows.sort(key=lambda r: -r[0])
    out = ["| file | row | column | v1 | v2 | Δ |", "|---|---|---|---|---|---|"]
    for _, name, i, col, x, y, d in rows[:top]:
        out.append(f"| {name} | {i} | {col} | {fmt(x)} | {fmt(y)} | {d:+.4g} |")
    return "\n".join(out), {"numeric_cells_compared": n_cells, "cells_changed": n_changed, "files_missing_in_v2": missing}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--v1", default=str(V1))
    ap.add_argument("--v2", required=True)
    ap.add_argument("--manuscript", required=True)
    ap.add_argument("--column-report", default=str(REPO / "phase1" / "data" / "frozen" / "frozen_matrix_v2_column_report.tsv"))
    ap.add_argument("--tests", default=None, help="path to a pytest output file to embed")
    ap.add_argument("--out", default="docs/v1-v2-delta.md")
    a = ap.parse_args()
    v1d, v2d = Path(a.v1), Path(a.v2)
    cmn = load_checker()
    cur, hl = curated(v1d, v2d)
    every, st = every_manuscript_number(v1d, v2d, Path(a.manuscript).expanduser(), cmn)
    cells, cst = all_cells(v1d, v2d)
    if Path(a.column_report).exists():
        cr = pd.read_csv(a.column_report, sep="\t")
        colrep = "\n".join(["| " + " | ".join(cr.columns) + " |", "|" + "---|" * len(cr.columns)]
                           + ["| " + " | ".join(fmt(v) if isinstance(v, float) else str(v) for v in r) + " |" for r in cr.itertuples(index=False)])
    else:
        colrep = "(column report not found)"
    tests = Path(a.tests).read_text() if a.tests and Path(a.tests).exists() else "(not run)"
    text = "\n".join([
        "# frozen-matrix-v1 → v2: what moves", "",
        f"v1 outputs: `{v1d}` (tracked). v2 outputs: `{v2d}`. Generated by `scripts/v1_v2_delta.py`; "
        "seeds, code and every setting identical between the two runs — only the analysis set differs "
        "(SpliceAI and Pangolin at full precision, Nucleotide Transformer re-scored under the pinned environment; "
        "TP53 scored for all ten predictors).", "",
        "## 0. Tie structure of the replaced columns", "", colrep, "",
        "## A. Conclusion-bearing quantities (sorted by |Δ|)", "", *[f"- {h}" for h in hl], "", cur, "",
        "## B. Every number printed in the manuscript, traced to its v1 cell and read back from v2", "",
        f"Tokens: {st['tokens']}; traced to a specific v1 cell: {st['traced']}; of those, the printed value moves under v2: "
        f"{st['changed']}; ambiguous (≥{cmn.WEAK_MATCH_MIN} candidate cells, not traced): {st['ambiguous']}; "
        f"no v1 cell within tolerance: {st['unmatched']}. Where a token matches several cells the one with the largest v2 change is shown.", "",
        every, "",
        "## C. Every numeric cell in every result table (top 150 by |Δ|)", "",
        f"{cst['numeric_cells_compared']:,} numeric cells compared; {cst['cells_changed']:,} differ. "
        f"Files absent from the v2 run: {cst['files_missing_in_v2'] or 'none'}.", "", cells, "",
        "## D. Test suite against the v2 outputs", "", "```", tests, "```", "",
    ])
    out = REPO / a.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text, encoding="utf8")
    print(f"wrote {out}; manuscript tokens traced {st['traced']}/{st['tokens']}, printed values that move: {st['changed']}; "
          f"cells changed {cst['cells_changed']}/{cst['numeric_cells_compared']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
