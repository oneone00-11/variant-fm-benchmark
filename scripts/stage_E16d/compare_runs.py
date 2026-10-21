"""E16d: what relabelling PALB2 and BARD1 changed in the evidence outputs, read from
two commits (or a commit and the working tree): the in-sample thresholds, the interval
likelihood ratio at them, the held-out folds that reach each tier, and the counts.
Commit 2 of E16d adds this script without running it; commit 3 runs it, comparing its
own outputs with commit 2's.

Writes
  phase1/reports/evidence/label_change_effect.csv: one row per predictor x stratum x
    tier x side (PP3, BP4), from evidence_thresholds.csv, the interval-LR curves and
    evidence_thresholds_logo_folds.csv:
      threshold_old / _new        the in-sample threshold;
      lr_at_threshold_old / _new  the curve's local LR at the grid point the threshold
                                  is (interval_lr_<tool>_<stratum>.csv);
      lr_bound_old / _new         the bound the tier is read from there (the gene
                                  bootstrap's lower bound for PP3, upper for BP4);
      folds_reaching_old / _new   held-out folds whose six training genes give a
                                  threshold for the tier;
      reached_old / _new, tier_change   in sample: gained, lost or same;
      n_pos / n_neg old and new, and n_labelled_between: the labelled variants of the
      new analysis set whose score lies between the old and the new threshold (NA for
      the fusion, which the stages compute and the set does not carry);
  docs/label-change-effect.md: every cell whose tier changed, in sample and held out;
    the count changes by gene and stratum; and the ten cells whose threshold moved
    across the most labelled variants.

It also holds the new counts to the audit's prediction (E16d.0 rule 3): every cell of
label_audit_impact.csv as the old commit has it (the counts each relabelled gene
would have under its closest rule) must equal the new depth_counts.csv / set_counts.csv
cell for cell, and the seven genes' totals must move by the predicted amounts. The
comparison is printed in the record page; any difference exits non-zero, after the
files are written, because it would mean something besides the labels changed.

Run from the repository root:
    .venv/bin/python scripts/stage_E16d/compare_runs.py --old <commit> [--new <commit>]
(without --new, the working tree is the new side)
"""
from __future__ import annotations

import argparse
import io
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel  # noqa: E402

sys.path.insert(0, str(REPO / "phase1"))
from src import evid_common as K  # noqa: E402

REP = "phase1/reports/evidence"
OUT_CSV = REPO / REP / "label_change_effect.csv"
OUT_MD = REPO / "docs" / "label-change-effect.md"
SET = "phase1/data/evidence/analysis_set_v1.parquet"
TIERS = ["supporting", "moderate", "strong", "very_strong"]


def _bytes(ref: str | None, path: str) -> bytes | None:
    if ref is None:
        p = REPO / path
        return p.read_bytes() if p.exists() else None
    r = subprocess.run(["git", "-C", str(REPO), "show", f"{ref}:{path}"], capture_output=True)
    return r.stdout if r.returncode == 0 else None


def _csv(ref: str | None, path: str) -> pd.DataFrame | None:
    b = _bytes(ref, path)
    return None if b is None else pd.read_csv(io.BytesIO(b))


def _folds(df: pd.DataFrame) -> dict:
    ok = df[np.isfinite(pd.to_numeric(df["threshold_from_6_genes"], errors="coerce"))]
    return ok.groupby(["tool", "stratum", "tier", "side"]).size().to_dict()


def _at(curve: pd.DataFrame | None, thr: float, side: str) -> tuple[float, float]:
    if curve is None or not np.isfinite(thr):
        return np.nan, np.nan
    hit = curve[curve["score"] == thr]
    if hit.empty:
        hit = curve[np.isclose(curve["score"], thr, rtol=1e-12, atol=0)]
    if hit.empty:
        return np.nan, np.nan
    r = hit.iloc[0]
    return float(r["local_lr"]), float(r["lr_lo_gene_boot" if side == "pp3" else "lr_hi_gene_boot"])


def _between(aset: pd.DataFrame, tool: str, stratum: str, a: float, b: float) -> float:
    if tool not in aset.columns or not (np.isfinite(a) and np.isfinite(b)):
        return np.nan
    sub = K.stratum_frame(aset, stratum).dropna(subset=["y_assay", tool])
    lo, hi = min(a, b), max(a, b)
    return int(((sub[tool] >= lo) & (sub[tool] < hi)).sum())


def compare(old: str, new: str | None) -> pd.DataFrame:
    t_old = _csv(old, f"{REP}/evidence_thresholds.csv")
    t_new = _csv(new, f"{REP}/evidence_thresholds.csv")
    f_old = _folds(_csv(old, f"{REP}/evidence_thresholds_logo_folds.csv"))
    f_new = _folds(_csv(new, f"{REP}/evidence_thresholds_logo_folds.csv"))
    aset = pd.read_parquet(io.BytesIO(_bytes(new, SET)))
    keys = ["tool", "stratum", "tier"]
    m = t_old.merge(t_new, on=keys, suffixes=("_old", "_new"), validate="one_to_one")
    curves = {}
    rows = []
    for r in m.to_dict("records"):
        tool, stratum, tier = r["tool"], r["stratum"], r["tier"]
        for ref, tag in ((old, "old"), (new, "new")):
            if (tool, stratum, tag) not in curves:
                curves[(tool, stratum, tag)] = _csv(ref, f"{REP}/interval_lr_{tool}_{stratum}.csv")
        for side in ("pp3", "bp4"):
            told = float(r[f"{side}_threshold_insample_old"])
            tnew = float(r[f"{side}_threshold_insample_new"])
            rold, rnew = bool(r[f"{side}_threshold_reachable_old"] is True or
                              str(r[f"{side}_threshold_reachable_old"]) == "True"), \
                bool(r[f"{side}_threshold_reachable_new"] is True or
                     str(r[f"{side}_threshold_reachable_new"]) == "True")
            lo_old, bd_old = _at(curves[(tool, stratum, "old")], told, side)
            lo_new, bd_new = _at(curves[(tool, stratum, "new")], tnew, side)
            rows.append({
                "tool": tool, "stratum": stratum, "tier": tier, "side": side,
                "status_old": r["status_old"], "status_new": r["status_new"],
                "cut_lr": r["path_cut_lr_new"] if side == "pp3" else r["ben_cut_lr_new"],
                "threshold_old": told, "threshold_new": tnew,
                "lr_at_threshold_old": lo_old, "lr_at_threshold_new": lo_new,
                "lr_bound_old": bd_old, "lr_bound_new": bd_new,
                "folds_reaching_old": f_old.get((tool, stratum, tier, side), 0),
                "folds_reaching_new": f_new.get((tool, stratum, tier, side), 0),
                "reached_old": rold, "reached_new": rnew,
                "tier_change": ("gained" if rnew and not rold else
                                "lost" if rold and not rnew else "same"),
                "n_pos_old": r["n_pos_old"], "n_pos_new": r["n_pos_new"],
                "n_neg_old": r["n_neg_old"], "n_neg_new": r["n_neg_new"],
                "n_labelled_between": _between(aset, tool, stratum, told, tnew),
            })
    out = pd.DataFrame(rows)
    out["n_labelled_between"] = out["n_labelled_between"].astype("Int64")
    return out


def _f(x) -> str:
    if x is None or (isinstance(x, float) and math.isnan(x)) or x is pd.NA:
        return "n/a"
    if isinstance(x, float):
        return f"{x:.6g}"
    return str(x)


def counts(old: str, new: str | None) -> pd.DataFrame:
    """The analysis set's counts by gene and stratum that changed, old and new."""
    d_old, d_new = _csv(old, f"{REP}/depth_counts.csv"), _csv(new, f"{REP}/depth_counts.csv")
    s_old, s_new = _csv(old, f"{REP}/set_counts.csv"), _csv(new, f"{REP}/set_counts.csv")
    rows = []
    for gene in sorted(set(d_new["gene"])):
        for b in ("3-10", "11-50", "3-50"):
            o = d_old[(d_old["gene"] == gene) & (d_old["depth_nt"] == b)]
            n = d_new[(d_new["gene"] == gene) & (d_new["depth_nt"] == b)]
            for col in ("labelled", "damaging", "normal"):
                if len(o) and len(n) and int(o[col].iloc[0]) != int(n[col].iloc[0]):
                    rows.append({"gene": gene, "stratum": b, "count": col,
                                 "old": int(o[col].iloc[0]), "new": int(n[col].iloc[0])})
    for gene in sorted(set(s_new["gene"])):
        o = s_old[(s_old["gene"] == gene) & (s_old["stratum"] == "pm12")]
        n = s_new[(s_new["gene"] == gene) & (s_new["stratum"] == "pm12")]
        for col in ("n_labelled", "n_damaging", "n_normal"):
            if int(o[col].sum()) != int(n[col].sum()):
                rows.append({"gene": gene, "stratum": "1-2", "count": col.removeprefix("n_"),
                             "old": int(o[col].sum()), "new": int(n[col].sum())})
    return pd.DataFrame(rows, columns=["gene", "stratum", "count", "old", "new"])


BINS = {"pm12": None, "s3_10": "3-10", "s11_50": "11-50", "s3_50": "3-50"}


def _counts_at(ref: str | None, gene: str, stratum: str) -> dict:
    if stratum == "pm12":
        s = _csv(ref, f"{REP}/set_counts.csv")
        s = s[(s["gene"] == gene) & (s["stratum"] == "pm12")]
        return {"n": int(s["n"].sum()), "labelled": int(s["n_labelled"].sum()),
                "damaging": int(s["n_damaging"].sum()), "normal": int(s["n_normal"].sum())}
    d = _csv(ref, f"{REP}/depth_counts.csv")
    r = d[(d["gene"] == gene) & (d["depth_nt"] == BINS[stratum])].iloc[0]
    return {"n": int(r["n"]), "labelled": int(r["labelled"]), "damaging": int(r["damaging"]),
            "normal": int(r["normal"])}


def predictions(old: str, new: str | None) -> pd.DataFrame:
    """Each cell of the old commit's label_audit_impact.csv against the new counts, and
    the seven genes' totals (old total plus the predicted change) against the new ones."""
    imp = _csv(old, f"{REP}/label_audit_impact.csv")
    rows = []
    for r in imp.itertuples():
        got = _counts_at(new, r.gene, r.stratum)
        for k in ("labelled", "damaging", "normal", "excluded"):
            want = getattr(r, f"n_{k}_rule")
            actual = got["n"] - got["labelled"] if k == "excluded" else got[k]
            rows.append({"gene": r.gene, "stratum": r.stratum, "count": k,
                         "now": int(getattr(r, f"n_{k}_now")), "predicted": int(want),
                         "actual": actual})
    for st in ("s3_10", "s11_50", "s3_50"):
        before, after = _counts_at(old, "all", st), _counts_at(new, "all", st)
        for k in ("labelled", "damaging", "normal", "excluded"):
            sub = imp[imp["stratum"] == st]
            delta = int((sub[f"n_{k}_rule"] - sub[f"n_{k}_now"]).sum())
            b = before["n"] - before["labelled"] if k == "excluded" else before[k]
            a = after["n"] - after["labelled"] if k == "excluded" else after[k]
            rows.append({"gene": "all seven", "stratum": st, "count": k, "now": b,
                         "predicted": b + delta, "actual": a})
    out = pd.DataFrame(rows)
    out["met"] = out["predicted"] == out["actual"]
    return out


def write_md(eff: pd.DataFrame, cnt: pd.DataFrame, old: str, new: str | None,
             pred: pd.DataFrame) -> None:
    L = ["# What relabelling PALB2 and BARD1 changed",
         "",
         f"Written by `scripts/stage_E16d/compare_runs.py`, comparing the evidence outputs of "
         f"commit {old} (old) with {new or 'the working tree of the next commit'} (new); do not "
         "edit by hand. The only change between the two is the label files of PALB2 and BARD1, "
         "now built from the investigator-provided ranges of the score sets Supplementary Table "
         "S1 cites (stage E16d). Full rows: `phase1/reports/evidence/label_change_effect.csv`.",
         "", "## Tier changes", ""]
    tc = eff[eff["tier_change"] != "same"]
    L.append(f"In sample, {len(tc)} cell(s) gained or lost a tier "
             "(predictor x stratum x tier x side):")
    L.append("")
    if len(tc):
        L += ["| side | predictor | stratum | tier | change | threshold old -> new | LR there old -> new | held-out folds old -> new |",
              "|---|---|---|---|---|---|---|---|"]
        for r in tc.itertuples():
            L.append(f"| {r.side.upper()} | {r.tool} | {r.stratum} | {r.tier} | {r.tier_change} | "
                     f"{_f(r.threshold_old)} -> {_f(r.threshold_new)} | {_f(r.lr_at_threshold_old)} "
                     f"-> {_f(r.lr_at_threshold_new)} | {r.folds_reaching_old} -> {r.folds_reaching_new} |")
        L.append("")
    fc = eff[eff["folds_reaching_old"] != eff["folds_reaching_new"]]
    L.append(f"Held out, {len(fc)} cell(s) changed the number of folds that reach the tier:")
    L.append("")
    if len(fc):
        L += ["| side | predictor | stratum | tier | folds old -> new |", "|---|---|---|---|---|"]
        for r in fc.itertuples():
            L.append(f"| {r.side.upper()} | {r.tool} | {r.stratum} | {r.tier} | "
                     f"{r.folds_reaching_old} -> {r.folds_reaching_new} |")
        L.append("")
    L += ["## The audit's prediction", "",
          f"Every cell of `label_audit_impact.csv` at {old} (the counts the closest rule would "
          "give), against the new counts; and the seven genes' totals. "
          f"{int(pred['met'].sum())} of {len(pred)} cells met.", "",
          "| gene | stratum | count | now | predicted | actual | met |", "|---|---|---|---|---|---|---|"]
    L += [f"| {r.gene} | {r.stratum} | {r.count} | {r.now:,} | {r.predicted:,} | {r.actual:,} | "
          f"{'yes' if r.met else '**no**'} |" for r in pred.itertuples()]
    L += ["", "## Count changes", "",
          "The analysis set's counts that changed, by gene and stratum "
          "(`depth_counts.csv`, `set_counts.csv`):", ""]
    if len(cnt):
        L += ["| gene | stratum | count | old | new |", "|---|---|---|---|---|"]
        L += [f"| {r.gene} | {r.stratum} | {r.count} | {r.old:,} | {r.new:,} |" for r in cnt.itertuples()]
    else:
        L.append("none")
    L += ["", "## The ten largest threshold moves", "",
          "Ranked by the labelled variants of the new analysis set that lie between the old and "
          "the new threshold (the fusion, which the set does not carry, is left out):", ""]
    mv = eff.dropna(subset=["n_labelled_between"])
    mv = mv[mv["n_labelled_between"] > 0].sort_values(
        ["n_labelled_between", "tool", "stratum", "tier", "side"], ascending=[False, True, True, True, True]).head(10)
    if len(mv):
        L += ["| side | predictor | stratum | tier | threshold old -> new | variants between |",
              "|---|---|---|---|---|---|"]
        L += [f"| {r.side.upper()} | {r.tool} | {r.stratum} | {r.tier} | {_f(r.threshold_old)} -> "
              f"{_f(r.threshold_new)} | {int(r.n_labelled_between)} |" for r in mv.itertuples()]
    else:
        L.append("No threshold moved across a labelled variant.")
    L.append("")
    OUT_MD.write_text("\n".join(L))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--old", required=True)
    ap.add_argument("--new", default=None)
    args = ap.parse_args()
    with RunLog(__file__) as log:
        old = subprocess.run(["git", "-C", str(REPO), "rev-parse", args.old], check=True,
                             capture_output=True, text=True).stdout.strip()
        new = (subprocess.run(["git", "-C", str(REPO), "rev-parse", args.new], check=True,
                              capture_output=True, text=True).stdout.strip() if args.new else None)
        log.note("old", old)
        log.note("new", new or "working tree")
        eff = compare(old, new)
        eff.to_csv(OUT_CSV, index=False)
        log.output(OUT_CSV)
        cnt = counts(old, new)
        pred = predictions(old, new)
        write_md(eff, cnt, old[:7], new[:7] if new else None, pred)
        log.output(OUT_MD)
        log.note("prediction_cells_met", f"{int(pred['met'].sum())}/{len(pred)}")
        log.note("tier_changes", int((eff["tier_change"] != "same").sum()))
        print(f"{len(eff)} rows; {int((eff['tier_change'] != 'same').sum())} tier change(s); "
              f"{len(cnt)} count change(s); prediction {int(pred['met'].sum())}/{len(pred)} "
              f"cells met; wrote {rel(OUT_CSV)}, {rel(OUT_MD)}")
        if not pred["met"].all():
            raise SystemExit("the counts differ from the audit's prediction; stopped (E16d.6)")


if __name__ == "__main__":
    main()
