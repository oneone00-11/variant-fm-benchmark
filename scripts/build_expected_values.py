"""Generate tests/fixtures/expected_values.json from the pipeline outputs.

The guardrail tests (tests/test_reported_numbers.py, tests/test_phase8_robustness.py)
pin the numbers the manuscript prints to the pipeline outputs that produced them.
The pins used to be literals typed into the test files; this script derives every one
of them from phase1/reports/phase1/ at the manuscript's printed precision, so the
expected values are never hand-written and can be regenerated with one command
whenever the analysis set changes:

    python scripts/build_expected_values.py            # writes the fixture
    python scripts/build_expected_values.py --check    # exit 1 if the fixture is stale

Each entry records its source (file, row selector, column), the printed precision, and
the value as printed. The tests compare the live output with the fixture at half a unit
of that precision, exactly as before; nothing here changes a tolerance.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
REPORTS = REPO / "phase1" / "reports" / "phase1"
OUT = REPO / "tests" / "fixtures" / "expected_values.json"

# kind: num (float at dp), pct (x100 at dp), int, str, frac (numerator of p over n_assignments),
#       ci_lo / ci_hi (from a "[lo, hi]" string), bool
S = []


def add(name, file, where, col, kind="num", dp=None, note=""):
    S.append(dict(name=name, file=file, where=where, col=col, kind=kind, dp=dp, note=note))


PRIMARY = {"set": "y_assay/BRCA1_included", "calib": "isotonic"}
# ---- H1 (Results 3.1, Abstract, Table 1)
for c in ("fusion_rho", "single_rho", "delta"):
    add(f"h1.overall.{c}", "phase2_H1_stratified.csv", {"stratum": "overall"}, c, dp=3)
add("h1.overall.ci_lo", "phase2_H1_stratified.csv", {"stratum": "overall"}, "ci95", "ci_lo", 3)
add("h1.overall.ci_hi", "phase2_H1_stratified.csv", {"stratum": "overall"}, "ci95", "ci_hi", 3)
for s in ("core_like", "region"):
    add(f"h1.{s}.delta", "phase2_H1_stratified.csv", {"stratum": s}, "delta", dp=4)
# ---- H3 (Results 3.3)
for cond in ("y_assay/BRCA1_included", "y_assay/BRCA1_excluded", "y_clinvar/BRCA1_included", "y_clinvar/BRCA1_excluded"):
    add(f"h3.dBrier.{cond}", "phase3_H3_headline.csv", {"set": cond, "calib": "isotonic"}, "dBrier", dp=4)
add("h3.yield.n_significant_isotonic", "phase3_H3_headline.csv", {"calib": "isotonic"}, "dYield_ci", "count_ci_lo_gt0")
for obj, tag in (("fusion_M1", "fusion"), ("single:pangolin", "pangolin")):
    add(f"h3.yield.{tag}.frac_pct", "phase3_calibration_summary.csv", {**PRIMARY, "model": obj}, "actionable_frac", "pct", 1)
    add(f"h3.yield.{tag}.acc_pct", "phase3_calibration_summary.csv", {**PRIMARY, "model": obj}, "actionable_acc", "pct", 1)
add("h3.yield.direction", "phase3_calibration_summary.csv", PRIMARY, "actionable_frac", "yield_direction",
    note="sign of fusion minus Pangolin high-confidence fraction in the primary condition")
# ---- TP53 (Results 3.4, Table 4)
add("tp53.fusion.brier", "phase4_tp53_external.csv", {"model": "startswith:fusion"}, "Brier", dp=3)
add("tp53.best.brier", "phase4_tp53_external.csv", {"model": "startswith:best_single"}, "Brier", dp=3)
add("tp53.fusion.row", "phase4_tp53_external.csv", {"model": "startswith:fusion"}, "model", "str")
add("tp53.best.row", "phase4_tp53_external.csv", {"model": "startswith:best_single"}, "model", "str")
# ---- sampling frame (Methods 2.4)
for c in ("w_min", "w_max", "ESS_over_n"):
    add(f"ipw.decile.{c}", "ipw_weight_diagnostics.csv", {"scheme": "decile"}, c, dp=2)
add("ipw.whole.observed", "ipw_frame_balance_wholegene.csv", {"gene": "POOLED"}, "mean_abs_z_observed", dp=2)
add("ipw.whole.remainder", "ipw_frame_balance_wholegene.csv", {"gene": "POOLED"}, "mean_abs_z_remainder", dp=2)
add("ipw.window.observed", "ipw_frame_balance.csv", {"gene": "POOLED"}, "mean_abs_z_observed", dp=2)
add("ipw.window.remainder", "ipw_frame_balance.csv", {"gene": "POOLED"}, "mean_abs_z_remainder", dp=2)
# ---- H4 (Results 3.3)
H4 = {"condition": "y_assay/BRCA1_included", "calibration": "raw", "cal_method": "isotonic"}
add("h4.raw.fusion.lr_plus", "phase5_lr_headline.csv", H4, "lr_plus_fusion", dp=1)
add("h4.raw.best.lr_plus", "phase5_lr_headline.csv", H4, "lr_plus_best_single", dp=1)
add("h4.raw.fusion.lr_plus_interp", "phase5_lr_headline.csv", H4, "lr_plus_fusion_interp", dp=1)
add("h4.raw.best.lr_plus_interp", "phase5_lr_headline.csv", H4, "lr_plus_best_single_interp", dp=1)
for cond in ("y_clinvar/BRCA1_included", "y_clinvar/BRCA1_excluded"):
    add(f"h4.clinvar_strong_raw_observed.{cond}", "phase5_likelihood_ratios.csv",
        {"condition": cond, "calibration": "raw", "cal_method": "isotonic"}, "acmg_tier_observed", "count_eq:Strong")
    add(f"h4.clinvar_strong_raw_interp.{cond}", "phase5_likelihood_ratios.csv",
        {"condition": cond, "calibration": "raw", "cal_method": "isotonic"}, "acmg_tier", "count_eq:Strong")
add("h4.assay_strong_interp_any_stage", "phase5_likelihood_ratios.csv",
    {"condition": "startswith:y_assay", "cal_method": "isotonic"}, "acmg_tier", "count_eq:Strong")
add("h4.assay_strong_observed_any_stage", "phase5_likelihood_ratios.csv",
    {"condition": "startswith:y_assay", "cal_method": "isotonic"}, "acmg_tier_observed", "count_eq:Strong")
add("h4.interp_resolves", "phase5_lr_headline.csv", {"cal_method": "isotonic"}, "fusion_better_interp", "rows_true",
    note="condition|calibration rows whose interpolated paired-difference interval excludes zero")
add("h4.observed_resolves", "phase5_lr_headline.csv", {"cal_method": "isotonic"}, "fusion_better", "rows_true")
add("h4.interp_tier_crossings", "phase5_lr_headline.csv", {"cal_method": "isotonic"}, "crosses_a_tier_boundary", "count_true")
add("h4.observed_tier_crossings", "phase5_lr_headline.csv", {"cal_method": "isotonic"}, "crosses_a_tier_boundary_observed", "count_true")
add("h4.cal_effect.median_abs_change_interp", "phase5_calibration_effect.csv", {"cal_method": "isotonic"}, "median_abs_change_lr_plus_interp", dp=2)
add("h4.cal_effect.cells_changing_tier", "phase5_calibration_effect.csv", {"cal_method": "isotonic"}, "cells_changing_tier", "int")
add("h4.cal_effect.cells_with_a_tier", "phase5_calibration_effect.csv", {"cal_method": "isotonic"}, "cells_with_a_tier", "int")
add("h5b.fusion.yield_conf", "phase5b_evidence_yield.csv", {"condition": "y_assay/BRCA1_included", "object": "fusion_M1"}, "yield_confidence_0.90_0.10", dp=3)
add("h5b.fusion.yield_moderate", "phase5b_evidence_yield.csv", {"condition": "y_assay/BRCA1_included", "object": "fusion_M1"}, "yield_acmg_moderate_or_above", dp=3)
# ---- weights / thinning (Discussion)
add("thin.mean6", "phase6_training_gene_summary.csv", {"n_training_genes": 6}, "mean", dp=4)
add("thin.mean3", "phase6_training_gene_summary.csv", {"n_training_genes": 3}, "mean", dp=4)
add("thin.n_genes_worse", "phase6_training_gene_summary.csv", {"n_training_genes": 3}, "n_genes_worse_than_6", "int")
for f in ("alphagenome", "pangolin", "spliceai"):
    add(f"weights.{f}.mean", "phase6_enet_weight_stability.csv", {"_index": f}, "mean", dp=3)
    add(f"weights.{f}.sd", "phase6_enet_weight_stability.csv", {"_index": f}, "sd", dp=3)
# ---- label contrast (3.3)
for c in ("assay_labelled", "clinvar_labelled", "both", "assay_only"):
    add(f"labels.{c}", "phase7_base.csv", {}, c, "int")
add("labels.agreement", "phase7_label_agreement.csv", {}, "agreement", dp=3)
add("labels.kappa", "phase7_label_agreement.csv", {}, "kappa", dp=2)
add("labels.pos_assay", "phase7_label_agreement.csv", {}, "positive_rate_assay", dp=2)
add("labels.pos_clinvar", "phase7_label_agreement.csv", {}, "positive_rate_clinvar", dp=2)
add("labels.strong_fun", "phase7_label_contrast.csv", {}, "tier_fun", "count_eq:Strong")
add("labels.strong_cv", "phase7_label_contrast.csv", {}, "tier_cv", "count_eq:Strong")
add("sel.fusion.recorded.sens", "phase7_selection_test.csv", {"object": "fusion_M1", "subset": "ClinVar-recorded"}, "sensitivity", dp=3)
add("sel.fusion.assay_only.sens", "phase7_selection_test.csv", {"object": "fusion_M1", "subset": "assay-only"}, "sensitivity", dp=3)
# ---- phase 8
SF = "phase8_sign_flip_exact.csv"
add("sf.rho.incl.p", SF, {"analysis": "delta_rho", "condition": "BRCA1_included"}, "p_exact", dp=4)
add("sf.rho.incl.num", SF, {"analysis": "delta_rho", "condition": "BRCA1_included"}, "p_exact", "frac")
add("sf.rho.excl.p", SF, {"analysis": "delta_rho", "condition": "BRCA1_excluded"}, "p_exact", dp=4)
add("sf.brier.primary.num", SF, {"analysis": "delta_brier", "condition": "isotonic/y_assay/BRCA1_included"}, "p_exact", "frac")
add("sf.brier.primary.p", SF, {"analysis": "delta_brier", "condition": "isotonic/y_assay/BRCA1_included"}, "p_exact", dp=4)
add("sf.brier.primary.T_obs", SF, {"analysis": "delta_brier", "condition": "isotonic/y_assay/BRCA1_included"}, "T_obs", dp=4)
for cond in ("isotonic/y_assay/BRCA1_included", "isotonic/y_assay/BRCA1_excluded", "isotonic/y_clinvar/BRCA1_included",
             "isotonic/y_clinvar/BRCA1_excluded", "platt/y_assay/BRCA1_included", "platt/y_assay/BRCA1_excluded",
             "platt/y_clinvar/BRCA1_included", "platt/y_clinvar/BRCA1_excluded"):
    add(f"sf.brier.{cond}.p", SF, {"analysis": "delta_brier", "condition": cond}, "p_exact", dp=4)
    add(f"sf.brier.{cond}.obs", SF, {"analysis": "delta_brier", "condition": cond}, "delta_brier_obs", dp=4)
add("sf.brier.primary.wild", SF, {"analysis": "delta_brier", "condition": "isotonic/y_assay/BRCA1_included"}, "p_wild_rademacher", dp=4)
add("sf.rho.incl.wild", SF, {"analysis": "delta_rho", "condition": "BRCA1_included"}, "p_wild_rademacher", dp=4)
L2 = "phase8_leave_two_genes_out.csv"
add("l2go.min.dropped", L2, {"_argmin": "delta_rho"}, "dropped", "str")
add("l2go.max.dropped", L2, {"_argmax": "delta_rho"}, "dropped", "str")
add("l2go.min", L2, {}, "delta_rho", "agg:min", 4)
add("l2go.median", L2, {}, "delta_rho", "agg:median", 4)
add("l2go.max", L2, {}, "delta_rho", "agg:max", 4)
HK = "phase8_leaderboard_hk.csv"
for c in ("pooled_rho", "dl_lo", "dl_hi", "hk_lo", "hk_hi"):
    add(f"hk.M1_enet.{c}", HK, {"model": "M1_enet"}, c, dp=3)
for c in ("dl_lo", "dl_hi", "hk_lo", "hk_hi"):
    add(f"hk.nt.{c}", HK, {"model": "single:nt"}, c, dp=3)
for c in ("hk_lo", "hk_hi"):
    add(f"hk.gpn_msa.{c}", HK, {"model": "single:gpn_msa"}, c, dp=3)
MU = {"condition": "isotonic/y_assay/BRCA1_included", "binning": "width10"}
for c in ("brier_fusion", "brier_single", "rel_fusion", "rel_single", "res_fusion", "res_single", "unc",
          "d_brier", "d_brier_lo", "d_brier_hi", "d_rel", "d_rel_lo", "d_rel_hi", "d_res", "d_res_lo", "d_res_hi"):
    add(f"murphy.{c}", "phase8_murphy_decomposition.csv", MU, c, dp=4)
SEL = "phase8_selection_test_stratified.csv"
for stratum, obj, subset in (("overall", "fusion_M1", "ClinVar-recorded"), ("overall", "fusion_M1", "assay-only"),
                             ("overall", "cadd", "ClinVar-recorded"), ("overall", "cadd", "assay-only"),
                             ("minus_BRCA1", "fusion_M1", "ClinVar-recorded"), ("minus_BRCA1", "fusion_M1", "assay-only"),
                             ("minus_BRCA1", "cadd", "ClinVar-recorded"), ("minus_BRCA1", "cadd", "assay-only"),
                             ("offset_core_like", "fusion_M1", "ClinVar-recorded"), ("offset_core_like", "fusion_M1", "assay-only"),
                             ("offset_core_like", "cadd", "ClinVar-recorded"), ("offset_core_like", "cadd", "assay-only"),
                             ("offset_region", "fusion_M1", "ClinVar-recorded"), ("offset_region", "fusion_M1", "assay-only"),
                             ("offset_region", "cadd", "ClinVar-recorded"), ("offset_region", "cadd", "assay-only"),
                             ("assay_only_VUS", "fusion_M1", "assay-only"), ("assay_only_VUS", "cadd", "assay-only"),
                             ("assay_only_Conflicting/Other", "fusion_M1", "assay-only"), ("assay_only_Conflicting/Other", "cadd", "assay-only")):
    add(f"sel8.{stratum}.{obj}.{subset}", SEL, {"stratum": stratum, "object": obj, "subset": subset}, "lr_plus", dp=1)
PG = "phase8_per_gene_deltas.csv"
add("pergene.rho_favours_n", PG, {}, "rho_favours_fusion", "count_true")
add("pergene.brier_favours_n", PG, {}, "brier_favours_fusion", "count_true")
add("pergene.rho_dissent", PG, {"rho_favours_fusion": False}, "gene", "str")
add("pergene.brier_dissent", PG, {"brier_favours_fusion": False}, "gene", "str")
add("pergene.best_single", PG, {"gene": "BAP1"}, "best_single", "str")
LD = "phase4_tp53_label_definitions.csv"
add("tp53.rules.fusion_better_n", LD, {}, "fusion_better", "count_true")
for rule, tag in (("control-anchored RFS>0", "control"), ("median split (top 50%)", "median"), ("NA mid-band |RFS|>=0.5", "midband")):
    add(f"tp53.rules.{tag}.dBrier", LD, {"label_rule": rule}, "dBrier", dp=3)
    add(f"tp53.rules.{tag}.ci_lo", LD, {"label_rule": rule}, "ci_lo", dp=3)
    add(f"tp53.rules.{tag}.ci_hi", LD, {"label_rule": rule}, "ci_hi", dp=3)
OP = "phase8_lr_plus_operating_points.csv"
for brca in ("included", "excluded"):
    for spec in (0.9, 0.95, 0.975, 0.99):
        add(f"op.fusion.{brca}.{spec}", OP, {"condition": f"y_assay/BRCA1_{brca}", "object": "fusion_M1", "target_spec": spec}, "lr_plus", dp=2)
        add(f"op.fusion.{brca}.{spec}.interp", OP, {"condition": f"y_assay/BRCA1_{brca}", "object": "fusion_M1", "target_spec": spec}, "lr_plus_interp", dp=2)
for obj, spec in (("pangolin", 0.95), ("spliceai", 0.95), ("cadd", 0.99), ("nt", 0.95), ("mean_M0b", 0.95)):
    add(f"op.{obj}.included.{spec}", OP, {"condition": "y_assay/BRCA1_included", "object": obj, "target_spec": spec}, "lr_plus", dp=2)
    add(f"op.{obj}.included.{spec}.interp", OP, {"condition": "y_assay/BRCA1_included", "object": obj, "target_spec": spec}, "lr_plus_interp", dp=2)


def _select(df: pd.DataFrame, where: dict) -> pd.DataFrame:
    if "_index" in where:
        return df[df.iloc[:, 0] == where["_index"]]
    if "_argmin" in where:
        return df.loc[[df[where["_argmin"]].idxmin()]]
    if "_argmax" in where:
        return df.loc[[df[where["_argmax"]].idxmax()]]
    m = np.ones(len(df), bool)
    for k, v in where.items():
        col = df[k]
        if isinstance(v, str) and v.startswith("startswith:"):
            m &= col.astype(str).str.startswith(v.split(":", 1)[1]).to_numpy()
        elif isinstance(v, float):
            m &= np.isclose(col.astype(float), v)
        else:
            m &= (col == v).to_numpy()
    return df[m]


def _fmt(x: float, dp: int) -> str:
    return f"{x:.{dp}f}"


def compute(entry: dict) -> dict:
    df = pd.read_csv(REPORTS / entry["file"])
    sub = _select(df, entry["where"])
    kind, col, dp = entry["kind"], entry["col"], entry["dp"]
    if kind.startswith("count_eq:"):
        val = int((sub[col] == kind.split(":", 1)[1]).sum())
    elif kind == "count_true":
        val = int(sub[col].astype(str).str.lower().eq("true").sum())
    elif kind == "rows_true":
        keep = sub[col].astype(str).str.lower().eq("true")
        val = sorted(f"{r.condition}|{r.calibration}" for r in sub[keep].itertuples())
    elif kind == "count_ci_lo_gt0":
        val = int(sum(float(re.findall(r"-?\d+\.?\d*", str(s))[0]) > 0 for s in sub[col]))
    elif kind == "yield_direction":
        f = float(sub[sub["model"] == "fusion_M1"][col].iloc[0]); p = float(sub[sub["model"] == "single:pangolin"][col].iloc[0])
        val = "fusion_higher" if f > p else ("pangolin_higher" if p > f else "equal")
    elif kind.startswith("agg:"):
        val = _fmt(float(getattr(sub[col], kind.split(":", 1)[1])()), dp)
    else:
        assert len(sub) == 1, f"{entry['name']}: selector matched {len(sub)} rows"
        raw = sub[col].iloc[0]
        if kind == "num":
            val = _fmt(float(raw), dp)
        elif kind == "pct":
            val = _fmt(float(raw) * 100, dp)
        elif kind == "int":
            val = int(raw)
        elif kind == "str":
            val = str(raw)
        elif kind == "frac":
            n = int(sub["n_assignments"].iloc[0])
            val = {"numerator": int(round(float(raw) * n)), "denominator": n}
        elif kind == "ci_lo":
            val = _fmt(float(re.findall(r"-?\d+\.?\d*", str(raw))[0]), dp)
        elif kind == "ci_hi":
            val = _fmt(float(re.findall(r"-?\d+\.?\d*", str(raw))[1]), dp)
        else:
            raise ValueError(kind)
    return {"value": val, "file": entry["file"], "where": entry["where"], "col": col, "kind": kind, "dp": dp,
            **({"note": entry["note"]} if entry["note"] else {})}


def build() -> dict:
    values = {e["name"]: compute(e) for e in S}
    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip()
    except Exception:
        commit = "unknown"
    return {"_generated_by": "scripts/build_expected_values.py", "_generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "_git_commit": commit, "_reports": str(REPORTS.relative_to(REPO)), "values": values}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    new = build()
    if a.check:
        old = json.loads(OUT.read_text())["values"] if OUT.exists() else {}
        stale = [k for k, v in new["values"].items() if old.get(k, {}).get("value") != v["value"]]
        missing = [k for k in old if k not in new["values"]]
        if stale or missing:
            print("STALE:", stale, "| MISSING:", missing)
            return 1
        print(f"fixture up to date ({len(new['values'])} values)")
        return 0
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(new, indent=1) + "\n")
    print(f"wrote {OUT.relative_to(REPO)} with {len(new['values'])} values from {REPORTS.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
