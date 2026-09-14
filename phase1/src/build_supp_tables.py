"""Regenerate the pipeline-derived supplementary tables from the report CSVs.

The supplement states that every table is compiled without recomputation from
the frozen analysis outputs. This script is that compilation: each supported
table is rebuilt cell by cell from the CSV named in its legend, using one
formatting rule per table, so the document cannot drift from the outputs.

    python -m src.build_supp_tables --in S.docx --out S.docx --reports reports/phase1
    python -m src.build_supp_tables --in S.docx --reports reports/phase1_v1 --verify

`--verify` regenerates every table in memory and compares it with what the
document holds, cell by cell, without writing anything; run against the v1
outputs it must reproduce the published (v1) supplement, which is how the
formatting rules were checked before the v2 tables were written.

Tables and sources
------------------
  S2   coverage_by_region.csv               S3a  phase3_H3_headline.csv
  S3b  phase3_calibration_summary.csv       S4   phase3_pertool_brier_ci.csv
  S5a  phase2_H2_ablation.csv               S5b  phase2_leaderboard.csv
  S6   directionality_check.csv             S8a  phase4_tp53_external*.csv
  S8b  phase4_tp53_label_definitions.csv    S10a-d  ipw_*.csv (src.build_supp_ipw)
  S11  phase5_likelihood_ratios.csv (+ interpolated columns)
  S12  phase6_enet_weight_stability.csv     S13  phase6_training_gene_summary.csv
  S14  phase7_label_contrast.csv + phase7_selection_test.csv
  S15  phase8_murphy_decomposition.csv      S16  phase8_sign_flip_exact / leave_two_genes_out / leaderboard_hk
  S17  phase8_selection_test_stratified.csv S19  phase8_lr_plus_operating_points.csv (+ interpolated)
  S20  phase8_lr_placement_decomposition.csv (new with frozen-matrix-v2)
S1, S7 and S18 are provenance tables with no pipeline output behind them and are
left as they are.
"""
from __future__ import annotations

import argparse
import copy
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import docx
from docx.oxml.ns import qn
from docx.table import Table
from docx.text.paragraph import Paragraph

from .build_supp_ipw import panels as ipw_panels, PRETTY, COND

MINUS = "−"


def f(x, nd):
    return "" if pd.isna(x) else f"{float(x):.{nd}f}"


def s(x, nd):
    """Signed, with a true minus sign."""
    return "" if pd.isna(x) else f"{float(x):+.{nd}f}".replace("-", MINUS)


def ne(x, nd):
    return "n/e" if pd.isna(x) else f"{float(x):.{nd}f}"


def na(x, nd):
    return "n/a" if pd.isna(x) else f"{float(x):.{nd}f}"


def _half_up(x, nd):
    from decimal import Decimal, ROUND_HALF_UP
    return str(Decimal(str(x)).quantize(Decimal(1).scaleb(-nd), rounding=ROUND_HALF_UP))


def col_digits(values, main_nd, at_least=None, main_scale=1.0):
    """Decimals to print a column at so that a reader who rounds the printed value
    (half-up) to any precision the main text uses gets the main text's cell.

    main_nd is one precision or several (Brier is quoted to three decimals in Table 3
    and four in the Figure 2 caption); the main text prints value * main_scale (a
    fraction as a percentage has main_scale 100). Starts one decimal above the finest
    main precision and adds digits until no printed value is a tie the unrounded value
    is not: 0.032495 printed at five decimals is 0.03250, which rounds to 0.033, but the
    value itself rounds to 0.032, so that column needs six. The arithmetic is decimal,
    as a reader's is: in binary floats 0.7865 * 100 is 78.649..., which hides exactly
    the tie this exists to catch.
    """
    from decimal import Decimal
    mains = [main_nd] if isinstance(main_nd, int) else list(main_nd)
    vals = [Decimal(repr(float(v))) for v in values if pd.notna(v)]
    scale = Decimal(repr(float(main_scale)))
    shift = len(str(int(main_scale))) - 1          # 100 -> two decimals of the fraction
    nd = max(max(mains) + shift + 1, at_least or 0)
    while True:
        ok = all(_half_up(Decimal(_half_up(v, nd)) * scale, m) == _half_up(v * scale, m)
                 for v in vals for m in mains)
        if ok:
            return nd
        nd += 1
        assert nd <= 10, "no printable precision reproduces the main-text rounding"


def fmt_col(values, main_nd, at_least=None, main_scale=1.0, blank=""):
    """Format a column at col_digits' precision; NaN prints as `blank`."""
    nd = col_digits(values, main_nd, at_least, main_scale)
    return [blank if pd.isna(v) else _half_up(float(v), nd) for v in values]


def is_number(s):
    try:
        float(str(s).replace("−", "-"))
        return True
    except ValueError:
        return False


def raw(df):
    """The CSV's own text for every cell (what the earlier tables were pasted from)."""
    return df.fillna("").astype(str)


# ---------------------------------------------------------------------------
# table builders: each returns (header, rows) as lists of strings
# ---------------------------------------------------------------------------
def t_s2(rep):
    d = pd.read_csv(rep / "coverage_by_region.csv")
    cols = [c for c in d.columns if c != "predictor"]
    return (["predictor", *cols],
            [[r["predictor"], *[f"{_half_up(100 * float(r[c]), 1)}%" for c in cols]] for _, r in d.iterrows()])


def t_s3a(rep):
    d = pd.read_csv(rep / "phase3_H3_headline.csv")
    # the main text prints these differences to four decimals; the interval strings are
    # already formatted once by phase 3 and are copied
    dece, dbri, dyld = (fmt_col(d[c], 4, at_least=5) for c in ("dECE", "dBrier", "dYield"))
    return (["Condition", "Calib", "ΔECE", "ΔECE 95% CI", "ΔBrier", "ΔBrier CI (gene-clust)",
             "ΔBrier CI (BCa)", "ΔBrier CI (variant)", "ΔYield", "ΔYield 95% CI", "Fusion better on"],
            [[COND[r["set"]], r["calib"], a, r["dECE_ci"], b, r["dBrier_ci_geneclust"],
              r["dBrier_ci_BCa"], r["dBrier_ci_variant(sens)"], c, r["dYield_ci"], r["fusion_better_on"]]
             for (_, r), a, b, c in zip(d.iterrows(), dece, dbri, dyld)])


def _summary_rows(d):
    """Table 3 prints ECE and Brier to three decimals and the two yields as percentages
    to one decimal; the fractions stay fractions here, printed at the precision
    col_digits finds for a reader who multiplies by 100 and rounds."""
    ece, bri = fmt_col(d["ECE"], 3), fmt_col(d["Brier"], (3, 4))   # Brier: Table 3 and the Figure 2 caption
    frac = fmt_col(d["actionable_frac"], 1, main_scale=100.0)
    acc = fmt_col(d["actionable_acc"], 1, main_scale=100.0)
    n = ["" if pd.isna(v) else str(int(v)) for v in d["n_actionable"]]
    return list(zip(ece, bri, frac, acc, n))


def t_s3b(rep):
    d = pd.read_csv(rep / "phase3_calibration_summary.csv")
    return (["Condition", "Calib", "Object", "ECE", "Brier", "High-conf. frac", "High-conf. acc", "n high-conf."],
            [[COND[r["set"]], r["calib"], r["model"], *vals]
             for (_, r), vals in zip(d.iterrows(), _summary_rows(d))])


def t_s4(rep):
    d = pd.read_csv(rep / "phase3_pertool_brier_ci.csv")
    dbri = fmt_col(d["dBrier_vs_fusion"], 4, at_least=5)
    return (["Condition", "Calib", "Single tool", "ΔBrier vs fusion", "gene-clust 95% CI", "Fusion lower"],
            [[COND[r["set"]], r["calib"], r["tool"], v, "" if pd.isna(r["ci95"]) else r["ci95"],
              "" if pd.isna(r["fusion_lower"]) else str(r["fusion_lower"])]
             for (_, r), v in zip(d.iterrows(), dbri)])


def t_s5a(rep):
    d = pd.read_csv(rep / "phase2_H2_ablation.csv")
    full, abl = fmt_col(d["full_rho"], 3, at_least=4), fmt_col(d["ablated_rho"], 3, at_least=4)
    delta = fmt_col(d["delta_full_minus_ablated"], (3, 4), at_least=4)
    return (["Model", "Ablation", "Full ρ", "Ablated ρ", "Δρ (full−ablated)", "95% CI", "Evo contributes"],
            [[r["model"], r["ablation"], a, b, c, r["ci95"], r["evo_contributes"]]
             for (_, r), a, b, c in zip(d.iterrows(), full, abl, delta)])


def t_s5b(rep):
    d = pd.read_csv(rep / "phase2_leaderboard.csv")
    # Table 1 prints rho and the interval to three decimals and I^2 as an integer
    # Table 1 prints rho and the interval to three decimals, the abstract rho to two
    nd_rho = col_digits(pd.concat([d["pooled_rho"], d.get("lo", d["pooled_rho"]), d.get("hi", d["pooled_rho"])]), (2, 3), at_least=4)
    rho = ["" if pd.isna(v) else _half_up(v, nd_rho) for v in d["pooled_rho"]]
    if "lo" in d.columns:
        ci = ["" if pd.isna(lo) else f"[{_half_up(lo, nd_rho)}, {_half_up(hi, nd_rho)}]" for lo, hi in zip(d["lo"], d["hi"])]
    else:
        ci = ["" if (pd.isna(v) or v == "n/a") else v for v in d["ci95"]]
    i2 = fmt_col(d["I2"], 0, at_least=1)
    return (["Model / tool", "Pooled ρ", "95% CI", "I² (%)", "k"],
            [[r["model"], a, b, c, str(int(r["k"]))] for (_, r), a, b, c in zip(d.iterrows(), rho, ci, i2)])


def t_s6(rep):
    d = raw(pd.read_csv(rep / "directionality_check.csv", dtype=str))
    return (["Gene", "n LoF", "n syn", "median path. LoF", "median path. syn", "Δ (LoF−syn)",
             "control AUROC", "Verdict"],
            [[r["gene"], r["n_lof"], r["n_syn"], r["median_path_lof"], r["median_path_syn"],
              r["delta_lof_minus_syn"], r["control_auroc"], r["verdict"]] for _, r in d.iterrows()])


def _tp53_n(rep):
    ld = pd.read_csv(rep / "phase4_tp53_label_definitions.csv")
    n_main = int(ld[ld.label_rule.str.startswith("control")].n.iloc[0])
    n_band = int(ld[ld.label_rule.str.startswith("NA")].n.iloc[0])
    return n_main, n_band


def t_s8a(rep):
    n_main, n_band = _tp53_n(rep)
    rows = []
    for fn, lab in (("phase4_tp53_external.csv", f"Control-anchored RFS>0 (n={n_main})"),
                    ("phase4_tp53_external_naband.csv", f"NA mid-band |RFS|≥0.5 (n={n_band})")):
        d = pd.read_csv(rep / fn)
        rows += [[lab, r["model"], *vals] for (_, r), vals in zip(d.iterrows(), _summary_rows(d))]
    return (["Label", "Object", "ECE", "Brier", "High-conf. frac", "High-conf. acc", "n high-conf."], rows)


LABEL_RULE = {"control-anchored RFS>0": "Control-anchored (RFS>0)",
              "median split (top 50%)": "Median split (top 50%)",
              "NA mid-band |RFS|>=0.5": "NA mid-band (|RFS|≥0.5)"}


def t_s8b(rep):
    d = pd.read_csv(rep / "phase4_tp53_label_definitions.csv")
    return (["Label definition", "n", "ΔBrier (best−fusion)", "variant-level 95% CI", "Favours fusion"],
            [[LABEL_RULE.get(r.label_rule, r.label_rule), str(int(r.n)), f(r.dBrier, 3),
              f"{f(r.ci_lo, 3)} – {f(r.ci_hi, 3)}", "yes" if bool(r.fusion_better) else "no"]
             for r in d.itertuples()])


def _obj(o):
    return o[len("single:"):] if o.startswith("single:") else o


def t_s11(rep):
    d = pd.read_csv(rep / "phase5_likelihood_ratios.csv")
    d = d[d.cal_method == "isotonic"]
    stage = {"raw": "raw", "calibrated_isotonic": "calibrated"}
    rows = []
    dual = "acmg_tier_observed" in d.columns   # the v1 tables predate the interpolated basis
    for _, r in d.iterrows():
        rows.append([r.condition, stage[r.calibration], _obj(r.object), str(int(r.n)), str(int(r.n_pos)),
                     str(int(r.n_neg)), na(r.lr_plus, 3), na(r.lr_plus_lo, 3), na(r.lr_plus_hi, 3),
                     na(r.sensitivity_at_that_point, 3), na(r.specificity_achieved, 3),
                     r.acmg_tier_observed if dual else r.acmg_tier,
                     na(r.lr_minus_at_spec95, 3), na(r.lr_minus_at_sens95, 3), r.benign_tier_at_sens95,
                     # frozen-matrix-v2: the interpolated basis
                     *([na(r.lr_plus_interp, 3), na(r.lr_plus_interp_lo, 3), na(r.lr_plus_interp_hi, 3),
                        na(r.sensitivity_interp, 3), r.acmg_tier] if dual else [""] * 5)])
    return (["Condition", "Calibration", "Object", "n", "n pos", "n neg", "LR+ (scanned)", "LR+ lo", "LR+ hi",
             "Sensitivity", "Specificity achieved", "Band (scanned)", "LR− at spec 95", "LR− at sens 95",
             "Benign band", "LR+ (interpolated)", "interp lo", "interp hi", "Sensitivity (interp.)",
             "Band (interpolated)"], rows)


def t_s12(rep):
    d = pd.read_csv(rep / "phase6_enet_weight_stability.csv")
    d = d[(d["mean"].abs() > 0) | (d["mean_abs"] > 0)]
    return (["Feature", "Mean", "SD", "Min", "Max", "Range", "Mean |coef|", "Folds non-zero",
             "Folds with sign flip", "Range / mean |coef|"],
            [[r[1], f(r.mean, 3), f(r.sd, 3), f(r.min, 3), f(r.max, 3), f(r.range, 3), f(r.mean_abs, 3),
              str(int(r.n_folds_nonzero)), str(int(r.n_folds_sign_flip)), f(r.range_over_mean_abs, 3)]
             for r in d.itertuples()])


def t_s13(rep):
    d = pd.read_csv(rep / "phase6_training_gene_summary.csv")
    return (["Training genes", "Subsets evaluated", "Mean Brier", "SD", "Min", "Max",
             "Held-out genes worse than at 6"],
            [[str(int(r["n_training_genes"])), str(int(r["count"])), f(r["mean"], 4), f(r["std"], 4),
              f(r["min"], 4), f(r["max"], 4), str(int(r.n_genes_worse_than_6))] for _, r in d.iterrows()])


def t_s14(rep):
    d = pd.read_csv(rep / "phase7_label_contrast.csv")
    st = pd.read_csv(rep / "phase7_selection_test.csv")
    rec = st[st.subset == "ClinVar-recorded"].set_index("object").lr_plus
    ao = st[st.subset == "assay-only"].set_index("object").lr_plus
    d = d.sort_values("delta_lr", ascending=False, na_position="last")
    rows = []
    for r in d.itertuples():
        rows.append([r.object, r.training_group, str(int(r.n_fun)), na(r.lr_fun, 3), na(r.sens_fun, 3),
                     r.tier_fun, na(r.lr_cv, 3), na(r.sens_cv, 3), r.tier_cv, na(r.delta_lr, 3),
                     na(r.delta_sens, 3), na(rec.get(r.object, np.nan), 3), na(ao.get(r.object, np.nan), 3)])
    return (["Object", "Training signal", "n", "LR+ (assay)", "Sens (assay)", "Band (assay)", "LR+ (ClinVar)",
             "Sens (ClinVar)", "Band (ClinVar)", "ΔLR+", "ΔSens", "LR+ recorded (assay labels)",
             "LR+ assay-only (assay labels)"], rows)


MURPHY_COND = {"isotonic/y_assay/BRCA1_included": "Assay, +BRCA1", "isotonic/y_assay/BRCA1_excluded": "Assay, −BRCA1",
               "isotonic/y_clinvar/BRCA1_included": "ClinVar, +BRCA1", "isotonic/y_clinvar/BRCA1_excluded": "ClinVar, −BRCA1"}
BINNING = {"width10": "width-10", "freq10": "freq-10", "width15": "width-15"}


def t_s15a(rep):
    d = pd.read_csv(rep / "phase8_murphy_decomposition.csv")
    r = d[(d.condition == "isotonic/y_assay/BRCA1_included") & (d.binning == "width10")].iloc[0]
    return (["Object", "Brier", "REL", "RES", "UNC"],
            [["Fusion (M1)", f(r.brier_fusion, 4), f(r.rel_fusion, 4), f(r.res_fusion, 4), f(r.unc, 4)],
             ["Pangolin", f(r.brier_single, 4), f(r.rel_single, 4), f(r.res_single, 4), f(r.unc, 4)]])


def t_s15b(rep):
    d = pd.read_csv(rep / "phase8_murphy_decomposition.csv")
    rows = []
    for cond in MURPHY_COND:
        for b in BINNING:
            r = d[(d.condition == cond) & (d.binning == b)].iloc[0]
            rows.append([MURPHY_COND[cond], BINNING[b], s(r.d_brier, 4),
                         f"{s(r.d_rel, 4)} [{s(r.d_rel_lo, 4)}, {s(r.d_rel_hi, 4)}]",
                         f"{s(r.d_res, 4)} [{s(r.d_res_lo, 4)}, {s(r.d_res_hi, 4)}]"])
    return (["Condition", "Binning", "ΔBrier", "ΔREL [95% CI]", "ΔRES [95% CI]"], rows)


def _sf_label(analysis, condition):
    if analysis == "delta_rho":
        return "Δρ", "functional, " + ("+BRCA1" if condition.endswith("included") else "−BRCA1")
    method, label, brca = condition.split("/")
    m = {"isotonic": "isotonic", "platt": "Platt"}[method]
    c = {"y_assay": "Assay", "y_clinvar": "ClinVar"}[label]
    return "ΔBrier", f"{m} × {c}, " + ("+BRCA1" if brca.endswith("included") else "−BRCA1")


def t_s16a(rep):
    d = pd.read_csv(rep / "phase8_sign_flip_exact.csv")
    rows = []
    for r in d.itertuples():
        q, c = _sf_label(r.analysis, r.condition)
        rows.append([q, c, f(r.p_exact, 4), f(r.p_wild_rademacher, 4)])
    return (["Quantity", "Condition", "Exact sign-flip p", "Wild bootstrap p"], rows)


def t_s16b(rep):
    d = pd.read_csv(rep / "phase8_leave_two_genes_out.csv").sort_values("dropped")
    return (["Dropped pair", "k remaining", "Fusion ρ", "Best-single ρ", "Δρ"],
            [[r.dropped, str(int(r.k_remaining)), f(r.fusion_rho, 3), f(r.single_rho, 3),
              f"{float(r.delta_rho):+.4f}".replace("-", MINUS)] for r in d.itertuples()])


HK_NAMES = {"M1_enet": PRETTY["fusion_M1"], "M2_gbt": "Fusion — gradient-boosted trees (M2)",
            "M0b_mean": PRETTY["mean_M0b"]}


def t_s16c(rep):
    d = pd.read_csv(rep / "phase8_leaderboard_hk.csv")
    rows = []
    for r in d.itertuples():
        name = HK_NAMES.get(r.model, PRETTY.get(r.model, r.model))
        if pd.isna(r.pooled_rho):
            rows.append([name, "—", "—", "—"])
        else:
            rows.append([name, f(r.pooled_rho, 3),
                         f"[{f(r.dl_lo, 3)}, {f(r.dl_hi, 3)}]", f"[{f(r.hk_lo, 3)}, {f(r.hk_hi, 3)}]"])
    return (["Model / tool", "Pooled ρ", "DerSimonian–Laird 95% CI", "Hartung–Knapp 95% CI"], rows)


OBJ_ORDER = ["fusion_M1", "mean_M0b", "pangolin", "spliceai", "alphagenome", "gpn_msa", "cadd",
             "phylop", "phastcons", "nt", "gnomad_af", "alphamissense"]
S17_COLS = [("overall", "ClinVar-recorded"), ("overall", "assay-only"),
            ("minus_BRCA1", "ClinVar-recorded"), ("minus_BRCA1", "assay-only"),
            ("offset_core_like", "ClinVar-recorded"), ("offset_core_like", "assay-only"),
            ("offset_region", "ClinVar-recorded"), ("offset_region", "assay-only"),
            ("assay_only_VUS", "assay-only"), ("assay_only_Conflicting/Other", "assay-only")]


def _pretty(o):
    return PRETTY.get(o, PRETTY.get("single:" + o, o))


def t_s17(rep):
    d = pd.read_csv(rep / "phase8_selection_test_stratified.csv").set_index(["stratum", "subset", "object"]).lr_plus
    rows = [[_pretty(o), *[ne(d.get((st, sub, o), np.nan), 2) for st, sub in S17_COLS]] for o in OBJ_ORDER]
    return (["Object", "Overall recorded", "Overall assay-only", "−BRCA1 recorded", "−BRCA1 assay-only",
             "Core recorded", "Core assay-only", "Region recorded", "Region assay-only", "Assay-only: VUS",
             "Assay-only: Conflicting/Other"], rows)


SPECS = [0.9, 0.95, 0.975, 0.99]


def t_s19(rep):
    d = pd.read_csv(rep / "phase8_lr_plus_operating_points.csv")
    rows = []
    for cond in ("y_assay/BRCA1_included", "y_assay/BRCA1_excluded"):
        for o in OBJ_ORDER:
            sub = d[(d.condition == cond) & (d.object == o)].set_index("target_spec")
            obs = [ne(sub.lr_plus.get(t, np.nan), 2) for t in SPECS]
            itp = ([ne(sub.lr_plus_interp.get(t, np.nan), 2) for t in SPECS]
                   if "lr_plus_interp" in sub.columns else [""] * len(SPECS))   # v1 predates the interpolated basis
            rows.append([_pretty(o), COND[cond], *obs, *itp])
    return (["Object", "Condition", "LR+ @90% spec", "LR+ @95% spec", "LR+ @97.5% spec", "LR+ @99% spec",
             "interp. @90%", "interp. @95%", "interp. @97.5%", "interp. @99%"], rows)


def t_s20(rep):
    d = pd.read_csv(rep / "phase8_lr_placement_decomposition.csv")
    stage = {"raw": "raw", "calibrated_isotonic": "calibrated"}
    rows = []
    for r in d.itertuples():
        rows.append([r.condition, stage[r.calibration], _obj(r.object),
                     ne(r.lr_plus_obs_v1, 2), ne(r.lr_plus_obs_v2, 2), ne(r.spec_achieved_v2, 3),
                     ne(r.lr_plus_interp_v1, 2), ne(r.lr_plus_interp_v2, 2),
                     s(r.delta_obs, 2) or "n/e", s(r.delta_interp_score_change, 2) or "n/e",
                     s(r.delta_placement, 2) or "n/e", s(r.placement_within_v2, 2) or "n/e",
                     f"{r.tier_obs_v1} → {r.tier_obs_v2}", f"{r.tier_interp_v1} → {r.tier_interp_v2}"])
    return (["Condition", "Calibration", "Object", "LR+ scanned, v1", "LR+ scanned, v2", "Specificity achieved, v2",
             "LR+ interp., v1", "LR+ interp., v2", "Δ scanned (v2−v1)", "of which score change (Δ interp.)",
             "of which threshold placement", "Placement within v2 (scanned − interp.)", "Band (scanned) v1 → v2",
             "Band (interp.) v1 → v2"], rows)


def t_ipw(k):
    def build(rep):
        title, caption, header, rows = ipw_panels(rep)[k]
        return header, rows
    return build


# caption prefix -> builder; each caption may be followed by several tables (S15, S16)
TABLES = [
    ("Supplementary Table S2.", [t_s2]),
    ("Supplementary Table S3a.", [t_s3a]),
    ("Supplementary Table S3b.", [t_s3b]),
    ("Supplementary Table S4.", [t_s4]),
    ("Supplementary Table S5a.", [t_s5a]),
    ("Supplementary Table S5b.", [t_s5b]),
    ("Supplementary Table S6.", [t_s6]),
    ("Supplementary Table S8a.", [t_s8a]),
    ("Supplementary Table S8b.", [t_s8b]),
    ("Supplementary Table S10a.", [t_ipw(0)]),
    ("Supplementary Table S10b.", [t_ipw(1)]),
    ("Supplementary Table S10c.", [t_ipw(2)]),
    ("Supplementary Table S10d.", [t_ipw(3)]),
    ("Supplementary Table S11.", [t_s11]),
    ("Supplementary Table S12.", [t_s12]),
    ("Supplementary Table S13.", [t_s13]),
    ("Supplementary Table S14.", [t_s14]),
    ("Supplementary Table S15.", [t_s15a, t_s15b]),
    ("Supplementary Table S16.", [t_s16a, t_s16b, t_s16c]),
    ("Supplementary Table S17.", [t_s17]),
    ("Supplementary Table S19.", [t_s19]),
    ("Supplementary Table S20.", [t_s20]),
]


# ---------------------------------------------------------------------------
# docx plumbing
# ---------------------------------------------------------------------------
def locate(doc):
    """caption prefix -> [Table, ...] in document order."""
    out, current = {}, None
    for child in doc.element.body.iterchildren():
        if child.tag == qn("w:p"):
            text = Paragraph(child, doc).text.strip()
            if text.startswith("Supplementary Table S") or text.startswith("Supplementary Note S"):
                current = (text.split(" ")[2].rstrip(".") if text.startswith("Supplementary Table")
                           else None)
                if current:
                    current = "Supplementary Table " + current
        elif child.tag == qn("w:tbl") and current:
            out.setdefault(current, []).append(Table(child, doc))
    return out


def _set_cell(cell, text):
    """Replace a cell's text keeping the first run's formatting."""
    p = cell.paragraphs[0]
    for extra in cell.paragraphs[1:]:
        extra._p.getparent().remove(extra._p)
    runs = p.runs
    if not runs:
        p.add_run(text)
        return
    runs[0].text = text
    for r in runs[1:]:
        r._r.getparent().remove(r._r)


def _grid_widths(table):
    return [gc.get(qn("w:w")) for gc in table._tbl.tblGrid.findall(qn("w:gridCol"))]


def _add_columns(table, n_new):
    """Append n_new columns, cloning the last cell of every row (formatting incl.)."""
    grid = table._tbl.tblGrid
    last = grid.findall(qn("w:gridCol"))[-1]
    for _ in range(n_new):
        grid.append(copy.deepcopy(last))
    for tr in table._tbl.findall(qn("w:tr")):
        tcs = tr.findall(qn("w:tc"))
        for _ in range(n_new):
            tr.append(copy.deepcopy(tcs[-1]))


def fill(table, header, rows):
    """Header row kept (text refreshed); data rows rebuilt from the first data row."""
    trs = table._tbl.findall(qn("w:tr"))
    assert len(trs) >= 2, "table needs a header and one data row to serve as template"
    n_have = len(trs[0].findall(qn("w:tc")))
    if len(header) > n_have:
        _add_columns(table, len(header) - n_have)
        trs = table._tbl.findall(qn("w:tr"))
    assert len(header) == len(trs[0].findall(qn("w:tc"))), (len(header), n_have)
    template = copy.deepcopy(trs[1])
    for tr in trs[1:]:
        table._tbl.remove(tr)
    for r in rows:
        assert len(r) == len(header), (len(r), len(header), r)
        table._tbl.append(copy.deepcopy(template))
    tbl = Table(table._tbl, table._parent)
    for c, h in zip(tbl.rows[0].cells, header):
        _set_cell(c, h)
    for row, vals in zip(tbl.rows[1:], rows):
        for c, v in zip(row.cells, vals):
            _set_cell(c, v)


def compare(table, header, rows, key):
    have = [[c.text.strip() for c in r.cells] for r in table.rows]
    want = [header, *rows]
    bad = []
    ncol = min(len(have[0]), len(header))
    if len(have) != len(want):
        bad.append(f"row count {len(have)} in document vs {len(want)} regenerated")
    for i, (h, w) in enumerate(zip(have, want)):
        for j in range(ncol):
            if i == 0:
                continue  # headers may be renamed deliberately
            if h[j] == w[j]:
                continue
            # the same value printed at another precision (the v1 tables were pasted as the
            # CSV's own strings, the v2 tables are formatted): equal as numbers is equal
            if is_number(h[j]) and is_number(w[j]) and abs(float(h[j].replace("−", "-")) - float(w[j].replace("−", "-"))) < 1e-9:
                continue
            bad.append(f"r{i}c{j}: doc {h[j]!r} vs {w[j]!r}")
    if len(header) > len(have[0]):
        bad.append(f"{len(header) - len(have[0])} new column(s) not in the document")
    return bad


def insert_table_after(doc, after_table, title, caption, header, rows, template_table):
    """A new captioned table cloned from template_table's formatting, placed after after_table."""
    tbl = copy.deepcopy(template_table._tbl)
    t = Table(tbl, doc)
    p_title = copy.deepcopy(doc.paragraphs[0]._p)
    p_cap = copy.deepcopy(doc.paragraphs[0]._p)
    after_table._tbl.addnext(p_title)
    p_title.addnext(p_cap)
    p_cap.addnext(tbl)
    for el, text in ((p_title, title), (p_cap, caption)):
        para = Paragraph(el, doc)
        for r in para.runs[1:]:
            r._r.getparent().remove(r._r)
        if para.runs:
            para.runs[0].text = text
        else:
            para.add_run(text)
    # shrink/grow the cloned grid to the header width before filling
    n_have = len(tbl.findall(qn("w:tr"))[0].findall(qn("w:tc")))
    if len(header) < n_have:
        grid = tbl.tblGrid
        for gc in grid.findall(qn("w:gridCol"))[len(header):]:
            grid.remove(gc)
        for tr in tbl.findall(qn("w:tr")):
            for tc in tr.findall(qn("w:tc"))[len(header):]:
                tr.remove(tc)
    fill(t, header, rows)
    return t


S20_TITLE = "Supplementary Table S20. LR+ at the 95%-specificity operating point, frozen-matrix-v1 → v2: score change versus threshold placement"
S20_CAPTION = ("For every object × condition × calibration stage, the LR+ at the scanned observed-value threshold "
               "under the published frozen-matrix-v1 and under frozen-matrix-v2, the specificity the v2 threshold "
               "achieves, the LR+ interpolated to exactly 95% specificity under each version, and the decomposition "
               "of the scanned-threshold change into score change (the interpolated Δ) and threshold placement "
               "(the remainder). Placement within v2 is the scanned value minus the interpolated value under v2: "
               "at exactly 95% specificity LR+ is bounded by 1/0.05 = 20, so a scanned value above 20 — AlphaGenome's "
               "24.34 with BRCA1 excluded after calibration, interpolated 14.03, placement +10.31 — is the achieved "
               "specificity running past the nominal one, not stronger evidence. Both matrices are run through the "
               "identical phase-5 code path. 'n/e' = not evaluable. Source: phase1/reports/phase1/"
               "phase8_lr_placement_decomposition.csv (src/phase8_lr_decomposition.py).")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="src", required=True)
    ap.add_argument("--out", dest="dst")
    ap.add_argument("--reports", default="reports/phase1")
    ap.add_argument("--verify", action="store_true")
    a = ap.parse_args()
    rep = Path(a.reports)
    doc = docx.Document(a.src)
    found = locate(doc)
    total_bad = 0
    for prefix, builders in TABLES:
        key = prefix.rstrip(".")
        tables = found.get(key, [])
        if not tables:
            if a.verify or prefix != "Supplementary Table S20.":
                print(f"[supp] {key}: not in document" + ("" if a.verify else " (skipped)"))
            if prefix == "Supplementary Table S20." and not a.verify:
                header, rows = t_s20(rep)
                after = found["Supplementary Table S19"][0]
                insert_table_after(doc, after, S20_TITLE, S20_CAPTION, header, rows, after)
                print(f"[supp] {key}: inserted after S19 ({len(rows)} rows)")
            continue
        for k, (build, table) in enumerate(zip(builders, tables)):
            header, rows = build(rep)
            tag = key + ("" if len(builders) == 1 else f" ({'abc'[k]})")
            if a.verify:
                bad = compare(table, header, rows, tag)
                total_bad += len(bad)
                print(f"[verify] {tag}: {len(rows)} rows, {len(bad)} mismatch(es)")
                for b in bad[:8]:
                    print("    " + b)
            else:
                fill(table, header, rows)
                print(f"[supp] {tag}: {len(rows)} rows written")
    if a.verify:
        print(f"[verify] total mismatches: {total_bad}")
        return 1 if total_bad else 0
    doc.save(a.dst or a.src)
    print(f"[supp] wrote {a.dst or a.src}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
