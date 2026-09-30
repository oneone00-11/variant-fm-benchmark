"""E15 -- practical metrics: sensitivity, specificity and ROC at the evidence thresholds.

The likelihood ratio is the quantity to compare predictors on and to assign an evidence
tier with, but it does not tell a user what to expect from a tool in the terms a
laboratory reads first: how many damaging variants a threshold catches, how many normal
ones it flags, how many variants get any evidence at all, and how the tool behaves
across every threshold. This stage puts those classical measures beside the tiers, at
the thresholds the study has already fitted. Nothing is fitted here. Every threshold is
read from E2's configuration or from E3's tables, and no output of E1-E14 changes.

Outputs, under reports/evidence/ unless named otherwise:

  practical_metrics_by_threshold.csv   one row per tool x stratum x threshold kind x
      threshold source: counts; sensitivity and specificity with Wilson 95% intervals;
      coverage; prevalence; the observed PPV; the band likelihood ratio; the posterior
      at a prior of 0.10; and two tiers, side by side
  practical_metrics_logo_by_gene.csv   the held-out rows, one row per fold
  practical_metrics_external.csv       DDX3X and TP53 at the thresholds E7 carried over
  practical_roc_curves.csv, practical_pr_curves.csv   every distinct score threshold
  practical_curve_auc.csv              pooled AUROC and AUPRC of each curve, beside
      the gene-pooled AUROC of fig_data/fig3_evidence_by_territory.csv
  figures/figure_practical.pdf, .png   ROC and precision-recall curves, tiers marked
  ../docs/practical-metrics-summary.md the one-page table, the sentences and the legend

The columns are E3's: the thirteen score columns E1 puts in the panel and the fusion.
alphagenome_v061 is left out, as in E2-E14: it is not a panel column, scores 3,242 of
the 8,853 variants, and E3 fits no threshold on it.

Definitions, fixed here so that no column can be read two ways:

  * Labels, strata and score orientation are evid_common's: y_assay (1 damaging,
    0 normal), and every column larger-is-damaging (feature_orientation.csv).
  * PP3 side: a variant is called damaging when its score is AT OR ABOVE the
    threshold, the direction E3's thresholds are defined in (Pejaver's Eq. 8).
    sensitivity = damaging variants called / damaging variants; specificity = normal
    variants not called / normal variants.
  * BP4 side: a variant is called benign when its score is AT OR BELOW the threshold.
    bp4_sensitivity = normal variants called benign / normal variants; bp4_specificity
    = damaging variants not called benign / damaging variants. The prefix keeps them
    apart from the PP3 columns. coverage, coverage_labelled, n_scored_in_band, the
    posterior and tier_by_band_lr serve both sides; `side` says which band they count.
  * n, n_pos and n_neg count labelled variants, and sensitivity, specificity,
    prevalence and PPV are over them. n_scored counts every scored variant in the
    stratum, labelled or not, and coverage is the share of them in the band, because a
    user applying a tool does not know the label: it is the share of the stratum's
    variants that receive the evidence at all. coverage_labelled is the same share
    among labelled variants, the denominator E2's frac_pp3 uses.
  * band_lr_plus is E2's band ratio, P(score >= t | damaging) / P(score >= t |
    normal) = sensitivity / (1 - specificity), with evid_common.band_lr's rules: the
    denominator bounded at 1/(n_neg + 1), and no ratio for a band holding fewer than
    ten labelled variants. bp4_band_lr is the ratio of the band at or below t; its tier
    is read off the rule-of-three bound where that band holds no damaging variant, as
    E2 reads it, while the posterior is E2's, from the point estimate.
  * tier_by_interval_rule is the tier E3's interval rule gives the threshold: for a
    fitted threshold, the tier it was fitted for; for a published cut point, the tier
    E3's in-sample thresholds give a variant scoring exactly the cut point.
    tier_by_band_lr is the tier the band's ratio reaches against the same cuts. The
    first is a bound on every score above the threshold, the second an average over
    the band, and they can differ; both are kept.
  * threshold_source is `published` for Walker's cut points (0.2 and 0.1, applied to
    the published-basis SpliceAI column only, as everywhere in this study), `in_sample`
    for E3's thresholds fitted on the whole stratum and `logo` for E3's
    leave-one-gene-out thresholds, fitted on the other genes (six, or five at 11-50 bp,
    where BRCA2 has no variant). A logo row applies each fold's threshold to that
    fold's held-out gene and ADDS the counts over the folds whose threshold reached the
    tier, so its sensitivity is a cross-validated proportion of variants, not a mean of
    proportions. Those folds include ones whose band holds too few variants for a
    likelihood ratio, which Table 3 counts as not evaluable; folds_with_band_lr says
    how many had one. The combined Atlas score's BRCA1 and RAD51C folds are left out,
    as everywhere held-out genes are counted (its model was selected on those assays),
    and the fusion's logo rows carry a note adapted from E3's leakage note. Where a
    logo row has no counts (status other than ok) its n columns are empty rather than
    the stratum's.
  * Wilson intervals treat variants as independent. The seven genes are clusters, so
    the intervals understate the spread between genes; the by-gene table shows it.

E3 writes per-fold thresholds for the PP3 side only, so the BP4 kinds have no held-out
numbers. Their logo rows are written with a status that says so, not dropped.

The curves pool the seven genes: a user applies one threshold to every gene, and the
pooled curve is that use. The gene-pooled AUROC of the study (per-gene AUROCs combined
on the logit scale with random effects; Methods, Supplementary Table S10), which
fig_data/fig3_evidence_by_territory.csv carries for all variants, leaves out how the
genes' scores and prevalences differ. The AUC table carries both, and splits the pooled
AUROC into its same-gene and cross-gene pairs so the difference can be read off rather
than argued. Each curve table starts at threshold +inf with nothing called: (0, 0) on
the ROC, and on the precision-recall curve the conventional (recall 0, precision 1),
which is not an observation and is not drawn.

Run (PYTHONPATH=phase1):  python -m src.evid_practical_metrics
"""
from __future__ import annotations

import math
from fractions import Fraction
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import rankdata
from sklearn.metrics import (average_precision_score, precision_recall_curve,
                             roc_auc_score, roc_curve)

from . import evid_common as K
from . import evid_external as X
from . import evid_figures as F

REPORT_DIR = Path("reports/evidence")
CONFIG_PATH = Path("config/walker2023.yaml")
OUT_TABLE = REPORT_DIR / "practical_metrics_by_threshold.csv"
OUT_BY_GENE = REPORT_DIR / "practical_metrics_logo_by_gene.csv"
OUT_EXTERNAL = REPORT_DIR / "practical_metrics_external.csv"
OUT_ROC = REPORT_DIR / "practical_roc_curves.csv"
OUT_PR = REPORT_DIR / "practical_pr_curves.csv"
OUT_AUC = REPORT_DIR / "practical_curve_auc.csv"
OUT_MD = Path("../docs/practical-metrics-summary.md")
FIGURE = "figure_practical"
GENE_POOLED = REPORT_DIR / "fig_data" / "fig3_evidence_by_territory.csv"

# the four columns of the main tables; every other E3 column is carried as well
MAIN_TOOLS = ["spliceai_walker", "pangolin", "alphagenome", "avi"]
CURVE_STRATA = ["s3_10", "s11_50", "s3_50"]
# pm12 is computed and not drawn: the canonical dinucleotides are outside the
# recommendation, and the tables mark the pool the way the change log does
TABLE_STRATA = CURVE_STRATA + ["pm12"]
OUT_OF_SCOPE = "out_of_scope_incl_pm12"
PP3_TIERS = ["supporting", "moderate", "strong", "very_strong"]
BP4_TIERS = ["supporting", "moderate"]
EXTERNAL_TIERS = ["supporting", "moderate", "strong"]        # the tiers E7 carries over
TIER_NAME = {"supporting": "Supporting", "moderate": "Moderate", "strong": "Strong",
             "very_strong": "Very strong"}
# DDX3X under the deposit's own classification, TP53 under its control-anchored
# construction: each gene's primary label, as in Figures 1 and 4
EXTERNAL = [("DDX3X", F.PRIMARY_LABEL["DDX3X"]), ("TP53", F.PRIMARY_LABEL["TP53"])]
PRIOR = 0.10                   # the prior Tavtigian's tier cuts are solved at
Z95 = 1.959963984540054        # two-sided 95% quantile of the standard normal
MATCHED_SPECIFICITY = (Fraction(90, 100), Fraction(95, 100))
BP4_NO_LOGO = ("no held-out threshold: E3 writes per-fold thresholds for the PP3 "
               "side only")
LEAKAGE = ("adapted from E3: fusion out-of-fold scores in the training folds were "
           "produced by models trained on the held-out gene; held-out values are "
           "optimistic")
AVI_EXCLUDED = ("the combined Atlas score's model was selected on this assay; left out "
                "of the held-out counts")

TABLE_COLUMNS = [
    "tool", "main_table", "stratum", "scope", "side", "threshold_kind",
    "threshold_source", "status", "note", "in_sample_reached", "threshold",
    "threshold_fold_min", "threshold_fold_max", "folds_total", "folds_excluded",
    "folds_used", "folds_with_band_lr", "tier_by_interval_rule", "tier_by_band_lr",
    "n_scored", "n", "n_pos", "n_neg", "prevalence", "n_scored_in_band", "coverage",
    "coverage_labelled", "n_pos_above", "n_neg_above", "n_pos_below", "n_neg_below",
    "sensitivity", "sensitivity_lo", "sensitivity_hi",
    "specificity", "specificity_lo", "specificity_hi", "ppv_observed", "band_lr_plus",
    "bp4_n_neg_at_or_below", "bp4_n_pos_at_or_below", "bp4_n_neg_above",
    "bp4_n_pos_above", "bp4_sensitivity", "bp4_sensitivity_lo", "bp4_sensitivity_hi",
    "bp4_specificity", "bp4_specificity_lo", "bp4_specificity_hi", "bp4_npv_observed",
    "bp4_band_lr", "posterior_at_prior_0.10",
]
PP3_METRICS = ["n_scored_in_band", "coverage", "coverage_labelled", "n_pos_above",
               "n_neg_above", "n_pos_below", "n_neg_below", "sensitivity",
               "sensitivity_lo", "sensitivity_hi", "specificity", "specificity_lo",
               "specificity_hi", "ppv_observed", "band_lr_plus", "posterior_at_prior_0.10",
               "tier_by_band_lr"]
N_COLUMNS = ["n_scored", "n", "n_pos", "n_neg", "prevalence"]


def _read(path: Path) -> pd.DataFrame:
    """A table as written. pandas' default float parser can land one binary digit away
    from the value written, which moves a threshold enough to drop the variants that
    sit exactly on it (every threshold here is an observed score) from a band that
    includes its boundary."""
    return pd.read_csv(path, float_precision="round_trip")


# ---------------------------------------------------------------------------
# counts and the quantities read off them
# ---------------------------------------------------------------------------
def _cell(frame: pd.DataFrame, tool: str, label: str, thr: float) -> dict:
    """Counts at one threshold on both sides: the PP3 band (score >= thr) and the BP4
    band (score <= thr), over labelled variants, and over every scored variant for
    coverage. The complement of each band is counted, not inferred."""
    s_all = frame[tool].to_numpy(dtype=float)
    y_all = frame[label].to_numpy(dtype=float)
    scored = np.isfinite(s_all)
    lab = scored & np.isfinite(y_all)
    s, y = s_all[lab], y_all[lab]
    pos, neg = s[y == 1], s[y == 0]
    sc = s_all[scored]
    return {
        "n_scored": int(scored.sum()), "n": int(lab.sum()),
        "n_pos": int(len(pos)), "n_neg": int(len(neg)),
        "ge_pos": int((pos >= thr).sum()), "lt_pos": int((pos < thr).sum()),
        "ge_neg": int((neg >= thr).sum()), "lt_neg": int((neg < thr).sum()),
        "le_pos": int((pos <= thr).sum()), "gt_pos": int((pos > thr).sum()),
        "le_neg": int((neg <= thr).sum()), "gt_neg": int((neg > thr).sum()),
        "ge_all": int((sc >= thr).sum()), "le_all": int((sc <= thr).sum()),
    }


def _add(cells: list[dict]) -> dict:
    return {k: int(sum(c[k] for c in cells)) for k in cells[0]}


def _describe(frame: pd.DataFrame, tool: str, label: str) -> dict:
    """A population's own numbers, which a row carries beside its metrics."""
    c = _cell(frame, tool, label, np.inf)
    return _n_columns(c)


def _n_columns(c: dict) -> dict:
    return {"n_scored": c["n_scored"], "n": c["n"], "n_pos": c["n_pos"],
            "n_neg": c["n_neg"], "prevalence": _ratio(c["n_pos"], c["n"])}


def _ratio(k: int, n: int) -> float:
    return k / n if n else np.nan


def _wilson(k: int, n: int) -> tuple[float, float]:
    """Wilson score interval, 95%. At k = 0 the lower end is 0 and at k = n the upper
    end is 1, exactly: the closed form gives them, the floating-point sum does not."""
    if n <= 0:
        return np.nan, np.nan
    p, z2 = k / n, Z95 * Z95
    denom = 1.0 + z2 / n
    centre = (p + z2 / (2 * n)) / denom
    half = Z95 * math.sqrt(p * (1.0 - p) / n + z2 / (4.0 * n * n)) / denom
    lo = 0.0 if k == 0 else max(0.0, centre - half)
    hi = 1.0 if k == n else min(1.0, centre + half)
    return lo, hi


def _band_lr(k_pos: int, n_pos: int, k_neg: int, n_neg: int) -> float:
    """evid_common.band_lr from counts: the same gates and the same bounded
    denominator, so a row here and the ratio E2 or E7 reports for the same band are
    the same number."""
    if n_pos < K.MIN_POS or n_neg < K.MIN_NEG or k_pos + k_neg < K.MIN_BAND:
        return np.nan
    p = k_pos / n_pos
    if p == 0.0:
        return 0.0
    return p / max(k_neg / n_neg, 1.0 / (n_neg + 1))


def _benign_bound(k_pos: int, n_pos: int, k_neg: int, n_neg: int) -> float:
    """evid_common.band_lr_benign_bound from counts: the value a BP4 tier is read off."""
    lr = _band_lr(k_pos, n_pos, k_neg, n_neg)
    if np.isfinite(lr) and k_pos == 0:
        return K.rule_of_three_lr(n_pos, k_neg / n_neg, n_neg)
    return lr


def _posterior(lr: float) -> float:
    """The posterior probability of pathogenicity at Tavtigian's prior, as E2 writes it."""
    if not np.isfinite(lr):
        return np.nan
    return float(lr * PRIOR / ((lr - 1.0) * PRIOR + 1.0))


def _pp3_values(c: dict, bands: dict) -> dict:
    lr = _band_lr(c["ge_pos"], c["n_pos"], c["ge_neg"], c["n_neg"])
    sens_lo, sens_hi = _wilson(c["ge_pos"], c["n_pos"])
    spec_lo, spec_hi = _wilson(c["lt_neg"], c["n_neg"])
    return {
        "n_scored_in_band": c["ge_all"], "coverage": _ratio(c["ge_all"], c["n_scored"]),
        "coverage_labelled": _ratio(c["ge_pos"] + c["ge_neg"], c["n"]),
        "n_pos_above": c["ge_pos"], "n_neg_above": c["ge_neg"],
        "n_pos_below": c["lt_pos"], "n_neg_below": c["lt_neg"],
        "sensitivity": _ratio(c["ge_pos"], c["n_pos"]),
        "sensitivity_lo": sens_lo, "sensitivity_hi": sens_hi,
        "specificity": _ratio(c["lt_neg"], c["n_neg"]),
        "specificity_lo": spec_lo, "specificity_hi": spec_hi,
        "ppv_observed": _ratio(c["ge_pos"], c["ge_pos"] + c["ge_neg"]),
        "band_lr_plus": lr, "posterior_at_prior_0.10": _posterior(lr),
        "tier_by_band_lr": K.tier_of(lr, bands, "pathogenic"),
    }


def _bp4_values(c: dict, bands: dict) -> dict:
    lr = _band_lr(c["le_pos"], c["n_pos"], c["le_neg"], c["n_neg"])
    bound = _benign_bound(c["le_pos"], c["n_pos"], c["le_neg"], c["n_neg"])
    sens_lo, sens_hi = _wilson(c["le_neg"], c["n_neg"])
    spec_lo, spec_hi = _wilson(c["gt_pos"], c["n_pos"])
    return {
        "n_scored_in_band": c["le_all"], "coverage": _ratio(c["le_all"], c["n_scored"]),
        "coverage_labelled": _ratio(c["le_pos"] + c["le_neg"], c["n"]),
        "bp4_n_neg_at_or_below": c["le_neg"], "bp4_n_pos_at_or_below": c["le_pos"],
        "bp4_n_neg_above": c["gt_neg"], "bp4_n_pos_above": c["gt_pos"],
        "bp4_sensitivity": _ratio(c["le_neg"], c["n_neg"]),
        "bp4_sensitivity_lo": sens_lo, "bp4_sensitivity_hi": sens_hi,
        "bp4_specificity": _ratio(c["gt_pos"], c["n_pos"]),
        "bp4_specificity_lo": spec_lo, "bp4_specificity_hi": spec_hi,
        "bp4_npv_observed": _ratio(c["le_neg"], c["le_neg"] + c["le_pos"]),
        "bp4_band_lr": lr, "posterior_at_prior_0.10": _posterior(lr),
        "tier_by_band_lr": K.tier_of(bound, bands, "benign"),
    }


def _values(c: dict, side: str, bands: tuple[dict, dict]) -> dict:
    return _pp3_values(c, bands[0]) if side == "pp3" else _bp4_values(c, bands[1])


def _interval_tier(e: pd.DataFrame, side: str, score: float) -> str:
    """The tier E3's in-sample thresholds give a variant scoring `score`: the highest
    tier whose threshold it clears, at or above it on the PP3 side and at or below it
    on the BP4 side."""
    col = "pp3_threshold_insample" if side == "pp3" else "bp4_threshold_insample"
    best = "below supporting"
    for tier in PP3_TIERS:
        if tier not in e.index or e.loc[tier, "status"] != "ok":
            return "not evaluable"
        tau = float(e.loc[tier, col])
        if np.isfinite(tau) and (score >= tau if side == "pp3" else score <= tau):
            best = TIER_NAME[tier]
    return best


# ---------------------------------------------------------------------------
# E4.1 and E4.2: the seven genes, at every threshold, in-sample and held out
# ---------------------------------------------------------------------------
def by_threshold(df: pd.DataFrame, ev: pd.DataFrame, folds: pd.DataFrame, cfg: dict,
                 tools: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    bands = K.acmg_bands(cfg)
    cuts = (("pp3", "walker_pp3", float(cfg["thresholds"]["pp3"]["value"])),
            ("bp4", "walker_bp4", float(cfg["thresholds"]["bp4"]["value"])))
    folds = folds[folds["side"] == "pp3"]
    no_counts = {c: np.nan for c in N_COLUMNS}
    rows, gene_rows = [], []
    for tool in tools:
        for stratum in TABLE_STRATA:
            frame = K.stratum_frame(df, stratum)
            e = ev[(ev.tool == tool) & (ev.stratum == stratum)].set_index("tier")
            base = {"tool": tool, "main_table": tool in MAIN_TOOLS, "stratum": stratum,
                    "scope": "in_scope" if stratum in K.IN_SCOPE_STRATA else OUT_OF_SCOPE}
            whole = _describe(frame, tool, "y_assay")
            if tool == "spliceai_walker":
                for side, kind, cut in cuts:
                    c = _cell(frame, tool, "y_assay", cut)
                    rows.append(base | whole | {
                        "side": side, "threshold_kind": kind,
                        "threshold_source": "published", "status": "ok", "note": "",
                        "threshold": cut,
                        "tier_by_interval_rule": _interval_tier(e, side, cut)}
                        | _values(c, side, bands))
            for side, tiers in (("pp3", PP3_TIERS), ("bp4", BP4_TIERS)):
                for tier in tiers:
                    kind = tier if side == "pp3" else f"bp4_{tier}"
                    evaluable = tier in e.index and e.loc[tier, "status"] == "ok"
                    reason = ("" if evaluable else "E3 has no row here"
                              if tier not in e.index else str(e.loc[tier, "reason"]))
                    tau = (float(e.loc[tier, f"{side}_threshold_insample"])
                           if evaluable else np.nan)
                    # carried on the held-out row too: a fold can reach a tier the whole
                    # stratum does not, and E7 and Table 3 do not use such a tier
                    head = base | {"side": side, "threshold_kind": kind,
                                   "in_sample_reached": bool(np.isfinite(tau))}
                    # in-sample: E3's threshold on the whole stratum
                    row = head | whole | {"threshold_source": "in_sample"}
                    if not evaluable:
                        row |= {"status": "not evaluable", "note": reason}
                    elif not np.isfinite(tau):
                        row |= {"status": "not reached", "note": ""}
                    else:
                        row |= {"status": "ok", "note": "", "threshold": tau,
                                "tier_by_interval_rule": TIER_NAME[tier]}
                        row |= _values(_cell(frame, tool, "y_assay", tau), side, bands)
                    rows.append(row)
                    # held out: each fold's threshold on its held-out gene
                    row = head | no_counts | {"threshold_source": "logo"}
                    if side == "bp4":
                        rows.append(row | {"status": "no held-out threshold",
                                           "note": BP4_NO_LOGO})
                        continue
                    if not evaluable:
                        rows.append(row | {"status": "not evaluable", "note": reason})
                        continue
                    f = folds[(folds.tool == tool) & (folds.stratum == stratum)
                              & (folds.tier == tier)].sort_values("heldout_gene")
                    excluded = sorted(g for g in f.heldout_gene
                                      if tool == "avi" and g in K.AVI_SEEN_IN_TRAINING)
                    used = []
                    for fr in f.itertuples():
                        tau_g = float(fr.threshold_from_6_genes)
                        gframe = frame[frame["gene"].astype(str) == fr.heldout_gene]
                        counted = fr.heldout_gene not in excluded and np.isfinite(tau_g)
                        g = {"tool": tool, "stratum": stratum, "scope": base["scope"],
                             "tier": tier, "heldout_gene": fr.heldout_gene,
                             "threshold_from_6_genes": tau_g,
                             "e3_heldout_lr": fr.heldout_lr,
                             "counted_in_logo_row": counted,
                             "status": ("ok" if np.isfinite(tau_g)
                                        else "not reached in this fold"),
                             "note": (AVI_EXCLUDED if fr.heldout_gene in excluded
                                      else LEAKAGE if tool == K.FUSION else "")}
                        g |= _describe(gframe, tool, "y_assay")
                        if np.isfinite(tau_g):
                            c = _cell(gframe, tool, "y_assay", tau_g)
                            g |= _pp3_values(c, bands[0])
                            if counted:
                                used.append((tau_g, c, g["band_lr_plus"]))
                        gene_rows.append(g)
                    row |= {"folds_total": len(f), "folds_excluded": "; ".join(excluded),
                            "folds_used": len(used),
                            "folds_with_band_lr": sum(np.isfinite(lr) for _, _, lr in used),
                            "note": LEAKAGE if tool == K.FUSION else ""}
                    if not len(f):
                        rows.append(row | {"status": "no folds",
                                           "note": "E3 ran no leave-one-gene-out fold here"})
                        continue
                    if not used:
                        rows.append(row | {"status": "not reached"})
                        continue
                    c = _add([c for _, c, _ in used])
                    taus = [t for t, _, _ in used]
                    # the held-out genes' own numbers, over the folds counted
                    row |= {"status": "ok", "threshold_fold_min": min(taus),
                            "threshold_fold_max": max(taus),
                            "tier_by_interval_rule": TIER_NAME[tier]} | _n_columns(c)
                    rows.append(row | _pp3_values(c, bands[0]))
    table = pd.DataFrame(rows).reindex(columns=TABLE_COLUMNS)
    lead = ["tool", "stratum", "scope", "tier", "heldout_gene", "threshold_from_6_genes",
            "e3_heldout_lr", "counted_in_logo_row", "status", "note"] + N_COLUMNS
    by_gene = pd.DataFrame(gene_rows).reindex(columns=lead + PP3_METRICS)
    return table, by_gene


# ---------------------------------------------------------------------------
# E4.3: the two genes outside the seven, at the thresholds carried over
# ---------------------------------------------------------------------------
def external(table: pd.DataFrame, cfg: dict, tools: list[str]) -> pd.DataFrame:
    """At the thresholds E7 applied, read from E7's own tables. Where E7 carried no
    threshold, the row takes E7's words: its tier label as the status and its note."""
    bands = K.acmg_bands(cfg)
    pp3_cut = float(cfg["thresholds"]["pp3"]["value"])
    seven = table.set_index(["tool", "stratum", "threshold_kind", "threshold_source"])
    rows = []
    for gene, label in EXTERNAL:
        data, _ = X._load_external(gene)
        ext = _read(REPORT_DIR / f"external_{gene.lower()}.csv")
        ext = ext[ext.label_definition == label]
        basis = _read(REPORT_DIR / f"external_{gene.lower()}_column_basis.csv")
        basis = basis.set_index("tool")
        for tool in [t for t in tools if t in set(ext.tool)]:
            for stratum in TABLE_STRATA:
                r = ext[(ext.tool == tool) & (ext.stratum == stratum)]
                if not len(r):
                    continue
                r = r.iloc[0]
                frame = K.stratum_frame(data, stratum)
                base = {"gene": gene, "label_definition": label, "tool": tool,
                        "main_table": tool in MAIN_TOOLS, "stratum": stratum,
                        "scope": ("in_scope" if stratum in K.IN_SCOPE_STRATA
                                  else OUT_OF_SCOPE),
                        "external_test": (gene, tool) not in F.NOT_EXTERNAL}
                base |= _describe(frame, tool, label)
                kinds = ([("walker_pp3", "published", pp3_cut)]
                         if tool == "spliceai_walker" else [])
                kinds += [(t, "logo_median", r.get(f"e3_{t}_threshold", np.nan))
                          for t in EXTERNAL_TIERS]
                for kind, source, thr in kinds:
                    row = base | {"threshold_kind": kind, "threshold_source": source}
                    note = ("" if row["external_test"] else
                            "not an external test: the combined Atlas score's model was "
                            "selected on the DDX3X assay")
                    if r.status != "ok":
                        rows.append(row | {"status": r.status, "note": r.reason})
                        continue
                    if source == "logo_median" and not bool(r.column_basis_comparable):
                        rows.append(row | {"status": "not carried over",
                                           "note": "column not comparable: "
                                                   + str(basis.loc[tool, "reason"])})
                        continue
                    thr = float(thr) if thr is not None else np.nan
                    if not np.isfinite(thr):
                        label_here = r.get(f"e3_{kind}_tier_here", np.nan)
                        why = r.get(f"e3_{kind}_note", np.nan)
                        rows.append(row | {
                            "status": (label_here if isinstance(label_here, str)
                                       else "not carried over"),
                            "note": why if isinstance(why, str)
                            else "no threshold carried over"})
                        continue
                    c = _cell(frame, tool, label, thr)
                    row |= {"status": "ok", "note": note, "threshold": thr}
                    row |= _pp3_values(c, bands[0])
                    # the same tool, band and threshold kind on the seven genes: held
                    # out for a fitted threshold, the pooled value for the cut point
                    key = (tool, stratum, kind, "published" if source == "published"
                           else "logo")
                    if key in seven.index and seven.loc[key, "status"] == "ok":
                        s7 = seven.loc[key]
                        row |= {"seven_gene_source": key[3],
                                "seven_gene_sensitivity": s7["sensitivity"],
                                "seven_gene_specificity": s7["specificity"],
                                "sensitivity_minus_seven_gene":
                                    row["sensitivity"] - s7["sensitivity"],
                                "specificity_minus_seven_gene":
                                    row["specificity"] - s7["specificity"]}
                    rows.append(row)
    lead = ["gene", "label_definition", "tool", "main_table", "stratum", "scope",
            "external_test", "threshold_kind", "threshold_source", "status", "note",
            "threshold"] + N_COLUMNS
    tail = ["seven_gene_source", "seven_gene_sensitivity", "seven_gene_specificity",
            "sensitivity_minus_seven_gene", "specificity_minus_seven_gene"]
    return pd.DataFrame(rows).reindex(columns=lead + PP3_METRICS + tail)


# ---------------------------------------------------------------------------
# E4.4: the curves, the seven genes pooled
# ---------------------------------------------------------------------------
def _mann_whitney(y: np.ndarray, s: np.ndarray) -> tuple[float, int]:
    """U over the damaging-normal pairs, ties counted half, and the number of pairs."""
    r = rankdata(s, method="average")
    n_pos, n_neg = int((y == 1).sum()), int((y == 0).sum())
    return float(r[y == 1].sum() - n_pos * (n_pos + 1) / 2.0), n_pos * n_neg


def _sens_at_specificity(y: np.ndarray, s: np.ndarray, q: Fraction) -> tuple[float, float]:
    """(sensitivity, threshold) at the most sensitive threshold whose specificity is at
    least q, over the distinct scores; the comparison is made in integers."""
    order = np.argsort(-s, kind="mergesort")
    ss, yy = s[order], y[order]
    last = np.append(np.flatnonzero(np.diff(ss)), len(ss) - 1)
    tp, fp = np.cumsum(yy == 1)[last], np.cumsum(yy == 0)[last]
    n_pos, n_neg = int(tp[-1]), int(fp[-1])
    ok = fp * q.denominator <= (q.denominator - q.numerator) * n_neg
    if not ok.any():
        return 0.0, np.inf
    j = int(np.flatnonzero(ok)[-1])
    return tp[j] / n_pos, float(ss[last[j]])


def curves(df: pd.DataFrame, tools: list[str]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    pooled = _read(GENE_POOLED)
    pooled = pooled[pooled.arm == "all"].set_index(["tool", "stratum"])
    roc_rows, pr_rows, auc_rows = [], [], []
    for tool in tools:
        for stratum in CURVE_STRATA:
            sub = K.stratum_frame(df, stratum).dropna(subset=["y_assay", tool])
            y = sub["y_assay"].to_numpy(dtype=float).astype(int)
            s = sub[tool].to_numpy(dtype=float)
            genes = sub["gene"].astype(str).to_numpy()
            n_pos, n_neg = int((y == 1).sum()), int((y == 0).sum())
            # one point per distinct score, from "nothing called" (threshold +inf)
            fpr, tpr, thr = roc_curve(y, s, drop_intermediate=False)
            roc_rows.append(pd.DataFrame({"tool": tool, "stratum": stratum,
                                          "threshold": thr, "fpr": fpr, "tpr": tpr,
                                          "n_pos": n_pos, "n_neg": n_neg}))
            # sklearn returns these in increasing threshold, closed by the conventional
            # (recall 0, precision 1) point; written in the ROC table's order instead,
            # that point first at threshold +inf
            prec, rec, thr_pr = precision_recall_curve(y, s, drop_intermediate=False)
            pr_rows.append(pd.DataFrame({
                "tool": tool, "stratum": stratum,
                "threshold": np.concatenate([[np.inf], thr_pr[::-1]]),
                "recall": np.concatenate([[rec[-1]], rec[:-1][::-1]]),
                "precision": np.concatenate([[prec[-1]], prec[:-1][::-1]]),
                "n_pos": n_pos, "n_neg": n_neg}))
            # The pooled AUROC split into its damaging-normal pairs: pairs within one
            # gene, whose AUROC is the per-gene AUROCs averaged with weights n_pos x
            # n_neg, and pairs across genes. The same-gene part is also given over the
            # genes the gene-pooled estimate keeps, so pooled minus gene-pooled reads as
            # three terms: cross-gene pairs, the genes that estimate leaves out, and its
            # own weighting (inverse variance, random effects, logit scale).
            u_all, pairs_all = _mann_whitney(y, s)
            u_same = pairs_same = u_kept = pairs_kept = 0
            dropped, dropped_pos, dropped_neg = [], 0, 0
            for g in sorted(set(genes)):
                m = genes == g
                u_g, p_g = _mann_whitney(y[m], s[m])
                u_same += u_g
                pairs_same += p_g
                gp, gn = int((y[m] == 1).sum()), int((y[m] == 0).sum())
                if gp < K.MIN_POS or gn < K.MIN_NEG:
                    dropped.append(g)
                    dropped_pos += gp
                    dropped_neg += gn
                else:
                    u_kept += u_g
                    pairs_kept += p_g
            row = {"tool": tool, "stratum": stratum, "n_pos": n_pos, "n_neg": n_neg,
                   "prevalence": n_pos / (n_pos + n_neg),
                   "auroc_pooled_curve": roc_auc_score(y, s),
                   "auprc_pooled_curve": average_precision_score(y, s),
                   "auroc_same_gene_pairs": _ratio(u_same, pairs_same),
                   "auroc_same_gene_pairs_gene_pool_genes": _ratio(u_kept, pairs_kept),
                   "auroc_cross_gene_pairs": _ratio(u_all - u_same, pairs_all - pairs_same),
                   "share_cross_gene_pairs": 1.0 - pairs_same / pairs_all}
            for q in MATCHED_SPECIFICITY:
                sens, t = _sens_at_specificity(y, s, q)
                row[f"sensitivity_at_specificity_{float(q):.2f}"] = sens
                row[f"threshold_at_specificity_{float(q):.2f}"] = t
            g = pooled.loc[(tool, stratum)] if (tool, stratum) in pooled.index else None
            row |= {
                "auroc_gene_pooled": g["auroc"] if g is not None else np.nan,
                "auroc_gene_pooled_lo": g["auroc_lo"] if g is not None else np.nan,
                "auroc_gene_pooled_hi": g["auroc_hi"] if g is not None else np.nan,
                "auprc_median_per_gene": g["prauc"] if g is not None else np.nan,
                "k_genes_gene_pooled": g["k_genes"] if g is not None else np.nan,
                "n_pos_gene_pooled": g["n_pos"] if g is not None else np.nan,
                "n_neg_gene_pooled": g["n_neg"] if g is not None else np.nan,
                # genes the gene-pooled estimate leaves out (fewer than ten on a side);
                # the pooled curve keeps them
                "genes_not_in_gene_pool": "; ".join(dropped),
                "n_pos_not_in_gene_pool": dropped_pos,
                "n_neg_not_in_gene_pool": dropped_neg,
            }
            row["auroc_pooled_minus_gene_pooled"] = (row["auroc_pooled_curve"]
                                                     - row["auroc_gene_pooled"])
            auc_rows.append(row)
    return (pd.concat(roc_rows, ignore_index=True), pd.concat(pr_rows, ignore_index=True),
            pd.DataFrame(auc_rows))


# ---------------------------------------------------------------------------
# E4.5: the figure
# ---------------------------------------------------------------------------
# Supplementary Figure S2's colours for these four columns: the only place the study
# gives each tool its own colour (Figure 3 draws every tool in one colour). The markers
# encode the tier here, so each curve takes a dash pattern as its second channel.
TOOL_COLOUR = {t: c for t, c, _ in F.DILUTION_TOOLS}
TOOL_DASH = {"spliceai_walker": "-", "pangolin": (0, (4.0, 1.4)),
             "alphagenome": (0, (1.1, 1.1)), "avi": (0, (4.0, 1.2, 1.1, 1.2))}
SHORT = {"spliceai_walker": "SpliceAI", "pangolin": "Pangolin",
         "alphagenome": "AlphaGenome", "avi": "Atlas combined"}
TIER_MARKER = {"supporting": ("o", 3.6), "moderate": ("s", 3.4), "strong": ("^", 4.0)}
# Two tools' thresholds can land on the same point; each tool's markers are drawn a
# fixed distance (in points) off the curve vertex so that neither hides the other
MARKER_NUDGE = {"spliceai_walker": (0.0, 1.3), "pangolin": (1.3, 0.0),
                "alphagenome": (0.0, -1.3), "avi": (-1.3, 0.0)}
# Every fitted threshold sits at a false-positive rate of a few per cent, where the
# full ROC cannot separate them; an inset draws that end of the curves, its range
# set from the thresholds themselves in steps of these sizes
INSET_X_STEP, INSET_Y_STEP = 0.02, 0.2
INSET_BOUNDS = [0.45, 0.42, 0.52, 0.30]      # axes fraction: above the AUROC key


def _point(roc: pd.DataFrame, pr: pd.DataFrame, tool: str, stratum: str,
           thr: float) -> tuple[float, float, float, float]:
    """(fpr, tpr, recall, precision) of "score >= thr" on the drawn curves: the vertex
    at the lowest curve threshold still at or above thr, which is thr itself when thr
    is an observed score."""
    r = roc[(roc.tool == tool) & (roc.stratum == stratum)]
    p = pr[(pr.tool == tool) & (pr.stratum == stratum)]
    j = int(np.flatnonzero(r.threshold.to_numpy() >= thr)[-1])
    k = int(np.flatnonzero(p.threshold.to_numpy() >= thr)[-1])
    return (float(r.fpr.iloc[j]), float(r.tpr.iloc[j]),
            float(p.recall.iloc[k]), float(p.precision.iloc[k]))


def _tier_marks(table: pd.DataFrame) -> pd.DataFrame:
    """The in-sample Supporting, Moderate and Strong thresholds of the drawn tools."""
    return table[(table.threshold_source == "in_sample") & (table.side == "pp3")
                 & (table.status == "ok") & table.tool.isin(MAIN_TOOLS)
                 & table.stratum.isin(CURVE_STRATA)
                 & table.threshold_kind.isin(list(TIER_MARKER))]


def _inset_range(roc: pd.DataFrame, marks: pd.DataFrame) -> tuple[float, float]:
    """One inset range for every panel, so the insets compare: the largest false-
    positive rate of a drawn threshold, and the largest sensitivity the drawn curves
    reach within it, each rounded up to its step."""
    xmax = math.ceil(float((1.0 - marks.specificity).max()) / INSET_X_STEP) * INSET_X_STEP
    r = roc[roc.tool.isin(MAIN_TOOLS) & roc.stratum.isin(CURVE_STRATA)
            & (roc.fpr <= xmax)]
    ymax = math.ceil(float(r.tpr.max()) / INSET_Y_STEP) * INSET_Y_STEP
    return round(xmax, 10), round(min(ymax, 1.0), 10)


def figure(roc: pd.DataFrame, pr: pd.DataFrame, auc: pd.DataFrame,
           table: pd.DataFrame, cfg: dict) -> None:
    import matplotlib.pyplot as plt
    import matplotlib.transforms as mtransforms
    from matplotlib.lines import Line2D
    F._style()
    pp3_cut = float(cfg["thresholds"]["pp3"]["value"])
    marks = _tier_marks(table)
    xmax, ymax = _inset_range(roc, marks)
    unit = [0, 0.25, 0.5, 0.75, 1]
    unit_labels = ["0", "0.25", "0.5", "0.75", "1"]
    fig, axes = plt.subplots(2, 3, figsize=(F.DOUBLE, 128 * F.MM),
                             gridspec_kw={"hspace": 0.42, "wspace": 0.28})

    def nudged(ax, tool):
        dx, dy = MARKER_NUDGE[tool]
        return mtransforms.offset_copy(ax.transData, fig=fig, x=dx, y=dy, units="points")

    for j, st in enumerate(CURVE_STRATA):
        top, bottom = axes[0, j], axes[1, j]
        inset = top.inset_axes(INSET_BOUNDS)
        top.plot([0, 1], [0, 1], color=F.GRID, lw=0.6, ls=(0, (1.5, 1.5)), zorder=0)
        a = auc[auc.stratum == st].set_index("tool")
        prevalence = float(a["prevalence"].iloc[0])
        bottom.axhline(prevalence, color=F.INK2, lw=0.7, ls=(0, (3, 1.5)), zorder=0)
        bottom.text(0.04, prevalence, f"prevalence {prevalence:.2f}", fontsize=5.8,
                    color=F.INK2, ha="left", va="bottom")
        handles = []
        for tool in MAIN_TOOLS:
            colour, dash = TOOL_COLOUR[tool], TOOL_DASH[tool]
            r = roc[(roc.tool == tool) & (roc.stratum == st)]
            p = pr[(pr.tool == tool) & (pr.stratum == st)]
            top.plot(r.fpr, r.tpr, color=colour, lw=1.0, ls=dash, zorder=2)
            inset.plot(r.fpr, r.tpr, color=colour, lw=0.8, ls=dash, zorder=2)
            # precision as the step function average precision sums; the (recall 0,
            # precision 1) convention in the table's first row is not an observation
            bottom.plot(p.recall.iloc[1:], p.precision.iloc[1:], color=colour, lw=1.0,
                        ls=dash, drawstyle="steps-pre", zorder=2)
            for m in marks[(marks.tool == tool) & (marks.stratum == st)].itertuples():
                x, yv, rc, pc = _point(roc, pr, tool, st, m.threshold)
                shape, size = TIER_MARKER[m.threshold_kind]
                style = dict(marker=shape, ms=size, color=colour, mec="white", mew=0.4,
                             ls="", zorder=4, clip_on=False)
                top.plot(x, yv, transform=nudged(top, tool), **style)
                inset.plot(x, yv, transform=nudged(inset, tool), **style)
                bottom.plot(rc, pc, transform=nudged(bottom, tool), **style)
            if tool == "spliceai_walker":
                x, yv, rc, pc = _point(roc, pr, tool, st, pp3_cut)
                style = dict(marker="D", ms=3.8, mfc="white", mec=colour, mew=0.9, ls="",
                             zorder=5)
                top.plot(x, yv, **style)
                if x <= xmax:
                    inset.plot(x, yv, **style)
                bottom.plot(rc, pc, **style)
            handles.append(Line2D([], [], color=colour, lw=1.0, ls=dash,
                                  label=f"{SHORT[tool]}  "
                                        f"{a.loc[tool, 'auroc_pooled_curve']:.2f}"))
        top.legend(handles=handles, loc="lower right", title="Pooled AUROC",
                   alignment="left", fontsize=5.8, title_fontsize=5.8, handlelength=2.2,
                   handletextpad=0.4, borderaxespad=0.2, labelspacing=0.25)
        inset.set_xlim(0, xmax)
        inset.set_ylim(0, ymax)
        xt = [round(v, 10) for v in np.arange(0, xmax + INSET_X_STEP / 2, INSET_X_STEP)]
        inset.set_xticks(xt)
        inset.set_xticklabels([f"{v:g}" for v in xt], fontsize=5.5)
        yt = [0, round(ymax / 2, 10), ymax]
        inset.set_yticks(yt)
        inset.set_yticklabels([f"{v:g}" for v in yt], fontsize=5.5)
        inset.tick_params(length=1.5, pad=1)
        for spine in inset.spines.values():
            spine.set_visible(True)
            spine.set_color(F.MUTED)
            spine.set_linewidth(0.5)
        for ax in (top, bottom):
            ax.set_xlim(-0.01, 1.01)
            ax.set_ylim(-0.01, 1.01)
            ax.set_xticks(unit)
            ax.set_yticks(unit)
            ax.set_xticklabels(unit_labels)
            ax.set_yticklabels(unit_labels)
            ax.set_aspect("equal")
        top.set_title(F.STRATUM_LABEL[st], fontweight="bold")
        top.set_xlabel("1 − specificity")
        bottom.set_xlabel("Sensitivity (recall)")
        if j == 0:
            top.set_ylabel("Sensitivity")
            bottom.set_ylabel("Precision")
    F._panel_letter(axes[0, 0], "a", x=-0.2)
    F._panel_letter(axes[1, 0], "b", x=-0.2)
    tools = [Line2D([], [], color=TOOL_COLOUR[t], lw=1.0, ls=TOOL_DASH[t],
                    label=F.NAME[t] + ("\u00a0†" if t == "avi" else ""))
             for t in MAIN_TOOLS]
    tier_marks = [Line2D([], [], marker=TIER_MARKER[t][0], ms=TIER_MARKER[t][1], ls="",
                         color=F.INK2, mec="white", mew=0.4,
                         label=f"Fitted {TIER_NAME[t]} threshold")
                  for t in ("supporting", "moderate", "strong")]
    tier_marks += [Line2D([], [], marker="D", ms=3.8, ls="", mfc="white", mec=F.INK2,
                          mew=0.9, label=f"Published cut point, {pp3_cut:g}")]
    # the legend fills column by column: tools on the first row, markers on the second
    handles = [h for pair in zip(tools, tier_marks) for h in pair]
    handles.append(Line2D([], [], color=F.INK2, lw=0.7, ls=(0, (3, 1.5)),
                          label="Prevalence"))
    fig.legend(handles=handles, loc="lower center", ncol=5, bbox_to_anchor=(0.5, 0.0),
               handletextpad=0.4, columnspacing=1.0, handlelength=2.4)
    fig.subplots_adjust(left=0.07, right=0.99, top=0.95, bottom=0.17)
    F._save(fig, FIGURE)


# ---------------------------------------------------------------------------
# E4.6: the one-page summary, written from the tables above
# ---------------------------------------------------------------------------
# a predictor as a sentence names it, lower case where it is not at the start
IN_SENTENCE = {"spliceai_walker": "SpliceAI", "pangolin": "Pangolin",
               "alphagenome": "the AlphaGenome splice score",
               "avi": "the Atlas combined score"}


def _two(x) -> str:
    """Two decimals; a proportion within rounding of 1 or 0 is not printed as one."""
    if x is None or not np.isfinite(x):
        return "–"
    if 0.995 <= x < 1.0:
        return ">0.99"
    if 0.0 < x < 0.005:
        return "<0.01"
    return f"{x:.2f}"


def _up(x: float, places: int = 2) -> float:
    """Rounded up, for a sentence that states a bound ("never more than")."""
    f = 10 ** places
    return math.ceil(round(x * f, 9)) / f


def _pct(x: float) -> str:
    return f"{100 * x:.0f}%"


def _row(table: pd.DataFrame, tool: str, stratum: str, kind: str, source: str):
    r = table[(table.tool == tool) & (table.stratum == stratum)
              & (table.threshold_kind == kind) & (table.threshold_source == source)]
    return r.iloc[0] if len(r) else None


def _held_out_folds(lg) -> int:
    """The folds a held-out row could have counted: those E3 ran, less any left out."""
    if lg is None or not np.isfinite(lg.folds_total):
        return 0
    excluded = [g for g in str(lg.folds_excluded).split("; ") if g and g != "nan"]
    return int(lg.folds_total) - len(excluded)


def _shown(ins, lg) -> bool:
    """A held-out value is shown when the tier is reached in-sample, the rule E7 applies
    before carrying a threshold anywhere, and when at least half the held-out folds
    reach it, the majority rule of the figure tables; otherwise it rests on one or two
    genes and only the fold count is shown."""
    if ins is None or ins.status != "ok" or lg is None or lg.status != "ok":
        return False
    return 2 * int(lg.folds_used) >= _held_out_folds(lg)


def _direction(v: pd.Series) -> str | None:
    v = v.dropna()
    if v.empty:
        return None
    if (v > 0).all():
        return "rises"
    if (v < 0).all():
        return "falls"
    if (v == 0).all():
        return "is unchanged"
    return "moves both ways"


def _distal_crossings(roc: pd.DataFrame) -> tuple[float, float]:
    """Where some tool's 11-50 bp curve lies above its 3-10 bp curve, both read as step
    functions: the largest such false-positive rate below 0.5 (0 if none) and the
    smallest at or above 0.5 (1 if none). At the start a few top-scoring variants
    decide the curve; at the end the band whose damaging variants run out first
    reaches full sensitivity first."""
    low, high = 0.0, 1.0
    for tool in MAIN_TOOLS:
        a = roc[(roc.tool == tool) & (roc.stratum == "s11_50")]
        b = roc[(roc.tool == tool) & (roc.stratum == "s3_10")]
        grid = np.union1d(a.fpr.to_numpy(), b.fpr.to_numpy())

        def tpr_at(c, x):
            f, t = c.fpr.to_numpy(), c.tpr.to_numpy()
            return np.array([t[f <= v].max() for v in x])
        above = grid[tpr_at(a, grid) > tpr_at(b, grid)]
        if (above < 0.5).any():
            low = max(low, float(above[above < 0.5].max()))
        if (above >= 0.5).any():
            high = min(high, float(above[above >= 0.5].min()))
    return low, high


def _sentences(table: pd.DataFrame, auc: pd.DataFrame, ext: pd.DataFrame,
               roc: pd.DataFrame) -> list[str]:
    """Five sentences for the slides, each carrying at most two numbers. Every claim is
    computed from the tables, and each has a branch for data that would not support
    its wording."""
    main = auc[auc.tool.isin(MAIN_TOOLS)]
    out = []
    # 1. the distal band on the ROC
    best_distal = float(main[main.stratum == "s11_50"].auroc_pooled_curve.max())
    low_proximal = float(main[main.stratum == "s3_10"].auroc_pooled_curve.min())
    low, high = _distal_crossings(roc)
    if best_distal < low_proximal and low < 0.01 and high > 0.9:
        out.append("At 11–50 bp each ROC curve lies below the same tool's 3–10 bp curve "
                   "except near its two ends; the best pooled AUROC there is "
                   f"{best_distal:.2f}, below the lowest at 3–10 bp ({low_proximal:.2f}).")
    else:
        out.append(f"At 11–50 bp the best pooled AUROC is {best_distal:.2f}, against "
                   f"{low_proximal:.2f} for the lowest at 3–10 bp.")
    # 2 and 3. Strong at 3-10 bp, held out, where the table shows a held-out value
    strong = []
    for tool in MAIN_TOOLS:
        ins, lg = (_row(table, tool, "s3_10", "strong", s) for s in ("in_sample", "logo"))
        if _shown(ins, lg):
            strong.append(lg)
    if strong:
        sens = [float(r.sensitivity) for r in strong]
        out.append(f"Held out, a Strong threshold at 3–10 bp catches {_pct(min(sens))} to "
                   f"{_pct(max(sens))} of damaging variants, depending on the tool.")
        cover = max(float(r.coverage) for r in strong)
        out.append("The cost is coverage: Strong evidence reaches no more than "
                   f"{_pct(_up(cover))} of the variants in that band.")
    else:
        out.append("Held out, no tool's Strong threshold at 3–10 bp holds in most genes.")
    # 4. matched specificity
    out.append(_sentence_matched(main))
    # 5. the external genes against the held-out seven
    out.append(_sentence_external(ext))
    return out


def _sentence_matched(main: pd.DataFrame) -> str:
    """Which tool is most sensitive at a specificity of 0.95, and by how much."""
    col = "sensitivity_at_specificity_0.95"
    leaders, margin, spread = [], 0.0, 0.0
    for st in CURVE_STRATA:
        v = main[main.stratum == st].set_index("tool")[col].sort_values(ascending=False)
        leaders.append(set(v.index[v == v.iloc[0]]))
        below = v[v < v.iloc[0]]
        margin = max(margin, float(v.iloc[0] - below.iloc[0]) if len(below) else 0.0)
        spread = max(spread, float(v.iloc[0] - v.iloc[-1]))
    common = set.intersection(*leaders)
    if len(common) == 1:
        (lead,) = common
        tied = any(len(s) > 1 for s in leaders)
        text = (f"At a specificity of 0.95, {IN_SENTENCE[lead]} is the most sensitive"
                f"{', or tied for it,' if tied else ''} in all three bands, though by no "
                f"more than {_up(margin):.2f}.")
    elif len(common) > 1:
        names = " and ".join(IN_SENTENCE[t] for t in MAIN_TOOLS if t in common)
        text = (f"At a specificity of 0.95, {names} are tied for the most sensitive in "
                "all three bands.")
    else:
        text = ("At a specificity of 0.95 no tool is the most sensitive in every band, and "
                f"the four are never more than {_up(spread):.2f} apart.")
    return text[0].upper() + text[1:]


def _sentence_external(ext: pd.DataFrame) -> str:
    """The external genes at the carried thresholds against the held-out seven genes."""
    e = ext[(ext.status == "ok") & ext.external_test.astype(bool) & ext.main_table
            & ext.stratum.isin(CURVE_STRATA) & ext.threshold_kind.isin(EXTERNAL_TIERS)
            & ext.sensitivity_minus_seven_gene.notna()]
    dd, tp = e[e.gene == "DDX3X"], e[e.gene == "TP53"]
    d_sens = _direction(dd.sensitivity_minus_seven_gene)
    d_spec = _direction(dd.specificity_minus_seven_gene)
    if d_sens is None:
        return "No threshold carried to DDX3X can be set against the seven genes."
    words = {"rises": "more", "falls": "less"}
    if d_sens in words and d_spec in words:
        head = (f"Carried to DDX3X, the fitted thresholds are {words[d_sens]} sensitive "
                "than in the held-out seven genes, by as much as "
                f"{dd.sensitivity_minus_seven_gene.abs().max():.2f}, and {words[d_spec]} "
                "specific, by as much as "
                f"{dd.specificity_minus_seven_gene.abs().max():.2f}")
    else:
        head = (f"Carried to DDX3X, sensitivity {d_sens} against the held-out seven genes, "
                f"by as much as {dd.sensitivity_minus_seven_gene.abs().max():.2f}, and "
                f"specificity {d_spec}")
    t_sens = _direction(tp.sensitivity_minus_seven_gene)
    if t_sens is None:
        tail = ""
    elif t_sens == "moves both ways":
        # name the exceptions to the majority direction when one tool holds them all
        up = tp[tp.sensitivity_minus_seven_gene > 0]
        down = tp[tp.sensitivity_minus_seven_gene < 0]
        minority = down if len(down) < len(up) else up
        majority_word = "rises" if len(down) < len(up) else "falls"
        if len(up) != len(down) and minority.tool.nunique() == 1:
            tiers = [TIER_NAME[k] for k in EXTERNAL_TIERS
                     if k in set(minority.threshold_kind)]
            tail = (f"; on TP53 sensitivity {majority_word} too, except for "
                    f"{IN_SENTENCE[minority.tool.iloc[0]]}'s {' and '.join(tiers)}")
        else:
            tail = "; on TP53 sensitivity moves both ways"
    else:
        tail = f"; on TP53 sensitivity {t_sens}" + (" too" if t_sens == d_sens else "")
    return head + tail + "."


def summary(table: pd.DataFrame, auc: pd.DataFrame, ext: pd.DataFrame,
            roc: pd.DataFrame, xmax: float) -> str:
    tiers = ("supporting", "moderate", "strong")
    metrics = ("sensitivity", "specificity", "coverage")
    head = (["Predictor", "Band"]
            + [f"{TIER_NAME[t]}: {m}" for t in tiers for m in metrics]
            + ["Held-out folds reaching S, M, St", "Pooled AUROC", "Pooled AUPRC"])
    body = []
    au = auc.set_index(["tool", "stratum"])
    for tool in MAIN_TOOLS:
        for st in CURVE_STRATA:
            cells = [F.NAME[tool] + (" †" if tool == "avi" else ""), F.STRATUM_LABEL[st]]
            folds, held_out = [], 0
            for t in tiers:
                ins = _row(table, tool, st, t, "in_sample")
                lg = _row(table, tool, st, t, "logo")
                reached = ins is not None and ins.status == "ok"
                held = lg if _shown(ins, lg) else None
                for m in metrics:
                    cells.append(f"{_two(held[m]) if held is not None else '–'} "
                                 f"({_two(ins[m]) if reached else '–'})")
                held_out = max(held_out, _held_out_folds(lg))
                folds.append(str(int(lg.folds_used)) if lg is not None
                             and np.isfinite(lg.folds_used) else "–")
            cells.append(", ".join(folds) + (f" of {held_out}" if held_out else ""))
            cells.append(_two(au.loc[(tool, st), "auroc_pooled_curve"]))
            cells.append(_two(au.loc[(tool, st), "auprc_pooled_curve"]))
            body.append(cells)

    md = ["# Practical metrics at the evidence thresholds", "",
          "Written by `phase1/src/evid_practical_metrics.py` (stage E15 of "
          "`scripts/reproduce_evidence.py`) from the tables it writes to "
          "`phase1/reports/evidence/`; do not edit by hand. The thresholds are the "
          "study's own: nothing here is fitted.", "",
          "## Sensitivity, specificity and coverage at the fitted thresholds", "",
          "Each cell gives the held-out value first and the in-sample value in "
          "brackets, for the seven genes.", "",
          "| " + " | ".join(head) + " |", "|" + "|".join("---" for _ in head) + "|"]
    md += ["| " + " | ".join(r) + " |" for r in body]
    md += ["",
           "- Held out: each fold's threshold, fitted on the other genes (six, or five at "
           "11–50 bp, where BRCA2 has no variant), applied to the gene left out; the "
           "counts are added over the folds whose threshold reached the tier, not "
           "averaged. In brackets: the threshold fitted on all the genes, applied to "
           "them.",
           "- Held-out folds reaching S, M, St: for each tier, the folds whose threshold "
           "reached it, of the genes held out. Those folds include ones whose band holds "
           "too few variants for a likelihood ratio, which Table 3 counts as not "
           "evaluable, so the count is not Table 3's k/n.",
           "- A held-out value is shown when the tier is reached in-sample (the rule E7 "
           "applies before carrying a threshold to another gene) and in at least half "
           "the held-out folds; otherwise it would rest on one or two genes, and only "
           "the fold count is given. A dash in brackets: the tier is not reached "
           "in-sample.",
           "- The Atlas combined score (†): its model was selected on the BRCA1 and RAD51C "
           "assays, so its held-out values leave those two genes out. Its in-sample "
           "values, in brackets, include them.",
           "- Sensitivity: the share of damaging variants at or above the threshold. "
           "Specificity: the share of normal variants below it. Coverage: the share of "
           "all scored variants in the band, labelled or not, at or above it; these are "
           "the variants a user would get the evidence for. The share among labelled "
           "variants only, the denominator of the manuscript's '14.8% of the variants' at "
           "0.2, is `coverage_labelled`.",
           "- Pooled AUROC and AUPRC: the seven genes pooled into one curve, the way one "
           "threshold is used across genes. The manuscript's AUROC (Methods) is "
           "gene-pooled instead: per-gene AUROCs combined on the logit scale with "
           "random effects. For all variants it is `auroc_gene_pooled` in "
           "`practical_curve_auc.csv`, beside the pooled value and a split of the "
           "difference.",
           "- Two decimals; >0.99 marks a value from 0.995 up to, not including, 1.",
           "- Wilson intervals, the BP4 side, the other ten columns and the fusion are "
           "in `practical_metrics_by_threshold.csv`; each fold in "
           "`practical_metrics_logo_by_gene.csv`; DDX3X and TP53 in "
           "`practical_metrics_external.csv`. The curve tables start at threshold "
           "+inf, where nothing is called; the precision-recall table's first row is "
           "the conventional (recall 0, precision 1), not an observation.",
           "", "## In five sentences", ""]
    md += [f"{i}. {s}" for i, s in enumerate(_sentences(table, auc, ext, roc), 1)]
    md += ["", "## Figure legend", "",
           f"`phase1/reports/evidence/figures/{FIGURE}.pdf`. ROC curves (a) and "
           "precision-recall curves (b) of four predictors, the seven genes pooled, by "
           "distance from the exon boundary. Colours are those of Supplementary Figure "
           "S2, and each predictor also has its own line pattern. Filled markers are the "
           "fitted in-sample thresholds: circle Supporting, square Moderate, triangle "
           "Strong; a tier that is not reached is not drawn, and markers of different "
           "predictors are drawn slightly apart so that coinciding ones stay visible. "
           "The open diamond is the published cut point, 0.2, on SpliceAI (published "
           "basis). Insets in (a) enlarge the high-specificity end, 1 − specificity up to "
           f"{xmax:g} against sensitivity, where every fitted threshold lies. The key in "
           "each ROC panel gives the pooled AUROC; the dashed line in (b) is the band's "
           "prevalence. † The Atlas combined score's model was selected on the BRCA1 and "
           "RAD51C assays, which these curves include.", ""]
    return "\n".join(md)


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------
def load() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict, list[str]]:
    cfg = yaml.safe_load(CONFIG_PATH.read_text())
    ev = _read(REPORT_DIR / "evidence_thresholds.csv")
    folds = _read(REPORT_DIR / "evidence_thresholds_logo_folds.csv")
    df = K.load_set()
    have = set(ev["tool"])
    tools = [c for c in F.COLUMNS if c in have] + ([K.FUSION] if K.FUSION in have else [])
    if K.FUSION in have:
        # the out-of-fold fusion E3 fitted its thresholds on, rebuilt the same way; the
        # maths library raises floating-point flags inside sklearn's matmul that do not
        # reach the scores. A variant without a score is left out the way E3 leaves it.
        with np.errstate(all="ignore"):
            df[K.FUSION] = K.logo_fusion(
                df, [t for t in K.FUSION_FEATURES if t in df.columns])
    return df, ev, folds, cfg, tools


def main() -> None:
    df, ev, folds, cfg, tools = load()
    table, by_gene = by_threshold(df, ev, folds, cfg, tools)
    ext = external(table, cfg, tools)
    roc, pr, auc = curves(df, tools)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    for path, d in ((OUT_TABLE, table), (OUT_BY_GENE, by_gene), (OUT_EXTERNAL, ext),
                    (OUT_ROC, roc), (OUT_PR, pr), (OUT_AUC, auc)):
        d.to_csv(path, index=False)
        print(f"[E15] wrote {path} ({len(d):,} rows)")
    figure(roc, pr, auc, table, cfg)
    xmax, _ = _inset_range(roc, _tier_marks(table))
    OUT_MD.write_text(summary(table, auc, ext, roc, xmax), encoding="utf-8")
    print(f"[E15] wrote {OUT_MD}")


if __name__ == "__main__":
    main()
