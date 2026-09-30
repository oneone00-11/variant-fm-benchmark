"""A threshold read back from a CSV is the threshold that was written, and a printed
threshold selects the variants the fitted one selects.

The defect this exists for: evid_arms read E3's thresholds with pandas' default CSV
float parser, which is not correctly rounded. It returned the combined Atlas score's
Moderate threshold at 3-10 bp, written as 1.0710335969924927, one binary digit high,
and the one damaging PALB2 variant scoring exactly that value fell out of the band
"score >= threshold": Supplementary Table S10 printed 340 variants in band,
sensitivity 0.516 and a band ratio of 28.1 where the fitted threshold selects 341,
0.518 and 28.2. The same misreading made the supplementary tables choose the digits of
several printed thresholds against the wrong value. Every stage now reads thresholds
through evid_common.read_back.
"""
from __future__ import annotations

import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

PHASE1 = Path(__file__).resolve().parents[1]
if str(PHASE1) not in sys.path:
    sys.path.insert(0, str(PHASE1))
EV = PHASE1 / "reports/evidence"
SUPP = EV / "supplement"
THRESHOLDS = EV / "evidence_thresholds.csv"
FOLDS = EV / "evidence_thresholds_logo_folds.csv"
ARMS = EV / "arms_at_tool_threshold.csv"
SET = PHASE1 / "data/evidence/analysis_set_v1.parquet"
# the damaging PALB2 variant whose combined Atlas score is exactly E3's Moderate
# threshold for that score at 3-10 bp
PALB2 = "urn:mavedb:00001259-a-2#9971"

need = pytest.mark.skipif(not (THRESHOLDS.exists() and ARMS.exists()),
                          reason="E3 or E2.3 not run")
need_supp = pytest.mark.skipif(not ((SUPP / "tableS6_thresholds_pathogenic.csv").exists()
                                    and SET.exists()),
                               reason="supplementary tables not built, or no analysis set")

# Every default-parser read left in the evidence stages, and why none carries a
# threshold. A new read, or a change to one of these, fails the test below until it
# is either routed through read_back or added here with its reason.
DEFAULT_READS = {
    "evid_build_set.py: lab = pd.read_csv(path, sep=None, engine=\"python\")":
        "the deposits' label files, before any threshold exists",
    "evid_delta.py: return pd.read_csv(p) if p.exists() else pd.DataFrame()":
        "the published study's tables (read twice), compared to four decimals",
    "evid_diagnostics.py: c = pd.read_csv(path).dropna(subset=[\"local_lr\"])":
        "local likelihood-ratio curves, for the monotonicity counts",
    "evid_external.py: cls = pd.read_csv(cpath, usecols=[\"hgvs_nt\", \"score\", DDX3X_CLASS_COLUMN])":
        "the DDX3X deposit's classification",
    "evid_fig_data.py: w = pd.read_csv(REPORT_DIR / \"walker_thresholds.csv\")":
        "ratios at Walker's cut points, which come from the configuration",
    "evid_fig_data.py: a = pd.read_csv(arms)":
        "ratios without BRCA1 at Walker's cut point",
    "evid_fig_data.py: c = pd.read_csv(f)":
        "local likelihood-ratio curves for a draft figure",
    "evid_fig_data.py: tm = pd.read_csv(REPORT_DIR / \"territory_metrics.csv\")":
        "AUROC columns only; its threshold column is not used",
    "evid_figures.py: d = pd.read_csv(FIG_DATA / \"fig1_walker_cutpoints.csv\")":
        "ratios at Walker's cut points",
    "evid_figures.py: c = pd.read_csv(p).sort_values(\"score\")":
        "local likelihood-ratio curves drawn in Figure 2",
    "evid_figures.py: basis = {g: pd.read_csv(REPORT_DIR / f\"external_{g.lower()}_column_basis.csv\")":
        "the external genes' column-basis flags",
    "evid_figures.py: a = pd.read_csv(REPORT_DIR / \"territory_metrics_arms_noBRCA1.csv\")":
        "AUROCs by ClinVar arm",
    "evid_figures.py: d = pd.read_csv(p)":
        "shares of dilution subsets for Supplementary Figure S2",
    "evid_supp_tables.py: return pd.read_csv(_need(p), **kw)":
        "the helper behind the _csv reads listed here",
    "evid_supp_tables.py: comp = (pd.read_csv(ATLAS_COMPOSITION, sep=\"\\t\").set_index(\"gene\")":
        "the atlas's composition table",
    "evid_supp_tables.py: sc = _csv(REPORTS / \"set_counts.csv\")":
        "counts (read twice)",
    "evid_supp_tables.py: dep = _csv(_need(TP53_DEPOSIT), sep=\"\\t\", usecols=[\"accession\", \"HGVS(cDNA)\"])":
        "the TP53 deposit's identifiers",
    "evid_supp_tables.py: train = _csv(REPORTS / \"predictor_training_provenance.csv\").set_index(\"score_column\")":
        "training provenance, text",
    "evid_supp_tables.py: orient = _csv(REPORTS / \"feature_orientation.csv\").set_index(\"feature\")":
        "orientation correlations",
    "evid_supp_tables.py: cc = _csv(REPORTS / \"column_concordance.csv\")":
        "column concordance",
    "evid_supp_tables.py: w = _csv(REPORTS / \"walker_thresholds.csv\")":
        "ratios at Walker's cut points, which come from the configuration",
    "evid_supp_tables.py: t = _csv(REPORTS / \"territory_metrics_arms_noBRCA1.csv\",":
        "AUROCs and ratios at Walker's cut point by arm",
    "evid_supp_tables.py: w = _csv(REPORTS / \"arms_within_gene.csv\")":
        "within-gene arm comparisons",
    "evid_supp_tables.py: basis = _csv(bp).set_index(\"tool\")":
        "the external genes' column-basis flags",
    "evid_supp_tables.py: sm = _csv(path)":
        "in-frame attribution summaries",
    "evid_supp_tables.py: \"Rows\": str(len(pd.read_csv(path))),":
        "a row count",
    "evid_tables.py: sc = pd.read_csv(REPORT_DIR / \"set_counts.csv\")":
        "counts",
    "evid_tables.py: t = pd.read_csv(REPORT_DIR / \"predictor_training_provenance.csv\")":
        "training provenance, text",
    "evid_tables.py: am = pd.read_csv(REPORT_DIR / \"territory_metrics_arms_noBRCA1.csv\")":
        "AUROCs by ClinVar arm",
    "evid_tables.py: print(f\"[tables] wrote {f.name} ({len(pd.read_csv(f))} rows)\")":
        "a row count",
    "evid_tp53_extend.py: return pd.read_csv(buf, sep=\"\\t\")":
        "the TP53 deposit, checked against its manifest",
}
TWICE = {"evid_delta.py: return pd.read_csv(p) if p.exists() else pd.DataFrame()",
         "evid_supp_tables.py: sc = _csv(REPORTS / \"set_counts.csv\")"}


def test_every_default_read_is_one_that_carries_no_threshold():
    """Each read in the evidence stages that does not go through read_back is on the
    list above, with the reason it needs no exact read; paths held in variables are
    caught too, because the check is the read itself, not the file name in it."""
    found = Counter()
    for p in sorted((PHASE1 / "src").glob("evid_*.py")):
        for line in p.read_text().splitlines():
            code = line.strip()
            if (re.search(r"\bread_csv\(|\b_csv\(", code) and "read_back" not in code
                    and "to_csv" not in code and not code.startswith("def ")
                    and "float_precision=\"round_trip\"" not in code):
                found[f"{p.name}: {code}"] += 1
    expected = Counter({k: (2 if k in TWICE else 1) for k in DEFAULT_READS})
    assert found == expected, (
        "default-parser reads not on the list: "
        + "; ".join(sorted(found - expected)) + " | listed but gone: "
        + "; ".join(sorted(expected - found)))


@need
def test_every_stored_threshold_reads_back_as_written():
    """Each threshold cell, and each copy of one, read through read_back is float() of
    its own text."""
    from src import evid_common as K
    columns = {
        THRESHOLDS: ["pp3_threshold_insample", "bp4_threshold_insample",
                     "pp3_threshold_variant_boot", "bp4_threshold_variant_boot",
                     "pp3_logo_threshold_median"],
        FOLDS: ["threshold_from_6_genes"],
        ARMS: ["threshold"],
        EV / "tier_logo_table.csv": ["in_sample_threshold"],
        EV / "territory_metrics.csv": ["e3_threshold_at_best_tier"],
        EV / "threshold_dilution_summary.csv": ["threshold_median", "threshold_all_genes",
                                               "threshold_q25", "threshold_q75"],
        EV / "external_ddx3x.csv": [f"e3_{t}_threshold" for t in ("supporting", "moderate", "strong")],
        EV / "external_tp53.csv": [f"e3_{t}_threshold" for t in ("supporting", "moderate", "strong")],
    }
    checked = 0
    for path, cols in columns.items():
        text = pd.read_csv(path, dtype=str, keep_default_na=False)
        exact = K.read_back(path)
        for col in cols:
            for t, v in zip(text[col], exact[col]):
                if t:
                    assert v == float(t), (path.name, col, t, v)
                    checked += 1
    assert checked > 1000


@need
def test_the_palb2_variant_on_the_threshold_is_in_its_band():
    """The regression itself: the threshold read back equals the value written, the
    variant scoring exactly that value is in the band, and so the arm counts are the
    fitted threshold's: 341 in band and 295 of 570 damaging over all arms, 120 in the
    recorded-unclassified arm."""
    from src import evid_common as K
    text = pd.read_csv(THRESHOLDS, dtype=str)
    t = text[(text.tool == "avi") & (text.stratum == "s3_10") & (text.tier == "moderate")]
    written = t.pp3_threshold_insample.iloc[0]
    assert written == "1.0710335969924927"
    ev = K.read_back(THRESHOLDS)
    tau = ev[(ev.tool == "avi") & (ev.stratum == "s3_10")
             & (ev.tier == "moderate")].pp3_threshold_insample.iloc[0]
    assert tau == float(written)

    if SET.exists():
        df = pd.read_parquet(SET)
        v = df[df.variant_id == PALB2].iloc[0]
        assert (v.gene, v.stratum, v.y_assay, v.clinvar_arm) == (
            "PALB2", "s3_10", 1.0, "recorded_unclassified")
        assert v.avi == tau                        # on the threshold, so in the band

    arms = K.read_back(ARMS)
    a = arms[(arms.tool == "avi") & (arms.stratum == "s3_10")
             & (arms.threshold_basis == "fitted moderate") & (arms.gene_set == "all_genes")]
    assert (a.threshold == tau).all()
    whole = a[a.clinvar_arm == "all"].iloc[0]
    assert (whole.n_in_band, whole.n_pos) == (341, 570)
    assert whole.sens == 295 / 570
    assert a[a.clinvar_arm == "recorded_unclassified"].n_in_band.iloc[0] == 120

    s10 = pd.read_csv(SUPP / "tableS10_clinvar_status.csv", dtype=str)
    row = s10[(s10.iloc[:, 0] == "Seven genes") & (s10.iloc[:, 2] == "All variants")
              & (s10.iloc[:, 3] == "Atlas combined score") & (s10.iloc[:, 1] == "3–10")]
    assert list(row.iloc[0, [5, 9, 10, 12]]) == ["1.071", "341", "0.518", "28.2"]


@need_supp
def test_every_printed_threshold_selects_the_fitted_band():
    """S6, S7, S8, S10 and S12: each printed threshold is fmt_thr of the exact value it
    prints, and selects exactly the variants that value selects (at or above it on the
    PP3 side, at or below it on the BP4 side)."""
    from src import evid_common as K
    from src import evid_supp_tables as S
    df = pd.read_parquet(SET)
    tools = {v: k for k, v in S.tool_labels().items()}
    strata = {v: k for k, v in S.STRATUM_LABEL.items()}
    bad, checked = [], 0

    def check(where, printed, exact, tool, side):
        nonlocal checked
        if printed == "" and not np.isfinite(exact):
            return
        s = df[tool].dropna().to_numpy(dtype=float)
        count = (lambda v: int((s >= v).sum())) if side == "upper" else (lambda v: int((s <= v).sum()))
        checked += 1
        if printed != S.fmt_thr(exact, tool, side) or count(float(printed)) != count(exact):
            bad.append((where, printed, repr(exact)))

    def text(name):
        return pd.read_csv(SUPP / name, dtype=str, keep_default_na=False)

    ev = K.read_back(THRESHOLDS).set_index(["tool", "stratum", "tier"])
    for name, side, cols in (
            ("tableS6_thresholds_pathogenic.csv", "upper",
             [("In-sample threshold (gene-clustered bootstrap)", "pp3_threshold_insample"),
              ("In-sample threshold (variant-level bootstrap)", "pp3_threshold_variant_boot")]),
            ("tableS7_thresholds_benign.csv", "lower",
             [("In-sample threshold (gene-clustered bootstrap)", "bp4_threshold_insample"),
              ("In-sample threshold (variant-level bootstrap)", "bp4_threshold_variant_boot")])):
        for _, r in text(name).iterrows():
            key = (tools[r["Tool"]], strata[r["Stratum (intronic offset, bp)"]], r["Tier"].lower())
            for c, col in cols:
                if c in r:
                    check((name, key, c), r[c], float(ev.loc[key, col]), key[0], side)
    dil = K.read_back(EV / "threshold_dilution_summary.csv")
    for _, r in text("tableS8_dilution.csv").iterrows():
        tool, st = tools[r["Tool"]], strata[r["Stratum (intronic offset, bp)"]]
        d = dil[(dil.tool == tool) & (dil.stratum == st) & (dil.tier == r["Tier"].lower())
                & (dil.k == int(r["Training genes"]))].iloc[0]
        check(("S8", tool, st, r["Tier"], "median"), r["Threshold, median"],
              float(d.threshold_median), tool, "upper")
        check(("S8", tool, st, r["Tier"], "all"), r["Threshold on all genes in the stratum"],
              float(d.threshold_all_genes), tool, "upper")
        if r["Threshold, interquartile range"]:
            lo, hi = r["Threshold, interquartile range"].split(" to ")
            check(("S8", tool, st, "q25"), lo, float(d.threshold_q25), tool, "upper")
            check(("S8", tool, st, "q75"), hi, float(d.threshold_q75), tool, "upper")
    arms = K.read_back(ARMS)
    fitted = arms[arms.threshold_basis.str.startswith("fitted ")]
    genes = {"Seven genes": "all_genes", "Without BRCA1": "no_BRCA1"}
    arm = {"All variants": "all", "Classified (P/LP or B/LB)": "classified",
           "Recorded, unclassified": "recorded_unclassified", "Unrecorded": "unrecorded"}
    for _, r in text("tableS10_clinvar_status.csv").iterrows():
        if r["Fitted threshold"] == "":
            continue
        tool, st = tools[r["Tool"]], strata[r["Stratum (intronic offset, bp)"]]
        f = fitted[(fitted.tool == tool) & (fitted.stratum == st)
                   & (fitted.gene_set == genes[r["Genes"]])
                   & (fitted.clinvar_arm == arm[r["ClinVar record status"]])].iloc[0]
        check(("S10", r["Genes"], tool, st, r["ClinVar record status"]),
              r["Fitted threshold"], float(f.threshold), tool, "upper")
    ext = {g: K.read_back(EV / f"external_{g.lower()}.csv") for g in ("DDX3X", "TP53")}
    for _, r in text("tableS12_external_genes.csv").iterrows():
        tool = tools[r["Tool"]]
        for tier in ("Supporting", "Moderate", "Strong"):
            printed = r[f"{tier}: threshold (median of leave-one-gene-out fits)"]
            if not printed:
                continue
            carried = ext[r["Gene"]]
            carried = carried[carried.tool == tool][f"e3_{tier.lower()}_threshold"].dropna().unique()
            match = [c for c in carried if S.fmt_thr(float(c), tool, "upper") == printed]
            checked += 1
            if len(match) != 1:
                bad.append(("S12", r["Gene"], tool, tier, printed, list(map(repr, carried))))
    assert checked > 1000
    assert not bad, bad[:10]


def test_printed_ratios_round_exact_ties_half_up():
    """A ratio of counts on a decimal tie is computed one binary digit off it; the
    printed value is the tie's half-up rounding, not the stored double's."""
    from src import evid_supp_tables as S
    assert (22 / 88) / (5 / 291) == 14.549999999999999          # 291/20, one digit low
    assert S.fmt_lr((22 / 88) / (5 / 291)) == "14.6"
    assert S.fmt_lr(14.549) == "14.5" and S.fmt_lr(14.551) == "14.6"
    assert S.fmt_lr(243 / 4) == "60.8" and S.fmt_lr(float("nan")) == ""


def test_read_back_round_trips_what_to_csv_writes(tmp_path):
    """On values the default parser misreads, including the threshold of the defect
    and a small one it truncates by 1,787 units in the last place."""
    from src import evid_common as K
    rng = np.random.default_rng(3)
    values = np.concatenate([[1.0710335969924927, 0.37547463178634644,
                              0.026346683502197266, -0.09308166056871414,
                              0.0003352165222167969],
                             rng.random(2000) * 3 - 0.5])
    pd.DataFrame({"x": values}).to_csv(tmp_path / "t.csv", index=False)
    assert np.array_equal(K.read_back(tmp_path / "t.csv").x.to_numpy(), values)
