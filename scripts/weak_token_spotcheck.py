"""Hand spot-check of the number check's WEAK matches: 20 tokens, seeded.

A WEAK match is a manuscript token whose value lies within tolerance of at least 25
pipeline values, so the forward check (phase1/src/check_manuscript_numbers.py) accepts
it without being able to say which value it is. The v2 review found 161 of them on
calibration_draft.docx sha256 893076f1... and asked for a seeded sample of 20 to be
traced to the cell each should come from.

The token list is frozen in docs/weak-token-spotcheck.tokens.json so the sample can
be redrawn after the text moves on. The source of each sampled token is recorded by
hand in SPEC below -- that is the manual part -- but every source VALUE is read from
the outputs or the pipeline's constants when this script runs, and agreement is
judged at the token's printed precision with the forward check's own tolerance.

    python scripts/weak_token_spotcheck.py      # writes docs/weak-token-spotcheck.md
"""
from __future__ import annotations

import json
import sys
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "phase1"))
from src.check_manuscript_numbers import _tolerance, pipeline_values, RESULTS  # noqa: E402
import src.phase3_calibration as P3  # noqa: E402
import src.phase5_likelihood_ratios as P5  # noqa: E402

REP = REPO / "phase1" / "reports" / "phase1"
TOKENS = REPO / "docs" / "weak-token-spotcheck.tokens.json"
OUT = REPO / "docs" / "weak-token-spotcheck.md"
SEED = 20260914
EXPECTED_SAMPLE = [0, 29, 37, 46, 47, 51, 52, 66, 85, 91, 92, 95, 103, 110, 112, 124, 126, 130, 137, 157]


def csv(name):
    return pd.read_csv(REP / name)


def cell(name, where, col):
    d = csv(name)
    for k, v in where.items():
        d = d[d[k] == v]
    assert len(d) == 1, (name, where, len(d))
    return d.iloc[0][col]


def ci_part(s, k):
    import re
    return float(re.findall(r"[-+]?\d*\.?\d+", s)[k])


LR = ("phase5_likelihood_ratios.csv", {"condition": "y_assay/BRCA1_included", "calibration": "calibrated_isotonic",
                                       "object": "single:pangolin", "cal_method": "isotonic"}, "lr_plus")
H2 = ("phase2_H2_ablation.csv", {"model": "M1_enet", "ablation": "drop_conservation"})
TP = ("phase4_tp53_label_definitions.csv", {"label_rule": "median split (top 50%)"})

# (block, token) -> (kind, source described, value in the token's units or None)
SPEC = {
    ("P3", "11"): ("not a quantity", "digits of the repository name in github.com/oneone00-11/variant-fm-benchmark", None),
    ("P26", "2.2"): ("not a quantity", "cross-reference to Methods 2.2", None),
    ("P29", "19.9"): ("measurement", "phase5_likelihood_ratios.csv · Pangolin, functional standard +BRCA1, isotonic-calibrated · lr_plus",
                      lambda: cell(*LR)),
    ("P39", "2.1"): ("not a quantity", "cross-reference to Methods 2.1", None),
    ("P41", "1"): ("not a quantity", "caption label 'Figure 1'", None),
    ("P43", "0"): ("measurement", "coverage_by_region.csv · AlphaMissense · largest of splice_core / splice_region, as % scored",
                   lambda: 100 * max(cell("coverage_by_region.csv", {"predictor": "alphamissense"}, "splice_core"),
                                     cell("coverage_by_region.csv", {"predictor": "alphamissense"}, "splice_region"))),
    ("P45", "+0.001"): ("measurement", "phase2_H2_ablation.csv · elastic net, drop conservation · delta_full_minus_ablated",
                        lambda: cell(*H2, "delta_full_minus_ablated")),
    ("P47", "7"): ("derived", "rank of phyloP by mean |coefficient| in phase6_enet_weight_stability.csv",
                   lambda: 1 + list(csv("phase6_enet_weight_stability.csv").sort_values("mean_abs", ascending=False)
                                    .iloc[:, 0]).index("phylop")),
    ("P59", "19.9"): ("measurement", "same cell as P29 · 19.9", lambda: cell(*LR)),
    ("P63", "2.4"): ("not a quantity", "cross-reference to Methods 2.4", None),
    ("P63", "93"): ("derived", "Strong threshold 18.7 × (1 − TARGET_SPEC): the sensitivity Strong requires, in % ('more than 93%')",
                    lambda: 100 * 18.7 * (1 - P5.TARGET_SPEC)),
    ("P64", "0.10"): ("constant", "phase3_calibration.LO, the lower edge of the 0.90/0.10 confidence band", lambda: P3.LO),
    ("P68", "−0.014"): ("measurement", "phase4_tp53_label_definitions.csv · median split · ci_lo", lambda: cell(*TP, "ci_lo")),
    ("P74", "95"): ("constant", "phase5_likelihood_ratios.TARGET_SPEC, in %", lambda: 100 * P5.TARGET_SPEC),
    ("P75", "95"): ("constant", "phase5_likelihood_ratios.TARGET_SPEC, in %", lambda: 100 * P5.TARGET_SPEC),
    ("P86", "0"): ("not a quantity", "digit inside MaveDB accession 00000673-0-1", None),
    ("P86", "0.10"): ("constant", "phase5_likelihood_ratios.PRIOR, the Tavtigian prior", lambda: P5.PRIOR),
    ("T1r8c4", "0"): ("measurement", "phase2_leaderboard.csv · GPN-MSA · I2", lambda: cell("phase2_leaderboard.csv", {"model": "single:gpn_msa"}, "I2")),
    ("T2r1c5", "−0.001"): ("measurement", "phase2_H2_ablation.csv · elastic net, drop conservation · ci95 lower bound",
                           lambda: ci_part(cell(*H2, "ci95"), 0)),
    ("T4r2c2", "0.005"): ("measurement", "phase4_tp53_label_definitions.csv · median split · dBrier", lambda: cell(*TP, "dBrier")),
}


def agrees(kind, token, value, tol):
    printed = float(token.replace("−", "-").replace("+", ""))
    if kind == "derived" and token == "93":           # stated as a lower bound: "more than 93%"
        return printed <= value < printed + 1
    return abs(value - printed) <= tol


def main():
    frozen = json.loads(TOKENS.read_text())
    weak = frozen["weak"]
    pick = sorted(np.random.default_rng(SEED).choice(len(weak), size=20, replace=False).tolist())
    assert pick == EXPECTED_SAMPLE, pick
    pool = pipeline_values(RESULTS)
    rows, kinds, errors = [], {}, 0
    for n in pick:
        t = weak[n]
        kind, source, fn = SPEC[(t["uid"], t["token"])]
        kinds[kind] = kinds.get(kind, 0) + 1
        tol = _tolerance(t["token"])
        v = abs(t["value"])
        near = pool[np.argmin(np.minimum(np.abs(pool - v), np.abs(pool - v / 100) * 100))]
        if fn is None:
            val, verdict = "—", "n/a (not a quantity)"
        else:
            x = float(fn())
            ok = agrees(kind, t["token"], x, tol)
            errors += 0 if ok else 1
            val, verdict = f"{x:.6g}", "yes" if ok else "**NO**"
        ctx = t["context"].strip()
        ctx = (ctx[:72] + "…") if len(ctx) > 72 else ctx
        rows.append(f"| {n} | {t['uid']} | {t['token']} | …{ctx} | {t['n_pool_matches']} within tolerance; nearest {near:.6g} "
                    f"| {kind}: {source} | {val} | {verdict} |")

    def once(x, nd):
        return str(Decimal(repr(float(x))).quantize(Decimal(1).scaleb(-nd), rounding=ROUND_HALF_UP))

    lb = csv("phase2_leaderboard.csv").set_index("model")
    pool = lb.rename(columns={"pooled_rho": "rho"})
    pr = csv("phase3_calibration_summary.csv")
    pr = pr[(pr.set == "y_assay/BRCA1_included") & (pr.calib == "isotonic")].set_index("model")
    v1 = pd.read_csv(REPO / "phase1" / "reports" / "phase1_v1" / "phase3_calibration_summary.csv")
    v1 = v1[(v1.set == "y_assay/BRCA1_included") & (v1.calib == "isotonic")].set_index("model")
    tp = csv("phase4_tp53_external.csv").set_index("model")
    best_tp = [i for i in tp.index if i.startswith("best_single")][0]
    gate = csv("phase4_tp53_orientation.csv").set_index("orientation")
    gate_run = float(gate.loc[[i for i in gate.index if i.startswith("as run")][0], "control_auroc"])
    # the values as the tables stored them (four decimals; one for I^2) before they went unrounded
    stored = {"I2": 59.5, "pangolin_ece": 0.0325, "spliceai_frac": 0.7045, "nt_ece": 0.0475,
              "nt_frac": 0.0866, "tp53_share": 0.7865}
    # (where, cell, printed before the fix, unrounded value, decimals, cause)
    found = [
        ("Table 1", "Fusion (M1) · I²", "60", pool.loc["M1_enet", "I2"], 0,
         f"second rounding: the leaderboard stored I² to one decimal ({stored['I2']}) and the v1→v2 pass kept 60"),
        ("Table 3", "Pangolin · ECE", "0.033", pr.loc["single:pangolin", "ECE"], 3,
         f"second rounding of the stored {stored['pangolin_ece']}, which did not move between v1 and v2, so the pass left the cell"),
        ("Figure 2 caption", "Pangolin · ECE", "0.033", pr.loc["single:pangolin", "ECE"], 3, "the same value"),
        ("Table 3", "SpliceAI · high-confidence %", "70.5", 100 * pr.loc["single:spliceai", "actionable_frac"], 1,
         f"second rounding: the pass wrote 70.5 in from the stored {stored['spliceai_frac']}"),
        ("Table 3", "Nucleotide Transformer · ECE", "0.048", pr.loc["single:nt", "ECE"], 3,
         f"second rounding: the pass wrote 0.048 in from the stored {stored['nt_ece']}"),
        ("Table 3", "Nucleotide Transformer · high-confidence %", "8.8", 100 * pr.loc["single:nt", "actionable_frac"], 1,
         f"stale v1 value (v1 stored {v1.loc['single:nt', 'actionable_frac']}): the pass classed 8.8 as a section "
         "number and never replaced it"),
        ("3.4", "Pangolin on TP53 · high-confidence %", "78.7", 100 * tp.loc[best_tp, "actionable_frac"], 1,
         f"second rounding: stored {stored['tp53_share']}; found when phase 4 went unrounded (P8)"),
        ("3.4", "TP53 post-orientation control AUROC", "0.998", gate_run, 3,
         "never persisted: the gate was printed once in July and quoted from a comment; the number check matched 0.998 "
         "to BRCA2's row of directionality_check.csv. The supplement said 0.997. Now written to phase4_tp53_orientation.csv"),
    ]
    out_rows = []
    for where, what, printed, unrounded, nd, cause in found:
        correct = once(unrounded, nd)
        assert correct != printed, f"{where} {what}: no longer a misprint"
        out_rows.append(f"| {where} | {what} | {printed} | **{correct}** | {float(unrounded):.6g} | {cause} |")

    # every token the v1 -> v2 pass exempted through the section-number pattern
    import re
    rep = json.loads((REPO / "docs" / "manuscript-replacements.json").read_text())
    sect = re.compile(r"^(?:[1-9]|1[0-9])\.[0-9]$")
    exempted = [t for t in rep if t["status"] == "whitelisted" and sect.match(t["token"].lstrip("+-−"))]
    sel = csv("phase8_selection_test_stratified.csv").set_index(["stratum", "subset", "object"]).lr_plus
    ratios = [("fusion, BRCA1 excluded, recorded", "20.2", sel[("minus_BRCA1", "ClinVar-recorded", "fusion_M1")]),
              ("fusion, BRCA1 excluded, assay-only", "13.7", sel[("minus_BRCA1", "assay-only", "fusion_M1")]),
              ("CADD, BRCA1 excluded, recorded", "19.1", sel[("minus_BRCA1", "ClinVar-recorded", "cadd")]),
              ("CADD, BRCA1 excluded, assay-only", "9.9", sel[("minus_BRCA1", "assay-only", "cadd")]),
              ("fusion, assay-only VUS", "12.1", sel[("assay_only_VUS", "assay-only", "fusion_M1")])]
    ratio_text = "; ".join(f"{lab} {printed} (unrounded {v:.4f}, {'agrees' if once(v, 1) == printed else 'DISAGREES'})"
                           for lab, printed, v in ratios)
    assert all(once(v, 1) == printed for _, printed, v in ratios)

    lines = [
        "# WEAK-token spot-check",
        "",
        f"Manuscript: `calibration_draft.docx` sha256 `{frozen['manuscript_sha256_prefix']}…` (the version the v2 review examined). "
        f"Forward-check counts at that version: matched {frozen['counts']['matched']}, WEAK {frozen['counts']['weak']}, "
        f"whitelisted {frozen['counts']['whitelisted']}, no source {frozen['counts']['no_source']}, "
        f"scoped mismatch {frozen['counts']['scoped_mismatch']}.",
        "",
        f"Sample: `numpy.random.default_rng({SEED}).choice(161, 20, replace=False)`, sorted → {pick}. "
        "Token list frozen in `docs/weak-token-spotcheck.tokens.json`; this page is written by `scripts/weak_token_spotcheck.py`.",
        "",
        "A WEAK token is one the forward check accepts because at least 25 pipeline values lie within its tolerance; the check "
        "picks none of them, so the \"gate\" column gives that count and the nearest pool value rather than a matched cell. "
        "The correct source is recorded by hand; its value is read from the outputs or the pipeline's constants when the "
        "script runs, and agreement uses the forward check's own tolerance (half a unit in the last printed digit).",
        "",
        "| # | block | token | context | gate | correct source | source value | agrees |",
        "|---|---|---|---|---|---|---|---|",
        *rows,
        "",
        f"**Result.** {20 - errors} of 20 agree; {errors} disagree. By kind: "
        + ", ".join(f"{k} {v}" for k, v in sorted(kinds.items())) + ". "
        "No sampled token needed a fix, so no claims anchor was added for the sample.",
        "",
        "## Found outside the sample: eight misprints, three causes",
        "",
        "Redrawing Figure 1 from the tables, and tracing the cells behind it, found misprints the sample could not "
        "reach. The forward check had accepted every one, because each lies within tolerance of some pipeline value. "
        "Every correct value below is recomputed from an unrounded output when this page is written.",
        "",
        "| where | cell | printed | correct | unrounded | cause |",
        "|---|---|---|---|---|---|",
        *out_rows,
        "",
        "**Second rounding.** These tables stored four decimals (I² one) and the manuscript prints three (one). Rounding "
        "the stored value again misprints the cell whenever it lands on a tie. The tables the manuscript prints from "
        "now store their point estimates unrounded (`phase2_leaderboard.csv`, `phase3_calibration_summary.csv`, "
        "`phase4_tp53_external.csv` and the others listed in `tests/test_reported_numbers.py::FULL_PRECISION`); the "
        "supplement prints them at a precision a reader can round back to the main text "
        "(`build_supp_tables.col_digits`), and `tests/test_cross_document_rounding.py` checks that boundary. Table "
        "cells cannot carry claims anchors, which bind a unique substring of a paragraph, so the guards are tests: "
        "`test_table_1_rounds_the_unrounded_pooled_values_once` and `test_table_3_rounds_the_calibration_summary_once` "
        "in `tests/test_reported_numbers.py` compare every cell of the two tables with the unrounded values rounded once.",
        "",
        "**A number with no persisted source.** The TP53 orientation gate on the full SGE was computed once in July, "
        "printed, and quoted from a comment; the number check accepted 0.998 because BRCA2's control AUROC (0.9975) "
        "lies within tolerance of it. `phase4_external_tp53.py` now writes the gate to `phase4_tp53_orientation.csv` "
        "on every run (0.9969 as run, 0.0031 under the default orientation), and 3.4 and Supplementary Table S6's "
        "note both quote 0.997 from it.",
        "",
        f"**A measurement exempted as a section number.** The whitelist exempts section numbers by pattern (a digit or "
        f"two, a point, a digit), and a pattern cannot tell a cross-reference from a measurement. The v1→v2 pass "
        f"exempted {len(exempted)} tokens this way. Besides the Table 3 cell, the measurements among them are the "
        "Tavtigian threshold 18.7 (a literature constant), numbers in two sentences rewritten in full during the v2 "
        "revision, and five likelihood ratios in the selection-test paragraph, rechecked here against "
        f"`phase8_selection_test_stratified.csv`: {ratio_text}. The exemption now holds only for numbers that open a "
        "heading, in the forward check and in the replacement tool, and "
        "`tests/test_analysis_claims.py::test_the_section_number_exemption_covers_headings_not_measurements` pins it.",
        "",
        "**Checked and clean.** Every other cell of Tables 1 and 3; all of Table 2, including its two ties (the elastic "
        "net's drop-all-evolution lower bound, stored 0.0005, is 0.0005103 in a replay of phase 2's bootstrap, so "
        "+0.001 stands; the gradient-boosted drop-conservation ablated ρ, stored 0.7595, recomputes to 0.7594899, so "
        "0.759 stands); Table 4, whose source stores unrounded values; the pooled ρ quoted in 3.1 against "
        "`phase2_pooled_rho.csv`; and the other numbers of the Figure 2 caption.",
    ]
    OUT.write_text("\n".join(lines) + "\n")
    print(f"{20 - errors}/20 agree ({errors} disagree); kinds {kinds} -> {OUT.relative_to(REPO)}")


if __name__ == "__main__":
    main()
