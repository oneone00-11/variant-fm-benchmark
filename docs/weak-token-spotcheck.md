# WEAK-token spot-check

Manuscript: `calibration_draft.docx` sha256 `893076f1ca2bccff…` (the version the v2 review examined). Forward-check counts at that version: matched 378, WEAK 161, whitelisted 4, no source 0, scoped mismatch 0.

Sample: `numpy.random.default_rng(20260914).choice(161, 20, replace=False)`, sorted → [0, 29, 37, 46, 47, 51, 52, 66, 85, 91, 92, 95, 103, 110, 112, 124, 126, 130, 137, 157]. Token list frozen in `docs/weak-token-spotcheck.tokens.json`; this page is written by `scripts/weak_token_spotcheck.py`.

A WEAK token is one the forward check accepts because at least 25 pipeline values lie within its tolerance; the check picks none of them, so the "gate" column gives that count and the nearest pool value rather than a matched cell. The correct source is recorded by hand; its value is read from the outputs or the pipeline's constants when the script runs, and agreement uses the forward check's own tolerance (half a unit in the last printed digit).

| # | block | token | context | gate | correct source | source value | agrees |
|---|---|---|---|---|---|---|---|
| 0 | P3 | 11 | …Repository: github.com/oneone00-11/variant-fm-benchmark · Frozen anal | 96 within tolerance; nearest 11 | not a quantity: digits of the repository name in github.com/oneone00-11/variant-fm-benchmark | — | n/a (not a quantity) |
| 29 | P26 | 2.2 | …hich is why everything is rank-normalised within gene (2.2, 2.3). They a… | 33 within tolerance; nearest 0.0220372 | not a quantity: cross-reference to Methods 2.2 | — | n/a (not a quantity) |
| 37 | P29 | 19.9 | …standard (AlphaGenome 24.3 at 97.2% realised, Pangolin 19.9 at 96.0%) pa… | 38 within tolerance; nearest 19.8993 | measurement: phase5_likelihood_ratios.csv · Pangolin, functional standard +BRCA1, isotonic-calibrated · lr_plus | 19.8861 | yes |
| 46 | P39 | 2.1 | …six audited UTR-intron and seven exon-side variants of 2.1), eleven of t… | 33 within tolerance; nearest 2.0966 | not a quantity: cross-reference to Methods 2.1 | — | n/a (not a quantity) |
| 47 | P41 | 1 | …Figure 1. Ranking among top splice-aware to | 2117 within tolerance; nearest 0.01 | not a quantity: caption label 'Figure 1' | — | n/a (not a quantity) |
| 51 | P43 | 0 | …ed trees, M0b equal-weight mean. AlphaMissense scores ≈0 splice variants… | 5922 within tolerance; nearest 0 | measurement: coverage_by_region.csv · AlphaMissense · largest of splice_core / splice_region, as % scored | 0 | yes |
| 52 | P45 | +0.001 | …barely moved ranking: dropping conservation gave Δρ = +0.001 (95% CI [−0… | 40 within tolerance; nearest 0.0009758 | measurement: phase2_H2_ablation.csv · elastic net, drop conservation · delta_full_minus_ablated | 0.0007 | yes |
| 66 | P47 | 7 | …109 (rank 4 of 20), phastCons 0.050 (5), phyloP 0.018 (7). | 140 within tolerance; nearest 7 | derived: rank of phyloP by mean |coefficient| in phase6_enet_weight_stability.csv | 7 | yes |
| 85 | P59 | 19.9 | …thresholds two point estimates cross it — Pangolin at 19.9 (achieved spe… | 38 within tolerance; nearest 19.8993 | measurement: same cell as P29 · 19.9 | 19.8861 | yes |
| 91 | P63 | 2.4 | …ct the sampling-frame check finds on the ranking axis (2.4), reappearing… | 34 within tolerance; nearest 0.024 | not a quantity: cross-reference to Methods 2.4 | — | n/a (not a quantity) |
| 92 | P63 | 93 | …Strong at this operating point means calling more than 93% of true posit… | 63 within tolerance; nearest 0.93 | derived: Strong threshold 18.7 × (1 − TARGET_SPEC): the sensitivity Strong requires, in % ('more than 93%') | 93.5 | yes |
| 95 | P64 | 0.10 | …the functional standard, against 74.9% under the 0.90/0.10 confidence ba… | 75 within tolerance; nearest 0.1 | constant: phase3_calibration.LO, the lower edge of the 0.90/0.10 confidence band | 0.1 | yes |
| 103 | P68 | −0.014 | …mall and its interval spans zero (ΔBrier 0.005, 95% CI −0.014 to 0.024) … | 35 within tolerance; nearest 0.0140392 | measurement: phase4_tp53_label_definitions.csv · median split · ci_lo | -0.0143973 | yes |
| 110 | P74 | 95 | …oth sides, with no ClinVar label involved anywhere. At 95% specificity t… | 100 within tolerance; nearest 0.95 | constant: phase5_likelihood_ratios.TARGET_SPEC, in % | 95 | yes |
| 112 | P75 | 95 | …y. The operating point also bounds the conclusion — at 95% specificity L… | 100 within tolerance; nearest 0.95 | constant: phase5_likelihood_ratios.TARGET_SPEC, in % | 95 | yes |
| 124 | P86 | 0 | …, 00001250-a-2 (BARD1), 00001259-a-2 (PALB2), 00000673-0-1 (RAD51C), 000… | 5922 within tolerance; nearest 0 | not a quantity: digit inside MaveDB accession 00000673-0-1 | — | n/a (not a quantity) |
| 126 | P86 | 0.10 | …work (95% specificity operating point, Tavtigian prior 0.10) was coordin… | 75 within tolerance; nearest 0.1 | constant: phase5_likelihood_ratios.PRIOR, the Tavtigian prior | 0.1 | yes |
| 130 | T1r8c4 | 0 | …0 | 5922 within tolerance; nearest 0 | measurement: phase2_pooled_rho.csv · GPN-MSA · I2 | 0 | yes |
| 137 | T2r1c5 | −0.001 | …−0.001 to +0.002 | 40 within tolerance; nearest 0.0009758 | measurement: phase2_H2_ablation.csv · elastic net, drop conservation · ci95 lower bound | -0.001 | yes |
| 157 | T4r2c2 | 0.005 | …0.005 | 26 within tolerance; nearest 0.005 | measurement: phase4_tp53_label_definitions.csv · median split · dBrier | 0.00511928 | yes |

**Result.** 20 of 20 agree; 0 disagree. By kind: constant 4, derived 2, measurement 8, not a quantity 6. No sampled token needed a fix, so no claims anchor was added for the sample.

## Found outside the sample: six misprinted cells, two causes

Redrawing Figure 1 from the tables, and tracing the cells behind it, found misprints the sample could not reach. The forward check had accepted every one, because each lies within tolerance of some pipeline value. Every correct value below is recomputed from an unrounded output when this page is written.

| where | cell | printed | correct | unrounded | cause |
|---|---|---|---|---|---|
| Table 1 | Fusion (M1) · I² | 60 | **59** | 59.4887 | second rounding: the v1→v2 pass read the one-decimal 59.5 and kept 60 |
| Table 3 | Pangolin · ECE | 0.033 | **0.032** | 0.0324953 | second rounding of the stored 0.0325, which did not move between v1 and v2, so the pass left the cell |
| Figure 2 caption | Pangolin · ECE | 0.033 | **0.032** | 0.0324953 | the same value |
| Table 3 | SpliceAI · high-confidence % | 70.5 | **70.4** | 70.4478 | second rounding: the pass wrote 70.5 in from the stored 0.7045 |
| Table 3 | Nucleotide Transformer · ECE | 0.048 | **0.047** | 0.0474878 | second rounding: the pass wrote 0.048 in from the stored 0.0475 |
| Table 3 | Nucleotide Transformer · high-confidence % | 8.8 | **8.7** | 8.65672 | stale v1 value (v1 stored 0.0884): the pass classed 8.8 as a section number and never replaced it |

**Second rounding.** These tables store four decimals (I² one) and the manuscript prints three (one). Rounding the stored value again misprints the cell whenever it lands on a tie. Phase 2 now also writes `phase2_pooled_rho.csv` and phase 3 `phase3_calibration_summary_precise.csv`, the same estimates unrounded, beside the published tables, which are unchanged. Table cells cannot carry claims anchors, which bind a unique substring of a paragraph, so the guards are tests: `test_table_1_rounds_the_unrounded_pooled_values_once` and `test_table_3_rounds_the_precise_calibration_values_once` in `tests/test_reported_numbers.py` compare every cell of the two tables with the unrounded values rounded once.

**A measurement exempted as a section number.** The whitelist exempts section numbers by pattern (a digit or two, a point, a digit), and a pattern cannot tell a cross-reference from a measurement. The v1→v2 pass exempted 44 tokens this way. Besides the Table 3 cell, the measurements among them are the Tavtigian threshold 18.7 (a literature constant), numbers in two sentences rewritten in full during the v2 revision, and five likelihood ratios in the selection-test paragraph, rechecked here against `phase8_selection_test_stratified.csv`: fusion, BRCA1 excluded, recorded 20.2 (unrounded 20.2279, agrees); fusion, BRCA1 excluded, assay-only 13.7 (unrounded 13.6559, agrees); CADD, BRCA1 excluded, recorded 19.1 (unrounded 19.0921, agrees); CADD, BRCA1 excluded, assay-only 9.9 (unrounded 9.8547, agrees); fusion, assay-only VUS 12.1 (unrounded 12.0896, agrees). The exemption now holds only for numbers that open a heading, in the forward check and in the replacement tool, and `tests/test_analysis_claims.py::test_the_section_number_exemption_covers_headings_not_measurements` pins it.

**Checked and clean.** Every other cell of Tables 1 and 3; all of Table 2, including its two ties (the elastic net's drop-all-evolution lower bound, stored 0.0005, is 0.0005103 in a replay of phase 2's bootstrap, so +0.001 stands; the gradient-boosted drop-conservation ablated ρ, stored 0.7595, recomputes to 0.7594899, so 0.759 stands); Table 4, whose source stores unrounded values; the pooled ρ quoted in 3.1 against `phase2_pooled_rho.csv`; and the other numbers of the Figure 2 caption.
