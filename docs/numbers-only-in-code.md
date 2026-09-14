# Measured numbers that no output carries, and copies of them in code

Manuscript `calibration_draft.docx`, 413 measured numbers. Written by `python scripts/numbers_only_in_code.py <docx>`. The number check's own record of every token is `docs/manuscript-binding.md`.

## 1. Measured numbers with no persisted source: 0

A measured number is persisted when the number check binds it to a cell of an output its paragraph or table declares, or to a derivation the check computes from those outputs. A pool match, a whitelisted value, or a derivation declared in the claims file without being computed has no output that carries it.

None.

## 2. Found and resolved in the binding round

| printed | what it is | where it lived | resolution | binding now |
|---|---|---|---|---|
| 0.997 | TP53 control AUROC after orientation; printed 0.998 until the rounding round | a config comment; the value check matched 0.998 to BRCA2's 0.9975 | phase4_external_tp53.py writes it to phase4_tp53_orientation.csv (0.996925) on every run and the text prints 0.997; a cell now binds only if it rounds to the print, so 0.9975 no longer does | P69: one of 13 cells |
| 21,409 | SpliceAI values that re-round to the v1 print (Methods 2.5) | an assertion in tests/test_frozen_matrix_v2.py, and the whitelist | phase1_build_frozen_matrix_v2.py writes it to frozen_matrix_v2_column_report.tsv (n_rerounded_to_v1_print) | P33: single cell |
| 0.2149994 | full-precision score of the one re-rounding exception | the same test, and the whitelist | written to frozen_matrix_v2_column_report.tsv (exception_v2_value) | P33: single cell |
| 0.9997 | rank correlation of the re-scored Nucleotide Transformer column with v1 | frozen_matrix_v2_column_report.tsv, outside the outputs the check read; whitelisted | the paragraph declares the column report; bound to its nt row | P33: one of 2 cells; P36: single cell |
| 0.672 | agreement of AlphaGenome client v0.6.1 with the companion's v0.7.0 | docs/column-provenance.md, the page a one-off audit script writes | scripts/column_provenance_audit.py also writes data/rescore/column_provenance_shared.tsv (alphagenome, spearman) | P33: single cell |
| 72 | fusion-versus-single-tool comparisons whose interval excludes zero | no cell; the value check matched it to unrelated 0.72 values read as a percentage | a percentage reading now needs a printed %; the count is computed from phase3_pertool_brier_ci.csv | P53: computed derivation |
| 30 | disagreements between the two label sets | the claims file, as typed arithmetic (18 + 12) | computed from phase7_label_confusion.csv (its off-diagonal cells) | P65: computed derivation |
| 1.44 | smallest LR+ ratio, ClinVar-recorded over assay-only | the claims file, as typed arithmetic (20.333 / 14.141) | computed from phase7_selection_test.csv (min_ratio) | P65: computed derivation |
| 4.20 | largest such ratio | the claims file, as typed arithmetic (12.518 / 2.979) | computed from phase7_selection_test.csv (max_ratio) | P65: computed derivation |
| 46,392 | SNVs in the companion atlas's matrix | the whitelist | a property of the companion study's data, not this pipeline's: set aside as cited, with its reference | P18: not a measurement (cited) |

## 3. Distinctive measured numbers with a copy outside the outputs: 145 of 256

The 256 measured numbers with three or more significant digits, or of 100 or more, were searched for in 162 tracked source, configuration, test, README and hand-written documentation files; pipeline outputs and generated pages are excluded. A shorter literal such as 7 or 0.5 occurs in hundreds of unrelated lines, so a copy of it says nothing about where a number came from; section 1 covers those numbers too.

| block | number | binding | copies (first three) | context |
|---|---|---|---|---|
| P7 | 1,781 | one of 3 cells | README.md:38, README.md:73, README.md:85 (+47) | …ds. Ten predictors and three fusions were evaluated on 1,781 splice-region varia… |
| P8 | 74.9 | single cell | phase1/config/analysis_claims.json:130 | …our; in the primary condition the interval spans zero (74.9% against 73.6%). Acc… |
| P9 | 1.61 | one of 2 cells | README.md:46 | …fusion's advantage narrows to ΔLR+ 1.0 (95% CI 0.04 to 1.61) — one unit inside t… |
| P18 | 21,410 | one of 2 cells | README.md:9, README.md:211, docs/NOTE_S9.md:18 (+45) | …The analysis reuses a frozen matrix of 21,410 GRCh38 single-nucleotide variants… |
| P18 | 21,394 | single cell | docs/column-provenance.md:3, docs/column-provenance.md:26, docs/column-provenance.md:27 (+13) | …ncertain, conflicting and other. Of these 21,410 SNVs, 21,394 are shared with th… |
| P18 | 0.868 | single cell | docs/prose-density.md:158 | …against synonymous/nonsense controls (control AUROC ≥ 0.868 in every gene) (Supp… |
| P18 | 583 | one of 2 cells | docs/prose-density.md:42, docs/prose-density.md:107, phase1/src/check_manuscript_numbers.py:113 (+2) | …ron-side splice region: splice-core (\|offset\| ≤ 2; n = 583, including six audite… |
| P18 | 1,191 | one of 3 cells | docs/prose-density.md:42, docs/prose-density.md:107, phase1/src/check_manuscript_numbers.py:113 (+2) | …on classifier had missed) plus splice-region (3–8; n = 1,191) plus seven exon-si… |
| P18 | 1,781 | one of 4 cells | README.md:38, README.md:73, README.md:85 (+47) | …s them as coding/UTR — and are grouped with the core — 1,781 in all, so the regi… |
| P18 | 590 | single cell | docs/prose-density.md:42, docs/prose-density.md:107 | …region-stratified H1 analysis contrasts core-like (n = 590) with splice-region (… |
| P19 | 1,675 | single cell | docs/lr-stability-check.md:114, docs/lr-stability-check.md:115, docs/lr-stability-check.md:118 (+11) | …ing = 1, normal = 0, indeterminate = NA). This covered 1,675 of 1,781 splice var… |
| P19 | 1,781 | one of 3 cells | README.md:38, README.md:73, README.md:85 (+47) | …normal = 0, indeterminate = NA). This covered 1,675 of 1,781 splice variants (94… |
| P29 | 95.0 | one of 38 cells | docs/lr-stability-check.md:313, docs/lr-stability-check.md:408, phase1/src/phase5_likelihood_ratios.py:24 | …bserved score value the realised specificity runs from 95.0% to 98.0%, so ratios… |
| P29 | 24.3 | one of 2 cells | docs/prose-density.md:53, docs/prose-density.md:147, docs/prose-density.md:164 (+1) | …, so ratios a little above 20 do occur: the largest is 24.3 on the functional st… |
| P29 | 97.2 | one of 2 cells | docs/prose-density.md:147 | …ratios on the functional standard (AlphaGenome 24.3 at 97.2% realised, Pangolin … |
| P29 | 19.9 | one of 39 cells | docs/lr-stability-check.md:3, phase1/src/phase8_lr_decomposition.py:5, scripts/weak_token_spotcheck.py:66 (+1) | …standard (AlphaGenome 24.3 at 97.2% realised, Pangolin 19.9 at 96.0%) partly ref… |
| P30 | 1.26 | one of 2 cells | docs/prose-density.md:181, phase1/config/analysis_claims.json:46 | …ction is strongly enriched for large effects (mean \|z\| 1.26 against 0.62 for var… |
| P30 | 1.24 | one of 2 cells | docs/lr-stability-check.md:54, docs/lr-stability-check.md:358, docs/prose-density.md:37 (+1) | …\|z\| decile distribution gives weights between 0.72 and 1.24 and a Kish effective… |
| P30 | 1,768 | one of 7 cells | docs/prose-density.md:37, docs/prose-density.md:50, docs/prose-density.md:135 (+6) | …0.0005, +0.0226] — an interval that spans zero on this 1,768-variant baseline, a… |
| P30 | 1,781 | single cell | README.md:38, README.md:73, README.md:85 (+47) | …as it does unweighted ([−0.0004, +0.0218]), where the 1,781-variant headline of … |
| P33 | 21,410 | one of 10 cells | README.md:9, README.md:211, docs/NOTE_S9.md:18 (+45) | …eproduces the published frozen-matrix-v1 column at all 21,410 variants, and the … |
| P33 | 21,409 | single cell | phase1/config/manuscript_number_whitelist.json:2, phase1/config/manuscript_number_whitelist.json:63, phase1/config/manuscript_number_whitelist.json:65 | …all 21,410 variants, and the pinned SpliceAI column at 21,409 of them. The excep… |
| P33 | 0.2149994 | single cell | phase1/config/manuscript_number_whitelist.json:2, phase1/config/manuscript_number_whitelist.json:63, phase1/config/manuscript_number_whitelist.json:66 (+1) | …ue (chr17:58696697 C>T), whose full-precision score of 0.2149994 falls within 1 … |
| P33 | 0.9997 | one of 2 cells | README.md:198, docs/column-provenance.md:13, docs/column-provenance.md:30 (+7) | …olumn agrees with the original at rank correlation ρ = 0.9997. The other predict… |
| P36 | 0.868 | single cell | docs/prose-density.md:158 | …ses in every gene (synonymous/nonsense control AUROC ≥ 0.868). The SpliceAI, Pan… |
| P36 | 0.9997 | single cell | README.md:198, docs/column-provenance.md:13, docs/column-provenance.md:30 (+7) | …eement to the printed precision for the first two, ρ = 0.9997 for the third. The… |
| P37 | 413 | one of 2 cells | phase1/config/pipeline_facts.json:64, phase1/config/pipeline_facts.json:65 | …t, rank or ratio computed from them, rounds to it. All 413 measured numbers are … |
| P37 | 183 | single cell | phase1/config/pipeline_facts.json:66 | …hem, rounds to it. All 413 measured numbers are bound: 183 to a single cell and … |
| P37 | 230 | single cell | phase1/config/pipeline_facts.json:67 | …3 measured numbers are bound: 183 to a single cell and 230 to one of several cel… |
| P37 | 148 | single cell | docs/lr-stability-check.md:254, docs/lr-stability-check.md:255, docs/milestone3_analysis_ready.md:64 (+1) | …firms the value but not which cell is meant. The other 148 numeric tokens — sect… |
| P41 | 0.761 | single cell | README.md:44, docs/milestone8_results.md:21, docs/milestone8b_results.md:14 (+2) | …t, high band (Table 1): pooled per-gene Spearman ρ was 0.761 (Pangolin), 0.752 (… |
| P41 | 0.750 | single cell | docs/milestone8_results.md:23, docs/milestone8b_results.md:16, docs/prose-density.md:84 (+1) | …Spearman ρ was 0.761 (Pangolin), 0.752 (SpliceAI) and 0.750 (AlphaGenome), then … |
| P41 | 0.694 | single cell | docs/milestone8_results.md:24, docs/milestone8b_results.md:17, docs/prose-density.md:84 (+1) | …, 0.752 (SpliceAI) and 0.750 (AlphaGenome), then CADD (0.694), GPN-MSA (0.665, o… |
| P41 | 0.665 | single cell | docs/prose-density.md:84, docs/prose-density.md:106 | …) and 0.750 (AlphaGenome), then CADD (0.694), GPN-MSA (0.665, oriented), phyloP … |
| P41 | 0.639 | one of 3 cells | docs/milestone8_results.md:25, docs/milestone8_results.md:26, docs/milestone8b_results.md:18 (+3) | …then CADD (0.694), GPN-MSA (0.665, oriented), phyloP (0.639), phastCons (0.612),… |
| P41 | 0.612 | single cell | docs/milestone8_results.md:27, docs/milestone8b_results.md:20, docs/prose-density.md:84 (+1) | …GPN-MSA (0.665, oriented), phyloP (0.639), phastCons (0.612), Nucleotide Transfo… |
| P41 | 0.186 | single cell | docs/milestone8_results.md:28, docs/milestone8b_results.md:22, docs/prose-density.md:84 (+1) | …(0.612), Nucleotide Transformer (0.492) and gnomAD AF (0.186). AlphaMissense cou… |
| P42 | 0.774 | one of 3 cells | README.md:44, docs/lr-stability-check.md:49 | …The elastic-net fusion reached ρ = 0.774, exceeding the best single tool (P… |
| P42 | 128 | one of 5 cells | docs/column-provenance.md:13, docs/lr-stability-check.md:117, docs/lr-stability-check.md:422 (+11) | …s point the same way. An exact sign-flip test over all 128 assignments of the se… |
| P42 | +0.0217 | single cell | docs/prose-density.md:67, docs/prose-density.md:129 | …splice core (+0.0337) than in the splice-region band (+0.0217), but neither inte… |
| P42 | +0.118 | single cell | docs/prose-density.md:153 | …), but neither interval excludes zero (95% CI [−0.035, +0.118] and [−0.0056, +0.… |
| P44 | 1,781 | one of 3 cells | README.md:38, README.md:73, README.md:85 (+47) | …) between prediction and functional pathogenicity over 1,781 splice-region varia… |
| P46 | 0.730 | single cell | docs/lr-stability-check.md:20 | …s widens the fusion's interval from [0.741, 0.803] to [0.730, 0.811] without cha… |
| P46 | 0.811 | single cell | docs/lr-stability-check.md:37 | …s the fusion's interval from [0.741, 0.803] to [0.730, 0.811] without changing t… |
| P53 | +0.0103 | one of 4 cells | README.md:45 | …ary condition, the functional standard with BRCA1 (Δ = +0.0103, 95% CI [+0.0047,… |
| P53 | +0.0182 | single cell | README.md:45 | …nal standard with BRCA1 (Δ = +0.0103, 95% CI [+0.0047, +0.0182]); six of the sev… |
| P53 | +0.0117 | one of 5 cells | docs/prose-density.md:66, docs/prose-density.md:169 | …the interval likewise excluded zero without BRCA1 (Δ = +0.0117, [+0.0042, +0.018… |
| P53 | +0.0187 | one of 2 cells | docs/prose-density.md:121 | …se excluded zero without BRCA1 (Δ = +0.0117, [+0.0042, +0.0187]) and on ClinVar … |
| P53 | 1,768 | one of 6 cells | docs/prose-density.md:37, docs/prose-density.md:50, docs/prose-density.md:135 (+6) | …t sensitivity behind the modest difference between the 1,768-variant reweighting… |
| P53 | 1,781 | single cell | README.md:38, README.md:73, README.md:85 (+47) | …between the 1,768-variant reweighting baseline and the 1,781-variant headline in… |
| P53 | 74.9 | single cell | phase1/config/analysis_claims.json:130 | …ional standard including BRCA1 — the margin was small: 74.9% high-confidence at … |
| P53 | 96.7 | single cell | docs/prose-density.md:183, phase1/config/analysis_claims.json:163 | …BRCA1 — the margin was small: 74.9% high-confidence at 96.7% accuracy versus Pan… |
| P53 | 96.3 | one of 3 cells | docs/prose-density.md:184 | …onfidence at 96.7% accuracy versus Pangolin's 73.6% at 96.3% (1,254 against 1,23… |
| P54 | 0.0630 | one of 3 cells | README.md:45 | …reliability, resolution and uncertainty. The fusion's 0.0630 against Pangolin's … |
| P54 | 0.0733 | one of 3 cells | README.md:45, tests/test_analysis_claims.py:271 | …nd uncertainty. The fusion's 0.0630 against Pangolin's 0.0733 — the headline ΔBr… |
| P54 | +0.0103 | one of 4 cells | README.md:45 | …630 against Pangolin's 0.0733 — the headline ΔBrier is +0.0103 [+0.0047, +0.0182… |
| P54 | +0.0182 | one of 3 cells | README.md:45 | …in's 0.0733 — the headline ΔBrier is +0.0103 [+0.0047, +0.0182] — divides as fol… |
| P54 | 0.2499 | one of 6 cells | docs/prose-density.md:77, docs/prose-density.md:111 | …divides as follows: the uncertainty term is identical (0.2499); the reliability … |
| P54 | 0.0103 | one of 4 cells | README.md:45 | …lementary Table S15. Rounded table entries subtract to 0.0103; the binned reliab… |
| P55 | 0.118 | single cell | docs/prose-density.md:153 | …icing sharpness — GPN-MSA reaches ECE 0.025 with Brier 0.118 and 47% high-confid… |
| P57 | 1,675 | single cell | docs/lr-stability-check.md:114, docs/lr-stability-check.md:115, docs/lr-stability-check.md:118 (+11) | …ins (mean predicted versus observed fraction damaging; 1,675 variants carrying a… |
| P57 | 0.0630 | single cell | README.md:45 | …imary condition the fusion attains ECE 0.025 and Brier 0.0630 (Pangolin 0.032 an… |
| P57 | 0.0733 | single cell | README.md:45, tests/test_analysis_claims.py:271 | …attains ECE 0.025 and Brier 0.0630 (Pangolin 0.032 and 0.0733; Table 3). The two… |
| P57 | 0.975 | one of 2 cells | phase1/src/build_supp_tables.py:378, phase1/src/phase2_model.py:151, phase1/src/phase3_calibration.py:180 (+6) | …gonal closely (0.026 predicted against 0.027 observed; 0.975 against 0.962).… |
| P57 | 0.962 | one of 2 cells | docs/prose-density.md:193 | …(0.026 predicted against 0.027 observed; 0.975 against 0.962).… |
| P59 | 1,675 | single cell | docs/lr-stability-check.md:114, docs/lr-stability-check.md:115, docs/lr-stability-check.md:118 (+11) | …the seven genes), across the functional standard (n = 1,675 with BRCA1, 1,168 wi… |
| P59 | 1,168 | single cell | docs/lr-stability-check.md:162, docs/lr-stability-check.md:163, docs/lr-stability-check.md:166 (+4) | …across the functional standard (n = 1,675 with BRCA1, 1,168 without) and ClinVar… |
| P59 | 946 | single cell | docs/lr-stability-check.md:210, docs/lr-stability-check.md:211, docs/lr-stability-check.md:214 (+6) | …(n = 1,675 with BRCA1, 1,168 without) and ClinVar (n = 946 and 734). Positive fa… |
| P59 | 734 | single cell | docs/lr-stability-check.md:258, docs/lr-stability-check.md:259, docs/lr-stability-check.md:262 (+4) | …75 with BRCA1, 1,168 without) and ClinVar (n = 946 and 734). Positive favours th… |
| P62 | 19.9 | one of 53 cells | docs/lr-stability-check.md:3, phase1/src/phase8_lr_decomposition.py:5, scripts/weak_token_spotcheck.py:66 (+1) | …thresholds two point estimates cross it — Pangolin at 19.9 (achieved specificity… |
| P62 | 24.3 | one of 4 cells | docs/prose-density.md:53, docs/prose-density.md:147, docs/prose-density.md:164 (+1) | …nic calibration; interpolated 15.9) and AlphaGenome at 24.3 (97.2%; BRCA1 exclud… |
| P62 | 97.2 | one of 6 cells | docs/prose-density.md:147 | …libration; interpolated 15.9) and AlphaGenome at 24.3 (97.2%; BRCA1 excluded, ca… |
| P62 | 17.78 | one of 5 cells | docs/lr-stability-check.md:319 | …s interpolated LR+ in the primary condition runs 9.35, 17.78, 31.33 and 63.82, c… |
| P62 | 17.8 | one of 11 cells | README.md:46 | …hable (2.4). At 95% specificity the fusion attains LR+ 17.8 (95% CI 16.3–18.6) a… |
| P62 | 0.889 | one of 6 cells | docs/lr-stability-check.md:11, docs/prose-density.md:41, docs/prose-density.md:64 (+2) | …ttains LR+ 17.8 (95% CI 16.3–18.6) at a sensitivity of 0.889, the best single to… |
| P62 | 16.8 | one of 11 cells | README.md:46, docs/prose-density.md:64, docs/prose-density.md:141 | …–18.6) at a sensitivity of 0.889, the best single tool 16.8 (15.6–18.5), and ten… |
| P63 | 1.02 | one of 5 cells | docs/prose-density.md:44, docs/prose-density.md:161 | …one unit on the functional standard — with BRCA1, ΔLR+ 1.02 (95% CI 0.04 to 1.61… |
| P63 | 1.61 | one of 2 cells | README.md:46 | …ional standard — with BRCA1, ΔLR+ 1.02 (95% CI 0.04 to 1.61) on raw scores and 1… |
| P63 | 1.96 | single cell | phase1/config/manuscript_number_whitelist.json:15, phase1/src/phase2_model.py:128, phase1/src/phase2_model.py:154 (+5) | …(95% CI 0.04 to 1.61) on raw scores and 1.64 (0.07 to 1.96) after isotonic calib… |
| P63 | 19.9 | one of 21 cells | docs/lr-stability-check.md:3, phase1/src/phase8_lr_decomposition.py:5, scripts/weak_token_spotcheck.py:66 (+1) | …angolin's threshold lands at 96.0% specificity and its 19.9 crosses Strong while… |
| P63 | 17.6 | one of 5 cells | docs/prose-density.md:62, docs/prose-density.md:114 | …, does not; the same cell interpolates to 15.9 against 17.6, and the difference … |
| P64 | 1,675 | single cell | docs/lr-stability-check.md:114, docs/lr-stability-check.md:115, docs/lr-stability-check.md:118 (+11) | …ndard none does — but those are not the same variants: 1,675 carry an assay labe… |
| P64 | 946 | single cell | docs/lr-stability-check.md:210, docs/lr-stability-check.md:211, docs/lr-stability-check.md:214 (+6) | …are not the same variants: 1,675 carry an assay label, 946 a ClinVar label, and … |
| P64 | 907 | one of 19 cells | docs/prose-density.md:94, docs/prose-density.md:183, phase1/config/analysis_claims.json:121 (+3) | …: 1,675 carry an assay label, 946 a ClinVar label, and 907 carry both. Restricte… |
| P65 | 907 | one of 37 cells | docs/prose-density.md:94, docs/prose-density.md:183, phase1/config/analysis_claims.json:121 (+3) | …the two label sets is excluded by measurement: on the 907 shared variants they a… |
| P65 | 96.7 | one of 3 cells | docs/prose-density.md:183, phase1/config/analysis_claims.json:163 | …measurement: on the 907 shared variants they agree on 96.7% of calls (Cohen's κ … |
| P65 | 1.44 | computed derivation | docs/prose-density.md:41, docs/prose-density.md:136, phase1/config/analysis_claims.json:192 | …of VUS, conflicting and other records — by ratios from 1.44 to 4.20, and the fus… |
| P65 | 4.20 | computed derivation | docs/milestone3_analysis_ready.md:31, docs/milestone6_spliceai_pangolin.md:62, docs/milestone7_dnalm_alphagenome.md:69 (+2) | …conflicting and other records — by ratios from 1.44 to 4.20, and the fusion's se… |
| P65 | 0.971 | one of 13 cells | docs/lr-stability-check.md:73, docs/lr-stability-check.md:98, docs/prose-density.md:41 (+1) | …the fusion's sensitivity at 95% specificity falls from 0.971 to 0.681, CADD's fr… |
| P65 | 0.681 | one of 3 cells | docs/prose-density.md:41, docs/prose-density.md:136 | …n's sensitivity at 95% specificity falls from 0.971 to 0.681, CADD's from 0.908 … |
| P65 | 0.908 | one of 4 cells | docs/prose-density.md:41, docs/prose-density.md:136 | …95% specificity falls from 0.971 to 0.681, CADD's from 0.908 to 0.327 (operating… |
| P65 | 0.327 | one of 5 cells | docs/prose-density.md:41, docs/prose-density.md:136 | …ficity falls from 0.971 to 0.681, CADD's from 0.908 to 0.327 (operating points a… |
| P65 | 0.889 | one of 4 cells | docs/lr-stability-check.md:11, docs/prose-density.md:41, docs/prose-density.md:64 (+2) | …he subset sensitivities need not average to the pooled 0.889). The gap survives … |
| P65 | 20.2 | one of 6 cells | docs/prose-density.md:200, scripts/weak_token_spotcheck.py:176 | …led 0.889). The gap survives BRCA1 removal (fusion LR+ 20.2 against 13.7; CADD 1… |
| P65 | 13.7 | one of 6 cells | docs/prose-density.md:200, scripts/weak_token_spotcheck.py:177 | …he gap survives BRCA1 removal (fusion LR+ 20.2 against 13.7; CADD 19.1 against 9… |
| P65 | 19.1 | one of 17 cells | docs/prose-density.md:200, scripts/weak_token_spotcheck.py:178 | …ives BRCA1 removal (fusion LR+ 20.2 against 13.7; CADD 19.1 against 9.9). Two st… |
| P65 | 12.1 | one of 8 cells | docs/prose-density.md:62, docs/prose-density.md:114, scripts/weak_token_spotcheck.py:180 | …d within the assay-only arm the VUS subset (fusion LR+ 12.1) and the conflicting… |
| P67 | 74.9 | single cell | phase1/config/analysis_claims.json:130 | …Moderate or above on the functional standard, against 74.9% under the 0.90/0.10 … |
| P67 | 21.0 | single cell | docs/prose-density.md:78, environment.yml:17, requirements.txt:11 | …rp scores: in that same condition phastCons moves from 21.0% to 75.0%, because a… |
| P67 | 75.0 | single cell | docs/prose-density.md:78 | …: in that same condition phastCons moves from 21.0% to 75.0%, because a confiden… |
| P69 | 192 | one of 7 cells | README.md:14, docs/frozen-matrix-v2.md:33, docs/milestone4_candidate_genes.md:38 (+18) | …seven genes, then applied without refitting to TP53's 192 intron-side splice SNV… |
| P69 | 0.997 | one of 13 cells | docs/prose-density.md:198, phase1/src/check_manuscript_numbers.py:571, phase1/src/config.py:108 (+6) | …d gate used throughout (post-orientation control AUROC 0.997) (Supplementary Tab… |
| P71 | 0.118 | single cell | docs/prose-density.md:153 | …r than the best single tool's (Pangolin): 0.077 versus 0.118 (ΔBrier 0.041; vari… |
| P71 | 0.122 | single cell | docs/prose-density.md:89, docs/prose-density.md:119 | …95% CI 0.025–0.059), with lower calibration error (ECE 0.122 versus 0.157) and h… |
| P71 | 0.968 | single cell | docs/lr-stability-check.md:87 | …and higher accuracy within the high-confidence subset (0.968 versus 0.940), whic… |
| P72 | 192 | one of 2 cells | README.md:14, docs/frozen-matrix-v2.md:33, docs/milestone4_candidate_genes.md:38 (+18) | …gene rather than a second significance claim: with n = 192 from one gene there i… |
| P79 | 0.218 | one of 2 cells | docs/milestone8b_results.md:38 | …d SpliceAI carry the weight in every fold (mean 0.224, 0.218 and 0.186, SD 0.019… |
| P79 | 0.186 | one of 2 cells | docs/milestone8_results.md:28, docs/milestone8b_results.md:22, docs/prose-density.md:84 (+1) | …carry the weight in every fold (mean 0.224, 0.218 and 0.186, SD 0.019, 0.029 and… |
| P79 | 0.0596 | one of 2 cells | phase1/config/analysis_claims.json:152 | …s the calibration only slightly: mean Brier moves from 0.0596 with six training … |
| T1r1c2 | 0.774 | single cell | README.md:44, docs/lr-stability-check.md:49 | …0.774… |
| T1r2c2 | 0.761 | single cell | README.md:44, docs/milestone8_results.md:21, docs/milestone8b_results.md:14 (+2) | …0.761… |
| T1r2c3 | 0.731 | one of 4 cells | docs/milestone8_results.md:24, docs/milestone8b_results.md:17 | …0.731 – 0.788… |
| T1r2c3 | 0.788 | one of 2 cells | docs/lr-stability-check.md:38 | …0.731 – 0.788… |
| T1r4c3 | 0.783 | one of 2 cells | docs/lr-stability-check.md:39 | …0.716 – 0.783… |
| T1r6c2 | 0.750 | single cell | docs/milestone8_results.md:23, docs/milestone8b_results.md:16, docs/prose-density.md:84 (+1) | …0.750… |
| T1r6c3 | 0.714 | one of 2 cells | docs/milestone8_results.md:23, docs/milestone8b_results.md:16 | …0.714 – 0.782… |
| T1r6c3 | 0.782 | one of 2 cells | docs/milestone8_results.md:23, docs/milestone8b_results.md:16 | …0.714 – 0.782… |
| T1r7c2 | 0.694 | single cell | docs/milestone8_results.md:24, docs/milestone8b_results.md:17, docs/prose-density.md:84 (+1) | …0.694… |
| T1r7c3 | 0.654 | one of 2 cells | docs/milestone8_results.md:24, docs/milestone8b_results.md:17 | …0.654 – 0.731… |
| T1r7c3 | 0.731 | one of 4 cells | docs/milestone8_results.md:24, docs/milestone8b_results.md:17 | …0.654 – 0.731… |
| T1r8c2 | 0.665 | single cell | docs/prose-density.md:84, docs/prose-density.md:106 | …0.665… |
| T1r8c3 | 0.639 | one of 3 cells | docs/milestone8_results.md:25, docs/milestone8_results.md:26, docs/milestone8b_results.md:18 (+3) | …0.639 – 0.691… |
| T1r8c3 | 0.691 | one of 2 cells | docs/milestone8_results.md:25, docs/milestone8b_results.md:18 | …0.639 – 0.691… |
| T1r9c2 | 0.639 | one of 3 cells | docs/milestone8_results.md:25, docs/milestone8_results.md:26, docs/milestone8b_results.md:18 (+3) | …0.639… |
| T1r9c3 | 0.607 | one of 2 cells | docs/milestone8_results.md:26, docs/milestone8b_results.md:19 | …0.607 – 0.669… |
| T1r9c3 | 0.669 | one of 2 cells | docs/lr-stability-check.md:31, docs/milestone8_results.md:26, docs/milestone8b_results.md:19 | …0.607 – 0.669… |
| T1r10c2 | 0.612 | single cell | docs/milestone8_results.md:27, docs/milestone8b_results.md:20, docs/prose-density.md:84 (+1) | …0.612… |
| T1r10c3 | 0.563 | one of 2 cells | docs/milestone8_results.md:27, docs/milestone8b_results.md:20 | …0.563 – 0.657… |
| T1r10c3 | 0.657 | one of 2 cells | docs/milestone8_results.md:27, docs/milestone8b_results.md:20 | …0.563 – 0.657… |
| T1r12c2 | 0.186 | single cell | docs/milestone8_results.md:28, docs/milestone8b_results.md:22, docs/prose-density.md:84 (+1) | …0.186… |
| T1r12c3 | 0.298 | one of 2 cells | docs/milestone8_results.md:28, docs/milestone8b_results.md:22 | …0.069 – 0.298… |
| T2r1c2 | 0.774 | one of 2 cells | README.md:44, docs/lr-stability-check.md:49 | …0.774… |
| T2r1c3 | 0.773 | single cell | docs/lr-stability-check.md:19, docs/prose-density.md:59, docs/prose-density.md:117 | …0.773… |
| T2r2c2 | 0.774 | one of 2 cells | README.md:44, docs/lr-stability-check.md:49 | …0.774… |
| T2r2c3 | 0.771 | single cell | docs/lr-stability-check.md:50 | …0.771… |
| T2r3c3 | 0.759 | single cell | docs/lr-stability-check.md:42, docs/lr-stability-check.md:311, scripts/weak_token_spotcheck.py:247 | …0.759… |
| T3r1c4 | 74.9 | single cell | phase1/config/analysis_claims.json:130 | …74.9… |
| T3r1c5 | 96.7 | single cell | docs/prose-density.md:183, phase1/config/analysis_claims.json:163 | …96.7… |
| T3r2c5 | 96.3 | one of 3 cells | docs/prose-density.md:184 | …96.3… |
| T3r5c5 | 96.3 | one of 3 cells | docs/prose-density.md:184 | …96.3… |
| T3r7c3 | 0.118 | single cell | docs/prose-density.md:153 | …0.118… |
| T3r9c4 | 21.0 | single cell | docs/prose-density.md:78, environment.yml:17, requirements.txt:11 | …21.0… |
| T4r1c1 | 192 | one of 2 cells | README.md:14, docs/frozen-matrix-v2.md:33, docs/milestone4_candidate_genes.md:38 (+18) | …192… |
| T4r2c1 | 192 | one of 2 cells | README.md:14, docs/frozen-matrix-v2.md:33, docs/milestone4_candidate_genes.md:38 (+18) | …192… |
| T4r3c1 | 181 | single cell | docs/lr-stability-check.md:206, docs/lr-stability-check.md:207 | …181… |
