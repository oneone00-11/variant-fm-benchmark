# Measured numbers that no output carries, and copies of them in code

Manuscript `calibration_draft.docx`, 410 measured numbers. Written by `python scripts/numbers_only_in_code.py <docx>`. The number check's own record of every token is `docs/manuscript-binding.md`.

## 1. Measured numbers with no persisted source: 0

A measured number is persisted when the number check binds it to a cell of an output its paragraph or table declares, or to a derivation the check computes from those outputs. A pool match, a whitelisted value, or a derivation declared in the claims file without being computed has no output that carries it.

None.

## 2. Found and resolved in the binding round

| printed | what it is | where it lived | resolution | binding now |
|---|---|---|---|---|
| 0.997 | TP53 control AUROC after orientation; printed 0.998 until the rounding round | a config comment; the value check matched 0.998 to BRCA2's 0.9975 | phase4_external_tp53.py writes it to phase4_tp53_orientation.csv (0.996925) on every run and the text prints 0.997; a cell now binds only if it rounds to the print, so 0.9975 no longer does, and the paragraph's claim pins the number to the orientation cell | P75: single cell |
| 21,409 | SpliceAI values that re-round to the v1 print (Methods 2.5) | an assertion in tests/test_frozen_matrix_v2.py, and the whitelist | phase1_build_frozen_matrix_v2.py writes it to frozen_matrix_v2_column_report.tsv (n_rerounded_to_v1_print) | P39: single cell |
| 0.2149994 | full-precision score of the one re-rounding exception | the same test, and the whitelist | written to frozen_matrix_v2_column_report.tsv (exception_v2_value) | P39: single cell |
| 0.9997 | rank correlation of the re-scored Nucleotide Transformer column with v1 | frozen_matrix_v2_column_report.tsv, outside the outputs the check read; whitelisted | the paragraph declares the column report; bound to its nt row | P39: one of 2 cells; P42: single cell |
| 0.672 | agreement of AlphaGenome client v0.6.1 with the companion's v0.7.0 | docs/column-provenance.md, the page a one-off audit script writes | scripts/column_provenance_audit.py also writes data/rescore/column_provenance_shared.tsv (alphagenome, spearman) | P39: single cell |
| 72 | fusion-versus-single-tool comparisons whose interval excludes zero | no cell; the value check matched it to unrelated 0.72 values read as a percentage | a percentage reading now needs a printed %; the count is computed from phase3_pertool_brier_ci.csv | P59: computed derivation |
| 30 | disagreements between the two label sets | the claims file, as typed arithmetic (18 + 12) | computed from phase7_label_confusion.csv (its off-diagonal cells) | P71: computed derivation |
| 1.44 | smallest LR+ ratio, ClinVar-recorded over assay-only | the claims file, as typed arithmetic (20.333 / 14.141) | computed from phase7_selection_test.csv (min_ratio) | P71: computed derivation |
| 4.20 | largest such ratio | the claims file, as typed arithmetic (12.518 / 2.979) | computed from phase7_selection_test.csv (max_ratio) | P71: computed derivation |
| 46,392 | SNVs in the companion atlas's matrix | the whitelist | a property of the companion study's data, not this pipeline's: set aside as cited, with its reference | P19: not a measurement (cited) |

## 3. Distinctive measured numbers with a copy outside the outputs: 117 of 253

The 253 measured numbers with three or more significant digits, or of 100 or more, were searched for in 164 tracked source, configuration, test, README and hand-written documentation files; pipeline outputs and generated pages are excluded. A shorter literal such as 7 or 0.5 occurs in hundreds of unrelated lines, so a copy of it says nothing about where a number came from; section 1 covers those numbers too.

| block | number | binding | copies (first three) | context |
|---|---|---|---|---|
| P7 | 1,781 | one of 3 cells | README.md:38, README.md:73, README.md:85 (+42) | …ds. Ten predictors and three fusions were evaluated on 1,781 splice-region varia… |
| P8 | 74.9 | single cell | phase1/config/analysis_claims.json:133 | …our; in the primary condition the interval spans zero (74.9% against 73.6%). Acc… |
| P9 | 1.61 | one of 2 cells | README.md:46 | …fusion's advantage narrows to ΔLR+ 1.0 (95% CI 0.04 to 1.61) — one unit inside t… |
| P18 | 21,410 | one of 2 cells | README.md:9, README.md:211, docs/NOTE_S9.md:10 (+42) | …The analysis reuses a frozen matrix of 21,410 GRCh38 single-nucleotide variants… |
| P19 | 21,410 | one of 2 cells | README.md:9, README.md:211, docs/NOTE_S9.md:10 (+42) | …Of these 21,410 SNVs, 21,394 are shared with the c… |
| P19 | 21,394 | single cell | docs/column-provenance.md:3, docs/column-provenance.md:26, docs/column-provenance.md:27 (+12) | …Of these 21,410 SNVs, 21,394 are shared with the companion atla… |
| P21 | 1,781 | one of 3 cells | README.md:38, README.md:73, README.md:85 (+42) | …The primary subset is the intron-side splice region, 1,781 variants in all. It c… |
| P21 | 583 | one of 2 cells | phase1/src/check_manuscript_numbers.py:113, tests/test_frozen_matrix.py:53, tests/test_frozen_matrix.py:55 | …n all. It comprises the splice core (\|offset\| ≤ 2; n = 583), the splice-region b… |
| P21 | 1,191 | one of 3 cells | phase1/src/check_manuscript_numbers.py:113, tests/test_frozen_matrix.py:53, tests/test_frozen_matrix.py:56 | …et\| ≤ 2; n = 583), the splice-region band (3–8 bp; n = 1,191), and seven exon-si… |
| P24 | 1,675 | single cell | docs/lr-stability-check.md:114, docs/lr-stability-check.md:115, docs/lr-stability-check.md:118 (+6) | …ing = 1, normal = 0, indeterminate = NA). This covered 1,675 of 1,781 splice var… |
| P24 | 1,781 | one of 3 cells | README.md:38, README.md:73, README.md:85 (+42) | …normal = 0, indeterminate = NA). This covered 1,675 of 1,781 splice variants (94… |
| P35 | 95.0 | one of 38 cells | docs/lr-stability-check.md:313, docs/lr-stability-check.md:408, phase1/src/phase5_likelihood_ratios.py:24 | …bserved score value the realised specificity runs from 95.0% to 98.0%, so ratios… |
| P35 | 24.3 | one of 2 cells | phase1/src/phase8_lr_decomposition.py:14 | …, so ratios a little above 20 do occur: the largest is 24.3 on the functional st… |
| P35 | 19.9 | one of 39 cells | docs/lr-stability-check.md:3, phase1/src/phase8_lr_decomposition.py:5, scripts/weak_token_spotcheck.py:66 (+1) | …standard (AlphaGenome 24.3 at 97.2% realised, Pangolin 19.9 at 96.0%) partly ref… |
| P36 | 1.26 | one of 2 cells | phase1/config/analysis_claims.json:46 | …ction is strongly enriched for large effects (mean \|z\| 1.26 against 0.62 for var… |
| P36 | 1.24 | one of 2 cells | docs/lr-stability-check.md:54, docs/lr-stability-check.md:358 | …\|z\| decile distribution gives weights between 0.72 and 1.24 and a Kish effective… |
| P36 | 1,768 | one of 7 cells | docs/NOTE_S9.md:16, phase1/config/analysis_claims.json:272, phase1/src/build_supp_ipw.py:102 (+3) | …0.0005, +0.0226] — an interval that spans zero on this 1,768-variant baseline, a… |
| P36 | 1,781 | single cell | README.md:38, README.md:73, README.md:85 (+42) | …as it does unweighted ([−0.0004, +0.0218]), where the 1,781-variant headline of … |
| P39 | 21,410 | one of 10 cells | README.md:9, README.md:211, docs/NOTE_S9.md:10 (+42) | …eproduces the published frozen-matrix-v1 column at all 21,410 variants, and the … |
| P39 | 21,409 | single cell | phase1/config/manuscript_number_whitelist.json:2, phase1/config/manuscript_number_whitelist.json:63, phase1/config/manuscript_number_whitelist.json:65 (+1) | …all 21,410 variants, and the pinned SpliceAI column at 21,409 of them. The excep… |
| P39 | 0.2149994 | single cell | phase1/config/manuscript_number_whitelist.json:2, phase1/config/manuscript_number_whitelist.json:63, phase1/config/manuscript_number_whitelist.json:66 (+2) | …ue (chr17:58696697 C>T), whose full-precision score of 0.2149994 falls within 1 … |
| P39 | 0.9997 | one of 2 cells | README.md:198, docs/column-provenance.md:13, docs/column-provenance.md:30 (+8) | …olumn agrees with the original at rank correlation ρ = 0.9997. The other predict… |
| P39 | 0.672 | single cell | scripts/numbers_only_in_code.py:55 | …ore definitions, and the two columns agree at only ρ = 0.672 (companion Note S11… |
| P42 | 0.9997 | single cell | README.md:198, docs/column-provenance.md:13, docs/column-provenance.md:30 (+8) | …eement to the printed precision for the first two, ρ = 0.9997 for the third. The… |
| P47 | 0.761 | single cell | README.md:44, docs/milestone8_results.md:21, docs/milestone8b_results.md:14 | …t, high band (Table 1): pooled per-gene Spearman ρ was 0.761 (Pangolin), 0.752 (… |
| P47 | 0.750 | single cell | docs/milestone8_results.md:23, docs/milestone8b_results.md:16 | …Spearman ρ was 0.761 (Pangolin), 0.752 (SpliceAI) and 0.750 (AlphaGenome), then … |
| P47 | 0.694 | single cell | docs/milestone8_results.md:24, docs/milestone8b_results.md:17 | …, 0.752 (SpliceAI) and 0.750 (AlphaGenome), then CADD (0.694), GPN-MSA (0.665, o… |
| P47 | 0.639 | one of 3 cells | docs/milestone8_results.md:25, docs/milestone8_results.md:26, docs/milestone8b_results.md:18 (+1) | …then CADD (0.694), GPN-MSA (0.665, oriented), phyloP (0.639), phastCons (0.612),… |
| P47 | 0.612 | single cell | docs/milestone8_results.md:27, docs/milestone8b_results.md:20 | …GPN-MSA (0.665, oriented), phyloP (0.639), phastCons (0.612), Nucleotide Transfo… |
| P47 | 0.186 | single cell | docs/milestone8_results.md:28, docs/milestone8b_results.md:22 | …(0.612), Nucleotide Transformer (0.492) and gnomAD AF (0.186). AlphaMissense cou… |
| P48 | 0.774 | one of 3 cells | README.md:44, docs/lr-stability-check.md:49 | …The elastic-net fusion reached ρ = 0.774, exceeding the best single tool (P… |
| P48 | 128 | one of 5 cells | docs/column-provenance.md:13, docs/lr-stability-check.md:117, docs/lr-stability-check.md:422 (+11) | …s point the same way. An exact sign-flip test over all 128 assignments of the se… |
| P50 | 1,781 | one of 3 cells | README.md:38, README.md:73, README.md:85 (+42) | …) between prediction and functional pathogenicity over 1,781 splice-region varia… |
| P52 | 0.730 | single cell | docs/lr-stability-check.md:20 | …s widens the fusion's interval from [0.741, 0.803] to [0.730, 0.811] without cha… |
| P52 | 0.811 | single cell | docs/lr-stability-check.md:37 | …s the fusion's interval from [0.741, 0.803] to [0.730, 0.811] without changing t… |
| P59 | +0.0103 | one of 4 cells | README.md:45 | …ary condition, the functional standard with BRCA1 (Δ = +0.0103, 95% CI [+0.0047,… |
| P59 | +0.0182 | single cell | README.md:45 | …nal standard with BRCA1 (Δ = +0.0103, 95% CI [+0.0047, +0.0182]); six of the sev… |
| P59 | 1,768 | one of 6 cells | docs/NOTE_S9.md:16, phase1/config/analysis_claims.json:272, phase1/src/build_supp_ipw.py:102 (+3) | …t sensitivity behind the modest difference between the 1,768-variant reweighting… |
| P59 | 1,781 | single cell | README.md:38, README.md:73, README.md:85 (+42) | …between the 1,768-variant reweighting baseline and the 1,781-variant headline in… |
| P59 | 74.9 | single cell | phase1/config/analysis_claims.json:133 | …ional standard including BRCA1 — the margin was small: 74.9% high-confidence at … |
| P59 | 96.7 | single cell | phase1/config/analysis_claims.json:166 | …BRCA1 — the margin was small: 74.9% high-confidence at 96.7% accuracy versus Pan… |
| P60 | 0.0630 | one of 3 cells | README.md:45 | …reliability, resolution and uncertainty. The fusion's 0.0630 against Pangolin's … |
| P60 | 0.0733 | one of 3 cells | README.md:45, tests/test_analysis_claims.py:271 | …nd uncertainty. The fusion's 0.0630 against Pangolin's 0.0733 — the headline ΔBr… |
| P60 | +0.0103 | one of 4 cells | README.md:45 | …630 against Pangolin's 0.0733 — the headline ΔBrier is +0.0103 [+0.0047, +0.0182… |
| P60 | +0.0182 | one of 3 cells | README.md:45 | …in's 0.0733 — the headline ΔBrier is +0.0103 [+0.0047, +0.0182] — divides as fol… |
| P60 | 0.0103 | one of 4 cells | README.md:45 | …lementary Table S15. Rounded table entries subtract to 0.0103; the binned reliab… |
| P63 | 1,675 | single cell | docs/lr-stability-check.md:114, docs/lr-stability-check.md:115, docs/lr-stability-check.md:118 (+6) | …ins (mean predicted versus observed fraction damaging; 1,675 variants carrying a… |
| P63 | 0.0630 | single cell | README.md:45 | …imary condition the fusion attains ECE 0.025 and Brier 0.0630 (Pangolin 0.032 an… |
| P63 | 0.0733 | single cell | README.md:45, tests/test_analysis_claims.py:271 | …attains ECE 0.025 and Brier 0.0630 (Pangolin 0.032 and 0.0733; Table 3). The two… |
| P63 | 0.975 | one of 2 cells | phase1/src/build_supp_tables.py:378, phase1/src/phase2_model.py:151, phase1/src/phase3_calibration.py:180 (+6) | …gonal closely (0.026 predicted against 0.027 observed; 0.975 against 0.962).… |
| P65 | 1,675 | single cell | docs/lr-stability-check.md:114, docs/lr-stability-check.md:115, docs/lr-stability-check.md:118 (+6) | …the seven genes), across the functional standard (n = 1,675 with BRCA1, 1,168 wi… |
| P65 | 1,168 | single cell | docs/lr-stability-check.md:162, docs/lr-stability-check.md:163, docs/lr-stability-check.md:166 (+2) | …across the functional standard (n = 1,675 with BRCA1, 1,168 without) and ClinVar… |
| P65 | 946 | single cell | docs/lr-stability-check.md:210, docs/lr-stability-check.md:211, docs/lr-stability-check.md:214 (+3) | …(n = 1,675 with BRCA1, 1,168 without) and ClinVar (n = 946 and 734). Positive fa… |
| P65 | 734 | single cell | docs/lr-stability-check.md:258, docs/lr-stability-check.md:259, docs/lr-stability-check.md:262 (+2) | …75 with BRCA1, 1,168 without) and ClinVar (n = 946 and 734). Positive favours th… |
| P68 | 19.9 | one of 53 cells | docs/lr-stability-check.md:3, phase1/src/phase8_lr_decomposition.py:5, scripts/weak_token_spotcheck.py:66 (+1) | …thresholds two point estimates cross it — Pangolin at 19.9 (achieved specificity… |
| P68 | 24.3 | one of 4 cells | phase1/src/phase8_lr_decomposition.py:14 | …nic calibration; interpolated 15.9) and AlphaGenome at 24.3 (97.2%; BRCA1 exclud… |
| P68 | 17.78 | one of 5 cells | docs/lr-stability-check.md:319 | …s interpolated LR+ in the primary condition runs 9.35, 17.78, 31.33 and 63.82, c… |
| P68 | 17.8 | one of 11 cells | README.md:46 | …hable (2.4). At 95% specificity the fusion attains LR+ 17.8 (95% CI 16.3–18.6) a… |
| P68 | 0.889 | one of 6 cells | docs/lr-stability-check.md:11 | …ttains LR+ 17.8 (95% CI 16.3–18.6) at a sensitivity of 0.889, the best single to… |
| P68 | 16.8 | one of 11 cells | README.md:46 | …–18.6) at a sensitivity of 0.889, the best single tool 16.8 (15.6–18.5), and ten… |
| P69 | 1.61 | one of 2 cells | README.md:46 | …ional standard — with BRCA1, ΔLR+ 1.02 (95% CI 0.04 to 1.61) on raw scores and 1… |
| P69 | 1.96 | single cell | phase1/config/manuscript_number_whitelist.json:15, phase1/src/phase2_model.py:128, phase1/src/phase2_model.py:154 (+5) | …(95% CI 0.04 to 1.61) on raw scores and 1.64 (0.07 to 1.96) after isotonic calib… |
| P69 | 19.9 | one of 21 cells | docs/lr-stability-check.md:3, phase1/src/phase8_lr_decomposition.py:5, scripts/weak_token_spotcheck.py:66 (+1) | …angolin's threshold lands at 96.0% specificity and its 19.9 crosses Strong while… |
| P70 | 1,675 | single cell | docs/lr-stability-check.md:114, docs/lr-stability-check.md:115, docs/lr-stability-check.md:118 (+6) | …ndard none does — but those are not the same variants: 1,675 carry an assay labe… |
| P70 | 946 | single cell | docs/lr-stability-check.md:210, docs/lr-stability-check.md:211, docs/lr-stability-check.md:214 (+3) | …are not the same variants: 1,675 carry an assay label, 946 a ClinVar label, and … |
| P70 | 907 | one of 19 cells | phase1/config/analysis_claims.json:124, phase1/src/phase7_label_contrast.py:12, phase1/src/phase7_label_contrast.py:13 (+1) | …: 1,675 carry an assay label, 946 a ClinVar label, and 907 carry both. Restricte… |
| P71 | 907 | one of 37 cells | phase1/config/analysis_claims.json:124, phase1/src/phase7_label_contrast.py:12, phase1/src/phase7_label_contrast.py:13 (+1) | …the two label sets is excluded by measurement: on the 907 shared variants they a… |
| P71 | 96.7 | one of 3 cells | phase1/config/analysis_claims.json:166 | …measurement: on the 907 shared variants they agree on 96.7% of calls (Cohen's κ … |
| P71 | 1.44 | computed derivation | phase1/config/analysis_claims.json:195, scripts/numbers_only_in_code.py:63 | …of VUS, conflicting and other records — by ratios from 1.44 to 4.20, and the fus… |
| P71 | 4.20 | computed derivation | docs/milestone3_analysis_ready.md:31, docs/milestone6_spliceai_pangolin.md:62, docs/milestone7_dnalm_alphagenome.md:69 (+1) | …conflicting and other records — by ratios from 1.44 to 4.20, and the fusion's se… |
| P71 | 0.971 | one of 13 cells | docs/lr-stability-check.md:73, docs/lr-stability-check.md:98 | …the fusion's sensitivity at 95% specificity falls from 0.971 to 0.681, CADD's fr… |
| P71 | 0.889 | one of 4 cells | docs/lr-stability-check.md:11 | …he subset sensitivities need not average to the pooled 0.889). The gap survives … |
| P71 | 20.2 | one of 6 cells | scripts/weak_token_spotcheck.py:176 | …led 0.889). The gap survives BRCA1 removal (fusion LR+ 20.2 against 13.7; CADD 1… |
| P71 | 13.7 | one of 6 cells | scripts/weak_token_spotcheck.py:177 | …he gap survives BRCA1 removal (fusion LR+ 20.2 against 13.7; CADD 19.1 against 9… |
| P71 | 19.1 | one of 17 cells | scripts/weak_token_spotcheck.py:178 | …ives BRCA1 removal (fusion LR+ 20.2 against 13.7; CADD 19.1 against 9.9). Two st… |
| P71 | 12.1 | one of 8 cells | scripts/weak_token_spotcheck.py:180 | …d within the assay-only arm the VUS subset (fusion LR+ 12.1) and the conflicting… |
| P73 | 74.9 | single cell | phase1/config/analysis_claims.json:133 | …Moderate or above on the functional standard, against 74.9% under the 0.90/0.10 … |
| P73 | 21.0 | single cell | environment.yml:17, requirements.txt:11 | …rp scores: in that same condition phastCons moves from 21.0% to 75.0%, because a… |
| P75 | 192 | one of 7 cells | README.md:14, docs/frozen-matrix-v2.md:33, docs/milestone4_candidate_genes.md:38 (+15) | …seven genes, then applied without refitting to TP53's 192 intron-side splice SNV… |
| P75 | 0.997 | single cell | phase1/config/analysis_claims.json:70, phase1/src/check_manuscript_numbers.py:341, phase1/src/check_manuscript_numbers.py:342 (+17) | …d gate used throughout (post-orientation control AUROC 0.997) (Supplementary Tab… |
| P77 | 0.968 | single cell | docs/lr-stability-check.md:87 | …and higher accuracy within the high-confidence subset (0.968 versus 0.940), whic… |
| P78 | 192 | one of 2 cells | README.md:14, docs/frozen-matrix-v2.md:33, docs/milestone4_candidate_genes.md:38 (+15) | …gene rather than a second significance claim: with n = 192 from one gene there i… |
| P85 | 0.218 | one of 2 cells | docs/milestone8b_results.md:38 | …d SpliceAI carry the weight in every fold (mean 0.224, 0.218 and 0.186, SD 0.019… |
| P85 | 0.186 | one of 2 cells | docs/milestone8_results.md:28, docs/milestone8b_results.md:22 | …carry the weight in every fold (mean 0.224, 0.218 and 0.186, SD 0.019, 0.029 and… |
| P85 | 0.0596 | one of 2 cells | phase1/config/analysis_claims.json:155 | …s the calibration only slightly: mean Brier moves from 0.0596 with six training … |
| T1r1c2 | 0.774 | single cell | README.md:44, docs/lr-stability-check.md:49 | …0.774… |
| T1r2c2 | 0.761 | single cell | README.md:44, docs/milestone8_results.md:21, docs/milestone8b_results.md:14 | …0.761… |
| T1r2c3 | 0.731 | one of 4 cells | docs/milestone8_results.md:24, docs/milestone8b_results.md:17 | …0.731 – 0.788… |
| T1r2c3 | 0.788 | one of 2 cells | docs/lr-stability-check.md:38 | …0.731 – 0.788… |
| T1r4c3 | 0.783 | one of 2 cells | docs/lr-stability-check.md:39 | …0.716 – 0.783… |
| T1r6c2 | 0.750 | single cell | docs/milestone8_results.md:23, docs/milestone8b_results.md:16 | …0.750… |
| T1r6c3 | 0.714 | one of 2 cells | docs/milestone8_results.md:23, docs/milestone8b_results.md:16 | …0.714 – 0.782… |
| T1r6c3 | 0.782 | one of 2 cells | docs/milestone8_results.md:23, docs/milestone8b_results.md:16 | …0.714 – 0.782… |
| T1r7c2 | 0.694 | single cell | docs/milestone8_results.md:24, docs/milestone8b_results.md:17 | …0.694… |
| T1r7c3 | 0.654 | one of 2 cells | docs/milestone8_results.md:24, docs/milestone8b_results.md:17 | …0.654 – 0.731… |
| T1r7c3 | 0.731 | one of 4 cells | docs/milestone8_results.md:24, docs/milestone8b_results.md:17 | …0.654 – 0.731… |
| T1r8c3 | 0.639 | one of 3 cells | docs/milestone8_results.md:25, docs/milestone8_results.md:26, docs/milestone8b_results.md:18 (+1) | …0.639 – 0.691… |
| T1r8c3 | 0.691 | one of 2 cells | docs/milestone8_results.md:25, docs/milestone8b_results.md:18 | …0.639 – 0.691… |
| T1r9c2 | 0.639 | one of 3 cells | docs/milestone8_results.md:25, docs/milestone8_results.md:26, docs/milestone8b_results.md:18 (+1) | …0.639… |
| T1r9c3 | 0.607 | one of 2 cells | docs/milestone8_results.md:26, docs/milestone8b_results.md:19 | …0.607 – 0.669… |
| T1r9c3 | 0.669 | one of 2 cells | docs/lr-stability-check.md:31, docs/milestone8_results.md:26, docs/milestone8b_results.md:19 | …0.607 – 0.669… |
| T1r10c2 | 0.612 | single cell | docs/milestone8_results.md:27, docs/milestone8b_results.md:20 | …0.612… |
| T1r10c3 | 0.563 | one of 2 cells | docs/milestone8_results.md:27, docs/milestone8b_results.md:20 | …0.563 – 0.657… |
| T1r10c3 | 0.657 | one of 2 cells | docs/milestone8_results.md:27, docs/milestone8b_results.md:20 | …0.563 – 0.657… |
| T1r12c2 | 0.186 | single cell | docs/milestone8_results.md:28, docs/milestone8b_results.md:22 | …0.186… |
| T1r12c3 | 0.298 | one of 2 cells | docs/milestone8_results.md:28, docs/milestone8b_results.md:22 | …0.069 – 0.298… |
| T2r1c2 | 0.774 | one of 2 cells | README.md:44, docs/lr-stability-check.md:49 | …0.774… |
| T2r1c3 | 0.773 | single cell | docs/lr-stability-check.md:19 | …0.773… |
| T2r2c2 | 0.774 | one of 2 cells | README.md:44, docs/lr-stability-check.md:49 | …0.774… |
| T2r2c3 | 0.771 | single cell | docs/lr-stability-check.md:50 | …0.771… |
| T2r3c3 | 0.759 | single cell | docs/lr-stability-check.md:42, docs/lr-stability-check.md:311, scripts/weak_token_spotcheck.py:247 | …0.759… |
| T3r1c4 | 74.9 | single cell | phase1/config/analysis_claims.json:133 | …74.9… |
| T3r1c5 | 96.7 | single cell | phase1/config/analysis_claims.json:166 | …96.7… |
| T3r9c4 | 21.0 | single cell | environment.yml:17, requirements.txt:11 | …21.0… |
| T4r1c1 | 192 | one of 2 cells | README.md:14, docs/frozen-matrix-v2.md:33, docs/milestone4_candidate_genes.md:38 (+15) | …192… |
| T4r2c1 | 192 | one of 2 cells | README.md:14, docs/frozen-matrix-v2.md:33, docs/milestone4_candidate_genes.md:38 (+15) | …192… |
| T4r3c1 | 181 | single cell | docs/lr-stability-check.md:206, docs/lr-stability-check.md:207 | …181… |
