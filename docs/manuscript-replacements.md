# Manuscript number replacements: `calibration_draft.docx`

v1 outputs `/Users/cliffzhang/work/variant-fm-benchmark/phase1/reports/phase1_v1` → v2 outputs `/Users/cliffzhang/work/variant-fm-benchmark/phase1/reports/phase1`. Tokens: 506; by status: context 59, decided 45, kept_integer 31, manual 7, no_v1_source 19, no_v2_counterpart 1, out_of_scope 12, unambiguous 242, whitelisted 90. Tokens whose printed value changes: 174.

## Replacements (position, old, new)

| block | old | new | status | cell(s) (v1→v2) | context |
|---|---|---|---|---|---|
| P8 | +0.012 | +0.013 | decided | phase2_H1_stratified.csv[0].delta (0.0117→0.0129) | es Spearman ρ 0.76 and an elastic-net fusion adds only +0.012, not a robust advantage — saturati |
| P8 | 75.1 | 74.9 | context | phase3_calibration_summary.csv[0].actionable_frac (75.1→74.87); phase5b_evidence_yield.csv[0].yield_confidence_0.90_0.10 (75.1045→74.8657) | ervals excluding zero, a tie in the primary condition (75.1% vs 75.2%), with no accuracy diffe |
| P8 | 75.2 | 73.6 | unambiguous | phase3_calibration_summary.csv[3].actionable_frac (75.22→73.55); phase5b_evidence_yield.csv[3].yield_confidence_0.90_0.10 (75.2239→73.5522) | cluding zero, a tie in the primary condition (75.1% vs 75.2%), with no accuracy difference in |
| P9 | 31.75 | 31.85 | unambiguous | phase8_lr_plus_operating_points.csv[2].lr_plus (31.7507→31.8498) | sion and several single tools reach Strong (fusion LR+ 31.75 at 97.5%). Predictors do better on |
| P29 | 21.1 | 19.6 | unambiguous | phase5_likelihood_ratios.csv[79].lr_plus (21.0588→19.6416); phase5_lr_headline.csv[7].lr_plus_best_single (21.0588→19.6416) | wo largest ratios on the functional standard (Pangolin 21.1 at 96.5% realised, AlphaGenome 24. |
| P30 | +0.0095 | +0.0104 | decided | ipw_ranking.csv[11].rho (0.00946136→0.0103785) | s quintile, decile and twentile binning, Δρ moves from +0.0095 to +0.0098 [+0.0006, +0.0205] and |
| P30 | +0.0098 | +0.0105 | decided | ipw_ranking.csv[35].rho (0.00984553→0.010486) | decile and twentile binning, Δρ moves from +0.0095 to +0.0098 [+0.0006, +0.0205] and ΔBrier from |
| P30 | +0.0006 | −0.0005 | unambiguous | ipw_ranking.csv[35].ci_lo (0.000603014→-0.000468309); ipw_ranking.csv[59].ci_lo (0.000603014→-0.000468309) | nd twentile binning, Δρ moves from +0.0095 to +0.0098 [+0.0006, +0.0205] and ΔBrier from 0.0082 t |
| P30 | +0.0205 | +0.0226 | unambiguous | ipw_ranking.csv[35].ci_hi (0.0205399→0.022601); ipw_ranking.csv[59].ci_hi (0.0205399→0.022601) | le binning, Δρ moves from +0.0095 to +0.0098 [+0.0006, +0.0205] and ΔBrier from 0.0082 to 0.0083 |
| P30 | 0.0082 | 0.0098 | context | ipw_headline.csv[0].dBrier (0.00823781→0.00981805) | +0.0095 to +0.0098 [+0.0006, +0.0205] and ΔBrier from 0.0082 to 0.0083 [+0.0020, +0.0178], sign |
| P30 | 0.0083 | 0.0099 | unambiguous | ipw_headline.csv[2].dBrier (0.00828659→0.00985881); ipw_headline.csv[4].dBrier (0.00828659→0.00985881) | o +0.0098 [+0.0006, +0.0205] and ΔBrier from 0.0082 to 0.0083 [+0.0020, +0.0178], significant in |
| P30 | +0.0020 | +0.0031 | unambiguous | ipw_headline.csv[2].dBrier_lo (0.00195041→0.00308802); ipw_headline.csv[4].dBrier_lo (0.00195041→0.00308802) | 8 [+0.0006, +0.0205] and ΔBrier from 0.0082 to 0.0083 [+0.0020, +0.0178], significant in the same |
| P30 | +0.0178 | +0.0194 | unambiguous | ipw_headline.csv[2].dBrier_hi (0.0177781→0.019403); ipw_headline.csv[4].dBrier_hi (0.0177781→0.019403) | 6, +0.0205] and ΔBrier from 0.0082 to 0.0083 [+0.0020, +0.0178], significant in the same three of |
| P38 | 0.753 | 0.752 | unambiguous | phase2_leaderboard.csv[3].pooled_rho (0.753→0.7517); phase8_leaderboard_hk.csv[3].pooled_rho (0.752989→0.751658) | e 1): pooled per-gene Spearman ρ was 0.761 (Pangolin), 0.753 (SpliceAI) and 0.750 (AlphaGenome) |
| P38 | 0.494 | 0.492 | unambiguous | phase2_leaderboard.csv[10].pooled_rho (0.4942→0.4923); phase8_leaderboard_hk.csv[10].pooled_rho (0.494233→0.492287) | oP (0.639), phastCons (0.612), Nucleotide Transformer (0.494) and gnomAD AF (0.186). AlphaMisse |
| P39 | 0.773 | 0.774 | context | phase2_H1_stratified.csv[0].fusion_rho (0.773→0.7738); phase2_leaderboard.csv[0].pooled_rho (0.773→0.7738); phase8_leaderboard_hk.csv[0].pooled_rho (0.773→0.773833) | The elastic-net fusion reached ρ = 0.773, exceeding the best single tool (P |
| P39 | +0.012 | +0.013 | decided | phase2_H1_stratified.csv[0].delta (0.0117→0.0129) | 773, exceeding the best single tool (Pangolin) by Δρ = +0.012 (95% CI [+0.002, +0.021]) — a posi |
| P39 | +0.021 | +0.024 | unambiguous | phase2_H1_stratified.csv[0].ci95#2 (0.0214→0.0239); phase8_sign_flip_exact.csv[1].T_obs (0.0207941→0.0242656) | single tool (Pangolin) by Δρ = +0.012 (95% CI [+0.002, +0.021]) — a positive but marginal margin |
| P39 | +0.0113 | +0.0138 | unambiguous | phase8_leave_two_genes_out.csv[12].delta_rho (0.0113312→0.0137831) | all 21 gene pairs keep Δρ positive throughout — median +0.0113, from +0.0061 (dropping BARD1 and |
| P39 | +0.0061 | +0.0064 | unambiguous | phase8_leave_two_genes_out.csv[6].delta_rho (0.00608977→0.00636961) | irs keep Δρ positive throughout — median +0.0113, from +0.0061 (dropping BARD1 and BRCA1) to +0.0 |
| P39 | +0.0187 | +0.0207 | unambiguous | phase8_leave_two_genes_out.csv[2].delta_rho (0.0186941→0.0207216) | an +0.0113, from +0.0061 (dropping BARD1 and BRCA1) to +0.0187 (dropping BAP1 and BRCA2) — so the |
| P39 | +0.0095 | +0.0104 | decided | phase2_no_offset_drop.csv[1].value (0.00946136→0.0103785) | disproportionate share, and dropping them moves Δρ to +0.0095, below all two hundred random thir |
| P39 | +0.0117 | +0.0129 | unambiguous | phase2_H1_stratified.csv[0].delta (0.0117→0.0129); phase2_no_offset_drop.csv[0].value (0.0117366→0.0128749); phase2_no_offset_drop.csv[5].value (0.0117231→0.0128624) | , below all two hundred random thirteen-variant drops (+0.0117 ± 0.0004). The interval still excl |
| P39 | +0.0005 | −0.0003 | unambiguous | phase2_no_offset_drop.csv[2].value (0.000496937→-0.000253412) | 2,000 resamples but its lower bound falls 4.6-fold (to +0.0005; the per-draw values are archived |
| P39 | +0.0372 | +0.0337 | unambiguous | phase2_H1_stratified.csv[1].delta (0.0372→0.0337) | point estimates the gain is larger in the splice core (+0.0372) than in the splice-region band (+ |
| P39 | −0.034 | −0.035 | unambiguous | phase2_H1_stratified.csv[1].ci95#1 (-0.0336→-0.0345) | (+0.0217), but neither interval excludes zero (95% CI [−0.034, +0.119] and [−0.0001, +0.045]) an |
| P39 | +0.119 | +0.118 | unambiguous | phase2_H1_stratified.csv[1].ci95#2 (0.1194→0.1178) | ), but neither interval excludes zero (95% CI [−0.034, +0.119] and [−0.0001, +0.045]) and absolu |
| P39 | −0.0001 | −0.0056 | unambiguous | phase2_H1_stratified.csv[2].ci95#1 (-0.0001→-0.0056) | r interval excludes zero (95% CI [−0.034, +0.119] and [−0.0001, +0.045]) and absolute performance |
| P39 | +0.045 | +0.051 | unambiguous | phase2_H1_stratified.csv[2].ci95#2 (0.0446→0.0506) | l excludes zero (95% CI [−0.034, +0.119] and [−0.0001, +0.045]) and absolute performance in the |
| P39 | 0.198 | 0.203 | unambiguous | phase2_H1_stratified.csv[1].single_rho (0.198→0.2032) | 0.045]) and absolute performance in the core is low (ρ 0.198–0.235), so the larger gain there i |
| P39 | 0.235 | 0.237 | unambiguous | phase2_H1_stratified.csv[1].fusion_rho (0.2353→0.2369) | ) and absolute performance in the core is low (ρ 0.198–0.235), so the larger gain there is meas |
| P43 | 0.740 | 0.741 | unambiguous | phase2_leaderboard.csv[0].ci95#1 (0.74→0.741); phase8_leaderboard_hk.csv[0].dl_lo (0.73975→0.740855) | ensitivity analysis widens the fusion's interval from [0.740, 0.802] to [0.729, 0.810] without |
| P43 | 0.802 | 0.803 | unambiguous | phase2_leaderboard.csv[0].ci95#2 (0.802→0.803); phase8_leaderboard_hk.csv[0].dl_hi (0.802485→0.80309) | ity analysis widens the fusion's interval from [0.740, 0.802] to [0.729, 0.810] without changin |
| P43 | 0.729 | 0.730 | context | phase8_leaderboard_hk.csv[0].hk_lo (0.729227→0.73037) | s widens the fusion's interval from [0.740, 0.802] to [0.729, 0.810] without changing the order |
| P43 | 0.810 | 0.811 | unambiguous | phase8_leaderboard_hk.csv[0].hk_hi (0.810472→0.811056) | s the fusion's interval from [0.740, 0.802] to [0.729, 0.810] without changing the ordering (fu |
| P45 | −0.001 | 0.001 | unambiguous | phase2_H2_ablation.csv[0].delta_full_minus_ablated (-0.0009→0.0007); phase2_H2_ablation.csv[1].ci95#1 (-0.0014→0.0005) | d not degrade ranking: dropping conservation gave Δρ = −0.001 (95% CI [−0.003, +0.001]) and drop |
| P45 | −0.003 | −0.001 | unambiguous | phase2_H2_ablation.csv[0].ci95#1 (-0.0025→-0.001) | nking: dropping conservation gave Δρ = −0.001 (95% CI [−0.003, +0.001]) and dropping conservatio |
| P45 | +0.001 | +0.002 | decided | phase2_H2_ablation.csv[0].ci95#2 (0.0009→0.0021) | ropping conservation gave Δρ = −0.001 (95% CI [−0.003, +0.001]) and dropping conservation plus t |
| P45 | +0.002 | +0.003 | decided | phase2_H2_ablation.csv[1].delta_full_minus_ablated (0.0022→0.0032) | vation plus the alignment-conditioned signal gave Δρ = +0.002 ([−0.001, +0.007]); the gradient-b |
| P45 | +0.007 | +0.006 | decided | phase2_H2_ablation.csv[1].ci95#2 (0.0066→0.0064) | lignment-conditioned signal gave Δρ = +0.002 ([−0.001, +0.007]); the gradient-boosted check agre |
| P47 | 0.111 | 0.109 | unambiguous | phase6_enet_weight_stability.csv[3].mean (0.110624→0.109354); phase6_enet_weight_stability.csv[3].mean_abs (0.110624→0.109354) | Elastic-net /coefficients/ in the full model: GPN-MSA 0.111 (rank 4 of 20), phastCons 0.051 (5 |
| P47 | 0.051 | 0.050 | context | phase6_enet_weight_stability.csv[4].mean (0.0508968→0.0502411); phase6_enet_weight_stability.csv[4].mean_abs (0.0508968→0.0502411) | he full model: GPN-MSA 0.111 (rank 4 of 20), phastCons 0.051 (5), phyloP 0.017 (7). |
| P47 | 0.017 | 0.018 | decided | phase6_enet_weight_stability.csv[6].mean (0.0172281→0.018472) | -MSA 0.111 (rank 4 of 20), phastCons 0.051 (5), phyloP 0.017 (7). |
| P50 | +0.0085 | +0.0103 | decided | phase3_H3_headline.csv[0].dBrier (0.0085→0.0103) | ary condition, the functional standard with BRCA1 (Δ = +0.0085, 95% CI [+0.0032, +0.0165]); six o |
| P50 | +0.0032 | +0.0047 | context | phase3_H3_headline.csv[0].dBrier_ci_geneclust#1 (0.0032→0.0047) | e functional standard with BRCA1 (Δ = +0.0085, 95% CI [+0.0032, +0.0165]); six of the seven genes |
| P50 | +0.0165 | +0.0182 | decided | phase3_H3_headline.csv[0].dBrier_ci_geneclust#2 (0.0165→0.0182) | nal standard with BRCA1 (Δ = +0.0085, 95% CI [+0.0032, +0.0165]); six of the seven genes favour t |
| P50 | +0.0116 | +0.0117 | decided | phase3_H3_headline.csv[1].dBrier (0.0116→0.0117) | the interval likewise excluded zero without BRCA1 (Δ = +0.0116, [+0.0040, +0.0185]) and on ClinVa |
| P50 | +0.0040 | +0.0042 | unambiguous | phase3_H3_headline.csv[1].dBrier_ci_geneclust#1 (0.004→0.0042) | al likewise excluded zero without BRCA1 (Δ = +0.0116, [+0.0040, +0.0185]) and on ClinVar (+0.0119 |
| P50 | +0.0185 | +0.0187 | unambiguous | phase3_H3_headline.csv[1].dBrier_ci_geneclust#2 (0.0185→0.0187) | se excluded zero without BRCA1 (Δ = +0.0116, [+0.0040, +0.0185]) and on ClinVar (+0.0119, [+0.000 |
| P50 | +0.0119 | +0.0124 | decided | phase3_H3_headline.csv[2].dBrier (0.0119→0.0124) | RCA1 (Δ = +0.0116, [+0.0040, +0.0185]) and on ClinVar (+0.0119, [+0.0004, +0.0250]; +0.0192, [+0. |
| P50 | +0.0004 | +0.0009 | decided | phase3_H3_headline.csv[2].dBrier_ci_geneclust#1 (0.0004→0.0009) | +0.0116, [+0.0040, +0.0185]) and on ClinVar (+0.0119, [+0.0004, +0.0250]; +0.0192, [+0.0028, +0.0 |
| P50 | +0.0250 | +0.0257 | decided | phase3_H3_headline.csv[2].dBrier_ci_geneclust#2 (0.025→0.0257) | [+0.0040, +0.0185]) and on ClinVar (+0.0119, [+0.0004, +0.0250]; +0.0192, [+0.0028, +0.0361]). Th |
| P50 | +0.0192 | +0.0196 | unambiguous | phase3_H3_headline.csv[3].dBrier (0.0192→0.0196); phase3_pertool_brier_ci.csv[31].dBrier_vs_fusion (0.0192→0.0196); phase8_murphy_decomposition.csv[9].d_brier (0.0192207→0.0196039) | +0.0185]) and on ClinVar (+0.0119, [+0.0004, +0.0250]; +0.0192, [+0.0028, +0.0361]). The fusion h |
| P50 | +0.0028 | +0.0042 | context | phase8_murphy_decomposition.csv[6].rel_fusion (0.00284893→0.00418837) | and on ClinVar (+0.0119, [+0.0004, +0.0250]; +0.0192, [+0.0028, +0.0361]). The fusion had the low |
| P50 | +0.0361 | +0.0349 | decided | phase3_H3_headline.csv[3].dBrier_ci_geneclust#2 (0.0361→0.0349) | inVar (+0.0119, [+0.0004, +0.0250]; +0.0192, [+0.0028, +0.0361]). The fusion had the lowest Brier |
| P50 | 0.0469 | 0.0312 | unambiguous | phase8_sign_flip_exact.csv[2].p_exact (0.046875→0.03125) | on the primary-condition ΔBrier gives a two-sided p = 0.0469 — significant at 0.05 but not at 0 |
| P50 | 72 | 74 | context | phase8_lr_plus_operating_points.csv[51].threshold (72.4745→73.7411) | e S16); the interval excluded zero in every one of the 72 fusion-versus-single-tool comparis |
| P50 | −0.0001 | 0.0011 | context | ipw_headline.csv[10].dBrier_lo (-0.000111195→0.00113042) | 000 gene resamples the ClinVar +BRCA1 interval reaches −0.0001 in the unweighted baseline and −0. |
| P50 | 75.1 | 74.9 | context | phase3_calibration_summary.csv[0].actionable_frac (75.1→74.87); phase5b_evidence_yield.csv[0].yield_confidence_0.90_0.10 (75.1045→74.8657) | dard including BRCA1 — the two were indistinguishable: 75.1% high-confidence at 96.8% accuracy |
| P50 | 96.8 | 96.7 | decided | phase3_calibration_summary.csv[0].actionable_acc (96.82→96.73) | e two were indistinguishable: 75.1% high-confidence at 96.8% accuracy versus Pangolin's 75.2% |
| P50 | 75.2 | 73.6 | context | phase3_calibration_summary.csv[3].actionable_frac (75.22→73.55); phase5b_evidence_yield.csv[3].yield_confidence_0.90_0.10 (75.2239→73.5522) | 1% high-confidence at 96.8% accuracy versus Pangolin's 75.2% at 96.3%, a difference of two var |
| P51 | 0.0631 | 0.0630 | unambiguous | phase3_calibration_summary.csv[0].Brier (0.0631→0.063); phase8_murphy_decomposition.csv[0].brier_fusion (0.0631095→0.0630194); phase8_murphy_decomposition.csv[1].brier_fusion (0.0631095→0.0630194) | reliability, resolution and uncertainty. The fusion's 0.0631 against Pangolin's 0.0717 — the he |
| P51 | 0.0717 | 0.0733 | context | phase3_calibration_summary.csv[3].Brier (0.0717→0.0733) | nd uncertainty. The fusion's 0.0631 against Pangolin's 0.0717 — the headline ΔBrier is +0.0085 [ |
| P51 | +0.0085 | +0.0103 | context | phase3_pertool_brier_ci.csv[1].dBrier_vs_fusion (0.0085→0.0103) | 631 against Pangolin's 0.0717 — the headline ΔBrier is +0.0085 [+0.0032, +0.0165] — divides as fo |
| P51 | +0.0032 | +0.0047 | context | phase3_H3_headline.csv[0].dBrier_ci_geneclust#1 (0.0032→0.0047) | st Pangolin's 0.0717 — the headline ΔBrier is +0.0085 [+0.0032, +0.0165] — divides as follows: th |
| P51 | +0.0165 | +0.0182 | decided | phase3_H3_headline.csv[0].dBrier_ci_geneclust#2 (0.0165→0.0182) | in's 0.0717 — the headline ΔBrier is +0.0085 [+0.0032, +0.0165] — divides as follows: the uncerta |
| P51 | 0.0047 | 0.0037 | decided | phase8_murphy_decomposition.csv[0].rel_fusion (0.00468447→0.00369326) | (0.2499); the reliability terms are indistinguishable (0.0047 versus 0.0048, ΔREL +0.0001 [−0.00 |
| P51 | 0.0048 | 0.0044 | context | phase8_murphy_decomposition.csv[0].rel_single (0.00481779→0.00442177) | reliability terms are indistinguishable (0.0047 versus 0.0048, ΔREL +0.0001 [−0.0008, +0.0067]); |
| P51 | +0.0001 | +0.0007 | decided | phase8_murphy_decomposition.csv[0].d_rel (0.000133317→0.000728508) | erms are indistinguishable (0.0047 versus 0.0048, ΔREL +0.0001 [−0.0008, +0.0067]); and essential |
| P51 | −0.0008 | −0.0011 | unambiguous | phase8_murphy_decomposition.csv[0].d_rel_lo (-0.000822836→-0.00114282) | indistinguishable (0.0047 versus 0.0048, ΔREL +0.0001 [−0.0008, +0.0067]); and essentially the wh |
| P51 | +0.0067 | +0.0141 | decided | phase8_murphy_decomposition.csv[0].d_rel_hi (0.0067169→0.0140778) | uishable (0.0047 versus 0.0048, ΔREL +0.0001 [−0.0008, +0.0067]); and essentially the whole gain |
| P51 | 0.1911 | 0.1899 | unambiguous | phase8_murphy_decomposition.csv[0].res_fusion (0.191122→0.189885) | ]); and essentially the whole gain sits in resolution (0.1911 versus 0.1822, ΔRES +0.0090 [+0.00 |
| P51 | 0.1822 | 0.1807 | unambiguous | phase8_murphy_decomposition.csv[0].res_single (0.182167→0.18075) | ially the whole gain sits in resolution (0.1911 versus 0.1822, ΔRES +0.0090 [+0.0006, +0.0128]). |
| P51 | +0.0090 | +0.0091 | unambiguous | phase8_murphy_decomposition.csv[0].d_res (0.00895552→0.009135) | le gain sits in resolution (0.1911 versus 0.1822, ΔRES +0.0090 [+0.0006, +0.0128]). Fusion theref |
| P51 | +0.0006 | +0.0002 | decided | phase8_murphy_decomposition.csv[0].d_res_lo (0.000550274→0.000151505) | its in resolution (0.1911 versus 0.1822, ΔRES +0.0090 [+0.0006, +0.0128]). Fusion therefore buys |
| P51 | +0.0128 | +0.0121 | decided | phase8_murphy_decomposition.csv[0].d_res_hi (0.0127932→0.0121457) | solution (0.1911 versus 0.1822, ΔRES +0.0090 [+0.0006, +0.0128]). Fusion therefore buys a better |
| P51 | 0.0085 | 0.0103 | decided | phase3_H3_headline.csv[0].dBrier (0.0085→0.0103) | e binned terms sum to ≈0.009; the unbinned headline is 0.0085 — the binned decomposition approxi |
| P52 | +0.0065 | +0.0052 | context | phase3_H3_headline.csv[2].dECE (0.0065→0.0052) | ntage on ClinVar carried intervals excluding zero (Δ = +0.0065, [+0.0002, +0.0247] with BRCA1; +0 |
| P52 | +0.0002 | +0.0008 | decided | phase3_H3_headline.csv[2].dECE_ci#1 (0.0002→0.0008) | linVar carried intervals excluding zero (Δ = +0.0065, [+0.0002, +0.0247] with BRCA1; +0.0157, [+0 |
| P52 | +0.0247 | +0.0215 | context | phase3_H3_headline.csv[2].dECE_ci#2 (0.0247→0.0215) | rried intervals excluding zero (Δ = +0.0065, [+0.0002, +0.0247] with BRCA1; +0.0157, [+0.0063, +0 |
| P52 | +0.0157 | +0.0177 | context | phase3_pertool_brier_ci.csv[1].ci95#2 (0.0157→0.0177) | ding zero (Δ = +0.0065, [+0.0002, +0.0247] with BRCA1; +0.0157, [+0.0063, +0.0396] without), and |
| P52 | +0.0063 | +0.0042 | context | phase5b_evidence_yield.csv[24].delta (0.00634249→0.00422833) | (Δ = +0.0065, [+0.0002, +0.0247] with BRCA1; +0.0157, [+0.0063, +0.0396] without), and did not on |
| P52 | +0.0396 | +0.0451 | context | phase3_H3_headline.csv[2].dYield_ci#2 (0.0396→0.0451) | 065, [+0.0002, +0.0247] with BRCA1; +0.0157, [+0.0063, +0.0396] without), and did not on the func |
| P54 | 0.026 | 0.025 | decided | phase3_calibration_summary.csv[0].ECE (0.026→0.0249) | area; in the primary condition the fusion attains ECE 0.026 and Brier 0.0631 (Pangolin 0.033 a |
| P54 | 0.0631 | 0.0630 | unambiguous | phase3_calibration_summary.csv[0].Brier (0.0631→0.063); phase8_murphy_decomposition.csv[0].brier_fusion (0.0631095→0.0630194); phase8_murphy_decomposition.csv[1].brier_fusion (0.0631095→0.0630194) | imary condition the fusion attains ECE 0.026 and Brier 0.0631 (Pangolin 0.033 and 0.0717; Table |
| P54 | 0.0717 | 0.0733 | context | phase3_calibration_summary.csv[3].Brier (0.0717→0.0733) | attains ECE 0.026 and Brier 0.0631 (Pangolin 0.033 and 0.0717; Table 3). The two extreme bins ho |
| P54 | 1,258 | 1,254 | decided | phase3_calibration_summary.csv[0].n_actionable (1258→1254) | 0.033 and 0.0717; Table 3). The two extreme bins hold 1,258 of 1,675 variants and track the di |
| P54 | 0.974 | 0.975 | decided | phase3_reliability_fusion.csv[9].mean_pred (0.974175→0.975369) | gonal closely (0.026 predicted against 0.026 observed; 0.974 against 0.962). |
| P56 | 1,168 | 1,173 | context | phase3_calibration_summary.csv[1].n_actionable (1168→1173) | across the functional standard (n = 1,675 with BRCA1, 1,168 without) and ClinVar (n = 946 and |
| P59 | 21.1 | 19.6 | context | phase5_likelihood_ratios.csv[79].lr_plus (21.0588→19.6416); phase5_lr_headline.csv[7].lr_plus_best_single (21.0588→19.6416) | that cross the 18.7 threshold elsewhere — Pangolin at 21.1 and AlphaGenome at 24.3, both with |
| P59 | 31.75 | 31.85 | unambiguous | phase8_lr_plus_operating_points.csv[2].lr_plus (31.7507→31.8498) | usion's LR+ in the primary condition runs 9.40, 18.08, 31.75 and 67.48, crossing Strong from 97 |
| P59 | 67.48 | 68.13 | context | phase8_lr_plus_operating_points.csv[3].lr_plus (67.4826→68.1328) | + in the primary condition runs 9.40, 18.08, 31.75 and 67.48, crossing Strong from 97.5% up — a |
| P59 | 30.26 | 30.22 | unambiguous | phase8_lr_plus_operating_points.csv[14].lr_plus (30.2647→30.2152) | he other ten evaluable objects cross as well (Pangolin 30.26, SpliceAI 29.28, the equal-weight |
| P59 | 29.28 | 28.73 | unambiguous | phase8_lr_plus_operating_points.csv[10].lr_plus (29.2815→28.7292) | luable objects cross as well (Pangolin 30.26, SpliceAI 29.28, the equal-weight mean 29.22, Alph |
| P60 | 1.02 | 1.04 | context | phase5_lr_headline.csv[0].delta_lr_plus (1.01543→1.04019); phase5_lr_headline.csv[8].delta_lr_plus (1.01543→1.04019) | fter calibration (functional standard with BRCA1: ΔLR+ 1.02, 95% CI −0.73 to 1.62), and in sev |
| P60 | −0.73 | −0.36 | context | phase5_lr_headline.csv[0].lo (-0.730301→-0.358564) | ion (functional standard with BRCA1: ΔLR+ 1.02, 95% CI −0.73 to 1.62), and in seven of the eigh |
| P60 | 0.93 | 0.78 | context | phase5_calibration_effect.csv[0].median_abs_change_lr_plus (0.930527→0.777965) | introduces: the median change in LR+ on calibrating is 0.93, and three of 48 model × condition |
| P64 | 86.8 | 87.0 | unambiguous | phase5b_evidence_yield.csv[0].yield_acmg_moderate_or_above (86.806→86.9851) | -based definition of clinical yield, the fusion places 86.8% of variants at Moderate or above |
| P66 | 0.997 | 0.998 | unambiguous | directionality_check.csv[3].control_auroc (0.9975→0.9975) | d gate used throughout (post-orientation control AUROC 0.997) (Supplementary Table S6). |
| P68 | 0.079 | 0.077 | unambiguous | phase4_tp53_external.csv[0].Brier (0.0785→0.0768) | TP53 was lower than the best single tool's (SpliceAI): 0.079 versus 0.122 (ΔBrier 0.043; varian |
| P68 | 0.122 | 0.118 | context | phase4_tp53_external.csv[1].Brier (0.1216→0.1178) | r than the best single tool's (SpliceAI): 0.079 versus 0.122 (ΔBrier 0.043; variant-level 95% C |
| P68 | 0.043 | 0.041 | unambiguous | phase4_tp53_label_definitions.csv[0].dBrier (0.0430949→0.0410328) | t single tool's (SpliceAI): 0.079 versus 0.122 (ΔBrier 0.043; variant-level 95% CI 0.016–0.071) |
| P68 | 0.016 | 0.025 | decided | phase4_tp53_label_definitions.csv[0].ci_lo (0.0155127→0.0245144) | 0.079 versus 0.122 (ΔBrier 0.043; variant-level 95% CI 0.016–0.071), with lower calibration err |
| P68 | 0.071 | 0.059 | unambiguous | phase4_tp53_label_definitions.csv[0].ci_hi (0.0707184→0.0589604) | versus 0.122 (ΔBrier 0.043; variant-level 95% CI 0.016–0.071), with lower calibration error (EC |
| P68 | 0.120 | 0.122 | unambiguous | phase4_tp53_external.csv[0].ECE (0.12→0.1225) | 95% CI 0.016–0.071), with lower calibration error (ECE 0.120 versus 0.154) and higher accuracy |
| P68 | 0.154 | 0.157 | unambiguous | phase4_tp53_external.csv[1].ECE (0.1539→0.1569) | 0.071), with lower calibration error (ECE 0.120 versus 0.154) and higher accuracy within the hi |
| P68 | 0.978 | 0.968 | unambiguous | phase4_tp53_external.csv[0].actionable_acc (0.9781→0.9677) | and higher accuracy within the high-confidence subset (0.978 versus 0.912). The control-anchore |
| P68 | 0.912 | 0.940 | unambiguous | phase4_tp53_external.csv[1].actionable_acc (0.9118→0.9404) | curacy within the high-confidence subset (0.978 versus 0.912). The control-anchored, median-spl |
| P76 | 0.222 | 0.224 | context | phase6_enet_weight_stability.csv[0].mean (0.222161→0.224002); phase6_enet_weight_stability.csv[0].mean_abs (0.222161→0.224002) | olin and SpliceAI carry the weight in every fold (mean 0.222, 0.210 and 0.198, SD 0.019, 0.029 |
| P76 | 0.210 | 0.218 | unambiguous | phase6_enet_weight_stability.csv[1].mean (0.209755→0.218342); phase6_enet_weight_stability.csv[1].mean_abs (0.209755→0.218342) | d SpliceAI carry the weight in every fold (mean 0.222, 0.210 and 0.198, SD 0.019, 0.029 and 0.0 |
| P76 | 0.198 | 0.186 | context | phase6_enet_weight_stability.csv[2].mean (0.198191→0.185773); phase6_enet_weight_stability.csv[2].mean_abs (0.198191→0.185773) | carry the weight in every fold (mean 0.222, 0.210 and 0.198, SD 0.019, 0.029 and 0.022), with |
| P76 | 0.022 | 0.021 | decided | phase6_enet_weight_stability.csv[2].sd (0.0217264→0.0206481) | fold (mean 0.222, 0.210 and 0.198, SD 0.019, 0.029 and 0.022), with no sign change in any featu |
| P76 | 0.0589 | 0.0596 | context | phase6_training_gene_summary.csv[3].mean (0.0589067→0.0596437); phase6_training_gene_summary.csv[3].mean_paired_over_held_out_genes (0.0589067→0.0596437) | s the calibration only slightly: mean Brier moves from 0.0589 with six training genes to 0.0608 |
| P76 | 0.0608 | 0.0609 | unambiguous | phase6_training_gene_summary.csv[0].mean (0.0607935→0.0608979); phase6_training_gene_summary.csv[0].mean_paired_over_held_out_genes (0.0607935→0.0608979) | ean Brier moves from 0.0589 with six training genes to 0.0608 with three, worse in all seven hel |
| T1r1c2 | 0.773 | 0.774 | unambiguous | phase2_leaderboard.csv[0].pooled_rho (0.773→0.7738); phase8_leaderboard_hk.csv[0].pooled_rho (0.773→0.773833) | 0.773 |
| T1r1c3 | 0.740 | 0.741 | unambiguous | phase2_leaderboard.csv[0].ci95#1 (0.74→0.741); phase8_leaderboard_hk.csv[0].dl_lo (0.73975→0.740855) | 0.740 – 0.802 |
| T1r1c3 | 0.802 | 0.803 | unambiguous | phase2_leaderboard.csv[0].ci95#2 (0.802→0.803); phase8_leaderboard_hk.csv[0].dl_hi (0.802485→0.80309) | 0.740 – 0.802 |
| T1r2c3 | 0.729 | 0.731 | unambiguous | phase2_leaderboard.csv[1].ci95#1 (0.729→0.731); phase8_leaderboard_hk.csv[1].dl_lo (0.729192→0.731216) | 0.729 – 0.790 |
| T1r2c3 | 0.790 | 0.788 | unambiguous | phase2_leaderboard.csv[1].ci95#2 (0.79→0.788); phase8_leaderboard_hk.csv[1].dl_hi (0.789998→0.787813) | 0.729 – 0.790 |
| T1r2c4 | 54 | 47 | decided | phase8_leaderboard_hk.csv[1].I2 (53.5236→46.5256) | 54 |
| T1r3c2 | 0.759 | 0.758 | unambiguous | phase2_leaderboard.csv[2].pooled_rho (0.7589→0.7577); phase8_leaderboard_hk.csv[2].pooled_rho (0.758887→0.757702) | 0.759 |
| T1r3c3 | 0.738 | 0.737 | unambiguous | phase2_leaderboard.csv[2].ci95#1 (0.738→0.737); phase8_leaderboard_hk.csv[2].dl_lo (0.73836→0.73709) | 0.738 – 0.778 |
| T1r3c3 | 0.778 | 0.777 | unambiguous | phase2_leaderboard.csv[2].ci95#2 (0.778→0.777); phase8_leaderboard_hk.csv[2].dl_hi (0.778008→0.776906) | 0.738 – 0.778 |
| T1r4c2 | 0.753 | 0.752 | unambiguous | phase2_leaderboard.csv[3].pooled_rho (0.753→0.7517); phase8_leaderboard_hk.csv[3].pooled_rho (0.752989→0.751658) | 0.753 |
| T1r4c3 | 0.718 | 0.716 | unambiguous | phase2_leaderboard.csv[3].ci95#1 (0.718→0.716); phase8_leaderboard_hk.csv[3].dl_lo (0.718227→0.716248) | 0.718 – 0.784 |
| T1r4c3 | 0.784 | 0.783 | unambiguous | phase2_leaderboard.csv[3].ci95#2 (0.784→0.783); phase8_leaderboard_hk.csv[3].dl_hi (0.784002→0.783206) | 0.718 – 0.784 |
| T1r4c4 | 58 | 59 | unambiguous | phase2_leaderboard.csv[3].I2 (57.8→58.9); phase8_leaderboard_hk.csv[3].I2 (57.7916→58.8765) | 58 |
| T1r5c3 | 0.722 | 0.723 | unambiguous | phase2_leaderboard.csv[4].ci95#1 (0.722→0.723); phase8_leaderboard_hk.csv[4].dl_lo (0.72229→0.723071) | 0.722 – 0.776 |
| T1r5c4 | 38 | 35 | unambiguous | phase2_leaderboard.csv[4].I2 (37.8→35.1); phase8_leaderboard_hk.csv[4].I2 (37.8127→35.0584) | 38 |
| T1r11c2 | 0.494 | 0.492 | unambiguous | phase2_leaderboard.csv[10].pooled_rho (0.4942→0.4923); phase8_leaderboard_hk.csv[10].pooled_rho (0.494233→0.492287) | 0.494 |
| T1r11c3 | 0.374 | 0.371 | unambiguous | phase2_leaderboard.csv[10].ci95#1 (0.374→0.371); phase8_leaderboard_hk.csv[10].dl_lo (0.373562→0.371312) | 0.374 – 0.598 |
| T1r11c3 | 0.598 | 0.597 | unambiguous | phase2_leaderboard.csv[10].ci95#2 (0.598→0.597); phase8_leaderboard_hk.csv[10].dl_hi (0.598454→0.596827) | 0.374 – 0.598 |
| T2r1c2 | 0.773 | 0.774 | unambiguous | phase2_H2_ablation.csv[0].full_rho (0.773→0.7738) | 0.773 |
| T2r1c3 | 0.774 | 0.773 | unambiguous | phase2_H2_ablation.csv[0].ablated_rho (0.7739→0.7732) | 0.774 |
| T2r1c4 | -0.001 | 0.001 | unambiguous | phase2_H2_ablation.csv[0].delta_full_minus_ablated (-0.0009→0.0007) | -0.001 |
| T2r1c5 | -0.003 | −0.001 | unambiguous | phase2_H2_ablation.csv[0].ci95#1 (-0.0025→-0.001) | -0.003 to +0.001 |
| T2r1c5 | +0.001 | +0.002 | unambiguous | phase2_H2_ablation.csv[0].ci95#2 (0.0009→0.0021) | -0.003 to +0.001 |
| T2r2c2 | 0.773 | 0.774 | unambiguous | phase2_H2_ablation.csv[1].full_rho (0.773→0.7738) | 0.773 |
| T2r2c4 | +0.002 | +0.003 | unambiguous | phase2_H2_ablation.csv[1].delta_full_minus_ablated (0.0022→0.0032) | +0.002 |
| T2r2c5 | -0.001 | 0.001 | unambiguous | phase2_H2_ablation.csv[1].ci95#1 (-0.0014→0.0005) | -0.001 to +0.007 |
| T2r2c5 | +0.007 | +0.006 | unambiguous | phase2_H2_ablation.csv[1].ci95#2 (0.0066→0.0064) | -0.001 to +0.007 |
| T2r3c2 | 0.759 | 0.758 | unambiguous | phase2_H2_ablation.csv[2].full_rho (0.7589→0.7577) | 0.759 |
| T2r3c3 | 0.758 | 0.759 | unambiguous | phase2_H2_ablation.csv[2].ablated_rho (0.7583→0.7595) | 0.758 |
| T2r3c4 | +0.001 | −0.002 | unambiguous | phase2_H2_ablation.csv[2].delta_full_minus_ablated (0.0005→-0.0018) | +0.001 |
| T2r3c5 | -0.006 | −0.011 | unambiguous | phase2_H2_ablation.csv[2].ci95#1 (-0.0056→-0.0108) | -0.006 to +0.004 |
| T2r4c2 | 0.759 | 0.758 | unambiguous | phase2_H2_ablation.csv[3].full_rho (0.7589→0.7577) | 0.759 |
| T2r4c3 | 0.754 | 0.757 | unambiguous | phase2_H2_ablation.csv[3].ablated_rho (0.7543→0.7573) | 0.754 |
| T2r4c4 | +0.005 | +0.000 | unambiguous | phase2_H2_ablation.csv[3].delta_full_minus_ablated (0.0046→0.0004) | +0.005 |
| T2r4c5 | -0.011 | −0.013 | unambiguous | phase2_H2_ablation.csv[3].ci95#1 (-0.0109→-0.0127) | -0.011 to +0.015 |
| T2r4c5 | +0.015 | +0.010 | unambiguous | phase2_H2_ablation.csv[3].ci95#2 (0.0154→0.01) | -0.011 to +0.015 |
| T3r1c2 | 0.026 | 0.025 | unambiguous | phase3_calibration_summary.csv[0].ECE (0.026→0.0249) | 0.026 |
| T3r1c4 | 75.1 | 74.9 | unambiguous | phase3_calibration_summary.csv[0].actionable_frac (75.1→74.87) | 75.1 |
| T3r1c5 | 96.8 | 96.7 | unambiguous | phase3_calibration_summary.csv[0].actionable_acc (96.82→96.73) | 96.8 |
| T3r2c3 | 0.072 | 0.073 | unambiguous | phase3_calibration_summary.csv[3].Brier (0.0717→0.0733) | 0.072 |
| T3r2c4 | 75.2 | 73.6 | unambiguous | phase3_calibration_summary.csv[3].actionable_frac (75.22→73.55) | 75.2 |
| T3r3c2 | 0.012 | 0.015 | unambiguous | phase3_calibration_summary.csv[1].ECE (0.0118→0.0155) | 0.012 |
| T3r3c3 | 0.077 | 0.078 | unambiguous | phase3_calibration_summary.csv[1].Brier (0.0775→0.0776) | 0.077 |
| T3r3c4 | 69.7 | 70.0 | unambiguous | phase3_calibration_summary.csv[1].actionable_frac (69.73→70.03) | 69.7 |
| T3r5c2 | 0.047 | 0.050 | unambiguous | phase3_calibration_summary.csv[2].ECE (0.047→0.0495) | 0.047 |
| T3r5c3 | 0.081 | 0.080 | unambiguous | phase3_calibration_summary.csv[2].Brier (0.0807→0.0804) | 0.081 |
| T3r5c4 | 74.3 | 70.5 | unambiguous | phase3_calibration_summary.csv[2].actionable_frac (74.33→70.45) | 74.3 |
| T3r5c5 | 94.5 | 94.4 | unambiguous | phase3_calibration_summary.csv[2].actionable_acc (94.46→94.41) | 94.5 |
| T3r10c2 | 0.038 | 0.048 | unambiguous | phase3_calibration_summary.csv[6].ECE (0.0385→0.0475) | 0.038 |
| T3r10c3 | 0.198 | 0.199 | unambiguous | phase3_calibration_summary.csv[6].Brier (0.1981→0.1986) | 0.198 |
| T3r10c5 | 79.0 | 78.6 | unambiguous | phase3_calibration_summary.csv[6].actionable_acc (79.05→78.62) | 79.0 |
| T4r1c2 | 0.043 | 0.041 | unambiguous | phase4_tp53_label_definitions.csv[0].dBrier (0.0430949→0.0410328) | 0.043 |
| T4r1c3 | 0.016 | 0.025 | unambiguous | phase4_tp53_label_definitions.csv[0].ci_lo (0.0155127→0.0245144) | 0.016 – 0.071 |
| T4r1c3 | 0.071 | 0.059 | unambiguous | phase4_tp53_label_definitions.csv[0].ci_hi (0.0707184→0.0589604) | 0.016 – 0.071 |
| T4r2c2 | 0.038 | 0.005 | unambiguous | phase4_tp53_label_definitions.csv[1].dBrier (0.0380774→0.00511928) | 0.038 |
| T4r2c3 | 0.010 | −0.014 | unambiguous | phase4_tp53_label_definitions.csv[1].ci_lo (0.0102309→-0.0143973) | 0.010 – 0.068 |
| T4r2c3 | 0.068 | 0.024 | unambiguous | phase4_tp53_label_definitions.csv[1].ci_hi (0.0680576→0.0237496) | 0.010 – 0.068 |
| T4r3c2 | 0.044 | 0.041 | unambiguous | phase4_tp53_label_definitions.csv[2].dBrier (0.0439214→0.0412353) | 0.044 |
| T4r3c3 | 0.016 | 0.024 | decided | phase4_tp53_label_definitions.csv[2].ci_lo (0.0159163→0.0238157) | 0.016 – 0.075 |
| T4r3c3 | 0.075 | 0.061 | unambiguous | phase4_tp53_label_definitions.csv[2].ci_hi (0.0754446→0.0605688) | 0.016 – 0.075 |

## Ambiguous tokens needing a decision (0)

| block | token | candidate new values | cells | context |
|---|---|---|---|---|

## Tokens with no v1 cell within tolerance in their section's outputs, or no v2 counterpart (20)

| block | token | status | context |
|---|---|---|---|
| P18 | 21,410 | no_v1_source | The analysis reuses a frozen matrix of 21,410 GRCh38 single-nucleotide variants |
| P18 | 15 | no_v1_source | nVar records were taken from the GRCh38 VCF release of 15 June 2026 ; records of every revi |
| P18 | 583 | no_v1_source | ron-side splice region: splice-core (/offset/ ≤ 2; n = 583, including six audited UTR-intron |
| P18 | 1,191 | no_v1_source | on classifier had missed) plus splice-region (3–8; n = 1,191) plus seven exon-side variants at |
| P18 | 1,781 | no_v1_source | s them as coding/UTR — and are grouped with the core — 1,781 in all, so the region-stratified H |
| P18 | 590 | no_v1_source | region-stratified H1 analysis contrasts core-like (n = 590) with splice-region (n = 1,191). T |
| P18 | 11 | no_v1_source | the region band 3–8, where the atlas stratifies 3–10, 11–50 and >50 bp. Of the nine UTR-int |
| P18 | 50 | no_v1_source | e region band 3–8, where the atlas stratifies 3–10, 11–50 and >50 bp. Of the nine UTR-intron |
| P18 | −2 | no_v1_source | identified in the shared audit, the six at offsets −1/−2 are corrected into the splice core |
| P18 | −3 | no_v1_source | −1/−2 are corrected into the splice core; the three at −3 carry no usable intron offset in t |
| P19 | 1,781 | no_v1_source | normal = 0, indeterminate = NA). This covered 1,675 of 1,781 splice variants (94%) at a near-ba |
| P19 | 94 | no_v1_source | te = NA). This covered 1,675 of 1,781 splice variants (94%) at a near-balanced positive frac |
| P29 | 96.5 | no_v2_counterpart | st ratios on the functional standard (Pangolin 21.1 at 96.5% realised, AlphaGenome 24.3 at 97. |
| P39 | 21 | no_v1_source | that test either. Leave-two-genes-out refits over all 21 gene pairs keep Δρ positive throug |
| P50 | 1,768 | no_v1_source | t sensitivity behind the modest difference between the 1,768-variant reweighting baseline and t |
| P50 | 1,781 | no_v1_source | between the 1,768-variant reweighting baseline and the 1,781-variant headline in 2.4. The fusio |
| P62 | 4.20 | no_v1_source | conflicting and other records — by ratios from 1.44 to 4.20, and the fusion's sensitivity at 9 |
| P78 | 11 | no_v1_source | g is not established; the window ends at 8 bp, and the 11–50 bp band in which the companion |
| P78 | 50 | no_v1_source | s not established; the window ends at 8 bp, and the 11–50 bp band in which the companion atl |
| P78 | 6 | no_v1_source | nd reporting results without BRCA1. Labels abstain on ~6% of splice variants, excluded from |

## Applied: 174 tokens replaced in place; not located in text: 0

