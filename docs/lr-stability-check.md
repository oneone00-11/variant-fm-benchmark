# LR+ stability check: v1 versus v2 (read-only)

Purpose: explain the two evidence-tier changes the v1→v2 delta showed (Pangolin, isotonic-calibrated, primary condition 16.4 → 19.9; BRCA1-excluded 21.1 → 17.4) by separating *score change* from *threshold placement*. Score vectors were reconstructed with the phase-5 code path itself (same features, LOGO fusion, isotonic calibration, min–max scaling, seed); the reconstruction reproduces every `lr_plus` in `phase5_likelihood_ratios.csv` to 1e-9 for both versions (96/96 rows each). Nothing in the repository was modified; this file is the only output.

Notation: *chosen* = the most sensitive observed score value whose FPR ≤ 5% (the scan `lr_at_specificity` uses); *below* = the next lower unique value, i.e. the candidate the scan rejected because its FPR exceeded 5%; *above* = the next higher unique value. Tie sizes count variants at exactly that value (all / negatives / positives). Conditions are the four label sets; stage is raw or isotonic-calibrated (the `cal_method = isotonic` rows of the CSV).

## 1. Achieved specificity at the scanned threshold, v1 versus v2

| condition | stage | object | thr v1 | spec v1 | sens v1 | LR+ v1 (tier) | thr v2 | spec v2 | sens v2 | LR+ v2 (tier) | ΔLR+ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| y_assay/BRCA1_included | raw | fusion_M1 | 0.5644 | 95.08% | 0.889 | 18.08 (Moderate) | 0.5765 | 95.08% | 0.889 | 18.08 (Moderate) | +0.00 |
| y_assay/BRCA1_included | calibrated_isotonic | fusion_M1 | 0.7333 | 95.08% | 0.881 | 17.91 (Moderate) | 0.7612 | 95.08% | 0.877 | 17.83 (Moderate) | -0.07 |
| y_assay/BRCA1_included | raw | mean_M0b | 0.6154 | 95.08% | 0.815 | 16.57 (Moderate) | 0.6152 | 95.08% | 0.816 | 16.59 (Moderate) | +0.02 |
| y_assay/BRCA1_included | calibrated_isotonic | mean_M0b | 0.7778 | 95.20% | 0.787 | 16.39 (Moderate) | 0.7500 | 95.20% | 0.784 | 16.34 (Moderate) | -0.05 |
| y_assay/BRCA1_included | raw | single:spliceai | 0.5671 | 95.08% | 0.832 | 16.92 (Moderate) | 0.5811 | 95.08% | 0.831 | 16.89 (Moderate) | -0.02 |
| y_assay/BRCA1_included | calibrated_isotonic | single:spliceai | 0.8571 | 95.67% | 0.739 | 17.06 (Moderate) | 0.8389 | 95.43% | 0.743 | 16.27 (Moderate) | -0.80 |
| y_assay/BRCA1_included | raw | single:pangolin | 0.5485 | 95.08% | 0.839 | 17.06 (Moderate) | 0.5788 | 95.08% | 0.838 | 17.04 (Moderate) | -0.02 |
| y_assay/BRCA1_included | calibrated_isotonic | single:pangolin | 0.7882 | 95.08% | 0.806 | **16.40 (Moderate)** | 0.7986 | 96.02% | 0.792 | **19.89 (Strong)** | +3.49 |
| y_assay/BRCA1_included | raw | single:alphagenome | 0.6099 | 95.08% | 0.773 | 15.73 (Moderate) | 0.6099 | 95.08% | 0.773 | 15.73 (Moderate) | +0.00 |
| y_assay/BRCA1_included | calibrated_isotonic | single:alphagenome | 0.8205 | 95.08% | 0.730 | 14.84 (Moderate) | 0.8205 | 95.08% | 0.730 | 14.84 (Moderate) | +0.00 |
| y_assay/BRCA1_included | raw | single:gpn_msa | 0.6728 | 95.08% | 0.636 | 12.93 (Moderate) | 0.6728 | 95.08% | 0.636 | 12.93 (Moderate) | +0.00 |
| y_assay/BRCA1_included | calibrated_isotonic | single:gpn_msa | 0.8514 | 96.25% | 0.590 | 15.73 (Moderate) | 0.8514 | 96.25% | 0.590 | 15.73 (Moderate) | +0.00 |
| y_assay/BRCA1_included | raw | single:nt | 0.8236 | 95.08% | 0.318 | 6.46 (Moderate) | 0.8231 | 95.08% | 0.318 | 6.46 (Moderate) | +0.00 |
| y_assay/BRCA1_included | calibrated_isotonic | single:nt | 0.8167 | 96.02% | 0.290 | 7.28 (Moderate) | 0.8214 | 96.14% | 0.253 | 6.56 (Moderate) | -0.72 |
| y_assay/BRCA1_included | raw | single:cadd | 0.6730 | 95.10% | 0.634 | 12.95 (Moderate) | 0.6730 | 95.10% | 0.634 | 12.95 (Moderate) | +0.00 |
| y_assay/BRCA1_included | calibrated_isotonic | single:cadd | 0.8790 | 96.65% | 0.538 | 16.07 (Moderate) | 0.8790 | 96.65% | 0.538 | 16.07 (Moderate) | +0.00 |
| y_assay/BRCA1_included | raw | single:alphamissense | — | — | — | not evaluable | — | — | — | not evaluable | |
| y_assay/BRCA1_included | calibrated_isotonic | single:alphamissense | — | — | — | not evaluable | — | — | — | not evaluable | |
| y_assay/BRCA1_included | raw | single:gnomad_af | 0.9624 | 96.26% | 0.051 | 1.35 (below supporting) | 0.9624 | 96.26% | 0.051 | 1.35 (below supporting) | +0.00 |
| y_assay/BRCA1_included | calibrated_isotonic | single:gnomad_af | 0.5854 | 96.26% | 0.045 | 1.21 (below supporting) | 0.5854 | 96.26% | 0.045 | 1.21 (below supporting) | +0.00 |
| y_assay/BRCA1_included | raw | single:phylop | 0.6610 | 95.08% | 0.669 | 13.60 (Moderate) | 0.6610 | 95.08% | 0.669 | 13.60 (Moderate) | +0.00 |
| y_assay/BRCA1_included | calibrated_isotonic | single:phylop | 0.7969 | 95.08% | 0.630 | 12.80 (Moderate) | 0.7969 | 95.08% | 0.630 | 12.80 (Moderate) | +0.00 |
| y_assay/BRCA1_included | raw | single:phastcons | 0.9885 | 98.01% | 0.132 | 6.61 (Moderate) | 0.9885 | 98.01% | 0.132 | 6.61 (Moderate) | +0.00 |
| y_assay/BRCA1_included | calibrated_isotonic | single:phastcons | 1.0000 | 98.01% | 0.104 | 5.20 (Moderate) | 1.0000 | 98.01% | 0.104 | 5.20 (Moderate) | +0.00 |
| y_assay/BRCA1_excluded | raw | fusion_M1 | 0.5586 | 95.12% | 0.875 | 17.95 (Moderate) | 0.5661 | 95.12% | 0.879 | 18.02 (Moderate) | +0.07 |
| y_assay/BRCA1_excluded | calibrated_isotonic | fusion_M1 | 0.7000 | 95.30% | 0.864 | 18.36 (Moderate) | 0.7000 | 95.12% | 0.864 | 17.70 (Moderate) | -0.66 |
| y_assay/BRCA1_excluded | raw | mean_M0b | 0.6150 | 95.12% | 0.813 | 16.67 (Moderate) | 0.6142 | 95.12% | 0.811 | 16.63 (Moderate) | -0.03 |
| y_assay/BRCA1_excluded | calibrated_isotonic | mean_M0b | 0.7857 | 95.30% | 0.790 | 16.79 (Moderate) | 0.7667 | 95.47% | 0.788 | 17.39 (Moderate) | +0.61 |
| y_assay/BRCA1_excluded | raw | single:spliceai | 0.5828 | 95.30% | 0.790 | 16.79 (Moderate) | 0.5937 | 95.12% | 0.783 | 16.05 (Moderate) | -0.74 |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:spliceai | 0.8902 | 95.12% | 0.682 | 13.98 (Moderate) | 0.8852 | 95.12% | 0.705 | 14.46 (Moderate) | +0.48 |
| y_assay/BRCA1_excluded | raw | single:pangolin | 0.5550 | 95.12% | 0.808 | 16.57 (Moderate) | 0.5820 | 95.12% | 0.813 | 16.67 (Moderate) | +0.10 |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:pangolin | 0.8095 | 96.52% | 0.736 | **21.11 (Strong)** | 0.8102 | 95.64% | 0.759 | **17.43 (Moderate)** | -3.68 |
| y_assay/BRCA1_excluded | raw | single:alphagenome | 0.6099 | 95.12% | 0.753 | 15.43 (Moderate) | 0.6099 | 95.12% | 0.753 | 15.43 (Moderate) | +0.00 |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:alphagenome | 0.8750 | 97.21% | 0.678 | 24.34 (Strong) | 0.8750 | 97.21% | 0.678 | 24.34 (Strong) | +0.00 |
| y_assay/BRCA1_excluded | raw | single:gpn_msa | 0.6604 | 95.12% | 0.652 | 13.36 (Moderate) | 0.6604 | 95.12% | 0.652 | 13.36 (Moderate) | +0.00 |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:gpn_msa | 0.8723 | 96.17% | 0.530 | 13.84 (Moderate) | 0.8723 | 96.17% | 0.530 | 13.84 (Moderate) | +0.00 |
| y_assay/BRCA1_excluded | raw | single:nt | 0.7485 | 95.12% | 0.466 | 9.56 (Moderate) | 0.7485 | 95.12% | 0.466 | 9.56 (Moderate) | +0.00 |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:nt | 0.8229 | 95.12% | 0.409 | 8.39 (Moderate) | 0.8276 | 95.64% | 0.380 | 8.74 (Moderate) | +0.35 |
| y_assay/BRCA1_excluded | raw | single:cadd | 0.5922 | 95.01% | 0.774 | 15.52 (Moderate) | 0.5922 | 95.01% | 0.774 | 15.52 (Moderate) | +0.00 |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:cadd | 0.6667 | 95.01% | 0.771 | 15.44 (Moderate) | 0.6667 | 95.01% | 0.771 | 15.44 (Moderate) | +0.00 |
| y_assay/BRCA1_excluded | raw | single:alphamissense | — | — | — | not evaluable | — | — | — | not evaluable | |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:alphamissense | — | — | — | not evaluable | — | — | — | not evaluable | |
| y_assay/BRCA1_excluded | raw | single:gnomad_af | 0.9600 | 95.10% | 0.074 | 1.51 (below supporting) | 0.9600 | 95.10% | 0.074 | 1.51 (below supporting) | +0.00 |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:gnomad_af | 0.6047 | 95.51% | 0.056 | 1.24 (below supporting) | 0.6047 | 95.51% | 0.056 | 1.24 (below supporting) | +0.00 |
| y_assay/BRCA1_excluded | raw | single:phylop | 0.6610 | 95.12% | 0.660 | 13.53 (Moderate) | 0.6610 | 95.12% | 0.660 | 13.53 (Moderate) | +0.00 |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:phylop | 0.8182 | 95.47% | 0.618 | 13.64 (Moderate) | 0.8182 | 95.47% | 0.618 | 13.64 (Moderate) | +0.00 |
| y_assay/BRCA1_excluded | raw | single:phastcons | 0.9885 | 97.04% | 0.182 | 6.14 (Moderate) | 0.9885 | 97.04% | 0.182 | 6.14 (Moderate) | +0.00 |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:phastcons | 1.0000 | 97.04% | 0.143 | 4.83 (Moderate) | 1.0000 | 97.04% | 0.143 | 4.83 (Moderate) | +0.00 |
| y_clinvar/BRCA1_included | raw | fusion_M1 | 0.4363 | 95.14% | 0.998 | 20.54 (Strong) | 0.4541 | 95.14% | 0.998 | 20.54 (Strong) | +0.00 |
| y_clinvar/BRCA1_included | calibrated_isotonic | fusion_M1 | 0.0365 | 96.16% | 0.993 | 25.88 (Strong) | 0.1462 | 96.42% | 0.993 | 27.73 (Strong) | +1.85 |
| y_clinvar/BRCA1_included | raw | mean_M0b | 0.4634 | 95.14% | 0.996 | 20.50 (Strong) | 0.4658 | 95.14% | 0.996 | 20.50 (Strong) | +0.00 |
| y_clinvar/BRCA1_included | calibrated_isotonic | mean_M0b | 0.1579 | 95.40% | 0.993 | 21.57 (Strong) | 0.2857 | 96.68% | 0.991 | 29.81 (Strong) | +8.24 |
| y_clinvar/BRCA1_included | raw | single:spliceai | 0.5054 | 95.14% | 0.977 | 20.10 (Strong) | 0.5165 | 95.14% | 0.977 | 20.10 (Strong) | +0.00 |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:spliceai | 0.5000 | 95.14% | 0.959 | 19.73 (Strong) | 0.5678 | 95.40% | 0.953 | 20.70 (Strong) | +0.98 |
| y_clinvar/BRCA1_included | raw | single:pangolin | 0.4478 | 95.40% | 0.996 | 21.64 (Strong) | 0.4762 | 95.14% | 0.996 | 20.50 (Strong) | -1.14 |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:pangolin | 0.4131 | 95.40% | 0.978 | 21.25 (Strong) | 0.2900 | 95.14% | 0.978 | 20.13 (Strong) | -1.12 |
| y_clinvar/BRCA1_included | raw | single:alphagenome | 0.5054 | 95.14% | 0.960 | 19.76 (Strong) | 0.5054 | 95.14% | 0.960 | 19.76 (Strong) | +0.00 |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:alphagenome | 0.7013 | 96.16% | 0.919 | 23.95 (Strong) | 0.7013 | 96.16% | 0.919 | 23.95 (Strong) | +0.00 |
| y_clinvar/BRCA1_included | raw | single:gpn_msa | 0.5420 | 95.14% | 0.903 | 18.58 (Moderate) | 0.5420 | 95.14% | 0.903 | 18.58 (Moderate) | +0.00 |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:gpn_msa | 0.6923 | 95.14% | 0.899 | 18.50 (Moderate) | 0.6923 | 95.14% | 0.899 | 18.50 (Moderate) | +0.00 |
| y_clinvar/BRCA1_included | raw | single:nt | 0.7008 | 95.14% | 0.582 | 11.98 (Moderate) | 0.7029 | 95.14% | 0.575 | 11.83 (Moderate) | -0.15 |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:nt | 0.8307 | 95.40% | 0.580 | 12.60 (Moderate) | 0.8285 | 95.40% | 0.571 | 12.41 (Moderate) | -0.20 |
| y_clinvar/BRCA1_included | raw | single:cadd | 0.4699 | 95.03% | 0.971 | 19.53 (Strong) | 0.4699 | 95.03% | 0.971 | 19.53 (Strong) | +0.00 |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:cadd | 0.6000 | 95.58% | 0.956 | 21.62 (Strong) | 0.6000 | 95.58% | 0.956 | 21.62 (Strong) | +0.00 |
| y_clinvar/BRCA1_included | raw | single:alphamissense | — | — | — | not evaluable | — | — | — | not evaluable | |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:alphamissense | — | — | — | not evaluable | — | — | — | not evaluable | |
| y_clinvar/BRCA1_included | raw | single:gnomad_af | 0.9247 | 95.03% | 0.165 | 3.33 (Supporting) | 0.9247 | 95.03% | 0.165 | 3.33 (Supporting) | +0.00 |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:gnomad_af | 0.6923 | 95.65% | 0.102 | 2.35 (Supporting) | 0.6923 | 95.65% | 0.102 | 2.35 (Supporting) | +0.00 |
| y_clinvar/BRCA1_included | raw | single:phylop | 0.5481 | 95.14% | 0.924 | 19.02 (Strong) | 0.5481 | 95.14% | 0.924 | 19.02 (Strong) | +0.00 |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:phylop | 0.6178 | 96.68% | 0.919 | 27.64 (Strong) | 0.6178 | 96.68% | 0.919 | 27.64 (Strong) | +0.00 |
| y_clinvar/BRCA1_included | raw | single:phastcons | 0.9102 | 95.14% | 0.791 | 16.28 (Moderate) | 0.9102 | 95.14% | 0.791 | 16.28 (Moderate) | +0.00 |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:phastcons | 0.9332 | 95.14% | 0.791 | 16.28 (Moderate) | 0.9332 | 95.14% | 0.791 | 16.28 (Moderate) | +0.00 |
| y_clinvar/BRCA1_excluded | raw | fusion_M1 | 0.4381 | 95.12% | 0.998 | 20.45 (Strong) | 0.4564 | 95.12% | 0.998 | 20.45 (Strong) | +0.00 |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | fusion_M1 | 0.0416 | 96.04% | 0.988 | 24.92 (Strong) | 0.1767 | 96.34% | 0.988 | 27.00 (Strong) | +2.08 |
| y_clinvar/BRCA1_excluded | raw | mean_M0b | 0.4658 | 95.12% | 0.995 | 20.40 (Strong) | 0.4677 | 95.12% | 0.995 | 20.40 (Strong) | +0.00 |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | mean_M0b | 0.1755 | 95.43% | 0.988 | 21.60 (Strong) | 0.2671 | 96.65% | 0.985 | 29.38 (Strong) | +7.78 |
| y_clinvar/BRCA1_excluded | raw | single:spliceai | 0.5167 | 95.43% | 0.968 | 21.17 (Strong) | 0.5198 | 95.12% | 0.970 | 19.89 (Strong) | -1.27 |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:spliceai | 0.6154 | 95.43% | 0.924 | 20.20 (Strong) | 0.6137 | 95.12% | 0.936 | 19.19 (Strong) | -1.01 |
| y_clinvar/BRCA1_excluded | raw | single:pangolin | 0.4478 | 95.12% | 0.995 | 20.40 (Strong) | 0.4861 | 95.12% | 0.995 | 20.40 (Strong) | +0.00 |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:pangolin | 0.6856 | 95.43% | 0.963 | 21.06 (Strong) | 0.6154 | 95.12% | 0.958 | 19.64 (Strong) | -1.42 |
| y_clinvar/BRCA1_excluded | raw | single:alphagenome | 0.5018 | 95.12% | 0.953 | 19.54 (Strong) | 0.5018 | 95.12% | 0.953 | 19.54 (Strong) | +0.00 |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:alphagenome | 0.7442 | 95.73% | 0.906 | 21.24 (Strong) | 0.7442 | 95.73% | 0.906 | 21.24 (Strong) | +0.00 |
| y_clinvar/BRCA1_excluded | raw | single:gpn_msa | 0.5450 | 95.12% | 0.906 | 18.58 (Moderate) | 0.5450 | 95.12% | 0.906 | 18.58 (Moderate) | +0.00 |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:gpn_msa | 0.6667 | 95.73% | 0.884 | 20.72 (Strong) | 0.6667 | 95.73% | 0.884 | 20.72 (Strong) | +0.00 |
| y_clinvar/BRCA1_excluded | raw | single:nt | 0.6956 | 95.12% | 0.623 | 12.77 (Moderate) | 0.6888 | 95.12% | 0.635 | 13.03 (Moderate) | +0.25 |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:nt | 0.7826 | 95.43% | 0.631 | 13.79 (Moderate) | 0.7857 | 95.43% | 0.621 | 13.57 (Moderate) | -0.22 |
| y_clinvar/BRCA1_excluded | raw | single:cadd | 0.4699 | 95.03% | 0.979 | 19.70 (Strong) | 0.4699 | 95.03% | 0.979 | 19.70 (Strong) | +0.00 |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:cadd | 0.4355 | 95.03% | 0.971 | 19.54 (Strong) | 0.4355 | 95.03% | 0.971 | 19.54 (Strong) | +0.00 |
| y_clinvar/BRCA1_excluded | raw | single:alphamissense | — | — | — | not evaluable | — | — | — | not evaluable | |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:alphamissense | — | — | — | not evaluable | — | — | — | not evaluable | |
| y_clinvar/BRCA1_excluded | raw | single:gnomad_af | 0.9594 | 95.52% | 0.090 | 2.01 (below supporting) | 0.9594 | 95.52% | 0.090 | 2.01 (below supporting) | +0.00 |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:gnomad_af | 0.7500 | 96.27% | 0.030 | 0.80 (below supporting) | 0.7500 | 96.27% | 0.030 | 0.80 (below supporting) | +0.00 |
| y_clinvar/BRCA1_excluded | raw | single:phylop | 0.5553 | 95.12% | 0.929 | 19.04 (Strong) | 0.5553 | 95.12% | 0.929 | 19.04 (Strong) | +0.00 |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:phylop | 0.5474 | 96.65% | 0.929 | 27.69 (Strong) | 0.5474 | 96.65% | 0.929 | 27.69 (Strong) | +0.00 |
| y_clinvar/BRCA1_excluded | raw | single:phastcons | 0.9102 | 96.04% | 0.778 | 19.64 (Strong) | 0.9102 | 96.04% | 0.778 | 19.64 (Strong) | +0.00 |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:phastcons | 0.9332 | 96.04% | 0.778 | 19.64 (Strong) | 0.9332 | 96.04% | 0.778 | 19.64 (Strong) | +0.00 |

Achieved specificity runs from 95.01% to 98.01% in v1 and from 95.01% to 98.01% in v2. Rows where it moves by more than 0.5 percentage points: y_assay/BRCA1_included / calibrated_isotonic / single:pangolin (95.08% → 96.02%); y_assay/BRCA1_excluded / calibrated_isotonic / single:pangolin (96.52% → 95.64%); y_assay/BRCA1_excluded / calibrated_isotonic / single:nt (95.12% → 95.64%); y_clinvar/BRCA1_included / calibrated_isotonic / mean_M0b (95.40% → 96.68%); y_clinvar/BRCA1_excluded / calibrated_isotonic / mean_M0b (95.43% → 96.65%). Bold marks the only two cells whose ACMG tier changes.

## 2. Tie groups at and around the scanned threshold

| condition | stage | object | ver | distinct | below: value / tie all,neg,pos / spec | chosen: value / tie all,neg,pos / spec | above: value / tie all,neg,pos / spec |
|---|---|---|---|---|---|---|---|
| y_assay/BRCA1_included | raw | fusion_M1 | v1 | 1675 | 0.5641 / 1,1,0 / 94.96% | 0.5644 / 1,0,1 / 95.08% | 0.5664 / 1,0,1 / 95.08% |
| y_assay/BRCA1_included | raw | fusion_M1 | v2 | 1675 | 0.5741 / 1,1,0 / 94.96% | 0.5765 / 1,0,1 / 95.08% | 0.5775 / 1,0,1 / 95.08% |
| y_assay/BRCA1_included | calibrated_isotonic | fusion_M1 | v1 | 131 | 0.7273 / 7,2,5 / 94.85% | 0.7333 / 3,1,2 / 95.08% | 0.7647 / 5,0,5 / 95.20% |
| y_assay/BRCA1_included | calibrated_isotonic | fusion_M1 | v2 | 128 | 0.7273 / 7,2,5 / 94.85% | 0.7612 / 10,0,10 / 95.08% | 0.7692 / 3,2,1 / 95.08% |
| y_assay/BRCA1_included | raw | mean_M0b | v1 | 1675 | 0.6150 / 1,1,0 / 94.96% | 0.6154 / 1,0,1 / 95.08% | 0.6157 / 1,0,1 / 95.08% |
| y_assay/BRCA1_included | raw | mean_M0b | v2 | 1673 | 0.6142 / 1,1,0 / 94.96% | 0.6152 / 1,0,1 / 95.08% | 0.6180 / 1,0,1 / 95.08% |
| y_assay/BRCA1_included | calibrated_isotonic | mean_M0b | v1 | 165 | 0.7647 / 3,2,1 / 94.96% | 0.7778 / 6,2,4 / 95.20% | 0.7868 / 1,0,1 / 95.43% |
| y_assay/BRCA1_included | calibrated_isotonic | mean_M0b | v2 | 137 | 0.7368 / 3,2,1 / 94.96% | 0.7500 / 2,0,2 / 95.20% | 0.7637 / 1,0,1 / 95.20% |
| y_assay/BRCA1_included | raw | single:spliceai | v1 | 374 | 0.5652 / 1,1,0 / 94.96% | 0.5671 / 1,0,1 / 95.08% | 0.5686 / 2,0,2 / 95.08% |
| y_assay/BRCA1_included | raw | single:spliceai | v2 | 1670 | 0.5806 / 1,1,0 / 94.96% | 0.5811 / 1,0,1 / 95.08% | 0.5815 / 1,0,1 / 95.08% |
| y_assay/BRCA1_included | calibrated_isotonic | single:spliceai | v1 | 108 | 0.8529 / 8,6,2 / 94.96% | 0.8571 / 4,3,1 / 95.67% | 0.8656 / 1,0,1 / 96.02% |
| y_assay/BRCA1_included | calibrated_isotonic | single:spliceai | v2 | 133 | 0.8387 / 6,4,2 / 94.96% | 0.8389 / 1,0,1 / 95.43% | 0.8438 / 4,0,4 / 95.43% |
| y_assay/BRCA1_included | raw | single:pangolin | v1 | 383 | 0.5477 / 2,2,0 / 94.85% | 0.5485 / 3,0,3 / 95.08% | 0.5537 / 1,0,1 / 95.08% |
| y_assay/BRCA1_included | raw | single:pangolin | v2 | 1672 | 0.5786 / 1,1,0 / 94.96% | 0.5788 / 1,0,1 / 95.08% | 0.5806 / 1,0,1 / 95.08% |
| y_assay/BRCA1_included | calibrated_isotonic | single:pangolin | v1 | 115 | 0.7800 / 6,4,2 / 94.61% | 0.7882 / 8,0,8 / 95.08% | 0.7895 / 1,1,0 / 95.08% |
| y_assay/BRCA1_included | calibrated_isotonic | single:pangolin | v2 | 129 | 0.7973 / 15,10,5 / 94.85% | 0.7986 / 1,0,1 / 96.02% | 0.8000 / 14,1,13 / 96.02% |
| y_assay/BRCA1_included | raw | single:alphagenome | v1 | 1611 | 0.6093 / 1,1,0 / 94.96% | 0.6099 / 1,1,0 / 95.08% | 0.6099 / 1,1,0 / 95.20% |
| y_assay/BRCA1_included | raw | single:alphagenome | v2 | 1611 | 0.6093 / 1,1,0 / 94.96% | 0.6099 / 1,1,0 / 95.08% | 0.6099 / 1,1,0 / 95.20% |
| y_assay/BRCA1_included | calibrated_isotonic | single:alphagenome | v1 | 134 | 0.8148 / 6,2,4 / 94.85% | 0.8205 / 5,1,4 / 95.08% | 0.8235 / 2,0,2 / 95.20% |
| y_assay/BRCA1_included | calibrated_isotonic | single:alphagenome | v2 | 134 | 0.8148 / 6,2,4 / 94.85% | 0.8205 / 5,1,4 / 95.08% | 0.8235 / 2,0,2 / 95.20% |
| y_assay/BRCA1_included | raw | single:gpn_msa | v1 | 1467 | 0.6722 / 1,1,0 / 94.96% | 0.6728 / 1,1,0 / 95.08% | 0.6729 / 1,0,1 / 95.20% |
| y_assay/BRCA1_included | raw | single:gpn_msa | v2 | 1467 | 0.6722 / 1,1,0 / 94.96% | 0.6728 / 1,1,0 / 95.08% | 0.6729 / 1,0,1 / 95.20% |
| y_assay/BRCA1_included | calibrated_isotonic | single:gpn_msa | v1 | 153 | 0.8219 / 31,11,20 / 94.96% | 0.8514 / 8,0,8 / 96.25% | 0.8571 / 12,1,11 / 96.25% |
| y_assay/BRCA1_included | calibrated_isotonic | single:gpn_msa | v2 | 153 | 0.8219 / 31,11,20 / 94.96% | 0.8514 / 8,0,8 / 96.25% | 0.8571 / 12,1,11 / 96.25% |
| y_assay/BRCA1_included | raw | single:nt | v1 | 1671 | 0.8231 / 1,1,0 / 94.96% | 0.8236 / 1,0,1 / 95.08% | 0.8240 / 1,0,1 / 95.08% |
| y_assay/BRCA1_included | raw | single:nt | v2 | 1671 | 0.8226 / 1,1,0 / 94.96% | 0.8231 / 1,1,0 / 95.08% | 0.8236 / 1,1,0 / 95.20% |
| y_assay/BRCA1_included | calibrated_isotonic | single:nt | v1 | 153 | 0.8145 / 50,17,33 / 94.03% | 0.8167 / 12,1,11 / 96.02% | 0.8214 / 8,0,8 / 96.14% |
| y_assay/BRCA1_included | calibrated_isotonic | single:nt | v2 | 155 | 0.8197 / 51,17,34 / 94.15% | 0.8214 / 11,1,10 / 96.14% | 0.8219 / 10,1,9 / 96.25% |
| y_assay/BRCA1_included | raw | single:cadd | v1 | 1120 | 0.6716 / 3,1,2 / 94.97% | 0.6730 / 1,0,1 / 95.10% | 0.6733 / 2,0,2 / 95.10% |
| y_assay/BRCA1_included | raw | single:cadd | v2 | 1120 | 0.6716 / 3,1,2 / 94.97% | 0.6730 / 1,0,1 / 95.10% | 0.6733 / 2,0,2 / 95.10% |
| y_assay/BRCA1_included | calibrated_isotonic | single:cadd | v1 | 113 | 0.8713 / 46,15,31 / 94.72% | 0.8790 / 2,2,0 / 96.65% | 0.9385 / 9,0,9 / 96.91% |
| y_assay/BRCA1_included | calibrated_isotonic | single:cadd | v2 | 113 | 0.8713 / 46,15,31 / 94.72% | 0.8790 / 2,2,0 / 96.65% | 0.9385 / 9,0,9 / 96.91% |
| y_assay/BRCA1_included | raw | single:alphamissense | v1 | — | — | — | — |
| y_assay/BRCA1_included | raw | single:alphamissense | v2 | — | — | — | — |
| y_assay/BRCA1_included | calibrated_isotonic | single:alphamissense | v1 | — | — | — | — |
| y_assay/BRCA1_included | calibrated_isotonic | single:alphamissense | v2 | — | — | — | — |
| y_assay/BRCA1_included | raw | single:gnomad_af | v1 | 371 | 0.9600 / 7,4,3 / 94.90% | 0.9624 / 1,1,0 / 96.26% | 0.9628 / 6,2,4 / 96.60% |
| y_assay/BRCA1_included | raw | single:gnomad_af | v2 | 371 | 0.9600 / 7,4,3 / 94.90% | 0.9624 / 1,1,0 / 96.26% | 0.9628 / 6,2,4 / 96.60% |
| y_assay/BRCA1_included | calibrated_isotonic | single:gnomad_af | v1 | 56 | 0.5818 / 11,7,4 / 93.88% | 0.5854 / 8,4,4 / 96.26% | 0.5952 / 5,3,2 / 97.62% |
| y_assay/BRCA1_included | calibrated_isotonic | single:gnomad_af | v2 | 56 | 0.5818 / 11,7,4 / 93.88% | 0.5854 / 8,4,4 / 96.26% | 0.5952 / 5,3,2 / 97.62% |
| y_assay/BRCA1_included | raw | single:phylop | v1 | 749 | 0.6592 / 6,1,5 / 94.96% | 0.6610 / 2,2,0 / 95.08% | 0.6616 / 2,1,1 / 95.32% |
| y_assay/BRCA1_included | raw | single:phylop | v2 | 749 | 0.6592 / 6,1,5 / 94.96% | 0.6610 / 2,2,0 / 95.08% | 0.6616 / 2,1,1 / 95.32% |
| y_assay/BRCA1_included | calibrated_isotonic | single:phylop | v1 | 133 | 0.7959 / 17,3,14 / 94.73% | 0.7969 / 7,2,5 / 95.08% | 0.8000 / 4,0,4 / 95.32% |
| y_assay/BRCA1_included | calibrated_isotonic | single:phylop | v2 | 133 | 0.7959 / 17,3,14 / 94.73% | 0.7969 / 7,2,5 / 95.08% | 0.8000 / 4,0,4 / 95.32% |
| y_assay/BRCA1_included | raw | single:phastcons | v1 | 222 | 0.9829 / 207,39,168 / 93.44% | 0.9885 / 23,0,23 / 98.01% | 1.0000 / 102,17,85 / 98.01% |
| y_assay/BRCA1_included | raw | single:phastcons | v2 | 222 | 0.9829 / 207,39,168 / 93.44% | 0.9885 / 23,0,23 / 98.01% | 1.0000 / 102,17,85 / 98.01% |
| y_assay/BRCA1_included | calibrated_isotonic | single:phastcons | v1 | 54 | 0.8830 / 207,39,168 / 93.44% | 1.0000 / 102,17,85 / 98.01% | — |
| y_assay/BRCA1_included | calibrated_isotonic | single:phastcons | v2 | 54 | 0.8830 / 207,39,168 / 93.44% | 1.0000 / 102,17,85 / 98.01% | — |
| y_assay/BRCA1_excluded | raw | fusion_M1 | v1 | 1168 | 0.5558 / 1,1,0 / 94.95% | 0.5586 / 1,1,0 / 95.12% | 0.5592 / 1,0,1 / 95.30% |
| y_assay/BRCA1_excluded | raw | fusion_M1 | v2 | 1168 | 0.5659 / 1,1,0 / 94.95% | 0.5661 / 1,0,1 / 95.12% | 0.5666 / 1,0,1 / 95.12% |
| y_assay/BRCA1_excluded | calibrated_isotonic | fusion_M1 | v1 | 106 | 0.6818 / 9,4,5 / 94.60% | 0.7000 / 1,0,1 / 95.30% | 0.7500 / 13,1,12 / 95.30% |
| y_assay/BRCA1_excluded | calibrated_isotonic | fusion_M1 | v2 | 99 | 0.6923 / 10,5,5 / 94.25% | 0.7000 / 1,0,1 / 95.12% | 0.7270 / 1,1,0 / 95.12% |
| y_assay/BRCA1_excluded | raw | mean_M0b | v1 | 1168 | 0.6109 / 1,1,0 / 94.95% | 0.6150 / 1,1,0 / 95.12% | 0.6154 / 1,0,1 / 95.30% |
| y_assay/BRCA1_excluded | raw | mean_M0b | v2 | 1166 | 0.6129 / 1,1,0 / 94.95% | 0.6142 / 1,1,0 / 95.12% | 0.6184 / 1,0,1 / 95.30% |
| y_assay/BRCA1_excluded | calibrated_isotonic | mean_M0b | v1 | 116 | 0.7500 / 3,2,1 / 94.95% | 0.7857 / 10,3,7 / 95.30% | 0.7931 / 3,0,3 / 95.82% |
| y_assay/BRCA1_excluded | calibrated_isotonic | mean_M0b | v2 | 111 | 0.7500 / 7,4,3 / 94.77% | 0.7667 / 5,0,5 / 95.47% | 0.7812 / 3,0,3 / 95.47% |
| y_assay/BRCA1_excluded | raw | single:spliceai | v1 | 292 | 0.5803 / 4,3,1 / 94.77% | 0.5828 / 1,1,0 / 95.30% | 0.5865 / 7,0,7 / 95.47% |
| y_assay/BRCA1_excluded | raw | single:spliceai | v2 | 1165 | 0.5932 / 1,1,0 / 94.95% | 0.5937 / 1,1,0 / 95.12% | 0.5955 / 1,0,1 / 95.30% |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:spliceai | v1 | 90 | 0.8889 / 3,2,1 / 94.77% | 0.8902 / 18,0,18 / 95.12% | 0.8947 / 6,0,6 / 95.12% |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:spliceai | v2 | 99 | 0.8571 / 7,3,4 / 94.60% | 0.8852 / 16,0,16 / 95.12% | 0.8889 / 3,2,1 / 95.12% |
| y_assay/BRCA1_excluded | raw | single:pangolin | v1 | 302 | 0.5541 / 3,1,2 / 94.95% | 0.5550 / 2,2,0 / 95.12% | 0.5565 / 1,0,1 / 95.47% |
| y_assay/BRCA1_excluded | raw | single:pangolin | v2 | 1166 | 0.5815 / 1,1,0 / 94.95% | 0.5820 / 1,0,1 / 95.12% | 0.5823 / 1,0,1 / 95.12% |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:pangolin | v1 | 88 | 0.8085 / 15,9,6 / 94.95% | 0.8095 / 3,0,3 / 96.52% | 0.8113 / 1,1,0 / 96.52% |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:pangolin | v2 | 98 | 0.8000 / 6,5,1 / 94.77% | 0.8102 / 1,0,1 / 95.64% | 0.8148 / 8,4,4 / 95.64% |
| y_assay/BRCA1_excluded | raw | single:alphagenome | v1 | 1129 | 0.6093 / 1,1,0 / 94.95% | 0.6099 / 1,1,0 / 95.12% | 0.6115 / 1,0,1 / 95.30% |
| y_assay/BRCA1_excluded | raw | single:alphagenome | v2 | 1129 | 0.6093 / 1,1,0 / 94.95% | 0.6099 / 1,1,0 / 95.12% | 0.6115 / 1,0,1 / 95.30% |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:alphagenome | v1 | 92 | 0.8673 / 29,14,15 / 94.77% | 0.8750 / 14,1,13 / 97.21% | 0.8763 / 1,0,1 / 97.39% |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:alphagenome | v2 | 92 | 0.8673 / 29,14,15 / 94.77% | 0.8750 / 14,1,13 / 97.21% | 0.8763 / 1,0,1 / 97.39% |
| y_assay/BRCA1_excluded | raw | single:gpn_msa | v1 | 1065 | 0.6572 / 2,2,0 / 94.77% | 0.6604 / 1,0,1 / 95.12% | 0.6610 / 1,0,1 / 95.12% |
| y_assay/BRCA1_excluded | raw | single:gpn_msa | v2 | 1065 | 0.6572 / 2,2,0 / 94.77% | 0.6604 / 1,0,1 / 95.12% | 0.6610 / 1,0,1 / 95.12% |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:gpn_msa | v1 | 116 | 0.8696 / 33,8,25 / 94.77% | 0.8723 / 18,5,13 / 96.17% | 0.8765 / 4,1,3 / 97.04% |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:gpn_msa | v2 | 116 | 0.8696 / 33,8,25 / 94.77% | 0.8723 / 18,5,13 / 96.17% | 0.8765 / 4,1,3 / 97.04% |
| y_assay/BRCA1_excluded | raw | single:nt | v1 | 1165 | 0.7480 / 1,1,0 / 94.95% | 0.7485 / 1,0,1 / 95.12% | 0.7505 / 2,0,2 / 95.12% |
| y_assay/BRCA1_excluded | raw | single:nt | v2 | 1165 | 0.7480 / 1,1,0 / 94.95% | 0.7485 / 1,0,1 / 95.12% | 0.7505 / 2,0,2 / 95.12% |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:nt | v1 | 106 | 0.8220 / 6,2,4 / 94.77% | 0.8229 / 28,6,22 / 95.12% | 0.8318 / 32,9,23 / 96.17% |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:nt | v2 | 121 | 0.8242 / 27,5,22 / 94.77% | 0.8276 / 6,2,4 / 95.64% | 0.8333 / 2,0,2 / 95.99% |
| y_assay/BRCA1_excluded | raw | single:cadd | v1 | 770 | 0.5918 / 1,1,0 / 94.82% | 0.5922 / 2,0,2 / 95.01% | 0.5937 / 1,1,0 / 95.01% |
| y_assay/BRCA1_excluded | raw | single:cadd | v2 | 770 | 0.5918 / 1,1,0 / 94.82% | 0.5922 / 2,0,2 / 95.01% | 0.5937 / 1,1,0 / 95.01% |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:cadd | v1 | 93 | 0.6250 / 4,4,0 / 94.24% | 0.6667 / 3,0,3 / 95.01% | 0.7143 / 1,0,1 / 95.01% |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:cadd | v2 | 93 | 0.6250 / 4,4,0 / 94.24% | 0.6667 / 3,0,3 / 95.01% | 0.7143 / 1,0,1 / 95.01% |
| y_assay/BRCA1_excluded | raw | single:alphamissense | v1 | — | — | — | — |
| y_assay/BRCA1_excluded | raw | single:alphamissense | v2 | — | — | — | — |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:alphamissense | v1 | — | — | — | — |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:alphamissense | v2 | — | — | — | — |
| y_assay/BRCA1_excluded | raw | single:gnomad_af | v1 | 347 | 0.9594 / 2,1,1 / 94.69% | 0.9600 / 7,4,3 / 95.10% | 0.9628 / 6,2,4 / 96.73% |
| y_assay/BRCA1_excluded | raw | single:gnomad_af | v2 | 347 | 0.9594 / 2,1,1 / 94.69% | 0.9600 / 7,4,3 / 95.10% | 0.9628 / 6,2,4 / 96.73% |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:gnomad_af | v1 | 52 | 0.6000 / 5,3,2 / 94.29% | 0.6047 / 11,6,5 / 95.51% | 0.6098 / 9,5,4 / 97.96% |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:gnomad_af | v2 | 52 | 0.6000 / 5,3,2 / 94.29% | 0.6047 / 11,6,5 / 95.51% | 0.6098 / 9,5,4 / 97.96% |
| y_assay/BRCA1_excluded | raw | single:phylop | v1 | 604 | 0.6575 / 2,1,1 / 94.95% | 0.6610 / 2,2,0 / 95.12% | 0.6616 / 2,1,1 / 95.47% |
| y_assay/BRCA1_excluded | raw | single:phylop | v2 | 604 | 0.6575 / 2,1,1 / 94.95% | 0.6610 / 2,2,0 / 95.12% | 0.6616 / 2,1,1 / 95.47% |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:phylop | v1 | 107 | 0.7895 / 6,4,2 / 94.77% | 0.8182 / 3,0,3 / 95.47% | 0.8204 / 4,1,3 / 95.47% |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:phylop | v2 | 107 | 0.7895 / 6,4,2 / 94.77% | 0.8182 / 3,0,3 / 95.47% | 0.8204 / 4,1,3 / 95.47% |
| y_assay/BRCA1_excluded | raw | single:phastcons | v1 | 181 | 0.9472 / 123,13,110 / 94.77% | 0.9885 / 23,0,23 / 97.04% | 1.0000 / 102,17,85 / 97.04% |
| y_assay/BRCA1_excluded | raw | single:phastcons | v2 | 181 | 0.9472 / 123,13,110 / 94.77% | 0.9885 / 23,0,23 / 97.04% | 1.0000 / 102,17,85 / 97.04% |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:phastcons | v1 | 40 | 0.8771 / 123,13,110 / 94.77% | 1.0000 / 102,17,85 / 97.04% | — |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:phastcons | v2 | 40 | 0.8771 / 123,13,110 / 94.77% | 1.0000 / 102,17,85 / 97.04% | — |
| y_clinvar/BRCA1_included | raw | fusion_M1 | v1 | 946 | 0.4355 / 1,1,0 / 94.88% | 0.4363 / 1,1,0 / 95.14% | 0.4381 / 1,1,0 / 95.40% |
| y_clinvar/BRCA1_included | raw | fusion_M1 | v2 | 946 | 0.4483 / 1,1,0 / 94.88% | 0.4541 / 1,1,0 / 95.14% | 0.4564 / 1,1,0 / 95.40% |
| y_clinvar/BRCA1_included | calibrated_isotonic | fusion_M1 | v1 | 25 | 0.0333 / 12,12,0 / 93.09% | 0.0365 / 1,1,0 / 96.16% | 0.1383 / 1,1,0 / 96.42% |
| y_clinvar/BRCA1_included | calibrated_isotonic | fusion_M1 | v2 | 23 | 0.0345 / 14,14,0 / 92.84% | 0.1462 / 1,0,1 / 96.42% | 0.1732 / 1,1,0 / 96.42% |
| y_clinvar/BRCA1_included | raw | mean_M0b | v1 | 946 | 0.4625 / 1,1,0 / 94.88% | 0.4634 / 1,1,0 / 95.14% | 0.4658 / 1,1,0 / 95.40% |
| y_clinvar/BRCA1_included | raw | mean_M0b | v2 | 944 | 0.4656 / 1,1,0 / 94.88% | 0.4658 / 1,1,0 / 95.14% | 0.4677 / 1,1,0 / 95.40% |
| y_clinvar/BRCA1_included | calibrated_isotonic | mean_M0b | v1 | 37 | 0.1429 / 4,4,0 / 94.37% | 0.1579 / 1,1,0 / 95.40% | 0.1970 / 1,1,0 / 95.65% |
| y_clinvar/BRCA1_included | calibrated_isotonic | mean_M0b | v2 | 35 | 0.2000 / 8,7,1 / 94.88% | 0.2857 / 1,0,1 / 96.68% | 0.2866 / 1,1,0 / 96.68% |
| y_clinvar/BRCA1_included | raw | single:spliceai | v1 | 205 | 0.4925 / 1,1,0 / 94.88% | 0.5054 / 1,0,1 / 95.14% | 0.5076 / 2,2,0 / 95.14% |
| y_clinvar/BRCA1_included | raw | single:spliceai | v2 | 945 | 0.5128 / 1,1,0 / 94.88% | 0.5165 / 1,1,0 / 95.14% | 0.5198 / 1,0,1 / 95.40% |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:spliceai | v1 | 35 | 0.4545 / 1,1,0 / 94.88% | 0.5000 / 3,1,2 / 95.14% | 0.5006 / 1,1,0 / 95.40% |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:spliceai | v2 | 43 | 0.5625 / 2,2,0 / 94.88% | 0.5678 / 1,1,0 / 95.40% | 0.5714 / 2,1,1 / 95.65% |
| y_clinvar/BRCA1_included | raw | single:pangolin | v1 | 226 | 0.4413 / 2,2,0 / 94.88% | 0.4478 / 1,0,1 / 95.40% | 0.4512 / 1,1,0 / 95.40% |
| y_clinvar/BRCA1_included | raw | single:pangolin | v2 | 945 | 0.4744 / 1,1,0 / 94.88% | 0.4762 / 1,1,0 / 95.14% | 0.4861 / 1,0,1 / 95.40% |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:pangolin | v1 | 34 | 0.3333 / 3,3,0 / 94.63% | 0.4131 / 2,2,0 / 95.40% | 0.5000 / 6,3,3 / 95.91% |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:pangolin | v2 | 35 | 0.2500 / 8,3,5 / 94.37% | 0.2900 / 1,1,0 / 95.14% | 0.4706 / 4,1,3 / 95.40% |
| y_clinvar/BRCA1_included | raw | single:alphagenome | v1 | 917 | 0.5037 / 1,1,0 / 94.88% | 0.5054 / 1,0,1 / 95.14% | 0.5055 / 1,1,0 / 95.14% |
| y_clinvar/BRCA1_included | raw | single:alphagenome | v2 | 917 | 0.5037 / 1,1,0 / 94.88% | 0.5054 / 1,0,1 / 95.14% | 0.5055 / 1,1,0 / 95.14% |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:alphagenome | v1 | 42 | 0.6957 / 6,6,0 / 94.63% | 0.7013 / 1,1,0 / 96.16% | 0.8462 / 4,4,0 / 96.42% |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:alphagenome | v2 | 42 | 0.6957 / 6,6,0 / 94.63% | 0.7013 / 1,1,0 / 96.16% | 0.8462 / 4,4,0 / 96.42% |
| y_clinvar/BRCA1_included | raw | single:gpn_msa | v1 | 890 | 0.5417 / 1,1,0 / 94.88% | 0.5420 / 1,1,0 / 95.14% | 0.5450 / 1,0,1 / 95.40% |
| y_clinvar/BRCA1_included | raw | single:gpn_msa | v2 | 890 | 0.5417 / 1,1,0 / 94.88% | 0.5420 / 1,1,0 / 95.14% | 0.5450 / 1,0,1 / 95.40% |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:gpn_msa | v1 | 77 | 0.6667 / 2,2,0 / 94.63% | 0.6923 / 5,0,5 / 95.14% | 0.7143 / 5,0,5 / 95.14% |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:gpn_msa | v2 | 77 | 0.6667 / 2,2,0 / 94.63% | 0.6923 / 5,0,5 / 95.14% | 0.7143 / 5,0,5 / 95.14% |
| y_clinvar/BRCA1_included | raw | single:nt | v1 | 945 | 0.7000 / 1,1,0 / 94.88% | 0.7008 / 1,1,0 / 95.14% | 0.7029 / 1,0,1 / 95.40% |
| y_clinvar/BRCA1_included | raw | single:nt | v2 | 946 | 0.7008 / 1,1,0 / 94.88% | 0.7029 / 1,1,0 / 95.14% | 0.7038 / 1,0,1 / 95.40% |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:nt | v1 | 104 | 0.7778 / 8,3,5 / 94.63% | 0.8307 / 1,0,1 / 95.40% | 0.8571 / 9,0,9 / 95.40% |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:nt | v2 | 107 | 0.8182 / 4,3,1 / 94.63% | 0.8285 / 1,0,1 / 95.40% | 0.8421 / 7,0,7 / 95.40% |
| y_clinvar/BRCA1_included | raw | single:cadd | v1 | 575 | 0.4677 / 1,1,0 / 94.75% | 0.4699 / 1,1,0 / 95.03% | 0.4729 / 1,0,1 / 95.30% |
| y_clinvar/BRCA1_included | raw | single:cadd | v2 | 575 | 0.4677 / 1,1,0 / 94.75% | 0.4699 / 1,1,0 / 95.03% | 0.4729 / 1,0,1 / 95.30% |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:cadd | v1 | 58 | 0.5714 / 3,3,0 / 94.75% | 0.6000 / 3,0,3 / 95.58% | 0.6006 / 1,0,1 / 95.58% |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:cadd | v2 | 58 | 0.5714 / 3,3,0 / 94.75% | 0.6000 / 3,0,3 / 95.58% | 0.6006 / 1,0,1 / 95.58% |
| y_clinvar/BRCA1_included | raw | single:alphamissense | v1 | — | — | — | — |
| y_clinvar/BRCA1_included | raw | single:alphamissense | v2 | — | — | — | — |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:alphamissense | v1 | — | — | — | — |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:alphamissense | v2 | — | — | — | — |
| y_clinvar/BRCA1_included | raw | single:gnomad_af | v1 | 235 | 0.9082 / 1,1,0 / 94.41% | 0.9247 / 1,0,1 / 95.03% | 0.9324 / 1,0,1 / 95.03% |
| y_clinvar/BRCA1_included | raw | single:gnomad_af | v2 | 235 | 0.9082 / 1,1,0 / 94.41% | 0.9247 / 1,0,1 / 95.03% | 0.9324 / 1,0,1 / 95.03% |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:gnomad_af | v1 | 60 | 0.6774 / 7,4,3 / 93.17% | 0.6923 / 3,0,3 / 95.65% | 0.7065 / 1,0,1 / 95.65% |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:gnomad_af | v2 | 60 | 0.6774 / 7,4,3 / 93.17% | 0.6923 / 3,0,3 / 95.65% | 0.7065 / 1,0,1 / 95.65% |
| y_clinvar/BRCA1_included | raw | single:phylop | v1 | 489 | 0.5443 / 2,2,0 / 94.63% | 0.5481 / 1,0,1 / 95.14% | 0.5482 / 1,0,1 / 95.14% |
| y_clinvar/BRCA1_included | raw | single:phylop | v2 | 489 | 0.5443 / 2,2,0 / 94.63% | 0.5481 / 1,0,1 / 95.14% | 0.5482 / 1,0,1 / 95.14% |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:phylop | v1 | 72 | 0.6111 / 8,8,0 / 94.63% | 0.6178 / 1,1,0 / 96.68% | 0.6429 / 9,0,9 / 96.93% |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:phylop | v2 | 72 | 0.6111 / 8,8,0 / 94.63% | 0.6178 / 1,1,0 / 96.68% | 0.6429 / 9,0,9 / 96.93% |
| y_clinvar/BRCA1_included | raw | single:phastcons | v1 | 148 | 0.8767 / 56,5,51 / 93.86% | 0.9102 / 48,7,41 / 95.14% | 0.9323 / 77,4,73 / 96.93% |
| y_clinvar/BRCA1_included | raw | single:phastcons | v2 | 148 | 0.8767 / 56,5,51 / 93.86% | 0.9102 / 48,7,41 / 95.14% | 0.9323 / 77,4,73 / 96.93% |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:phastcons | v1 | 53 | 0.8689 / 6,6,0 / 93.61% | 0.9332 / 48,7,41 / 95.14% | 0.9386 / 77,4,73 / 96.93% |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:phastcons | v2 | 53 | 0.8689 / 6,6,0 / 93.61% | 0.9332 / 48,7,41 / 95.14% | 0.9386 / 77,4,73 / 96.93% |
| y_clinvar/BRCA1_excluded | raw | fusion_M1 | v1 | 734 | 0.4363 / 1,1,0 / 94.82% | 0.4381 / 1,1,0 / 95.12% | 0.4432 / 1,1,0 / 95.43% |
| y_clinvar/BRCA1_excluded | raw | fusion_M1 | v2 | 734 | 0.4541 / 1,1,0 / 94.82% | 0.4564 / 1,1,0 / 95.12% | 0.4569 / 1,1,0 / 95.43% |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | fusion_M1 | v1 | 19 | 0.0385 / 12,12,0 / 92.38% | 0.0416 / 1,1,0 / 96.04% | 0.1416 / 1,1,0 / 96.34% |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | fusion_M1 | v2 | 17 | 0.0400 / 14,14,0 / 92.07% | 0.1767 / 1,1,0 / 96.34% | 0.2108 / 1,1,0 / 96.65% |
| y_clinvar/BRCA1_excluded | raw | mean_M0b | v1 | 734 | 0.4634 / 1,1,0 / 94.82% | 0.4658 / 1,1,0 / 95.12% | 0.4660 / 1,0,1 / 95.43% |
| y_clinvar/BRCA1_excluded | raw | mean_M0b | v2 | 732 | 0.4658 / 1,1,0 / 94.82% | 0.4677 / 1,1,0 / 95.12% | 0.4683 / 1,0,1 / 95.43% |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | mean_M0b | v1 | 30 | 0.1667 / 4,4,0 / 94.21% | 0.1755 / 1,1,0 / 95.43% | 0.1983 / 1,1,0 / 95.73% |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | mean_M0b | v2 | 25 | 0.2222 / 6,6,0 / 94.82% | 0.2671 / 1,1,0 / 96.65% | 0.2841 / 1,1,0 / 96.95% |
| y_clinvar/BRCA1_excluded | raw | single:spliceai | v1 | 156 | 0.5076 / 2,2,0 / 94.82% | 0.5167 / 1,0,1 / 95.43% | 0.5205 / 1,0,1 / 95.43% |
| y_clinvar/BRCA1_excluded | raw | single:spliceai | v2 | 734 | 0.5165 / 1,1,0 / 94.82% | 0.5198 / 1,0,1 / 95.12% | 0.5274 / 1,1,0 / 95.12% |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:spliceai | v1 | 26 | 0.6111 / 9,2,7 / 94.82% | 0.6154 / 10,0,10 / 95.43% | 0.7826 / 2,0,2 / 95.43% |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:spliceai | v2 | 35 | 0.6111 / 4,2,2 / 94.51% | 0.6137 / 1,0,1 / 95.12% | 0.6466 / 1,0,1 / 95.12% |
| y_clinvar/BRCA1_excluded | raw | single:pangolin | v1 | 177 | 0.4413 / 2,2,0 / 94.51% | 0.4478 / 1,0,1 / 95.12% | 0.4514 / 1,0,1 / 95.12% |
| y_clinvar/BRCA1_excluded | raw | single:pangolin | v2 | 733 | 0.4762 / 1,1,0 / 94.82% | 0.4861 / 1,0,1 / 95.12% | 0.4872 / 1,1,0 / 95.12% |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:pangolin | v1 | 29 | 0.6471 / 2,2,0 / 94.82% | 0.6856 / 1,0,1 / 95.43% | 0.7500 / 10,5,5 / 95.43% |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:pangolin | v2 | 33 | 0.5815 / 1,1,0 / 94.82% | 0.6154 / 1,1,0 / 95.12% | 0.6271 / 1,1,0 / 95.43% |
| y_clinvar/BRCA1_excluded | raw | single:alphagenome | v1 | 713 | 0.5000 / 1,1,0 / 94.82% | 0.5018 / 1,1,0 / 95.12% | 0.5025 / 1,0,1 / 95.43% |
| y_clinvar/BRCA1_excluded | raw | single:alphagenome | v2 | 713 | 0.5000 / 1,1,0 / 94.82% | 0.5018 / 1,1,0 / 95.12% | 0.5025 / 1,0,1 / 95.43% |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:alphagenome | v1 | 41 | 0.7000 / 4,4,0 / 94.51% | 0.7442 / 1,0,1 / 95.73% | 0.7500 / 2,2,0 / 95.73% |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:alphagenome | v2 | 41 | 0.7000 / 4,4,0 / 94.51% | 0.7442 / 1,0,1 / 95.73% | 0.7500 / 2,2,0 / 95.73% |
| y_clinvar/BRCA1_excluded | raw | single:gpn_msa | v1 | 698 | 0.5420 / 1,1,0 / 94.82% | 0.5450 / 1,0,1 / 95.12% | 0.5494 / 1,1,0 / 95.12% |
| y_clinvar/BRCA1_excluded | raw | single:gpn_msa | v2 | 698 | 0.5420 / 1,1,0 / 94.82% | 0.5450 / 1,0,1 / 95.12% | 0.5494 / 1,1,0 / 95.12% |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:gpn_msa | v1 | 62 | 0.6364 / 6,5,1 / 94.21% | 0.6667 / 10,1,9 / 95.73% | 0.7143 / 5,1,4 / 96.04% |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:gpn_msa | v2 | 62 | 0.6364 / 6,5,1 / 94.21% | 0.6667 / 10,1,9 / 95.73% | 0.7143 / 5,1,4 / 96.04% |
| y_clinvar/BRCA1_excluded | raw | single:nt | v1 | 733 | 0.6928 / 1,1,0 / 94.82% | 0.6956 / 1,0,1 / 95.12% | 0.6968 / 1,1,0 / 95.12% |
| y_clinvar/BRCA1_excluded | raw | single:nt | v2 | 734 | 0.6878 / 1,1,0 / 94.82% | 0.6888 / 1,0,1 / 95.12% | 0.6892 / 1,0,1 / 95.12% |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:nt | v1 | 78 | 0.7812 / 6,3,3 / 94.51% | 0.7826 / 8,3,5 / 95.43% | 0.8788 / 9,0,9 / 96.34% |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:nt | v2 | 90 | 0.7778 / 6,2,4 / 94.82% | 0.7857 / 3,0,3 / 95.43% | 0.8096 / 1,0,1 / 95.43% |
| y_clinvar/BRCA1_excluded | raw | single:cadd | v1 | 441 | 0.4677 / 1,1,0 / 94.70% | 0.4699 / 1,1,0 / 95.03% | 0.4760 / 1,1,0 / 95.36% |
| y_clinvar/BRCA1_excluded | raw | single:cadd | v2 | 441 | 0.4677 / 1,1,0 / 94.70% | 0.4699 / 1,1,0 / 95.03% | 0.4760 / 1,1,0 / 95.36% |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:cadd | v1 | 26 | 0.4286 / 1,1,0 / 94.70% | 0.4355 / 1,0,1 / 95.03% | 0.4444 / 3,0,3 / 95.03% |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:cadd | v2 | 26 | 0.4286 / 1,1,0 / 94.70% | 0.4355 / 1,0,1 / 95.03% | 0.4444 / 3,0,3 / 95.03% |
| y_clinvar/BRCA1_excluded | raw | single:alphamissense | v1 | — | — | — | — |
| y_clinvar/BRCA1_excluded | raw | single:alphamissense | v2 | — | — | — | — |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:alphamissense | v1 | — | — | — | — |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:alphamissense | v2 | — | — | — | — |
| y_clinvar/BRCA1_excluded | raw | single:gnomad_af | v1 | 204 | 0.9355 / 9,1,8 / 94.78% | 0.9594 / 1,0,1 / 95.52% | 0.9600 / 3,2,1 / 95.52% |
| y_clinvar/BRCA1_excluded | raw | single:gnomad_af | v2 | 204 | 0.9355 / 9,1,8 / 94.78% | 0.9594 / 1,0,1 / 95.52% | 0.9600 / 3,2,1 / 95.52% |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:gnomad_af | v1 | 43 | 0.7368 / 5,2,3 / 94.78% | 0.7500 / 1,0,1 / 96.27% | 0.7502 / 1,1,0 / 96.27% |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:gnomad_af | v2 | 43 | 0.7368 / 5,2,3 / 94.78% | 0.7500 / 1,0,1 / 96.27% | 0.7502 / 1,1,0 / 96.27% |
| y_clinvar/BRCA1_excluded | raw | single:phylop | v1 | 397 | 0.5527 / 2,1,1 / 94.82% | 0.5553 / 2,2,0 / 95.12% | 0.5599 / 1,1,0 / 95.73% |
| y_clinvar/BRCA1_excluded | raw | single:phylop | v2 | 397 | 0.5527 / 2,1,1 / 94.82% | 0.5553 / 2,2,0 / 95.12% | 0.5599 / 1,1,0 / 95.73% |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:phylop | v1 | 49 | 0.5385 / 7,7,0 / 94.51% | 0.5474 / 1,1,0 / 96.65% | 0.6111 / 10,0,10 / 96.95% |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:phylop | v2 | 49 | 0.5385 / 7,7,0 / 94.51% | 0.5474 / 1,1,0 / 96.65% | 0.6111 / 10,0,10 / 96.95% |
| y_clinvar/BRCA1_excluded | raw | single:phastcons | v1 | 124 | 0.8767 / 56,5,51 / 94.51% | 0.9102 / 48,7,41 / 96.04% | 0.9323 / 77,4,73 / 98.17% |
| y_clinvar/BRCA1_excluded | raw | single:phastcons | v2 | 124 | 0.8767 / 56,5,51 / 94.51% | 0.9102 / 48,7,41 / 96.04% | 0.9323 / 77,4,73 / 98.17% |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:phastcons | v1 | 42 | 0.7789 / 56,5,51 / 94.51% | 0.9332 / 48,7,41 / 96.04% | 0.9475 / 77,4,73 / 98.17% |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:phastcons | v2 | 42 | 0.7789 / 56,5,51 / 94.51% | 0.9332 / 48,7,41 / 96.04% | 0.9475 / 77,4,73 / 98.17% |

### The two cells that change tier

**Pangolin, isotonic-calibrated, primary condition (n = 1675, 854 negatives; one negative = 0.117 pp of FPR).** In v1 the scan stops at calibrated value 0.7882, a tie group of 8 variants (all positives), with FPR = 4.92% (42 negatives above it): achieved specificity 95.08%, LR+ 16.40. In v2 the value just below the 5% line, 0.7973, is a tie group of 15 (10 negatives, 5 positives) that would put FPR at 5.15% — over the line by 2 negative(s) — so the scan takes the next value up, 0.7986 (a single positive), where FPR is 3.98% (34 negatives) and specificity 96.02%. Sensitivity actually *falls* (0.806 → 0.792); the LR+ rises from 16.40 to 19.89 only because the denominator FPR shrinks from 0.0492 to 0.0398. Crossing 18.7 here is a placement effect of a 10-negative tie group created by isotonic calibration, not a discrimination gain.

**Pangolin, isotonic-calibrated, BRCA1 excluded (n = 1168, 574 negatives; one negative = 0.174 pp).** The same mechanism in reverse. In v1 the value below the line, 0.8085, is a tie of 15 (9 negatives) giving FPR 5.05%, so the scan jumps to 0.8095 at specificity 96.52% (20 negatives) and LR+ 21.11 — the v1 'Strong' rested on landing 1.5 points above 95%. In v2 the tie below is smaller (6, 5 negatives, FPR 5.23%), the chosen value 0.8102 lands at 95.64%, sensitivity is *higher* than in v1 (0.736 → 0.759), and LR+ is 17.43. The v2 score is at least as discriminating; the tier drop is the loss of a favourable placement.

## 3. Sensitivity analysis: LR+ at exactly 95.0% specificity by linear interpolation of the ROC

Between the two adjacent unique-threshold operating points that bracket FPR = 0.05 (point a: the chosen threshold, FPR ≤ 0.05; point b: the next lower value, FPR > 0.05), sensitivity is interpolated linearly to FPR = 0.05 and LR+ = sens / 0.05. This removes the dependence on where an observed value happens to fall; it is not the paper's estimator (which is defined on observed values so that a tie group cannot drive the FPR towards 1), it is a check on it.

| condition | stage | object | LR+ obs v1 | LR+ interp v1 | LR+ obs v2 | LR+ interp v2 | Δ obs | Δ interp | tier interp v1 → v2 |
|---|---|---|---|---|---|---|---|---|---|
| y_assay/BRCA1_included | raw | fusion_M1 | 18.08 | 17.78 | 18.08 | 17.78 | +0.00 | +0.00 | Moderate → Moderate |
| y_assay/BRCA1_included | calibrated_isotonic | fusion_M1 | 17.91 | 17.66 | 17.83 | 17.58 | -0.07 | -0.07 | Moderate → Moderate |
| y_assay/BRCA1_included | raw | mean_M0b | 16.57 | 16.30 | 16.59 | 16.32 | +0.02 | +0.02 | Moderate → Moderate |
| y_assay/BRCA1_included | calibrated_isotonic | mean_M0b | 16.39 | 15.76 | 16.34 | 15.71 | -0.05 | -0.05 | Moderate → Moderate |
| y_assay/BRCA1_included | raw | single:spliceai | 16.92 | 16.64 | 16.89 | 16.61 | -0.02 | -0.02 | Moderate → Moderate |
| y_assay/BRCA1_included | calibrated_isotonic | single:spliceai | 17.06 | 14.83 | 16.27 | 14.90 | -0.80 | +0.07 | Moderate → Moderate |
| y_assay/BRCA1_included | raw | single:pangolin | 17.06 | 16.78 | 17.04 | 16.76 | -0.02 | -0.02 | Moderate → Moderate |
| y_assay/BRCA1_included | calibrated_isotonic | single:pangolin | 16.40 | 16.14 | 19.89 | 15.94 | +3.49 | -0.19 | Moderate → Moderate |
| y_assay/BRCA1_included | raw | single:alphagenome | 15.73 | 15.47 | 15.73 | 15.47 | +0.00 | +0.00 | Moderate → Moderate |
| y_assay/BRCA1_included | calibrated_isotonic | single:alphagenome | 14.84 | 14.63 | 14.84 | 14.63 | +0.00 | +0.00 | Moderate → Moderate |
| y_assay/BRCA1_included | raw | single:gpn_msa | 12.93 | 12.72 | 12.93 | 12.72 | +0.00 | +0.00 | Moderate → Moderate |
| y_assay/BRCA1_included | calibrated_isotonic | single:gpn_msa | 15.73 | 12.26 | 15.73 | 12.26 | +0.00 | +0.00 | Moderate → Moderate |
| y_assay/BRCA1_included | raw | single:nt | 6.46 | 6.36 | 6.46 | 6.36 | +0.00 | +0.00 | Moderate → Moderate |
| y_assay/BRCA1_included | calibrated_isotonic | single:nt | 7.28 | 6.21 | 6.56 | 5.54 | -0.72 | -0.67 | Moderate → Moderate |
| y_assay/BRCA1_included | raw | single:cadd | 12.95 | 12.73 | 12.95 | 12.73 | +0.00 | +0.00 | Moderate → Moderate |
| y_assay/BRCA1_included | calibrated_isotonic | single:cadd | 16.07 | 11.45 | 16.07 | 11.45 | +0.00 | +0.00 | Moderate → Moderate |
| y_assay/BRCA1_included | raw | single:gnomad_af | 1.35 | 1.29 | 1.35 | 1.29 | +0.00 | +0.00 | below supporting → below supporting |
| y_assay/BRCA1_included | calibrated_isotonic | single:gnomad_af | 1.21 | 1.12 | 1.21 | 1.12 | +0.00 | +0.00 | below supporting → below supporting |
| y_assay/BRCA1_included | raw | single:phylop | 13.60 | 13.46 | 13.60 | 13.46 | +0.00 | +0.00 | Moderate → Moderate |
| y_assay/BRCA1_included | calibrated_isotonic | single:phylop | 12.80 | 12.67 | 12.80 | 12.67 | +0.00 | +0.00 | Moderate → Moderate |
| y_assay/BRCA1_included | raw | single:phastcons | 6.61 | 5.33 | 6.61 | 5.33 | +0.00 | +0.00 | Moderate → Moderate |
| y_assay/BRCA1_included | calibrated_isotonic | single:phastcons | 5.20 | 4.77 | 5.20 | 4.77 | +0.00 | +0.00 | Moderate → Moderate |
| y_assay/BRCA1_excluded | raw | fusion_M1 | 17.95 | 17.51 | 18.02 | 17.58 | +0.07 | +0.07 | Moderate → Moderate |
| y_assay/BRCA1_excluded | calibrated_isotonic | fusion_M1 | 18.36 | 17.34 | 17.70 | 17.30 | -0.66 | -0.05 | Moderate → Moderate |
| y_assay/BRCA1_excluded | raw | mean_M0b | 16.67 | 16.26 | 16.63 | 16.23 | -0.03 | -0.03 | Moderate → Moderate |
| y_assay/BRCA1_excluded | calibrated_isotonic | mean_M0b | 16.79 | 15.82 | 17.39 | 15.83 | +0.61 | +0.01 | Moderate → Moderate |
| y_assay/BRCA1_excluded | raw | single:spliceai | 16.79 | 15.81 | 16.05 | 15.66 | -0.74 | -0.15 | Moderate → Moderate |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:spliceai | 13.98 | 13.65 | 14.46 | 14.14 | +0.48 | +0.49 | Moderate → Moderate |
| y_assay/BRCA1_excluded | raw | single:pangolin | 16.57 | 16.21 | 16.67 | 16.26 | +0.10 | +0.05 | Moderate → Moderate |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:pangolin | 21.11 | 14.91 | 17.43 | 15.21 | -3.68 | +0.30 | Moderate → Moderate |
| y_assay/BRCA1_excluded | raw | single:alphagenome | 15.43 | 15.05 | 15.43 | 15.05 | +0.00 | +0.00 | Moderate → Moderate |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:alphagenome | 24.34 | 14.03 | 24.34 | 14.03 | +0.00 | +0.00 | Moderate → Moderate |
| y_assay/BRCA1_excluded | raw | single:gpn_msa | 13.36 | 13.03 | 13.36 | 13.03 | +0.00 | +0.00 | Moderate → Moderate |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:gpn_msa | 13.84 | 11.31 | 13.84 | 11.31 | +0.00 | +0.00 | Moderate → Moderate |
| y_assay/BRCA1_excluded | raw | single:nt | 9.56 | 9.33 | 9.56 | 9.33 | +0.00 | +0.00 | Moderate → Moderate |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:nt | 8.39 | 8.23 | 8.74 | 8.16 | +0.35 | -0.07 | Moderate → Moderate |
| y_assay/BRCA1_excluded | raw | single:cadd | 15.52 | 15.49 | 15.52 | 15.49 | +0.00 | +0.00 | Moderate → Moderate |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:cadd | 15.44 | 15.42 | 15.44 | 15.42 | +0.00 | +0.00 | Moderate → Moderate |
| y_assay/BRCA1_excluded | raw | single:gnomad_af | 1.51 | 1.51 | 1.51 | 1.51 | +0.00 | +0.00 | below supporting → below supporting |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:gnomad_af | 1.24 | 1.21 | 1.24 | 1.21 | +0.00 | +0.00 | below supporting → below supporting |
| y_assay/BRCA1_excluded | raw | single:phylop | 13.53 | 13.22 | 13.53 | 13.22 | +0.00 | +0.00 | Moderate → Moderate |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:phylop | 13.64 | 12.40 | 13.64 | 12.40 | +0.00 | +0.00 | Moderate → Moderate |
| y_assay/BRCA1_excluded | raw | single:phastcons | 6.14 | 6.97 | 6.14 | 6.97 | +0.00 | +0.00 | Moderate → Moderate |
| y_assay/BRCA1_excluded | calibrated_isotonic | single:phastcons | 4.83 | 6.20 | 4.83 | 6.20 | +0.00 | +0.00 | Moderate → Moderate |
| y_clinvar/BRCA1_included | raw | fusion_M1 | 20.54 | 19.96 | 20.54 | 19.96 | +0.00 | +0.00 | Strong → Strong |
| y_clinvar/BRCA1_included | calibrated_isotonic | fusion_M1 | 25.88 | 19.86 | 27.73 | 19.86 | +1.85 | +0.00 | Strong → Strong |
| y_clinvar/BRCA1_included | raw | mean_M0b | 20.50 | 19.93 | 20.50 | 19.93 | +0.00 | +0.00 | Strong → Strong |
| y_clinvar/BRCA1_included | calibrated_isotonic | mean_M0b | 21.57 | 19.86 | 29.81 | 19.85 | +8.24 | -0.00 | Strong → Strong |
| y_clinvar/BRCA1_included | raw | single:spliceai | 20.10 | 19.53 | 20.10 | 19.53 | +0.00 | +0.00 | Strong → Strong |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:spliceai | 19.73 | 19.17 | 20.70 | 19.06 | +0.98 | -0.11 | Strong → Strong |
| y_clinvar/BRCA1_included | raw | single:pangolin | 21.64 | 19.93 | 20.50 | 19.93 | -1.14 | +0.00 | Strong → Strong |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:pangolin | 21.25 | 19.57 | 20.13 | 19.60 | -1.12 | +0.03 | Strong → Strong |
| y_clinvar/BRCA1_included | raw | single:alphagenome | 19.76 | 19.21 | 19.76 | 19.21 | +0.00 | +0.00 | Strong → Strong |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:alphagenome | 23.95 | 18.38 | 23.95 | 18.38 | +0.00 | +0.00 | Moderate → Moderate |
| y_clinvar/BRCA1_included | raw | single:gpn_msa | 18.58 | 18.05 | 18.58 | 18.05 | +0.00 | +0.00 | Moderate → Moderate |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:gpn_msa | 18.50 | 17.98 | 18.50 | 17.98 | +0.00 | +0.00 | Moderate → Moderate |
| y_clinvar/BRCA1_included | raw | single:nt | 11.98 | 11.64 | 11.83 | 11.50 | -0.15 | -0.14 | Moderate → Moderate |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:nt | 12.60 | 11.70 | 12.41 | 11.44 | -0.20 | -0.25 | Moderate → Moderate |
| y_clinvar/BRCA1_included | raw | single:cadd | 19.53 | 19.42 | 19.53 | 19.42 | +0.00 | +0.00 | Strong → Strong |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:cadd | 21.62 | 19.11 | 21.62 | 19.11 | +0.00 | +0.00 | Strong → Strong |
| y_clinvar/BRCA1_included | raw | single:gnomad_af | 3.33 | 3.31 | 3.33 | 3.31 | +0.00 | +0.00 | Supporting → Supporting |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:gnomad_af | 2.35 | 2.17 | 2.35 | 2.17 | +0.00 | +0.00 | Supporting → Supporting |
| y_clinvar/BRCA1_included | raw | single:phylop | 19.02 | 18.49 | 19.02 | 18.49 | +0.00 | +0.00 | Moderate → Moderate |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:phylop | 27.64 | 18.38 | 27.64 | 18.38 | +0.00 | +0.00 | Moderate → Moderate |
| y_clinvar/BRCA1_included | raw | single:phastcons | 16.28 | 16.02 | 16.28 | 16.02 | +0.00 | +0.00 | Moderate → Moderate |
| y_clinvar/BRCA1_included | calibrated_isotonic | single:phastcons | 16.28 | 15.82 | 16.28 | 15.82 | +0.00 | +0.00 | Moderate → Moderate |
| y_clinvar/BRCA1_excluded | raw | fusion_M1 | 20.45 | 19.95 | 20.45 | 19.95 | +0.00 | +0.00 | Strong → Strong |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | fusion_M1 | 24.92 | 19.75 | 27.00 | 19.75 | +2.08 | +0.00 | Strong → Strong |
| y_clinvar/BRCA1_excluded | raw | mean_M0b | 20.40 | 19.90 | 20.40 | 19.90 | +0.00 | +0.00 | Strong → Strong |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | mean_M0b | 21.60 | 19.75 | 29.38 | 19.70 | +7.78 | -0.05 | Strong → Strong |
| y_clinvar/BRCA1_excluded | raw | single:spliceai | 21.17 | 19.36 | 19.89 | 19.41 | -1.27 | +0.05 | Strong → Strong |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:spliceai | 20.20 | 18.71 | 19.19 | 18.74 | -1.01 | +0.02 | Strong → Strong |
| y_clinvar/BRCA1_excluded | raw | single:pangolin | 20.40 | 19.90 | 20.40 | 19.90 | +0.00 | +0.00 | Strong → Strong |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:pangolin | 21.06 | 19.26 | 19.64 | 19.16 | -1.42 | -0.10 | Strong → Strong |
| y_clinvar/BRCA1_excluded | raw | single:alphagenome | 19.54 | 19.06 | 19.54 | 19.06 | +0.00 | +0.00 | Strong → Strong |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:alphagenome | 21.24 | 18.13 | 21.24 | 18.13 | +0.00 | +0.00 | Moderate → Moderate |
| y_clinvar/BRCA1_excluded | raw | single:gpn_msa | 18.58 | 18.13 | 18.58 | 18.13 | +0.00 | +0.00 | Moderate → Moderate |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:gpn_msa | 20.72 | 17.71 | 20.72 | 17.71 | +0.00 | +0.00 | Moderate → Moderate |
| y_clinvar/BRCA1_excluded | raw | single:nt | 12.77 | 12.46 | 13.03 | 12.71 | +0.25 | +0.25 | Moderate → Moderate |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:nt | 13.79 | 12.68 | 13.57 | 12.55 | -0.22 | -0.13 | Moderate → Moderate |
| y_clinvar/BRCA1_excluded | raw | single:cadd | 19.70 | 19.57 | 19.70 | 19.57 | +0.00 | +0.00 | Strong → Strong |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:cadd | 19.54 | 19.41 | 19.54 | 19.41 | +0.00 | +0.00 | Strong → Strong |
| y_clinvar/BRCA1_excluded | raw | single:gnomad_af | 2.01 | 2.92 | 2.01 | 2.92 | +0.00 | +0.00 | Supporting → Supporting |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:gnomad_af | 0.80 | 1.11 | 0.80 | 1.11 | +0.00 | +0.00 | below supporting → below supporting |
| y_clinvar/BRCA1_excluded | raw | single:phylop | 19.04 | 18.59 | 19.04 | 18.59 | +0.00 | +0.00 | Moderate → Moderate |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:phylop | 27.69 | 18.57 | 27.69 | 18.57 | +0.00 | +0.00 | Moderate → Moderate |
| y_clinvar/BRCA1_excluded | raw | single:phastcons | 19.64 | 17.27 | 19.64 | 17.27 | +0.00 | +0.00 | Moderate → Moderate |
| y_clinvar/BRCA1_excluded | calibrated_isotonic | single:phastcons | 19.64 | 17.27 | 19.64 | 17.27 | +0.00 | +0.00 | Moderate → Moderate |

At the observed thresholds the largest |ΔLR+| across the 88 evaluable cells is 8.24 and two cells change tier; at interpolated 95.0% the largest |ΔLR+| is 0.67 and **no cell changes tier**. Decomposition for the two Pangolin cells:

| cell | observed Δ | of which score change (interpolated Δ) | of which threshold placement |
|---|---|---|---|
| primary, calibrated | +3.49 | -0.19 (16.14 → 15.94) | +3.69 |
| BRCA1 excluded, calibrated | -3.68 | +0.30 (14.91 → 15.21) | -3.98 |

Read: the scores did not get better or worse in any way the 95% operating point can see; the tier movements are where a calibrated tie group of ten to fifteen variants sits relative to the 5% line. This is the tie artefact the phase-5 docstring already anticipates for isotonic calibration ("any difference found is a tie artefact, not an improvement in evidence"), now visible between two versions of the inputs rather than between raw and calibrated. The raw-score cells, which carry the paper's H4 statements for the fusion, move by at most 0.03 (fusion 18.08 → 18.08; Pangolin 17.06 → 17.04).

## 4. Sign-flip exact test, ΔBrier, primary condition (isotonic / y_assay / BRCA1 included)

| | v1 | v2 |
|---|---|---|
| T_obs (n-weighted per-gene ΔBrier) | 0.008541 | 0.010300 |
| k genes / assignments | 7 / 128 | 7 / 128 |
| exact two-sided p | 0.046875 = **6/128** = 3 × (2/128) | 0.03125 = **4/128** = 2 × (2/128) |
| wild-cluster (Rademacher) p | 0.056 | 0.037 |

With seven clusters the two-sided exact p can only take the values 2j/128; the minimum is 2/128 = 0.0156. v1 sat at 6/128 (the observed statistic is exceeded or matched by six of the 128 sign assignments), v2 at 4/128. The v2 value is 'significant at 0.05 but not at 0.01' exactly as the v1 sentence says; the finest resolution of the test is unchanged. Source: `phase8_sign_flip_exact.csv`, v1 and v2 rows `delta_brier, isotonic/y_assay/BRCA1_included`.

## 5. Test-count reconciliation

`pytest.ini` sets `testpaths = tests phase1/tests`, so a bare `pytest` from the repository root collects both directories (`phase1/tests/test_ipw_reweight.py` holds 11 tests). Earlier figures quoted in this session as '61' and '19 failed / 41 passed / 1 skipped' came from `pytest tests/`, which omits those 11; the root-level numbers below are the ones that correspond to the manuscript.

| checkout | invocation | collected | passed | failed | skipped | skips are |
|---|---|---|---|---|---|---|
| 5b79395 (= v1.0.0-submission + 2; the state the paper describes), fresh git clone, manuscript absent | `pytest` (root) | 66 | 63 | 0 | 3 | the three tests that read the manuscript file |
| same, manuscript present | `pytest` (root) | 66 | 66 | 0 | 0 | — |
| 5b79395 exported with `git archive` (no `.git`), manuscript present — measured this session | `pytest` (root) | 66 | 65 | 0 | 1 | `test_reproduction_coverage` needs `git ls-files` |
| same, manuscript hidden — measured this session | `pytest` (root) | 66 | 62 | 0 | 4 | the three manuscript tests + the coverage test |
| HEAD a5e94f5, v1 default (`FROZEN_VERSION` unset), manuscript present — measured | `pytest` (root) | **72** | 72 | 0 | 0 | — |
| HEAD a5e94f5, v1 default, manuscript hidden — measured | `pytest` (root) | 72 | 69 | 0 | 3 | the three manuscript tests |
| HEAD a5e94f5, `reports/phase1` replaced by the v2 outputs (rsync copy without `.git`) — measured | `pytest` (root) | 72 | 52 | **19** | 1 | coverage test (no `.git` in the copy) |
| HEAD, `pytest tests/` only (what `~/work/tools/run_tests_on_v2.sh` ran) | `pytest tests/` | 61 | 41 (v2) / 61 (v1) | 19 (v2) / 0 | 1 / 0 | — |

Where the six P2 structural tests (`tests/test_frozen_matrix_v2.py`) fall: they read the two parquet files, not the report tables, so they pass on **both** sides — 6 of the 72 passed under v1 defaults and 6 of the 52 passed with the v2 outputs swapped in. The five `tests/test_frozen_matrix.py` tests (v1 hash and shape) likewise pass on both sides. The 19 failures under v2 are 11 in `test_phase8_robustness.py` and 8 in `test_reported_numbers.py`, every one an expected value pinned to a published v1 number (one, `test_tp53_external_validation`, fails on the row name `fusion_8feat`, which the ten-predictor v2 run writes as `fusion_10feat`).

The manuscript's 'Sixty-three automated tests pass from a fresh clone … sixty-six collected; the three that read the manuscript file skip in its absence … where the manuscript file is present, all sixty-six pass' therefore corresponds to the first two rows: commit 5b79395, a git clone (not an archive extract), run from the repository root. At HEAD the same sentence would read 69 / 72 / 3 and 72, because the six v2 structural tests were added; nothing else in the suite changed.
