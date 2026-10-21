# Rebaseline in the reference image

Written by `scripts/stage_B6/rebaseline_report.py` from the per-file report of `scripts/reproduce_evidence.py --verify` (`reproduce_report.txt`), after every stage ran in the reference image (`container/run.sh`: inputs fetched online, stages and `--verify` with no network), against the committed outputs of commit e12220c; do not edit by hand. The comparator alone grades the files; the differences below are computed here and listed, not judged.

## The comparator's verdict

Quoted as it wrote them:

> NOT REPRODUCED: every file not marked identical
>
> differs: 2; differs (figure): 13; differs (numerically equal): 21; identical: 429; identical (build time aside): 6

## Files not identical

| status | file |
|---|---|
| differs | `phase1/data/evidence/external/tp53_splice_scored_12nt.manifest.json` |
| differs (numerically equal) | `phase1/data/evidence/external/tp53_splice_scored_12nt.parquet` |
| differs (numerically equal) | `phase1/reports/evidence/evid_delta.csv` |
| differs (numerically equal) | `phase1/reports/evidence/evidence_thresholds.csv` |
| differs (numerically equal) | `phase1/reports/evidence/evidence_thresholds_logo_folds.csv` |
| differs (figure) | `phase1/reports/evidence/fig_data/fig1_draft.png` |
| differs (numerically equal) | `phase1/reports/evidence/fig_data/fig1_walker_cutpoints.csv` |
| differs (figure) | `phase1/reports/evidence/fig_data/fig2_draft.png` |
| differs (numerically equal) | `phase1/reports/evidence/fig_data/fig2_interval_lr.csv` |
| differs (figure) | `phase1/reports/evidence/fig_data/fig3_draft.png` |
| differs (numerically equal) | `phase1/reports/evidence/fig_data/fig3_evidence_by_territory.csv` |
| differs (figure) | `phase1/reports/evidence/fig_data/fig4_draft.png` |
| differs (figure) | `phase1/reports/evidence/figures/figure1.png` |
| differs (figure) | `phase1/reports/evidence/figures/figure2.png` |
| differs (figure) | `phase1/reports/evidence/figures/figure3.png` |
| differs (figure) | `phase1/reports/evidence/figures/figure4.png` |
| differs (figure) | `phase1/reports/evidence/figures/figureS1.png` |
| differs (figure) | `phase1/reports/evidence/figures/figureS2.png` |
| differs (figure) | `phase1/reports/evidence/figures/figureS3.png` |
| differs (figure) | `phase1/reports/evidence/figures/figure_practical.png` |
| differs (figure) | `phase1/reports/evidence/figures/figure_practical_slides.png` |
| differs (numerically equal) | `phase1/reports/evidence/fusion_fold_coefficients.csv` |
| differs (numerically equal) | `phase1/reports/evidence/interval_lr_fusion_enet_all_1_50.csv` |
| differs (numerically equal) | `phase1/reports/evidence/interval_lr_fusion_enet_pm12.csv` |
| differs (numerically equal) | `phase1/reports/evidence/interval_lr_fusion_enet_s11_50.csv` |
| differs (numerically equal) | `phase1/reports/evidence/interval_lr_fusion_enet_s3_10.csv` |
| differs (numerically equal) | `phase1/reports/evidence/interval_lr_fusion_enet_s3_50.csv` |
| differs (numerically equal) | `phase1/reports/evidence/local_lr_monotonicity.txt` |
| differs (numerically equal) | `phase1/reports/evidence/practical_curve_auc.csv` |
| differs (numerically equal) | `phase1/reports/evidence/practical_metrics_by_threshold.csv` |
| differs (numerically equal) | `phase1/reports/evidence/practical_metrics_logo_by_gene.csv` |
| differs (numerically equal) | `phase1/reports/evidence/practical_pr_curves.csv` |
| differs (numerically equal) | `phase1/reports/evidence/practical_roc_curves.csv` |
| differs | `phase1/reports/evidence/supplement/supplement_tables_manifest.csv` |
| differs (numerically equal) | `phase1/reports/evidence/territory_metrics.csv` |
| differs (numerically equal) | `phase1/reports/evidence/territory_metrics_arms_noBRCA1.csv` |

## How far the numerically equal files moved

Per file, the columns (for JSON, the keys; `[]` stands for any list position) whose numbers moved, the cells that moved, and the largest absolute and relative difference between the committed file and the image's (relative: the difference over the larger magnitude).

| file | column | cells moved | largest absolute | largest relative |
|---|---|---|---|---|
| `phase1/data/evidence/external/tp53_splice_scored_12nt.parquet` | `nt` | 5 | 8.88e-16 | 2.32e-16 |
| `phase1/reports/evidence/evid_delta.csv` | `evid_set_8853` | 6 | 1.11e-16 | 2.21e-16 |
| `phase1/reports/evidence/evidence_thresholds.csv` | `bp4_threshold_insample` | 3 | 2.22e-16 | 5.06e-16 |
| `phase1/reports/evidence/evidence_thresholds.csv` | `bp4_threshold_variant_boot` | 4 | 2.22e-16 | 5.06e-16 |
| `phase1/reports/evidence/evidence_thresholds.csv` | `pp3_logo_threshold_median` | 5 | 2.22e-16 | 3.11e-16 |
| `phase1/reports/evidence/evidence_thresholds.csv` | `pp3_logo_threshold_spread` | 2 | 2.22e-16 | 9.47e-15 |
| `phase1/reports/evidence/evidence_thresholds.csv` | `pp3_threshold_insample` | 2 | 2.22e-16 | 3.12e-16 |
| `phase1/reports/evidence/evidence_thresholds.csv` | `pp3_threshold_variant_boot` | 2 | 1.11e-16 | 1.77e-16 |
| `phase1/reports/evidence/evidence_thresholds_logo_folds.csv` | `threshold_from_6_genes` | 23 | 2.22e-16 | 3.37e-16 |
| `phase1/reports/evidence/fig_data/fig1_walker_cutpoints.csv` | `lr_pp3` | 1 | 8.88e-16 | 1.73e-16 |
| `phase1/reports/evidence/fig_data/fig1_walker_cutpoints.csv` | `lr_pp3_hi` | 1 | 3.55e-15 | 1.52e-16 |
| `phase1/reports/evidence/fig_data/fig1_walker_cutpoints.csv` | `lr_pp3_lo` | 1 | 8.88e-16 | 1.48e-16 |
| `phase1/reports/evidence/fig_data/fig2_interval_lr.csv` | `hi` | 2972 | 1.14e-13 | 3.83e-16 |
| `phase1/reports/evidence/fig_data/fig2_interval_lr.csv` | `lo` | 1379 | 1.42e-14 | 3.91e-16 |
| `phase1/reports/evidence/fig_data/fig2_interval_lr.csv` | `local_lr` | 1271 | 1.42e-14 | 3.8e-16 |
| `phase1/reports/evidence/fig_data/fig2_interval_lr.csv` | `score` | 8094 | 8.88e-16 | 9.02e-16 |
| `phase1/reports/evidence/fig_data/fig3_evidence_by_territory.csv` | `auroc` | 1 | 1.11e-16 | 1.61e-16 |
| `phase1/reports/evidence/fig_data/fig3_evidence_by_territory.csv` | `auroc_lo` | 1 | 1.11e-16 | 1.75e-16 |
| `phase1/reports/evidence/fusion_fold_coefficients.csv` | `coef_heldout_BAP1` | 5 | 1.11e-16 | 1.16e-15 |
| `phase1/reports/evidence/fusion_fold_coefficients.csv` | `coef_heldout_BARD1` | 6 | 2.5e-16 | 5.11e-15 |
| `phase1/reports/evidence/fusion_fold_coefficients.csv` | `coef_heldout_BRCA1` | 6 | 1.11e-16 | 1.11e-15 |
| `phase1/reports/evidence/fusion_fold_coefficients.csv` | `coef_heldout_BRCA2` | 5 | 1.67e-16 | 1.41e-15 |
| `phase1/reports/evidence/fusion_fold_coefficients.csv` | `coef_heldout_PALB2` | 5 | 1.11e-16 | 3.76e-15 |
| `phase1/reports/evidence/fusion_fold_coefficients.csv` | `coef_heldout_RAD51C` | 6 | 3.05e-16 | 5.82e-15 |
| `phase1/reports/evidence/fusion_fold_coefficients.csv` | `coef_heldout_VHL` | 6 | 1.11e-16 | 4.08e-15 |
| `phase1/reports/evidence/fusion_fold_coefficients.csv` | `max` | 6 | 2.5e-16 | 1.13e-15 |
| `phase1/reports/evidence/fusion_fold_coefficients.csv` | `mean` | 5 | 5.55e-17 | 8.93e-16 |
| `phase1/reports/evidence/fusion_fold_coefficients.csv` | `mean_abs` | 5 | 5.55e-17 | 8.93e-16 |
| `phase1/reports/evidence/fusion_fold_coefficients.csv` | `min` | 4 | 6.59e-17 | 5.11e-15 |
| `phase1/reports/evidence/fusion_fold_coefficients.csv` | `sd` | 6 | 6.94e-17 | 4.28e-15 |
| `phase1/reports/evidence/interval_lr_fusion_enet_all_1_50.csv` | `eps` | 5374 | 4.44e-16 | 1.63e-13 |
| `phase1/reports/evidence/interval_lr_fusion_enet_all_1_50.csv` | `score` | 4064 | 2.22e-16 | 6.65e-16 |
| `phase1/reports/evidence/interval_lr_fusion_enet_pm12.csv` | `eps` | 410 | 4.44e-16 | 1.39e-13 |
| `phase1/reports/evidence/interval_lr_fusion_enet_pm12.csv` | `score` | 259 | 2.22e-16 | 3.55e-16 |
| `phase1/reports/evidence/interval_lr_fusion_enet_s11_50.csv` | `eps` | 3126 | 4.44e-16 | 1.21e-13 |
| `phase1/reports/evidence/interval_lr_fusion_enet_s11_50.csv` | `score` | 2341 | 2.22e-16 | 6.65e-16 |
| `phase1/reports/evidence/interval_lr_fusion_enet_s3_10.csv` | `eps` | 1984 | 4.44e-16 | 6.87e-14 |
| `phase1/reports/evidence/interval_lr_fusion_enet_s3_10.csv` | `score` | 1464 | 2.22e-16 | 6.29e-16 |
| `phase1/reports/evidence/interval_lr_fusion_enet_s3_50.csv` | `eps` | 4995 | 4.44e-16 | 1.63e-13 |
| `phase1/reports/evidence/interval_lr_fusion_enet_s3_50.csv` | `score` | 3805 | 2.22e-16 | 6.65e-16 |
| `phase1/reports/evidence/local_lr_monotonicity.txt` | `(words of the text)` | 6 | 5.68e-14 | 3.15e-16 |
| `phase1/reports/evidence/practical_curve_auc.csv` | `auroc_gene_pooled` | 1 | 1.11e-16 | 1.61e-16 |
| `phase1/reports/evidence/practical_curve_auc.csv` | `auroc_gene_pooled_lo` | 1 | 1.11e-16 | 1.75e-16 |
| `phase1/reports/evidence/practical_curve_auc.csv` | `auroc_pooled_minus_gene_pooled` | 1 | 1.11e-16 | 1.47e-14 |
| `phase1/reports/evidence/practical_curve_auc.csv` | `threshold_at_specificity_0.90` | 2 | 1.11e-16 | 1.9e-16 |
| `phase1/reports/evidence/practical_metrics_by_threshold.csv` | `threshold` | 3 | 2.22e-16 | 5.06e-16 |
| `phase1/reports/evidence/practical_metrics_by_threshold.csv` | `threshold_fold_max` | 1 | 2.22e-16 | 3.09e-16 |
| `phase1/reports/evidence/practical_metrics_by_threshold.csv` | `threshold_fold_min` | 2 | 2.22e-16 | 3.09e-16 |
| `phase1/reports/evidence/practical_metrics_logo_by_gene.csv` | `threshold_from_6_genes` | 16 | 2.22e-16 | 3.09e-16 |
| `phase1/reports/evidence/practical_pr_curves.csv` | `threshold` | 7610 | 2.22e-16 | 6.65e-16 |
| `phase1/reports/evidence/practical_roc_curves.csv` | `threshold` | 7610 | 2.22e-16 | 6.65e-16 |
| `phase1/reports/evidence/territory_metrics.csv` | `auroc` | 1 | 1.11e-16 | 1.61e-16 |
| `phase1/reports/evidence/territory_metrics.csv` | `auroc_i2_pct` | 2 | 1.42e-14 | 1.9e-16 |
| `phase1/reports/evidence/territory_metrics.csv` | `auroc_lo` | 1 | 1.11e-16 | 1.75e-16 |
| `phase1/reports/evidence/territory_metrics.csv` | `e3_threshold_at_best_tier` | 4 | 2.22e-16 | 3.12e-16 |
| `phase1/reports/evidence/territory_metrics_arms_noBRCA1.csv` | `auroc_lo` | 1 | 1.11e-16 | 1.59e-16 |

Largest relative difference over all of them: 1.63e-13.

## Files that differ in more than their numbers' last digits

The comparator called these "differs". Their changed lines, HEAD's first (`-`), then the image's (`+`); a sha256 that is the hash of a file above is named after it.

`phase1/data/evidence/external/tp53_splice_scored_12nt.manifest.json`, 12 changed lines:

```
-      "n_identical": 24,
-      "max_abs_diff": 0.0,
-      "differing": []
+      "n_identical": 23,
+      "max_abs_diff": 4.440892098500626e-16,
+      "differing": [
+        "17-7673696-C-A"
+      ]
-  "sha256": "<sha256 of phase1/data/evidence/external/tp53_splice_scored_12nt.parquet>",
+  "sha256": "<sha256 of phase1/data/evidence/external/tp53_splice_scored_12nt.parquet>",
-  "parquet_sha256": "<sha256 of phase1/data/evidence/external/tp53_splice_scored_12nt.parquet>"
+  "parquet_sha256": "<sha256 of phase1/data/evidence/external/tp53_splice_scored_12nt.parquet>"
```

`phase1/reports/evidence/supplement/supplement_tables_manifest.csv`, 4 changed cells:

```
row S1, Source sha256: eb3c0af7abf1...; b178232474d1...; bfbc8bab6d93...; n/a (code); n/a (code); cf51a401813a...; a5f1ccdd9d52...; <sha256 of phase1/data/evidence/external/tp53_splice_scored_12nt.parquet>; <sha256 of phase1/data/evidence/external/tp53_splice_scored_12nt.manifest.json>; b0fc17b4ffad... -> eb3c0af7abf1...; b178232474d1...; bfbc8bab6d93...; n/a (code); n/a (code); cf51a401813a...; a ...
row S6, Source sha256: <sha256 of phase1/reports/evidence/evidence_thresholds.csv>; b8ba8ae0d90b...; <sha256 of phase1/reports/evidence/evidence_thresholds_logo_folds.csv> -> <sha256 of phase1/reports/evidence/evidence_thresholds.csv>; b8ba8ae0d90b...; <sha256 of phase1/reports/evidence/evidence_thresholds_logo_folds.csv>
row S7, Source sha256: <sha256 of phase1/reports/evidence/evidence_thresholds.csv> -> <sha256 of phase1/reports/evidence/evidence_thresholds.csv>
row S10, Source sha256: b99e15e13428...; <sha256 of phase1/reports/evidence/territory_metrics_arms_noBRCA1.csv>; afdcb4e2f57b... -> b99e15e13428...; <sha256 of phase1/reports/evidence/territory_metrics_arms_noBRCA1.csv>; afdcb4e2f57b...
```

## The image and the machine

- Image: `ghcr.io/oneone00-11/variant-fm-benchmark-reference@sha256:0dc03436cd94cf46effe25160efa3015ab9808852179d58908c72e16e0537776` (image ID, the config digest: `sha256:8ecf15fc8944ac940876fa1a9a5820767753b6342b40dc4c6ec064a70444c78a`), built at commit 9b15c34 by https://github.com/oneone00-11/variant-fm-benchmark/actions/runs/39621528080 (`container/ENVIRONMENT.json`, section image).
- CPU: VirtualApple @ 2.50GHz (VirtualApple), 6 processors; flags: aes apic avx avx2 bmi1 bmi2 clflush cmov cpuid cx16 cx8 de f16c fma fpu fxsr lahf_lm lm mmx movbe nopl nx pat pclmulqdq pni popcnt pse36 rdrand rep_good sep sse sse2 sse4_1 sse4_2 ssse3 syscall tsc (`container/cpu_probe/mac-rosetta-pinned.json`).
- OpenBLAS kernels: libscipy_openblas 0.3.27: Nehalem; libopenblas 0.3.27: Nehalem.
- numpy dispatch in use: SSSE3 SSE41 POPCNT SSE42.
- Settings in force: OPENBLAS_CORETYPE=Nehalem, NPY_DISABLE_CPU_FEATURES=AVX F16C FMA3 AVX2 AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL, OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1, MKL_NUM_THREADS=1, VECLIB_MAXIMUM_THREADS=1, NUMEXPR_NUM_THREADS=1, PYTHONHASHSEED=0, TZ=UTC, LC_ALL=C.UTF-8.
