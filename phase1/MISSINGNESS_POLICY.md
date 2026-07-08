# Phase 1 — missingness policy

This records the per-predictor coverage facts on the frozen matrix (`frozen_matrix_v1.parquet`)
and the handling convention Phase 2 will follow. No imputation is baked into the frozen matrix;
all imputation is deferred to Phase 2 and fit **inside each training fold** to prevent leakage.

## (a) Coverage facts

Source: `reports/phase1/coverage_by_region.csv` (fraction of variants with a non-NA score,
by `analysis_region`). Splice-set regions are `splice_core`, `splice_region`, `splice_exon_edge`.

| predictor | splice_core | splice_region | splice_exon_edge | intronic | missense | notes |
|---|---|---|---|---|---|---|
| spliceai | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | complete |
| pangolin | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | complete |
| alphagenome | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | complete |
| gpn_msa | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | complete |
| nt | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | complete |
| phylop | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | complete |
| phastcons | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | complete |
| cadd | 0.935 | 0.917 | 1.00 | 0.920 | 0.932 | ~92–93% on splice |
| gnomad_af | 0.232 | 0.328 | 0.00 | 0.391 | 0.365 | sparse (only common variants) |
| **alphamissense** | **0.00** | **0.00** | **0.143** | **0.00** | 0.994 | **≈0 on splice — protein model is blind here** |

**Headline fact:** **AlphaMissense coverage on the splice region is ≈0** — 0.00 on both `splice_core`
and `splice_region`, and only 1/7 `splice_exon_edge` variants (0.143) because those are exon-side
(missense-scorable). Over the full 1,781-variant splice set this is ~0.06% (1 variant). This is the
central motivation of the study: the protein-level predictor cannot score splice/intronic variants,
which is exactly where the fusion must rely on the splice-aware and evolution/alignment predictors.

Secondary facts: `gnomad_af` is sparse everywhere (~23–39% on splice; only variants observed in
gnomAD get an allele frequency — absence is informative, not random). `cadd` is ~92–93% complete on
splice. The seven splice-aware / DNA-LM / conservation predictors are 100% complete on every region.

## (b) Phase 2 handling convention

Missingness is a first-class feature, never silently dropped. Handling depends on the model family:

- **Tree models (LightGBM / XGBoost — the M2 workhorse):** use **native NA handling**. The learner
  sends missing values down its own default branch; no imputation. The `*_isna` indicator columns are
  available if a split on missingness is useful but are not required.
- **Linear / MLP models (M1 elastic net, M3 shallow MLP):** use the **`*_isna` indicator columns
  (already frozen in) + median imputation fit INSIDE each training fold**. The median is computed on
  the training genes of that fold only and applied to the held-out gene, so no cross-gene statistic
  leaks. The indicator lets the model learn a distinct effect for "missing".

- **No global imputation.** No median/mean is computed over the whole matrix and frozen in. Every
  fold-dependent statistic (imputation values, rank-normalisation) is fit within the fold in Phase 2.
  The frozen matrix carries only raw scores, labels, and `*_isna` flags.

- **No variant is dropped for missing a non-required score.** A variant missing AlphaMissense or
  gnomAD still participates; its other predictors and the continuous `func_pathogenicity` target
  remain valid.

### Rationale for the `*_isna` flags being in the frozen matrix

They are a pure per-row property (no cross-gene statistic), so freezing them in cannot leak. They let
the linear/MLP branch model informative missingness (especially gnomAD absence ⇒ rare ⇒ weakly
enriched for pathogenic) without any fold-dependent computation at freeze time.
