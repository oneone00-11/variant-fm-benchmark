# Milestone 8 — Results (Phase 3)

Final performance evaluation on `score_matrix_final.tsv` (7 genes, 21,410 SNVs).
All models oriented so **higher = more pathogenic**; functional pathogenicity =
`-functional_score`, so a **positive Spearman = model agrees with pathogenicity**.
Every estimate carries a 95% CI. NT/Evo2/ESM are **pending** (no GPU) and excluded.

**Main result = correlation with the continuous, independent functional gold
standard (per gene, then random-effects meta).** ClinVar binary metrics are
**auxiliary** (and circular for ClinVar-trained models — see (c)).

---

## (a) Which models are reliable on splice variants, and by how much

Splice subset (core+region; N=1,781), **meta Spearman ρ vs functional score**
(DerSimonian-Laird across 7 genes):

| model | meta ρ | 95% CI | I² |
|---|--:|---|--:|
| Pangolin | **0.761** | 0.729–0.790 | 54% |
| SpliceAI | **0.753** | 0.718–0.784 | 58% |
| AlphaGenome | **0.750** | 0.714–0.782 | 59% |
| CADD | 0.694 | 0.654–0.731 | 53% |
| GPN-MSA | 0.666 | 0.639–0.691 | 0% |
| phyloP100way | 0.639 | 0.607–0.669 | 15% |
| phastCons100way | 0.612 | 0.563–0.657 | 57% |
| gnomAD AF | 0.186 | 0.069–0.298 | 44% |
| **AlphaMissense** | **— (n=1)** | — | — |

**Pairwise (paired bootstrap on the same splice variants, Δ mean per-gene ρ):**
SpliceAI ≈ Pangolin ≈ AlphaGenome — **all three differences non-significant**
(SpliceAI−Pangolin −0.006 [−0.020,+0.010]; SpliceAI−AlphaGenome 0.000
[−0.019,+0.019]; Pangolin−AlphaGenome +0.006 [−0.015,+0.027]). All three
**significantly** beat GPN-MSA (Δ≈+0.08, CI excludes 0) and CADD (Δ≈+0.06).
GPN-MSA significantly below CADD (−0.024 [−0.049,−0.000]).

**Read:** the three splice specialists (SpliceAI, Pangolin, AlphaGenome) are the
top tier and **statistically indistinguishable** from each other on splice — we do
**not** claim a single winner. The DNA language model GPN-MSA and the
genome-wide CADD are clearly behind but still strong (ρ≈0.67–0.69). Conservation
moderate; allele frequency weak. (`spearman_meta.tsv`, `pairwise_splice.tsv`.)

## (b) The coverage gap — and whether DNA models fill it

**AlphaMissense (protein/missense) scores 1 of 1,781 splice variants** — it is
structurally blind to splice/non-coding variation. In contrast, **all four
DNA/splice-reading models score all 1,781** splice variants, at ρ≈0.67–0.76.

So on the project's banner region, the gap left by the protein model is **fully
covered** by the DNA/splice models, with strong agreement to the functional gold
standard. By region (meta ρ), the picture inverts as expected:

| region | best DNA/splice | AlphaMissense |
|---|--:|--:|
| splice | 0.75 (SpliceAI/Pangolin/AlphaGenome) | n=1 (no coverage) |
| coding | 0.13 (splice models ≈ noise) | **0.471** |
| missense | 0.10 | **0.471** |

i.e. splice models and the protein model are **complementary**, each strong only
in its own regime; CADD and GPN-MSA are the most consistent genome-wide
(coding ρ≈0.41–0.46, splice ρ≈0.67–0.69). This complementarity is the core
finding. (Figures: `fig1_splice_performance.png`, `fig3_region_heatmap.png`.)

## (c) Circularity sensitivity (with vs without BRCA1)

ClinVar binary AUROC, **splice** subset, both versions are essentially identical:

| model | incl BRCA1 | excl BRCA1 |
|---|--:|--:|
| Pangolin | 0.9993 | 0.9992 |
| SpliceAI | 0.9975 | 0.9974 |
| AlphaGenome | 0.9970 | 0.9972 |
| GPN-MSA | 0.9764 | 0.9774 |

Removing BRCA1 (the only `circularity_flag=TRUE` gene) **does not change the
splice conclusions**. Two caveats on the ClinVar binary metric (why it is
auxiliary, not main):
1. On splice the binary task is near-trivial (canonical sites are P/LP, deep
   intronic B/LB) → AUROC saturates near 1.0 for every splice-aware model; it does
   **not** discriminate them. The **continuous functional correlation does** —
   hence it is the main result.
2. **CADD AUROC ≈ 0.995 on the full clean_PB set is inflated** — CADD (and
   AlphaMissense) are trained on ClinVar-like labels, so ClinVar AUROC is circular
   at the *model-training* level, not just for BRCA1. The functional gold standard
   avoids this entirely. (`auroc_clinvar.tsv`.)

---

## Auxiliary / clinical (exploratory — do not over-read)

- **Sensitivity @ 95% specificity** (clean_PB, excl BRCA1, all regions): CADD 0.985
  [0.978–0.991] (inflated, see above), AlphaMissense 0.90 (missense, n_pos=301),
  GPN-MSA 0.831 [0.801–0.852], conservation 0.65–0.69; splice models 0.39–0.48
  (low because evaluated genome-wide — they fire only on splice). (`sens_at_95spec.tsv`.)
- **VUS reclassification (EXPLORATORY)**: per-gene functional threshold (Youden-J
  on clean_PB) labelled VUS functionally abnormal/normal; model AUROC to separate
  them — AlphaMissense 0.841 [0.829–0.853] (missense VUS), GPN-MSA 0.763, CADD
  0.781, AlphaGenome 0.641, SpliceAI/Pangolin ~0.61 (most VUS are coding, so splice
  models help little); gnomAD ~0.50. n_abnormal≈1,267 / n_normal≈6,584. The
  functional threshold is derived from ClinVar-labelled variants, so this is
  **exploratory**, not a primary claim. (`vus_reclass.tsv`.)

## Statistics & robustness
- Per-gene Spearman with Fisher-z CIs; cross-gene pooling by DerSimonian-Laird
  random effects (I², τ² reported) — never mixing gene-specific functional scales.
- AUROC/AUPRC: stratified bootstrap (2000×, resampled within gene), percentile CIs.
- Pairwise: paired bootstrap on the same variant set → Δρ + CI (no point-estimate-
  only claims).
- Missingness: each cell uses only variants the model scored (n reported);
  gnomAD-absent treated as **missing, not AF=0**. No imputation/extrapolation on NA.
- `splice_core` binary AUROC skipped (0 benign — biological, per M3); it appears in
  the continuous result only.

## Deliverables
- `results/tables/`: `spearman_by_gene.tsv`, `spearman_meta.tsv`,
  `auroc_clinvar.tsv`, `pairwise_splice.tsv`, `vus_reclass.tsv`, `sens_at_95spec.tsv`
- `results/figures/`: `fig1_splice_performance.png`, `fig2_forest_spliceai_splice.png`,
  `fig3_region_heatmap.png`

## Bottom line
On the splice variants where an independent functional gold standard exists,
**SpliceAI, Pangolin and AlphaGenome are the top, mutually-tied tier (ρ≈0.75)**;
GPN-MSA and CADD follow (ρ≈0.67–0.69). The protein model AlphaMissense **cannot
score splice variants at all (n=1/1781)** — the DNA/splice models fill exactly that
gap. Conclusions hold with or without BRCA1 and rest on the continuous functional
standard, not the (circular, saturated) ClinVar binary metric.

## Pending (future)
- Nucleotide Transformer, Evo2, ESM — need GPU; would add more DNA-LM / protein-LM
  points (NA in the matrix now).
