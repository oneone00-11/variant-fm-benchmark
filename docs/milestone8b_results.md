# Milestone 8b — Phase 3 re-run with Nucleotide Transformer (NT)

Re-ran the full Milestone 8 evaluation with NT added (cloud-GPU scored,
masked REF/ALT log-likelihood ratio, oriented higher=pathogenic). Same methods
and main/auxiliary hierarchy as M8 — only NT is new. NT covers **21,410/21,410**.
Evo2/ESM still NA (pending). This documents the **changes vs M8**.

## (a) NT on splice — meta ρ and where it ranks

Splice subset (N=1,781), meta Spearman ρ vs functional gold standard:

| rank | model | meta ρ | 95% CI | I² |
|---|---|--:|---|--:|
| 1 | Pangolin | 0.761 | 0.729–0.790 | 54% |
| 2 | SpliceAI | 0.753 | 0.718–0.784 | 58% |
| 3 | AlphaGenome | 0.750 | 0.714–0.782 | 59% |
| 4 | CADD | 0.694 | 0.654–0.731 | 53% |
| 5 | GPN-MSA | 0.666 | 0.639–0.691 | 0% |
| 6 | phyloP | 0.639 | 0.607–0.669 | 15% |
| 7 | phastCons | 0.612 | 0.563–0.657 | 57% |
| **8** | **NT** | **0.494** | **0.374–0.599** | **89%** |
| 9 | gnomAD AF | 0.186 | 0.069–0.298 | 44% |
| — | AlphaMissense | — (n=1) | — | — |

**NT lands 8th** — below every other DNA-reading model and conservation, but
clearly above the allele-frequency baseline. Its **I²=89%** means NT's splice
performance is highly variable across the 7 genes (least consistent model).

## (b) The DNA-language-model line (NT + GPN-MSA) vs splice specialists

Pairwise paired-bootstrap Δ(mean per-gene ρ) on the splice subset:

| comparison | Δρ | 95% CI | significant |
|---|--:|---|---|
| GPN-MSA − NT | **+0.183** | +0.140, +0.221 | **yes** |
| SpliceAI − NT | +0.262 | +0.222, +0.301 | yes |
| Pangolin − NT | +0.268 | +0.230, +0.307 | yes |
| AlphaGenome − NT | +0.262 | +0.218, +0.305 | yes |
| CADD − NT | +0.207 | +0.167, +0.248 | yes |
| Pangolin − GPN-MSA | +0.085 | +0.061, +0.111 | yes |

**The two DNA language models are not equal:** the alignment/MSA-based **GPN-MSA
significantly outperforms the single-sequence NT** on splice (Δ +0.183, CI excludes
0). The ordering on splice is:

> supervised splice specialists (≈0.75) > CADD (0.69) ≈ GPN-MSA (0.67) > NT (0.49) ≫ gnomAD (0.19)

So the zero-shot DNA-LM line **can read splice** (both NT and GPN-MSA beat the
frequency baseline), but neither zero-shot DNA-LM yet matches the supervised
splice models; and within DNA-LMs the MSA-based model is much stronger than the
single-sequence one. Both DNA-LMs trail the splice specialists significantly
(GPN-MSA gap ≈ −0.08; NT gap ≈ −0.26).

By region, NT mirrors GPN-MSA's pattern but weaker: splice 0.49, intronic 0.03,
coding 0.07, missense ≈ 0.00 — i.e. NT's masked-LLR captures **splice-proximal**
signal but essentially **no missense** signal (unlike AlphaMissense 0.47).

## (c) Does adding NT change the main conclusion? No.

- **Main result unchanged.** The top tier on splice is still SpliceAI ≈ Pangolin ≈
  AlphaGenome (≈0.75, mutually non-significant); NT is an **additional, weaker**
  DNA-LM data point (0.49), not a new leader. The continuous-functional,
  excl-BRCA1 main conclusion stands exactly as in M8.
- **Auxiliary (ClinVar AUROC, splice):** NT 0.872 [0.848–0.895] excl-BRCA1 vs
  0.828 [0.802–0.854] incl-BRCA1 — below all other splice-aware models (≈0.98–1.0);
  the BRCA1 sensitivity is small and does not affect ranking.
- **Clinical (exploratory):** NT VUS-reclassification AUROC 0.527 [0.509–0.545]
  (≈chance — most VUS are coding, where NT is weak); sensitivity@95%-spec 0.538
  [0.513–0.565] (above splice specialists on the all-variant set, below GPN-MSA/CADD).
  Exploratory; small effect.

## Updated deliverables (overwrite M8 same-named outputs)
- `results/tables/`: `spearman_meta.tsv`, `spearman_by_gene.tsv`, `auroc_clinvar.tsv`,
  `pairwise_splice.tsv` (now incl. NT), `vus_reclass.tsv`, `sens_at_95spec.tsv`,
  `full_model_coverage.tsv` (NT = 21,410)
- `results/figures/`: `fig1_splice_performance.png` (incl. NT), forest, heatmap
- `data/processed/score_matrix_final.tsv` now carries the populated NT column.

## Bottom line (vs M8)
Adding NT **strengthens, not changes** the story: it is a second DNA language model
that demonstrably reads splice signal (ρ=0.49 ≫ frequency), but in zero-shot it is
the weakest DNA-reading model and is **significantly outperformed by the MSA-based
GPN-MSA (Δ +0.18)** and by the supervised splice specialists (Δ +0.26). The headline
— DNA/splice models cover and predict the splice gap that AlphaMissense cannot
(n=1/1781), with the supervised splice specialists on top — is unchanged.

## Pending
- Evo2, ESM — still NA (need GPU). Adding them would re-run this same pipeline.
