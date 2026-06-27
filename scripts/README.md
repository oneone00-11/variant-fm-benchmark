# scripts/ — pipeline & analysis code, by milestone

Scripts keep a **flat layout** with numeric prefixes that encode phase order
(34 of them do `from config import …` / `from phase3_lib import …`, so they must
share one directory). The numbering, not subfolders, is the organization.

**To reproduce the results you only need the Phase-3 block** (100–103), wrapped by
`reproduce_main.py`. Everything ≤ 96 is the data-assembly / scoring provenance and
requires the raw inputs (and, for some, a GPU or API key) that are **not** in the repo.

## Shared libraries
| File | Role |
|---|---|
| `config.py` | repo-relative paths, 7-gene list, ClinVar CLNSIG/CLNREVSTAT maps |
| `phase3_lib.py` | model orientation (`ORIENT`), Fisher-z CIs, DerSimonian–Laird meta, stratified bootstrap |
| `hgvs_region.py` | HGVS c.-token parsing → region/intron-offset class |
| `mutalyzer_map.py` | transcript-HGVS → GRCh38 mapping (splice/intronic), cached |

## Milestone 1 — ClinVar inventory
`01_download_clinvar.py` · `02_clinvar_inventory.py`

## Milestone 2 — Functional gold standard (MaveDB SGE) + coordinate mapping
Probes `10`–`13`; catalog/build `20_catalog_mavedb.py`, `21_build_functional_master.py`;
mapping/validation `22`–`28`; `30_proteingym.py`; `40_coverage_matrix.py`

## Milestone 3 — Analysis-ready set + power analysis
`50_inspect_for_merge.py` · `51_clean_and_merge.py` · `52_precision.py`

## Milestone 4 — Panel expansion (VHL, BAP1)
`60_catalog_candidates.py` · `60b_check_assays.py` · `61_build_candidate_master.py`
· `62_clinvar_candidates.py` · `63_candidate_coverage.py`

## Milestone 5 — Precomputed scores (CADD, AlphaMissense, phyloP/phastCons, gnomAD)
`70_merge_v2.py` (7-gene merge) · `71_probe_scores.py` · `72_gnomad.py` ·
`73_conservation.py` · `74_alphamissense.py` · `75_cadd.py` · `76_assemble_score_matrix.py`

## Milestone 6 — Local SpliceAI + Pangolin
`80_make_vcf.py` · `82_assemble_v2.py`  (the SpliceAI/Pangolin runs themselves are
the WSL/conda shell scripts in `refs/*.sh` — reference provenance)

## Milestone 7 — GPN-MSA + AlphaGenome  *(re-scoring; needs network / API key)*
`90_score_gpn.py` (remote tabix; needs `pysam` + GPN `.tbi`) ·
`91_score_alphagenome.py` (needs `ALPHAGENOME_API_KEY` env var) · `94_assemble_final.py`

## Milestone 8b / 11 — Nucleotide Transformer  *(reproducible; needs GPU)*
`92_score_nt.py` — **archived NT scoring procedure** (InstaDeepAI/nucleotide-
transformer-v2-500m-multi-species; masked 6-mer REF/ALT log-likelihood ratio,
higher=pathogenic). Regenerates the git-ignored `nt_cache2.tsv`. The original run
was on a cloud GPU; this restores reproducibility — see the file header and README.
Requires a GPU env (`torch`, `transformers`, `pyfaidx`) and the Ensembl release-112
reference FASTA. · `95_check_nt.py` (verify the NT-augmented matrix) · `96_update_coverage.py`

## Milestone 8 — Phase-3 evaluation  ← reproduce these
| Script | Output |
|---|---|
| `100_spearman.py` | `spearman_meta.tsv`, `spearman_by_gene.tsv` — **PRIMARY metric** |
| `101_auroc.py` | `auroc_clinvar.tsv` — auxiliary, with/without BRCA1 |
| `102_clinical_pairwise.py` | `pairwise_splice.tsv`, `vus_reclass.tsv`, `sens_at_95spec.tsv` |
| `103_figures.py` | `fig1_splice_performance.png`, `fig2_forest_*.png`, `fig3_region_heatmap.png` |

## Milestone 9 — one-click reproduction
`reproduce_main.py` — runs 100→103 from `data/processed/score_matrix_final.tsv`.
