# Milestone 6 — Local SpliceAI + Pangolin (splice specialists)

Filled the splice/non-coding banner gap left by Milestone 5 (the Broad API was
GFW-blocked) by running **SpliceAI** and **Pangolin** locally. No external API.
No model scoring beyond these two; no performance metrics (that is Phase 3).

## Environment (Step 0)

Native Windows could not run either tool: **`pysam` has no Windows build** and
both CLIs depend on it. Resolved by running in **WSL Ubuntu + conda (miniforge)**:

- `spliceai` env: **SpliceAI 1.3.1**, **tensorflow 2.19.1 (CPU)**, pyfaidx, pysam, samtools (bioconda). Needed `setuptools<80` (SpliceAI imports `pkg_resources`, removed in setuptools ≥81).
- `pangolin` env: **Pangolin** (git `tkzeng/Pangolin`), **torch 2.12.0 (CPU)**, gffutils, pyfastx. Pangolin's VCF writer is incompatible with the installed PyVCF3 (`_Info` signature) → used **CSV input mode** instead.
- CPU-only, no GPU. 16 cores, 7 GB RAM (the binding constraint).

## Reference files (cached on D:, `D:\variant-fm-benchmark\refs`)

- **Ensembl GRCh38 release-112** per-chromosome FASTA for **chr 2, 3, 13, 16, 17**
  (our 7 genes' chromosomes), concatenated → `grch38_subset.fa` + `.fai`.
- **Ensembl GRCh38.112 GTF** → gffutils DB (`grch38_subset.gtf.db`, subset to those
  chroms) for Pangolin's annotation.
- SpliceAI annotation: its bundled **grch38** table.

**Chromosome naming — verified consistent (no `chr` prefix)** across all three:
Ensembl FASTA sequence names (`2`,`3`,…,`17`), SpliceAI's grch38 annotation
(`1`,`2`,…), and our VCF (built from analysis_ready_v2, `17`/`3`). No renaming needed.

## Input & parameters

- Input VCF: all **21,410** unique SNVs from `analysis_ready_v2` (all 7 genes, all
  regions — not only splice), GRCh38, sorted.
- **SpliceAI**: default distance **-D 50**. Summary = **max of the 4 deltas**
  (DS_AG/DS_AL/DS_DG/DS_DL); all 4 components retained.
- **Pangolin**: default distance **-d 50**, default mask. Summary = **max(splice
  gain, |splice loss|)**; gain and loss retained.
- **Directionality (both): higher = more splice-altering = more likely pathogenic.**

## Coverage (`results/tables/splice_model_coverage.tsv`)

Both tools scored **all 21,410 variants (NA = 0)** — every variant lies within an
annotated gene/transcript window. Per the task's rule, a **near-0 score is
"scored-low", not NA** — distinguished below (threshold < 0.01):

| region | n | SpliceAI scored / of which low | Pangolin scored / of which low |
|---|--:|--:|--:|
| splice | 1781 | 1781 / 228 | 1781 / 440 |
| intronic | 1581 | 1581 / 653 | 1581 / 974 |
| missense | 11979 | 11979 / 5433 | 11979 / 8499 |
| synonymous | 4915 | 4915 / 2260 | 4915 / 3564 |
| nonsense | 956 | 956 / 319 | 956 / 352 |
| utr | 180 | 180 / 84 | 180 / 79 |

Interpretation: in the **splice** subset most variants get a *meaningful* (non-low)
score (SpliceAI 1553, Pangolin 1341); in coding subsets most are low (expected — a
missense/synonymous change usually does not alter splicing). The low scores are
genuine model output, not missing data.

## Spot-check (known pathogenic splice → both tools high, concordant with function)

| variant | functional score | SpliceAI | Pangolin | ClinVar |
|---|--:|--:|--:|---|
| BRCA1 c.5194-1G>C | −4.20 | 1.00 | 0.87 | P/LP |
| BRCA1 c.4986+1G>A | −3.35 | 0.99 | 0.83 | P/LP |
| BRCA1 c.5074+1G>A | −2.30 | 0.74 | 0.70 | P/LP |
| VHL c.340+1G>T | −3.21 | 0.98 | 0.86 | P/LP |
| VHL c.341-2A>T | −2.67 | 0.93 | 0.78 | P/LP |

Canonical splice-site P/LP variants get near-maximal SpliceAI **and** Pangolin
scores, concordant with their strongly-damaging functional scores. Ref bases match
the minus-strand genes (alignment correct).

## Alignment & ref-base check

Scores aligned to `score_matrix_v2.tsv` by exact GRCh38 (chrom,pos,ref,alt) — so
every matched score inherently agrees on the reference allele (21,410/21,410 matched
for both tools, independently confirming our ref alleles).

## Problems encountered & how resolved (full reproducibility)
1. **pysam not installable on Windows** → moved to WSL conda.
2. **SpliceAI `pkg_resources` ImportError** → pinned `setuptools<80`.
3. **Pangolin VCF writer crash (PyVCF3 `_Info`)** → CSV input/output mode.
4. **Pangolin OOM at startup** (8 simultaneous torch loads on 7 GB) → 5 processes
   with **staggered (25 s) startup**.
5. **Laptop hibernation / battery-death (×3) interrupting the run** → (a) disabled
   AC+DC sleep/hibernate and lid action via `powercfg`; (b) **cumulative banking +
   self-resuming loop** (`scored_so_far.csv` persisted every pass — progress can't
   be lost); (c) **detached run (`setsid`)** so it survives the session closing.
6. **WSL vhdx on C: growing from conda envs** (C: low) → `conda clean`; outputs all
   write to D:.

## Deliverables
- `data/processed/score_matrix_v2.tsv` — v1 columns (CADD/AlphaMissense/phyloP/
  phastCons/gnomAD) **+ spliceai_ds (+4 components) + pangolin_score (+gain/loss)**
  + coverage flags. AlphaGenome column still NA (batch B, needs Google API key).
- `results/tables/splice_model_coverage.tsv`
- raw outputs cached: `data/raw/scores/spliceai/spliceai_out.vcf`,
  `data/raw/scores/pangolin/pangolin_out.csv`; run scripts in `refs/`.

## Still pending (next batch / Phase 3)
- **AlphaGenome** — needs your Google AlphaGenome API key.
- Self-hosted DNA-LMs (Nucleotide Transformer, Evo2, ESM, GPN-MSA).
- No metrics computed yet (Spearman/AUROC = Phase 3).
