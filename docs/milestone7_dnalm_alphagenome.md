# Milestone 7 — DNA language models + AlphaGenome

Batch 2B: scored the remaining models onto the 7-gene set (21,410 SNVs). No
performance metrics (Phase 3). Credentials: AlphaGenome key read **only** from
`ALPHAGENOME_API_KEY` env, forwarded into WSL via `WSLENV` — never written to
code, logs, cache, or this report (verified, see end).

## Step 0 — availability & resources (checked before running)

| Model | Status | Resource / method |
|---|---|---|
| **GPN-MSA** | ✅ run | Precomputed genome-wide scores on HuggingFace `songlab/gpn-msa-hg38-scores` (tabix TSV, 81 GB). **Remote-queried** the 7 gene regions via pysam HTTP range (downloaded only the 2.65 MB `.tbi` index) — no big download. |
| **AlphaGenome** | ✅ run | `alphagenome` 0.6.1 client; key via env; 1-variant test call passed (auth + network + splice scorers OK). |
| **Nucleotide Transformer** | ⛔ deferred | Needs GPU for practicality; CPU inference on 21,410×2 sequences infeasible on this box (no GPU, 7 GB RAM). |
| **Evo2** | ⛔ deferred | Requires GPU + large weights/transformer-engine; cannot run on CPU/7 GB. |
| **ESM** | ⛔ deferred | Protein/missense-only (overlaps AlphaMissense); CPU inference slow; low priority. |

Network from this China host was reachable for HuggingFace and googleapis at run time.

## GPN-MSA (done: 21,410 / 21,410)

- Source: `songlab/gpn-msa-hg38-scores` (Benegas et al., *Nat Biotechnol* 2025), hg38,
  columns chrom/pos/ref/alt/score, 1-based, ~9 B SNV predictions.
- Method: per-gene tabix region fetch (pysam, direct HF-CDN URL + local index),
  exact (chrom,pos,ref,alt) match. The HF Xet CDN was flaky from China
  (intermittent timeouts/aborts) → handled with **per-region retry + URL
  re-resolve** (all 7 genes recovered: BRCA1 3893, BARD1 3403, BRCA2 3757,
  PALB2 4339, RAD51C 1880, VHL 1103, BAP1 3035).
- **Directionality: LOWER (more negative) = more deleterious / pathogenic.**

## AlphaGenome (done: 21,410 / 21,410)

- Client `alphagenome` **0.6.1**, recommended splice scorers **SPLICE_SITES +
  SPLICE_SITE_USAGE + SPLICE_JUNCTIONS** (GeneMaskSplicingScorer), interval
  **1 MB** (`SEQUENCE_LENGTH_1MB`).
- Per-variant aggregate = **max |raw_score|** across all splice tracks (the
  REF/ALT effect at splice sites/junctions), via `variant_scorers.tidy_scores`.
- **Directionality: HIGHER = more splice-altering / pathogenic.**
- ~21,410 API calls: **6 concurrent workers**, exponential-backoff retry,
  **incremental cache + resumable** (`ag_cache.tsv`, skips done), run **detached
  (setsid)** so it survived session/sleep. Well under the API's quota.
- AlphaGenome uses `chr`-prefixed chromosomes → added `chr` to our (no-chr) coords
  for the calls; results joined back on the no-chr key.

## Coverage (`results/tables/full_model_coverage.tsv`) — all models

| model | scored / 21,410 | direction |
|---|--:|---|
| phyloP100way | 21,410 | higher=path |
| phastCons100way | 21,410 | higher=path |
| SpliceAI | 21,410 | higher=path |
| Pangolin | 21,410 | higher=path |
| **GPN-MSA** | **21,410** | **lower=path** |
| **AlphaGenome** | **21,410** | **higher=path** |
| CADD | 19,897 | higher=path |
| AlphaMissense | 11,926 (missense) | higher=path |
| gnomAD AF (global/popmax) | 7,708 | higher=benign |
| Nucleotide Transformer / Evo2 / ESM | 0 (deferred) | — |

NA is preserved (never 0-filled); for splice-effect models a near-0 value is
"scored-low", distinguished from NA in the coverage table. **GPN-MSA and
AlphaGenome give full (21,410) coverage** — the two DNA-reading models that, with
SpliceAI/Pangolin, complete the splice/non-coding banner.

## Spot-check — known pathogenic splice (4-way concordant)

| variant | functional | GPN-MSA | SpliceAI | AlphaGenome | ClinVar |
|---|--:|--:|--:|--:|---|
| BRCA1 c.5194-1G>C | −4.20 | −9.49 | 1.00 | 7.38 | P/LP |
| BRCA1 c.4986+1G>A | −3.35 | −7.85 | 0.99 | 3.83 | P/LP |
| VHL c.340+1G>T | −3.21 | −10.73 | 0.98 | 5.32 | P/LP |
| VHL c.341-2A>T | −2.67 | −12.86 | 0.93 | 7.79 | P/LP |

All four DNA/splice models flag these pathogenic splice variants in the correct
direction (GPN-MSA strongly negative; AlphaGenome/SpliceAI high), concordant with
the strongly-damaging functional scores. Aligned on exact GRCh38 (chrom,pos,ref,alt)
→ ref alleles verified by construction.

## Credential-leak check
`grep` for the key value across scripts, docs, `data/raw/scores/`, and all logs
returned **no matches** — the key appears nowhere on disk. `ag_cache.tsv` columns
are only chrom/pos/ref/alt/alphagenome_splice.

## Deliverables
- `data/processed/score_matrix_final.tsv` — all models: CADD, AlphaMissense,
  phyloP, phastCons, gnomAD(×2), SpliceAI(+4 comp), Pangolin(+gain/loss),
  **GPN-MSA, AlphaGenome**, + NA columns for NT/Evo2/ESM, each with `*_covered`.
- `results/tables/full_model_coverage.tsv` — model × region scored/low/NA.
- Caches: `data/raw/scores/gpn/gpn.tsv`, `data/raw/scores/alphagenome/ag_cache.tsv`;
  GPN index `refs/gpn_scores.tsv.bgz.tbi`.

## Pending (future / Phase 3)
- NT / Evo2 / ESM — need GPU (or a GPU cloud env); deferred, NA in matrix.
- Phase 3 = performance metrics (Spearman vs functional score; AUROC vs ClinVar;
  stratified by region/gene; CIs). **No metrics computed in Phase 2.**
