# Milestone 5 — Panel merge (7 genes) + first-batch model scoring

First batch: **ready-made / precomputed scores only**. No self-hosted DNA-LM
inference (NT/Evo2/ESM = batch 2). **No performance metrics** (that is Phase 3).
This milestone produces the aligned score matrix + coverage audit.

## A. Merge → analysis_ready_v2 (7 genes)

BAP1 + VHL merged into the same column structure & region rule as analysis_ready.

- **analysis_ready_v2.tsv: 21,410 rows, 7 genes** (existing 5 unchanged at 17,272;
  BAP1 +3,035, VHL +1,103).
- Completeness check (matches M3+M4): region_class × gene, label_set per gene, and
  circularity_flag (only **BRCA1 = TRUE**; BAP1/VHL/others = FALSE) all consistent.
- label_set totals: clean_PB 8,373 · VUS 7,851 · conflicting 3,120 · other 2,066.
- Spot-check passed: VHL splice_core P/LP strongly damaging (c.340+1G>T −3.21,
  c.463+1G>C −2.84); BAP1 splice_core P/LP at the bottom of its compressed scale
  (−0.20..−0.21). Ref bases consistent with gene strand.

## B. First-batch scores (aligned on exact GRCh38 (chrom,pos,ref,alt))

`data/processed/score_matrix.tsv` — one row per variant, one column per model +
a `*_covered` flag. **Coverage is left NA where a model does not score a variant
— no 0-fill, no extrapolation.** Because alignment uses the exact 4-tuple
(chrom,pos,ref,alt), every matched score inherently agrees on the reference
allele (ref-base verified by construction; e.g. 19,897 variants independently
found in CADD at the same ref/alt confirm our GRCh38 ref alleles).

| Model | Source / version (fetched 2026-06-17) | Directionality | Scored / 21,410 | By-design scope |
|---|---|---|--:|---|
| **CADD** PHRED | API `cadd.gs.washington.edu` **GRCh38-v1.7** | higher = pathogenic | 19,897 (93%) | all regions (key baseline) |
| **AlphaMissense** | Zenodo rec. 8208688, `AlphaMissense_hg38.tsv.gz` | higher = pathogenic | 11,926 | **missense only** |
| **phyloP100way** | UCSC track API, hg38 | higher = pathogenic | 21,410 (100%) | all (conservation) |
| **phastCons100way** | UCSC track API, hg38 | higher = pathogenic | 21,410 (100%) | all (conservation) |
| **gnomAD AF** (global+popmax) | GraphQL API, **gnomad_r4** (v4) | higher = **benign** | 7,708 (36%) | all (naïve freq baseline) |
| SpliceAI Δ | — | (higher = pathogenic) | **0 — not run** | splice-proximal |
| Pangolin | — | (higher = pathogenic) | **0 — not run** | splice-proximal |
| AlphaGenome splice | — | (higher = pathogenic) | **0 — not run** | splice |

Full model × region coverage: `results/tables/score_coverage_report.tsv`.

### Coverage notes (the coverage IS a result)
- **gnomAD 36%**: the other 64% are simply absent from gnomAD v4 (rare/unobserved).
  NA here means "not in gnomAD", NOT AF = 0 — left NA per guardrail; the eval stage
  decides how to treat gnomAD-absent variants.
- **AlphaMissense** is missense-by-design; it scored 11,906 of our missense variants
  plus ~20 variants we labelled nonsense/other (cross-transcript consequence edge
  cases where AM's canonical protein call is missense). Non-missense = NA, not extrapolated.
- **CADD 93%**: ~1,500 SNVs returned an empty CADD API record (left NA).

## ⚠️ Key limitation — the banner region is under-covered this batch

For the **splice subset (1,781 variants)** the only scores so far are generic:
CADD 1,644, phyloP/phastCons 1,781, gnomAD 526 — while the **splice-specialist
models (SpliceAI, Pangolin, AlphaGenome) are 0**. Those three are exactly the
models that matter most for the project's splice/non-coding banner, and all three
were blocked this batch:

- **SpliceAI & Pangolin** — the Broad lookup API (`spliceailookup-api.broadinstitute.org`)
  is **SSL-reset from this network** (host-specific GFW block; note gnomAD on the
  same broadinstitute.org domain works). Options for batch 2: (a) run SpliceAI &
  Pangolin locally (needs the reference FASTA + the annotation GTF; CPU is fine for
  ~21k variants, no GPU required for SpliceAI raw); (b) download Illumina's
  precomputed SpliceAI SNV VCF; (c) reach the Broad API from a different network/VPN.
- **AlphaGenome** — needs a **Google AlphaGenome API key** (not available here) and
  the `alphagenome` Python client; the Google endpoint is also likely proxy-gated
  from this network. Requires you to provide an API key.

These belong to batch 2 alongside the self-hosted DNA-LMs.

## Deliverables
- `data/processed/analysis_ready_v2.tsv` (7 genes, 21,410)
- `data/processed/score_matrix.tsv` (variant × model scores + coverage flags +
  gene/region_class/splice_class/functional_score/clnsig_class/label_set/circularity_flag)
- `results/tables/score_coverage_report.tsv` (model × region scored/total/NA)
- raw caches in `data/raw/scores/` (cadd_cache.json, AlphaMissense_hg38.tsv.gz,
  gnomad_af.tsv, conservation.tsv, *_source.txt) — on D: drive, reproducible.

## Not done (by design / blocked) — for the next batch
- Self-hosted DNA-LMs (Nucleotide Transformer, Evo2, ESM, GPN-MSA…) — batch 2.
- SpliceAI, Pangolin, AlphaGenome — blocked this batch (see above).
- GERP — not natively served for hg38 by the UCSC API; phyloP100way + phastCons100way
  used as the conservation baselines instead (GERP hg38 bigWig deferred if needed).
- No Spearman/AUROC/any performance metric — Phase 3.
- No GRCh37 liftover needed (all sources provided GRCh38).
