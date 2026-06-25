# Variant-FM-Benchmark

**Functional-evidence benchmarking of splice and sequence-based variant-effect
predictors in hereditary-cancer genes.**

We benchmark language models and classical predictors for the clinical
interpretation of cancer-gene variants, using **Saturation Genome Editing (SGE) /
MAVE functional assays as the primary, independent gold standard** (not ClinVar
labels — to avoid circularity). The signature question: **can DNA-reading models
score the splice / non-coding variants that protein models (AlphaMissense) cannot?**

- **Gene panel (7):** BRCA1, BRCA2, BARD1, PALB2, RAD51C, VHL, BAP1
- **Evaluation set:** 21,410 SNVs (GRCh38) that have **both** a functional score
  **and** a ClinVar record — `data/processed/score_matrix_final.tsv`
- **Primary metric:** per-gene Spearman of model score vs. continuous functional
  score, pooled across genes by DerSimonian–Laird random-effects meta-analysis
- **Auxiliary metric:** ClinVar binary AUROC, always reported with/without BRCA1
  (BRCA1's ClinVar labels are circular with its SGE)

## Headline result

On the splice subset (N = 1,781), meta Spearman ρ vs. the functional gold standard:

| Model | meta ρ | Type |
|---|--:|---|
| Pangolin | **0.76** | splice specialist |
| SpliceAI | **0.75** | splice specialist |
| AlphaGenome | **0.75** | splice specialist (DNA) |
| CADD | 0.69 | classical |
| GPN-MSA | 0.67 | DNA language model (MSA) |
| phyloP / phastCons | 0.61–0.64 | conservation |
| Nucleotide Transformer | 0.49 | DNA language model (single-seq) |
| gnomAD AF | 0.19 | frequency baseline |
| **AlphaMissense** | **— (n=1/1781)** | protein/missense — **cannot score splice** |

The three splice specialists are the top, **statistically tied** tier; the two DNA
language models follow (MSA-based GPN-MSA > single-sequence NT); the protein model
is structurally blind to splice. Models are **complementary** — each strong only in
its own regime. Full results: [docs/milestone8_results.md](docs/milestone8_results.md)
and [docs/milestone8b_results.md](docs/milestone8b_results.md).

## Repository layout

```
variant-fm-benchmark/
├── README.md                  this file
├── LICENSE                    MIT (code); third-party data licenses noted within
├── requirements.txt           pip deps (CPU-only)  /  environment.yml (conda)
├── data/
│   └── processed/             evaluation tables (IN repo).  score_matrix_final.tsv = FINAL matrix
│       (data/raw/ is git-ignored — large, re-downloadable; see "Data acquisition")
├── scripts/                   pipeline + analysis code, numbered by phase (see scripts/README.md)
│   ├── config.py              shared paths, gene list, ClinVar maps
│   ├── phase3_lib.py          orientation, Fisher-z CIs, DL meta, bootstrap
│   ├── reproduce_main.py      >>> one-command reproduction of all main results <<<
│   └── 0x–10x_*.py            data assembly → model scoring → Phase-3 evaluation
├── docs/                      milestone reports 1–9 (the written record / Methods+Results)
├── results/
│   ├── tables/                all result TSVs (spearman_meta, auroc_clinvar, …)
│   └── figures/               fig1 splice perf, fig2 forest, fig3 region heatmap
└── refs/                      scoring provenance scripts (FASTA/GTF binaries git-ignored)
```

## Installation (CPU-only — no GPU)

```bash
git clone <repo-url> && cd variant-fm-benchmark
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
# (conda alternative:  conda env create -f environment.yml)
```
Python 3.9; tested with pandas 2.3.3 / numpy 2.0.2 / scipy 1.13.1 /
scikit-learn 1.6.1 / matplotlib 3.9.4.

## Reproduce the main results (one command)

```bash
python scripts/reproduce_main.py
```
Starting from `data/processed/score_matrix_final.tsv` (no model re-scoring), this
regenerates every table in `results/tables/` and all three figures in
`results/figures/`. All randomness is seeded (`np.random.default_rng(20260619)`),
so output is deterministic. The script prints the headline splice ranking at the end.

## The final evaluation matrix — `data/processed/score_matrix_final.tsv`

One row per SNV (21,410). Key columns:
- **Variant:** `gene, chrom, pos, ref, alt, region_class, splice_class, intron_offset, hgvs_nt`
- **Gold standards:** `functional_score` (continuous; **more negative = more
  damaging**, per-gene scale), `clnsig_class`, `path_binary` (1=P/LP, 0=B/LB),
  `label_set`, `circularity_flag` (TRUE only for BRCA1), `source_urn`, `assay_type`
- **Model scores** (+ `*_covered` flags): `cadd_phred, alphamissense, phylop100way,
  phastcons100way, gnomad_af_global, gnomad_af_popmax, spliceai_ds, pangolin_score,
  gpn_msa_score, alphagenome_splice, nucleotide_transformer` (`evo2, esm` reserved, NA)

**Directionality** (see `scripts/phase3_lib.py` `ORIENT`): all oriented so higher =
more pathogenic. GPN-MSA and gnomAD AF are negated; functional pathogenicity = `-functional_score`.

## Data acquisition (inputs NOT committed — large / re-downloadable)

`data/raw/` and the large `refs/` binaries are git-ignored. To rebuild the pipeline
from scratch (not needed to reproduce results — the final matrix ships in the repo):

| Input | Source | Notes |
|---|---|---|
| ClinVar VCF (GRCh38) | `ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz` | release 2026-06-15; ≥1★ kept |
| gnomAD allele frequencies | gnomAD v4 (Broad) | benign proxy + popmax |
| Reference genome / GTF | Ensembl **release-112** GRCh38 (FASTA + GTF) | for SpliceAI/Pangolin local scoring |
| AlphaMissense (hg38) | Zenodo `AlphaMissense_hg38.tsv.gz` | missense only; CC BY-NC-SA 4.0 |
| GPN-MSA scores | HuggingFace `songlab/gpn-msa-hg38-scores` (remote tabix) | only the small `.tbi` index is local |
| MaveDB SGE score sets | MaveDB API `api.mavedb.org/api/v1` (URNs below) | the functional gold standard |

China-network note: pip via Tsinghua mirror + `--timeout 180 --retries 10`; HuggingFace
via `hf-mirror.com` if slow. (No proxy is required on a normal connection.)

## Models / tools evaluated

| Model | Kind | Version / source |
|---|---|---|
| SpliceAI | supervised splice | Illumina, default models |
| Pangolin | supervised splice | Tiger Genomics |
| AlphaGenome | DNA model (splice tracks) | Google AlphaGenome API |
| GPN-MSA | DNA LM (multiple-sequence alignment) | Song Lab, hg38 precomputed |
| Nucleotide Transformer | DNA LM (single sequence) | InstaDeep; masked REF/ALT log-likelihood ratio |
| AlphaMissense | protein/missense LM | Google DeepMind, hg38 precomputed |
| CADD | classical ensemble | CADD v1.6+ (GRCh38) |
| phyloP / phastCons | conservation | UCSC 100-way |
| gnomAD AF | frequency baseline | gnomAD v4 |

## Functional gold-standard datasets (cite & attribute)

Primary SGE/MAVE score set per gene (the independent gold standard):

| Gene | MaveDB URN | Transcript | License | Reference |
|---|---|---|---|---|
| BRCA1 | urn:mavedb:00000097-0-2 | NM_007294.3 | CC0 | Findlay 2018 (PMID 30209399) |
| BRCA2 | urn:mavedb:00001225-a-1 | ENST00000380152.8 | CC0 | 2025 (PMID 39779857) |
| BARD1 | urn:mavedb:00001250-a-2 | NM_000465.4 | CC0 | — |
| PALB2 | urn:mavedb:00001259-a-2 | NM_024675.4 | CC0 | — |
| RAD51C | urn:mavedb:00000673-0-1 | ENST00000337432.9 | CC BY 4.0 | 2024 (PMID 39299233) |
| VHL | urn:mavedb:00000675-a-1 | ENST00000256474.3 | CC BY 4.0 | Buckley 2023 (PMID 38969834) |
| BAP1 | urn:mavedb:00000662-0-1 | ENST00000460680.6 | CC BY 4.0 | Waters 2024 (PMID 38969833) |

ENST-based sets (BRCA2, RAD51C, VHL, BAP1) were coordinate-mapped to GRCh38 via
their MANE-Select RefSeq equivalents and validated 12/12 against MaveDB's own
post-mapping (see [docs/milestone2_goldstandard.md](docs/milestone2_goldstandard.md),
[docs/milestone4_candidate_genes.md](docs/milestone4_candidate_genes.md)).

## Credentials

No credentials are needed to reproduce the results. **Re-scoring** AlphaGenome
requires a Google AlphaGenome API key, read **only** from the environment variable
`ALPHAGENOME_API_KEY` (never hard-coded, cached, or committed):
```bash
export ALPHAGENOME_API_KEY=...   # only for scripts/91_score_alphagenome.py
```

## License

Code: MIT (see [LICENSE](LICENSE)). Third-party data and model scores retain their
own licenses (CC0 / CC BY 4.0 for MaveDB sets; AlphaMissense CC BY-NC-SA 4.0,
non-commercial) — see the table above and the data-licensing note in `LICENSE`.

## Status

All milestones complete (see `docs/`). Phase-3 evaluation is fully reproducible on
CPU. Pending future work: Evo2 / ESM (need GPU; columns reserved as NA).
