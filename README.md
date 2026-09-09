# Variant-FM-Benchmark

**Functional-evidence benchmarking of splice and sequence-based variant-effect
predictors in hereditary-cancer genes**, using **Saturation Genome Editing (SGE) /
MAVE functional assays as the primary, independent gold standard** (not ClinVar
labels — to avoid circularity).

- **Gene panel (7):** BRCA1, BRCA2, BARD1, PALB2, RAD51C, VHL, BAP1
- **Analysis set:** 21,410 SNVs (GRCh38) that have **both** a functional score
  **and** a ClinVar record. Content-hash-pinned as **`frozen-matrix-v1`**
  (`phase1/data/frozen/`); the flat TSV view is `data/processed/score_matrix_final.tsv`
- **External held-out gene:** TP53 (192 splice SNVs), never used for fitting

## What this repository contains

This repository carries **two related analyses** over the same frozen matrix. They
answer different questions and have separate entry points.

### Study A — calibration and fusion (the current manuscript)

*"Ranking is saturated, calibration is not, and neither changes the evidence
strength."* Given that the top splice predictors are statistically tied on
**ranking**, the question becomes whether combining them buys anything ranking
cannot show — namely **calibrated probabilities**, and then whether that gain
survives translation into ACMG evidence. It does not: at a fixed specificity the
evidence strength is a function of the ranking, and the ranking is saturated. Elastic-net fusion, out-of-gene isotonic/Platt calibration,
ECE / Brier / clinical yield, all under **leave-one-gene-out (LOGO)**, plus
external validation on the held-out gene **TP53**.

- Code: [`phase1/src/`](phase1/src/) (`phase2_model.py`, `phase3_calibration.py`,
  `phase4_external_tp53.py`)
- Reproduce: **`python scripts/reproduce_calibration.py`** (see below)
- Results: `phase1/reports/phase1/`
- Method detail: [`docs/NOTE_S9.md`](docs/NOTE_S9.md) (= Supplementary Note S9)

Headline numbers (frozen-matrix-v1, splice subset N = 1,781, LOGO):

| | fusion (elastic net) | best single | Δ (fusion − single) |
|---|--:|--:|---|
| pooled Spearman ρ | 0.773 | 0.761 (Pangolin) | **+0.012** `[0.002, 0.021]` |
| Brier (isotonic, `y_assay`) | 0.0631 | 0.0717 (Pangolin) | **−0.0085** `[0.0032, 0.0165]` |
| Brier — TP53, fully held out | 0.0785 | 0.1216 (SpliceAI) | **−0.0431** (directional replication) |

Ranking is effectively saturated (Δρ ≈ 0.01); calibration is not — and the
calibration gap replicates on a gene the model has never seen.

### Study B — coverage and model complementarity

*Can DNA-reading models score the splice / non-coding variants that protein models
(AlphaMissense) cannot?* Per-gene Spearman vs. the continuous functional score,
pooled by DerSimonian–Laird random-effects meta-analysis; ClinVar binary AUROC as
an auxiliary metric, always reported with/without BRCA1 (BRCA1's ClinVar labels are
circular with its SGE).

- Code: [`scripts/`](scripts/) (numbered by phase — see `scripts/README.md`)
- Reproduce: **`python scripts/reproduce_main.py`**
- Results: `results/tables/`, `results/figures/`
- Method/results detail: `docs/milestone1`–`milestone8b`

## Study B headline result

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
├── scripts/                   Study B pipeline + analysis, numbered by phase (see scripts/README.md)
│   ├── config.py              shared paths, gene list, ClinVar maps
│   ├── phase3_lib.py          orientation, Fisher-z CIs, DL meta, bootstrap
│   ├── reproduce_main.py      >>> one-command reproduction of Study B <<<
│   └── reproduce_calibration.py  >>> one-command reproduction of Study A <<<
├── phase1/                    Study A (fusion + calibration + TP53)
│   ├── src/                   config.py, phase1_* (build/freeze), phase2_model.py
│   │                          (fusion), phase3_calibration.py, phase4_external_tp53.py
│   ├── data/frozen/           frozen-matrix-v1 + manifest (sha256-pinned)
│   └── reports/phase1/        Study A result CSVs + figures
├── docs/                      milestone reports 1–8b (Methods/Results record)
│   └── NOTE_S9.md             Supplementary Note S9 (Study A reproduction)
├── results/
│   ├── tables/                all Study B result TSVs (spearman_meta, auroc_clinvar, …)
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

## Reproduce the results

Clean-clone reproduction is two commands per study — install, then run. Neither
re-scores any model: both start from the frozen matrix that ships in the repo.
(Note S9's "two-command frozen-matrix rebuild" is a different command pair —
`make_raw` + `phase1_build_frozen_matrix`, for rebuilding the frozen matrix
itself; see [docs/NOTE_S9.md](docs/NOTE_S9.md).)

**Study A — fusion & calibration** (the current manuscript):
```bash
pip install -r requirements.txt
python scripts/reproduce_calibration.py
```
Runs, in order: build input → frozen-matrix integrity check (sha256 against
`phase1/data/frozen/manifest_v1.json`) → directionality gate → H1/H2 fusion →
H3 calibration → TP53 external validation → H4 likelihood ratios and Tavtigian
bands → evidence yield → weight stability and training-gene thinning → label
contrast and selection test → intronic-offset drop sensitivity → sampling-frame
reweighting (run and balance diagnostics). Thirteen stages. Writes
`phase1/reports/phase1/` and prints the headline Δρ / ΔBrier table. Seeded
(`RANDOM_SEED` in `phase1/src/config.py`); deterministic.
Checkout `frozen-matrix-v1` to pin the exact matrix the manuscript used:
```bash
git checkout frozen-matrix-v1
```

**Study A — reviewer robustness analyses** (Phase 8: exact sign-flip tests,
leave-two-genes-out, Hartung–Knapp CIs, Murphy decomposition, stratified
selection test, LR+ operating points):
```bash
pip install -r requirements.txt
python scripts/reproduce_robustness.py
```
Same frozen matrix and hash check as above; writes `phase1/reports/phase1/phase8_*.csv`.

**Between the two scripts, every tracked file in `phase1/reports/phase1/` is
regenerated** — `reproduce_calibration.py` writes all of it except the
`phase8_*` tables, which `reproduce_robustness.py` writes.
`tests/test_reproduction_coverage.py` fails if any tracked result is produced by
a module that neither script runs, so the entry points cannot fall behind the
pipeline again.

**Study B — coverage & complementarity:**
```bash
pip install -r requirements.txt
python scripts/reproduce_main.py
```
Starting from `data/processed/score_matrix_final.tsv` (no model re-scoring), this
regenerates every table in `results/tables/` and all three figures in
`results/figures/`. All randomness is seeded (`np.random.default_rng(20260619)`),
so output is deterministic. The script prints the headline splice ranking at the end.

### Reproducibility scope — what is and is not reproducible

**The analyses above are fully reproducible on CPU from the committed matrix.**
**The upstream model *scores* are not all reproducible**, and two are known not to be:

| Model | Score reproducible? | Why |
|---|---|---|
| Nucleotide Transformer | **No** | The original scoring was a one-off cloud-GPU run; that script and its environment were not committed. `scripts/92_score_nt.py` is an *archival reconstruction* of the procedure, not the original code, and the exact context window of the original run is not confirmed. Re-running it needs a GPU and will not be byte-identical. |
| Pangolin | **No** | Installed from an unpinned git revision; the exact revision was not recorded, so the scoring environment cannot be reconstructed byte-identically. |
| All others | Yes | Deterministic given the pinned inputs in "Data acquisition". |

Both models are therefore **excluded from the TP53 external validation** (Study A
uses an 8-feature fusion there; see `phase1/src/phase4_external_tp53.py`). This
matches the manuscript's Methods §2.5 — the limitation is declared, not worked
around.

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

## Models / tools evaluated

| Model | Kind | Version / source |
|---|---|---|
| SpliceAI | supervised splice | Illumina, default models |
| Pangolin | supervised splice | Tiger Genomics |
| AlphaGenome | DNA model (splice tracks) | Google AlphaGenome API |
| GPN-MSA | DNA LM (multiple-sequence alignment) | Song Lab, hg38 precomputed |
| Nucleotide Transformer | DNA LM (single sequence) | InstaDeepAI/nucleotide-transformer-v2-500m-multi-species; masked 6-mer REF/ALT log-likelihood ratio. Reproducible via `scripts/92_score_nt.py` (needs GPU; regenerates the git-ignored `nt_cache2.tsv`) |
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

## Author

**Ningyi Zhang** — Department of Biological Sciences, National University of
Singapore · cliffzhang@u.nus.edu ·
[ORCID 0009-0004-3020-4044](https://orcid.org/0009-0004-3020-4044)

Development was AI-assisted (Claude, Anthropic) under the author's direction;
the author designed the study, verified all reported numbers, and takes full
responsibility for the content.

## License

Code: MIT (see [LICENSE](LICENSE)). Third-party data and model scores retain their
own licenses (CC0 / CC BY 4.0 for MaveDB sets; AlphaMissense CC BY-NC-SA 4.0,
non-commercial) — see the table above and the data-licensing note in `LICENSE`.

## Status

Both analyses are complete and reproducible on CPU from the committed frozen matrix.
Model *scoring* is not uniformly reproducible: the Nucleotide Transformer and
Pangolin scores cannot be regenerated byte-identically (see "Reproducibility scope"
above) — this is a declared limitation, and both are excluded from the external
validation. Pending future work: Evo2 / ESM (need GPU; columns reserved as NA).
