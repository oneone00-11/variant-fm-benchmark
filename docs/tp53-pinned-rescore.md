# P1 run record — TP53 and sanity scoring with the atlas's pinned scorers

Driver: `scripts/phase4b_score_pinned.py` (imports the atlas scorer modules by path; contains no scoring logic). Reference: Ensembl GRCh38 release-112 subset FASTA (chr 2/3/13/16/17) and, for Pangolin, the gffutils database built from the release-112 GTF — the same files the atlas and this repository's seven-gene runs used.

## Provenance

| model | set | n | scored | tool / version | pin | key parameters | device | runtime (s) | input sha256 |
|---|---|---|---|---|---|---|---|---|---|
| pangolin | tp53 | 192 | 192 | Pangolin 1.0.2 torch 2.13.0 | `5cf94b8db938c658391b4305cd7ce33297d44ff7` | distance=50; mask=True (default); score_cutoff=None; score_exons=False | cpu | 821.3 | `2410f56e6d0c…` |
| pangolin | sanity | 200 | 200 | Pangolin 1.0.2 torch 2.13.0 | `5cf94b8db938c658391b4305cd7ce33297d44ff7` | distance=50; mask=True (default); score_cutoff=None; score_exons=False | cpu | 617.5 | `64e31fb1dced…` |
| pangolin | frozen_only16 | 16 | 16 | Pangolin 1.0.2 torch 2.13.0 | `5cf94b8db938c658391b4305cd7ce33297d44ff7` | distance=50; mask=True (default); score_cutoff=None; score_exons=False | cpu | 44.3 | `64e31fb1dced…` |
| spliceai | tp53 | 192 | 192 | SpliceAI 1.3.1 TF 2.19.1 | `1.3.1` | annotation=grch38; distance=50; mask=0 | cpu | 205.3 | `2410f56e6d0c…` |
| spliceai | sanity | 200 | 200 | SpliceAI 1.3.1 TF 2.19.1 | `1.3.1` | annotation=grch38; distance=50; mask=0 | cpu | 240.2 | `64e31fb1dced…` |
| spliceai | frozen_only16 | 16 | 16 | SpliceAI 1.3.1 TF 2.19.1 | `1.3.1` | annotation=grch38; distance=50; mask=0 | cpu | 26.9 | `64e31fb1dced…` |
| nt | tp53 | 192 | 192 | Nucleotide Transformer v2 500M multi-species  torch 2.14.0 | `06615c1660c892fc199840c18123f8385b3542a8` | window_bp=6000; kmer=6 | cpu | 494.4 | `2410f56e6d0c…` |
| nt | sanity | 200 | 200 | Nucleotide Transformer v2 500M multi-species  torch 2.14.0 | `06615c1660c892fc199840c18123f8385b3542a8` | window_bp=6000; kmer=6 | cpu | 384.7 | `64e31fb1dced…` |
| nt | frozen_only16 | 16 | 16 | Nucleotide Transformer v2 500M multi-species  torch 2.14.0 | `06615c1660c892fc199840c18123f8385b3542a8` | window_bp=6000; kmer=6 | cpu | 33.5 | `64e31fb1dced…` |

Full provenance (reference and annotation hashes, scorer-module hashes, environment versions) is in the `*.provenance.json` next to each table.

## Sanity check — 200 splice variants of frozen-matrix-v1

| model | n | Spearman ρ vs v1 column | Pearson r | identical after rounding to 2 d.p. | max |Δ| | expected |
|---|---|---|---|---|---|---|
| pangolin | 200 | 0.9974 | 1.0000 | 1.000 | 2.861e-08 | ρ near 1.000 (bounded by the ties in the rounded v1 column); re-rounded values identical |
| spliceai | 200 | 0.9981 | 1.0000 | 1.000 | 0 | ρ near 1.000 (bounded by the ties in the rounded v1 column); re-rounded values identical |
| nt | 200 | 0.9997 | 0.9995 |  | 0.3714 | ρ ≈ 0.9997 (the atlas's concordance against this column) |

TP53 SpliceAI: the full-precision scores re-round to the existing TP53 column exactly (max |round(v2,2) − v1| = 0 over 192 variants; distinct values 192 against 41).

## Outputs

```
data/tp53/nt_tp53.provenance.json
data/tp53/nt_tp53.tsv
data/tp53/pangolin_tp53.provenance.json
data/tp53/pangolin_tp53.tsv
data/tp53/spliceai_tp53.provenance.json
data/tp53/spliceai_tp53.tsv
data/rescore/nt_frozen_only16.provenance.json
data/rescore/nt_frozen_only16.tsv
data/rescore/nt_sanity.provenance.json
data/rescore/nt_sanity.tsv
data/rescore/pangolin_frozen_only16.provenance.json
data/rescore/pangolin_frozen_only16.tsv
data/rescore/pangolin_sanity.provenance.json
data/rescore/pangolin_sanity.tsv
data/rescore/spliceai_frozen_only16.provenance.json
data/rescore/spliceai_frozen_only16.tsv
data/rescore/spliceai_sanity.provenance.json
data/rescore/spliceai_sanity.tsv
```
