# frozen-matrix-v2

`frozen-matrix-v2` is `frozen-matrix-v1` with three predictor columns replaced and nothing
else changed. It exists so that the effect of two scoring decisions can be measured on the
same 21,410 variants with the same code and seeds:

| column | v1 (published) | v2 |
|---|---|---|
| `spliceai` | SpliceAI 1.3.1 CLI output, printed to 2 d.p. (`-A grch38 -D 50`, max of four deltas) | same tool, annotation and distance; the CLI's `{:.2f}` output formatting removed, so the float the model computes is kept (`functional-standard-atlas/models/spliceai/score_fullprec.py`) |
| `pangolin` | Pangolin from an unpinned `tkzeng/Pangolin` checkout, `-d 50`, `max(gain, |loss|)` printed to 2 d.p. | Pangolin pinned to git `5cf94b8db938c658391b4305cd7ce33297d44ff7`, same distance and definition, `round(…, 2)` formatting removed (`models/pangolin/score_fullprec.py`) |
| `nt` | one-off cloud-GPU run of NT-v2-500M-multi-species whose script was not committed and whose context window was not recorded | the atlas's pinned scorer (`models/nt/score.py`: 6,000-bp window, variant aligned inside one 6-mer, masked-token log-likelihood ratio), run on an NVIDIA A6000 for the 21,394 shared SNVs and on CPU here for the 16 RAD51C SNVs the atlas lacks |

The v1 SpliceAI and Pangolin values are the 2-d.p. roundings of the v2 values (checked row by
row; see `tests/test_frozen_matrix_v2.py`), so the two versions differ only in output precision
for those two columns. For NT the two versions are the same model and score definition under
two runs; they agree at Spearman ρ = 0.9997 over the 21,410 variants.

## Build

```
cd phase1
python -m src.phase1_build_frozen_matrix_v2     # -> data/frozen/frozen_matrix_v2.parquet, manifest_v2.json
```

Inputs, all tracked:

- `phase1/data/frozen/frozen_matrix_v1.parquet` — verified against the v1 pin (`2a0e249b…3199`) before anything is replaced;
- `phase1/data/rescore/atlas_columns_v2.tsv` — the three columns for the 21,394 shared SNVs, extracted from
  `functional-standard-atlas/results/score_matrix_atlas_v2.parquet` (Zenodo 10.5281/zenodo.22671905, release v2.4.4);
  provenance in `atlas_columns_v2.provenance.json`;
- `data/rescore/{spliceai,pangolin,nt}_frozen_only16.tsv` — the sixteen RAD51C SNVs absent from the atlas,
  scored with the same pinned scorers by `scripts/phase4b_score_pinned.py` (provenance JSON beside each);
- `data/tp53/{spliceai,pangolin,nt}_tp53.tsv` — the 192 TP53 splice SNVs, same scorers, which the build also
  folds into `phase1/data/external/tp53_splice_scored_v2.parquet` so the external validation can use all ten predictors.

The build asserts that the variant set and row order are identical to v1 and that every column other
than the three replaced ones (and their `_isna` flags) is value-identical. Those assertions are also
tests (`tests/test_frozen_matrix_v2.py`), together with the content hash.

| | v1 | v2 |
|---|---|---|
| rows / genes / splice | 21,410 / 7 / 1,781 | 21,410 / 7 / 1,781 |
| content sha256 (canonical CSV, rows sorted by `variant_id`) | `2a0e249b44906f11ed91ce4746aba7389d8fff70d791a17710a3e26ec66d3199` | `8666d0258e3078d05ff2af131a962f6f18c647559308650b090341440ed54752` |
| git tag | `frozen-matrix-v1` | `frozen-matrix-v2` |

## What changed in the three columns

From `phase1/data/frozen/frozen_matrix_v2_column_report.tsv`:

| column | Spearman ρ v1 vs v2 | distinct values v1 → v2 (all 21,410) | distinct values v1 → v2 (1,781 splice) | tie-imposed ρ ceiling on the splice subset v1 → v2 |
|---|---|---|---|---|
| spliceai | 0.9599 | 101 → 21,344 | 98 → 1,779 | 0.9982 → 1.000 |
| pangolin | 0.9362 | 92 → 21,373 | 92 → 1,781 | 0.9966 → 1.000 |
| nt | 0.9997 | 21,377 → 21,383 | 1,780 → 1,781 | 1.000 → 1.000 |

The ceiling is the largest Spearman ρ any target could reach against the column given its ties
(the atlas's `max_spearman_given_ties`). On the splice subset the rounding cost was small but real:
0.2–0.3% of the ceiling overall, and per gene between 0.1% (RAD51C, SpliceAI) and 0.9% (BAP1, Pangolin) —
see `docs/column-provenance.md`.

## Running the pipeline on v2

```
FROZEN_VERSION=v2 PHASE1_REPORT_DIR=reports/phase1_v2 python scripts/reproduce_calibration.py
FROZEN_VERSION=v2 PHASE1_REPORT_DIR=reports/phase1_v2 python scripts/reproduce_robustness.py
```

Both scripts rebuild and verify v1, then build v2 from it and verify v2's hash against the pin in
`scripts/reproduce_calibration.py`, exactly as they do for v1; with the variables unset they behave
as before and write to `reports/phase1/`. The TP53 stage reads `tp53_splice_scored_v2.parquet` and
uses all ten predictors when `FROZEN_VERSION=v2` (eight, as published, otherwise).
