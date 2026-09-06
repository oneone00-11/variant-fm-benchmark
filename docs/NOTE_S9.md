# Supplementary Note S9 — Reproducibility

This file is the repository-side copy of **Supplementary Note S9** of the
manuscript. It is kept consistent with the submitted supplement.

The entire analysis is pinned to a **content-hash-verified frozen feature matrix**
and is regenerable from tracked sources by a **two-command frozen-matrix rebuild**
(the `make_raw` + `phase1_build_frozen_matrix` pair below — distinct from the
README's "install, then run" pair, which reproduces the *analysis* from the
already-frozen matrix).

## The frozen analysis set

| | |
|---|---|
| **Frozen matrix** | `phase1/data/frozen/frozen_matrix_v1.parquet` |
| **SHA-256** | `2a0e249b44906f11ed91ce4746aba7389d8fff70d791a17710a3e26ec66d3199` |
| **Content** | 21,410 GRCh38 SNVs across 7 genes; 1,781 splice-region variants |
| **Frozen** | 2026-07-08 (UTC) |
| **Git tag** | `frozen-matrix-v1` (commit `2c78945a41b37fb919d8ae8d94936fe84ccf60d2`) |
| **Build commit** | `10d94bbe01d9834165dcfad31cf48a6bb37ce648` |
| **Random seed** | `20260708` (`RANDOM_SEED` in `phase1/src/config.py`) |
| **Manifest** | `phase1/data/frozen/manifest_v1.json` |

The SHA-256 is taken over a **canonical CSV serialisation** of the matrix with rows
sorted by `variant_id` (see `phase1/src/phase1_build_frozen_matrix.freeze()`), not
over the parquet bytes — parquet is not byte-reproducible across writer versions.
The rebuild below reproduces the SHA-256 above.

## Two-command frozen-matrix rebuild

From `phase1/`:

```bash
python -m src.make_raw                    # regenerate the build input (deterministic)
python -m src.phase1_build_frozen_matrix  # rebuild the frozen matrix + manifest
```

`data/raw/variant_scores.parquet` is intentionally **not** committed; it is
regenerated deterministically from the tracked upstream matrix
`data/processed/score_matrix_final.tsv`, so a clean clone can rebuild
`frozen_matrix_v1` without it.

## Running the full analysis

To rebuild the frozen matrix **and** run every analysis stage behind the
manuscript's results, from the repository root:

```bash
pip install -r requirements.txt
python scripts/reproduce_calibration.py
```

This runs, in order, and verifies the frozen-matrix hash before any analysis:

| Stage | Module | Produces |
|---|---|---|
| 0 | `src.make_raw` | Phase-1 build input |
| 1 | `src.phase1_build_frozen_matrix` | frozen matrix + manifest (hash-verified) |
| 2 | `src.phase1_directionality_check` | control-anchored orientation gate |
| 3 | `src.phase2_model` | H1 fusion vs best single; H2 evolution-axis ablation |
| 4 | `src.phase3_calibration` | H3 calibration — ECE / Brier / clinical yield |
| 5 | `src.phase4_external_tp53` | external validation on held-out TP53 |

Outputs land in `phase1/reports/phase1/`. All randomness is seeded; output is
deterministic.

## Reproducibility scope — an explicit limitation

The **analysis** is fully reproducible on CPU from the frozen matrix. The upstream
model **scores** are not uniformly reproducible:

- **Nucleotide Transformer** — the original scoring was a one-off cloud-GPU run
  whose script and environment were not committed. `scripts/92_score_nt.py` is an
  archival reconstruction of the procedure, not the original code; the exact
  context window used is not confirmed. Re-running it will not be byte-identical.
- **Pangolin** — installed from an unpinned git revision that was not recorded, so
  the scoring environment cannot be reconstructed byte-identically.

Both are therefore **excluded from the TP53 external validation**, which uses an
8-feature fusion (`phase1/src/phase4_external_tp53.py`). This matches Methods §2.5:
the limitation is declared rather than worked around.

## Availability

- **Code + frozen dataset:** `github.com/oneone00-11/variant-fm-benchmark`
- **Functional data:** MaveDB and the cited primary publications (per-gene URNs,
  transcripts, licences and references are tabulated in the top-level `README.md`)
