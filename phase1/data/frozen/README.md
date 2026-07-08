# Frozen feature matrix — v1

- **`frozen_matrix_v1.parquet`** — the Phase-1 frozen, versioned feature matrix: 21,410 GRCh38 SNVs × (7 genes) with raw predictor scores, splice sub-classes, two label tracks (`y_clinvar` Method A; official `y_assay` Method B), `*_isna` missingness flags, and provenance columns. Built by `phase1/src/phase1_build_frozen_matrix.py` (`METHOD_B="official"`); no imputation/rank-norm baked in (see `phase1/MISSINGNESS_POLICY.md`).
- **sha256 = `2a0e249b44906f11ed91ce4746aba7389d8fff70d791a17710a3e26ec66d3199`** (over the canonical CSV serialisation; recorded in `manifest_v1.json`). Regenerate identically with `python -m src.phase1_build_frozen_matrix` from `phase1/`.
- **source commit = `6bd274f`** (repo HEAD at freeze time; upstream benchmark scoring `data/processed/score_matrix_final.tsv`). Tagged `frozen-matrix-v1`.

**Clean-clone rebuild** (the build input `data/raw/variant_scores.parquet` is *not* committed — it is deterministically regenerated from the tracked upstream TSV):

```
cd phase1
python -m src.make_raw                        # TSV -> data/raw/variant_scores.parquet
python -m src.phase1_build_frozen_matrix       # -> frozen_matrix_v1.parquet, sha256 2a0e249b…3199
```

Verified: regenerating the input via `make_raw` produces a byte-logically identical input and the build reproduces sha256 `2a0e249b…3199` exactly.
