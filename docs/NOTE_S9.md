# Supplementary Note S9 — Reproducibility

This file is the repository-side copy of **Supplementary Note S9** of the
manuscript. It is kept consistent with the submitted supplement: the text below is the
note as the supplement prints it.

The entire analysis is pinned to a content-hash-verified frozen feature matrix and is regenerable from tracked sources by a two-command frozen-matrix rebuild (the make_raw + phase1_build_frozen_matrix pair below) — distinct from the install-then-run pair (pip install -r requirements.txt; python scripts/reproduce_calibration.py), which reproduces the analysis from the already-frozen matrix.

- **Frozen matrix:** frozen_matrix_v2.parquet  —  SHA-256 = 8666d0258e3078d05ff2af131a962f6f18c647559308650b090341440ed54752 (built from frozen_matrix_v1.parquet, SHA-256 = 2a0e249b44906f11ed91ce4746aba7389d8fff70d791a17710a3e26ec66d3199, by replacing the SpliceAI, Pangolin and Nucleotide Transformer columns)
- **Content:** 21,410 GRCh38 SNVs across 7 genes; 1,781 splice-region variants; v1 frozen 2026-07-08 (UTC); v2 built 2026-09-13 (UTC).
- **Git tags:** frozen-matrix-v2 (commit 0fff9715d30bcedb0469738515eaf5213129fd67; build commit 5ea00510277a3cab188b124767d1d49eb52eb75e) and frozen-matrix-v1 (commit 2c78945a41b37fb919d8ae8d94936fe84ccf60d2; build commit 10d94bbe01d9834165dcfad31cf48a6bb37ce648).
- **Rebuild (from repo root, phase1/):**   python -m src.make_raw   then   python -m src.phase1_build_frozen_matrix (rebuilds v1)   then   python -m src.phase1_build_frozen_matrix_v2 (builds v2 from v1 and the tracked pinned-scorer columns)
- Random seed fixed (20260708) and reported; the rebuild reproduces both SHA-256 values above.
- **Supplementary tables:** every pipeline-derived table is written by phase1/src/build_supp_tables.py from the report tables. Before the frozen-matrix-v2 tables were written, the generator was run with --verify against the frozen-matrix-v1 outputs and the supplement as it stood before the revision, and reproduced every existing cell with no numeric mismatch (a cell whose printed precision changed is compared as a number; the columns this revision adds have no earlier cell to compare). Only then was it used to write the v2 tables, which pass the same check against the v2 outputs. Both runs are recorded in docs/supplement-verify.md.
- **Code + frozen dataset:** github.com/oneone00-11/variant-fm-benchmark. Functional data: MaveDB and the cited primary publications.
- **Sampling-frame baseline:** the reweighted analysis (Supplementary Tables S10a–S10d) is evaluated on the 1,768 splice variants that can be placed in the sampling frame — the thirteen without a usable intron offset are excluded at that stage — so its unweighted baseline differs modestly from the headline estimates on all 1,781; both pipelines share the same frozen out-of-fold predictions and the same isotonic calibrator.
