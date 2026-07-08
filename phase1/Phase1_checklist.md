# Phase 1 — build & freeze the feature matrix

**Goal of this phase:** produce one *frozen, versioned, audited* raw feature matrix with splice sub-classes, two sets of functional-class labels, and a full missingness audit — so that Phase 2 (the fusion model) starts from a fixed, traceable target instead of a moving one.

**Estimated:** 5–8 focused working days.

**What this phase deliberately does NOT do:** no imputation, no rank-normalisation, no model. Those are fold-dependent and are fit *inside* each training fold in Phase 2. Phase 1 only builds, labels, audits, freezes.

---

## Setup (30 min)

```
your-repo/
└── phase1/
    ├── src/
    │   ├── config.py                     # <- you edit this
    │   └── phase1_build_frozen_matrix.py # <- you run this
    └── Phase1_checklist.md
```

Run from the `phase1/` folder: `python -m src.phase1_build_frozen_matrix`
Requires: pandas, numpy, scikit-learn, pyarrow.

---

## Day 1 — wire it to your real data
- [ ] In `config.py`, set `RAW_MATRIX_PATH` to your actual scoring matrix.
- [ ] Fill in the **right-hand side** of every entry in `COLUMNS` to match your column names. Set any predictor you don't yet have to `None`.
- [ ] Run the pipeline. It will tell you exactly which mapped columns are missing — fix `config.py` until it runs.
- **Acceptance:** `manifest_v1.json` is written and `n_rows` ≈ your expected total (~21,410).

## Day 2 — get the splice subset right (this is the highest-stakes step)
- [ ] If you already have a trusted region/consequence label from VEP, set `COLUMNS["region_label"]` to it — the pipeline will honour it verbatim. **Prefer this** over HGVS parsing.
- [ ] If not, the pipeline derives the offset from `hgvs_c`. **Spot-check 20 variants by hand**: does `intron_offset` match the `+N/-N` in the HGVS string? Does `region` (core/region/intronic) look right?
- [ ] Confirm your splice-subset size lands near the benchmark's ~1,781. A large discrepancy means the region rule or the offset parse is off.
- **Acceptance:** `subset_counts.csv` per-gene splice counts are biologically plausible; splice total is in the right ballpark.

## Day 3 — Method A labels (ClinVar clean_PB)
- [ ] Check that your ClinVar strings actually contain "Pathogenic"/"Benign"/"Conflicting" text. If your column uses codes/enums instead, adjust the matching in `classes_from_clinvar`.
- [ ] Review `class_balance.csv`: are the `y_clinvar` positive/negative counts sane, and does the with/without-BRCA1 split behave as expected?
- **Acceptance:** clean P/LP and B/LB counts match a manual ClinVar filter you trust.

## Day 4–5 — Method B labels (assay-intrinsic) — the fiddly one
- [ ] The pipeline auto-splits each gene's functional score into damaging/normal via a 2-component mixture. **This is a placeholder.** For every gene where the original SGE/MAVE paper reports an official cutoff (control-based, or a stated bimodal threshold), **override `y_assay` with that published threshold.** Hand-tuning here is worth the time — it's what makes "independent functional standard" defensible.
- [ ] Sanity-check orientation: plot `func_pathogenicity` per gene; the damaging mode must be the *higher* one. If a gene is flipped, add it to `FLIP_GENES` in `config.py` and re-run.
- [ ] Compare `y_assay` vs `y_clinvar` agreement on the clean_PB subset — high agreement is a good sign; systematic disagreement in one gene flags an orientation or threshold problem.
- **Acceptance:** per-gene damaging/normal split matches the assay paper's own classification where one exists.

## Day 6 — missingness audit & decision
- [ ] Read `coverage_by_region.csv`. Confirm the expected pattern: AlphaMissense ≈ 0 on the splice rows; splice tools and conservation ≈ complete.
- [ ] **Decide and write down** the Phase-2 handling per predictor: native NA (trees) vs indicator + median-impute-in-fold (linear/MLP). The `*_isna` flags are already in the frozen matrix for this.
- **Acceptance:** a one-paragraph note in your repo stating the coverage facts and your handling decision.

## Day 7–8 — freeze, verify, tag
- [ ] Re-run end-to-end. Record the `sha256` from `manifest_v1.json`.
- [ ] Commit the frozen matrix + manifest, then **git-tag** the commit (e.g. `frozen-matrix-v1`). This is the version Phase 2 pins to.
- [ ] Write a 3-line README next to the frozen file: what it is, the sha256, and the source commit.
- **Acceptance:** you can regenerate the identical `sha256` from a clean checkout. If the upstream benchmark scoring later changes (e.g. after reviewer requests), bump `FROZEN_VERSION` to `v2` and re-freeze — never edit `v1` in place.

---

## Guardrails baked into the code (don't undo these)
- **Gene is never a feature.** It's used only for labelling and (in Phase 2) for the leave-one-gene-out folds. Adding gene identity as a model input would break the generalisation test.
- **No cross-gene statistics are frozen in.** Only raw scores, labels, and is-NA flags are frozen. Imputation/rank-norm stay in Phase 2 folds.
- **BRCA1 circularity is handled by reporting, not deletion.** Every ClinVar-referenced and calibration table is produced with *and* without BRCA1. Keep it that way.
- **Two label definitions travel together.** `y_clinvar` and `y_assay` are both carried forward so the calibration step (Phase 3) can show results are not an artefact of ClinVar dependence.

## Definition of done for Phase 1
A single `frozen_matrix_v1.parquet` + `manifest_v1.json`, git-tagged, with: correct splice sub-classes, two label columns, per-predictor missingness flags, and three diagnostic CSVs whose numbers you have personally eyeballed and trust. When that's true, Phase 2 can begin.
