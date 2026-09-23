# Evidence-strength rework: change log by round

A dated log, kept by hand. It records what changed between rounds and why; the
numbers in it are as they stood when each round was written. The current value of
every quantity is in `phase1/reports/evidence/`, and `docs/evidence-delta.md` (written
by `src/evid_delta.py`) is the generated comparison with the published analysis set.

## Round 2 (2026-09-20): what moved relative to the first round's report

A list, no interpretation.

### Decisions that changed which number is quoted

| | first round | now |
|---|---|---|
| the pool every "overall" number is read off | `all_1_50` | `s3_50` (3-50 bp, pm12 excluded) |
| `all_1_50` in the tables | the pooled row | kept, column labelled `out_of_scope_incl_pm12` |
| the 588 ClinVar records with no clinical assertion | inside `recorded_unclassified`, not visible | still inside it; `clinvar_arm_strict` breaks them out beside the three arms, and `clinvar_arm_no_assertion.csv` is the crosstab |
| tools in the panel | 8 + spliceai_walker | + `avi`, `avi_splice_sites`, `avi_splice_site_usage`, `avi_splice_junctions` |

### Numbers that moved because the pool changed

| quantity | `all_1_50` | `s3_50` |
|---|---|---|
| Walker cut point, LR+ (all arms, spliceai_walker) | 9.55 | 8.28 |
| Walker cut point, LR(<=0.1) | 0.126 (Moderate) | 0.252 (Supporting) |
| sensitivity at >= 0.2 | 0.849 | 0.696 |
| ClinVar arms on AUROC: tools where `classified` is highest | 13/13 | 3/13 |
| the same, BRCA1 dropped | 13/13 | 2/13 |
| ClinVar arms on Walker LR+: tools where `classified` is highest | 13/13 | 11/13 |

### New quantities

- `avi` against functional pathogenicity: rho 0.393 (0.424 at 3-10 bp, 0.098 at 11-50 bp);
  `avi_splice_sites` 0.405. Against `alphagenome` (v0.7.0): 0.729 and 0.872.
- Strong under leave-one-gene-out at 3-50 bp, folds whose held-out ratio clears 18.7:
  avi 7/7, alphagenome 6/7, pangolin 6/7, spliceai 5/7, spliceai_walker 5/7.
- Strong under leave-one-gene-out at 11-50 bp: 0/6 for every tool, including the three
  that reach it in-sample.
- In-frame share among the variants a tool calls and the assay does not, 3-10 bp,
  SpliceAI: seven genes 40.1%, DDX3X 26.8%, TP53 9.5%.
- Pangolin on DDX3X under E3's leave-one-gene-out thresholds: Moderate in five of six
  evaluable cells, Strong in one (11-50 bp, FDR labels, LR 29.6).

### Corrections carried from the first round's internal review

Listed in that round's report; the tables here are the corrected ones. The only
correction whose direction is worth restating: the local likelihood ratio's numerator
is no longer floored, so BP4 Strong and Very strong are no longer unreachable by
construction, and BP4 Moderate is reachable at 3-10 bp.

### Correction to the round-2 entry above (2026-09-23)

The row "ClinVar arms on Walker LR+: tools where `classified` is highest" compared
arms by applying SpliceAI's 0.2 cut point to every column's raw score. That cut
point is defined on SpliceAI's delta score only; for AlphaGenome (range about 0.8
to 2.2), CADD or phyloP it selects nearly every variant and the ratio means
nothing. The comparison is withdrawn. It is replaced by `arms_at_tool_threshold.csv`,
which compares arms at each column's own fitted threshold, and at the published cut
point for the two SpliceAI columns only.

## Round 3 (2026-09-23)

- ClinVar arms compared at each column's own fitted threshold
  (`src/evid_arms.py`, `arms_at_tool_threshold.csv`); the direction differs by
  distance band, which the earlier pooled comparison could not show.
- The external gene DDX3X now carries the four AlphaGenome Atlas columns, merged
  from the same Atlas pass as the analysis set; the DDX3X merge and both external
  in-frame attributions are stages of `scripts/reproduce_evidence.py`.
- New stages: training-signal table (E9), fusion coefficients across folds (E10),
  thresholds refitted on every subset of training genes (E11), printed main tables
  (E12), printed supplementary tables (E13), figures (E14).
- The DDX3X assay was earlier described in a docstring as a different kind of
  readout. It is not: like most of the seven, it scores depletion from HAP1 cells.
- Printed thresholds keep as many figures as it takes to select the same variants as
  the fitted value; AlphaGenome's sit just under its ceiling of 2.2.
- The companion atlas is located through `evid_common.ATLAS_REPO`, defaulting to a
  checkout beside this repository.

## Round 4 (2026-09-23)

Changes that alter results, each with the reason:

- **DDX3X labels.** The deposit does publish a per-variant classification, in a
  separate score set (urn:mavedb:00000658-0-1), a random-forest call on the assay's
  fold-changes; it is now the primary label, as the seven genes use their deposits'
  own classifications. The control-anchored label is kept as a sensitivity label,
  ungated and gated: in score set q-1 (the final exon) the nonsense controls are not
  depleted (control AUROC 0.51), and that set supplied 37 of the label's 71 in-scope
  damaging calls. The earlier finding that the Strong thresholds did not transfer to
  DDX3X came from that set.
- **External intervals.** Every external likelihood ratio now carries a
  variant-level bootstrap interval within the gene and the damaging and normal counts
  in its band, and a tier is reported from the lower bound as well as the point
  estimate (`evid_external._band_detail`).
- **TP53 to |offset| 12.** The archived TP53 table stopped at 8 because of the
  earlier study's splice-window constant; `src/evid_tp53_extend.py` scores the 96
  deposited SNVs at offsets 9-12 with the same scorers and carries the 192 archived
  rows unchanged. It runs offline from cached scores as a stage of the entry point.
- **Fusion without the AlphaGenome splice score.** The AlphaGenome terms of service
  bar the use of its outputs to train other models (LICENSE-DATA), so the elastic net
  is refitted on the other seven panel columns (`evid_common.FUSION_FEATURES`). The
  single-column thresholds were rerun with it and came back identical.
- **Atlas combined score provenance.** It is a supervised model trained on gnomAD
  allele frequency, with AlphaMissense and conservation inputs, whose checkpoint was
  selected on four saturation genome editing datasets including the BRCA1, RAD51C and
  DDX3X assays used here. Its provenance row is corrected, and Table 3 counts its
  held-out folds over the other five genes.
- **Walker configuration.** `config/walker2023.yaml` now records the applied weight
  (Supporting, p. 1056) beside the calibrated strength (Moderate, p. 1051), and the
  BRCA1 check at a score of 0.5 (p. 1055); an earlier note told readers not to carry
  "supporting", which was wrong.
- New outputs: `depth_counts.csv` (labelled and damaging variants by intronic depth),
  terms-of-use columns in `predictor_training_provenance.csv` and Table 2, and an
  in-scope split in `clinvar_arm_no_assertion.csv`.
