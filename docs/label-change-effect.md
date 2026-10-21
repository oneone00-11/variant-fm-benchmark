# What relabelling PALB2 and BARD1 changed

Written by `scripts/stage_E16d/compare_runs.py`, comparing the evidence outputs of commit aab88bc (old) with the working tree of the next commit (new); do not edit by hand. The only change between the two is the label files of PALB2 and BARD1, now built from the investigator-provided ranges of the score sets Supplementary Table S1 cites (stage E16d). Full rows: `phase1/reports/evidence/label_change_effect.csv`.

## Tier changes

In sample, 3 cell(s) gained or lost a tier (predictor x stratum x tier x side):

| side | predictor | stratum | tier | change | threshold old -> new | LR there old -> new | held-out folds old -> new |
|---|---|---|---|---|---|---|---|
| PP3 | cadd | s11_50 | supporting | lost | 22.6 -> n/a | 5.09435 -> n/a | 3 -> 3 |
| BP4 | avi_splice_sites | s3_50 | supporting | gained | n/a -> 0.0102539 | n/a -> 0.164075 | 0 -> 0 |
| BP4 | avi_splice_site_usage | s3_10 | moderate | lost | 0.0136719 -> n/a | 0.136908 -> n/a | 0 -> 0 |

Held out, 1 cell(s) changed the number of folds that reach the tier:

| side | predictor | stratum | tier | folds old -> new |
|---|---|---|---|---|
| PP3 | nt | s3_50 | supporting | 0 -> 1 |

## The audit's prediction

Every cell of `label_audit_impact.csv` at aab88bc (the counts the closest rule would give), against the new counts; and the seven genes' totals. 44 of 44 cells met.

| gene | stratum | count | now | predicted | actual | met |
|---|---|---|---|---|---|---|
| BARD1 | pm12 | labelled | 119 | 119 | 119 | yes |
| BARD1 | pm12 | damaging | 119 | 119 | 119 | yes |
| BARD1 | pm12 | normal | 0 | 0 | 0 | yes |
| BARD1 | pm12 | excluded | 0 | 0 | 0 | yes |
| BARD1 | s3_10 | labelled | 459 | 459 | 459 | yes |
| BARD1 | s3_10 | damaging | 66 | 66 | 66 | yes |
| BARD1 | s3_10 | normal | 393 | 393 | 393 | yes |
| BARD1 | s3_10 | excluded | 14 | 14 | 14 | yes |
| BARD1 | s11_50 | labelled | 1,435 | 1,436 | 1,436 | yes |
| BARD1 | s11_50 | damaging | 21 | 19 | 19 | yes |
| BARD1 | s11_50 | normal | 1,414 | 1,417 | 1,417 | yes |
| BARD1 | s11_50 | excluded | 23 | 22 | 22 | yes |
| BARD1 | s3_50 | labelled | 1,894 | 1,895 | 1,895 | yes |
| BARD1 | s3_50 | damaging | 87 | 85 | 85 | yes |
| BARD1 | s3_50 | normal | 1,807 | 1,810 | 1,810 | yes |
| BARD1 | s3_50 | excluded | 37 | 36 | 36 | yes |
| PALB2 | pm12 | labelled | 141 | 141 | 141 | yes |
| PALB2 | pm12 | damaging | 139 | 139 | 139 | yes |
| PALB2 | pm12 | normal | 2 | 2 | 2 | yes |
| PALB2 | pm12 | excluded | 3 | 3 | 3 | yes |
| PALB2 | s3_10 | labelled | 411 | 445 | 445 | yes |
| PALB2 | s3_10 | damaging | 87 | 87 | 87 | yes |
| PALB2 | s3_10 | normal | 324 | 358 | 358 | yes |
| PALB2 | s3_10 | excluded | 90 | 56 | 56 | yes |
| PALB2 | s11_50 | labelled | 464 | 477 | 477 | yes |
| PALB2 | s11_50 | damaging | 15 | 15 | 15 | yes |
| PALB2 | s11_50 | normal | 449 | 462 | 462 | yes |
| PALB2 | s11_50 | excluded | 76 | 63 | 63 | yes |
| PALB2 | s3_50 | labelled | 875 | 922 | 922 | yes |
| PALB2 | s3_50 | damaging | 102 | 102 | 102 | yes |
| PALB2 | s3_50 | normal | 773 | 820 | 820 | yes |
| PALB2 | s3_50 | excluded | 166 | 119 | 119 | yes |
| all seven | s3_10 | labelled | 3,075 | 3,109 | 3,109 | yes |
| all seven | s3_10 | damaging | 570 | 570 | 570 | yes |
| all seven | s3_10 | normal | 2,505 | 2,539 | 2,539 | yes |
| all seven | s3_10 | excluded | 215 | 181 | 181 | yes |
| all seven | s11_50 | labelled | 4,537 | 4,551 | 4,551 | yes |
| all seven | s11_50 | damaging | 225 | 223 | 223 | yes |
| all seven | s11_50 | normal | 4,312 | 4,328 | 4,328 | yes |
| all seven | s11_50 | excluded | 170 | 156 | 156 | yes |
| all seven | s3_50 | labelled | 7,612 | 7,660 | 7,660 | yes |
| all seven | s3_50 | damaging | 795 | 793 | 793 | yes |
| all seven | s3_50 | normal | 6,817 | 6,867 | 6,867 | yes |
| all seven | s3_50 | excluded | 385 | 337 | 337 | yes |

## Count changes

The analysis set's counts that changed, by gene and stratum (`depth_counts.csv`, `set_counts.csv`):

| gene | stratum | count | old | new |
|---|---|---|---|---|
| BARD1 | 11-50 | labelled | 1,435 | 1,436 |
| BARD1 | 11-50 | damaging | 21 | 19 |
| BARD1 | 11-50 | normal | 1,414 | 1,417 |
| BARD1 | 3-50 | labelled | 1,894 | 1,895 |
| BARD1 | 3-50 | damaging | 87 | 85 |
| BARD1 | 3-50 | normal | 1,807 | 1,810 |
| PALB2 | 3-10 | labelled | 411 | 445 |
| PALB2 | 3-10 | normal | 324 | 358 |
| PALB2 | 11-50 | labelled | 464 | 477 |
| PALB2 | 11-50 | normal | 449 | 462 |
| PALB2 | 3-50 | labelled | 875 | 922 |
| PALB2 | 3-50 | normal | 773 | 820 |
| all | 3-10 | labelled | 3,075 | 3,109 |
| all | 3-10 | normal | 2,505 | 2,539 |
| all | 11-50 | labelled | 4,537 | 4,551 |
| all | 11-50 | damaging | 225 | 223 |
| all | 11-50 | normal | 4,312 | 4,328 |
| all | 3-50 | labelled | 7,612 | 7,660 |
| all | 3-50 | damaging | 795 | 793 |
| all | 3-50 | normal | 6,817 | 6,867 |

## The ten largest threshold moves

Ranked by the labelled variants of the new analysis set that lie between the old and the new threshold (the fusion, which the set does not carry, is left out):

| side | predictor | stratum | tier | threshold old -> new | variants between |
|---|---|---|---|---|---|
| BP4 | avi | s3_10 | supporting | 0.041665 -> 0.191895 | 642 |
| BP4 | avi | s3_50 | supporting | -0.0669537 -> -0.0517293 | 129 |
| BP4 | spliceai | all_1_50 | supporting | 0.00184184 -> 0.00170457 | 72 |
| PP3 | avi | s3_50 | moderate | 0.824646 -> 0.892125 | 67 |
| BP4 | spliceai | all_1_50 | moderate | 0.000335217 -> 0.00031347 | 24 |
| PP3 | spliceai | s3_10 | supporting | 0.375475 -> 0.397732 | 24 |
| PP3 | avi_splice_site_usage | s3_10 | moderate | 0.678711 -> 0.656784 | 14 |
| PP3 | alphagenome | s3_10 | moderate | 2.19952 -> 2.19955 | 9 |
| PP3 | avi_splice_site_usage | s3_10 | supporting | 0.539062 -> 0.568359 | 9 |
| BP4 | alphagenome | s3_50 | supporting | 1.24515 -> 1.2469 | 8 |
