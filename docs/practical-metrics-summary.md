# Practical metrics at the evidence thresholds

Written by `phase1/src/evid_practical_metrics.py` (stage E15 of `scripts/reproduce_evidence.py`) from the tables it writes to `phase1/reports/evidence/`; do not edit by hand. The thresholds are the study's own: nothing here is fitted.

## Sensitivity, specificity and coverage at the fitted thresholds

Each cell gives the held-out value first and the in-sample value in brackets, for the seven genes.

| Predictor | Band | Supporting: sensitivity | Supporting: specificity | Supporting: coverage | Moderate: sensitivity | Moderate: specificity | Moderate: coverage | Strong: sensitivity | Strong: specificity | Strong: coverage | Held-out folds reaching S, M, St | Pooled AUROC | Pooled AUPRC |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| SpliceAI (published basis) | 3–10 bp | 0.61 (0.67) | 0.95 (0.95) | 0.15 (0.16) | 0.56 (0.57) | 0.97 (0.97) | 0.13 (0.13) | 0.21 (0.08) | 0.99 (>0.99) | 0.04 (0.02) | 7, 7, 7 of 7 | 0.90 | 0.73 |
| SpliceAI (published basis) | 11–50 bp | 0.40 (0.41) | 0.97 (0.97) | 0.05 (0.05) | 0.35 (0.37) | 0.98 (0.98) | 0.04 (0.04) | – (–) | – (–) | – (–) | 6, 6, 1 of 6 | 0.79 | 0.36 |
| SpliceAI (published basis) | 3–50 bp | 0.61 (0.63) | 0.96 (0.95) | 0.10 (0.11) | 0.54 (0.55) | 0.97 (0.97) | 0.08 (0.08) | 0.19 (0.20) | >0.99 (>0.99) | 0.02 (0.02) | 7, 7, 7 of 7 | 0.88 | 0.63 |
| Pangolin | 3–10 bp | 0.67 (0.69) | 0.95 (0.95) | 0.16 (0.17) | 0.54 (0.58) | 0.97 (0.97) | 0.12 (0.13) | 0.38 (0.39) | 0.99 (0.99) | 0.07 (0.07) | 7, 7, 7 of 7 | 0.91 | 0.76 |
| Pangolin | 11–50 bp | 0.41 (0.39) | 0.96 (0.97) | 0.06 (0.05) | 0.27 (0.26) | 0.99 (0.99) | 0.02 (0.02) | – (–) | – (–) | – (–) | 6, 6, 0 of 6 | 0.81 | 0.39 |
| Pangolin | 3–50 bp | 0.59 (0.60) | 0.97 (0.97) | 0.09 (0.09) | 0.55 (0.56) | 0.98 (0.98) | 0.08 (0.08) | 0.34 (0.33) | >0.99 (>0.99) | 0.04 (0.04) | 7, 7, 7 of 7 | 0.89 | 0.66 |
| AlphaGenome splice score | 3–10 bp | 0.67 (0.68) | 0.95 (0.95) | 0.16 (0.16) | 0.59 (0.59) | 0.97 (0.97) | 0.13 (0.13) | 0.39 (0.41) | 0.99 (0.99) | 0.08 (0.08) | 7, 7, 7 of 7 | 0.90 | 0.75 |
| AlphaGenome splice score | 11–50 bp | 0.37 (0.41) | 0.97 (0.97) | 0.04 (0.05) | 0.36 (0.31) | 0.98 (0.99) | 0.04 (0.03) | – (0.25) | – (>0.99) | – (0.02) | 6, 6, 1 of 6 | 0.81 | 0.41 |
| AlphaGenome splice score | 3–50 bp | 0.58 (0.59) | 0.97 (0.97) | 0.09 (0.09) | 0.54 (0.53) | 0.98 (0.98) | 0.07 (0.07) | 0.35 (0.35) | 0.99 (0.99) | 0.04 (0.04) | 7, 7, 7 of 7 | 0.89 | 0.66 |
| Atlas combined score † | 3–10 bp | 0.61 (0.68) | 0.95 (0.95) | 0.14 (0.17) | 0.46 (0.52) | 0.98 (0.98) | 0.09 (0.11) | 0.24 (0.27) | >0.99 (>0.99) | 0.04 (0.05) | 5, 5, 5 of 5 | 0.91 | 0.79 |
| Atlas combined score † | 11–50 bp | 0.36 (0.35) | 0.97 (0.98) | 0.05 (0.04) | 0.29 (0.28) | 0.98 (0.99) | 0.03 (0.03) | – (0.11) | – (>0.99) | – (0.01) | 4, 4, 1 of 4 | 0.79 | 0.37 |
| Atlas combined score † | 3–50 bp | 0.55 (0.62) | 0.97 (0.96) | 0.09 (0.10) | 0.48 (0.55) | 0.98 (0.98) | 0.07 (0.08) | 0.28 (0.35) | >0.99 (>0.99) | 0.03 (0.04) | 5, 5, 5 of 5 | 0.88 | 0.68 |

- Held out: each fold's threshold, fitted on the other genes (six, or five at 11–50 bp, where BRCA2 has no variant), applied to the gene left out; the counts are added over the folds whose threshold reached the tier, not averaged. In brackets: the threshold fitted on all the genes, applied to them.
- Held-out folds reaching S, M, St: for each tier, the folds whose threshold reached it, of the genes held out. Those folds include ones whose band holds too few variants for a likelihood ratio, which Table 3 counts as not evaluable, so the count is not Table 3's k/n.
- A held-out value is shown when the tier is reached in-sample (the rule E7 applies before carrying a threshold to another gene) and in at least half the held-out folds; otherwise it would rest on one or two genes, and only the fold count is given. A dash in brackets: the tier is not reached in-sample.
- The Atlas combined score (†): its model was selected on the BRCA1 and RAD51C assays, so its held-out values leave those two genes out. Its in-sample values, in brackets, include them.
- Sensitivity: the share of damaging variants at or above the threshold. Specificity: the share of normal variants below it. Coverage: the share of all scored variants in the band, labelled or not, at or above it; these are the variants a user would get the evidence for. The share among labelled variants only, the denominator of the manuscript's '14.8% of the variants' at 0.2, is `coverage_labelled`.
- Pooled AUROC and AUPRC: the seven genes pooled into one curve, the way one threshold is used across genes. The manuscript's AUROC (Methods) is gene-pooled instead: per-gene AUROCs combined on the logit scale with random effects. For all variants it is `auroc_gene_pooled` in `practical_curve_auc.csv`, beside the pooled value and a split of the difference.
- Two decimals; >0.99 marks a value from 0.995 up to, not including, 1.
- Wilson intervals, the BP4 side, the other ten columns and the fusion are in `practical_metrics_by_threshold.csv`; each fold in `practical_metrics_logo_by_gene.csv`; DDX3X and TP53 in `practical_metrics_external.csv`. The curve tables start at threshold +inf, where nothing is called; the precision-recall table's first row is the conventional (recall 0, precision 1), not an observation.

## In five sentences

1. At 11–50 bp each ROC curve lies below the same tool's 3–10 bp curve except near its two ends; the best pooled AUROC there is 0.81, below the lowest at 3–10 bp (0.90).
2. Held out, a Strong threshold at 3–10 bp catches 21% to 39% of damaging variants, depending on the tool.
3. The cost is coverage: Strong evidence reaches no more than 8% of the variants in that band.
4. At a specificity of 0.95, the AlphaGenome splice score is the most sensitive, or tied for it, in all three bands, though by no more than 0.01.
5. Carried to DDX3X, the fitted thresholds are more sensitive than in the held-out seven genes, by as much as 0.38, and less specific, by as much as 0.05; on TP53 sensitivity rises too, except for the Atlas combined score's Moderate and Strong.

## Figure legend

`phase1/reports/evidence/figures/figure_practical.pdf`. ROC curves (a) and precision-recall curves (b) of four predictors, the seven genes pooled, by distance from the exon boundary. Colours are those of Supplementary Figure S2, and each predictor also has its own line pattern. Filled markers are the fitted in-sample thresholds: circle Supporting, square Moderate, triangle Strong; a tier that is not reached is not drawn, and markers of different predictors are drawn slightly apart so that coinciding ones stay visible. The open diamond is the published cut point, 0.2, on SpliceAI (published basis). Insets in (a) enlarge the high-specificity end, 1 − specificity up to 0.06 against sensitivity, where every fitted threshold lies. The key in each ROC panel gives the pooled AUROC; the dashed line in (b) is the band's prevalence. † The Atlas combined score's model was selected on the BRCA1 and RAD51C assays, which these curves include.
