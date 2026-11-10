"""E8.13: docs/manuscript_bindings_v6.tsv, the binding table of the v6 text, as
scripts/stage_E5/build_bindings_v5.py built v5's (which stays as it is).

The texts. Written here from the manuscript and the abstract of the F3 file (both read
only) into a folder outside the repository (--texts), each paragraph on one line, in the
form the number gate reads (scripts/check_evidence_manuscript_numbers.py):
- v6_main.md: "# Abstract" and the F3 file's abstract; the manuscript's main text, its
  figure legends among it, from the first heading after the title to "References"; then
  the tables' titles and notes (each "Table N." paragraph and the parts' titles); then
  "# References", where the gate stops reading;
- v6_supplement.md: the supplementary notes, from "Supplementary Notes" to "References
  cited in the Supplementary Notes";
- the same two with a source tag after each sentence that holds a measured number,
  naming the files its numbers are bound to (v6_main_keyed.md, v6_supplement_keyed.md),
  for the gate's bound mode and self-test. List items are written without their numbers
  (the gate does not read a list's numbering).

The table. One row per number the gate reads, in reading order:
- a measured number gets the file, rows, column and derivation of its cell, below
  (BINDINGS), each keyed by the sentence it stands in (a piece of it), the number as
  printed and its place among the numbers so printed in that sentence;
- a constant gets the reason the gate gives it, or, where the gate's reason does not fit
  the sentence or the gate does not declare the number, the reason below (REASONS);
- a number the gate declares that is a measurement in its sentence (the m of a fold
  count "k/n (m)", a count read as a band's bound) is bound, as any measured number.

Checks before writing; any failure stops it:
- every number of the two texts has one entry, and every entry one number;
- every excerpt is a piece of its sentence as the gate reads it;
- the texts are the ones the manuscript and the F3 file give (their sha256 recorded).

    .venv/bin/python scripts/stage_E8/build_bindings_v6.py \\
        --manuscript ~/Desktop/calibration/CALInewversion.docx \\
        --f3 ~/Desktop/calibration/F3_摘要标题与声明_草稿_20261108.docx \\
        --texts ~/Desktop/整理文件/E8b_工作文件/v6_texts
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog  # noqa: E402

sys.path.insert(0, str(REPO / "scripts"))
import check_evidence_manuscript_numbers as gate  # noqa: E402

OUT = REPO / "docs" / "manuscript_bindings_v6.tsv"
R = "phase1/reports/evidence/"
WT, SC, LC = R + "walker_thresholds.csv", R + "set_counts.csv", R + "walker_cutpoint_local_vs_band.csv"
DC, DS = R + "depth_counts.csv", R + "threshold_dilution_summary.csv"
C5 = R + "suppl_analyses/c5_fold_counts_point_vs_bound.csv"
C1, C1L = R + "suppl_analyses/c1_clinvar_reference_metrics.csv", R + "suppl_analyses/c1_low_circularity_metrics.csv"
C2 = R + "suppl_analyses/c2_clinvar_coverage_by_stratum.csv"
# two sources carry comment lines above their header, which the gate's reader does not
# skip: their numbers are bound to the supplementary tables the layout prints from them
S28, S31 = R + "layout/supplementary/inframe_true_positives.csv", R + "layout/supplementary/bp4_at_bp7_positions.csv"
PM, ES = R + "practical_metrics_by_threshold.csv", R + "external_summary.csv"
T4B, T4C = R + "tables/table4b_fitted_threshold_highest_arm.csv", R + "tables/table4c_auroc_highest_arm.csv"
S11 = R + "supplement/tableS11_within_gene_arms_summary.csv"
S3P = R + "supplement/tableS3_predictor_columns.csv"
XX2, CB = R + "external_xrcc2.csv", R + "external_ddx3x_column_basis.csv"
VL = R + "layout/variants/analysed_variants.csv"
YAML = "phase1/config/walker2023.yaml"
DDX = "phase1/data/evidence/external/ddx3x_splice.manifest.json"
ANCH = "phase1/data/evidence/external/e18_anchors_cadd.parquet.provenance.json"
AGD = "phase1/data/evidence/alphagenome_official/drift_2b_start_scores.tsv.provenance.json"
Q1 = '["urn:mavedb:00000658-q-1"]'


def b(file, column, rows="", derive="value", scale=1):
    return dict(file=file, rows=rows, column=column, derive=derive, scale=scale)


def wt(stratum, column, arm="all", scale=1):
    return b(WT, column, f"score_column=spliceai_walker; stratum={stratum}; clinvar_arm={arm}; scope=pooled",
             scale=scale)


def c5(tools, stratum, tier, column):
    """A fold count of Supplementary Table S20(b): k (point_k, bound_k), n (n_folds), m
    (point_m, bound_m); several predictors or bands at once, each the same."""
    tool = f"tool in {tools}" if "|" in tools else f"tool={tools}"
    st = f"stratum in {stratum}" if "|" in stratum else f"stratum={stratum}"
    return b(C5, column, f"{tool}; {st}; tier={tier}", derive="all" if "|" in tools + stratum else "value")


def es(stratum, tool, side, tier):
    """An external summary's cell "k/n (m)", read as text: the number is one it prints."""
    return b(ES, "External genes reaching the tier, k/n (m)",
             f"Gene set=Seven external genes; Stratum (intronic offset, bp)={stratum}; Tool={tool}; "
             f"Side={side}; Tier={tier}", derive="text")


def c1(tool, threshold_kind, column, stratum="s3_10", file=C1):
    return b(file, column, f"star_level=stars_ge_1; analysis=main; measure=threshold; tool={tool}; "
                           f"stratum={stratum}; threshold_kind={threshold_kind}")


def band(gene):
    """How deep an external gene's deposit reaches: its variants' largest signed offset in
    the list of analysed variants (on the donor side for SFPQ and TP53)."""
    return b(VL, "Signed intronic offset", f"Gene={gene}", derive="max")


FOLDS = "Supplementary Table S20(b)"
# Each measured number: (document, a piece of its sentence, the number as printed, its
# rank among the numbers so printed in that sentence) -> its cell.
BINDINGS: list[tuple[str, str, str, int, dict]] = [
    # the abstract
    ("main", "the published PP3 cut point was 8.27", "8.27", 1, wt("s3_50", "lr_pp3")),
    ("main", "the published PP3 cut point was 8.27", "6.36", 1, wt("s3_50", "lr_pp3_lo")),
    ("main", "the published PP3 cut point was 8.27", "11.1", 1, wt("s3_50", "lr_pp3_hi")),
    # introduction
    ("main", "a score of 0.5 or above gave Strong evidence", "0.5", 1,
     b(YAML, "thresholds.applied_weight.higher_score_check_cut_point")),
    ("main", "The analysis set comprises 8,853", "8,853", 1, b(SC, "n", derive="sum")),
    ("main", "The analysis set comprises 8,853", "13", 1, b(S3P, "Score column", derive="count")),
    # materials and methods
    ("main", "Its 8,853 variants were taken", "8,853", 1, b(SC, "n", derive="sum")),
    ("main", "We evaluated 13 predictor scores", "13", 1, b(S3P, "Score column", derive="count")),
    # results: the published cut points
    ("main", "gave a band likelihood ratio of 8.27", "8.27", 1, wt("s3_50", "lr_pp3")),
    ("main", "gave a band likelihood ratio of 8.27", "6.36", 1, wt("s3_50", "lr_pp3_lo")),
    ("main", "gave a band likelihood ratio of 8.27", "11.1", 1, wt("s3_50", "lr_pp3_hi")),
    ("main", "The ratio was 5.13 in the proximal band", "5.13", 1, wt("s3_10", "lr_pp3")),
    ("main", "The ratio was 5.13 in the proximal band", "10.4", 1, wt("s11_50", "lr_pp3")),
    ("main", "At the BP4 cut point the band ratio was 0.172", "0.172", 1, wt("s3_10", "lr_bp4")),
    ("main", "the ratio at the PP3 cut point was 1.03", "1.03", 1, wt("pm12", "lr_pp3")),
    ("main", "local likelihood ratio was 0.85 in the proximal band", "0.85", 1,
     b(LC, "local_lr_at_cut", "tool=spliceai_walker; stratum=s3_10")),
    ("main", "local likelihood ratio was 0.85 in the proximal band", "1.77", 1,
     b(LC, "local_lr_at_cut", "tool=spliceai_walker; stratum=s3_50")),
    # results: tiers held out (Supplementary Table S20(b): k/n (m) by the point estimate,
    # and by the lower bound)
    ("main", "Pangolin cleared it in 6/7 (1) folds and SpliceAI and the AlphaGenome", "6", 1,
     c5("pangolin", "s3_50", "strong", "point_k")),
    ("main", "Pangolin cleared it in 6/7 (1) folds and SpliceAI and the AlphaGenome", "7", 1,
     c5("pangolin", "s3_50", "strong", "n_folds")),
    ("main", "Pangolin cleared it in 6/7 (1) folds and SpliceAI and the AlphaGenome", "1", 1,
     c5("pangolin", "s3_50", "strong", "point_m")),
    ("main", "Pangolin cleared it in 6/7 (1) folds and SpliceAI and the AlphaGenome", "5", 1,
     c5("spliceai_walker|alphagenome", "s3_50", "strong", "point_k")),
    ("main", "Pangolin cleared it in 6/7 (1) folds and SpliceAI and the AlphaGenome", "7", 2,
     c5("spliceai_walker|alphagenome", "s3_50", "strong", "n_folds")),
    ("main", "Pangolin cleared it in 6/7 (1) folds and SpliceAI and the AlphaGenome", "1", 2,
     c5("spliceai_walker|alphagenome", "s3_50", "strong", "point_m")),
    ("main", "In the proximal band Pangolin and AlphaGenome cleared it", "6", 1,
     c5("pangolin|alphagenome", "s3_10", "strong", "point_k")),
    ("main", "In the proximal band Pangolin and AlphaGenome cleared it", "7", 1,
     c5("pangolin|alphagenome", "s3_10", "strong", "n_folds")),
    ("main", "In the proximal band Pangolin and AlphaGenome cleared it", "1", 1,
     c5("pangolin|alphagenome", "s3_10", "strong", "point_m")),
    ("main", "In the proximal band Pangolin and AlphaGenome cleared it", "4", 1,
     c5("spliceai_walker", "s3_10", "strong", "point_k")),
    ("main", "In the proximal band Pangolin and AlphaGenome cleared it", "7", 2,
     c5("spliceai_walker", "s3_10", "strong", "n_folds")),
    ("main", "In the proximal band Pangolin and AlphaGenome cleared it", "1", 2,
     c5("spliceai_walker", "s3_10", "strong", "point_m")),
    ("main", "cleared it in 4/6 (2) folds each", "4", 1, c5("spliceai_walker|pangolin", "s11_50", "moderate", "point_k")),
    ("main", "cleared it in 4/6 (2) folds each", "6", 1, c5("spliceai_walker|pangolin", "s11_50", "moderate", "n_folds")),
    ("main", "cleared it in 4/6 (2) folds each", "2", 1, c5("spliceai_walker|pangolin", "s11_50", "moderate", "point_m")),
    ("main", "Most of the distal band’s 223 damaging variants", "223", 1, wt("s11_50", "n_pos")),
    ("main", "SpliceAI’s Strong count fell to 3/7 (1)", "3", 1, c5("spliceai_walker", "s3_10|s3_50", "strong", "bound_k")),
    ("main", "SpliceAI’s Strong count fell to 3/7 (1)", "7", 1, c5("spliceai_walker", "s3_10|s3_50", "strong", "n_folds")),
    ("main", "SpliceAI’s Strong count fell to 3/7 (1)", "1", 1, c5("spliceai_walker", "s3_10|s3_50", "strong", "bound_m")),
    ("main", "SpliceAI’s Strong count fell to 3/7 (1)", "4", 1, c5("pangolin", "s3_10", "strong", "bound_k")),
    ("main", "SpliceAI’s Strong count fell to 3/7 (1)", "7", 2, c5("pangolin", "s3_10", "strong", "n_folds")),
    ("main", "SpliceAI’s Strong count fell to 3/7 (1)", "1", 2, c5("pangolin", "s3_10", "strong", "bound_m")),
    # the share of held-out pairs clearing Strong over subsets of the training genes
    ("main", "proximal Strong threshold cleared at least 95%", "95%", 1,
     b(DS, "frac_pairs_cleared", "tool=pangolin; stratum=s3_10; tier=strong", "min", 100)),
    ("main", "proximal Strong threshold cleared at least 95%", "55", 1,
     b(DS, "frac_pairs_cleared", "tool=spliceai_walker; stratum=s3_10; tier=strong", "min", 100)),
    ("main", "proximal Strong threshold cleared at least 95%", "67%", 1,
     b(DS, "frac_pairs_cleared", "tool=spliceai_walker; stratum=s3_10; tier=strong", "max", 100)),
    # held-out performance in the proximal band (Table 4(a), Supplementary Table S17)
    ("main", "caught about half of the damaging variants at a specificity of 0.97", "0.97", 1,
     b(PM, "specificity", "tool in spliceai_walker|pangolin; stratum=s3_10; threshold_kind=moderate; "
                          "threshold_source=logo", "all")),
    ("main", "The published PP3 cut point caught 0.79", "0.79", 1,
     b(PM, "sensitivity", "tool=spliceai_walker; stratum=s3_10; threshold_kind=walker_pp3")),
    ("main", "The published PP3 cut point caught 0.79", "0.85", 1,
     b(PM, "specificity", "tool=spliceai_walker; stratum=s3_10; threshold_kind=walker_pp3")),
    # results: ClinVar record status
    ("main", "ClinVar has no record of 52%", "52%", 1,
     b(C2, "share_unrecorded", "gene=Seven genes; stratum=s3_10", scale=100)),
    ("main", "ClinVar has no record of 52%", "77%", 1,
     b(C2, "share_unrecorded", "gene=Seven genes; stratum=s11_50", scale=100)),
    ("main", "the band ratio was 16.3 among variants", "16.3", 1, wt("s3_50", "lr_pp3", arm="classified")),
    ("main", "the band ratio was 16.3 among variants", "5.0", 1, wt("s3_50", "lr_pp3", arm="recorded_unclassified")),
    ("main", "the classified arm gave the highest ratio in 13/13 columns", "13", 1,
     b(T4B, "Highest: Classified", "Genes=Seven genes; Band=3–10 bp")),
    ("main", "the classified arm gave the highest ratio in 13/13 columns", "13", 2,
     b(T4B, "Columns with all three arms estimable", "Genes=Seven genes; Band=3–10 bp")),
    ("main", "the classified arm gave the highest ratio in 13/13 columns", "11", 1,
     b(T4B, "Highest: Classified", "Genes=Seven genes; Band=3–50 bp")),
    ("main", "the classified arm gave the highest ratio in 13/13 columns", "11", 2,
     b(T4B, "Columns with all three arms estimable", "Genes=Seven genes; Band=3–50 bp")),
    ("main", "it had the lowest point estimate in 10/13 columns", "10", 1,
     b(T4C, "Lowest AUROC: Classified", "Pool=3–50 bp; Genes=Seven genes")),
    ("main", "it had the lowest point estimate in 10/13 columns", "13", 1,
     b(T4C, "Columns", "Pool=3–50 bp; Genes=Seven genes")),
    ("main", "it had the lowest point estimate in 10/13 columns", "3", 1,
     b(T4C, "Highest AUROC: Classified", "Pool=3–50 bp; Genes=Seven genes")),
    ("main", "it had the lowest point estimate in 10/13 columns", "13", 2,
     b(T4C, "Columns", "Pool=3–50 bp; Genes=Seven genes")),
    ("main", "the unrecorded arm had the highest AUROC in 34/52", "34", 1,
     b(S11, "Highest in: Unrecorded", "Stratum (intronic offset, bp)=3–50; Gene startswith All 4 genes")),
    ("main", "the unrecorded arm had the highest AUROC in 34/52", "52", 1,
     b(S11, "Columns compared", "Stratum (intronic offset, bp)=3–50; Gene startswith All 4 genes")),
    ("main", "the unrecorded arm had the highest AUROC in 34/52", "8", 1,
     b(S11, "Highest in: Classified (P/LP or B/LB)", "Stratum (intronic offset, bp)=3–50; Gene startswith All 4 genes")),
    ("main", "the unrecorded arm had the highest AUROC in 34/52", "52", 2,
     b(S11, "Columns compared", "Stratum (intronic offset, bp)=3–50; Gene startswith All 4 genes")),
    ("main", "all 74 pathogenic or likely pathogenic variants in the proximal band", "74", 1,
     c1("spliceai_walker", "walker_pp3", "n_pos_above")),
    ("main", "The same held for all 12 whose assertion", "12", 1,
     c1("spliceai_walker", "walker_pp3", "n_pos_above", file=C1L)),
    ("main", "the share fell to 11/74 for SpliceAI", "11", 1, c1("spliceai_walker", "strong", "n_pos_above")),
    ("main", "the share fell to 11/74 for SpliceAI", "74", 1, c1("spliceai_walker", "strong", "n_plp")),
    ("main", "the share fell to 11/74 for SpliceAI", "52", 1, c1("pangolin", "strong", "n_pos_above")),
    ("main", "the share fell to 11/74 for SpliceAI", "74", 2, c1("pangolin", "strong", "n_plp")),
    # results: external genes
    ("main", "XRCC2, with 7 damaging variants there", "7", 1,
     b(XX2, "n_pos", "label_definition=y_investigator_range; tool=spliceai_walker; stratum=s3_10")),
    ("main", "their Strong thresholds cleared the Strong boundary in 6/7 (1) genes", "6", 1,
     es("3–50", "Pangolin", "PP3", "Strong")),
    ("main", "their Strong thresholds cleared the Strong boundary in 6/7 (1) genes", "7", 1,
     es("3–50", "Pangolin", "PP3", "Strong")),
    ("main", "their Strong thresholds cleared the Strong boundary in 6/7 (1) genes", "1", 1,
     es("3–50", "Pangolin", "PP3", "Strong")),
    ("main", "their Strong thresholds cleared the Strong boundary in 6/7 (1) genes", "4", 1,
     es("3–50", "SpliceAI (published basis)", "PP3", "Strong")),
    ("main", "their Strong thresholds cleared the Strong boundary in 6/7 (1) genes", "7", 2,
     es("3–50", "SpliceAI (published basis)", "PP3", "Strong")),
    ("main", "their Strong thresholds cleared the Strong boundary in 6/7 (1) genes", "3", 1,
     es("3–50", "SpliceAI (published basis)", "PP3", "Strong")),
    ("main", "the Strong boundary was cleared in 2/7 (1) genes for Pangolin", "2", 1,
     es("3–10", "Pangolin", "PP3", "Strong")),
    ("main", "the Strong boundary was cleared in 2/7 (1) genes for Pangolin", "7", 1,
     es("3–10", "Pangolin", "PP3", "Strong")),
    ("main", "the Strong boundary was cleared in 2/7 (1) genes for Pangolin", "1", 1,
     es("3–10", "Pangolin", "PP3", "Strong")),
    ("main", "the Strong boundary was cleared in 2/7 (1) genes for Pangolin", "1", 2,
     es("3–10", "SpliceAI (published basis)", "PP3", "Strong")),
    ("main", "the Strong boundary was cleared in 2/7 (1) genes for Pangolin", "7", 2,
     es("3–10", "SpliceAI (published basis)", "PP3", "Strong")),
    ("main", "the Strong boundary was cleared in 2/7 (1) genes for Pangolin", "6", 1,
     es("3–10", "SpliceAI (published basis)", "PP3", "Strong")),
    ("main", "these BP4 thresholds cleared the Supporting boundary in 4/7 (1)", "4", 1,
     es("3–10", "SpliceAI (published basis)", "BP4", "Supporting")),
    ("main", "these BP4 thresholds cleared the Supporting boundary in 4/7 (1)", "7", 1,
     es("3–10", "SpliceAI (published basis)", "BP4", "Supporting")),
    ("main", "these BP4 thresholds cleared the Supporting boundary in 4/7 (1)", "1", 1,
     es("3–10", "SpliceAI (published basis)", "BP4", "Supporting")),
    ("main", "these BP4 thresholds cleared the Supporting boundary in 4/7 (1)", "5", 1,
     es("3–10", "Pangolin", "BP4", "Supporting")),
    ("main", "these BP4 thresholds cleared the Supporting boundary in 4/7 (1)", "7", 2,
     es("3–10", "Pangolin", "BP4", "Supporting")),
    ("main", "these BP4 thresholds cleared the Supporting boundary in 4/7 (1)", "1", 2,
     es("3–10", "Pangolin", "BP4", "Supporting")),
    # discussion
    ("main", "lower than the 16.0 they report", "16.0", 1, b(YAML, "thresholds.pp3.reported_lr")),
    ("main", "with a ratio of 2.52 for scores of 0.20–0.22", "2.52", 1,
     b(YAML, "thresholds.applied_weight.narrow_band_lr")),
    ("main", "with a ratio of 2.52 for scores of 0.20–0.22", "0.20", 1,
     b(YAML, "thresholds.applied_weight.narrow_band_scores[0]")),
    ("main", "with a ratio of 2.52 for scores of 0.20–0.22", "0.22", 1,
     b(YAML, "thresholds.applied_weight.narrow_band_scores[1]")),
    ("main", "In the proximal band, 40.5% of the variants above the cut point", "40.5%", 1,
     b(S28, "In-frame share", "Gene set=Seven genes; Events read from=SpliceAI (published basis); "
                              "Stratum (intronic offset, bp)=3–10; True or false positive=False positives", scale=100)),
    ("main", "In the proximal band, 40.5% of the variants above the cut point", "29.8%", 1,
     b(S28, "In-frame share", "Gene set=Seven genes; Events read from=SpliceAI (published basis); "
                              "Stratum (intronic offset, bp)=3–10; True or false positive=True positives", scale=100)),
    ("main", "nonsense controls were not depleted (control AUROC 0.51)", "0.51", 1,
     b(DDX, f"controls.per_scoreset{Q1}.control_auroc")),
    ("main", "1.6% of the variants scoring at or below the BP4 cut point were damaging", "1.6%", 1,
     b(S31, "Damaging share of variants meeting the condition",
       "Condition startswith (a); Splice-site side=both sides; Gene=Seven genes", scale=100)),
    # the figures' legends: the classes in each band (Figure 1 and Figure 3)
    *[("main", f"Damaging and normal variants: 570 and 2,539{tail}", printed, 1, wt(st, col))
      for tail in (" at 3–10 bp, 223 and 4,328 at 11–50 bp, and 793 and 6,867 at 3–50 bp; the counts",
                   " at 3–10 bp, 223 and 4,328 at 11–50 bp, and 793 and 6,867 at 3–50 bp. †")
      for printed, st, col in (("570", "s3_10", "n_pos"), ("2,539", "s3_10", "n_neg"), ("223", "s11_50", "n_pos"),
                               ("4,328", "s11_50", "n_neg"), ("793", "s3_50", "n_pos"), ("6,867", "s3_50", "n_neg"))],
    ("main", "Each gene’s bands stop at the depth its deposit reaches: 25 bp in DDX3X", "25", 1,
     b(DDX, "max_intron_offset")),
    ("main", "Each gene’s bands stop at the depth its deposit reaches: 25 bp in DDX3X", "44", 1, band("SFPQ")),
    ("main", "Each gene’s bands stop at the depth its deposit reaches: 25 bp in DDX3X", "12", 1, band("TP53")),
    # the tables' notes
    ("main", "All 13 predictors scored every variant", "13", 1, b(S3P, "Score column", derive="count")),
    ("main", "The 856 variants at the canonical dinucleotides", "856", 1, b(SC, "n", "stratum=pm12", "sum")),
    ("main", "BRCA2’s deposit stops at 10 bp", "10", 1, b(DC, "max_offset_in_gene", "gene=BRCA2", "max")),
    ("main", "The distal band, with 3 P/LP variants, is left out", "3", 1,
     c1("spliceai_walker", "walker_pp3", "n_plp", stratum="s11_50")),
    # the supplementary notes
    ("supplement", "39 anchor variants for the locally run scores", "39", 1, b(ANCH, "same_variable.n_anchors")),
    ("supplement", "39 anchor variants for the locally run scores", "10", 1, b(AGD, "n_variants")),
    ("supplement", "In each distance band, the 13 predictors", "13", 1, b(S3P, "Score column", derive="count")),
    ("supplement", "nonsense controls were not depleted (control AUROC 0.51)", "0.51", 1,
     b(DDX, f"controls.per_scoreset{Q1}.control_auroc")),
    ("supplement", "which therefore uses 12 scores", "12", 1, b(CB, "tool", "comparable=True", "count")),
    ("supplement", "12 bp in TP53 and 25 bp in DDX3X", "12", 1, band("TP53")),
    ("supplement", "12 bp in TP53 and 25 bp in DDX3X", "25", 1, b(DDX, "max_intron_offset")),
]

# Constants the gate does not declare, or declares for a reason the sentence does not
# bear: (document, a piece of the sentence, the number as printed, rank) -> the reason
BAND = "a distance band's bound, fixed by the design (3-10, 11-50, 3-50 bp; config/walker2023.yaml strata)"
POSITION = "a position group's bound, fixed by the design (src/evid_suppl_c6_c7.py)"
BP7 = "the BP7 positions of the recommendation, donor +7 and beyond, acceptor -21 and beyond (Walker et al.)"
REASONS: list[tuple[str, str, str, int, str]] = [
    ("main", "Confidence intervals (95%) for the training genes", "95%", 1, "confidence level"),
    ("main", "The BP4 cut point gave Moderate evidence within 10 bp", "10", 1, BAND),
    ("main", "split at 10 bp into a proximal and a distal band", "10", 1, BAND),
    ("main", "A ratio this close to 1 at the cut point", "1", 1, "a likelihood ratio of 1, no evidence either way"),
    ("main", "We have shown that the published PP3 cut point holds at Moderate", "10", 1, BAND),
    ("main", "Beyond 10 bp, Moderate is the most", "10", 1, BAND),
    ("main", "at a prior of 0.05 SpliceAI and Pangolin lose Strong", "0.05", 1,
     "a prior of the sensitivity analysis, fixed in advance (src/evid_suppl_c4.py)"),
    ("supplement", "Ensembl release 112 annotation", "112", 1, "the Ensembl release of Pangolin's annotation (its provenance)"),
    ("supplement", "variant-level bootstrap of 10,000 draws", "10,000", 1,
     "variant-level bootstrap draws, Pejaver's count (src/evid_interval_lr.py N_BOOT_VARIANT)"),
    ("supplement", "Each variant is scored over a 1 Mb window", "1", 1,
     "the scoring window of the developers' notebook, 1 Mb (data/evidence/alphagenome_official/NOTEBOOK_SOURCE.json)"),
    ("supplement", "its normal proportion bounded at 1/(n + 1)", "1", 1, "the bound's formula, 1/(n + 1)"),
    ("supplement", "its normal proportion bounded at 1/(n + 1)", "1", 2, "the bound's formula, 1/(n + 1)"),
    *[("supplement", "three donor groups (+3 to +6", p, k, POSITION)
      for p, k in (("3", 1), ("6", 1), ("7", 1), ("10", 1), ("11", 1), ("50", 1), ("3", 2), ("10", 2),
                   ("11", 2), ("20", 1), ("21", 1), ("50", 2))],
    ("supplement", "at donor +7 and beyond and acceptor −21 and beyond", "7", 1, BP7),
    ("supplement", "at donor +7 and beyond and acceptor −21 and beyond", "21", 1, BP7),
    ("supplement", "recomputed at priors of 0.05 and 0.185", "0.05", 1,
     "a prior of the sensitivity analysis, fixed in advance (src/evid_suppl_c4.py)"),
    ("supplement", "recomputed at priors of 0.05 and 0.185", "0.185", 1,
     "a prior of the sensitivity analysis, fixed in advance (src/evid_suppl_c4.py)"),
    ("supplement", "fewer than 10 variants at or above the fold’s threshold", "10", 1,
     "minimum band occupancy (src/evid_common.py)"),
]


# ------------------------------------------------------------------ the texts
def _paras(doc):
    import docx
    from docx.oxml.ns import qn
    for el in doc.element.body:
        if el.tag == qn("w:p"):
            p = docx.text.paragraph.Paragraph(el, doc)
            yield p.style.name, p.text.strip()


def write_texts(manuscript: Path, f3: Path, folder: Path) -> dict[str, Path]:
    """The two texts as the gate reads them (see the module's docstring)."""
    import docx
    d, f = docx.Document(str(manuscript)), docx.Document(str(f3))
    boxes = [t.cell(0, 0).text.strip() for t in f.tables if len(t.rows) == 1 and len(t.columns) == 1]
    abstracts = [x for x in boxes if x.startswith("Computational predictions of splicing")]
    if len(abstracts) != 1:
        raise SystemExit(f"{f3.name}: {len(abstracts)} boxes hold the abstract, not one")
    main, tables, supp, where = ["# Abstract", "", abstracts[0], ""], [], [], "title"
    for style, t in _paras(d):
        heading = style.startswith("Heading")
        if heading and where == "title" and t in ("Introduction", "Abstract"):
            where = "main"
        if heading and t == "References" and where == "main":
            where = "references"
            continue
        if heading and t == "Tables":
            where = "tables"
            continue
        if heading and t == "Supplementary Notes":
            where = "notes"
        if heading and t == "References cited in the Supplementary Notes":
            supp += ["# References cited in the Supplementary Notes", ""]
            where = "end"
            continue
        if not t or where in ("title", "references", "end"):
            continue
        level = int(style.split()[-1]) if heading and style.split()[-1].isdigit() else 0
        line = ("#" * level + " " + t) if level else t
        {"main": main, "tables": tables, "notes": supp}[where].extend([line, ""])
    main += ["# Tables", ""] + tables + ["# References", ""]
    folder.mkdir(parents=True, exist_ok=True)
    out = {"main": folder / "v6_main.md", "supplement": folder / "v6_supplement.md"}
    out["main"].write_text("\n".join(main) + "\n")
    out["supplement"].write_text("\n".join(supp) + "\n")
    return out


def _key(entries, doc: str, tok, rank: int):
    hits = [e for e in entries if e[0] == doc and e[1] in tok.sentence and e[2] == tok.printed and e[3] == rank]
    return hits


def rows_for(doc: str, text: str) -> tuple[list[dict], list]:
    """Each number of a text with its entry; the entries matched, for the check that every
    entry is used once."""
    out, used = [], []
    ranks: dict[tuple[str, str], int] = {}
    for t in gate.read_tokens(text):
        k = (t.sentence, t.printed)
        ranks[k] = ranks.get(k, 0) + 1
        bound, why = _key(BINDINGS, doc, t, ranks[k]), _key(REASONS, doc, t, ranks[k])
        if len(bound) + len(why) > 1:
            raise SystemExit(f"{doc} line {t.line}: {t.printed} has {len(bound) + len(why)} entries")
        clean = t.excerpt_space
        at = clean.find(t.context.clean[t.context.start:t.context.end].strip())
        excerpt = clean[max(0, t.context.start - 30):t.context.end + 30].strip()
        if excerpt not in clean:
            raise SystemExit(f"{doc} line {t.line}: the excerpt is not a piece of the sentence")
        row = {"document": doc, "line": t.line, "occurrence": t.occurrence, "printed": t.printed,
               "excerpt": excerpt}
        if bound:
            used.append(bound[0])
            row |= {"kind": "measured", "reason": "", **bound[0][4]}
        elif why:
            used.append(why[0])
            row |= {"kind": "constant", "reason": why[0][4], "file": "", "rows": "", "column": "",
                    "derive": "", "scale": ""}
        elif t.reason:
            row |= {"kind": "constant", "reason": t.reason, "file": "", "rows": "", "column": "",
                    "derive": "", "scale": ""}
        else:
            raise SystemExit(f"{doc} line {t.line}: {t.printed} is measured and has no binding "
                             f"({t.sentence[:120]})")
        out.append(row)
    return out, used


def keyed(text: str, rows: list[dict]) -> str:
    """The text with a source tag after each sentence that holds a measured number: the
    files its numbers are bound to, each named from the folder the gate searches. The
    sentences are cut as the gate cuts them (gate._pieces), each paragraph a line."""
    dirs = [Path(d).resolve() for d in gate.SEARCH_DIRS]
    files: dict[tuple[int, str], list[str]] = {}
    for t, r in zip(gate.read_tokens(text), rows):
        if r["kind"] != "measured":
            continue
        path = (REPO / r["file"]).resolve()
        name = min((str(path.relative_to(d)) for d in dirs if path.is_relative_to(d)), key=len)
        names = files.setdefault((t.line, t.sentence), [])
        if name not in names:
            names.append(name)
    lines = text.split("\n")
    for line in sorted({k[0] for k in files}):
        src, out, at = lines[line - 1], [], 0
        for off, sent in gate._pieces(src, gate.SENTENCE_BREAK):
            end = off + len(sent.rstrip())
            names = files.get((line, " ".join(sent.split())))
            out.append(src[at:end] + (" [" + "; ".join(names) + "]" if names else ""))
            at = end
        out.append(src[at:])
        lines[line - 1] = "".join(out)
    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manuscript", required=True, type=Path)
    ap.add_argument("--f3", required=True, type=Path)
    ap.add_argument("--texts", required=True, type=Path)
    args = ap.parse_args()
    with RunLog(__file__) as log:
        ms, f3 = args.manuscript.expanduser().resolve(), args.f3.expanduser().resolve()
        log.note("manuscript_sha256", hashlib.sha256(ms.read_bytes()).hexdigest())
        log.note("f3_sha256", hashlib.sha256(f3.read_bytes()).hexdigest())
        texts = write_texts(ms, f3, args.texts.expanduser().resolve())
        rows, used = [], []
        for doc in ("main", "supplement"):
            text = texts[doc].read_text()
            r, u = rows_for(doc, text)
            rows += r
            used += u
            (texts[doc].parent / f"v6_{doc}_keyed.md").write_text(keyed(text, r))
            log.note(f"{doc}_sha256", hashlib.sha256(texts[doc].read_bytes()).hexdigest())
        unused = [e[:4] for e in BINDINGS + REASONS if e not in used]
        twice = [e[:4] for e in used if used.count(e) > 1]
        if unused or twice:
            raise SystemExit(f"entries matched no number: {unused[:5]}; matched twice: {twice[:5]}")
        gate.write_bindings(OUT, rows)
        n = {k: sum(1 for r in rows if r["kind"] == k) for k in ("measured", "constant")}
        log.note("rows", n | {"total": len(rows)})
        print(f"{OUT.relative_to(REPO)}: {len(rows)} numbers, {n['measured']} measured, "
              f"{n['constant']} constants; texts in {texts['main'].parent}")


if __name__ == "__main__":
    main()
