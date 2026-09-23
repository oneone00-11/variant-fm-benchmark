"""E9 -- where each predictor's training signal comes from, and what it overlaps.

The question is whether the evaluated predictors carry information that overlaps
the data they are scored against, by two routes in particular: evolutionary
conservation, and functional annotations shared with the training sets. The question cannot be answered with a correlation, so it is answered with
provenance: for every score column in the panel, what supervised the model, whether
any clinical variant classification entered that supervision, and whether the
training data include measurements of the same kind as this study's standard.

Three dependencies are kept apart, because only the first is usually discussed.

  clinical circularity  a predictor trained on clinical assertions, then scored
                        against clinical assertions. No panel column is trained
                        on clinical assertions, and ClinVar status is a
                        stratifying variable rather than a label. The DDX3X
                        deposit's own classification, the primary label there,
                        has a decision boundary trained on clinically classified
                        variants, which is why DDX3X takes no part in the ClinVar
                        comparison.
  assay overlap         a predictor trained on multiplexed assay measurements,
                        then scored against measurements of the same kind.
  shared conservation   a predictor whose training signal is cross-species
                        constraint, scored against variants whose functional
                        effect is itself correlated with constraint. This is not
                        leakage of the labels; it is a shared cause, and it cannot
                        be removed by holding genes out.

Rows for the eight panel columns and the Walker-basis SpliceAI re-score are taken
from the companion atlas's curated table rather than restated here, so the two
papers cannot disagree about what a model was trained on. That table holds a
documentation URL per row and is the single place either paper edits.

The AlphaGenome Atlas columns are not in it -- the Atlas was released after it was
built -- so their rows are defined below, with the same standard of evidence: what
the project's own documentation states. An earlier version recorded the combined
Atlas score (AVI) as undocumented and as an aggregation of AlphaGenome tracks; the
Atlas paper documents it, and it is neither. AVI is a supervised model: a neural
network over AlphaGenome features, AlphaMissense, three protein-termination
features and two conservation scores, trained to separate gnomAD v4.1 variants
above and below a filtering allele frequency of 0.1%, with four saturation genome
editing studies -- including the BRCA1, RAD51C and DDX3X assays used here -- in the
validation set that selected its checkpoint. Its row says so, and it is the one
column whose training overlaps this study's standard.

Run (PYTHONPATH=phase1):  python -m src.evid_training_provenance
"""
from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path

import pandas as pd

from .evid_common import ATLAS_REPO as _atlas_repo_default  # noqa: E402
ATLAS_REPO = _atlas_repo_default  # EVID_ATLAS_REPO overrides; see evid_common
REPORT_DIR = Path("reports/evidence")
OUT = REPORT_DIR / "predictor_training_provenance.csv"

# score column -> the predictor row in the atlas's curated table
FROM_ATLAS = {
    "spliceai": "SpliceAI",
    "spliceai_walker": "SpliceAI",
    "pangolin": "Pangolin",
    "alphagenome": "AlphaGenome",
    "cadd": "CADD",
    "phylop": "phyloP-100way",
    "phastcons": "phastCons-100way",
    "gpn_msa": "GPN-MSA",
    "nt": "NT-v2-500M",
}

# What the training signal IS, in the terms this paper's question needs. A column
# is assigned to exactly one class; the classes are the routes by which a score
# could share information with a splicing readout.
SIGNAL_CLASS = {
    "spliceai": "Splicing (annotated junctions)",
    "spliceai_walker": "Splicing (annotated junctions)",
    "pangolin": "Splicing (measured site usage)",
    "alphagenome": "Functional genomics tracks, including splicing",
    "avi": "Supervised on population allele frequency, with AlphaGenome, "
           "AlphaMissense and conservation inputs",
    "avi_splice_sites": "Functional genomics tracks, including splicing",
    "avi_splice_site_usage": "Functional genomics tracks, including splicing",
    "avi_splice_junctions": "Functional genomics tracks, including splicing",
    "cadd": "Proxy contrast of simulated against fixed derived variants, over "
            "annotations that include SpliceAI and MMSplice predictions",
    "phylop": "Cross-species sequence constraint",
    "phastcons": "Cross-species sequence constraint",
    "gpn_msa": "Cross-species sequence constraint (self-supervised)",
    "nt": "Genomic sequence (self-supervised)",
    "fusion_enet": "Within-gene rank of the seven assays' continuous scores, "
                   "refitted per fold",
}

_AG_TRAINING = ("Multi-task supervision on functional genomics tracks "
                "(expression, splicing, chromatin) from ENCODE/GTEx-class data; "
                "no variant labels")
_ATLAS_PAPER = ("AlphaGenome Atlas team, AlphaGenome Atlas: in silico mutagenesis "
                "of the entire human genome improves prioritization and "
                "interpretation of non-coding variants (Google DeepMind, 2026), "
                "Methods")
_NONE_DOC = "None documented"

# The Atlas columns, defined here to the same standard as the atlas table's rows.
ATLAS_COLUMNS = {
    "avi": (
        "Neural network over precomputed AlphaGenome predictions (ten modality "
        "groups, including splicing, each the per-variant maximum across tissues), "
        "AlphaMissense, three protein-termination features and two conservation "
        "scores (PhastCons 470-way, Cactus 241-way), trained to separate gnomAD v4.1 "
        "variants above and below a filtering allele frequency of 0.1%",
        "No",
        "Yes: its checkpoint selection and early stopping used a validation set of "
        "complex-trait variants and four genome editing assays, the BRCA1 (Findlay "
        "2018), RAD51C (Olvera-León 2024) and DDX3X (Radford 2023) saturation genome "
        "editing assays and an ATM prime-editing screen (Lee 2025); BRCA2, BARD1, "
        "PALB2, VHL and BAP1 were among its test sets",
        f"{_ATLAS_PAPER}; https://alphagenome.google/atlas"),
    "avi_splice_sites": (
        f"{_AG_TRAINING}. Precomputed SPLICE_SITES scorer, retrieved by lookup",
        "No", _NONE_DOC,
        "https://www.alphagenomedocs.com/"),
    "avi_splice_site_usage": (
        f"{_AG_TRAINING}. Precomputed SPLICE_SITE_USAGE scorer, retrieved by lookup",
        "No", _NONE_DOC,
        "https://www.alphagenomedocs.com/"),
    "avi_splice_junctions": (
        f"{_AG_TRAINING}. Precomputed SPLICE_JUNCTIONS scorer, retrieved by lookup",
        "No", _NONE_DOC,
        "https://www.alphagenomedocs.com/"),
    "fusion_enet": (
        "Elastic net over the eight panel columns, refitted within each "
        "leave-one-gene-out fold on the training genes' functional labels",
        "No",
        "By construction: the fitted weights come from this study's own "
        "saturation genome editing labels",
        "phase1/src/evid_common.py (logo_fusion)"),
}

# Whether the column's own training signal overlaps the readout it is scored
# against here -- a cell survival or depletion measurement of gene function.
READOUT_OVERLAP = {
    "Splicing (annotated junctions)":
        "Indirect: predicts a splice event, not gene function; the assay reads "
        "function, so agreement requires the event to change function",
    "Splicing (measured site usage)":
        "Indirect: predicts a splice event, not gene function; the assay reads "
        "function, so agreement requires the event to change function",
    "Functional genomics tracks, including splicing":
        "Indirect: predicts molecular effects including splicing, not gene "
        "function; no multiplexed assay measurement documented in training",
    "Proxy contrast of simulated against fixed derived variants, over "
    "annotations that include SpliceAI and MMSplice predictions":
        "Shared cause: the proxy contrast is itself shaped by selection, which "
        "also shapes which variants damage function; and a shared input, since "
        "SpliceAI scores are among its annotations",
    "Cross-species sequence constraint":
        "Shared cause: constraint and functional damage are both consequences of "
        "selection; holding genes out does not remove it",
    "Cross-species sequence constraint (self-supervised)":
        "Shared cause: constraint and functional damage are both consequences of "
        "selection; holding genes out does not remove it",
    "Supervised on population allele frequency, with AlphaGenome, AlphaMissense "
    "and conservation inputs":
        "Direct for BRCA1 and RAD51C here and for the external gene DDX3X, whose "
        "assays selected its checkpoint; shared cause otherwise, since allele "
        "frequency and constraint are both shaped by selection",
    "Genomic sequence (self-supervised)":
        "Shared cause, indirectly: sequence patterns learned across species "
        "reflect selection; no labels of any kind enter training",
    "Within-gene rank of the seven assays' continuous scores, refitted per fold":
        "Direct: fitted on the standard itself, which is why it is evaluated "
        "leave-one-gene-out and reported in the supplement only",
}


# Terms of use per column, as recorded in this repository's LICENSE-DATA (the Atlas
# columns and, by the same terms of service, the AlphaGenome splice score) and in
# the companion atlas's LICENSE-DATA audit (every other column). Three questions a
# laboratory needs answered: may it be used commercially, may it inform a clinical
# decision, and may its output train another model.
_AG_TOS = ("AlphaGenome Terms of Service: non-commercial, research only; outputs not "
           "for clinical decision-making and not to train other models")
TERMS = {
    "spliceai": ("Trained models CC BY-NC 4.0 (Illumina); commercial use needs a "
                 "licence", "No", "No"),
    "spliceai_walker": ("Trained models CC BY-NC 4.0 (Illumina); commercial use needs "
                        "a licence", "No", "No"),
    "pangolin": ("Software GPL-3.0; its output is not encumbered by that licence",
                 "No", "No"),
    "alphagenome": (_AG_TOS, "Yes", "Yes"),
    "avi": (_AG_TOS, "Yes", "Yes"),
    "avi_splice_sites": (_AG_TOS, "Yes", "Yes"),
    "avi_splice_site_usage": (_AG_TOS, "Yes", "Yes"),
    "avi_splice_junctions": (_AG_TOS, "Yes", "Yes"),
    "cadd": ("Non-commercial; commercial licence from the University of Washington",
             "No", "No"),
    "phylop": ("UCSC track data, free for commercial use", "No", "No"),
    "phastcons": ("UCSC track data, free for commercial use", "No", "No"),
    "gpn_msa": ("MIT", "No", "No"),
    "nt": ("CC BY-NC-SA 4.0", "No", "No"),
    "fusion_enet": ("This study's model, under the terms of its inputs; the AlphaGenome "
                    "splice score is excluded from its features for that reason",
                    "No", "No"),
}


def _atlas_training() -> dict:
    spec = importlib.util.spec_from_file_location(
        "atlas_predictor_resources", ATLAS_REPO / "src/atlas/predictor_resources.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("atlas_predictor_resources", mod)
    spec.loader.exec_module(mod)
    return mod.TRAINING


def build() -> pd.DataFrame:
    """One row per score column, in the order the panel is reported."""
    training = _atlas_training()
    order = list(SIGNAL_CLASS)
    rows = []
    for col in order:
        if col in FROM_ATLAS:
            name = FROM_ATLAS[col]
            data, clinical, mave, src = training[name]
            origin = f"companion atlas predictor_resources.TRAINING['{name}']"
        else:
            data, clinical, mave, src = ATLAS_COLUMNS[col]
            origin = "this module"
        signal = SIGNAL_CLASS[col]
        rows.append({
            "score_column": col,
            "training_signal_class": signal,
            "training_data": data,
            "contains_clinical_labels": clinical,
            "documented_mave_or_sge_in_training": mave,
            "overlap_with_the_functional_readout": READOUT_OVERLAP[signal],
            "training_source": src,
            "row_origin": origin,
            "terms_of_use": TERMS[col][0],
            "terms_bar_clinical_decision_making": TERMS[col][1],
            "terms_bar_training_other_models": TERMS[col][2],
        })
    return pd.DataFrame(rows)


def main() -> None:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    out = build()
    out.to_csv(OUT, index=False)
    n_clin = int((out["contains_clinical_labels"] != "No").sum())
    mave = out.loc[out["documented_mave_or_sge_in_training"] != "None documented",
                   "score_column"].tolist()
    print(f"[E9] wrote {OUT} ({len(out)} score columns)")
    print(f"     columns whose training carries clinical classifications: {n_clin}")
    print("     columns with any documented multiplexed-assay component: "
          f"{len(mave)}{' -- ' + ', '.join(mave) if mave else ''}")
    print(out[["score_column", "training_signal_class",
               "contains_clinical_labels"]].to_string(index=False))


if __name__ == "__main__":
    main()
