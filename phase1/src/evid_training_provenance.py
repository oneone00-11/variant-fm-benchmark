"""E9 -- where each predictor's training signal comes from, and what it overlaps.

The question is whether the evaluated predictors carry information that overlaps
the data they are scored against, by two routes in particular: evolutionary
conservation, and functional annotations shared with the training sets. The question cannot be answered with a correlation, so it is answered with
provenance: for every score column in the panel, what supervised the model, whether
any clinical variant classification entered that supervision, and whether the
training data include measurements of the same kind as this study's standard.

Three dependencies are kept apart, because only the first is usually discussed.

  clinical circularity  a predictor trained on clinical assertions, then scored
                        against clinical assertions. This study cannot have it:
                        the standard is an assay measurement, and ClinVar status
                        is a stratifying variable rather than a label.
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
the project's own documentation states, and an explicit "not documented in the
sources checked" where it states nothing.

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
    "avi": "Functional genomics tracks, including splicing",
    "avi_splice_sites": "Functional genomics tracks, including splicing",
    "avi_splice_site_usage": "Functional genomics tracks, including splicing",
    "avi_splice_junctions": "Functional genomics tracks, including splicing",
    "cadd": "Proxy contrast of simulated against fixed derived variants",
    "phylop": "Cross-species sequence constraint",
    "phastcons": "Cross-species sequence constraint",
    "gpn_msa": "Cross-species sequence constraint (self-supervised)",
    "nt": "Genomic sequence (self-supervised)",
    "fusion_enet": "This study's own functional labels, refitted per fold",
}

_AG_TRAINING = ("Multi-task supervision on functional genomics tracks "
                "(expression, splicing, chromatin) from ENCODE/GTEx-class data; "
                "no variant labels")
_ATLAS_AGG = ("Precomputed AlphaGenome predictions aggregated by the Atlas; the "
              "combination behind the AVI score is not documented in the sources "
              "checked (product page and API documentation, accessed 2026-09-20)")
_NONE_DOC = "None documented"

# The Atlas columns, defined here to the same standard as the atlas table's rows.
ATLAS_COLUMNS = {
    "avi": (
        f"{_AG_TRAINING}. {_ATLAS_AGG}",
        "No", _NONE_DOC,
        "https://www.alphagenomedocs.com/ ; "
        "https://deepmind.google/science/alphagenome/"),
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
    "Proxy contrast of simulated against fixed derived variants":
        "Shared cause: the proxy contrast is itself shaped by selection, which "
        "also shapes which variants damage function",
    "Cross-species sequence constraint":
        "Shared cause: constraint and functional damage are both consequences of "
        "selection; holding genes out does not remove it",
    "Cross-species sequence constraint (self-supervised)":
        "Shared cause: constraint and functional damage are both consequences of "
        "selection; holding genes out does not remove it",
    "Genomic sequence (self-supervised)":
        "None documented: no labels of any kind enter training",
    "This study's own functional labels, refitted per fold":
        "Direct: fitted on the standard itself, which is why it is evaluated "
        "leave-one-gene-out and reported in the supplement only",
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
