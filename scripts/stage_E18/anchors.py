"""E18.0 rule 4(a): the anchor variants every external scorer re-scores.

An external column is the analysis set's variable only if the scorer that produces it,
in the same run, reproduces the analysis set's stored values bit for bit on anchor
variants drawn from the analysis set. The anchors: for each of the seven training genes
and each of the strata s3_10 and s11_50, three variants (all of them where a cell holds
fewer), drawn without replacement with numpy's default generator seeded with
config.RANDOM_SEED, cell by cell in a fixed order (genes alphabetically, then s3_10,
s11_50), from the cell's variant_ids sorted. AlphaGenome's anchors are the pilot's ten
of the determinism check instead (E18.1).

The list, with the analysis set's stored value of every predictor column, goes to
phase1/data/evidence/external/e18_anchors.parquet; scorers read it from there.

Run from the repository root:
    .venv/bin/python scripts/stage_E18/anchors.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel  # noqa: E402

sys.path.insert(0, str(REPO / "phase1"))
from src import config as C  # noqa: E402

SET = REPO / "phase1/data/evidence/analysis_set_v1.parquet"
OUT = REPO / "phase1/data/evidence/external/e18_anchors.parquet"
STRATA = ["s3_10", "s11_50"]
PER_CELL = 3
# the analysis set's predictor columns, the thirteen E18 makes for the external genes
# (alphagenome's anchors are the drift checks; its stored value is kept for the record)
COLUMNS = ["spliceai", "spliceai_walker", "pangolin", "alphagenome", "nt", "cadd",
           "phylop", "phastcons", "gpn_msa", "avi", "avi_splice_sites",
           "avi_splice_site_usage", "avi_splice_junctions"]


def draw(df: pd.DataFrame) -> pd.DataFrame:
    rng = np.random.default_rng(C.RANDOM_SEED)
    picked = []
    for gene in sorted(C.ASSAY_LABELS):
        for st in STRATA:
            ids = np.array(sorted(df.loc[(df.gene == gene) & (df.stratum == st), "variant_id"]))
            if len(ids):
                picked += list(rng.choice(ids, size=min(PER_CELL, len(ids)), replace=False))
    return df.set_index("variant_id").loc[picked].reset_index()


def main() -> None:
    with RunLog(__file__) as log:
        df = pd.read_parquet(SET)
        log.input(SET)
        a = draw(df)
        keep = ["variant_id", "gene", "stratum", "chrom", "pos", "ref", "alt"] + COLUMNS
        a = a[keep]
        a.to_parquet(OUT, index=False)
        log.output(OUT)
        n = a.groupby(["gene", "stratum"]).size().to_dict()
        log.note("cells", {f"{g}/{s}": int(k) for (g, s), k in n.items()})
        log.note("n", len(a))
        print(f"{len(a)} anchors -> {rel(OUT)}")
        print(a.groupby(["gene", "stratum"]).size().unstack(fill_value=0).to_string())


if __name__ == "__main__":
    main()
