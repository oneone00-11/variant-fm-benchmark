"""
Deterministically regenerate the Phase-1 raw build input from the tracked
upstream scoring matrix, so a clean clone can rebuild frozen_matrix_v1 without
`variant_scores.parquet` itself being committed.

This reproduces EXACTLY the conversion originally used to create the frozen-v1
input: read the tab-separated `data/processed/score_matrix_final.tsv` as strings,
add the `variant_id` key, and coerce the same numeric columns. `load_raw()` reads
parquet or comma-CSV (not tab-TSV), hence the parquet materialisation.

Run from phase1/:  python -m src.make_raw
Then:              python -m src.phase1_build_frozen_matrix
must reproduce manifest sha256 = 2a0e249b44906f11ed91ce4746aba7389d8fff70d791a17710a3e26ec66d3199
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd
from . import config as C

_HERE = Path(__file__).resolve()
PHASE1_DIR = _HERE.parents[1]                 # .../phase1
REPO_ROOT  = _HERE.parents[2]                 # repo root
SRC_TSV    = REPO_ROOT / "data" / "processed" / "score_matrix_final.tsv"

# Columns coerced to numeric in the original conversion (order/list is load-bearing
# for hash reproducibility -- do NOT add/remove without re-verifying the sha256).
NUMERIC_COLS = [
    "pos", "functional_score", "cadd_phred", "cadd_raw", "alphamissense",
    "gnomad_af_global", "gnomad_af_popmax", "phylop100way", "phastcons100way",
    "spliceai_ds", "pangolin_score", "alphagenome_splice", "gpn_msa_score",
    "nucleotide_transformer", "intron_offset", "stars",
]


def make_raw() -> Path:
    if not SRC_TSV.exists():
        raise SystemExit(f"[make_raw] upstream TSV not found: {SRC_TSV}")
    df = pd.read_csv(SRC_TSV, sep="\t", dtype=str)
    df["variant_id"] = (df["chrom"].astype(str) + "-" + df["pos"].astype(str) + "-"
                        + df["ref"].astype(str) + "-" + df["alt"].astype(str))
    for c in NUMERIC_COLS:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    out = C.RAW_MATRIX_PATH
    if not out.is_absolute():
        out = PHASE1_DIR / out            # resolve relative to phase1/, CWD-independent
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out, index=False)
    print(f"[make_raw] wrote {out}  shape={df.shape}  (from {SRC_TSV.name})")
    return out


if __name__ == "__main__":
    make_raw()
