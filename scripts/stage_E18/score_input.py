"""E18.2: the variants every commit-2 scorer scores, in one run each.

phase1/data/evidence/external/e18_score_input.parquet: the 39 anchors of
e18_anchors.parquet (set "anchor"; scripts/stage_E18/anchors.py), then the external
genes' variants (set "external"), in E18.1's order and by E18.1's rule: DDX3X, every
row of ddx3x_splice.parquet; TP53, every row of tp53_splice_scored_12nt.parquet; CTCF,
SFPQ, RAD51D, XRCC2 and SBDS, the rows of <gene>_variants.parquet at
3 <= |offset| <= 50. Columns: set, gene, variant_id, chrom, pos, ref, alt. Each scorer
scores the anchors and the external variants it is asked for in the same run, and its
anchors must equal the analysis set's stored values bit for bit (E18.0 rule 4).

Run from the repository root:
    .venv/bin/python scripts/stage_E18/score_input.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel  # noqa: E402

EXT = REPO / "phase1/data/evidence/external"
ANCHORS = EXT / "e18_anchors.parquet"
OUT = EXT / "e18_score_input.parquet"
NEW_GENES = ["CTCF", "SFPQ", "RAD51D", "XRCC2", "SBDS"]
COLS = ["set", "gene", "variant_id", "chrom", "pos", "ref", "alt"]


def main() -> None:
    with RunLog(__file__) as log:
        a = pd.read_parquet(ANCHORS)
        log.input(ANCHORS)
        parts = [a.assign(set="anchor")[COLS]]
        sources = [("DDX3X", EXT / "ddx3x_splice.parquet", False),
                   ("TP53", EXT / "tp53_splice_scored_12nt.parquet", False)] + [
            (g, EXT / f"{g.lower()}_variants.parquet", True) for g in NEW_GENES]
        for g, path, window in sources:
            v = pd.read_parquet(path)
            log.input(path)
            if window:
                v = v[v["offset"].abs().between(3, 50)]
            if set(v["gene"]) != {g}:
                raise SystemExit(f"{rel(path)} names genes {set(v['gene'])}")
            parts.append(v.assign(set="external")[COLS])
        df = pd.concat(parts, ignore_index=True)
        df["chrom"], df["pos"] = df["chrom"].astype(str), df["pos"].astype(int)
        n = df[df.set == "external"].gene.value_counts().to_dict()
        if (len(a), n["DDX3X"], n["TP53"], sum(n[g] for g in NEW_GENES), n["SBDS"]) != (
                39, 1857, 288, 4522, 447) or df["variant_id"].duplicated().any():
            raise SystemExit(f"not E18's variants: {len(a)} anchors, {n}")
        df.to_parquet(OUT, index=False)
        log.output(OUT)
        log.note("counts", {"anchor": len(a)} | n)
        print(f"{len(df)} variants ({len(a)} anchors, {len(df) - len(a)} external) -> {rel(OUT)}")


if __name__ == "__main__":
    main()
