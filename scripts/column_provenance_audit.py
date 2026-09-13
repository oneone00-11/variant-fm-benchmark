"""Column-provenance audit: frozen-matrix-v1 against the companion atlas (read-only).

For every predictor column the two repositories share, on the variants they share:
Spearman rho between the two columns, whether the values are identical, the distinct
value counts, and the arithmetic ceiling Spearman rho can reach given the ties in each
column (the atlas's own `max_spearman_given_ties`, reused rather than re-implemented).
For SpliceAI and Pangolin the 1,781-variant splice subset (|offset| <= 8, Methods 2.1)
is audited separately, rounded and full-precision, overall and per gene.

    python scripts/column_provenance_audit.py --atlas-repo ../functional-standard-atlas \
        --out docs/column-provenance.md

Writes nothing except the report.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parents[1]
FROZEN = REPO / "phase1" / "data" / "frozen" / "frozen_matrix_v1.parquet"

# canonical frozen-matrix column -> atlas column
SHARED = {
    "spliceai": "spliceai_ds",
    "pangolin": "pangolin_score",
    "alphagenome": "alphagenome",          # different definition (v0.6.1 vs v0.7.0)
    "gpn_msa": "gpn_msa",                  # atlas stores the sign-flipped score
    "nt": "nucleotide_transformer",
    "cadd": "cadd",
    "alphamissense": "alphamissense",
    "phylop": "phylop100way",
    "phastcons": "phastcons100way",
    "gnomad_af": "gnomad_af_global",       # atlas stores the sign-flipped AF
}

# Provenance as recorded in each repository (Supplementary Table S1 and scripts/ here;
# models/*/NOTES.md, results/predictor_resources_v1.tsv and Methods there).
PROVENANCE = {
    "spliceai":      ("SpliceAI 1.3.1 CLI (TF 2.19.1 CPU, WSL), -A grch38 -D 50, max of 4 deltas as "
                      "printed by the CLI (2 d.p.); scripts/80_make_vcf.py -> 82_assemble_v2.py",
                      "SpliceAI 1.3.1, same annotation/distance; models/spliceai/score_fullprec.py "
                      "patches the CLI's {:.2f} to {:.17g} (--validate: re-rounded output reproduces "
                      "the stock CLI bit for bit); spliceai_ds = full precision, spliceai_ds_cli_rounded kept"),
    "pangolin":      ("Pangolin, git tkzeng/Pangolin main (~2026-06, unpinned; later shown to equal "
                      "commit 5cf94b8), torch 2.12.0 CPU, -d 50 CSV mode, max(gain,|loss|) as printed (2 d.p.)",
                      "Pangolin pinned to 5cf94b8db938c658391b4305cd7ce33297d44ff7, torch 2.13.0 CPU, -d 50; "
                      "models/pangolin/score_fullprec.py removes the two round(...,2) format calls; "
                      "pangolin_score = full precision, pangolin_score_cli_rounded kept"),
    "nt":            ("InstaDeepAI/nucleotide-transformer-v2-500m-multi-species, masked 6-mer LLR, one-off "
                      "cloud-GPU run whose script was not committed; window unrecorded; scripts/92_score_nt.py "
                      "is an archival reconstruction (WINDOW_BP = 6000 assumed)",
                      "same checkpoint, models/nt/score.py: 6,000 bp window, variant aligned inside one 6-mer, "
                      "RunPod A6000, torch 2.11.0+cu128, transformers 4.46.3 (pipfreeze_runpod.txt); "
                      "validated against this repo's column at rho = 0.9997"),
    "alphagenome":   ("AlphaGenome API client 0.6.1, 1-Mb interval, recommended SPLICE_* scorers, "
                      "aggregate = max |raw score| (scripts/91_score_alphagenome.py)",
                      "API client 0.7.0, merged-quantile splice score over a 16-kb window "
                      "(models/alphagenome/score.py; definition sweep in sweep_definitions.py); the atlas "
                      "also stores this repo's column as data/external/companion_alphagenome_v061.parquet"),
    "gpn_msa":       ("songlab/gpn-msa-hg38-scores, revision cf1718a9, tabix lookup (scripts/90_score_gpn.py); "
                      "raw sign (negative = more deleterious), reoriented in phase 2",
                      "same precomputed scores, sign-flipped at scoring time (models/gpn_msa/score.py)"),
    "cadd":          ("CADD GRCh38-v1.7 PHRED via REST API (scripts/75_cadd.py)",
                      "CADD GRCh38-v1.7 precomputed scores via remote tabix (models/cadd/score.py)"),
    "alphamissense": ("AlphaMissense_hg38.tsv.gz (Zenodo 8208688) lookup (scripts/74_alphamissense.py)",
                      "same file (models/alphamissense/score.py)"),
    "phylop":        ("UCSC phyloP100way via track API (scripts/73_conservation.py)",
                      "same track API (models/conservation/score.py)"),
    "phastcons":     ("UCSC phastCons100way via track API (scripts/73_conservation.py)",
                      "same track API (models/conservation/score.py)"),
    "gnomad_af":     ("gnomAD v4 global AF via GraphQL (scripts/72_gnomad.py); raw AF, reoriented in phase 2",
                      "gnomAD v4.1 GraphQL, sign-flipped at scoring time (models/gnomad/score.py)"),
}


def load_atlas(atlas: Path) -> pd.DataFrame:
    m = pd.read_parquet(atlas / "results" / "score_matrix_atlas_v2.parquet")
    m["variant_id_cpra"] = (m["chrom"].astype(str) + "-" + m["pos"].astype(int).astype(str)
                            + "-" + m["ref"] + "-" + m["alt"])
    return m.set_index("variant_id_cpra")


def tie_stats(v: np.ndarray, ceiling) -> dict:
    v = v[~np.isnan(v)]
    vals, counts = np.unique(v, return_counts=True)
    return {"n": len(v), "distinct": int(len(vals)), "max_tie": int(counts.max()) if len(counts) else 0,
            "rho_ceiling": ceiling(v) if len(v) >= 3 else float("nan")}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--atlas-repo", required=True)
    ap.add_argument("--out", default="docs/column-provenance.md")
    a = ap.parse_args()
    atlas = Path(a.atlas_repo).resolve()
    sys.path.insert(0, str(atlas / "src"))
    from atlas.robustness import max_spearman_given_ties as ceiling  # reuse, do not reimplement

    fz = pd.read_parquet(FROZEN)
    at = load_atlas(atlas)
    shared_ids = fz["variant_id"][fz["variant_id"].isin(at.index)]
    only = fz[~fz["variant_id"].isin(at.index)]
    fs = fz.set_index("variant_id").loc[shared_ids]
    asub = at.loc[shared_ids]

    rows = []
    for col, acol in SHARED.items():
        x = fs[col].astype(float).to_numpy()
        y = asub[acol].astype(float).to_numpy()
        ok = ~(np.isnan(x) | np.isnan(y))
        rho = stats.spearmanr(x[ok], y[ok]).statistic if ok.sum() > 2 else float("nan")
        eq = float(np.mean(np.isclose(x[ok], y[ok], rtol=0, atol=1e-12)))
        eq_neg = float(np.mean(np.isclose(x[ok], -y[ok], rtol=0, atol=1e-12)))
        row = {
            "column": col, "atlas_column": acol, "n_shared_scored": int(ok.sum()),
            "n_frozen_scored": int((~np.isnan(x)).sum()), "n_atlas_scored": int((~np.isnan(y)).sum()),
            "spearman": rho, "identical_frac": eq, "identical_up_to_sign_frac": eq_neg,
            "distinct_frozen": int(pd.Series(x[ok]).nunique()), "distinct_atlas": int(pd.Series(y[ok]).nunique()),
            "ceiling_frozen": ceiling(x[ok]), "ceiling_atlas": ceiling(y[ok]),
        }
        if col in ("spliceai", "pangolin"):
            r = asub[acol + "_cli_rounded"].astype(float).to_numpy()
            okr = ok & ~np.isnan(r)
            row["identical_to_atlas_cli_rounded_frac"] = float(np.mean(np.isclose(x[okr], r[okr], atol=1e-12)))
            row["rounded_fullprec_reround_frac"] = float(np.mean(np.isclose(np.round(y[okr], 2), r[okr], atol=1e-9)))
        rows.append(row)
    tab = pd.DataFrame(rows)

    # splice subset audit, rounded (frozen) vs full precision (atlas)
    sp = fz[fz["is_splice"]].copy()
    sp_shared = sp[sp["variant_id"].isin(at.index)]
    sp_rows, gene_rows = [], []
    for col, acol in (("spliceai", "spliceai_ds"), ("pangolin", "pangolin_score")):
        r_all = sp[col].astype(float).to_numpy()
        f_all = at.loc[sp_shared["variant_id"], acol].astype(float).to_numpy()
        sp_rows.append({"column": col, "version": "rounded (frozen-matrix-v1, n=%d)" % len(sp), **tie_stats(r_all, ceiling)})
        sp_rows.append({"column": col, "version": "full precision (atlas, n=%d shared)" % len(sp_shared), **tie_stats(f_all, ceiling)})
        for g, gdf in sp_shared.groupby("gene", observed=True):
            rr = gdf[col].astype(float).to_numpy()
            ff = at.loc[gdf["variant_id"], acol].astype(float).to_numpy()
            gene_rows.append({"column": col, "gene": g, "n": len(gdf),
                              "distinct_rounded": tie_stats(rr, ceiling)["distinct"],
                              "ceiling_rounded": tie_stats(rr, ceiling)["rho_ceiling"],
                              "distinct_fullprec": tie_stats(ff, ceiling)["distinct"],
                              "ceiling_fullprec": tie_stats(ff, ceiling)["rho_ceiling"],
                              "rho_rounded_vs_fullprec": stats.spearmanr(rr, ff).statistic})
    sp_tab, gene_tab = pd.DataFrame(sp_rows), pd.DataFrame(gene_rows)

    def md(df: pd.DataFrame, floatfmt="{:.4f}") -> str:
        cols = list(df.columns)
        out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
        for _, r in df.iterrows():
            cells = []
            for c in cols:
                v = r[c]
                cells.append(floatfmt.format(v) if isinstance(v, float) and not pd.isna(v) else str(v))
            out.append("| " + " | ".join(cells) + " |")
        return "\n".join(out)

    prov = ["| column | this repository (frozen-matrix-v1) | companion atlas | full precision here / there |",
            "|---|---|---|---|"]
    for col in SHARED:
        here, there = PROVENANCE[col]
        fp = {"spliceai": "no / yes", "pangolin": "no / yes"}.get(col, "n/a (no output rounding involved)")
        prov.append(f"| {col} | {here} | {there} | {fp} |")

    only_txt = ", ".join(only["variant_id"]) if len(only) else "none"
    note = (
        "**Reading the table.** Eight columns are the same measurements in both repositories: CADD, "
        "AlphaMissense, phyloP, phastCons, gnomAD AF and GPN-MSA are identical value for value (the last two "
        "up to the atlas's sign flip), and the frozen SpliceAI and Pangolin columns are identical to the "
        "atlas's *rounded* CLI columns, whose full-precision counterparts re-round to them exactly. "
        "The two therefore differ only in output precision, not in model, version, reference or distance. "
        "Nucleotide Transformer is the same checkpoint and score definition re-run under a pinned "
        "environment; it is not value-identical (rho as tabulated) because the original run's context "
        "window was never recorded. AlphaGenome is a genuinely different score definition (client 0.6.1 "
        "max-|raw| versus 0.7.0 merged-quantile, 16-kb window) and is out of scope for a precision swap. "
        f"Sixteen frozen variants ({only_txt}) are absent from the atlas's RAD51C deposit and have no "
        "atlas column at all."
    )

    text = "\n".join([
        "# Column provenance: frozen-matrix-v1 versus the companion atlas", "",
        f"Generated by `scripts/column_provenance_audit.py` (read-only). Frozen matrix: {len(fz):,} SNVs; "
        f"atlas SNV rows matched by chrom-pos-ref-alt: {len(shared_ids):,}; frozen-only: {len(only)}.", "",
        "## 1. How each column was produced", "", *prov, "",
        "## 2. Agreement on the shared variants", "",
        "`spearman` is between the two columns; `identical_frac` is the share of shared scored variants whose "
        "values agree to 1e-12 (`_up_to_sign` allows the atlas's sign flip); `ceiling_*` is the largest Spearman "
        "rho any target could reach against that column given its ties (atlas `max_spearman_given_ties`). "
        "For SpliceAI and Pangolin, `identical_to_atlas_cli_rounded_frac` compares the frozen column with the "
        "atlas's rounded CLI column and `rounded_fullprec_reround_frac` checks that the atlas's full-precision "
        "column re-rounded to 2 d.p. reproduces it.", "",
        md(tab), "",
        "## 3. The 1,781-variant splice subset: tie structure, rounded versus full precision", "",
        md(sp_tab), "", "Per gene (the unit the per-gene Spearman is computed in):", "", md(gene_tab), "",
        "## 4. Summary", "", note, "",
    ])
    out = REPO / a.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text, encoding="utf8")
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
