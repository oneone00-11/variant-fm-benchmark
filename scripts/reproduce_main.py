#!/usr/bin/env python3
"""One-click reproduction of the main results (Milestone 8 / 8b).

Starts from the final, NT-included evaluation matrix
(`data/processed/score_matrix_final.tsv`) and regenerates every result table and
all three main figures. No model re-scoring is performed and no score is modified
here — this only re-runs the Phase-3 evaluation on the already-computed scores.

Reproduces:
  results/tables/spearman_meta.tsv      <- PRIMARY metric (meta Spearman vs functional GS)
  results/tables/spearman_by_gene.tsv
  results/tables/auroc_clinvar.tsv      <- auxiliary (ClinVar binary, with/without BRCA1)
  results/tables/pairwise_splice.tsv
  results/tables/vus_reclass.tsv
  results/tables/sens_at_95spec.tsv
  results/figures/fig1_splice_performance.png
  results/figures/fig2_forest_spliceai_splice.png
  results/figures/fig3_region_heatmap.png

All randomness is seeded (np.random.default_rng(20260619)), so outputs are
deterministic across machines.

Usage:
    python scripts/reproduce_main.py
"""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
MATRIX = ROOT / "data" / "processed" / "score_matrix_final.tsv"

# Ordered Phase-3 pipeline. 103 (figures) depends on 100's tables.
STEPS = [
    ("100_spearman.py",          "Per-gene Spearman + DerSimonian-Laird meta  [PRIMARY]"),
    ("101_auroc.py",             "ClinVar binary AUROC (auxiliary, +/- BRCA1)"),
    ("102_clinical_pairwise.py", "Paired-bootstrap deltas, VUS reclass, sens@95%spec"),
    ("103_figures.py",           "Figures 1-3"),
]


def main():
    if not MATRIX.exists():
        sys.exit(
            f"ERROR: {MATRIX} not found.\n"
            "This file (the final evaluation matrix) ships in the repo under "
            "data/processed/. If missing, see README 'Data acquisition'."
        )

    print("=" * 72)
    print("Reproducing main results from:", MATRIX.relative_to(ROOT))
    print("=" * 72)
    for i, (script, desc) in enumerate(STEPS, 1):
        print(f"\n[{i}/{len(STEPS)}] {script}  —  {desc}")
        # Run each script standalone (it puts scripts/ on sys.path for `from config import`).
        subprocess.run([sys.executable, str(SCRIPTS / script)], check=True)

    _print_headline()


def _print_headline():
    """Print the headline splice ranking from the regenerated spearman_meta.tsv."""
    import pandas as pd
    meta = pd.read_csv(ROOT / "results" / "tables" / "spearman_meta.tsv", sep="\t")
    sp = meta[meta.region == "splice"].copy()
    sp["meta_rho"] = pd.to_numeric(sp["meta_rho"], errors="coerce")
    sp = sp.sort_values("meta_rho", ascending=False, na_position="last")
    print("\n" + "=" * 72)
    print("HEADLINE — splice subset, meta Spearman rho vs functional gold standard")
    print("=" * 72)
    for r in sp.itertuples():
        rho = "n=1 (no coverage)" if r.meta_rho != r.meta_rho else f"{r.meta_rho:.3f}"
        print(f"  {r.model:24s} {rho:>17s}   n={r.total_n}")
    print("\nExpected top tier (statistically tied): Pangolin ~0.76, SpliceAI ~0.75,")
    print("AlphaGenome ~0.75 > CADD ~0.69 ~ GPN-MSA ~0.67 > NT ~0.49 > gnomAD ~0.19.")
    print("AlphaMissense scores only 1/1781 splice variants — the gap DNA models fill.")
    print("\nAll tables -> results/tables/ ; all figures -> results/figures/")


if __name__ == "__main__":
    main()
