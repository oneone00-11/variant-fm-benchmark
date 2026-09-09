#!/usr/bin/env python3
"""
Study A -- one-command reproduction of every analysis the main manuscript
reports: elastic-net fusion, out-of-gene isotonic/Platt calibration, ECE/Brier/
yield under leave-one-gene-out, external validation on the held-out gene TP53,
the likelihood ratios and their Tavtigian bands, the evidence yield, the
weight-stability and training-gene-thinning checks, the ClinVar/assay label
contrast and its selection test, the intronic-offset drop sensitivity, and the
sampling-frame reweighting.

The reviewer-requested Phase-8 robustness analyses have their own entry point,
`scripts/reproduce_robustness.py`. Between them the two scripts regenerate every
tracked file in `phase1/reports/phase1/`, which
tests/test_reproduction_coverage.py enforces.

Clean-clone reproduction is two commands:

    pip install -r requirements.txt
    python scripts/reproduce_calibration.py

This script re-scores NO model. It starts from the tracked upstream matrix
`data/processed/score_matrix_final.tsv`, regenerates the Phase-1 build input,
rebuilds the frozen matrix, and VERIFIES its sha256 against the pinned manifest
`phase1/data/frozen/manifest_v1.json` (git tag `frozen-matrix-v1`) before any
analysis runs. If the hash does not match, it stops -- a silent input change would
invalidate every downstream number.

Stages
------
   1  make_raw                   regenerate the Phase-1 build input (deterministic)
   2  phase1_build_frozen_matrix rebuild + sha256-verify frozen-matrix-v1
   3  phase1_directionality_check control-anchored orientation gate (all genes)
   4  phase2_model               H1 fusion vs best single, H2 evolution-axis ablation
   5  phase3_calibration         H3 calibration (ECE / Brier / clinical yield)
   6  phase4_external_tp53       external validation on TP53 (fully held out)
   7  phase5_likelihood_ratios   H4 LR+/LR- at 95% specificity, Tavtigian bands
   8  phase5b_evidence_yield     fraction reaching Moderate evidence or above
   9  phase6_sensitivity         cross-fold weight stability, training-gene thinning
  10  phase7_label_contrast      ClinVar/assay label divergence, selection test
  11  phase2b_no_offset_drop     intronic-offset drop sensitivity
  12  ipw_reweight               sampling-frame reweighting
  13  ipw_reweight balance       sampling-frame balance diagnostics

Outputs land in `phase1/reports/phase1/`. Everything is seeded
(`RANDOM_SEED` in `phase1/src/config.py`); output is deterministic.

Stages 7 to 13 carry the gene-clustered bootstraps and the leave-one-gene-out
refits, and dominate the runtime; stages 1 to 6 finish in about two minutes, so
the Delta-rho / Delta-Brier headline appears long before the run ends.

Note on reproducibility scope: the analysis is fully reproducible on CPU, but the
upstream Nucleotide Transformer and Pangolin *scores* are not (see the
"Reproducibility scope" table in README.md). Both are excluded from stage 5.
"""
from __future__ import annotations
import hashlib
import subprocess
import sys
from pathlib import Path

REPO   = Path(__file__).resolve().parents[1]
PHASE1 = REPO / "phase1"

# Every stage the main manuscript reports, in order. Stages 7 to 13 were absent
# from this list while the sections that rest on them -- the likelihood ratios,
# the evidence yield, the weight-stability and training-gene sensitivity checks,
# the label contrast and the sampling-frame reweighting -- were already in the
# paper. Their modules were committed and their outputs are tracked, so the
# numbers were never in doubt; what did not exist was a documented command that
# regenerates them. tests/test_reproduction_coverage.py now fails if any tracked
# output has no entry point, so this list cannot fall behind the pipeline again.
STAGES = [
    ("src.make_raw",                   "regenerate Phase-1 build input", []),
    ("src.phase1_build_frozen_matrix", "rebuild + verify frozen-matrix-v1", []),
    ("src.phase1_directionality_check","directionality gate", []),
    ("src.phase2_model",               "H1 fusion / H2 ablation", []),
    ("src.phase3_calibration",         "H3 calibration", []),
    ("src.phase4_external_tp53",       "TP53 external validation", []),
    ("src.phase5_likelihood_ratios",   "H4 likelihood ratios / Tavtigian bands", []),
    ("src.phase5b_evidence_yield",     "evidence yield at Moderate and above", []),
    ("src.phase6_sensitivity",         "weight stability / training-gene thinning", []),
    ("src.phase7_label_contrast",      "label divergence and the selection test", []),
    ("src.phase2b_no_offset_drop",     "intronic-offset drop sensitivity", []),
    ("src.ipw_reweight",               "sampling-frame reweighting", []),
    ("src.ipw_reweight",               "sampling-frame balance diagnostics", ["balance"]),
]


def run_stage(module: str, label: str, n: int, total: int,
              args: list | None = None) -> None:
    print(f"\n{'='*74}\n[{n}/{total}] {module}  --  {label}\n{'='*74}", flush=True)
    r = subprocess.run([sys.executable, "-m", module, *(args or [])], cwd=PHASE1)
    if r.returncode != 0:
        sys.exit(f"\nFAILED at stage {n}/{total} ({module}); exit code {r.returncode}")


def verify_frozen() -> None:
    """Fail loudly if the rebuilt matrix is not the one the manuscript used.

    The pin below is the published content hash of frozen-matrix-v1 (git tag
    `frozen-matrix-v1`, and `phase1/data/frozen/manifest_v1.json` as committed).
    It is hardcoded on purpose: the build step REWRITES the manifest, so checking
    the rebuild against the freshly written manifest would be circular.

    The hash is taken over a canonical CSV serialisation (rows sorted by
    `variant_id`), matching `phase1/src/phase1_build_frozen_matrix.freeze()` --
    parquet itself is not byte-reproducible across writer versions.
    """
    import pandas as pd

    EXPECTED_SHA256 = "2a0e249b44906f11ed91ce4746aba7389d8fff70d791a17710a3e26ec66d3199"
    EXPECTED_ROWS, EXPECTED_GENES, EXPECTED_SPLICE = 21410, 7, 1781

    built = PHASE1 / "data" / "frozen" / "frozen_matrix_v1.parquet"
    df = pd.read_parquet(built).sort_values("variant_id").reset_index(drop=True)
    got = hashlib.sha256(df.to_csv(index=False).encode()).hexdigest()

    print("\n--- frozen-matrix-v1 integrity check ---")
    print(f"  expected sha256 : {EXPECTED_SHA256}")
    print(f"  rebuilt  sha256 : {got}")
    if got != EXPECTED_SHA256:
        sys.exit(
            "\nFROZEN MATRIX MISMATCH -- the rebuilt matrix differs from the one the\n"
            "manuscript used. Refusing to continue: downstream numbers would not be\n"
            "comparable to the published results. Check out tag `frozen-matrix-v1`\n"
            "and re-run, or investigate the upstream matrix."
        )
    n_rows, n_genes = len(df), df["gene"].nunique()
    n_splice = int(df["is_splice"].sum())
    assert (n_rows, n_genes, n_splice) == (EXPECTED_ROWS, EXPECTED_GENES, EXPECTED_SPLICE), (
        f"shape drift: {n_rows}/{n_genes}/{n_splice} != "
        f"{EXPECTED_ROWS}/{EXPECTED_GENES}/{EXPECTED_SPLICE}")
    print(f"  OK -- {n_rows:,} rows, {n_genes} genes, {n_splice:,} splice variants\n")


def headline() -> None:
    """Print the two numbers the manuscript's claim rests on."""
    import pandas as pd
    rep = PHASE1 / "reports" / "phase1"
    print(f"\n{'='*74}\nHEADLINE -- Study A\n{'='*74}")
    try:
        h1 = pd.read_csv(rep / "phase2_H1_stratified.csv")
        row = h1[h1.stratum == "overall"].iloc[0]
        print(f"\nH1  ranking      delta-rho = {row.delta:+.4f}  {row.ci95}   [{row.H1}]")
    except Exception as e:
        print(f"  (H1 table unavailable: {e})")
    try:
        h3 = pd.read_csv(rep / "phase3_H3_headline.csv")
        r = h3[(h3.set == "y_assay/BRCA1_included") & (h3.calib == "isotonic")].iloc[0]
        print(f"H3  calibration  delta-Brier = {r.dBrier:+.4f}  {r.dBrier_ci_geneclust}")
    except Exception as e:
        print(f"  (H3 table unavailable: {e})")
    try:
        tp = pd.read_csv(rep / "phase4_tp53_external.csv").set_index("model")
        f = tp.loc["fusion_8feat", "Brier"]
        s = [i for i in tp.index if i.startswith("best_single")][0]
        print(f"TP53 external    Brier fusion = {f:.4f}  vs  {s} = {tp.loc[s,'Brier']:.4f}")
    except Exception as e:
        print(f"  (TP53 table unavailable: {e})")
    print(f"\nFull result tables: phase1/reports/phase1/\n")


def main() -> None:
    if not PHASE1.is_dir():
        sys.exit(f"phase1/ not found under {REPO}")
    total = len(STAGES)
    for i, (mod, label, args) in enumerate(STAGES, 1):
        run_stage(mod, label, i, total, args)
        if mod == "src.phase1_build_frozen_matrix":
            verify_frozen()
    headline()


if __name__ == "__main__":
    main()
