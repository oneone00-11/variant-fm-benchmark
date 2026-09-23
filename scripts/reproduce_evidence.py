#!/usr/bin/env python3
"""
evidence-strength reframe -- one command for every analysis stage of the reworked study.

What this reproduces, and what it does not. The analysis stages below run from
tracked inputs and are deterministic under the seed in `phase1/src/config.py`.
Three inputs are NOT rebuilt here because they are third-party downloads or model
re-scores, and each is pinned by sha256 in a manifest instead:

  * the ClinVar GRCh38 VCF of 15 June 2026 (192 MB; `phase1/data/evidence/clinvar/`,
    gitignored, md5 checked against NCBI's own file on download);
  * the `spliceai_walker` column, re-scored at Walker's -D 4999 in the atlas's
    pinned SpliceAI 1.3.1 environment (`src/evid_score_spliceai_walker.py`);
  * the external-gene columns (`src/evid_score_spliceai_events.py` and the E7
    scoring), likewise;
  * the AlphaGenome Atlas columns (`src/evid_score_avi.py`), which need an API key
    outside the tree and a separate venv carrying alphagenome>=0.9.0 -- the pinned
    .venv keeps 0.7.0 because that is the provenance of the alphagenome column.

Each of those prints the command that produces it when its input is missing.

Stages
------
  1  evid_build_set          rebuild the 1 <= |offset| <= 50 analysis set (E1)
  2  evid_walker_thresholds  the ClinGen SVI fixed cut points, by stratum and arm (E2)
  3  evid_interval_lr        score-to-evidence intervals, Pejaver's procedure (E3)
  4  evid_territory_metrics  AUROC / PR-AUC by territory and ClinVar arm (E4)
  5  evid_inframe            what the false positives are predicting (E6), then the
                             same attribution on the two external genes (E2.7)
  6  evid_external TP53      fixed and fitted thresholds on the held-out gene (E7)
  7  evid_external DDX3X     merge the scored columns, then the same on a gene
                             outside the seven (E7)
  8  evid_diagnostics        cut-point and monotonicity diagnostic tables (E8)
  9  evid_tier_logo          in-sample tier vs held-out ratio, per fold (E2.2)
 10  evid_arms               ClinVar arms without BRCA1, and within gene (E2.3/E2.4)
 11  evid_fig_data           the figure tables and draft PNGs (E2.9)
 12  evid_training_provenance where each predictor's training signal comes from (E9)
 13  evid_fusion_stability   the fusion's coefficients in every fold (E10)
 14  evid_dilution           thresholds refitted on every subset of training genes (E11)
 15  evid_tables             the four main tables as printed (E12)
 16  evid_supp_tables        the supplementary tables as printed (E13)
 17  evid_figures            main and supplementary figures, PDF and PNG (E14)
 18  evid_delta              every published quantity with a counterpart (E8)

Outputs land in `phase1/reports/evidence/`.

    python scripts/reproduce_evidence.py
    python scripts/reproduce_evidence.py --from 3      # resume at a stage
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PHASE1 = REPO / "phase1"

STAGES = [
    ("src.evid_build_set",          "E1  rebuild the analysis set", []),
    ("src.evid_walker_thresholds",  "E2  Walker fixed cut points", []),
    ("src.evid_interval_lr",        "E3  score-to-evidence intervals", []),
    ("src.evid_territory_metrics",  "E4  territory and ClinVar-arm metrics", []),
    ("src.evid_inframe",            "E6  in-frame attribution", ["--attribute"]),
    ("src.evid_inframe",            "E2.7 in-frame attribution: DDX3X", ["--external-attribute", "ddx3x"]),
    ("src.evid_inframe",            "E2.7 in-frame attribution: TP53", ["--external-attribute", "tp53"]),
    ("src.evid_external",           "E7  external gene: TP53", ["--apply", "TP53"]),
    ("src.evid_external",           "E7  external gene: DDX3X, merge scored columns",
     ["--merge-scores", "ddx3x"]),
    ("src.evid_external",           "E7  external gene: DDX3X", ["--apply", "DDX3X"]),
    ("src.evid_diagnostics",        "E8  cut-point and monotonicity diagnostics", []),
    ("src.evid_tier_logo",          "E2.2 in-sample tier against held-out ratio", []),
    ("src.evid_arms",               "E2.3/E2.4 ClinVar arms without the gene confound", []),
    ("src.evid_fig_data",           "E2.9 one tidy table per figure, plus draft PNGs", []),
    ("src.evid_training_provenance", "E9  predictor training signals and overlaps", []),
    ("src.evid_fusion_stability",   "E10 fusion coefficients across folds", []),
    ("src.evid_dilution",           "E11 thresholds refitted on subsets of training genes", []),
    ("src.evid_tables",             "E12 the four main tables, as printed", []),
    ("src.evid_supp_tables",        "E13 the supplementary tables, as printed", []),
    ("src.evid_figures",            "E14 main and supplementary figures", []),
    ("src.evid_delta",              "E8  old/new quantity list -> docs/evidence-delta.md", []),
]

# Stages that need an input this script does not produce, and what produces it.
NEEDS = {
    "src.evid_tier_logo": [
        ("phase1/reports/evidence/evidence_thresholds_logo_folds.csv",
         "python -m src.evid_interval_lr (stage 3)"),
    ],
    "src.evid_build_set": [
        ("phase1/data/evidence/clinvar/clinvar_20260615.vcf.gz",
         "curl -O https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/"
         "clinvar_20260615.vcf.gz   (into phase1/data/evidence/clinvar/)"),
    ],
    "src.evid_inframe": [
        ("phase1/data/evidence/inframe_events.parquet",
         "python -m src.evid_inframe --subset, then the SpliceAI-environment command "
         "it prints"),
    ],
}


def run_stage(module: str, label: str, n: int, total: int, args: list) -> bool:
    """Run one stage. Returns False if it was skipped for a missing input.

    A missing input skips its own stage and lets the rest run. Exiting the whole
    script would mean an un-scored model column silently prevents the stages after
    it from regenerating at all, which is how outputs drift from the code that is
    supposed to produce them.
    """
    print(f"\n{'=' * 74}\n[{n}/{total}] {module}  --  {label}\n{'=' * 74}", flush=True)
    missing = [(rel, how) for rel, how in NEEDS.get(module, [])
               if not (REPO / rel).exists()]
    if missing:
        for rel, how in missing:
            print(f"  SKIPPED -- missing input {rel}\n    produce it with:  {how}")
        return False
    r = subprocess.run([sys.executable, "-m", module, *args], cwd=PHASE1)
    if r.returncode != 0:
        sys.exit(f"\nFAILED at stage {n}/{total} ({module}); exit code {r.returncode}")
    return True


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--from", dest="start", type=int, default=1)
    args = ap.parse_args()
    total = len(STAGES)
    skipped = []
    for i, (mod, label, extra) in enumerate(STAGES, 1):
        if i < args.start:
            print(f"[{i}/{total}] {mod} -- skipped (--from)")
            continue
        if not run_stage(mod, label, i, total, extra):
            skipped.append(f"{i}. {mod} -- {label}")
    print(f"\nTables: phase1/reports/evidence/")
    if skipped:
        print("\nStages skipped for a missing input -- their outputs are whatever "
              "the last successful run left:")
        for s in skipped:
            print(f"  {s}")
        sys.exit(1)
    print("All stages ran.\n")


if __name__ == "__main__":
    main()
