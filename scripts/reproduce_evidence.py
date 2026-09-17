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
    scoring), likewise.

Each of those prints the command that produces it when its input is missing.

Stages
------
  1  evid_build_set          rebuild the 1 <= |offset| <= 50 analysis set (E1)
  2  evid_walker_thresholds  the ClinGen SVI fixed cut points, by stratum and arm (E2)
  3  evid_interval_lr        score-to-evidence intervals, Pejaver's procedure (E3)
  4  evid_territory_metrics  AUROC / PR-AUC by territory and ClinVar arm (E4)
  5  evid_inframe            what the false positives are predicting (E6)
  6  evid_external TP53      fixed and fitted thresholds on the held-out gene (E7)
  7  evid_external DDX3X     the same, on a gene outside the seven (E7)

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
    ("src.evid_external",           "E7  external gene: TP53", ["--apply", "TP53"]),
    ("src.evid_external",           "E7  external gene: DDX3X", ["--apply", "DDX3X"]),
]

# Stages that need an input this script does not produce, and what produces it.
NEEDS = {
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


def run_stage(module: str, label: str, n: int, total: int, args: list) -> None:
    print(f"\n{'=' * 74}\n[{n}/{total}] {module}  --  {label}\n{'=' * 74}", flush=True)
    for rel, how in NEEDS.get(module, []):
        if not (REPO / rel).exists():
            sys.exit(f"\nMISSING INPUT for stage {n}: {rel}\n  produce it with:  {how}")
    r = subprocess.run([sys.executable, "-m", module, *args], cwd=PHASE1)
    if r.returncode != 0:
        sys.exit(f"\nFAILED at stage {n}/{total} ({module}); exit code {r.returncode}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--from", dest="start", type=int, default=1)
    args = ap.parse_args()
    total = len(STAGES)
    for i, (mod, label, extra) in enumerate(STAGES, 1):
        if i < args.start:
            print(f"[{i}/{total}] {mod} -- skipped")
            continue
        run_stage(mod, label, i, total, extra)
    print(f"\nAll stages done. Tables: phase1/reports/evidence/\n")


if __name__ == "__main__":
    main()
