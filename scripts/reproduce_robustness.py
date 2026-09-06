#!/usr/bin/env python3
"""
Study A -- one-command reproduction of the reviewer-requested robustness
analyses (Phase 8). This is the robustness companion to
`scripts/reproduce_calibration.py`: same frozen matrix, same seed, same hash
verification, different stage.

Clean-clone reproduction is two commands:

    pip install -r requirements.txt
    python scripts/reproduce_robustness.py

Stages
------
  0  make_raw                    regenerate the Phase-1 build input (deterministic)
  1  phase1_build_frozen_matrix  rebuild + sha256-verify frozen-matrix-v1
  2  phase8_robustness           the six robustness analyses:
                                 (a) exact sign-flip tests (2^k assignments) +
                                     wild cluster bootstrap, delta-rho and
                                     delta-Brier (all 8 conditions)
                                 (b) leave-two-genes-out delta-rho (21 pairs)
                                 (c) Hartung-Knapp CIs for the Table-1 leaderboard
                                 (d) Murphy Brier decomposition (3 binnings x
                                     4 conditions, cluster-bootstrap CIs)
                                 (e) stratified selection test (+/-BRCA1,
                                     offset strata, assay-only record types)
                                 (f) LR+ at specificity 0.90/0.95/0.975/0.99

Outputs land in `phase1/reports/phase1/` as phase8_*.csv. Everything is seeded
(`RANDOM_SEED` in `phase1/src/config.py`); output is deterministic. Runtime is
dominated by the leave-one-gene-out elastic-net refits (a few minutes on CPU).
"""
from __future__ import annotations
import subprocess
import sys
from pathlib import Path

REPO   = Path(__file__).resolve().parents[1]
PHASE1 = REPO / "phase1"
sys.path.insert(0, str(Path(__file__).resolve().parent))
from reproduce_calibration import verify_frozen  # noqa: E402  (shared hash pin)

STAGES = [
    ("src.make_raw",                   "regenerate Phase-1 build input"),
    ("src.phase1_build_frozen_matrix", "rebuild + verify frozen-matrix-v1"),
    ("src.phase8_robustness",          "reviewer robustness analyses"),
]


def run_stage(module: str, label: str, n: int, total: int) -> None:
    print(f"\n{'='*74}\n[{n}/{total}] {module}  --  {label}\n{'='*74}", flush=True)
    r = subprocess.run([sys.executable, "-m", module], cwd=PHASE1)
    if r.returncode != 0:
        sys.exit(f"\nFAILED at stage {n}/{total} ({module}); exit code {r.returncode}")


def headline() -> None:
    """Print the robustness numbers the revision's response letter rests on."""
    import pandas as pd
    rep = PHASE1 / "reports" / "phase1"
    print(f"\n{'='*74}\nHEADLINE -- Study A robustness (Phase 8)\n{'='*74}")
    try:
        sf = pd.read_csv(rep / "phase8_sign_flip_exact.csv")
        r = sf[(sf.analysis == "delta_rho") & (sf.condition == "BRCA1_included")].iloc[0]
        print(f"\n(a) sign-flip delta-rho   exact two-sided p = {r.p_exact:.4f} "
              f"(2^{r.k} = {r.n_assignments} assignments; wild p = {r.p_wild_rademacher:.4f})")
        r = sf[(sf.analysis == "delta_brier") &
               (sf.condition == "isotonic/y_assay/BRCA1_included")].iloc[0]
        print(f"    sign-flip delta-Brier exact two-sided p = {r.p_exact:.4f} "
              f"(primary condition; wild p = {r.p_wild_rademacher:.4f})")
    except Exception as e:
        print(f"  (sign-flip table unavailable: {e})")
    try:
        l2 = pd.read_csv(rep / "phase8_leave_two_genes_out.csv")
        print(f"(b) leave-two-genes-out   min {l2.delta_rho.min():+.4f} / "
              f"median {l2.delta_rho.median():+.4f} / max {l2.delta_rho.max():+.4f} "
              f"over {len(l2)} pairs (all positive: {bool((l2.delta_rho > 0).all())})")
    except Exception as e:
        print(f"  (leave-two-genes-out table unavailable: {e})")
    try:
        hk = pd.read_csv(rep / "phase8_leaderboard_hk.csv").set_index("model")
        f = hk.loc["M1_enet"]
        print(f"(c) fusion rho {f.pooled_rho:.3f}  DL [{f.dl_lo:.3f},{f.dl_hi:.3f}]  "
              f"HK [{f.hk_lo:.3f},{f.hk_hi:.3f}]")
    except Exception as e:
        print(f"  (HK leaderboard unavailable: {e})")
    try:
        m = pd.read_csv(rep / "phase8_murphy_decomposition.csv")
        r = m[(m.condition == "isotonic/y_assay/BRCA1_included") & (m.binning == "width10")].iloc[0]
        print(f"(d) Murphy (primary)      dBrier {r.d_brier:+.4f} [{r.d_brier_lo:.4f},{r.d_brier_hi:.4f}]  "
              f"dREL {r.d_rel:+.4f} [{r.d_rel_lo:.4f},{r.d_rel_hi:.4f}]  "
              f"dRES {r.d_res:+.4f} [{r.d_res_lo:.4f},{r.d_res_hi:.4f}]")
    except Exception as e:
        print(f"  (Murphy table unavailable: {e})")
    print(f"\nFull result tables: phase1/reports/phase1/phase8_*.csv\n")


def main() -> None:
    if not PHASE1.is_dir():
        sys.exit(f"phase1/ not found under {REPO}")
    total = len(STAGES)
    for i, (mod, label) in enumerate(STAGES, 1):
        run_stage(mod, label, i, total)
        if mod == "src.phase1_build_frozen_matrix":
            verify_frozen()
    headline()


if __name__ == "__main__":
    main()
