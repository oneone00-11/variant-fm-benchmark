"""Every tracked result must be reachable from a documented command.

The defect this exists for: `phase5_likelihood_ratios`, `phase5b_evidence_yield`,
`phase6_sensitivity`, `phase7_label_contrast`, `phase2b_no_offset_drop` and
`ipw_reweight` produced twenty-six of the forty-eight tracked report files, and
appeared in no reproduction script and nowhere in README, docs/ or scripts/.
Their outputs were committed, so every number in the manuscript could be checked
-- but a reader following the README could not regenerate them, and nothing in
the repository noticed. Writing the missing stages into the scripts fixes the
instance; this fixes the class.

The check is deliberately not a hand-maintained list of "modules that should be
covered": such a list is one more thing to forget, and forgetting is the failure
mode. It starts from the outputs that actually ship and works backwards.
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SRC = REPO / "phase1" / "src"
REPORTS = REPO / "phase1" / "reports" / "phase1"

# A stage entry looks like ("src.phase2_model", "H1 fusion / H2 ablation", []).
STAGE = re.compile(r'\(\s*"src\.([A-Za-z_0-9]+)"')


def _tracked(paths: str) -> list[str]:
    r = subprocess.run(["git", "ls-files", paths], cwd=REPO,
                       capture_output=True, text=True)
    if r.returncode != 0:
        return []
    return [ln.strip() for ln in r.stdout.splitlines() if ln.strip()]


def entry_points() -> dict[str, set[str]]:
    """script name -> the modules its stage table runs."""
    out = {}
    for script in sorted((REPO / "scripts").glob("reproduce_*.py")):
        text = script.read_text(errors="ignore")
        block = re.search(r"STAGES\s*=\s*\[(.*?)\n\]", text, re.S)
        if not block:
            continue
        out[script.name] = set(STAGE.findall(block.group(1)))
    return out


def producers(name: str) -> set[str]:
    """Modules under phase1/src that WRITE this output.

    Writers put the filename and `to_csv` in one statement:
        base.to_csv(C.REPORT_DIR / "phase7_base.csv", index=False)
    Merely naming the file is not producing it -- `phase8_robustness` imports
    from `phase5_likelihood_ratios`, and `phase7` prints the names it wrote.
    An earlier draft of this check matched any mention and therefore reported
    phase5_likelihood_ratios.csv as covered by phase8, which does not write it.

    Matched across the whole call rather than one line, because the writes wrap:

        frame_balance(frame, set(df["key"])).to_csv(
            C.REPORT_DIR / "ipw_frame_balance.csv", index=False)
    """
    out = set()
    for p in SRC.glob("*.py"):
        flat = re.sub(r"\s+", " ", p.read_text(errors="ignore"))
        for m in re.finditer(r"to_csv\(", flat):
            if name in flat[m.end():m.end() + 100]:
                out.add(p.stem)
                break
    return out


def test_every_reproduce_script_exposes_a_stage_table():
    eps = entry_points()
    assert eps, "no scripts/reproduce_*.py with a STAGES table was found"
    for name, mods in eps.items():
        assert mods, f"{name} has a STAGES block but no src.* modules in it"


def test_every_tracked_result_has_a_documented_entry_point():
    tracked = [p for p in _tracked("phase1/reports/phase1") if p.endswith(".csv")]
    if not tracked:
        pytest.skip("no tracked report files in this checkout")

    covered = set().union(*entry_points().values())
    orphans, unknown = [], []
    for rel in tracked:
        stem = Path(rel).name
        who = producers(stem)
        if not who:
            unknown.append(rel)
        elif not (who & covered):
            orphans.append(f"{Path(rel).name}  (written by {', '.join(sorted(who))})")

    assert not unknown, (
        "tracked results whose producing module could not be identified:\n  "
        + "\n  ".join(unknown))
    assert not orphans, (
        f"{len(orphans)} tracked result(s) are produced by modules that no "
        "reproduction script runs, so a reader following the README cannot "
        "regenerate them:\n  " + "\n  ".join(sorted(orphans)))


def test_the_readme_names_every_entry_point():
    """A stage table nobody is told to run is not a documented path."""
    readme = (REPO / "README.md").read_text(errors="ignore")
    missing = [s for s in entry_points() if s not in readme]
    assert not missing, (
        "reproduction scripts absent from README.md: " + ", ".join(sorted(missing)))
