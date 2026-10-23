"""E18.0 rule 1: the freeze guard.

The thresholds of the seven training genes were frozen at the tag frozen-thresholds-v1.
This records, for every file they rest on or are read from, its sha256 as the tag has
it (`git show frozen-thresholds-v1:<path>`), in phase1/data/evidence/FROZEN_AT_TAG.json:

  * the analysis set and its manifest;
  * the AlphaGenome scores of the seven genes, and every other record of
    phase1/data/evidence/alphagenome_official/ as tagged (E18.1: the seven genes'
    score files do not change by a byte);
  * the thresholds: evidence_thresholds.csv, evidence_thresholds_logo_folds.csv and
    tier_logo_table.csv;
  * every table3* file under tables/, and Supplementary Tables S6 and S7;
  * every interval_lr_*.csv;
  * every output of the fusion model (its coefficients, Table S9 and its curves).

phase1/tests/test_evid_frozen_guard.py then checks that HEAD holds each of these files
byte for byte as recorded. The tag object's own name, commit and time are recorded
beside them.

Run from the repository root:
    .venv/bin/python scripts/stage_E18/freeze_guard.py
"""
from __future__ import annotations

import fnmatch
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel  # noqa: E402

TAG = "frozen-thresholds-v1"
OUT = REPO / "phase1" / "data" / "evidence" / "FROZEN_AT_TAG.json"
EV = "phase1/reports/evidence/"
# (what the files are, the glob patterns that pick them out of the tag's tree)
GROUPS = [
    ("the analysis set and its manifest",
     ["phase1/data/evidence/analysis_set_v1.parquet",
      "phase1/data/evidence/analysis_set_v1.manifest.json"]),
    ("the seven genes' AlphaGenome scores and their records",
     ["phase1/data/evidence/alphagenome_official/*"]),
    ("the thresholds", [EV + "evidence_thresholds.csv", EV + "evidence_thresholds_logo_folds.csv",
                        EV + "tier_logo_table.csv"]),
    ("Table 3", [EV + "tables/table3*"]),
    ("Supplementary Tables S6 and S7", [EV + "supplement/tableS6_*", EV + "supplement/tableS7_*"]),
    ("the interval likelihood-ratio curves", [EV + "interval_lr_*.csv"]),
    ("the fusion model's outputs", [EV + "fusion_fold_coefficients.csv",
                                    EV + "supplement/tableS9_fusion.csv",
                                    EV + "interval_lr_fusion_enet_*.csv"]),
]


def _git(*args: str) -> bytes:
    return subprocess.run(["git", "-C", str(REPO), *args], check=True,
                          capture_output=True).stdout


def main() -> None:
    with RunLog(__file__) as log:
        obj = _git("cat-file", "-p", TAG).decode()
        head = dict(line.split(" ", 1) for line in obj.split("\n\n")[0].splitlines())
        if head.get("type") != "commit" or head.get("tag") != TAG:
            raise SystemExit(f"{TAG} is not an annotated tag of a commit")
        stamp = int(head["tagger"].rsplit(" ", 2)[1])
        tagged = datetime.fromtimestamp(stamp, timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        tree = _git("ls-tree", "-r", "--name-only", TAG).decode().splitlines()
        files, groups = {}, {}
        for what, patterns in GROUPS:
            hits = sorted({p for p in tree for pat in patterns if fnmatch.fnmatchcase(p, pat)})
            if not hits:
                raise SystemExit(f"{TAG} holds no file for {what}: {patterns}")
            groups[what] = hits
            for p in hits:
                blob = _git("show", f"{TAG}:{p}")
                files[p] = {"sha256": hashlib.sha256(blob).hexdigest(), "bytes": len(blob)}
        record = {
            "written_by": rel(Path(__file__).resolve()),
            "tag": TAG,
            "tag_object": _git("rev-parse", TAG).decode().strip(),
            "commit": head["object"],
            "tagged_utc": tagged,
            "how": f"sha256 of `git show {TAG}:<path>`",
            "groups": groups,
            "files": dict(sorted(files.items())),
        }
        OUT.write_text(json.dumps(record, indent=1, ensure_ascii=False) + "\n")
        log.output(OUT)
        log.note("tag", {k: record[k] for k in ("tag", "tag_object", "commit", "tagged_utc")})
        log.note("files", len(files))
        print(f"{TAG} ({record['tag_object'][:7]} -> {record['commit'][:7]}, {tagged}): "
              f"{len(files)} files -> {rel(OUT)}")
        for what, hits in groups.items():
            print(f"  {len(hits):3d}  {what}")


if __name__ == "__main__":
    main()
