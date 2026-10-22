"""E17.1 (read-only): are the four Atlas columns raw scores or quantiles?

avi, avi_splice_sites, avi_splice_site_usage and avi_splice_junctions came from the
AlphaGenome Atlas API with client 0.9.0 (phase1/data/evidence/avi.parquet.provenance.json).
phase1/src/evid_score_avi._aggregate takes, per variant, the maximum of the returned
AnnData's X. What X holds is the client's choice: this downloads the 0.9.0 wheel, reads
alphagenome/atlas/atlas.py, and checks the lines that fill X (the service's scores) and
the 'quantiles' layer (its calibrated scores). Nothing is changed; the answer and the
lines it rests on, with the sha256 of each file read, go to
phase1/data/evidence/alphagenome_official/ATLAS_COLUMNS_CHECK.json.

Run from the repository root:
    .venv-ag-official/bin/python scripts/stage_E17/atlas_columns_check.py
"""
from __future__ import annotations

import hashlib
import inspect
import json
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel, sha256_of  # noqa: E402

PHASE1 = REPO / "phase1"
sys.path.insert(0, str(PHASE1))
OUT = PHASE1 / "data" / "evidence" / "alphagenome_official" / "ATLAS_COLUMNS_CHECK.json"
PROV = PHASE1 / "data" / "evidence" / "avi.parquet.provenance.json"
CLIENT = "0.9.0"
# (what the line shows, the text it holds)
ATLAS_LINES = [
    ("the service's scores, as float32", "x = np.frombuffer(score.scores, dtype=np.float32)"),
    ("its calibrated scores", "score.calibrated_scores, dtype=np.float32"),
    ("the calibrated scores go to a layer named 'quantiles'",
     "layers = {'quantiles': np.concatenate(calibrated_scores)}"),
    ("X is the service's scores", "X=scores, obs=obs, var=var, layers=layers"),
]


def main() -> None:
    with RunLog(__file__) as log:
        prov = json.loads(PROV.read_text())
        log.input(PROV)
        if prov["client_version"] != CLIENT:
            raise SystemExit(f"the Atlas columns came from client {prov['client_version']}")
        with tempfile.TemporaryDirectory() as d:
            subprocess.run([sys.executable, "-m", "pip", "download", "--no-deps", "--quiet",
                            "--disable-pip-version-check", f"alphagenome=={CLIENT}",
                            "--dest", d], check=True)
            wheel = next(Path(d).glob("*.whl"))
            wheel_sha = hashlib.sha256(wheel.read_bytes()).hexdigest()
            with zipfile.ZipFile(wheel) as z:
                atlas = z.read("alphagenome/atlas/atlas.py").decode()
        lines = atlas.splitlines()
        found = []
        for what, text in ATLAS_LINES:
            hits = [i + 1 for i, line in enumerate(lines) if text in line]
            if not hits:
                raise SystemExit(f"atlas.py of {CLIENT} has no line holding {text!r}")
            found.append({"shows": what, "line": hits[0], "text": lines[hits[0] - 1].strip()})
        from src import evid_score_avi as A
        agg = inspect.getsource(A._aggregate)
        if "np.asarray(ad.X, dtype=float)" not in agg:
            raise SystemExit("evid_score_avi._aggregate does not read ad.X")
        record = {
            "written_by": rel(Path(__file__).resolve()),
            "columns": prov["columns"],
            "provenance": rel(PROV), "client_version": prov["client_version"],
            "scorers_requested": prov["scorers_requested"],
            "aggregation_recorded": prov["aggregation"],
            "our_code": {"file": "phase1/src/evid_score_avi.py", "function": "_aggregate",
                         "reads": "ad.X", "sha256": sha256_of(PHASE1 / "src" / "evid_score_avi.py")},
            "client_source": {"wheel": wheel.name, "wheel_sha256": wheel_sha,
                              "file": "alphagenome/atlas/atlas.py",
                              "sha256": hashlib.sha256(atlas.encode()).hexdigest(),
                              "lines": found},
            "answer": "raw: each column is the maximum, over tracks and rows, of X, which "
                      "client 0.9.0 fills with the service's scores; the calibrated scores "
                      "(quantiles) sit in the 'quantiles' layer, which the columns do not "
                      "read",
        }
        OUT.parent.mkdir(parents=True, exist_ok=True)
        OUT.write_text(json.dumps(record, indent=1, ensure_ascii=False) + "\n")
        log.output(OUT)
        log.note("answer", record["answer"])
        print(record["answer"])
        for f in found:
            print(f"  atlas.py:{f['line']}  {f['text']}")


if __name__ == "__main__":
    main()
