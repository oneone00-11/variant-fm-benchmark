"""E20b.1 (commit 1): the outputs downstream of E6's re-taken subsets, against HEAD.

  * inframe_attribution.csv: every row of HEAD's (tool + variant_id) unchanged in every
    field; the new rows exactly the variants E20b added (write_e6_subsets.py's record);
  * inframe_attribution_summary.csv: only the 3-10 nt cells change, each by the frames of
    the added variants;
  * E2.7's files (DDX3X, TP53), c13_vcep_band_lr.csv: byte for byte as at HEAD;
  * c12_inframe_tp_fp.csv: the true-positive rows unchanged; the false-positive rows of the
    seven genes at 3-10 nt the only ones that change; the header's one sentence;
  * Supplementary Table S13 and the supplement's manifest: what changed, row by row.
Anything else changing stops (and check_scope.py checks every output against HEAD).

Run in the reference image (修改单 §0 rule 9), offline, after the pipeline from stage 5:
    docker run ... IMAGE bash scripts/stage_E20b/in_image.sh scripts/stage_E20b/check_e6_downstream.py
"""
from __future__ import annotations

import glob
import io
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog  # noqa: E402

REP = "phase1/reports/evidence"
SAME = [f"{REP}/inframe_attribution_ddx3x.csv", f"{REP}/inframe_attribution_ddx3x_summary.csv",
        f"{REP}/inframe_attribution_tp53.csv", f"{REP}/inframe_attribution_tp53_summary.csv",
        f"{REP}/suppl_analyses/c13_vcep_band_lr.csv"]
FRAMES = ["in_frame", "out_of_frame", "undetermined"]


def head(path: str) -> bytes:
    return subprocess.run(["git", "-C", str(REPO), "show", f"HEAD:{path}"], capture_output=True,
                          check=True).stdout


def table(data: bytes, comment: str | None = None) -> pd.DataFrame:
    return pd.read_csv(io.BytesIO(data), comment=comment, keep_default_na=False, dtype=str)


def main() -> None:
    failures, notes = [], {}
    rec = json.loads(Path(sorted(glob.glob(str(REPO / "scripts/stage_E20b/run_log/write_e6_subsets_*.json")))[-1]).read_text())
    added = {t: set(rec["notes"]["checks"][t]["added"]) for t in ("spliceai", "pangolin")}
    with RunLog(__file__) as log:
        # E6's per-variant attribution
        path = f"{REP}/inframe_attribution.csv"
        log.input(REPO / path)
        old, new = table(head(path)), table((REPO / path).read_bytes())
        key = ["tool", "variant_id"]
        o, n = old.set_index(key), new.set_index(key)
        new_only = n.loc[~n.index.isin(o.index)]
        unchanged = n.loc[o.index].equals(o) and list(n.columns) == list(o.columns)
        got_added = {t: set(new_only.loc[t].index) if t in new_only.index.get_level_values(0) else set()
                     for t in added}
        notes["attribution"] = {"head_rows": len(old), "new_rows": len(new),
                                "head_rows_unchanged": bool(unchanged),
                                "new_rows_are_the_added_variants": got_added == added,
                                "added_frames": {t: new_only.loc[t]["frame"].value_counts().to_dict()
                                                 for t in added if t in new_only.index.get_level_values(0)}}
        if not unchanged or got_added != added or (o.index.difference(n.index).size):
            failures.append("inframe_attribution.csv: not HEAD's rows plus the added variants")
        # the summary
        path = f"{REP}/inframe_attribution_summary.csv"
        log.input(REPO / path)
        so, sn = table(head(path)).set_index(["tool", "stratum"]), table((REPO / path).read_bytes()).set_index(["tool", "stratum"])
        changed = [i for i in sn.index if not sn.loc[i].equals(so.loc[i])]
        notes["summary_changed_cells"] = {f"{t} {s}": {"head": so.loc[(t, s)].to_dict(), "new": sn.loc[(t, s)].to_dict()}
                                          for t, s in changed}
        if set(s for _, s in changed) - {"s3_10"} or set(so.index) != set(sn.index):
            failures.append("the summary changes outside the 3-10 nt cells")
        for t, s in changed:
            frames = notes["attribution"]["added_frames"].get(t, {})
            for f in FRAMES:
                if int(sn.loc[(t, s), f]) - int(so.loc[(t, s), f]) != int(frames.get(f, 0)):
                    failures.append(f"summary {t} {s} {f}: the change is not the added variants' count")
        # E2.7 and C13: byte for byte
        notes["unchanged_files"] = {}
        for path in SAME:
            log.input(REPO / path)
            same = head(path) == (REPO / path).read_bytes()
            notes["unchanged_files"][path] = same
            if not same:
                failures.append(f"{path} changed")
        # C12: true-positive rows unchanged, false positives only at seven genes 3-10 nt
        path = f"{REP}/suppl_analyses/c12_inframe_tp_fp.csv"
        log.input(REPO / path)
        co, cn = table(head(path), "#"), table((REPO / path).read_bytes(), "#")
        k = ["gene_set", "label_definition", "tool", "stratum", "group"]
        co, cn = co.set_index(k), cn.set_index(k)
        diff = [i for i in cn.index if not cn.loc[i].equals(co.loc[i])]
        notes["c12_changed_rows"] = {" ".join(i): {"head": co.loc[i].to_dict(), "new": cn.loc[i].to_dict()} for i in diff}
        if any(i[4] != "false_positive" or i[0] != "seven genes" or i[3] != "s3_10" for i in diff) \
                or set(co.index) != set(cn.index):
            failures.append("c12: a row changed beyond the seven genes' false positives at 3-10 nt")
        hl = [x for x in head(path).decode().splitlines() if x.startswith("#")]
        nl = [x for x in (REPO / path).read_text().splitlines() if x.startswith("#")]
        notes["c12_header_lines_changed"] = [{"head": a, "new": b} for a, b in zip(hl, nl) if a != b]
        # S13 and the supplement manifest
        for path in (f"{REP}/supplement/tableS13_inframe_attribution.csv",
                     f"{REP}/supplement/supplement_tables_manifest.csv"):
            log.input(REPO / path)
            a, b = table(head(path)), table((REPO / path).read_bytes())
            rows = [{"head": a.iloc[i].to_dict(), "new": b.iloc[i].to_dict()}
                    for i in range(min(len(a), len(b))) if not a.iloc[i].equals(b.iloc[i])]
            notes[f"changed_rows {Path(path).name}"] = rows
            if len(a) != len(b):
                failures.append(f"{path}: the number of rows changed")
        log.note("checks", notes)
        log.note("failures", failures)
        print(json.dumps(notes, indent=1, ensure_ascii=False)[:6000])
        if failures:
            for f in failures:
                print("FAIL", f)
            raise SystemExit("unexpected changes downstream of E6: stopped")
        print("downstream of E6: only the expected changes")


if __name__ == "__main__":
    main()
