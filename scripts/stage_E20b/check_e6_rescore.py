"""E20b.1 (commit 1): the two anchors of E6's re-scored event records, and, only when both
hold, the records put in place of the committed ones.

Anchor 1, the committed records. Every variant of E6's committed event files (SpliceAI 609,
Pangolin 524) is in the new run with every field of its record identical bit for bit: the
four delta scores and their four positions, the largest and the gene symbol (SpliceAI);
the gain and loss scores, their positions and the largest (Pangolin), and the CLI's raw
record line, read from the committed work folder and the new one.

Anchor 2, the score columns. Every variant of the new subsets stands in E20.1's relation
to the analysis set: SpliceAI's largest delta equals spliceai_walker bit for bit;
Pangolin's largest value equals pangolin rounded as Pangolin's CLI prints it (round() on
the float32 value to two decimals).

With --install, and only if both hold, the new event files, their provenance and
environment records and Pangolin's work folder replace the committed ones (the work folder
whole), and phase1/data/e20b_staging/ is removed; the sha256 of every file before and
after goes into the run record. A failure stops before anything is replaced.

Run in the reference image (修改单 §0 rule 9), offline:
    docker run ... IMAGE bash scripts/stage_E20b/in_image.sh scripts/stage_E20b/check_e6_rescore.py --install
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel, sha256_of  # noqa: E402

PHASE1 = REPO / "phase1"
sys.path.insert(0, str(PHASE1))
os.chdir(PHASE1)
from src import evid_common as K  # noqa: E402

EV = PHASE1 / "data/evidence"
STAGING = PHASE1 / "data/e20b_staging"
SAI = ["variant_id", "gene", "chrom", "pos", "ref", "alt", "spliceai_event_max", "spliceai_event_symbol",
       "spliceai_ds_ag", "spliceai_ds_al", "spliceai_ds_dg", "spliceai_ds_dl",
       "spliceai_dp_ag", "spliceai_dp_al", "spliceai_dp_dg", "spliceai_dp_dl"]
PANG = ["variant_id", "gene", "chrom", "pos", "ref", "alt", "pangolin_event_max", "pangolin_gain",
        "pangolin_gain_pos", "pangolin_loss", "pangolin_loss_pos"]
RUNS = {"spliceai": ("inframe_events.parquet", "inframe_subset.parquet", SAI, "spliceai_walker"),
        "pangolin": ("inframe_events_pangolin.parquet", "inframe_subset_pangolin.parquet", PANG, "pangolin")}
WORK = "pangolin_events_work"


def same_bits(a: pd.Series, b: pd.Series) -> bool:
    if a.dtype.kind == "f" or b.dtype.kind == "f":
        x, y = a.to_numpy(dtype=float), b.to_numpy(dtype=float)
        return bool(np.array_equal(x.view(np.int64), y.view(np.int64)))
    return bool((a.astype(str).to_numpy() == b.astype(str).to_numpy()).all())


def cli_lines(folder: Path) -> dict:
    out = {}
    for p in sorted((folder / "chunks").glob("chunk_*.out.csv")):
        for line in p.read_text().splitlines()[1:]:
            chrom, pos, ref, alt, _ = line.split(",", 4)
            out[(chrom, int(pos), ref, alt)] = line
    return out


def cli_round(x: float) -> float:
    return float(round(np.float32(x), 2))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--install", action="store_true")
    install = ap.parse_args().install
    failures, record = [], {}
    with RunLog(__file__) as log:
        log.input(K.SET_PATH)
        df = K.load_set().set_index("variant_id")
        for tool, (events, subset, fields, col) in RUNS.items():
            old_p, new_p, sub_p = EV / events, STAGING / events, EV / subset
            for p in (old_p, new_p, sub_p):
                log.input(p)
            old, new, sub = pd.read_parquet(old_p), pd.read_parquet(new_p), pd.read_parquet(sub_p)
            rec = {"committed_rows": int(len(old)), "new_rows": int(len(new)), "subset_rows": int(len(sub))}
            if sorted(new["variant_id"]) != sorted(sub["variant_id"]) or new["variant_id"].duplicated().any():
                failures.append(f"{tool}: the new events are not the new subset, once each")
            # anchor 1: every committed record, every field
            keep = new.set_index("variant_id").loc[old["variant_id"]].reset_index()
            rec["anchor1_fields_identical"] = {f: same_bits(old[f], keep[f]) for f in fields}
            if not all(rec["anchor1_fields_identical"].values()):
                failures.append(f"{tool}: a committed record differs in the new run")
            if tool == "pangolin":
                lo, ln = cli_lines(EV / WORK), cli_lines(STAGING / WORK)
                keys = [(str(r.chrom), int(r.pos), r.ref, r.alt) for r in old.itertuples()]
                rec["anchor1_cli_lines_identical"] = sum(
                    1 for k in keys if ln.get(k) is not None and lo.get(k) == ln.get(k))
                if rec["anchor1_cli_lines_identical"] != len(keys):
                    failures.append("pangolin: a committed CLI line differs in the new run")
            # anchor 2: every variant of the new subset against the score column
            stored = df.loc[new["variant_id"], col].to_numpy(dtype=float)
            got = new[fields[6]].to_numpy(dtype=float)
            want = stored if tool == "spliceai" else np.array([cli_round(x) for x in stored])
            ok = (got == want) & np.isfinite(got)
            rec["anchor2_holding"] = int(ok.sum())
            rec["anchor2_failing"] = new.loc[~ok, "variant_id"].tolist()
            if not ok.all():
                failures.append(f"{tool}: {int((~ok).sum())} new rows break the relation to {col}")
            record[tool] = rec
        log.note("anchors", record)
        log.note("failures", failures)
        if failures:
            for f in failures:
                print("FAIL", f)
            raise SystemExit("an anchor does not hold: nothing is replaced")
        print("both anchors hold")
        if not install:
            return
        moved = []
        for tool, (events, *_rest) in RUNS.items():
            for suffix in ("", ".provenance.json", ".environment.json"):
                src, dst = STAGING / (events + suffix), EV / (events + suffix)
                before = sha256_of(dst) if dst.exists() else None
                shutil.copyfile(src, dst)
                moved.append({"file": rel(dst), "sha256_before": before, "sha256_after": sha256_of(dst)})
                log.output(dst)
        old_work = sorted(p for p in (EV / WORK).rglob("*") if p.is_file())
        before = {rel(p): sha256_of(p) for p in old_work}
        shutil.rmtree(EV / WORK)
        shutil.copytree(STAGING / WORK, EV / WORK)
        after = {rel(p): sha256_of(p) for p in sorted((EV / WORK).rglob("*")) if p.is_file()}
        for p in sorted((EV / WORK).rglob("*")):
            if p.is_file():
                log.output(p)
        log.note("moved", moved)
        log.note("work_folder", {"before": before, "after": after})
        shutil.rmtree(STAGING)
        print(f"installed; {rel(STAGING)} removed")


if __name__ == "__main__":
    main()
