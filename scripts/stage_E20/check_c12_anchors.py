"""E20.1 (C12): the two anchors of the true positives' event records. Any failure stops.

Anchor 1, the score column. The relation between an event record's largest score and the
stored score column is read first on the existing records (E6's and E2.7's), then
required of every row of the new ones:
  * SpliceAI: spliceai_event_max equals the column spliceai_walker bit for bit (the
    analysis set's for the seven genes, the scored table's for DDX3X and for TP53);
  * Pangolin: pangolin_event_max equals the analysis set's pangolin rounded as Pangolin's
    CLI prints it, round() on the float32 value to two decimals (the scorer reads the
    CLI's output; the column was kept at full precision) -- output precision only.

Anchor 2, the existing events. The 10 SpliceAI and 10 Pangolin variants
write_c12_subsets.py drew from E6's event files, scored in the same runs as the true
positives, must give E6's records back in every field, bit for bit: the four delta scores
and their four positions, the largest and the gene symbol (SpliceAI); the gain and loss
scores, their positions and the largest (Pangolin), and the CLI's raw record line, read
from E6's work folder and the new one.

Run in the reference image (修改单 §0 rule 9), offline:
    docker run ... IMAGE bash scripts/stage_E20/in_image.sh scripts/stage_E20/check_c12_anchors.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog  # noqa: E402

PHASE1 = REPO / "phase1"
sys.path.insert(0, str(PHASE1))
os.chdir(PHASE1)
from src import evid_common as K  # noqa: E402

EV = Path("data/evidence")
EXT = EV / "external"
SAI_FIELDS = ["spliceai_event_max", "spliceai_event_symbol", "spliceai_ds_ag", "spliceai_ds_al",
              "spliceai_ds_dg", "spliceai_ds_dl", "spliceai_dp_ag", "spliceai_dp_al",
              "spliceai_dp_dg", "spliceai_dp_dl"]
PANG_FIELDS = ["pangolin_event_max", "pangolin_gain", "pangolin_gain_pos", "pangolin_loss",
               "pangolin_loss_pos"]
ID_FIELDS = ["variant_id", "gene", "chrom", "pos", "ref", "alt"]
# (existing events, new events, score input, scored table, column, kind)
SETS = {
    "seven_genes_spliceai": (EV / "inframe_events.parquet", EV / "inframe_events_tp.parquet",
                             EV / "inframe_score_input_tp.parquet", None, "spliceai_walker", "bits"),
    "seven_genes_pangolin": (EV / "inframe_events_pangolin.parquet",
                             EV / "inframe_events_pangolin_tp.parquet",
                             EV / "inframe_score_input_pangolin_tp.parquet", None, "pangolin", "cli"),
    "ddx3x_spliceai": (EXT / "ddx3x_inframe_events.parquet", EXT / "ddx3x_inframe_events_tp.parquet",
                       EXT / "ddx3x_inframe_subset_tp.parquet", EXT / "ddx3x_scored.parquet",
                       "spliceai_walker", "bits"),
    "tp53_spliceai": (EXT / "tp53_inframe_events_12nt.parquet",
                      EXT / "tp53_inframe_events_12nt_tp.parquet",
                      EXT / "tp53_inframe_subset_12nt_tp.parquet",
                      EXT / "tp53_splice_scored_12nt.parquet", "spliceai_walker", "bits"),
}
WORK = {"old": EV / "pangolin_events_work", "new": EV / "pangolin_events_work_tp"}


def cli_round(x: float) -> float:
    """Pangolin's CLI: round(gain[g], 2) on a float32, printed and read back as a float."""
    return float(round(np.float32(x), 2))


def holds(events: pd.DataFrame, table: pd.DataFrame, col: str, kind: str) -> dict:
    m = events.merge(table[["variant_id", col]], on="variant_id", how="left", validate="1:1")
    ev_col = "spliceai_event_max" if kind == "bits" else "pangolin_event_max"
    a = m[ev_col].to_numpy(dtype=float)
    b = m[col].to_numpy(dtype=float)
    want = b if kind == "bits" else np.array([cli_round(x) for x in b])
    ok = (a == want) & np.isfinite(a)
    return {"n": int(len(m)), "n_holding": int(ok.sum()),
            "failing": m.loc[~ok, "variant_id"].tolist(),
            "max_abs_difference_from_column": float(np.nanmax(np.abs(a - b))) if len(m) else 0.0}


def same_bits(a: pd.Series, b: pd.Series) -> bool:
    if a.dtype.kind == "f" or b.dtype.kind == "f":
        x, y = a.to_numpy(dtype=float), b.to_numpy(dtype=float)
        return bool(np.array_equal(x.view(np.int64), y.view(np.int64)))
    return bool((a.astype(str).to_numpy() == b.astype(str).to_numpy()).all())


def raw_lines(folder: Path) -> dict:
    out = {}
    for p in sorted((folder / "chunks").glob("chunk_*.out.csv")):
        for line in p.read_text().splitlines()[1:]:
            chrom, pos, ref, alt, _ = line.split(",", 4)
            out[(chrom, int(pos), ref, alt)] = line
    return out


def main() -> None:
    failures = []
    with RunLog(__file__) as log:
        log.input(K.SET_PATH)
        df = K.load_set()
        for name, (old, new, inp, scored, col, kind) in SETS.items():
            for p in (old, new, inp) + ((scored,) if scored else ()):
                log.input(p)
            table = df if scored is None else pd.read_parquet(scored)
            ev_old, ev_new = pd.read_parquet(old), pd.read_parquet(new)
            rec_old = holds(ev_old, table, col, kind)
            rec_new = holds(ev_new, table, col, kind)
            if rec_old["n_holding"] != rec_old["n"]:
                failures.append(f"{name}: the existing records do not settle the relation")
            if rec_new["n_holding"] != rec_new["n"]:
                failures.append(f"{name}: {len(rec_new['failing'])} new rows break the relation")
            log.note(f"anchor1_{name}", {"relation": "bit for bit" if kind == "bits" else
                                         "the column rounded as Pangolin's CLI prints it",
                                         "existing": rec_old, "new": rec_new})
            roles = pd.read_parquet(inp)
            if "role" not in roles.columns:
                continue
            anc = roles.loc[roles["role"] == "anchor", "variant_id"]
            a_old = ev_old.set_index("variant_id").loc[anc].reset_index()
            a_new = ev_new.set_index("variant_id").loc[anc].reset_index()
            fields = ID_FIELDS + (SAI_FIELDS if kind == "bits" else PANG_FIELDS)
            per_field = {f: same_bits(a_old[f], a_new[f]) for f in fields}
            rec = {"n": int(len(anc)), "fields_identical": per_field}
            if kind == "cli":
                for k, folder in WORK.items():
                    for p in sorted(folder.rglob("*")):
                        if p.is_file():
                            log.input(p)
                lo, ln = raw_lines(WORK["old"]), raw_lines(WORK["new"])
                keys = [(str(r.chrom), int(r.pos), r.ref, r.alt) for r in a_old.itertuples()]
                rec["cli_lines_identical"] = sum(1 for k in keys
                                                 if ln.get(k) is not None and lo.get(k) == ln.get(k))
                if rec["cli_lines_identical"] != len(keys):
                    failures.append(f"{name}: an anchor's CLI line differs from E6's")
            if not all(per_field.values()):
                failures.append(f"{name}: an anchor's event record differs from E6's")
            log.note(f"anchor2_{name}", rec)
        log.note("failures", failures)
        for f in failures:
            print("FAIL", f)
        if failures:
            raise SystemExit("an anchor does not hold: stopped")
        print("both anchors hold on every event file")


if __name__ == "__main__":
    main()
