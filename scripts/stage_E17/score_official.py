"""E17.3, E17.4: AlphaGenome's splicing score by its official definition, for the seven
training genes of the analysis set only.

The scoring code is the notebook's own, as quoted in
phase1/data/evidence/alphagenome_official/NOTEBOOK_SOURCE.json, and it is executed from
there: cell 5 (the imports) and cell 9 (the three scorers) as written; per variant, the
lines of cell 15 that score one variant (an interval of dna_client.SEQUENCE_LENGTH_1MB
centred on it, dna_model.score_variant with the three scorers, Homo sapiens), taken out
of the notebook's loop so that two threads can each score one; then
variant_scorers.tidy_scores and cell 13's compute_merged_splicing_score, as written. The
model is created as in cell 7 except for the key, which comes from
phase1/src/evid_score_avi.read_key() (E17.0 rule 1). The key is never printed, logged or
recorded: only whether one was found and where it came from.

Modes
  --examples  the notebook's three example variants (cell 15). The notebook prints no
              output, so there is nothing to compare them with; they show that the
              quoted code runs as written on this client.
  --pilot     120 variants of the analysis set, an even share in each non-empty gene x
              stratum cell (pm12, s3_10, s11_50), drawn with config.RANDOM_SEED; and,
              drawn from those, the 10 the determinism check scores again.
  --repeat    those 10 again, later: they must equal the pilot's scores bit for bit.
  --all       every row of the analysis set (commit 2, after gate 1).

Scoring is gentle: --workers threads at most (2 by default), a pause and exponential
backoff on errors, a stop on an authentication error or on quota errors that do not
clear. Each scored variant is appended to a cache in phase1/data/alphagenome_official_cache/
(gitignored), so a run can be resumed. Outputs, in alphagenome_official/: the scores
(<mode>_scores.tsv: variant_id, the three components, the composite, and per component
the rows scored and their smallest raw score), a provenance record beside it, and
ENVIRONMENT.json.

Run from the repository root:
    .venv-ag-official/bin/python scripts/stage_E17/score_official.py --pilot
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel, sha256_of, utc_now  # noqa: E402

PHASE1 = REPO / "phase1"
sys.path.insert(0, str(PHASE1))
from src import config as C                                # noqa: E402
from src.evid_score_avi import KEY_FILE, read_key          # noqa: E402

DIR = PHASE1 / "data" / "evidence" / "alphagenome_official"
CACHE = PHASE1 / "data" / "alphagenome_official_cache"
SOURCE = DIR / "NOTEBOOK_SOURCE.json"
SET = PHASE1 / "data" / "evidence" / "analysis_set_v1.parquet"
TRAINING = sorted(C.ASSAY_LABELS)                          # the seven genes
STRATA = ["pm12", "s3_10", "s11_50"]
PILOT_N, REPEAT_N = 120, 10
COMPONENTS = ["SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS"]
COMPOSITE = "alphagenome_splicing"
MAX_ATTEMPTS, BACKOFF = 6, 5.0
FATAL = ("UNAUTHENTICATED", "PERMISSION_DENIED", "API key", "api key")


def key_report() -> dict:
    key = read_key()
    found = bool(key)
    del key
    return {"api_key_present": found,
            "api_key_source": ("ALPHAGENOME_API_KEY" if os.environ.get("ALPHAGENOME_API_KEY")
                               else "~/" + str(KEY_FILE.relative_to(Path.home())) if found
                               else None)}


def _clean(msg: str) -> str:
    """An error message with the key, should it ever appear, cut out."""
    key = read_key()
    if key and key in msg:
        msg = msg.replace(key, "<key>")
    del key
    return msg[:300]


class Notebook:
    """The quoted cells, executed as written; the model created with read_key()."""

    def __init__(self):
        src = json.loads(SOURCE.read_text())
        self.source_sha256 = sha256_of(SOURCE)
        cells = {c["cell"]: c["source"] for c in src["quoted_cells"]}
        self.ns: dict = {}
        exec(cells[5], self.ns)                                    # imports
        key = read_key()
        if not key:
            raise SystemExit("no AlphaGenome API key (ALPHAGENOME_API_KEY or "
                             "~/.config/alphagenome/key): stopping (E17.0 rule 4)")
        self.ns["dna_model"] = self.ns["dna_client"].create(key)   # cell 7, the key aside
        del key
        exec(cells[9], self.ns)                                    # the three scorers
        lines = cells[15].split("\n")
        start = next(i for i, l in enumerate(lines) if l.startswith("for variant in variants:"))
        end = next(i for i, l in enumerate(lines) if l.strip() == "all_scores.append(scores)")
        body = lines[start + 1:end]                                # interval = ...; scores = ...
        self.body = "\n".join(body)
        exec("def score_one(variant):\n" + self.body + "\n  return scores\n", self.ns)
        merged = cells[13].split("\n\n\nmerged_scores")[0]         # the function, not its call
        exec(merged, self.ns)

    def variant(self, chrom: str, pos: int, ref: str, alt: str):
        return self.ns["genome"].Variant(chromosome=f"chr{chrom}", position=int(pos),
                                         reference_bases=ref, alternate_bases=alt)

    def score(self, variant) -> dict:
        scores = self.ns["score_one"](variant)
        df = self.ns["variant_scorers"].tidy_scores([scores])
        merged = self.ns["compute_merged_splicing_score"](df.copy())
        if len(merged) != 1:
            raise RuntimeError(f"{len(merged)} merged rows for one variant")
        row = merged.iloc[0]
        out = {c: float(row[c]) if c in merged.columns else 0.0 for c in COMPONENTS}
        out[COMPOSITE] = float(row[COMPOSITE])
        for c in COMPONENTS:
            part = df.loc[df["output_type"] == c, "raw_score"]
            out[f"rows_{c}"] = int(len(part))
            out[f"nan_{c}"] = int(part.isna().sum())
            out[f"min_raw_{c}"] = float(part.min()) if part.notna().any() else None
        return out


def pilot_sample(df):
    """An even share of PILOT_N in each non-empty gene x stratum cell, drawn with the
    project seed in a fixed cell order; a cell short of its share gives the rest to the
    others. Then REPEAT_N of those, for the determinism check."""
    import numpy as np
    rng = np.random.default_rng(C.RANDOM_SEED)
    cells = [(g, s) for g in TRAINING for s in STRATA
             if ((df.gene == g) & (df.stratum == s)).any()]
    size = {c: int(((df.gene == c[0]) & (df.stratum == c[1])).sum()) for c in cells}
    share, left = {c: 0 for c in cells}, PILOT_N
    while left:
        open_ = [c for c in cells if share[c] < size[c]]
        if not open_:
            raise SystemExit("the analysis set holds fewer than the pilot's variants")
        each = max(1, left // len(open_))
        for c in open_:
            take = min(each, size[c] - share[c], left)
            share[c] += take
            left -= take
            if not left:
                break
    picked = []
    for c in cells:
        ids = df.loc[(df.gene == c[0]) & (df.stratum == c[1]), "variant_id"].to_numpy()
        picked += list(rng.choice(ids, size=share[c], replace=False))
    repeat = sorted(rng.choice(np.array(sorted(picked)), size=REPEAT_N, replace=False))
    return sorted(picked), list(repeat), {f"{g}/{s}": n for (g, s), n in share.items()}


def environment(nb: Notebook, log) -> None:
    import alphagenome
    freeze = sorted(subprocess.run([sys.executable, "-m", "pip", "freeze"], check=True,
                                   capture_output=True, text=True).stdout.split())
    env = {
        "written_by": rel(Path(__file__).resolve()),
        "written_utc": utc_now(),
        "python": sys.version,
        "platform": platform.platform(),
        "pip_freeze": freeze,
        "requirements": "scripts/stage_E17/requirements.txt",
        "requirements_sha256": sha256_of(REPO / "scripts/stage_E17/requirements.txt"),
        "client": {"package": "alphagenome", "version": alphagenome.__version__},
        "model": {"model_version_requested": None,
                  "model_version_requested_note": "dna_client.create(...) as in the "
                  "notebook's cell 7, with no model_version: the service's default",
                  "model_version_in_response": "not_provided",
                  "response_metadata_keys": nb.uns_keys},
    } | key_report()
    (DIR / "ENVIRONMENT.json").write_text(json.dumps(env, indent=1, ensure_ascii=False) + "\n")
    log.output(DIR / "ENVIRONMENT.json")


def main() -> None:
    ap = argparse.ArgumentParser()
    mode = ap.add_mutually_exclusive_group(required=True)
    for m in ("examples", "pilot", "repeat", "all"):
        mode.add_argument(f"--{m}", dest="mode", action="store_const", const=m)
    ap.add_argument("--workers", type=int, default=2)
    args = ap.parse_args()
    import pandas as pd
    with RunLog(__file__) as log:
        log.note("key", key_report())
        DIR.mkdir(parents=True, exist_ok=True)
        CACHE.mkdir(parents=True, exist_ok=True)
        nb = Notebook()
        log.input(SOURCE)
        if args.mode == "examples":
            g = dict(nb.ns)
            src = {c["cell"]: c["source"] for c in json.loads(SOURCE.read_text())["quoted_cells"]}
            exec(src[15].split("\n\n# Score each variant")[0], g)      # the variants, as written
            todo = pd.DataFrame([{"variant_id": str(v), "gene": "", "stratum": "",
                                  "chrom": v.chromosome.removeprefix("chr"), "pos": v.position,
                                  "ref": v.reference_bases, "alt": v.alternate_bases}
                                 for v in g["variants"]])
        else:
            df = pd.read_parquet(SET, columns=["variant_id", "gene", "stratum", "chrom",
                                                "pos", "ref", "alt"])
            log.input(SET)
            if set(df.gene) - set(TRAINING):
                raise SystemExit(f"genes outside the seven: {set(df.gene) - set(TRAINING)}")
            pilot, repeat, share = pilot_sample(df)
            log.note("pilot_share", share)
            ids = {"pilot": pilot, "repeat": repeat, "all": df.variant_id.tolist()}[args.mode]
            todo = df.set_index("variant_id").loc[ids].reset_index()
        cache = CACHE / f"{args.mode}.jsonl"
        done = {}
        if cache.exists():
            for line in cache.read_text().splitlines():
                r = json.loads(line)
                done[r["variant_id"]] = r
        todo_rows = [r for r in todo.itertuples() if r.variant_id not in done]
        print(f"[{args.mode}] {len(todo)} variants, {len(done)} cached, "
              f"{len(todo_rows)} to score with {args.workers} worker(s)", flush=True)
        nb.uns_keys = None
        lock, stats = threading.Lock(), {"retries": 0, "failures": []}
        start = utc_now()
        t0 = time.time()

        def one(r):
            v = nb.variant(r.chrom, r.pos, r.ref, r.alt)
            for attempt in range(1, MAX_ATTEMPTS + 1):
                try:
                    out = nb.score(v)
                    break
                except Exception as e:                          # noqa: BLE001
                    msg = _clean(f"{type(e).__name__}: {e}")
                    if any(f in msg for f in FATAL):
                        raise SystemExit(f"authentication error: stopping (E17.0 rule 4): "
                                         f"{msg}")
                    if attempt == MAX_ATTEMPTS:
                        with lock:
                            stats["failures"].append({"variant_id": r.variant_id,
                                                      "error": msg})
                        return None
                    with lock:
                        stats["retries"] += 1
                    time.sleep(BACKOFF * 2 ** (attempt - 1))
            rec = {"variant_id": r.variant_id, "scored_utc": utc_now()} | out
            with lock:
                with open(cache, "a") as f:
                    f.write(json.dumps(rec) + "\n")
                done[r.variant_id] = rec
                if len(done) % 20 == 0:
                    print(f"  {len(done)} scored, {time.time() - t0:.0f} s", flush=True)
            return rec

        # the response's own metadata, once, for ENVIRONMENT.json
        probe = nb.ns["score_one"](nb.variant(*todo.iloc[0][["chrom", "pos", "ref", "alt"]]))
        nb.uns_keys = sorted({k for ad in probe for k in ad.uns.keys()})
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            for fut in as_completed([pool.submit(one, r) for r in todo_rows]):
                fut.result()
        elapsed = time.time() - t0
        if stats["failures"]:
            raise SystemExit(f"{len(stats['failures'])} variants failed after "
                             f"{MAX_ATTEMPTS} attempts: {stats['failures'][:3]}")
        rows = [done[v] for v in todo.variant_id]
        table = todo[["variant_id", "gene", "stratum"]].merge(pd.DataFrame(rows),
                                                              on="variant_id")
        cols = (["variant_id", "gene", "stratum"] + COMPONENTS + [COMPOSITE]
                + [f"{p}_{c}" for c in COMPONENTS for p in ("rows", "nan", "min_raw")]
                + ["scored_utc"])
        out = DIR / f"{args.mode}_scores.tsv"
        table[cols].to_csv(out, sep="\t", index=False,          # exact, and read back as such
                           float_format=lambda x: repr(float(x)))
        log.output(out)
        environment(nb, log)
        prov = {
            "written_by": rel(Path(__file__).resolve()),
            "mode": args.mode,
            "definition": rel(SOURCE), "definition_sha256": nb.source_sha256,
            "code_run": {"cells_as_written": [5, 9, 13, "15 (the lines scoring one "
                                              "variant, out of the loop)"],
                         "cell_7": "dna_client.create(read_key()): the key aside, as written",
                         "per_variant_lines": nb.body},
            "client_version": __import__("alphagenome").__version__,
            "model_version_requested": None, "model_version_in_response": "not_provided",
            "start_utc": start, "end_utc": utc_now(),
            "seconds_scoring": round(elapsed, 1),
            "n_variants": len(table), "n_scored_this_run": len(todo_rows),
            "n_from_cache": len(todo) - len(todo_rows),
            "retries": stats["retries"], "failures": stats["failures"],
            "workers": args.workers,
            "sample": ({"seed": C.RANDOM_SEED, "pilot_n": PILOT_N, "repeat_n": REPEAT_N,
                        "share": log.notes.get("pilot_share")}
                       if args.mode in ("pilot", "repeat") else None),
            "genes": sorted(set(table.gene) - {""}),
            "output": rel(out), "output_sha256": sha256_of(out),
        } | key_report()
        (DIR / f"{args.mode}_scores.tsv.provenance.json").write_text(
            json.dumps(prov, indent=1, ensure_ascii=False) + "\n")
        log.output(DIR / f"{args.mode}_scores.tsv.provenance.json")
        log.note("seconds_scoring", round(elapsed, 1))
        log.note("retries", stats["retries"])
        print(f"[{args.mode}] wrote {rel(out)}: {len(table)} variants, "
              f"{elapsed:.0f} s for {len(todo_rows)}, {stats['retries']} retries")


if __name__ == "__main__":
    main()
