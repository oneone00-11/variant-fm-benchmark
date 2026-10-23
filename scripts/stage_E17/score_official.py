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
  --all       every row of the analysis set (commit 2, after gate 1). Its own cache
              means the pilot's 120 are scored again; they are compared with the
              pilot's scores bit for bit (gate 1 review, R4).
  --drift     the determinism check's 10 once more, after --all, from an empty cache
              of its own, so that each is a fresh call: they must equal the pilot's
              scores bit for bit, or the service's model changed during the full run
              (E17.4; gate 1 review, R1).
After the thresholds were frozen (tag frozen-thresholds-v1; stage E18.1):
  --drift-2b-start  the same 10, from an empty cache of their own, before any external
              gene is scored: they must equal the pilot's scores bit for bit, or the
              service's model has changed since E17 and the external scores could not
              go with the frozen thresholds.
  --external  the external genes' variants (external_variants), only after the start
              check found its 10 identical.
  --drift-2b-end  the same 10 once more, after --external, from an empty cache: equal
              to the pilot's, bit for bit, and each scored after the external run ended.
These write ENVIRONMENT_2b.json, so that E17's ENVIRONMENT.json stays as tagged, and
every score they write must post-date the tag (its time is read from
phase1/data/evidence/FROZEN_AT_TAG.json).
A comparison that finds a difference, or a raw score below zero (the scorers are
unsigned, so the maximum would no longer be the maximum absolute score), is written to
the record and then stops the script (E17.10; E18.8).

Scoring is gentle: --workers threads at most (2 by default), and a pause before each new
attempt, doubling each time. The run stops (gate 1 review, R3) on an authentication
error (E17.0 rule 4), when three variants in a row have used all their attempts whatever
the error (E17.4: errors that do not clear), and on Ctrl+C: the variants not yet started
are cancelled, and a flag every worker checks before each call lets the calls in flight
end and no other begin. Each scored variant is appended to a cache in
phase1/data/alphagenome_official_cache/ (gitignored, one file per mode), so a stopped
run can be resumed. Outputs, in alphagenome_official/: the scores (<mode>_scores.tsv:
variant_id, the three components, the composite, and per component the rows scored and
their smallest raw score), a provenance record beside it, and ENVIRONMENT.json. Both
records quote the client's line naming the model a request with no model_version gets
(R6a); the request, as in the notebook, names none.

Run from the repository root:
    .venv-ag-official/bin/python scripts/stage_E17/score_official.py --pilot
"""
from __future__ import annotations

import argparse
import base64
import hashlib
import json
import os
import platform
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
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
REQ = REPO / "scripts" / "stage_E17" / "requirements.txt"
TRAINING = sorted(C.ASSAY_LABELS)                          # the seven genes
STRATA = ["pm12", "s3_10", "s11_50"]
PILOT_N, REPEAT_N = 120, 10
FULL = "all_scores.tsv"
# E18.1: after the tag, the external genes
EXT = PHASE1 / "data" / "evidence" / "external"
FROZEN = PHASE1 / "data" / "evidence" / "FROZEN_AT_TAG.json"
EXTERNAL = "external_scores.tsv"
NEW_GENES = ["CTCF", "SFPQ", "RAD51D", "XRCC2", "SBDS"]
E18_MODES = ("drift_2b_start", "external", "drift_2b_end")
COMPONENTS = ["SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS"]
COMPOSITE = "alphagenome_splicing"
MAX_ATTEMPTS, BACKOFF = 6, 5.0
IN_A_ROW = 3              # variants in a row that use all their attempts: the run stops
FATAL = ("UNAUTHENTICATED", "PERMISSION_DENIED", "API key", "api key")
DNA_MODEL = "alphagenome/models/v1/dna_model.py"
DEFAULT_LINE = "ALL_FOLDS = enum.auto()  # Default model version if None explicitly set."


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


class Stopped(SystemExit):
    """The run stopped before its end; .state holds what it had done."""

    def __init__(self, why: str, state: dict):
        super().__init__(why)
        self.state = state


def run(rows, score, record, *, workers: int, attempts: int = MAX_ATTEMPTS,
        backoff: float = BACKOFF, in_a_row: int = IN_A_ROW,
        stop: threading.Event | None = None) -> dict:
    """score(row) for every row, `workers` at a time, and record(row, result) for each
    result. An error is tried again after a pause (backoff seconds, doubling), up to
    `attempts` calls; a variant that uses them all is a failure, and the caller stops
    after the run.

    The run stops early, raising Stopped, on an authentication error (E17.0 rule 4),
    when `in_a_row` variants in a row have used all their attempts, whatever the error
    (E17.4: errors that do not clear), and on Ctrl+C. The variants not yet started are
    then cancelled; `stop`, which every worker checks before each call and which ends any
    pause, lets no new call begin, so only the calls already in flight end, at most one
    per worker. Returns the calls, retries and failures."""
    stop = stop or threading.Event()
    lock = threading.Lock()
    state = {"calls": 0, "retries": 0, "failures": [], "in_a_row": 0, "why": None}

    def halt(why: str) -> None:                # with the lock held
        state["why"] = state["why"] or why
        stop.set()

    def one(r) -> None:
        for attempt in range(1, attempts + 1):
            if stop.is_set():
                return
            with lock:
                state["calls"] += 1
            try:
                out = score(r)
            except Exception as e:                              # noqa: BLE001
                msg = _clean(f"{type(e).__name__}: {e}")
                if any(f in msg for f in FATAL):
                    with lock:
                        halt(f"authentication error: stopping (E17.0 rule 4): {msg}")
                    return
                if attempt < attempts:
                    with lock:
                        state["retries"] += 1
                    stop.wait(backoff * 2 ** (attempt - 1))
                    continue
                with lock:
                    state["failures"].append({"variant_id": r.variant_id, "error": msg})
                    state["in_a_row"] += 1
                    if state["in_a_row"] >= in_a_row:
                        halt(f"{in_a_row} variants in a row failed after {attempts} "
                             f"attempts each: stopping (E17.4); the cache keeps every "
                             f"variant scored")
                return
            with lock:
                state["in_a_row"] = 0
            record(r, out)
            return

    pool = ThreadPoolExecutor(max_workers=workers)
    try:
        for fut in as_completed([pool.submit(one, r) for r in rows]):
            fut.result()                       # an error of this script's own ends the run
            if stop.is_set():
                break
    except KeyboardInterrupt:
        with lock:
            halt("interrupted (Ctrl+C): stopping; the cache keeps every variant scored")
    finally:
        stop.set()                             # no call begins after this
        pool.shutdown(wait=True, cancel_futures=True)
    if state["why"]:
        raise Stopped(state["why"], state)
    return state


def differing(path: Path, ref: Path, ids: list[str]) -> list[str]:
    """The variants among ids whose components or composite, as written, are not those
    of ref to the last character. Both files are written with repr, so to the last bit."""
    import pandas as pd
    cols = COMPONENTS + [COMPOSITE]
    a, b = (pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
            .set_index("variant_id")[cols] for p in (path, ref))
    return [v for v in ids if list(a.loc[v]) != list(b.loc[v])]


def model_default() -> dict:
    """Which model a request with no model_version gets, in the client's own words (gate
    1 review, R6a): the line of dna_model.py that says so, with its number; the file's
    sha256, which must be the one the wheel's RECORD lists for it; and the wheel's sha256
    as requirements.txt pins it."""
    from importlib import metadata
    dist = metadata.distribution("alphagenome")
    entry = next(f for f in dist.files if str(f) == DNA_MODEL)
    data = Path(entry.locate()).read_bytes()
    digest = hashlib.sha256(data)
    if (entry.hash.mode, entry.hash.value) != (
            "sha256", base64.urlsafe_b64encode(digest.digest()).rstrip(b"=").decode()):
        raise SystemExit(f"{DNA_MODEL} is not the file the alphagenome wheel installed")
    lines = data.decode().split("\n")
    hits = [n for n, line in enumerate(lines, 1) if DEFAULT_LINE in line]
    pins = [l for l in REQ.read_text().splitlines()
            if l.startswith(f"alphagenome=={dist.version} --hash=sha256:")]
    if len(hits) != 1 or len(pins) != 1:
        raise SystemExit(f"{len(hits)} lines of {DNA_MODEL} hold {DEFAULT_LINE!r}, and "
                         f"{len(pins)} lines of {rel(REQ)} pin alphagenome {dist.version}")
    return {"file": DNA_MODEL, "line": hits[0], "text": lines[hits[0] - 1],
            "file_sha256": digest.hexdigest(), "file_sha256_in_wheel_record": True,
            "wheel": f"alphagenome=={dist.version}",
            "wheel_sha256": pins[0].split("--hash=sha256:")[1],
            "wheel_sha256_from": rel(REQ)}


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


def external_variants():
    """E18.1: DDX3X, every row of ddx3x_splice.parquet; TP53, every row of
    tp53_splice_scored_12nt.parquet, the 288 E7 reads; the five further genes, the rows
    of <gene>_variants.parquet at 3 <= |offset| <= 50. The counts must be the
    instruction's: 1,857; 288; 4,522, of which 4,075 for CTCF, SFPQ, RAD51D and XRCC2
    and 447 for SBDS."""
    import pandas as pd
    cols = ["gene", "variant_id", "chrom", "pos", "ref", "alt"]
    parts = [("DDX3X", pd.read_parquet(EXT / "ddx3x_splice.parquet")),
             ("TP53", pd.read_parquet(EXT / "tp53_splice_scored_12nt.parquet"))]
    for g in NEW_GENES:
        v = pd.read_parquet(EXT / f"{g.lower()}_variants.parquet")
        parts.append((g, v[v["offset"].abs().between(3, 50)]))
    for g, d in parts:
        if set(d["gene"]) != {g}:
            raise SystemExit(f"the {g} rows name genes {set(d['gene'])}")
    df = pd.concat([d[cols] for _, d in parts], ignore_index=True)
    df["chrom"], df["pos"] = df["chrom"].astype(str), df["pos"].astype(int)
    n = df["gene"].value_counts().to_dict()
    if (n.get("DDX3X"), n.get("TP53"), sum(n.get(g, 0) for g in NEW_GENES),
            sum(n.get(g, 0) for g in NEW_GENES[:4]), n.get("SBDS")) != (1857, 288, 4522, 4075, 447):
        raise SystemExit(f"the external variants are not the instruction's: {n}")
    if df["variant_id"].duplicated().any():
        raise SystemExit("a variant_id comes twice among the external variants")
    return df


def tagged_utc() -> str:
    """The tag's time, as the freeze guard recorded it from the tag object."""
    return json.loads(FROZEN.read_text())["tagged_utc"]


def environment(nb: Notebook, log, name: str = "ENVIRONMENT.json") -> None:
    import alphagenome
    freeze = sorted(subprocess.run([sys.executable, "-m", "pip", "freeze"], check=True,
                                   capture_output=True, text=True).stdout.split())
    env = {
        "written_by": rel(Path(__file__).resolve()),
        "written_utc": utc_now(),
        "python": sys.version,
        "platform": platform.platform(),
        "pip_freeze": freeze,
        "requirements": rel(REQ),
        "requirements_sha256": sha256_of(REQ),
        "client": {"package": "alphagenome", "version": alphagenome.__version__},
        "model": {"model_version_requested": None,
                  "model_version_requested_note": "dna_client.create(...) as in the "
                  "notebook's cell 7, with no model_version: the service's default",
                  "default_when_not_requested": model_default(),
                  "model_version_in_response": "not_provided",
                  "response_metadata_keys": nb.uns_keys},
    } | key_report()
    (DIR / name).write_text(json.dumps(env, indent=1, ensure_ascii=False) + "\n")
    log.output(DIR / name)


def main() -> None:
    ap = argparse.ArgumentParser()
    mode = ap.add_mutually_exclusive_group(required=True)
    for m in ("examples", "pilot", "repeat", "all", "drift",
              "drift-2b-start", "external", "drift-2b-end"):
        mode.add_argument(f"--{m}", dest="mode", action="store_const", const=m.replace("-", "_"))
    ap.add_argument("--workers", type=int, default=2)
    args = ap.parse_args()
    import pandas as pd
    e18 = args.mode in E18_MODES
    with RunLog(__file__) as log:
        log.note("key", key_report())
        DIR.mkdir(parents=True, exist_ok=True)
        CACHE.mkdir(parents=True, exist_ok=True)
        cache = CACHE / f"{args.mode}.jsonl"
        after, tag_utc = None, (tagged_utc() if e18 else None)
        if e18:
            log.input(FROZEN)
        # what each run must come after, read from the record of the run before it
        before = {"drift": (FULL, "--all"), "external": ("drift_2b_start_scores.tsv", "--drift-2b-start"),
                  "drift_2b_end": (EXTERNAL, "--external")}.get(args.mode)
        if before:
            prev = DIR / f"{before[0]}.provenance.json"
            if not prev.exists():
                raise SystemExit(f"--{args.mode.replace('_', '-')} runs after {before[1]}, and "
                                 f"{rel(prev)} is not there")
            log.input(prev)
            p = json.loads(prev.read_text())
            after = {"provenance": rel(prev), "sha256": sha256_of(prev), "end_utc": p["end_utc"]}
            if args.mode == "external":
                d = p["drift_check"]
                if d["identical"] != d["variants"] or d["variants"] != REPEAT_N:
                    raise SystemExit("the start drift check did not find its variants identical: "
                                     "no external gene is scored (E18.1)")
                after["drift_check"] = f"{d['identical']}/{d['variants']} identical"
        if args.mode.startswith("drift") and cache.exists() and cache.read_text().strip():
            raise SystemExit(f"{rel(cache)} is not empty: a drift check scores its "
                             f"{REPEAT_N} variants afresh and does not resume")
        nb = Notebook()
        log.input(SOURCE)
        key = ["variant_id", "gene", "stratum"]
        if args.mode == "examples":
            g = dict(nb.ns)
            src = {c["cell"]: c["source"] for c in json.loads(SOURCE.read_text())["quoted_cells"]}
            exec(src[15].split("\n\n# Score each variant")[0], g)      # the variants, as written
            todo = pd.DataFrame([{"variant_id": str(v), "gene": "", "stratum": "",
                                  "chrom": v.chromosome.removeprefix("chr"), "pos": v.position,
                                  "ref": v.reference_bases, "alt": v.alternate_bases}
                                 for v in g["variants"]])
        elif args.mode == "external":
            todo = external_variants()
            for f in ["ddx3x_splice.parquet", "tp53_splice_scored_12nt.parquet"] + [
                    f"{g.lower()}_variants.parquet" for g in NEW_GENES]:
                log.input(EXT / f)
            key = ["gene", "variant_id", "chrom", "pos", "ref", "alt"]
        else:
            df = pd.read_parquet(SET, columns=["variant_id", "gene", "stratum", "chrom",
                                                "pos", "ref", "alt"])
            log.input(SET)
            if set(df.gene) - set(TRAINING):
                raise SystemExit(f"genes outside the seven: {set(df.gene) - set(TRAINING)}")
            pilot, repeat, share = pilot_sample(df)
            log.note("pilot_share", share)
            ids = {"pilot": pilot, "repeat": repeat, "drift": repeat, "drift_2b_start": repeat,
                   "drift_2b_end": repeat, "all": df.variant_id.tolist()}[args.mode]
            todo = df.set_index("variant_id").loc[ids].reset_index()
        done = {}
        if cache.exists():
            for line in cache.read_text().splitlines():
                r = json.loads(line)
                done[r["variant_id"]] = r
        todo_rows = [r for r in todo.itertuples() if r.variant_id not in done]
        print(f"[{args.mode}] {len(todo)} variants, {len(done)} cached, "
              f"{len(todo_rows)} to score with {args.workers} worker(s)", flush=True)
        nb.uns_keys = None
        lock = threading.Lock()
        start = utc_now()
        t0 = time.time()

        def score(r):
            return nb.score(nb.variant(r.chrom, r.pos, r.ref, r.alt))

        def record(r, out):
            rec = {"variant_id": r.variant_id, "scored_utc": utc_now()} | out
            with lock:
                with open(cache, "a") as f:
                    f.write(json.dumps(rec) + "\n")
                done[r.variant_id] = rec
                if len(done) % 20 == 0:
                    print(f"  {len(done)} scored, {time.time() - t0:.0f} s", flush=True)

        try:                                     # the response's own metadata, once
            probe = nb.ns["score_one"](nb.variant(*todo.iloc[0][["chrom", "pos", "ref", "alt"]]))
        except Exception as e:                                  # noqa: BLE001
            raise SystemExit("the first call failed: stopping (E17.0 rule 4, if it is an "
                             f"authentication error): {_clean(f'{type(e).__name__}: {e}')}"
                             ) from None
        nb.uns_keys = sorted({k for ad in probe for k in ad.uns.keys()})
        try:
            stats = run(todo_rows, score, record, workers=args.workers)
        except Stopped as s:
            log.note("stopped", {"why": str(s), "calls": s.state["calls"],
                                 "retries": s.state["retries"],
                                 "failures": s.state["failures"], "in_cache": len(done)})
            raise
        elapsed = time.time() - t0
        if stats["failures"]:
            raise SystemExit(f"{len(stats['failures'])} variants failed after "
                             f"{MAX_ATTEMPTS} attempts: {stats['failures'][:3]}")
        rows = [done[v] for v in todo.variant_id]
        table = todo[key].merge(pd.DataFrame(rows), on="variant_id")
        cols = (key + COMPONENTS + [COMPOSITE]
                + [f"{p}_{c}" for c in COMPONENTS for p in ("rows", "nan", "min_raw")]
                + ["scored_utc"])
        out = DIR / f"{args.mode}_scores.tsv"
        table[cols].to_csv(out, sep="\t", index=False,          # exact, and read back as such
                           float_format=lambda x: repr(float(x)))
        log.output(out)
        environment(nb, log, "ENVIRONMENT_2b.json" if e18 else "ENVIRONMENT.json")
        checks, stops = {}, []
        summary = {c: {"variants_with_no_row": int((table[f"rows_{c}"] == 0).sum()),
                       "nan_raw_scores": int(table[f"nan_{c}"].sum()),
                       "min_raw": float(table[f"min_raw_{c}"].min())} for c in COMPONENTS}
        if [c for c in COMPONENTS if summary[c]["min_raw"] < 0]:
            stops.append("a raw score below zero")
        if args.mode in ("all", "drift", "drift_2b_start", "drift_2b_end"):
            ref = DIR / "pilot_scores.tsv"
            ids = pilot if args.mode == "all" else repeat
            diff = differing(out, ref, ids)
            checks["pilot_cross_check" if args.mode == "all" else "drift_check"] = {
                "against": rel(ref), "against_sha256": sha256_of(ref),
                "columns": COMPONENTS + [COMPOSITE], "variants": len(ids),
                "identical": len(ids) - len(diff), "differing": diff}
            if diff:
                stops.append(f"{len(diff)} of {len(ids)} variants differ from the pilot")
        if args.mode in ("drift", "drift_2b_end"):
            end = datetime.fromisoformat(after["end_utc"])
            early = [r["variant_id"] for r in rows if datetime.fromisoformat(r["scored_utc"]) <= end]
            checks["after_full_run" if args.mode == "drift" else "after_external"] = (
                after | {"every_scored_utc_later": not early})
            if early or len(todo_rows) != REPEAT_N:
                stops.append("the drift check did not score each of its variants afresh "
                             f"after {before[1]}")
        if args.mode == "drift_2b_start" and len(todo_rows) != REPEAT_N:
            stops.append("the drift check did not score each of its variants afresh")
        if args.mode == "external":
            checks["after_drift_2b_start"] = after
        if e18:
            when = datetime.fromisoformat(tag_utc.replace("Z", "+00:00"))
            early = [r["variant_id"] for r in rows if datetime.fromisoformat(r["scored_utc"]) <= when]
            checks["after_tag"] = {"tag": "frozen-thresholds-v1", "tagged_utc": tag_utc,
                                   "record": rel(FROZEN), "every_scored_utc_later": not early}
            if early:
                stops.append("a score is not later than the tag")
        extra = {}
        if e18:
            extra["code_sha256"] = {rel(p): sha256_of(p) for p in (
                Path(__file__).resolve(), Path(__file__).resolve().parent / "_runlog.py")}
        if args.mode == "external":
            extra["n_by_gene"] = {g: int(n) for g, n in table["gene"].value_counts().items()}
            extra["selection"] = ("DDX3X: every row of ddx3x_splice.parquet; TP53: every row of "
                                  "tp53_splice_scored_12nt.parquet; CTCF, SFPQ, RAD51D, XRCC2, SBDS: "
                                  "the rows of <gene>_variants.parquet at 3 <= |offset| <= 50")
            extra["sources"] = {rel(EXT / f): sha256_of(EXT / f) for f in (
                ["ddx3x_splice.parquet", "tp53_splice_scored_12nt.parquet"]
                + [f"{g.lower()}_variants.parquet" for g in NEW_GENES])}
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
            "model_version_default": model_default(),
            "start_utc": start, "end_utc": utc_now(),
            "seconds_scoring": round(elapsed, 1),
            "n_variants": len(table), "n_scored_this_run": len(todo_rows),
            "n_from_cache": len(todo) - len(todo_rows),
            "scored_utc_range": [min(table.scored_utc), max(table.scored_utc)],
            "calls": stats["calls"], "retries": stats["retries"],
            "failures": stats["failures"], "workers": args.workers,
            "sample": ({"seed": C.RANDOM_SEED, "pilot_n": PILOT_N, "repeat_n": REPEAT_N,
                        "share": log.notes.get("pilot_share")}
                       if args.mode in ("pilot", "repeat", "drift", "drift_2b_start",
                                        "drift_2b_end") else None),
            "genes": sorted(set(table.gene) - {""}),
            "summary": summary,
        } | extra | checks | {"output": rel(out), "output_sha256": sha256_of(out)} | key_report()
        (DIR / f"{args.mode}_scores.tsv.provenance.json").write_text(
            json.dumps(prov, indent=1, ensure_ascii=False) + "\n")
        log.output(DIR / f"{args.mode}_scores.tsv.provenance.json")
        log.note("seconds_scoring", round(elapsed, 1))
        log.note("retries", stats["retries"])
        for k, v in checks.items():
            log.note(k, v)
        print(f"[{args.mode}] wrote {rel(out)}: {len(table)} variants, "
              f"{elapsed:.0f} s for {len(todo_rows)}, {stats['retries']} retries")
        for k, v in checks.items():
            print(f"  {k}: " + (f"{v['identical']}/{v['variants']} identical"
                                if "identical" in v else json.dumps(v)))
        if stops:
            raise SystemExit("; ".join(stops) + ": stopping (E17.10; E18.8)")


if __name__ == "__main__":
    main()
