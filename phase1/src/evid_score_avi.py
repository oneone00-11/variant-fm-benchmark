"""E5 -- AlphaGenome Variant Impact (AVI) for the analysis set.

AVI is the single combined coding/non-coding score released with AlphaGenome Atlas
on 8 September 2026. It reaches the analysis set through the Atlas API, which is a
different service from the per-variant scoring API the `alphagenome` column came
from, and which appeared in client 0.9.0 (`alphagenome.atlas`, CHANGELOG 0.9.0:
"Add AlphaGenome Atlas API for programmatic access to Atlas scores").

Two prerequisites, neither of which this machine has:

  * `alphagenome >= 0.9.0`. The pinned environment used for the atlas columns is
    0.7.0, which has no `alphagenome.atlas` module at all.
  * `ALPHAGENOME_API_KEY`. The key is issued per person, free for non-commercial
    use, at https://deepmind.google.com/science/alphagenome. It is not on this
    machine and cannot be obtained without the account holder.

So this module is written to run, not run yet. It refuses loudly rather than
writing a half-scored column, and it prints the scorer names the service actually
offers before picking one, because the AVI scorer's key is served by
`scorer_metadata()` at runtime and is not a constant in the client.

Licence, to be recorded in docs/column-provenance.md and LICENSE-DATA once a
column exists: outputs of the AlphaGenome API and information in AlphaGenome Atlas
are for NON-COMMERCIAL use only and must not be used to train other models, except
for those artefacts the AlphaGenome Terms of Service designate as Permissive Use
Downloadable Artifacts. That is more restrictive than the CC BY 4.0 the frozen data
carries, so AVI must be declared per column, as the companion atlas already does
for its eight restricted score columns.

Run:  python -m src.evid_score_avi --check      # environment only, no network
      python -m src.evid_score_avi --list       # what the service offers
      python -m src.evid_score_avi --run
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import time

import numpy as np
import pandas as pd

SET_PATH = Path("data/evidence/analysis_set_v1.parquet")
OUT = Path("data/evidence/avi_scores.parquet")
CACHE = Path("data/evidence/avi_cache.jsonl")

MIN_CLIENT = (0, 9, 0)
AVI_HINTS = ("avi", "variant_impact", "variant impact")


KEY_FILE = Path.home() / ".config" / "alphagenome" / "key"


def read_key() -> str | None:
    """The key, from the environment or from a file OUTSIDE the repository.

    Never from the tree: tests/test_no_secrets.py greps git ls-files for
    provider-anchored patterns and fails on a match, so a key committed here would
    break the build before it could be pushed.
    """
    env = os.environ.get("ALPHAGENOME_API_KEY")
    if env:
        return env.strip()
    if KEY_FILE.exists():
        return KEY_FILE.read_text().strip() or None
    return None


def environment_report() -> dict:
    key = read_key()
    rep = {"api_key_present": bool(key),
           "api_key_source": ("ALPHAGENOME_API_KEY" if os.environ.get("ALPHAGENOME_API_KEY")
                              else (str(KEY_FILE) if key else None))}
    try:
        import alphagenome
        rep["client_version"] = getattr(alphagenome, "__version__", "unknown")
    except ImportError:
        rep["client_version"] = None
    try:
        from alphagenome.atlas import atlas  # noqa: F401
        rep["atlas_module"] = True
    except Exception:
        rep["atlas_module"] = False
    v = rep.get("client_version")
    try:
        rep["client_new_enough"] = (
            tuple(int(x) for x in str(v).split(".")[:3]) >= MIN_CLIENT)
    except Exception:
        rep["client_new_enough"] = False
    rep["ready"] = bool(rep["api_key_present"] and rep["atlas_module"]
                        and rep["client_new_enough"])
    return rep


def _blocked(rep: dict) -> int:
    print("[E5] BLOCKED -- AVI cannot be scored here.")
    print(json.dumps(rep, indent=2))
    print("\nTo unblock, in this order:")
    print("  1. obtain a free non-commercial API key at")
    print("     https://deepmind.google.com/science/alphagenome")
    print("     and export ALPHAGENOME_API_KEY")
    print("  2. create a separate environment with alphagenome>=0.9.0 -- do NOT")
    print("     upgrade the 0.7.0 environment the published alphagenome column was")
    print("     scored in, or that column stops being reproducible")
    print("  3. re-run: python -m src.evid_score_avi --list, then --run")
    return 2


def _client():
    from alphagenome.atlas import atlas
    key = read_key()
    if not key:
        raise SystemExit("[E2.1] no API key: set ALPHAGENOME_API_KEY or write "
                         f"{KEY_FILE}")
    return atlas.create(key)


def _describe(meta) -> dict:
    """Whatever the service says about a scorer, as plain data."""
    out = {}
    for attr in ("name", "description", "version", "output_type", "organism",
                 "scorer_type", "requires_gene", "ontology_terms", "num_tracks"):
        v = getattr(meta, attr, None)
        if v is None:
            continue
        try:
            out[attr] = v if isinstance(v, (str, int, float, bool)) else str(v)
        except Exception:                                    # noqa: BLE001
            pass
    if not out:
        out["repr"] = str(meta)[:400]
    return out


def list_scorers() -> int:
    """List every scorer the service offers and record it verbatim.

    The AVI scorer's key is served at runtime, not a constant in the client, so the
    list is what identifies it -- and the list is what a reader needs in order to
    check that the right one was chosen.
    """
    rep = environment_report()
    if not rep["ready"]:
        return _blocked(rep)
    md = _client().scorer_metadata()
    record = {
        "listed_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "client_version": rep["client_version"],
        "n_scorers": len(md),
        "scorers": {name: _describe(meta) for name, meta in sorted(md.items())},
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d")
    path = OUT.parent / f"avi_scorers_{stamp}.json"
    path.write_text(json.dumps(record, indent=2, default=str) + "\n")
    print(f"[E2.1] the Atlas service offers {len(md)} scorers; written to {path}")
    for name in sorted(md):
        mark = "  <-- AVI candidate" if any(h in name.lower() for h in AVI_HINTS) else ""
        print(f"  {name}{mark}")
    return 0


SPLICE_SCORERS = ["SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS"]
AVI_SCORER = "AVI_SCORE"

TARGETS = {
    "analysis": Path("data/evidence/analysis_set_v1.parquet"),
    "ddx3x": Path("data/evidence/external/ddx3x_splice.parquet"),
    "tp53": Path("data/external/tp53_splice_scored_v2.parquet"),
}


def _aggregate(ad, variants_order: list[str]) -> np.ndarray:
    """One value per variant: the maximum over tracks, and over rows where a
    scorer returns several (SPLICE_JUNCTIONS emits one row per junction).

    The maximum is the atlas's own aggregation over splice tracks, so the columns
    are comparable in construction even though the score type is not the same one
    the v0.7.0 column was built on.
    """
    import pandas as _pd
    x = np.asarray(ad.X, dtype=float)
    per_row = np.nanmax(x, axis=1) if x.ndim == 2 and x.shape[1] > 1 else x.ravel()
    key = ad.obs["variant"].astype(str).to_numpy()
    s = _pd.Series(per_row).groupby(key).max()
    return s.reindex(variants_order).to_numpy(dtype=float)


MAX_RETRIES = 6
BACKOFF_BASE = 20.0          # seconds; the quota is per minute


def run(scorer: str | None, targets: list[str], with_splice: bool,
        batch: int, out_path: Path, pause: float = 0.0) -> int:
    rep = environment_report()
    if not rep["ready"]:
        return _blocked(rep)
    from alphagenome.data import genome

    client = _client()
    md = client.scorer_metadata()
    if scorer is None:
        scorer = AVI_SCORER if AVI_SCORER in md else None
        if scorer is None:
            hits = [n for n in md if any(h in n.lower() for h in AVI_HINTS)]
            print(f"[E2.1] cannot pick the AVI scorer automatically; candidates: {hits}")
            return 3
    wanted = [scorer] + ([s for s in SPLICE_SCORERS if s in md] if with_splice else [])
    colname = {scorer: "avi"}
    for s in SPLICE_SCORERS:
        colname[s] = "avi_splice_" + s.replace("SPLICE_", "").lower()
    print(f"[E2.1] scorers: {wanted}")

    frames = []
    for name in targets:
        src = TARGETS[name]
        if not src.exists():
            print(f"[E2.1] {name}: {src} missing, skipped")
            continue
        df = pd.read_parquet(src, columns=["variant_id", "chrom", "pos", "ref", "alt"])
        df = df.drop_duplicates("variant_id").reset_index(drop=True)
        df["set"] = name
        frames.append(df)
    todo_all = pd.concat(frames, ignore_index=True)
    cols = [colname[s] for s in wanted]

    # Resume from the checkpoint if one survived, otherwise from a previous output:
    # a run that hit the service's per-minute quota still wrote what it had, and
    # re-querying thousands of variants that are already answered wastes the quota
    # the remaining ones need.
    ckpt = Path(str(out_path) + ".partial.parquet")
    done = {}
    for src in (ckpt, out_path):
        if not src.exists():
            continue
        prev = pd.read_parquet(src)
        if not all(c in prev.columns for c in cols):
            continue
        for r in prev[["set", "variant_id"] + cols].itertuples(index=False):
            if not any(pd.isna(x) for x in r[2:]):
                done[(r[0], r[1])] = tuple(r[2:])
        break
    if done:
        print(f"[resume] {len(done):,} of {len(todo_all):,} already scored")

    failures = []
    t0 = time.time()
    pending = todo_all[[(s, v) not in done
                        for s, v in zip(todo_all["set"], todo_all["variant_id"])]]
    n_requests = 0
    for start in range(0, len(pending), batch):
        chunk = pending.iloc[start:start + batch]
        vs = [genome.Variant(chromosome=f"chr{r.chrom}", position=int(r.pos),
                             reference_bases=r.ref, alternate_bases=r.alt)
              for r in chunk.itertuples()]
        order = [f"chr{r.chrom}:{int(r.pos)}:{r.ref}>{r.alt}" for r in chunk.itertuples()]
        # The service enforces a per-minute request quota and answers a burst past
        # it with RESOURCE_EXHAUSTED. Back off and retry rather than dropping the
        # batch: a dropped batch is a hole in a column, and a hole that only shows
        # up as NaN later is worse than waiting.
        res = None
        for attempt in range(MAX_RETRIES):
            try:
                res = client.query_variants(vs, requested_scorers=wanted,
                                            progress_bar=False)
                n_requests += len(vs)
                break
            except Exception as e:                            # noqa: BLE001
                msg = str(e)
                exhausted = "RESOURCE_EXHAUSTED" in msg or "Quota exceeded" in msg
                if attempt == MAX_RETRIES - 1:
                    failures.append({"first_variant_id": chunk["variant_id"].iloc[0],
                                     "n": len(chunk), "attempts": MAX_RETRIES,
                                     "error": msg[:300]})
                    print(f"  batch at {start} FAILED after {MAX_RETRIES}: "
                          f"{msg[:100]}", flush=True)
                    break
                wait = BACKOFF_BASE * (2 ** attempt) if exhausted else 2.0
                print(f"  batch at {start} {'quota' if exhausted else 'error'}, "
                      f"retry {attempt + 1}/{MAX_RETRIES - 1} in {wait:.0f}s",
                      flush=True)
                time.sleep(wait)
        if res is None:
            continue
        vals = {s: _aggregate(res[s], order) for s in wanted}
        for k, r in enumerate(chunk.itertuples()):
            done[(r.set, r.variant_id)] = tuple(float(vals[s][k]) for s in wanted)
        frame = todo_all.copy()
        arr = np.asarray([done.get((s, v), (np.nan,) * len(cols))
                          for s, v in zip(frame["set"], frame["variant_id"])], dtype=float)
        for c_i, c in enumerate(cols):
            frame[c] = arr[:, c_i]
        keep = [(s, v) in done for s, v in zip(frame["set"], frame["variant_id"])]
        frame[keep].to_parquet(ckpt, index=False)
        el = time.time() - t0
        n_left = len(pending) - min(start + batch, len(pending))
        print(f"  {len(done):>6}/{len(todo_all)}  {el/60:5.1f} min, "
              f"~{n_left * el / max(start + batch, 1) / 60:5.1f} min left", flush=True)
        if n_left and pause:
            time.sleep(pause)

    frame = todo_all.copy()
    arr = np.asarray([done.get((s, v), (np.nan,) * len(cols))
                      for s, v in zip(frame["set"], frame["variant_id"])], dtype=float)
    for c_i, c in enumerate(cols):
        frame[c] = arr[:, c_i]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(out_path, index=False)
    if not failures and not frame[cols].isna().any().any():
        ckpt.unlink(missing_ok=True)
    else:
        print(f"[E2.1] {int(frame[cols[0]].isna().sum())} variants unscored; "
              f"checkpoint kept at {ckpt.name} -- re-run to continue")

    prov = {
        "columns": cols,
        "product": "AlphaGenome Atlas -- AlphaGenome Variant Impact (AVI)",
        "atlas_released": "2026-09-08",
        "client_version": rep["client_version"],
        "scorers_requested": wanted,
        "scorer_list_recorded_in": sorted(
            str(p) for p in OUT.parent.glob("avi_scorers_*.json")),
        "aggregation": "per variant, the maximum over tracks and over rows; "
                       "SPLICE_JUNCTIONS returns one row per junction",
        "accessed_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "n_variants": int(len(frame)),
        "n_scored": {c: int(frame[c].notna().sum()) for c in cols},
        "by_set": {str(name): {"n": int(len(g)),
                               **{c: int(g[c].notna().sum()) for c in cols}}
                   for name, g in frame.groupby("set")},
        "n_requests": n_requests,
        "failed_batches": failures,
        "runtime_minutes": round((time.time() - t0) / 60, 1),
        "licence": "AlphaGenome Terms of Service -- NON-COMMERCIAL use only; outputs "
                   "must not be used to train other models",
        "licence_source": "https://deepmind.google.com/science/alphagenome/terms",
        "output_sha256": hashlib.sha256(out_path.read_bytes()).hexdigest(),
    }
    Path(str(out_path) + ".provenance.json").write_text(
        json.dumps(prov, indent=2, default=str) + "\n")
    print(json.dumps({k: v for k, v in prov.items() if k != "scorer_list_recorded_in"},
                     indent=2, default=str))
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--scorer", default=None)
    ap.add_argument("--targets", default="analysis,ddx3x,tp53")
    ap.add_argument("--no-splice", action="store_true",
                    help="score only the AVI score, not the splice sub-scores")
    ap.add_argument("--batch", type=int, default=250)
    ap.add_argument("--pause", type=float, default=0.0,
                    help="seconds between batches, to stay under the request quota")
    ap.add_argument("--out", default="data/evidence/avi.parquet")
    args = ap.parse_args()
    if args.check:
        rep = environment_report()
        print(json.dumps(rep, indent=2))
        return 0 if rep["ready"] else _blocked(rep)
    if args.list:
        return list_scorers()
    if args.run:
        return run(args.scorer, [x for x in args.targets.split(",") if x],
                   not args.no_splice, args.batch, Path(args.out), args.pause)
    ap.error("pass --check, --list or --run")


if __name__ == "__main__":
    sys.exit(main())
