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

import numpy as np
import pandas as pd

SET_PATH = Path("data/evidence/analysis_set_v1.parquet")
OUT = Path("data/evidence/avi_scores.parquet")
CACHE = Path("data/evidence/avi_cache.jsonl")

MIN_CLIENT = (0, 9, 0)
AVI_HINTS = ("avi", "variant_impact", "variant impact")


def environment_report() -> dict:
    rep = {"api_key_present": bool(os.environ.get("ALPHAGENOME_API_KEY"))}
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
    return atlas.create(os.environ["ALPHAGENOME_API_KEY"])


def list_scorers() -> int:
    rep = environment_report()
    if not rep["ready"]:
        return _blocked(rep)
    md = _client().scorer_metadata()
    print(f"[E5] the Atlas service offers {len(md)} scorers:")
    for name in sorted(md):
        mark = "  <-- AVI?" if any(h in name.lower() for h in AVI_HINTS) else ""
        print(f"  {name}{mark}")
    return 0


def run(scorer: str | None) -> int:
    rep = environment_report()
    if not rep["ready"]:
        return _blocked(rep)
    from alphagenome.data import genome
    from alphagenome.atlas import atlas  # noqa: F401

    client = _client()
    md = client.scorer_metadata()
    if scorer is None:
        hits = [n for n in md if any(h in n.lower() for h in AVI_HINTS)]
        if len(hits) != 1:
            print(f"[E5] cannot pick the AVI scorer automatically; candidates: {hits}")
            print("     re-run with --scorer <name> after checking --list")
            return 3
        scorer = hits[0]
    print(f"[E5] scorer: {scorer}")

    df = pd.read_parquet(SET_PATH, columns=["variant_id", "chrom", "pos", "ref", "alt"])
    variants = [genome.Variant(chromosome=f"chr{r.chrom}", position=int(r.pos),
                               reference_bases=r.ref, alternate_bases=r.alt)
                for r in df.itertuples()]
    scores = client.query_variants(variants, requested_scorers=[scorer])

    ad = scores[scorer]
    values = np.asarray(ad.X).reshape(len(variants), -1)
    out = df[["variant_id"]].copy()
    out["avi"] = values.max(axis=1) if values.shape[1] > 1 else values.ravel()
    if values.shape[1] > 1:
        out["avi_n_tracks"] = values.shape[1]

    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(OUT, index=False)
    prov = {
        "column": "avi",
        "product": "AlphaGenome Atlas -- AlphaGenome Variant Impact (AVI)",
        "released": "2026-09-08",
        "client_version": rep["client_version"],
        "scorer": scorer,
        "accessed_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "n_variants": int(len(out)),
        "n_scored": int(out["avi"].notna().sum()),
        "licence": "AlphaGenome Terms of Service -- non-commercial use only; "
                   "outputs must not be used to train other models",
        "licence_source": "https://deepmind.google.com/science/alphagenome/terms",
        "output_sha256": hashlib.sha256(OUT.read_bytes()).hexdigest(),
    }
    Path(str(OUT) + ".provenance.json").write_text(json.dumps(prov, indent=2) + "\n")
    print(json.dumps(prov, indent=2))
    print("\n[E5] now re-run E1 (to merge the column), then E2, E3 and E4.")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--scorer", default=None)
    args = ap.parse_args()
    if args.check:
        rep = environment_report()
        print(json.dumps(rep, indent=2))
        return 0 if rep["ready"] else _blocked(rep)
    if args.list:
        return list_scorers()
    if args.run:
        return run(args.scorer)
    ap.error("pass --check, --list or --run")


if __name__ == "__main__":
    sys.exit(main())
