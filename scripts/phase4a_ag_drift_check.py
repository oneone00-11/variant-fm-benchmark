"""
phase4a_ag_drift_check.py -- AlphaGenome drift gate BEFORE scoring TP53.

Why: the 7-gene AlphaGenome scores were computed 2026-06-19 with client
`alphagenome==0.6.1`. The served model is UNVERSIONED (the API returns no model
version/timestamp; ag_cache.tsv stores only chrom/pos/ref/alt/score), and client
0.7.0 shipped 2026-06-26 (7 days after that run). So the only way to know TP53
would be scored by the SAME model as the 7 genes is to re-score a few 7-gene
variants today and check they still reproduce the cached values.

This script re-scores 14 known 7-gene variants (4 high-signal splice-core anchors
+ 10 sampled across the cached-score range) and compares to `ag_cache.tsv`.

  ALL within tolerance  -> NO DRIFT: safe to score TP53 by the same API path.
  ANY out of tolerance  -> DRIFT:    STOP; TP53 and the 7 genes must be
                                     re-scored together for comparability.

The scoring code path (score_one) is COPIED VERBATIM from 91_score_alphagenome.py
so the check exercises the identical logic (SPLICE scorers, 1 MB interval,
max|raw_score| via tidy_scores, chr-prefixed coords, 5x backoff). Do NOT edit it
here without editing it there.

Run (in the same env as the 7-gene run; key from env only, never printed):
    pip install 'alphagenome==0.6.1'
    export ALPHAGENOME_API_KEY=...      # your key; this script never logs it
    python scripts/phase4a_ag_drift_check.py
Paste the printed VERDICT + table back.
"""
from __future__ import annotations
import csv
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "data" / "raw" / "scores" / "alphagenome" / "ag_cache.tsv"
OUT = ROOT / "data" / "external" / "ag_drift_check_result.tsv"

# --- tolerance -------------------------------------------------------------
# AlphaGenome splice raw_score is a float16-quantized model output (cached values
# fall on ~1/256 grid), so an UNCHANGED, deterministic model reproduces to ~1e-3.
# PASS if within EITHER 0.01 absolute or 1% relative; anything larger = drift.
TOL_ABS = 0.01
TOL_REL = 0.01

# --- the 4 high-signal anchor variants (pathogenic splice-core; milestone7) ---
# (gene, hgvs, chrom, pos, ref, alt, cached_alphagenome_splice)
ANCHORS = [
    ("BRCA1", "NM_007294.3:c.5194-1G>C", "17", 43057136, "C", "G", 7.378906),
    ("BRCA1", "NM_007294.3:c.4986+1G>A", "17", 43070927, "C", "T", 3.832031),
    ("VHL",   "ENST00000256474.3:c.340+1G>T", "3", 10142188, "G", "T", 5.324219),
    ("VHL",   "ENST00000256474.3:c.341-2A>T", "3", 10146512, "A", "T", 7.789062),
]

N_SAMPLE = 10   # additional variants sampled deterministically across the score range


def require_env_and_pin():
    if not os.environ.get("ALPHAGENOME_API_KEY"):
        sys.exit("ALPHAGENOME_API_KEY not set - aborting (key must come from env).")
    try:
        import alphagenome
    except ImportError:
        sys.exit("alphagenome not installed - run: pip install 'alphagenome==0.6.1'")
    ver = getattr(alphagenome, "__version__", "?")
    if ver != "0.6.1":
        sys.exit(f"alphagenome=={ver} installed, but the 7-gene run used 0.6.1. "
                 f"Pin it: pip install 'alphagenome==0.6.1'  (drift check must match the run's client).")
    print(f"[drift] alphagenome client {ver} (pinned to the 7-gene run) OK", flush=True)


def sample_from_cache(exclude_keys):
    """Deterministically pick N_SAMPLE cached variants spanning min->max cached
    score (evenly spaced by rank), skipping any coord already in `exclude_keys`."""
    rows = []
    with open(CACHE) as f:
        for d in csv.DictReader(f, delimiter="\t"):
            try:
                s = float(d["alphagenome_splice"])
            except (ValueError, KeyError):
                continue
            key = (d["chrom"], int(d["pos"]), d["ref"], d["alt"])
            if key in exclude_keys:
                continue
            rows.append((s, key))
    rows.sort(key=lambda x: (x[0], x[1]))          # by score, then coord -> deterministic
    n = len(rows)
    picks, seen = [], set()
    for k in range(N_SAMPLE):
        idx = round(k * (n - 1) / (N_SAMPLE - 1)) if N_SAMPLE > 1 else 0
        # nudge off duplicates
        while idx in seen and idx < n - 1:
            idx += 1
        seen.add(idx)
        s, key = rows[idx]
        picks.append((key[0], key[1], key[2], key[3], s))
    return picks


# ===========================================================================
# score_one -- COPIED VERBATIM from 91_score_alphagenome.py (identical logic)
# ===========================================================================
def build_scorer():
    from alphagenome.models import dna_client, variant_scorers
    model = dna_client.create(os.environ["ALPHAGENOME_API_KEY"])
    recs = variant_scorers.RECOMMENDED_VARIANT_SCORERS
    SPLICE = [recs[n] for n in recs if "SPLICE" in str(n).upper()]
    print(f"[drift] {len(SPLICE)} splice scorers loaded (same as the run)", flush=True)
    return model, SPLICE, dna_client, variant_scorers


def score_one(v, model, SPLICE, dna_client, variant_scorers, probe=None):
    from alphagenome.data import genome
    chrom, pos, ref, alt = v
    var = genome.Variant(chromosome="chr" + chrom, position=pos,
                         reference_bases=ref, alternate_bases=alt)
    iv = var.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
    for attempt in range(5):
        try:
            out = model.score_variant(interval=iv, variant=var, variant_scorers=SPLICE)
            if probe is not None and not probe:      # probe the response metadata once
                probe.update(_probe_meta(out))
            df = variant_scorers.tidy_scores(out)
            if df is None or len(df) == 0 or "raw_score" not in df.columns:
                return ""
            return f"{df['raw_score'].abs().max():.6f}"
        except Exception as e:
            time.sleep(2 * (attempt + 1))
            last = str(e)[:80]
    return None  # give up -> None (report as ERROR)


def _probe_meta(out):
    """Opportunistically capture any model version/timestamp the API exposes
    (we established the cache stored none; this checks the live response)."""
    info = {}
    try:
        items = out if isinstance(out, (list, tuple)) else [out]
        a = items[0]
        for attr in ("uns", "attrs", "metadata"):
            d = getattr(a, attr, None)
            if isinstance(d, dict):
                for k, val in d.items():
                    if any(t in str(k).lower() for t in ("version", "model", "date", "time", "build", "commit")):
                        info[f"{attr}.{k}"] = str(val)[:80]
    except Exception as e:
        info["_probe_error"] = str(e)[:80]
    return info or {"_note": "no version/model metadata found on API response"}


def verdict(cached, rescored):
    if rescored is None:
        return "ERROR", None, None
    if rescored == "":
        return ("PASS" if (cached is None) else "FAIL(empty)"), None, None
    r = float(rescored)
    ab = abs(r - cached)
    rel = ab / abs(cached) if cached else (0.0 if ab == 0 else 1.0)
    return ("PASS" if (ab <= TOL_ABS or rel <= TOL_REL) else "FAIL"), ab, rel


def main():
    require_env_and_pin()
    model, SPLICE, dna_client, variant_scorers = build_scorer()

    anchor_keys = {(c, p, r, a) for (_, _, c, p, r, a, _) in ANCHORS}
    checks = [(g, h, c, p, r, a, cv) for (g, h, c, p, r, a, cv) in ANCHORS]
    for (c, p, r, a, cv) in sample_from_cache(anchor_keys):
        checks.append(("(cache)", f"{c}:{p}:{r}:{a}", c, p, r, a, cv))

    print(f"[drift] re-scoring {len(checks)} variants "
          f"({len(ANCHORS)} anchors + {len(checks)-len(ANCHORS)} sampled)...\n", flush=True)

    probe, results, n_pass, n_fail, n_err = {}, [], 0, 0, 0
    hdr = f"{'gene':7} {'hgvs/key':30} {'chrom:pos:ref>alt':22} {'cached':>10} {'rescored':>10} {'absΔ':>8} {'relΔ':>7}  verdict"
    print(hdr); print("-" * len(hdr))
    for (g, h, c, p, r, a, cv) in checks:
        rs = score_one((c, p, r, a), model, SPLICE, dna_client, variant_scorers,
                       probe=probe)
        vd, ab, rel = verdict(cv, rs)
        n_pass += vd == "PASS"; n_fail += vd.startswith("FAIL"); n_err += vd == "ERROR"
        coord = f"{c}:{p}:{r}>{a}"
        rss = rs if rs not in (None, "") else (rs or "NA")
        abs_s = f"{ab:.4f}" if ab is not None else "-"
        rel_s = f"{rel*100:.2f}%" if rel is not None else "-"
        print(f"{g:7} {h[:30]:30} {coord:22} {cv:>10.4f} {str(rss):>10} {abs_s:>8} {rel_s:>7}  {vd}")
        results.append((g, h, coord, cv, rss, abs_s, rel_s, vd))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        f.write("gene\thgvs_or_key\tcoord\tcached\trescored\tabs_delta\trel_delta\tverdict\n")
        for row in results:
            f.write("\t".join(str(x) for x in row) + "\n")

    print("\n[drift] API response metadata probe:", probe)
    print(f"[drift] wrote {OUT}")
    overall = "NO DRIFT" if (n_fail == 0 and n_err == 0) else (
              "INCONCLUSIVE (errors)" if n_fail == 0 else "DRIFT DETECTED")
    print("\n" + "=" * 62)
    print(f"VERDICT: {overall}   (pass={n_pass}  fail={n_fail}  error={n_err}  of {len(checks)})")
    if overall == "NO DRIFT":
        print("-> Served model reproduces the 7-gene AlphaGenome scores.")
        print("-> SAFE to score TP53 via the API (9-feature set; NT excluded; GPN pinned to cf1718a9).")
    elif overall.startswith("DRIFT"):
        print("-> Served model has CHANGED since 2026-06-19.")
        print("-> STOP: do NOT score TP53 alone. TP53 AND the 7 genes must be re-scored together.")
    else:
        print("-> API errors prevented a clean check; re-run (network/quota), then re-read the verdict.")
    print("=" * 62)


if __name__ == "__main__":
    main()
