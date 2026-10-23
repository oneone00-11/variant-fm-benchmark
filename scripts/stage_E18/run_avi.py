"""E18.2: the four Atlas columns for TP53 and the five further genes, by the analysis
set's own scorer.

The analysis set's avi, avi_splice_sites, avi_splice_site_usage and
avi_splice_junctions come from phase1/src/evid_score_avi.py's run() (the AlphaGenome
Atlas, client 0.9.0, scorers AVI_SCORE, SPLICE_SITES, SPLICE_SITE_USAGE,
SPLICE_JUNCTIONS; the maximum over tracks and rows), in one run that also scored DDX3X
and TP53's first 192 rows (avi.parquet). TP53's other 96 rows came from a later run that
re-scored TP53 rows but no analysis-set variant, so TP53 is scored again here, all 288,
with the five further genes' 4,522 and the 39 anchors, in one run of the same run(),
unchanged, pointed at these targets as evid_tp53_extend points it at TP53, from phase1/
as it expects. DDX3X keeps its values from the analysis set's own run.

run() retries a failed batch and, when every attempt fails, goes on to the next one. The
client it gets is wrapped so that E18.0 rule 12 holds around it: an authentication
error stops the run at once (E17.0 rule 4), any other error backs off exponentially
before run() retries it, and a batch whose every attempt fails stops the run; run()'s
checkpoint keeps every batch already scored, so a rerun resumes. Error messages are
recorded with the key, should it ever appear in one, cut out; the key is read only by
evid_score_avi.read_key() and never written.

Outputs (phase1/data/e18_cache/, gitignored): avi_e18.parquet and the provenance run()
writes beside it, avi_e18.environment.json, and avi_e18.segments.json, the start and
end of every invocation. --trial scores the anchors and the first 50 of these external
variants into avi_e18_trial.parquet.

Run in the Atlas environment (alphagenome 0.9.0):
    ../envs/alphagenome-atlas/bin/python scripts/stage_E18/run_avi.py
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel, utc_now  # noqa: E402
import _e18 as E  # noqa: E402

PHASE1 = REPO / "phase1"
GENES = ["TP53", "CTCF", "SFPQ", "RAD51D", "XRCC2", "SBDS"]
COLUMNS = ["avi", "avi_splice_sites", "avi_splice_site_usage", "avi_splice_junctions"]
FATAL = ("UNAUTHENTICATED", "PERMISSION_DENIED", "API key", "api key")   # as E17's
BACKOFF = 5.0


class Stop(SystemExit):
    """Not an Exception, so run()'s retry loop lets it through."""


class Guarded:
    def __init__(self, client, A):
        self._c, self._A = client, A
        self.in_a_row, self.errors, self.scorers = 0, [], None

    def _clean(self, msg: str) -> str:
        key = self._A.read_key()
        if key and key in msg:
            msg = msg.replace(key, "<key>")
        del key
        return msg[:300]

    def scorer_metadata(self):
        md = self._c.scorer_metadata()
        self.scorers = sorted(md)
        return md

    def query_variants(self, *args, **kwargs):
        try:
            out = self._c.query_variants(*args, **kwargs)
        except Exception as e:                                  # noqa: BLE001
            msg = self._clean(f"{type(e).__name__}: {e}")
            if any(f in msg for f in FATAL):
                raise Stop(f"authentication error: stopping (E17.0 rule 4): {msg}")
            self.in_a_row += 1
            self.errors.append({"utc": utc_now(), "error": msg})
            if self.in_a_row >= self._A.MAX_RETRIES:
                raise Stop(f"a batch failed all {self._A.MAX_RETRIES} attempts: stopping "
                           "(E18.0 rule 12); the checkpoint keeps every batch scored")
            if "RESOURCE_EXHAUSTED" not in msg and "Quota exceeded" not in msg:
                time.sleep(BACKOFF * 2 ** (self.in_a_row - 1))  # run() adds its own 2 s
            raise
        self.in_a_row = 0
        return out

    def __getattr__(self, name):
        return getattr(self._c, name)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--trial", action="store_true",
                    help=f"the anchors and the first {E.TRIAL_N} external variants only")
    args = ap.parse_args()
    os.chdir(PHASE1)                       # run() reads data/evidence/ relative to phase1/
    sys.path.insert(0, str(PHASE1))
    from src import evid_score_avi as A
    out = E.cache_path("avi_e18", args.trial, ".parquet")
    segf = E.cache_path("avi_e18", args.trial, ".segments.json")
    with RunLog(__file__) as log:
        log.input(E.INPUT)
        df = E.load_input(args.trial, keep=lambda x: x["gene"].isin(GENES))
        E.CACHE.mkdir(parents=True, exist_ok=True)
        targets, names = {}, []
        for name, part in [("anchor", df[df["set"] == "anchor"])] + [
                (g.lower(), df[(df["set"] == "external") & (df["gene"] == g)]) for g in GENES]:
            if not len(part):
                continue
            p = E.cache_path(f"avi_target_{name}", args.trial, ".parquet")
            part[["variant_id", "chrom", "pos", "ref", "alt"]].to_parquet(p, index=False)
            targets[name] = p
            names.append(name)
            log.output(p)
        A.TARGETS = dict(A.TARGETS, **targets)
        real = A._client
        guard = {}

        def guarded():
            guard["g"] = Guarded(real(), A)
            return guard["g"]

        A._client = guarded
        segments = json.loads(segf.read_text()) if segf.exists() else []
        segments.append({"start_utc": utc_now()})
        E.write_json(segf, segments)
        clock = E.Clock()
        try:
            rc = A.run(None, names, True, 250, out)
        finally:
            g = guard.get("g")
            segments[-1].update({"end_utc": utc_now(),
                                 "errors": g.errors if g else [],
                                 "scorers_offered": g.scorers if g else None})
            E.write_json(segf, segments)
            log.output(segf)
            log.note("segment", segments[-1])
        if rc:
            raise SystemExit(f"run() returned {rc}")
        seconds = clock.seconds()
        log.output(out)
        log.output(Path(str(out) + ".provenance.json"))
        import pandas as pd
        got = pd.read_parquet(out)
        anchors = got[got["set"] == "anchor"]
        checks = {c: E.anchor_check(c, dict(zip(anchors["variant_id"], anchors[c]))) for c in COLUMNS}
        for c in COLUMNS:
            E.report(log, checks[c], None)
        log.note("anchors", checks)
        if args.trial:
            E.report(log, None, E.throughput(
                seconds, len(df), len(E.load_input(False, keep=lambda x: x["gene"].isin(GENES)))))
        freeze = subprocess.run([sys.executable, "-m", "pip", "freeze"], check=True,
                                capture_output=True, text=True).stdout.split()
        env = {"scorer": "phase1/src/evid_score_avi.py::run", "targets": names,
               "interpreter": sys.executable, "python": platform.python_version(),
               "platform": platform.platform(), "pip_freeze": freeze,
               "segments": segments} | A.environment_report()
        envf = E.cache_path("avi_e18", args.trial, ".environment.json")
        E.write_json(envf, env)
        log.output(envf)
        print(f"done: {rel(out)}")


if __name__ == "__main__":
    main()
