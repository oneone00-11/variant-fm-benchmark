"""E9.1: how often each source of the label audit's uncommitted files serves the file.

Each source scripts/reproduce_sources.json lists is asked for its file by the entry
point's own single attempt (scripts/reproduce_evidence.py's _take: one download, no
retry, the file counted only with the sha256 its manifest records), into a temporary
folder, so the cache the stages read is not touched. The BRCA2 paper's Supplementary
Table 3, whose first source has failed runs before, is asked of each of its sources once
a round, the rounds --gap-minutes apart; the other four files' sources in the rounds
--others-in names. No two requests are less than a minute apart, and the order of the
sources turns by one each round. An attempt during which this machine slept (its
wall-clock time longer than its running time by more than five seconds: a laptop's lid
closed) says nothing of the source: it is recorded as void, not counted, and made again.
The table's order is then set from what this records, by hand.

    .venv/bin/python scripts/stage_E9/measure_label_audit_sources.py \\
        --rounds 10 --gap-minutes 20 --others-in 1,5,9 --journal <a file outside the repository>

Every attempt is appended to --journal as it ends (a line of JSON); the run record
holds them all, and for each source the attempts, the successes, the seconds and each
way it failed.

A run that ended before its last round (its process stopped with the session that
started it, so that it wrote no run record) is continued with --first-round: the rounds
before it are not made again but read from the journal, where every attempt was written
as it ended, and counted with this run's in the record, which says which run made which.

    .venv/bin/python scripts/stage_E9/measure_label_audit_sources.py \\
        --rounds 10 --first-round 5 --gap-minutes 20 --others-in 1,5,9 --journal <the same file>
"""
from __future__ import annotations

import argparse
import json
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog  # noqa: E402

sys.path.insert(0, str(REPO / "scripts"))
import reproduce_evidence as R  # noqa: E402

OFTEN = "phase1/data/evidence/label_audit_cache/brca2/41586_2024_8388_MOESM3_ESM.xlsx"
APART = 60.0        # seconds between two requests, at least
SLEPT = 5.0         # an attempt whose wall-clock time exceeds its running time by more is void


def cached_entries() -> list[dict]:
    out = []
    for man in (sorted(R.LABEL_AUDIT.glob("*_raw/manifest.json"))
                + sorted(R.SECOND_SET.glob("*_raw/label_source_manifest.json"))):
        out += [e for e in json.loads(man.read_text())["files"] if e.get("cached")]
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rounds", type=int, default=10)
    ap.add_argument("--first-round", type=int, default=1,
                    help="continue a run that ended early: the rounds before this one are read from --journal")
    ap.add_argument("--gap-minutes", type=float, default=20)
    ap.add_argument("--others-in", default="1,5,9")
    ap.add_argument("--journal", required=True, type=Path)
    args = ap.parse_args()
    others_in = {int(x) for x in args.others_in.split(",") if x}
    journal = args.journal.expanduser()
    journal.parent.mkdir(parents=True, exist_ok=True)
    earlier = []
    if args.first_round > 1:
        made = [json.loads(ln) for ln in journal.read_text().splitlines() if ln.strip()]
        earlier = [a for a in made if a["round"] < args.first_round]
        have = sorted({a["round"] for a in earlier})
        if have != list(range(1, args.first_round)) or any(a["round"] >= args.first_round for a in made):
            sys.exit(f"{journal}: holds rounds {sorted({a['round'] for a in made})}; to continue from round "
                     f"{args.first_round} it must hold rounds 1 to {args.first_round - 1} and no later one")
    with RunLog(__file__) as log, tempfile.TemporaryDirectory(prefix="measure_sources_") as tmp:
        log.input(R.FETCH_SOURCES)
        table = json.loads(R.FETCH_SOURCES.read_text())
        entries = {e["path"]: e for e in cached_entries()}
        attempts = [a for a in earlier if not a.get("void")]
        voided = [a for a in earlier if a.get("void")]
        last = None
        began = time.monotonic()
        for rnd in range(args.first_round, args.rounds + 1):
            wait = began + (rnd - args.first_round) * args.gap_minutes * 60 - time.monotonic()
            if wait > 0:
                time.sleep(wait)
            files = [OFTEN] + ([p for p in entries if p != OFTEN] if rnd in others_in else [])
            for path in files:
                sources = R.label_audit_sources(entries[path], table)
                turn = (rnd - 1) % len(sources)
                for source in sources[turn:] + sources[:turn]:
                    while True:
                        if last is not None and time.monotonic() - last < APART:
                            time.sleep(APART - (time.monotonic() - last))
                        at = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
                        wall, start = time.time(), time.monotonic()
                        why, _ = R._take(source, Path(tmp) / Path(path).name, entries[path]["sha256"])
                        seconds = round(time.monotonic() - start, 1)
                        asleep = round(time.time() - wall - (time.monotonic() - start), 1)
                        last = time.monotonic()
                        for left in Path(tmp).iterdir():        # a .part of a download cut short
                            left.unlink()
                        rec = {"utc": at, "round": rnd, "file": path, "url": source["url"], "what": source["what"],
                               "zip_member": source["zip_member"], "served_the_file": why is None,
                               "failure": why, "seconds": seconds}
                        if asleep > SLEPT:
                            rec |= {"void": True, "why_void": f"this machine slept {asleep} s during the attempt"}
                            voided.append(rec)
                        else:
                            attempts.append(rec)
                        with open(journal, "a", encoding="utf-8") as f:
                            f.write(json.dumps(rec, ensure_ascii=False) + "\n")
                        print(f"{at} round {rnd}: {'VOID (' + rec['why_void'] + '); ' if rec.get('void') else ''}"
                              f"{'the file' if why is None else 'FAILED: ' + why} in {seconds} s  {source['url'][:90]}", flush=True)
                        if not rec.get("void"):
                            break
        summary = {}
        for a in attempts:
            s = summary.setdefault(a["file"] + " <- " + a["url"], {"what": a["what"], "attempts": 0, "served_the_file": 0, "seconds_when_served": [], "failures": {}, "seconds_when_failed": []})
            s["attempts"] += 1
            if a["served_the_file"]:
                s["served_the_file"] += 1
                s["seconds_when_served"].append(a["seconds"])
            else:
                key = a["failure"].split(" after 1 attempt")[-1].strip(" ()")[:160] if "gave up on" in a["failure"] else a["failure"][:160]
                s["failures"][key] = s["failures"].get(key, 0) + 1
                s["seconds_when_failed"].append(a["seconds"])
        log.note("rounds", args.rounds)
        log.note("first_round_of_this_run", args.first_round)
        if earlier:
            log.note("rounds_read_from_the_journal", {
                "rounds": sorted({a["round"] for a in earlier}), "attempts": len(earlier),
                "first_utc": min(a["utc"] for a in earlier), "last_utc": max(a["utc"] for a in earlier),
                "why": "made by an earlier run of this script, whose process ended with the session that started "
                       "it before its last round, so that it wrote no run record; each attempt was written to the "
                       "journal as it ended, and is counted here as written"})
        log.note("gap_minutes", args.gap_minutes)
        log.note("others_in_rounds", sorted(others_in))
        log.note("attempts", attempts)
        log.note("void_attempts", voided)
        log.note("by_source", summary)
        for k, s in summary.items():
            print(f"{s['served_the_file']}/{s['attempts']}  {k.split(' <- ')[1][:100]}  failures: {s['failures']}")


if __name__ == "__main__":
    main()
