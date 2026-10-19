"""E16 pre-check of the four deposits as fetched, before any table is built.

Offline: reads only phase1/data/evidence/external/<gene>_raw/. For each gene it
counts the intron-side SNVs in the score table by distance into the intron, checks
how many of them have a GRCh38 postMapped record in variant_details.ndjson, applies
the investigator-provided ranges mechanically, and looks for ranges that overlap,
leave gaps, or have a score sitting exactly on a shared boundary. Everything it
finds goes into its run record; nothing is written elsewhere.

Run from the repository root:
    .venv/bin/python scripts/stage_E16/probe_raw.py
"""
from __future__ import annotations

import csv
import json
import math
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog  # noqa: E402

sys.path.insert(0, str(REPO / "phase1"))
from src import evid_external as E  # noqa: E402
from src import evid_newgenes as N  # noqa: E402

EXT = REPO / "phase1" / "data" / "evidence" / "external"
GENES = N.GENES
HGVS_NT = N.HGVS_NT
matches = N._matches          # the stage's own reading of a range and its bounds


def bin_of(depth: int) -> str:
    """The stage's layers: evid_common's names, and out_of_scope_deep past 50 nt."""
    return E.stratum_of(depth) if depth <= N.MAX_DEPTH else N.DEEP


def overlaps_and_gaps(classes: list[dict]) -> dict:
    """Order the intervals by lower end and compare neighbours."""
    key = lambda c: -math.inf if c["range"][0] is None else c["range"][0]   # noqa: E731
    cs = sorted(classes, key=key)
    found = {"overlaps": [], "gaps": [], "unbounded_ends": []}
    if cs and cs[0]["range"][0] is not None:
        found["unbounded_ends"].append(f"lowest class {cs[0]['label']!r} has a lower bound")
    if cs and cs[-1]["range"][1] is not None:
        found["unbounded_ends"].append(f"highest class {cs[-1]['label']!r} has an upper bound")
    for a, b in zip(cs, cs[1:]):
        hi_a, lo_b = a["range"][1], b["range"][0]
        if hi_a is None or lo_b is None:
            found["overlaps"].append(f"{a['label']!r} and {b['label']!r}: an open end "
                                     "inside the order")
            continue
        if hi_a > lo_b or (hi_a == lo_b and a["inclusiveUpperBound"]
                           and b["inclusiveLowerBound"]):
            found["overlaps"].append({"classes": [a["label"], b["label"]], "at": hi_a,
                                      "upper_inclusive": a["inclusiveUpperBound"],
                                      "lower_inclusive": b["inclusiveLowerBound"]})
        elif hi_a < lo_b or not (a["inclusiveUpperBound"] or b["inclusiveLowerBound"]):
            found["gaps"].append({"between": [a["label"], b["label"]], "from": hi_a,
                                  "to": lo_b})
    return found


def probe(gene: str) -> dict:
    raw = EXT / f"{gene.lower()}_raw"
    man = json.loads((raw / "manifest.json").read_text())
    meta = json.loads((raw / "scoreset_meta.json").read_text())
    target = meta["targetGenes"][0]["targetAccession"]["accession"]
    csv_path = raw / f"urn_mavedb_{man['urn'].split(':')[-1]}.csv"
    rows = list(csv.DictReader(csv_path.open(newline="")))
    details = {}
    with (raw / "variant_details.ndjson").open() as f:
        for line in f:
            if line.strip():
                r = json.loads(line)
                details[r["urn"]] = r
    ranges = json.loads((raw / "score_ranges.json").read_text())
    inv = ranges["investigator_provided"]
    out = {"urn": man["urn"], "n_rows": len(rows), "n_variant_details": len(details),
           "target_accession": target, "n_investigator_calibrations": len(inv)}

    intronic, other_tx = [], 0
    for r in rows:
        m = HGVS_NT.search(r["hgvs_nt"] or "")
        if not m:
            continue
        if m.group(1) != target:
            other_tx += 1
        intronic.append((r, int(m.group(3))))
    out["n_intronic_snv"] = len(intronic)
    out["n_intronic_on_other_transcript"] = other_tx
    out["deepest_offset"] = max((d for _, d in intronic), default=None)
    out["intronic_by_depth"] = dict(Counter(bin_of(d) for _, d in intronic))

    mapped = [(r, d) for r, d in intronic if (details.get(r["accession"]) or {}).get("postMapped")]
    out["intronic_with_postMapped"] = len(mapped)
    out["intronic_without_postMapped"] = len(intronic) - len(mapped)
    unmapped_kinds = Counter()
    for r in rows:
        if r["accession"] in details:
            continue
        h = r["hgvs_nt"] or ""
        kind = ("intronic SNV" if HGVS_NT.search(h) else "deletion" if "del" in h
                else "other")
        unmapped_kinds[kind] += 1
    out["rows_without_a_variant_details_record"] = dict(unmapped_kinds)
    out["postMapped_reference_accessions"] = dict(Counter(
        e["value"].split(":")[0] for d in details.values() if d.get("postMapped")
        for e in d["postMapped"].get("expressions", []) if e.get("syntax") == "hgvs.g"))

    missing_score = sum(1 for r, _ in intronic if r["score"] in ("", "NA", None))
    out["intronic_missing_score"] = missing_score
    if "rna_score" in rows[0]:
        out["intronic_with_rna_score"] = sum(1 for r, _ in intronic
                                             if r["rna_score"] not in ("", "NA", None))
    else:
        out["intronic_with_rna_score"] = "no rna_score column"

    if len(inv) != 1:
        out["ranges_problem"] = f"{len(inv)} investigator-provided calibrations"
        return out
    classes = inv[0]["functionalClassifications"]
    out["classes"] = [{"label": c["label"], "functionalClassification":
                       c["functionalClassification"], "range": c["range"],
                       "inclusiveLowerBound": c["inclusiveLowerBound"],
                       "inclusiveUpperBound": c["inclusiveUpperBound"]} for c in classes]
    out["range_structure"] = overlaps_and_gaps(classes)
    bounds = {b for c in classes for b in c["range"] if b is not None}
    on_bound = Counter()
    by_class = Counter()
    multi = none = 0
    abnormal_by_bin = Counter()
    for r, d in intronic:
        if r["score"] in ("", "NA", None):
            continue
        x = float(r["score"])
        if x in bounds:
            on_bound[x] += 1
        hit = [c for c in classes if matches(x, c)]
        if len(hit) > 1:
            multi += 1
        elif not hit:
            none += 1
        else:
            by_class[hit[0]["label"]] += 1
            if hit[0]["functionalClassification"] == "abnormal":
                abnormal_by_bin[bin_of(d)] += 1
    out["intronic_scores_exactly_on_a_boundary"] = {str(k): v for k, v in on_bound.items()}
    out["intronic_in_more_than_one_class"] = multi
    out["intronic_in_no_class"] = none
    out["intronic_by_class"] = dict(by_class)
    out["intronic_abnormal_by_depth"] = dict(abnormal_by_bin)
    return out


def main() -> None:
    with RunLog(__file__) as log:
        results = {}
        for gene in GENES:
            raw = EXT / f"{gene.lower()}_raw"
            for p in sorted(raw.iterdir()):
                log.input(p)
            results[gene] = probe(gene)
        log.note("per_gene", results)
        print(json.dumps(results, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
