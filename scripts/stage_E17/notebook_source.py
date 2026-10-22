"""E17.1: where the official definition of the AlphaGenome splicing score comes from.

The definition is the one in AlphaGenome's own notebook for scoring variants on their
splicing effect, in google-deepmind/alphagenome at a pinned commit. This fetches the
notebook and the repository's LICENSE at that commit, checks the notebook against the
sha256 recorded when it was first read, quotes verbatim the cells that compute the
score (with their cell numbers), and writes, item by item, what the quoted code fixes:
the client call, the scorers, the input window, the organism, tissue or ontology, the
score field, the aggregation, the combination, missing values, the strand and the
variant. Each item names the cells and lines it rests on, and this script checks that
each named line holds the code it is said to hold. The client's own values (the window
constant, the scorers as the client builds them, whether they are signed) are read here
from the installed client, so it runs in .venv-ag-official/. The whole notebook is not
committed; phase1/data/evidence/alphagenome_official/NOTEBOOK_SOURCE.json is the record.

Run from the repository root:
    .venv-ag-official/bin/python scripts/stage_E17/notebook_source.py
"""
from __future__ import annotations

import hashlib
import json
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel, utc_now  # noqa: E402

OUT = REPO / "phase1" / "data" / "evidence" / "alphagenome_official" / "NOTEBOOK_SOURCE.json"
REPOSITORY = "google-deepmind/alphagenome"
COMMIT = "55943112b83439e0c2512110f1ae607d4b5b492a"      # main on 2026-10-22
PATH = "colabs/v1/splicing_variant_scoring.ipynb"
SHA256 = "c39e6a6a09bc62512107c106cd43cabf7a016cbd073474ff872c669a1ad838df"
QUOTED_CODE = [3, 5, 7, 9, 11, 13, 15]       # the cells that install, load and score
QUOTED_MARKDOWN = [12]                       # the merged score's wording, for the sign

# item, value, and the (cell, line, text the line holds) it rests on
DEFINITION = [
    ("client", "alphagenome, unpinned in the notebook: the newest release on PyPI on the "
               "day of scoring",
     [(3, 5, "! pip install alphagenome")]),
    ("client call", "dna_client.create(api key), then dna_model.score_variant(interval, "
                    "variant, variant_scorers, organism) per variant. The notebook reads "
                    "the key with colab_utils.get_api_key(); here it comes from "
                    "phase1/src/evid_score_avi.read_key() (E17.0 rule 1), the only line "
                    "of the quoted code that is not run as written",
     [(7, 1, "dna_model = dna_client.create("), (11, 13, "dna_model.score_variant("),
      (15, 27, "dna_model.score_variant(")]),
    ("scorers", "variant_scorers.RECOMMENDED_VARIANT_SCORERS['SPLICE_SITES'], "
                "['SPLICE_SITE_USAGE'] and ['SPLICE_JUNCTIONS']",
     [(9, 2, "RECOMMENDED_VARIANT_SCORERS['SPLICE_SITES']"),
      (9, 3, "RECOMMENDED_VARIANT_SCORERS['SPLICE_SITE_USAGE']"),
      (9, 4, "RECOMMENDED_VARIANT_SCORERS['SPLICE_JUNCTIONS']")]),
    ("input window", "dna_client.SEQUENCE_LENGTH_1MB nucleotides, an interval of that "
                     "length centred on the variant (variant.reference_interval.resize)",
     [(11, 10, "variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)"),
      (15, 26, "variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)")]),
    ("organism", "dna_client.Organism.HOMO_SAPIENS",
     [(11, 17, "organism=dna_client.Organism.HOMO_SAPIENS"),
      (15, 31, "organism=dna_client.Organism.HOMO_SAPIENS")]),
    ("tissue or ontology", "none: score_variant is called with no ontology selection, "
                           "and the maximum runs over every track",
     [(11, 16, "variant_scorers=splicing_scorers,"),
      (13, 21, "df.groupby(['variant_id', 'output_type'])['raw_score']")]),
    ("score field", "raw_score, not quantile_score",
     [(13, 21, "['raw_score']")]),
    ("aggregation", "for each variant and scorer, the maximum of raw_score over all its "
                    "rows: every gene, every track and, for the junction scorer, every "
                    "junction. The three scorers are unsigned (is_signed False), so this "
                    "maximum is the maximum absolute score the notebook's text names",
     [(13, 21, "df.groupby(['variant_id', 'output_type'])['raw_score']"),
      (13, 22, ".max()"), (9, 8, "(signed={scorer.is_signed})"),
      (12, 4, "scorer, we take the maximum absolute score across all genes and tracks")]),
    ("combination", "SPLICE_SITES + SPLICE_SITE_USAGE + SPLICE_JUNCTIONS / 5.0",
     [(13, 30, "max_scores.get('SPLICE_SITES', 0.0)"),
      (13, 31, "+ max_scores.get('SPLICE_SITE_USAGE', 0.0)"),
      (13, 32, "+ max_scores.get('SPLICE_JUNCTIONS', 0.0) / 5.0")]),
    ("missing values", "a NaN raw_score is skipped by the maximum (pandas); a scorer with "
                       "no row for a variant counts 0.0 (fillna), and one with no row for "
                       "any variant 0.0 (the .get default)",
     [(13, 22, ".max()"), (13, 25, ".fillna(0.0)"),
      (13, 30, "max_scores.get('SPLICE_SITES', 0.0)")]),
    ("strand", "none chosen: rows of genes on either strand all enter the maximum",
     [(13, 21, "df.groupby(['variant_id', 'output_type'])['raw_score']")]),
    ("variant", "genome.Variant(chromosome 'chrN', position 1-based, reference_bases, "
                "alternate_bases)",
     [(11, 2, "variant = genome.Variant("), (11, 3, "chromosome='chr13',"),
      (11, 4, "position=32316462,")]),
]


def _get(url: str, accept: str | None = None) -> bytes:
    headers = {"User-Agent": "variant-fm-benchmark E17"}
    if accept:
        headers["Accept"] = accept
    with urllib.request.urlopen(urllib.request.Request(url, headers=headers), timeout=60) as r:
        return r.read()


def main() -> None:
    with RunLog(__file__) as log:
        raw = f"https://raw.githubusercontent.com/{REPOSITORY}/{COMMIT}/{PATH}"
        blob = _get(raw)
        got = hashlib.sha256(blob).hexdigest()
        if got != SHA256:
            raise SystemExit(f"{PATH} at {COMMIT[:12]} hashes to {got}, not {SHA256}")
        licence = _get(f"https://raw.githubusercontent.com/{REPOSITORY}/{COMMIT}/LICENSE")
        history = json.loads(_get(f"https://api.github.com/repos/{REPOSITORY}/commits"
                                  f"?path={PATH}&sha={COMMIT}&per_page=1",
                                  "application/vnd.github+json"))
        nb = json.loads(blob)
        cells = nb["cells"]
        text = {i: "".join(c.get("source", [])) for i, c in enumerate(cells)}
        printed = sum(len(c.get("outputs", [])) for c in cells)
        definition = []
        for item, value, refs in DEFINITION:
            for cell, line, holds in refs:
                lines = text[cell].split("\n")
                if holds not in lines[line - 1]:
                    raise SystemExit(f"{item}: cell {cell} line {line} does not hold "
                                     f"{holds!r}: {lines[line - 1]!r}")
            definition.append({"item": item, "value": value,
                               "quoted": [{"cell": c, "line": n, "holds": h}
                                          for c, n, h in refs]})
        # what the installed client makes of the quoted code (cells 5 and 9, as written)
        ns: dict = {}
        exec(text[5], ns)
        exec(text[9], ns)
        dna_client = ns["dna_client"]
        import alphagenome
        client = {"package": "alphagenome", "version": alphagenome.__version__,
                  "SEQUENCE_LENGTH_1MB": dna_client.SEQUENCE_LENGTH_1MB,
                  "organism": str(dna_client.Organism.HOMO_SAPIENS),
                  "scorers": [{"repr": repr(s), "name": s.name, "is_signed": s.is_signed}
                              for s in ns["splicing_scorers"]]}
        if any(s["is_signed"] for s in client["scorers"]):
            raise SystemExit("a splicing scorer is signed: the maximum of raw_score and the "
                             "maximum absolute score differ (E17.0 rule 8)")
        record = {
            "repository": f"https://github.com/{REPOSITORY}",
            "commit": COMMIT,
            "path": PATH,
            "raw_url": raw,
            "sha256": got,
            "bytes": len(blob),
            "last_commit_touching_the_file": {
                "sha": history[0]["sha"], "date": history[0]["commit"]["committer"]["date"],
                "message": history[0]["commit"]["message"].split("\n")[0]},
            "licence": {"file": "LICENSE", "sha256": hashlib.sha256(licence).hexdigest(),
                        "first_lines": [l.strip() for l in licence.decode().splitlines()
                                        if l.strip()][:2],
                        "spdx_by_github": json.loads(_get(
                            f"https://api.github.com/repos/{REPOSITORY}/license",
                            "application/vnd.github+json"))["license"]["spdx_id"]},
            "accessed_utc": utc_now(),
            "printed_outputs": printed,
            "printed_outputs_note": ("none: every cell's outputs are empty, so the "
                                     "notebook's examples have no printed reference"
                                     if printed == 0 else "see the notebook"),
            "quoted_cells": [{"cell": i, "type": cells[i]["cell_type"], "source": text[i]}
                             for i in sorted(QUOTED_CODE + QUOTED_MARKDOWN)],
            "definition": definition,
            "client": client,
        }
        OUT.parent.mkdir(parents=True, exist_ok=True)
        OUT.write_text(json.dumps(record, indent=1, ensure_ascii=False) + "\n")
        log.output(OUT)
        log.note("sha256", got)
        log.note("client", client)
        print(f"{PATH} @ {COMMIT[:12]}: sha256 ok, {len(definition)} items, "
              f"{printed} printed outputs -> {rel(OUT)}")
        for s in client["scorers"]:
            print(f"  {s['repr']} signed={s['is_signed']}")


if __name__ == "__main__":
    main()
