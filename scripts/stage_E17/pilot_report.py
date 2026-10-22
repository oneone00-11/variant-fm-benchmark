"""E17.3: the pilot of AlphaGenome's official splicing score, as docs/alphagenome-official-pilot.md.

Reads what the pilot wrote (phase1/data/evidence/alphagenome_official/: the definition's
record, the environment, the notebook's examples, the 120 pilot scores, the 10 scored
again, and the read-only check of the Atlas columns) and the analysis set's two other
AlphaGenome columns. Reports the definition item by item with the notebook lines it
rests on, the examples, the determinism check, missing values, the distribution, the
Spearman correlation with the other two columns (values only, no threshold), throughput
and the projection for the whole analysis set. It computes no likelihood ratio,
threshold or tier (E17.0 rule 9).

Run from the repository root:  .venv/bin/python scripts/stage_E17/pilot_report.py
"""
from __future__ import annotations

import json
import sys
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel  # noqa: E402

DIR = REPO / "phase1" / "data" / "evidence" / "alphagenome_official"
SET = REPO / "phase1" / "data" / "evidence" / "analysis_set_v1.parquet"
OUT = REPO / "docs" / "alphagenome-official-pilot.md"
COMPONENTS = ["SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS"]
COMPOSITE = "alphagenome_splicing"
VALUES = COMPONENTS + [COMPOSITE]


def _read(name: str):
    import pandas as pd
    return pd.read_csv(DIR / name, sep="\t", dtype=str, keep_default_na=False)


def _num(frame, cols):
    out = frame.copy()
    for c in cols:
        out[c] = out[c].astype(float)
    return out


def g(x: float) -> str:
    return f"{x:.6g}"


def main() -> None:
    import pandas as pd
    from scipy.stats import spearmanr
    with RunLog(__file__) as log:
        names = ["NOTEBOOK_SOURCE.json", "ENVIRONMENT.json", "ATLAS_COLUMNS_CHECK.json",
                 "examples_scores.tsv", "pilot_scores.tsv", "repeat_scores.tsv",
                 "pilot_scores.tsv.provenance.json", "repeat_scores.tsv.provenance.json",
                 "examples_scores.tsv.provenance.json"]
        for n in names:
            log.input(DIR / n)
        src = json.loads((DIR / "NOTEBOOK_SOURCE.json").read_text())
        env = json.loads((DIR / "ENVIRONMENT.json").read_text())
        atlas = json.loads((DIR / "ATLAS_COLUMNS_CHECK.json").read_text())
        p_prov = json.loads((DIR / "pilot_scores.tsv.provenance.json").read_text())
        r_prov = json.loads((DIR / "repeat_scores.tsv.provenance.json").read_text())
        ex_text, pilot_text, rep_text = (_read(n) for n in
                                         ("examples_scores.tsv", "pilot_scores.tsv",
                                          "repeat_scores.tsv"))
        pilot = _num(pilot_text, VALUES)
        aset = pd.read_parquet(SET, columns=["variant_id", "gene", "stratum", "alphagenome",
                                              "alphagenome_v061"])
        log.input(SET)
        md = ["# AlphaGenome by its official definition: the pilot", "",
              "Written by `scripts/stage_E17/pilot_report.py` from "
              "`phase1/data/evidence/alphagenome_official/`; do not edit by hand. Stage E17, "
              "commit 1. It computes no likelihood ratio, threshold or tier (E17.0 rule 9): "
              "it shows whether the implementation is the notebook's, not how the score does.",
              ""]
        # ---------------------------------------------------------------- definition
        lic = src["licence"]
        md += ["## The definition and where it comes from", "",
               f"- Notebook: `{src['path']}` in {src['repository']}, commit "
               f"`{src['commit']}` (the last commit touching the file: "
               f"`{src['last_commit_touching_the_file']['sha'][:12]}`, "
               f"{src['last_commit_touching_the_file']['date']}, \"{src['last_commit_touching_the_file']['message']}\").",
               f"- sha256 `{src['sha256']}` ({src['bytes']:,} bytes); read "
               f"{src['accessed_utc']}. The repository's licence: `LICENSE`, "
               f"\"{' '.join(lic['first_lines'])}\" ({lic['spdx_by_github']}), sha256 "
               f"`{lic['sha256']}`.",
               f"- Printed outputs in the notebook: {src['printed_outputs']} "
               "(none: there is no reference output to reproduce).",
               f"- The quoted cells, verbatim, with their numbers: "
               f"`phase1/data/evidence/alphagenome_official/NOTEBOOK_SOURCE.json` "
               f"(cells {', '.join(str(c['cell']) for c in src['quoted_cells'])}).", "",
               "| item | what the quoted code fixes | cell:line |", "|---|---|---|"]
        for d in src["definition"]:
            refs = ", ".join(f"{q['cell']}:{q['line']}" for q in d["quoted"])
            md.append(f"| {d['item']} | {d['value']} | {refs} |")
        cl = src["client"]
        md += ["", f"Read from the installed client ({cl['package']} {cl['version']}): "
               f"`SEQUENCE_LENGTH_1MB` = {cl['SEQUENCE_LENGTH_1MB']:,}; organism "
               f"`{cl['organism']}`; the scorers "
               + "; ".join(f"`{s['repr']}` (signed: {s['is_signed']})" for s in cl["scorers"])
               + ".", "",
               "Arithmetic: the client returns raw scores as float32, so cell 13 takes the "
               "maxima and adds and divides in float32. The components and the composite "
               "are float32 values, written exactly (each number reads back to the same "
               "bits), and the composite is (sites + usage) + junctions / 5 in float32 to the "
               "last bit (`phase1/tests/test_evid_alphagenome_official.py`). Read them back "
               "exactly: pandas' default float parser is not round-trip exact.", ""]
        # ---------------------------------------------------------------- environment
        md += ["## The scoring environment", "",
               f"- Client: {env['client']['package']} {env['client']['version']}. The "
               "notebook installs it unpinned (cell 3), so it is the newest release on PyPI on "
               "the day of scoring; `scripts/stage_E17/requirements.txt` pins it and its "
               f"{len(env['pip_freeze'])} packages, each with the sha256 of its file "
               f"(sha256 of the file: `{env['requirements_sha256']}`), and "
               "`.venv-ag-official/` (gitignored) is built from it.",
               f"- Python: {env['python'].split()[0]}; {env['platform']}.",
               f"- Model version: requested {env['model']['model_version_requested']} "
               f"({env['model']['model_version_requested_note']}); in the response: "
               f"{env['model']['model_version_in_response']} (the response's metadata keys: "
               f"{', '.join(env['model']['response_metadata_keys'])}).",
               f"- API key: present {env['api_key_present']}, from "
               f"`{env['api_key_source']}`; its value is recorded nowhere.", ""]
        # ---------------------------------------------------------------- examples
        md += ["## The notebook's own examples", "",
               "The notebook prints no output (无参照输出), so its three examples have no "
               "reference to match. Scored here with the quoted code as written (cell 7's key "
               "aside), they show that code running on this client:", "",
               "| variant | " + " | ".join(VALUES) + " |", "|---|" + "---|" * len(VALUES)]
        for r in ex_text.itertuples():
            md.append(f"| `{r.variant_id}` | " + " | ".join(getattr(r, c) for c in VALUES) + " |")
        md.append("")
        # ---------------------------------------------------------------- sample
        share = p_prov["sample"]["share"]
        md += ["## The pilot sample", "",
               f"{p_prov['n_variants']} variants of `analysis_set_v1`, an even share in each "
               f"non-empty gene x stratum cell, drawn with `config.RANDOM_SEED` "
               f"({p_prov['sample']['seed']}) by `scripts/stage_E17/score_official.py`:", "",
               "| gene | pm12 | s3_10 | s11_50 |", "|---|---|---|---|"]
        for gene in sorted({k.split("/")[0] for k in share}):
            md.append(f"| {gene} | " + " | ".join(str(share.get(f"{gene}/{s}", 0))
                                                  for s in ("pm12", "s3_10", "s11_50")) + " |")
        md.append("")
        # ---------------------------------------------------------------- determinism
        merged = pilot_text.merge(rep_text, on="variant_id", suffixes=("_pilot", "_repeat"))
        same = merged.apply(lambda r: all(r[f"{c}_pilot"] == r[f"{c}_repeat"] for c in VALUES),
                            axis=1)
        gap = (datetime.fromisoformat(r_prov["start_utc"])
               - datetime.fromisoformat(p_prov["end_utc"])).total_seconds() / 60
        md += ["## Determinism", "",
               f"{len(merged)} of the pilot's variants (drawn with the same seed) were scored "
               f"again {gap:.0f} minutes after the pilot ended ({p_prov['end_utc']} -> "
               f"{r_prov['start_utc']}). Identical to the last bit in all four values: "
               f"**{int(same.sum())} of {len(merged)}**.", "",
               "| variant | " + " | ".join(VALUES) + " | identical |",
               "|---|" + "---|" * (len(VALUES) + 1)]
        for r, s in zip(merged.itertuples(index=False), same):
            row = r._asdict()
            md.append(f"| `{row['variant_id']}` | "
                      + " | ".join(row[f"{c}_pilot"] for c in VALUES)
                      + f" | {'yes' if s else 'NO'} |")
        md.append("")
        # ---------------------------------------------------------------- missing
        md += ["## Missing values", "",
               "| component | variants with no row | NaN raw scores | smallest raw score |",
               "|---|---|---|---|"]
        for c in COMPONENTS:
            rows = pilot_text[f"rows_{c}"].astype(int)
            nans = pilot_text[f"nan_{c}"].astype(int)
            mins = pd.to_numeric(pilot_text[f"min_raw_{c}"], errors="coerce")
            md.append(f"| {c} | {int((rows == 0).sum())} | {int(nans.sum())} | {g(mins.min())} |")
        md += ["", f"Composite missing: {int(pilot[COMPOSITE].isna().sum())} of {len(pilot)}.",
               ""]
        # ---------------------------------------------------------------- distribution
        md += ["## Distribution", "",
               "| value | min | 25% | median | 75% | max |", "|---|---|---|---|---|---|"]
        for c in VALUES:
            q = pilot[c].quantile([0, .25, .5, .75, 1]).tolist()
            md.append(f"| {c} | " + " | ".join(g(x) for x in q) + " |")
        md += ["", "| stratum | n | median composite | max composite |", "|---|---|---|---|"]
        for s in ("pm12", "s3_10", "s11_50"):
            sub = pilot[pilot.stratum == s][COMPOSITE]
            md.append(f"| {s} | {len(sub)} | {g(sub.median())} | {g(sub.max())} |")
        md.append("")
        # ---------------------------------------------------------------- correlation
        both = pilot.merge(aset[["variant_id", "alphagenome", "alphagenome_v061"]],
                           on="variant_id")
        md += ["## Spearman correlation with the other two AlphaGenome columns", "",
               "Over the pilot's variants, values only, no threshold:", "",
               "| with | n | Spearman rho |", "|---|---|---|"]
        for other, label in (("alphagenome", "`alphagenome` (the atlas quantile composite; "
                                             "`alphagenome_atlas_quantile` from commit 3)"),
                             ("alphagenome_v061", "`alphagenome_v061` (client 0.6.1)")):
            ok = both[[COMPOSITE, other]].dropna()
            rho = spearmanr(ok[COMPOSITE], ok[other]).statistic
            md.append(f"| {label} | {len(ok)} | {rho:.3f} |")
            log.note(f"spearman_{other}", {"n": len(ok), "rho": float(rho)})
        md.append("")
        # ---------------------------------------------------------------- throughput
        rate = p_prov["n_scored_this_run"] / p_prov["seconds_scoring"]
        n_all = len(aset)
        md += ["## Throughput", "",
               f"The pilot scored {p_prov['n_scored_this_run']} variants in "
               f"{p_prov['seconds_scoring']:.0f} s with {p_prov['workers']} workers "
               f"({rate:.2f} variants a second), with {p_prov['retries']} retries and "
               f"{len(p_prov['failures'])} failures. At that rate the whole analysis set, "
               f"{n_all:,} variants, takes about {n_all / rate / 60:.0f} minutes and "
               f"{n_all:,} calls of `score_variant`, one per variant, before any retry.", ""]
        # ---------------------------------------------------------------- atlas columns
        md += ["## The four Atlas columns: raw scores or quantiles (read-only)", "",
               f"{', '.join('`' + c + '`' for c in atlas['columns'])}: **{atlas['answer']}**.",
               f"- The record: `{atlas['provenance']}` (client {atlas['client_version']}, "
               f"scorers {', '.join(atlas['scorers_requested'])}; aggregation: "
               f"\"{atlas['aggregation_recorded']}\").",
               f"- Our code: `{atlas['our_code']['file']}`, `{atlas['our_code']['function']}`, "
               f"reads `{atlas['our_code']['reads']}`.",
               f"- The client's code: `{atlas['client_source']['file']}` in "
               f"`{atlas['client_source']['wheel']}` (sha256 "
               f"`{atlas['client_source']['wheel_sha256']}`):"]
        for ln in atlas["client_source"]["lines"]:
            md.append(f"  - line {ln['line']}, {ln['shows']}: `{ln['text']}`")
        md.append("")
        # ---------------------------------------------------------------- per variant
        md += ["## Every pilot variant", "",
               "| variant | gene | stratum | " + " | ".join(VALUES) + " |",
               "|---|---|---|" + "---|" * len(VALUES)]
        for r in pilot_text.itertuples():
            md.append(f"| `{r.variant_id}` | {r.gene} | {r.stratum} | "
                      + " | ".join(g(float(getattr(r, c))) for c in VALUES) + " |")
        md.append("")
        OUT.write_text("\n".join(md))
        log.output(OUT)
        log.note("determinism", f"{int(same.sum())} of {len(merged)}")
        print(f"wrote {rel(OUT)}: determinism {int(same.sum())}/{len(merged)}, "
              f"{n_all / rate / 60:.0f} min projected")


if __name__ == "__main__":
    main()
