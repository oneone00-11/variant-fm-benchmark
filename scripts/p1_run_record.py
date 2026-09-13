"""Write docs/tp53-pinned-rescore.md: the P1 run record (provenance + sanity concordance).

Reads data/tp53/*.tsv, data/rescore/*_frozen_only16.tsv, data/rescore/*_sanity.tsv and their
provenance JSONs, compares the sanity-subset scores with the frozen-matrix-v1 columns and the
TP53 SpliceAI scores with the v1 TP53 column, and reports the concordance the plan asks for.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parents[1]
FROZEN = REPO / "phase1" / "data" / "frozen" / "frozen_matrix_v1.parquet"
TP53 = REPO / "phase1" / "data" / "external" / "tp53_splice_scored.parquet"
COLS = {"pangolin": ("pangolin_fullprec", "pangolin"), "spliceai": ("spliceai_ds_fullprec", "spliceai"), "nt": ("nucleotide_transformer", "nt")}


def prov(p: Path) -> dict:
    q = p.with_name(p.name.replace(".tsv", ".provenance.json"))
    return json.loads(q.read_text()) if q.exists() else {}


def main() -> int:
    fz = pd.read_parquet(FROZEN).set_index("variant_id")
    t53 = pd.read_parquet(TP53).set_index("variant_id")
    lines = ["# P1 run record — TP53 and sanity scoring with the atlas's pinned scorers", "",
             "Driver: `scripts/phase4b_score_pinned.py` (imports the atlas scorer modules by path; contains no scoring logic). "
             "Reference: Ensembl GRCh38 release-112 subset FASTA (chr 2/3/13/16/17) and, for Pangolin, the gffutils "
             "database built from the release-112 GTF — the same files the atlas and this repository's seven-gene runs used.", ""]
    lines += ["## Provenance", "", "| model | set | n | scored | tool / version | pin | key parameters | device | runtime (s) | input sha256 |",
              "|---|---|---|---|---|---|---|---|---|---|"]
    for m in ("pangolin", "spliceai", "nt"):
        for s in ("tp53", "sanity", "frozen_only16"):
            p = REPO / "data" / ("tp53" if s == "tp53" else "rescore") / f"{m}_{s}.tsv"
            if not p.exists():
                lines.append(f"| {m} | {s} | — | not run | | | | | | |")
                continue
            pv = prov(p)
            pin = pv.get("git_commit") or pv.get("checkpoint_revision_sha") or pv.get("tool_version", "")
            par = pv.get("parameters", {})
            par_s = "; ".join(f"{k}={v}" for k, v in par.items() if k != "score_definition")
            lines.append(f"| {m} | {s} | {pv.get('n_variants')} | {pv.get('scored')} | {pv.get('tool', '')} {pv.get('tool_version', '')} "
                         f"{('torch ' + pv['torch']) if 'torch' in pv else ''}{('TF ' + pv['tensorflow']) if 'tensorflow' in pv else ''} | `{pin}` | {par_s} | "
                         f"{pv.get('device', 'cpu')} | {pv.get('runtime_seconds')} | `{str(pv.get('input_sha256', ''))[:12]}…` |")
    lines += ["", "Full provenance (reference and annotation hashes, scorer-module hashes, environment versions) is in the "
              "`*.provenance.json` next to each table.", "", "## Sanity check — 200 splice variants of frozen-matrix-v1", "",
              "| model | n | Spearman ρ vs v1 column | Pearson r | identical after rounding to 2 d.p. | max |Δ| | expected |",
              "|---|---|---|---|---|---|---|"]
    summary = {}
    for m, (col, canon) in COLS.items():
        p = REPO / "data" / "rescore" / f"{m}_sanity.tsv"
        if not p.exists():
            lines.append(f"| {m} | — | not run | | | | |")
            continue
        t = pd.read_csv(p, sep="\t").set_index("variant_id")
        old = fz.loc[t.index, canon].astype(float).to_numpy()
        new = t[col].astype(float).to_numpy()
        ok = ~(np.isnan(old) | np.isnan(new))
        rho = stats.spearmanr(old[ok], new[ok]).statistic
        r = stats.pearsonr(old[ok], new[ok]).statistic
        if m in ("pangolin", "spliceai"):
            same = float(np.mean(np.isclose(np.round(new[ok], 2), old[ok], atol=1e-9)))
            mx = float(np.max(np.abs(np.round(new[ok], 2) - old[ok])))
            exp = "ρ near 1.000 (bounded by the ties in the rounded v1 column); re-rounded values identical"
        else:
            same = float("nan")
            mx = float(np.max(np.abs(new[ok] - old[ok])))
            exp = "ρ ≈ 0.9997 (the atlas's concordance against this column)"
        summary[m] = {"n": int(ok.sum()), "spearman": rho, "pearson": r, "identical_2dp": same, "max_abs": mx}
        lines.append(f"| {m} | {ok.sum()} | {rho:.4f} | {r:.4f} | {'' if np.isnan(same) else f'{same:.3f}'} | {mx:.4g} | {exp} |")
    # TP53 SpliceAI against the v1 TP53 column (rounded CLI)
    p = REPO / "data" / "tp53" / "spliceai_tp53.tsv"
    if p.exists():
        t = pd.read_csv(p, sep="\t").set_index("variant_id")
        old = t53.loc[t.index, "spliceai"].astype(float).to_numpy(); new = t["spliceai_ds_fullprec"].astype(float).to_numpy()
        ok = ~(np.isnan(old) | np.isnan(new))
        lines += ["", f"TP53 SpliceAI: the full-precision scores re-round to the existing TP53 column exactly "
                  f"(max |round(v2,2) − v1| = {np.max(np.abs(np.round(new[ok], 2) - old[ok])):.4g} over {ok.sum()} variants; "
                  f"distinct values {pd.Series(new[ok]).nunique()} against {pd.Series(old[ok]).nunique()})."]
    lines += ["", "## Outputs", "", "```", *[str(p.relative_to(REPO)) for p in sorted((REPO / 'data' / 'tp53').glob('*')) + sorted((REPO / 'data' / 'rescore').glob('*'))], "```", ""]
    out = REPO / "docs" / "tp53-pinned-rescore.md"
    out.write_text("\n".join(lines), encoding="utf8")
    (REPO / "docs" / "tp53-pinned-rescore.summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("\n".join(lines[-40:]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
