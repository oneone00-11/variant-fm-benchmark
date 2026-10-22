"""E17.5: when the pipeline switched to AlphaGenome's own definition, did only what
depends on AlphaGenome change?

Every output of this tree (the comparator's own list: the tracked files under
scripts/reproduce_evidence.py's output folders, and any new file there that git does
not ignore) is compared with BASE, the commit before the switch, which scored the
analysis set but did not read the scores:

  * a table whose rows name a predictor (PREDICTOR_COLUMNS): every row that is not
    about AlphaGenome must be the same, byte for byte, in every column that is not
    AlphaGenome's own, and in the same order; AlphaGenome's own rows may differ, in
    number too (a curve has a point per distinct score). A text report or page, line by
    line, the same way;
  * the fusion's and `alphagenome_v061`'s outputs must be the same, byte for byte:
    their rows by the rule above, and the files that hold nothing else (MUST_HOLD);
  * a manifest may change only in the keys JSON_KEYS names, each with its reason; the
    supplementary tables' manifest only where a table it lists, or a source it hashes,
    changed for AlphaGenome;
  * any other change (a cell outside AlphaGenome's rows and columns, a page's line that
    does not name AlphaGenome) must carry the reason it depends on AlphaGenome
    (REASONS); a change without one fails;
  * a figure must be the same, byte for byte, unless it draws AlphaGenome (AG_FIGURES).

Then the comparison table's old definition, `alphagenome_atlas_quantile`, against
BASE's outputs for `alphagenome`, which was that column then: float for float against
E3's threshold table, E2.2's held-out table, E4's territory metrics and E15's practical
metrics; and as printed, against Table 3 and Supplementary Tables S6 and S7 (E17.5;
the author's decision of 2026-10-22).

The verdicts go to the run log (scripts/stage_E17/run_log/) and to stdout; the exit
status is 1 if anything fails. Run from the repository root, on the outputs written in
the reference image:
    .venv/bin/python scripts/stage_E17/check_alphagenome_only.py
"""
from __future__ import annotations

import csv
import difflib
import importlib.util
import io
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog  # noqa: E402

PHASE1 = REPO / "phase1"
sys.path.insert(0, str(PHASE1))

BASE = "34a6c8c"                 # E17's full scoring commit: the scores, not yet read
OLD = "alphagenome_atlas_quantile"
NAME = "AlphaGenome splice score"
AG = {"alphagenome", OLD, NAME}
PREDICTOR_COLUMNS = {"tool", "feature", "column_a", "column_b", "object", "panel_column",
                     "Tool", "Score column", "Column A", "Column B", "Predictor"}
EV = "phase1/reports/evidence/"
COMPARISON = EV + "alphagenome_definition_comparison.csv"
NEW = {COMPARISON: "E17's comparison of the two definitions (src.evid_alphagenome_definitions)"}
# outputs that hold the fusion or alphagenome_v061 and nothing else
MUST_HOLD = ([EV + "fusion_fold_coefficients.csv", EV + "supplement/tableS9_fusion.csv"]
             + [EV + f"interval_lr_fusion_enet_{s}.csv"
                for s in ("pm12", "s3_10", "s11_50", "s3_50", "all_1_50")]
             + ["phase1/data/evidence/external/tp53_12nt_alphagenome_v061.parquet"])
# the figures that draw AlphaGenome, and where (figure 1 draws SpliceAI alone; figure
# 4's AlphaGenome row is empty, since DDX3X has no AlphaGenome column and TP53's is not
# comparable; neither may change)
_MAIN = "the four main predictors (src.evid_practical_metrics.MAIN_TOOLS)"
AG_FIGURES = {
    EV + "figures/figure2.pdf": "interval curves of four predictors, AlphaGenome's panel among them",
    EV + "figures/figure2.png": "the same figure as a PNG",
    EV + "figures/figureS1.pdf": "interval curves of every column (src.evid_figures.COLUMNS)",
    EV + "figures/figureS1.png": "the same figure as a PNG",
    EV + "figures/figure3.pdf": "held-out ratios, a row per column, AlphaGenome's among them",
    EV + "figures/figure3.png": "the same figure as a PNG",
    EV + "figures/figureS2.pdf": "thresholds refitted on subsets of genes, an AlphaGenome line",
    EV + "figures/figureS2.png": "the same figure as a PNG",
    EV + "figures/figureS3.pdf": "a row per column, AlphaGenome's among them",
    EV + "figures/figureS3.png": "the same figure as a PNG",
    EV + "figures/figure_practical.pdf": f"ROC curves and tier marks of {_MAIN}",
    EV + "figures/figure_practical.png": "the same figure as a PNG",
    EV + "figures/figure_practical_slides.png": f"the slide version, {_MAIN}",
    EV + "fig_data/fig2_draft.png": "the draft of figure 2, AlphaGenome's curve among four",
    EV + "fig_data/fig3_draft.png": "pooled AUROC of every column by territory",
}
# reasons for changes outside AlphaGenome's own rows, columns and lines, each checked
# against the run (E17 report, section 6)
REASONS = {
    "docs/practical-metrics-summary.md":
        "numbered statements 2 to 4 are taken over the four main predictors, AlphaGenome "
        "among them: held out, the highest Strong sensitivity at 3-10 bp (39%) and the "
        "highest Strong coverage there (8.1%) were AlphaGenome's and are now Pangolin's "
        "(38%, 7.5%); at a specificity of 0.95 AlphaGenome was the most sensitive, or "
        "tied for it, in all three bands, and now no predictor is",
    EV + "supplement/tableS11_within_gene_arms_summary.csv":
        "counts, per gene, of the columns whose AUROC is highest in each ClinVar arm; at "
        "3-50 bp in BAP1 the AlphaGenome splice score's highest arm moved from unrecorded "
        "(0.818) to recorded unclassified (0.811, a margin below the printed precision), "
        "so one column moves from one count to the other, in the BAP1 row and in the "
        "row of the four genes",
}
# manifests: the keys that may change, and why
_TP53_WORDING = ("the AlphaGenome wording E17 corrected (src.evid_tp53_extend, "
                 "src.evid_score_avi); the manifest hashes its code")
JSON_KEYS = {
    "phase1/data/evidence/analysis_set_v1.manifest.json": {
        "built_utc": "build time",
        "sha256": "the analysis set's AlphaGenome columns (checked above)",
        "columns": "alphagenome_atlas_quantile, after alphagenome",
        "sources.alphagenome_official": "the score file the AlphaGenome column now comes from",
    },
    "phase1/data/evidence/external/tp53_splice_scored_12nt.manifest.json": {
        "columns.alphagenome.definition": "names the analysis-set column the TP53 column is "
                                          "not comparable with, AlphaGenome's own definition",
        "inputs_sha256.phase1/src/evid_tp53_extend.py": _TP53_WORDING,
        "inputs_sha256.phase1/src/evid_score_avi.py": _TP53_WORDING,
    },
}
# files this commit changes for AlphaGenome that the supplementary tables read
SOURCES_CHANGED = {"LICENSE-DATA", "docs/column-provenance.md",
                   "phase1/data/evidence/alphagenome_official/NOTEBOOK_SOURCE.json",
                   "phase1/data/evidence/alphagenome_official/all_scores.tsv.provenance.json"}


def _reproduce():
    spec = importlib.util.spec_from_file_location("reproduce_evidence",
                                                  REPO / "scripts/reproduce_evidence.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _git(*args: str, data: bytes | None = None) -> bytes:
    return subprocess.run(["git", "-C", str(REPO), *args], input=data, check=True,
                          capture_output=True).stdout


def _names(*args: str, roots: list[str]) -> list[str]:
    return [n for n in _git(*args, "--", *roots).decode().split("\n") if n]


def base_blobs(paths: list[str]) -> dict[str, bytes]:
    """Each path as BASE has it, through one `git cat-file --batch`."""
    out = _git("cat-file", "--batch", data="".join(f"{BASE}:{p}\n" for p in paths).encode())
    blobs, pos = {}, 0
    for p in paths:
        nl = out.index(b"\n", pos)
        size = int(out[pos:nl].split()[2])
        blobs[p] = out[nl + 1: nl + 1 + size]
        pos = nl + 1 + size + 1
    return blobs


def _ag_line(line: str) -> bool:
    low = line.lower()
    return "alphagenome" in low and "v061" not in low and "0.6.1" not in low


def compare_table(old: bytes, new: bytes) -> list[str]:
    """The rows not about AlphaGenome, in order, outside AlphaGenome's own columns:
    what differs among them ([] if nothing)."""
    a = list(csv.reader(io.StringIO(old.decode("utf-8"))))
    b = list(csv.reader(io.StringIO(new.decode("utf-8"))))
    if not a or not b or a[0] != b[0]:
        return ["the header changed"]
    head = a[0]
    pred = [i for i, h in enumerate(head) if h in PREDICTOR_COLUMNS]
    keep = [i for i, h in enumerate(head)
            if not (h.strip() in AG or (h.startswith("alphagenome_") and "v061" not in h))]

    def others(rows):
        return [(k, [r[i] for i in keep]) for k, r in enumerate(rows[1:], start=2)
                if not any(r[i].strip() in AG for i in pred)]
    oa, ob = others(a), others(b)
    if len(oa) != len(ob):
        return [f"{len(oa)} rows not about AlphaGenome, now {len(ob)}"]
    return [f"line {ka}: " + ", ".join(f"{head[keep[i]]}: {x!r} -> {y!r}"
                                       for i, (x, y) in enumerate(zip(ra, rb)) if x != y)[:300]
            for (ka, ra), (_, rb) in zip(oa, ob) if ra != rb]


def compare_text(old: bytes, new: bytes) -> list[str]:
    """The lines that do not name AlphaGenome, in order: what differs among them, as
    the lines removed (-) and added (+)."""
    a = [l for l in old.decode("utf-8").splitlines() if not _ag_line(l)]
    b = [l for l in new.decode("utf-8").splitlines() if not _ag_line(l)]
    out = []
    for tag, i1, i2, j1, j2 in difflib.SequenceMatcher(None, a, b, autojunk=False).get_opcodes():
        if tag != "equal":
            out += [f"- {x[:160]}" for x in a[i1:i2]] + [f"+ {y[:160]}" for y in b[j1:j2]]
    return out


def _flat(x, p: str = "") -> dict:
    if not isinstance(x, dict):
        return {p: x}
    out = {}
    for k, v in x.items():
        out.update(_flat(v, f"{p}.{k}" if p else k))
    return out


def compare_json(path: str, old: bytes, new: bytes) -> list[str]:
    """Every changed key must be one JSON_KEYS allows for this file."""
    allowed = JSON_KEYS.get(path, {})
    a, b = _flat(json.loads(old)), _flat(json.loads(new))
    bad = []
    for k in sorted(set(a) | set(b)):
        if a.get(k) == b.get(k):
            continue
        if not any(k == x or k.startswith(x + ".") for x in allowed):
            bad.append(f"{k} changed")
    if path.endswith("analysis_set_v1.manifest.json"):
        cols = list(a["columns"])
        cols.insert(cols.index("alphagenome") + 1, OLD)
        if b["columns"] != cols:
            bad.append("the column list changed beyond alphagenome_atlas_quantile")
    return bad


def check_analysis_set(old: bytes) -> list[str]:
    """`alphagenome` holds the new scores and the old ones sit beside it, value for
    value, as alphagenome_atlas_quantile; every other column as it was."""
    a = pd.read_parquet(io.BytesIO(old))
    b = pd.read_parquet(PHASE1 / "data/evidence/analysis_set_v1.parquet")
    cols = list(a.columns)
    cols.insert(cols.index("alphagenome") + 1, OLD)
    bad = []
    if list(b.columns) != cols or len(a) != len(b):
        bad.append("the columns or rows changed beyond alphagenome_atlas_quantile")
    else:
        bad += [f"{c} changed" for c in a.columns if c != "alphagenome" and not a[c].equals(b[c])]
        if not np.array_equal(a["alphagenome"].to_numpy(dtype=float),
                              b[OLD].to_numpy(dtype=float), equal_nan=True):
            bad.append("alphagenome_atlas_quantile is not BASE's alphagenome")
    return bad


def explain_supp_manifest(old: bytes, new: bytes, accepted: set[str]) -> tuple[list[str], list[str]]:
    """Row by row and field by field: every source whose sha256 moved, and every table
    whose own sha256 moved, must be a change accepted here or a file SOURCES_CHANGED
    names. Returns (what changed, what is not explained)."""
    a = pd.read_csv(io.BytesIO(old), dtype=str, keep_default_na=False).set_index("File")
    b = pd.read_csv(io.BytesIO(new), dtype=str, keep_default_na=False).set_index("File")
    if list(a.index) != list(b.index) or list(a.columns) != list(b.columns):
        return [], ["the manifest's tables or columns changed"]
    notes, bad = [], []
    for t in a.index:
        ra, rb = a.loc[t], b.loc[t]
        fields = [c for c in a.columns if ra[c] != rb[c]]
        if not fields:
            continue
        src_a = dict(zip(ra["Source files"].split("; "), ra["Source sha256"].split("; ")))
        src_b = dict(zip(rb["Source files"].split("; "), rb["Source sha256"].split("; ")))
        moved = [s for s in src_b if src_a.get(s) != src_b[s]]
        unexplained = [s for s in moved if s not in accepted | SOURCES_CHANGED]
        if "File sha256" in fields and EV + "supplement/" + t not in accepted:
            unexplained.append(t)
        bad += [f"{t}: {s}" for s in unexplained]
        notes.append(f"{t}: {', '.join(fields)}" + (f"; sources: {', '.join(moved)}" if moved else ""))
    return notes, bad


def check_old_definition(blob) -> list[str]:
    """The comparison table's old definition against BASE's outputs for `alphagenome`:
    float for float against E3, E2.2, E4 and E15's files, and as printed in Table 3 and
    Supplementary Tables S6 and S7. Returns the disagreements."""
    from src import evid_supp_tables as S
    from src.evid_practical_metrics import BP4_TIERS
    from src.evid_tables import _folds, STRATUM_LABEL as T3_LABEL

    def frame(rel: str, **kw) -> pd.DataFrame:
        return pd.read_csv(io.BytesIO(blob(EV + rel)), float_precision="round_trip", **kw)

    ev = frame("evidence_thresholds.csv")
    tl = frame("tier_logo_table.csv")
    tm = frame("territory_metrics.csv")
    pm = frame("practical_metrics_by_threshold.csv")
    text = {"dtype": str, "keep_default_na": False}
    t3 = frame("tables/table3_evidence_tiers.csv", **text)
    s6 = frame("supplement/tableS6_thresholds_pathogenic.csv", **text)
    s7 = frame("supplement/tableS7_thresholds_benign.csv", **text)
    t = pd.read_csv(REPO / COMPARISON, float_precision="round_trip")
    bad = []

    def same(x, y) -> bool:
        return (pd.isna(x) and pd.isna(y)) or x == y

    for r in t.itertuples():
        side, st, tier = r.side.lower(), r.stratum, r.tier
        got = {c: getattr(r, f"{OLD}_{c}") for c in (
            "reached_in_sample", "in_sample_threshold", "band_lr_at_threshold",
            "heldout_k_of_n", "heldout_lr_median", "pooled_auroc")}
        got["reached_in_sample"] = str(got["reached_in_sample"]) == "True"
        e = ev[(ev.tool == "alphagenome") & (ev.stratum == st) & (ev.tier == tier)].iloc[0]
        want = {"reached_in_sample": bool(e[f"{side}_threshold_reachable"]),
                "in_sample_threshold": e[f"{side}_threshold_insample"],
                "pooled_auroc": tm[(tm.tool == "alphagenome") & (tm.stratum == st)
                                   & (tm.clinvar_arm == "all")]["auroc"].iloc[0]}
        if side == "pp3":
            want["heldout_k_of_n"] = _folds(
                tl[(tl.tool == "alphagenome") & (tl.stratum == st) & (tl.tier == tier)],
                "alphagenome", st, tier)
        elif want["reached_in_sample"]:
            want["heldout_lr_median"] = e["bp4_logo_heldout_lr_median"]
        if want["reached_in_sample"] and (side == "pp3" or tier in BP4_TIERS):
            kind = tier if side == "pp3" else f"bp4_{tier}"
            p = pm[(pm.tool == "alphagenome") & (pm.stratum == st) & (pm.side == side)
                   & (pm.threshold_kind == kind) & (pm.threshold_source == "in_sample")]
            want["band_lr_at_threshold"] = p["band_lr_plus" if side == "pp3"
                                             else "bp4_band_lr"].iloc[0]
        bad += [f"{st} {tier} {side} {c}: {got[c]!r} in the table, {w!r} at BASE"
                for c, w in want.items() if not same(got[c], w)]
        # as printed
        tab, n = (s6, 6) if side == "pp3" else (s7, 7)
        row = tab[(tab.Tool == NAME) & (tab["Stratum (intronic offset, bp)"] == S.STRATUM_LABEL[st])
                  & (tab.Tier == tier.capitalize())].iloc[0]
        printed = S.fmt_thr(got["in_sample_threshold"], OLD, "upper" if side == "pp3" else "lower")
        if row["In-sample threshold (gene-clustered bootstrap)"] != printed:
            bad.append(f"{st} {tier} {side}: prints {printed!r}, BASE's S{n} "
                       f"{row['In-sample threshold (gene-clustered bootstrap)']!r}")
        if side == "bp4" and pd.notna(got["heldout_lr_median"]) and \
                row["Held-out LR, median over folds"] != S.fmt_lr(got["heldout_lr_median"]):
            bad.append(f"{st} {tier} bp4: held-out median prints "
                       f"{S.fmt_lr(got['heldout_lr_median'])!r}, BASE's S7 "
                       f"{row['Held-out LR, median over folds']!r}")
        if side == "pp3" and tier in ("moderate", "strong"):
            cell = t3.loc[t3.Predictor == NAME, f"{T3_LABEL[st]}: {tier.capitalize()}, held out"].iloc[0]
            if cell != got["heldout_k_of_n"]:
                bad.append(f"{st} {tier}: k/n {got['heldout_k_of_n']!r}, BASE's Table 3 {cell!r}")
    return bad


def main() -> int:
    R = _reproduce()
    roots = [str(r.relative_to(REPO)) for r in R.OUTPUT_ROOTS]
    keep = lambda names: [n for n in names if R._in_roots(REPO / n)]  # noqa: E731
    base = base_blobs(keep(_names("ls-tree", "-r", "--name-only", BASE, roots=roots)))
    now = set(keep(_names("ls-files", roots=roots)
                   + _names("ls-files", "--others", "--exclude-standard", roots=roots)))
    with RunLog(__file__) as log:
        log.note("base", _git("rev-parse", BASE).decode().strip())
        verdict: dict[str, tuple[str, str]] = {}
        for p in sorted(set(base) | now):
            if p not in now:
                verdict[p] = ("FAIL", "removed")
            elif p not in base:
                verdict[p] = (("ok", f"new: {NEW[p]}") if p in NEW
                              else ("FAIL", "new, not expected"))
            elif base[p] == (REPO / p).read_bytes():
                verdict[p] = ("identical", "")
        changed = [p for p in sorted(set(base) & now) if p not in verdict]
        manifest = EV + "supplement/supplement_tables_manifest.csv"
        for p in changed:
            if p == manifest:
                continue                                     # after the files it lists
            new, bad = (REPO / p).read_bytes(), None
            if p in MUST_HOLD:
                bad, why = ["the fusion's or v0.6.1's own output changed"], ""
            elif p.endswith((".png", ".pdf")):
                bad, why = ([] if p in AG_FIGURES else ["a figure without AlphaGenome changed"],
                            f"draws AlphaGenome: {AG_FIGURES.get(p, '')}")
            elif p.endswith(".parquet") and p.endswith("analysis_set_v1.parquet"):
                bad, why = check_analysis_set(base[p]), ("AlphaGenome's columns only: the "
                           "official scores in its place, the composite beside it")
            elif p.endswith(".json"):
                bad = compare_json(p, base[p], new)
                why = "; ".join(f"{k}: {v}" for k, v in JSON_KEYS.get(p, {}).items())
            elif p.endswith((".csv", ".tsv")):
                bad, why = compare_table(base[p], new), "AlphaGenome's rows or columns only"
            elif p.endswith((".txt", ".md")):
                bad, why = compare_text(base[p], new), "AlphaGenome's lines only"
            else:
                bad, why = ["changed"], ""
            if bad and p in REASONS:
                bad, why = [], (f"{why}; besides, {len(bad)} other change(s) "
                                f"({'; '.join(bad)[:300]}): {REASONS[p]}")
            verdict[p] = ("FAIL", "; ".join(bad)[:500]) if bad else ("ok", why)
        if manifest in changed:
            accepted = {p for p, (v, _) in verdict.items() if v in ("ok", "identical")}
            notes, bad = explain_supp_manifest(base[manifest], (REPO / manifest).read_bytes(),
                                               accepted)
            verdict[manifest] = (("FAIL", "; ".join(bad)) if bad else
                                 ("ok", "only the tables and sources that changed for "
                                        "AlphaGenome: " + " | ".join(notes)))
        old_bad = check_old_definition(lambda rel: base[rel])
        counts = pd.Series([v for v, _ in verdict.values()]).value_counts().to_dict()
        log.note("counts", counts)
        log.note("changed", {p: list(v) for p, v in verdict.items() if v[0] != "identical"})
        log.note("old_definition_disagreements", old_bad)
        print(f"outputs compared with {BASE}: {counts}")
        for p, (v, why) in verdict.items():
            if v != "identical":
                print(f"  {v:4s} {p}  {why}")
        print(f"the old definition against {BASE}'s outputs, float for float and as "
              f"printed: {'every number agrees' if not old_bad else f'{len(old_bad)} disagree'}")
        for b in old_bad:
            print("  ", b)
        failed = [p for p, (v, _) in verdict.items() if v == "FAIL"]
        log.note("failed", failed)
        return 1 if failed or old_bad else 0


if __name__ == "__main__":
    sys.exit(main())
