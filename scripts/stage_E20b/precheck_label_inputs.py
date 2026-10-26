"""E20b.0: the committed inputs a one-off script wrote while reading a label column, derived
again from HEAD's analysis set and label tables by the same code, and compared with the
committed files variant by variant. Nothing is written into the repository but this run's
record; the derived files go to a temporary folder.

The inputs (scripts/reproduce_inputs.txt) a one-off script wrote while reading a label
column (y_assay, or an external gene's label columns), and the code that wrote each:

  inframe_subset.parquet, inframe_subset_pangolin.parquet          E6
      evid_inframe.write_subset (--subset; y_assay), called with its two paths redirected
  pangolin_events_work/variants.csv
      evid_score_pangolin_events.emit_csv, from the Pangolin subset just derived
  external/ddx3x_inframe_subset.parquet                            E2.7
      evid_inframe.external_subset("ddx3x") (--external-subset; DDX3X's four labels),
      called with its path redirected
  external/tp53_inframe_subset.parquet                             E2.7, the 10 nt subset
      evid_tp53_extend.inframe_subset, which derives it again from the 12 nt table on
      every run of its stage and checks it against this file; called on the stage's table
  inframe_subset_tp.parquet, inframe_subset_pangolin_tp.parquet,   E20.1
  inframe_score_input_tp.parquet, inframe_score_input_pangolin_tp.parquet,
  external/ddx3x_inframe_subset_tp.parquet, external/tp53_inframe_subset_12nt_tp.parquet
      scripts/stage_E20/write_c12_subsets.py's own functions (seven_genes, external, anchors)
  pangolin_events_work_tp/variants.csv
      emit_csv, from the Pangolin score input just derived
  external/newgenes_scored.parquet                                 E18.2
      scripts/stage_E18/build_columns.assemble: the five genes' rows with their labels

Reviewed and left out, with the reason (written into the run record):
  * the event files: the scorers' output on the subsets, derived again only by scoring
    (E6's: E20b commit 1);
  * external/ddx3x_splice.parquet: evid_external --prepare ddx3x builds DDX3X's labels
    from the deposit; it reads no label column;
  * the other one-off inputs (E16-E19's score inputs, anchors, records and extracts): their
    scripts read no label column;
  * files a stage writes on every run (tp53_inframe_subset_12nt.parquet, the five genes'
    <gene>_variants.parquet) are not inputs, and the comparison report
    label_change_effect.csv is not one either.

Run in the reference image (修改单 §0 rule 9), offline:
    docker run ... IMAGE bash scripts/stage_E20b/in_image.sh scripts/stage_E20b/precheck_label_inputs.py
"""
from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel  # noqa: E402

PHASE1 = REPO / "phase1"
sys.path.insert(0, str(PHASE1))
os.chdir(PHASE1)
from src import evid_common as K  # noqa: E402
from src import evid_inframe as F  # noqa: E402
from src import evid_score_pangolin_events as PG  # noqa: E402
from src import evid_tp53_extend as T  # noqa: E402

sys.path.insert(1, str(REPO / "scripts/stage_E20"))
import write_c12_subsets as W  # noqa: E402
sys.path.insert(1, str(REPO / "scripts/stage_E18"))
import build_columns as B  # noqa: E402

EV = PHASE1 / "data/evidence"
LEFT_OUT = {
    "event files (inframe_events*.parquet, external/*_inframe_events*.parquet)":
        "the scorers' output on the subsets, derived again only by scoring them (E6's: E20b commit 1)",
    "external/ddx3x_splice.parquet":
        "evid_external --prepare ddx3x builds DDX3X's labels from the deposit; it reads no label column",
    "E16-E19's other one-off inputs (score inputs, anchors, records, extracts)":
        "their scripts read no label column (searched: y_assay, y_*, label_cols, PRIMARY_LABEL)",
    "tp53_inframe_subset_12nt.parquet, external/<gene>_variants.parquet":
        "written by a stage on every run, so not inputs",
    "label_change_effect.csv": "a comparison report (COMPARISON_REPORTS), not an input",
}


def compare(name: str, committed: pd.DataFrame, derived: pd.DataFrame, key: str = "variant_id") -> dict:
    a, b = committed.set_index(key), derived.set_index(key)
    only_c = sorted(set(a.index) - set(b.index))
    only_d = sorted(set(b.index) - set(a.index))
    common = sorted(set(a.index) & set(b.index))
    cols_c, cols_d = list(a.columns), list(b.columns)
    differing = {}
    for c in [c for c in cols_c if c in cols_d]:
        x = a.loc[common, c].to_numpy(dtype=object)
        y = b.loc[common, c].to_numpy(dtype=object)
        same = (pd.isna(x) & pd.isna(y)) | (x == y)
        if not same.all():
            differing[c] = [common[i] for i in np.where(~same)[0]][:20]
    order_same = (list(committed[key]) == [v for v in derived[key] if v in set(committed[key])])
    rec = {"file": name, "committed_rows": int(len(committed)), "derived_rows": int(len(derived)),
           "columns_same": cols_c == cols_d, "only_committed": only_c, "only_derived": only_d,
           "values_differ": differing, "order_of_common_rows_same": bool(order_same)}
    rec["identical"] = (not only_c and not only_d and not differing and rec["columns_same"]
                        and len(committed) == len(derived) and order_same)
    return rec


def csv_lines(path: Path) -> list[str]:
    return path.read_text().splitlines()


def main() -> None:
    tmp = Path(tempfile.mkdtemp(prefix="e20b_precheck_"))
    results = []
    with RunLog(__file__) as log:
        df = K.load_set()
        log.input(K.SET_PATH)
        # E6, by write_subset itself
        F.SUBSET, F.PANG_SUBSET = tmp / "inframe_subset.parquet", tmp / "inframe_subset_pangolin.parquet"
        F.write_subset()
        for name in ("inframe_subset.parquet", "inframe_subset_pangolin.parquet"):
            log.input(EV / name)
            results.append(compare(name, pd.read_parquet(EV / name), pd.read_parquet(tmp / name)))
        # the Pangolin work folder's variant list, from the subset just derived
        PG.emit_csv(pd.read_parquet(tmp / "inframe_subset_pangolin.parquet"), tmp / "variants.csv")
        old, new = csv_lines(EV / "pangolin_events_work/variants.csv"), csv_lines(tmp / "variants.csv")
        log.input(EV / "pangolin_events_work/variants.csv")
        results.append({"file": "pangolin_events_work/variants.csv", "committed_rows": len(old) - 1,
                        "derived_rows": len(new) - 1, "only_committed": sorted(set(old) - set(new)),
                        "only_derived": sorted(set(new) - set(old)), "identical": old == new})
        # E2.7, DDX3X, by external_subset itself
        spec = dict(F.EXTERNAL["ddx3x"])
        F.EXTERNAL["ddx3x"] = spec | {"subset": tmp / "ddx3x_inframe_subset.parquet"}
        F.external_subset("ddx3x")
        log.input(EV / "external/ddx3x_inframe_subset.parquet")
        results.append(compare("external/ddx3x_inframe_subset.parquet",
                               pd.read_parquet(EV / "external/ddx3x_inframe_subset.parquet"),
                               pd.read_parquet(tmp / "ddx3x_inframe_subset.parquet")))
        # E2.7, TP53's 10 nt subset, by the stage's own check
        table = pd.read_parquet(EV / "external/tp53_splice_scored_12nt.parquet")
        sub = T.inframe_subset(table)           # asserts the frozen part equals the file
        frozen = pd.read_parquet(T.FROZEN_SUBSET)
        mine = sub[sub["variant_id"].isin(pd.read_parquet(T.FROZEN_TABLE)["variant_id"])]
        log.input(T.FROZEN_SUBSET)
        results.append(compare("external/tp53_inframe_subset.parquet",
                               frozen, mine[list(frozen.columns)].reset_index(drop=True)))
        # E20.1, by write_c12_subsets' own functions
        for tool, col in (("spliceai", "spliceai_walker"), ("pangolin", "pangolin")):
            tp = W.seven_genes(df, col, 1, W.BANDS)
            name = W.OUT[tool].name
            log.input(EV / name)
            results.append(compare(name, pd.read_parquet(EV / name), tp))
            inp = pd.concat([tp[W.SCORER_COLS].assign(role="true_positive"), W.anchors(tool)],
                            ignore_index=True)
            name = W.OUT[f"{tool}_input"].name
            log.input(EV / name)
            results.append(compare(name, pd.read_parquet(EV / name), inp))
            if tool == "pangolin":
                PG.emit_csv(inp, tmp / "variants_tp.csv")
                old = csv_lines(EV / "pangolin_events_work_tp/variants.csv")
                new = csv_lines(tmp / "variants_tp.csv")
                log.input(EV / "pangolin_events_work_tp/variants.csv")
                results.append({"file": "pangolin_events_work_tp/variants.csv",
                                "committed_rows": len(old) - 1, "derived_rows": len(new) - 1,
                                "only_committed": sorted(set(old) - set(new)),
                                "only_derived": sorted(set(new) - set(old)), "identical": old == new})
        for gene in ("ddx3x", "tp53"):
            name = "external/" + W.OUT[gene].name
            log.input(EV / name)
            results.append(compare(name, pd.read_parquet(EV / name), W.external(gene, 1, W.BANDS)))
        # E18.2, the five genes with their labels, by build_columns.assemble
        built, _ = B.assemble()
        log.input(EV / "external/newgenes_scored.parquet")
        results.append(compare("external/newgenes_scored.parquet",
                               pd.read_parquet(EV / "external/newgenes_scored.parquet"), built))
        log.note("compared", results)
        log.note("left_out", LEFT_OUT)
        differ = [r["file"] for r in results if not r["identical"]]
        log.note("differing", differ)
        for r in results:
            extra = "" if r["identical"] else (
                f"  only committed {len(r['only_committed'])}, only derived {len(r['only_derived'])}, "
                f"values differ in {sorted(r.get('values_differ', {}))}")
            print(f"{'identical' if r['identical'] else 'DIFFERS  '}  {r['file']}{extra}")


if __name__ == "__main__":
    main()
