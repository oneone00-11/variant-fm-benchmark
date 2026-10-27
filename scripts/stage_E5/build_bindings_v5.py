"""E5.1 (commit 1): docs/manuscript_bindings_v5.tsv, the binding table of the v5 text, from
the E5.0 number table.

E5.0, the read-only check before the revision (2026-10-27), read every number of the v5
main text and of its supplementary notes against the cell its sentence means, and wrote
one row per number, in reading order, to its number table: where the number stands, as
printed, its category (constant, matches, does not match, not found by the bound mode),
the cell, and the value there at 72dbcb1. This writes the same numbers as a binding
table for the gate's cells mode (scripts/check_evidence_manuscript_numbers.py --cells):

- a measured number gets the file, rows, column and derivation of the cell E5.0 names,
  written out once below (BINDINGS, keyed by E5.0's number id: m = main text, s =
  supplementary notes, then the rank in reading order);
- a constant gets the reason the gate gives it, or, for the five numbers the gate does
  not declare and E5.0 found to be constants, the reason below (REASONS).

It also writes tests/fixtures/manuscript_bindings_v5_expected.tsv: E5.0's verdict on
every number, and the verdict the cells mode must give at 72dbcb1 (E5.0 §2): E5.0's
"matches" and constants as they are; E5.0's "does not match", and its "not found" where
the value there differs from the print, as a mismatch with E5.0's value; and the five
"not found" numbers of E5.0 §2.3 whose value is the printed one as what they are,
constants or (the derived 71) a match.

Checks before writing; any failure stops it:
- the gate reads the same numbers from the two texts, in the same order and on the same
  lines, as E5.0's table lists them;
- every excerpt is a piece of its sentence as the gate reads it;
- every measured number has one binding and every binding one number.

    .venv/bin/python scripts/stage_E5/build_bindings_v5.py \\
        --table ~/Desktop/E5.0_数字表_20261027.tsv \\
        --main ~/Desktop/整理文件/calibration_manuscript_v5_round2_keyed.md \\
        --supplement ~/Desktop/整理文件/calibration_supplement_notes_v5_keyed.md
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog  # noqa: E402

sys.path.insert(0, str(REPO / "scripts"))
import check_evidence_manuscript_numbers as gate  # noqa: E402

OUT = REPO / "docs" / "manuscript_bindings_v5.tsv"
EXPECTED = REPO / "tests" / "fixtures" / "manuscript_bindings_v5_expected.tsv"
EXPECTED_COLUMNS = ["document", "line", "occurrence", "printed", "e5_0_category", "e5_0_now",
                    "expected_status", "expected_now"]
DOCS = {"正文": "main", "补充说明": "supplement"}
CATEGORY = {"常数": "constant", "对上": "match", "对不上": "mismatch", "找不到": "not_found"}

R = "phase1/reports/evidence/"
WT, SC, DC = R + "walker_thresholds.csv", R + "set_counts.csv", R + "depth_counts.csv"
NA, LC, TL = R + "clinvar_arm_no_assertion.csv", R + "walker_cutpoint_local_vs_band.csv", R + "tier_logo_table.csv"
DS, AR = R + "threshold_dilution_summary.csv", R + "arms_at_tool_threshold.csv"
XD, XT, ET = R + "external_ddx3x.csv", R + "external_tp53.csv", R + "evidence_thresholds.csv"
IF, IFD, IFT = (R + "inframe_attribution_summary.csv", R + "inframe_attribution_ddx3x_summary.csv",
                R + "inframe_attribution_tp53_summary.csv")
CC, TB, LF = R + "column_concordance.csv", R + "external_tp53_column_basis.csv", R + "evidence_thresholds_logo_folds.csv"
S2 = R + "supplement/tableS2_analysis_set.csv"
YAML = "phase1/config/walker2023.yaml"
DDX = "phase1/data/evidence/external/ddx3x_splice.manifest.json"
TPM = "phase1/data/evidence/external/tp53_splice_scored_12nt.manifest.json"
AVI = "phase1/data/evidence/avi_scorers_20260920.json"
AGP = "phase1/data/evidence/alphagenome_official/all_scores.tsv.provenance.json"
Q1 = '["urn:mavedb:00000658-q-1"]'
THREE = "tool=alphagenome; stratum in s3_10|s11_50|s3_50; tier in supporting|moderate|strong"


def b(file, column, rows="", derive="value", scale=1):
    return dict(file=file, rows=rows, column=column, derive=derive, scale=scale)


def wt(stratum, column, scale=1):
    return b(WT, column, f"score_column=spliceai_walker; stratum={stratum}; clinvar_arm=all; scope=pooled",
             scale=scale)


def arm(arm_, column, scale=1):
    return b(AR, column, "tool=spliceai_walker; threshold_basis=published cut point; gene_set=all_genes; "
                         f"stratum=s3_50; clinvar_arm={arm_}", scale=scale)


def xd(label, stratum, column, tool="spliceai_walker"):
    return b(XD, column, f"label_definition={label}; stratum={stratum}; tool={tool}")


def xt(label, stratum, column):
    return b(XT, column, f"label_definition={label}; stratum={stratum}; tool=spliceai_walker")


# E5.0's cell for every measured number, as a binding (E5.0 tools/claims.py, one by one)
BINDINGS = {
    # main text
    "m001": wt("s3_50", "n"),
    "m002": wt("s3_50", "lr_pp3"), "m037": wt("s3_50", "lr_pp3"),
    "m004": wt("s3_50", "lr_pp3_lo"), "m039": wt("s3_50", "lr_pp3_lo"),
    "m005": wt("s3_50", "lr_pp3_hi"), "m040": wt("s3_50", "lr_pp3_hi"),
    "m009": b(YAML, "thresholds.applied_weight.higher_score_check_cut_point"),
    "m013": b(SC, "n", derive="sum"), "m014": b(SC, "n", derive="sum"),
    "m015": b(SC, "n_labelled", derive="sum"),
    "m016": wt("pm12", "n_pos"), "m017": wt("pm12", "n"),
    "m018": b(DC, "damaging", "gene=all; depth_nt=11-30"),
    "m019": b(DC, "damaging", "gene=all; depth_nt=11-50"),
    "m028": b(NA, "n", "window_nt=3-50; clinvar_arm_strict=recorded_no_assertion; gene=BRCA1"),
    "m029": b(NA, "n", "window_nt=3-50; clinvar_arm_strict=recorded_no_assertion; gene=(all)"),
    "m030": b(DDX, "n_intronic_snv_mapped"),
    "m031": b(DDX, "max_intron_offset"),
    "m032": b(DDX, f"controls.per_scoreset{Q1}.control_auroc"),
    "m036": wt("s3_50", "frac_pp3", scale=100),
    "m041": b(YAML, "thresholds.pp3.reported_lr"),
    "m042": wt("s3_10", "lr_pp3"), "m043": wt("s11_50", "lr_pp3"),
    "m044": wt("s3_10", "lr_bp4"), "m045": wt("s11_50", "lr_bp4"), "m046": wt("s3_50", "lr_bp4"),
    "m047": b(LC, "local_lr_at_cut", "tool=spliceai_walker; stratum=s3_10"),
    "m048": b(LC, "local_lr_at_cut", "tool=spliceai_walker; stratum=s3_50"),
    "m049": b(TL, "heldout_lr_median", "tool=spliceai_walker; stratum=s11_50; tier=strong"),
    "m050": b(DS, "frac_pairs_cleared", "tier=moderate; k>=3", "min", 100),
    "m051": b(DS, "frac_pairs_cleared", "tool=pangolin; tier=strong; stratum!=s11_50", "min", 100),
    "m052": b(DS, "frac_pairs_cleared", "tool=spliceai_walker; tier=strong; stratum!=s11_50", "min", 100),
    "m053": b(DS, "frac_pairs_cleared", "tool=spliceai_walker; tier=strong; stratum!=s11_50", "max", 100),
    "m054": arm("classified", "band_lr"), "m055": arm("recorded_unclassified", "band_lr"),
    "m056": arm("unrecorded", "band_lr"),
    "m057": arm("classified", "spec", 100), "m058": arm("unrecorded", "spec", 100),
    "m059": xd("y_deposit", "s3_50", "walker_lr_pp3"),
    "m061": xd("y_deposit", "s3_50", "walker_lr_pp3_lo"),
    "m062": xd("y_deposit", "s3_50", "walker_lr_pp3_hi"),
    "m063": xd("y_deposit", "s3_10", "walker_lr_pp3"),
    "m064": xd("y_deposit", "s11_50", "walker_lr_pp3"),
    # the median of the folds' thresholds: e3_strong_lr_here in v3.0.1, renamed
    # pp3_logo_median_strong_lr_here by E18 (889e400) with the same cells
    "m065": xd("y_deposit", "s11_50", "pp3_logo_median_strong_lr_here", tool="avi_splice_site_usage"),
    "m066": xd("y_control_anchored", "s3_50", "walker_lr_pp3"),
    "m067": xd("y_control_anchored_gated", "s3_50", "walker_lr_pp3"),
    "m068": xt("y_control_anchored", "s3_10", "n_pos"), "m115": xt("y_control_anchored", "s3_10", "n_pos"),
    "m069": xt("y_control_anchored", "s3_10", "n_neg"), "m116": xt("y_control_anchored", "s3_10", "n_neg"),
    "m070": xt("y_control_anchored", "s3_10", "walker_lr_pp3"),
    "m072": xt("y_control_anchored", "s3_10", "walker_lr_pp3_lo"),
    "m073": xt("y_control_anchored", "s3_10", "walker_lr_pp3_hi"),
    "m074": xt("y_control_anchored", "s11_50", "n"),
    "m076": b(IF, "in_frame_share", "tool=spliceai; stratum=s3_10", scale=100),
    "m078": b(DC, "labelled", "gene=all; depth_nt=31-50"),
    "m079": b(ET, "pp3_threshold_insample", THREE, "max"),
    "m103": xd("y_deposit", "s3_10", "n_pos"), "m104": xd("y_deposit", "s3_10", "n_neg"),
    "m107": xd("y_deposit", "s11_50", "n_pos"), "m108": xd("y_deposit", "s11_50", "n_neg"),
    "m111": xd("y_deposit", "s3_50", "n_pos"), "m112": xd("y_deposit", "s3_50", "n_neg"),
    "m110": b(DDX, "max_intron_offset"), "m114": b(DDX, "max_intron_offset"),
    "m119": b(S2, "Variants", "Stratum (intronic offset, bp) startswith Canonical; Gene=All genes; "
                              "ClinVar record status=All variants"),
    # supplementary notes
    "s003": b(YAML, "thresholds.pp3.reported_lr"),
    "s004": b(YAML, "thresholds.bp4.reported_lr"),
    "s005": b(YAML, "thresholds.applied_weight.narrow_band_scores[0]"),
    "s006": b(YAML, "thresholds.applied_weight.narrow_band_scores[1]"),
    "s007": b(YAML, "thresholds.applied_weight.narrow_band_lr"),
    "s009": b(YAML, "thresholds.applied_weight.higher_score_check_cut_point"),
    "s012": b(CC, "spearman_rho", "column_a=spliceai; column_b=spliceai_walker"),
    "s013": b(CC, "n", "column_a=spliceai; column_b=spliceai_walker"),
    "s014": b(CC, "n_a_below_b_at_or_above", "column_a=spliceai; column_b=spliceai_walker"),
    # the interval the AlphaGenome column is scored over: recorded as the scoring code
    # in the column's provenance, so the number must appear in that text
    "s015": b(AGP, "code_run.per_variant_lines", derive="text"),
    "s016": b(TB, "analysis_set_max", "tool=alphagenome"),
    "s018": b(AVI, "n_scorers"),
    "s020": b(CC, "spearman_rho", "column_a=avi; column_b=alphagenome"),
    "s021": b(CC, "n", "column_a=avi; column_b=alphagenome"),
    "s022": b(CC, "spearman_rho", "column_a=avi_splice_sites; column_b=alphagenome"),
    "s023": b(ET, "pp3_threshold_insample", THREE, "all"),
    "s024": b(ET, "pp3_threshold_insample", "tool=alphagenome; stratum=s3_50; tier=supporting"),
    "s025": b(ET, "pp3_threshold_insample", "tool=alphagenome; stratum=s3_50; tier=strong"),
    "s033": b(LF, "heldout_lr", "tool=fusion_enet; stratum in s3_10|s11_50|s3_50; "
                                "tier in supporting|moderate|strong", "count_nonnull"),
    "s034": b(LF, "heldout_lr", "tool=fusion_enet; stratum in s3_10|s11_50|s3_50; "
                                "tier in supporting|moderate|strong", "count"),
    "s035": b(DDX, "n_intronic_snv_mapped"), "s042": b(DDX, "n_intronic_snv_mapped"),
    "s036": b(DDX, "max_intron_offset"), "s037": b(DDX, "max_intron_offset"),
    "s038": b(DDX, "strata.s3_10"), "s039": b(DDX, "strata.s11_50"),
    "s041": b(DDX, "labels.n_damaging_deposit"),
    "s043": b(DDX, f"controls.per_scoreset{Q1}.control_auroc"),
    # "0.99 to 1.0 in fifteen of the other sixteen": the sixteen score sets other than
    # q-1, of which p-1 is the one left out
    "s044": b(DDX, "controls.per_scoreset.*.control_auroc",
              "key not in urn:mavedb:00000658-q-1|urn:mavedb:00000658-p-1", "min"),
    "s045": b(DDX, "controls.per_scoreset.*.control_auroc", "key!=urn:mavedb:00000658-q-1", "max"),
    "s046": b(DDX, f"controls.per_scoreset{Q1}.n_damaging_control_anchored_in_scope"),
    "s047": b(DDX, "controls.per_scoreset.*.n_damaging_control_anchored_in_scope", derive="sum"),
    "s049": b(TPM, "n_rows"), "s050": b(TPM, "n_frozen_rows"), "s051": b(TPM, "n_new_rows"),
    "s055": xt("y_control_anchored", "s11_50", "n"),
    "s057": b(IF, "n", "tool=spliceai; stratum=s3_10"),
    "s058": b(IFD, "n", "label_definition=y_deposit; tool=spliceai; stratum=s3_10"),
    "s059": b(IFT, "n", "tool=spliceai; stratum=s3_10"),
}

# The numbers E5.0 found to be constants that the gate does not declare (E5.0 §2.3 and
# its number table): a parameter in the code, a value quoted from the literature, the
# deposit's own cut, and a digit of a tolerance written as 1e-9.
REASONS = {
    "s029": "variant-level bootstrap resamples, Pejaver et al.'s 10,000 (src/evid_interval_lr.py, N_BOOT_VARIANT)",
    "s031": "prior of Pejaver et al.'s boundaries, quoted and not used (src/evid_interval_lr.py)",
    "s040": "the deposit's own cut on its random-forest posterior (Radford et al.; src/evid_external.py)",
    "s054": "margin around zero that the third TP53 label sets aside (src/evid_external.py)",
    "s061": "the exponent in 1e-9, a relative tolerance, not a measured number",
}
POSITION = re.compile(r"^(正文|补充说明) 第 (\d+) 行「…(.*)…」$")


def read_table(path: Path) -> list[dict]:
    rows = []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            m = POSITION.match(r["位置"])
            if not m:
                raise SystemExit(f"cannot read the position {r['位置']!r}")
            rows.append(dict(document=DOCS[m.group(1)], line=int(m.group(2)), excerpt=m.group(3),
                             printed=r["印的数"], category=CATEGORY[r["类别"]], source=r["出处"],
                             now=r["现值"], reason=r["原因"]))
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", type=Path, required=True)
    ap.add_argument("--main", type=Path, required=True)
    ap.add_argument("--supplement", type=Path, required=True)
    args = ap.parse_args()
    with RunLog(__file__) as log:
        log.input(REPO / "scripts/check_evidence_manuscript_numbers.py")
        for p in (args.table, args.main, args.supplement):
            log.input(p.expanduser())
        table = read_table(args.table.expanduser())
        out, expected, used = [], [], set()
        for doc, md in (("main", args.main), ("supplement", args.supplement)):
            tokens = gate.read_tokens(md.expanduser().read_text())
            rows = [r for r in table if r["document"] == doc]
            if [(t.line, t.printed) for t in tokens] != [(r["line"], r["printed"]) for r in rows]:
                raise SystemExit(f"{doc}: the gate reads other numbers than E5.0's table lists")
            for i, (t, r) in enumerate(zip(tokens, rows)):
                nid = f"{doc[0] if doc == 'main' else 's'}{i:03d}"
                if r["excerpt"] not in t.excerpt_space:
                    raise SystemExit(f"{nid} line {t.line}: the excerpt is not in the sentence")
                row = dict(document=doc, line=t.line, occurrence=t.occurrence, printed=t.printed,
                           excerpt=r["excerpt"])
                consistent_not_found = r["category"] == "not_found" and nid in REASONS | {"s047": ""}
                if nid in BINDINGS:
                    if r["category"] == "constant":
                        raise SystemExit(f"{nid}: E5.0 calls it a constant, yet it has a binding")
                    row.update(kind="measured", reason="", **BINDINGS[nid])
                    used.add(nid)
                elif r["category"] == "constant" or nid in REASONS:
                    reason = REASONS.get(nid) or t.reason
                    if not reason:
                        raise SystemExit(f"{nid}: a constant with no reason")
                    row.update(kind="constant", reason=reason, file="", rows="", column="",
                               derive="", scale="")
                else:
                    raise SystemExit(f"{nid} line {t.line} {t.printed}: measured in E5.0, no binding")
                out.append(row)
                if r["category"] in ("constant", "match"):
                    status, now = r["category"], r["now"]
                elif consistent_not_found:
                    status, now = ("constant" if nid in REASONS else "match"), r["now"]
                else:
                    status, now = "mismatch", r["now"]
                expected.append(dict(document=doc, line=t.line, occurrence=t.occurrence,
                                     printed=t.printed, e5_0_category=r["category"], e5_0_now=r["now"],
                                     expected_status=status, expected_now=now))
        left = sorted(set(BINDINGS) - used)
        if left:
            raise SystemExit(f"bindings with no number: {left}")
        gate.write_bindings(OUT, out)
        with open(EXPECTED, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=EXPECTED_COLUMNS, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(expected)
        log.output(OUT)
        log.output(EXPECTED)
        counts = {}
        for e in expected:
            k = (e["document"], e["e5_0_category"], e["expected_status"])
            counts[" / ".join(k)] = counts.get(" / ".join(k), 0) + 1
        log.note("rows", len(out))
        log.note("e5_0_category_to_expected_status", counts)
        for k, v in sorted(counts.items()):
            print(f"{v:4d}  {k}")
        print(f"wrote {len(out)} rows to {OUT.relative_to(REPO)} and {EXPECTED.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
