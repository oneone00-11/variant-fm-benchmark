"""E1 -- rebuild the analysis set for the evidence-strength reframe.

What changes from the set the published calibration study used:

  * the window is the intron-side splice region at 1 <= |offset| <= 50, not
    |offset| <= 8, and it is stratified pm12 / s3_10 / s11_50 -- the strata of the
    companion atlas, so the two papers cut the territory the same way;
  * a ClinVar record is no longer required for entry, so the set can carry a
    `unrecorded` arm at all;
  * a complete predictor row is no longer required for entry;
  * scores come from the atlas matrix, which is the wider of the two freezes.

Nothing here fits anything. It selects, joins labels, annotates ClinVar, checks
orientation, and writes a manifest.

Run (PYTHONPATH=phase1):  python -m src.evid_build_set
"""
from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from . import config as C
from .phase1_method_b_labels import _norm_hgvs

ATLAS_REPO = Path(os.environ.get(
    "evidence-strength_ATLAS_REPO", "/Users/cliffzhang/work/functional-standard-atlas"))
ATLAS_MATRIX = ATLAS_REPO / "results" / "score_matrix_atlas_v2.parquet"
ATLAS_V061 = ATLAS_REPO / "results" / "alphagenome_v061_scores.parquet"

OUT_DIR = Path("data/evidence")
REPORT_DIR = Path("reports/evidence")
OUT_PARQUET = OUT_DIR / "analysis_set_v1.parquet"
OUT_MANIFEST = OUT_DIR / "analysis_set_v1.manifest.json"

CLINVAR_VCF = OUT_DIR / "clinvar" / "clinvar_20260615.vcf.gz"

# atlas region bucket -> our stratum name
STRATA = {"splice_1_2": "pm12", "splice_3_10": "s3_10", "splice_11_50": "s11_50"}
STRATUM_ORDER = ["pm12", "s3_10", "s11_50"]

# canonical name -> atlas column. AlphaMissense (zero splice coverage) and the two
# gnomAD allele-frequency columns (independent ACMG evidence: BA1/BS1/PM2, not
# PP3) are deliberately out of the panel; see the brief, E1.
PANEL = {
    "spliceai": "spliceai_ds",
    "pangolin": "pangolin_score",
    "alphagenome": "alphagenome",
    "cadd": "cadd",
    "phylop": "phylop100way",
    "phastcons": "phastcons100way",
    "gpn_msa": "gpn_msa",
    "nt": "nucleotide_transformer",
}
# scored later, merged if the file is present
OPTIONAL_COLUMNS = {
    "spliceai_walker": OUT_DIR / "spliceai_walker.parquet",
    "avi": OUT_DIR / "avi_scores.parquet",
}

CLINVAR_ARMS = ["classified", "recorded_unclassified", "unrecorded"]


# ---------------------------------------------------------------------------
# provenance helpers
# ---------------------------------------------------------------------------
def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def git_head(repo: Path) -> str:
    try:
        return subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return "unknown"


def _load_atlas_module(name: str, relpath: str):
    """Import an atlas module by path, as phase4b_score_pinned does."""
    spec = importlib.util.spec_from_file_location(name, ATLAS_REPO / relpath)
    mod = importlib.util.module_from_spec(spec)
    sys.modules.setdefault(name, mod)
    spec.loader.exec_module(mod)
    return mod


def _load_clnsig_classifier():
    """The project's own CLNSIG collapse, so the new annotation uses the rule the
    published frozen matrix was built with."""
    spec = importlib.util.spec_from_file_location(
        "_legacy_scripts_config", Path("../scripts/config.py").resolve())
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.classify_clnsig


# ---------------------------------------------------------------------------
# ClinVar
# ---------------------------------------------------------------------------
def annotate_clinvar(keys: set[tuple[str, int, str, str]], vcf: Path) -> pd.DataFrame:
    """One pass over the ClinVar VCF, keeping the records our variants match.

    A variant is `recorded` if the release carries a record at its exact
    chrom/pos/ref/alt. The clinical-significance field is collapsed by the
    project's own classify_clnsig, so `classified` means P/LP or B/LB and
    `recorded_unclassified` means VUS, Conflicting or Other.
    """
    classify = _load_clnsig_classifier()
    wanted_chroms = {k[0] for k in keys}
    rows, file_date = [], None
    with gzip.open(vcf, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("##"):
                if line.startswith("##fileDate"):
                    file_date = line.strip().split("=", 1)[1]
                continue
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t", 8)
            if p[0] not in wanted_chroms:
                continue
            key = (p[0], int(p[1]), p[3], p[4])
            if key not in keys:
                continue
            info = dict(kv.split("=", 1) for kv in p[7].split(";") if "=" in kv)
            rows.append({
                "chrom": p[0], "pos": int(p[1]), "ref": p[3], "alt": p[4],
                "variation_id": p[2],
                "clnsig_raw": info.get("CLNSIG", ""),
                "clnsig_class": classify(info.get("CLNSIG", "")),
                "clnrevstat": info.get("CLNREVSTAT", ""),
            })
    out = pd.DataFrame(rows)
    out.attrs["file_date"] = file_date
    return out


def clinvar_arm(cls: str | float) -> str:
    if not isinstance(cls, str) or cls == "":
        return "unrecorded"
    return "classified" if cls in ("P/LP", "B/LB") else "recorded_unclassified"


def clinvar_arm_strict(cls: str | float, clnsig_raw: str | float) -> str:
    """The same three arms, except that a record carrying NO clinical assertion at
    all is separated out.

    classify_clnsig maps an empty CLNSIG to "Other", which clinvar_arm then files
    under recorded_unclassified -- so an allele-only ClinVar entry, which nobody has
    ever assessed, sits in the same arm as a variant ClinVar looked at and could not
    call. On this set 588 rows are such records and 576 of them are BRCA1, which is
    85% of BRCA1's middle arm; separating them moves that arm's likelihood ratio by
    a full evidence tier. The three arms the brief specifies stay as they are and
    this runs beside them.
    """
    if not isinstance(cls, str) or cls == "":
        return "unrecorded"
    if cls in ("P/LP", "B/LB"):
        return "classified"
    if not isinstance(clnsig_raw, str) or clnsig_raw.strip() == "":
        return "recorded_no_assertion"
    return "recorded_unclassified"


# ---------------------------------------------------------------------------
# build
# ---------------------------------------------------------------------------
def load_atlas() -> pd.DataFrame:
    ev = _load_atlas_module("atlas_evaluate", "src/atlas/evaluate.py")
    df = pd.read_parquet(ATLAS_MATRIX)
    df["region_atlas"] = df["hgvs_c"].map(ev.classify_region)
    df["intron_offset_abs"] = df["hgvs_c"].map(ev.splice_offset)
    return df


def attach_labels(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Official per-assay labels, joined on normalised transcript c. notation.

    The join key is the same one phase1_method_b_labels uses. It is valid on the
    atlas matrix because the atlas's MANE-Select c. strings agree with the frozen
    matrix's hgvs_nt wherever the two products overlap. The manifest records the
    agreement and the size of that overlap: it is the analysis set's intersection
    with the frozen matrix, not all 21,394 rows of the frozen matrix.
    """
    df = df.copy()
    df["_join"] = df["hgvs_c"].map(_norm_hgvs)
    y = pd.Series(pd.NA, index=df.index, dtype="Float64")
    src = pd.Series("unlabelled", index=df.index, dtype=object)

    cov = []
    for gene, spec in C.ASSAY_LABELS.items():
        path = C.ASSAY_LABEL_DIR / spec["file"]
        mask = df["gene"] == gene
        if not path.exists():
            cov.append(dict(gene=gene, label_file="MISSING", n=int(mask.sum()),
                            n_labelled=0, coverage=np.nan))
            continue
        lab = pd.read_csv(path, sep=None, engine="python")
        lut = dict(zip(lab[spec["key"]].map(_norm_hgvs),
                       lab[spec["class_col"]].map(spec["map"])))
        mapped = df.loc[mask, "_join"].map(lut)
        y.loc[mask] = mapped.astype("Float64")
        # a row present in the label file but mapped to None is indeterminate --
        # it stays in the set and out of the labelled analyses, per the brief
        in_file = df.loc[mask, "_join"].isin(lut)
        src.loc[mask] = np.where(mapped.notna(), "official",
                                 np.where(in_file, "indeterminate", "unlabelled"))
        cov.append(dict(gene=gene, label_file=spec["file"], n=int(mask.sum()),
                        n_labelled=int(mapped.notna().sum()),
                        coverage=round(float(mapped.notna().mean()), 4)))
    df["y_assay"] = y
    df["y_assay_source"] = src
    return df.drop(columns=["_join"]), pd.DataFrame(cov)


def directionality(atlas: pd.DataFrame) -> pd.DataFrame:
    """Control-anchored orientation gate, per gene, on the whole atlas rows for
    that gene -- the controls are coding, so they are not in the splice set."""
    from .phase1_directionality_check import check
    frame = atlas.rename(columns={"consequence": "region",
                                  "functional_pathogenicity": "func_pathogenicity"})
    return check(frame[["gene", "region", "func_pathogenicity"]])


def orientation_check(df: pd.DataFrame) -> pd.DataFrame:
    """Sign of each panel column against functional pathogenicity, measured.

    The atlas already orients gpn_msa and the allele-frequency columns; this
    records the measured sign rather than trusting either repo's convention.
    """
    rows = []
    extra = [c for c in list(OPTIONAL_COLUMNS) + ["alphagenome_v061"]
             if c in df.columns]
    for name in list(PANEL) + extra:
        s = df[[name, "func_pathogenicity"]].dropna()
        rho = (float(s[name].corr(s["func_pathogenicity"], method="spearman"))
               if len(s) > 30 else np.nan)
        rows.append({"feature": name, "n": len(s), "spearman_vs_pathogenicity": rho,
                     "oriented_larger_is_damaging": bool(rho > 0) if np.isfinite(rho) else None})
    return pd.DataFrame(rows)


def old_set_disposition(new: pd.DataFrame) -> pd.DataFrame:
    """Where each of the published study's 1,781 splice variants went."""
    fm = pd.read_parquet(C.OUTPUT_DIR / "frozen_matrix_v2.parquet")
    old = fm[fm["is_splice"]].copy()
    atlas_keys = set(new["key"])
    ev = _load_atlas_module("atlas_evaluate", "src/atlas/evaluate.py")

    atlas_all = pd.read_parquet(ATLAS_MATRIX, columns=["chrom", "pos", "ref", "alt", "hgvs_c"])
    atlas_all["key"] = (atlas_all.chrom.astype(str) + "-" + atlas_all.pos.astype(str)
                        + "-" + atlas_all.ref + "-" + atlas_all.alt)
    atlas_region = dict(zip(atlas_all["key"], atlas_all["hgvs_c"].map(ev.classify_region)))

    def disposition(row):
        k = row["variant_id"]
        if k in atlas_keys:
            return "in_new_set", ""
        if k not in atlas_region:
            return "absent_from_atlas", ("not in the atlas matrix: RAD51C rows the "
                                         "atlas deposit does not carry")
        r = atlas_region[k]
        if r == "coding_or_utr":
            return "reclassified_coding_or_utr", (
                "no intron offset in the atlas c. notation -- exon-side last-base "
                "donor variants and the 3 UTR-intron variants at offset -3; the "
                "atlas classifier files them as coding/UTR")
        if r == "splice_deep":
            return "outside_window", "|offset| > 50"
        return "other", r

    disp = old["variant_id"].to_frame()
    disp[["disposition", "reason"]] = old.apply(
        lambda r: pd.Series(disposition(r)), axis=1)
    disp["gene"] = old["gene"].to_numpy()
    disp["old_splice_subclass"] = old["splice_subclass"].to_numpy()
    return disp


def atlas_absent_rows() -> pd.DataFrame:
    """The frozen-matrix rows the atlas does not carry, and where they fall.

    The published study named sixteen RAD51C SNVs absent from the atlas deposit.
    They are accounted for here whether or not they are splice variants, because a
    reader checking the two products against each other will look for them.
    """
    fm = pd.read_parquet(C.OUTPUT_DIR / "frozen_matrix_v2.parquet")
    atlas_all = pd.read_parquet(ATLAS_MATRIX, columns=["chrom", "pos", "ref", "alt"])
    keys = set(atlas_all.chrom.astype(str) + "-" + atlas_all.pos.astype(str)
               + "-" + atlas_all.ref + "-" + atlas_all.alt)
    miss = fm[~fm["variant_id"].isin(keys)].copy()
    return miss[["variant_id", "gene", "region_label", "intron_offset",
                 "is_splice", "splice_subclass"]]


def build() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    atlas = load_atlas()
    is_snv = (atlas["ref"].str.len() == 1) & (atlas["alt"].str.len() == 1)
    sel = atlas[is_snv & atlas["region_atlas"].isin(STRATA)].copy()
    sel["stratum"] = sel["region_atlas"].map(STRATA)
    sel["key"] = (sel.chrom.astype(str) + "-" + sel.pos.astype(str)
                  + "-" + sel.ref + "-" + sel.alt)

    # panel
    keep = {v: k for k, v in PANEL.items()}
    sel = sel.rename(columns=keep)
    sel["func_pathogenicity"] = sel["functional_pathogenicity"]

    # AlphaGenome v0.6.1, the published study's definition, as a sensitivity column
    if ATLAS_V061.exists():
        v061 = pd.read_parquet(ATLAS_V061)
        # Named explicitly. This file is a full matrix, not a two-column score file,
        # so "the first column that is not variant_id" picks `gene` and the column
        # silently becomes a gene symbol -- which is what it was until this was
        # caught, and nothing downstream noticed because the v0.6.1 column is a
        # sensitivity column outside the panel and outside the orientation check.
        src_col = "alphagenome_splice"
        if src_col not in v061.columns:
            raise SystemExit(f"[E1] {ATLAS_V061} has no {src_col!r}; "
                             f"columns are {list(v061.columns)}")
        sel = sel.merge(v061[["variant_id", src_col]].rename(
            columns={src_col: "alphagenome_v061"}), on="variant_id", how="left")

    for name, path in OPTIONAL_COLUMNS.items():
        if path.exists():
            extra = pd.read_parquet(path)
            cols = [c for c in extra.columns
                    if c == name or c.startswith(name + "_")]
            sel = sel.merge(extra[["variant_id"] + cols], on="variant_id", how="left")

    sel, label_cov = attach_labels(sel)

    # ClinVar, three arms
    keys = set(zip(sel.chrom.astype(str), sel.pos.astype(int), sel.ref, sel.alt))
    cv = annotate_clinvar(keys, CLINVAR_VCF)
    file_date = cv.attrs.get("file_date")
    sel = sel.merge(cv, on=["chrom", "pos", "ref", "alt"], how="left")
    sel["clinvar_arm"] = sel["clnsig_class"].map(clinvar_arm)
    sel["clinvar_arm_strict"] = [
        clinvar_arm_strict(c, r) for c, r in zip(sel["clnsig_class"], sel["clnsig_raw"])]
    sel["clinvar_has_assertion"] = sel["clinvar_arm_strict"] != "recorded_no_assertion"

    dirn = directionality(atlas)
    orient = orientation_check(sel)

    # ---- reports -----------------------------------------------------------
    counts = (sel.assign(labelled=sel.y_assay.notna(),
                         pos=(sel.y_assay == 1))
                 .groupby(["gene", "stratum", "clinvar_arm", "clinvar_arm_strict"],
                          observed=True)
                 .agg(n=("key", "size"), n_labelled=("labelled", "sum"),
                      n_damaging=("pos", "sum"))
                 .reset_index())
    counts["n_normal"] = counts["n_labelled"] - counts["n_damaging"]
    counts["positive_rate"] = (counts["n_damaging"] / counts["n_labelled"]).round(4)
    counts.to_csv(REPORT_DIR / "set_counts.csv", index=False)

    disp = old_set_disposition(sel)
    disp.to_csv(REPORT_DIR / "old_set_disposition.csv", index=False)
    disp_summary = disp.groupby(["disposition", "reason"], dropna=False).size() \
                       .reset_index(name="n")
    disp_summary.to_csv(REPORT_DIR / "old_set_disposition_summary.csv", index=False)

    absent = atlas_absent_rows()
    absent.to_csv(REPORT_DIR / "frozen_rows_absent_from_atlas.csv", index=False)

    label_cov.to_csv(REPORT_DIR / "label_coverage.csv", index=False)
    dirn.to_csv(REPORT_DIR / "directionality_check.csv", index=False)
    orient.to_csv(REPORT_DIR / "feature_orientation.csv", index=False)

    drop = ["region_atlas", "mapping_status"]
    sel = sel.drop(columns=[c for c in drop if c in sel.columns])
    sel.to_parquet(OUT_PARQUET, index=False)

    # the scoring input for the model re-scores, written here so a clean clone has
    # it without an ad-hoc step
    sel[["variant_id", "gene", "chrom", "pos", "ref", "alt"]].to_parquet(
        OUT_DIR / "analysis_set_variants.parquet", index=False)

    # ---- manifest ----------------------------------------------------------
    fm = pd.read_parquet(C.OUTPUT_DIR / "frozen_matrix_v2.parquet",
                         columns=["variant_id", "hgvs_nt"])
    shared = sel.merge(fm, left_on="key", right_on="variant_id", suffixes=("", "_fm"))
    hgvs_agree = float((shared["hgvs_c"].map(_norm_hgvs)
                        == shared["hgvs_nt"].map(_norm_hgvs)).mean()) if len(shared) else np.nan

    manifest = {
        "product": "evidence analysis_set_v1",
        "built_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "n_rows": int(len(sel)),
        "sha256": sha256_of(OUT_PARQUET),
        "selection": {
            "variant_type": "SNV",
            "region": "intron-side splice, 1 <= |intron_offset| <= 50",
            "strata": STRATA,
            "requires_clinvar_record": False,
            "requires_complete_predictor_row": False,
        },
        "sources": {
            "atlas_matrix": {
                "path": str(ATLAS_MATRIX),
                "sha256": sha256_of(ATLAS_MATRIX),
                "atlas_git_head": git_head(ATLAS_REPO),
            },
            "clinvar_vcf": {
                "path": str(CLINVAR_VCF),
                "sha256": sha256_of(CLINVAR_VCF),
                "file_date": file_date,
                "release": "GRCh38 weekly, 15 June 2026",
                "review_status_filter": "none -- records of every review status kept",
            },
            "assay_labels": {
                gene: {"file": spec["file"],
                       "sha256": sha256_of(C.ASSAY_LABEL_DIR / spec["file"])}
                for gene, spec in C.ASSAY_LABELS.items()
                if (C.ASSAY_LABEL_DIR / spec["file"]).exists()
            },
            "frozen_matrix_v2": {
                "path": str(C.OUTPUT_DIR / "frozen_matrix_v2.parquet"),
                "sha256": sha256_of(C.OUTPUT_DIR / "frozen_matrix_v2.parquet"),
            },
        },
        "checks": {
            "hgvs_agreement_with_frozen_v2_on_shared": hgvs_agree,
            "n_shared_with_frozen_v2": int(len(shared)),
            "old_splice_set_total": int(len(disp)),
            "old_splice_set_in_new": int((disp.disposition == "in_new_set").sum()),
            "old_splice_set_explained": int((disp.disposition != "in_new_set").sum()),
        },
        "columns": list(sel.columns),
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")

    # ---- console -----------------------------------------------------------
    print(f"[E1] analysis_set_v1: {len(sel):,} SNVs  sha256 {manifest['sha256'][:16]}...")
    print(sel.groupby(["stratum", "clinvar_arm"], observed=True).size()
             .unstack(fill_value=0).reindex(STRATUM_ORDER).to_string())
    print("\n[E1] labelled:", int(sel.y_assay.notna().sum()),
          " damaging:", int((sel.y_assay == 1).sum()),
          " normal:", int((sel.y_assay == 0).sum()),
          " indeterminate:", int((sel.y_assay_source == "indeterminate").sum()))
    print("\n[E1] old-set disposition (must sum to 1,781):")
    print(disp_summary.to_string(index=False))
    print(f"\n[E1] frozen-matrix rows absent from the atlas: {len(absent)}"
          f"  (splice: {int(absent.is_splice.sum())})")
    if len(absent):
        print(absent.groupby(["gene", "region_label", "is_splice"], observed=True)
                    .size().reset_index(name="n").to_string(index=False))
    print(f"\n[E1] ClinVar fileDate: {file_date}")
    strict = sel["clinvar_arm_strict"].value_counts().to_dict()
    print(f"[E1] arms with no-assertion records separated: {strict}")
    noass = sel[sel["clinvar_arm_strict"] == "recorded_no_assertion"]
    if len(noass):
        print(f"     of the {len(noass)} records carrying no clinical assertion, "
              f"by gene: {noass.gene.value_counts().to_dict()}")
    print("\n[E1] directionality gate:")
    print(dirn.to_string(index=False))
    print("\n[E1] measured feature orientation:")
    print(orient.to_string(index=False))


if __name__ == "__main__":
    build()
