"""E7 -- external genes: do the fixed and the fitted thresholds transfer?

Two genes, and they answer different questions.

TP53 is the gene the published study already held out. Its archived table stopped
at |offset| 8 because of that study's splice-window constant, although the deposit
reaches 12; evid_tp53_extend carries the 192 archived rows unchanged and scores the
96 deposited SNVs at offsets 9-12 with the same scorers, so TP53 now populates pm12,
s3_10 and eleven to twelve nucleotides of s11_50. The deposit publishes no
per-variant classification, so its labels are constructed (below). Its 11-12 nt
rows hold no damaging variant under the control-anchored label, so that band is
reported by count only.

DDX3X is new, and is outside the gene class everything else here sits in: it is
X-linked, and its disease association is neurodevelopmental rather than cancer
predisposition. Its assay is NOT a different kind of readout, and an earlier version
of this docstring implied it was: like all seven, the DDX3X deposit scores
variants by depletion from HAP1 cells over a time course, where the gene is
essential (Radford et al. 2023). What differs is the gene, not the measurement. Its saturation
editing deposit (MaveDB urn:mavedb:00000658, seventeen exon-level score sets, CC0)
carries 1,857 intron-side SNVs at 1 <= |offset| <= 50 -- more splice-region variants
than the published study's whole analysis set -- with 192 / 768 / 897 across the
three strata. That is what makes it worth the work: it is the only candidate found
that populates the distal band, to 25 nt, the territory the companion atlas locates the real failure
in, on a gene set that shares nothing with the seven.

Labels follow the rule used for the seven genes: the assay authors' own
per-variant classification. The DDX3X deposit publishes one in a separate score
set (urn:mavedb:00000658-0-1): a random-forest call whose only inputs are the
assay's combined log fold-changes at days 7, 11 and 15, with a decision boundary
trained on clinically classified variants for the neurodevelopmental context,
'abnormal' above a posterior of 0.5 (Radford et al. 2023). An earlier version of
this module said the deposit published no classification; that was wrong, and the
constructions it used instead are kept below as sensitivity labels:

  deposit                   : the published call (primary)
  control_anchored_gated    : damaging if the score falls below the midpoint of
                              the score set's synonymous and nonsense control
                              medians, in score sets whose controls separate
                              (control AUROC >= 0.80); unlabelled elsewhere
  control_anchored          : the same without the gate
  fdr                       : damaging if the trend FDR clears 0.01 and the score
                              sits below the synonymous median

The gate exists because one score set, urn:mavedb:00000658-q-1 (the final-exon
tile), has nonsense controls that are not depleted -- control AUROC 0.51 against
0.99 to 1.0 in fifteen of the other sixteen -- so its midpoint carries no
information and the ungated label is close to a coin flip there. Without the gate
that one tile supplies over half of the in-scope control-anchored damaging calls.
Its behaviour is consistent with nonsense-mediated decay escape in the last exon.
The ungated label is kept so the effect of the gate is visible.

Since E18.3 every column E7 tests on either gene is the analysis set's own variable:
the same scorer, version, precision and direction, shown by anchors scored in the
same run or by a committed value-for-value record (docs/external-column-identity.md).
TP53's AlphaGenome is the official definition, its client-0.6.1 column kept beside it
as alphagenome_v061 and not tested; its SpliceAI, Pangolin, NT, Walker-basis and
Atlas columns were scored again by the analysis set's scorers; its GPN-MSA is turned
to the analysis set's direction. DDX3X gains SpliceAI at distance 50, AlphaGenome and
NT, and its Pangolin is full precision.

Run (PYTHONPATH=phase1):
    python -m src.evid_external --prepare ddx3x
    python -m src.evid_external --apply TP53
    python -m src.evid_external --apply DDX3X
"""
from __future__ import annotations

import argparse
import csv
import glob
import hashlib
import json
import os
import re
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from scipy.stats import spearmanr

from . import config as C
from . import evid_common as K

from .evid_common import ATLAS_REPO as _atlas_repo_default  # noqa: E402
ATLAS_REPO = _atlas_repo_default  # EVID_ATLAS_REPO overrides; see evid_common
# Transcript models and reference sequence the mapping reads, tracked in this
# repository (a copy of the atlas cache, including the models fetched after the
# atlas archive was made), so these stages need neither that cache nor the
# Mutalyzer and Ensembl services.
REF_CACHE = Path("data/evidence/reference_cache")
EXT_DIR = Path("data/evidence/external")
REPORT_DIR = Path("reports/evidence")
CONFIG_PATH = Path("config/walker2023.yaml")

DDX3X_EXPERIMENT = "urn:mavedb:00000658"
DDX3X_SCORESETS = [f"{DDX3X_EXPERIMENT}-{c}-1" for c in "abcdefghijklmnopq"]
DDX3X_TRANSCRIPT = "ENST00000644876.2"
# the deposit's own per-variant classification, published as its own score set
DDX3X_CLASSIFICATION_SCORESET = f"{DDX3X_EXPERIMENT}-0-1"
DDX3X_CLASS_COLUMN = "SGE_prediction_of_variant_function_in_NDD_context"
DDX3X_CLASS_MAP = {"abnormal": 1.0, "normal": 0.0}
# per-score-set control separation below which the control-anchored label is not
# assigned in that set; the same 0.80 the pooled orientation gate uses
CONTROL_GATE = 0.80
# primary first; the rest are sensitivity labels
DDX3X_LABELS = ["y_deposit", "y_control_anchored_gated", "y_control_anchored", "y_fdr"]
N_BOOT = 2000
MAVEDB_API = "https://api.mavedb.org/api/v1/score-sets"

def intronic_snv_regex(transcript: str = r"ENST[\d.]+") -> re.Pattern:
    """A transcript-level HGVS of an intron-side SNV. Groups: the transcript, the c.
    string, the distance into the intron. DDX3X's hgvs_splice column names Ensembl
    transcripts; evid_newgenes reads RefSeq ones from hgvs_nt with the same pattern."""
    return re.compile(rf"({transcript}):(c\.[-*]?\d+[+-](\d+)[ACGT]>[ACGT])$")


_SPLICE = intronic_snv_regex()
_PRO_SYN = re.compile(r"p\.[A-Za-z]{3}\d+=")
_PRO_NON = re.compile(r"p\.[A-Za-z]{3}\d+(Ter|\*)")

FDR_CUT = 0.01
# TP53 to |offset| 12: the archived 192 rows plus the 96 scored by evid_tp53_extend
TP53_TABLE = EXT_DIR / "tp53_splice_scored_12nt.parquet"
# The names every output gives the external genes' primary labels, from this one place:
# each is named for its source as the source's record calls it, with the URN of the
# score set it comes from (the author's decision of 2026-10-24). A label read from a
# score set's investigator-provided ranges is an "Investigator-provided range"; one
# copied from a deposit's per-variant class column a "Deposit classification". The
# figures write the short forms, and the URN goes into the caption.
#   TP53 (E18.4): its construction in _load_external (relative fitness score > 0) is,
# variant by variant, the classification by the investigator-provided ranges of
# urn:mavedb:00001277-a-1, the calibration its record titles "Investigator-provided
# functional ranges" (investigatorProvided: true; the record's two other calibrations,
# IGVF's and ExCALIBR's, are not used): all 288 rows E7 reads, at their own coordinates
# (scripts/stage_E18/check_tp53_ranges.py; tp53_raw/). Its column keeps the name
# y_control_anchored.
#   DDX3X: the deposit's own call, the class column of the classification score set.
LABEL_SOURCE_NAME = {"ranges": "Investigator-provided range",
                     "classification": "Deposit classification"}
PRIMARY_LABEL_SOURCE = {"DDX3X": ("classification", DDX3X_CLASSIFICATION_SCORESET),
                        "TP53": ("ranges", "urn:mavedb:00001277-a-1")}


def primary_label_name(gene: str) -> str:
    """The tables' name for a gene's primary label: its source and the score set's URN."""
    kind, urn = PRIMARY_LABEL_SOURCE[gene]
    return f"{LABEL_SOURCE_NAME[kind]} ({urn})"


def primary_label_short(gene: str) -> str:
    """The figures' name for it, without the URN."""
    return LABEL_SOURCE_NAME[PRIMARY_LABEL_SOURCE[gene][0]].lower()


TP53_PRIMARY_LABEL_NAME = primary_label_name("TP53")
TP53_PRIMARY_LABEL_SHORT = primary_label_short("TP53")
DDX3X_PRIMARY_LABEL_NAME = primary_label_name("DDX3X")
DDX3X_PRIMARY_LABEL_SHORT = primary_label_short("DDX3X")


def stratum_of(offset: int) -> str:
    if offset <= 2:
        return "pm12"
    if offset <= 10:
        return "s3_10"
    return "s11_50"


# ---------------------------------------------------------------------------
# DDX3X preparation
# ---------------------------------------------------------------------------
def _fetch(url: str) -> bytes:
    with urllib.request.urlopen(url, timeout=120) as r:
        return r.read()


def prepare_ddx3x() -> None:
    sys.path.insert(0, str(ATLAS_REPO / "src"))
    from atlas import mapping as M
    # atlas.mapping resolves its reference cache RELATIVE to the working directory
    # ("data/raw/reference"). These stages run from phase1/, so every cached
    # Mutalyzer model and Ensembl lookup would be missed and re-fetched -- slowly,
    # and with a second copy of the cache appearing under phase1/. Point it at the
    # atlas's own cache.
    M.REF_CACHE = (REF_CACHE if REF_CACHE.exists()
                   else ATLAS_REPO / "data/raw/reference")

    raw_dir = EXT_DIR / "ddx3x_raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    tables = {}
    for urn in DDX3X_SCORESETS:
        path = raw_dir / f"{urn.replace(':', '_')}.csv"
        if not path.exists():
            path.write_bytes(_fetch(f"{MAVEDB_API}/{urn}/scores"))
        tables[urn] = list(csv.DictReader(path.open()))

    mane = M.load_mane_records(["DDX3X"])["DDX3X"]
    tmap, prov = M.build_transcript_map(DDX3X_TRANSCRIPT, mane)
    lo = min(e.g_start for e in tmap.exons) - 5000
    hi = max(e.g_end for e in tmap.exons) + 5000
    # MANE reports DDX3X's chromosome as '23'; the Ensembl sequence endpoint
    # names it 'X'. Every other gene in this project is autosomal, which is why
    # the atlas mapping code never had to make the distinction.
    s0, _, seq = M.fetch_ensembl_region("X", lo, hi)
    rseq = M.RegionSeq(s0, seq)

    rows, controls, failures = [], [], {}
    for urn, table in tables.items():
        for r in table:
            score = r.get("score")
            if score in (None, "", "NA"):
                continue
            pro = (r.get("hgvs_pro") or "").strip()
            if _PRO_SYN.search(pro):
                controls.append({"scoreset": urn, "class": "synonymous",
                                 "score": float(score)})
            elif _PRO_NON.search(pro):
                controls.append({"scoreset": urn, "class": "nonsense",
                                 "score": float(score)})
            m = _SPLICE.search((r.get("hgvs_splice") or "").strip())
            if not m:
                continue
            try:
                pos, ref, alt = M.map_variant(tmap, rseq, m.group(2))
            except Exception as e:                       # noqa: BLE001
                failures[str(e)[:80]] = failures.get(str(e)[:80], 0) + 1
                continue
            fdr = r.get("cLFC_trend_BH_FDR")
            rows.append({
                "variant_id": f"X-{pos}-{ref}-{alt}", "gene": "DDX3X",
                "scoreset": urn, "mavedb_accession": r["accession"],
                "hgvs_c": m.group(2), "chrom": "X", "pos": pos, "ref": ref, "alt": alt,
                "intron_offset_abs": int(m.group(3)),
                "stratum": stratum_of(int(m.group(3))),
                "assay_score": float(score),
                "fdr": float(fdr) if fdr not in (None, "", "NA") else np.nan,
            })

    df = pd.DataFrame(rows).drop_duplicates("variant_id").reset_index(drop=True)
    ctl = pd.DataFrame(controls)
    med_syn = float(ctl.loc[ctl["class"] == "synonymous", "score"].median())
    med_non = float(ctl.loc[ctl["class"] == "nonsense", "score"].median())

    # orientation gate, control-anchored, exactly as phase1_directionality_check
    from sklearn.metrics import roc_auc_score
    y = (ctl["class"] == "nonsense").astype(int).to_numpy()
    # damage is depletion here, so pathogenicity is the negated score
    auroc = float(roc_auc_score(y, -ctl["score"].to_numpy()))
    df["func_pathogenicity"] = -df["assay_score"]

    # The seventeen score sets are separate experiments on separate exon tiles and
    # they are NOT on a common scale: the per-set nonsense median runs from -0.0009
    # to -0.2700. A single pooled midpoint therefore sits outside several sets'
    # whole score range, and those sets contribute no damaging calls at all -- one
    # of them labels a canonical donor variant with a trend FDR of 1e-26 as normal.
    # The control anchoring is done per score set; the pooled values are kept only
    # for the record.
    per_set = (ctl.pivot_table(index="scoreset", columns="class", values="score",
                               aggfunc="median")
                  .rename(columns={"synonymous": "med_syn", "nonsense": "med_non"}))
    per_set["midpoint"] = (per_set["med_syn"] + per_set["med_non"]) / 2.0
    midpoint = (med_syn + med_non) / 2.0          # pooled, reported not used

    mid = df["scoreset"].map(per_set["midpoint"])
    syn = df["scoreset"].map(per_set["med_syn"])
    missing = mid.isna()
    if missing.any():                              # a set with no controls falls back
        mid = mid.fillna(midpoint)
        syn = syn.fillna(med_syn)
    df["control_midpoint_used"] = mid
    df["y_control_anchored"] = (df["assay_score"] < mid).astype(float)
    df["y_fdr"] = ((df["fdr"] < FDR_CUT) & (df["assay_score"] < syn)).astype(float)

    # The pooled gate above can pass while one score set's controls do not separate
    # at all, so separation is also checked per set, and the control-anchored label
    # is withheld in a set that fails it.
    set_auroc = {}
    for urn, g in ctl.groupby("scoreset"):
        yy = (g["class"] == "nonsense").astype(int).to_numpy()
        set_auroc[urn] = (float(roc_auc_score(yy, -g["score"].to_numpy()))
                          if 0 < yy.sum() < len(yy) else np.nan)
    per_set["control_auroc"] = pd.Series(set_auroc)
    # how much of the ungated label each set supplies, in the in-scope window
    ins = df[df["stratum"] != "pm12"]
    per_set["n_in_scope"] = ins.groupby("scoreset").size()
    per_set["n_damaging_control_anchored_in_scope"] = (
        ins.groupby("scoreset")["y_control_anchored"].sum())
    per_set = per_set.fillna({"n_in_scope": 0, "n_damaging_control_anchored_in_scope": 0})
    gated = sorted(u for u, a in set_auroc.items() if not (a >= CONTROL_GATE))
    df["y_control_anchored_gated"] = np.where(df["scoreset"].isin(gated), np.nan,
                                              df["y_control_anchored"])

    # The deposit's own call, joined on the transcript-level HGVS.
    cpath = raw_dir / f"{DDX3X_CLASSIFICATION_SCORESET.replace(':', '_')}.csv"
    if not cpath.exists():
        cpath.write_bytes(_fetch(f"{MAVEDB_API}/{DDX3X_CLASSIFICATION_SCORESET}/scores"))
    cls = pd.read_csv(cpath, usecols=["hgvs_nt", "score", DDX3X_CLASS_COLUMN])
    cls = cls[cls["hgvs_nt"].str.startswith(DDX3X_TRANSCRIPT + ":", na=False)]
    cls["hgvs_c"] = cls["hgvs_nt"].str.split(":", n=1).str[1]
    if cls["hgvs_c"].duplicated().any():
        raise SystemExit("[E7] DDX3X classification: duplicated HGVS in the deposit")
    cls = cls.set_index("hgvs_c")
    df["deposit_class"] = df["hgvs_c"].map(cls[DDX3X_CLASS_COLUMN])
    df["deposit_posterior"] = df["hgvs_c"].map(cls["score"])
    df["y_deposit"] = df["deposit_class"].map(DDX3X_CLASS_MAP)
    n_unmatched = int(df["deposit_class"].isna().sum())

    EXT_DIR.mkdir(parents=True, exist_ok=True)
    out = EXT_DIR / "ddx3x_splice.parquet"
    df.to_parquet(out, index=False)
    manifest = {
        "gene": "DDX3X", "source": "MaveDB", "experiment": DDX3X_EXPERIMENT,
        "score_sets": DDX3X_SCORESETS, "licence": "CC0",
        "transcript": DDX3X_TRANSCRIPT, "mane_select": True,
        "assembly": "GRCh38", "chrom": "X",
        "n_intronic_snv_mapped": int(len(df)),
        "n_mapping_failures": int(sum(failures.values())),
        "mapping_failures": failures,
        "strata": df["stratum"].value_counts().to_dict(),
        "n_in_scope_3_to_max": int((df["stratum"] != "pm12").sum()),
        "max_intron_offset": int(df["intron_offset_abs"].max()),
        "controls": {"n_synonymous": int((ctl["class"] == "synonymous").sum()),
                     "n_nonsense": int((ctl["class"] == "nonsense").sum()),
                     "median_synonymous_pooled": med_syn,
                     "median_nonsense_pooled": med_non,
                     "midpoint_pooled_not_used": midpoint,
                     "anchoring": "per score set; the seventeen tiles are separate "
                                  "experiments and are not on a common scale",
                     "per_scoreset": per_set.round(6).to_dict(orient="index"),
                     "n_scoresets_without_controls": int(missing.sum()),
                     "control_auroc": auroc,
                     "gate": "OK" if auroc >= 0.80 else
                             ("WEAK_SEPARATION" if auroc >= 0.65 else "FLIP_NEEDED"),
                     "per_scoreset_gate": CONTROL_GATE,
                     "scoresets_gated_out": gated},
        "labels": {
            "primary": "deposit",
            "deposit": (f"{DDX3X_CLASSIFICATION_SCORESET} column {DDX3X_CLASS_COLUMN}: "
                        "random-forest call on the assay's combined log fold-changes, "
                        "decision boundary trained on clinically classified variants; "
                        "abnormal -> 1, normal -> 0"),
            "control_anchored_gated": ("assay_score < midpoint of the score set's "
                                       "control medians, in score sets with control "
                                       f"AUROC >= {CONTROL_GATE}; unlabelled elsewhere"),
            "control_anchored": "assay_score < midpoint of control medians",
            "fdr": f"trend BH-FDR < {FDR_CUT} and assay_score < synonymous median",
            "note": "the deposit label is the assay authors' own call, as for the "
                    "seven genes; the other three are constructed here and reported "
                    "as sensitivity labels",
            "n_unmatched_to_deposit_call": n_unmatched,
            "n_damaging_deposit": int(df["y_deposit"].sum()),
            "n_damaging_control_anchored_gated": int(df["y_control_anchored_gated"].sum()),
            "n_damaging_control_anchored": int(df["y_control_anchored"].sum()),
            "n_damaging_fdr": int(df["y_fdr"].sum()),
        },
        "transcript_map_provenance": prov,
        "sha256": hashlib.sha256(out.read_bytes()).hexdigest(),
        "built_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    (EXT_DIR / "ddx3x_splice.manifest.json").write_text(json.dumps(manifest, indent=2,
                                                                  default=str) + "\n")
    print(f"[E7] DDX3X: {len(df):,} intron-side SNVs mapped, "
          f"{sum(failures.values())} failures")
    print(f"     strata {manifest['strata']}")
    print(f"     control gate AUROC {auroc:.4f} ({manifest['controls']['gate']}); "
          f"pooled syn median {med_syn:.4f}, nonsense median {med_non:.4f}")
    print(f"     anchored per score set; per-set midpoints "
          f"{per_set['midpoint'].min():.4f} to {per_set['midpoint'].max():.4f}")
    lab = manifest["labels"]
    print(f"     damaging: deposit call {lab['n_damaging_deposit']} (primary; "
          f"{lab['n_unmatched_to_deposit_call']} unmatched), control-anchored gated "
          f"{lab['n_damaging_control_anchored_gated']}, ungated "
          f"{lab['n_damaging_control_anchored']}, FDR {lab['n_damaging_fdr']}")
    print(f"     score sets gated out of the control-anchored label: {gated or 'none'}")
    print(f"     wrote {out}")
    print("\n[E7] next: score the panel for these variants, then --apply DDX3X")


# ---------------------------------------------------------------------------
# apply the thresholds
# ---------------------------------------------------------------------------
# The external genes' columns (E18.3). Each is the analysis set's own variable: scored
# after the thresholds were frozen by the scorer of the analysis set's column, which
# in the same run reproduced the analysis set's stored values of 39 anchor variants
# bit for bit (scripts/stage_E18/build_columns.py), or shown so by a committed
# value-for-value record; docs/external-column-identity.md gives the evidence cell by
# cell. The files E7 read before (DDX3X's ddx3x_cadd.parquet and the rest, Pangolin at
# the CLI's two decimals among them) stay in place, unread.
DDX3X_SCORE_FILES = {
    "ddx3x_spliceai.parquet": ["spliceai"],
    "ddx3x_spliceai_walker_rescored.parquet": ["spliceai_walker"],
    "ddx3x_pangolin_fullprec.parquet": ["pangolin"],
    "ddx3x_cadd_rescored.parquet": ["cadd"],
    "ddx3x_conservation_rescored.parquet": ["phylop", "phastcons"],
    "ddx3x_gpn_msa_rescored.parquet": ["gpn_msa"],
    "ddx3x_nt_a6000.parquet": ["nt"],
}
# TP53: the 288-row table's columns that are not the analysis set's variables,
# replaced by E18.2's. Its CADD, phyloP, phastCons and GPN-MSA are kept: the published
# study's scorers give the analysis set's values (docs/column-provenance.md, section 2),
# GPN-MSA on the raw sign, which _load_external turns to the analysis set's (E18.0
# rule 5). Its AlphaGenome is the client-0.6.1 definition, kept as alphagenome_v061,
# which has no frozen threshold and is not tested.
TP53_SCORE_FILES = {
    "tp53_spliceai_rescored.parquet": ["spliceai"],
    "tp53_spliceai_walker_rescored.parquet": ["spliceai_walker", "spliceai_walker_ds_ag",
                                              "spliceai_walker_ds_al", "spliceai_walker_ds_dg",
                                              "spliceai_walker_ds_dl"],
    "tp53_pangolin_rescored.parquet": ["pangolin"],
    "tp53_nt_rescored.parquet": ["nt"],
    "tp53_avi_rescored.parquet": ["avi", "avi_splice_sites", "avi_splice_site_usage",
                                  "avi_splice_junctions"],
}
# AlphaGenome by its official definition (E17), scored for the external genes after
# the tag (E18.1) and read as the analysis set reads its own (evid_build_set:
# the round-trip parser)
AG_EXTERNAL = Path("data/evidence/alphagenome_official/external_scores.tsv")
# NT, scored again on the original GPU, did not reproduce its anchors, so by the rule the
# author fixed before the run (2026-10-23) no external gene has an NT column of the
# analysis set's variable and NT is not tested; TP53's published-study NT column goes
# with it. The column check names the column all the same.
NT_RECORD = EXT_DIR / "nt_gpu_run" / "record.json"


def _no_same_variable(tool: str) -> str | None:
    """Why an external gene has no column of the analysis set's variable, if recorded."""
    if tool != "nt" or not NT_RECORD.exists():
        return None
    r = json.loads(NT_RECORD.read_text())["anchors"]
    return (f"no column of the same variable: NT on the original GPU reproduced "
            f"{r['n_identical']} of the {r['n']} anchors bit for bit; by the rule fixed "
            "before the run (the author's decision of 2026-10-23) it does not enter the "
            "external test")


def _alphagenome_for(gene: str) -> pd.DataFrame:
    ag = pd.read_csv(AG_EXTERNAL, sep="\t", usecols=["gene", "variant_id", "alphagenome_splicing"],
                     float_precision="round_trip")
    ag = ag[ag["gene"] == gene]
    return ag[["variant_id", "alphagenome_splicing"]].rename(
        columns={"alphagenome_splicing": "alphagenome"})


# The AlphaGenome Atlas columns are not in the registry above because they do not
# arrive as a per-gene file: one Atlas pass covered the analysis set and both
# external genes, and wrote a single table carrying a `set` label. They are taken
# from it here, filtered to this gene's rows, so the external gene is scored on the
# same lookup the seven genes are.
AVI_PATH = Path("data/evidence/avi.parquet")
AVI_COLUMNS = ["avi", "avi_splice_sites", "avi_splice_site_usage",
               "avi_splice_junctions"]


def _avi_for(set_name: str) -> pd.DataFrame | None:
    if not AVI_PATH.exists():
        return None
    a = pd.read_parquet(AVI_PATH)
    if "set" not in a.columns:
        return None
    a = a[a["set"] == set_name]
    take = [c for c in AVI_COLUMNS if c in a.columns]
    return a[["variant_id"] + take] if len(a) and take else None


def merge_ddx3x_scores() -> None:
    """Collect whatever has been scored into one table the apply step reads."""
    base = pd.read_parquet(EXT_DIR / "ddx3x_splice.parquet")[["variant_id"]]
    present, missing = [], []
    for fname, cols in DDX3X_SCORE_FILES.items():
        path = EXT_DIR / fname
        if not path.exists():
            missing.append(fname)
            continue
        d = pd.read_parquet(path)
        take = [c for c in cols if c in d.columns]
        if not take:
            missing.append(f"{fname} (no expected column)")
            continue
        base = base.merge(d[["variant_id"] + take], on="variant_id", how="left")
        present += take
    base = base.merge(_alphagenome_for("DDX3X"), on="variant_id", how="left")
    present.append("alphagenome")
    avi = _avi_for("ddx3x")
    if avi is not None:
        base = base.merge(avi, on="variant_id", how="left")
        present += [c for c in avi.columns if c != "variant_id"]
    else:
        missing.append("avi.parquet (Atlas columns)")
    out = EXT_DIR / "ddx3x_scored.parquet"
    base.to_parquet(out, index=False)
    cov = {c: int(base[c].notna().sum()) for c in base.columns if c != "variant_id"}
    print(f"[E7] merged {len(cov)} scored columns for {len(base):,} DDX3X variants")
    for c, n in cov.items():
        print(f"     {c:20s} {n:>5}/{len(base)}")
    if missing:
        print("     not scored: " + ", ".join(missing))
    print(f"     wrote {out}")


def _load_external(gene: str) -> tuple[pd.DataFrame, list[str]]:
    if gene.upper() == "TP53":
        df = pd.read_parquet(TP53_TABLE)
        df = df[df["is_splice"]].copy()
        df["intron_offset_abs"] = df["intron_offset"].abs()
        df["stratum"] = df["intron_offset_abs"].map(
            lambda o: stratum_of(int(o)) if np.isfinite(o) else "pm12")
        # The deposit carries no per-variant classification, so the published study
        # constructed one. Its three definitions are reproduced here from
        # phase4_external_tp53.py rather than re-invented: relative fitness score is
        # anchored at -1 for synonymous and +1 for nonsense, so the control-anchored
        # midpoint is RFS > 0; the median split and the mid-band exclusion are the
        # two sensitivity definitions that study reports beside it.
        rfs = df["func_pathogenicity"].to_numpy(dtype=float)
        df["y_control_anchored"] = (rfs > 0).astype(float)
        df["y_median_split"] = (rfs > np.median(rfs)).astype(float)
        mid = np.abs(rfs) >= 0.5
        df["y_mid_band_excluded"] = np.where(mid, (rfs > 0).astype(float), np.nan)
        # E18.3: the analysis set's variables (TP53_SCORE_FILES)
        df = df.rename(columns={"alphagenome": "alphagenome_v061"})
        for fname, cols in TP53_SCORE_FILES.items():
            df = df.drop(columns=[c for c in cols if c in df.columns])
            if not (EXT_DIR / fname).exists():
                # NT only: by the author's decision of 2026-10-23 it is scored on the
                # original GPU and counts only if its anchors hold; without it TP53 has
                # no NT column of the analysis set's variable, so none is tested
                if cols != ["nt"]:
                    raise SystemExit(f"[E7] {EXT_DIR / fname} is missing")
                continue
            d = pd.read_parquet(EXT_DIR / fname)[["variant_id"] + cols]
            df = df.merge(d, on="variant_id", how="left")
        df = df.merge(_alphagenome_for("TP53"), on="variant_id", how="left")
        df["gpn_msa"] = -df["gpn_msa"]       # the raw sign, turned (E18.0 rule 5)
        return df, ["y_control_anchored", "y_median_split", "y_mid_band_excluded"]
    if gene.upper() == "DDX3X":
        df = pd.read_parquet(EXT_DIR / "ddx3x_splice.parquet")
        scored = EXT_DIR / "ddx3x_scored.parquet"
        if scored.exists():
            extra = pd.read_parquet(scored)
            df = df.merge(extra, on="variant_id", how="left", suffixes=("", "_dup"))
        return df, list(DDX3X_LABELS)
    raise SystemExit(f"[E7] no external gene '{gene}'")


# Definition mismatches that no statistical check can be relied on to find, because
# the two columns can overlap in range and still be different variables. Recorded
# per external gene and per column, with the source of the claim. Empty since E18.3:
# TP53's AlphaGenome entry went when E7 began reading the official definition, scored
# for TP53 after the tag (E18.1); the client-0.6.1 column it described is kept as
# alphagenome_v061, outside the test.
DECLARED_MISMATCH: dict[tuple[str, str], str] = {}


def column_basis_check(df: pd.DataFrame, tools: list[str],
                       gene_name: str = "") -> pd.DataFrame:
    """Is each external column the same variable as the analysis-set column of that
    name? A threshold fitted on one is meaningless on the other.

    Since E18.3 every column is the analysis set's variable by construction
    (TP53_SCORE_FILES, DDX3X_SCORE_FILES; docs/external-column-identity.md), and this
    check stays as the independent guard it was: a wrong orientation, a range that
    misses the fitted column, or a declared definition change still disqualifies.
    """
    base = K.load_set()
    rows = []
    for tool in tools:
        if tool not in base.columns:
            rows.append({"tool": tool, "comparable": False,
                         "reason": "no column of this name in the analysis set"})
            continue
        e = df[tool].dropna().to_numpy(dtype=float)
        a = base[tool].dropna().to_numpy(dtype=float)
        lab = df.get("func_pathogenicity")
        rho = np.nan
        if lab is not None:
            m = df[tool].notna() & lab.notna()
            if m.sum() > 30:
                rho = float(spearmanr(df.loc[m, tool], lab[m]).statistic)
        oriented = (not np.isfinite(rho)) or rho > 0
        # a range that sits largely outside the fitted column's range means a
        # threshold from that column lands somewhere else on this one
        overlap = (min(e.max(), a.max()) - max(e.min(), a.min())) / \
                  max(a.max() - a.min(), 1e-12)
        # Where the fitted column's 90th percentile -- roughly where a PP3 threshold
        # sits -- lands in this gene's distribution. This is REPORTED, not used to
        # disqualify: a threshold landing at a different quantile on a new gene is
        # the phenomenon E7 exists to measure, not a reason to refuse the test. Only
        # a different VARIABLE disqualifies, and that shows up as a wrong
        # orientation or as a declared definition change.
        q90 = float(np.quantile(a, 0.90))
        q90_here = float((e <= q90).mean())
        declared = DECLARED_MISMATCH.get((gene_name, tool))
        ok = bool(oriented and overlap > 0.5 and not declared)
        reason = ""
        if declared:
            reason = f"declared definition mismatch: {declared}"
        elif not oriented:
            reason = "anti-correlated on this gene: different orientation convention"
        elif overlap <= 0.5:
            reason = "range barely overlaps the fitted column"
        rows.append({
            "tool": tool, "comparable": ok,
            "external_min": float(e.min()), "external_max": float(e.max()),
            "analysis_set_min": float(a.min()), "analysis_set_max": float(a.max()),
            "spearman_vs_pathogenicity": rho, "range_overlap": float(overlap),
            "fitted_q90": q90, "fitted_q90_percentile_here": q90_here,
            "reason": reason,
        })
    return pd.DataFrame(rows)


def _rng_for(*key) -> np.random.Generator:
    """Seeded from the cell's identity, as evid_walker_thresholds does, so one row
    can be reproduced on its own and adding a row changes no other."""
    h = hashlib.sha256(("|".join(str(k) for k in key)).encode()).digest()
    return np.random.default_rng([C.RANDOM_SEED, int.from_bytes(h[:8], "big")])


def _band_detail(y: np.ndarray, v: np.ndarray, thr: float, side: str,
                 bands: dict, direction: str, cell: tuple) -> dict:
    """Counts in the band, the ratio, and a variant-level bootstrap interval.

    An external gene is one gene, so there is no cluster structure to resample and
    the gene-clustered interval of the seven-gene tables does not exist here. The
    interval below resamples variants within each class (2,000 draws, class sizes
    fixed): it describes how precisely this gene's own variants estimate the ratio,
    not how the ratio would vary between genes. A tier is reported both from the
    point estimate and from the lower bound, so a transfer that rests on the point
    estimate alone is visible as such.
    """
    in_band = (v >= thr) if side == "upper" else (v <= thr)
    out = {"n_pos_band": int((in_band & (y == 1)).sum()),
           "n_neg_band": int((in_band & (y == 0)).sum())}
    lr = K.band_lr(y, v, thr, side)
    out["lr"] = lr
    out["lr_lo"] = out["lr_hi"] = np.nan
    if np.isfinite(lr):
        rng = _rng_for(*cell)
        ip, ineg = np.where(y == 1)[0], np.where(y == 0)[0]
        draws = np.empty(N_BOOT)
        for b in range(N_BOOT):
            idx = np.concatenate([rng.choice(ip, ip.size), rng.choice(ineg, ineg.size)])
            draws[b] = K.band_lr(y[idx], v[idx], thr, side)
        draws = draws[np.isfinite(draws)]
        if draws.size >= N_BOOT * 0.5:
            out["lr_lo"], out["lr_hi"] = (float(q) for q in np.percentile(draws, [2.5, 97.5]))
    out["tier"] = K.tier_of(lr, bands, direction)
    # the conservative end of the interval: the lower bound on the pathogenic side,
    # the upper bound on the benign side
    bound = out["lr_lo"] if direction == "pathogenic" else out["lr_hi"]
    out["tier_at_bound"] = K.tier_of(bound, bands, direction)
    return out


def apply_thresholds(gene: str) -> None:
    cfg = yaml.safe_load(CONFIG_PATH.read_text())
    pp3 = cfg["thresholds"]["pp3"]["value"]
    bp4 = cfg["thresholds"]["bp4"]["value"]
    path_bands, ben_bands = K.acmg_bands(cfg)

    df, label_cols = _load_external(gene)
    ev_path = REPORT_DIR / "evidence_thresholds.csv"
    ev = K.read_back(ev_path) if ev_path.exists() else None

    rows = []
    tools = [t for t in K.PANEL + K.OPTIONAL if t in df.columns]
    basis = column_basis_check(df, tools, gene.upper())
    # every panel column is accounted for: tested, or named with the reason it is not
    absent = [t for t in K.PANEL + K.OPTIONAL if t not in df.columns]
    if any(_no_same_variable(t) is None for t in absent):
        raise SystemExit(f"[E7] {gene.upper()}: no column for "
                         f"{[t for t in absent if _no_same_variable(t) is None]} and no recorded reason")
    basis = pd.concat([basis, pd.DataFrame([{"tool": t, "comparable": False,
                                             "reason": _no_same_variable(t)} for t in absent])],
                      ignore_index=True)
    order = {t: i for i, t in enumerate(K.PANEL + K.OPTIONAL)}
    basis = basis.sort_values("tool", key=lambda c: c.map(order), kind="mergesort").reset_index(drop=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    basis.to_csv(REPORT_DIR / f"external_{gene.lower()}_column_basis.csv", index=False)
    comparable = set(basis.loc[basis["comparable"].astype(bool), "tool"])
    if set(tools) - comparable:
        print(f"[E7] {gene.upper()}: E3 thresholds NOT carried onto "
              f"{sorted(set(tools) - comparable)} -- not the same variable as the "
              "column they were fitted on; see the column-basis table")
    for label in label_cols:
        for stratum in K.STRATA:
            sub = (df if stratum == "all_1_50"
                   else df[df["stratum"] != "pm12"] if stratum == "s3_50"
                   else df[df["stratum"] == stratum])
            for tool in tools:
                s = sub.dropna(subset=[label, tool])
                y = s[label].to_numpy(dtype=float)
                v = s[tool].to_numpy(dtype=float)
                base = {"gene": gene.upper(), "label_definition": label,
                        "stratum": stratum, "tool": tool, "n": int(len(s)),
                        "n_pos": int((y == 1).sum()), "n_neg": int((y == 0).sum()),
                        # pm12 and any pool containing it sit outside the ClinGen
                        # SVI recommendation; without this flag a reader takes a
                        # pm12-inclusive row as a threshold result
                        "walker_in_scope": stratum in K.IN_SCOPE_STRATA}
                if base["n_pos"] < K.MIN_POS or base["n_neg"] < K.MIN_NEG:
                    rows.append(base | {"status": "not evaluable",
                                        "reason": "fewer than 10 labelled on one side"})
                    continue
                row = base | {"status": "ok", "reason": ""}
                # --- Walker's fixed cut points, applied unchanged
                if tool in ("spliceai", "spliceai_walker"):
                    row["walker_frac_pp3"] = float((v >= pp3).mean())
                    row["walker_frac_grey"] = float(((v > bp4) & (v < pp3)).mean())
                    row["walker_sens_pp3"] = float((v[y == 1] >= pp3).mean())
                    row["walker_spec_pp3"] = float((v[y == 0] < pp3).mean())
                    row["walker_lr_pp3"] = K.band_lr(y, v, pp3, "upper")
                    row["walker_lr_bp4"] = K.band_lr(y, v, bp4, "lower")
                    row["walker_tier_pp3"] = K.tier_of(row["walker_lr_pp3"],
                                                      path_bands, "pathogenic")
                    row["walker_tier_bp4"] = K.tier_of(
                        K.band_lr_benign_bound(y, v, bp4), ben_bands, "benign")
                    d = _band_detail(y, v, pp3, "upper", path_bands, "pathogenic",
                                     (gene.upper(), label, stratum, tool, "walker_pp3"))
                    row |= {"walker_pp3_n_pos_band": d["n_pos_band"],
                            "walker_pp3_n_neg_band": d["n_neg_band"],
                            "walker_lr_pp3_lo": d["lr_lo"], "walker_lr_pp3_hi": d["lr_hi"],
                            "walker_tier_pp3_at_bound": d["tier_at_bound"]}
                    d = _band_detail(y, v, bp4, "lower", ben_bands, "benign",
                                     (gene.upper(), label, stratum, tool, "walker_bp4"))
                    row |= {"walker_bp4_n_pos_band": d["n_pos_band"],
                            "walker_bp4_n_neg_band": d["n_neg_band"],
                            "walker_lr_bp4_lo": d["lr_lo"], "walker_lr_bp4_hi": d["lr_hi"]}
                row["column_basis_comparable"] = tool in comparable
                # --- E3's thresholds, applied unchanged (E18.5). The main test uses the
                # in-sample thresholds Table 3(b) recommends, read at full precision from
                # the frozen evidence_thresholds.csv: PP3 at Supporting, Moderate and Strong
                # and BP4 at Supporting and Moderate, each only where E3 reached the tier
                # in-sample; a tier not reached carries no threshold. A PP3 tier is read on
                # the upper side by its point estimate, with the tier at the interval's
                # lower bound beside it; a BP4 tier on the lower side by
                # evid_common.band_lr_benign_bound, as the published BP4 cut point is. The
                # medians of the leave-one-gene-out folds' PP3 thresholds, the external
                # test's thresholds before E18, stay as a sensitivity analysis, computed as
                # before (the same cells, the same draws) in columns of their own.
                if ev is not None and tool in comparable:
                    e = ev[(ev.tool == tool) & (ev.stratum == stratum)
                           & (ev.status == "ok")]
                    if not len(e):
                        row["e3_thresholds"] = (
                            f"none: E3 has no rows for tool={tool}, stratum={stratum}")
                    for tier in ["supporting", "moderate", "strong"]:
                        t = e[e.tier == tier]
                        if not len(t):
                            continue
                        reachable = bool(t["pp3_threshold_reachable"].iloc[0])
                        p = f"pp3_insample_{tier}"
                        row[f"{p}_reachable"] = reachable
                        if reachable:
                            tau = float(t["pp3_threshold_insample"].iloc[0])
                            d = _band_detail(y, v, tau, "upper", path_bands, "pathogenic",
                                             (gene.upper(), label, stratum, tool,
                                              f"insample_{tier}"))
                            row |= {f"{p}_threshold": tau, f"{p}_lr_here": d["lr"],
                                    f"{p}_tier_here": d["tier"],
                                    f"{p}_n_pos_band": d["n_pos_band"],
                                    f"{p}_n_neg_band": d["n_neg_band"],
                                    f"{p}_lr_lo": d["lr_lo"], f"{p}_lr_hi": d["lr_hi"],
                                    f"{p}_tier_at_bound": d["tier_at_bound"]}
                        else:
                            row |= {f"{p}_threshold": np.nan, f"{p}_lr_here": np.nan,
                                    f"{p}_tier_here": "not evaluable",
                                    f"{p}_note": "tier not reached in-sample: no threshold"}
                        # the sensitivity analysis: the LOGO median. A tier E3 could not
                        # reach in-sample still has a finite LOGO median, because that
                        # median is taken over only the folds where a threshold was found;
                        # carrying it over would publish a "moderate" cut looser than the
                        # "supporting" cut.
                        q = f"pp3_logo_median_{tier}"
                        folds = float(t["pp3_logo_folds_reaching"].iloc[0] or 0)
                        row[f"{q}_reachable_in_sample"] = reachable
                        row[f"{q}_logo_folds"] = folds
                        if not reachable:
                            row[f"{q}_threshold"] = np.nan
                            row[f"{q}_lr_here"] = np.nan
                            row[f"{q}_tier_here"] = "not evaluable"
                            row[f"{q}_note"] = (
                                "tier unreachable in-sample; LOGO median not carried over")
                            continue
                        tau = float(t["pp3_logo_threshold_median"].iloc[0])
                        row[f"{q}_threshold"] = tau
                        d = _band_detail(y, v, tau, "upper", path_bands, "pathogenic",
                                         (gene.upper(), label, stratum, tool, tier))
                        row[f"{q}_lr_here"] = d["lr"]
                        row[f"{q}_tier_here"] = d["tier"]
                        row[f"{q}_n_pos_band"] = d["n_pos_band"]
                        row[f"{q}_n_neg_band"] = d["n_neg_band"]
                        row[f"{q}_lr_lo"] = d["lr_lo"]
                        row[f"{q}_lr_hi"] = d["lr_hi"]
                        row[f"{q}_tier_at_bound"] = d["tier_at_bound"]
                    for tier in ["supporting", "moderate"]:
                        t = e[e.tier == tier]
                        if not len(t):
                            continue
                        reachable = bool(t["bp4_threshold_reachable"].iloc[0])
                        p = f"bp4_insample_{tier}"
                        row[f"{p}_reachable"] = reachable
                        if not reachable:
                            row |= {f"{p}_threshold": np.nan, f"{p}_lr_here": np.nan,
                                    f"{p}_tier_here": "not evaluable",
                                    f"{p}_note": "tier not reached in-sample: no threshold"}
                            continue
                        tau = float(t["bp4_threshold_insample"].iloc[0])
                        d = _band_detail(y, v, tau, "lower", ben_bands, "benign",
                                         (gene.upper(), label, stratum, tool,
                                          f"insample_bp4_{tier}"))
                        row |= {f"{p}_threshold": tau, f"{p}_lr_here": d["lr"],
                                f"{p}_tier_here": K.tier_of(K.band_lr_benign_bound(y, v, tau),
                                                            ben_bands, "benign"),
                                f"{p}_n_pos_band": d["n_pos_band"],
                                f"{p}_n_neg_band": d["n_neg_band"],
                                f"{p}_lr_lo": d["lr_lo"], f"{p}_lr_hi": d["lr_hi"]}
                rows.append(row)

    out = pd.DataFrame(rows)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    path = REPORT_DIR / f"external_{gene.lower()}.csv"
    out.to_csv(path, index=False)
    print(f"[E7] wrote {path} ({len(out)} rows)")
    ok = out[out.status == "ok"]
    cols = [c for c in ["label_definition", "stratum", "tool", "n", "n_pos", "n_neg",
                        "walker_sens_pp3", "walker_spec_pp3", "walker_lr_pp3",
                        "walker_tier_pp3", "walker_lr_bp4", "walker_tier_bp4"]
            if c in ok.columns]
    sa = ok[ok.tool.isin(["spliceai", "spliceai_walker"])]
    if len(sa):
        print("\n--- Walker cut points on this gene ---")
        print(sa[cols].round(4).to_string(index=False))
    ne = out[out.status != "ok"]
    print(f"\n[E7] not evaluable: {len(ne)} cells")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--prepare", choices=["ddx3x"])
    ap.add_argument("--merge-scores", choices=["ddx3x"])
    ap.add_argument("--apply")
    args = ap.parse_args()
    if args.prepare == "ddx3x":
        prepare_ddx3x()
    elif args.merge_scores == "ddx3x":
        merge_ddx3x_scores()
    elif args.apply:
        apply_thresholds(args.apply)
    else:
        ap.error("pass --prepare ddx3x, --merge-scores ddx3x or --apply <gene>")


if __name__ == "__main__":
    main()
