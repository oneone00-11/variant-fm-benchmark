"""Supplementary tables -- the print versions of the evidence-strength outputs.

What this stage answers: what a reader of the supplement sees, cell for cell, and
which file each cell came from. It writes one CSV per supplementary table into
reports/evidence/supplement/, and a manifest giving each table's title, its source
files with their sha256, its row count and its own sha256.

It fits nothing. Every value is read from a file an earlier stage wrote or from a
tracked input, and the only arithmetic done here is arithmetic a reader could redo
from the same files: totals of counts, and a rate or share from its numerator and
denominator. The choices that change what is printed, and why:

  * Rounding happens once, here, from the unrounded stored value: likelihood
    ratios to three significant figures, proportions and AUROC to three decimals,
    thresholds to four significant figures. It is half-up on the value's shortest
    decimal representation, which is the digit string the source CSV itself holds,
    so a reader rounding the CSV by hand gets the printed digit. Python's format()
    and round() work on the binary value instead, which for a stored 2.675 lies
    just under the tie, and print 2.67 where the reader writes 2.68.
  * Where a source stores a quantity already rounded, it is recomputed from its
    parts rather than rounded again, because rounding a rounded number can move the
    last printed digit. The in-frame shares (stored at four decimals) are rebuilt
    from their counts, and the per-fold held-out ratios of S6 from the per-fold
    file rather than parsed out of the four-figure string in tier_logo_table.csv.
  * Where two upstream files carry the same quantity, both are read and the stage
    stops if they disagree: set_counts.csv against a recount of the analysis set,
    tier_logo_table.csv against evidence_thresholds.csv and the per-fold file, and
    the external-gene tables against the figure table built from them. A supplement
    that silently disagrees with a figure is the failure this is here to prevent.
  * S10 carries none of the walker_* columns of territory_metrics_arms_noBRCA1.csv.
    Those apply SpliceAI's 0.2 cut point to every column's raw scale so that the
    table stays rectangular, and for AlphaGenome, CADD or phyloP the result has no
    meaning. The published-cut-point ratio in S10 comes from
    arms_at_tool_threshold.csv, which computes it for the two SpliceAI columns only.
  * The constructed labels of the external genes are described from the
    implementation rather than restated. TP53's cut values are read out of
    evid_external's source and re-applied to the data, and the stage stops if they
    do not reproduce the labels that module produces; DDX3X's come from the
    manifest its preparation step wrote, with the FDR cut checked against the
    module constant.
  * Inputs whose manifests record a hash are checked against it before use.
  * The output is deterministic and machine-independent: no timestamp, no absolute
    path (atlas files are named relative to the atlas repository), every cell
    written as text so that pandas' float printing cannot vary it, and rows in a
    fixed order.

S8 and S9 are written by their own stages; the numbering here is fixed and leaves
them their places. The manifest lists the tables this stage writes.

Run (PYTHONPATH=phase1):  python -m src.evid_supp_tables
"""
from __future__ import annotations

import argparse
import contextlib
import hashlib
import importlib.util
import inspect
import json
import os
import re
import sys
from decimal import ROUND_HALF_UP, Decimal
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from . import config as C
from . import evid_common as K
from . import evid_external as EXT
from . import evid_inframe as INF
from . import evid_walker_thresholds as WT
from .evid_training_provenance import FROM_ATLAS
from .phase1_build_frozen_matrix_v2 import canonical_sha256

PHASE1 = Path(__file__).resolve().parents[1]
REPO = PHASE1.parent
from .evid_common import ATLAS_REPO as _atlas_repo_default  # noqa: E402
ATLAS_REPO = _atlas_repo_default  # EVID_ATLAS_REPO overrides; see evid_common
REPORTS = PHASE1 / "reports/evidence"
OUT_DIR = REPORTS / "supplement"
MANIFEST = "supplement_tables_manifest.csv"

SET_PATH = PHASE1 / "data/evidence/analysis_set_v1.parquet"
SET_MANIFEST = PHASE1 / "data/evidence/analysis_set_v1.manifest.json"
DDX3X_PARQUET = PHASE1 / "data/evidence/external/ddx3x_splice.parquet"
DDX3X_MANIFEST = PHASE1 / "data/evidence/external/ddx3x_splice.manifest.json"
TP53_PARQUET = PHASE1 / "data/external/tp53_splice_scored_v2.parquet"
TP53_MANIFEST = PHASE1 / "data/external/tp53_splice_scored_v2.manifest.json"
TP53_DEPOSIT = REPO / "data/external/tp53_mavedb_scores.tsv"
WALKER_CFG = PHASE1 / "config/walker2023.yaml"
WALKER_PROV = PHASE1 / "data/evidence/spliceai_walker.parquet.provenance.json"
AVI_PROV = PHASE1 / "data/evidence/avi.parquet.provenance.json"
LICENSE_DATA = REPO / "LICENSE-DATA"
COLUMN_PROV = REPO / "docs/column-provenance.md"
ATLAS_COMPOSITION = ATLAS_REPO / "results/table1_atlas_composition.tsv"
ATLAS_RESOURCES = ATLAS_REPO / "src/atlas/predictor_resources.py"

SEVEN = ["BAP1", "BARD1", "BRCA1", "BRCA2", "PALB2", "RAD51C", "VHL"]
EXTERNAL = ["DDX3X", "TP53"]
# Short names of the assay publications, as the manuscript cites them.
PUBLICATION = {
    "BRCA1": "Findlay 2018", "BRCA2": "Huang 2025", "BARD1": "Woo 2025",
    "PALB2": "Boonen 2026", "RAD51C": "Olvera-Leon 2024", "VHL": "Buckley 2024",
    "BAP1": "Waters 2024", "DDX3X": "Radford 2023", "TP53": "Funk 2025",
}

TOOLS = K.PANEL + K.OPTIONAL
IN_SCOPE = ["s3_10", "s11_50", "s3_50"]
TIERS = ["supporting", "moderate", "strong"]
ARMS = ["classified", "recorded_unclassified", "unrecorded"]
STRICT_ARMS = ["classified", "recorded_unclassified", "recorded_no_assertion",
               "unrecorded"]
WALKER_TOOLS = ("spliceai", "spliceai_walker")

STRATUM_LABEL = {
    "pm12": "±1,2 (out of scope)", "s3_10": "3–10", "s11_50": "11–50",
    "s3_50": "3–50", "all_1_50": "1–50 incl. ±1,2, out of scope",
}
ARM_LABEL = {
    "all": "All variants", "classified": "Classified (P/LP or B/LB)",
    "recorded_unclassified": "Recorded, not classified",
    "unrecorded": "Not in ClinVar",
}
STRICT_ARM_LABEL = {
    "classified": "Classified (P/LP or B/LB)",
    "recorded_unclassified": "Recorded: VUS, conflicting or other",
    "recorded_no_assertion": "Recorded: no clinical assertion",
    "unrecorded": "Not in ClinVar",
}
GENE_SET_LABEL = {"all_genes": "All genes", "no_BRCA1": "BRCA1 excluded"}
LABEL_NAME = {"y_control_anchored": "Control-anchored", "y_fdr": "FDR",
              "y_median_split": "Median split",
              "y_mid_band_excluded": "Mid-band excluded"}

# What each panel column is, in words. Numbers belong to the configuration column,
# which reads them from the provenance files.
WHAT = {
    "spliceai": "Largest of SpliceAI's four delta scores (acceptor and donor, gain "
                "and loss) at the companion atlas's scoring window; the panel's "
                "ranking column",
    "spliceai_walker": "The same SpliceAI statistic re-scored at the window the "
                       "published cut points were calibrated on; the column those "
                       "cut points are applied to",
    "pangolin": "Pangolin's largest predicted change in splice-site usage, gain or "
                "loss",
    "alphagenome": "AlphaGenome splice score, merged-quantile definition, computed "
                   "by the model at query time",
    "cadd": "CADD deleteriousness score, PHRED-scaled",
    "phylop": "Cross-species conservation at the variant base (phyloP)",
    "phastcons": "Conserved-element posterior probability at the variant base "
                 "(phastCons)",
    "gpn_msa": "GPN-MSA alignment-conditioned language-model score, sign-flipped so "
               "that larger is more damaging",
    "nt": "Nucleotide Transformer masked-token log-likelihood ratio",
}


# ---------------------------------------------------------------------------
# formatting: every printed number passes through one of these
# ---------------------------------------------------------------------------
def _finite(x) -> bool:
    try:
        return x is not None and not pd.isna(x) and np.isfinite(float(x))
    except (TypeError, ValueError):
        return False


def _plain(q: Decimal) -> str:
    s = format(q, "f")
    return s[1:] if s.startswith("-") and Decimal(s) == 0 else s


def fmt_sig(x, sig: int) -> str:
    """Half-up to `sig` significant figures, on the shortest decimal repr."""
    if not _finite(x):
        return ""
    x = float(x)
    if x == 0:
        return "0"
    d = Decimal(repr(x))
    exp = d.adjusted()
    q = d.quantize(Decimal(1).scaleb(exp - sig + 1), rounding=ROUND_HALF_UP)
    if q.adjusted() > exp:                  # 9.996 -> 10.00: one digit too many
        q = q.quantize(Decimal(1).scaleb(exp - sig + 2), rounding=ROUND_HALF_UP)
    return _plain(q)


def fmt_dp(x, dp: int = 3) -> str:
    """Half-up to `dp` decimals, on the shortest decimal repr."""
    if not _finite(x):
        return ""
    q = Decimal(repr(float(x))).quantize(Decimal(1).scaleb(-dp),
                                        rounding=ROUND_HALF_UP)
    return _plain(q)


def fmt_lr(x) -> str:
    return fmt_sig(x, 3)


_SCORES: dict[str, np.ndarray] = {}


def _scores(tool: str) -> np.ndarray:
    if tool not in _SCORES:
        df = pd.read_parquet(SET_PATH)
        _SCORES[tool] = (np.sort(df[tool].dropna().to_numpy(dtype=float))
                         if tool in df.columns else np.array([]))
    return _SCORES[tool]


def fmt_thr(x, tool: str | None = None, side: str = "upper") -> str:
    """Four significant figures, or more where four would move a variant across
    the threshold. AlphaGenome's thresholds sit just under its ceiling of 2.2, and
    at four figures a printed 2.200 selects no variant at all; the printed value
    must select exactly the band the fitted value selects."""
    if not _finite(x):
        return ""
    s = _scores(tool) if tool else np.array([])
    for sig in range(4, 12):
        txt = fmt_sig(x, sig)
        if not len(s):
            return txt
        v = float(txt)
        if side == "upper":
            same = int((s >= v).sum()) == int((s >= float(x)).sum())
        else:
            same = int((s <= v).sum()) == int((s <= float(x)).sum())
        if same:
            return txt
    return fmt_sig(x, 12)


def fmt_int(x) -> str:
    return str(int(round(float(x)))) if _finite(x) else ""


def fmt_ci(v, lo, hi, fmt=fmt_lr) -> str:
    if not _finite(v):
        return ""
    if _finite(lo) and _finite(hi):
        return f"{fmt(v)} ({fmt(lo)}–{fmt(hi)})"
    return f"{fmt(v)} (no interval)"


def fmt_const(x) -> str:
    """A configuration constant, printed as the configuration writes it."""
    return repr(float(x)).rstrip("0").rstrip(".") if float(x) != int(x) else str(int(x))


def yes_no(v) -> str:
    if v is None or (not isinstance(v, str) and pd.isna(v)):
        return ""
    if isinstance(v, str):
        v = v.strip().lower() in ("true", "1", "yes")
    return "Yes" if bool(v) else "No"


def tier_text(t) -> str:
    if not isinstance(t, str) or not t:
        return ""
    return {"below supporting": "Below Supporting"}.get(t, t)


def _text(v) -> str:
    return "" if v is None or (not isinstance(v, str) and pd.isna(v)) else str(v)


def _n_variants(n) -> str:
    n = int(round(float(n)))
    return f"{n} variant" if n == 1 else f"{n} variants"


def _plain_reason(v) -> str:
    """Upstream reasons name the calibration stage by its internal label."""
    return re.sub(r"\bE3\b", "the interval calibration", _text(v))


# ---------------------------------------------------------------------------
# inputs
# ---------------------------------------------------------------------------
@contextlib.contextmanager
def _in_phase1():
    """The stage modules resolve their paths relative to phase1/."""
    here = Path.cwd()
    os.chdir(PHASE1)
    try:
        yield
    finally:
        os.chdir(here)


def _sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _display(p: Path) -> str:
    """A path a reader on another machine can resolve: repo-relative, or relative
    to the atlas repository."""
    p = Path(p).resolve()
    for root, prefix in ((REPO.resolve(), ""), (ATLAS_REPO.resolve(), "atlas:")):
        try:
            return prefix + p.relative_to(root).as_posix()
        except ValueError:
            continue
    raise ValueError(f"source outside both repositories: {p}")


def _need(p: Path) -> Path:
    if not p.exists():
        hint = (" (set EVID_ATLAS_REPO to the companion atlas checkout)"
                if str(p).startswith(str(ATLAS_REPO)) else "")
        raise SystemExit(f"[supp] missing input {p}{hint}")
    return p


def _csv(p: Path, **kw) -> pd.DataFrame:
    return pd.read_csv(_need(p), **kw)


@lru_cache(maxsize=1)
def _walker_cfg() -> dict:
    return yaml.safe_load(_need(WALKER_CFG).read_text())


@lru_cache(maxsize=1)
def _analysis_set() -> pd.DataFrame:
    man = json.loads(_need(SET_MANIFEST).read_text())
    if _sha256(_need(SET_PATH)) != man["sha256"]:
        raise SystemExit(f"[supp] {SET_PATH.name} does not match the sha256 in "
                         f"{SET_MANIFEST.name}; rebuild the set or its manifest")
    return K.load_set(SET_PATH)


@lru_cache(maxsize=1)
def _atlas_resources():
    name = "atlas_predictor_resources"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _need(ATLAS_RESOURCES))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def _resource_rows() -> dict:
    return {r[0]: dict(zip(["predictor", "version", "source", "accessed", "licence",
                            "licence_source"], r))
            for r in _atlas_resources().ROWS}


@lru_cache(maxsize=None)
def _json(p: Path) -> dict:
    return json.loads(_need(p).read_text())


def _distance_label(d: int) -> str:
    return f"±{d:,} nt"


@lru_cache(maxsize=1)
def tool_labels() -> dict:
    """Printed name of every score column, built from the provenance sources."""
    rows = _resource_rows()
    cfg = _walker_cfg()
    walker = _json(WALKER_PROV)
    avi = _json(AVI_PROV)
    lab = {}
    for col, name in FROM_ATLAS.items():
        lab[col] = rows[name]["predictor"]
    d50 = int(cfg["column_check"]["atlas_definition"]["distance"])
    lab["spliceai"] = f"SpliceAI ({_distance_label(d50)})"
    lab["spliceai_walker"] = (f"SpliceAI ({_distance_label(int(walker['distance']))}, "
                              "Walker basis)")
    if len(avi["columns"]) != len(avi["scorers_requested"]):
        raise SystemExit("[supp] AVI provenance: columns and scorers do not pair")
    for col, scorer in zip(avi["columns"], avi["scorers_requested"]):
        name = "AVI" if scorer == "AVI_SCORE" else scorer.replace("_", " ").lower()
        lab[col] = f"AlphaGenome Atlas {name}"
    v = re.fullmatch(r"alphagenome_v(\d)(\d)(\d)", "alphagenome_v061")
    lab["alphagenome_v061"] = f"AlphaGenome (client {'.'.join(v.groups())} definition)"
    missing = [t for t in TOOLS if t not in lab]
    if missing:
        raise SystemExit(f"[supp] no printed name for {missing}")
    if len(set(lab.values())) != len(lab):
        raise SystemExit("[supp] two score columns share a printed name")
    return lab


def _tool(t: str) -> str:
    return tool_labels()[t]


# ---------------------------------------------------------------------------
# S1 -- functional standards
# ---------------------------------------------------------------------------
def _tp53_rules() -> dict:
    """TP53's constructed labels, read from the code that constructs them."""
    src = inspect.getsource(EXT._load_external)
    ca = re.search(r'df\["y_control_anchored"\] = \(rfs > ([-0-9.]+)\)', src)
    mb = re.search(r'mid = np\.abs\(rfs\) >= ([0-9.]+)', src)
    ms = re.search(r'df\["y_median_split"\] = \(rfs > np\.median\(rfs\)\)', src)
    if not (ca and mb and ms):
        raise SystemExit("[supp] TP53 label construction in evid_external changed; "
                         "the S1 wording must be revisited")
    return {"control": ca.group(1), "midband": mb.group(1)}


def _external_labels(gene: str) -> tuple[pd.DataFrame, list[str]]:
    with _in_phase1():
        return EXT._load_external(gene)


def _check_tp53_rules(df: pd.DataFrame, rules: dict) -> float:
    rfs = df["func_pathogenicity"].to_numpy(dtype=float)
    c0, c1 = float(rules["control"]), float(rules["midband"])
    med = float(np.median(rfs))
    want = {
        "y_control_anchored": (rfs > c0).astype(float),
        "y_median_split": (rfs > med).astype(float),
        "y_mid_band_excluded": np.where(np.abs(rfs) >= c1, (rfs > c0).astype(float),
                                        np.nan),
    }
    for col, w in want.items():
        got = df[col].to_numpy(dtype=float)
        if not np.array_equal(got, w, equal_nan=True):
            raise SystemExit(f"[supp] TP53 {col}: the rule read from evid_external "
                             "does not reproduce its labels")
    return med


def _label_counts(df: pd.DataFrame, col: str) -> dict:
    y = df[col]
    return {"labelled": int(y.notna().sum()), "damaging": int((y == 1).sum()),
            "normal": int((y == 0).sum()), "missing": int(y.isna().sum())}


def _offsets(s: pd.Series) -> str:
    return f"{fmt_int(s.min())}–{fmt_int(s.max())}"


def s1() -> tuple[pd.DataFrame, list[Path]]:
    df = _analysis_set()
    comp = (pd.read_csv(ATLAS_COMPOSITION, sep="\t").set_index("gene")
            if ATLAS_COMPOSITION.exists() else None)
    sc = _csv(REPORTS / "set_counts.csv")
    rows = []
    for g in SEVEN:
        spec = C.ASSAY_LABELS[g]
        sub = df[df["gene"] == g]
        # The analysis set carries one accession and one transcript per gene; it
        # is the tracked source. The atlas composition table is gitignored in the
        # atlas repository, so it is used only to cross-check when present.
        if sub["urn"].nunique() != 1 or sub["transcript"].nunique() != 1:
            raise SystemExit(f"[supp] {g}: more than one accession or transcript")
        urn, tx = sub["urn"].iloc[0], sub["transcript"].iloc[0]
        if comp is not None and (comp.loc[g, "mavedb_urn"] != urn
                                 or comp.loc[g, "transcript"] != tx):
            raise SystemExit(f"[supp] {g}: analysis set and atlas composition table "
                             "disagree on accession or transcript")
        y = sub["y_assay"]
        n, lab = len(sub), int(y.notna().sum())
        dam, nor = int((y == 1).sum()), int((y == 0).sum())
        c = sc[sc["gene"] == g][["n", "n_labelled", "n_damaging", "n_normal"]].sum()
        if (int(c.n), int(c.n_labelled), int(c.n_damaging), int(c.n_normal)) != \
                (n, lab, dam, nor):
            raise SystemExit(f"[supp] {g}: set_counts.csv disagrees with the set")
        src = sub["y_assay_source"].value_counts()
        cls = spec["map"]
        rows.append({
            "Gene": g, "Role": "Analysis set",
            "Assay publication": PUBLICATION[g],
            "MaveDB accession": urn, "Transcript": tx,
            "Published class field": spec["class_col"],
            "Classes mapped to damaging": "; ".join(k for k, v in cls.items() if v == 1.0),
            "Classes mapped to normal": "; ".join(k for k, v in cls.items() if v == 0.0),
            "Classes mapped to missing": "; ".join(k for k, v in cls.items() if v is None),
            "Intronic offset (bp)": _offsets(sub["intron_offset_abs"]),
            "Intron-side SNVs": str(n), "Labelled": str(lab),
            "Damaging": str(dam), "Normal": str(nor), "Not labelled": str(n - lab),
            "Note": (f"Not labelled: {int(src.get('indeterminate', 0))} in a class "
                     f"mapped to missing, {int(src.get('unlabelled', 0))} absent from "
                     "the label file"),
        })

    # DDX3X: the deposit publishes scores only; both label definitions are built
    # by evid_external's preparation step and recorded in its manifest
    man = _json(DDX3X_MANIFEST)
    if _sha256(_need(DDX3X_PARQUET)) != man["sha256"]:
        raise SystemExit("[supp] ddx3x_splice.parquet does not match its manifest")
    fdr_text = man["labels"]["fdr"]
    if f"< {EXT.FDR_CUT}" not in fdr_text:
        raise SystemExit("[supp] DDX3X manifest FDR cut differs from evid_external")
    d, dcols = _external_labels("DDX3X")
    ctrl = man["controls"]
    n_sets = len(man["score_sets"])
    if ctrl.get("n_scoresets_without_controls", 0):
        raise SystemExit("[supp] a DDX3X score set lacks controls; S1 wording assumes "
                         "every set is anchored on its own controls")
    cnt = {c: _label_counts(d, c) for c in dcols}
    miss = sum(v["missing"] for v in cnt.values())
    rows.append({
        "Gene": "DDX3X", "Role": "External gene",
        "Assay publication": PUBLICATION["DDX3X"],
        "MaveDB accession": (f"{man['score_sets'][0]} to {man['score_sets'][-1]} "
                             f"({n_sets} score sets)"),
        "Transcript": man["transcript"],
        "Published class field": "None published; labels constructed from the assay score",
        "Classes mapped to damaging": (
            f"{LABEL_NAME['y_control_anchored']} (primary): score below the midpoint "
            f"of the synonymous and nonsense control medians of its own score set, "
            f"each of the {n_sets} sets anchored separately. "
            f"{LABEL_NAME['y_fdr']}: trend BH-FDR < {EXT.FDR_CUT} and score below "
            "the synonymous median"),
        "Classes mapped to normal": "Every other scored variant, under each definition",
        "Classes mapped to missing": (
            "No variant: every scored variant is labelled under both definitions"
            if miss == 0 else f"{miss} variants without a label"),
        "Intronic offset (bp)": _offsets(d["intron_offset_abs"]),
        "Intron-side SNVs": str(len(d)),
        "Labelled": str(cnt["y_control_anchored"]["labelled"]),
        "Damaging": str(cnt["y_control_anchored"]["damaging"]),
        "Normal": str(cnt["y_control_anchored"]["normal"]),
        "Not labelled": str(cnt["y_control_anchored"]["missing"]),
        "Note": "Counts are for the primary definition. " + "; ".join(
            f"{LABEL_NAME[c]}: {v['damaging']} damaging, {v['normal']} normal"
            for c, v in cnt.items()) + f". {man['labels']['note'].capitalize()}.",
    })

    # TP53: the published study's constructions, carried unchanged by evid_external
    tman = _json(TP53_MANIFEST)
    t, tcols = _external_labels("TP53")
    if canonical_sha256(pd.read_parquet(TP53_PARQUET)) != tman["sha256"]:
        raise SystemExit("[supp] tp53_splice_scored_v2.parquet does not match its "
                         "manifest's content hash")
    rules = _tp53_rules()
    med = _check_tp53_rules(t, rules)
    dep = _csv(_need(TP53_DEPOSIT), sep="\t", usecols=["accession", "HGVS(cDNA)"])
    acc = dep["accession"].str.split("#").str[0].unique()
    tx = t["hgvs_nt"].astype(str).str.split(":").str[0].unique()
    dep_tx = dep["HGVS(cDNA)"].astype(str).str.split(":").str[0].unique()
    if len(acc) != 1 or len(tx) != 1 or set(dep_tx) != set(tx):
        raise SystemExit("[supp] TP53 accession or transcript is not unique")
    cnt = {c: _label_counts(t, c) for c in tcols}
    c0, c1 = rules["control"], rules["midband"]
    rows.append({
        "Gene": "TP53", "Role": "External gene",
        "Assay publication": PUBLICATION["TP53"],
        "MaveDB accession": acc[0], "Transcript": tx[0],
        "Published class field": ("None published; labels constructed from the "
                                  "relative fitness score (RFS)"),
        "Classes mapped to damaging": (
            f"{LABEL_NAME['y_control_anchored']} (primary): RFS > {c0}. "
            f"{LABEL_NAME['y_median_split']}: RFS above the median of the analysed "
            f"variants ({fmt_sig(med, 3)}). {LABEL_NAME['y_mid_band_excluded']}: "
            f"RFS > {c0} with |RFS| ≥ {c1}"),
        "Classes mapped to normal": (
            f"{LABEL_NAME['y_control_anchored']}: RFS ≤ {c0}. "
            f"{LABEL_NAME['y_median_split']}: RFS at or below the median. "
            f"{LABEL_NAME['y_mid_band_excluded']}: RFS ≤ {c0} with |RFS| ≥ {c1}"),
        "Classes mapped to missing": (
            f"{LABEL_NAME['y_mid_band_excluded']} only: |RFS| < {c1} "
            f"({cnt['y_mid_band_excluded']['missing']} variants)"),
        "Intronic offset (bp)": _offsets(t["intron_offset_abs"]),
        "Intron-side SNVs": str(len(t)),
        "Labelled": str(cnt["y_control_anchored"]["labelled"]),
        "Damaging": str(cnt["y_control_anchored"]["damaging"]),
        "Normal": str(cnt["y_control_anchored"]["normal"]),
        "Not labelled": str(cnt["y_control_anchored"]["missing"]),
        "Note": "Counts are for the primary definition. " + "; ".join(
            f"{LABEL_NAME[c]}: {v['damaging']} damaging, {v['normal']} normal"
            for c, v in cnt.items()) + ".",
    })
    srcs = [SET_PATH, SET_MANIFEST, REPORTS / "set_counts.csv",
            PHASE1 / "src/config.py", PHASE1 / "src/evid_external.py",
            DDX3X_PARQUET, DDX3X_MANIFEST, TP53_PARQUET, TP53_MANIFEST, TP53_DEPOSIT]
    return pd.DataFrame(rows), srcs


# ---------------------------------------------------------------------------
# S2 -- analysis-set composition
# ---------------------------------------------------------------------------
def s2() -> tuple[pd.DataFrame, list[Path]]:
    sc = _csv(REPORTS / "set_counts.csv")
    df = _analysis_set()
    # set_counts.csv must be the current set's counts, not an earlier build's
    recount = (df.assign(lab=df["y_assay"].notna(), dam=(df["y_assay"] == 1))
             .groupby(["gene", "stratum", "clinvar_arm_strict"])
             .agg(n=("variant_id", "size"), n_labelled=("lab", "sum"),
                  n_damaging=("dam", "sum")).reset_index())
    mine = sc.groupby(["gene", "stratum", "clinvar_arm_strict"])[
        ["n", "n_labelled", "n_damaging"]].sum().reset_index()
    chk = mine.merge(recount, on=["gene", "stratum", "clinvar_arm_strict"],
                     how="outer", suffixes=("", "_set"))
    for c in ("n", "n_labelled", "n_damaging"):
        if not (chk[c].fillna(-1).astype(int) == chk[f"{c}_set"].fillna(-1).astype(int)).all():
            raise SystemExit(f"[supp] set_counts.csv column {c} is stale against the set")

    cols = ["n", "n_labelled", "n_damaging", "n_normal"]

    def row(stratum_label, gene, arm, v):
        lab = int(v["n_labelled"])
        return {"Stratum (intronic offset, bp)": stratum_label, "Gene": gene,
                "ClinVar record status": arm,
                "Variants": str(int(v["n"])), "Labelled": str(lab),
                "Damaging": str(int(v["n_damaging"])),
                "Normal": str(int(v["n_normal"])),
                "Positive rate (damaging / labelled)":
                    fmt_dp(v["n_damaging"] / lab) if lab else ""}

    rows = []
    strata = [s for s in K.STRATA if s in set(sc["stratum"])]
    for st in strata:
        s = sc[sc["stratum"] == st]
        g = s.groupby(["gene", "clinvar_arm_strict"])[cols].sum()
        for gene in SEVEN:
            for arm in STRICT_ARMS:
                if (gene, arm) in g.index:
                    rows.append(row(STRATUM_LABEL[st], gene, STRICT_ARM_LABEL[arm],
                                    g.loc[(gene, arm)]))
        by_arm = s.groupby("clinvar_arm_strict")[cols].sum()
        for arm in STRICT_ARMS:
            if arm in by_arm.index:
                rows.append(row(STRATUM_LABEL[st], "All genes", STRICT_ARM_LABEL[arm],
                                by_arm.loc[arm]))
        rows.append(row(STRATUM_LABEL[st], "All genes", ARM_LABEL["all"], s[cols].sum()))
    for pool, keep, label in (("s3_50", lambda x: x != "pm12", STRATUM_LABEL["s3_50"]),
                              ("all_1_50", lambda x: True, "1–50 incl. ±1,2 (out of scope)")):
        s = sc[sc["stratum"].map(keep)]
        by_arm = s.groupby("clinvar_arm_strict")[cols].sum()
        for arm in STRICT_ARMS:
            if arm in by_arm.index:
                rows.append(row(label, "All genes", STRICT_ARM_LABEL[arm], by_arm.loc[arm]))
        rows.append(row(label, "All genes", ARM_LABEL["all"], s[cols].sum()))
    out = pd.DataFrame(rows)
    total = int(out.iloc[-1]["Variants"])
    if total != len(df):
        raise SystemExit(f"[supp] S2 total {total} != analysis set {len(df)}")
    return out, [REPORTS / "set_counts.csv", SET_PATH]


# ---------------------------------------------------------------------------
# S3 -- predictor columns
# ---------------------------------------------------------------------------
def _licence_sections() -> dict:
    text = _need(LICENSE_DATA).read_text()
    heads = [(m.group(1), m.group(2), m.start())
             for m in re.finditer(r"^(\d+)\. (.+)$", text, flags=re.M)]
    out = {}
    for i, (num, title, start) in enumerate(heads):
        end = heads[i + 1][2] if i + 1 < len(heads) else len(text)
        out[num] = (title, text[start:end])
    return out


def _avi_definitions() -> dict:
    text = _need(COLUMN_PROV).read_text()
    lines = text.splitlines()
    head = "| column | source | definition | terms |"
    i = lines.index(head)
    out = {}
    for ln in lines[i + 2:]:
        if not ln.startswith("|"):
            break
        cells = [c.strip() for c in ln.strip().strip("|").split("|")]
        if len(cells) != 4:
            raise SystemExit(f"[supp] unexpected AVI row in {COLUMN_PROV.name}: {ln}")
        out[cells[0].strip("`")] = cells[2]
    return out


def _pangolin_distance() -> int:
    for ln in _need(COLUMN_PROV).read_text().splitlines():
        if ln.startswith("| pangolin |"):
            found = set(re.findall(r"-d (\d+)", ln))
            if len(found) == 1:
                return int(found.pop())
    raise SystemExit(f"[supp] no single Pangolin distance in {COLUMN_PROV.name}")


def s3() -> tuple[pd.DataFrame, list[Path]]:
    train = _csv(REPORTS / "predictor_training_provenance.csv").set_index("score_column")
    df = _analysis_set()
    rows_by = _resource_rows()
    cfg = _walker_cfg()
    walker = _json(WALKER_PROV)
    avi = _json(AVI_PROV)
    orient = _csv(REPORTS / "feature_orientation.csv").set_index("feature")
    lic = _licence_sections()
    avi_sec = [n for n, (t, _) in lic.items() if "AVI" in t]
    set_sec = [n for n, (t, _) in lic.items() if "frozen analysis set" in t]
    if len(avi_sec) != 1 or len(set_sec) != 1:
        raise SystemExit("[supp] LICENSE-DATA sections not found")
    avi_body = lic[avi_sec[0]][1]
    for col in avi["columns"]:
        if f"`{col}`" not in avi_body or "NON-COMMERCIAL" not in avi_body:
            raise SystemExit(f"[supp] LICENSE-DATA section {avi_sec[0]} does not "
                             f"declare {col}")
    avi_def = _avi_definitions()
    scorer = dict(zip(avi["columns"], avi["scorers_requested"]))
    ad = cfg["column_check"]["atlas_definition"]
    pang_d = _pangolin_distance()

    rows = []
    for col in TOOLS:
        if col in avi["columns"]:
            what = avi_def[col].replace(" AND ", " and ")
            what = (what[0].upper() + what[1:] if scorer[col] == "AVI_SCORE"
                    else f"Atlas {scorer[col]} scorer; {what}")
            version = (f"AlphaGenome API client {avi['client_version']}; Atlas "
                       f"released {avi['atlas_released']}; scorer {scorer[col]}; "
                       "precomputed values retrieved by lookup")
            source = avi["product"]
            accessed = avi["accessed_utc"][:10]
            licence, lic_src = avi["licence"], avi["licence_source"]
            audit = f"LICENSE-DATA section {avi_sec[0]} (this repository)"
        else:
            r = rows_by[FROM_ATLAS[col]]
            what = WHAT[col]
            source, licence, lic_src = r["source"], r["licence"], r["licence_source"]
            audit = (f"Companion atlas audit (predictor_resources.ROWS), to which "
                     f"LICENSE-DATA section {set_sec[0]} defers")
            if col == "spliceai_walker":
                version = (f"{walker['package_version']}; {walker['statistic']}; "
                           f"distance {walker['distance']}; "
                           f"{'masked' if walker['mask'] else 'unmasked'}; "
                           f"annotation {walker['annotation']}; basis "
                           f"{walker['basis'].split(';')[0]}")
                accessed = walker["scored_utc"][:10]
            else:
                version = r["version"]
                accessed = r["accessed"]
                if col == "spliceai":
                    version += (f"; distance {ad['distance']}; "
                                f"{'masked' if ad['masked'] else 'unmasked'}")
                elif col == "pangolin":
                    version += f"; distance {pang_d}"
        rho = orient.loc[col, "spearman_vs_pathogenicity"] if col in orient.index else np.nan
        tr = train.loc[col] if col in train.index else None
        rows.append({
            "Tool": _tool(col), "Score column": col,
            "Role": "Fusion panel" if col in K.PANEL else "Evaluated beside the panel",
            "What the column is": what,
            "Version, resource and scoring configuration": version,
            "Source": source, "Accessed": accessed,
            "Licence": licence, "Licence source": lic_src, "Licence audit": audit,
            "Spearman rho with functional pathogenicity": fmt_dp(rho),
            "Variants scored": str(int(df[col].notna().sum())),
            "Analysis-set variants": str(len(df)),
            # the documentation behind Table 2, one row per column
            "Training signal": tr["training_signal_class"] if tr is not None else "",
            "Training data": tr["training_data"] if tr is not None else "",
            "Clinical classifications in training":
                tr["contains_clinical_labels"] if tr is not None else "",
            "Multiplexed-assay measurements in training":
                tr["documented_mave_or_sge_in_training"] if tr is not None else "",
            "Training documentation": tr["training_source"] if tr is not None else "",
        })
    srcs = [ATLAS_RESOURCES, WALKER_PROV, AVI_PROV, WALKER_CFG, LICENSE_DATA,
            COLUMN_PROV, REPORTS / "feature_orientation.csv", SET_PATH,
            REPORTS / "predictor_training_provenance.csv",
            PHASE1 / "src/evid_common.py"]
    return pd.DataFrame(rows), srcs


# ---------------------------------------------------------------------------
# S4 -- agreement between related columns
# ---------------------------------------------------------------------------
def s4() -> tuple[pd.DataFrame, list[Path]]:
    cc = _csv(REPORTS / "column_concordance.csv")
    rows = []
    for r in cc.itertuples():
        has_cut = _finite(r.cut_point)
        rows.append({
            "Column A": _tool(r.column_a), "Column B": _tool(r.column_b),
            "Variants compared": fmt_int(r.n),
            "Spearman rho": fmt_dp(r.spearman_rho),
            "Cut point": fmt_const(r.cut_point) if has_cut else "",
            "Below the cut on A, at or above it on B":
                fmt_int(r.n_a_below_b_at_or_above) if has_cut else "",
            "At or above the cut on A, below it on B":
                fmt_int(r.n_a_at_or_above_b_below) if has_cut else "",
        })
    return pd.DataFrame(rows), [REPORTS / "column_concordance.csv"]


# ---------------------------------------------------------------------------
# S5 -- the published cut points
# ---------------------------------------------------------------------------
def s5() -> tuple[pd.DataFrame, list[Path]]:
    cfg = _walker_cfg()
    pp3 = fmt_const(cfg["thresholds"]["pp3"]["value"])
    bp4 = fmt_const(cfg["thresholds"]["bp4"]["value"])
    route = cfg["canonical_pm12"]["route"]
    w = _csv(REPORTS / "walker_thresholds.csv")
    w = w[w["score_column"] == "spliceai_walker"]
    if w.empty:
        raise SystemExit("[supp] walker_thresholds.csv carries no Walker-basis rows")

    def row(r, gene_label):
        ok = r.status == "ok"
        notes = []
        if r.stratum == "pm12":
            notes.append(f"Out of scope: the recommendation excludes ±1,2 variants, "
                         f"which go through the {route}")
        elif r.stratum == "all_1_50":
            notes.append("Out of scope: the pool includes the ±1,2 stratum")
        if not ok:
            notes.append(f"Not evaluable: {r.reason}")
        pp3_cell = bp4_cell = ""
        if ok:
            pp3_cell = (fmt_ci(r.lr_pp3, r.lr_pp3_lo, r.lr_pp3_hi)
                        if _finite(r.lr_pp3) else "not evaluable")
            bp4_cell = (fmt_ci(r.lr_bp4, r.lr_bp4_lo, r.lr_bp4_hi)
                        if _finite(r.lr_bp4) else "not evaluable")
            if r.scope == "gene":
                pp3_cell = fmt_lr(r.lr_pp3) if _finite(r.lr_pp3) else "not evaluable"
                bp4_cell = fmt_lr(r.lr_bp4) if _finite(r.lr_bp4) else "not evaluable"
            if not _finite(r.lr_bp4) and _finite(r.n_in_bp4_band) \
                    and r.n_in_bp4_band < K.MIN_BAND:
                notes.append(f"BP4 band holds {_n_variants(r.n_in_bp4_band)}, "
                             f"fewer than {K.MIN_BAND}")
            if _finite(r.lr_bp4) and float(r.lr_bp4) == 0 and _finite(r.lr_bp4_for_tier):
                notes.append("No damaging variant in the BP4 band; the tier is read "
                             "off the rule-of-three upper bound "
                             f"{fmt_lr(r.lr_bp4_for_tier)}")
        if r.scope == "pooled" and ok:
            basis = (f"Gene-clustered bootstrap over {fmt_int(r.n_genes)} genes; "
                     f"draws evaluable: PP3 LR {fmt_int(r.lr_pp3_n_boot_used)} of "
                     f"{WT.N_BOOT}, BP4 LR {fmt_int(r.lr_bp4_n_boot_used)} of {WT.N_BOOT}")
        else:
            basis = _text(r.ci_basis).capitalize()
        return {
            "Stratum (intronic offset, bp)": STRATUM_LABEL[r.stratum],
            "In scope of the recommendation": yes_no(r.walker_in_scope),
            "ClinVar record status": ARM_LABEL[r.clinvar_arm],
            "Gene": gene_label,
            "Damaging": fmt_int(r.n_pos), "Normal": fmt_int(r.n_neg),
            f"Fraction of labelled scoring ≥ {pp3}": fmt_dp(r.frac_pp3) if ok else "",
            "Sensitivity": fmt_dp(r.sens_pp3) if ok else "",
            "Specificity": fmt_dp(r.spec_pp3) if ok else "",
            f"PP3 band (≥ {pp3}) LR (95% CI)": pp3_cell,
            "PP3 tier": tier_text(r.tier_pp3) if ok else "not evaluable",
            f"BP4 band (≤ {bp4}) LR (95% CI)": bp4_cell,
            "BP4 tier": tier_text(r.tier_bp4) if ok else "not evaluable",
            "CI basis": basis, "Note": "; ".join(notes),
        }

    rows = []
    for st in K.STRATA:
        pooled = w[(w.stratum == st) & (w.scope == "pooled")]
        for arm in K.ARMS:
            for r in pooled[pooled.clinvar_arm == arm].itertuples():
                rows.append(row(r, "All genes (pooled)"))
        if st in IN_SCOPE:
            genes = w[(w.stratum == st) & (w.scope == "gene") & (w.clinvar_arm == "all")]
            for r in genes.sort_values("gene").itertuples():
                rows.append(row(r, r.gene))
    return pd.DataFrame(rows), [REPORTS / "walker_thresholds.csv", WALKER_CFG]


# ---------------------------------------------------------------------------
# S6 / S7 -- interval-calibration thresholds
# ---------------------------------------------------------------------------
def _evidence() -> pd.DataFrame:
    ev = _csv(REPORTS / "evidence_thresholds.csv")
    return ev[(ev.tool != K.FUSION) & ev.stratum.isin(IN_SCOPE) & ev.tier.isin(TIERS)]


def _fold_string(f: pd.DataFrame) -> str:
    parts = []
    for r in f.sort_values("heldout_gene").itertuples():
        if _finite(r.heldout_lr):
            parts.append(f"{r.heldout_gene}={fmt_lr(r.heldout_lr)}")
        elif not _finite(r.threshold_from_6_genes):
            parts.append(f"{r.heldout_gene}=no threshold")
        else:
            parts.append(f"{r.heldout_gene}=n.e.")
    return "; ".join(parts)


def s6() -> tuple[pd.DataFrame, list[Path]]:
    ev = _evidence().set_index(["tool", "stratum", "tier"])
    tl = _csv(REPORTS / "tier_logo_table.csv").set_index(["tool", "stratum", "tier"])
    folds = _csv(REPORTS / "evidence_thresholds_logo_folds.csv")
    folds = folds[folds.side == "pp3"]
    rows = []
    for tool in TOOLS:
        for st in IN_SCOPE:
            for tier in TIERS:
                key = (tool, st, tier)
                e, t = ev.loc[key], tl.loc[key]
                base = {"Tool": _tool(tool), "Stratum (intronic offset, bp)":
                        STRATUM_LABEL[st], "Tier": tier.capitalize(),
                        "Tier LR cut": fmt_lr(e.path_cut_lr)}
                if e.status != "ok" or t.status != "ok":
                    rows.append(base | {"Reached in-sample": "not evaluable"})
                    continue
                f = folds[(folds.tool == tool) & (folds.stratum == st)
                          & (folds.tier == tier)]
                lrs = f["heldout_lr"].astype(float)
                ev_lrs = lrs.dropna()
                med, lo, hi = ((ev_lrs.median(), ev_lrs.min(), ev_lrs.max())
                               if len(ev_lrs) else (np.nan, np.nan, np.nan))
                cut = float(e.path_cut_lr)
                # the two upstream tables and the per-fold file must tell one story
                same = [
                    (bool(e.pp3_threshold_reachable), bool(t.in_sample_tier_reached)),
                    (int(e.pp3_logo_folds_reaching), int(t.folds_reached)),
                    (int(e.pp3_logo_folds_with_heldout_lr), int(t.folds_with_heldout_lr)),
                    (len(f), int(t.n_folds)),
                    (int(f.threshold_from_6_genes.notna().sum()), int(t.folds_reached)),
                    (int((lrs >= cut).sum()), int(t.folds_heldout_lr_above_cut)),
                ]
                num = [(e.pp3_threshold_insample, t.in_sample_threshold),
                       (med, t.heldout_lr_median),
                       (lo, t.heldout_lr_min), (hi, t.heldout_lr_max)]
                if any(a != b for a, b in same) or not all(
                        (pd.isna(a) and pd.isna(b)) or np.isclose(a, b, rtol=0, atol=1e-12)
                        for a, b in num):
                    raise SystemExit(f"[supp] S6 {key}: tier_logo_table.csv, "
                                     "evidence_thresholds.csv and the per-fold file "
                                     "disagree; rerun the upstream stages")
                rows.append(base | {
                    "Reached in-sample": yes_no(e.pp3_threshold_reachable),
                    "In-sample threshold (gene-clustered bootstrap)":
                        fmt_thr(e.pp3_threshold_insample, tool),
                    "In-sample threshold (variant-level bootstrap)":
                        fmt_thr(e.pp3_threshold_variant_boot, tool),
                    "Held-out genes": str(len(f)),
                    "Folds with a fitted threshold": str(int(t.folds_reached)),
                    "Folds with an evaluable held-out LR": str(int(t.folds_with_heldout_lr)),
                    "Held-out genes clearing the cut": str(int(t.folds_heldout_lr_above_cut)),
                    "Held-out LR, median": fmt_lr(med),
                    "Held-out LR, minimum": fmt_lr(lo),
                    "Held-out LR, maximum": fmt_lr(hi),
                    "Held-out LR per fold": _fold_string(f),
                })
    srcs = [REPORTS / "evidence_thresholds.csv", REPORTS / "tier_logo_table.csv",
            REPORTS / "evidence_thresholds_logo_folds.csv"]
    return pd.DataFrame(rows), srcs


def s7() -> tuple[pd.DataFrame, list[Path]]:
    ev = _evidence().set_index(["tool", "stratum", "tier"])
    rows = []
    for tool in TOOLS:
        for st in IN_SCOPE:
            for tier in TIERS:
                e = ev.loc[(tool, st, tier)]
                base = {"Tool": _tool(tool), "Stratum (intronic offset, bp)":
                        STRATUM_LABEL[st], "Tier": tier.capitalize(),
                        "Benign tier LR cut": fmt_lr(e.ben_cut_lr)}
                if e.status != "ok":
                    rows.append(base | {"Reached in-sample": "not evaluable"})
                    continue
                rows.append(base | {
                    "Reached in-sample": yes_no(e.bp4_threshold_reachable),
                    "In-sample threshold (gene-clustered bootstrap)":
                        fmt_thr(e.bp4_threshold_insample, tool, "lower"),
                    "In-sample threshold (variant-level bootstrap)":
                        fmt_thr(e.bp4_threshold_variant_boot, tool, "lower"),
                    "Held-out LR, median over folds":
                        fmt_lr(e.bp4_logo_heldout_lr_median)
                        if bool(e.bp4_threshold_reachable) else "",
                })
    return pd.DataFrame(rows), [REPORTS / "evidence_thresholds.csv"]


# ---------------------------------------------------------------------------
# S8 -- thresholds refitted on subsets of training genes
# ---------------------------------------------------------------------------
def _pct(x) -> str:
    return f"{100 * float(x):.0f}%" if _finite(x) else ""


def s8() -> tuple[pd.DataFrame, list[Path]]:
    src = REPORTS / "threshold_dilution_summary.csv"
    d = _csv(src)
    rows = []
    order = {t: i for i, t in enumerate(TOOLS)}
    d = d.assign(_t=d.tool.map(order), _s=d.stratum.map({s_: i for i, s_ in enumerate(IN_SCOPE)}),
                 _r=d.tier.map({"moderate": 0, "strong": 1}))
    for r in d.sort_values(["_t", "_s", "_r", "k"], kind="mergesort").itertuples():
        iqr = (f"{fmt_thr(r.threshold_q25, r.tool)} to {fmt_thr(r.threshold_q75, r.tool)}"
               if _finite(r.threshold_q25) else "")
        lr_iqr = (f"{fmt_lr(r.heldout_lr_q25)} to {fmt_lr(r.heldout_lr_q75)}"
                  if _finite(r.heldout_lr_q25) else "")
        rows.append({
            "Tool": _tool(r.tool), "Stratum (intronic offset, bp)": STRATUM_LABEL[r.stratum],
            "Tier": r.tier.capitalize(), "Training genes": fmt_int(r.k),
            "Subsets": fmt_int(r.n_subsets),
            "Subsets reaching the tier": f"{fmt_int(r.n_subsets_reached)} ({_pct(r.frac_subsets_reached)})",
            "Threshold, median": fmt_thr(r.threshold_median, r.tool),
            "Threshold, interquartile range": iqr,
            "Threshold on all genes in the stratum": fmt_thr(r.threshold_all_genes, r.tool),
            "Held-out pairs evaluable": f"{fmt_int(r.n_pairs_evaluable)} of {fmt_int(r.n_pairs)}",
            "Held-out pairs clearing the cut":
                f"{fmt_int(r.n_pairs_cleared)} ({_pct(r.frac_pairs_cleared)})"
                if int(r.n_pairs_evaluable) else "",
            "Held-out LR, median": fmt_lr(r.heldout_lr_median),
            "Held-out LR, interquartile range": lr_iqr,
        })
    return pd.DataFrame(rows), [src]


# ---------------------------------------------------------------------------
# S10 -- ClinVar record status
# ---------------------------------------------------------------------------
def s10() -> tuple[pd.DataFrame, list[Path]]:
    pp3 = fmt_const(_walker_cfg()["thresholds"]["pp3"]["value"])
    a = _csv(REPORTS / "arms_at_tool_threshold.csv")
    pub = a[a.threshold_basis == "published cut point"]
    if set(pub.tool) - set(WALKER_TOOLS):
        raise SystemExit("[supp] a published-cut-point row for a non-SpliceAI column")
    fit = a[a.threshold_basis.str.startswith("fitted ")]
    key = ["gene_set", "stratum", "clinvar_arm", "tool"]
    if fit.duplicated(key).any() or pub.duplicated(key).any():
        raise SystemExit("[supp] arms_at_tool_threshold.csv has duplicate cells")
    fit, pub = fit.set_index(key), pub.set_index(key)
    # only the AUROC columns: the walker_* columns of this file are not read
    t = _csv(REPORTS / "territory_metrics_arms_noBRCA1.csv",
             usecols=key + ["k_genes", "auroc", "auroc_lo", "auroc_hi", "status"])
    t = t.set_index(key)

    def auroc(k):
        if k not in t.index:
            return {"AUROC (95% CI)": "", "Genes pooled for AUROC": ""}
        r = t.loc[k]
        if r.status != "ok":
            return {"AUROC (95% CI)": "not evaluable",
                    "Genes pooled for AUROC": fmt_int(r.k_genes)}
        return {"AUROC (95% CI)": fmt_ci(r.auroc, r.auroc_lo, r.auroc_hi, fmt_dp),
                "Genes pooled for AUROC": fmt_int(r.k_genes)}

    def lr(r):
        if r.status != "ok" or not _finite(r.band_lr):
            return "not evaluable"
        return fmt_lr(r.band_lr)

    rows = []
    for gs in ("all_genes", "no_BRCA1"):
        for st in IN_SCOPE + ["all_1_50"]:
            arms = (["all"] if st != "all_1_50" else []) + ARMS
            for arm in arms:
                for tool in TOOLS:
                    k = (gs, st, arm, tool)
                    row = {"Genes": GENE_SET_LABEL[gs],
                           "Stratum (intronic offset, bp)": STRATUM_LABEL[st],
                           "ClinVar record status": ARM_LABEL[arm], "Tool": _tool(tool)}
                    band = {c: "" for c in (
                        "Fitted threshold tier", "Fitted threshold",
                        "Genes contributing", "Damaging", "Normal", "Variants in band",
                        "Sensitivity", "Specificity", "Band LR at fitted threshold",
                        f"Band LR at published cut point (≥ {pp3})")}
                    if st != "all_1_50":
                        f = fit.loc[k] if k in fit.index else None
                        p = pub.loc[k] if k in pub.index else None
                        ref = f if f is not None else p
                        if f is not None and p is not None and \
                                (f.n_pos, f.n_neg) != (p.n_pos, p.n_neg):
                            raise SystemExit(f"[supp] S10 {k}: counts differ between "
                                             "the fitted and published rows")
                        if ref is not None:
                            band.update({"Genes contributing": fmt_int(ref.n_genes),
                                         "Damaging": fmt_int(ref.n_pos),
                                         "Normal": fmt_int(ref.n_neg)})
                        if f is None:
                            band["Fitted threshold tier"] = "none reached in-sample"
                        else:
                            tier = f.threshold_basis.split(" ", 1)[1]
                            band.update({
                                "Fitted threshold tier": (
                                    "Moderate" if tier == "moderate"
                                    else f"{tier.capitalize()} (Moderate not reached)"),
                                "Fitted threshold": fmt_thr(f.threshold, tool),
                                "Variants in band": fmt_int(f.n_in_band),
                                "Sensitivity": fmt_dp(f.sens) if f.status == "ok" else "",
                                "Specificity": fmt_dp(f.spec) if f.status == "ok" else "",
                                "Band LR at fitted threshold": lr(f)})
                        if p is not None:
                            band[f"Band LR at published cut point (≥ {pp3})"] = lr(p)
                    rows.append(row | band | auroc(k))
    srcs = [REPORTS / "arms_at_tool_threshold.csv",
            REPORTS / "territory_metrics_arms_noBRCA1.csv", WALKER_CFG]
    return pd.DataFrame(rows), srcs


# ---------------------------------------------------------------------------
# S11 -- arms within gene
# ---------------------------------------------------------------------------
S11_STRATA = ["s3_10", "s3_50"]


def _within() -> pd.DataFrame:
    w = _csv(REPORTS / "arms_within_gene.csv")
    return w[w.stratum.isin(S11_STRATA)]


def _highest(aucs: dict) -> str:
    """The arm with the highest AUROC, on the unrounded values. Where the runner-up
    prints the same three decimals, the table says so rather than let two equal
    printed numbers carry an unexplained winner."""
    ranked = sorted(ARMS, key=lambda a: (-aucs[a], ARMS.index(a)))
    best = aucs[ranked[0]]
    top = [a for a in ARMS if aucs[a] == best]
    if len(top) > 1:
        return "tie: " + " / ".join(ARM_LABEL[a] for a in top)
    out = ARM_LABEL[top[0]]
    if fmt_dp(aucs[ranked[1]]) == fmt_dp(best):
        out += " (margin below the printed precision)"
    return out


def _three_arm_genes(w: pd.DataFrame, st: str) -> list[str]:
    ok = w[(w.stratum == st)].assign(ok=lambda d: d.status == "ok")
    per = ok.groupby(["gene", "tool"])["ok"].sum().unstack()
    full = per.eq(len(ARMS))
    if not full.all(axis=1).eq(full.any(axis=1)).all():
        raise SystemExit(f"[supp] S11 {st}: evaluable arms differ between columns")
    return sorted(full.index[full.all(axis=1)])


def s11() -> tuple[pd.DataFrame, list[Path]]:
    w = _within()
    idx = w.set_index(["stratum", "gene", "tool", "clinvar_arm"])
    rows = []
    for st in S11_STRATA:
        for gene in sorted(w.loc[w.stratum == st, "gene"].unique()):
            for tool in TOOLS:
                row = {"Stratum (intronic offset, bp)": STRATUM_LABEL[st],
                       "Gene": gene, "Tool": _tool(tool)}
                aucs = {}
                for arm in ARMS:
                    k = (st, gene, tool, arm)
                    if k in idx.index:
                        r = idx.loc[k]
                        row[f"{ARM_LABEL[arm]}: damaging/normal"] = \
                            f"{fmt_int(r.n_pos)}/{fmt_int(r.n_neg)}"
                        if r.status == "ok":
                            aucs[arm] = float(r.auroc)
                            row[f"{ARM_LABEL[arm]}: AUROC"] = fmt_dp(r.auroc)
                        else:
                            row[f"{ARM_LABEL[arm]}: AUROC"] = "not evaluable"
                    else:
                        row[f"{ARM_LABEL[arm]}: damaging/normal"] = "0/0"
                        row[f"{ARM_LABEL[arm]}: AUROC"] = "not evaluable"
                row["Arm with the highest AUROC"] = (_highest(aucs) if len(aucs) == len(ARMS)
                                                     else "")
                rows.append(row)
    return pd.DataFrame(rows), [REPORTS / "arms_within_gene.csv"]


def s11_summary() -> tuple[pd.DataFrame, list[Path]]:
    w = _within()
    rows = []
    for st in S11_STRATA:
        genes = _three_arm_genes(w, st)
        tot = {a: 0 for a in ARMS} | {"tie": 0, "n": 0}
        for gene in genes:
            c = {a: 0 for a in ARMS} | {"tie": 0, "n": 0}
            for tool in TOOLS:
                d = w[(w.stratum == st) & (w.gene == gene) & (w.tool == tool)]
                aucs = dict(zip(d.clinvar_arm, d.auroc.astype(float)))
                best = max(aucs.values())
                top = [a for a in ARMS if aucs[a] == best]
                c["tie" if len(top) > 1 else top[0]] += 1
                c["n"] += 1
            for k in tot:
                tot[k] += c[k]
            rows.append((st, gene, c))
        rows.append((st, f"All {len(genes)} genes with three evaluable arms "
                         f"({', '.join(genes)})", tot))
    out = []
    for st, gene, c in rows:
        top = max(ARMS, key=lambda a: (c[a], -ARMS.index(a)))
        tied = [a for a in ARMS if c[a] == c[top]]
        out.append({
            "Stratum (intronic offset, bp)": STRATUM_LABEL[st], "Gene": gene,
            "Columns compared": str(c["n"]),
            **{f"Highest in: {ARM_LABEL[a]}": str(c[a]) for a in ARMS},
            "Ties": str(c["tie"]),
            "Arm most often highest": (ARM_LABEL[top] if len(tied) == 1 else
                                       "tie: " + " / ".join(ARM_LABEL[a] for a in tied)),
        })
    return pd.DataFrame(out), [REPORTS / "arms_within_gene.csv"]


# ---------------------------------------------------------------------------
# S12 -- external genes
# ---------------------------------------------------------------------------
def s12() -> tuple[pd.DataFrame, list[Path]]:
    pp3 = fmt_const(_walker_cfg()["thresholds"]["pp3"]["value"])
    fig = _csv(REPORTS / "fig_data/fig4_external.csv")
    rows, srcs = [], []
    for gene in EXTERNAL:
        ep = REPORTS / f"external_{gene.lower()}.csv"
        bp = REPORTS / f"external_{gene.lower()}_column_basis.csv"
        srcs += [ep, bp]
        e = _csv(ep)
        basis = _csv(bp).set_index("tool")
        e = e[e.stratum.isin(IN_SCOPE) & (e.n > 0)]
        labels = list(dict.fromkeys(e.label_definition))
        for lab in labels:
            for st in IN_SCOPE:
                for tool in TOOLS:
                    sel = e[(e.label_definition == lab) & (e.stratum == st) & (e.tool == tool)]
                    if sel.empty:
                        continue
                    r = sel.iloc[0]
                    ok = r.status == "ok"
                    b = basis.loc[tool] if tool in basis.index else None
                    comparable = bool(b.comparable) if b is not None else False
                    notes = []
                    if not ok:
                        notes.append(f"Not evaluable: {r.reason}")
                    row = {"Gene": gene, "Label definition": LABEL_NAME[lab],
                           "Stratum (intronic offset, bp)": STRATUM_LABEL[st],
                           "Tool": _tool(tool),
                           "Damaging": fmt_int(r.n_pos), "Normal": fmt_int(r.n_neg),
                           "Evaluated": yes_no(ok),
                           "Same variable as the fitted column": yes_no(comparable),
                           "Why not comparable": ("" if comparable or b is None
                                                  else _plain_reason(b.reason))}
                    wl = f"LR at published cut point (≥ {pp3})"
                    row[wl] = row["Tier at published cut point"] = ""
                    if ok and tool in WALKER_TOOLS:
                        row[wl] = fmt_lr(r.walker_lr_pp3) if _finite(r.walker_lr_pp3) \
                            else "not evaluable"
                        row["Tier at published cut point"] = tier_text(r.walker_tier_pp3)
                    for tier in TIERS:
                        name = tier.capitalize()
                        tcol, lcol, gcol = (f"{name}: threshold (median of leave-one-gene-out fits)",
                                            f"{name}: LR in this gene",
                                            f"{name}: tier reached")
                        row[tcol] = row[lcol] = row[gcol] = ""
                        if not (ok and comparable):
                            continue
                        thr = r.get(f"e3_{tier}_threshold")
                        lr = r.get(f"e3_{tier}_lr_here")
                        note = r.get(f"e3_{tier}_note")
                        row[tcol] = fmt_thr(thr, tool)
                        row[lcol] = fmt_lr(lr) if _finite(lr) else (
                            "not evaluable" if _finite(thr) else "")
                        row[gcol] = tier_text(r.get(f"e3_{tier}_tier_here"))
                        if _text(note):
                            notes.append(f"{name}: {note}")
                        # the figure table is built from this file; it must agree
                        fr = fig[(fig.gene == gene) & (fig.label_definition == lab)
                                 & (fig.stratum == st) & (fig.tool == tool)
                                 & (fig.tier == tier)]
                        if len(fr) == 1:
                            f0 = fr.iloc[0]
                            pairs = [(thr, f0.logo_threshold), (lr, f0.lr_at_logo_threshold)]
                            if tool in WALKER_TOOLS:
                                pairs.append((r.walker_lr_pp3, f0.lr_at_walker_cut))
                            for x, y in pairs:
                                if not ((pd.isna(x) and pd.isna(y)) or
                                        np.isclose(float(x), float(y), rtol=1e-9, atol=0)):
                                    raise SystemExit(
                                        f"[supp] S12 {gene}/{lab}/{st}/{tool}/{tier}: "
                                        "fig4_external.csv disagrees with "
                                        f"{ep.name}; rerun evid_fig_data")
                    row["Note"] = "; ".join(notes)
                    rows.append(row)
    srcs.append(REPORTS / "fig_data/fig4_external.csv")
    return pd.DataFrame(rows), srcs


# ---------------------------------------------------------------------------
# S13 -- in-frame attribution
# ---------------------------------------------------------------------------
def s13() -> tuple[pd.DataFrame, list[Path]]:
    pp3 = fmt_const(_walker_cfg()["thresholds"]["pp3"]["value"])
    sets = [
        ("Seven genes", REPORTS / "inframe_attribution_summary.csv",
         {"spliceai": PHASE1 / INF.SUBSET, "pangolin": PHASE1 / INF.PANG_SUBSET},
         "assay normal"),
    ]
    for g in ("ddx3x", "tp53"):
        spec = INF.EXTERNAL[g]
        normal = ("normal under either label definition" if len(spec["label_cols"]) > 1
                  else "normal under the control-anchored definition")
        # the external gene's Pangolin events are read on SpliceAI's subset
        sets.append((g.upper(), REPORTS / f"inframe_attribution_{g}_summary.csv",
                     {"spliceai": PHASE1 / spec["subset"],
                      "pangolin": PHASE1 / spec["subset"]}, normal))
    rows, srcs = [], []
    for name, path, subset, normal in sets:
        sm = _csv(path)
        srcs.append(path)
        for tool in ("spliceai", "pangolin"):
            part = sm[sm.tool == tool]
            if part.empty:
                continue
            sub = pd.read_parquet(_need(subset[tool]), columns=["score_column"])
            srcs.append(subset[tool])
            col = sub["score_column"].unique()
            if len(col) != 1:
                raise SystemExit(f"[supp] {subset[tool].name}: mixed score columns")
            for st in K.STRATA:
                for r in part[part.stratum == st].itertuples():
                    total = int(r.in_frame + r.out_of_frame + r.undetermined)
                    if total != int(r.n):
                        raise SystemExit(f"[supp] {path.name}: counts do not add up")
                    resolved = int(r.in_frame + r.out_of_frame)
                    rows.append({
                        "Gene set": name,
                        "Events read from": _tool(tool) if tool != "spliceai"
                        else _tool(col[0]),
                        "Variant subset": f"{_tool(col[0])} ≥ {pp3} and {normal}",
                        "Stratum (intronic offset, bp)": STRATUM_LABEL[st],
                        "In-frame": str(int(r.in_frame)),
                        "Out-of-frame": str(int(r.out_of_frame)),
                        "Undetermined": str(int(r.undetermined)),
                        "Total": str(total),
                        "In-frame share": fmt_dp(r.in_frame / total) if total else "",
                        "In-frame share among resolved":
                            fmt_dp(r.in_frame / resolved) if resolved else "",
                    })
    return pd.DataFrame(rows), list(dict.fromkeys(srcs))


# ---------------------------------------------------------------------------
# the table list, in the fixed numbering
# ---------------------------------------------------------------------------
TABLES = [
    ("S1", "tableS1_functional_standards.csv", s1,
     "Functional standards: assay, accession, transcript, published class labels "
     "and how each maps to damaging, normal or missing, with variant counts",
     "External genes publish no classification; their constructed definitions are "
     "stated as implemented, and counts are for the primary (control-anchored) "
     "definition."),
    ("S2", "tableS2_analysis_set.csv", s2,
     "Analysis-set composition by gene, intronic-offset stratum and ClinVar record "
     "status",
     "ClinVar status is the strict split, which separates records carrying no "
     "clinical assertion. The ±1,2 stratum is outside the scope of the splicing "
     "recommendation and is shown for completeness."),
    ("S3", "tableS3_predictor_columns.csv", s3,
     "Predictor score columns: definition, version and scoring configuration, "
     "source, access date, licence, orientation, coverage and training signal",
     "Spearman rho is against functional pathogenicity on the whole analysis set. "
     "Every column is positively oriented as stored, including GPN-MSA, whose sign "
     "was flipped when it was scored, so no column is flipped here."),
    ("S4", "tableS4_column_agreement.csv", s4,
     "Agreement between related score columns",
     "The cut point applies only to the two SpliceAI columns."),
    ("S5", "tableS5_published_cut_points.csv", s5,
     "The published SpliceAI cut points on the Walker-basis column, by stratum and "
     "ClinVar record status, with per-gene point estimates",
     "Band likelihood ratios are point estimates; pooled intervals are 95% "
     "percentile intervals of a gene-clustered bootstrap; per-gene rows carry no "
     "interval. The BP4 tier is read off the rule-of-three bound where the band "
     "holds no damaging variant."),
    ("S6", "tableS6_thresholds_pathogenic.csv", s6,
     "Interval-calibration thresholds, pathogenic side: tier reached in-sample, "
     "thresholds under the gene-clustered and variant-level bootstraps, and "
     "leave-one-gene-out transfer",
     "Held-out likelihood ratios are point estimates in the held-out gene at the "
     "threshold fitted on the other genes. n.e. = not evaluable: the training genes "
     "did not reach the tier, fewer than ten of the held-out gene's variants scored "
     "above the threshold, or the held-out gene had fewer than ten damaging or ten "
     "normal variants. Thresholds are printed at four significant figures, or more "
     "where four would move a variant across the threshold; AlphaGenome's lie just "
     "below its ceiling of 2.2. The fusion is excluded."),
    ("S7", "tableS7_thresholds_benign.csv", s7,
     "Interval-calibration thresholds, benign side",
     "A held-out likelihood ratio of 0 means the band held no damaging variant in "
     "the held-out gene; no tier is read off it. The held-out median is left blank "
     "for a tier not reached in-sample. The fusion is excluded."),
    ("S8", "tableS8_dilution.csv", s8,
     "Thresholds refitted on every subset of training genes: share of subsets "
     "reaching each tier, the threshold across subsets, and the share of held-out "
     "genes clearing the cut",
     "Four splice-aware columns, Moderate and Strong tiers, every training subset of "
     "two or more genes. A held-out pair is one training subset and one gene left "
     "out of it. With two or three training genes the gene-clustered bootstrap gives "
     "wide intervals, so a low share of subsets reaching Strong partly reflects that "
     "width. Thresholds are printed at four significant figures, or more where four "
     "would move a variant across the threshold."),
    ("S10", "tableS10_clinvar_status.csv", s10,
     "Band likelihood ratios and AUROC by ClinVar record status, with and without "
     "BRCA1",
     "Each column's threshold is its own in-sample Moderate threshold, or "
     "Supporting where Moderate is not reached, fitted on all arms together. The "
     "published cut point is reported for the two SpliceAI columns only. Rows "
     "without BRCA1 use the thresholds fitted with BRCA1 included. AUROC is "
     "pooled across genes by random-effects meta-analysis; it is not computed for "
     "the all-variants arm."),
    ("S11", "tableS11_within_gene_arms.csv", s11,
     "Within-gene comparison of ClinVar record status arms: AUROC per gene, arm "
     "and column",
     "Comparing arms within a gene removes the confound of gene composition, which "
     "differs between arms."),
    ("S11", "tableS11_within_gene_arms_summary.csv", s11_summary,
     "Within-gene comparison, summary: the arm with the highest AUROC, per gene, "
     "over all columns",
     "Genes contributing all three arms only."),
    ("S12", "tableS12_external_genes.csv", s12,
     "External genes: the published cut point and the seven-gene fitted thresholds "
     "applied to DDX3X and TP53",
     "In-scope strata only. Fitted thresholds are the median of the "
     "leave-one-gene-out (LOGO) fits, carried over only where the tier is reached "
     "in-sample and the external column is the same variable. Label definitions "
     "are stated in Table S1."),
    ("S13", "tableS13_inframe_attribution.csv", s13,
     "In-frame attribution of variants scored at or above the PP3 cut point and "
     "called normal by the assay",
     "An inference from each tool's own predicted events, not a transcript "
     "measurement. Shares are recomputed from the counts. The PP3 cut point of 0.2 "
     "is calibrated for SpliceAI only; the Pangolin rows use the same score as a "
     "comparison device, so that both tools' false positives are defined alike."),
]


def _check(df: pd.DataFrame, name: str) -> None:
    if df.empty:
        raise SystemExit(f"[supp] {name} is empty")
    bad = df.isin(["nan", "NaN", "None", "<NA>"]).any().any() or df.isna().any().any()
    if bad:
        raise SystemExit(f"[supp] {name} carries a missing-value token")


def write_all(out_dir: Path = OUT_DIR) -> pd.DataFrame:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = []
    for table, fname, build, title, notes in TABLES:
        df, srcs = build()
        df = df.fillna("").astype(str)
        _check(df, fname)
        path = out_dir / fname
        df.to_csv(path, index=False, lineterminator="\n", encoding="utf-8")
        srcs = list(dict.fromkeys(Path(s) for s in srcs))
        data_srcs = [s_ for s_ in srcs if s_.suffix != ".py"]
        manifest.append({
            "Table": table, "File": fname, "Title": title,
            "Source files": "; ".join(_display(s) for s in srcs),
            "Source sha256": "; ".join(_sha256(s_) for s_ in data_srcs),
            "Rows": str(len(df)), "File sha256": _sha256(path), "Notes": notes,
        })
    # S8 and S9 are written by their own stages into the same directory
    for table, fname, stage, title in (
            ("S9", "tableS9_fusion.csv", "src/evid_fusion_stability.py",
             "Elastic-net combination: tier, per-fold thresholds and coefficient stability"),):
        path = OUT_DIR / fname
        if path.exists():
            manifest.append({"Table": table, "File": fname, "Title": title,
                             "Source files": f"written by {stage}", "Source sha256": "",
                             "Rows": str(len(pd.read_csv(path))),
                             "File sha256": _sha256(path),
                             "Notes": "Written by its own stage; indexed here."})
    order = {f"S{i}": i for i in range(1, 20)}
    man = pd.DataFrame(manifest).sort_values("Table", key=lambda c: c.map(order),
                                              kind="stable")
    man.to_csv(out_dir / MANIFEST, index=False, lineterminator="\n", encoding="utf-8")
    return man


def main() -> None:
    ap = argparse.ArgumentParser(description="Write the supplementary tables.")
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR,
                    help="where to write (default reports/evidence/supplement); the "
                         "tests use it to rerun into a scratch directory")
    out = ap.parse_args().out_dir
    man = write_all(out)
    print(f"[supp] wrote {len(man)} tables and {MANIFEST} to {out}")
    for _, r in man.iterrows():
        print(f"  {r['Table']:4s} {r['File']:42s} {int(r['Rows']):>5} rows  "
              f"sha256 {r['File sha256']}")


if __name__ == "__main__":
    main()
