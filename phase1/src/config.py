"""
config.py -- Phase 1 configuration for the splice-variant fusion study.

This is the single source of truth. Everything downstream reads from here.
FILL IN every line marked  # >>> EDIT  to match your actual data, then run:

    python -m src.phase1_build_frozen_matrix

Nothing in Phase 1 imputes or fits any statistic that could leak across genes:
we only build, label, audit, and FREEZE the raw matrix. All fold-dependent
transforms (imputation, rank-normalisation) are deferred to Phase 2, where they
are fit inside each training fold. See NOTE in transforms section below.
"""

from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
# >>> EDIT: point these at your real files / folders
# Source matrix is score_matrix_final.tsv, converted to parquet (with an added
# variant_id key) because load_raw() reads parquet or comma-CSV, not tab-TSV.
# Relative (resolved from phase1/) so a clean clone rebuilds portably.
# Regenerate this file from the tracked upstream TSV with:  python -m src.make_raw
RAW_MATRIX_PATH = Path("data/raw/variant_scores.parquet")
OUTPUT_DIR      = Path("data/frozen")                       # frozen matrix + manifest go here
REPORT_DIR      = Path("reports/phase1")                    # diagnostic tables go here

FROZEN_VERSION  = "v1"  # bump this whenever the upstream benchmark scoring changes

# ---------------------------------------------------------------------------
# Column mapping -- map YOUR column names to the canonical names used downstream.
# Left  = canonical name (DO NOT change). Right = your raw column name.  # >>> EDIT
# Set a value to None if you don't have that column yet.
# ---------------------------------------------------------------------------
COLUMNS = {
    # identity / annotation
    "variant_id":   "variant_id",      # constructed as "{chrom}-{pos}-{ref}-{alt}" in the parquet
    "gene":         "gene",            # gene symbol
    "chrom":        "chrom",
    "pos":          "pos",
    "ref":          "ref",
    "alt":          "alt",
    "hgvs_c":       None,              # no separate coding-HGVS col; offset comes from precomputed intron_offset
    "hgvs_nt":      "hgvs_nt",         # transcript c.-notation (e.g. "NM_007294.3:c.5467+18A>T") -> Method B join key
    "consequence":  "mc_terms",        # VEP molecular consequence terms
    "region_label": "region_class",    # authoritative broad region -> defines is_splice
    "intron_offset": "intron_offset",  # precomputed signed offset (preferred over HGVS parse); None to derive from HGVS

    # gold standard + clinical
    "func_score":   "functional_score",  # continuous SGE/MAVE score, raw as published
    "clinvar":      "clnsig_class",       # ClinVar column (text, enum, or precomputed 0/1 -- see CLINVAR_MODE)

    # predictor scores (the 10). Set to None any you don't have yet.
    "spliceai":      "spliceai_ds",
    "pangolin":      "pangolin_score",
    "alphagenome":   "alphagenome_splice",
    "gpn_msa":       "gpn_msa_score",
    "nt":            "nucleotide_transformer",  # Nucleotide Transformer masked-token LLR
    "cadd":          "cadd_phred",
    "alphamissense": "alphamissense",
    "phylop":        "phylop100way",
    "phastcons":     "phastcons100way",
    "gnomad_af":     "gnomad_af_global",
}

# ---------------------------------------------------------------------------
# Feature groups -- kept separable so ablation in Phase 2 is clean.
# ---------------------------------------------------------------------------
PREDICTOR_FEATURES = ["spliceai", "pangolin", "alphagenome", "gpn_msa", "nt",
                      "cadd", "alphamissense", "gnomad_af"]
# Evolutionary / alignment axis -- this is the headline ablation group.
EVO_FEATURES       = ["phylop", "phastcons", "gpn_msa"]
# NOTE: gpn_msa is alignment-conditioned, so it legitimately sits in BOTH groups.
# In Phase 2 the evo-axis ablation must be run TWO ways: (a) dropping only the
# pure conservation features (phyloP/phastCons), and (b) dropping all three,
# so the alignment-vs-conservation contributions can be separated.

GENES = ["BRCA1", "BRCA2", "BARD1", "PALB2", "RAD51C", "VHL", "BAP1"]

# ---------------------------------------------------------------------------
# ClinVar (Method A) label source -- how to turn the `clinvar` column into y_clinvar
# ---------------------------------------------------------------------------
#   "enum"   -> exact-match an enum column using CLINVAR_ENUM_MAP (unmatched -> NA)
#   "binary" -> the column is already a clean 0/1 (P=1, B=0, VUS=NA); used as-is
#   "text"   -> substring-match free text for pathogenic/benign (original behaviour)
CLINVAR_MODE = "enum"
CLINVAR_ENUM_MAP = {"P/LP": 1.0, "B/LB": 0.0}   # everything else -> NA (VUS, Conflicting, ...)

# ---------------------------------------------------------------------------
# Functional-score orientation
# ---------------------------------------------------------------------------
# We define functional PATHOGENICITY = -func_score  (more damaging -> larger).
# If a gene's assay is oriented so HIGHER raw score = MORE damaging, list it here
# and it will be flipped so all genes share one convention.
FLIP_GENES = []  # >>> EDIT if needed, e.g. ["VHL"]

# ---------------------------------------------------------------------------
# Splice sub-classification (by |intron offset| to nearest splice site)
# ---------------------------------------------------------------------------
SPLICE_CORE_MAX   = 2   # |offset| <= 2 -> splice-core
SPLICE_REGION_MAX = 8   # 3..8 -> splice-region ; >8 -> intronic
# NOTE for THIS dataset: the functional standard covers the proximal window only;
# max |offset| = 8, so the intronic (>8) tier is empty. State this scope in the paper.

# ---------------------------------------------------------------------------
# Audited edge-case corrections (see Phase1_splice_nooffset_report). Applied
# reproducibly here rather than by hand-editing the frozen matrix.
# ---------------------------------------------------------------------------
# 6 UTR-intron variants (c.-19-1 / c.-19-2) whose upstream offset was NaN because
# the annotator missed the c.-N-M notation. True |offset| is 1/2 -> real core.
OFFSET_OVERRIDES = {
    "17-43124116-C-A": -1, "17-43124116-C-G": -1, "17-43124116-C-T": -1,  # BRCA1 c.-19-1
    "17-43124117-T-A": -2, "17-43124117-T-C": -2, "17-43124117-T-G": -2,  # BRCA1 c.-19-2
}
# 7 genuine exon-side splice variants (last exonic base, splice_donor) for which
# no intron offset exists. Kept in the splice set, labelled splice_exon_edge,
# grouped with core in the coarse tier. >>> confirm these match your variant_id format.
EXON_EDGE_VARIANTS = [
    "17-43104875-A-C", "17-43104875-A-G", "17-43104875-A-T",  # BRCA1 c.294
    "17-43106476-A-G",                                         # BRCA1 c.192
    "17-43115728-G-A",                                         # BRCA1 c.132
    "16-23641125-A-T", "16-23641126-C-T",                      # PALB2 c.33 / c.32
]

# ---------------------------------------------------------------------------
# Assay-intrinsic class definition (Method B)
# ---------------------------------------------------------------------------
# Minimum variants required in a gene to attempt a 2-component mixture split.
MIN_N_FOR_MIXTURE = 30
RANDOM_SEED = 20260708


# ---------------------------------------------------------------------------
# Method B -- official per-assay labels (replaces the placeholder GMM entirely).
# Join key is hgvs_nt (transcript c.-notation): assembly-independent, no liftover.
# For each gene: label file, its join column, its class column, and the collapse
# to the common binary (damaging=1 / normal=0 / ambiguous=NA).
# ---------------------------------------------------------------------------
from pathlib import Path as _Path
ASSAY_LABEL_DIR = _Path("data/assay_labels")   # one label file per gene

ASSAY_LABELS = {
    "BRCA1":  dict(file="brca1_findlay2018.tsv",  key="hgvs_nt", class_col="func.class",
                   map={"non-functional": 1.0, "functional": 0.0, "intermediate": None}),
    "BRCA2":  dict(file="brca2_huang2025.tsv",    key="hgvs_nt", class_col="coarse_class",
                   map={"abnormal": 1.0, "normal": 0.0, "intermediate": None}),
    "BARD1":  dict(file="bard1_woo2025.tsv",      key="hgvs_nt", class_col="class",
                   map={"functionally abnormal": 1.0, "normal": 0.0, "indeterminate": None}),
    "PALB2":  dict(file="palb2_mavedb.csv",       key="hgvs_nt", class_col="class",
                   map={"abnormal": 1.0, "normal": 0.0, "indeterminate": None}),
    "RAD51C": dict(file="rad51c_mavedb.csv",      key="hgvs_nt", class_col="functional_classification",
                   map={"Fast depleted": 1.0, "Slow depleted": 1.0, "Unchanged": 0.0, "Enriched": None}),
    "VHL":    dict(file="vhl_buckley2024.tsv",    key="hgvs_nt", class_col="class",
                   map={"LOF1": 1.0, "LOF2": 1.0, "neutral": 0.0, "intermediate": None}),
    "BAP1":   dict(file="bap1_waters2024.tsv",    key="hgvs_nt", class_col="class",
                   map={"strongly depleted": 1.0, "weakly depleted": 1.0,
                        "unchanged": 0.0, "enriched": None}),
}


# ---------------------------------------------------------------------------
# Method B mode: "official" (join published labels, the validated default) or
# "gmm" (the old unanchored placeholder, kept only for comparison/fallback).
# ---------------------------------------------------------------------------
METHOD_B = "official"
