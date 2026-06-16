"""Shared configuration: target genes, paths, ClinVar classification maps.

Keeping these in one place so every script agrees on gene set, file locations,
and how raw ClinVar fields collapse into clean categories. Edit the gene list
or add genes here only.
"""
from pathlib import Path

# --- Project paths -----------------------------------------------------------
ROOT = Path(__file__).resolve().parent.parent
DATA_RAW = ROOT / "data" / "raw"
DATA_INTERIM = ROOT / "data" / "interim"
DATA_PROCESSED = ROOT / "data" / "processed"
RESULTS_TABLES = ROOT / "results" / "tables"
RESULTS_FIGURES = ROOT / "results" / "figures"

for _p in (DATA_RAW, DATA_INTERIM, DATA_PROCESSED, RESULTS_TABLES, RESULTS_FIGURES):
    _p.mkdir(parents=True, exist_ok=True)

# --- Target genes ------------------------------------------------------------
# DNA-repair / hereditary breast-ovarian cancer genes with dense ClinVar
# annotation and SGE/MAVE functional gold standards.
TARGET_GENES = ["BRCA1", "BRCA2", "BARD1", "RAD51C", "PALB2", "ATM"]

# --- ClinVar data source (GRCh38) -------------------------------------------
CLINVAR_VCF_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
CLINVAR_VCF_GZ = DATA_RAW / "clinvar_GRCh38.vcf.gz"

# --- CLNREVSTAT -> review star rating ---------------------------------------
# ClinVar's 0-4 gold-star scale. We keep variants with >= 1 star for the
# primary clinical-label reference (the plan's "ClinVar >= 1 star" rule).
REVSTAT_STARS = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_single_submitter": 1,
    "criteria_provided,_conflicting_classifications": 1,
    "criteria_provided,_conflicting_interpretations": 1,
    "no_assertion_criteria_provided": 0,
    "no_assertion_provided": 0,
    "no_interpretation_for_the_single_variant": 0,
    "no_classification_for_the_single_variant": 0,
    "no_classifications_from_unflagged_records": 0,
}


def revstat_to_stars(revstat: str) -> int:
    """Map a CLNREVSTAT string to its 0-4 gold-star count (default 0)."""
    return REVSTAT_STARS.get(revstat, 0)


def classify_clnsig(clnsig: str) -> str:
    """Collapse a raw CLNSIG value into one of:
    P/LP, B/LB, VUS, Conflicting, Other.
    """
    if not clnsig:
        return "Other"
    s = clnsig.lower()
    # Conflicting first — it can co-occur with other tokens.
    if "conflicting" in s:
        return "Conflicting"
    if "pathogenic" in s and "benign" not in s:
        # Pathogenic, Likely_pathogenic, Pathogenic/Likely_pathogenic
        return "P/LP"
    if "benign" in s and "pathogenic" not in s:
        return "B/LB"
    if "uncertain_significance" in s or s == "uncertain_significance":
        return "VUS"
    return "Other"


# Sequence-Ontology consequence -> coarse variant class for stratified analysis.
# Used to find the splice / non-coding subset the project focuses on.
def classify_consequence(mc_terms) -> str:
    """Given a list of SO consequence terms (from the MC field), return a
    coarse class: missense / synonymous / nonsense / splice / noncoding / other.
    Priority order matters: splice beats missense if both are present.
    """
    if not mc_terms:
        return "unknown"
    terms = {t.lower() for t in mc_terms}
    if any("splice" in t for t in terms):
        return "splice"
    if any(t in terms for t in ("nonsense", "stop_gained")):
        return "nonsense"
    if "missense_variant" in terms:
        return "missense"
    if "synonymous_variant" in terms:
        return "synonymous"
    noncoding_markers = (
        "intron_variant", "5_prime_utr_variant", "3_prime_utr_variant",
        "upstream", "downstream", "non-coding", "non_coding",
        "genic_upstream", "genic_downstream", "intergenic",
    )
    if any(any(m in t for m in noncoding_markers) for t in terms):
        return "noncoding"
    if any(t in terms for t in ("frameshift_variant",)):
        return "frameshift"
    if any("inframe" in t for t in terms):
        return "inframe_indel"
    return "other"
