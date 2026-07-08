"""
Phase 1 -- build and freeze the raw feature matrix for the splice-fusion study.

Steps (each is a pure function so you can test/re-run in isolation):
  1. load_raw ................ read your matrix, rename to canonical columns
  2. orient_functional ....... define func_pathogenicity = -func_score (with per-gene flips)
  3. derive_intron_offset .... parse HGVS "+N/-N" (fallback if no precomputed offset)
  4. assign_region .......... splice-core / splice-region / intronic / non-splice
  5. classes_from_clinvar .... Method A binary labels (clean_PB)
  6. classes_from_assay ...... Method B binary labels (per-gene 2-component mixture)
  7. add_missingness_flags ... one is-missing indicator per predictor
  8. audit + freeze .......... coverage tables, class balance, sha256 manifest, parquet

Run:  python -m src.phase1_build_frozen_matrix
"""
from __future__ import annotations

import hashlib
import json
import re
import subprocess
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture

from . import config as C


# ---------------------------------------------------------------------------
# 1. Load + canonicalise
# ---------------------------------------------------------------------------
def load_raw() -> pd.DataFrame:
    p = C.RAW_MATRIX_PATH
    if not p.exists():
        sys.exit(f"[phase1] raw matrix not found at {p} -- edit RAW_MATRIX_PATH in config.py")
    df = pd.read_parquet(p) if p.suffix == ".parquet" else pd.read_csv(p)

    rename = {src: canon for canon, src in C.COLUMNS.items() if src is not None}
    missing = [src for src in rename if src not in df.columns]
    if missing:
        sys.exit(f"[phase1] these mapped columns are absent from your file: {missing}\n"
                 f"        -> fix the right-hand side of COLUMNS in config.py")
    df = df.rename(columns=rename)[list(rename.values())].copy()

    # keep only the seven target genes, in a fixed category order (stable sorting)
    df = df[df["gene"].isin(C.GENES)].copy()
    df["gene"] = pd.Categorical(df["gene"], categories=C.GENES, ordered=True)
    return df


# ---------------------------------------------------------------------------
# 2. Orient the functional score
# ---------------------------------------------------------------------------
def orient_functional(df: pd.DataFrame) -> pd.DataFrame:
    s = df["func_score"].astype(float).copy()
    if C.FLIP_GENES:
        flip = df["gene"].isin(C.FLIP_GENES)
        s = s.where(~flip, -s)
    # pathogenicity: larger = more damaging
    df["func_pathogenicity"] = -s
    return df


# ---------------------------------------------------------------------------
# 3. Intron offset  (prefer a precomputed value; parse HGVS only as fallback)
# ---------------------------------------------------------------------------
_HGVS_OFFSET = re.compile(r"[+-]\d+")

def _offset_from_hgvs(hgvs: str) -> float:
    """Return signed intron offset from coding HGVS, e.g. 'c.4096+3G>A' -> 3,
    'c.4097-2A>G' -> -2. Exonic variants (no +/-) return NaN."""
    if not isinstance(hgvs, str):
        return np.nan
    m = _HGVS_OFFSET.search(hgvs)
    return float(m.group()) if m else np.nan

def derive_intron_offset(df: pd.DataFrame) -> pd.DataFrame:
    if C.COLUMNS.get("intron_offset") and "intron_offset" in df.columns:
        df["intron_offset"] = pd.to_numeric(df["intron_offset"], errors="coerce")
        df["offset_source"] = np.where(df["intron_offset"].notna(), "upstream", "missing")
    elif "hgvs_c" in df.columns:
        df["intron_offset"] = df["hgvs_c"].map(_offset_from_hgvs)
        df["offset_source"] = np.where(df["intron_offset"].notna(), "hgvs", "missing")
    else:
        df["intron_offset"] = np.nan
        df["offset_source"] = "missing"

    # reproducible manual corrections for audited edge cases
    overrides = getattr(C, "OFFSET_OVERRIDES", {})
    if overrides:
        m = df["variant_id"].isin(overrides)
        df.loc[m, "intron_offset"] = df.loc[m, "variant_id"].map(overrides)
        df.loc[m, "offset_source"] = "override"
    # donor (+) vs acceptor (-) side; NaN for exonic
    df["splice_side"] = np.where(df["intron_offset"] > 0, "donor",
                          np.where(df["intron_offset"] < 0, "acceptor", pd.NA))
    return df


# ---------------------------------------------------------------------------
# 4. Region assignment
# ---------------------------------------------------------------------------
def assign_region(df: pd.DataFrame) -> pd.DataFrame:
    """Two granularities, decoupled:
       - `region` / `is_splice`  : broad membership. From the authoritative
         region_label if provided (so the splice count matches the benchmark's),
         else derived from offset.
       - `splice_subclass`       : fine core/region/intronic tier, ALWAYS derived
         from intron_offset, populated only within the splice set."""
    if C.COLUMNS.get("region_label"):
        df["region"] = df["region_label"].astype(str)
        r = df["region"].str.lower().str.strip()
        # contains "splice" but is not a "non_splice"/"non-splice" style negation
        df["is_splice"] = r.str.contains("splice", na=False) & ~r.str.startswith("non")
    else:
        d0 = df["intron_offset"].abs()
        df["region"] = np.where(d0.notna(), "splice", "non_splice")
        df["is_splice"] = d0.notna()

    d = df["intron_offset"].abs()
    sub = pd.Series(pd.NA, index=df.index, dtype=object)
    sub[df["is_splice"] & d.le(C.SPLICE_CORE_MAX)] = "splice_core"
    sub[df["is_splice"] & d.gt(C.SPLICE_CORE_MAX) & d.le(C.SPLICE_REGION_MAX)] = "splice_region"
    sub[df["is_splice"] & d.gt(C.SPLICE_REGION_MAX)] = "intronic"
    # genuine exon-side splice variants: no intron offset exists -> explicit label
    edge = df["variant_id"].isin(getattr(C, "EXON_EDGE_VARIANTS", []))
    sub[df["is_splice"] & edge] = "splice_exon_edge"
    # any remaining splice variant with no offset and not on the edge list
    sub[df["is_splice"] & d.isna() & sub.isna()] = "splice_nooffset"
    df["splice_subclass"] = sub

    # coarse 2-tier for THIS dataset (intronic tier is empty; edge folds into core-like)
    core_like = df["splice_subclass"].isin(["splice_core", "splice_exon_edge"])
    df["splice_bin"] = np.where(df["is_splice"],
                        np.where(core_like, "core_like", "region"), pd.NA)

    # analysis_region = fine tier within splice, broad label elsewhere (used by audits)
    df["analysis_region"] = np.where(df["is_splice"],
                                     df["splice_subclass"].astype(object),
                                     df["region"].astype(object))
    return df


# ---------------------------------------------------------------------------
# 5. Method A -- ClinVar clean_PB binary labels
# ---------------------------------------------------------------------------
_PATH = ("pathogenic",)          # matches 'Pathogenic' and 'Likely_pathogenic'
_BENIGN = ("benign",)            # matches 'Benign' and 'Likely_benign'

def classes_from_clinvar(df: pd.DataFrame) -> pd.DataFrame:
    mode = getattr(C, "CLINVAR_MODE", "text")
    label = pd.Series(pd.NA, index=df.index, dtype="Float64")

    if mode == "enum":
        # exact-match an enum column (e.g. "P/LP" -> 1, "B/LB" -> 0); unmatched -> NA
        mapped = df["clinvar"].map(C.CLINVAR_ENUM_MAP)
        label = mapped.astype("Float64")
    elif mode == "binary":
        # column is already a clean 0/1 (P=1, B=0, VUS=NA)
        label = pd.to_numeric(df["clinvar"], errors="coerce").astype("Float64")
    else:  # "text"
        s = df["clinvar"].astype(str).str.lower()
        is_conflict = s.str.contains("conflict")
        is_path = s.str.contains("|".join(_PATH)) & ~is_conflict
        is_ben  = s.str.contains("|".join(_BENIGN)) & ~is_conflict
        label[is_path] = 1.0
        label[is_ben]  = 0.0

    df["y_clinvar"] = label
    return df


# ---------------------------------------------------------------------------
# 6. Method B -- assay-intrinsic binary labels via per-gene 2-component mixture
# ---------------------------------------------------------------------------
def classes_from_assay(df: pd.DataFrame) -> pd.DataFrame:
    """Per-gene Gaussian mixture on func_pathogenicity: the component with the
    higher mean = 'damaging'. Emits a hard label + a posterior probability.
    Where the original assay paper gives official control-based thresholds,
    OVERRIDE this column -- see the checklist."""
    y = pd.Series(pd.NA, index=df.index, dtype="Float64")
    post = pd.Series(np.nan, index=df.index, dtype=float)

    for gene, idx in df.groupby("gene", observed=True).groups.items():
        vals = df.loc[idx, "func_pathogenicity"].astype(float)
        ok = vals.dropna()
        if len(ok) < C.MIN_N_FOR_MIXTURE:
            continue
        gm = GaussianMixture(n_components=2, random_state=C.RANDOM_SEED, n_init=5)
        gm.fit(ok.values.reshape(-1, 1))
        dmg_comp = int(np.argmax(gm.means_.ravel()))     # higher mean = more damaging
        p_dmg = gm.predict_proba(ok.values.reshape(-1, 1))[:, dmg_comp]
        post.loc[ok.index] = p_dmg
        y.loc[ok.index] = (p_dmg >= 0.5).astype(float)

    df["y_assay"] = y
    df["y_assay_posterior"] = post
    return df


# ---------------------------------------------------------------------------
# 7. Missingness flags  (safe to precompute: no cross-gene statistic)
# ---------------------------------------------------------------------------
def add_missingness_flags(df: pd.DataFrame) -> pd.DataFrame:
    present = [f for f in (C.PREDICTOR_FEATURES + C.EVO_FEATURES)
               if f in df.columns and C.COLUMNS.get(f)]
    for f in dict.fromkeys(present):            # dedupe, preserve order
        df[f"{f}_isna"] = df[f].isna().astype("int8")
    return df


# ---------------------------------------------------------------------------
# 8a. Audits
# ---------------------------------------------------------------------------
def coverage_by_region(df: pd.DataFrame) -> pd.DataFrame:
    feats = [f for f in dict.fromkeys(C.PREDICTOR_FEATURES + C.EVO_FEATURES)
             if f in df.columns and C.COLUMNS.get(f)]
    rows = []
    for region, g in df.groupby("analysis_region", observed=True):
        for f in feats:
            rows.append({"region": region, "predictor": f,
                         "n": len(g), "n_scored": int(g[f].notna().sum()),
                         "coverage": round(g[f].notna().mean(), 4)})
    return pd.DataFrame(rows).pivot(index="predictor", columns="region",
                                    values="coverage").fillna(0.0)

def class_balance(df: pd.DataFrame) -> pd.DataFrame:
    out = []
    for label_col in ["y_clinvar", "y_assay"]:
        for drop_brca1 in (False, True):
            sub = df if not drop_brca1 else df[df["gene"] != "BRCA1"]
            splice = sub[sub["is_splice"]]
            lab = splice[label_col].dropna()
            out.append({
                "label": label_col,
                "brca1": "excluded" if drop_brca1 else "included",
                "n_labelled": int(lab.shape[0]),
                "n_pos": int((lab == 1).sum()),
                "n_neg": int((lab == 0).sum()),
                "pos_frac": round(float((lab == 1).mean()), 4) if len(lab) else np.nan,
            })
    return pd.DataFrame(out)

def subset_counts(df: pd.DataFrame) -> pd.DataFrame:
    return (df.groupby(["gene", "analysis_region"], observed=True).size()
              .rename("n").reset_index()
              .pivot(index="gene", columns="analysis_region", values="n")
              .fillna(0).astype(int))


# ---------------------------------------------------------------------------
# 8b. Freeze
# ---------------------------------------------------------------------------
def _git_commit() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"],
                                       stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return "no-git"

def freeze(df: pd.DataFrame) -> dict:
    C.OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    df = df.sort_values("variant_id").reset_index(drop=True)

    # content hash over a canonical CSV serialisation (stable across machines)
    canonical = df.to_csv(index=False).encode()
    sha = hashlib.sha256(canonical).hexdigest()

    out_path = C.OUTPUT_DIR / f"frozen_matrix_{C.FROZEN_VERSION}.parquet"
    df.to_parquet(out_path, index=False)

    manifest = {
        "version": C.FROZEN_VERSION,
        "frozen_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "git_commit": _git_commit(),
        "source_path": str(C.RAW_MATRIX_PATH),
        "n_rows": int(df.shape[0]),
        "n_genes": int(df["gene"].nunique()),
        "n_splice": int(df["is_splice"].sum()),
        "sha256": sha,
        "output_path": str(out_path),
    }
    (C.OUTPUT_DIR / f"manifest_{C.FROZEN_VERSION}.json").write_text(json.dumps(manifest, indent=2))
    return manifest


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------
def build() -> tuple[pd.DataFrame, dict]:
    df = load_raw()
    df = orient_functional(df)
    df = derive_intron_offset(df)
    df = assign_region(df)
    df = classes_from_clinvar(df)

    C.REPORT_DIR.mkdir(parents=True, exist_ok=True)
    # Method B: official published labels (default) or the old GMM placeholder
    if getattr(C, "METHOD_B", "official") == "official":
        from .phase1_method_b_labels import ingest   # lazy import avoids circularity
        df, cov = ingest(df)
        cov.to_csv(C.REPORT_DIR / "methodB_coverage.csv", index=False)
    else:
        df = classes_from_assay(df)

    df = add_missingness_flags(df)

    coverage_by_region(df).to_csv(C.REPORT_DIR / "coverage_by_region.csv")
    class_balance(df).to_csv(C.REPORT_DIR / "class_balance.csv", index=False)
    subset_counts(df).to_csv(C.REPORT_DIR / "subset_counts.csv")

    manifest = freeze(df)   # single, correct manifest for the FINAL matrix
    return df, manifest


def main() -> None:
    df, manifest = build()
    print("[phase1] frozen matrix written")
    print(json.dumps(manifest, indent=2))
    print("\n[phase1] splice-subset subclass counts:")
    print(df[df["is_splice"]].groupby("splice_subclass", observed=True).size().to_string())
    print(f"\n[phase1] reports in: {C.REPORT_DIR}/")


if __name__ == "__main__":
    main()
