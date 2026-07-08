"""
Phase 1 -- Method B label ingestion (parallel track, needed before Phase 3).

Replaces the placeholder GMM y_assay with official published per-variant classes,
joined on transcript c.-notation (hgvs_nt) so assembly differences (BRCA1 hg19)
never matter. Produces:
  - y_assay        : official binary label (1 damaging / 0 normal / NA)
  - y_assay_source : "official" | "unlabelled"
and a per-gene join-coverage report.

Genes with no matched label stay NA -> those variants simply rely on Method A
(y_clinvar) and the continuous func_pathogenicity target. That is fine: Method B
only feeds Phase-3 calibration.

Run:  python -m src.phase1_method_b_labels
"""
from __future__ import annotations
import re
import sys
import numpy as np
import pandas as pd
from . import config as C

_CDOT = re.compile(r"c\.[^ \t]+$")

def _norm_hgvs(x) -> str | float:
    """Reduce a transcript HGVS to a comparable 'c.xxx' string:
    strip the accession prefix (NM_/ENST...:) and version, lowercase, no spaces."""
    if not isinstance(x, str):
        return np.nan
    s = x.strip()
    if ":" in s:
        s = s.split(":")[-1]           # drop 'NM_007294.3:' prefix
    m = _CDOT.search(s.replace(" ", ""))
    return m.group().lower() if m else s.replace(" ", "").lower()


def _load_gene_labels(gene: str, spec: dict) -> pd.DataFrame:
    path = C.ASSAY_LABEL_DIR / spec["file"]
    if not path.exists():
        return pd.DataFrame(columns=["_join", "y"])   # not fetched yet -> empty
    lab = pd.read_csv(path, sep=None, engine="python")
    if spec["key"] not in lab.columns or spec["class_col"] not in lab.columns:
        sys.exit(f"[methodB] {gene}: expected columns {spec['key']!r}/{spec['class_col']!r} "
                 f"not in {path.name}; columns present: {list(lab.columns)}")
    out = pd.DataFrame({
        "_join": lab[spec["key"]].map(_norm_hgvs),
        "y": lab[spec["class_col"]].map(spec["map"]),   # unmapped values -> NaN
    })
    return out.dropna(subset=["_join"]).drop_duplicates("_join")


def ingest(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    if "hgvs_nt" not in df.columns:
        sys.exit("[methodB] frozen matrix has no hgvs_nt column -- map it in config.COLUMNS")
    df = df.copy()
    df["_join"] = df["hgvs_nt"].map(_norm_hgvs)
    y = pd.Series(pd.NA, index=df.index, dtype="Float64")

    cov = []
    for gene, spec in C.ASSAY_LABELS.items():
        lab = _load_gene_labels(gene, spec)
        mask = df["gene"] == gene
        n_splice = int((mask & df["is_splice"]).sum())
        if lab.empty:
            cov.append(dict(gene=gene, label_file="MISSING", n_splice=n_splice,
                            n_labelled_splice=0, splice_coverage=0.0))
            continue
        lut = dict(zip(lab["_join"], lab["y"]))
        mapped = df.loc[mask, "_join"].map(lut)
        y.loc[mask] = mapped.astype("Float64")
        lab_splice = int((mask & df["is_splice"] & y.notna()).sum())
        cov.append(dict(gene=gene, label_file=spec["file"], n_splice=n_splice,
                        n_labelled_splice=lab_splice,
                        splice_coverage=round(lab_splice / n_splice, 4) if n_splice else np.nan))

    df["y_assay"] = y
    df["y_assay_source"] = np.where(y.notna(), "official", "unlabelled")
    df = df.drop(columns=["_join"])
    return df, pd.DataFrame(cov)


def main() -> None:
    p = C.OUTPUT_DIR / f"frozen_matrix_{C.FROZEN_VERSION}.parquet"
    df = pd.read_parquet(p)
    df, cov = ingest(df)
    C.REPORT_DIR.mkdir(parents=True, exist_ok=True)
    cov.to_csv(C.REPORT_DIR / "methodB_coverage.csv", index=False)
    from .phase1_build_frozen_matrix import freeze   # lazy import avoids circularity
    manifest = freeze(df)                              # re-writes parquet AND manifest together
    print("[methodB] official labels joined; y_assay replaced")
    print(cov.to_string(index=False))
    print(f"\n[methodB] re-froze matrix; new sha256 = {manifest['sha256']}")


if __name__ == "__main__":
    main()
