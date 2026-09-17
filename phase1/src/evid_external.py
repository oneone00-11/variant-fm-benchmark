"""E7 -- external genes: do the fixed and the fitted thresholds transfer?

Two genes, and they answer different questions.

TP53 is the gene the published study already held out. Its scores exist, so it is
carried forward unchanged and the Walker cut points and the E3 thresholds are
applied to it as they stand. Its window is the published one (|offset| <= 8), so it
populates pm12 and s3_10 and says nothing about 11-50 bp.

DDX3X is new, and is outside the gene class everything else here sits in: it is not
a tumour suppressor read out by cell survival, and it is X-linked. Its saturation
editing deposit (MaveDB urn:mavedb:00000658, seventeen exon-level score sets, CC0)
carries 1,857 intron-side SNVs at 1 <= |offset| <= 50 -- more splice-region variants
than the published study's whole analysis set -- with 192 / 768 / 897 across the
three strata. That is what makes it worth the work: it is the only candidate found
that populates 11-50 bp, the territory the companion atlas locates the real failure
in, on a gene set that shares nothing with the seven.

Labels are the honest weak point and are reported twice rather than once. The
deposit publishes no functional classification, so a binary has to be constructed,
exactly as the published study had to for TP53. Two constructions are carried side
by side and both are reported:

  control_anchored : damaging if the score falls below the midpoint of the
                     synonymous and nonsense control medians
  fdr              : damaging if the trend FDR clears 0.01 and the score sits
                     below the synonymous median

Neither is the assay authors' own call. Any statement about DDX3X has to survive
both or be reported as depending on the label definition.

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

from . import config as C
from . import evid_common as K

ATLAS_REPO = Path(os.environ.get(
    "evidence-strength_ATLAS_REPO", "/Users/cliffzhang/work/functional-standard-atlas"))
EXT_DIR = Path("data/evidence/external")
REPORT_DIR = Path("reports/evidence")
CONFIG_PATH = Path("config/walker2023.yaml")

DDX3X_EXPERIMENT = "urn:mavedb:00000658"
DDX3X_SCORESETS = [f"{DDX3X_EXPERIMENT}-{c}-1" for c in "abcdefghijklmnopq"]
DDX3X_TRANSCRIPT = "ENST00000644876.2"
MAVEDB_API = "https://api.mavedb.org/api/v1/score-sets"

_SPLICE = re.compile(r"(ENST[\d.]+):(c\.[-*]?\d+[+-](\d+)[ACGT]>[ACGT])$")
_PRO_SYN = re.compile(r"p\.[A-Za-z]{3}\d+=")
_PRO_NON = re.compile(r"p\.[A-Za-z]{3}\d+(Ter|\*)")

FDR_CUT = 0.01


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
                controls.append({"class": "synonymous", "score": float(score)})
            elif _PRO_NON.search(pro):
                controls.append({"class": "nonsense", "score": float(score)})
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

    midpoint = (med_syn + med_non) / 2.0
    df["y_control_anchored"] = (df["assay_score"] < midpoint).astype(float)
    df["y_fdr"] = ((df["fdr"] < FDR_CUT) & (df["assay_score"] < med_syn)).astype(float)

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
        "controls": {"n_synonymous": int((ctl["class"] == "synonymous").sum()),
                     "n_nonsense": int((ctl["class"] == "nonsense").sum()),
                     "median_synonymous": med_syn, "median_nonsense": med_non,
                     "midpoint": midpoint, "control_auroc": auroc,
                     "gate": "OK" if auroc >= 0.80 else
                             ("WEAK_SEPARATION" if auroc >= 0.65 else "FLIP_NEEDED")},
        "labels": {
            "control_anchored": "assay_score < midpoint of control medians",
            "fdr": f"trend BH-FDR < {FDR_CUT} and assay_score < synonymous median",
            "note": "the deposit publishes no functional classification; both "
                    "definitions are constructed here and neither is the assay "
                    "authors' own call",
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
          f"syn median {med_syn:.4f}, nonsense median {med_non:.4f}")
    print(f"     damaging: control-anchored {manifest['labels']['n_damaging_control_anchored']}, "
          f"FDR {manifest['labels']['n_damaging_fdr']}")
    print(f"     wrote {out}")
    print("\n[E7] next: score the panel for these variants, then --apply DDX3X")


# ---------------------------------------------------------------------------
# apply the thresholds
# ---------------------------------------------------------------------------
# Per-scorer outputs for DDX3X, written by the atlas's own model scripts, and the
# panel column each carries. Anything absent is simply not in the external table;
# the applied-threshold report then covers fewer tools and says so.
DDX3X_SCORE_FILES = {
    "ddx3x_spliceai_walker.parquet": ["spliceai_walker"],
    "ddx3x_spliceai_d50.parquet": ["spliceai"],
    "ddx3x_cadd.parquet": ["cadd"],
    "ddx3x_phyloP100way.parquet": ["phylop100way"],
    "ddx3x_phastCons100way.parquet": ["phastcons100way"],
    "ddx3x_gpn_msa.parquet": ["gpn_msa"],
    "ddx3x_pangolin.parquet": ["pangolin_score"],
    "ddx3x_nt.parquet": ["nucleotide_transformer"],
}
DDX3X_RENAME = {"phylop100way": "phylop", "phastcons100way": "phastcons",
                "pangolin_score": "pangolin", "nucleotide_transformer": "nt"}


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
    base = base.rename(columns=DDX3X_RENAME)
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
        df = pd.read_parquet("data/external/tp53_splice_scored_v2.parquet")
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
        walker = EXT_DIR / "tp53_spliceai_walker.parquet"
        if walker.exists():
            w = pd.read_parquet(walker)
            cols = [c for c in w.columns if c.startswith("spliceai_walker")]
            df = df.merge(w[["variant_id"] + cols], on="variant_id", how="left")
        return df, ["y_control_anchored", "y_median_split", "y_mid_band_excluded"]
    if gene.upper() == "DDX3X":
        df = pd.read_parquet(EXT_DIR / "ddx3x_splice.parquet")
        scored = EXT_DIR / "ddx3x_scored.parquet"
        if scored.exists():
            extra = pd.read_parquet(scored)
            df = df.merge(extra, on="variant_id", how="left", suffixes=("", "_dup"))
        return df, ["y_control_anchored", "y_fdr"]
    raise SystemExit(f"[E7] no external gene '{gene}'")


def apply_thresholds(gene: str) -> None:
    cfg = yaml.safe_load(CONFIG_PATH.read_text())
    pp3 = cfg["thresholds"]["pp3"]["value"]
    bp4 = cfg["thresholds"]["bp4"]["value"]
    path_bands, ben_bands = K.acmg_bands(cfg)

    df, label_cols = _load_external(gene)
    ev_path = REPORT_DIR / "evidence_thresholds.csv"
    ev = pd.read_csv(ev_path) if ev_path.exists() else None

    rows = []
    tools = [t for t in K.PANEL + K.OPTIONAL if t in df.columns]
    for label in label_cols:
        for stratum in ["pm12", "s3_10", "s11_50", "all_1_50"]:
            sub = df if stratum == "all_1_50" else df[df["stratum"] == stratum]
            for tool in tools:
                s = sub.dropna(subset=[label, tool])
                y = s[label].to_numpy(dtype=float)
                v = s[tool].to_numpy(dtype=float)
                base = {"gene": gene.upper(), "label_definition": label,
                        "stratum": stratum, "tool": tool, "n": int(len(s)),
                        "n_pos": int((y == 1).sum()), "n_neg": int((y == 0).sum())}
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
                    row["walker_tier_bp4"] = K.tier_of(row["walker_lr_bp4"],
                                                       ben_bands, "benign")
                # --- E3's LOGO thresholds, applied unchanged
                if ev is not None:
                    e = ev[(ev.tool == tool) & (ev.stratum == stratum)
                           & (ev.status == "ok")]
                    for tier in ["supporting", "moderate", "strong"]:
                        t = e[e.tier == tier]
                        if not len(t):
                            continue
                        tau = float(t["pp3_logo_threshold_median"].iloc[0])
                        row[f"e3_{tier}_threshold"] = tau
                        row[f"e3_{tier}_lr_here"] = K.band_lr(y, v, tau, "upper")
                        row[f"e3_{tier}_tier_here"] = K.tier_of(
                            row[f"e3_{tier}_lr_here"], path_bands, "pathogenic")
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
