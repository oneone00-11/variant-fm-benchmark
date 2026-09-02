"""Why does the same score reach a higher evidence tier against ClinVar labels?

Section 3.5 reports that on the functional standard no object reaches Strong,
while against ClinVar most do. Three explanations are possible and they are not
the same claim:

  (a) circularity     the predictor was trained or thresholded on clinical labels
  (b) label divergence  the two label sets are not judging the same thing
  (c) selection       variants that carry a ClinVar record are easier to call

The published contrast could not separate them, because it was not computed on
the same variants: 1,675 carry an assay label, 946 a ClinVar label, 907 carry
both. Everything here is computed on those 907, so (c) is controlled by
construction at the level of set membership, and the residual effect-size
imbalance inside the intersection is examined separately.

Usage (PYTHONPATH=phase1):  python -m src.phase7_label_contrast
"""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd

from . import config as C
from .phase2_model import (usable_features, rank_within_gene,
                           rank_target_within_gene, logo_oof)
from .phase5_likelihood_ratios import (lr_at_specificity, acmg_tier,
                                       TARGET_SPEC, N_BOOT, MIN_POS, MIN_NEG)

RNG = np.random.default_rng(C.RANDOM_SEED)

# A's ten predictors, grouped by whether their training signal contains clinical
# assertions. The grouping is the companion atlas paper's Supplemental Table S11.
# Group 1 is empty here: the five ClinVar/HGMD-trained meta-predictors (REVEL,
# BayesDel, ClinPred, MetaRNN, VEST4) are in that paper's nineteen and none of
# them is among these ten. That is a finding, not an oversight -- it means
# explanation (a) cannot account for an inflation seen in this panel.
CLINICAL_LABEL_TRAINED: set[str] = set()
NOT_CLINICAL_LABEL_TRAINED = {
    "spliceai", "pangolin", "alphagenome",     # splice supervision, no variant labels
    "gpn_msa", "nt",                            # self-supervised sequence
    "phylop", "phastcons",                      # alignment statistics, not trained
    "cadd",                                     # simulated vs fixed derived
    "alphamissense",                            # population-frequency proxy
    "gnomad_af",                                # an observed frequency, not a predictor
}


def build_scores(df):
    feats = usable_features(df)
    genes = df["gene"].astype(str)
    Xrn = rank_within_gene(df, feats)
    ytr = (rank_target_within_gene(df)
           if getattr(C, "TRAIN_TARGET", "rank_within_gene") == "rank_within_gene"
           else df["func_pathogenicity"])
    scores = {"fusion_M1": logo_oof(Xrn, ytr, genes, feats, "enet"),
              "mean_M0b": logo_oof(Xrn, ytr, genes, feats, "mean")}
    for f in feats:
        scores[f"single:{f}"] = Xrn[f].to_numpy()
    for k in scores:
        s = scores[k]
        lo, hi = np.nanmin(s), np.nanmax(s)
        scores[k] = (s - lo) / (hi - lo) if hi > lo else s
    return scores


def agreement(sub):
    """AE-2: do the two label sets judge the same thing on the same variants?"""
    a = sub["y_assay"].astype(int).to_numpy()
    c = sub["y_clinvar"].astype(int).to_numpy()
    tab = pd.crosstab(pd.Series(c, name="clinvar"), pd.Series(a, name="assay"))
    n = len(a)
    obs = float((a == c).mean())
    # chance agreement under independent marginals, for Cohen's kappa
    pe = float((a == 1).mean() * (c == 1).mean() + (a == 0).mean() * (c == 0).mean())
    kappa = (obs - pe) / (1 - pe) if pe < 1 else np.nan
    return tab, {"n": n, "agreement": obs, "kappa": kappa,
                 "positive_rate_assay": float(a.mean()),
                 "positive_rate_clinvar": float(c.mean())}


def contrast(sub, scores_sub, weights=None, tag="intersection"):
    """AE-3/AE-6: LR+ and sensitivity under both label sets, same variants."""
    genes = sub["gene"].astype(str).to_numpy()
    rows = []
    for key, s in scores_sub.items():
        name = key.replace("single:", "")
        row = {"tag": tag, "object": name,
               "training_group": ("clinical labels" if name in CLINICAL_LABEL_TRAINED
                                  else "no clinical labels" if name in NOT_CLINICAL_LABEL_TRAINED
                                  else "fusion of the above")}
        for lab in ("y_assay", "y_clinvar"):
            y = sub[lab].astype(float).to_numpy()
            m = ~np.isnan(y) & ~np.isnan(s)
            lrp, _, thr, tpr, spec = lr_at_specificity(y[m].astype(int), s[m])
            short = "fun" if lab == "y_assay" else "cv"
            row |= {f"n_{short}": int(m.sum()),
                    f"lr_{short}": lrp, f"sens_{short}": tpr, f"spec_{short}": spec,
                    f"tier_{short}": acmg_tier(lrp)}
        row["delta_lr"] = row["lr_cv"] - row["lr_fun"]
        row["delta_sens"] = row["sens_cv"] - row["sens_fun"]
        rows.append(row)
    return pd.DataFrame(rows)


def boot_delta_by_group(sub, scores_sub, groups):
    """AE-4: is the inflation different between training groups? Gene-clustered."""
    genes = sub["gene"].astype(str).to_numpy()
    ug = pd.unique(genes)
    ya = sub["y_assay"].astype(float).to_numpy()
    yc = sub["y_clinvar"].astype(float).to_numpy()

    def deltas(idx):
        out = {}
        for key, s in scores_sub.items():
            name = key.replace("single:", "")
            ma = idx[~np.isnan(ya[idx]) & ~np.isnan(s[idx])]
            mc = idx[~np.isnan(yc[idx]) & ~np.isnan(s[idx])]
            if len(ma) < MIN_POS + MIN_NEG or len(mc) < MIN_POS + MIN_NEG:
                continue
            a = lr_at_specificity(ya[ma].astype(int), s[ma])[0]
            c = lr_at_specificity(yc[mc].astype(int), s[mc])[0]
            if np.isfinite(a) and np.isfinite(c):
                out[name] = c - a
        return out

    obs = deltas(np.arange(len(genes)))
    boots = []
    for _ in range(N_BOOT):
        idx = np.concatenate([np.where(genes == g)[0]
                              for g in RNG.choice(ug, len(ug), replace=True)])
        boots.append(deltas(idx))

    res = []
    for gname, members in groups.items():
        have = [m for m in members if m in obs]
        if not have:
            res.append({"group": gname, "n_objects": 0, "median_delta_lr": np.nan,
                        "lo": np.nan, "hi": np.nan, "min": np.nan, "max": np.nan})
            continue
        vals = [obs[m] for m in have]
        bs = [np.median([b[m] for m in have if m in b])
              for b in boots if any(m in b for m in have)]
        lo, hi = (np.nanpercentile(bs, [2.5, 97.5]) if len(bs) > N_BOOT * 0.5
                  else (np.nan, np.nan))
        res.append({"group": gname, "n_objects": len(have),
                    "median_delta_lr": float(np.median(vals)),
                    "lo": float(lo), "hi": float(hi),
                    "min": float(np.min(vals)), "max": float(np.max(vals))})
    return pd.DataFrame(res)


def by_effect_size(sub, scores_sub, n_strata=3):
    """AE-5: does the inflation survive inside bands of functional effect size?

    If the ClinVar-recorded variants are simply more extreme, holding |z| fixed
    should remove the gap. |z| is the functional score standardised within gene,
    the same quantity the sampling-frame analysis uses.
    """
    s = sub.copy()
    g = s.groupby("gene")["func_pathogenicity"]
    s["z"] = (s["func_pathogenicity"] - g.transform("mean")) / g.transform("std")
    s["abs_z"] = s["z"].abs()
    s["band"] = pd.qcut(s["abs_z"], n_strata,
                        labels=[f"|z| tertile {i+1}" for i in range(n_strata)])
    out = []
    for band, idx in s.groupby("band", observed=True).groups.items():
        pos = s.index.get_indexer(idx)
        sub_b = s.loc[idx]
        sc = {k: v[pos] for k, v in scores_sub.items()}
        t = contrast(sub_b, sc, tag=str(band))
        t["mean_abs_z"] = float(sub_b["abs_z"].mean())
        out.append(t)
    return pd.concat(out, ignore_index=True)


def selection_test(df, scores):
    """The decisive test for (c): are ClinVar-recorded variants easier anyway?

    Both subsets are scored against the SAME functional labels, so the label set
    is held fixed and only membership varies. If the recorded subset gives a
    higher LR+ under the assay labels too, the evidence gap is a property of
    which variants carry a ClinVar record, not of what ClinVar says about them.
    """
    rec = (df["y_assay"].notna() & df["y_clinvar"].notna()).to_numpy()
    only = (df["y_assay"].notna() & df["y_clinvar"].isna()).to_numpy()
    y = df["y_assay"].astype(float).to_numpy()
    rows = []
    for label, mask in (("ClinVar-recorded", rec), ("assay-only", only)):
        for key, s in scores.items():
            name = key.replace("single:", "")
            m = mask & ~np.isnan(y) & ~np.isnan(s)
            lrp, _, _, tpr, spec = lr_at_specificity(y[m].astype(int), s[m])
            rows.append({"subset": label, "object": name, "n": int(m.sum()),
                         "n_pos": int((y[m] == 1).sum()), "lr_plus": lrp,
                         "sensitivity": tpr, "specificity": spec,
                         "tier": acmg_tier(lrp)})
    return pd.DataFrame(rows)


def run():
    p = C.OUTPUT_DIR / f"frozen_matrix_{C.FROZEN_VERSION}.parquet"
    if not p.exists():
        sys.exit(f"[phase7] frozen matrix not found at {p}")
    df = pd.read_parquet(p)
    df = df[df["is_splice"] & df["func_pathogenicity"].notna()].reset_index(drop=True)
    scores = build_scores(df)

    both = df["y_assay"].notna() & df["y_clinvar"].notna()
    sub = df[both].reset_index(drop=True)
    pos = np.where(both.to_numpy())[0]
    sc = {k: v[pos] for k, v in scores.items()}

    C.REPORT_DIR.mkdir(parents=True, exist_ok=True)
    print("=== AE-1: the base of comparison ===")
    base = pd.DataFrame([{
        "assay_labelled": int(df["y_assay"].notna().sum()),
        "clinvar_labelled": int(df["y_clinvar"].notna().sum()),
        "both": int(both.sum()),
        "frac_of_assay": float(both.sum() / df["y_assay"].notna().sum()),
        "frac_of_clinvar": float(both.sum() / df["y_clinvar"].notna().sum()),
        "assay_only": int((df["y_assay"].notna() & df["y_clinvar"].isna()).sum()),
        "clinvar_only": int((df["y_clinvar"].notna() & df["y_assay"].isna()).sum())}])
    base.to_csv(C.REPORT_DIR / "phase7_base.csv", index=False)
    print(base.round(4).to_string(index=False))

    print("\n=== AE-2: do the two label sets agree, on the same variants? ===")
    tab, stats = agreement(sub)
    print(tab.to_string())
    print(pd.DataFrame([stats]).round(4).to_string(index=False))
    tab.to_csv(C.REPORT_DIR / "phase7_label_confusion.csv")
    pd.DataFrame([stats]).to_csv(C.REPORT_DIR / "phase7_label_agreement.csv", index=False)

    print("\n=== AE-3 / AE-6: LR+ and sensitivity under both labels, same 907 variants ===")
    tbl = contrast(sub, sc)
    tbl.to_csv(C.REPORT_DIR / "phase7_label_contrast.csv", index=False)
    print(tbl[["object", "training_group", "n_fun", "lr_fun", "sens_fun", "tier_fun",
               "lr_cv", "sens_cv", "tier_cv", "delta_lr", "delta_sens"]]
          .sort_values("delta_lr", ascending=False).round(3).to_string(index=False))

    print("\n=== AE-4: inflation by training provenance ===")
    groups = {"clinical labels": sorted(CLINICAL_LABEL_TRAINED),
              "no clinical labels": sorted(NOT_CLINICAL_LABEL_TRAINED),
              "fusions": ["fusion_M1", "mean_M0b"]}
    gb = boot_delta_by_group(sub, sc, groups)
    gb.to_csv(C.REPORT_DIR / "phase7_delta_by_training_group.csv", index=False)
    print(gb.round(3).to_string(index=False))

    print("\n=== AE-5: inflation within bands of functional effect size ===")
    bands = by_effect_size(sub, sc)
    bands.to_csv(C.REPORT_DIR / "phase7_by_effect_size.csv", index=False)
    print(bands.groupby("tag", observed=True)
          .agg(mean_abs_z=("mean_abs_z", "first"),
               median_lr_fun=("lr_fun", "median"),
               median_lr_cv=("lr_cv", "median"),
               median_delta=("delta_lr", "median"),
               median_delta_sens=("delta_sens", "median")).round(3).to_string())
    print("\n=== AE-5b: the same functional labels, on recorded versus assay-only variants ===")
    sel = selection_test(df, scores)
    sel.to_csv(C.REPORT_DIR / "phase7_selection_test.csv", index=False)
    piv = sel.pivot_table(index="object", columns="subset",
                          values=["lr_plus", "sensitivity", "n"])
    print(piv.round(3).to_string())

    print(f"\n[phase7] wrote phase7_base.csv, phase7_label_confusion.csv, "
          f"phase7_label_agreement.csv, phase7_label_contrast.csv, "
          f"phase7_delta_by_training_group.csv, phase7_by_effect_size.csv, "
          f"phase7_selection_test.csv to {C.REPORT_DIR}")


if __name__ == "__main__":
    run()
