"""M5-B: assemble score_matrix.tsv by aligning every model's scores to
analysis_ready_v2 on (chrom,pos,ref,alt). Coverage left NA where a model does not
score that variant (no 0-fill, no extrapolation). Also writes the coverage report.

Directionality (documented; applied at eval time, not here):
  cadd_phred, alphamissense, phylop100way, phastcons100way : higher = more pathogenic
  gnomad_af_global, gnomad_af_popmax                       : higher = more BENIGN
"""
import pandas as pd

from config import DATA_PROCESSED, DATA_RAW, RESULTS_TABLES

SCORES = DATA_RAW / "scores"
KEY = ["chrom", "pos", "ref", "alt"]

DIRECTION = {
    "cadd_phred": "higher=pathogenic", "alphamissense": "higher=pathogenic",
    "phylop100way": "higher=pathogenic", "phastcons100way": "higher=pathogenic",
    "gnomad_af_global": "higher=benign", "gnomad_af_popmax": "higher=benign",
}
# which region_classes each model is BY DESIGN meant to score (for coverage interpretation)
APPLICABLE = {
    "cadd_phred": "all", "phylop100way": "all", "phastcons100way": "all",
    "gnomad_af_global": "all", "gnomad_af_popmax": "all",
    "alphamissense": "missense-only",
}


def load(name, cols):
    f = SCORES / name
    if not f.exists():
        print(f"WARNING: {name} missing — its columns will be all-NA")
        return None
    d = pd.read_csv(f, sep="\t", dtype={"chrom": str}, low_memory=False)
    d["pos"] = d["pos"].astype(int)
    return d[KEY + cols].drop_duplicates(KEY)


def main():
    ar = pd.read_csv(DATA_PROCESSED / "analysis_ready_v2.tsv", sep="\t",
                     dtype={"chrom": str}, low_memory=False)
    ar["pos"] = ar["pos"].astype(int)
    m = ar.copy()
    for name, cols in [("cadd.tsv", ["cadd_phred", "cadd_raw"]),
                       ("alphamissense.tsv", ["alphamissense", "am_class"]),
                       ("gnomad_af.tsv", ["gnomad_af_global", "gnomad_af_popmax"]),
                       ("conservation.tsv", ["phylop100way", "phastcons100way"])]:
        d = load(name, cols)
        if d is not None:
            before = len(m)
            m = m.merge(d, on=KEY, how="left")
            assert len(m) == before, "merge changed row count!"

    # models not done this batch -> explicit NA columns so schema is stable
    for blocked in ["spliceai_ds", "pangolin_score", "alphagenome_splice"]:
        m[blocked] = pd.NA

    SCORE_COLS = ["cadd_phred", "alphamissense", "phylop100way", "phastcons100way",
                  "gnomad_af_global", "gnomad_af_popmax",
                  "spliceai_ds", "pangolin_score", "alphagenome_splice"]
    for c in SCORE_COLS:
        m[c + "_covered"] = m[c].notna()

    m.to_csv(DATA_PROCESSED / "score_matrix.tsv", sep="\t", index=False)
    print(f"Wrote score_matrix.tsv: {len(m)} variants x {len(SCORE_COLS)} models")

    # --- coverage report: model x region ---
    rows = []
    for model in SCORE_COLS:
        for reg, g in m.groupby("region_class"):
            n_tot = len(g)
            n_sc = int(g[model].notna().sum())
            rows.append({"model": model, "region_class": reg,
                         "n_total": n_tot, "n_scored": n_sc, "n_NA": n_tot - n_sc,
                         "applicable": APPLICABLE.get(model, "n/a (not run this batch)")})
    cov = pd.DataFrame(rows)
    cov.to_csv(RESULTS_TABLES / "score_coverage_report.tsv", sep="\t", index=False)

    print("\n=== coverage by model (overall scored / total) ===")
    ov = m[SCORE_COLS].notna().sum()
    for c in SCORE_COLS:
        print(f"  {c:20s} {int(ov[c]):6d}/{len(m)}   dir={DIRECTION.get(c,'-')}  "
              f"applicable={APPLICABLE.get(c,'pending 2nd batch')}")
    print("\n=== coverage by model x region (pivot of n_scored) ===")
    print(cov.pivot_table(index="model", columns="region_class",
                          values="n_scored", fill_value=0).to_string())


if __name__ == "__main__":
    main()
