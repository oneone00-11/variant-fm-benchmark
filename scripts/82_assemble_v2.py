"""M6: parse SpliceAI VCF + Pangolin CSV outputs, align to score_matrix on
(chrom,pos,ref,alt), and write score_matrix_v2.tsv. Distinguishes
"scored-low" (tool returned a near-0 score) from "NA" (tool did not score it).

Directionality: SpliceAI_DS and Pangolin score: higher = more splice-altering /
more likely pathogenic.
"""
import re
import pandas as pd

from config import DATA_PROCESSED, DATA_RAW, RESULTS_TABLES

SCORES = DATA_RAW / "scores"
KEY = ["chrom", "pos", "ref", "alt"]
SPLICEAI_VCF = SCORES / "spliceai" / "spliceai_out.vcf"
PANGOLIN_CSV = SCORES / "pangolin" / "pangolin_out.csv"

PANG_RE = re.compile(r"\|(-?\d+):(-?[0-9.]+)\|(-?\d+):(-?[0-9.]+)\|Warnings")


def parse_spliceai():
    rows = []
    with open(SPLICEAI_VCF) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 8:
                continue
            chrom, pos, _id, ref, alt, _q, _f, info = c[:8]
            m = re.search(r"SpliceAI=([^;\t]+)", info)
            if not m:
                continue  # not scored -> stays NA after merge
            best = None; comps = None
            for entry in m.group(1).split(","):
                p = entry.split("|")
                if len(p) < 6:
                    continue
                try:
                    ds = [float(p[2]), float(p[3]), float(p[4]), float(p[5])]
                except ValueError:
                    continue
                mx = max(ds)
                if best is None or mx > best:
                    best = mx; comps = ds
            if best is not None:
                rows.append({"chrom": chrom, "pos": int(pos), "ref": ref, "alt": alt,
                             "spliceai_ds": best, "spliceai_ds_ag": comps[0],
                             "spliceai_ds_al": comps[1], "spliceai_ds_dg": comps[2],
                             "spliceai_ds_dl": comps[3]})
    return pd.DataFrame(rows)


def parse_pangolin():
    """Manual parse: the Pangolin field can contain commas (Warnings text), so
    split on only the first 4 commas."""
    rows = []
    with open(PANGOLIN_CSV, encoding="utf-8") as f:
        header = f.readline()
        for line in f:
            parts = line.rstrip("\n").split(",", 4)
            if len(parts) < 5:
                continue
            chrom, pos, ref, alt, val = parts
            best = gain_b = loss_b = None
            for g_pos, gain, l_pos, loss in PANG_RE.findall(val):
                gn, ls = float(gain), float(loss)
                agg = max(gn, -ls)
                if best is None or agg > best:
                    best, gain_b, loss_b = agg, gn, ls
            if best is not None:
                rows.append({"chrom": str(chrom), "pos": int(pos), "ref": ref,
                             "alt": alt, "pangolin_score": best,
                             "pangolin_gain": gain_b, "pangolin_loss": loss_b})
    return pd.DataFrame(rows)


def main():
    m = pd.read_csv(DATA_PROCESSED / "score_matrix.tsv", sep="\t",
                    dtype={"chrom": str}, low_memory=False)
    m["pos"] = m["pos"].astype(int)
    # drop the v1 placeholder columns we are now filling
    for col in ["spliceai_ds", "pangolin_score", "spliceai_ds_covered",
                "pangolin_score_covered", "alphagenome_splice", "alphagenome_splice_covered"]:
        if col in m.columns:
            m = m.drop(columns=col)

    n0 = len(m)
    if SPLICEAI_VCF.exists():
        sp = parse_spliceai()
        m = m.merge(sp, on=KEY, how="left")
        print(f"SpliceAI: parsed {len(sp)} scored; merged (rows {len(m)})")
    else:
        print("SpliceAI output missing!"); m["spliceai_ds"] = pd.NA
    if PANGOLIN_CSV.exists():
        pg = parse_pangolin()
        m = m.merge(pg, on=KEY, how="left")
        print(f"Pangolin: parsed {len(pg)} scored; merged (rows {len(m)})")
    else:
        print("Pangolin output missing!"); m["pangolin_score"] = pd.NA
    assert len(m) == n0, "merge changed row count"

    m["alphagenome_splice"] = pd.NA  # still batch-B / needs key
    for c in ["spliceai_ds", "pangolin_score", "alphagenome_splice"]:
        m[c + "_covered"] = m[c].notna()

    m.to_csv(DATA_PROCESSED / "score_matrix_v2.tsv", sep="\t", index=False)
    print(f"Wrote score_matrix_v2.tsv: {len(m)} rows")

    # coverage: scored vs scored-low vs NA, per model x region
    LOW = 0.01
    rows = []
    for model in ["spliceai_ds", "pangolin_score"]:
        for reg, g in m.groupby("region_class"):
            scored = g[model].notna()
            n_sc = int(scored.sum())
            n_low = int((scored & (g[model].fillna(9) < LOW)).sum())
            rows.append({"model": model, "region_class": reg, "n_total": len(g),
                         "n_scored": n_sc, "n_scored_low(<0.01)": n_low,
                         "n_NA": len(g) - n_sc})
    cov = pd.DataFrame(rows)
    cov.to_csv(RESULTS_TABLES / "splice_model_coverage.tsv", sep="\t", index=False)
    print("\n=== splice_model_coverage (scored / scored-low / NA) ===")
    print(cov.to_string(index=False))

    # spot-check known pathogenic splice variants
    print("\n=== SPOT-CHECK known pathogenic splice ===")
    checks = m[m.hgvs_nt.isin([
        "NM_007294.3:c.5194-1G>C", "NM_007294.3:c.4986+1G>A",
        "NM_007294.3:c.5074+1G>A", "ENST00000256474.3:c.340+1G>T"]) |
        ((m.gene == "VHL") & (m.splice_class == "splice_core") & (m.clnsig_class == "P/LP"))]
    cols = ["gene", "hgvs_nt", "region_class", "splice_class", "functional_score",
            "spliceai_ds", "pangolin_score", "clnsig_class"]
    print(checks[cols].head(12).to_string(index=False))


if __name__ == "__main__":
    main()
