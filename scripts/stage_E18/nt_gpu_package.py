"""E18.2, NT: the package for scoring Nucleotide Transformer on the original GPU setup.

The analysis set's `nt` column is the companion atlas's, scored on a rented NVIDIA RTX
A6000 (RunPod, 2026-08-01) with the environment recorded in the atlas's
models/nt/pipfreeze_runpod.txt (torch 2.11.0+cu128, transformers 4.46.3). On this Mac
the same scorer reproduces 4 of the 39 anchors bit for bit on the CPU and 2 on the GPU
(MPS), the largest difference 1.1e-05 (E18 report, section 4). The author's decision of
2026-10-23: score NT on an RTX A6000 in the recorded environment, the 39 anchors and the
6,667 external variants in one run; it counts only if all 39 anchors equal the analysis
set's values bit for bit. Otherwise NT stays out of the external test, with no
tolerance and no re-scoring of the seven genes.

This builds the package the author runs on the pod. It holds no key (NT needs none):
  variants.tsv          the 39 anchors (with the analysis set's stored value), then the
                        6,667 external variants, in that order
  regions.fa            reference slices around them: Ensembl release 112, from the
                        atlas's grch38_subset.fa for chromosomes 2, 3, 13, 16 and 17 and
                        from phase1/data/evidence/refs/ for 1, 7 and X; each slice extends
                        at least SLACK nt beyond its variants, and every variant's model
                        window, cut from its slice, is checked here to equal the window cut
                        from the whole chromosome
  atlas_nt_score.py     the atlas's models/nt/score.py at release v2.5.0, verbatim
  pipfreeze_runpod.txt  the atlas's record of the pod's environment, verbatim
  requirements_pod.txt  the same, less the one line conda's own packaging left in it
  run_nt_gpu.py         the run: it refuses to start unless the GPU, the CUDA build, the
                        packages and the model revision are the recorded ones
  README.md             the steps
  MANIFEST.sha256       every file's sha256

What the original run left unrecorded is taken from score.py's defaults and named in
the package: float32, one variant per forward pass, the 6,000-nt window, and the model
repository's default branch (it has pointed at 06615c16 since 2025-10-22, before the
run). The model loads as score.py loads it; it needs the repository's own code, which
transformers runs only after an interactive "y" or trust_remote_code=True, so the run
passes the latter, and the tokenizer loads with no argument, as in score.py.

Run from the repository root, in the atlas's NT environment (it has pyfaidx):
    ../envs/nt/bin/python scripts/stage_E18/nt_gpu_package.py --out ~/Desktop/E18_NT_GPU_package
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import os
import shutil
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog, rel, sha256_of  # noqa: E402

# the companion atlas checkout, by evid_common's rule (EVID_ATLAS_REPO, else the checkout
# beside this repository); not imported, because this runs in the atlas's NT environment,
# which has pyfaidx and not this repository's packages
ATLAS_REPO = Path(os.environ.get("EVID_ATLAS_REPO", str(REPO.parent / "functional-standard-atlas")))

EXT = REPO / "phase1/data/evidence/external"
ANCHORS = EXT / "e18_anchors.parquet"
ATLAS_SCORE = ATLAS_REPO / "models/nt/score.py"
ATLAS_FREEZE = ATLAS_REPO / "models/nt/pipfreeze_runpod.txt"
SCORE_SHA = "d69318b80ce001342ea493d4cef38ea007735746d5cbb1b5a6005ad99e692d01"
FREEZE_SHA = "21bff23888b101bd273f6fb4554b321ad1a5151088fafb51343f743453b1ced7"
SUBSET_FA = ATLAS_REPO / "data/refs/grch38_subset.fa"
REFS = REPO / "phase1/data/evidence/refs"
SLACK = 4000          # nt beyond the outermost variants of a slice (the window is 3,000 each way)
GAP = 20000           # variants further apart than this go to separate slices
NEW_GENES = ["CTCF", "SFPQ", "RAD51D", "XRCC2", "SBDS"]
RUN_SCRIPT = Path(__file__).resolve().parent / "nt_gpu_run.py"
README = Path(__file__).resolve().parent / "nt_gpu_README.md"


def external_variants() -> pd.DataFrame:
    """E18.1's set, in E18.1's order: DDX3X, TP53, then the five further genes at
    3 <= |offset| <= 50."""
    parts = [pd.read_parquet(EXT / "ddx3x_splice.parquet"),
             pd.read_parquet(EXT / "tp53_splice_scored_12nt.parquet")]
    for g in NEW_GENES:
        v = pd.read_parquet(EXT / f"{g.lower()}_variants.parquet")
        parts.append(v[v["offset"].abs().between(3, 50)])
    cols = ["gene", "variant_id", "chrom", "pos", "ref", "alt"]
    df = pd.concat([p[cols] for p in parts], ignore_index=True)
    df["chrom"], df["pos"] = df["chrom"].astype(str), df["pos"].astype(int)
    return df


def _fasta(chrom: str):
    from pyfaidx import Fasta
    if chrom in {"2", "3", "13", "16", "17"}:
        return Fasta(str(SUBSET_FA))[chrom]
    return Fasta(str(REFS / f"chr{chrom}.fa"))[chrom]


def _window(nt, seq: str, pos0: int) -> tuple[int, str]:
    """The window masked_llr cuts, by its own lines (start and the sequence)."""
    half = nt.WINDOW_BP // 2
    start = max(0, pos0 - half)
    end = min(len(seq), pos0 + half)
    var_off = pos0 - start
    start += var_off % nt.KMER
    return start, seq[start:end].upper()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--anchors", default=str(ANCHORS))
    args = ap.parse_args()
    out = Path(args.out).expanduser()
    anchors_path = Path(args.anchors).expanduser()
    with RunLog(__file__) as log:
        for p, want in ((ATLAS_SCORE, SCORE_SHA), (ATLAS_FREEZE, FREEZE_SHA)):
            if sha256_of(p) != want:
                raise SystemExit(f"{p} is not the atlas release's file")
            log.input(p)
        spec = importlib.util.spec_from_file_location("atlas_nt", ATLAS_SCORE)
        nt = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(nt)
        a = pd.read_parquet(anchors_path)
        log.input(anchors_path)
        anchors = a[["gene", "variant_id", "chrom", "pos", "ref", "alt"]].assign(
            set="anchor", stored_nt=[repr(float(x)) for x in a["nt"]])
        ext = external_variants().assign(set="external", stored_nt="")
        rows = pd.concat([anchors, ext], ignore_index=True)[
            ["set", "gene", "variant_id", "chrom", "pos", "ref", "alt", "stored_nt"]]
        if rows["variant_id"].duplicated().any() or len(ext) != 6667 or len(anchors) != 39:
            raise SystemExit("the variant list is not the 39 anchors and 6,667 external variants")
        if out.exists():
            raise SystemExit(f"{out} exists; the package is written once")
        out.mkdir(parents=True)
        # reference slices, and every window checked against the whole chromosome
        slices, checked = [], 0
        for chrom, g in rows.groupby("chrom", sort=True):
            whole = str(_fasta(chrom)[:])
            pos = sorted(g["pos"])
            groups, cur = [], [pos[0]]
            for p in pos[1:]:
                if p - cur[-1] > GAP:
                    groups.append(cur)
                    cur = [p]
                else:
                    cur.append(p)
            groups.append(cur)
            for grp in groups:
                lo, hi = max(0, grp[0] - 1 - SLACK), min(len(whole), grp[-1] + SLACK)
                seq = whole[lo:hi]
                slices.append((chrom, lo, seq))
                for r in g[g["pos"].between(grp[0], grp[-1])].itertuples():
                    s_full, w_full = _window(nt, whole, r.pos - 1)
                    s_cut, w_cut = _window(nt, seq, r.pos - 1 - lo)
                    if w_full != w_cut or s_full != s_cut + lo or whole[r.pos - 1].upper() != r.ref:
                        raise SystemExit(f"{r.variant_id}: the slice does not give the window")
                    checked += 1
        if checked != len(rows):
            raise SystemExit(f"{checked} windows checked, {len(rows)} variants")
        with open(out / "regions.fa", "w") as f:
            for chrom, lo, seq in slices:
                f.write(f">{chrom}:{lo}\n")
                for i in range(0, len(seq), 60):
                    f.write(seq[i:i + 60] + "\n")
        rows.to_csv(out / "variants.tsv", sep="\t", index=False)
        shutil.copyfile(ATLAS_SCORE, out / "atlas_nt_score.py")
        shutil.copyfile(ATLAS_FREEZE, out / "pipfreeze_runpod.txt")
        # what pip installs: the record, less the line conda's own packaging left in it
        (out / "requirements_pod.txt").write_text("".join(
            line + "\n" for line in ATLAS_FREEZE.read_text().splitlines()
            if line.strip() and not line.startswith("packaging @ ")))
        shutil.copyfile(RUN_SCRIPT, out / "run_nt_gpu.py")
        shutil.copyfile(README, out / "README.md")
        files = sorted(p for p in out.iterdir() if p.is_file())
        (out / "MANIFEST.sha256").write_text("".join(
            f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}\n" for p in files))
        for p in files + [out / "MANIFEST.sha256"]:
            log.output(p)
        log.note("variants", {"anchors": len(anchors), "external": len(ext)})
        log.note("slices", [{"chrom": c, "start0": lo, "length": len(s)} for c, lo, s in slices])
        log.note("windows_checked", checked)
        print(f"{rel(out)}: {len(anchors)} anchors + {len(ext)} external variants; "
              f"{len(slices)} reference slices; {checked} windows checked")


if __name__ == "__main__":
    main()
