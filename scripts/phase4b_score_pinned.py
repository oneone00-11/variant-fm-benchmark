"""Score variant sets with the companion atlas's pinned scorers (P1 of the rescore plan).

This script contains NO scoring logic. It imports the atlas's own scorer modules by
path and calls their functions, so every score it emits is definition-identical to the
atlas column of the same name:

  pangolin  models/pangolin/score_fullprec.py  (git 5cf94b8, distance 50, full precision)
  spliceai  models/spliceai/score_fullprec.py  (SpliceAI 1.3.1, -A grch38 -D 50, full precision)
  nt        models/nt/score.py                 (NT-v2-500M-multi-species, 6,000-bp masked 6-mer LLR)

Each model must run under ITS OWN environment (the atlas's model-local venvs for Pangolin
and SpliceAI; a torch+transformers env for NT):

  <pangolin venv>/bin/python scripts/phase4b_score_pinned.py --model pangolin --set tp53
  <spliceai venv>/bin/python scripts/phase4b_score_pinned.py --model spliceai --set tp53
  <nt env>/bin/python       scripts/phase4b_score_pinned.py --model nt --set tp53 [--device cpu|cuda|mps]

Sets: tp53 (the 192 intron-side TP53 splice SNVs of the external validation), sanity
(200 splice variants of frozen-matrix-v1, seed 20260913, for concordance against the
existing columns), frozen_only16 (the sixteen RAD51C SNVs absent from the atlas).

Outputs <out-dir>/<model>_<set>.tsv and <out-dir>/<model>_<set>.provenance.json.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
FROZEN = REPO / "phase1" / "data" / "frozen" / "frozen_matrix_v1.parquet"
TP53 = REPO / "phase1" / "data" / "external" / "tp53_splice_scored.parquet"
COLS = {"pangolin": "pangolin_fullprec", "spliceai": "spliceai_ds_fullprec", "nt": "nucleotide_transformer"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def load_by_path(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def variant_set(which: str, atlas: Path) -> pd.DataFrame:
    if which == "tp53":
        df = pd.read_parquet(TP53, columns=["variant_id", "gene", "chrom", "pos", "ref", "alt"])
        assert len(df) == 192, len(df)
        return df.reset_index(drop=True)
    fz = pd.read_parquet(FROZEN, columns=["variant_id", "gene", "chrom", "pos", "ref", "alt", "is_splice"])
    if which == "sanity":
        return fz[fz["is_splice"]].sample(n=200, random_state=20260913).drop(columns="is_splice").reset_index(drop=True)
    if which == "frozen_only16":
        at = pd.read_parquet(atlas / "results" / "score_matrix_atlas_v2.parquet",
                             columns=["chrom", "pos", "ref", "alt"])
        ids = set(at["chrom"].astype(str) + "-" + at["pos"].astype(int).astype(str) + "-" + at["ref"] + "-" + at["alt"])
        out = fz[~fz["variant_id"].isin(ids)].drop(columns="is_splice").reset_index(drop=True)
        assert len(out) == 16, len(out)
        return out
    raise SystemExit(f"unknown set {which}")


def git_commit(repo: Path) -> str:
    try:
        return subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    except Exception:
        return "unknown"


def pkg_version(name: str) -> str:
    try:
        from importlib.metadata import version
        return version(name)
    except Exception:
        return "unknown"


def score_pangolin(df: pd.DataFrame, atlas: Path, prov: dict):
    mod = load_by_path("atlas_pangolin_fullprec", atlas / "models" / "pangolin" / "score_fullprec.py")
    import pangolin.pangolin as PP
    import torch
    torch.set_num_threads(int(os.environ.get("OMP_NUM_THREADS", "3")))
    args = mod.Args(mod.DISTANCE)
    gtf, models = mod.setup(PP)
    fullprec = mod.build_fullprec(PP)
    out = []
    for r in df.itertuples():
        sc = fullprec(0, str(r.chrom), int(r.pos), r.ref, r.alt, gtf, models, args)
        out.append(mod.aggregate(sc) if sc != -1 else np.nan)
    direct = Path(PP.__file__).resolve().parents[1]
    try:
        import glob
        du = json.load(open(glob.glob(str(Path(PP.__file__).resolve().parents[1] / "pangolin-*.dist-info" / "direct_url.json"))[0]))
        commit = du["vcs_info"]["commit_id"]
    except Exception:
        commit = "see models/pangolin/requirements.txt"
    prov.update({
        "tool": "Pangolin", "tool_version": pkg_version("pangolin"), "git_commit": commit,
        "torch": pkg_version("torch"), "gffutils": pkg_version("gffutils"),
        "parameters": {"distance": mod.DISTANCE, "mask": "True (default)", "score_cutoff": None,
                       "score_exons": "False", "score_definition": "max(splice gain, |splice loss|) over transcript records"},
        "precision": "full float (round(...,2) formatting removed by score_fullprec.build_fullprec)",
        "reference_fasta": str(mod.REF_FASTA), "reference_fasta_sha256": sha256(mod.REF_FASTA),
        "annotation_db": str(mod.GTF_DB), "annotation_db_sha256": sha256(mod.GTF_DB),
        "scorer_module": str(atlas / "models/pangolin/score_fullprec.py"),
        "scorer_module_sha256": sha256(atlas / "models/pangolin/score_fullprec.py"),
    })
    return out


def score_spliceai(df: pd.DataFrame, atlas: Path, prov: dict):
    mod = load_by_path("atlas_spliceai_fullprec", atlas / "models" / "spliceai" / "score_fullprec.py")
    from spliceai.utils import Annotator
    ann = Annotator(str(mod.REF_FASTA), mod.ANNOTATION)
    fullprec = mod.build_fullprec()
    out = []
    for r in df.itertuples():
        rec = mod.Rec(str(r.chrom), int(r.pos), r.ref, r.alt)
        try:
            out.append(mod.aggregate(fullprec(rec, ann, mod.DISTANCE, mod.MASK)))
        except Exception:
            out.append(np.nan)
    prov.update({
        "tool": "SpliceAI", "tool_version": pkg_version("spliceai"), "tensorflow": pkg_version("tensorflow"),
        "parameters": {"annotation": mod.ANNOTATION, "distance": mod.DISTANCE, "mask": mod.MASK,
                       "score_definition": "max(DS_AG, DS_AL, DS_DG, DS_DL)"},
        "precision": "full float ({:.2f} formatting replaced by {:.17g} in score_fullprec.build_fullprec)",
        "reference_fasta": str(mod.REF_FASTA), "reference_fasta_sha256": sha256(mod.REF_FASTA),
        "scorer_module": str(atlas / "models/spliceai/score_fullprec.py"),
        "scorer_module_sha256": sha256(atlas / "models/spliceai/score_fullprec.py"),
    })
    return out


def score_nt(df: pd.DataFrame, atlas: Path, prov: dict, device: str):
    mod = load_by_path("atlas_nt", atlas / "models" / "nt" / "score.py")
    import torch
    from pyfaidx import Fasta
    from transformers import AutoModelForMaskedLM, AutoTokenizer
    torch.set_num_threads(int(os.environ.get("OMP_NUM_THREADS", "4")))
    tok = AutoTokenizer.from_pretrained(mod.CHECKPOINT, trust_remote_code=True)
    model = AutoModelForMaskedLM.from_pretrained(mod.CHECKPOINT, trust_remote_code=True)
    if device == "cuda" and not torch.cuda.is_available():
        raise SystemExit("cuda requested but unavailable")
    if device == "mps" and not torch.backends.mps.is_available():
        raise SystemExit("mps requested but unavailable")
    model.to(device).eval()
    fasta = Fasta(str(mod.REF_FASTA))
    seq_cache: dict = {}
    out = []
    for r in df.itertuples():
        chrom = str(r.chrom)
        if chrom not in seq_cache:
            seq_cache[chrom] = str(fasta[chrom])
        try:
            v = mod.masked_llr(tok, model, device, seq_cache[chrom], int(r.pos) - 1, r.ref, r.alt)
        except Exception:
            v = None
        out.append(np.nan if v is None else v)
    snap = Path(model.config._name_or_path) if False else None
    try:
        from huggingface_hub import snapshot_download
        snap = snapshot_download(mod.CHECKPOINT, local_files_only=True)
        revision = Path(snap).name
    except Exception:
        revision = "unknown"
    prov.update({
        "tool": "Nucleotide Transformer v2 500M multi-species", "checkpoint": mod.CHECKPOINT,
        "checkpoint_revision_sha": revision, "torch": pkg_version("torch"), "transformers": pkg_version("transformers"),
        "device": device,
        "parameters": {"window_bp": mod.WINDOW_BP, "kmer": mod.KMER,
                       "score_definition": "logP(ref 6-mer | masked) - logP(alt 6-mer | masked); higher = more pathogenic"},
        "reference_fasta": str(mod.REF_FASTA), "reference_fasta_sha256": sha256(mod.REF_FASTA),
        "scorer_module": str(atlas / "models/nt/score.py"), "scorer_module_sha256": sha256(atlas / "models/nt/score.py"),
        "trust_remote_code": True, "trust_remote_code_note": "the checkpoint repository now ships custom model code; "
                "transformers 4.46.3 refuses to load it without this flag (the atlas run predates that change)",
        "note": "the atlas column was produced on a CUDA A6000 (models/nt/pipfreeze_runpod.txt); this run's "
                "device is recorded above and the concordance against the atlas column is reported in the run record",
    })
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--model", required=True, choices=list(COLS))
    ap.add_argument("--set", required=True, choices=["tp53", "sanity", "frozen_only16"])
    ap.add_argument("--atlas-repo", default=str(REPO.parent / "functional-standard-atlas"))
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--device", default="cpu")
    a = ap.parse_args()
    atlas = Path(a.atlas_repo).resolve()
    out_dir = Path(a.out_dir) if a.out_dir else (REPO / "data" / ("tp53" if a.set == "tp53" else "rescore"))
    out_dir.mkdir(parents=True, exist_ok=True)

    df = variant_set(a.set, atlas)
    src = TP53 if a.set == "tp53" else FROZEN
    prov = {
        "model": a.model, "set": a.set, "n_variants": int(len(df)),
        "input_file": str(src), "input_sha256": sha256(src),
        "atlas_repo": str(atlas), "atlas_commit": git_commit(atlas), "this_repo_commit": git_commit(REPO),
        "python": platform.python_version(), "platform": platform.platform(),
        "started_utc": datetime.now(timezone.utc).isoformat(),
    }
    t0 = time.time()
    if a.model == "pangolin":
        vals = score_pangolin(df, atlas, prov)
    elif a.model == "spliceai":
        vals = score_spliceai(df, atlas, prov)
    else:
        vals = score_nt(df, atlas, prov, a.device)
    df[COLS[a.model]] = vals
    prov.update({"runtime_seconds": round(time.time() - t0, 1), "scored": int(df[COLS[a.model]].notna().sum()),
                 "unscored": int(df[COLS[a.model]].isna().sum()),
                 "finished_utc": datetime.now(timezone.utc).isoformat()})
    stem = out_dir / f"{a.model}_{a.set}"
    df.to_csv(stem.with_suffix(".tsv"), sep="\t", index=False)
    prov["output_sha256"] = sha256(stem.with_suffix(".tsv"))
    (out_dir / f"{a.model}_{a.set}.provenance.json").write_text(json.dumps(prov, indent=2) + "\n")
    print(f"{a.model}/{a.set}: scored {prov['scored']}/{len(df)} in {prov['runtime_seconds']} s -> {stem}.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
