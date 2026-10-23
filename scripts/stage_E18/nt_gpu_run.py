"""Score Nucleotide Transformer on the pod: the 39 anchors and the 6,667 external
variants in one run (E18.2; the author's decision of 2026-10-23).

It refuses to start unless:
  * every file of the package is the one MANIFEST.sha256 lists;
  * the GPU is an NVIDIA RTX A6000 (Ampere, compute capability 8.6), the card of the
    original run, and not an RTX 6000 Ada;
  * Python is 3.12, torch is 2.11.0+cu128 built for CUDA 12.8, transformers is 4.46.3,
    and `pip freeze --all` holds every package of pipfreeze_runpod.txt at its recorded
    version and nothing else but pip, wheel and packaging;
  * the model resolves to the revision its repository's default branch has named since
    2025-10-22, before the original run: 06615c1660c892fc199840c18123f8385b3542a8.

Then it scores every variant with the atlas's own masked_llr (atlas_nt_score.py, the
atlas's models/nt/score.py verbatim), as score.py does: the tokenizer and the model
from_pretrained(CHECKPOINT), the model with trust_remote_code=True (the interactive "y"
the original run had to give, since the model is the repository's own code), on "cuda",
in eval mode, float32, one variant per forward pass, in the order of variants.tsv,
anchors first. Nothing the original run did not record is tuned here.

It writes out/: results.tsv (every value as Python's repr, and for each anchor whether
it equals the analysis set's stored value bit for bit), environment.json (nvidia-smi, the
driver, CUDA and cuDNN, the GPU's name and capability, Python, pip freeze, the model's
revision, the times) and MANIFEST.sha256. No key is needed or read.

--trial first (E18.2: every scorer is first run on the anchors and 50 variants, to
measure throughput): the same checks, the 39 anchors and the first 50 external
variants, into out_trial/, and the full run's time projected from it. The trial's
values are not used; the full run scores the anchors again, with all 6,667.

    python run_nt_gpu.py --trial
    python run_nt_gpu.py
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import platform
import re
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
TRIAL_N = 50
GPU = "NVIDIA RTX A6000"
CAPABILITY = (8, 6)
TORCH, CUDA, TRANSFORMERS = "2.11.0+cu128", "12.8", "4.46.3"
REVISION = "06615c1660c892fc199840c18123f8385b3542a8"
ALLOWED_EXTRA = {"pip", "wheel", "packaging"}


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def _norm(name: str) -> str:
    return re.sub(r"[-_.]+", "-", name).lower()


def check_package() -> None:
    for line in (HERE / "MANIFEST.sha256").read_text().splitlines():
        sha, name = line.split("  ", 1)
        if _sha(HERE / name) != sha:
            raise SystemExit(f"{name} is not the package's file: stopping")


def check_environment() -> dict:
    import torch
    import transformers
    if sys.version_info[:2] != (3, 12):
        raise SystemExit(f"Python {platform.python_version()}, not 3.12: stopping")
    if not torch.cuda.is_available():
        raise SystemExit("no CUDA device: stopping")
    name, cap = torch.cuda.get_device_name(0), torch.cuda.get_device_capability(0)
    if name != GPU or tuple(cap) != CAPABILITY:
        raise SystemExit(f"the GPU is {name} (capability {cap}), not {GPU} {CAPABILITY}: stopping")
    found = (torch.__version__, torch.version.cuda, transformers.__version__)
    if found != (TORCH, CUDA, TRANSFORMERS):
        raise SystemExit(f"torch, CUDA, transformers are {found}, not "
                         f"{(TORCH, CUDA, TRANSFORMERS)}: stopping")
    recorded = {}
    for line in (HERE / "pipfreeze_runpod.txt").read_text().splitlines():
        if "==" in line:
            n, v = line.split("==", 1)
            recorded[_norm(n)] = v.strip()
    freeze = subprocess.run([sys.executable, "-m", "pip", "freeze", "--all"], check=True,
                            capture_output=True, text=True).stdout
    have = {}
    for line in freeze.splitlines():
        if "==" in line:
            n, v = line.split("==", 1)
            have[_norm(n)] = v.strip()
        elif " @ " in line:
            have[_norm(line.split(" @ ", 1)[0])] = line.split(" @ ", 1)[1].strip()
    wrong = {n: (v, have.get(n)) for n, v in recorded.items() if have.get(n) != v}
    extra = sorted(set(have) - set(recorded) - ALLOWED_EXTRA)
    if wrong or extra:
        raise SystemExit(f"the environment is not the recorded one: differing {wrong}; "
                         f"not recorded {extra}: stopping")
    smi = subprocess.run(["nvidia-smi"], capture_output=True, text=True).stdout
    query = subprocess.run(["nvidia-smi", "--query-gpu=name,driver_version,memory.total,"
                            "compute_cap", "--format=csv"], capture_output=True, text=True).stdout
    return {"gpu_name": name, "gpu_capability": list(cap), "torch": torch.__version__,
            "cuda": torch.version.cuda, "cudnn": torch.backends.cudnn.version(),
            "transformers": transformers.__version__, "python": sys.version,
            "platform": platform.platform(), "nvidia_smi": smi, "nvidia_smi_query": query,
            "pip_freeze_all": freeze.splitlines()}


def regions() -> list[tuple[str, int, str]]:
    out, head, seq = [], None, []
    for line in (HERE / "regions.fa").read_text().splitlines():
        if line.startswith(">"):
            if head:
                out.append((*head, "".join(seq)))
            chrom, start0 = line[1:].split(":")
            head, seq = (chrom, int(start0)), []
        else:
            seq.append(line)
    out.append((*head, "".join(seq)))
    return out


def main() -> int:
    if sys.argv[1:] not in ([], ["--trial"]):
        raise SystemExit("usage: python run_nt_gpu.py [--trial]")
    trial = sys.argv[1:] == ["--trial"]
    out_dir = HERE / ("out_trial" if trial else "out")
    check_package()
    env = check_environment()
    import pandas as pd
    import torch
    from transformers import AutoModelForMaskedLM, AutoTokenizer
    spec = importlib.util.spec_from_file_location("atlas_nt", HERE / "atlas_nt_score.py")
    nt = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(nt)
    variants = pd.read_csv(HERE / "variants.tsv", sep="\t", dtype=str, keep_default_na=False)
    n_full = len(variants)
    if trial:
        variants = pd.concat([variants[variants.set == "anchor"],
                              variants[variants.set == "external"].head(TRIAL_N)])
    slices = regions()
    start, t0 = _now(), time.time()
    tok = AutoTokenizer.from_pretrained(nt.CHECKPOINT)
    model = AutoModelForMaskedLM.from_pretrained(nt.CHECKPOINT, trust_remote_code=True)
    revision = getattr(model.config, "_commit_hash", None)
    if revision != REVISION or REVISION not in type(model).__module__:
        raise SystemExit(f"the model resolved to revision {revision} "
                         f"({type(model).__module__}), not {REVISION}: stopping")
    device = "cuda"
    model.to(device).eval()
    print(f"nt: {nt.CHECKPOINT} at {revision} on {env['gpu_name']}, {model.dtype}", flush=True)
    rows, t_loop = [], time.time()
    for i, r in enumerate(variants.itertuples(index=False), 1):
        pos0 = int(r.pos) - 1
        hit = [(s0, seq) for c, s0, seq in slices if c == r.chrom and s0 <= pos0 < s0 + len(seq)]
        if len(hit) != 1:
            raise SystemExit(f"{r.variant_id}: {len(hit)} reference slices hold it: stopping")
        s0, seq = hit[0]
        val = nt.masked_llr(tok, model, device, seq, pos0 - s0, r.ref, r.alt)
        got = "" if val is None else repr(float(val))
        same = (r.set == "anchor" and got != "" and float(got) == float(r.stored_nt))
        rows.append(r._asdict() | {"nt": got, "identical": same if r.set == "anchor" else ""})
        if i % 500 == 0:
            print(f"  {i}/{len(variants)} scored, {time.time() - t0:.0f} s", flush=True)
    loop_seconds = time.time() - t_loop
    out_dir.mkdir(exist_ok=True)
    res = pd.DataFrame(rows)
    res.to_csv(out_dir / "results.tsv", sep="\t", index=False)
    anc = res[res.set == "anchor"]
    n_same = int((anc["identical"] == True).sum())  # noqa: E712
    env |= {"checkpoint": nt.CHECKPOINT, "model_revision": revision,
            "model_class": f"{type(model).__module__}.{type(model).__name__}",
            "dtype": str(model.dtype), "device": device, "start_utc": start, "end_utc": _now(),
            "seconds": round(time.time() - t0, 1), "scoring_seconds": round(loop_seconds, 1),
            "trial": trial, "n_variants": len(res),
            "n_scored": int((res["nt"] != "").sum()), "anchors_identical": n_same,
            "anchors": len(anc), "package_manifest_sha256": _sha(HERE / "MANIFEST.sha256")}
    if trial:
        per = loop_seconds / len(res)
        env["projected_full_run_minutes"] = round(per * n_full / 60, 1)
    (out_dir / "environment.json").write_text(json.dumps(env, indent=1) + "\n")
    files = sorted(p for p in out_dir.iterdir() if p.is_file() and p.name != "MANIFEST.sha256")
    (out_dir / "MANIFEST.sha256").write_text("".join(f"{_sha(p)}  {p.name}\n" for p in files))
    print(f"anchors: {n_same}/{len(anc)} identical to the analysis set, bit for bit")
    if trial:
        print(f"throughput: {per:.3f} s per variant; the full run of {n_full} variants "
              f"projects to {env['projected_full_run_minutes']} min")
    print(f"wrote {out_dir}: results.tsv, environment.json, MANIFEST.sha256")
    return 0


if __name__ == "__main__":
    sys.exit(main())
