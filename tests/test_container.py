"""The reference image (stage B6) is defined by files that hold each other to account.

The image's lock has the evidence lock's packages at the same versions, each with the
sha256 of one wheel; container/Dockerfile pins what container/ENVIRONMENT.json recorded
from the registry and from the base image, and copies nothing of the repository but
that lock; container/run.sh runs with reproduce.sh's settings and fixes the rest; the
build workflow does nothing unless started by hand; and the entry point's --fetch-only
fetches and runs no stage.
"""
from __future__ import annotations

import importlib.util
import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
LOCK = REPO / "requirements-evidence.lock.txt"
IMAGE_LOCK = REPO / "container" / "requirements.linux-amd64.txt"
DOCKERFILE = REPO / "container" / "Dockerfile"
RUN_SH = REPO / "container" / "run.sh"
ENVIRONMENT = REPO / "container" / "ENVIRONMENT.json"
WORKFLOW = REPO / ".github" / "workflows" / "container.yml"
THREADS = re.compile(r"\b([A-Z_]+_(?:NUM_THREADS|MAXIMUM_THREADS))=(\d+)")


def _lines(path: Path) -> list[str]:
    return [line for line in path.read_text().splitlines()
            if line.strip() and not line.startswith("#")]


def _entry_point():
    spec = importlib.util.spec_from_file_location(
        "reproduce_evidence", REPO / "scripts" / "reproduce_evidence.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# --------------------------------------------------------------------------
# the lock
# --------------------------------------------------------------------------
def test_the_image_lock_has_the_evidence_locks_packages_at_the_same_versions():
    assert {line.split()[0] for line in _lines(IMAGE_LOCK)} == set(_lines(LOCK))


def test_every_package_of_the_image_lock_carries_the_sha256_of_one_wheel():
    for line in _lines(IMAGE_LOCK):
        assert re.fullmatch(r"[A-Za-z0-9_.\-]+==\S+ --hash=sha256:[0-9a-f]{64}", line), line


# --------------------------------------------------------------------------
# the Dockerfile and the record it is pinned to
# --------------------------------------------------------------------------
def test_the_base_image_is_the_recorded_python_3_12_13_by_digest():
    base = json.loads(ENVIRONMENT.read_text())["base_image"]
    assert base["config_env"]["PYTHON_VERSION"] == "3.12.13"
    assert (base["config_architecture"], base["config_os"]) == ("amd64", "linux")
    assert base["reference"] == f"python:{base['tag']}@{base['linux_amd64_manifest_digest']}"
    assert re.findall(r"^FROM\s+(\S+)", DOCKERFILE.read_text(), re.M) == [base["reference"]]


def test_every_debian_package_is_pinned_to_the_recorded_version():
    pins = json.loads(ENVIRONMENT.read_text())["apt_pins"]["packages"]
    installed = dict(re.findall(r"^\s+([a-z0-9.+-]+)=(\S+) \\$", DOCKERFILE.read_text(), re.M))
    assert installed == {pkg: p["version"] for pkg, p in pins.items()}
    assert "git" in installed                     # --verify reads HEAD with git


def test_the_image_installs_wheels_by_hash_and_copies_nothing_else():
    text = DOCKERFILE.read_text()
    assert "--require-hashes --no-deps --only-binary=:all:" in text
    assert "pip check" in text and "pip freeze" in text
    assert re.findall(r"^(?:COPY|ADD)\s+(.+)$", text, re.M) == [
        "requirements.linux-amd64.txt /opt/reference/requirements.linux-amd64.txt"]


def test_each_section_of_the_record_names_a_script_that_exists():
    for name, section in json.loads(ENVIRONMENT.read_text()).items():
        assert (REPO / section["written_by"]).is_file(), name


# --------------------------------------------------------------------------
# container/run.sh
# --------------------------------------------------------------------------
def test_run_sh_keeps_reproduce_sh_settings_and_fixes_the_rest():
    threads = set(THREADS.findall((REPO / "reproduce.sh").read_text()))
    run = RUN_SH.read_text()
    assert threads and {v for _, v in threads} == {"1"}
    assert set(THREADS.findall(run)) == threads
    for setting in ("PYTHONHASHSEED=0", "TZ=UTC", "LC_ALL=C.UTF-8", "unset MATPLOTLIBRC",
                    'git config --global --add safe.directory "$PWD"'):
        assert setting in run, setting
    assert "scripts/reproduce_evidence.py --fetch-only" in run
    assert "scripts/reproduce_evidence.py --verify" in run


def test_run_sh_holds_the_cpu_settings_in_one_place():
    """B6.4's interface: the two values, and the variables exported only from them."""
    run = RUN_SH.read_text()
    for var in ("OPENBLAS_CORETYPE", "NPY_DISABLE_CPU_FEATURES"):
        assert len(re.findall(rf"^{var}_PIN=", run, re.M)) == 1, var
        assert f'export {var}="${var}_PIN"' in run, var


# --------------------------------------------------------------------------
# the build workflow
# --------------------------------------------------------------------------
def test_the_container_workflow_does_nothing_unless_started_by_hand():
    import yaml
    wf = yaml.safe_load(WORKFLOW.read_text())
    on = wf[True]                       # YAML 1.1 reads the key "on" as true
    assert set(on) == {"push", "workflow_dispatch"}
    # the push only registers the workflow, and only a change to it starts one
    assert on["push"] == {"branches": ["evidence/rework"],
                          "paths": [".github/workflows/container.yml"]}
    for name, job in wf["jobs"].items():
        assert job["if"].startswith("github.event_name == 'workflow_dispatch' && "), name
    assert wf["permissions"] == {"contents": "read"}
    assert wf["jobs"]["build"]["permissions"] == {"contents": "read", "packages": "write"}


def test_the_image_is_named_after_the_owner_github_reports():
    text = WORKFLOW.read_text()
    assert 'gh repo view "$GITHUB_REPOSITORY" --json owner --jq .owner.login' in text
    assert "--platform linux/amd64" in text and "docker save" in text


# --------------------------------------------------------------------------
# --fetch-only
# --------------------------------------------------------------------------
def test_fetch_only_fetches_every_input_and_runs_no_stage(monkeypatch, capsys):
    E = _entry_point()
    calls = []
    for name in ("fetch_clinvar", "fetch_reference_chromosomes", "fetch_label_audit_cache"):
        monkeypatch.setattr(E, name, lambda name=name: calls.append(name))
    monkeypatch.setattr(E, "choose_atlas", lambda release_only, fetch: (
        calls.append(("atlas", release_only, fetch)), Path("/atlas"))[1])
    monkeypatch.setattr(E, "atlas_version", lambda atlas: "2.5.0")
    monkeypatch.setattr(E, "run_stage", lambda *a, **k: pytest.fail("a stage ran"))
    monkeypatch.setattr(E, "baseline", lambda *a, **k: pytest.fail("a baseline was taken"))
    monkeypatch.setattr(sys, "argv", ["reproduce_evidence.py", "--fetch-only"])
    E.main()
    assert calls == ["fetch_clinvar", "fetch_reference_chromosomes",
                     "fetch_label_audit_cache", ("atlas", True, True)]
    assert "no stage ran" in capsys.readouterr().out


@pytest.mark.parametrize("extra", [["--verify"], ["--from", "3"]])
def test_fetch_only_refuses_a_stage_option(monkeypatch, extra):
    E = _entry_point()
    monkeypatch.setattr(sys, "argv", ["reproduce_evidence.py", "--fetch-only", *extra])
    monkeypatch.setattr(E, "fetch_clinvar", lambda: pytest.fail("fetched"))
    with pytest.raises(SystemExit):
        E.main()


# --------------------------------------------------------------------------
# B6.4: the CPU settings and the probe records they rest on
# --------------------------------------------------------------------------
PROBES = REPO / "container" / "cpu_probe"
KERNEL_FLAGS = {"Nehalem": {"sse4_2", "popcnt"}}      # what the kernel type needs


def _pin(var: str) -> str:
    return re.search(rf'^{var}_PIN="([^"]*)"', RUN_SH.read_text(), re.M).group(1)


def test_the_cpu_settings_rest_on_the_probes_of_both_machines():
    core, disabled = _pin("OPENBLAS_CORETYPE"), _pin("NPY_DISABLE_CPU_FEATURES").split()
    probes = [json.loads((PROBES / f"{n}.json").read_text())
              for n in ("github-linux", "mac-rosetta")]
    assert {p["settings"]["OPENBLAS_CORETYPE"] for p in probes} == {None}   # before
    assert core in KERNEL_FLAGS
    for p in probes:
        assert KERNEL_FLAGS[core] <= set(p["cpuinfo"]["flags"]), p["name"]
        assert set(disabled) <= set(p["numpy_cpu"]["dispatch"]), p["name"]
        kept = [d for d in p["numpy_cpu"]["dispatch"] if d not in disabled]
        assert set(kept) <= set(p["numpy_cpu"]["found"]), p["name"]


@pytest.mark.parametrize("name", ["github-linux-pinned", "mac-rosetta-pinned"])
def test_the_pinned_probes_show_the_settings_in_force(name):
    """B6.4 step 3, once each machine has been probed again under the settings."""
    path = PROBES / f"{name}.json"
    if not path.exists():
        pytest.skip(f"{path.name} is written by the probe after the settings are pushed")
    p = json.loads(path.read_text())
    core, disabled = _pin("OPENBLAS_CORETYPE"), _pin("NPY_DISABLE_CPU_FEATURES")
    assert p["settings"]["OPENBLAS_CORETYPE"] == core
    assert p["settings"]["NPY_DISABLE_CPU_FEATURES"] == disabled
    blas = [t["architecture"] for t in p["threadpool_info"] if t.get("internal_api") == "openblas"]
    assert len(blas) == 2 and set(blas) == {core}
    assert p["numpy_cpu"]["dispatch_in_use"] == [
        d for d in p["numpy_cpu"]["dispatch"] if d not in disabled.split()
        and d in p["numpy_cpu"]["found"]]


# --------------------------------------------------------------------------
# B6.5: the required CI job runs in the recorded image
# --------------------------------------------------------------------------
BYTE_FOR_BYTE = "REPRODUCED: every output matches the committed one byte for byte."


def test_the_required_job_runs_the_recorded_image_and_demands_reproduced():
    import yaml
    jobs = yaml.safe_load((REPO / ".github/workflows/reproduce.yml").read_text())["jobs"]
    assert set(jobs) == {"image", "macos-arm64"}        # the native Linux job is gone
    image = jobs["image"]
    assert not image.get("continue-on-error")
    steps = "\n".join(s.get("run", "") for s in image["steps"])
    assert 'json.load(open("container/ENVIRONMENT.json"))["image"]["reference"]' in steps
    assert "docker pull \"$IMAGE\"" in steps
    assert "bash container/run.sh fetch" in steps
    assert "--network none" in steps and "bash container/run.sh verify" in steps
    assert f'"# {BYTE_FOR_BYTE}"' in steps                   # not REPRODUCED NUMERICALLY
    assert BYTE_FOR_BYTE in (REPO / "scripts" / "reproduce_evidence.py").read_text()


def test_the_native_macos_job_is_kept_for_the_record_only():
    import yaml
    mac = yaml.safe_load((REPO / ".github/workflows/reproduce.yml").read_text())["jobs"]["macos-arm64"]
    assert mac["continue-on-error"] is True                  # the author's decision
    assert any(s.get("run") == "bash reproduce.sh" for s in mac["steps"])


def test_the_readme_gives_the_recorded_image():
    ref = json.loads(ENVIRONMENT.read_text())["image"]["reference"]
    readme = (REPO / "README.md").read_text()
    assert f"IMAGE={ref}" in readme
    assert "bash container/run.sh fetch" in readme and "bash container/run.sh verify" in readme


# --------------------------------------------------------------------------
# PNG figures, pixel by pixel (B6, the author's wrap-up of 2026-10-22)
# --------------------------------------------------------------------------
def _pixels():
    spec = importlib.util.spec_from_file_location(
        "compare_png_pixels", REPO / "scripts" / "stage_B6" / "compare_png_pixels.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _png(array, level: int) -> bytes:
    import io
    from PIL import Image
    buf = io.BytesIO()
    Image.fromarray(array).save(buf, format="png", compress_level=level)
    return buf.getvalue()


def test_the_pixel_comparison_tells_the_three_cases_apart():
    import numpy as np
    P = _pixels()
    image = np.zeros((20, 30, 4), dtype=np.uint8)
    image[..., 3] = 255
    image[5:10, 5:20, 0] = 200
    a, b = _png(image, 9), _png(image, 1)               # the same pixels, other bytes
    other = image.copy()
    other[2, 3, 1] = 7
    c = _png(other, 9)
    assert a != b
    assert P.classify(a, a)["status"] == P.SAME_BYTES
    assert P.classify(a, b)["status"] == P.SAME_PIXELS
    r = P.classify(a, c)
    assert (r["status"], r["pixels_differing"], r["box"]) == (P.OTHER_PIXELS, 1, (3, 2, 4, 3))


def test_the_rebaselined_pngs_differ_in_their_bytes_only():
    P = _pixels()
    doc = (REPO / "docs" / "rebaseline-png-pixels.md").read_text()
    rows = dict(re.findall(r"^\| `([^`]+)` \| ([^|]+?) \|", doc, re.M))
    assert len(rows) == 13 and set(rows.values()) == {P.SAME_PIXELS}
    assert "old = commit e12220c, new = commit 7125246" in doc
    if subprocess.run(["git", "-C", str(REPO), "cat-file", "-e", "e12220c^{commit}"],
                      capture_output=True).returncode:
        pytest.skip("the two commits are not in this clone; the document is checked")
    assert sorted(rows) == P.pngs("e12220c") == P.pngs("7125246")
    for path, status in rows.items():
        assert P.classify(P.read("e12220c", path), P.read("7125246", path))["status"] == status


def test_reproduce_sh_names_the_reference_image_and_keeps_the_runs_exit_status():
    ref = json.loads(ENVIRONMENT.read_text())["image"]["reference"]
    text = (REPO / "reproduce.sh").read_text()
    for needed in (f"IMAGE={ref}", "bash container/run.sh fetch", "--network none",
                   "bash container/run.sh verify", "scripts/stage_B6/compare_png_pixels.py"):
        assert needed in text, needed
    run = '.venv/bin/python scripts/reproduce_evidence.py --fetch-inputs --verify "$@"'
    tail = text[text.index("status=0"):]
    assert run in tail
    for code in (0, 1):
        script = "set -euo pipefail\n" + tail.replace(run, f"(exit {code})")
        r = subprocess.run(["bash", "-c", script], capture_output=True, text=True)
        assert r.returncode == code and "container/run.sh verify" in r.stdout


def test_the_macos_job_compares_the_pngs_by_pixel_after_the_comparator():
    import yaml
    steps = yaml.safe_load((REPO / ".github/workflows/reproduce.yml").read_text()
                           )["jobs"]["macos-arm64"]["steps"]
    runs = [s.get("run", "") for s in steps]
    comparator = runs.index("bash reproduce.sh")
    pixels = next(i for i, r in enumerate(runs) if "compare_png_pixels.py" in r)
    assert pixels == comparator + 1
    assert steps[pixels]["if"] == "always()" and "GITHUB_STEP_SUMMARY" in runs[pixels]
