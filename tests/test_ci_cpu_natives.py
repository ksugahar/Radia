"""CPU natives cross from the CPU build job to the Windows GPU job only as recorded.

Guards: the Windows GPU wheel job stages exactly the file set the CPU build job
recorded, byte for byte, with a clean core receipt for the same commit, and
never accepts CUDA modules, private build records or unlisted files.
Assumes: tools/ci_cpu_natives.py is the only handoff used by build-test.yml.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
SHA = "a" * 40


def _tool():
    spec = importlib.util.spec_from_file_location("ci_cpu_natives", ROOT / "tools/ci_cpu_natives.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _build_tree(root: Path, sha=SHA, dirty=False):
    root.mkdir()
    core = root / "_radia_pybind.pyd"
    core.write_bytes(b"core")
    (root / "_radia_pybind.pyd.build.json").write_text(json.dumps({
        "source_commit": sha, "source_dirty": dirty,
        "binary_sha256": hashlib.sha256(b"core").hexdigest()}), encoding="utf-8")
    (root / "radia_motor_rom.dll").write_bytes(b"rom")
    (root / "peec_matrices.pyd").write_bytes(b"peec")
    (root / "_radia_pybind.locked-old.pyd").write_bytes(b"old")
    (root / "__init__.py").write_text("", encoding="utf-8")
    return root


def _recorded(tmp_path):
    tool = _tool()
    out = tmp_path / "artifact"
    tool.record(_build_tree(tmp_path / "src"), out, SHA)
    dest = tmp_path / "dest"
    dest.mkdir()
    return tool, out, dest


def test_record_then_stage_copies_exactly_the_listed_files(tmp_path):
    tool, out, dest = _recorded(tmp_path)
    names = {e["name"] for e in json.loads((out / tool.MANIFEST).read_text())["files"]}
    assert names == {"_radia_pybind.pyd", "_radia_pybind.pyd.build.json",
                     "radia_motor_rom.dll", "peec_matrices.pyd"}, names
    tool.stage(out, dest, SHA, tmp_path / "evidence.json")
    assert {p.name for p in dest.iterdir()} == names
    assert (dest / "peec_matrices.pyd").read_bytes() == b"peec"


@pytest.mark.parametrize("defect", ["unlisted", "missing", "bytes", "sha", "cuda", "private", "dest"])
def test_stage_rejects_any_deviation_from_the_record(tmp_path, defect):
    tool, out, dest = _recorded(tmp_path)
    sha = SHA
    if defect == "unlisted":
        (out / "extra.dll").write_bytes(b"x")
    elif defect == "missing":
        (out / "peec_matrices.pyd").unlink()
    elif defect == "bytes":
        (out / "peec_matrices.pyd").write_bytes(b"PEEC")
    elif defect == "sha":
        sha = "b" * 40
    elif defect == "cuda":
        (out / "_radia_hdiv_cuda.pyd").write_bytes(b"gpu")
        manifest = json.loads((out / tool.MANIFEST).read_text())
        manifest["files"].append({"name": "_radia_hdiv_cuda.pyd", "bytes": 3,
                                  "sha256": hashlib.sha256(b"gpu").hexdigest()})
        (out / tool.MANIFEST).write_text(json.dumps(manifest))
    elif defect == "private":
        (out / "x.pyd.private-build.json").write_text("{}")
    elif defect == "dest":
        (dest / "_radia_pybind.pyd").write_bytes(b"stale")
    with pytest.raises(ValueError):
        tool.stage(out, dest, sha, tmp_path / "evidence.json")


@pytest.mark.parametrize("defect", ["dirty", "other-commit", "no-rom"])
def test_record_requires_a_clean_core_receipt_and_required_members(tmp_path, defect):
    tool = _tool()
    src = _build_tree(tmp_path / "src", sha="c" * 40 if defect == "other-commit" else SHA,
                      dirty=defect == "dirty")
    if defect == "no-rom":
        (src / "radia_motor_rom.dll").unlink()
    with pytest.raises(ValueError):
        tool.record(src, tmp_path / "artifact", SHA)


def test_native_release_uses_the_recorded_handoff_on_both_sides():
    jobs = yaml.safe_load((ROOT / ".github/workflows/build-test.yml").read_text(encoding="utf-8"))["jobs"]
    cpu = "\n".join(step.get("run", "") for step in jobs["build-test"]["steps"])
    gpu = "\n".join(step.get("run", "") for step in jobs["windows-gpu-wheel"]["steps"])
    assert "ci_cpu_natives.py record" in cpu
    assert "ci_cpu_natives.py stage" in gpu
    assert "Copy-Item" not in "\n".join(
        step.get("run", "") for step in jobs["windows-gpu-wheel"]["steps"]
        if "radia-cpu-pyd" in step.get("run", ""))
    upload = next(step for step in jobs["build-test"]["steps"]
                  if step.get("with", {}).get("name") == "radia-cpu-pyd")
    assert upload["with"]["path"] == "${{ env.RADIA_CPU_NATIVES_DIR }}"
