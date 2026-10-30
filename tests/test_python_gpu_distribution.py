"""Guards: formal wheel GPU inventory/identity and explicit architecture compilation.

Assumes: synthetic native bytes are packaging fixtures, never GPU numerical evidence.
"""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import zipfile

import pytest

ROOT = Path(__file__).resolve().parents[1]


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "tools" / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def record(name, data, source="a" * 40):
    return dict(schema="radia.native-build-provenance.v1", source_commit=source,
                source_dirty=False, binary_name=name, binary_sha256=hashlib.sha256(data).hexdigest(),
                architectures=["75", "80"], ptx_architectures=[],
                fault_injection=False, device_fault_injection=False, lifecycle_diagnostic=False,
                source_hashes={"src/core/gpu/kernel.cu": "b" * 64},
                pybind_header_hashes={"pybind11/header.h": "c" * 64})


@pytest.mark.parametrize("mutation", ["none", "missing", "source", "hash", "architecture", "diagnostic", "private"])
def test_formal_wheel_requires_complete_same_source_gpu_payload(tmp_path, mutation):
    payload = load("python_gpu_payload")
    verifier = load("verify_cpu_wheel")
    wheel = tmp_path / "radia-5.4.1-cp312-cp312-win_amd64.whl"
    with zipfile.ZipFile(wheel, "w") as z:
        z.writestr("radia-5.4.1.dist-info/METADATA", "Name: radia\nVersion: 5.4.1\n")
        core = b"synthetic CPU"
        z.writestr("radia/_radia_pybind.pyd", core)
        z.writestr("radia/_radia_pybind.pyd.build.json", json.dumps(record("_radia_pybind.pyd", core)))
        for module in payload.MODULES:
            if module == "_hdiv_cuda" and mutation == "missing":
                continue
            name = module + ".cp312-win_amd64.pyd"
            data = b"synthetic non-executable GPU"
            receipt = record(name, data)
            if module == "_hdiv_cuda":
                if mutation == "source": receipt["source_commit"] = "d" * 40
                if mutation == "hash": receipt["binary_sha256"] = "0" * 64
                if mutation == "architecture": receipt["architectures"] = ["75"]
                if mutation == "diagnostic": receipt["fault_injection"] = True
            z.writestr("radia/" + name, data)
            z.writestr("radia/" + name + ".build.json", json.dumps(receipt))
        if mutation == "private": z.writestr("radia/commands.json", "private original")
    if mutation == "none":
        result = verifier.verify(wheel, require_python_gpu=True)
        assert len(result["python_gpu_binaries"]) == 4, result
        assert result["runtime_qualified"] is False, result
    else:
        with pytest.raises(ValueError): verifier.verify(wheel, require_python_gpu=True)


def test_cpu_only_wheel_cannot_satisfy_gpu_release(tmp_path):
    wheel = tmp_path / "radia-5.4.1-cp312-cp312-win_amd64.whl"
    with zipfile.ZipFile(wheel, "w") as z:
        z.writestr("radia-5.4.1.dist-info/METADATA", "Name: radia\nVersion: 5.4.1\n")
    with pytest.raises((ValueError, KeyError)):
        load("verify_cpu_wheel").verify(wheel, require_python_gpu=True)


def test_stage_validates_all_members_before_copying(tmp_path):
    payload = load("python_gpu_payload")
    out, package = tmp_path / "out", tmp_path / "package"
    out.mkdir(); package.mkdir()
    for module in payload.MODULES:
        binary = out / (module + ".cp312-win_amd64.pyd")
        binary.write_bytes(b"synthetic")
        receipt = record(binary.name, binary.read_bytes())
        if module == "_hdiv_cuda": receipt["source_commit"] = "b" * 40
        Path(str(binary) + ".build.json").write_text(json.dumps(receipt))
    with pytest.raises(ValueError): payload.stage(out, package, "a" * 40, ".cp312-win_amd64.pyd")
    assert list(package.iterdir()) == [], list(package.iterdir())


def test_builder_compiles_both_windows_architectures_without_ptx(tmp_path, monkeypatch):
    builder = load("build_python_gpu")
    monkeypatch.setattr(builder, "identity", lambda: "a" * 40)
    monkeypatch.setattr(builder, "admission", lambda *a: {"owner": "fixture"})
    monkeypatch.setattr(builder.platform, "system", lambda: "Windows")
    monkeypatch.setattr(builder, "check_toolchain", lambda *a: (tmp_path / "nvcc", "CUDA12.8", "MSVC14.44"))
    sdk = tmp_path / "headers"; (sdk / "pybind11").mkdir(parents=True)
    (sdk / "pybind11/header.h").write_text("fixture")
    monkeypatch.setattr(builder, "headers", lambda: sdk)
    commands = []
    def execute(cmd, out, label, saved):
        commands.append(list(map(str, cmd)))
        for word in cmd:
            if str(word).startswith("/OUT:"): Path(str(word)[5:]).write_bytes(b"not executable")
    monkeypatch.setattr(builder, "execute", execute)
    matlab = tmp_path / "matlab"
    monkeypatch.setattr(sys, "argv", ["builder", "--cuda-root", str(matlab / "sys/cuda/win64/cuda"),
        "--matlab-root", str(matlab), "--host-compiler", str(tmp_path / "cl"),
        "--architectures", "75", "80", "--output", str(tmp_path / "built"),
        "--admission", str(tmp_path / "admission"), "--admission-sha256", "a" * 64,
        "--components", "_gpu_csr"])
    builder.main()
    cuda = commands[0]
    assert "-gencode=arch=compute_75,code=sm_75" in cuda, cuda
    assert "-gencode=arch=compute_80,code=sm_80" in cuda, cuda
    assert not any("code=compute_" in arg for arg in cuda), cuda
    assert "/MD" in cuda, cuda


@pytest.mark.parametrize("architectures", [["75"], ["80"], ["75", "75"], ["80", "75"], ["75", "80", "90"]])
def test_windows_release_rejects_incomplete_or_implicit_architectures(tmp_path, monkeypatch, architectures):
    builder = load("build_python_gpu")
    monkeypatch.setattr(builder, "capture", lambda cmd: "Cuda compilation tools, release 12.8, V12.8")
    with pytest.raises(RuntimeError, match="Explicit architectures"):
        builder.check_toolchain(tmp_path, tmp_path / "cl", architectures, True)
