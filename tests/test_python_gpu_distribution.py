"""Guards: formal wheel GPU inventory/identity and explicit architecture compilation.

Assumes: synthetic native bytes are packaging fixtures, never GPU numerical evidence.
"""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from types import SimpleNamespace
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


@pytest.mark.parametrize("mutate_sdk", [False, True])
def test_builder_compiles_both_windows_architectures_without_ptx(tmp_path, monkeypatch, mutate_sdk):
    builder = load("build_python_gpu")
    monkeypatch.setattr(builder, "identity", lambda: "a" * 40)
    monkeypatch.setattr(builder, "admission", lambda *a: {"owner": "fixture"})
    monkeypatch.setattr(builder.platform, "system", lambda: "Windows")
    monkeypatch.setattr(builder, "check_toolchain", lambda root, *a: (root / "bin/nvcc.exe", "release 12.8, V12.8.61", "MSVC14.44"))
    sdk = tmp_path / "headers"; (sdk / "pybind11").mkdir(parents=True)
    (sdk / "pybind11/header.h").write_text("fixture")
    monkeypatch.setattr(builder, "headers", lambda: sdk)
    commands = []
    def execute(cmd, out, label, saved):
        commands.append(list(map(str, cmd)))
        if mutate_sdk and label.endswith("-cuda"):
            (matlab / "bin/win64/libdevice.10.bc").write_bytes(b"changed SDK during compile fixture")
        for word in cmd:
            if str(word).startswith("/OUT:"): Path(str(word)[5:]).write_bytes(b"not executable")
    monkeypatch.setattr(builder, "execute", execute)
    matlab = tmp_path / "matlab"
    for relative in ("sys/cuda/win64/cuda/bin/nvcc.exe", "sys/cuda/win64/cuda/bin/nvcc.profile",
                     "sys/cuda/win64/cuda/nvvm/bin/cicc.exe", "bin/win64/libdevice.10.bc",
                     "sys/cuda/win64/cuda/include/cuda.h"):
        fixture = matlab / relative
        fixture.parent.mkdir(parents=True, exist_ok=True)
        fixture.write_bytes(b"#define CUDA_VERSION 12080\n" if relative.endswith("cuda.h") else b"synthetic SDK input")
    monkeypatch.setattr(sys, "argv", ["builder", "--cuda-root", str(matlab / "sys/cuda/win64/cuda"),
        "--matlab-root", str(matlab), "--host-compiler", str(tmp_path / "cl"),
        "--architectures", "75", "80", "--output", str(tmp_path / "built"),
        "--admission", str(tmp_path / "admission"), "--admission-sha256", "a" * 64,
        "--components", "_gpu_csr"])
    if mutate_sdk:
        with pytest.raises(RuntimeError, match="SDK inputs changed during build"):
            builder.main()
        assert not list((tmp_path / "built").glob("*.build.json"))
        return
    builder.main()
    cuda = commands[0]
    assert "-gencode=arch=compute_75,code=sm_75" in cuda, cuda
    assert "-gencode=arch=compute_80,code=sm_80" in cuda, cuda
    assert not any("code=compute_" in arg for arg in cuda), cuda
    assert "/MD" in cuda, cuda
    assert cuda[0] == str(matlab / "sys/cuda/win64/cuda/bin/nvcc.exe")
    assert cuda[cuda.index("--libdevice-directory") + 1] == str(matlab / "bin/win64")
    assert "--dont-use-profile" not in cuda
    private = list((tmp_path / "built").glob("*.private-build.json"))
    assert len(private) == 1
    identity = json.loads(private[0].read_text())["cuda_sdk_identity"]
    assert identity["files"]["libdevice"]["sha256"] == builder.sha(matlab / "bin/win64/libdevice.10.bc")


@pytest.mark.parametrize("architectures", [["75"], ["80"], ["75", "75"], ["80", "75"], ["75", "80", "90"]])
def test_windows_release_rejects_incomplete_or_implicit_architectures(tmp_path, monkeypatch, architectures):
    builder = load("build_python_gpu")
    monkeypatch.setattr(builder, "capture", lambda cmd: "Cuda compilation tools, release 12.8, V12.8")
    with pytest.raises(RuntimeError, match="Explicit architectures"):
        builder.check_toolchain(tmp_path, tmp_path / "cl", architectures, True)


def test_hdiv_host_cmake_values_do_not_reinterpret_windows_path_escapes(tmp_path, monkeypatch):
    builder = load("build_python_gpu")
    compiler = tmp_path / "compiler with spaces" / "cl.exe"
    sdk = tmp_path / "installed SDK"
    sdk.mkdir()
    configs = {"ngsolve": sdk / "NGSolveConfig.cmake", "netgen-mesher": sdk / "NetgenConfig.cmake"}
    for path in configs.values():path.write_text("# fixture")
    monkeypatch.setattr(builder, "installed_config", lambda package, filename: configs[package])
    monkeypatch.setattr(builder.platform, "system", lambda: "Windows")
    import shutil
    monkeypatch.setattr(shutil, "which", lambda name: str(tmp_path / "cmake.exe"))
    out = tmp_path / "build with spaces"; out.mkdir()
    binary = out / "_hdiv_cuda.pyd"
    commands = []
    def execute(argv, output, label, saved):
        commands.append(argv)
        if label == "_hdiv_cuda-host":
            binary.write_bytes(b"not executable")
            (out / "radia_hdiv_kernels.lib").write_bytes(b"!<arch>\nnot a real archive member")
    monkeypatch.setattr(builder, "execute", execute)
    cuda = tmp_path / "cuda toolkit"
    (out / "hdiv.obj").write_bytes(b"not a real object")
    builder.build_hdiv_host(SimpleNamespace(host_compiler=compiler), binary,
                            out / "hdiv.obj", sdk, cuda, out, [])
    expected = {"CMAKE_CXX_COMPILER": compiler, "NGSolve_DIR": sdk, "Netgen_DIR": sdk,
                "Python3_EXECUTABLE": Path(sys.executable), "RADIA_PYBIND_INCLUDE": sdk,
                "RADIA_CUDA_OBJECT": out / "hdiv.obj", "RADIA_OUTPUT": out,
                **{"RADIA_"+name.upper():cuda / "lib/x64" / (name+".lib") for name in ("cudart","cublas","cusolver")}}
    values = dict(arg[2:].split("=",1) for arg in commands[0] if isinstance(arg,str) and arg.startswith("-D"))
    for key,path in expected.items():
        assert values[key] == path.as_posix() and "\\" not in values[key], (key,values[key],path)


@pytest.mark.parametrize("failure", ["same-type", "missing", "wrong-public-name"])
def test_coexistence_rejects_colliding_or_missing_public_types(failure):
    """Guards: coexistence cannot accept an aliased Resident type.

    Assumes: these are refusal fixtures; only actual native imports certify coexistence.
    """
    checker = load("check_python_gpu_coexistence")
    amg = type("Resident", (), {})
    ams = type("Resident", (), {})
    if failure == "same-type": ams = amg
    if failure == "missing": ams = None
    if failure == "wrong-public-name": ams = type("WrongName", (), {})
    modules = {"_gpu_resident_amg": SimpleNamespace(Resident=amg),
               "_gpu_ams": SimpleNamespace(Resident=ams)}
    with pytest.raises(RuntimeError, match="distinct Resident types"):
        checker.resident_types(modules)
