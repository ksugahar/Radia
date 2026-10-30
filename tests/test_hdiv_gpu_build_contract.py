"""Guards: Radia-owned HDiv build closure and host/device separation.

Assumes: command fixtures are not native/GPU numerical validation.
"""
from pathlib import Path
import importlib.util
import json
import sys
import types
import pytest
ROOT = Path(__file__).resolve().parents[1]


def builder():
    spec = importlib.util.spec_from_file_location("hdiv_builder_contract", ROOT / "tools/build_python_gpu.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_hdiv_host_build_uses_installed_ngsolve_and_exact_pip_cuda(tmp_path, monkeypatch):
    b = builder()
    configs = {}
    for package, name in [("ngsolve", "NGSolveConfig.cmake"), ("netgen-mesher", "NetgenConfig.cmake")]:
        p = tmp_path / package / name
        p.parent.mkdir(); p.write_text("# independent SDK fixture")
        configs[package] = p
    monkeypatch.setattr(b, "installed_config", lambda package, name: configs[package])
    monkeypatch.setattr(b.platform, "system", lambda: "Linux")
    monkeypatch.setattr(b, "cuda_libraries", lambda *a: {n: str(tmp_path / "pip" / (n + ".so")) for n in ("cudart", "cublas", "cusolver")})
    import shutil
    monkeypatch.setattr(shutil, "which", lambda name: "/fixture/cmake")
    out = tmp_path / "out"; out.mkdir(); binary = out / "_hdiv_cuda.so"
    commands = []
    def execute(command, output, label, recorded):
        recorded.append(list(map(str, command)))
        if "--build" in command:
            binary.write_bytes(b"synthetic artifact; never executed")
    monkeypatch.setattr(b, "execute", execute)
    b.build_hdiv_host(types.SimpleNamespace(host_compiler=Path("/explicit/g++")), binary,
                      out / "cuda.o", tmp_path / "headers", tmp_path / "toolkit", out, commands)
    configure, compile_host = commands
    assert "-DCMAKE_CXX_COMPILER=" + Path("/explicit/g++").as_posix() in configure
    assert "-DPython3_EXECUTABLE=" + Path(sys.executable).as_posix() in configure
    assert "-DNGSolve_DIR=" + configs["ngsolve"].parent.as_posix() in configure
    assert "-DNetgen_DIR=" + configs["netgen-mesher"].parent.as_posix() in configure
    assert all("-DRADIA_" + name.upper() + "=" + (tmp_path / "pip" / (name + ".so")).as_posix() in configure for name in ("cudart", "cublas", "cusolver"))
    assert compile_host[-2:] == ["--parallel", "2"]
    assert not any("matlab" in word.lower() or "radia-simulink" in word for command in commands for word in command)
    assert set(json.loads((out / "hdiv-sdk-private.json").read_text())) == {str(p) for p in configs.values()}


def test_missing_or_ambiguous_sdk_is_rejected(tmp_path, monkeypatch):
    b = builder()
    files = []
    monkeypatch.setattr(b.metadata, "distribution", lambda _: types.SimpleNamespace(files=files, locate_file=lambda x: tmp_path / x))
    with pytest.raises(RuntimeError, match="Missing or ambiguous"):
        b.installed_config("ngsolve", "NGSolveConfig.cmake")
    for directory in ("one", "two"):
        p = tmp_path / directory / "NGSolveConfig.cmake"; p.parent.mkdir(); p.write_text("fixture")
        files.append(p.relative_to(tmp_path))
    with pytest.raises(RuntimeError, match="Missing or ambiguous"):
        b.installed_config("ngsolve", "NGSolveConfig.cmake")


def test_python_entry_is_lazy_and_not_a_simulink_dependency():
    spec = importlib.util.spec_from_file_location("hdiv_cold_entry", ROOT / "src/radia/hdiv_gpu.py")
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    assert callable(module.PreparedHDivNonlinear)
    assert callable(module.prepare_hdiv_host_setup)
    with pytest.raises(ValueError, match="explicit backend"):
        module.PreparedHDivNonlinear("unused.vol", {}, backend="cpu", label_contract="unused.json", mesh_check_receipt="unused.receipt", mesh_check_receipt_sha256="a"*64, mesh_check_provenance={})
    with pytest.raises(ValueError, match="explicit backend"):
        module.prepare_hdiv_host_setup("unused.vol", backend="cpu", label_contract="unused.json", mesh_check_receipt="unused.receipt", mesh_check_receipt_sha256="a"*64, mesh_check_provenance={})


def test_cuda_translation_unit_does_not_pull_fe_or_mex_headers():
    import re
    directory = ROOT / "src/core/gpu/hdiv"
    pending = [directory / "hdiv_resident_python_cuda.cu"]; seen = set()
    while pending:
        p = pending.pop()
        if p in seen: continue
        seen.add(p); text = p.read_text()
        for include in re.findall(r'^#include [<"]([^">]+)[">]', text, re.M):
            assert include not in {"mex.h", "comp.hpp", "hdivfe.hpp", "meshaccess.hpp"}, (p, include)
            if (p.parent / include).is_file():pending.append(p.parent / include)
    assert len(seen) > 4
    assert not (ROOT / "packages/radia-simulink/native/hdiv_resident_python.cpp").exists()


def test_headers_follow_wheel_metadata_layout(tmp_path, monkeypatch):
    b=builder()
    sdk=tmp_path/"include/netgen"
    header=sdk/"pybind11/pybind11.h";header.parent.mkdir(parents=True);header.write_text("fixture")
    monkeypatch.setitem(sys.modules,"netgen",types.SimpleNamespace(__file__=str(tmp_path/"netgen/__init__.py")))
    monkeypatch.setitem(sys.modules,"ngsolve",types.SimpleNamespace())
    monkeypatch.setattr(b.metadata,"version",lambda _: "6.2.2608")
    monkeypatch.setattr(b.metadata,"distribution",lambda _:types.SimpleNamespace(files=[header.relative_to(tmp_path)],locate_file=lambda f:tmp_path/f))
    assert b.headers()==sdk.resolve()
