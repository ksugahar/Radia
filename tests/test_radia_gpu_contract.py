"""GPU adapter and runtime source contracts.

Guards: Explicit backend selection, exact pip versions/RECORD bytes, conflicting
ELF providers, and resident-handle close/reject behavior without GPU fallback.
Assumes: CUDA loading and installed metadata are mocked. Temporary files are
fixtures, not native libraries. No compilation, CUDA context or solve executes.
"""
import base64
import hashlib
import importlib.util
from pathlib import Path
from types import SimpleNamespace
import tomllib

import pytest

ROOT = Path(__file__).resolve().parents[1]


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "src/radia" / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_exact_gpu_pins_are_read_from_installed_metadata(monkeypatch):
    runtime = load("_cuda_runtime")
    requirements = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]["optional-dependencies"]["gpu"]
    monkeypatch.setattr(runtime.metadata, "requires", lambda _: requirements)
    pins = runtime._pins()
    assert set(pins) == {name for name, _ in runtime._LIBRARIES}
    monkeypatch.setattr(runtime.metadata, "requires", lambda _: requirements[:-1])
    with pytest.raises(ImportError, match="one exact pin"):
        runtime._pins()


def test_library_record_hash_and_version_are_enforced(tmp_path, monkeypatch):
    runtime = load("_cuda_runtime")
    path = tmp_path / "libtest.so.12"
    path.write_bytes(b"fixture library bytes")
    expected = base64.urlsafe_b64encode(hashlib.sha256(path.read_bytes()).digest()).decode().rstrip("=")
    class File:
        hash = SimpleNamespace(mode="sha256", value=expected)
        def __str__(self):
            return "libtest.so.12"
    dist = SimpleNamespace(version="12.8.1", files=[File()], locate_file=lambda _: path)
    monkeypatch.setattr(runtime, "_LIBRARIES", (("cuda-fixture", ("libtest.so.12",)),))
    monkeypatch.setattr(runtime.metadata, "distribution", lambda _: dist)
    assert runtime._library_paths({"cuda-fixture": "12.8.1"})[path.name]["path"] == str(path)
    dist.version = "13.3.0"
    with pytest.raises(ImportError, match="required 12.8.1, found 13.3.0"):
        runtime._library_paths({"cuda-fixture": "12.8.1"})
    dist.version = "12.8.1"
    path.write_bytes(b"tampered")
    with pytest.raises(ImportError, match="RECORD mismatch"):
        runtime._library_paths({"cuda-fixture": "12.8.1"})


def test_conflicting_site_elf_is_rejected_before_loading(monkeypatch):
    runtime = load("_cuda_runtime")
    expected = Path("/pip/libcudart.so.12")
    libraries = {expected.name: {"path": str(expected)}}
    monkeypatch.setattr(runtime, "_maps", lambda: {Path("/site/libcudart.so.13")})
    with pytest.raises(ImportError, match="Conflicting CUDA"):
        runtime._check_loaded(libraries)
    monkeypatch.setattr(runtime, "_maps", lambda: set())
    with pytest.raises(ImportError, match="could not be verified"):
        runtime._check_loaded(libraries, require_all=True)
    monkeypatch.setattr(runtime, "_maps", lambda: {expected})
    runtime._check_loaded(libraries, require_all=True)


def test_unsupported_backend_does_not_load_native(monkeypatch):
    gpu = load("gpu")
    monkeypatch.setattr(gpu, "_native", lambda _: pytest.fail("must not load GPU"))
    with pytest.raises(ValueError, match="explicitly cpu or gpu"):
        gpu.solve_spd(None, None, backend="auto")
    with pytest.raises(ValueError, match="explicit backend"):
        gpu.ResidentAMG([], None, backend="cpu")


def test_resident_handle_closes_once_and_never_falls_back(monkeypatch):
    gpu = load("gpu")
    calls = []
    class Native:
        def solve(self, *args):
            calls.append(("solve", args))
            return args[0].copy(), {"backend": "gpu"}
        def close(self):
            calls.append(("close",))
    owner = Native()
    monkeypatch.setattr(gpu, "_native", lambda _: SimpleNamespace(Resident=lambda *_: owner))
    import numpy as np
    levels = [{"a": {"n": 1, "m": 1, "row": np.array([0, 1], dtype=np.int32),
                       "col": np.array([0], dtype=np.int32), "val": np.array([1.])}}]
    with gpu.ResidentAMG(levels, None, backend="gpu") as context:
        x, info = context.solve(np.array([2.]), max_iterations=100)
        assert x.tolist() == [2.] and info["residual_checked"]
        with pytest.raises(ValueError, match="tolerance"):
            context.solve("rhs", tolerance=float("nan"))
    context.close()
    assert [item[0] for item in calls] == ["solve", "close"]
    with pytest.raises(RuntimeError, match="closed"):
        context.solve("rhs")


def test_build_entry_rejects_unreviewed_toolchain_and_admission(tmp_path, monkeypatch):
    spec = importlib.util.spec_from_file_location("gpu_builder", ROOT / "tools/build_python_gpu.py")
    builder = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(builder)
    monkeypatch.setattr(builder, "capture", lambda _: "Cuda compilation tools, release 13.3, V13.3.0")
    with pytest.raises(RuntimeError, match="12.8"):
        builder.check_toolchain(tmp_path, tmp_path / "g++", "100", False)
    record = tmp_path / "admission.json"
    record.write_text('{"source_sha":"wrong", "owner":"CodexFlow2", "gpu_build_approved":true}')
    with pytest.raises(RuntimeError, match="exact source"):
        builder.admission(record, builder.sha(record), "a" * 40)


def test_shared_source_closure_resolves_without_legacy_search():
    import re
    core = ROOT / "src/core/gpu"
    for source in core.iterdir():
        if source.is_file():
            for include in re.findall(r'^#include "([^"\n]+)"', source.read_text(), re.M):
                assert (source.parent / include).is_file(), (source, include)
    native = ROOT / "packages/radia-simulink/native"
    for name in ("radia_cuda_ams.cu", "radia_cuda_csr_pcg.cu"):
        source = native / name
        includes = re.findall(r'^#include "([^"\n]+)"', source.read_text(), re.M)
        moved = [name for name in includes if "src/core/gpu/" in name]
        assert moved and all((source.parent / name).is_file() for name in moved), source
    assert not (native / "ams_amg_kernels.cuh").exists()
    assert (core / "ams_amg_kernels.cuh").is_file()


def test_windows_rechecks_record_bytes_on_warm_calls(tmp_path, monkeypatch):
    import sys
    runtime, windows = load("_cuda_runtime"), load("_cuda_runtime_windows")
    windows.__package__ = "gpu_runtime_fixture"
    monkeypatch.setitem(sys.modules, "gpu_runtime_fixture._cuda_runtime", runtime)
    monkeypatch.setattr(runtime, "_pins", lambda: {"cuda-fixture": "12.8.1"})
    path = tmp_path / "cudart64_12.dll"
    path.write_bytes(b"synthetic DLL payload")
    expected = base64.urlsafe_b64encode(hashlib.sha256(path.read_bytes()).digest()).decode().rstrip("=")
    class File:
        hash = SimpleNamespace(mode="sha256", value=expected)
        def __str__(self):
            return path.name
    dist = SimpleNamespace(version="12.8.1", files=[File()], locate_file=lambda _: path)
    monkeypatch.setattr(runtime.metadata, "distribution", lambda _: dist)
    monkeypatch.setattr(windows.sys, "platform", "win32")
    monkeypatch.setattr(windows, "_loaded_modules", lambda: {path})
    monkeypatch.setattr(windows.os, "add_dll_directory", lambda _: SimpleNamespace(close=lambda: None), raising=False)
    windows.prepare_cuda_runtime()
    windows.prepare_cuda_runtime()
    assert len(windows._handles) == 1
    monkeypatch.setattr(windows, "_loaded_modules", lambda: {tmp_path / "cudart64_13.dll"})
    with pytest.raises(ImportError, match="Conflicting CUDA DLL"):
        windows.prepare_cuda_runtime()
    monkeypatch.setattr(windows, "_loaded_modules", lambda: {path})
    path.write_bytes(b"tampered after the first load")
    with pytest.raises(ImportError, match="RECORD mismatch"):
        windows.prepare_cuda_runtime()


def test_public_build_receipt_excludes_private_admission_and_commands(tmp_path, monkeypatch):
    import json
    import sys
    spec = importlib.util.spec_from_file_location("gpu_receipt_builder", ROOT / "tools/build_python_gpu.py")
    builder = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(builder)
    out, matlab = tmp_path / "output", tmp_path / "private-matlab-location"
    admission = tmp_path / "private-admission.json"
    admission.write_text('{"private_secret":"do-not-publish"}')
    sha = "a" * 40
    monkeypatch.setattr(builder, "identity", lambda: sha)
    monkeypatch.setattr(builder, "admission", lambda *_: {"private_secret": "do-not-publish"})
    monkeypatch.setattr(builder.platform, "system", lambda: "Windows")
    monkeypatch.setattr(builder, "check_toolchain", lambda *_: (tmp_path / "nvcc", "CUDA12.8", "MSVC14.44"))
    sdk = tmp_path / "sdk"
    (sdk / "pybind11").mkdir(parents=True)
    (sdk / "pybind11/header.h").write_text("synthetic header")
    monkeypatch.setattr(builder, "headers", lambda: sdk)
    def compile_fixture(command, output, label, commands):
        commands.append(list(map(str, command)))
        for item in command:
            if str(item).startswith("/OUT:"):
                Path(str(item)[5:]).write_bytes(b"synthetic non-executable payload")
    monkeypatch.setattr(builder, "execute", compile_fixture)
    monkeypatch.setattr(sys, "argv", ["build_python_gpu.py", "--cuda-root", str(matlab / "sys/cuda/win64/cuda"),
                        "--matlab-root", str(matlab), "--host-compiler", str(tmp_path / "cl.exe"),
                        "--architecture", "75", "--output", str(out), "--admission", str(admission),
                        "--admission-sha256", builder.sha(admission)])
    builder.main()
    receipts = list(out.glob("*.build.json"))
    assert len(receipts) == len(builder.COMPONENTS) == 3
    for receipt in receipts:
        record = json.loads(receipt.read_text())
        assert "admission" not in record and "commands" not in record
        assert "do-not-publish" not in receipt.read_text()
        assert str(matlab) not in receipt.read_text()
        private = Path(str(receipt).removesuffix(".build.json") + ".private-build.json")
        assert record["private_receipt_sha256"] == builder.sha(private)
        assert "do-not-publish" in private.read_text()
    assert "do-not-publish" in (out / "admission-private.json").read_text()


@pytest.mark.parametrize("result", [
    ([2.], {"backend": "cpu"}),
    ([99.], {"backend": "gpu", "true_relative_residual": 0.}),
])
def test_resident_amg_rejects_backend_and_false_native_residual(monkeypatch, result):
    import numpy as np
    gpu = load("gpu")
    class Native:
        def solve(self, *args):
            return np.array(result[0]), result[1]
        def close(self):
            pass
    monkeypatch.setattr(gpu, "_native", lambda _: SimpleNamespace(Resident=lambda *_: Native()))
    values = np.array([1.])
    levels = [{"a": {"n": 1, "m": 1, "row": np.array([0, 1], dtype=np.int32),
                      "col": np.array([0], dtype=np.int32), "val": values}}]
    with gpu.ResidentAMG(levels, None, backend="gpu") as context:
        values[:] = 2. / 99.  # Caller mutation must not change the reference equation.
        with pytest.raises(RuntimeError, match="backend|residual"):
            context.solve(np.array([2.]))


def test_native_core_receipt_missing_raises_import_error(tmp_path, monkeypatch):
    import sys
    import sysconfig
    gpu = load("gpu")
    monkeypatch.setattr(gpu, "__file__", str(tmp_path / "gpu.py"))
    binary = tmp_path / ("_gpu_csr" + sysconfig.get_config_var("EXT_SUFFIX"))
    binary.write_bytes(b"fixture")
    Path(str(binary) + ".build.json").write_text("{}")
    core = SimpleNamespace(__file__=str(tmp_path / "_radia_pybind.pyd"))
    monkeypatch.setitem(sys.modules, "radia", SimpleNamespace(_radia_pybind=core))
    monkeypatch.setitem(sys.modules, "radia._radia_pybind", core)
    with pytest.raises(ImportError, match="CPU artifact build provenance"):
        gpu._native("_gpu_csr")
