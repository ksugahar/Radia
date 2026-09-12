"""Source-level contracts for the MATLAB/MEX runtime boundary."""

from __future__ import annotations

import tomllib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _source(path: str) -> str:
    return (ROOT / path).read_text(encoding="utf-8")


def test_matlab_setup_checks_the_native_ngsolve_abi_and_mkl_dispatcher() -> None:
    setup = _source("matlab/+radia/setup.m")
    expected = _source("matlab/+radia/+internal/expectedNGSolveVersion.m")
    with (ROOT / "pyproject.toml").open("rb") as stream:
        dependencies = tomllib.load(stream)["project"]["dependencies"]

    assert 'value = "6.2.2606";' in expected
    assert "ngsolve==6.2.2606" in dependencies
    assert "netgen-mesher==6.2.2606" in dependencies
    assert 'setenv("MKL_THREADING_LAYER", "SEQUENTIAL")' in setup
    assert '"mkl_threading_layer_requested"' in setup
    assert '"mkl_threading_layer",' not in setup
    assert "RADIA_NGSOLVE_VERSION:" in setup
    assert "RADIA_NETGEN_VERSION:" in setup
    assert "radia.internal.expectedNGSolveVersion()" in setup
    assert "radia:setup:NGSolveABI" in setup


def test_matlab_release_uses_sequential_mkl_without_changing_python() -> None:
    build = _source("Build.ps1")
    cmake = _source("CMakeLists.txt")
    package = _source("tools/package_simulink_release.py")
    verifier = _source("tools/verify_simulink_release.py")

    assert build.index("$PythonLibrary") < build.index("$env:MKLROOT")
    assert 'if ($MatlabMexOnly)' in build
    assert '"mkl_sequential"' in build
    assert "mkl_intel_thread" in cmake
    assert "mkl_sequential" in cmake
    assert "IN LISTS MKL_MATLAB_RUNTIME_DLLS" in cmake
    assert "mkl_sequential.3.dll" in package
    assert "mkl_intel_thread.3.dll" not in package
    assert "+radia/+internal/expectedNGSolveVersion.m" in package
    assert "matlab/mkl_sequential.3.dll" in verifier
    assert "matlab/+radia/+internal/expectedNGSolveVersion.m" in verifier


def test_release_schema_preserves_historical_threaded_mkl_contracts() -> None:
    import importlib.util

    path = ROOT / "tools" / "verify_simulink_release.py"
    spec = importlib.util.spec_from_file_location("verify_runtime_contract", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)

    for historical in (
        module.PREVIEW_V2_REQUIRED_MEMBERS,
        module.FULL_REQUIRED_MEMBERS_V2,
        module.FULL_REQUIRED_MEMBERS_V3,
    ):
        assert "matlab/mkl_intel_thread.3.dll" in historical
        assert "matlab/mkl_sequential.3.dll" not in historical
        assert "matlab/+radia/+internal/expectedNGSolveVersion.m" not in historical
    for current in (module.REQUIRED_MEMBERS, module.FULL_REQUIRED_MEMBERS):
        assert "matlab/mkl_sequential.3.dll" in current
        assert "matlab/mkl_intel_thread.3.dll" not in current
        assert "matlab/+radia/+internal/expectedNGSolveVersion.m" in current
