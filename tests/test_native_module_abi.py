"""Regression tests for the source-tree native ABI probe."""

from __future__ import annotations

import _ctypes
import importlib.util
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CHECKER_PATH = ROOT / "tools" / "check_native_module_abi.py"


def _load_checker():
    spec = importlib.util.spec_from_file_location(
        "check_native_module_abi", CHECKER_PATH
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_probe_accepts_a_native_extension_from_the_requested_directory():
    checker = _load_checker()
    native_root = Path(_ctypes.__file__).resolve().parent

    result = checker.probe_module("_ctypes", native_root)

    assert Path(result["path"]).resolve() == Path(_ctypes.__file__).resolve()
    assert result["python"] == sys.version.split()[0]


def test_probe_rejects_a_python_shim_in_place_of_a_native_module(tmp_path):
    package = tmp_path / "demo"
    package.mkdir()
    (package / "__init__.py").write_text("", encoding="ascii")
    (package / "native.py").write_text("value = 1\n", encoding="ascii")

    completed = subprocess.run(
        [
            sys.executable,
            str(CHECKER_PATH),
            "--source-root",
            str(tmp_path),
            "--module",
            "demo.native",
        ],
        check=False,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )

    assert completed.returncode == 1
    assert "not a native extension" in completed.stderr


def test_build_contract_syncs_and_load_checks_peec_matrices():
    cmake = (ROOT / "CMakeLists.txt").read_text(encoding="utf-8")
    build = (ROOT / "Build.ps1").read_text(encoding="utf-8")

    peec_section = cmake[
        cmake.index("pybind11_add_module(peec_matrices") : cmake.index(
            'message(STATUS "peec_matrices module configured")'
        )
    ]
    assert '"${CMAKE_SOURCE_DIR}/src/radia/peec_matrices.pyd"' in peec_section
    assert "ERROR: peec_matrices build failed" in build
    assert 'dst = "peec_matrices.pyd"; required = $true' in build
    assert '"radia.peec_matrices"' in build
    assert "check_native_module_abi.py" in build


def test_native_release_keeps_generated_context_outside_source_tree():
    import yaml

    workflow = yaml.safe_load(
        (ROOT / ".github/workflows/build-test.yml").read_text(encoding="utf-8")
    )
    job = workflow["jobs"]["build-test"]
    steps = {step.get("name"): step for step in job["steps"]}
    context = steps["Record exact CI ref context"]["run"]
    assert "Join-Path $env:RADIA_CI_OUTPUT_DIR 'ci-release-context.json'" in context
    assert "Join-Path $env:GITHUB_WORKSPACE 'ci-release-context.json'" not in context
    assert steps["Upload exact CI ref context"]["with"]["path"] == (
        "${{ env.RADIA_CI_OUTPUT_DIR }}/ci-release-context.json"
    )
    build = steps["Build with MSVC + MKL"]["run"]
    assert "@rebuildFlag -RequireNativeProvenance" in build
    assert "memory-before-build.json" in build
    assert "memory-after-build.json" in build
    memory = steps["Upload build memory diagnostics"]
    assert memory["if"] == "always()"
    assert memory["with"]["path"] == (
        "${{ env.RADIA_CI_OUTPUT_DIR }}/memory-*-build.json"
    )
    for key, value in {
        "OMP_NUM_THREADS": "8", "NGS_NUM_THREADS": "8",
        "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
    }.items():
        assert job["env"][key] == value
