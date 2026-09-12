"""Contract checks for the persisted pybind11/MEX transfer-map study."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "validation_test" / "accelerator" / "benchmark_beam_transfer_mex_python.py"
SPEC = importlib.util.spec_from_file_location("beam_transfer_validation", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def _report(backend: str) -> dict:
    return {
        "schema": MODULE.RESULT_SCHEMA,
        "backend": backend,
        "case_id": "nonlinear-three-region-transfer-map",
        "case_sha256": "a" * 64,
        "radia_git_head": "b" * 40,
        "binary": {
            "path": f"C:/temp/{backend}.pyd",
            "bytes": 123,
            "sha256": "c" * 64,
            "source_commit": "b" * 40,
            "build_manifest": f"C:/temp/{backend}.pyd.build.json",
        },
        "repeats": 3,
        "first_s": 0.5,
        "median_s": 0.25,
        "min_s": 0.2,
        "observables": {"R_fro": 1.0, "T_211": 0.6, "U_3111": 1.08},
    }


def test_case_contract_builds_expected_sparse_tensors():
    case_path = ROOT / "validation_test" / "accelerator" / "beam_transfer_benchmark_case.json"
    raw, value, lengths, a, f2, f3 = MODULE._case(case_path)

    assert len(raw) > 0
    assert value["schema"] == MODULE.CASE_SCHEMA
    assert lengths.shape == (3,)
    assert a.shape == (3, 6, 6)
    assert f2[0, 1, 0, 0] == pytest.approx(2.0)
    assert f3[2, 3, 0, 0, 0] == pytest.approx(-0.7)


def test_fast_ci_declares_numpy_for_array_contracts():
    import yaml

    workflow = yaml.safe_load((ROOT / ".github/workflows/radia-fast.yml").read_text())
    steps = workflow["jobs"]["fast-contracts"]["steps"]
    setup = next(step for step in steps if step.get("name") == "Create isolated fast-CI environment")
    install = next(line for line in setup["run"].splitlines() if "pip install" in line)
    assert "numpy" in install.split()


def test_pybind_lane_fails_loudly_without_matching_source_artifact(tmp_path):
    with pytest.raises(RuntimeError, match="installed-wheel fallback is forbidden"):
        MODULE._load_source_backend(tmp_path)

    stale = tmp_path / "stale" / "src" / "radia" / "_radia_pybind.pyd"
    stale.parent.mkdir(parents=True)
    stale.write_bytes(b"old-binary-must-not-be-imported")
    with pytest.raises(RuntimeError, match="provenance manifest is missing"):
        MODULE._load_source_backend(tmp_path / "stale")


def test_build_manifest_must_match_binary_bytes_and_selected_commit(tmp_path, monkeypatch):
    source_root = tmp_path / "source"
    binary = source_root / "src" / "radia" / "_radia_pybind.pyd"
    binary.parent.mkdir(parents=True)
    binary.write_bytes(b"native-v1")
    manifest = {
        "schema": MODULE.PROVENANCE_SCHEMA,
        "source_commit": "a" * 40,
        "source_dirty": False,
        "source_change_fingerprint_sha256": "d" * 64,
        "binary_name": binary.name,
        "binary_bytes": binary.stat().st_size,
        "binary_sha256": MODULE._sha256(binary),
    }
    manifest_path = Path(f"{binary}.build.json")
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    monkeypatch.setattr(MODULE, "_git_head", lambda root=MODULE.ROOT: "a" * 40)

    assert MODULE._load_build_provenance(binary, source_root)["source_commit"] == "a" * 40
    binary.write_bytes(b"native-v2")
    with pytest.raises(ValueError, match="binary_sha256"):
        MODULE._load_build_provenance(binary, source_root)


def test_mocked_provenance_and_comparison_contract(tmp_path):
    python = _report("python-pybind11")
    matlab = _report("matlab-mex")
    matlab["binary"]["bytes"] = 123.0  # MATLAB JSON represents doubles as numbers.
    python_path = tmp_path / "python.json"
    matlab_path = tmp_path / "matlab.json"
    python_path.write_text(json.dumps(python), encoding="utf-8")
    matlab_path.write_text(json.dumps(matlab), encoding="utf-8")

    comparison = MODULE.compare(python_path, matlab_path, rtol=1e-11, atol=1e-12)
    assert comparison["pass"], comparison
    assert comparison["binary_provenance"]["python"]["sha256"] == "c" * 64

    matlab["binary"]["source_commit"] = "d" * 40
    matlab_path.write_text(json.dumps(matlab), encoding="utf-8")
    with pytest.raises(ValueError, match="binary.source_commit"):
        MODULE.compare(python_path, matlab_path, rtol=1e-11, atol=1e-12)

    matlab = _report("matlab-mex")
    matlab["repeats"] = 0
    matlab_path.write_text(json.dumps(matlab), encoding="utf-8")
    with pytest.raises(ValueError, match="repeats must be a positive integer"):
        MODULE.compare(python_path, matlab_path, rtol=1e-11, atol=1e-12)


def test_comparison_rejects_nonfinite_values_and_handles_zero_tolerances(tmp_path):
    python = _report("python-pybind11")
    matlab = _report("matlab-mex")
    python["observables"] = {"zero": 0.0}
    matlab["observables"] = {"zero": 0.0}
    python_path = tmp_path / "python.json"
    matlab_path = tmp_path / "matlab.json"
    python_path.write_text(json.dumps(python), encoding="utf-8")
    matlab_path.write_text(json.dumps(matlab), encoding="utf-8")

    assert MODULE.compare(python_path, matlab_path, rtol=0.0, atol=0.0)["pass"]
    with pytest.raises(ValueError, match="rtol must be finite"):
        MODULE.compare(python_path, matlab_path, rtol=float("inf"), atol=0.0)

    matlab["observables"]["zero"] = 1.0
    matlab_path.write_text(json.dumps(matlab), encoding="utf-8")
    comparison = MODULE.compare(python_path, matlab_path, rtol=0.0, atol=0.0)
    assert not comparison["pass"]
    assert comparison["checks"]["zero"]["relative_error"] is None
    MODULE._json_text(comparison)

    matlab_path.write_text(json.dumps(matlab).replace("1.0", "NaN", 1), encoding="utf-8")
    with pytest.raises(ValueError, match="non-finite JSON constant"):
        MODULE.compare(python_path, matlab_path, rtol=0.0, atol=0.0)
    with pytest.raises(ValueError, match="Out of range float values"):
        MODULE._json_text({"forbidden": float("nan")})

    python = _report("python-pybind11")
    matlab = _report("matlab-mex")
    python["observables"] = {"huge": 1.0e308}
    matlab["observables"] = {"huge": -1.0e308}
    python_path.write_text(json.dumps(python), encoding="utf-8")
    matlab_path.write_text(json.dumps(matlab), encoding="utf-8")
    overflow_case = MODULE.compare(python_path, matlab_path, rtol=1.9, atol=0.0)
    assert not overflow_case["pass"]
    MODULE._json_text(overflow_case)


def test_benchmark_function_rejects_invalid_repeats_before_native_import():
    with pytest.raises(ValueError, match="positive integer"):
        MODULE.benchmark(Path("unused.json"), repeats=0)
    with pytest.raises(ValueError, match="positive integer"):
        MODULE.benchmark(Path("unused.json"), repeats=True)


def test_matlab_runner_hashes_mex_and_avoids_simulink_setup_mutation():
    source = (
        ROOT / "validation_test" / "accelerator" / "benchmark_beam_transfer_mex_matlab.m"
    ).read_text(encoding="utf-8")
    assert "radia.beam.propagateVariationalMap" in source
    assert '"sha256", sha256File(mexPath)' in source
    assert '"source_commit", sourceCommit' in source
    assert "readBuildProvenance(buildManifestPath, mexPath, repoRoot)" in source
    assert "ConfigureSimulinkFileGeneration=false" in source
    assert "radia_mex resolved outside the selected MATLAB source tree" in source
    assert "requireCleanSource(repoRoot)" in source
    assert "Benchmark timings must be finite and positive" in source
    assert "Benchmark observables must be finite" in source
    assert '"radia_legacy_core_version", radia.UtiVer()' in source
    assert '"path", "matlab/" + string(mexFile.name)' in source
    assert "permute(double(caseData.A_per_m), [2, 3, 1])" in source
    assert "item(2)+1, item(3)+1, item(4)+1, item(1)+1" in source
    assert "item(2)+1, item(3)+1, item(4)+1, item(5)+1, item(1)+1" in source


def test_current_native_api_retains_canonical_h5_and_f4_outputs():
    python_api = (ROOT / "src" / "radia" / "beam.py").read_text(encoding="utf-8")
    matlab_api = (
        ROOT / "matlab" / "+radia" / "+beam" / "canonicalBodyHamiltonianJet.m"
    ).read_text(encoding="utf-8")
    mex_api = (ROOT / "src" / "matlab" / "radia_beam_mex_commands.cpp").read_text(
        encoding="utf-8"
    )

    assert "H2/H3/H4/H5" in python_api
    assert "A/F2/F3/F4" in python_api
    assert "H2/H3/H4/H5" in matlab_api
    assert "A/F2/F3/F4" in matlab_api
    assert '"H5_per_m"' in mex_api
    assert '"F4_per_m"' in mex_api


def test_build_writes_hash_bound_native_provenance_manifests():
    source = (ROOT / "Build.ps1").read_text(encoding="utf-8")
    helper = (ROOT / "tools" / "native_build_provenance.ps1").read_text(
        encoding="utf-8"
    )
    assert "Get-NativeBuildSourceIdentity" in source
    assert "Clear-NativeBuildProvenance" in source
    assert "$srcHash -eq $dstHash" in source
    assert "radia.native-build-provenance.v1" in helper
    assert "binary_sha256" in helper
    assert "source_dirty" in helper
    assert "source_change_fingerprint_sha256" in helper
