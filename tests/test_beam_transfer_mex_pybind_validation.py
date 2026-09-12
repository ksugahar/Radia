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
        },
        "median_s": 0.25,
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


def test_pybind_lane_fails_loudly_without_matching_source_artifact(tmp_path):
    with pytest.raises(RuntimeError, match="installed-wheel fallback is forbidden"):
        MODULE._load_source_backend(tmp_path)


def test_mocked_provenance_and_comparison_contract(tmp_path):
    python = _report("python-pybind11")
    matlab = _report("matlab-mex")
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


def test_matlab_runner_hashes_mex_and_avoids_simulink_setup_mutation():
    source = (
        ROOT / "validation_test" / "accelerator" / "benchmark_beam_transfer_mex_matlab.m"
    ).read_text(encoding="utf-8")
    assert "radia.beam.propagateVariationalMap" in source
    assert '"sha256", sha256File(mexPath)' in source
    assert '"source_commit", sourceCommit' in source
    assert "ConfigureSimulinkFileGeneration=false" in source
    assert "radia_mex resolved outside the selected MATLAB source tree" in source
    assert "requireCleanSource(repoRoot)" in source
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
