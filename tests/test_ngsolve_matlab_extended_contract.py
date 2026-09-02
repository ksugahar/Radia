"""Fast artifact checks for the manual NGSolve/MATLAB validation tiers."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULT_DIR = ROOT / "validation_test" / "ngsolve_matlab_parity"


def _load(name: str) -> dict:
    return json.loads((RESULT_DIR / name).read_text(encoding="utf-8"))


def test_breadth_result_covers_500_small_2d_and_3d_cases() -> None:
    result = _load("results_ngsolve_matlab_breadth_500.json")
    assert result["schema"] == "radia.ngsolve-python-matlab-breadth-500.v1"
    assert result["all_passed"] is True
    assert result["case_count"] == result["passed_count"] == 500
    assert result["failed_count"] == 0
    assert result["mesh_count"] == 10
    assert result["two_d_case_count"] == result["three_d_case_count"] == 250
    assert result["complex_case_count"] == 60
    assert result["spatial_coefficient_case_count"] == 190
    assert result["boundary_case_count"] == 40
    assert result["maximum_dofs"] < result["dof_limit"] == 1_000_000

    cases = result["case_results"]
    assert len({case["case_id"] for case in cases}) == 500
    assert {case["space"] for case in cases} == {"h1", "hcurl", "hdiv"}
    assert {case["form"] for case in cases} == {
        "boundary_mass",
        "curlcurl",
        "divdiv",
        "mass",
        "stiffness",
    }


def test_scale_result_covers_20_matrix_free_cases_above_10k_dofs() -> None:
    result = _load("results_ngsolve_matlab_scale_20.json")
    assert result["schema"] == "radia.ngsolve-python-matlab-scale-20.v1"
    assert result["all_passed"] is True
    assert result["case_count"] == result["passed_count"] == 20
    assert result["failed_count"] == 0
    assert result["minimum_dofs"] >= result["dof_floor"] == 10_000
    assert result["maximum_dofs"] < result["dof_limit"] == 1_000_000
    cases = result["case_results"]
    assert len({case["mesh_id"] for case in cases}) == 4
    assert len({case["case_id"] for case in cases}) == 20


def test_manufactured_result_covers_five_element_families() -> None:
    result = _load("results_ngsolve_matlab_manufactured.json")
    assert result["schema"] == "radia.ngsolve-matlab-manufactured-solution.v1"
    assert result["all_passed"] is True
    assert result["case_count"] == result["passed_count"] == 15
    assert result["failed_count"] == 0
    assert result["mesh_family_count"] == 5
    assert len({case["mesh_id"] for case in result["case_results"]}) == 5


def test_hibino_result_records_the_exact_native_runtime() -> None:
    result = _load("results_ngsolve_matlab_extended_hibino.json")
    assert result["schema"] == (
        "radia.ngsolve-python-matlab-extended-validation.v1"
    )
    assert result["all_passed"] is True
    assert result["host_name"].upper() == "HIBINO"
    assert result["matlab_release"] == "2026a"
    assert result["ngsolve_version"] == "6.2.2606"
    assert len(result["mex_sha256"]) == 64
    int(result["mex_sha256"], 16)
    assert result["native_handle_count_after"] == 0


def test_canonical_tier_results_are_the_recorded_hibino_results() -> None:
    for stem in (
        "results_ngsolve_matlab_breadth_500",
        "results_ngsolve_matlab_scale_20",
        "results_ngsolve_matlab_manufactured",
    ):
        assert _load(f"{stem}.json") == _load(f"{stem}_hibino.json")
