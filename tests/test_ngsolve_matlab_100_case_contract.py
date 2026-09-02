"""Fast catalog checks for the manual 100-case NGSolve/MATLAB validation."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CATALOG_PATH = (
    ROOT / "validation_test" / "ngsolve_matlab_parity" / "case_catalog.py"
)
RESULT_PATH = CATALOG_PATH.with_name(
    "results_ngsolve_matlab_100_case_parity.json"
)


def _catalog_module():
    spec = importlib.util.spec_from_file_location("ngsolve_matlab_case_catalog", CATALOG_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_catalog_has_100_bounded_cases_with_2d_and_3d() -> None:
    module = _catalog_module()
    cases = module.build_case_catalog()
    assert len(cases) == 100
    assert len({case["case_id"] for case in cases}) == 100
    assert module.DOF_LIMIT == 1_000_000
    assert sum(case["dimension"] == 2 for case in cases) == 50
    assert sum(case["dimension"] == 3 for case in cases) == 50
    assert sum(bool(case["solve"]) for case in cases) == 72


def test_catalog_covers_linear_ngsolve_space_and_operator_families() -> None:
    cases = _catalog_module().build_case_catalog()
    combinations = {(case["space"], case["form"]) for case in cases}
    assert combinations == {
        ("h1", "mass"),
        ("h1", "stiffness"),
        ("hcurl", "mass"),
        ("hcurl", "curlcurl"),
        ("hdiv", "mass"),
        ("hdiv", "divdiv"),
    }
    assert max(case["order"] for case in cases if case["space"] == "h1") == 4
    assert all(case["weight"] > 0.0 for case in cases)
    assert all(
        case["dirichlet"] == ".*"
        for case in cases
        if case["space"] == "h1" and case["form"] == "stiffness"
    )


def test_committed_validation_result_records_100_passes() -> None:
    result = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    assert result["schema"] == (
        "radia.ngsolve-python-matlab-100-case-validation.v1"
    )
    assert result["all_passed"] is True
    assert result["case_count"] == result["passed_count"] == 100
    assert result["two_d_case_count"] == result["three_d_case_count"] == 50
    assert result["maximum_dofs"] < result["dof_limit"] == 1_000_000
    assert result["native_handle_count_after"] == 0
