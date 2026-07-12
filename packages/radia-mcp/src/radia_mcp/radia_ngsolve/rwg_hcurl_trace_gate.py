"""Solver-neutral RWG-to-HCurl trace and discrete de Rham gate."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any


_POLICY_LIMITS = {
    "topology": 1.0e-14,
    "trace": 1.0e-12,
    "de_rham": 1.0e-12,
    "reference": 1.0e-10,
}


def _finite(value: object, name: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(parsed):
        raise ValueError(f"{name} must be finite")
    return parsed


def _nonnegative(value: object, name: str) -> float:
    parsed = _finite(value, name)
    if parsed < 0.0:
        raise ValueError(f"{name} must be nonnegative")
    return parsed


def _integer(value: object, name: str) -> int:
    parsed = _finite(value, name)
    if not parsed.is_integer():
        raise ValueError(f"{name} must be an integer")
    return int(parsed)


def rwg_hcurl_trace_consistency_gate(
    summary: Mapping[str, object],
) -> dict[str, Any]:
    """Gate first-order RWG/HCurl topology, trace, Gram, and reference matrices."""

    if not isinstance(summary, Mapping):
        raise ValueError("summary must be an object")
    cases_value = summary.get("cases")
    if not isinstance(cases_value, Sequence) or isinstance(cases_value, (str, bytes)):
        raise ValueError("cases must be an array")
    cases = list(cases_value)
    if any(not isinstance(row, Mapping) for row in cases):
        raise ValueError("cases must contain objects")
    trace = summary.get("trace_diagnostics")
    capabilities = summary.get("capabilities")
    tests = summary.get("matlab_tests")
    if not isinstance(trace, Mapping):
        raise ValueError("trace_diagnostics must be an object")
    if not isinstance(capabilities, Mapping):
        raise ValueError("capabilities must be an object")
    if not isinstance(tests, Mapping):
        raise ValueError("matlab_tests must be an object")
    tolerances = summary.get("tolerances", {})
    if not isinstance(tolerances, Mapping):
        raise ValueError("tolerances must be an object")
    limits = {
        "topology": _nonnegative(
            tolerances.get("maximum_topology_error", _POLICY_LIMITS["topology"]),
            "maximum_topology_error",
        ),
        "trace": _nonnegative(
            tolerances.get("maximum_trace_error", _POLICY_LIMITS["trace"]),
            "maximum_trace_error",
        ),
        "de_rham": _nonnegative(
            tolerances.get("maximum_curl_grad_error", _POLICY_LIMITS["de_rham"]),
            "maximum_curl_grad_error",
        ),
        "reference": _nonnegative(
            tolerances.get("maximum_ngsolve_matrix_error", _POLICY_LIMITS["reference"]),
            "maximum_ngsolve_matrix_error",
        ),
    }
    if any(limits[name] > _POLICY_LIMITS[name] for name in limits):
        raise ValueError("tolerances cannot exceed the policy limits")

    case_ids = [str(row.get("id") or "") for row in cases]
    expected_ids = [f"GYP-{index:03d}" for index in range(81, 91)]
    edge_ids_value = trace.get("rwg_to_hcurl_edge_ids")
    edge_ids = (
        [_integer(value, "rwg_to_hcurl_edge_ids") for value in edge_ids_value]
        if isinstance(edge_ids_value, Sequence) and not isinstance(edge_ids_value, (str, bytes))
        else []
    )
    rwg_dofs = _integer(trace.get("rwg_dofs"), "rwg_dofs")
    hcurl_dofs = _integer(trace.get("hcurl_dofs"), "hcurl_dofs")
    selector_error = _nonnegative(
        trace.get("selector_orthogonality_error"), "selector_orthogonality_error"
    )
    rotated_error = _nonnegative(
        trace.get("rotated_trace_magnitude_error"), "rotated_trace_magnitude_error"
    )
    gram_symmetry = _nonnegative(
        trace.get("rwg_gram_symmetry_error"), "rwg_gram_symmetry_error"
    )
    gram_minimum = _finite(trace.get("rwg_gram_min_eigenvalue"), "rwg_gram_min_eigenvalue")
    mass_symmetry = _nonnegative(
        trace.get("hcurl_mass_symmetry_error"), "hcurl_mass_symmetry_error"
    )
    curlcurl_symmetry = _nonnegative(
        trace.get("hcurl_curlcurl_symmetry_error"), "hcurl_curlcurl_symmetry_error"
    )
    curl_grad = _nonnegative(trace.get("curl_grad_error"), "curl_grad_error")
    ngsolve_mass = _nonnegative(
        trace.get("ngsolve_hcurl_mass_error"), "ngsolve_hcurl_mass_error"
    )
    ngsolve_curlcurl = _nonnegative(
        trace.get("ngsolve_hcurl_curlcurl_error"), "ngsolve_hcurl_curlcurl_error"
    )
    selector_nonzeros = _integer(trace.get("selector_nonzeros"), "selector_nonzeros")
    rotated_nonzeros = _integer(
        trace.get("rotated_trace_nonzeros"), "rotated_trace_nonzeros"
    )
    rotated_rank = _integer(trace.get("rotated_trace_rank"), "rotated_trace_rank")
    quadrature_order = _integer(trace.get("rwg_quadrature_order"), "rwg_quadrature_order")
    test_count = _integer(tests.get("test_count"), "matlab_tests.test_count")
    passed_count = _integer(tests.get("passed"), "matlab_tests.passed")
    failed_count = _integer(tests.get("failed"), "matlab_tests.failed")
    incomplete_count = _integer(tests.get("incomplete"), "matlab_tests.incomplete")
    checks = {
        "ten_declared_source_cases_pass": case_ids == expected_ids
        and all(row.get("passed") is True for row in cases),
        "unit_tetra_rwg_hcurl_dof_counts_match": rwg_dofs == 6
        and hcurl_dofs == 6
        and _integer(capabilities.get("hcurl_dofs", -1), "capabilities.hcurl_dofs") == 6,
        "rwg_edges_map_one_to_one_into_hcurl": edge_ids == list(range(1, 7))
        and trace.get("edge_ids_unique") is True
        and selector_nonzeros == rwg_dofs
        and selector_error <= limits["topology"],
        "rotated_trace_is_full_rank_with_inverse_edge_scaling": rotated_nonzeros == rwg_dofs
        and rotated_rank == rwg_dofs
        and rotated_error <= limits["trace"],
        "rwg_gram_is_symmetric_positive_definite": gram_symmetry <= limits["trace"]
        and gram_minimum > 0.0
        and quadrature_order >= 3,
        "hcurl_mass_and_curlcurl_are_symmetric": max(mass_symmetry, curlcurl_symmetry)
        <= limits["trace"],
        "discrete_curl_grad_identity_closes": curl_grad <= limits["de_rham"],
        "independent_hcurl_reference_matrices_match": capabilities.get("ok") is True
        and bool(str(capabilities.get("ngsolve_version") or "").strip())
        and max(ngsolve_mass, ngsolve_curlcurl) <= limits["reference"],
        "maxwell_bem_capability_is_not_overclaimed": capabilities.get("has_maxwell_sl") is True
        or summary.get("maxwell_bem_claimed") is False,
        "related_matlab_tests_pass": passed_count == test_count
        and test_count >= 1
        and failed_count == 0
        and incomplete_count == 0,
    }
    accepted = all(checks.values())
    return {
        "policy": "rwg_hcurl_trace_consistency_gate_v1",
        "status": "ok" if accepted else "needs_attention",
        "solver_ready": accepted,
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "limits": limits,
        "metrics": {
            "case_count": len(cases),
            "rwg_dofs": rwg_dofs,
            "hcurl_dofs": hcurl_dofs,
            "selector_orthogonality_error": selector_error,
            "rotated_trace_magnitude_error": rotated_error,
            "rwg_gram_min_eigenvalue": gram_minimum,
            "curl_grad_error": curl_grad,
            "ngsolve_hcurl_mass_error": ngsolve_mass,
            "ngsolve_hcurl_curlcurl_error": ngsolve_curlcurl,
            "related_test_count": test_count,
        },
        "lesson": (
            "A readable first-order Maxwell trace must preserve one-to-one boundary-edge identity, "
            "the rotated Nedelec-to-RWG edge-length scaling, a positive RWG Gram matrix, the discrete "
            "curl-grad identity, and an independent HCurl matrix reference without claiming unavailable Maxwell BEM."
        ),
    }
