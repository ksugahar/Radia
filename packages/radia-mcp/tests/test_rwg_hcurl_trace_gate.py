import copy
import json

import pytest

from radia_mcp.radia_ngsolve.rwg_hcurl_trace_gate import (
    rwg_hcurl_trace_consistency_gate as build_gate,
)
from radia_mcp.radia_ngsolve.server import rwg_hcurl_trace_consistency_gate


def summary() -> dict:
    cases = [{"id": f"GYP-{index:03d}", "passed": True} for index in range(81, 91)]
    return {
        "cases": cases,
        "trace_diagnostics": {
            "rwg_dofs": 6,
            "hcurl_dofs": 6,
            "rwg_to_hcurl_edge_ids": [1, 2, 3, 4, 5, 6],
            "edge_ids_unique": True,
            "selector_nonzeros": 6,
            "selector_orthogonality_error": 0.0,
            "rotated_trace_nonzeros": 6,
            "rotated_trace_rank": 6,
            "rotated_trace_magnitude_error": 0.0,
            "rwg_gram_symmetry_error": 4.0e-17,
            "rwg_gram_min_eigenvalue": 0.5,
            "rwg_quadrature_order": 3,
            "hcurl_mass_symmetry_error": 0.0,
            "hcurl_curlcurl_symmetry_error": 0.0,
            "curl_grad_error": 0.0,
            "ngsolve_hcurl_mass_error": 8.0e-16,
            "ngsolve_hcurl_curlcurl_error": 0.0,
        },
        "capabilities": {
            "ok": True,
            "ngsolve_version": "6.2.2604",
            "hcurl_dofs": 6,
            "has_maxwell_sl": False,
        },
        "maxwell_bem_claimed": False,
        "matlab_tests": {"test_count": 18, "passed": 18, "failed": 0, "incomplete": 0},
    }


def test_accepts_rwg_hcurl_trace_and_de_rham_evidence():
    result = build_gate(summary())
    assert result["status"] == "ok"
    assert result["solver_ready"] is True
    assert json.loads(rwg_hcurl_trace_consistency_gate(json.dumps(summary())))["status"] == "ok"


def test_rejects_edge_alias_and_curl_grad_drift():
    bad = summary()
    bad["trace_diagnostics"]["rwg_to_hcurl_edge_ids"][-1] = 5
    bad["trace_diagnostics"]["edge_ids_unique"] = False
    bad["trace_diagnostics"]["curl_grad_error"] = 1.0e-3
    result = build_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["rwg_edges_map_one_to_one_into_hcurl"] is False
    assert result["checks"]["discrete_curl_grad_identity_closes"] is False


def test_rejects_false_maxwell_bem_claim_and_reference_drift():
    bad = copy.deepcopy(summary())
    bad["maxwell_bem_claimed"] = True
    bad["trace_diagnostics"]["ngsolve_hcurl_mass_error"] = 1.0e-4
    bad["matlab_tests"]["failed"] = 1
    bad["matlab_tests"]["passed"] = 17
    result = build_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["maxwell_bem_capability_is_not_overclaimed"] is False
    assert result["checks"]["independent_hcurl_reference_matrices_match"] is False
    assert result["checks"]["related_matlab_tests_pass"] is False


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("rwg_dofs", 6.5),
        ("selector_nonzeros", 6.5),
        ("selector_orthogonality_error", -1.0e-16),
    ],
)
def test_rejects_nonintegral_inventory_and_negative_error_norms(field, value):
    bad = summary()
    bad["trace_diagnostics"][field] = value
    with pytest.raises(ValueError):
        build_gate(bad)


def test_rejects_tolerance_relaxation_and_missing_reference_version():
    relaxed = summary()
    relaxed["tolerances"] = {"maximum_ngsolve_matrix_error": 1.0}
    with pytest.raises(ValueError, match="policy limits"):
        build_gate(relaxed)

    missing_version = summary()
    missing_version["capabilities"]["ngsolve_version"] = ""
    result = build_gate(missing_version)
    assert result["checks"]["independent_hcurl_reference_matrices_match"] is False
