from __future__ import annotations

import math

from radia_mcp.radia_ngsolve.femm_v46_identity import validate_public_identity


FORCE = "v46_public_axisymmetric_force_weighted_stress_contour_open_boundary_partial_mismatch"
BH = "v46_public_nonlinear_bh_curve_unit_scale_convergence_status_nan_mismatch"


def _identity() -> dict[str, object]:
    generation = "test-force-open-v46"
    force = {
        "generation": generation,
        **{key: generation for key in ("force_generation", "contour_generation", "open_boundary_generation", "solve_generation", "result_generation")},
        "force_method": "weighted_stress_tensor", "result_force_method": "weighted_stress_tensor",
        "open_boundary": "outer_kelvin", "result_open_boundary": "outer_kelvin",
        "solve_state": "complete", "result_solve_state": "complete",
        "partial_solve_status": "none", "result_partial_solve_status": "none",
        "force_n": [12.0, -3.0], "result_force_n": [12.0, -3.0],
        "axisymmetric_factor": 2.0 * math.pi, "result_axisymmetric_factor": 2.0 * math.pi,
        "contour_owner": "contour:test-force-open-v46", "result_contour_owner": "contour:test-force-open-v46",
        "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64,
    }
    generation = "test-bh-unit-v46"
    bh = {
        "generation": generation,
        **{key: generation for key in ("bh_curve_generation", "unit_generation", "convergence_generation", "result_generation")},
        "bh_unit": "tesla_ampere_per_meter", "result_bh_unit": "tesla_ampere_per_meter",
        "convergence_status": "converged", "result_convergence_status": "converged",
        "nonfinite_status": "none", "result_nonfinite_status": "none",
        "bh_points": [[0.0, 0.0], [0.8, 500.0]], "result_bh_points": [[0.0, 0.0], [0.8, 500.0]],
        "material_owner": "material:test-bh-unit-v46", "result_material_owner": "material:test-bh-unit-v46",
        "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
    }
    return {FORCE: force, BH: bh}


def test_v46_public_femm_identity_accepts_closed_artifacts():
    checks = validate_public_identity(_identity())
    assert checks and all(checks.values())


def test_v46_public_femm_identity_rejects_force_partial_nan_mutation():
    identity = _identity()
    identity[FORCE]["result_partial_solve_status"] = "partial"
    identity[FORCE]["result_force_n"] = [float("nan"), -3.0]
    checks = validate_public_identity(identity)
    assert checks and not all(checks.values())


def test_v46_public_femm_identity_rejects_bh_unit_convergence_mutation():
    identity = _identity()
    identity[BH]["result_bh_unit"] = "gauss_oersted"
    identity[BH]["result_convergence_status"] = "not_converged"
    checks = validate_public_identity(identity)
    assert checks and not all(checks.values())


def test_v46_case_ids_are_frozen():
    assert FORCE.startswith("v46_public_")
    assert BH.startswith("v46_public_")
