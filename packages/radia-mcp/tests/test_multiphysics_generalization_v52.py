from copy import deepcopy

from radia_mcp.radia_ngsolve.adjoint_weakform_identity_v52 import ADJOINT, WEAK_FORM, validate_public_v52_identity


CASE_IDS = {
    "v52_public_adjoint_sensitivity_objective_scaling_complex_conjugation_design_owner_mismatch",
    "v52_public_weakform_testfunction_sign_boundary_orientation_measure_owner_mismatch",
}


def _generation(prefix: str, names: tuple[str, ...]) -> dict[str, str]:
    return {"generation": prefix, **{name: prefix for name in names}}


def _records():
    adjoint = {
        **_generation("adj-v52-test", ("objective_generation", "scaling_generation", "conjugation_generation", "design_generation", "gradient_generation", "owner_generation", "result_generation")),
        "objective_tag": "obj_loss", "result_objective_tag": "obj_loss",
        "objective_scale": 0.01, "result_objective_scale": 0.01,
        "complex_adjoint_convention": "hermitian_conjugate", "result_complex_adjoint_convention": "hermitian_conjugate",
        "design_variable_order": ["x", "y"], "result_design_variable_order": ["x", "y"],
        "scaled_gradient": [1.0, -2.0], "result_scaled_gradient": [1.0, -2.0],
        "solution_owner": "solution:adj-v52-test", "result_solution_owner": "solution:adj-v52-test",
        "result_sha256": "a" * 64, "accepted_result_sha256": "a" * 64,
    }
    weak = {
        **_generation("weak-v52-test", ("testfunction_generation", "sign_generation", "orientation_generation", "measure_generation", "term_generation", "owner_generation", "result_generation")),
        "test_function": "test(u)", "result_test_function": "test(u)",
        "residual_sign": "lhs_minus_rhs", "result_residual_sign": "lhs_minus_rhs",
        "boundary_orientation": "outward_normal", "result_boundary_orientation": "outward_normal",
        "integration_measure": "surface_jacobian", "result_integration_measure": "surface_jacobian",
        "weak_terms": ["test(u)*u", "dot(grad(test(u)),grad(u))"],
        "result_weak_terms": ["test(u)*u", "dot(grad(test(u)),grad(u))"],
        "form_owner": "weak-form:v52-test", "result_form_owner": "weak-form:v52-test",
        "result_sha256": "b" * 64, "accepted_result_sha256": "b" * 64,
    }
    return {ADJOINT: adjoint, WEAK_FORM: weak}


def test_v52_public_positive_replay_is_accepted():
    assert validate_public_v52_identity(_records())["status"] == "ok"


def test_v52_public_mixed_adjoint_identity_is_rejected():
    value = deepcopy(_records())
    value[ADJOINT]["result_complex_adjoint_convention"] = "transpose_without_conjugation"
    assert validate_public_v52_identity(value)["status"] == "needs_attention"


def test_v52_public_mixed_weakform_identity_is_rejected():
    value = deepcopy(_records())
    value[WEAK_FORM]["result_boundary_orientation"] = "inward_normal"
    assert validate_public_v52_identity(value)["status"] == "needs_attention"
