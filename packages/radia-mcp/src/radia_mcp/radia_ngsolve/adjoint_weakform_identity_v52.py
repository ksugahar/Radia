"""Adjoint-sensitivity and weak-form identity checks for v52 artifacts."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

from .coupled_periodic_identity_v53 import validate_public_v53_identity


ADJOINT = "adjoint_objective_scaling_conjugation_design_owner_identity"
WEAK_FORM = "weakform_testfunction_sign_boundary_orientation_measure_owner_identity"


def _digest(value: object) -> bool:
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generation(row: Mapping[str, object], *names: str) -> bool:
    value = str(row.get("generation") or "")
    return bool(value) and all(row.get(name) == value for name in names)


def _result_identity(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _finite_list(value: object, *, length: int | None = None) -> bool:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        return False
    values = list(value)
    return (length is None or len(values) == length) and bool(values) and all(
        isinstance(item, (int, float))
        and not isinstance(item, bool)
        and math.isfinite(float(item))
        for item in values
    )


def _adjoint_ok(row: Mapping[str, object]) -> bool:
    design = row.get("design_variable_order")
    gradient = row.get("scaled_gradient")
    return (
        _generation(
            row,
            "objective_generation", "scaling_generation", "conjugation_generation",
            "design_generation", "gradient_generation", "owner_generation", "result_generation",
        )
        and str(row.get("objective_tag") or "").startswith("obj_")
        and row.get("result_objective_tag") == row.get("objective_tag")
        and isinstance(row.get("objective_scale"), (int, float))
        and not isinstance(row.get("objective_scale"), bool)
        and math.isfinite(float(row["objective_scale"]))
        and float(row["objective_scale"]) > 0.0
        and row.get("result_objective_scale") == row.get("objective_scale")
        and row.get("complex_adjoint_convention") == "hermitian_conjugate"
        and row.get("result_complex_adjoint_convention") == row.get("complex_adjoint_convention")
        and isinstance(design, Sequence)
        and not isinstance(design, (str, bytes))
        and bool(design)
        and len(design) == len(set(design))
        and all(isinstance(name, str) and name for name in design)
        and row.get("result_design_variable_order") == design
        and _finite_list(gradient, length=len(design))
        and row.get("result_scaled_gradient") == gradient
        and str(row.get("solution_owner") or "").startswith("solution:")
        and row.get("result_solution_owner") == row.get("solution_owner")
        and _result_identity(row)
    )


def _weak_form_ok(row: Mapping[str, object]) -> bool:
    terms = row.get("weak_terms")
    return (
        _generation(
            row,
            "testfunction_generation", "sign_generation", "orientation_generation",
            "measure_generation", "term_generation", "owner_generation", "result_generation",
        )
        and str(row.get("test_function") or "").startswith("test(")
        and row.get("result_test_function") == row.get("test_function")
        and row.get("residual_sign") == "lhs_minus_rhs"
        and row.get("result_residual_sign") == row.get("residual_sign")
        and row.get("boundary_orientation") == "outward_normal"
        and row.get("result_boundary_orientation") == row.get("boundary_orientation")
        and row.get("integration_measure") == "surface_jacobian"
        and row.get("result_integration_measure") == row.get("integration_measure")
        and isinstance(terms, Sequence)
        and not isinstance(terms, (str, bytes))
        and len(terms) >= 2
        and all(isinstance(term, str) and term for term in terms)
        and row.get("result_weak_terms") == terms
        and str(row.get("form_owner") or "").startswith("weak-form:")
        and row.get("result_form_owner") == row.get("form_owner")
        and _result_identity(row)
    )


def validate_public_v52_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping):
        return {}
    checks: dict[str, bool] = {}
    v53 = validate_public_v53_identity(payload)
    if v53:
        checks.update(v53["checks"])
    adjoint = payload.get(ADJOINT)
    weak_form = payload.get(WEAK_FORM)
    if adjoint is not None:
        checks["v52_adjoint_objective_scaling_conjugation_design_owner"] = isinstance(adjoint, Mapping) and _adjoint_ok(adjoint)
    if weak_form is not None:
        checks["v52_weakform_testfunction_sign_orientation_measure_owner"] = isinstance(weak_form, Mapping) and _weak_form_ok(weak_form)
    if not checks:
        return {}
    return {
        "policy": "adjoint_weakform_identity_v52",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, accepted in checks.items() if not accepted],
    }
