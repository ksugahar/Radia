"""Fail-close structural geometry checks for Radia validation summaries."""

from __future__ import annotations

import json
import math
from collections.abc import Mapping

SCHEMA = "radia-visual-geometry-validation/v1"
NORMALIZATION = "dimensionless_relative_residuals"

_PROFILES = {
    "de_rham": ("de_rham", "hodge"),
    "mapped_em": ("de_rham", "pullback", "hodge", "maxwell"),
    "surface": ("hodge", "connection", "surface", "holonomy"),
    "full": (
        "de_rham",
        "pullback",
        "hodge",
        "connection",
        "surface",
        "holonomy",
        "maxwell",
    ),
}


def _number(value: object, name: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(parsed):
        raise ValueError(f"{name} must be finite")
    return parsed


def _nonnegative(value: object, name: str) -> float:
    parsed = _number(value, name)
    if parsed < 0.0:
        raise ValueError(f"{name} must be nonnegative")
    return parsed


def _positive(value: object, name: str) -> float:
    parsed = _number(value, name)
    if parsed <= 0.0:
        raise ValueError(f"{name} must be positive")
    return parsed


def _section(summary: Mapping[str, object], name: str) -> Mapping[str, object]:
    value = summary.get(name)
    if not isinstance(value, Mapping):
        raise TypeError(f"{name} must be an object for the selected profile")
    return value


def _tolerance(tolerances: Mapping[str, object], name: str) -> float:
    return _nonnegative(tolerances.get(name), f"tolerances.{name}")


def _wrapped_angle_error(left: float, right: float) -> float:
    return abs((left - right + math.pi) % (2.0 * math.pi) - math.pi)


def evaluate_visual_geometry(summary: Mapping[str, object]) -> dict:
    """Evaluate one visual-geometry validation profile."""

    if not isinstance(summary, Mapping):
        raise TypeError("summary must be an object")
    if summary.get("schema") != SCHEMA:
        raise ValueError(f"schema must be '{SCHEMA}'")
    if summary.get("normalization_policy") != NORMALIZATION:
        raise ValueError(f"normalization_policy must be '{NORMALIZATION}'")

    profile = str(summary.get("profile", "")).strip().lower()
    if profile not in _PROFILES:
        raise ValueError(
            "profile must be one of: " + ", ".join(sorted(_PROFILES))
        )
    tolerances = summary.get("tolerances")
    if not isinstance(tolerances, Mapping):
        raise TypeError("tolerances must be an object")

    required = _PROFILES[profile]
    for name in required:
        _section(summary, name)

    checks: dict[str, bool] = {}
    metrics: dict[str, float | int | str] = {}

    if "de_rham" in required:
        section = _section(summary, "de_rham")
        curl_grad = _nonnegative(
            section.get("curl_grad_relative_residual"),
            "de_rham.curl_grad_relative_residual",
        )
        div_curl = _nonnegative(
            section.get("div_curl_relative_residual"),
            "de_rham.div_curl_relative_residual",
        )
        limit = _tolerance(tolerances, "de_rham_relative")
        checks["curl_grad_is_zero"] = curl_grad <= limit
        checks["div_curl_is_zero"] = div_curl <= limit
        metrics["curl_grad_relative_residual"] = curl_grad
        metrics["div_curl_relative_residual"] = div_curl

    if "pullback" in required:
        section = _section(summary, "pullback")
        error = _nonnegative(
            section.get("exterior_derivative_commutator_relative_residual"),
            "pullback.exterior_derivative_commutator_relative_residual",
        )
        limit = _tolerance(tolerances, "pullback_relative")
        checks["pullback_commutes_with_exterior_derivative"] = error <= limit
        metrics["pullback_commutator_relative_residual"] = error

    if "hodge" in required:
        section = _section(summary, "hodge")
        symmetry = _nonnegative(
            section.get("symmetry_relative_error"),
            "hodge.symmetry_relative_error",
        )
        minimum_eigenvalue = _number(
            section.get("minimum_eigenvalue"), "hodge.minimum_eigenvalue"
        )
        symmetry_limit = _tolerance(tolerances, "hodge_symmetry_relative")
        positive_minimum = _positive(
            tolerances.get("hodge_positive_min"),
            "tolerances.hodge_positive_min",
        )
        checks["hodge_is_symmetric"] = symmetry <= symmetry_limit
        checks["hodge_is_positive"] = minimum_eigenvalue >= positive_minimum
        metrics["hodge_symmetry_relative_error"] = symmetry
        metrics["hodge_minimum_eigenvalue"] = minimum_eigenvalue

    if "connection" in required:
        convention = summary.get("cartan_sign_convention")
        if not isinstance(convention, str) or not convention.strip():
            raise ValueError(
                "cartan_sign_convention must declare the implemented convention"
            )
        section = _section(summary, "connection")
        metric_error = _nonnegative(
            section.get("metric_compatibility_relative_residual"),
            "connection.metric_compatibility_relative_residual",
        )
        torsion_error = _nonnegative(
            section.get("torsion_relative_residual"),
            "connection.torsion_relative_residual",
        )
        curvature_error = _nonnegative(
            section.get("curvature_form_relative_residual"),
            "connection.curvature_form_relative_residual",
        )
        bianchi_error = _nonnegative(
            section.get("bianchi_relative_residual"),
            "connection.bianchi_relative_residual",
        )
        limit = _tolerance(tolerances, "connection_relative")
        checks["connection_is_metric_compatible"] = metric_error <= limit
        checks["cartan_first_structure_equation_is_torsion_free"] = (
            torsion_error <= limit
        )
        checks["cartan_second_structure_equation_holds"] = curvature_error <= limit
        checks["bianchi_identity_holds"] = bianchi_error <= limit
        metrics["connection_metric_relative_residual"] = metric_error
        metrics["cartan_torsion_relative_residual"] = torsion_error
        metrics["cartan_curvature_relative_residual"] = curvature_error
        metrics["bianchi_relative_residual"] = bianchi_error
        metrics["cartan_sign_convention"] = convention.strip()

    surface_curvature = None
    if "surface" in required:
        orientation = summary.get("orientation_convention")
        if not isinstance(orientation, str) or not orientation.strip():
            raise ValueError(
                "orientation_convention must be declared for surface profiles"
            )
        section = _section(summary, "surface")
        surface_curvature = _number(
            section.get("curvature_integral_rad"),
            "surface.curvature_integral_rad",
        )
        boundary = _number(
            section.get("boundary_geodesic_curvature_rad", 0.0),
            "surface.boundary_geodesic_curvature_rad",
        )
        corners = _number(
            section.get("corner_turning_rad", 0.0),
            "surface.corner_turning_rad",
        )
        euler = _number(
            section.get("euler_characteristic"),
            "surface.euler_characteristic",
        )
        if abs(euler - round(euler)) > 1.0e-12:
            raise ValueError("surface.euler_characteristic must be an integer")
        residual = abs(
            surface_curvature + boundary + corners - 2.0 * math.pi * euler
        )
        limit = _tolerance(tolerances, "gauss_bonnet_rad")
        checks["gauss_bonnet_closes"] = residual <= limit
        metrics["gauss_bonnet_residual_rad"] = residual
        metrics["euler_characteristic"] = round(euler)
        metrics["orientation_convention"] = orientation.strip()

    if "holonomy" in required:
        section = _section(summary, "holonomy")
        measured = _number(section.get("angle_rad"), "holonomy.angle_rad")
        curvature = _number(
            section.get("curvature_integral_rad", surface_curvature),
            "holonomy.curvature_integral_rad",
        )
        orientation_sign = _number(
            section.get("orientation_sign"), "holonomy.orientation_sign"
        )
        if orientation_sign not in (-1.0, 1.0):
            raise ValueError("holonomy.orientation_sign must be +1 or -1")
        error = _wrapped_angle_error(measured, orientation_sign * curvature)
        limit = _tolerance(tolerances, "holonomy_rad")
        checks["holonomy_matches_integrated_curvature"] = error <= limit
        metrics["holonomy_residual_rad"] = error

    if "maxwell" in required:
        section = _section(summary, "maxwell")
        error = _nonnegative(
            section.get("dF_relative_residual"),
            "maxwell.dF_relative_residual",
        )
        limit = _tolerance(tolerances, "maxwell_relative")
        checks["homogeneous_maxwell_form_is_closed"] = error <= limit
        metrics["dF_relative_residual"] = error

    issues = [name for name, ok in checks.items() if not ok]
    return {
        "schema": SCHEMA,
        "profile": profile,
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "metrics": metrics,
        "notes": [
            "Exterior-derivative and pullback residuals must be normalized before this gate.",
            "The Hodge check covers geometry/material; it cannot repair missing cohomology.",
            "Holonomy signs depend on the declared orientation and connection convention.",
        ],
    }


def visual_geometry_gate(summary_json: str) -> str:
    """Parse JSON and return a stable, machine-readable gate result."""

    try:
        summary = json.loads(summary_json)
    except json.JSONDecodeError as exc:
        raise ValueError(f"summary_json must be valid JSON: {exc.msg}") from exc
    return json.dumps(evaluate_visual_geometry(summary), indent=2, sort_keys=True)
