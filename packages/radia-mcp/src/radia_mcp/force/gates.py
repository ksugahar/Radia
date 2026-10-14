"""Application-neutral force/torque consistency gates."""

from __future__ import annotations

import math
from typing import Any


def _vector3(values, name: str) -> list[float]:
    try:
        vector = [float(value) for value in values]
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a numeric three-vector") from exc
    if len(vector) != 3 or not all(math.isfinite(value) for value in vector):
        raise ValueError(f"{name} must be a finite numeric three-vector")
    return vector


def _norm(vector: list[float]) -> float:
    return math.sqrt(sum(value * value for value in vector))


def _difference(left: list[float], right: list[float]) -> list[float]:
    return [a - b for a, b in zip(left, right)]


def _sum(left: list[float], right: list[float]) -> list[float]:
    return [a + b for a, b in zip(left, right)]


def _relative(vector: list[float], left: list[float], right: list[float]) -> float:
    return _norm(vector) / max(_norm(left), _norm(right), 1.0e-300)


def electromagnetic_force_method_selection_gate(
    target_kind: str,
    requested_method: str,
    *,
    relative_permeability: float = 1.0,
    weighted_stress_available: bool = False,
    virtual_work_samples_available: bool = False,
    contour_clearance_mesh_layers: int = 0,
) -> dict[str, Any]:
    """Select and gate the robust primary force method for a target."""

    kind = str(target_kind or "").strip().lower().replace("-", "_").replace(" ", "_")
    method = str(requested_method or "").strip().lower().replace("-", "_").replace(" ", "_")
    permeability = float(relative_permeability)
    if not math.isfinite(permeability):
        raise ValueError("relative_permeability must be finite")
    layers = int(contour_clearance_mesh_layers)
    conductors = {"conductor", "current_conductor", "coil", "busbar"}
    magnetic_bodies = {"magnetic_body", "ferromagnetic_body", "magnet", "iron"}
    if kind in conductors and abs(permeability - 1.0) <= 1.0e-12:
        recommended = "lorentz_body_force"
        reason = "unit-permeability conductor supports direct J cross B integration"
    elif kind in magnetic_bodies and weighted_stress_available:
        recommended = "weighted_stress_volume"
        reason = "magnetic-body force should use mesh-robust weighted stress"
    elif kind in magnetic_bodies and virtual_work_samples_available:
        recommended = "coenergy_virtual_work"
        reason = "matched displacement samples support a coenergy derivative"
    else:
        recommended = "none"
        reason = "no robust primary force evidence is available"
    contour = method in {
        "contour_maxwell_stress",
        "maxwell_stress_contour",
        "line_maxwell_stress",
    }
    checks = {
        "target_kind_supported": kind in conductors | magnetic_bodies,
        "relative_permeability_positive": permeability > 0.0,
        "robust_primary_method_available": recommended != "none",
        "requested_method_matches_recommendation": method == recommended,
        "contour_not_used_as_primary": not contour,
        "contour_clearance_recorded_when_requested": not contour or layers > 0,
    }
    return {
        "policy": "radia.force-method-selection/v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "target_kind": kind,
        "requested_method": method,
        "recommended_method": recommended,
        "recommendation_reason": reason,
        "relative_permeability": permeability,
        "checks": checks,
    }


def force_torque_method_agreement_gate(
    primary: dict[str, Any],
    independent: dict[str, Any],
    *,
    maximum_force_relative_difference: float = 0.05,
    maximum_torque_relative_difference: float = 0.05,
) -> dict[str, Any]:
    """Compare two shared force-result records, including frame and pivot."""

    if not isinstance(primary, dict) or not isinstance(independent, dict):
        raise TypeError("primary and independent must be force-result mappings")
    force_tolerance = float(maximum_force_relative_difference)
    torque_tolerance = float(maximum_torque_relative_difference)
    if (
        not math.isfinite(force_tolerance)
        or not math.isfinite(torque_tolerance)
        or force_tolerance < 0.0
        or torque_tolerance < 0.0
    ):
        raise ValueError("relative-difference tolerances must be finite and >= 0")
    force_comparable = (
        primary.get("force_N") is not None and independent.get("force_N") is not None
    )
    torque_comparable = (
        primary.get("torque_Nm") is not None
        and independent.get("torque_Nm") is not None
    )
    checks: dict[str, bool] = {
        "at_least_one_comparable_resultant": force_comparable or torque_comparable,
        "different_methods": primary.get("method") != independent.get("method"),
        "same_frame": primary.get("frame") == independent.get("frame"),
        "same_pivot": primary.get("pivot_m") == independent.get("pivot_m"),
        "same_dimensionality": primary.get("dimensionality") == independent.get("dimensionality"),
        "same_per_unit_depth": primary.get("per_unit_depth") == independent.get("per_unit_depth"),
        "same_field_convention": primary.get("field_convention") == independent.get("field_convention"),
    }
    force_difference = None
    force_relative = None
    if force_comparable:
        force_primary = _vector3(primary["force_N"], "primary.force_N")
        force_independent = _vector3(independent["force_N"], "independent.force_N")
        force_difference = _difference(force_primary, force_independent)
        force_relative = _relative(force_difference, force_primary, force_independent)
        checks["force_agreement"] = force_relative <= force_tolerance
    else:
        checks["force_presence_matches"] = (
            primary.get("force_N") is None and independent.get("force_N") is None
        )
    torque_difference = None
    torque_relative = None
    if torque_comparable:
        torque_primary = _vector3(primary["torque_Nm"], "primary.torque_Nm")
        torque_independent = _vector3(independent["torque_Nm"], "independent.torque_Nm")
        torque_difference = _difference(torque_primary, torque_independent)
        torque_relative = _relative(torque_difference, torque_primary, torque_independent)
        checks["torque_agreement"] = torque_relative <= torque_tolerance
    else:
        checks["torque_presence_matches"] = (
            primary.get("torque_Nm") is None and independent.get("torque_Nm") is None
        )
    return {
        "policy": "radia.force-torque-method-agreement/v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "primary_method": primary.get("method"),
        "independent_method": independent.get("method"),
        "force_difference_N": force_difference,
        "force_relative_difference": force_relative,
        "torque_difference_Nm": torque_difference,
        "torque_relative_difference": torque_relative,
        "maximum_force_relative_difference": force_tolerance,
        "maximum_torque_relative_difference": torque_tolerance,
        "checks": checks,
    }


def force_action_reaction_gate(
    force_a_N,
    force_b_N,
    *,
    torque_a_Nm=None,
    torque_b_Nm=None,
    maximum_force_relative_residual: float = 0.01,
    maximum_torque_relative_residual: float = 0.01,
) -> dict[str, Any]:
    """Gate Newton-pair force closure and optional common-pivot torque closure."""

    force_a = _vector3(force_a_N, "force_a_N")
    force_b = _vector3(force_b_N, "force_b_N")
    force_residual = _sum(force_a, force_b)
    force_relative = _relative(force_residual, force_a, force_b)
    force_tolerance = float(maximum_force_relative_residual)
    torque_tolerance = float(maximum_torque_relative_residual)
    if (
        not math.isfinite(force_tolerance)
        or not math.isfinite(torque_tolerance)
        or force_tolerance < 0.0
        or torque_tolerance < 0.0
    ):
        raise ValueError("relative-residual tolerances must be finite and >= 0")
    checks = {"force_action_reaction": force_relative <= force_tolerance}
    torque_residual = None
    torque_relative = None
    if torque_a_Nm is not None or torque_b_Nm is not None:
        if torque_a_Nm is None or torque_b_Nm is None:
            raise ValueError("torque_a_Nm and torque_b_Nm must be supplied together")
        torque_a = _vector3(torque_a_Nm, "torque_a_Nm")
        torque_b = _vector3(torque_b_Nm, "torque_b_Nm")
        torque_residual = _sum(torque_a, torque_b)
        torque_relative = _relative(torque_residual, torque_a, torque_b)
        checks["common_pivot_torque_action_reaction"] = torque_relative <= torque_tolerance
    return {
        "policy": "radia.force-action-reaction/v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "force_residual_N": force_residual,
        "force_relative_residual": force_relative,
        "torque_residual_Nm": torque_residual,
        "torque_relative_residual": torque_relative,
        "checks": checks,
    }


def force_weight_equilibrium_gate(
    force_N,
    mass_kg: float,
    *,
    lift_axis: int = 2,
    gravity_m_per_s2: float = 9.80665,
    maximum_relative_residual: float = 0.02,
) -> dict[str, Any]:
    """Gate MagLev/static-bearing lift against weight along one axis."""

    force = _vector3(force_N, "force_N")
    mass = float(mass_kg)
    gravity = float(gravity_m_per_s2)
    tolerance = float(maximum_relative_residual)
    axis = int(lift_axis)
    if axis not in (0, 1, 2):
        raise ValueError("lift_axis must be 0, 1, or 2")
    if not math.isfinite(mass) or mass < 0.0:
        raise ValueError("mass_kg must be finite and >= 0")
    if not math.isfinite(gravity) or gravity <= 0.0:
        raise ValueError("gravity_m_per_s2 must be finite and > 0")
    if not math.isfinite(tolerance) or tolerance < 0.0:
        raise ValueError("maximum_relative_residual must be finite and >= 0")
    weight = mass * gravity
    residual = force[axis] - weight
    relative = abs(residual) / max(weight, abs(force[axis]), 1.0e-300)
    checks = {
        "lift_nonnegative": force[axis] >= 0.0,
        "lift_weight_equilibrium": relative <= tolerance,
    }
    return {
        "policy": "radia.force-weight-equilibrium/v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "lift_axis": axis,
        "lift_force_N": force[axis],
        "weight_N": weight,
        "residual_N": residual,
        "relative_residual": relative,
        "checks": checks,
    }
