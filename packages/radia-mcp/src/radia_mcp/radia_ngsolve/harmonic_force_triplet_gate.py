"""Solver-neutral closure gate for three harmonic magnetic-force methods."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


METHOD_ROLES = {
    "material_surface": "body",
    "maxwell_stress": "body",
    "coil_lorentz": "source",
}


def _vector(value: object, name: str) -> tuple[float, float, float]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or len(value) != 3:
        raise ValueError(f"{name} must contain three components")
    result = tuple(float(component) for component in value)
    if not all(math.isfinite(component) for component in result):
        raise ValueError(f"{name} must contain finite components")
    return result


def _norm(vector: Sequence[float]) -> float:
    return math.sqrt(sum(component * component for component in vector))


def harmonic_magnetic_force_triplet_closure_gate(
    summary: Mapping[str, object],
    *,
    maximum_body_method_relative_difference: float = 0.05,
    maximum_action_reaction_relative_residual: float = 0.01,
    maximum_transverse_relative: float = 1.0e-8,
) -> dict[str, object]:
    """Gate body-force agreement and source/body action-reaction closure."""
    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")

    tolerances = {
        "maximum_body_method_relative_difference": float(
            maximum_body_method_relative_difference
        ),
        "maximum_action_reaction_relative_residual": float(
            maximum_action_reaction_relative_residual
        ),
        "maximum_transverse_relative": float(maximum_transverse_relative),
    }
    if not all(math.isfinite(value) and value >= 0.0 for value in tolerances.values()):
        raise ValueError("all tolerances must be finite and nonnegative")

    records = summary.get("methods")
    if not isinstance(records, list) or len(records) != 3:
        raise ValueError("methods must contain exactly three records")
    indexed = {
        str(record.get("method") or ""): record
        for record in records
        if isinstance(record, Mapping)
    }
    if set(indexed) != set(METHOD_ROLES):
        raise ValueError(f"methods must be exactly {sorted(METHOD_ROLES)}")

    vectors = {
        method: _vector(record.get("force"), f"methods[{method}].force")
        for method, record in indexed.items()
    }
    norms = {method: _norm(vector) for method, vector in vectors.items()}
    material = vectors["material_surface"]
    stress = vectors["maxwell_stress"]
    lorentz = vectors["coil_lorentz"]
    if norms["maxwell_stress"] == 0.0:
        raise ValueError("maxwell_stress force must be nonzero")

    dominant_axis = max(range(3), key=lambda index: abs(stress[index]))
    body_denominator = 0.5 * (
        norms["material_surface"] + norms["maxwell_stress"]
    )
    reaction_denominator = 0.5 * (
        norms["maxwell_stress"] + norms["coil_lorentz"]
    )
    body_difference = _norm(
        tuple(material[index] - stress[index] for index in range(3))
    ) / max(body_denominator, 1.0e-300)
    action_reaction_residual = _norm(
        tuple(stress[index] + lorentz[index] for index in range(3))
    ) / max(reaction_denominator, 1.0e-300)
    transverse_relative = {
        method: _norm(
            tuple(vector[index] for index in range(3) if index != dominant_axis)
        )
        / max(abs(vector[dominant_axis]), 1.0e-300)
        for method, vector in vectors.items()
    }

    dimension = str(summary.get("quantity_dimension") or "")
    unit = str(summary.get("force_unit") or "")
    checks = {
        "harmonic_frequency_positive": math.isfinite(float(summary.get("frequency_hz", 0.0)))
        and float(summary.get("frequency_hz", 0.0)) > 0.0,
        "force_dimension_and_unit_consistent": (dimension, unit)
        in {("3d_total", "N"), ("2d_per_length", "N/m")},
        "global_component_frame_recorded": summary.get("component_frame") == "global",
        "method_roles_are_distinct": all(
            indexed[method].get("role") == role for method, role in METHOD_ROLES.items()
        ),
        "all_force_vectors_nonzero": all(value > 0.0 for value in norms.values()),
        "body_methods_share_dominant_axis": max(
            range(3), key=lambda index: abs(material[index])
        )
        == dominant_axis,
        "body_methods_share_dominant_sign": material[dominant_axis]
        * stress[dominant_axis]
        > 0.0,
        "source_force_opposes_body_force": lorentz[dominant_axis]
        * stress[dominant_axis]
        < 0.0,
        "body_method_difference_within_tolerance": body_difference
        <= tolerances["maximum_body_method_relative_difference"],
        "action_reaction_closure_within_tolerance": action_reaction_residual
        <= tolerances["maximum_action_reaction_relative_residual"],
        "transverse_components_within_tolerance": max(transverse_relative.values())
        <= tolerances["maximum_transverse_relative"],
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "harmonic_magnetic_force_triplet_closure_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "metrics": {
            "dominant_axis": "xyz"[dominant_axis],
            "body_method_relative_difference": body_difference,
            "action_reaction_relative_residual": action_reaction_residual,
            "maximum_transverse_relative": max(transverse_relative.values()),
        },
        "tolerances": tolerances,
        "lesson": (
            "A material-surface force and a closed-surface Maxwell-stress force are "
            "independent body-force estimates. The source-region Lorentz force should "
            "oppose the body force and close the action-reaction residual; method names, "
            "roles, units, frame, and dimensional convention must remain explicit."
        ),
    }
