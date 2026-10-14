"""Solver-neutral scaling contracts for linear magnetostatic sources."""

from __future__ import annotations

import math
from typing import Any


def _finite(value: Any, name: str, *, positive: bool = False) -> float:
    result = float(value)
    if not math.isfinite(result) or (positive and result <= 0.0):
        qualifier = "finite and positive" if positive else "finite"
        raise ValueError(f"{name} must be {qualifier}")
    return result


def _relative_error(actual: float, expected: float) -> float:
    return abs(actual - expected) / max(abs(actual), abs(expected), 1.0e-300)


def linear_magnetization_scaling_gate(
    summary: dict[str, Any],
    *,
    max_field_relative_error: float = 1.0e-10,
    max_energy_relative_error: float = 1.0e-10,
    max_independent_observable_relative_error: float = 3.0e-2,
    max_independent_refinement_relative_change: float = 2.0e-2,
    max_independent_linear_residual: float = 1.0e-10,
) -> dict[str, Any]:
    """Gate linear field and quadratic energy response to source scaling.

    For an unchanged linear magnetostatic problem, scaling every prescribed
    magnetization source by ``alpha`` must scale ``A`` and ``B`` by ``alpha``
    and field energy/coenergy by ``alpha**2``.  Compact artifacts should also
    carry vector-relative errors; scalar maxima alone cannot prove fieldwise
    scaling.
    """
    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")
    contract = summary.get("model_contract")
    full = summary.get("full_scale")
    scaled = summary.get("scaled")
    errors = summary.get("fieldwise_errors")
    independent = summary.get("independent_reference")
    units = summary.get("units")
    if not all(
        isinstance(value, dict)
        for value in (contract, full, scaled, errors, independent, units)
    ):
        raise ValueError(
            "model_contract, full_scale, scaled, fieldwise_errors, "
            "independent_reference and units must be mappings"
        )

    alpha = _finite(scaled.get("source_scale"), "scaled.source_scale", positive=True)
    if alpha >= 1.0:
        raise ValueError("scaled.source_scale must be between zero and one")
    limits = {
        "max_field_relative_error": _finite(
            max_field_relative_error, "max_field_relative_error"
        ),
        "max_energy_relative_error": _finite(
            max_energy_relative_error, "max_energy_relative_error"
        ),
        "max_independent_observable_relative_error": _finite(
            max_independent_observable_relative_error,
            "max_independent_observable_relative_error",
        ),
        "max_independent_refinement_relative_change": _finite(
            max_independent_refinement_relative_change,
            "max_independent_refinement_relative_change",
        ),
        "max_independent_linear_residual": _finite(
            max_independent_linear_residual, "max_independent_linear_residual"
        ),
    }
    if min(limits.values()) < 0.0:
        raise ValueError("tolerances must be nonnegative")

    full_a = _finite(full.get("max_abs_a"), "full_scale.max_abs_a", positive=True)
    full_b = _finite(full.get("max_b"), "full_scale.max_b", positive=True)
    full_energy = _finite(full.get("energy"), "full_scale.energy", positive=True)
    full_coenergy = _finite(full.get("coenergy"), "full_scale.coenergy", positive=True)
    scaled_a = _finite(scaled.get("max_abs_a"), "scaled.max_abs_a", positive=True)
    scaled_b = _finite(scaled.get("max_b"), "scaled.max_b", positive=True)
    scaled_energy = _finite(scaled.get("energy"), "scaled.energy", positive=True)
    scaled_coenergy = _finite(scaled.get("coenergy"), "scaled.coenergy", positive=True)
    a_field_error = _finite(errors.get("a_relative"), "fieldwise_errors.a_relative")
    b_field_error = _finite(errors.get("b_relative"), "fieldwise_errors.b_relative")
    independent_energy_error = _finite(
        independent.get("energy_relative_error"),
        "independent_reference.energy_relative_error",
    )
    independent_a_error = _finite(
        independent.get("max_abs_a_relative_error"),
        "independent_reference.max_abs_a_relative_error",
    )
    independent_refinement = _finite(
        independent.get("maximum_refinement_relative_change"),
        "independent_reference.maximum_refinement_relative_change",
    )
    independent_residual = _finite(
        independent.get("maximum_linear_residual"),
        "independent_reference.maximum_linear_residual",
    )

    metrics = {
        "source_scale": alpha,
        "a_amplitude_ratio": scaled_a / full_a,
        "b_amplitude_ratio": scaled_b / full_b,
        "energy_ratio": scaled_energy / full_energy,
        "coenergy_ratio": scaled_coenergy / full_coenergy,
        "a_fieldwise_relative_error": a_field_error,
        "b_fieldwise_relative_error": b_field_error,
        "full_energy_coenergy_relative_error": _relative_error(full_energy, full_coenergy),
        "scaled_energy_coenergy_relative_error": _relative_error(
            scaled_energy, scaled_coenergy
        ),
        "independent_energy_relative_error": independent_energy_error,
        "independent_max_abs_a_relative_error": independent_a_error,
        "independent_maximum_refinement_relative_change": independent_refinement,
        "independent_maximum_linear_residual": independent_residual,
    }
    checks = {
        "linear_material_contract": contract.get("all_materials_linear") is True,
        "same_mesh_and_boundary_contract": contract.get("same_mesh") is True
        and contract.get("same_boundary_conditions") is True
        and contract.get("only_source_scale_changed") is True,
        "magnetostatic_source_contract": str(contract.get("physics") or "").lower()
        == "magnetostatic"
        and str(contract.get("source") or "").lower() == "prescribed_magnetization",
        "si_observable_units_explicit": units
        == {
            "magnetic_vector_potential": "Wb/m",
            "magnetic_flux_density": "T",
            "energy": "J",
            "coenergy": "J",
        },
        "a_field_scales_linearly": a_field_error <= limits["max_field_relative_error"]
        and _relative_error(metrics["a_amplitude_ratio"], alpha)
        <= limits["max_field_relative_error"],
        "b_field_scales_linearly": b_field_error <= limits["max_field_relative_error"]
        and _relative_error(metrics["b_amplitude_ratio"], alpha)
        <= limits["max_field_relative_error"],
        "energy_scales_quadratically": _relative_error(
            metrics["energy_ratio"], alpha * alpha
        )
        <= limits["max_energy_relative_error"],
        "coenergy_scales_quadratically": _relative_error(
            metrics["coenergy_ratio"], alpha * alpha
        )
        <= limits["max_energy_relative_error"],
        "energy_equals_coenergy_for_linear_problem": metrics[
            "full_energy_coenergy_relative_error"
        ]
        <= limits["max_energy_relative_error"]
        and metrics["scaled_energy_coenergy_relative_error"]
        <= limits["max_energy_relative_error"],
        "independent_first_order_fem_recorded": str(
            independent.get("solver_family") or ""
        ).lower()
        == "independent_fem"
        and str(independent.get("space") or "").upper() == "H1 P1"
        and independent.get("separate_implementation") is True,
        "independent_observables_agree": independent_energy_error
        <= limits["max_independent_observable_relative_error"]
        and independent_a_error
        <= limits["max_independent_observable_relative_error"],
        "independent_refinement_is_stable": independent_refinement
        <= limits["max_independent_refinement_relative_change"],
        "independent_linear_residual_closes": independent_residual
        <= limits["max_independent_linear_residual"],
    }
    return {
        "policy": "linear_magnetization_scaling_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": metrics,
        "tolerances": limits,
        "lesson": (
            "On one unchanged linear magnetostatic discretization, prescribed "
            "magnetization scaling is linear in A and B and quadratic in energy; "
            "energy equals coenergy, and a separately assembled first-order FEM "
            "must agree after its own refinement and residual checks."
        ),
    }
