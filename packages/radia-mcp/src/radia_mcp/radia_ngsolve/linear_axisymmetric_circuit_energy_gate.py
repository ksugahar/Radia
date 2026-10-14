"""Solver-neutral identities for a linear axisymmetric winding sweep."""

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


def _relative_span(values: list[float]) -> float:
    return (max(values) - min(values)) / max(
        max(abs(value) for value in values), 1.0e-300
    )


def linear_axisymmetric_circuit_energy_gate(
    summary: dict[str, Any],
    *,
    maximum_scaling_relative_span: float = 1.0e-10,
    maximum_energy_identity_relative_error: float = 5.0e-4,
    maximum_energy_coenergy_relative_error: float = 1.0e-10,
    maximum_current_relative_error: float = 1.0e-10,
) -> dict[str, Any]:
    """Gate current, flux, field, and energy identities on one fixed mesh.

    A positive current sweep through an unchanged linear axisymmetric winding
    must keep ``lambda/I`` and ``B/I`` constant, keep ``W/I**2`` constant,
    and satisfy ``W = W'``.  The independently accumulated field-energy
    integral must also approach ``0.5*I*lambda`` within a mesh-discrete
    tolerance rather than being treated as an exact scalar identity.
    """
    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")
    contract = summary.get("model_contract")
    units = summary.get("units")
    rows = summary.get("rows")
    if not isinstance(contract, dict) or not isinstance(units, dict):
        raise ValueError("model_contract and units must be mappings")
    if not isinstance(rows, list) or len(rows) < 3:
        raise ValueError("rows must contain at least three current levels")

    limits = {
        "maximum_scaling_relative_span": _finite(
            maximum_scaling_relative_span, "maximum_scaling_relative_span"
        ),
        "maximum_energy_identity_relative_error": _finite(
            maximum_energy_identity_relative_error,
            "maximum_energy_identity_relative_error",
        ),
        "maximum_energy_coenergy_relative_error": _finite(
            maximum_energy_coenergy_relative_error,
            "maximum_energy_coenergy_relative_error",
        ),
        "maximum_current_relative_error": _finite(
            maximum_current_relative_error, "maximum_current_relative_error"
        ),
    }
    if min(limits.values()) < 0.0:
        raise ValueError("tolerances must be nonnegative")

    parsed: list[dict[str, float | int]] = []
    for index, row in enumerate(rows):
        if not isinstance(row, dict):
            raise ValueError(f"rows[{index}] must be a mapping")
        parsed.append(
            {
                "current_a": _finite(row.get("current_a"), f"rows[{index}].current_a", positive=True),
                "circuit_current_a": _finite(
                    row.get("circuit_current_a"),
                    f"rows[{index}].circuit_current_a",
                    positive=True,
                ),
                "flux_linkage_wb_turn": _finite(
                    row.get("flux_linkage_wb_turn"),
                    f"rows[{index}].flux_linkage_wb_turn",
                    positive=True,
                ),
                "magnetic_flux_density_t": _finite(
                    row.get("magnetic_flux_density_t"),
                    f"rows[{index}].magnetic_flux_density_t",
                    positive=True,
                ),
                "energy_j": _finite(row.get("energy_j"), f"rows[{index}].energy_j", positive=True),
                "coenergy_j": _finite(
                    row.get("coenergy_j"), f"rows[{index}].coenergy_j", positive=True
                ),
                "node_count": int(row.get("node_count")),
                "element_count": int(row.get("element_count")),
            }
        )
    currents = [float(row["current_a"]) for row in parsed]
    solved_currents = [float(row["circuit_current_a"]) for row in parsed]
    flux_per_current = [
        float(row["flux_linkage_wb_turn"]) / solved
        for row, solved in zip(parsed, solved_currents, strict=True)
    ]
    field_per_current = [
        float(row["magnetic_flux_density_t"]) / solved
        for row, solved in zip(parsed, solved_currents, strict=True)
    ]
    energy_per_current_squared = [
        float(row["energy_j"]) / solved**2
        for row, solved in zip(parsed, solved_currents, strict=True)
    ]
    current_errors = [
        _relative_error(solved, command)
        for command, solved in zip(currents, solved_currents, strict=True)
    ]
    energy_identity_errors = [
        _relative_error(
            float(row["energy_j"]),
            0.5 * solved * float(row["flux_linkage_wb_turn"]),
        )
        for row, solved in zip(parsed, solved_currents, strict=True)
    ]
    energy_coenergy_errors = [
        _relative_error(float(row["energy_j"]), float(row["coenergy_j"]))
        for row in parsed
    ]
    metrics = {
        "current_relative_error_max": max(current_errors),
        "flux_per_current_relative_span": _relative_span(flux_per_current),
        "field_per_current_relative_span": _relative_span(field_per_current),
        "energy_per_current_squared_relative_span": _relative_span(
            energy_per_current_squared
        ),
        "energy_identity_relative_error_max": max(energy_identity_errors),
        "energy_coenergy_relative_error_max": max(energy_coenergy_errors),
        "reference_inductance_h": flux_per_current[len(flux_per_current) // 2],
    }
    checks = {
        "linear_axisymmetric_winding_contract": contract
        == {
            "physics": "magnetostatic",
            "coordinate_system": "axisymmetric",
            "all_materials_linear": True,
            "same_mesh": True,
            "same_boundary_conditions": True,
            "only_circuit_current_changed": True,
        },
        "si_units_explicit": units
        == {
            "current": "A",
            "flux_linkage": "Wb-turn",
            "magnetic_flux_density": "T",
            "energy": "J",
            "coenergy": "J",
        },
        "positive_strictly_increasing_current_sweep": all(
            right > left for left, right in zip(currents, currents[1:])
        ),
        "circuit_current_matches_command": metrics["current_relative_error_max"]
        <= limits["maximum_current_relative_error"],
        "flux_linkage_scales_with_current": metrics[
            "flux_per_current_relative_span"
        ]
        <= limits["maximum_scaling_relative_span"],
        "field_scales_with_current": metrics["field_per_current_relative_span"]
        <= limits["maximum_scaling_relative_span"],
        "energy_scales_with_current_squared": metrics[
            "energy_per_current_squared_relative_span"
        ]
        <= limits["maximum_scaling_relative_span"],
        "energy_matches_half_i_lambda_within_mesh_tolerance": metrics[
            "energy_identity_relative_error_max"
        ]
        <= limits["maximum_energy_identity_relative_error"],
        "energy_equals_coenergy": metrics["energy_coenergy_relative_error_max"]
        <= limits["maximum_energy_coenergy_relative_error"],
        "fixed_positive_mesh_inventory": all(
            int(row["node_count"]) > 0 and int(row["element_count"]) > 0
            for row in parsed
        )
        and len({int(row["node_count"]) for row in parsed}) == 1
        and len({int(row["element_count"]) for row in parsed}) == 1,
    }
    return {
        "policy": "linear_axisymmetric_circuit_energy_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": metrics,
        "tolerances": limits,
        "lesson": (
            "On a fixed linear axisymmetric winding mesh, lambda and B scale "
            "with current, energy scales with current squared, energy equals "
            "coenergy, and the field-energy integral approaches 0.5 I lambda "
            "within an explicit mesh-discrete tolerance."
        ),
    }
