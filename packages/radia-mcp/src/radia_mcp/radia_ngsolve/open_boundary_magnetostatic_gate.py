"""Gauge-invariant equivalence gate for magnetostatic open-boundary formulations."""

from __future__ import annotations

import math
from typing import Any


def _finite_list(value: Any, name: str) -> list[float]:
    try:
        parsed = [float(item) for item in value]
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a numeric sequence") from exc
    if not parsed or not all(math.isfinite(item) for item in parsed):
        raise ValueError(f"{name} must contain finite values")
    return parsed


def _relative_error(left: float, right: float, floor: float = 1.0e-30) -> float:
    return abs(left - right) / max(abs(left), abs(right), floor)


def magnetostatic_open_boundary_equivalence_gate(
    formulation_rows: list[dict[str, Any]],
    *,
    physics_regime: str = "magnetostatic_open_boundary",
    axis_sample_indices: list[int] | None = None,
    max_dominant_b_relative_error: float = 0.01,
    max_axis_transverse_b_residual: float = 0.015,
    max_a_offset_relative_spread: float = 0.005,
    max_energy_coenergy_relative_error: float = 0.001,
    max_dominant_force_relative_error: float = 0.001,
    max_force_balance_relative: float = 0.002,
    max_transverse_force_difference_relative: float = 0.001,
) -> dict[str, Any]:
    """Compare two open-boundary solutions using gauge-invariant observables.

    ``A`` is intentionally treated only through the consistency of its additive
    offset. Direct relative error in ``A`` is reported as a diagnostic and never
    used as a pass/fail criterion.
    """

    if physics_regime != "magnetostatic_open_boundary":
        raise ValueError(
            "this gate is only for magnetostatic open boundaries; it does not "
            "select or validate radiation boundaries for wave problems"
        )
    if not isinstance(formulation_rows, list) or len(formulation_rows) != 2:
        raise ValueError("exactly two formulation rows are required")
    tolerances = [
        max_dominant_b_relative_error,
        max_axis_transverse_b_residual,
        max_a_offset_relative_spread,
        max_energy_coenergy_relative_error,
        max_dominant_force_relative_error,
        max_force_balance_relative,
        max_transverse_force_difference_relative,
    ]
    if any(not math.isfinite(float(value)) or float(value) < 0.0 for value in tolerances):
        raise ValueError("tolerances must be finite and nonnegative")

    rows: list[dict[str, Any]] = []
    for index, row in enumerate(formulation_rows):
        try:
            identifier = str(row["id"]).strip()
            a_values = _finite_list(row["a"], f"row {index}.a")
            bx_values = _finite_list(row["bx"], f"row {index}.bx")
            by_values = _finite_list(row["by"], f"row {index}.by")
            energy = _finite_list(row["energy"], f"row {index}.energy")
            coenergy = _finite_list(row["coenergy"], f"row {index}.coenergy")
            forces = [
                {
                    "x": float(force["x"]),
                    "y": float(force["y"]),
                }
                for force in row["forces"]
            ]
            mesh_nodes = int(row["mesh"]["nodes"])
            mesh_elements = int(row["mesh"]["elements"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"formulation row {index} is malformed") from exc
        if not identifier:
            raise ValueError(f"formulation row {index} has an empty id")
        sample_count = len(a_values)
        if sample_count < 3 or len(bx_values) != sample_count or len(by_values) != sample_count:
            raise ValueError("A, Bx, and By arrays must have the same length of at least three")
        if len(forces) < 2 or not all(
            math.isfinite(force[component]) for force in forces for component in ("x", "y")
        ):
            raise ValueError("each formulation needs at least two finite force vectors")
        rows.append(
            {
                "id": identifier,
                "a": a_values,
                "bx": bx_values,
                "by": by_values,
                "energy": energy,
                "coenergy": coenergy,
                "forces": forces,
                "mesh": {"nodes": mesh_nodes, "elements": mesh_elements},
            }
        )

    if rows[0]["id"] == rows[1]["id"]:
        raise ValueError("formulation ids must be unique")
    if len(rows[0]["energy"]) != len(rows[1]["energy"]) or len(rows[0]["coenergy"]) != len(rows[1]["coenergy"]):
        raise ValueError("energy and coenergy arrays must match across formulations")
    if len(rows[0]["forces"]) != len(rows[1]["forces"]):
        raise ValueError("force-vector counts must match across formulations")

    sample_count = len(rows[0]["a"])
    if len(rows[1]["a"]) != sample_count:
        raise ValueError("field sample counts must match across formulations")
    axes = [0, 1, 2] if axis_sample_indices is None else [int(value) for value in axis_sample_indices]
    if not axes or any(value < 0 or value >= sample_count for value in axes):
        raise ValueError("axis_sample_indices must select valid field samples")

    left, right = rows
    a_offsets = [right["a"][i] - left["a"][i] for i in range(sample_count)]
    a_offset_mean = sum(a_offsets) / len(a_offsets)
    a_offset_span = max(a_offsets) - min(a_offsets)
    if a_offset_span == 0.0:
        a_offset_relative_spread = 0.0
    else:
        a_offset_relative_spread = a_offset_span / max(abs(a_offset_mean), 1.0e-30)
    direct_a_relative_errors = [
        _relative_error(left["a"][i], right["a"][i], 1.0e-30)
        for i in range(sample_count)
    ]

    dominant_b_errors: list[float] = []
    dominant_b_components: list[str] = []
    for index in range(sample_count):
        bx_scale = max(abs(left["bx"][index]), abs(right["bx"][index]))
        by_scale = max(abs(left["by"][index]), abs(right["by"][index]))
        component = "bx" if bx_scale >= by_scale else "by"
        dominant_b_components.append(component)
        dominant_b_errors.append(
            _relative_error(left[component][index], right[component][index], 1.0e-30)
        )

    axis_residuals = [
        abs(row["by"][index])
        / max(abs(row["bx"][index]), abs(row["by"][index]), 1.0e-30)
        for row in rows
        for index in axes
    ]
    energy_errors = [
        _relative_error(a, b) for a, b in zip(left["energy"], right["energy"])
    ]
    coenergy_errors = [
        _relative_error(a, b) for a, b in zip(left["coenergy"], right["coenergy"])
    ]
    dominant_force_errors = [
        _relative_error(a["y"], b["y"])
        for a, b in zip(left["forces"], right["forces"])
    ]
    dominant_force_scale = max(
        abs(force["y"]) for row in rows for force in row["forces"]
    )
    force_balance = [
        math.hypot(
            sum(force["x"] for force in row["forces"]),
            sum(force["y"] for force in row["forces"]),
        )
        / max(max(abs(force["y"]) for force in row["forces"]), 1.0e-30)
        for row in rows
    ]
    transverse_difference = max(
        abs(a["x"] - b["x"])
        for a, b in zip(left["forces"], right["forces"])
    ) / max(dominant_force_scale, 1.0e-30)

    checks = {
        "two_distinct_formulations": left["id"] != right["id"],
        "positive_distinct_mesh_inventories": all(
            row["mesh"]["nodes"] > 0 and row["mesh"]["elements"] > 0 for row in rows
        )
        and left["mesh"] != right["mesh"],
        "vector_potential_additive_offset_consistent": a_offset_relative_spread
        <= max_a_offset_relative_spread,
        "dominant_flux_density_components_agree": max(dominant_b_errors)
        <= max_dominant_b_relative_error,
        "axis_transverse_flux_density_is_bounded": max(axis_residuals)
        <= max_axis_transverse_b_residual,
        "energy_and_coenergy_agree": max(energy_errors + coenergy_errors)
        <= max_energy_coenergy_relative_error,
        "dominant_force_components_agree": max(dominant_force_errors)
        <= max_dominant_force_relative_error,
        "whole_model_force_balance_is_bounded": max(force_balance)
        <= max_force_balance_relative,
        "small_transverse_force_is_scaled_by_dominant_force": transverse_difference
        <= max_transverse_force_difference_relative,
    }
    return {
        "policy": "magnetostatic_open_boundary_equivalence_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "physics_scope": "magnetostatic_open_boundary",
        "wave_boundary_policy": "not_applicable_do_not_infer_from_this_gate",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "sample_count": sample_count,
            "axis_sample_indices": axes,
            "dominant_b_components": dominant_b_components,
            "max_direct_a_relative_error_diagnostic_only": max(direct_a_relative_errors),
            "a_offset_mean": a_offset_mean,
            "a_offset_relative_spread": a_offset_relative_spread,
            "max_dominant_b_relative_error": max(dominant_b_errors),
            "max_axis_transverse_b_residual": max(axis_residuals),
            "max_energy_coenergy_relative_error": max(energy_errors + coenergy_errors),
            "max_dominant_force_relative_error": max(dominant_force_errors),
            "max_force_balance_relative": max(force_balance),
            "transverse_force_difference_relative": transverse_difference,
            "mesh_inventories": {row["id"]: row["mesh"] for row in rows},
        },
        "lesson": (
            "Open-boundary magnetostatic formulations may choose different additive gauges "
            "for vector potential. Compare the offset consistency, B, energy/coenergy, and "
            "properly scaled force balances; never reject the pair from direct A error alone."
        ),
    }
