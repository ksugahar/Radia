"""Solver-neutral validation for linear harmonic eddy-current levitation."""

from __future__ import annotations

import math
from typing import Any


def _finite_float(row: dict[str, Any], key: str) -> float:
    value = float(row[key])
    if not math.isfinite(value):
        raise ValueError(f"{key} must be finite")
    return value


def _relative_span(values: list[float]) -> float:
    scale = max((abs(value) for value in values), default=0.0)
    return (max(values) - min(values)) / scale if values and scale > 0.0 else math.inf


def linear_eddy_levitation_force_gate(
    summary: dict[str, Any],
    *,
    max_i2_coefficient_relative_span: float = 1.0e-6,
    max_force_method_relative_difference: float = 0.05,
    max_mesh_count_relative_span: float = 0.0,
    min_sample_count: int = 3,
) -> dict[str, Any]:
    """Validate a fixed-frequency current sweep for a linear conducting body.

    In a linear harmonic model, field amplitude scales with current, while the
    time-average Lorentz force, resistive loss, and stress-derived force scale
    with current squared. The Lorentz and stress routes are independent force
    extractions; agreement is a useful discretization and selection check.
    """
    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")
    contract = summary.get("contract") or {}
    if not isinstance(contract, dict):
        raise ValueError("contract must be a mapping")
    rows = summary.get("rows") or []
    if not isinstance(rows, list):
        raise ValueError("rows must be a list")

    parsed: list[dict[str, float]] = []
    parse_errors: list[str] = []
    for index, row in enumerate(rows):
        if not isinstance(row, dict):
            parse_errors.append(f"row {index} is not a mapping")
            continue
        try:
            parsed.append(
                {
                    "current_a": _finite_float(row, "current_a"),
                    "lorentz_dc_z_n": _finite_float(row, "lorentz_dc_z_n"),
                    "weighted_stress_dc_z_n": _finite_float(
                        row, "weighted_stress_dc_z_n"
                    ),
                    "resistive_loss_w": _finite_float(row, "resistive_loss_w"),
                    "node_count": _finite_float(row, "node_count"),
                    "element_count": _finite_float(row, "element_count"),
                }
            )
        except (KeyError, TypeError, ValueError) as exc:
            parse_errors.append(f"row {index}: {exc}")

    currents = [row["current_a"] for row in parsed]
    unique_increasing = bool(currents) and all(
        right > left > 0.0 for left, right in zip(currents, currents[1:])
    )
    coefficients: dict[str, list[float]] = {
        "lorentz_force_n_per_a2": [],
        "weighted_stress_force_n_per_a2": [],
        "loss_w_per_a2": [],
    }
    force_method_differences: list[float] = []
    for row in parsed:
        current_sq = row["current_a"] ** 2
        coefficients["lorentz_force_n_per_a2"].append(
            row["lorentz_dc_z_n"] / current_sq
        )
        coefficients["weighted_stress_force_n_per_a2"].append(
            row["weighted_stress_dc_z_n"] / current_sq
        )
        coefficients["loss_w_per_a2"].append(row["resistive_loss_w"] / current_sq)
        force_method_differences.append(
            abs(row["lorentz_dc_z_n"] - row["weighted_stress_dc_z_n"])
            / max(
                abs(row["lorentz_dc_z_n"]),
                abs(row["weighted_stress_dc_z_n"]),
                1.0e-30,
            )
        )

    coefficient_spans = {
        name: _relative_span(values) for name, values in coefficients.items()
    }
    node_counts = [row["node_count"] for row in parsed]
    element_counts = [row["element_count"] for row in parsed]
    node_span = _relative_span(node_counts)
    element_span = _relative_span(element_counts)
    force_signs_agree = bool(parsed) and all(
        row["lorentz_dc_z_n"] * row["weighted_stress_dc_z_n"] > 0.0
        for row in parsed
    )

    checks = {
        "rows_parsed_and_finite": not parse_errors and len(parsed) == len(rows),
        "sample_count_sufficient": len(parsed) >= int(min_sample_count),
        "positive_currents_strictly_increase": unique_increasing,
        "fixed_positive_frequency": math.isfinite(
            float(contract.get("frequency_hz", math.nan))
        )
        and float(contract.get("frequency_hz", math.nan)) > 0.0,
        "linear_material_contract_recorded": contract.get("linear_materials") is True,
        "dc_time_average_force_component": contract.get("force_component")
        == "time_average_dc",
        "conducting_target_recorded": contract.get("target_kind")
        == "conducting_body",
        "weighted_stress_mask_valid": contract.get("weighted_stress_mask")
        == "target_surrounded_by_air",
        "positive_force_and_loss": bool(parsed)
        and all(
            row["lorentz_dc_z_n"] > 0.0
            and row["weighted_stress_dc_z_n"] > 0.0
            and row["resistive_loss_w"] > 0.0
            for row in parsed
        ),
        "force_methods_have_same_sign": force_signs_agree,
        "lorentz_force_obeys_i2": coefficient_spans[
            "lorentz_force_n_per_a2"
        ]
        <= float(max_i2_coefficient_relative_span),
        "weighted_stress_force_obeys_i2": coefficient_spans[
            "weighted_stress_force_n_per_a2"
        ]
        <= float(max_i2_coefficient_relative_span),
        "resistive_loss_obeys_i2": coefficient_spans["loss_w_per_a2"]
        <= float(max_i2_coefficient_relative_span),
        "force_methods_agree": bool(force_method_differences)
        and max(force_method_differences)
        <= float(max_force_method_relative_difference),
        "mesh_inventory_positive_and_stable": bool(parsed)
        and min(node_counts) > 0.0
        and min(element_counts) > 0.0
        and node_span <= float(max_mesh_count_relative_span)
        and element_span <= float(max_mesh_count_relative_span),
    }
    issues = parse_errors + [name for name, ok in checks.items() if not ok]
    return {
        "policy": "linear_eddy_levitation_force_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "metrics": {
            "sample_count": len(parsed),
            "i2_coefficient_relative_spans": coefficient_spans,
            "max_force_method_relative_difference": max(
                force_method_differences, default=math.inf
            ),
            "node_count_relative_span": node_span,
            "element_count_relative_span": element_span,
        },
        "tolerances": {
            "max_i2_coefficient_relative_span": float(
                max_i2_coefficient_relative_span
            ),
            "max_force_method_relative_difference": float(
                max_force_method_relative_difference
            ),
            "max_mesh_count_relative_span": float(max_mesh_count_relative_span),
            "min_sample_count": int(min_sample_count),
        },
        "notes": [
            "Compare like force components: a harmonic time-average DC force is not a 2x force phasor.",
            "I-squared scaling is expected only for a fixed-frequency linear-material sweep.",
            "A weighted-stress body mask is valid when the selected target is surrounded by air; otherwise prefer virtual work or a Lorentz-volume route.",
        ],
    }
