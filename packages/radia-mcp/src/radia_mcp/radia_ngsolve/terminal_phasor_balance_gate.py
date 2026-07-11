"""Solver-neutral balance checks for cyclic multi-terminal phasors."""

from __future__ import annotations

import cmath
import math
from typing import Any


def _phasor(value: Any, name: str) -> complex:
    if isinstance(value, dict):
        if "real" not in value or "imag" not in value:
            raise ValueError(f"{name} must contain real and imag")
        result = complex(float(value["real"]), float(value["imag"]))
    elif isinstance(value, (list, tuple)) and len(value) == 2:
        result = complex(float(value[0]), float(value[1]))
    elif isinstance(value, (int, float)) and not isinstance(value, bool):
        result = complex(float(value), 0.0)
    else:
        raise ValueError(f"{name} must be a real/imag mapping or pair")
    if not (math.isfinite(result.real) and math.isfinite(result.imag)):
        raise ValueError(f"{name} must be finite")
    return result


def _relative_spread(values: list[complex]) -> float:
    magnitudes = [abs(value) for value in values]
    scale = max(magnitudes, default=0.0)
    if scale <= 0.0:
        return math.inf
    return (max(magnitudes) - min(magnitudes)) / scale


def _zero_sequence_residual(values: list[complex]) -> float:
    scale = sum(abs(value) for value in values)
    return abs(sum(values)) / scale if scale > 0.0 else math.inf


def _wrapped_angle_error_deg(actual: float, expected: float) -> float:
    return abs((actual - expected + 180.0) % 360.0 - 180.0)


def _phase_step_errors_deg(values: list[complex], expected_step_deg: float) -> list[float]:
    errors = []
    for left, right in zip(values, values[1:]):
        if abs(left) <= 0.0 or abs(right) <= 0.0:
            return [math.inf]
        actual = math.degrees(cmath.phase(right / left))
        errors.append(_wrapped_angle_error_deg(actual, expected_step_deg))
    return errors


def cyclic_terminal_phasor_balance_gate(
    summary: dict[str, Any],
    max_magnitude_relative_spread: float = 1.0e-5,
    max_phase_step_error_deg: float = 1.0e-2,
    max_zero_sequence_residual: float = 1.0e-5,
    max_terminal_kcl_residual: float = 1.0e-5,
    max_reference_current_relative_error: float = 2.0e-2,
) -> dict[str, Any]:
    """Gate cyclic voltage/current triplets and the all-terminal KCL residual.

    ``summary["groups"]`` contains at least two three-terminal groups.  Each
    group records ``voltage_phasors`` and outward-positive
    ``current_phasors`` as JSON-safe ``[real, imag]`` pairs or mappings.  The
    gate checks magnitude balance, the cyclic phase step, each triplet's
    zero-sequence residual, and cancellation of current over all terminals.
    """
    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")
    groups = summary.get("groups")
    if not isinstance(groups, list) or len(groups) < 2:
        raise ValueError("at least two terminal groups are required")

    tolerances = {
        "max_magnitude_relative_spread": float(max_magnitude_relative_spread),
        "max_phase_step_error_deg": float(max_phase_step_error_deg),
        "max_zero_sequence_residual": float(max_zero_sequence_residual),
        "max_terminal_kcl_residual": float(max_terminal_kcl_residual),
        "max_reference_current_relative_error": float(
            max_reference_current_relative_error
        ),
    }
    if any(not math.isfinite(value) or value < 0.0 for value in tolerances.values()):
        raise ValueError("all tolerances must be finite and nonnegative")

    expected_step = float(summary.get("expected_phase_step_deg", -120.0))
    if not math.isfinite(expected_step):
        raise ValueError("expected_phase_step_deg must be finite")
    voltage_unit = str(summary.get("voltage_unit") or "").strip()
    current_unit = str(summary.get("current_unit") or "").strip()

    checks: dict[str, bool] = {
        "voltage_unit_recorded": bool(voltage_unit),
        "current_unit_recorded": bool(current_unit),
    }
    labels: list[str] = []
    all_currents: list[complex] = []
    group_metrics: list[dict[str, Any]] = []

    for index, row in enumerate(groups, start=1):
        if not isinstance(row, dict):
            raise ValueError(f"group {index} must be a mapping")
        label = str(row.get("label") or "").strip()
        labels.append(label)
        voltages = [
            _phasor(value, f"group {index} voltage {item}")
            for item, value in enumerate(row.get("voltage_phasors") or [], start=1)
        ]
        currents = [
            _phasor(value, f"group {index} current {item}")
            for item, value in enumerate(row.get("current_phasors") or [], start=1)
        ]
        if len(voltages) != 3 or len(currents) != 3:
            raise ValueError(f"group {index} must contain three voltages and currents")
        all_currents.extend(currents)

        voltage_spread = _relative_spread(voltages)
        current_spread = _relative_spread(currents)
        voltage_phase_errors = _phase_step_errors_deg(voltages, expected_step)
        current_phase_errors = _phase_step_errors_deg(currents, expected_step)
        voltage_zero = _zero_sequence_residual(voltages)
        current_zero = _zero_sequence_residual(currents)
        mean_current = sum(abs(value) for value in currents) / 3.0

        prefix = f"group_{index}"
        checks[f"{prefix}_label_recorded"] = bool(label)
        checks[f"{prefix}_voltage_magnitudes_balanced"] = (
            voltage_spread <= tolerances["max_magnitude_relative_spread"]
        )
        checks[f"{prefix}_current_magnitudes_balanced"] = (
            current_spread <= tolerances["max_magnitude_relative_spread"]
        )
        checks[f"{prefix}_voltage_phase_sequence"] = (
            max(voltage_phase_errors) <= tolerances["max_phase_step_error_deg"]
        )
        checks[f"{prefix}_current_phase_sequence"] = (
            max(current_phase_errors) <= tolerances["max_phase_step_error_deg"]
        )
        checks[f"{prefix}_voltage_zero_sequence_small"] = (
            voltage_zero <= tolerances["max_zero_sequence_residual"]
        )
        checks[f"{prefix}_current_zero_sequence_small"] = (
            current_zero <= tolerances["max_zero_sequence_residual"]
        )

        reference_error = None
        if "reference_current_magnitude" in row:
            reference_unit = str(row.get("reference_current_unit") or "").strip()
            if reference_unit != current_unit:
                raise ValueError(
                    f"group {index} reference current unit must match current_unit"
                )
            reference = float(row["reference_current_magnitude"])
            if not math.isfinite(reference) or reference <= 0.0:
                raise ValueError(f"group {index} reference current must be positive")
            reference_error = abs(mean_current - reference) / reference
            checks[f"{prefix}_reference_current_matches"] = (
                reference_error
                <= tolerances["max_reference_current_relative_error"]
            )

        group_metrics.append(
            {
                "label": label,
                "voltage_magnitude_mean": sum(abs(value) for value in voltages) / 3.0,
                "current_magnitude_mean": mean_current,
                "voltage_magnitude_relative_spread": voltage_spread,
                "current_magnitude_relative_spread": current_spread,
                "max_voltage_phase_step_error_deg": max(voltage_phase_errors),
                "max_current_phase_step_error_deg": max(current_phase_errors),
                "voltage_zero_sequence_residual": voltage_zero,
                "current_zero_sequence_residual": current_zero,
                "reference_current_relative_error": reference_error,
            }
        )

    checks["group_labels_distinct"] = len(set(labels)) == len(labels) and all(labels)
    terminal_scale = sum(abs(value) for value in all_currents)
    terminal_kcl = (
        abs(sum(all_currents)) / terminal_scale
        if terminal_scale > 0.0
        else math.inf
    )
    checks["all_terminal_kcl_closes"] = (
        terminal_kcl <= tolerances["max_terminal_kcl_residual"]
    )

    return {
        "policy": "cyclic_terminal_phasor_balance_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "groups": group_metrics,
            "all_terminal_kcl_residual": terminal_kcl,
            "terminal_count": len(all_currents),
        },
        "units": {"voltage": voltage_unit, "current": current_unit},
        "expected_phase_step_deg": expected_step,
        "tolerances": tolerances,
        "notes": [
            "Terminal currents must use one consistent outward-positive convention.",
            "Cyclic symmetry and zero-sequence checks diagnose each triplet; all-terminal KCL diagnoses coupling between groups.",
            "An optional analytic current magnitude is a secondary approximation check, not a replacement for KCL.",
            "Never compare a per-length current with a terminal current unless an explicit extrusion depth converts the units.",
        ],
    }
