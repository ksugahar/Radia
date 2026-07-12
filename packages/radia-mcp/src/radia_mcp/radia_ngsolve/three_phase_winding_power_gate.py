"""Solver-neutral power and balance checks for coupled three-phase windings."""

from __future__ import annotations

import cmath
import math
from typing import Any


def _phasor(value: Any, name: str) -> complex:
    if isinstance(value, dict) and {"real", "imag"} <= set(value):
        result = complex(float(value["real"]), float(value["imag"]))
    elif isinstance(value, (list, tuple)) and len(value) == 2:
        result = complex(float(value[0]), float(value[1]))
    else:
        raise ValueError(f"{name} must be a real/imag mapping or pair")
    if not (math.isfinite(result.real) and math.isfinite(result.imag)):
        raise ValueError(f"{name} must be finite")
    return result


def _triplet(values: Any, name: str) -> list[complex]:
    if not isinstance(values, list) or len(values) != 3:
        raise ValueError(f"{name} must contain exactly three phasors")
    return [_phasor(value, f"{name}[{index}]") for index, value in enumerate(values)]


def _relative_spread(values: list[complex]) -> float:
    magnitudes = [abs(value) for value in values]
    scale = max(magnitudes, default=0.0)
    return (max(magnitudes) - min(magnitudes)) / scale if scale > 0.0 else math.inf


def _zero_sequence_residual(values: list[complex]) -> float:
    scale = sum(abs(value) for value in values)
    return abs(sum(values)) / scale if scale > 0.0 else math.inf


def _phase_step_error_deg(values: list[complex], expected_step_deg: float) -> float:
    errors = []
    for left, right in zip(values, values[1:]):
        if abs(left) == 0.0 or abs(right) == 0.0:
            return math.inf
        actual = math.degrees(cmath.phase(right / left))
        errors.append(abs((actual - expected_step_deg + 180.0) % 360.0 - 180.0))
    return max(errors, default=math.inf)


def _positive_resistance(value: Any, name: str) -> float:
    resistance = float(value)
    if not math.isfinite(resistance) or resistance <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return resistance


def three_phase_winding_power_balance_gate(
    summary: dict[str, Any],
    *,
    max_voltage_relative_spread: float = 1.0e-5,
    max_current_relative_spread: float = 5.0e-3,
    max_phase_step_error_deg: float = 1.0,
    max_star_kcl_residual: float = 1.0e-5,
    max_active_power_relative_residual: float = 1.0e-3,
) -> dict[str, Any]:
    """Gate a balanced source winding and copper loss in all coupled windings.

    The source-current direction is positive *into* the coupled field/circuit
    system.  ``source_winding`` contains three voltage and current phasors plus
    one phase resistance.  ``passive_windings`` contains zero or more current
    triplets and their phase resistances.  Peak and RMS phasors are both
    supported; their common power factor cancels in the relative residual.
    """
    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")
    convention = str(summary.get("phasor_convention") or "").strip().lower()
    if convention not in {"peak", "rms"}:
        raise ValueError("phasor_convention must be 'peak' or 'rms'")
    power_factor = 0.5 if convention == "peak" else 1.0
    expected_step = float(summary.get("expected_phase_step_deg", -120.0))
    if not math.isfinite(expected_step):
        raise ValueError("expected_phase_step_deg must be finite")

    tolerance_values = {
        "max_voltage_relative_spread": float(max_voltage_relative_spread),
        "max_current_relative_spread": float(max_current_relative_spread),
        "max_phase_step_error_deg": float(max_phase_step_error_deg),
        "max_star_kcl_residual": float(max_star_kcl_residual),
        "max_active_power_relative_residual": float(max_active_power_relative_residual),
    }
    if any(not math.isfinite(value) or value < 0.0 for value in tolerance_values.values()):
        raise ValueError("all tolerances must be finite and nonnegative")

    source = summary.get("source_winding")
    if not isinstance(source, dict):
        raise ValueError("source_winding must be a mapping")
    source_label = str(source.get("label") or "").strip()
    voltages = _triplet(source.get("voltage_phasors"), "source voltage_phasors")
    source_currents = _triplet(source.get("current_phasors"), "source current_phasors")
    source_resistance = _positive_resistance(
        source.get("phase_resistance_ohm"), "source phase_resistance_ohm"
    )

    passive = summary.get("passive_windings", [])
    if not isinstance(passive, list):
        raise ValueError("passive_windings must be a list")
    winding_rows = []
    copper_loss = source_resistance * sum(abs(value) ** 2 for value in source_currents)
    for index, row in enumerate(passive, start=1):
        if not isinstance(row, dict):
            raise ValueError(f"passive_windings[{index}] must be a mapping")
        label = str(row.get("label") or "").strip()
        currents = _triplet(row.get("current_phasors"), f"passive_windings[{index}] currents")
        resistance = _positive_resistance(
            row.get("phase_resistance_ohm"),
            f"passive_windings[{index}] phase_resistance_ohm",
        )
        loss = resistance * sum(abs(value) ** 2 for value in currents)
        copper_loss += loss
        winding_rows.append(
            {
                "label": label,
                "phase_resistance_ohm": resistance,
                "copper_loss_w": power_factor * loss,
                "current_magnitude_relative_spread": _relative_spread(currents),
            }
        )

    complex_input_power = power_factor * sum(
        voltage * current.conjugate()
        for voltage, current in zip(voltages, source_currents)
    )
    copper_loss *= power_factor
    power_scale = max(abs(complex_input_power.real), abs(copper_loss), 1.0e-300)
    active_power_residual = complex_input_power.real - copper_loss
    active_power_relative_residual = abs(active_power_residual) / power_scale
    voltage_spread = _relative_spread(voltages)
    current_spread = _relative_spread(source_currents)
    voltage_phase_error = _phase_step_error_deg(voltages, expected_step)
    current_phase_error = _phase_step_error_deg(source_currents, expected_step)
    voltage_zero = _zero_sequence_residual(voltages)
    current_zero = _zero_sequence_residual(source_currents)

    checks = {
        "source_label_recorded": bool(source_label),
        "voltage_unit_recorded": str(summary.get("voltage_unit") or "").strip() == "V",
        "current_unit_recorded": str(summary.get("current_unit") or "").strip() == "A",
        "resistance_unit_recorded": str(summary.get("resistance_unit") or "").strip() == "ohm",
        "source_voltage_magnitudes_balanced": voltage_spread
        <= tolerance_values["max_voltage_relative_spread"],
        "source_current_magnitudes_balanced": current_spread
        <= tolerance_values["max_current_relative_spread"],
        "source_voltage_phase_sequence": voltage_phase_error
        <= tolerance_values["max_phase_step_error_deg"],
        "source_current_phase_sequence": current_phase_error
        <= tolerance_values["max_phase_step_error_deg"],
        "source_voltage_zero_sequence_small": voltage_zero
        <= tolerance_values["max_star_kcl_residual"],
        "source_star_current_kcl_closes": current_zero
        <= tolerance_values["max_star_kcl_residual"],
        "passive_labels_recorded": all(row["label"] for row in winding_rows),
        "active_input_power_positive": complex_input_power.real > 0.0,
        "copper_loss_positive": copper_loss > 0.0,
        "active_power_closes_to_copper_loss": active_power_relative_residual
        <= tolerance_values["max_active_power_relative_residual"],
    }
    return {
        "policy": "three_phase_winding_power_balance_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "complex_input_power_va": {
                "real": complex_input_power.real,
                "imag": complex_input_power.imag,
            },
            "copper_loss_w": copper_loss,
            "active_power_residual_w": active_power_residual,
            "active_power_relative_residual": active_power_relative_residual,
            "source_voltage_magnitude_relative_spread": voltage_spread,
            "source_current_magnitude_relative_spread": current_spread,
            "source_voltage_zero_sequence_residual": voltage_zero,
            "source_current_zero_sequence_residual": current_zero,
            "max_source_voltage_phase_step_error_deg": voltage_phase_error,
            "max_source_current_phase_step_error_deg": current_phase_error,
            "passive_windings": winding_rows,
        },
        "phasor_convention": convention,
        "phasor_power_factor": power_factor,
        "expected_phase_step_deg": expected_step,
        "tolerances": tolerance_values,
        "notes": [
            "Normalize source currents positive into the coupled system before applying the gate.",
            "Apply the same peak/RMS convention to terminal power and every I^2R loss term.",
            "A passive winding need not satisfy zero-sequence KCL unless its terminals are explicitly connected in a star.",
        ],
    }
