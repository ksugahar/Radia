"""Solver-independent power-wave closure gate for a passive one-port sweep."""
from __future__ import annotations

import json
import math


def one_port_power_balance_gate(
    summary_json: str,
    max_power_relative_residual: float = 1.0e-9,
    max_balance_abs_residual: float = 1.0e-9,
    max_reference_impedance_relative_drift: float = 1.0e-9,
) -> dict:
    """Check accepted power against ``Pstim * (1 - |S11|**2)``."""

    tolerances = (
        max_power_relative_residual,
        max_balance_abs_residual,
        max_reference_impedance_relative_drift,
    )
    if any(value < 0.0 for value in tolerances):
        raise ValueError("tolerances must be nonnegative")
    try:
        summary = json.loads(summary_json)
    except json.JSONDecodeError as exc:
        raise ValueError(f"summary_json must be valid JSON: {exc.msg}") from exc
    if not isinstance(summary, dict):
        raise ValueError("summary_json must decode to an object")
    rows = summary.get("rows")
    if not isinstance(rows, list) or len(rows) < 2:
        raise ValueError("rows must contain at least two frequency samples")

    frequencies: list[float] = []
    power_residuals: list[float] = []
    balance_residuals: list[float] = []
    zrefs: list[complex] = []
    passive_rows: list[bool] = []
    accepted_bounds: list[bool] = []
    for index, row in enumerate(rows):
        if not isinstance(row, dict):
            raise ValueError(f"rows[{index}] must be an object")
        values = {
            "frequency": float(row["frequency"]),
            "s11_real": float(row["s11_real"]),
            "s11_imag": float(row["s11_imag"]),
            "stimulated_power_w": float(row["stimulated_power_w"]),
            "accepted_power_w": float(row["accepted_power_w"]),
            "balance_magnitude": float(row["balance_magnitude"]),
            "zref_real_ohm": float(row["zref_real_ohm"]),
            "zref_imag_ohm": float(row["zref_imag_ohm"]),
        }
        if not all(math.isfinite(value) for value in values.values()):
            raise ValueError("all row values must be finite")
        frequency = values["frequency"]
        s11 = complex(values["s11_real"], values["s11_imag"])
        stimulated = values["stimulated_power_w"]
        accepted = values["accepted_power_w"]
        if stimulated <= 0.0:
            raise ValueError("stimulated_power_w must be positive")
        expected_accepted = stimulated * (1.0 - abs(s11) ** 2)
        frequencies.append(frequency)
        power_residuals.append(abs(accepted - expected_accepted) / stimulated)
        balance_residuals.append(abs(values["balance_magnitude"] - abs(s11)))
        zrefs.append(complex(values["zref_real_ohm"], values["zref_imag_ohm"]))
        passive_rows.append(abs(s11) <= 1.0 + max_power_relative_residual)
        accepted_bounds.append(-max_power_relative_residual * stimulated <= accepted <= stimulated * (1.0 + max_power_relative_residual))

    zref_scale = max(abs(zrefs[0]), 1.0e-300)
    zref_drift = max(abs(value - zrefs[0]) for value in zrefs) / zref_scale
    checks = {
        "frequency_unit_recorded": bool(str(summary.get("frequency_unit", "")).strip()),
        "power_unit_is_W": summary.get("power_unit") == "W",
        "reference_impedance_unit_is_ohm": summary.get("reference_impedance_unit") == "ohm",
        "sparameter_basis_is_power_wave": summary.get("sparameter_basis") == "power_wave",
        "frequencies_strictly_increase": all(a < b for a, b in zip(frequencies, frequencies[1:])),
        "all_rows_passive": all(passive_rows),
        "accepted_power_is_bounded": all(accepted_bounds),
        "accepted_power_closes_from_s11": max(power_residuals) <= max_power_relative_residual,
        "reported_balance_matches_s11": max(balance_residuals) <= max_balance_abs_residual,
        "reference_impedance_real_positive": zrefs[0].real > 0.0,
        "reference_impedance_is_real": max(abs(value.imag) for value in zrefs) <= max_reference_impedance_relative_drift * zref_scale,
        "reference_impedance_is_stable": zref_drift <= max_reference_impedance_relative_drift,
    }
    return {
        "policy": "one_port_power_balance_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "row_count": len(rows),
            "frequency_min": min(frequencies),
            "frequency_max": max(frequencies),
            "maximum_s11_magnitude": max(
                abs(complex(float(row["s11_real"]), float(row["s11_imag"]))) for row in rows
            ),
            "maximum_power_relative_residual": max(power_residuals),
            "maximum_balance_abs_residual": max(balance_residuals),
            "reference_impedance_ohm": {"real": zrefs[0].real, "imag": zrefs[0].imag},
            "maximum_reference_impedance_relative_drift": zref_drift,
        },
        "notes": [
            "accepted power must be compared on the same sparse frequencies as the power result rows",
            "a dense S11 sweep may be interpolated only when the interpolation method and row identity are recorded",
        ],
    }
