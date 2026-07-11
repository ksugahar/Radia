"""Solver-neutral gates for linear-induction frequency sweeps."""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping


def linear_induction_frequency_sweep_gate(
    rows: Iterable[Mapping[str, object]],
    *,
    thrust_abs_tol_n: float = 0.75,
    thrust_rel_tol: float = 2.0e-3,
    phase_balance_atol_a: float = 1.0e-9,
) -> dict[str, object]:
    """Gate thrust, loss, power, and force-method agreement over frequency."""

    data = list(rows)
    if len(data) < 5 or any(not isinstance(row, Mapping) for row in data):
        raise ValueError("rows must contain at least five mappings")
    tolerances = (thrust_abs_tol_n, thrust_rel_tol, phase_balance_atol_a)
    if any(not math.isfinite(float(value)) or float(value) < 0.0 for value in tolerances):
        raise ValueError("tolerances must be finite and non-negative")

    def pair_magnitude(value: object, field: str) -> float:
        if not isinstance(value, Mapping):
            raise ValueError(f"{field} must be a real/imag mapping")
        real = float(value.get("real"))
        imag = float(value.get("imag"))
        if not math.isfinite(real) or not math.isfinite(imag):
            raise ValueError(f"{field} must be finite")
        return math.hypot(real, imag)

    try:
        frequencies = [float(row["frequency_hz"]) for row in data]
        lorentz = [float(row["lorentz_thrust_n"]) for row in data]
        stress = [float(row["weighted_stress_thrust_n"]) for row in data]
        losses = [float(row["plate_resistive_loss_w"]) for row in data]
        phase_imbalance = [pair_magnitude(row["phase_current_sum_a"], "phase_current_sum_a") for row in data]
        real_power = [float(row["three_phase_complex_power_va"]["real"]) for row in data]
        nodes = [int(row["node_count"]) for row in data]
        elements = [int(row["element_count"]) for row in data]
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"invalid linear-induction row: {exc}") from exc
    scalars = frequencies + lorentz + stress + losses + real_power
    if not all(math.isfinite(value) for value in scalars):
        raise ValueError("frequency-sweep scalars must be finite")

    method_errors = [abs(a - b) for a, b in zip(lorentz, stress)]
    method_limits = [
        float(thrust_abs_tol_n) + float(thrust_rel_tol) * max(abs(a), abs(b))
        for a, b in zip(lorentz, stress)
    ]
    peak_index = max(range(len(lorentz)), key=lambda index: lorentz[index])
    peak_thrust = lorentz[peak_index]
    checks = {
        "frequency_strictly_increasing_positive": all(
            value > 0.0 for value in frequencies
        ) and all(a < b for a, b in zip(frequencies, frequencies[1:])),
        "three_phase_currents_balanced": max(phase_imbalance) <= float(phase_balance_atol_a),
        "mesh_inventory_invariant": len(set(nodes)) == len(set(elements)) == 1,
        "lorentz_and_weighted_thrust_agree": all(
            error <= limit for error, limit in zip(method_errors, method_limits)
        ),
        "plate_loss_positive_increasing": all(value > 0.0 for value in losses)
        and all(a < b for a, b in zip(losses, losses[1:])),
        "real_input_power_positive": all(value > 0.0 for value in real_power),
        "plate_loss_below_real_input_power": all(loss < power for loss, power in zip(losses, real_power)),
        "positive_thrust_peak_is_interior": peak_thrust > 0.0 and 0 < peak_index < len(data) - 1,
        "thrust_reversal_observed": any(value > 0.0 for value in lorentz) and any(value < 0.0 for value in lorentz),
        "high_frequency_thrust_collapses": abs(lorentz[-1]) <= 0.01 * peak_thrust,
    }
    return {
        "policy": "linear_induction_frequency_sweep_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "peak_frequency_hz": frequencies[peak_index],
            "peak_lorentz_thrust_n": peak_thrust,
            "maximum_thrust_method_absolute_difference_n": max(method_errors),
            "maximum_phase_current_imbalance_a": max(phase_imbalance),
            "loss_growth_ratio": losses[-1] / losses[0],
            "node_count": nodes[0],
            "element_count": elements[0],
        },
        "lesson": (
            "Compare Lorentz and weighted-stress force only along the intended travel direction. "
            "Use an absolute-plus-relative tolerance near thrust zero crossings, and retain phase balance, "
            "loss, input-power, and mesh-inventory evidence with the frequency response."
        ),
    }
