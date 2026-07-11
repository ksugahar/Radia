"""Causality and geometric-arrival gate for CQ scattering responses."""
from __future__ import annotations

import math


def cq_scattering_arrival_gate(
    *,
    time_step_s: float,
    geometric_arrival_s: float,
    measured_peak_s: float,
    max_relative_residual: float,
    finite_response: bool,
    real_time_response: bool,
    max_peak_lag_steps: float = 3.0,
    max_residual: float = 1.0e-6,
) -> dict:
    """Require a real, converged scattered peak in its causal ray window."""

    values = [time_step_s, geometric_arrival_s, measured_peak_s, max_relative_residual, max_peak_lag_steps, max_residual]
    if not all(math.isfinite(float(value)) for value in values):
        raise ValueError("all scalar inputs must be finite")
    if time_step_s <= 0.0 or geometric_arrival_s < 0.0 or measured_peak_s < 0.0:
        raise ValueError("time_step_s must be positive and arrival times nonnegative")
    if max_peak_lag_steps < 0.0 or max_residual < 0.0 or max_relative_residual < 0.0:
        raise ValueError("tolerances and residuals must be nonnegative")

    peak_offset_steps = (measured_peak_s - geometric_arrival_s) / time_step_s
    checks = {
        "response_finite": finite_response is True,
        "inverse_transform_real": real_time_response is True,
        "linear_solves_converged": max_relative_residual <= max_residual,
        "peak_not_acausal": peak_offset_steps >= -1.0,
        "peak_near_geometric_arrival": peak_offset_steps <= max_peak_lag_steps,
    }
    return {
        "policy": "cq_scattering_arrival_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {"peak_offset_steps": peak_offset_steps, "max_relative_residual": max_relative_residual},
        "notes": [
            "gate the scattered-field peak, not a low-amplitude pulse tail, against the geometric ray arrival",
            "one time step of lead is allowed for discrete peak sampling; the declared lag window handles extended scatterers",
        ],
    }
