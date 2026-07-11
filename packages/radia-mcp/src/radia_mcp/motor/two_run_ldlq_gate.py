"""Solver-independent two-run d/q inductance extraction gate."""
from __future__ import annotations

import json
import math
from collections.abc import Mapping, Sequence


def _numeric_vector(values: object, name: str) -> list[float]:
    if not isinstance(values, Sequence) or isinstance(values, (str, bytes)):
        raise ValueError(f"{name} must be an array")
    try:
        result = [float(value) for value in values]
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must contain numeric values") from exc
    if not result or not all(math.isfinite(value) for value in result):
        raise ValueError(f"{name} must contain finite values")
    return result


def _three_component_rows(values: object, name: str) -> list[tuple[float, float, float]]:
    if not isinstance(values, Sequence) or isinstance(values, (str, bytes)):
        raise ValueError(f"{name} must be an array")
    result: list[tuple[float, float, float]] = []
    for index, row in enumerate(values):
        if not isinstance(row, Sequence) or isinstance(row, (str, bytes)) or len(row) != 3:
            raise ValueError(f"{name}[{index}] must contain three components")
        try:
            triple = tuple(float(value) for value in row)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{name}[{index}] must be numeric") from exc
        if not all(math.isfinite(value) for value in triple):
            raise ValueError(f"{name}[{index}] must be finite")
        result.append(triple)  # type: ignore[arg-type]
    return result


def _park(rows: Sequence[tuple[float, float, float]], angles_deg: Sequence[float], pole_pairs: int) -> tuple[list[float], list[float]]:
    d_values: list[float] = []
    q_values: list[float] = []
    for (a, b, c), angle_deg in zip(rows, angles_deg):
        alpha = (2.0 / 3.0) * (a - 0.5 * b - 0.5 * c)
        beta = (2.0 / 3.0) * ((math.sqrt(3.0) / 2.0) * (b - c))
        theta = math.radians(pole_pairs * angle_deg)
        d_values.append(alpha * math.cos(theta) + beta * math.sin(theta))
        q_values.append(-alpha * math.sin(theta) + beta * math.cos(theta))
    return d_values, q_values


def _mean(values: Sequence[float]) -> float:
    return sum(values) / len(values)


def _relative_pkpk(values: Sequence[float], scale: float) -> float:
    return (max(values) - min(values)) / max(abs(scale), 1.0e-30)


def evaluate_ipm_two_run_ldlq(
    summary: Mapping[str, object],
    angle_tolerance_deg: float = 1.0e-9,
    current_balance_tolerance_a: float = 1.0e-8,
    dq_current_ripple_relative_tolerance: float = 1.0e-6,
    pm_q_relative_tolerance: float = 3.0e-2,
    incremental_flux_ripple_relative_tolerance: float = 5.0e-2,
) -> dict:
    """Validate PM-only subtraction and extract two-run ``Ld`` and ``Lq``.

    The PM-only and current-on runs must expose their angle grids and phase
    order independently.  Subtracting their phase fluxes before the Park
    transform removes the permanent-magnet contribution.  The resulting
    finite-current values are secant perturbation inductances unless a
    frozen-permeability or small-signal workflow is documented separately.
    """
    if not isinstance(summary, Mapping):
        raise ValueError("summary must be an object")
    if any(value <= 0.0 for value in (
        angle_tolerance_deg,
        current_balance_tolerance_a,
        dq_current_ripple_relative_tolerance,
        pm_q_relative_tolerance,
        incremental_flux_ripple_relative_tolerance,
    )):
        raise ValueError("tolerances must be positive")
    try:
        pole_pairs = int(summary["pole_pairs"])
        pm_run = summary["pm_only"]
        on_run = summary["current_on"]
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("pole_pairs, pm_only and current_on are required") from exc
    if pole_pairs <= 0:
        raise ValueError("pole_pairs must be positive")
    if not isinstance(pm_run, Mapping) or not isinstance(on_run, Mapping):
        raise ValueError("pm_only and current_on must be objects")

    pm_angles = _numeric_vector(pm_run.get("mechanical_angles_deg"), "pm_only.mechanical_angles_deg")
    on_angles = _numeric_vector(on_run.get("mechanical_angles_deg"), "current_on.mechanical_angles_deg")
    pm_flux = _three_component_rows(pm_run.get("phase_flux_wb"), "pm_only.phase_flux_wb")
    on_flux = _three_component_rows(on_run.get("phase_flux_wb"), "current_on.phase_flux_wb")
    currents = _three_component_rows(on_run.get("phase_currents_a"), "current_on.phase_currents_a")
    lengths = {len(pm_angles), len(on_angles), len(pm_flux), len(on_flux), len(currents)}
    if len(lengths) != 1 or len(pm_angles) < 12:
        raise ValueError("both runs need matching arrays with at least 12 samples")

    pm_order = pm_run.get("phase_order")
    on_order = on_run.get("phase_order")
    angle_error = max(abs(a - b) for a, b in zip(pm_angles, on_angles))
    monotonic = all(b > a for a, b in zip(pm_angles, pm_angles[1:]))
    electrical_span = pole_pairs * (pm_angles[-1] - pm_angles[0])

    current_d, current_q = _park(currents, on_angles, pole_pairs)
    pm_d, pm_q = _park(pm_flux, pm_angles, pole_pairs)
    delta_flux = [tuple(on - pm for on, pm in zip(on_row, pm_row)) for on_row, pm_row in zip(on_flux, pm_flux)]
    delta_d, delta_q = _park(delta_flux, on_angles, pole_pairs)
    id_mean = _mean(current_d)
    iq_mean = _mean(current_q)
    delta_d_mean = _mean(delta_d)
    delta_q_mean = _mean(delta_q)
    pm_d_mean = _mean(pm_d)
    ld_h = delta_d_mean / id_mean if abs(id_mean) > 1.0e-30 else math.nan
    lq_h = delta_q_mean / iq_mean if abs(iq_mean) > 1.0e-30 else math.nan
    current_balance_max = max(abs(sum(row)) for row in currents)
    id_ripple = _relative_pkpk(current_d, id_mean)
    iq_ripple = _relative_pkpk(current_q, iq_mean)
    pm_q_max_relative = max(abs(value) for value in pm_q) / max(abs(pm_d_mean), 1.0e-30)
    delta_d_ripple = _relative_pkpk(delta_d, delta_d_mean)
    delta_q_ripple = _relative_pkpk(delta_q, delta_q_mean)

    expected_saliency = str(summary.get("expected_saliency", "none"))
    if expected_saliency not in {"none", "Lq_gt_Ld", "Ld_gt_Lq"}:
        raise ValueError("expected_saliency must be none, Lq_gt_Ld or Ld_gt_Lq")
    saliency_ok = (
        expected_saliency == "none"
        or (expected_saliency == "Lq_gt_Ld" and lq_h > ld_h)
        or (expected_saliency == "Ld_gt_Lq" and ld_h > lq_h)
    )
    checks = {
        "canonical_phase_order": pm_order == ["A", "B", "C"] and on_order == ["A", "B", "C"],
        "same_angle_grid": angle_error <= angle_tolerance_deg,
        "mechanical_angle_monotonic": monotonic,
        "complete_electrical_cycle": electrical_span >= 360.0 - angle_tolerance_deg,
        "phase_current_balance": current_balance_max <= current_balance_tolerance_a,
        "nonzero_dq_current": abs(id_mean) > 1.0e-30 and abs(iq_mean) > 1.0e-30,
        "dq_current_constant": id_ripple <= dq_current_ripple_relative_tolerance and iq_ripple <= dq_current_ripple_relative_tolerance,
        "pm_flux_d_axis_aligned": abs(pm_d_mean) > 1.0e-30 and pm_q_max_relative <= pm_q_relative_tolerance,
        "pm_subtraction_nonzero": max(abs(value) for value in (*delta_d, *delta_q)) > 1.0e-30,
        "incremental_flux_nearly_constant": delta_d_ripple <= incremental_flux_ripple_relative_tolerance and delta_q_ripple <= incremental_flux_ripple_relative_tolerance,
        "positive_finite_ld_lq": math.isfinite(ld_h) and math.isfinite(lq_h) and ld_h > 0.0 and lq_h > 0.0,
        "expected_saliency": saliency_ok,
    }
    return {
        "schema": "radia-motor-ipm-two-run-ldlq/v1",
        "policy": "same_angle_same_phase_pm_only_subtraction_before_park",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "sample_count": len(pm_angles),
        "pole_pairs": pole_pairs,
        "electrical_span_deg": electrical_span,
        "angle_grid_max_error_deg": angle_error,
        "current_balance_max_abs_a": current_balance_max,
        "current_d_mean_a": id_mean,
        "current_q_mean_a": iq_mean,
        "pm_flux_d_mean_wb": pm_d_mean,
        "pm_flux_q_max_relative": pm_q_max_relative,
        "delta_flux_d_mean_wb": delta_d_mean,
        "delta_flux_q_mean_wb": delta_q_mean,
        "ld_h": ld_h,
        "lq_h": lq_h,
        "saliency_ratio_lq_over_ld": lq_h / ld_h if ld_h > 0.0 else math.nan,
        "ripple_relative": {
            "current_d": id_ripple,
            "current_q": iq_ripple,
            "delta_flux_d": delta_d_ripple,
            "delta_flux_q": delta_q_ripple,
        },
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "interpretation": "finite-current secant perturbation; use frozen permeability or a small-signal pair for differential inductance",
    }


def ipm_two_run_ldlq_gate(summary_json: str) -> str:
    try:
        summary = json.loads(summary_json)
    except json.JSONDecodeError as exc:
        raise ValueError(f"summary_json must be valid JSON: {exc.msg}") from exc
    return json.dumps(evaluate_ipm_two_run_ldlq(summary), indent=2, sort_keys=True)
