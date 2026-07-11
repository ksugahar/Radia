"""Three-phase flux-linkage Park-axis alignment gate."""
from __future__ import annotations

import json
import math
from collections.abc import Sequence


def evaluate_phase_flux_park_alignment(
    mechanical_angles_deg: Sequence[float],
    phase_flux_wb: Sequence[Sequence[float]],
    pole_pairs: int,
    q_relative_tolerance: float = 3.0e-2,
    d_ripple_relative_tolerance: float = 2.0e-2,
) -> dict:
    """Project a PM-only abc flux sweep and gate its d/q alignment.

    A complete electrical cycle should become approximately constant ``d`` and
    near-zero ``q`` when the recorded mechanical angle is multiplied by the
    explicit pole-pair count. The transform is amplitude invariant.
    """
    if pole_pairs <= 0:
        raise ValueError("pole_pairs must be positive")
    if q_relative_tolerance <= 0 or d_ripple_relative_tolerance <= 0:
        raise ValueError("relative tolerances must be positive")
    try:
        angles = [float(value) for value in mechanical_angles_deg]
    except (TypeError, ValueError) as exc:
        raise ValueError("mechanical_angles_deg must be numeric") from exc
    if len(angles) < 12 or len(angles) != len(phase_flux_wb):
        raise ValueError("angles and phase_flux_wb need the same length of at least 12")
    if not all(math.isfinite(value) for value in angles):
        raise ValueError("angles must be finite")

    d_values: list[float] = []
    q_values: list[float] = []
    for index, row in enumerate(phase_flux_wb):
        if not isinstance(row, Sequence) or isinstance(row, (str, bytes)) or len(row) != 3:
            raise ValueError(f"phase_flux_wb[{index}] must contain three components")
        try:
            a, b, c = (float(value) for value in row)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"phase_flux_wb[{index}] must be numeric") from exc
        if not all(math.isfinite(value) for value in (a, b, c)):
            raise ValueError(f"phase_flux_wb[{index}] must be finite")
        alpha = (2.0 / 3.0) * (a - 0.5 * b - 0.5 * c)
        beta = (2.0 / 3.0) * ((math.sqrt(3.0) / 2.0) * (b - c))
        theta = math.radians(pole_pairs * angles[index])
        d_values.append(alpha * math.cos(theta) + beta * math.sin(theta))
        q_values.append(-alpha * math.sin(theta) + beta * math.cos(theta))

    d_mean = sum(d_values) / len(d_values)
    q_mean = sum(q_values) / len(q_values)
    scale = max(abs(d_mean), 1.0e-30)
    q_abs_max_relative = max(abs(value) for value in q_values) / scale
    d_pkpk_relative = (max(d_values) - min(d_values)) / scale
    electrical_span_deg = pole_pairs * (max(angles) - min(angles))
    checks = {
        "mechanical_angle_monotonic": all(b > a for a, b in zip(angles, angles[1:])),
        "complete_electrical_cycle": electrical_span_deg >= 360.0 - 1.0e-9,
        "d_axis_nonzero": abs(d_mean) > 1.0e-30,
        "q_axis_near_zero": q_abs_max_relative <= q_relative_tolerance,
        "d_axis_nearly_constant": d_pkpk_relative <= d_ripple_relative_tolerance,
    }
    return {
        "schema": "radia-motor-phase-flux-park-alignment/v1",
        "policy": "pm_only_phase_flux_requires_explicit_angle_basis_and_dq_alignment",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "sample_count": len(angles),
        "pole_pairs": pole_pairs,
        "electrical_span_deg": electrical_span_deg,
        "d_mean_Wb": d_mean,
        "q_mean_Wb": q_mean,
        "q_abs_max_relative": q_abs_max_relative,
        "d_pkpk_relative": d_pkpk_relative,
        "q_relative_tolerance": q_relative_tolerance,
        "d_ripple_relative_tolerance": d_ripple_relative_tolerance,
        "checks": checks,
    }


def phase_flux_park_alignment_gate(
    mechanical_angles_deg_json: str,
    phase_flux_wb_json: str,
    pole_pairs: int,
    q_relative_tolerance: float = 3.0e-2,
    d_ripple_relative_tolerance: float = 2.0e-2,
) -> str:
    try:
        angles = json.loads(mechanical_angles_deg_json)
        flux = json.loads(phase_flux_wb_json)
    except json.JSONDecodeError as exc:
        raise ValueError(f"angle and flux inputs must be valid JSON: {exc.msg}") from exc
    return json.dumps(
        evaluate_phase_flux_park_alignment(
            angles,
            flux,
            pole_pairs,
            q_relative_tolerance,
            d_ripple_relative_tolerance,
        ),
        indent=2,
        sort_keys=True,
    )
