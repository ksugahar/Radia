"""Solver-neutral gates for PWM current control and harmonic-loss result tables."""

from __future__ import annotations

import math
from typing import Any


def _vector(value: Any, name: str, *, minimum: int = 1) -> list[float]:
    try:
        parsed = [float(item) for item in value]
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a numeric sequence") from exc
    if len(parsed) < minimum or not all(math.isfinite(item) for item in parsed):
        raise ValueError(f"{name} must contain at least {minimum} finite values")
    return parsed


def _matrix(
    value: Any,
    name: str,
    *,
    rows: int,
    columns: int | None = None,
) -> list[list[float]]:
    try:
        parsed = [[float(item) for item in row] for row in value]
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a numeric matrix") from exc
    if len(parsed) != rows or not parsed:
        raise ValueError(f"{name} must have {rows} rows")
    width = len(parsed[0])
    if width == 0 or (columns is not None and width != columns):
        raise ValueError(f"{name} has an invalid column count")
    if any(len(row) != width for row in parsed):
        raise ValueError(f"{name} rows must have equal width")
    if not all(math.isfinite(item) for row in parsed for item in row):
        raise ValueError(f"{name} must contain finite values")
    return parsed


def _relative(residual: float, scale: float) -> float:
    return abs(residual) / max(abs(scale), 1.0e-30)


def pwm_controlled_motor_loss_gate(
    payload: dict[str, Any],
    *,
    max_three_phase_kcl_relative_error: float = 1.0e-3,
    max_tail_control_tracking_rms_relative_error: float = 0.05,
    max_angle_speed_integral_relative_error: float = 2.0e-4,
    max_power_sum_relative_error: float = 1.0e-10,
    max_loss_identity_relative_error: float = 1.0e-10,
    max_frequency_step_relative_span: float = 1.0e-10,
) -> dict[str, Any]:
    """Validate a PWM-controlled motor time series and harmonic-loss package.

    The loss spectrum uses a common table convention: row zero stores aggregate
    values, while positive-frequency rows store harmonic bins. Eddy loss is
    reconstructable from those bins; hysteresis and combined iron loss may be
    aggregate-only. Ratio-derived R/L values are excluded when current is near
    zero because they are ill-conditioned rather than physical outliers.
    """

    if not isinstance(payload, dict):
        raise ValueError("payload must be an object")
    tolerances = [
        max_three_phase_kcl_relative_error,
        max_tail_control_tracking_rms_relative_error,
        max_angle_speed_integral_relative_error,
        max_power_sum_relative_error,
        max_loss_identity_relative_error,
        max_frequency_step_relative_span,
    ]
    if any(not math.isfinite(float(value)) or float(value) < 0.0 for value in tolerances):
        raise ValueError("tolerances must be finite and nonnegative")

    time_series = payload.get("time_series")
    loss_spectrum = payload.get("loss_spectrum")
    if not isinstance(time_series, dict) or not isinstance(loss_spectrum, dict):
        raise ValueError("time_series and loss_spectrum objects are required")

    identity_value = payload.get("artifact_identity")
    identity_present = isinstance(identity_value, dict)
    cycle_generation_ok = True
    restart_phase_origin_ok = True
    if identity_value is not None and not identity_present:
        cycle_generation_ok = False
        restart_phase_origin_ok = False
    elif identity_present:
        torque_generation = str(identity_value.get("torque_cycle_generation", ""))
        loss_generation = str(identity_value.get("loss_cycle_generation", ""))
        cycle_generation_ok = (
            bool(torque_generation) and torque_generation == loss_generation
        )
        segments = identity_value.get("waveform_segments")
        if not isinstance(segments, list) or not segments:
            restart_phase_origin_ok = False
        else:
            phase_origins = []
            previous_end = -math.inf
            for segment in segments:
                if not isinstance(segment, dict):
                    restart_phase_origin_ok = False
                    break
                try:
                    phase_origin = float(segment["phase_origin_deg"])
                    start = float(segment["start_time_s"])
                    end = float(segment["end_time_s"])
                except (KeyError, TypeError, ValueError):
                    restart_phase_origin_ok = False
                    break
                if not all(math.isfinite(value) for value in (phase_origin, start, end)):
                    restart_phase_origin_ok = False
                    break
                if not str(segment.get("segment_id", "")) or end < start or start < previous_end:
                    restart_phase_origin_ok = False
                    break
                phase_origins.append(phase_origin)
                previous_end = end
            restart_phase_origin_ok = (
                restart_phase_origin_ok
                and bool(phase_origins)
                and len(set(phase_origins)) == 1
            )

    time_s = _vector(time_series.get("time_s"), "time_series.time_s", minimum=5)
    count = len(time_s)
    angle_deg = _vector(time_series.get("angle_deg"), "time_series.angle_deg", minimum=count)
    speed_rpm = _vector(time_series.get("speed_rpm"), "time_series.speed_rpm", minimum=count)
    torque_nm = _vector(time_series.get("torque_nm"), "time_series.torque_nm", minimum=count)
    total_power = _vector(
        time_series.get("reported_total_power_w"),
        "time_series.reported_total_power_w",
        minimum=count,
    )
    if any(len(values) != count for values in (angle_deg, speed_rpm, torque_nm, total_power)):
        raise ValueError("time-series vectors must have equal lengths")
    phase_currents = _matrix(
        time_series.get("phase_currents_a"),
        "time_series.phase_currents_a",
        rows=count,
        columns=3,
    )
    commands = _matrix(
        time_series.get("control_command_a"),
        "time_series.control_command_a",
        rows=count,
        columns=2,
    )
    feedback = _matrix(
        time_series.get("control_feedback_a"),
        "time_series.control_feedback_a",
        rows=count,
        columns=2,
    )
    power_components = _matrix(
        time_series.get("power_components_w"),
        "time_series.power_components_w",
        rows=count,
    )
    component_time_value = time_series.get("component_time_s")
    component_time_axes: list[list[float]] | None = None
    component_time_alignment_error = 0.0
    if component_time_value is not None:
        component_time_axes = _matrix(
            component_time_value,
            "time_series.component_time_s",
            rows=len(power_components[0]),
            columns=count,
        )
        component_time_alignment_error = max(
            abs(component_time_axes[component][index] - time_s[index])
            for component in range(len(component_time_axes))
            for index in range(count)
        )

    frequencies = _vector(
        loss_spectrum.get("frequency_hz"),
        "loss_spectrum.frequency_hz",
        minimum=3,
    )
    loss_count = len(frequencies)
    eddy = _matrix(
        loss_spectrum.get("eddy_components_w"),
        "loss_spectrum.eddy_components_w",
        rows=loss_count,
    )
    hysteresis = _matrix(
        loss_spectrum.get("hysteresis_components_w"),
        "loss_spectrum.hysteresis_components_w",
        rows=loss_count,
        columns=len(eddy[0]),
    )
    iron = _matrix(
        loss_spectrum.get("iron_components_w"),
        "loss_spectrum.iron_components_w",
        rows=loss_count,
        columns=len(eddy[0]),
    )
    reported_eddy = _vector(
        loss_spectrum.get("reported_eddy_total_w"),
        "loss_spectrum.reported_eddy_total_w",
        minimum=loss_count,
    )
    reported_hysteresis = _vector(
        loss_spectrum.get("reported_hysteresis_total_w"),
        "loss_spectrum.reported_hysteresis_total_w",
        minimum=loss_count,
    )
    reported_iron = _vector(
        loss_spectrum.get("reported_iron_total_w"),
        "loss_spectrum.reported_iron_total_w",
        minimum=loss_count,
    )
    if any(len(values) != loss_count for values in (reported_eddy, reported_hysteresis, reported_iron)):
        raise ValueError("loss total vectors must match the frequency axis")

    time_deltas = [right - left for left, right in zip(time_s, time_s[1:])]
    frequency_deltas = [right - left for left, right in zip(frequencies, frequencies[1:])]
    phase_scale = max(abs(value) for row in phase_currents for value in row)
    kcl_relative = max(abs(sum(row)) for row in phase_currents) / max(phase_scale, 1.0e-30)

    tail_start = max(1, int(0.75 * count))
    tracking_errors = []
    for component in range(2):
        command_rms = math.sqrt(
            sum(commands[index][component] ** 2 for index in range(tail_start, count))
            / (count - tail_start)
        )
        error_rms = math.sqrt(
            sum(
                (feedback[index][component] - commands[index][component]) ** 2
                for index in range(tail_start, count)
            )
            / (count - tail_start)
        )
        tracking_errors.append(error_rms / max(command_rms, 1.0e-30))

    integrated_angle = [angle_deg[0]]
    for index, delta_t in enumerate(time_deltas):
        integrated_angle.append(
            integrated_angle[-1]
            + 3.0 * (speed_rpm[index] + speed_rpm[index + 1]) * delta_t
        )
    angle_span = max(angle_deg) - min(angle_deg)
    angle_integral_relative = max(
        abs(actual - expected) for actual, expected in zip(angle_deg, integrated_angle)
    ) / max(angle_span, 1.0e-30)

    power_scale = max(abs(value) for value in total_power)
    power_sum_relative = max(
        abs(sum(components) - reported)
        for components, reported in zip(power_components, total_power)
    ) / max(power_scale, 1.0e-30)

    loss_scale = max(
        abs(value)
        for matrix in (eddy, hysteresis, iron)
        for row in matrix
        for value in row
    )
    loss_total_residual = max(
        abs(sum(matrix[index]) - reported[index])
        for matrix, reported in (
            (eddy, reported_eddy),
            (hysteresis, reported_hysteresis),
            (iron, reported_iron),
        )
        for index in range(loss_count)
    ) / max(loss_scale, 1.0e-30)
    eddy_reconstruction = max(
        [
            abs(eddy[0][component] - sum(row[component] for row in eddy[1:]))
            for component in range(len(eddy[0]))
        ]
        + [abs(reported_eddy[0] - sum(reported_eddy[1:]))]
    ) / max(loss_scale, 1.0e-30)
    aggregate_iron_decomposition = max(
        [
            abs(iron[0][component] - eddy[0][component] - hysteresis[0][component])
            for component in range(len(eddy[0]))
        ]
        + [abs(reported_iron[0] - reported_eddy[0] - reported_hysteresis[0])]
    ) / max(loss_scale, 1.0e-30)
    aggregate_only_bin_residual = max(
        [abs(value) for row in hysteresis[1:] for value in row]
        + [abs(value) for row in iron[1:] for value in row]
        + [abs(value) for value in reported_hysteresis[1:]]
        + [abs(value) for value in reported_iron[1:]]
    ) / max(loss_scale, 1.0e-30)
    frequency_step_span = (
        (max(frequency_deltas) - min(frequency_deltas))
        / max(abs(sum(frequency_deltas) / len(frequency_deltas)), 1.0e-30)
    )

    checks = {
        "time_axis_strictly_increases": all(delta > 0.0 for delta in time_deltas),
        "frequency_axis_has_zero_summary_then_uniform_positive_bins": frequencies[0] == 0.0
        and all(delta > 0.0 for delta in frequency_deltas)
        and frequency_step_span <= max_frequency_step_relative_span,
        "three_phase_current_satisfies_kcl": kcl_relative
        <= max_three_phase_kcl_relative_error,
        "tail_current_control_tracks_commands": max(tracking_errors)
        <= max_tail_control_tracking_rms_relative_error,
        "angle_matches_integrated_speed": angle_integral_relative
        <= max_angle_speed_integral_relative_error,
        "circuit_power_total_matches_component_sum": power_sum_relative
        <= max_power_sum_relative_error,
        "power_component_time_axes_match_common_axis_knotwise": (
            component_time_axes is None or component_time_alignment_error <= 1.0e-12
        ),
        "loss_total_columns_match_component_sums": loss_total_residual
        <= max_loss_identity_relative_error,
        "eddy_aggregate_matches_harmonic_bin_sum": eddy_reconstruction
        <= max_loss_identity_relative_error,
        "aggregate_iron_equals_eddy_plus_hysteresis": aggregate_iron_decomposition
        <= max_loss_identity_relative_error,
        "hysteresis_and_iron_are_aggregate_only": aggregate_only_bin_residual
        <= max_loss_identity_relative_error,
        "near_zero_current_ratio_outputs_are_diagnostic_only": payload.get(
            "ratio_diagnostic_policy"
        )
        == "exclude_ratio_when_denominator_current_is_below_floor",
        "torque_and_loss_share_periodic_cycle_generation": cycle_generation_ok,
        "restart_segments_preserve_phase_origin": restart_phase_origin_ok,
    }
    tail_torque = torque_nm[tail_start:]
    tail_torque_mean = sum(tail_torque) / len(tail_torque)
    tail_torque_ripple = (max(tail_torque) - min(tail_torque)) / max(
        abs(tail_torque_mean), 1.0e-30
    )
    return {
        "policy": "pwm_controlled_motor_loss_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "warnings": [] if identity_present else ["artifact_identity_not_recorded"],
        "metrics": {
            "time_row_count": count,
            "loss_row_count": loss_count,
            "loss_component_count": len(eddy[0]),
            "three_phase_kcl_global_relative_error": kcl_relative,
            "tail_control_tracking_rms_relative_errors": tracking_errors,
            "angle_speed_integral_relative_error": angle_integral_relative,
            "circuit_power_sum_global_relative_error": power_sum_relative,
            "power_component_time_axis_maximum_absolute_error_s": component_time_alignment_error,
            "loss_total_column_relative_error": loss_total_residual,
            "eddy_harmonic_reconstruction_relative_error": eddy_reconstruction,
            "aggregate_iron_decomposition_relative_error": aggregate_iron_decomposition,
            "aggregate_only_bin_relative_residual": aggregate_only_bin_residual,
            "frequency_step_relative_span": frequency_step_span,
            "tail_torque_mean_nm": tail_torque_mean,
            "tail_torque_peak_to_peak_relative": tail_torque_ripple,
        },
        "lesson": (
            "Gate PWM motor results with three-phase KCL, tail current-command tracking, "
            "integrated speed/angle, and power-column closure. In harmonic-loss tables, "
            "bind every time-domain power/loss component to the common time axis knot by knot; "
            "matching only an integrated loss cannot prove sample identity. "
            "treat row zero as the aggregate summary: eddy aggregate reconstructs from "
            "positive-frequency bins, while hysteresis and combined iron loss can remain "
            "aggregate-only. Exclude R/L ratios when their denominator current is near zero."
        ),
    }
