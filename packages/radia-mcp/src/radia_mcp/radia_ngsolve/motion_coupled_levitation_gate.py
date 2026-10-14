"""Solver-neutral validation for motion-coupled eddy-current levitation."""

from __future__ import annotations

import bisect
import math
from typing import Any


_UNITS = {
    "time": "s",
    "displacement": "m",
    "velocity": "m/s",
    "acceleration": "m/s^2",
    "force": "N",
    "mass": "kg",
}


def _series(row: dict[str, Any], key: str) -> list[float]:
    value = row.get(key)
    if not isinstance(value, list):
        raise ValueError(f"{key} must be a list")
    parsed = [float(item) for item in value]
    if not parsed or not all(math.isfinite(item) for item in parsed):
        raise ValueError(f"{key} must contain finite values")
    return parsed


def _strictly_increasing(values: list[float]) -> bool:
    return len(values) >= 2 and all(
        right > left for left, right in zip(values, values[1:])
    )


def _mean(values: list[float]) -> float:
    return sum(values) / len(values)


def _rms(values: list[float]) -> float:
    return math.sqrt(_mean([value * value for value in values]))


def _span(values: list[float]) -> float:
    return max(values) - min(values)


def _maximum_span_error(
    left: list[float], right: list[float], reference: list[float]
) -> float:
    scale = _span(reference)
    if len(left) != len(right) or scale <= 0.0:
        return math.inf
    return max(abs(a - b) for a, b in zip(left, right, strict=True)) / scale


def _relative_step_spread(values: list[float]) -> tuple[float, float]:
    if not _strictly_increasing(values):
        return math.inf, math.inf
    steps = [right - left for left, right in zip(values, values[1:])]
    mean_step = _mean(steps)
    spread = (max(steps) - min(steps)) / mean_step
    return mean_step, spread


def _interpolate(times: list[float], values: list[float], query: float) -> float:
    if len(times) != len(values) or not times or query < times[0] or query > times[-1]:
        raise ValueError("interpolation query is outside the sampled interval")
    index = bisect.bisect_left(times, query)
    if index == 0:
        return values[0]
    if index == len(times):
        return values[-1]
    if times[index] == query:
        return values[index]
    left = index - 1
    weight = (query - times[left]) / (times[index] - times[left])
    return values[left] + weight * (values[index] - values[left])


def _correlation(left: list[float], right: list[float]) -> float:
    if len(left) != len(right) or len(left) < 2:
        return -math.inf
    left_mean = _mean(left)
    right_mean = _mean(right)
    covariance = sum(
        (a - left_mean) * (b - right_mean)
        for a, b in zip(left, right, strict=True)
    )
    left_norm = math.sqrt(sum((value - left_mean) ** 2 for value in left))
    right_norm = math.sqrt(sum((value - right_mean) ** 2 for value in right))
    if left_norm <= 0.0 or right_norm <= 0.0:
        return -math.inf
    return covariance / (left_norm * right_norm)


def motion_coupled_eddy_levitation_transient_gate(
    summary: dict[str, Any],
    *,
    max_equation_residual_over_g: float = 0.01,
    max_displacement_replay_error_over_span: float = 0.005,
    max_force_replay_error_over_span: float = 0.02,
    max_experiment_rmse_over_span: float = 0.10,
    min_experiment_correlation: float = 0.80,
    min_probe_samples_per_force_period: float = 10.0,
    max_aliased_output_to_probe_force_span_ratio: float = 0.05,
) -> dict[str, Any]:
    """Gate a freely moving conducting body without hiding force aliasing.

    A sinusoidal magnetic excitation commonly gives a force component at twice
    the drive frequency. Requested output times can therefore land at one force
    phase and make a transient force look nearly constant. When that alias risk
    is present, the gate requires an adaptive-step force history interpolated to
    a common fixed grid. Adaptive row counts themselves are deliberately not a
    replay identity.
    """
    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")
    contract = summary.get("contract")
    units = summary.get("units")
    replays = summary.get("replays")
    experiment = summary.get("experiment")
    timing = summary.get("timing_breakdown_s")
    if not isinstance(contract, dict):
        raise ValueError("contract must be a mapping")
    if not isinstance(units, dict):
        raise ValueError("units must be a mapping")
    if not isinstance(replays, list) or len(replays) < 2:
        raise ValueError("at least two replays are required")
    if not isinstance(experiment, dict):
        raise ValueError("experiment must be a mapping")

    mass_kg = float(summary.get("mass_kg", math.nan))
    gravity_m_s2 = float(summary.get("gravity_m_s2", math.nan))
    drive_frequency_hz = float(contract.get("drive_frequency_hz", math.nan))
    force_frequency_hz = float(
        contract.get("expected_force_harmonic_hz", math.nan)
    )
    if not all(
        math.isfinite(value) and value > 0.0
        for value in (mass_kg, gravity_m_s2, drive_frequency_hz, force_frequency_hz)
    ):
        raise ValueError("mass, gravity, drive frequency, and force frequency must be positive")

    parsed: list[dict[str, Any]] = []
    equation_errors: list[float] = []
    output_steps: list[float] = []
    output_step_spreads: list[float] = []
    probe_steps: list[float] = []
    probe_step_spreads: list[float] = []
    adaptive_rows: list[int] = []
    adaptive_samples_per_period: list[float] = []
    force_coverage_ratios: list[float] = []
    all_cardinalities_match = True
    all_times_increase = True
    for replay in replays:
        if not isinstance(replay, dict):
            raise ValueError("replay entries must be mappings")
        output = {
            key: _series(replay, key)
            for key in (
                "output_time_s",
                "output_displacement_m",
                "output_velocity_m_s",
                "output_acceleration_m_s2",
                "output_lift_force_n",
                "output_gravity_force_n",
            )
        }
        probe = {
            key: _series(replay, key)
            for key in (
                "probe_time_s",
                "probe_displacement_m",
                "probe_lift_force_n",
            )
        }
        all_cardinalities_match &= len({len(value) for value in output.values()}) == 1
        all_cardinalities_match &= len({len(value) for value in probe.values()}) == 1
        all_times_increase &= _strictly_increasing(output["output_time_s"])
        all_times_increase &= _strictly_increasing(probe["probe_time_s"])
        output_step, output_spread = _relative_step_spread(output["output_time_s"])
        probe_step, probe_spread = _relative_step_spread(probe["probe_time_s"])
        output_steps.append(output_step)
        output_step_spreads.append(output_spread)
        probe_steps.append(probe_step)
        probe_step_spreads.append(probe_spread)
        adaptive_rows.append(int(replay.get("adaptive_probe_row_count", 0)))
        adaptive_samples_per_period.append(
            float(replay.get("adaptive_probe_median_samples_per_force_period", 0.0))
        )
        acceleration = output["output_acceleration_m_s2"]
        lift = output["output_lift_force_n"]
        gravity_force = output["output_gravity_force_n"]
        equation_errors.extend(
            abs(a - (force - gravity_force_n) / mass_kg) / gravity_m_s2
            for a, force, gravity_force_n in zip(
                acceleration, lift, gravity_force, strict=True
            )
        )
        probe_force_span = _span(probe["probe_lift_force_n"])
        force_coverage_ratios.append(
            _span(lift) / probe_force_span if probe_force_span > 0.0 else math.inf
        )
        parsed.append({"output": output, "probe": probe})

    reference = parsed[0]
    displacement_replay_errors: list[float] = []
    force_replay_errors: list[float] = []
    output_time_replay_errors: list[float] = []
    probe_time_replay_errors: list[float] = []
    for replay in parsed[1:]:
        displacement_replay_errors.append(
            _maximum_span_error(
                reference["probe"]["probe_displacement_m"],
                replay["probe"]["probe_displacement_m"],
                reference["probe"]["probe_displacement_m"],
            )
        )
        force_replay_errors.append(
            _maximum_span_error(
                reference["probe"]["probe_lift_force_n"],
                replay["probe"]["probe_lift_force_n"],
                reference["probe"]["probe_lift_force_n"],
            )
        )
        output_time_replay_errors.append(
            max(
                (
                    abs(left - right)
                    for left, right in zip(
                        reference["output"]["output_time_s"],
                        replay["output"]["output_time_s"],
                        strict=True,
                    )
                ),
                default=math.inf,
            )
            if len(reference["output"]["output_time_s"])
            == len(replay["output"]["output_time_s"])
            else math.inf
        )
        probe_time_replay_errors.append(
            max(
                (
                    abs(left - right)
                    for left, right in zip(
                        reference["probe"]["probe_time_s"],
                        replay["probe"]["probe_time_s"],
                        strict=True,
                    )
                ),
                default=math.inf,
            )
            if len(reference["probe"]["probe_time_s"])
            == len(replay["probe"]["probe_time_s"])
            else math.inf
        )

    experiment_time = _series(experiment, "time_s")
    experiment_displacement = _series(experiment, "displacement_m")
    if len(experiment_time) != len(experiment_displacement):
        raise ValueError("experiment time and displacement lengths must match")
    probe_time = reference["probe"]["probe_time_s"]
    probe_displacement = reference["probe"]["probe_displacement_m"]
    overlap = [
        (time_s, displacement_m)
        for time_s, displacement_m in zip(
            experiment_time, experiment_displacement, strict=True
        )
        if probe_time[0] <= time_s <= probe_time[-1]
    ]
    simulated = [
        _interpolate(probe_time, probe_displacement, time_s)
        for time_s, _ in overlap
    ]
    measured = [value for _, value in overlap]
    experiment_span = _span(measured) if measured else 0.0
    experiment_errors = [
        actual - expected for actual, expected in zip(simulated, measured, strict=True)
    ]
    experiment_rmse_over_span = (
        _rms(experiment_errors) / experiment_span
        if experiment_errors and experiment_span > 0.0
        else math.inf
    )
    experiment_correlation = _correlation(simulated, measured)

    output_step = max(output_steps, default=math.inf)
    force_phase_advance = output_step * force_frequency_hz
    output_alias_risk = (
        math.isfinite(force_phase_advance)
        and abs(force_phase_advance - round(force_phase_advance)) <= 1.0e-9
    )
    fine_probe_required = 1.0 / (
        force_frequency_hz * float(min_probe_samples_per_force_period)
    )
    force_sampling_strategy_ok = not output_alias_risk or (
        contract.get("force_observation")
        == "adaptive_internal_steps_interpolated_to_fixed_grid"
        and max(probe_steps, default=math.inf) <= fine_probe_required * (1.0 + 1.0e-9)
        and min(adaptive_rows, default=0)
        > max(len(item["output"]["output_time_s"]) for item in parsed)
        and min(adaptive_samples_per_period, default=0.0)
        >= float(min_probe_samples_per_force_period)
    )
    aliased_output_not_substituted = not output_alias_risk or max(
        force_coverage_ratios, default=math.inf
    ) <= float(max_aliased_output_to_probe_force_span_ratio)

    timing_ok = False
    if isinstance(timing, dict) and len(timing) == 4:
        try:
            timing_ok = all(
                math.isfinite(float(value)) and float(value) >= 0.0
                for value in timing.values()
            )
        except (TypeError, ValueError):
            timing_ok = False

    checks = {
        "si_units_explicit": units == _UNITS,
        "force_harmonic_is_twice_drive": abs(
            force_frequency_hz - 2.0 * drive_frequency_hz
        )
        <= 1.0e-12 * force_frequency_hz,
        "motion_equation_contract_recorded": contract.get("motion_equation")
        == "mass_times_acceleration_equals_lift_minus_gravity",
        "experimental_displacement_contract_recorded": contract.get(
            "experimental_comparison"
        )
        == "displacement_time_history",
        "series_cardinalities_match": all_cardinalities_match,
        "time_axes_strictly_increase": all_times_increase,
        "requested_output_grid_is_uniform": max(
            output_step_spreads, default=math.inf
        )
        <= 1.0e-10,
        "fixed_probe_grid_is_uniform": max(probe_step_spreads, default=math.inf)
        <= 1.0e-10,
        "requested_output_alias_risk_characterized": output_alias_risk,
        "force_sampling_strategy_resolves_alias": force_sampling_strategy_ok,
        "aliased_output_force_not_substituted_for_probe": aliased_output_not_substituted,
        "motion_equation_closes": max(equation_errors, default=math.inf)
        <= float(max_equation_residual_over_g),
        "fixed_grid_time_axes_replay": max(
            output_time_replay_errors + probe_time_replay_errors, default=math.inf
        )
        <= 1.0e-12,
        "displacement_history_replays": max(
            displacement_replay_errors, default=math.inf
        )
        <= float(max_displacement_replay_error_over_span),
        "force_history_replays": max(force_replay_errors, default=math.inf)
        <= float(max_force_replay_error_over_span),
        "experimental_overlap_is_sufficient": len(overlap) >= 100,
        "experimental_rmse_is_bounded": experiment_rmse_over_span
        <= float(max_experiment_rmse_over_span),
        "experimental_shape_is_correlated": experiment_correlation
        >= float(min_experiment_correlation),
        "exactly_four_timing_stages": timing_ok,
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "motion_coupled_eddy_levitation_transient_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "metrics": {
            "replay_count": len(parsed),
            "requested_output_step_s": output_step,
            "force_phase_advance_per_requested_output_cycles": force_phase_advance,
            "requested_output_alias_risk": output_alias_risk,
            "fixed_probe_step_s": max(probe_steps, default=math.inf),
            "adaptive_probe_row_counts": adaptive_rows,
            "maximum_motion_equation_residual_over_g": max(
                equation_errors, default=math.inf
            ),
            "maximum_displacement_replay_error_over_span": max(
                displacement_replay_errors, default=math.inf
            ),
            "maximum_force_replay_error_over_span": max(
                force_replay_errors, default=math.inf
            ),
            "maximum_requested_output_to_probe_force_span_ratio": max(
                force_coverage_ratios, default=math.inf
            ),
            "experimental_overlap_count": len(overlap),
            "experimental_rmse_over_span": experiment_rmse_over_span,
            "experimental_correlation": experiment_correlation,
        },
        "tolerances": {
            "max_equation_residual_over_g": float(max_equation_residual_over_g),
            "max_displacement_replay_error_over_span": float(
                max_displacement_replay_error_over_span
            ),
            "max_force_replay_error_over_span": float(
                max_force_replay_error_over_span
            ),
            "max_experiment_rmse_over_span": float(max_experiment_rmse_over_span),
            "min_experiment_correlation": float(min_experiment_correlation),
            "min_probe_samples_per_force_period": float(
                min_probe_samples_per_force_period
            ),
            "max_aliased_output_to_probe_force_span_ratio": float(
                max_aliased_output_to_probe_force_span_ratio
            ),
        },
        "lesson": (
            "A requested output grid can sample the same phase of a twice-drive "
            "force harmonic and hide almost the entire force waveform. Preserve "
            "adaptive internal-step probes, interpolate observables to a common "
            "fixed grid, and compare histories rather than adaptive row counts."
        ),
    }
