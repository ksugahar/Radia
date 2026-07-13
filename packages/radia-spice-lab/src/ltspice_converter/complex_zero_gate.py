"""Solver-neutral gate for a second-order complex-zero transfer function."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any


def _finite(value: object, name: str, *, positive: bool = False) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(parsed) or (positive and parsed <= 0.0):
        requirement = "positive and finite" if positive else "finite"
        raise ValueError(f"{name} must be {requirement}")
    return parsed


def _complex_pair(value: object, name: str) -> complex:
    if (
        not isinstance(value, Sequence)
        or isinstance(value, (str, bytes))
        or len(value) != 2
    ):
        raise ValueError(f"{name} must contain [real, imag]")
    return complex(_finite(value[0], name), _finite(value[1], name))


def _recover_pair(values: list[complex]) -> tuple[float, float, float]:
    omega = math.sqrt(abs(values[0] * values[1]))
    damping = -(values[0] + values[1]).real
    quality = omega / max(damping, 1.0e-300)
    conjugacy = abs(values[0] - values[1].conjugate()) / max(omega, 1.0e-300)
    return omega / (2.0 * math.pi), quality, conjugacy


def second_order_complex_zero_gate(summary: Mapping[str, object]) -> dict[str, Any]:
    """Gate two-pole/two-complex-zero response, roots, asymptotes, and replay."""
    if not isinstance(summary, Mapping):
        raise ValueError("summary must be an object")
    model = summary.get("model_contract")
    runs = summary.get("runs")
    pole_zero = summary.get("pole_zero")
    metrics = summary.get("metrics")
    timing = summary.get("timing_breakdown_s")
    if not all(isinstance(value, Mapping) for value in (model, pole_zero, metrics)):
        raise ValueError("model_contract, pole_zero, and metrics must be objects")
    if (
        not isinstance(runs, Sequence)
        or isinstance(runs, (str, bytes))
        or len(runs) != 2
        or any(not isinstance(run, Mapping) for run in runs)
    ):
        raise ValueError("runs must contain exactly two objects")

    pole_frequency = _finite(
        model.get("pole_natural_frequency_hz"),
        "pole_natural_frequency_hz",
        positive=True,
    )
    pole_quality = _finite(
        model.get("pole_quality_factor"), "pole_quality_factor", positive=True
    )
    zero_frequency = _finite(
        model.get("zero_natural_frequency_hz"),
        "zero_natural_frequency_hz",
        positive=True,
    )
    zero_quality = _finite(
        model.get("zero_quality_factor"), "zero_quality_factor", positive=True
    )
    dc_gain = _finite(model.get("dc_gain"), "dc_gain", positive=True)
    high_gain = _finite(
        model.get("high_frequency_gain"), "high_frequency_gain", positive=True
    )
    expected_high_gain = dc_gain * (pole_frequency / zero_frequency) ** 2

    normalized_runs = []
    for index, run in enumerate(runs):
        normalized_runs.append(
            {
                "point_count": int(
                    _finite(
                        run.get("point_count"),
                        f"runs[{index}].point_count",
                        positive=True,
                    )
                ),
                "frequency_min_hz": _finite(
                    run.get("frequency_min_hz"),
                    f"runs[{index}].frequency_min_hz",
                    positive=True,
                ),
                "frequency_max_hz": _finite(
                    run.get("frequency_max_hz"),
                    f"runs[{index}].frequency_max_hz",
                    positive=True,
                ),
                "minimum_magnitude": _finite(
                    run.get("minimum_magnitude"),
                    f"runs[{index}].minimum_magnitude",
                    positive=True,
                ),
                "minimum_magnitude_frequency_hz": _finite(
                    run.get("minimum_magnitude_frequency_hz"),
                    f"runs[{index}].minimum_magnitude_frequency_hz",
                    positive=True,
                ),
                "analytic_minimum_magnitude": _finite(
                    run.get("analytic_minimum_magnitude"),
                    f"runs[{index}].analytic_minimum_magnitude",
                    positive=True,
                ),
                "analytic_minimum_frequency_hz": _finite(
                    run.get("analytic_minimum_frequency_hz"),
                    f"runs[{index}].analytic_minimum_frequency_hz",
                    positive=True,
                ),
            }
        )

    poles_raw = pole_zero.get("poles")
    zeros_raw = pole_zero.get("zeros")
    if (
        not isinstance(poles_raw, Sequence)
        or not isinstance(zeros_raw, Sequence)
        or len(poles_raw) != 2
        or len(zeros_raw) != 2
    ):
        raise ValueError("pole_zero must contain two poles and two zeros")
    poles = [_complex_pair(value, "pole") for value in poles_raw]
    zeros = [_complex_pair(value, "zero") for value in zeros_raw]
    recovered_pf, recovered_pq, pole_conjugacy = _recover_pair(poles)
    recovered_zf, recovered_zq, zero_conjugacy = _recover_pair(zeros)

    analytic_l2 = _finite(
        metrics.get("maximum_analytic_complex_relative_l2"),
        "maximum_analytic_complex_relative_l2",
    )
    analytic_point = _finite(
        metrics.get("maximum_analytic_complex_point_relative_error"),
        "maximum_analytic_complex_point_relative_error",
    )
    replay_error = _finite(
        metrics.get("maximum_complex_replay_relative_error"),
        "maximum_complex_replay_relative_error",
    )

    timing_ok = False
    if isinstance(timing, Mapping) and len(timing) == 4:
        try:
            timing_ok = all(_finite(value, "timing") >= 0.0 for value in timing.values())
        except ValueError:
            timing_ok = False

    checks = {
        "second_order_complex_zero_model_contract": model.get("topology")
        == "second_order_complex_zero"
        and zero_frequency > pole_frequency,
        "dc_and_high_frequency_gain_follow_frequency_ratio": abs(
            high_gain - expected_high_gain
        )
        / expected_high_gain
        <= 1.0e-12,
        "two_dense_ac_replays_bracket_poles_and_zeros": all(
            run["point_count"] >= 200
            and run["frequency_min_hz"] <= pole_frequency / 10.0
            and run["frequency_max_hz"] >= zero_frequency * 5.0
            for run in normalized_runs
        ),
        "complex_transfer_matches_two_pole_two_zero_identity": analytic_l2
        <= 5.0e-6
        and analytic_point <= 5.0e-6,
        "stable_conjugate_poles_recover_f0_and_q": all(
            value.real < 0.0 for value in poles
        )
        and pole_conjugacy <= 1.0e-12
        and abs(recovered_pf - pole_frequency) / pole_frequency <= 1.0e-12
        and abs(recovered_pq - pole_quality) / pole_quality <= 1.0e-12,
        "minimum_phase_conjugate_zeros_recover_fn_and_qn": all(
            value.real < 0.0 for value in zeros
        )
        and zero_conjugacy <= 1.0e-12
        and abs(recovered_zf - zero_frequency) / zero_frequency <= 1.0e-12
        and abs(recovered_zq - zero_quality) / zero_quality <= 1.0e-12,
        "complex_zero_has_finite_analytic_real_axis_dip": all(
            run["minimum_magnitude"] > 0.0
            and abs(
                run["minimum_magnitude"] - run["analytic_minimum_magnitude"]
            )
            / run["analytic_minimum_magnitude"]
            <= 5.0e-6
            and abs(
                run["minimum_magnitude_frequency_hz"]
                - run["analytic_minimum_frequency_hz"]
            )
            / run["analytic_minimum_frequency_hz"]
            <= 1.0e-12
            for run in normalized_runs
        ),
        "complex_observable_replay_is_deterministic": replay_error <= 1.0e-12,
        "exactly_four_timing_stages": timing_ok,
    }
    return {
        "policy": "second_order_complex_zero_transfer_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "expected_high_frequency_gain": expected_high_gain,
            "recovered_pole_natural_frequency_hz": recovered_pf,
            "recovered_pole_quality_factor": recovered_pq,
            "recovered_zero_natural_frequency_hz": recovered_zf,
            "recovered_zero_quality_factor": recovered_zq,
            "maximum_analytic_complex_relative_l2": analytic_l2,
            "maximum_analytic_complex_point_relative_error": analytic_point,
            "maximum_complex_replay_relative_error": replay_error,
        },
        "lesson": (
            "A complex-zero section is a full pole-zero transfer problem, not "
            "a real-axis notch. Reconstruct pole and zero frequencies and Q "
            "from the roots, enforce the DC/high-frequency gain ratio, compare "
            "the full complex response, retain the finite dip, and replay it."
        ),
    }
