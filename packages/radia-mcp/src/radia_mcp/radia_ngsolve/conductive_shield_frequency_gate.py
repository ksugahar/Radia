"""Solver-neutral dual-regime magnetic and conductive shield sweep gate."""

from __future__ import annotations

import math
from typing import Any


def _complex_pair(value: object, name: str) -> complex:
    if not isinstance(value, list) or len(value) != 2:
        raise ValueError(f"{name} must contain [real, imag]")
    parsed = complex(float(value[0]), float(value[1]))
    if not math.isfinite(parsed.real) or not math.isfinite(parsed.imag):
        raise ValueError(f"{name} must be finite")
    return parsed


def _relative_error(left: complex | float, right: complex | float) -> float:
    return abs(left - right) / max(abs(left), abs(right), 1.0e-30)


def magnetic_conductive_shield_frequency_gate(
    summary: dict[str, Any],
    *,
    faraday_rtol: float = 1.0e-4,
    replay_rtol: float = 1.0e-12,
    primary_invariance_rtol: float = 5.0e-4,
    low_frequency_gain_minimum: float = 1.02,
    high_frequency_ratio_maximum: float = 0.80,
    response_flux_ratio_rtol: float = 5.0e-4,
) -> dict[str, Any]:
    """Gate low-frequency magnetic loading and high-frequency eddy shielding."""

    if not isinstance(summary, dict):
        raise ValueError("summary must be an object")
    models = summary.get("models")
    if not isinstance(models, list) or len(models) != 2:
        raise ValueError("models must contain baseline and shielded entries")

    parsed_models = []
    all_faraday_errors: list[float] = []
    replay_errors: list[float] = []
    for model in models:
        replays = model.get("replays") if isinstance(model, dict) else None
        if not isinstance(replays, list) or len(replays) != 2:
            raise ValueError("each model must contain exactly two replays")
        parsed_replays = []
        for replay in replays:
            rows = replay.get("rows") if isinstance(replay, dict) else None
            if not isinstance(rows, list) or len(rows) < 4:
                raise ValueError("each replay must contain at least four frequency rows")
            parsed_rows = []
            for row in rows:
                frequency = float(row.get("frequency_hz", math.nan))
                faraday_error = float(row.get("maximum_faraday_relative_error", math.nan))
                if not math.isfinite(frequency) or frequency <= 0.0:
                    raise ValueError("frequencies must be finite and positive")
                if not math.isfinite(faraday_error) or faraday_error < 0.0:
                    raise ValueError("Faraday errors must be finite and nonnegative")
                all_faraday_errors.append(faraday_error)
                parsed_rows.append(
                    {
                        "frequency_hz": frequency,
                        "primary_response": _complex_pair(row.get("primary_response"), "primary_response"),
                        "secondary_response": _complex_pair(row.get("secondary_response"), "secondary_response"),
                        "secondary_flux_linkage": _complex_pair(row.get("secondary_flux_linkage"), "secondary_flux_linkage"),
                        "maximum_faraday_relative_error": faraday_error,
                    }
                )
            frequencies = [row["frequency_hz"] for row in parsed_rows]
            if not all(right > left for left, right in zip(frequencies, frequencies[1:])):
                raise ValueError("frequency rows must be strictly increasing")
            parsed_replays.append(parsed_rows)
        for left, right in zip(parsed_replays[0], parsed_replays[1], strict=True):
            replay_errors.append(_relative_error(left["frequency_hz"], right["frequency_hz"]))
            for key in ("primary_response", "secondary_response", "secondary_flux_linkage"):
                replay_errors.append(_relative_error(left[key], right[key]))
        parsed_models.append(parsed_replays)

    baseline = parsed_models[0][0]
    shielded = parsed_models[1][0]
    baseline_frequencies = [row["frequency_hz"] for row in baseline]
    shielded_frequencies = [row["frequency_hz"] for row in shielded]
    grid_errors = [
        _relative_error(left, right)
        for left, right in zip(baseline_frequencies, shielded_frequencies, strict=True)
    ]
    primary_errors = []
    secondary_ratios = []
    flux_ratios = []
    ratio_gaps = []
    for open_row, shield_row in zip(baseline, shielded, strict=True):
        primary_errors.append(
            _relative_error(open_row["primary_response"], shield_row["primary_response"])
        )
        secondary_ratio = abs(shield_row["secondary_response"]) / max(
            abs(open_row["secondary_response"]), 1.0e-30
        )
        flux_ratio = abs(shield_row["secondary_flux_linkage"]) / max(
            abs(open_row["secondary_flux_linkage"]), 1.0e-30
        )
        secondary_ratios.append(secondary_ratio)
        flux_ratios.append(flux_ratio)
        ratio_gaps.append(_relative_error(secondary_ratio, flux_ratio))

    crossover_indices = [
        index
        for index, (left, right) in enumerate(zip(secondary_ratios, secondary_ratios[1:]))
        if (left - 1.0) * (right - 1.0) < 0.0
    ]
    crossover_bracket = None
    if len(crossover_indices) == 1:
        index = crossover_indices[0]
        crossover_bracket = [baseline_frequencies[index], baseline_frequencies[index + 1]]
    post_crossover_nonincreasing = bool(crossover_indices) and all(
        right <= left * (1.0 + 1.0e-9)
        for left, right in zip(
            secondary_ratios[crossover_indices[0] :],
            secondary_ratios[crossover_indices[0] + 1 :],
        )
    )

    timing = summary.get("timing_breakdown_s")
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
        "frequency_grids_match": len(baseline) == len(shielded)
        and max(grid_errors, default=math.inf) <= replay_rtol,
        "two_replays_are_deterministic_per_model": max(replay_errors, default=math.inf)
        <= replay_rtol,
        "faraday_identity_is_closed": max(all_faraday_errors, default=math.inf)
        <= faraday_rtol,
        "primary_response_is_nearly_invariant": max(primary_errors, default=math.inf)
        <= primary_invariance_rtol,
        "low_frequency_magnetic_loading_increases_coupling": secondary_ratios[0]
        >= low_frequency_gain_minimum,
        "high_frequency_conductive_shield_attenuates_coupling": secondary_ratios[-1]
        <= high_frequency_ratio_maximum,
        "single_gain_to_attenuation_crossover": len(crossover_indices) == 1,
        "post_crossover_attenuation_is_nonincreasing": post_crossover_nonincreasing,
        "response_and_flux_ratios_agree": max(ratio_gaps, default=math.inf)
        <= response_flux_ratio_rtol,
        "exactly_four_timing_stages": timing_ok,
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "magnetic_conductive_shield_frequency_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "metrics": {
            "frequency_count": len(baseline),
            "maximum_faraday_relative_error": max(all_faraday_errors),
            "maximum_replay_relative_error": max(replay_errors),
            "maximum_primary_response_relative_change": max(primary_errors),
            "low_frequency_secondary_coupling_ratio": secondary_ratios[0],
            "high_frequency_secondary_coupling_ratio": secondary_ratios[-1],
            "crossover_frequency_bracket_hz": crossover_bracket,
            "maximum_response_flux_ratio_relative_gap": max(ratio_gaps),
        },
        "lesson": (
            "A magnetically permeable conductive shield can increase coupling at low "
            "frequency and attenuate it at high frequency. Validate both regimes, the "
            "single crossover, Faraday closure, primary invariance, and independent replay."
        ),
    }
