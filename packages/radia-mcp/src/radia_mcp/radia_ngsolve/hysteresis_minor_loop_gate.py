"""Solver-neutral replay gate for stateful magnetic minor loops."""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Mapping, Sequence


def _mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be an object")
    return value


def _series(container: Mapping[str, object], key: str, *, minimum: int = 2) -> list[float]:
    value = container.get(key)
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError(f"{key} must be an array")
    rows = []
    for item in value:
        try:
            number = float(item)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{key} must contain finite numbers") from exc
        if not math.isfinite(number):
            raise ValueError(f"{key} must contain finite numbers")
        rows.append(number)
    if len(rows) < minimum:
        raise ValueError(f"{key} must contain at least {minimum} values")
    return rows


def _same_length(name: str, *values: Sequence[object]) -> None:
    if len({len(value) for value in values}) != 1:
        raise ValueError(f"{name} arrays must have equal lengths")


def _strictly_increasing(values: Sequence[float]) -> bool:
    return all(right > left for left, right in zip(values, values[1:]))


def _maximum_difference(left: Sequence[float], right: Sequence[float]) -> float:
    if len(left) != len(right):
        return math.inf
    return max((abs(a - b) for a, b in zip(left, right)), default=0.0)


def _reversal_count(values: Sequence[float], tolerance: float = 1.0e-12) -> int:
    signs = []
    for left, right in zip(values, values[1:]):
        difference = right - left
        if abs(difference) > tolerance:
            signs.append(1 if difference > 0.0 else -1)
    return sum(left != right for left, right in zip(signs, signs[1:]))


def _single_knot_alignment(
    saved_rows: Sequence[tuple[float, float, float]],
    fresh_rows: Sequence[tuple[float, float, float]],
) -> dict[str, object]:
    if len(saved_rows) == len(fresh_rows):
        difference = max(
            abs(a - b)
            for saved, fresh in zip(saved_rows, fresh_rows)
            for a, b in zip(saved, fresh)
        )
        return {"matched": True, "dropped_fresh_row": None, "maximum_absolute_error": difference}
    if len(fresh_rows) != len(saved_rows) + 1:
        return {"matched": False, "dropped_fresh_row": None, "maximum_absolute_error": math.inf}
    best_error = math.inf
    best_drop = None
    for drop in range(len(fresh_rows)):
        error = max(
            abs(a - b)
            for index, saved in enumerate(saved_rows)
            for a, b in zip(saved, fresh_rows[index + (index >= drop)])
        )
        if error < best_error:
            best_error = error
            best_drop = drop
    return {"matched": True, "dropped_fresh_row": best_drop, "maximum_absolute_error": best_error}


def _integral(time_s: Sequence[float], values: Sequence[float]) -> float:
    return sum(
        0.5 * (values[index] + values[index + 1]) * (time_s[index + 1] - time_s[index])
        for index in range(len(values) - 1)
    )


def hysteresis_minor_loop_replay_gate(
    summary: Mapping[str, object],
    *,
    minimum_reversal_count: int = 20,
    minimum_branch_spread: float = 1.0e-3,
    minimum_path_integral: float = 1.0e-3,
    maximum_saved_alignment_error: float = 1.0e-12,
    maximum_repeat_error: float = 1.0e-15,
    maximum_reference_relative_error: float = 2.0e-4,
    maximum_joule_power: float = 1.0e-15,
    maximum_loss_identity_error: float = 1.0e-15,
) -> dict[str, object]:
    """Validate a solved minor loop without assuming pointwise-positive loss power.

    A fresh result may expose one initial or repeated knot that a serialized
    result omits. The gate permits exactly one such row and discovers its
    position; broader row loss, reordering, or path drift is rejected.
    """

    if not isinstance(summary, Mapping):
        raise TypeError("summary must be an object")
    fresh = _mapping(summary.get("fresh"), "fresh")
    repeat = _mapping(summary.get("repeat"), "repeat")
    saved = _mapping(summary.get("saved_reference"), "saved_reference")
    reference = _mapping(summary.get("historical_reference"), "historical_reference")
    losses = _mapping(summary.get("losses"), "losses")

    time_s = _series(fresh, "time_s", minimum=9)
    drive = _series(fresh, "drive", minimum=9)
    response = _series(fresh, "response", minimum=9)
    _same_length("fresh", time_s, drive, response)
    repeat_time = _series(repeat, "time_s", minimum=9)
    repeat_drive = _series(repeat, "drive", minimum=9)
    repeat_response = _series(repeat, "response", minimum=9)
    _same_length("repeat", repeat_time, repeat_drive, repeat_response)
    saved_time = _series(saved, "time_s", minimum=8)
    saved_drive = _series(saved, "drive", minimum=8)
    saved_response = _series(saved, "response", minimum=8)
    _same_length("saved_reference", saved_time, saved_drive, saved_response)

    saved_alignment = _single_knot_alignment(
        list(zip(saved_time, saved_drive, saved_response)),
        list(zip(time_s, drive, response)),
    )
    repeat_error = max(
        _maximum_difference(time_s, repeat_time),
        _maximum_difference(drive, repeat_drive),
        _maximum_difference(response, repeat_response),
    )

    branch_values: dict[float, list[float]] = defaultdict(list)
    for input_value, output_value in zip(drive, response):
        branch_values[round(input_value, 9)].append(output_value)
    branch_spreads = [
        max(values) - min(values) for values in branch_values.values() if len(values) > 1
    ]
    maximum_branch_spread = max(branch_spreads, default=0.0)
    path_integral = sum(
        0.5 * (response[index] + response[index + 1]) * (drive[index + 1] - drive[index])
        for index in range(len(drive) - 1)
    )

    reference_time = _series(reference, "time_s", minimum=8)
    reference_magnitude = _series(reference, "response_magnitude", minimum=8)
    _same_length("historical_reference", reference_time, reference_magnitude)
    fresh_by_time = {round(value, 12): abs(response[index]) for index, value in enumerate(time_s)}
    matched = [
        (expected, fresh_by_time[round(time_value, 12)])
        for time_value, expected in zip(reference_time, reference_magnitude)
        if round(time_value, 12) in fresh_by_time
    ]
    reference_maximum_error = max(
        (abs(expected - observed) for expected, observed in matched), default=math.inf
    )
    reference_scale = max((abs(observed) for _, observed in matched), default=0.0)
    reference_relative_error = reference_maximum_error / max(reference_scale, 1.0e-300)

    loss_time = _series(losses, "time_s", minimum=2)
    joule_power = _series(losses, "joule_power_W", minimum=2)
    hysteresis_power = _series(losses, "hysteresis_power_W", minimum=2)
    iron_power = _series(losses, "iron_power_W", minimum=2)
    _same_length("losses", loss_time, joule_power, hysteresis_power, iron_power)
    maximum_joule = max(abs(value) for value in joule_power)
    maximum_loss_difference = _maximum_difference(hysteresis_power, iron_power)
    hysteresis_energy = _integral(loss_time, hysteresis_power)
    iron_energy = _integral(loss_time, iron_power)

    checks = {
        "fresh_time_is_strictly_increasing": _strictly_increasing(time_s),
        "signed_drive_spans_both_polarities": min(drive) < 0.0 < max(drive),
        "minor_loop_has_many_reversals": _reversal_count(drive) >= int(minimum_reversal_count),
        "response_is_multivalued_at_repeated_drive": maximum_branch_spread >= float(minimum_branch_spread),
        "path_integral_is_nonzero": abs(path_integral) >= float(minimum_path_integral),
        "saved_result_matches_after_one_knot_normalization": bool(saved_alignment["matched"])
        and float(saved_alignment["maximum_absolute_error"]) <= float(maximum_saved_alignment_error),
        "fresh_repeat_is_exact": repeat_error <= float(maximum_repeat_error),
        "historical_reference_time_grid_is_owned": len(matched) == len(reference_time),
        "historical_response_magnitude_agrees": reference_relative_error <= float(maximum_reference_relative_error),
        "eddy_free_joule_power_is_zero": maximum_joule <= float(maximum_joule_power),
        "iron_power_equals_hysteresis_power": maximum_loss_difference <= float(maximum_loss_identity_error),
        "signed_hysteresis_power_integrates_to_positive_loss": min(hysteresis_power) < 0.0 < max(hysteresis_power)
        and hysteresis_energy > 0.0
        and iron_energy > 0.0,
    }
    issues = [name for name, passed in checks.items() if not passed]
    return {
        "policy": "hysteresis_minor_loop_replay_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "metrics": {
            "sample_count": len(time_s),
            "drive_reversal_count": _reversal_count(drive),
            "repeated_drive_level_count": len(branch_spreads),
            "maximum_response_branch_spread": maximum_branch_spread,
            "path_integral_response_d_drive": path_integral,
            "saved_alignment": saved_alignment,
            "repeat_maximum_absolute_error": repeat_error,
            "historical_reference_maximum_absolute_error": reference_maximum_error,
            "historical_reference_relative_error": reference_relative_error,
            "maximum_joule_power_W": maximum_joule,
            "maximum_iron_hysteresis_power_error_W": maximum_loss_difference,
            "hysteresis_energy_J": hysteresis_energy,
            "iron_energy_J": iron_energy,
        },
        "notes": [
            "A minor loop must be multivalued at repeated drive levels; a single-valued B-H interpolation is not a hysteresis replay.",
            "Normalize at most one serialization-omitted initial or repeated knot, and reject broader row loss or reordering.",
            "Instantaneous hysteresis power may be signed. Gate positive integrated loss, not pointwise nonnegativity.",
        ],
    }
