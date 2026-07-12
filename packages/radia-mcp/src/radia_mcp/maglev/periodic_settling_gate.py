"""Periodic settling gate for rotating-conductor eddy-current responses."""

from __future__ import annotations

import math
from collections.abc import Iterable


def _finite(values: Iterable[float]) -> list[float]:
    result = [float(value) for value in values]
    if not result or not all(math.isfinite(value) for value in result):
        raise ValueError("response values must be finite and nonempty")
    return result


def _relative_l2(left: list[float], right: list[float]) -> float:
    residual = math.sqrt(sum((a - b) ** 2 for a, b in zip(left, right, strict=True)))
    scale = max(
        math.sqrt(sum(value * value for value in left)),
        math.sqrt(sum(value * value for value in right)),
        1.0e-30,
    )
    return residual / scale


def rotating_conductor_periodic_settling_gate(
    response: Iterable[float],
    *,
    steps_per_period: int,
    angle_step_deg: float,
    reference_response: Iterable[float] | None = None,
    maximum_final_period_relative_l2: float = 2.0e-3,
    maximum_contraction_factor: float = 0.35,
    reference_relative_l2_tolerance: float = 1.0e-8,
) -> dict[str, object]:
    """Validate periodic settling without mistaking the first turn for steady state."""

    values = _finite(response)
    steps = int(steps_per_period)
    if steps < 4 or len(values) < 3 * steps or len(values) % steps:
        raise ValueError("response must contain at least three complete equal periods")
    angle_step = float(angle_step_deg)
    if not math.isfinite(angle_step) or angle_step <= 0.0:
        raise ValueError("angle_step_deg must be finite and positive")
    periods = [values[index : index + steps] for index in range(0, len(values), steps)]
    transitions = [_relative_l2(left, right) for left, right in zip(periods, periods[1:])]
    contractions = [
        right / left if left > 0.0 else math.inf
        for left, right in zip(transitions, transitions[1:])
    ]
    reference_error = None
    if reference_response is not None:
        reference = _finite(reference_response)
        reference_error = _relative_l2(values, reference) if len(reference) == len(values) else math.inf

    checks = {
        "one_period_covers_full_rotation": abs(steps * angle_step - 360.0) <= 1.0e-10,
        "response_is_nonnegative": all(value >= 0.0 for value in values),
        "successive_period_error_strictly_decreases": all(
            right < left for left, right in zip(transitions, transitions[1:])
        ),
        "period_error_contracts_fast_enough": bool(contractions)
        and all(value <= float(maximum_contraction_factor) for value in contractions),
        "final_period_is_settled": transitions[-1] <= float(maximum_final_period_relative_l2),
        "reference_replay_matches": reference_error is None
        or reference_error <= float(reference_relative_l2_tolerance),
    }
    return {
        "policy": "rotating_conductor_periodic_settling_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "period_count": len(periods),
            "steps_per_period": steps,
            "successive_period_relative_l2": transitions,
            "contraction_factors": contractions,
            "final_period_relative_l2": transitions[-1],
            "reference_response_relative_l2": reference_error,
        },
        "lesson": (
            "A rotating-conductor transient is periodic-steady only after consecutive full-turn waveforms agree. "
            "Track contraction across several turns; do not accept the first revolution as the final response."
        ),
    }
