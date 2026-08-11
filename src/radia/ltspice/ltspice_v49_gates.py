"""Monte-Carlo and switch-state artifact identity gates."""

from __future__ import annotations

import math
from collections.abc import Mapping


MONTE_CARLO = "monte_carlo_seed_distribution_sample_parameter_run_order_measure_owner_identity"
SWITCH = "switch_hysteresis_state_timestep_breakpoint_raw_trace_owner_identity"


def _digest(value: object) -> bool:
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _generation(contract: Mapping[str, object], *names: str) -> bool:
    value = str(contract.get("generation_id") or "")
    return bool(value) and all(contract.get(name) == value for name in names)


def _finite(values: object, *, length: int | None = None, minimum: int = 1) -> bool:
    return isinstance(values, list) and len(values) >= minimum and (length is None or len(values) == length) and all(isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value)) for value in values)


def _monte_carlo_ok(contract: Mapping[str, object]) -> bool:
    distribution = contract.get("distribution")
    samples = contract.get("sample_values")
    count = len(samples) if isinstance(samples, list) else 0
    parameters = contract.get("parameter_values")
    order = contract.get("run_order")
    measures = contract.get("measure_values")
    return (
        _generation(contract, "seed_generation_id", "distribution_generation_id", "sample_generation_id", "parameter_generation_id", "run_generation_id", "measure_generation_id", "result_generation_id")
        and isinstance(contract.get("random_seed"), int) and not isinstance(contract.get("random_seed"), bool) and contract["random_seed"] >= 0 and contract.get("result_random_seed") == contract.get("random_seed")
        and isinstance(distribution, Mapping) and distribution.get("kind") == "normal" and isinstance(distribution.get("mean"), (int, float)) and math.isfinite(float(distribution["mean"])) and isinstance(distribution.get("stddev"), (int, float)) and math.isfinite(float(distribution["stddev"])) and float(distribution["stddev"]) > 0.0 and contract.get("result_distribution") == distribution
        and _finite(samples, minimum=3) and contract.get("result_sample_values") == samples
        and _finite(parameters, length=count) and contract.get("result_parameter_values") == parameters
        and isinstance(order, list) and order == [f"run:{index}" for index in range(count)] and contract.get("result_run_order") == order
        and _finite(measures, length=count) and contract.get("result_measure_values") == measures
        and str(contract.get("measure_owner") or "").startswith("measure:") and contract.get("result_measure_owner") == contract.get("measure_owner")
        and _digest(contract.get("result_sha256")) and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def _switch_ok(contract: Mapping[str, object]) -> bool:
    states = contract.get("hysteresis_state")
    times = contract.get("accepted_timesteps_s")
    breakpoints = contract.get("breakpoints_s")
    trace = contract.get("raw_trace_rows")
    count = len(times) if isinstance(times, list) else 0
    transitions = [times[index] for index in range(1, count) if states[index] != states[index - 1]] if isinstance(states, list) and len(states) == count else []
    trace_ok = isinstance(trace, list) and len(trace) == count and all(isinstance(row, list) and len(row) == 2 and row[0] == time and all(isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value)) for value in row) for row, time in zip(trace, times or []))
    return (
        _generation(contract, "hysteresis_generation_id", "state_generation_id", "timestep_generation_id", "breakpoint_generation_id", "raw_generation_id", "trace_generation_id", "owner_generation_id", "result_generation_id")
        and isinstance(states, list) and len(states) >= 3 and set(states) == {"off", "on"} and contract.get("result_hysteresis_state") == states
        and _finite(times, minimum=3) and all(float(left) < float(right) for left, right in zip(times, times[1:])) and contract.get("result_accepted_timesteps_s") == times
        and _finite(breakpoints) and breakpoints == transitions and contract.get("result_breakpoints_s") == breakpoints
        and trace_ok and contract.get("result_raw_trace_rows") == trace
        and str(contract.get("waveform_owner") or "").startswith("waveform:") and contract.get("result_waveform_owner") == contract.get("waveform_owner")
        and _digest(contract.get("result_sha256")) and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def validate_ltspice_v49_identity(positive: Mapping[str, object]) -> bool:
    if not isinstance(positive, Mapping):
        return False
    monte_carlo = positive.get(MONTE_CARLO); switch = positive.get(SWITCH)
    if monte_carlo is None and switch is None:
        return True
    return isinstance(monte_carlo, Mapping) and isinstance(switch, Mapping) and _monte_carlo_ok(monte_carlo) and _switch_ok(switch)
