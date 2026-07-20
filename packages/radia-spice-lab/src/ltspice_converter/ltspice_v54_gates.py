"""Monte Carlo and loop-gain replay identity gates."""

from __future__ import annotations

import math
from collections.abc import Mapping


MONTE = "montecarlo_seed_distribution_parameter_yield_sample_owner_identity"
LOOP = "loopgain_injection_breakpoint_sign_crossover_phasemargin_owner_identity"


def _digest(value: object) -> bool:
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _number(value: object, *, positive: bool = False) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value)) and (not positive or float(value) > 0.0)


def _generation(contract: Mapping[str, object], *fields: str) -> bool:
    generation = str(contract.get("generation_id") or "")
    return bool(generation) and all(contract.get(field) == generation for field in fields)


def _distribution_ok(value: object) -> bool:
    if not isinstance(value, Mapping):
        return False
    kind = value.get("distribution")
    if kind == "normal":
        return set(value) == {"distribution", "mean", "stddev"} and _number(value.get("mean")) and _number(value.get("stddev"), positive=True)
    if kind == "uniform":
        return set(value) == {"distribution", "lower", "upper"} and _number(value.get("lower")) and _number(value.get("upper")) and float(value["lower"]) < float(value["upper"])
    return False


def _monte_ok(contract: Mapping[str, object]) -> bool:
    seed = contract.get("rng_seed")
    distributions = contract.get("parameter_distributions")
    samples = contract.get("sampled_values")
    distributions_ok = (
        isinstance(distributions, Mapping) and bool(distributions)
        and all(isinstance(name, str) and name and _distribution_ok(value) for name, value in distributions.items())
    )
    samples_ok = isinstance(samples, list) and bool(samples)
    passed_count = 0
    if samples_ok:
        for expected_id, sample in enumerate(samples):
            if not isinstance(sample, Mapping) or set(sample) != {"sample_id", "parameters", "passed"}:
                samples_ok = False
                break
            parameters = sample["parameters"]
            if not (
                sample["sample_id"] == expected_id and isinstance(sample["passed"], bool)
                and isinstance(parameters, Mapping) and distributions_ok and set(parameters) == set(distributions)
                and all(_number(value) for value in parameters.values())
            ):
                samples_ok = False
                break
            for name, value in parameters.items():
                distribution = distributions[name]
                if distribution["distribution"] == "uniform" and not (float(distribution["lower"]) <= float(value) <= float(distribution["upper"])):
                    samples_ok = False
                    break
            if not samples_ok:
                break
            passed_count += int(sample["passed"])
    yield_count = contract.get("yield_count")
    return (
        _generation(contract, "seed_generation_id", "distribution_generation_id", "sample_generation_id", "yield_generation_id", "owner_generation_id", "result_generation_id")
        and isinstance(seed, int) and not isinstance(seed, bool) and seed >= 0 and contract.get("result_rng_seed") == seed
        and distributions_ok and contract.get("result_parameter_distributions") == distributions
        and samples_ok and contract.get("result_sampled_values") == samples
        and isinstance(yield_count, int) and not isinstance(yield_count, bool) and yield_count == passed_count and contract.get("result_yield_count") == yield_count
        and str(contract.get("run_owner") or "").startswith("run:") and contract.get("result_run_owner") == contract.get("run_owner")
        and _digest(contract.get("result_sha256")) and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def _loop_ok(contract: Mapping[str, object]) -> bool:
    rows = contract.get("loop_gain_rows")
    rows_ok = isinstance(rows, list) and len(rows) >= 3
    frequencies: list[float] = []
    crossover_rows: list[Mapping[str, object]] = []
    if rows_ok:
        for row in rows:
            if not isinstance(row, Mapping) or set(row) != {"frequency_hz", "real", "imag"} or not all(_number(row.get(field), positive=(field == "frequency_hz")) for field in ("frequency_hz", "real", "imag")):
                rows_ok = False
                break
            frequencies.append(float(row["frequency_hz"]))
            if math.isclose(math.hypot(float(row["real"]), float(row["imag"])), 1.0, rel_tol=1.0e-10, abs_tol=1.0e-12):
                crossover_rows.append(row)
        rows_ok = rows_ok and all(left < right for left, right in zip(frequencies, frequencies[1:])) and len(crossover_rows) == 1
    crossover = contract.get("crossover_frequency_hz")
    margin = contract.get("phase_margin_deg")
    if rows_ok:
        crossover_row = crossover_rows[0]
        expected_margin = 180.0 + math.degrees(math.atan2(float(crossover_row["imag"]), float(crossover_row["real"])))
    else:
        crossover_row = {}
        expected_margin = math.nan
    return (
        _generation(contract, "injection_generation_id", "sign_generation_id", "crossover_generation_id", "margin_generation_id", "owner_generation_id", "result_generation_id")
        and str(contract.get("injection_point") or "").startswith("node:") and contract.get("result_injection_point") == contract.get("injection_point")
        and contract.get("sign_convention") == "negative_feedback_return_ratio" and contract.get("result_sign_convention") == contract.get("sign_convention")
        and rows_ok and contract.get("result_loop_gain_rows") == rows
        and _number(crossover, positive=True) and math.isclose(float(crossover), float(crossover_row["frequency_hz"]), rel_tol=1.0e-12, abs_tol=1.0e-12) and contract.get("result_crossover_frequency_hz") == crossover
        and _number(margin) and math.isclose(float(margin), expected_margin, rel_tol=1.0e-10, abs_tol=1.0e-10) and 0.0 < float(margin) < 180.0 and contract.get("result_phase_margin_deg") == margin
        and str(contract.get("trace_owner") or "").startswith("trace:") and contract.get("result_trace_owner") == contract.get("trace_owner")
        and _digest(contract.get("result_sha256")) and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def validate_ltspice_v54_identity(positive: Mapping[str, object]) -> bool:
    """Validate optional v54 identities, requiring the pair when either is present."""
    if not isinstance(positive, Mapping):
        return False
    monte = positive.get(MONTE)
    loop = positive.get(LOOP)
    if monte is None and loop is None:
        return True
    return isinstance(monte, Mapping) and isinstance(loop, Mapping) and _monte_ok(monte) and _loop_ok(loop)
