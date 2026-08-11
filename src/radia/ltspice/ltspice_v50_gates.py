"""Stepped-analysis and noise-artifact identity gates."""

from __future__ import annotations

import itertools
import math
from collections.abc import Mapping

from .ltspice_v51_gates import validate_ltspice_v51_identity


STEP = "step_parameter_cartesian_nested_order_measure_row_owner_identity"
NOISE = "noise_input_output_source_contribution_bandwidth_integration_owner_identity"


def _digest(value: object) -> bool:
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _generation(contract: Mapping[str, object], *names: str) -> bool:
    value = str(contract.get("generation_id") or "")
    return bool(value) and all(contract.get(name) == value for name in names)


def _finite_list(
    values: object,
    *,
    length: int | None = None,
    positive: bool = False,
    nonnegative: bool = False,
) -> bool:
    if not isinstance(values, list) or not values or (length is not None and len(values) != length):
        return False
    for value in values:
        if not isinstance(value, (int, float)) or isinstance(value, bool) or not math.isfinite(float(value)):
            return False
        if positive and float(value) <= 0.0:
            return False
        if nonnegative and float(value) < 0.0:
            return False
    return True


def _row_key(row: list[object]) -> str:
    return ";".join(
        f"{row[index]}={format(float(row[index + 1]), '.15g')}"
        for index in range(0, len(row), 2)
    )


def _step_ok(contract: Mapping[str, object]) -> bool:
    names = contract.get("parameter_names")
    grid = contract.get("parameter_value_grid")
    nesting = contract.get("nesting_order")
    rows = contract.get("cartesian_step_rows")
    keys = contract.get("measure_row_keys")
    measures = contract.get("measure_values")
    if not (
        isinstance(names, list)
        and names
        and len(names) == len(set(names))
        and all(isinstance(name, str) and name for name in names)
        and isinstance(grid, Mapping)
        and set(grid) == set(names)
        and all(_finite_list(grid[name]) and len(grid[name]) == len(set(grid[name])) for name in names)
        and isinstance(nesting, list)
        and nesting == names
    ):
        return False
    expected_rows = [
        [item for pair in zip(nesting, values) for item in pair]
        for values in itertools.product(*(grid[name] for name in nesting))
    ]
    expected_keys = [_row_key(row) for row in expected_rows]
    return (
        _generation(
            contract,
            "parameter_generation_id",
            "cartesian_generation_id",
            "nesting_generation_id",
            "order_generation_id",
            "measure_generation_id",
            "owner_generation_id",
            "result_generation_id",
        )
        and rows == expected_rows
        and keys == expected_keys
        and _finite_list(measures, length=len(expected_rows))
        and contract.get("result_parameter_names") == names
        and contract.get("result_parameter_value_grid") == grid
        and contract.get("result_nesting_order") == nesting
        and contract.get("result_cartesian_step_rows") == rows
        and contract.get("result_measure_row_keys") == keys
        and contract.get("result_measure_values") == measures
        and str(contract.get("sweep_owner") or "").startswith("sweep:")
        and contract.get("result_sweep_owner") == contract.get("sweep_owner")
        and _digest(contract.get("result_sha256"))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def _noise_ok(contract: Mapping[str, object]) -> bool:
    input_reference = str(contract.get("input_reference") or "")
    output_reference = str(contract.get("output_reference") or "")
    frequencies = contract.get("frequency_hz")
    sources = contract.get("noise_source_order")
    contributions = contract.get("source_contribution_v2_per_hz")
    band = contract.get("integration_band_hz")
    integrated = contract.get("integrated_output_noise_v_rms")
    count = len(frequencies) if isinstance(frequencies, list) else 0
    traces_ok = (
        isinstance(sources, list)
        and bool(sources)
        and len(sources) == len(set(sources))
        and all(isinstance(source, str) and source for source in sources)
        and isinstance(contributions, Mapping)
        and set(contributions) == set(sources)
        and all(_finite_list(contributions[source], length=count, nonnegative=True) for source in sources)
    )
    integrated_ok = (
        isinstance(integrated, (int, float))
        and not isinstance(integrated, bool)
        and math.isfinite(float(integrated))
        and float(integrated) > 0.0
    )
    return (
        _generation(
            contract,
            "reference_generation_id",
            "source_generation_id",
            "frequency_generation_id",
            "bandwidth_generation_id",
            "integration_generation_id",
            "owner_generation_id",
            "result_generation_id",
        )
        and bool(input_reference)
        and bool(output_reference)
        and input_reference != output_reference
        and contract.get("result_input_reference") == input_reference
        and contract.get("result_output_reference") == output_reference
        and _finite_list(frequencies, positive=True)
        and count >= 2
        and all(float(left) < float(right) for left, right in zip(frequencies, frequencies[1:]))
        and contract.get("result_frequency_hz") == frequencies
        and traces_ok
        and contract.get("result_noise_source_order") == sources
        and contract.get("result_source_contribution_v2_per_hz") == contributions
        and _finite_list(band, length=2, positive=True)
        and float(band[0]) < float(band[1])
        and float(frequencies[0]) <= float(band[0])
        and float(band[1]) <= float(frequencies[-1])
        and contract.get("result_integration_band_hz") == band
        and integrated_ok
        and contract.get("result_integrated_output_noise_v_rms") == integrated
        and str(contract.get("trace_owner") or "").startswith("trace:")
        and contract.get("result_trace_owner") == contract.get("trace_owner")
        and _digest(contract.get("result_sha256"))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def validate_ltspice_v50_identity(positive: Mapping[str, object]) -> bool:
    if not isinstance(positive, Mapping):
        return False
    if not validate_ltspice_v51_identity(positive):
        return False
    step = positive.get(STEP)
    noise = positive.get(NOISE)
    if step is None and noise is None:
        return True
    return (
        isinstance(step, Mapping)
        and isinstance(noise, Mapping)
        and _step_ok(step)
        and _noise_ok(noise)
    )
