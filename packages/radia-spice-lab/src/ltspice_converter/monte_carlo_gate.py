"""Solver-neutral tolerance Monte Carlo family gate."""

from __future__ import annotations

import math
from collections.abc import Mapping


def _number(value: object, name: str, *, positive: bool = False) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(parsed):
        raise ValueError(f"{name} must be finite")
    if positive and parsed <= 0.0:
        raise ValueError(f"{name} must be positive")
    return parsed


def _statistics(group: Mapping[str, object], name: str) -> dict[str, float]:
    if not isinstance(group, Mapping):
        raise ValueError(f"{name} must be an object")
    return {
        "nominal": _number(group.get("nominal"), f"{name}.nominal", positive=True),
        "mean": _number(group.get("mean"), f"{name}.mean", positive=True),
        "standard_deviation": _number(
            group.get("standard_deviation"),
            f"{name}.standard_deviation",
            positive=True,
        ),
    }


def monte_carlo_tolerance_family_gate(summary: Mapping[str, object]) -> dict:
    """Gate uniform component tolerance statistics and root-N averaging.

    Equal, independent component deviations uniformly distributed over
    ``[-tol, +tol]`` have relative standard deviation ``tol/sqrt(3)``.
    Series/parallel combinations of ``N`` equal parts reduce the first-order
    relative spread by ``sqrt(N)``.  A symmetric two-arm divider has relative
    output spread ``tol/sqrt(6*N)`` when each arm contains ``N`` equal parts.
    """

    if not isinstance(summary, Mapping):
        raise ValueError("summary must be an object")
    units = summary.get("units")
    resistor = summary.get("resistor_equivalent")
    divider = summary.get("divider")
    tolerances = summary.get("gate_tolerances")
    if not all(isinstance(item, Mapping) for item in (units, resistor, divider, tolerances)):
        raise ValueError("units, resistor_equivalent, divider, and gate_tolerances must be objects")

    sample_count = int(_number(summary.get("sample_count"), "sample_count", positive=True))
    part_count = int(
        _number(
            summary.get("independent_parts_per_equivalent"),
            "independent_parts_per_equivalent",
            positive=True,
        )
    )
    tolerance_fraction = _number(
        summary.get("tolerance_fraction"), "tolerance_fraction", positive=True
    )
    if part_count < 2:
        raise ValueError("independent_parts_per_equivalent must be at least two")
    if tolerance_fraction >= 1.0:
        raise ValueError("tolerance_fraction must be less than one")

    rows = {
        "single_resistor": _statistics(resistor.get("single"), "resistor_equivalent.single"),
        "series_resistor": _statistics(resistor.get("series"), "resistor_equivalent.series"),
        "parallel_resistor": _statistics(resistor.get("parallel"), "resistor_equivalent.parallel"),
        "single_divider": _statistics(divider.get("single_per_arm"), "divider.single_per_arm"),
        "multi_divider": _statistics(divider.get("multi_per_arm"), "divider.multi_per_arm"),
    }
    relative_sigma = {
        name: row["standard_deviation"] / row["mean"] for name, row in rows.items()
    }
    expected = {
        "single_resistor": tolerance_fraction / math.sqrt(3.0),
        "series_resistor": tolerance_fraction / math.sqrt(3.0 * part_count),
        "parallel_resistor": tolerance_fraction / math.sqrt(3.0 * part_count),
        "single_divider": tolerance_fraction / math.sqrt(6.0),
        "multi_divider": tolerance_fraction / math.sqrt(6.0 * part_count),
    }
    theory_error = {
        name: abs(relative_sigma[name] - expected[name]) / expected[name]
        for name in rows
    }
    mean_error = {
        name: abs(row["mean"] - row["nominal"]) / row["nominal"]
        for name, row in rows.items()
    }
    combined_resistor_sigma = 0.5 * (
        relative_sigma["series_resistor"] + relative_sigma["parallel_resistor"]
    )
    series_parallel_difference = abs(
        relative_sigma["series_resistor"] - relative_sigma["parallel_resistor"]
    ) / combined_resistor_sigma
    resistor_reduction = combined_resistor_sigma / relative_sigma["single_resistor"]
    divider_reduction = relative_sigma["multi_divider"] / relative_sigma["single_divider"]
    expected_reduction = 1.0 / math.sqrt(part_count)

    minimum_samples = int(
        _number(tolerances.get("minimum_samples"), "gate_tolerances.minimum_samples", positive=True)
    )
    max_mean_error = _number(
        tolerances.get("mean_relative_error"),
        "gate_tolerances.mean_relative_error",
        positive=True,
    )
    max_theory_error = _number(
        tolerances.get("relative_sigma_theory_error"),
        "gate_tolerances.relative_sigma_theory_error",
        positive=True,
    )
    max_series_parallel_difference = _number(
        tolerances.get("series_parallel_relative_sigma_difference"),
        "gate_tolerances.series_parallel_relative_sigma_difference",
        positive=True,
    )
    max_reduction_error = _number(
        tolerances.get("root_n_reduction_relative_error"),
        "gate_tolerances.root_n_reduction_relative_error",
        positive=True,
    )

    checks = {
        "uniform_independent_tolerance_declared": summary.get("distribution")
        == "independent_uniform_symmetric",
        "units_explicit": units.get("output") == "V" and units.get("tolerance") == "fraction",
        "sample_count_sufficient": sample_count >= minimum_samples,
        "means_remain_nominal": max(mean_error.values()) <= max_mean_error,
        "relative_sigmas_match_uniform_theory": max(theory_error.values()) <= max_theory_error,
        "series_parallel_relative_sigmas_agree": (
            series_parallel_difference <= max_series_parallel_difference
        ),
        "resistor_family_follows_root_n_reduction": (
            abs(resistor_reduction - expected_reduction) / expected_reduction
            <= max_reduction_error
        ),
        "divider_family_follows_root_n_reduction": (
            abs(divider_reduction - expected_reduction) / expected_reduction
            <= max_reduction_error
        ),
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "schema": "radia-spice-monte-carlo-tolerance-family/v1",
        "policy": "uniform_tolerance_root_n_statistics_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "metrics": {
            "relative_standard_deviation": relative_sigma,
            "expected_relative_standard_deviation": expected,
            "relative_sigma_theory_error": theory_error,
            "mean_relative_error": mean_error,
            "series_parallel_relative_sigma_difference": series_parallel_difference,
            "resistor_root_n_reduction": resistor_reduction,
            "divider_root_n_reduction": divider_reduction,
            "expected_root_n_reduction": expected_reduction,
        },
        "notes": [
            "A tolerance function defines a distribution; a nominal tolerance alone is not a standard deviation.",
            "Root-N reduction assumes independent equal deviations and is a first-order result for parallel and divider combinations.",
            "Preserve sample count, distribution, seed policy, and raw step axis with every statistical result artifact.",
        ],
    }
