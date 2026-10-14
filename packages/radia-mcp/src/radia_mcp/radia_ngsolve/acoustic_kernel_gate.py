"""Solver-neutral gates for readable low-frequency acoustic kernels."""

from __future__ import annotations

import math
from typing import Any


def helmholtz_double_layer_low_frequency_gate(
    summary: dict[str, Any],
    *,
    min_positive_rows: int = 5,
    max_kr_abs: float = 1.0e-3,
    max_split_direct_relative_error: float = 1.0e-12,
    max_series_relative_error: float = 1.0e-6,
    max_quadratic_scaling_relative_error: float = 1.0e-6,
    min_smallest_cancellation_ratio: float = 1.0e8,
) -> dict[str, Any]:
    """Gate the low-frequency split of a source-normal Helmholtz double layer.

    For ``z = i*k*r``, the regular correction factor is
    ``exp(z)*(1-z)-1 = -z**2/2-z**3/3+...``.  The correction therefore starts
    quadratically in ``k*r``; validating only the total kernel can hide its
    cancellation error.
    """

    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")
    if int(min_positive_rows) < 3:
        raise ValueError("min_positive_rows must be at least 3")
    tolerances = {
        "max_kr_abs": float(max_kr_abs),
        "max_split_direct_relative_error": float(max_split_direct_relative_error),
        "max_series_relative_error": float(max_series_relative_error),
        "max_quadratic_scaling_relative_error": float(max_quadratic_scaling_relative_error),
        "min_smallest_cancellation_ratio": float(min_smallest_cancellation_ratio),
    }
    if any(not math.isfinite(value) or value < 0.0 for value in tolerances.values()):
        raise ValueError("tolerances must be finite and nonnegative")

    distance = float(summary.get("distance_m", math.nan))
    normal_dot = float(summary.get("normal_dot_m", math.nan))
    rows = summary.get("rows") or []
    if not isinstance(rows, list) or not all(isinstance(row, dict) for row in rows):
        raise ValueError("rows must be a list of mappings")

    parsed: list[dict[str, float]] = []
    required = (
        "wavenumber_per_m",
        "kr_abs",
        "laplace_real",
        "laplace_imag",
        "correction_real",
        "correction_imag",
        "split_direct_relative_error",
        "correction_series_relative_error",
    )
    for row in rows:
        values = {name: float(row.get(name, math.nan)) for name in required}
        if not all(math.isfinite(value) for value in values.values()):
            raise ValueError("required row values must be finite")
        cancellation = row.get("cancellation_ratio")
        values["cancellation_ratio"] = (
            math.inf if cancellation is None else float(cancellation)
        )
        parsed.append(values)

    positive = [row for row in parsed if row["wavenumber_per_m"] > 0.0]
    expected_laplace = (
        normal_dot / (4.0 * math.pi * distance**3)
        if math.isfinite(distance) and distance > 0.0 and math.isfinite(normal_dot)
        else math.nan
    )
    quadratic_reference = (
        expected_laplace * distance**2 / 2.0
        if math.isfinite(expected_laplace)
        else math.nan
    )
    kr_errors = [
        abs(row["kr_abs"] - row["wavenumber_per_m"] * distance)
        for row in parsed
    ] if math.isfinite(distance) else [math.inf]
    laplace_errors = [
        abs(complex(row["laplace_real"], row["laplace_imag"]) - expected_laplace)
        for row in parsed
    ] if math.isfinite(expected_laplace) else [math.inf]
    split_errors = [row["split_direct_relative_error"] for row in parsed]
    series_errors = [row["correction_series_relative_error"] for row in positive]
    quadratic_errors = [
        abs(row["correction_real"] / row["wavenumber_per_m"]**2 - quadratic_reference)
        / max(abs(quadratic_reference), 1.0e-300)
        for row in positive
    ] if math.isfinite(quadratic_reference) else [math.inf]
    positive_cancellation = [row["cancellation_ratio"] for row in positive]
    wavenumbers = [row["wavenumber_per_m"] for row in parsed]
    correction_magnitudes = [
        abs(complex(row["correction_real"], row["correction_imag"]))
        for row in positive
    ]
    scale = max(abs(expected_laplace), 1.0)
    kr_tolerance = max(1.0e-15, tolerances["max_kr_abs"] * 1.0e-12)

    checks = {
        "kernel_family_recorded": summary.get("kernel_family") == "helmholtz_source_normal_double_layer",
        "time_convention_recorded": summary.get("time_convention") == "exp(+i*k*r)",
        "distance_positive": math.isfinite(distance) and distance > 0.0,
        "normal_projection_physical": math.isfinite(normal_dot) and abs(normal_dot) <= distance,
        "row_count_sufficient": len(positive) >= int(min_positive_rows),
        "zero_wavenumber_row_present": bool(parsed) and wavenumbers[0] == 0.0,
        "wavenumbers_strictly_increasing": all(
            wavenumbers[index + 1] > wavenumbers[index]
            for index in range(max(0, len(wavenumbers) - 1))
        ),
        "kr_within_low_frequency_limit": bool(positive) and max(row["kr_abs"] for row in positive) <= tolerances["max_kr_abs"],
        "kr_matches_wavenumber_distance": max(kr_errors, default=math.inf) <= kr_tolerance,
        "laplace_part_matches_geometry": max(laplace_errors, default=math.inf) <= 1.0e-12 * scale,
        "split_matches_direct_kernel": max(split_errors, default=math.inf) <= tolerances["max_split_direct_relative_error"],
        "correction_matches_series": max(series_errors, default=math.inf) <= tolerances["max_series_relative_error"],
        "correction_starts_quadratic": max(quadratic_errors, default=math.inf) <= tolerances["max_quadratic_scaling_relative_error"],
        "correction_magnitude_nondecreasing": all(
            correction_magnitudes[index + 1] >= correction_magnitudes[index]
            for index in range(max(0, len(correction_magnitudes) - 1))
        ),
        "smallest_positive_row_exposes_cancellation": (
            bool(positive_cancellation)
            and positive_cancellation[0] >= tolerances["min_smallest_cancellation_ratio"]
        ),
    }
    return {
        "policy": "helmholtz_double_layer_low_frequency_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "row_count": len(parsed),
            "positive_row_count": len(positive),
            "maximum_kr_abs": max((row["kr_abs"] for row in positive), default=None),
            "maximum_split_direct_relative_error": max(split_errors, default=None),
            "maximum_series_relative_error": max(series_errors, default=None),
            "maximum_quadratic_scaling_relative_error": max(quadratic_errors, default=None),
            "smallest_positive_cancellation_ratio": positive_cancellation[0] if positive_cancellation else None,
            "expected_laplace_value": expected_laplace,
            "quadratic_reference": quadratic_reference,
        },
        "tolerances": {"min_positive_rows": int(min_positive_rows), **tolerances},
        "notes": [
            "the source-normal double-layer correction begins at order (k*r)^2",
            "check the regular correction separately because total-kernel agreement can hide cancellation",
            "the sign of the Laplace part depends on the source-normal orientation",
        ],
    }
