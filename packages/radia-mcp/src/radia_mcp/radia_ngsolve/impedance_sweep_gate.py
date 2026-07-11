"""Solver-neutral validation for multi-port complex impedance sweeps."""

from __future__ import annotations

import math
from typing import Any


def multiport_impedance_sweep_gate(
    frequency_rows: list[list[float]],
    impedance_real_rows: list[list[float]],
    impedance_imag_rows: list[list[float]],
    *,
    port_ids: list[str] | None = None,
    min_sample_count: int = 5,
    min_frequency_decades: float = 1.0,
    min_impedance_relative_span: float = 1.0e-6,
    passive_real_tolerance: float = 1.0e-9,
    common_grid_relative_tolerance: float = 1.0e-12,
) -> dict[str, Any]:
    """Require common frequency identity, positive-realness and nontrivial data."""

    frequencies = [[float(value) for value in row] for row in frequency_rows]
    real_rows = [[float(value) for value in row] for row in impedance_real_rows]
    imag_rows = [[float(value) for value in row] for row in impedance_imag_rows]
    port_count = len(frequencies)
    ids = list(port_ids) if port_ids is not None else [f"port_{index}" for index in range(port_count)]
    tolerances = (
        float(min_frequency_decades),
        float(min_impedance_relative_span),
        float(passive_real_tolerance),
        float(common_grid_relative_tolerance),
    )
    if int(min_sample_count) < 2:
        raise ValueError("min_sample_count must be at least 2")
    if any(not math.isfinite(value) or value < 0.0 for value in tolerances):
        raise ValueError("tolerances must be finite and nonnegative")
    if port_count == 0:
        raise ValueError("at least one impedance profile is required")
    if len(real_rows) != port_count or len(imag_rows) != port_count or len(ids) != port_count:
        raise ValueError("frequency, impedance and port rows must have equal profile counts")
    if len(set(ids)) != len(ids) or any(not str(value).strip() for value in ids):
        raise ValueError("port_ids must be nonempty and unique")
    for index, (freq, real, imag) in enumerate(zip(frequencies, real_rows, imag_rows)):
        if len(freq) != len(real) or len(freq) != len(imag):
            raise ValueError(f"profile {index} frequency and impedance rows must have equal length")
        if not all(math.isfinite(value) for value in freq + real + imag):
            raise ValueError(f"profile {index} contains a non-finite value")

    reference_grid = frequencies[0]
    sample_counts = [len(row) for row in frequencies]
    grids_match = all(
        len(row) == len(reference_grid)
        and all(
            abs(value - reference) <= common_grid_relative_tolerance * max(abs(reference), 1.0)
            for value, reference in zip(row, reference_grid)
        )
        for row in frequencies[1:]
    )
    strictly_increasing = [
        all(a < b for a, b in zip(row, row[1:])) for row in frequencies
    ]
    positive_frequency = [bool(row) and row[0] > 0.0 for row in frequencies]
    frequency_decades = [
        math.log10(row[-1] / row[0]) if len(row) >= 2 and row[0] > 0.0 else 0.0
        for row in frequencies
    ]
    min_real = [min(row, default=math.inf) for row in real_rows]
    magnitudes = [
        [math.hypot(real, imag) for real, imag in zip(real_row, imag_row)]
        for real_row, imag_row in zip(real_rows, imag_rows)
    ]
    relative_spans = []
    for row in magnitudes:
        scale = max(row, default=0.0)
        relative_spans.append((max(row) - min(row)) / scale if scale > 0.0 else 0.0)

    checks = {
        "port_ids_unique": len(set(ids)) == port_count,
        "sample_count_sufficient": all(count >= int(min_sample_count) for count in sample_counts),
        "common_frequency_grid": grids_match,
        "frequencies_positive": all(positive_frequency),
        "frequencies_strictly_increasing": all(strictly_increasing),
        "frequency_span_sufficient": all(value >= min_frequency_decades for value in frequency_decades),
        "impedance_positive_real_with_tolerance": all(value >= -passive_real_tolerance for value in min_real),
        "impedance_profiles_nonzero": all(max(row, default=0.0) > 0.0 for row in magnitudes),
        "impedance_profiles_nontrivial": all(value >= min_impedance_relative_span for value in relative_spans),
    }
    return {
        "policy": "multiport_impedance_sweep_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "port_count": port_count,
            "port_ids": ids,
            "sample_counts": sample_counts,
            "frequency_decades": frequency_decades,
            "minimum_real_impedance": min_real,
            "impedance_relative_spans": relative_spans,
        },
        "notes": [
            "cross-port comparison requires the same frequency grid before row-wise interpretation",
            "a passive impedance sweep should be positive real within a declared numerical tolerance",
        ],
    }
