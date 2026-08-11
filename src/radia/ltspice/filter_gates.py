"""Public-safe transfer-function gates for second-order active filters."""

from __future__ import annotations

import math
from typing import Any


def sallen_key_filter_family_gate(
    rows: list[dict[str, Any]],
    *,
    max_dc_gain_error: float = 0.01,
    max_minus3db_relative_error: float = 0.01,
    max_ideal_complex_relative_l2: float = 0.01,
    max_high_frequency_gain: float = 0.01,
) -> dict[str, Any]:
    """Gate unity-gain low-pass rows against the ideal two-pole identity."""

    if not isinstance(rows, list) or len(rows) < 2:
        raise ValueError("rows must contain at least two filter variants")
    tolerances = (
        float(max_dc_gain_error),
        float(max_minus3db_relative_error),
        float(max_ideal_complex_relative_l2),
        float(max_high_frequency_gain),
    )
    if any(not math.isfinite(value) or value < 0.0 for value in tolerances):
        raise ValueError("tolerances must be finite and nonnegative")

    normalized = []
    for index, row in enumerate(rows):
        try:
            identifier = str(row["id"]).strip()
            r1, r2 = float(row["R1_ohm"]), float(row["R2_ohm"])
            c1, c2 = float(row["C1_F"]), float(row["C2_F"])
            dc_gain = float(row["dc_gain"])
            peak_gain = float(row["peak_gain"])
            measured_cutoff = float(row["minus3dB_frequency_Hz"])
            high_gain = float(row["gain_at_100k"])
            complex_error = float(row["ideal_complex_relative_l2_to_20k"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"row {index} is missing a numeric filter field") from exc
        values = (r1, r2, c1, c2, dc_gain, peak_gain, measured_cutoff, high_gain, complex_error)
        if not identifier or any(not math.isfinite(value) for value in values):
            raise ValueError(f"row {index} contains an empty id or non-finite value")
        if min(r1, r2, c1, c2, measured_cutoff) <= 0.0:
            raise ValueError(f"row {index} components and cutoff must be positive")
        omega0 = 1.0 / math.sqrt(r1 * r2 * c1 * c2)
        natural_frequency = omega0 / (2.0 * math.pi)
        quality_factor = math.sqrt(r1 * r2 * c1 * c2) / (c2 * (r1 + r2))
        coefficient = 1.0 / (quality_factor * quality_factor) - 2.0
        ideal_cutoff = natural_frequency * math.sqrt(
            (-coefficient + math.sqrt(coefficient * coefficient + 4.0)) / 2.0
        )
        normalized.append({
            "id": identifier,
            "quality_factor": quality_factor,
            "natural_frequency_Hz": natural_frequency,
            "ideal_minus3dB_frequency_Hz": ideal_cutoff,
            "minus3dB_relative_error": abs(measured_cutoff - ideal_cutoff) / ideal_cutoff,
            "dc_gain_error": abs(dc_gain - 1.0),
            "peak_gain": peak_gain,
            "high_frequency_gain": high_gain,
            "ideal_complex_relative_l2": complex_error,
        })

    ids = [row["id"] for row in normalized]
    q_sorted = sorted(normalized, key=lambda row: row["quality_factor"])
    peaking_threshold = 1.0 / math.sqrt(2.0)
    peaking_matches_q = all(
        (row["peak_gain"] > 1.001) if row["quality_factor"] > peaking_threshold + 1.0e-3
        else (row["peak_gain"] <= 1.01)
        for row in normalized
    )
    checks = {
        "variant_ids_unique": len(set(ids)) == len(ids),
        "dc_gain_near_unity": all(row["dc_gain_error"] <= max_dc_gain_error for row in normalized),
        "minus3db_matches_two_pole_identity": all(
            row["minus3dB_relative_error"] <= max_minus3db_relative_error for row in normalized
        ),
        "complex_response_matches_ideal": all(
            row["ideal_complex_relative_l2"] <= max_ideal_complex_relative_l2 for row in normalized
        ),
        "high_frequency_response_is_attenuated": all(
            row["high_frequency_gain"] <= max_high_frequency_gain for row in normalized
        ),
        "peaking_behavior_matches_quality_factor": peaking_matches_q,
        "higher_quality_factor_does_not_reduce_peak": all(
            left["peak_gain"] <= right["peak_gain"] + 1.0e-6
            for left, right in zip(q_sorted, q_sorted[1:])
        ),
    }
    return {
        "schema": "radia-spice-lab.sallen-key-filter-family.v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "ok": all(checks.values()),
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "rows": normalized,
        "notes": [
            "the -3 dB frequency is not generally the natural frequency when Q is not Butterworth",
            "complex-response agreement is stronger than matching one cutoff scalar",
        ],
    }
