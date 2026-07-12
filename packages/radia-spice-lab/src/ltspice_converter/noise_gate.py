"""Solver-neutral thermal-noise identities for an RC low-pass pair."""

from __future__ import annotations

import math
from collections.abc import Mapping


def _relative(left: float, right: float) -> float:
    return abs(left - right) / max(abs(left), abs(right), 1.0e-300)


def rc_thermal_noise_psd_gate(
    summary: Mapping[str, object],
    *,
    maximum_density_relative_error: float = 1.0e-3,
    maximum_rms_relative_error: float = 1.0e-3,
    maximum_measure_relative_error: float = 5.0e-3,
    maximum_cross_case_relative_error: float = 1.0e-3,
) -> dict[str, object]:
    """Gate spectral density, finite-band RMS, and noise-measure semantics."""
    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    tolerances = {
        "maximum_density_relative_error": float(maximum_density_relative_error),
        "maximum_rms_relative_error": float(maximum_rms_relative_error),
        "maximum_measure_relative_error": float(maximum_measure_relative_error),
        "maximum_cross_case_relative_error": float(maximum_cross_case_relative_error),
    }
    if not all(math.isfinite(value) and value >= 0.0 for value in tolerances.values()):
        raise ValueError("tolerances must be finite and nonnegative")

    resistance = float(summary.get("resistance_ohm", math.nan))
    temperature = float(summary.get("temperature_k", math.nan))
    frequency_start = float(summary.get("frequency_start_hz", math.nan))
    frequency_stop = float(summary.get("frequency_stop_hz", math.nan))
    if not (
        math.isfinite(resistance)
        and resistance > 0.0
        and math.isfinite(temperature)
        and temperature > 0.0
        and math.isfinite(frequency_start)
        and math.isfinite(frequency_stop)
        and 0.0 <= frequency_start < frequency_stop
    ):
        raise ValueError("R, T, and frequency band must be finite and physical")

    cases = summary.get("cases")
    if not isinstance(cases, list) or len(cases) != 2:
        raise ValueError("cases must contain exactly two records")
    parsed = []
    for index, case in enumerate(cases):
        if not isinstance(case, Mapping):
            raise ValueError(f"cases[{index}] must be a mapping")
        row = {
            "capacitance": float(case.get("capacitance_f", math.nan)),
            "points": int(case.get("point_count", 0)),
            "density": float(case.get("low_frequency_density", math.nan)),
            "analytic_density": float(case.get("analytic_low_frequency_density", math.nan)),
            "numeric_rms": float(case.get("numeric_psd_integrated_rms", math.nan)),
            "analytic_rms": float(case.get("analytic_finite_band_rms", math.nan)),
            "noise_measure": float(case.get("noise_measure_integrated", math.nan)),
            "direct_density_integral": float(case.get("direct_density_integral", math.nan)),
        }
        if not all(math.isfinite(value) for value in row.values()):
            raise ValueError(f"cases[{index}] contains non-finite values")
        if row["capacitance"] <= 0.0 or row["points"] < 100:
            raise ValueError(f"cases[{index}] requires positive C and a resolved sweep")
        row["density_error"] = _relative(row["density"], row["analytic_density"])
        row["rms_error"] = _relative(row["numeric_rms"], row["analytic_rms"])
        row["measure_error"] = _relative(row["noise_measure"], row["numeric_rms"])
        row["direct_to_rms_ratio"] = row["direct_density_integral"] / max(
            row["numeric_rms"], 1.0e-300
        )
        parsed.append(row)

    numeric_ratio = parsed[1]["numeric_rms"] / parsed[0]["numeric_rms"]
    analytic_ratio = parsed[1]["analytic_rms"] / parsed[0]["analytic_rms"]
    cross_error = _relative(numeric_ratio, analytic_ratio)
    units = summary.get("units")
    if not isinstance(units, Mapping):
        raise ValueError("units must be a mapping")
    checks = {
        "noise_analysis_context_explicit": summary.get("analysis") == "noise",
        "psd_square_integration_explicit": summary.get("rms_integration")
        == "sqrt_integral_density_squared_df",
        "noise_measure_context_semantics_explicit": summary.get("measure_integ_semantics")
        == "noise_rms_not_ordinary_integral",
        "spectral_and_integrated_units_distinct": units.get("spectral_density")
        == "V/sqrt(Hz)"
        and units.get("rms_noise") == "V"
        and units.get("ordinary_density_integral") == "V*sqrt(Hz)",
        "density_matches_rc_thermal_identity": max(row["density_error"] for row in parsed)
        <= tolerances["maximum_density_relative_error"],
        "finite_band_rms_matches_analytic": max(row["rms_error"] for row in parsed)
        <= tolerances["maximum_rms_relative_error"],
        "noise_measure_matches_psd_rms": max(row["measure_error"] for row in parsed)
        <= tolerances["maximum_measure_relative_error"],
        "ordinary_density_integral_is_not_rms": min(
            row["direct_to_rms_ratio"] for row in parsed
        )
        >= 10.0,
        "capacitance_pair_ratio_closes": cross_error
        <= tolerances["maximum_cross_case_relative_error"],
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "rc_thermal_noise_psd_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "metrics": {
            "case_errors": [
                {
                    "density_relative_error": row["density_error"],
                    "rms_relative_error": row["rms_error"],
                    "measure_relative_error": row["measure_error"],
                    "ordinary_integral_to_rms_ratio": row["direct_to_rms_ratio"],
                }
                for row in parsed
            ],
            "numeric_rms_ratio": numeric_ratio,
            "analytic_rms_ratio": analytic_ratio,
            "cross_case_relative_error": cross_error,
        },
        "tolerances": tolerances,
        "lesson": (
            "Noise density has units V/sqrt(Hz), so RMS noise is sqrt(integral e_n^2 df). "
            "A noise-analysis INTEG measure may implement this RMS operation rather than "
            "ordinary waveform integration; preserve the analysis context and units."
        ),
    }
