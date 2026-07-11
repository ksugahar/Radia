"""One-port consistency gates for S, V/I, impedance, and power exports."""

from __future__ import annotations

import math
from bisect import bisect_left
from typing import Any


def _complex(value: Any, name: str) -> complex:
    if isinstance(value, dict):
        parsed = complex(float(value["real"]), float(value["imag"]))
    elif isinstance(value, (list, tuple)) and len(value) == 2:
        parsed = complex(float(value[0]), float(value[1]))
    else:
        parsed = complex(value)
    if not math.isfinite(parsed.real) or not math.isfinite(parsed.imag):
        raise ValueError(f"{name} must be finite")
    return parsed


def _relative_error(actual: complex | float, expected: complex | float) -> float:
    return abs(actual - expected) / max(abs(actual), abs(expected), 1.0e-30)


def one_port_vi_s_impedance_gate(
    summary: dict[str, Any],
    *,
    max_impedance_relative_error: float = 1.0e-12,
    max_power_relative_error: float = 1.0e-6,
    max_passivity_excess: float = 2.0e-4,
    max_frequency_match_relative_error: float = 1.0e-12,
    min_sample_count: int = 5,
) -> dict[str, Any]:
    """Validate one-port wave and terminal-variable representations."""
    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")
    rows = summary.get("rows") or []
    power_rows = summary.get("power_rows") or []
    if not isinstance(rows, list) or not isinstance(power_rows, list):
        raise ValueError("rows and power_rows must be lists")
    stimulated_power = float(summary.get("stimulated_power_w", math.nan))
    if not math.isfinite(stimulated_power) or stimulated_power <= 0.0:
        raise ValueError("stimulated_power_w must be finite and positive")

    parsed = []
    parse_errors = []
    for index, row in enumerate(rows):
        try:
            frequency = float(row["frequency_hz"])
            s11 = _complex(row["s11"], "s11")
            zref = _complex(row["zref_ohm"], "zref_ohm")
            voltage = _complex(row["voltage_v"], "voltage_v")
            current = _complex(row["current_a"], "current_a")
            if not math.isfinite(frequency) or frequency <= 0.0 or abs(current) == 0.0:
                raise ValueError("invalid frequency or zero current")
            if abs(1.0 - s11) == 0.0:
                raise ValueError("S11=1 makes impedance transform singular")
            z_vi = voltage / current
            z_s = zref * (1.0 + s11) / (1.0 - s11)
            terminal_power = 0.5 * (voltage * current.conjugate()).real
            wave_power = stimulated_power * (1.0 - abs(s11) ** 2)
            parsed.append(
                {
                    "frequency_hz": frequency,
                    "s11": s11,
                    "zref": zref,
                    "z_vi": z_vi,
                    "z_s": z_s,
                    "terminal_power": terminal_power,
                    "wave_power": wave_power,
                }
            )
        except (KeyError, TypeError, ValueError) as exc:
            parse_errors.append(f"row {index}: {exc}")

    frequencies = [row["frequency_hz"] for row in parsed]
    increasing = len(parsed) == len(rows) and all(
        right > left for left, right in zip(frequencies, frequencies[1:])
    )
    impedance_errors = [
        _relative_error(row["z_vi"], row["z_s"]) for row in parsed
    ]
    power_errors = [
        _relative_error(row["terminal_power"], row["wave_power"]) for row in parsed
    ]

    dense_frequencies = [row["frequency_hz"] for row in parsed]
    accepted_errors = []
    balance_errors = []
    frequency_match_errors = []
    sparse_errors = []
    for index, row in enumerate(power_rows):
        try:
            frequency = float(row["frequency_hz"])
            accepted = float(row["accepted_power_w"])
            stimulated = float(row["stimulated_power_w"])
            balance = float(row["balance_magnitude"])
            insertion = bisect_left(dense_frequencies, frequency)
            candidates = [
                position
                for position in (insertion - 1, insertion)
                if 0 <= position < len(parsed)
            ]
            if not candidates:
                raise ValueError("no dense frequency candidate")
            position = min(
                candidates,
                key=lambda candidate: abs(dense_frequencies[candidate] - frequency),
            )
            dense = parsed[position]
            match_error = abs(dense["frequency_hz"] - frequency) / max(
                abs(frequency), 1.0
            )
            if match_error > float(max_frequency_match_relative_error):
                raise ValueError("sparse frequency is off the dense grid")
            frequency_match_errors.append(match_error)
            accepted_errors.append(
                _relative_error(
                    accepted,
                    stimulated * (1.0 - abs(dense["s11"]) ** 2),
                )
            )
            balance_errors.append(abs(balance - abs(dense["s11"])))
            if abs(stimulated - stimulated_power) > 1.0e-12:
                sparse_errors.append(f"power row {index} stimulated power drift")
        except (KeyError, TypeError, ValueError):
            sparse_errors.append(f"power row {index} is incomplete or off-grid")

    checks = {
        "rows_parsed": not parse_errors and len(parsed) == len(rows),
        "sample_count_sufficient": len(rows) >= int(min_sample_count),
        "frequency_strictly_increases": increasing,
        "reference_impedance_positive_real": bool(parsed)
        and all(row["zref"].real > 0.0 and abs(row["zref"].imag) <= 1.0e-12 for row in parsed),
        "vi_matches_s_impedance_transform": bool(impedance_errors)
        and max(impedance_errors) <= float(max_impedance_relative_error),
        "terminal_power_matches_wave_power": bool(power_errors)
        and max(power_errors) <= float(max_power_relative_error),
        "sparse_power_rows_parsed_and_aligned": bool(power_rows)
        and not sparse_errors
        and len(accepted_errors) == len(power_rows),
        "accepted_power_matches_s11": bool(accepted_errors)
        and max(accepted_errors) <= float(max_power_relative_error),
        "balance_matches_s11_magnitude": bool(balance_errors)
        and max(balance_errors) <= 1.0e-12,
        "passivity_with_numerical_slack": bool(parsed)
        and max(abs(row["s11"]) for row in parsed) <= 1.0 + float(max_passivity_excess),
    }
    issues = parse_errors + sparse_errors + [name for name, ok in checks.items() if not ok]
    return {
        "policy": "one_port_vi_s_impedance_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "metrics": {
            "sample_count": len(rows),
            "power_sample_count": len(power_rows),
            "maximum_impedance_relative_error": max(impedance_errors, default=math.inf),
            "maximum_terminal_wave_power_relative_error": max(power_errors, default=math.inf),
            "maximum_accepted_power_relative_error": max(accepted_errors, default=math.inf),
            "maximum_balance_absolute_error": max(balance_errors, default=math.inf),
            "maximum_frequency_match_relative_error": max(
                frequency_match_errors, default=math.inf
            ),
            "maximum_s11_magnitude": max((abs(row["s11"]) for row in parsed), default=math.inf),
        },
        "tolerances": {
            "max_impedance_relative_error": float(max_impedance_relative_error),
            "max_power_relative_error": float(max_power_relative_error),
            "max_passivity_excess": float(max_passivity_excess),
            "max_frequency_match_relative_error": float(
                max_frequency_match_relative_error
            ),
            "min_sample_count": int(min_sample_count),
        },
        "notes": [
            "Use peak phasors: terminal real power is 0.5*Re(V*conj(I)).",
            "Small negative accepted power can arise when a numerically lossless one-port has |S11| slightly above one; gate the excess explicitly.",
            "Sparse power rows may differ from the dense S/V/I grid by floating-point roundoff; match by a recorded relative frequency tolerance.",
        ],
    }
