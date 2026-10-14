"""Solver-neutral Faraday gate for two-winding harmonic sweeps."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any


def _finite(value: object, name: str, *, positive: bool = False) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(parsed) or (positive and parsed <= 0.0):
        requirement = "positive and finite" if positive else "finite"
        raise ValueError(f"{name} must be {requirement}")
    return parsed


def _complex_pair(value: object, name: str) -> complex:
    if (
        not isinstance(value, Sequence)
        or isinstance(value, (str, bytes))
        or len(value) != 2
    ):
        raise ValueError(f"{name} must be [real, imag]")
    return complex(
        _finite(value[0], f"{name}[0]"),
        _finite(value[1], f"{name}[1]"),
    )


def _rows(value: object) -> list[Mapping[str, object]]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError("rows must be an array")
    rows = list(value)
    if len(rows) < 3 or any(not isinstance(row, Mapping) for row in rows):
        raise ValueError("rows must contain at least three objects")
    return rows


def two_winding_frequency_faraday_gate(
    summary: Mapping[str, object],
    *,
    maximum_faraday_relative_error: float = 1.0e-3,
    maximum_linkage_per_turn_relative_gap: float = 0.15,
) -> dict[str, Any]:
    """Validate ``R response = -j omega flux_linkage`` over a sweep.

    Each row contains two windings. ``response`` is the complex winding
    response normalized by resistance and ``flux_linkage_Wb_turn`` is the
    complex linked flux reported with the winding-turn factor applied.
    """
    if not isinstance(summary, Mapping):
        raise ValueError("summary must be an object")
    contract = summary.get("model_contract")
    if not isinstance(contract, Mapping):
        raise ValueError("model_contract must be an object")
    rows = _rows(summary.get("rows"))
    faraday_limit = _finite(
        maximum_faraday_relative_error,
        "maximum_faraday_relative_error",
        positive=True,
    )
    linkage_limit = _finite(
        maximum_linkage_per_turn_relative_gap,
        "maximum_linkage_per_turn_relative_gap",
        positive=True,
    )

    frequencies: list[float] = []
    faraday_errors: list[float] = []
    linkage_gaps: list[float] = []
    parsed_rows: list[dict[str, object]] = []
    for row_index, row in enumerate(rows):
        frequency = _finite(
            row.get("frequency_hz"), f"rows[{row_index}].frequency_hz", positive=True
        )
        windings = row.get("windings")
        if (
            not isinstance(windings, Sequence)
            or isinstance(windings, (str, bytes))
            or len(windings) != 2
            or any(not isinstance(winding, Mapping) for winding in windings)
        ):
            raise ValueError(f"rows[{row_index}].windings must contain two objects")
        frequencies.append(frequency)
        per_turn: list[complex] = []
        winding_metrics = []
        for winding_index, winding in enumerate(windings):
            prefix = f"rows[{row_index}].windings[{winding_index}]"
            turns = _finite(winding.get("turns"), f"{prefix}.turns", positive=True)
            resistance = _finite(
                winding.get("resistance_ohm"),
                f"{prefix}.resistance_ohm",
                positive=True,
            )
            response = _complex_pair(winding.get("response"), f"{prefix}.response")
            flux = _complex_pair(
                winding.get("flux_linkage_Wb_turn"),
                f"{prefix}.flux_linkage_Wb_turn",
            )
            voltage = resistance * response
            expected = -1j * 2.0 * math.pi * frequency * flux
            error = abs(voltage - expected) / max(abs(expected), 1.0e-300)
            faraday_errors.append(error)
            per_turn.append(flux / turns)
            winding_metrics.append({
                "turns": turns,
                "resistance_ohm": resistance,
                "faraday_relative_error": error,
            })
        linkage_gap = abs(per_turn[0] - per_turn[1]) / max(
            0.5 * (abs(per_turn[0]) + abs(per_turn[1])), 1.0e-300
        )
        linkage_gaps.append(linkage_gap)
        parsed_rows.append({
            "frequency_hz": frequency,
            "windings": winding_metrics,
            "linkage_per_turn_relative_gap": linkage_gap,
        })

    checks = {
        "harmonic_two_winding_contract": contract
        == {
            "physics": "harmonic_magnetics",
            "two_windings": True,
            "passive_secondary": True,
            "complex_phasors": True,
        },
        "positive_strictly_increasing_frequency_axis": all(
            right > left for left, right in zip(frequencies, frequencies[1:])
        ),
        "faraday_identity_holds_for_both_windings": max(faraday_errors)
        <= faraday_limit,
        "linkage_per_turn_is_consistent": max(linkage_gaps) <= linkage_limit,
    }
    return {
        "policy": "two_winding_frequency_faraday_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "frequency_count": len(frequencies),
            "frequency_min_hz": min(frequencies),
            "frequency_max_hz": max(frequencies),
            "maximum_faraday_relative_error": max(faraday_errors),
            "maximum_linkage_per_turn_relative_gap": max(linkage_gaps),
            "rows": parsed_rows,
        },
        "tolerances": {
            "maximum_faraday_relative_error": faraday_limit,
            "maximum_linkage_per_turn_relative_gap": linkage_limit,
        },
        "lesson": (
            "Cross-check complex winding response against linked flux at every "
            "frequency; a plausible magnitude alone does not establish phasor "
            "sign, resistance normalization, or Faraday consistency."
        ),
    }
