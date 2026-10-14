"""Cross-formulation gate for source-normalized single-loop fields."""

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


def _complex(value: object, name: str) -> complex:
    if (
        not isinstance(value, Sequence)
        or isinstance(value, (str, bytes))
        or len(value) != 2
    ):
        raise ValueError(f"{name} must be [real, imag]")
    return complex(_finite(value[0], f"{name}[0]"), _finite(value[1], f"{name}[1]"))


def _records(value: object, name: str, *, minimum: int = 3) -> list[Mapping[str, object]]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError(f"{name} must be an array")
    records = list(value)
    if len(records) < minimum or any(not isinstance(row, Mapping) for row in records):
        raise ValueError(f"{name} must contain at least {minimum} objects")
    return records


def _relative(first: float, second: float) -> float:
    return abs(first - second) / max(0.5 * (abs(first) + abs(second)), 1.0e-300)


def single_loop_source_normalized_field_gate(
    summary: Mapping[str, object],
    *,
    maximum_transfer_relative_gap: float = 0.05,
    maximum_component_reconstruction_error: float = 1.0e-10,
    maximum_transverse_to_axial_ratio: float = 5.0e-3,
    maximum_power_error_over_stimulated: float = 1.0e-8,
) -> dict[str, Any]:
    """Compare loop fields only after normalizing by solved source current."""
    if not isinstance(summary, Mapping):
        raise ValueError("summary must be an object")
    contract = summary.get("model_contract")
    if not isinstance(contract, Mapping):
        raise ValueError("model_contract must be an object")
    routes = _records(summary.get("routes"), "routes", minimum=2)
    if len(routes) != 2:
        raise ValueError("routes must contain exactly two formulations")
    limits = {
        "maximum_transfer_relative_gap": _finite(
            maximum_transfer_relative_gap,
            "maximum_transfer_relative_gap",
            positive=True,
        ),
        "maximum_component_reconstruction_error": _finite(
            maximum_component_reconstruction_error,
            "maximum_component_reconstruction_error",
            positive=True,
        ),
        "maximum_transverse_to_axial_ratio": _finite(
            maximum_transverse_to_axial_ratio,
            "maximum_transverse_to_axial_ratio",
            positive=True,
        ),
        "maximum_power_error_over_stimulated": _finite(
            maximum_power_error_over_stimulated,
            "maximum_power_error_over_stimulated",
            positive=True,
        ),
    }

    route_metrics = []
    frequency_axes = []
    reconstruction_errors = []
    transverse_ratios = []
    power_errors = []
    for route_index, route in enumerate(routes):
        name = str(route.get("name") or "").strip()
        if not name:
            raise ValueError(f"routes[{route_index}].name is required")
        field_rows = _records(
            route.get("field_rows"), f"routes[{route_index}].field_rows"
        )
        power_rows = _records(
            route.get("power_rows"), f"routes[{route_index}].power_rows"
        )
        frequencies = []
        transfers = []
        parsed_fields = []
        for row_index, row in enumerate(field_rows):
            prefix = f"routes[{route_index}].field_rows[{row_index}]"
            frequency = _finite(row.get("frequency_hz"), f"{prefix}.frequency_hz", positive=True)
            current = _complex(row.get("current_a"), f"{prefix}.current_a")
            h_components = row.get("h_components_a_per_m")
            if (
                not isinstance(h_components, Sequence)
                or isinstance(h_components, (str, bytes))
                or len(h_components) != 3
            ):
                raise ValueError(f"{prefix}.h_components_a_per_m must contain x, y, z")
            components = [
                _complex(value, f"{prefix}.h_components_a_per_m[{index}]")
                for index, value in enumerate(h_components)
            ]
            h_magnitude = _finite(
                row.get("h_magnitude_a_per_m"),
                f"{prefix}.h_magnitude_a_per_m",
                positive=True,
            )
            current_magnitude = abs(current)
            if current_magnitude <= 0.0:
                raise ValueError(f"{prefix}.current_a must be nonzero")
            reconstructed = math.sqrt(sum(abs(value) ** 2 for value in components))
            reconstruction_error = _relative(h_magnitude, reconstructed)
            transverse = math.sqrt(abs(components[0]) ** 2 + abs(components[1]) ** 2)
            transverse_ratio = transverse / max(abs(components[2]), 1.0e-300)
            transfer = h_magnitude / current_magnitude
            frequencies.append(frequency)
            transfers.append(transfer)
            reconstruction_errors.append(reconstruction_error)
            transverse_ratios.append(transverse_ratio)
            parsed_fields.append({
                "frequency_hz": frequency,
                "source_normalized_h_per_m": transfer,
                "component_reconstruction_relative_error": reconstruction_error,
                "transverse_to_axial_ratio": transverse_ratio,
            })
        parsed_power = []
        for row_index, row in enumerate(power_rows):
            prefix = f"routes[{route_index}].power_rows[{row_index}]"
            frequency = _finite(row.get("frequency_hz"), f"{prefix}.frequency_hz", positive=True)
            s11 = _complex(row.get("s11"), f"{prefix}.s11")
            stimulated = _finite(
                row.get("stimulated_power_w"), f"{prefix}.stimulated_power_w", positive=True
            )
            accepted = _finite(row.get("accepted_power_w"), f"{prefix}.accepted_power_w")
            expected = stimulated * (1.0 - abs(s11) ** 2)
            error = abs(accepted - expected) / stimulated
            power_errors.append(error)
            parsed_power.append({
                "frequency_hz": frequency,
                "one_port_power_error_over_stimulated": error,
            })
        frequency_axes.append(frequencies)
        route_metrics.append({
            "name": name,
            "formulation": str(route.get("formulation") or ""),
            "field_rows": parsed_fields,
            "power_rows": parsed_power,
            "transfers": transfers,
            "field_frequency_increasing": all(
                right > left for left, right in zip(frequencies, frequencies[1:])
            ),
        })

    transfer_gaps = [
        _relative(first, second)
        for first, second in zip(
            route_metrics[0]["transfers"], route_metrics[1]["transfers"]
        )
    ]
    checks = {
        "single_loop_cross_formulation_contract": contract
        == {
            "physics": "harmonic_maxwell",
            "single_turn_loop": True,
            "one_port_per_route": True,
            "same_probe_location": True,
            "raw_port_phase_comparable": False,
        },
        "distinct_named_formulations": len({row["name"] for row in route_metrics}) == 2
        and all(row["formulation"] for row in route_metrics),
        "same_strictly_increasing_field_frequency_axis": frequency_axes[0]
        == frequency_axes[1]
        and all(row["field_frequency_increasing"] for row in route_metrics),
        "field_components_reconstruct_magnitude": max(reconstruction_errors)
        <= limits["maximum_component_reconstruction_error"],
        "center_field_is_axially_dominant": max(transverse_ratios)
        <= limits["maximum_transverse_to_axial_ratio"],
        "source_normalized_field_transfer_agrees": max(transfer_gaps)
        <= limits["maximum_transfer_relative_gap"],
        "one_port_power_identity_holds": max(power_errors)
        <= limits["maximum_power_error_over_stimulated"],
        "raw_port_phase_is_not_used_for_cross_formulation_claim": contract.get(
            "raw_port_phase_comparable"
        )
        is False,
    }
    for row in route_metrics:
        row.pop("transfers")
        row.pop("field_frequency_increasing")
    return {
        "policy": "single_loop_source_normalized_field_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "field_frequency_count": len(frequency_axes[0]),
            "maximum_source_normalized_transfer_relative_gap": max(transfer_gaps),
            "maximum_component_reconstruction_relative_error": max(reconstruction_errors),
            "maximum_transverse_to_axial_ratio": max(transverse_ratios),
            "maximum_one_port_power_error_over_stimulated": max(power_errors),
            "routes": route_metrics,
        },
        "tolerances": limits,
        "lesson": (
            "When port conventions differ between formulations, compare the field "
            "transfer normalized by the solved source current. Raw reflection phase "
            "is not a cross-formulation observable until port calibration is aligned."
        ),
    }
