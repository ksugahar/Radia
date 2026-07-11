"""Solver-neutral checks for a constrained harmonic circuit artifact."""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any


def _finite_float(value: Any, name: str) -> float:
    parsed = float(value)
    if not math.isfinite(parsed):
        raise ValueError(f"{name} must be finite")
    return parsed


def _finite_complex(value: Any, name: str) -> complex:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a real/imag mapping")
    parsed = complex(
        _finite_float(value["real"], f"{name}.real"),
        _finite_float(value["imag"], f"{name}.imag"),
    )
    return parsed


def _relative_span(values: list[float]) -> float:
    scale = max((abs(value) for value in values), default=0.0)
    return (max(values) - min(values)) / scale if values and scale > 0.0 else math.inf


def harmonic_zero_net_circuit_gate(
    summary: dict[str, Any],
    *,
    max_faraday_relative_error: float = 1.0e-10,
    max_zero_net_current_ratio: float = 1.0e-10,
    max_source_integral_relative_span: float = 1.0e-9,
    max_loss_imaginary_fraction: float = 1.0e-10,
    max_mesh_count_relative_span: float = 0.0,
    min_sample_count: int = 3,
) -> dict[str, Any]:
    """Validate phasor identities without imposing a false square law.

    The constrained circuit carries zero net current but may have induced
    voltage and flux linkage. With the ``exp(+j*omega*t)`` convention, these
    quantities obey ``V = +j*omega*lambda``. Nonlinear magnetic material can
    make losses and forces depart strongly from amplitude-squared scaling, so
    those coefficients are reported as diagnostics rather than pass criteria.
    """
    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")
    contract = summary.get("contract") or {}
    rows = summary.get("rows") or []
    if not isinstance(contract, dict):
        raise ValueError("contract must be a mapping")
    if not isinstance(rows, list):
        raise ValueError("rows must be a list")

    frequency_hz = _finite_float(contract.get("frequency_hz", math.nan), "frequency_hz")
    omega = 2.0 * math.pi * frequency_hz
    parsed: list[dict[str, Any]] = []
    parse_errors: list[str] = []
    for index, row in enumerate(rows):
        if not isinstance(row, dict):
            parse_errors.append(f"row {index} is not a mapping")
            continue
        try:
            parsed.append(
                {
                    "source_amplitude": _finite_float(
                        row["source_amplitude"], f"row {index}.source_amplitude"
                    ),
                    "source_current_a": _finite_float(
                        row["source_current_a"], f"row {index}.source_current_a"
                    ),
                    "constrained_current": _finite_complex(
                        row["constrained_circuit_current_a"],
                        f"row {index}.constrained_circuit_current_a",
                    ),
                    "voltage": _finite_complex(
                        row["circuit_voltage_v"], f"row {index}.circuit_voltage_v"
                    ),
                    "flux": _finite_complex(
                        row["circuit_flux_linkage_wb_turn"],
                        f"row {index}.circuit_flux_linkage_wb_turn",
                    ),
                    "losses": [
                        _finite_complex(value, f"row {index}.conductive_losses_w")
                        for value in row["conductive_losses_w"]
                    ],
                    "dc_forces": [
                        _finite_float(value, f"row {index}.dc_force_components_n")
                        for value in row["dc_force_components_n"]
                    ],
                    "two_x_forces": [
                        _finite_complex(value, f"row {index}.two_x_force_phasors_n")
                        for value in row["two_x_force_phasors_n"]
                    ],
                    "node_count": _finite_float(
                        row["node_count"], f"row {index}.node_count"
                    ),
                    "element_count": _finite_float(
                        row["element_count"], f"row {index}.element_count"
                    ),
                }
            )
        except (KeyError, TypeError, ValueError) as exc:
            parse_errors.append(f"row {index}: {exc}")

    amplitudes = [row["source_amplitude"] for row in parsed]
    source_currents = [row["source_current_a"] for row in parsed]
    source_coefficients = [
        current / amplitude
        for current, amplitude in zip(source_currents, amplitudes, strict=True)
        if amplitude > 0.0
    ]
    faraday_errors: list[float] = []
    zero_net_ratios: list[float] = []
    loss_imaginary_fractions: list[float] = []
    for row in parsed:
        expected_voltage = 1j * omega * row["flux"]
        faraday_errors.append(
            abs(row["voltage"] - expected_voltage)
            / max(abs(row["voltage"]), abs(expected_voltage), 1.0e-30)
        )
        zero_net_ratios.append(
            abs(row["constrained_current"]) / max(abs(row["source_current_a"]), 1.0e-30)
        )
        for loss in row["losses"]:
            loss_imaginary_fractions.append(abs(loss.imag) / max(abs(loss.real), 1.0e-30))

    nodes = [row["node_count"] for row in parsed]
    elements = [row["element_count"] for row in parsed]
    component_counts_stable = bool(parsed) and all(
        len(row["losses"]) == len(parsed[0]["losses"])
        and len(row["dc_forces"]) == len(parsed[0]["dc_forces"])
        and len(row["two_x_forces"]) == len(parsed[0]["two_x_forces"])
        for row in parsed
    )
    loss_square_spans = (
        [
            _relative_span(
                [
                    row["losses"][component].real / row["source_amplitude"] ** 2
                    for row in parsed
                ]
            )
            for component in range(len(parsed[0]["losses"]))
        ]
        if component_counts_stable
        else []
    )
    dc_force_square_spans = (
        [
            _relative_span(
                [
                    abs(row["dc_forces"][component]) / row["source_amplitude"] ** 2
                    for row in parsed
                ]
            )
            for component in range(len(parsed[0]["dc_forces"]))
        ]
        if component_counts_stable
        else []
    )
    two_x_force_square_spans = (
        [
            _relative_span(
                [
                    abs(row["two_x_forces"][component])
                    / row["source_amplitude"] ** 2
                    for row in parsed
                ]
            )
            for component in range(len(parsed[0]["two_x_forces"]))
        ]
        if component_counts_stable
        else []
    )
    force_components = contract.get("force_components") or []
    checks = {
        "rows_parsed_and_finite": not parse_errors and len(parsed) == len(rows),
        "sample_count_sufficient": len(parsed) >= int(min_sample_count),
        "positive_source_amplitudes_strictly_increase": bool(amplitudes)
        and all(right > left > 0.0 for left, right in zip(amplitudes, amplitudes[1:])),
        "fixed_positive_frequency": frequency_hz > 0.0,
        "phasor_and_faraday_sign_recorded": contract.get("phasor_convention")
        == "exp(+j*omega*t)"
        and contract.get("faraday_identity") == "V=+j*omega*flux_linkage",
        "zero_net_circuit_contract_recorded": contract.get("circuit_constraint")
        == "zero_net_current",
        "nonlinear_scaling_not_overclaimed": contract.get("material_response")
        == "nonlinear"
        and contract.get("scaling_interpretation") == "diagnostic_only",
        "dc_and_two_x_force_components_separated": set(force_components)
        == {"dc_time_average", "two_x_phasor"}
        and bool(parsed)
        and component_counts_stable
        and all(row["dc_forces"] and row["two_x_forces"] for row in parsed),
        "force_scope_not_overclaimed": contract.get("force_method")
        == "lorentz_volume_current_density"
        and contract.get("force_scope")
        == "current_density_contribution_not_total_ferromagnetic_force",
        "faraday_voltage_flux_identity": bool(faraday_errors)
        and max(faraday_errors) <= float(max_faraday_relative_error),
        "zero_net_current_constraint": bool(zero_net_ratios)
        and max(zero_net_ratios) <= float(max_zero_net_current_ratio),
        "source_current_integral_scales_with_source": bool(source_coefficients)
        and _relative_span(source_coefficients) <= float(max_source_integral_relative_span),
        "conductive_losses_positive_and_real": bool(loss_imaginary_fractions)
        and all(loss.real > 0.0 for row in parsed for loss in row["losses"])
        and max(loss_imaginary_fractions) <= float(max_loss_imaginary_fraction),
        "mesh_inventory_positive_and_stable": bool(parsed)
        and min(nodes) > 0.0
        and min(elements) > 0.0
        and _relative_span(nodes) <= float(max_mesh_count_relative_span)
        and _relative_span(elements) <= float(max_mesh_count_relative_span),
    }
    issues = parse_errors + [name for name, ok in checks.items() if not ok]
    return {
        "policy": "harmonic_zero_net_circuit_faraday_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "metrics": {
            "sample_count": len(parsed),
            "max_faraday_relative_error": max(faraday_errors, default=math.inf),
            "max_zero_net_current_ratio": max(zero_net_ratios, default=math.inf),
            "source_current_per_amplitude_relative_span": _relative_span(source_coefficients),
            "max_loss_imaginary_fraction": max(loss_imaginary_fractions, default=math.inf),
            "loss_per_amplitude_squared_relative_spans": loss_square_spans,
            "max_loss_per_amplitude_squared_relative_span": max(
                loss_square_spans, default=math.inf
            ),
            "dc_force_per_amplitude_squared_relative_spans": dc_force_square_spans,
            "max_dc_force_per_amplitude_squared_relative_span": max(
                dc_force_square_spans, default=math.inf
            ),
            "two_x_force_per_amplitude_squared_relative_spans": two_x_force_square_spans,
            "max_two_x_force_per_amplitude_squared_relative_span": max(
                two_x_force_square_spans, default=math.inf
            ),
            "node_count_relative_span": _relative_span(nodes),
            "element_count_relative_span": _relative_span(elements),
        },
        "tolerances": {
            "max_faraday_relative_error": float(max_faraday_relative_error),
            "max_zero_net_current_ratio": float(max_zero_net_current_ratio),
            "max_source_integral_relative_span": float(max_source_integral_relative_span),
            "max_loss_imaginary_fraction": float(max_loss_imaginary_fraction),
            "max_mesh_count_relative_span": float(max_mesh_count_relative_span),
            "min_sample_count": int(min_sample_count),
        },
        "notes": [
            "A zero-net-current conductor can still have induced voltage and flux linkage.",
            "Match the Faraday sign to the recorded phasor convention before comparing solvers.",
            "Keep harmonic DC/time-average force separate from the complex two-times-frequency force.",
            "A Lorentz volume integral in conducting magnetic material is a current-density contribution, not the total ferromagnetic-body force.",
            "For nonlinear magnetic material, amplitude-squared loss and force coefficients are diagnostics, not acceptance laws.",
        ],
    }
