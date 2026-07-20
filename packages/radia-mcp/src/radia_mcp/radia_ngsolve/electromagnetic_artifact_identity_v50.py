"""Harmonic-circuit and iron-loss artifact identity checks for v50."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


CIRCUIT = "complex_current_phasor_circuit_series_parallel_turns_depth_owner_identity"
IRON = "iron_loss_steinmetz_frequency_flux_waveform_component_owner_identity"


def _digest(value: object) -> bool:
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generations(row: Mapping[str, object], *fields: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(field) == generation for field in fields)


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _finite_vector(value: object, length: int | None = None) -> bool:
    return (
        isinstance(value, Sequence)
        and not isinstance(value, (str, bytes))
        and (length is None or len(value) == length)
        and bool(value)
        and all(isinstance(item, (int, float)) and not isinstance(item, bool) and math.isfinite(float(item)) for item in value)
    )


def _circuit_ok(row: Mapping[str, object]) -> bool:
    current = row.get("current_phasor_a")
    paths = row.get("parallel_paths")
    turns = row.get("turns")
    depth = row.get("depth_m")
    return (
        _generations(
            row,
            "current_generation",
            "connection_generation",
            "turn_generation",
            "depth_generation",
            "owner_generation",
            "result_generation",
        )
        and _finite_vector(current, 2)
        and row.get("result_current_phasor_a") == current
        and row.get("winding_connection") in {"series", "parallel"}
        and row.get("result_winding_connection") == row.get("winding_connection")
        and isinstance(paths, int)
        and not isinstance(paths, bool)
        and paths > 0
        and row.get("result_parallel_paths") == paths
        and isinstance(turns, int)
        and not isinstance(turns, bool)
        and turns > 0
        and row.get("result_turns") == turns
        and isinstance(depth, (int, float))
        and not isinstance(depth, bool)
        and math.isfinite(float(depth))
        and float(depth) > 0.0
        and row.get("result_depth_m") == depth
        and str(row.get("circuit_owner") or "").startswith("circuit:")
        and row.get("result_circuit_owner") == row.get("circuit_owner")
        and _result(row)
    )


def _iron_ok(row: Mapping[str, object]) -> bool:
    coefficients = row.get("steinmetz_coefficients")
    frequency = row.get("frequency_hz")
    waveform = row.get("flux_density_waveform_t")
    components = row.get("loss_components")
    coefficient_ok = (
        isinstance(coefficients, Mapping)
        and set(coefficients) == {"kh", "ke", "alpha"}
        and all(isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value)) and float(value) > 0.0 for value in coefficients.values())
        and float(coefficients["alpha"]) > 1.0
    )
    component_ok = (
        isinstance(components, Mapping)
        and set(components) == {"hysteresis_w", "eddy_w", "excess_w", "total_w"}
        and all(isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value)) and float(value) >= 0.0 for value in components.values())
        and math.isclose(
            float(components["total_w"]),
            float(components["hysteresis_w"]) + float(components["eddy_w"]) + float(components["excess_w"]),
            rel_tol=1e-12,
            abs_tol=1e-12,
        )
    )
    waveform_ok = (
        _finite_vector(waveform)
        and len(waveform) >= 5
        and math.isclose(float(waveform[0]), float(waveform[-1]), abs_tol=1e-12)
        and min(float(value) for value in waveform) < 0.0 < max(float(value) for value in waveform)
    )
    return (
        _generations(
            row,
            "coefficient_generation",
            "frequency_generation",
            "waveform_generation",
            "component_generation",
            "owner_generation",
            "result_generation",
        )
        and coefficient_ok
        and row.get("result_steinmetz_coefficients") == coefficients
        and isinstance(frequency, (int, float))
        and math.isfinite(float(frequency))
        and float(frequency) > 0.0
        and row.get("result_frequency_hz") == frequency
        and waveform_ok
        and row.get("result_flux_density_waveform_t") == waveform
        and component_ok
        and row.get("result_loss_components") == components
        and str(row.get("material_owner") or "").startswith("material:")
        and row.get("result_material_owner") == row.get("material_owner")
        and _result(row)
    )


def validate_public_identity(identity: object) -> dict[str, bool]:
    if not isinstance(identity, Mapping):
        return {}
    checks: dict[str, bool] = {}
    circuit = identity.get(CIRCUIT)
    iron = identity.get(IRON)
    if circuit is not None:
        checks["v50_complex_current_connection_turns_depth_circuit_owner"] = isinstance(circuit, Mapping) and _circuit_ok(circuit)
    if iron is not None:
        checks["v50_steinmetz_frequency_waveform_components_material_owner"] = isinstance(iron, Mapping) and _iron_ok(iron)
    return checks
