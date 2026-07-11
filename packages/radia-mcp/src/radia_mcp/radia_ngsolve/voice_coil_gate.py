"""Solver-neutral force/current/flux gate for permanent-magnet voice coils."""

from __future__ import annotations

import math
from typing import Any, Iterable, Mapping


def voice_coil_force_flux_sweep_gate(
    rows: Iterable[Mapping[str, Any]],
    *,
    max_zero_force_relative: float = 0.01,
    max_odd_residual_relative: float = 0.10,
    max_force_constant_relative_span: float = 0.12,
) -> dict[str, Any]:
    """Gate a symmetric DC current sweep without requiring an ideal linear material."""

    records = [dict(row) for row in rows]
    if len(records) < 5 or len(records) % 2 == 0:
        raise ValueError("rows must contain an odd number of at least five samples")
    normalized = []
    for index, row in enumerate(records):
        try:
            current = float(row["current_a"])
            circuit_current = float(row["circuit_current_a"])
            force = float(row["axial_force_n"])
            flux_raw = row["flux_linkage_wb_turn"]
            flux = float(flux_raw[0] if isinstance(flux_raw, (list, tuple)) else flux_raw)
            nodes = int(row["node_count"])
            elements = int(row["element_count"])
        except (KeyError, TypeError, ValueError, IndexError) as exc:
            raise ValueError(f"row {index} is malformed") from exc
        values = (current, circuit_current, force, flux)
        if not all(math.isfinite(value) for value in values) or nodes <= 0 or elements <= 0:
            raise ValueError(f"row {index} contains invalid values")
        normalized.append({
            "current_a": current,
            "circuit_current_a": circuit_current,
            "axial_force_n": force,
            "flux_linkage_wb_turn": flux,
            "node_count": nodes,
            "element_count": elements,
        })

    currents = [row["current_a"] for row in normalized]
    forces = [row["axial_force_n"] for row in normalized]
    fluxes = [row["flux_linkage_wb_turn"] for row in normalized]
    midpoint = len(normalized) // 2
    zero_force = forces[midpoint]
    peak_force = max(abs(value) for value in forces)
    force_scale = max(peak_force, math.ulp(1.0))
    force_constants = [
        (force - zero_force) / current
        for current, force in zip(currents, forces)
        if current != 0.0
    ]
    mean_force_constant = sum(force_constants) / len(force_constants)
    force_constant_span = (
        (max(force_constants) - min(force_constants)) / abs(mean_force_constant)
        if mean_force_constant != 0.0 else math.inf
    )
    odd_residuals = [
        abs(forces[left] + forces[-left - 1] - 2.0 * zero_force) / force_scale
        for left in range(midpoint)
    ]
    incremental_flux_slopes = [
        (right_flux - left_flux) / (right_current - left_current)
        for left_current, right_current, left_flux, right_flux in zip(
            currents, currents[1:], fluxes, fluxes[1:]
        )
    ]
    checks = {
        "currents_strictly_increase": all(a < b for a, b in zip(currents, currents[1:])),
        "current_grid_is_symmetric_with_zero": (
            currents[midpoint] == 0.0
            and all(abs(currents[left] + currents[-left - 1]) <= 1.0e-12 for left in range(midpoint))
        ),
        "circuit_current_matches_requested": all(
            abs(row["circuit_current_a"] - row["current_a"]) <= 1.0e-12 for row in normalized
        ),
        "axial_force_strictly_increases": all(a < b for a, b in zip(forces, forces[1:])),
        "zero_current_force_is_small": abs(zero_force) / force_scale <= float(max_zero_force_relative),
        "force_odd_symmetry_envelope_ok": max(odd_residuals) <= float(max_odd_residual_relative),
        "force_constants_positive": all(value > 0.0 for value in force_constants),
        "force_constant_span_ok": force_constant_span <= float(max_force_constant_relative_span),
        "incremental_flux_linkage_positive": all(value > 0.0 for value in incremental_flux_slopes),
        "mesh_inventory_invariant": (
            len({row["node_count"] for row in normalized}) == 1
            and len({row["element_count"] for row in normalized}) == 1
        ),
    }
    return {
        "policy": "voice_coil_force_flux_sweep_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "zero_current_force_n": zero_force,
            "zero_force_relative": abs(zero_force) / force_scale,
            "maximum_odd_residual_relative": max(odd_residuals),
            "mean_force_constant_n_per_a": mean_force_constant,
            "force_constant_relative_span": force_constant_span,
            "minimum_incremental_flux_linkage_h": min(incremental_flux_slopes),
            "node_count": normalized[0]["node_count"],
            "element_count": normalized[0]["element_count"],
        },
        "lesson": (
            "A permanent-magnet voice coil with nonlinear iron need not be perfectly odd-linear. Gate a "
            "bounded odd-symmetry envelope, positive force constant, positive incremental flux linkage, "
            "small zero-current offset, and invariant mesh together."
        ),
    }
