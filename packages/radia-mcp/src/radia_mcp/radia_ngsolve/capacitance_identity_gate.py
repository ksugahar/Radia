"""Solver-neutral charge/energy identity for a two-conductor capacitor."""
from __future__ import annotations

import math


def two_conductor_capacitance_identity_gate(
    conductor_voltages_v,
    conductor_charges_c,
    stored_energy_j,
    *,
    driven_conductor_index: int = 1,
    planar_depth_m: float | None = None,
    max_capacitance_relative_error: float = 1.0e-5,
    max_charge_balance_relative_error: float = 1.0e-5,
):
    """Check ``C=Q/dV=2W/dV^2`` and equal/opposite conductor charge.

    ``planar_depth_m`` is optional for genuine 3D models and mandatory metadata
    for extruded 2D results.  Passing it makes the returned capacitance an
    absolute farad value rather than an ambiguous per-length quantity.
    """

    voltages = [float(value) for value in conductor_voltages_v]
    charges = [float(value) for value in conductor_charges_c]
    energy = float(stored_energy_j)
    if len(voltages) != 2 or len(charges) != 2:
        raise ValueError("exactly two conductor voltages and charges are required")
    if driven_conductor_index not in (0, 1):
        raise ValueError("driven_conductor_index must be 0 or 1")
    cap_tol = float(max_capacitance_relative_error)
    balance_tol = float(max_charge_balance_relative_error)
    if cap_tol < 0.0 or balance_tol < 0.0:
        raise ValueError("relative tolerances must be nonnegative")
    depth = None if planar_depth_m is None else float(planar_depth_m)
    if depth is not None and (not math.isfinite(depth) or depth <= 0.0):
        raise ValueError("planar_depth_m must be finite and positive when supplied")

    reference_index = 1 - driven_conductor_index
    delta_v = voltages[driven_conductor_index] - voltages[reference_index]
    driven_charge = charges[driven_conductor_index]
    finite = all(math.isfinite(value) for value in voltages + charges + [energy])
    nonzero_voltage = finite and abs(delta_v) > 0.0
    if nonzero_voltage:
        capacitance_charge = driven_charge / delta_v
        capacitance_energy = 2.0 * energy / (delta_v * delta_v)
        cap_rel_error = abs(capacitance_charge - capacitance_energy) / max(
            abs(capacitance_charge), abs(capacitance_energy), 1.0e-300
        )
    else:
        capacitance_charge = math.nan
        capacitance_energy = math.nan
        cap_rel_error = math.inf
    charge_balance_rel_error = abs(sum(charges)) / max(
        abs(charges[0]), abs(charges[1]), 1.0e-300
    )
    checks = {
        "all_finite": finite,
        "nonzero_potential_difference": nonzero_voltage,
        "stored_energy_positive": energy > 0.0,
        "capacitance_from_charge_positive": capacitance_charge > 0.0,
        "capacitance_from_energy_positive": capacitance_energy > 0.0,
        "charge_energy_capacitance_agree": cap_rel_error <= cap_tol,
        "conductor_charge_balance_ok": charge_balance_rel_error <= balance_tol,
    }
    return {
        "policy": "two_conductor_capacitance_identity_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "driven_conductor_index": driven_conductor_index,
        "reference_conductor_index": reference_index,
        "conductor_voltages_V": voltages,
        "conductor_charges_C": charges,
        "potential_difference_V": delta_v,
        "stored_energy_J": energy,
        "planar_depth_m": depth,
        "capacitance_from_charge_F": capacitance_charge,
        "capacitance_from_energy_F": capacitance_energy,
        "capacitance_relative_error": cap_rel_error,
        "charge_balance_relative_error": charge_balance_rel_error,
        "checks": checks,
        "lesson": (
            "A capacitance result is stronger when terminal charge, field energy, "
            "and equal/opposite conductor charge agree; an extruded 2D result must "
            "also carry its physical depth to distinguish F from F/m."
        ),
    }
