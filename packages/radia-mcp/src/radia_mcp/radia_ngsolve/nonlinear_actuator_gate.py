"""Solver-neutral saturation gate for nonlinear axisymmetric actuators."""
from __future__ import annotations
import math
from collections.abc import Mapping


def nonlinear_actuator_saturation_knee_gate(summary: Mapping[str, object]) -> dict[str, object]:
    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    rows = list(summary.get("rows") or [])
    if len(rows) < 5 or not all(isinstance(row, Mapping) for row in rows):
        raise ValueError("rows must contain at least five mappings")
    currents = [float(row["current_a"]) for row in rows]
    flux = [float(row["flux_linkage_wb_turn"]) for row in rows]
    inductance = [float(row["secant_inductance_h"]) for row in rows]
    axial = [float(row["weighted_stress_axial_force_n"]) for row in rows]
    radial = [float(row["weighted_stress_radial_force_n"]) for row in rows]
    force_i2 = [float(row["axial_force_per_i2_n_per_a2"]) for row in rows]
    fixed_b = [float(row["fixed_iron_b_t"]) for row in rows]
    plunger_b = [float(row["plunger_b_t"]) for row in rows]
    energy = [float(row["magnetic_energy_j"]) for row in rows]
    coenergy = [float(row["magnetic_coenergy_j"]) for row in rows]
    nodes = [int(row["node_count"]) for row in rows]
    elements = [int(row["element_count"]) for row in rows]
    values = currents + flux + inductance + axial + radial + force_i2 + fixed_b + plunger_b + energy + coenergy
    l_peak = max(range(len(rows)), key=inductance.__getitem__)
    f_peak = max(range(len(rows)), key=lambda index: abs(force_i2[index]))
    l_tail = inductance[l_peak:]
    f_tail = force_i2[f_peak:]
    energy_gap = abs(coenergy[-1] - energy[-1]) / max(abs(coenergy[-1]), 1.0e-300)
    checks = {
        "all_values_finite": all(math.isfinite(value) for value in values),
        "positive_uniformly_increasing_current": min(currents) > 0.0 and all(b > a for a, b in zip(currents, currents[1:])),
        "flux_and_iron_b_increase": all(b > a for a, b in zip(flux, flux[1:])) and all(b > a for a, b in zip(fixed_b, fixed_b[1:])) and all(b > a for a, b in zip(plunger_b, plunger_b[1:])),
        "shared_nonendpoint_knee": 0 < l_peak == f_peak < len(rows) - 2,
        "secant_inductance_decreases_after_knee": all(b < a for a, b in zip(l_tail, l_tail[1:])),
        "force_per_i2_decreases_after_knee": all(abs(b) < abs(a) for a, b in zip(f_tail, f_tail[1:])),
        "material_saturation_is_substantial": inductance[-1] / inductance[l_peak] < 0.6 and abs(force_i2[-1] / force_i2[f_peak]) < 0.3,
        "force_is_axial_and_single_signed": min(axial) * max(axial) > 0.0 and max(abs(value) for value in radial) <= 1.0e-9 * max(abs(value) for value in axial),
        "nonlinear_energy_coenergy_separate": min(energy) > 0.0 and min(coenergy) > 0.0 and energy_gap > 0.2,
        "mesh_inventory_stable": min(nodes) > 0 and min(elements) > 0 and len(set(nodes)) == 1 and len(set(elements)) == 1,
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "nonlinear_actuator_saturation_knee_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "metrics": {
            "knee_current_a": currents[l_peak],
            "secant_inductance_high_to_peak_ratio": inductance[-1] / inductance[l_peak],
            "force_per_i2_high_to_peak_ratio": abs(force_i2[-1] / force_i2[f_peak]),
            "high_current_energy_coenergy_relative_gap": energy_gap,
        },
        "lesson": "A nonlinear B-H curve may raise permeability before saturation. Detect a shared interior knee in lambda/I and F/I^2, then require both tails to decline; do not demand monotone decline from the first current sample.",
    }
