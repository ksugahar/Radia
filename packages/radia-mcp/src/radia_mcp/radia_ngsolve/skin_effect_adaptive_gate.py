"""Solver-neutral skin-effect port and adaptive-convergence validation."""
from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any


def _complex(value: Mapping[str, Any], name: str) -> complex:
    try:
        result = complex(float(value["real"]), float(value["imag"]))
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"{name} must provide finite real and imag values") from exc
    if not math.isfinite(result.real) or not math.isfinite(result.imag):
        raise ValueError(f"{name} must be finite")
    return result


def skin_effect_adaptive_energy_loss_gate(
    *, frequency_hz: float, current: Mapping[str, Any], voltage: Mapping[str, Any],
    impedance: Mapping[str, Any], power: Mapping[str, Any], flux_linkage: Mapping[str, Any],
    total_energy_j: float, total_loss_w: float, adaptive_rows: Sequence[Mapping[str, Any]],
    identity_rtol: float = 1.0e-4, flux_voltage_rtol: float = 1.0e-3,
    energy_tail_rtol: float = 2.0e-3, final_loss_change_rtol: float = 6.0e-2,
    final_adaptive_error_max: float = 5.0e-2,
) -> dict[str, object]:
    """Gate current-port identities and the slower convergence of skin loss."""
    f = float(frequency_hz); energy = float(total_energy_j); loss = float(total_loss_w)
    if not all(math.isfinite(v) and v > 0 for v in (f, energy, loss)):
        raise ValueError("frequency, energy, and loss must be finite and positive")
    i = _complex(current, "current"); v = _complex(voltage, "voltage")
    z = _complex(impedance, "impedance"); p = _complex(power, "power")
    psi = _complex(flux_linkage, "flux_linkage")
    if abs(i) == 0 or abs(v) == 0:
        raise ValueError("current and voltage must be nonzero")
    rows = [dict(row) for row in adaptive_rows]
    if len(rows) < 4:
        raise ValueError("adaptive_rows must contain at least four passes")
    cells = [int(row["mesh_cells"]) for row in rows]
    energies = [float(row["energy_j"]) for row in rows]
    losses = [float(row["loss_w"]) for row in rows]
    errors = [float(row["adaptive_error"]) for row in rows]
    if not all(math.isfinite(x) and x >= 0 for x in energies + losses + errors):
        raise ValueError("adaptive observables must be finite and nonnegative")
    omega = 2 * math.pi * f
    rel = lambda a, b: abs(a - b) / max(abs(b), 1.0e-300)
    energy_from_reactance = abs(z.imag) * abs(i) ** 2 / (2 * omega)
    energy_tail_spread = (max(energies[-3:]) - min(energies[-3:])) / energies[-1]
    final_loss_change = rel(losses[-1], losses[-2])
    checks = {
        "voltage_equals_impedance_current": rel(v, z * i) <= identity_rtol,
        "power_equals_voltage_current_conjugate": rel(p, v * i.conjugate()) <= identity_rtol,
        "loss_matches_oriented_real_power": rel(loss, abs(p.real)) <= identity_rtol,
        "energy_matches_reactive_impedance": rel(energy, energy_from_reactance) <= identity_rtol,
        "flux_linkage_matches_port_voltage": rel(v, -1j * omega * psi) <= flux_voltage_rtol,
        "mesh_cells_strictly_increase": all(a < b for a, b in zip(cells, cells[1:])),
        "energy_tail_stable": energy_tail_spread <= energy_tail_rtol,
        "skin_loss_change_resolved": final_loss_change <= final_loss_change_rtol,
        "final_loss_matches_adaptive_history": rel(loss, losses[-1]) <= identity_rtol,
        "final_adaptive_error_bounded": errors[-1] <= final_adaptive_error_max,
        "final_error_below_first_reported": errors[-1] < errors[0],
    }
    return {
        "policy": "skin_effect_adaptive_energy_loss_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks, "issues": [k for k, ok in checks.items() if not ok],
        "metrics": {"adaptive_pass_count": len(rows), "final_mesh_cells": cells[-1],
                    "energy_tail_relative_spread": energy_tail_spread,
                    "final_loss_relative_change": final_loss_change,
                    "final_adaptive_error": errors[-1]},
        "lesson": "Magnetic energy may converge before conductor surface loss. Gate port identities and both adaptive histories; solver completion alone is not skin-effect validation.",
    }
