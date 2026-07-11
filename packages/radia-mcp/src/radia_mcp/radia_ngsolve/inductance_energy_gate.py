"""Energy and analytic-reference gate for a coupled-coil excitation."""
from __future__ import annotations

import math


def inductance_energy_mutual_gate(
    *,
    self_inductance: float,
    energy_inductance: float,
    mutual_inductance: float,
    analytic_mutual_inductance: float,
    inductance_unit: str,
    max_energy_relative_error: float = 0.01,
    max_mutual_relative_error: float = 0.05,
) -> dict:
    """Gate L=2W/I^2 and one-direction mutual-inductance evidence."""

    values = (
        self_inductance,
        energy_inductance,
        mutual_inductance,
        analytic_mutual_inductance,
        max_energy_relative_error,
        max_mutual_relative_error,
    )
    if not all(math.isfinite(float(value)) for value in values):
        raise ValueError("all inductances and tolerances must be finite")
    if any(float(value) <= 0.0 for value in values[:4]):
        raise ValueError("all inductance measurements must be positive")
    if any(float(value) < 0.0 for value in values[4:]):
        raise ValueError("tolerances must be nonnegative")

    unit = str(inductance_unit).strip()
    recognized_unit = unit in {"H", "mH", "uH", "nH", "pH"}
    energy_error = abs(self_inductance - energy_inductance) / energy_inductance
    mutual_error = abs(mutual_inductance - analytic_mutual_inductance) / analytic_mutual_inductance
    checks = {
        "inductance_unit_explicit": recognized_unit,
        "self_inductance_positive": self_inductance > 0.0,
        "mutual_inductance_positive": mutual_inductance > 0.0,
        "energy_identity_closes": energy_error <= max_energy_relative_error,
        "mutual_analytic_reference_closes": mutual_error <= max_mutual_relative_error,
        "mutual_below_excited_self_inductance": abs(mutual_inductance) < self_inductance,
    }
    return {
        "policy": "inductance_energy_mutual_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "energy_relative_error": energy_error,
            "mutual_analytic_relative_error": mutual_error,
            "mutual_to_self_ratio": abs(mutual_inductance) / self_inductance,
        },
        "conventions": {
            "inductance_unit": unit,
            "energy_identity": "L=2*Wmag/I^2",
            "reciprocity_status": "not_tested_by_one_direction_gate",
        },
        "notes": [
            "This one-direction gate does not prove M12=M21; run two independent excitation solves for reciprocity.",
            "Evaluate each expected result expression after the study rebuild before accepting derived values.",
        ],
    }
