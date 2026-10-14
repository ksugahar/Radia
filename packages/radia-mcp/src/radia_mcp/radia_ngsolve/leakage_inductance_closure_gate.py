"""Cross-formulation gate for two-winding leakage inductance."""
from __future__ import annotations

import math


def _pair(value, name: str) -> tuple[float, float]:
    if not isinstance(value, list) or len(value) != 2:
        raise ValueError(f"{name} must contain two values")
    pair = tuple(float(item) for item in value)
    if not all(math.isfinite(item) for item in pair):
        raise ValueError(f"{name} must contain finite values")
    return pair


def leakage_inductance_closure_gate(
    summary: dict,
    *,
    max_reciprocity_relative_error: float = 0.01,
    max_closure_relative_error: float = 0.01,
    max_ampere_turn_relative_error: float = 1.0e-12,
    max_replay_abs_Wb: float = 1.0e-12,
) -> dict:
    """Compare compensated-current energy with a two-source L matrix.

    The matrix convention is ``lambda = L I``.  A compensated current pair
    satisfies ``N1*I1 + N2*I2 = 0``.  Its energy-derived leakage inductance
    must close against the quadratic form of the independently assembled
    inductance matrix.
    """
    if not isinstance(summary, dict):
        raise ValueError("summary must be an object")
    tolerances = (
        max_reciprocity_relative_error,
        max_closure_relative_error,
        max_ampere_turn_relative_error,
        max_replay_abs_Wb,
    )
    if any(value < 0.0 or not math.isfinite(value) for value in tolerances):
        raise ValueError("tolerances must be finite and nonnegative")

    currents = _pair(summary.get("compensated_currents_A"), "compensated_currents_A")
    turns = _pair(summary.get("turns"), "turns")
    flux = _pair(summary.get("compensated_flux_linkage_Wb"), "compensated_flux_linkage_Wb")
    matrix = summary.get("matrix_H")
    if not isinstance(matrix, list) or len(matrix) != 2:
        raise ValueError("matrix_H must be 2x2")
    if any(not isinstance(row, list) or len(row) != 2 for row in matrix):
        raise ValueError("matrix_H must be 2x2")
    l11, m12 = (float(value) for value in matrix[0])
    m21, l22 = (float(value) for value in matrix[1])
    values = (*currents, *turns, *flux, l11, m12, m21, l22)
    if not all(math.isfinite(value) for value in values):
        raise ValueError("all physical inputs must be finite")
    i1, i2 = currents
    n1, n2 = turns
    if i1 == 0.0 or n1 <= 0.0 or n2 <= 0.0 or l22 == 0.0:
        raise ValueError("I1, turns and L22 must be nonzero and turns positive")

    ampere_turn_sum = n1 * i1 + n2 * i2
    ampere_turn_scale = max(abs(n1 * i1), abs(n2 * i2), 1.0e-300)
    ampere_turn_error = abs(ampere_turn_sum) / ampere_turn_scale
    energy = 0.5 * (i1 * flux[0] + i2 * flux[1])
    direct_leakage = 2.0 * energy / (i1 * i1)
    matrix_leakage = (
        l11 * i1 * i1
        + (m12 + m21) * i1 * i2
        + l22 * i2 * i2
    ) / (i1 * i1)
    closure_error = abs(direct_leakage - matrix_leakage) / max(
        abs(direct_leakage), abs(matrix_leakage), 1.0e-300
    )
    reciprocity_error = abs(m12 - m21) / max(abs(m12), abs(m21), 1.0e-300)
    k2 = m12 * m21 / (l11 * l22) if l11 > 0.0 and l22 > 0.0 else math.inf
    short_circuit = l11 - m12 * m21 / l22
    replay = float(summary.get("replay_max_abs_Wb", math.inf))

    checks = {
        "positive_self_inductances": l11 > 0.0 and l22 > 0.0,
        "mutual_terms_have_consistent_sign": m12 == 0.0
        or m21 == 0.0
        or m12 * m21 > 0.0,
        "mutual_reciprocity_within_tolerance": reciprocity_error
        <= max_reciprocity_relative_error,
        "compensated_ampere_turns_close": ampere_turn_error
        <= max_ampere_turn_relative_error,
        "direct_leakage_inductance_positive": direct_leakage > 0.0,
        "matrix_leakage_inductance_positive": matrix_leakage > 0.0,
        "direct_and_matrix_leakage_close": closure_error
        <= max_closure_relative_error,
        "coupling_squared_physical": 0.0 <= k2 < 1.0,
        "short_circuit_inductance_positive_and_below_self": 0.0
        < short_circuit
        < l11,
        "fresh_replay_within_tolerance": math.isfinite(replay)
        and replay <= max_replay_abs_Wb,
    }
    return {
        "policy": "leakage_inductance_closure_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "observables": {
            "compensated_energy_J": energy,
            "leakage_inductance_direct_H": direct_leakage,
            "leakage_inductance_matrix_H": matrix_leakage,
            "short_circuit_inductance_H": short_circuit,
            "coupling_squared": k2,
            "reciprocity_relative_error": reciprocity_error,
            "ampere_turn_relative_error": ampere_turn_error,
            "direct_matrix_relative_error": closure_error,
            "replay_max_abs_Wb": replay,
        },
        "lesson": (
            "Validate leakage inductance by two independent routes: a direct "
            "compensated-current energy solve and a unit-current inductance "
            "matrix. Distinguish this leakage value from the Schur-complement "
            "short-circuit inductance, and gate reciprocity, coupling and replay."
        ),
    }
