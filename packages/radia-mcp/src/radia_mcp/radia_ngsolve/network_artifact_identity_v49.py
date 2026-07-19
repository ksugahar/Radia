"""Solver-neutral modal-port and particle-wake artifact identity checks."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


MODAL = "modal_port_mode_index_degeneracy_polarization_phase_impedance_mesh_owner_identity"
WAKE = "particle_wakefield_bunch_charge_time_reference_monitor_normalization_owner_identity"


def _sha(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _generation(row: Mapping[str, object], names: tuple[str, ...]) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(name) == generation for name in names)


def _digest(row: Mapping[str, object]) -> bool:
    return _sha(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _finite(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value))


def _finite_sequence(value: object, *, minimum: int = 1) -> bool:
    return (
        isinstance(value, Sequence)
        and not isinstance(value, (str, bytes))
        and len(value) >= minimum
        and all(_finite(item) for item in value)
    )


def _strictly_increasing(value: object) -> bool:
    return _finite_sequence(value, minimum=2) and all(
        float(left) < float(right) for left, right in zip(value, value[1:])
    )


def _orthonormal_2d_basis(value: object) -> bool:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or len(value) != 2:
        return False
    if not all(_finite_sequence(vector, minimum=2) and len(vector) == 2 for vector in value):
        return False
    first, second = value
    first_norm = sum(float(component) ** 2 for component in first)
    second_norm = sum(float(component) ** 2 for component in second)
    dot = sum(float(left) * float(right) for left, right in zip(first, second))
    return math.isclose(first_norm, 1.0, abs_tol=1.0e-12) and math.isclose(
        second_norm, 1.0, abs_tol=1.0e-12
    ) and math.isclose(dot, 0.0, abs_tol=1.0e-12)


def _modal_ok(row: Mapping[str, object]) -> bool:
    mode_index = row.get("mode_index")
    phase = row.get("polarization_phase_deg")
    impedance = row.get("reference_impedance_ohm")
    basis = row.get("degeneracy_basis")
    return (
        _generation(
            row,
            (
                "mode_generation",
                "degeneracy_generation",
                "polarization_generation",
                "impedance_generation",
                "mesh_generation",
                "port_generation",
                "result_generation",
            ),
        )
        and isinstance(mode_index, int)
        and not isinstance(mode_index, bool)
        and mode_index > 0
        and row.get("result_mode_index") == mode_index
        and _orthonormal_2d_basis(basis)
        and row.get("result_degeneracy_basis") == basis
        and _finite(phase)
        and -180.0 <= float(phase) <= 180.0
        and row.get("result_polarization_phase_deg") == phase
        and _finite(impedance)
        and float(impedance) > 0.0
        and row.get("result_reference_impedance_ohm") == impedance
        and _sha(row.get("mesh_sha256"))
        and row.get("result_mesh_sha256") == row.get("mesh_sha256")
        and str(row.get("port_owner") or "").startswith("port:")
        and row.get("result_port_owner") == row.get("port_owner")
        and _digest(row)
    )


def _wake_ok(row: Mapping[str, object]) -> bool:
    charge = row.get("bunch_charge_c")
    times = row.get("monitor_time_s")
    wake = row.get("wake_v_per_c")
    return (
        _generation(
            row,
            (
                "charge_generation",
                "time_generation",
                "monitor_generation",
                "normalization_generation",
                "result_generation",
            ),
        )
        and _finite(charge)
        and float(charge) > 0.0
        and row.get("result_bunch_charge_c") == charge
        and row.get("time_reference") == row.get("result_time_reference") == "bunch_center"
        and _strictly_increasing(times)
        and row.get("result_monitor_time_s") == times
        and _finite_sequence(wake, minimum=len(times) if isinstance(times, Sequence) else 1)
        and len(wake) == len(times)
        and row.get("result_wake_v_per_c") == wake
        and row.get("normalization") == row.get("result_normalization") == "per_coulomb"
        and str(row.get("result_owner") or "").startswith("result:")
        and row.get("accepted_result_owner") == row.get("result_owner")
        and _digest(row)
    )


def validate_public_v49_identity(payload: object) -> dict[str, bool]:
    if not isinstance(payload, Mapping):
        return {}
    rows = [row for row in (payload.get("runs") or []) if isinstance(row, Mapping)]
    checks: dict[str, bool] = {}
    modal = [row[MODAL] for row in rows if MODAL in row]
    wake = [row[WAKE] for row in rows if WAKE in row]
    if modal:
        checks["network_v49_modal_port_mode_polarization_mesh_owner"] = len(modal) == len(rows) and all(
            isinstance(row, Mapping) and _modal_ok(row) for row in modal
        )
    if wake:
        checks["network_v49_particle_wake_time_normalization_owner"] = len(wake) == len(rows) and all(
            isinstance(row, Mapping) and _wake_ok(row) for row in wake
        )
    return checks
