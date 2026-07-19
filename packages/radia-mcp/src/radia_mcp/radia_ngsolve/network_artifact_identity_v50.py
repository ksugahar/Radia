"""Solver-neutral Floquet and time-domain port identity checks for v50."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


FLOQUET = "periodic_floquet_phase_lattice_vector_mode_normalization_boundary_owner_identity"
TD_PORT = "time_domain_port_pulse_bandwidth_deembedding_reference_plane_owner_identity"


def _sha(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generation(row: Mapping[str, object], names: tuple[str, ...]) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(name) == generation for name in names)


def _digest(row: Mapping[str, object]) -> bool:
    return _sha(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _finite(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value))


def _vector(value: object, size: int) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes)) and len(value) == size and all(_finite(item) for item in value)


def _lattice(value: object) -> bool:
    if not isinstance(value, list) or len(value) != 2 or not all(_vector(item, 3) for item in value):
        return False
    first, second = value
    cross = [
        float(first[1]) * float(second[2]) - float(first[2]) * float(second[1]),
        float(first[2]) * float(second[0]) - float(first[0]) * float(second[2]),
        float(first[0]) * float(second[1]) - float(first[1]) * float(second[0]),
    ]
    return math.sqrt(sum(item * item for item in cross)) > 0.0


def _normalization(value: object) -> bool:
    return (
        isinstance(value, Mapping)
        and value.get("kind") == "unit-power"
        and _finite(value.get("value_w"))
        and math.isclose(float(value["value_w"]), 1.0, rel_tol=0.0, abs_tol=1.0e-12)
    )


def _floquet_ok(row: Mapping[str, object]) -> bool:
    phase = row.get("floquet_phase_deg")
    lattice = row.get("lattice_vectors_m")
    normalization = row.get("mode_normalization")
    owner = str(row.get("boundary_owner") or "")
    return (
        _generation(row, ("phase_generation", "lattice_generation", "mode_generation", "boundary_generation", "result_generation"))
        and _vector(phase, 2)
        and all(-180.0 <= float(item) <= 180.0 for item in phase)
        and row.get("result_floquet_phase_deg") == phase
        and _lattice(lattice)
        and row.get("result_lattice_vectors_m") == lattice
        and _normalization(normalization)
        and row.get("result_mode_normalization") == normalization
        and owner.startswith("boundary:")
        and row.get("result_boundary_owner") == owner
        and _digest(row)
    )


def _bandwidth(value: object) -> bool:
    return _vector(value, 2) and 0.0 < float(value[0]) < float(value[1])


def _td_port_ok(row: Mapping[str, object]) -> bool:
    bandwidth = row.get("pulse_bandwidth_hz")
    distance = row.get("deembedding_distance_m")
    plane = str(row.get("reference_plane") or "")
    owner = str(row.get("port_owner") or "")
    return (
        _generation(row, ("pulse_generation", "deembedding_generation", "reference_generation", "port_generation", "result_generation"))
        and _bandwidth(bandwidth)
        and row.get("result_pulse_bandwidth_hz") == bandwidth
        and _finite(distance)
        and float(distance) >= 0.0
        and row.get("result_deembedding_distance_m") == distance
        and plane.startswith("port-plane:")
        and row.get("result_reference_plane") == plane
        and owner.startswith("port:")
        and row.get("result_port_owner") == owner
        and _digest(row)
    )


def validate_public_v50_identity(payload: object) -> dict[str, bool]:
    """Validate optional v50 periodic and time-domain port records."""
    if not isinstance(payload, Mapping):
        return {}
    rows = [row for row in (payload.get("runs") or []) if isinstance(row, Mapping)]
    checks: dict[str, bool] = {}
    floquet = [row[FLOQUET] for row in rows if FLOQUET in row]
    ports = [row[TD_PORT] for row in rows if TD_PORT in row]
    if floquet:
        checks["network_v50_floquet_phase_lattice_normalization_boundary_owner"] = len(floquet) == len(rows) and all(isinstance(row, Mapping) and _floquet_ok(row) for row in floquet)
    if ports:
        checks["network_v50_time_port_bandwidth_deembedding_reference_owner"] = len(ports) == len(rows) and all(isinstance(row, Mapping) and _td_port_ok(row) for row in ports)
    return checks
