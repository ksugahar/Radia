"""Solver-neutral network and field-monitor artifact lineage checks."""

from __future__ import annotations

import math
from collections.abc import Mapping


SMATRIX = "v47_public_smatrix_port_order_reference_plane_network_owner_mismatch"
FIELD = "v47_public_field_energy_loss_q_monitor_frequency_row_key_mismatch"


def _sha(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _generation(row: Mapping[str, object], names: tuple[str, ...]) -> bool:
    value = str(row.get("generation") or "")
    return bool(value) and all(row.get(name) == value for name in names)


def _digest(row: Mapping[str, object]) -> bool:
    return _sha(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _finite_sequence(value: object, length: int) -> bool:
    return (
        isinstance(value, list)
        and len(value) == length
        and all(isinstance(item, (int, float)) and math.isfinite(float(item)) for item in value)
    )


def _smatrix_ok(row: Mapping[str, object]) -> bool:
    ports = row.get("port_order")
    planes = row.get("reference_plane_m")
    normalization = row.get("network_normalization_ohm")
    return (
        _generation(
            row,
            ("port_generation", "reference_plane_generation", "normalization_generation", "network_generation", "result_generation"),
        )
        and isinstance(ports, list)
        and len(ports) >= 1
        and all(isinstance(port, str) and port for port in ports)
        and len(ports) == len(set(ports))
        and row.get("result_port_order") == ports
        and isinstance(planes, Mapping)
        and set(planes) == set(ports)
        and all(isinstance(value, (int, float)) and math.isfinite(float(value)) for value in planes.values())
        and row.get("result_reference_plane_m") == planes
        and isinstance(normalization, (int, float))
        and math.isfinite(float(normalization))
        and float(normalization) > 0.0
        and row.get("result_network_normalization_ohm") == normalization
        and str(row.get("network_owner") or "").startswith("network:")
        and row.get("result_network_owner") == row.get("network_owner")
        and _digest(row)
    )


def _field_ok(row: Mapping[str, object]) -> bool:
    keys = row.get("frequency_row_keys")
    if not isinstance(keys, list) or not keys or len(keys) != len(set(keys)):
        return False
    count = len(keys)
    energy = row.get("field_energy_j")
    loss = row.get("loss_w")
    q_factor = row.get("q_factor")
    return (
        _generation(
            row,
            ("monitor_generation", "frequency_generation", "energy_generation", "loss_generation", "q_generation", "result_generation"),
        )
        and all(isinstance(key, str) and key for key in keys)
        and row.get("result_frequency_row_keys") == keys
        and str(row.get("monitor_identity") or "").startswith("monitor:")
        and row.get("result_monitor_identity") == row.get("monitor_identity")
        and _finite_sequence(energy, count)
        and row.get("result_field_energy_j") == energy
        and _finite_sequence(loss, count)
        and all(float(value) >= 0.0 for value in loss)
        and row.get("result_loss_w") == loss
        and _finite_sequence(q_factor, count)
        and all(float(value) >= 0.0 for value in q_factor)
        and row.get("result_q_factor") == q_factor
        and str(row.get("result_owner") or "").startswith("result:")
        and row.get("result_result_owner") == row.get("result_owner")
        and _digest(row)
    )


def validate_public_v47_identity(payload: object) -> dict[str, bool]:
    if not isinstance(payload, Mapping):
        return {}
    rows = [row for row in (payload.get("runs") or []) if isinstance(row, Mapping)]
    checks: dict[str, bool] = {}
    smatrices = [row[SMATRIX] for row in rows if SMATRIX in row]
    fields = [row[FIELD] for row in rows if FIELD in row]
    if smatrices:
        checks["network_v47_smatrix_port_plane_normalization_owner"] = (
            len(smatrices) == len(rows) and all(isinstance(row, Mapping) and _smatrix_ok(row) for row in smatrices)
        )
    if fields:
        checks["network_v47_field_energy_loss_q_monitor_rows"] = (
            len(fields) == len(rows) and all(isinstance(row, Mapping) and _field_ok(row) for row in fields)
        )
    return checks
