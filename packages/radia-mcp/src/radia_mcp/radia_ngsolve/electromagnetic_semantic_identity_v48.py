"""Semantic ownership checks for incremental magnetic and electrostatic artifacts."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


INCREMENTAL = "incremental_permeability_frozen_bias_harmonic_phasor_operating_point_owner_identity"
ELECTROSTATIC = "electrostatic_capacitance_charge_energy_voltage_sweep_identity"


def _digest(value: object) -> bool:
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _generations_ok(row: Mapping[str, object], fields: tuple[str, ...]) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(field) == generation for field in fields)


def _finite_vector(value: object, length: int | None = None) -> bool:
    return (
        isinstance(value, Sequence)
        and not isinstance(value, (str, bytes))
        and (length is None or len(value) == length)
        and bool(value)
        and all(isinstance(item, (int, float)) and math.isfinite(float(item)) for item in value)
    )


def _result_ok(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _tangent_ok(value: object) -> bool:
    if not isinstance(value, Mapping) or not value:
        return False
    for matrix in value.values():
        if not isinstance(matrix, list) or len(matrix) != 2 or not all(_finite_vector(row, 2) for row in matrix):
            return False
        a, b = (float(item) for item in matrix[0])
        c, d = (float(item) for item in matrix[1])
        if not math.isclose(b, c, rel_tol=1.0e-12, abs_tol=1.0e-12) or a <= 0.0 or a * d - b * c <= 0.0:
            return False
    return True


def _incremental_ok(row: Mapping[str, object]) -> bool:
    phasor = row.get("harmonic_phasor_a")
    tangent = row.get("material_tangent_h_per_m")
    owner = str(row.get("operating_point_owner") or "")
    return (
        _generations_ok(
            row,
            ("bias_generation", "phasor_generation", "tangent_generation", "operating_point_generation", "result_generation"),
        )
        and _digest(row.get("frozen_bias_sha256"))
        and row.get("result_frozen_bias_sha256") == row.get("frozen_bias_sha256")
        and _finite_vector(phasor, 2)
        and row.get("result_harmonic_phasor_a") == phasor
        and _tangent_ok(tangent)
        and row.get("result_material_tangent_h_per_m") == tangent
        and owner.startswith("operating-point:")
        and row.get("result_operating_point_owner") == owner
        and _result_ok(row)
    )


def _electrostatic_ok(row: Mapping[str, object]) -> bool:
    voltage = row.get("voltage_v")
    charge = row.get("charge_c")
    energy = row.get("field_energy_j")
    capacitance = row.get("capacitance_f")
    owners = row.get("conductor_owner_rows")
    valid_rows = (
        _finite_vector(voltage)
        and _finite_vector(charge, len(voltage))
        and _finite_vector(energy, len(voltage))
        and isinstance(owners, list)
        and len(owners) == len(voltage)
        and bool(owners)
        and len(set(owners)) == 1
        and all(str(owner).startswith("conductor:") for owner in owners)
    )
    closure = False
    if valid_rows and isinstance(capacitance, (int, float)) and math.isfinite(float(capacitance)) and float(capacitance) > 0.0:
        closure = all(
            math.isclose(float(q), float(capacitance) * float(v), rel_tol=1.0e-9, abs_tol=1.0e-18)
            and math.isclose(float(w), 0.5 * float(q) * float(v), rel_tol=1.0e-9, abs_tol=1.0e-18)
            for v, q, w in zip(voltage, charge, energy, strict=True)
        )
    return (
        _generations_ok(
            row,
            ("voltage_generation", "charge_generation", "energy_generation", "conductor_generation", "result_generation"),
        )
        and valid_rows
        and closure
        and row.get("result_voltage_v") == voltage
        and row.get("result_charge_c") == charge
        and row.get("result_field_energy_j") == energy
        and row.get("result_capacitance_f") == capacitance
        and row.get("result_conductor_owner_rows") == owners
        and _result_ok(row)
    )


def validate_public_identity(identity: object) -> dict[str, bool]:
    if not isinstance(identity, Mapping):
        return {}
    checks: dict[str, bool] = {}
    incremental = identity.get(INCREMENTAL)
    electrostatic = identity.get(ELECTROSTATIC)
    if incremental is not None:
        checks["v48_incremental_bias_phasor_tangent_owner"] = isinstance(incremental, Mapping) and _incremental_ok(incremental)
    if electrostatic is not None:
        checks["v48_electrostatic_charge_energy_voltage_owner_closure"] = isinstance(electrostatic, Mapping) and _electrostatic_ok(electrostatic)
    return checks
