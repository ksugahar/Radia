"""Induction and electrostatic artifact identity checks for v55."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


INDUCTION = "induction_skin_depth_jouleloss_complexfield_frequency_conductor_owner_identity"
CAPACITANCE = "electrostatic_capacitance_charge_voltage_energy_symmetry_owner_identity"


def _digest(value: object) -> bool:
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generations(row: Mapping[str, object], *fields: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(field) == generation for field in fields)


def _number(value: object, *, positive: bool = False, nonnegative: bool = False) -> bool:
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        return False
    number = float(value)
    return (
        math.isfinite(number)
        and (not positive or number > 0.0)
        and (not nonnegative or number >= 0.0)
    )


def _close(left: object, right: object) -> bool:
    return _number(left) and _number(right) and math.isclose(
        float(left), float(right), rel_tol=1.0e-10, abs_tol=1.0e-12
    )


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get(
        "accepted_result_sha256"
    ) == row.get("result_sha256")


def _induction_ok(row: Mapping[str, object]) -> bool:
    frequency = row.get("frequency_hz")
    conductivity = row.get("conductivity_s_m")
    permeability = row.get("permeability_h_m")
    skin_depth = row.get("skin_depth_m")
    field = row.get("complex_magnetic_field_a_m")
    expected_depth = (
        math.sqrt(
            2.0
            / (
                2.0
                * math.pi
                * float(frequency)
                * float(permeability)
                * float(conductivity)
            )
        )
        if all(_number(value, positive=True) for value in (frequency, conductivity, permeability))
        else math.nan
    )
    field_ok = (
        isinstance(field, Mapping)
        and set(field) == {"real", "imag"}
        and all(_number(value) for value in field.values())
    )
    names = (
        "frequency_hz",
        "conductivity_s_m",
        "permeability_h_m",
        "skin_depth_m",
        "complex_magnetic_field_a_m",
        "joule_loss_w",
        "conductor_owner",
    )
    return (
        _generations(
            row,
            "skin_generation",
            "field_generation",
            "loss_generation",
            "frequency_generation",
            "material_generation",
            "owner_generation",
            "result_generation",
        )
        and _number(skin_depth, positive=True)
        and _close(skin_depth, expected_depth)
        and field_ok
        and _number(row.get("joule_loss_w"), nonnegative=True)
        and str(row.get("conductor_owner") or "").startswith("conductor:")
        and all(row.get("result_" + name) == row.get(name) for name in names)
        and _result(row)
    )


def _vector(value: object, length: int) -> list[float] | None:
    if (
        not isinstance(value, Sequence)
        or isinstance(value, (str, bytes))
        or len(value) != length
        or not all(_number(item) for item in value)
    ):
        return None
    return [float(item) for item in value]


def _matrix(value: object, size: int) -> list[list[float]] | None:
    if (
        not isinstance(value, Sequence)
        or isinstance(value, (str, bytes))
        or len(value) != size
    ):
        return None
    rows = [_vector(row, size) for row in value]
    return None if any(row is None for row in rows) else rows


def _capacitance_ok(row: Mapping[str, object]) -> bool:
    order = row.get("conductor_order")
    if (
        not isinstance(order, Sequence)
        or isinstance(order, (str, bytes))
        or len(order) < 2
        or len(order) != len(set(order))
        or not all(isinstance(item, str) and item.startswith("conductor:") for item in order)
    ):
        return False
    size = len(order)
    matrix = _matrix(row.get("capacitance_matrix_f"), size)
    voltage = _vector(row.get("voltage_v"), size)
    charge = _vector(row.get("charge_c"), size)
    if matrix is None or voltage is None or charge is None:
        return False
    symmetry_ok = all(
        math.isclose(matrix[i][j], matrix[j][i], rel_tol=1.0e-10, abs_tol=1.0e-15)
        for i in range(size)
        for j in range(size)
    )
    expected_charge = [
        sum(matrix[i][j] * voltage[j] for j in range(size)) for i in range(size)
    ]
    charge_ok = all(_close(observed, expected) for observed, expected in zip(charge, expected_charge))
    expected_energy = 0.5 * sum(voltage[i] * charge[i] for i in range(size))
    names = (
        "conductor_order",
        "capacitance_matrix_f",
        "voltage_v",
        "charge_c",
        "stored_energy_j",
        "solution_owner",
    )
    return (
        _generations(
            row,
            "capacitance_generation",
            "charge_generation",
            "voltage_generation",
            "energy_generation",
            "symmetry_generation",
            "owner_generation",
            "result_generation",
        )
        and symmetry_ok
        and charge_ok
        and _number(row.get("stored_energy_j"), nonnegative=True)
        and _close(row.get("stored_energy_j"), expected_energy)
        and str(row.get("solution_owner") or "").startswith("solution:")
        and all(row.get("result_" + name) == row.get(name) for name in names)
        and _result(row)
    )


def validate_public_identity(identity: object) -> dict[str, bool]:
    if not isinstance(identity, Mapping):
        return {}
    checks: dict[str, bool] = {}
    induction = identity.get(INDUCTION)
    capacitance = identity.get(CAPACITANCE)
    if induction is not None:
        checks["v55_induction_skin_field_loss_frequency_conductor_owner"] = (
            isinstance(induction, Mapping) and _induction_ok(induction)
        )
    if capacitance is not None:
        checks["v55_electrostatic_capacitance_charge_energy_symmetry_owner"] = (
            isinstance(capacitance, Mapping) and _capacitance_ok(capacitance)
        )
    return checks
