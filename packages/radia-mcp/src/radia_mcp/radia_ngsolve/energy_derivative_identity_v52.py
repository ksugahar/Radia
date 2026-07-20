"""Virtual-displacement force and angular-energy torque checks for v52."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

from .energy_derivative_identity_v53 import validate_public_identity as validate_public_v53_identity


VIRTUAL_FORCE = "virtual_displacement_mesh_energy_sign_solution_owner_identity"
MAGNET_TORQUE = "magnet_torque_angular_energy_periodic_unwrap_derivative_owner_identity"


def _digest(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generations(row: Mapping[str, object], *fields: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(field) == generation for field in fields)


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _finite(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value))


def _finite_vector(value: object, length: int | None = None) -> bool:
    return (
        isinstance(value, Sequence)
        and not isinstance(value, (str, bytes))
        and bool(value)
        and (length is None or len(value) == length)
        and all(_finite(item) for item in value)
    )


def _virtual_force_ok(row: Mapping[str, object]) -> bool:
    axis = row.get("displacement_axis")
    step = row.get("displacement_step_m")
    energy_minus = row.get("energy_minus_j")
    energy_plus = row.get("energy_plus_j")
    force = row.get("force_n")
    axis_ok = _finite_vector(axis, 3) and math.isclose(
        sum(float(value) ** 2 for value in axis), 1.0, rel_tol=1.0e-12, abs_tol=1.0e-12
    )
    numeric_ok = all(_finite(value) for value in (step, energy_minus, energy_plus))
    expected_force = -((float(energy_plus) - float(energy_minus)) / (2.0 * float(step))) if numeric_ok and float(step) > 0.0 else None
    force_ok = (
        _finite_vector(force, 3)
        and expected_force is not None
        and all(
            math.isclose(float(component), expected_force * float(direction), rel_tol=1.0e-10, abs_tol=1.0e-10)
            for component, direction in zip(force, axis)
        )
    ) if axis_ok else False
    return (
        _generations(row, "mesh_generation", "energy_generation", "displacement_generation", "force_generation", "owner_generation", "result_generation")
        and _digest(row.get("minus_mesh_sha256"))
        and _digest(row.get("plus_mesh_sha256"))
        and row.get("minus_mesh_sha256") != row.get("plus_mesh_sha256")
        and row.get("result_minus_mesh_sha256") == row.get("minus_mesh_sha256")
        and row.get("result_plus_mesh_sha256") == row.get("plus_mesh_sha256")
        and axis_ok
        and row.get("result_displacement_axis") == axis
        and numeric_ok
        and float(step) > 0.0
        and row.get("result_displacement_step_m") == step
        and row.get("result_energy_minus_j") == energy_minus
        and row.get("result_energy_plus_j") == energy_plus
        and force_ok
        and row.get("result_force_n") == force
        and row.get("force_sign_convention") == "negative_energy_gradient"
        and row.get("result_force_sign_convention") == row.get("force_sign_convention")
        and str(row.get("solution_owner") or "").startswith("solution:")
        and row.get("result_solution_owner") == row.get("solution_owner")
        and _result(row)
    )


def _magnet_torque_ok(row: Mapping[str, object]) -> bool:
    wrapped = row.get("angles_wrapped_deg")
    unwrapped = row.get("angles_unwrapped_deg")
    energy = row.get("angular_energy_j")
    torque = row.get("torque_at_center_nm")
    vectors_ok = _finite_vector(wrapped) and _finite_vector(unwrapped) and _finite_vector(energy)
    samples_ok = (
        vectors_ok
        and len(wrapped) == len(unwrapped) == len(energy)
        and len(wrapped) >= 5
        and len(wrapped) % 2 == 1
        and all(float(left) < float(right) for left, right in zip(unwrapped, unwrapped[1:]))
        and all(math.isclose(float(wrapped_value) % 360.0, float(unwrapped_value) % 360.0, abs_tol=1.0e-12) for wrapped_value, unwrapped_value in zip(wrapped, unwrapped))
    )
    center = len(energy) // 2 if samples_ok else 0
    delta_angle = math.radians(float(unwrapped[center + 1]) - float(unwrapped[center - 1])) if samples_ok else 0.0
    expected_torque = -((float(energy[center + 1]) - float(energy[center - 1])) / delta_angle) if delta_angle > 0.0 else None
    return (
        _generations(row, "angle_generation", "unwrap_generation", "energy_generation", "derivative_generation", "owner_generation", "result_generation")
        and samples_ok
        and row.get("result_angles_wrapped_deg") == wrapped
        and row.get("result_angles_unwrapped_deg") == unwrapped
        and row.get("result_angular_energy_j") == energy
        and _finite(torque)
        and expected_torque is not None
        and math.isclose(float(torque), expected_torque, rel_tol=1.0e-10, abs_tol=1.0e-10)
        and row.get("result_torque_at_center_nm") == torque
        and row.get("derivative_sign_convention") == "negative_energy_gradient"
        and row.get("result_derivative_sign_convention") == row.get("derivative_sign_convention")
        and str(row.get("magnet_owner") or "").startswith("magnet:")
        and row.get("result_magnet_owner") == row.get("magnet_owner")
        and _result(row)
    )


def validate_public_identity(identity: object) -> dict[str, bool]:
    if not isinstance(identity, Mapping):
        return {}
    checks = validate_public_v53_identity(identity)
    virtual_force = identity.get(VIRTUAL_FORCE)
    magnet_torque = identity.get(MAGNET_TORQUE)
    if virtual_force is not None:
        checks["magnetic_force_v52_virtual_displacement_mesh_energy_sign_owner"] = (
            isinstance(virtual_force, Mapping) and _virtual_force_ok(virtual_force)
        )
    if magnet_torque is not None:
        checks["magnetic_force_v52_angular_energy_unwrap_derivative_owner"] = (
            isinstance(magnet_torque, Mapping) and _magnet_torque_ok(magnet_torque)
        )
    return checks
