"""Solver-neutral demagnetization and virtual-work artifact checks for v49."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


DEMAG = "nonlinear_demag_recoil_branch_temperature_loadstep_result_owner_identity"
VIRTUAL_WORK = "virtual_work_displacement_frame_mesh_state_energy_owner_identity"


def _digest(value: object) -> bool:
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generations(row: Mapping[str, object], *fields: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(field) == generation for field in fields)


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _finite_vector(value: object, length: int | None = None) -> bool:
    return (
        isinstance(value, Sequence)
        and not isinstance(value, (str, bytes))
        and (length is None or len(value) == length)
        and bool(value)
        and all(isinstance(item, (int, float)) and math.isfinite(float(item)) for item in value)
    )


def _recoil_branch(value: object) -> bool:
    if not isinstance(value, list) or len(value) < 3 or not all(_finite_vector(row, 2) for row in value):
        return False
    flux_density = [float(row[0]) for row in value]
    field = [float(row[1]) for row in value]
    return all(right > left for left, right in zip(flux_density, flux_density[1:])) and all(
        right > left for left, right in zip(field, field[1:])
    )


def _demag_ok(row: Mapping[str, object]) -> bool:
    branch = row.get("recoil_branch_t_a_per_m")
    temperature = row.get("temperature_c")
    load_step = row.get("load_step")
    owner = str(row.get("magnet_owner") or "")
    return (
        _generations(row, "branch_generation", "temperature_generation", "loadstep_generation", "result_generation")
        and _recoil_branch(branch)
        and row.get("result_recoil_branch_t_a_per_m") == branch
        and isinstance(temperature, (int, float))
        and math.isfinite(float(temperature))
        and -273.15 < float(temperature) <= 500.0
        and row.get("result_temperature_c") == temperature
        and isinstance(load_step, int)
        and not isinstance(load_step, bool)
        and load_step >= 0
        and row.get("result_load_step") == load_step
        and owner.startswith("magnet:")
        and row.get("result_magnet_owner") == owner
        and _result(row)
    )


def _virtual_work_ok(row: Mapping[str, object]) -> bool:
    displacement = row.get("displacement_m")
    frame = str(row.get("displacement_frame") or "")
    energy = row.get("energy_j")
    force_owner = str(row.get("force_owner") or "")
    symmetric_displacement = (
        _finite_vector(displacement, 3)
        and float(displacement[0]) < float(displacement[1]) < float(displacement[2])
        and math.isclose(float(displacement[1]), 0.0, rel_tol=0.0, abs_tol=1.0e-15)
        and math.isclose(float(displacement[0]), -float(displacement[2]), rel_tol=1.0e-12, abs_tol=1.0e-15)
    )
    return (
        _generations(row, "displacement_generation", "frame_generation", "mesh_generation", "energy_generation", "result_generation")
        and symmetric_displacement
        and row.get("result_displacement_m") == displacement
        and frame.startswith("frame:")
        and row.get("result_displacement_frame") == frame
        and _digest(row.get("mesh_state_sha256"))
        and row.get("result_mesh_state_sha256") == row.get("mesh_state_sha256")
        and _finite_vector(energy, 3)
        and row.get("result_energy_j") == energy
        and force_owner.startswith("force:")
        and row.get("result_force_owner") == force_owner
        and _result(row)
    )


def validate_public_identity(identity: object) -> dict[str, bool]:
    """Validate optional demagnetization and virtual-work identity records."""
    if not isinstance(identity, Mapping):
        return {}
    checks: dict[str, bool] = {}
    demag = identity.get(DEMAG)
    virtual_work = identity.get(VIRTUAL_WORK)
    if demag is not None:
        checks["demag_v49_recoil_temperature_loadstep_digest_owner"] = isinstance(demag, Mapping) and _demag_ok(demag)
    if virtual_work is not None:
        checks["force_v49_virtual_work_displacement_frame_mesh_energy_owner"] = isinstance(virtual_work, Mapping) and _virtual_work_ok(virtual_work)
    return checks
