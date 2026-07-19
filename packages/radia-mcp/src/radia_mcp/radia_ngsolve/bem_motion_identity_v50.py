"""BEM singular-integration and motional-EMF identity checks for v50."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


BEM = "bem_singular_quadrature_self_panel_nearfield_regularization_mesh_owner_identity"
MOTION = "motion_emf_velocity_frame_conductor_path_direction_result_owner_identity"


def _digest(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generations(row: Mapping[str, object], *fields: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(field) == generation for field in fields)


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _finite_vector(value: object, length: int) -> bool:
    return (
        isinstance(value, Sequence)
        and not isinstance(value, (str, bytes))
        and len(value) == length
        and all(isinstance(item, (int, float)) and math.isfinite(float(item)) for item in value)
    )


def _nearfield(value: object) -> bool:
    return (
        isinstance(value, Mapping)
        and set(value) == {"distance_ratio", "regularization", "max_depth"}
        and isinstance(value.get("distance_ratio"), (int, float))
        and math.isfinite(float(value["distance_ratio"]))
        and 0.0 < float(value["distance_ratio"]) <= 1.0
        and value.get("regularization") in {"adaptive-subdivision", "singularity-subtraction"}
        and isinstance(value.get("max_depth"), int)
        and 1 <= int(value["max_depth"]) <= 32
    )


def _bem_ok(row: Mapping[str, object]) -> bool:
    quadrature = str(row.get("singular_quadrature") or "")
    self_panel = str(row.get("self_panel_treatment") or "")
    nearfield = row.get("nearfield_regularization")
    owner = str(row.get("mesh_owner") or "")
    return (
        _generations(row, "quadrature_generation", "self_panel_generation", "nearfield_generation", "mesh_generation", "result_generation")
        and quadrature.startswith("duffy-order-")
        and row.get("result_singular_quadrature") == quadrature
        and self_panel in {"analytic-solid-angle", "principal-value"}
        and row.get("result_self_panel_treatment") == self_panel
        and _nearfield(nearfield)
        and row.get("result_nearfield_regularization") == nearfield
        and owner.startswith("mesh:")
        and row.get("result_mesh_owner") == owner
        and _result(row)
    )


def _path(value: object) -> bool:
    if not isinstance(value, list) or len(value) < 2 or not all(_finite_vector(point, 3) for point in value):
        return False
    return all(any(not math.isclose(float(a), float(b), abs_tol=0.0) for a, b in zip(left, right)) for left, right in zip(value, value[1:]))


def _motion_ok(row: Mapping[str, object]) -> bool:
    velocity = row.get("velocity_m_s")
    frame = str(row.get("velocity_frame") or "")
    path = row.get("conductor_path_m")
    direction = row.get("integration_direction")
    owner = str(row.get("emf_result_owner") or "")
    return (
        _generations(row, "velocity_generation", "frame_generation", "path_generation", "direction_generation", "result_generation")
        and _finite_vector(velocity, 3)
        and any(not math.isclose(float(item), 0.0, abs_tol=0.0) for item in velocity)
        and row.get("result_velocity_m_s") == velocity
        and frame.startswith("frame:")
        and row.get("result_velocity_frame") == frame
        and _path(path)
        and row.get("result_conductor_path_m") == path
        and direction in {"path-forward", "path-reverse"}
        and row.get("result_integration_direction") == direction
        and owner.startswith("emf:")
        and row.get("result_emf_owner") == owner
        and _result(row)
    )


def validate_public_identity(identity: object) -> dict[str, bool]:
    """Validate optional v50 BEM and motional-EMF identity records."""
    if not isinstance(identity, Mapping):
        return {}
    checks: dict[str, bool] = {}
    bem = identity.get(BEM)
    motion = identity.get(MOTION)
    if bem is not None:
        checks["magnetic_force_v50_bem_singular_self_nearfield_mesh_owner"] = isinstance(bem, Mapping) and _bem_ok(bem)
    if motion is not None:
        checks["magnetic_force_v50_motion_emf_velocity_frame_path_direction_owner"] = isinstance(motion, Mapping) and _motion_ok(motion)
    return checks
