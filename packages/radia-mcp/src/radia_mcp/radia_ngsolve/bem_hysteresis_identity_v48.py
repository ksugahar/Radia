"""Solver-neutral BEM discretization and hysteresis history identity checks."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


BEM = "bem_near_far_panel_quadrature_normal_solid_angle_mesh_revision_identity"
HYSTERESIS = "hysteresis_minor_loop_return_point_state_temperature_frequency_owner_identity"


def _digest(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _generation(row: Mapping[str, object], fields: tuple[str, ...]) -> bool:
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


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _bem_ok(row: Mapping[str, object]) -> bool:
    panels = row.get("panel_ids")
    near = row.get("near_quadrature_order")
    far = row.get("far_quadrature_order")
    normals = row.get("panel_normals")
    angles = row.get("panel_solid_angles_sr")
    count = len(panels) if isinstance(panels, list) else 0
    return (
        _generation(row, ("quadrature_generation", "normal_generation", "solid_angle_generation", "mesh_generation", "result_generation"))
        and count > 0
        and len(set(panels)) == count
        and isinstance(near, list)
        and isinstance(far, list)
        and len(near) == len(far) == count
        and all(isinstance(n, int) and isinstance(f, int) and n >= f > 0 for n, f in zip(near, far, strict=True))
        and isinstance(normals, list)
        and len(normals) == count
        and all(_finite_vector(normal, 3) and math.isclose(sum(float(value) ** 2 for value in normal), 1.0, rel_tol=1.0e-12, abs_tol=1.0e-12) for normal in normals)
        and _finite_vector(angles, count)
        and all(float(angle) > 0.0 for angle in angles)
        and math.isclose(sum(float(angle) for angle in angles), 4.0 * math.pi, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and row.get("result_panel_ids") == panels
        and row.get("result_near_quadrature_order") == near
        and row.get("result_far_quadrature_order") == far
        and row.get("result_panel_normals") == normals
        and row.get("result_panel_solid_angles_sr") == angles
        and str(row.get("mesh_revision") or "").startswith("mesh:")
        and row.get("result_mesh_revision") == row.get("mesh_revision")
        and _result(row)
    )


def _hysteresis_ok(row: Mapping[str, object]) -> bool:
    points = row.get("return_points")
    state = row.get("internal_state")
    temperature = row.get("temperature_k")
    frequency = row.get("frequency_hz")
    return (
        _generation(row, ("return_point_generation", "state_generation", "environment_generation", "material_generation", "result_generation"))
        and isinstance(points, list)
        and len(points) >= 2
        and all(_finite_vector(point, 2) for point in points)
        and row.get("result_return_points") == points
        and isinstance(state, Mapping)
        and state.get("branch") in {"ascending", "descending"}
        and isinstance(state.get("memory_depth"), int)
        and int(state["memory_depth"]) >= 0
        and row.get("result_internal_state") == state
        and isinstance(temperature, (int, float))
        and math.isfinite(float(temperature))
        and float(temperature) > 0.0
        and row.get("result_temperature_k") == temperature
        and isinstance(frequency, (int, float))
        and math.isfinite(float(frequency))
        and float(frequency) >= 0.0
        and row.get("result_frequency_hz") == frequency
        and str(row.get("material_owner") or "").startswith("material:")
        and row.get("result_material_owner") == row.get("material_owner")
        and _result(row)
    )


def validate_public_identity(identity: object) -> dict[str, bool]:
    if not isinstance(identity, Mapping):
        return {}
    checks: dict[str, bool] = {}
    bem = identity.get(BEM)
    hysteresis = identity.get(HYSTERESIS)
    if bem is not None:
        checks["bem_v48_panel_quadrature_normal_solid_angle_mesh"] = isinstance(bem, Mapping) and _bem_ok(bem)
    if hysteresis is not None:
        checks["hysteresis_v48_return_state_environment_material_owner"] = isinstance(hysteresis, Mapping) and _hysteresis_ok(hysteresis)
    return checks
