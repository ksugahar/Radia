"""Solver-neutral dual-lane motor and magnetic-force lineage checks."""

from __future__ import annotations

import math
from collections.abc import Mapping


MOTOR = "v47_public_motor_dual_lane_geometry_material_excitation_operating_point_identity_mismatch"
FORCE = "v47_public_force_coenergy_displacement_pair_body_owner_aggregation_mismatch"
MOTOR_LANES = ["ngsolve_age", "hdiv_mmm_hcurl_eddy_bubble"]


def _sha(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _generation(row: Mapping[str, object], names: tuple[str, ...]) -> bool:
    value = str(row.get("generation") or "")
    return bool(value) and all(row.get(name) == value for name in names)


def _digest(row: Mapping[str, object]) -> bool:
    return _sha(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _finite_mapping(value: object) -> bool:
    return (
        isinstance(value, Mapping)
        and bool(value)
        and all(isinstance(item, (int, float)) and math.isfinite(float(item)) for item in value.values())
    )


def _excitation(value: object) -> bool:
    if not isinstance(value, Mapping):
        return False
    phases = value.get("phase_order")
    currents = value.get("current_a")
    return (
        phases == ["A", "B", "C"]
        and isinstance(currents, list)
        and len(currents) == 3
        and all(isinstance(item, (int, float)) and math.isfinite(float(item)) for item in currents)
    )


def _motor_ok(row: Mapping[str, object]) -> bool:
    excitation = row.get("excitation_identity")
    geometry = row.get("geometry_identity_sha256")
    material = row.get("material_identity")
    operating_point = row.get("operating_point_key")
    return (
        _generation(
            row,
            (
                "geometry_generation",
                "material_generation",
                "excitation_generation",
                "operating_point_generation",
                "lane_a_generation",
                "lane_b_generation",
                "result_generation",
            ),
        )
        and row.get("lane_ids") == row.get("result_lane_ids") == MOTOR_LANES
        and _sha(geometry)
        and row.get("lane_a_geometry_identity_sha256") == geometry
        and row.get("lane_b_geometry_identity_sha256") == geometry
        and str(material or "").startswith("material:")
        and row.get("lane_a_material_identity") == material
        and row.get("lane_b_material_identity") == material
        and _excitation(excitation)
        and row.get("lane_a_excitation_identity") == excitation
        and row.get("lane_b_excitation_identity") == excitation
        and isinstance(operating_point, str)
        and bool(operating_point)
        and row.get("lane_a_operating_point_key") == operating_point
        and row.get("lane_b_operating_point_key") == operating_point
        and _digest(row)
    )


def _force_ok(row: Mapping[str, object]) -> bool:
    displacement = row.get("displacement_pair_m")
    coenergy = row.get("coenergy_pair_j")
    components = row.get("component_force_n")
    aggregate = row.get("aggregated_force_n")
    return (
        _generation(
            row,
            (
                "coenergy_generation",
                "displacement_generation",
                "body_owner_generation",
                "aggregation_generation",
                "result_generation",
            ),
        )
        and isinstance(displacement, list)
        and len(displacement) == 2
        and all(isinstance(value, (int, float)) and math.isfinite(float(value)) for value in displacement)
        and float(displacement[0]) < float(displacement[1])
        and row.get("result_displacement_pair_m") == displacement
        and isinstance(coenergy, list)
        and len(coenergy) == 2
        and all(isinstance(value, (int, float)) and math.isfinite(float(value)) for value in coenergy)
        and row.get("result_coenergy_pair_j") == coenergy
        and str(row.get("body_owner") or "").startswith("body:")
        and row.get("result_body_owner") == row.get("body_owner")
        and _finite_mapping(components)
        and row.get("result_component_force_n") == components
        and isinstance(aggregate, (int, float))
        and math.isfinite(float(aggregate))
        and math.isclose(sum(float(value) for value in components.values()), float(aggregate), rel_tol=1e-12, abs_tol=1e-12)
        and row.get("result_aggregated_force_n") == aggregate
        and _digest(row)
    )


def validate_public_v47_identity(identity: object) -> dict[str, bool]:
    if not isinstance(identity, Mapping):
        return {}
    checks: dict[str, bool] = {}
    motor = identity.get(MOTOR)
    force = identity.get(FORCE)
    if motor is not None:
        checks["motor_v47_dual_lane_shared_physics_identity"] = isinstance(motor, Mapping) and _motor_ok(motor)
    if force is not None:
        checks["magnetic_force_v47_coenergy_body_aggregation_identity"] = isinstance(force, Mapping) and _force_ok(force)
    return checks
