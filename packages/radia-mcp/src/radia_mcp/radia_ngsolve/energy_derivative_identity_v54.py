"""Magnetic-charge and directional-stiffness identity checks for v54."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


CHARGE = "magneticcharge_neutrality_surface_normal_material_region_owner_identity"
STIFFNESS = "maglev_stiffness_forcegradient_coordinate_loadpoint_owner_identity"


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


def _unit_vector(value: object, length: int = 3) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes)) and len(value) == length and all(_finite(item) for item in value) and math.isclose(sum(float(item) ** 2 for item in value), 1.0, rel_tol=0.0, abs_tol=1.0e-12)


def _charge_ok(row: Mapping[str, object]) -> bool:
    charges = row.get("surface_charges")
    charges_ok = isinstance(charges, Sequence) and not isinstance(charges, (str, bytes)) and len(charges) >= 2
    panel_charges: dict[int, float] = {}
    if charges_ok:
        for item in charges:
            if not isinstance(item, Mapping) or set(item) != {"panel", "magnetic_charge_a_m"}:
                charges_ok = False
                break
            panel = item["panel"]
            charge = item["magnetic_charge_a_m"]
            if not (isinstance(panel, int) and not isinstance(panel, bool) and panel > 0 and panel not in panel_charges and _finite(charge)):
                charges_ok = False
                break
            panel_charges[panel] = float(charge)
        charges_ok = charges_ok and math.isclose(sum(panel_charges.values()), 0.0, rel_tol=0.0, abs_tol=1.0e-12)
    panel_keys = {str(panel) for panel in panel_charges}
    normals = row.get("surface_normals")
    regions = row.get("material_region_map")
    orientations = row.get("boundary_orientation")
    normals_ok = isinstance(normals, Mapping) and set(normals) == panel_keys and all(_unit_vector(normal) for normal in normals.values())
    regions_ok = isinstance(regions, Mapping) and set(regions) == panel_keys and all(isinstance(region, str) and region.startswith("region:") for region in regions.values())
    orientations_ok = isinstance(orientations, Mapping) and set(orientations) == panel_keys
    if orientations_ok:
        orientations_ok = all(
            isinstance(orientation, int)
            and not isinstance(orientation, bool)
            and orientation in {-1, 1}
            and panel_charges[int(panel)] * orientation >= 0.0
            for panel, orientation in orientations.items()
        )
    return (
        _generations(row, "charge_generation", "normal_generation", "region_generation", "orientation_generation", "owner_generation", "result_generation")
        and charges_ok and normals_ok and regions_ok and orientations_ok
        and row.get("result_surface_charges") == charges
        and row.get("result_surface_normals") == normals
        and row.get("result_material_region_map") == regions
        and row.get("result_boundary_orientation") == orientations
        and str(row.get("solution_owner") or "").startswith("solution:")
        and row.get("result_solution_owner") == row.get("solution_owner")
        and _result(row)
    )


def _stiffness_ok(row: Mapping[str, object]) -> bool:
    gradient = row.get("force_gradient_n_per_m")
    stiffness = row.get("stiffness_n_per_m")
    direction = row.get("coordinate_direction")
    load_point = row.get("load_point_m")
    increment = row.get("displacement_increment_m")
    return (
        _generations(row, "gradient_generation", "coordinate_generation", "loadpoint_generation", "increment_generation", "owner_generation", "result_generation")
        and _finite(gradient) and _finite(stiffness)
        and float(stiffness) > 0.0
        and math.isclose(float(stiffness), -float(gradient), rel_tol=1.0e-12, abs_tol=1.0e-12)
        and row.get("result_force_gradient_n_per_m") == gradient
        and row.get("result_stiffness_n_per_m") == stiffness
        and _unit_vector(direction)
        and row.get("result_coordinate_direction") == direction
        and isinstance(load_point, Sequence) and not isinstance(load_point, (str, bytes)) and len(load_point) == 3 and all(_finite(value) for value in load_point)
        and row.get("result_load_point_m") == load_point
        and _finite(increment) and float(increment) > 0.0
        and row.get("result_displacement_increment_m") == increment
        and str(row.get("body_owner") or "").startswith("body:")
        and row.get("result_body_owner") == row.get("body_owner")
        and _result(row)
    )


def validate_public_identity(identity: object) -> dict[str, bool]:
    if not isinstance(identity, Mapping):
        return {}
    checks: dict[str, bool] = {}
    charge = identity.get(CHARGE)
    stiffness = identity.get(STIFFNESS)
    if charge is not None:
        checks["magnetic_force_v54_charge_neutrality_normal_region_orientation_owner"] = isinstance(charge, Mapping) and _charge_ok(charge)
    if stiffness is not None:
        checks["magnetic_force_v54_stiffness_gradient_coordinate_loadpoint_owner"] = isinstance(stiffness, Mapping) and _stiffness_ok(stiffness)
    return checks
