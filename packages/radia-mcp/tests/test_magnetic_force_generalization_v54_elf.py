from copy import deepcopy

from radia_mcp.radia_ngsolve.energy_derivative_identity_v54 import CHARGE, STIFFNESS, validate_public_identity


CASE_IDS = {
    "v54_public_magneticcharge_neutrality_surface_normal_material_region_owner_mismatch",
    "v54_public_maglev_stiffness_forcegradient_coordinate_loadpoint_owner_mismatch",
}


def _generations(generation: str, fields: tuple[str, ...]) -> dict[str, str]:
    return {"generation": generation, **{field: generation for field in fields}}


def _payload():
    charges = [{"panel": 11, "magnetic_charge_a_m": 0.25}, {"panel": 12, "magnetic_charge_a_m": -0.25}]; normals = {"11": [0.0, 0.0, 1.0], "12": [0.0, 0.0, -1.0]}; regions = {"11": "region:magnet", "12": "region:air"}; orientations = {"11": 1, "12": -1}
    charge = {**_generations("charge-v54", ("charge_generation", "normal_generation", "region_generation", "orientation_generation", "owner_generation", "result_generation")), "surface_charges": charges, "result_surface_charges": charges, "surface_normals": normals, "result_surface_normals": normals, "material_region_map": regions, "result_material_region_map": regions, "boundary_orientation": orientations, "result_boundary_orientation": orientations, "solution_owner": "solution:v54", "result_solution_owner": "solution:v54", "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64}
    direction = [0.0, 0.0, 1.0]; load_point = [0.0, 0.0, 0.025]
    stiffness = {**_generations("stiffness-v54", ("gradient_generation", "coordinate_generation", "loadpoint_generation", "increment_generation", "owner_generation", "result_generation")), "force_gradient_n_per_m": -1200.0, "result_force_gradient_n_per_m": -1200.0, "stiffness_n_per_m": 1200.0, "result_stiffness_n_per_m": 1200.0, "coordinate_direction": direction, "result_coordinate_direction": direction, "load_point_m": load_point, "result_load_point_m": load_point, "displacement_increment_m": 1.0e-5, "result_displacement_increment_m": 1.0e-5, "body_owner": "body:v54", "result_body_owner": "body:v54", "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64}
    return {CHARGE: charge, STIFFNESS: stiffness}


def test_v54_positive_identities_are_accepted():
    assert all(validate_public_identity(_payload()).values())


def test_v54_frozen_mutations_are_rejected():
    payload = deepcopy(_payload())
    payload[CHARGE]["result_surface_charges"] = [{"panel": 11, "magnetic_charge_a_m": 0.25}]
    payload[STIFFNESS]["result_coordinate_direction"] = [1.0, 0.0, 0.0]
    assert not all(validate_public_identity(payload).values())


def test_v54_self_consistent_nonphysical_records_are_rejected():
    payload = deepcopy(_payload())
    bad_charges = [{"panel": 11, "magnetic_charge_a_m": 0.25}, {"panel": 12, "magnetic_charge_a_m": 0.1}]
    payload[CHARGE]["surface_charges"] = payload[CHARGE]["result_surface_charges"] = bad_charges
    payload[STIFFNESS]["force_gradient_n_per_m"] = payload[STIFFNESS]["result_force_gradient_n_per_m"] = 1200.0
    assert not all(validate_public_identity(payload).values())


def test_v54_malformed_values_reject_without_raising():
    payload = deepcopy(_payload())
    payload[CHARGE]["boundary_orientation"] = {"11": [1], "12": -1}
    payload[STIFFNESS]["coordinate_direction"] = [[0.0], 0.0, 1.0]
    assert not all(validate_public_identity(payload).values())
