from copy import deepcopy

from radia_mcp.radia_ngsolve.magnetostatic_energy_identity_v55 import BEARING, DEMAG, validate_public_identity


CASE_IDS = {"v55_public_demag_energy_hb_curvebranch_volume_material_owner_mismatch", "v55_public_magneticbearing_crossstiffness_matrix_symmetry_coordinate_owner_mismatch"}


def _payload():
    gen = lambda name, fields: {"generation": name, **{field: name for field in fields}}
    samples = [{"h_a_per_m": 0.0, "b_t": 1.2}, {"h_a_per_m": -1.0e5, "b_t": 0.8}, {"h_a_per_m": -3.0e5, "b_t": 0.4}]
    demag = {**gen("demag-v55", ("hb_generation", "branch_generation", "volume_generation", "energy_generation", "material_generation", "owner_generation", "result_generation")), "hb_samples": samples, "result_hb_samples": samples, "curve_branch": "descending_demag", "result_curve_branch": "descending_demag", "material_volume_m3": 1.0e-5, "result_material_volume_m3": 1.0e-5, "demag_energy_j": 1.0, "result_demag_energy_j": 1.0, "material_revision": "magnet-v55-r4", "result_material_revision": "magnet-v55-r4", "material_owner": "material:magnet-v55", "result_material_owner": "material:magnet-v55", "solution_owner": "solution:demag-v55", "result_solution_owner": "solution:demag-v55", "result_sha256": "9" * 64, "accepted_result_sha256": "9" * 64}
    matrix = [[1200.0, 50.0], [50.0, 900.0]]; basis = {"x": [1.0, 0.0, 0.0], "y": [0.0, 1.0, 0.0]}; load = [0.0, 0.0, 0.025]
    bearing = {**gen("bearing-v55", ("matrix_generation", "basis_generation", "reciprocity_generation", "loadpoint_generation", "owner_generation", "result_generation")), "cross_stiffness_n_per_m": matrix, "result_cross_stiffness_n_per_m": matrix, "coordinate_basis": basis, "result_coordinate_basis": basis, "reciprocity_tolerance": 1.0e-10, "result_reciprocity_tolerance": 1.0e-10, "load_point_m": load, "result_load_point_m": load, "body_owner": "body:bearing-v55", "result_body_owner": "body:bearing-v55", "result_sha256": "a" * 64, "accepted_result_sha256": "a" * 64}
    return {DEMAG: demag, BEARING: bearing}


def test_v55_positive_identities_are_accepted():
    assert all(validate_public_identity(_payload()).values())


def test_v55_frozen_mutations_are_rejected():
    payload = deepcopy(_payload()); payload[DEMAG]["result_solution_owner"] = "solution:stale"; payload[BEARING]["result_body_owner"] = "body:stale"
    assert not all(validate_public_identity(payload).values())


def test_v55_self_consistent_nonphysical_records_are_rejected():
    payload = deepcopy(_payload()); payload[DEMAG]["demag_energy_j"] = payload[DEMAG]["result_demag_energy_j"] = -1.0
    payload[BEARING]["cross_stiffness_n_per_m"] = payload[BEARING]["result_cross_stiffness_n_per_m"] = [[1200.0, 500.0], [-50.0, 900.0]]
    assert not all(validate_public_identity(payload).values())


def test_v55_malformed_values_reject_without_raising():
    payload = deepcopy(_payload()); payload[DEMAG]["hb_samples"] = [{"h_a_per_m": [0.0], "b_t": 1.2}]; payload[BEARING]["coordinate_basis"] = {"x": [[1.0], 0.0, 0.0]}
    assert not all(validate_public_identity(payload).values())
