from __future__ import annotations

from test_magnetic_force_generalization_v31 import _gate
from test_magnetic_force_generalization_v33 import _identity_v33


_PROMOTED_CASE_IDS = (
    "v34_public_kelvin_transform_open_boundary_radius_energy_flux_far_field_mismatch",
    "v34_public_force_contour_path_airgap_stress_fourier_virtual_work_symmetry_mismatch",
)


def _identity_v34():
    identity = _identity_v33()
    generation = "kelvin-open-boundary-211"
    identity[
        "kelvin_transform_radius_permeability_jacobian_interface_energy_flux_far_field_mesh_owner_result_identity"
    ] = {
        "kelvin_generation": generation,
        **{key: generation for key in (
            "radius_generation", "mapping_generation", "jacobian_generation",
            "interface_generation", "energy_generation", "flux_generation",
            "far_field_generation", "mesh_generation", "owner_generation", "result_generation")},
        "kelvin_radius_m": 1.0, "result_kelvin_radius_m": 1.0,
        "mapped_permeability_relative": [1.0, 4.0, 16.0],
        "result_mapped_permeability_relative": [1.0, 4.0, 16.0],
        "mapping_jacobian_determinants": [1.0, 0.25, 0.0625],
        "result_mapping_jacobian_determinants": [1.0, 0.25, 0.0625],
        "interface_potential_jump": 0.0, "result_interface_potential_jump": 0.0,
        "interface_normal_flux_jump_wb": 0.0, "result_interface_normal_flux_jump_wb": 0.0,
        "magnetic_energy_j": 0.012, "result_magnetic_energy_j": 0.012,
        "inner_flux_wb": 0.02, "result_inner_flux_wb": 0.02,
        "outer_flux_wb": -0.02, "result_outer_flux_wb": -0.02,
        "far_field_radius_m": [2.0, 4.0, 8.0], "result_far_field_radius_m": [2.0, 4.0, 8.0],
        "far_field_potential": [0.5, 0.25, 0.125], "result_far_field_potential": [0.5, 0.25, 0.125],
        "kelvin_mesh_sha256": "1" * 64, "result_kelvin_mesh_sha256": "1" * 64,
        "kelvin_result_owner": "open-boundary/kelvin-211",
        "accepted_kelvin_result_owner": "open-boundary/kelvin-211",
        "kelvin_result_sha256": "2" * 64, "accepted_kelvin_result_sha256": "2" * 64,
    }
    generation = "airgap-force-211"
    identity[
        "force_airgap_contour_stress_fourier_virtual_work_symmetry_torque_origin_field_result_identity"
    ] = {
        "force_generation": generation,
        **{key: generation for key in (
            "contour_generation", "stress_generation", "fourier_generation",
            "virtual_work_generation", "symmetry_generation", "torque_generation",
            "field_generation", "owner_generation", "result_generation")},
        "airgap_contour_ids": [1, 2], "result_airgap_contour_ids": [1, 2],
        "contour_forces_n": [[10.0, 0.0], [10.00000001, 0.0]],
        "result_contour_forces_n": [[10.0, 0.0], [10.00000001, 0.0]],
        "fourier_stress_harmonics_n": [[0, 10.000000005], [1, 0.0], [2, 0.0]],
        "result_fourier_stress_harmonics_n": [[0, 10.000000005], [1, 0.0], [2, 0.0]],
        "virtual_work_force_n": 10.000000005, "result_virtual_work_force_n": 10.000000005,
        "virtual_work_displacement_m": 1.0e-5, "result_virtual_work_displacement_m": 1.0e-5,
        "torque_origin_m": [0.0, 0.0], "result_torque_origin_m": [0.0, 0.0],
        "torque_nm": 0.0, "result_torque_nm": 0.0,
        "force_symmetry": "mirror_y", "result_force_symmetry": "mirror_y",
        "force_field_sha256": "3" * 64, "result_force_field_sha256": "3" * 64,
        "force_result_owner": "force/airgap-211", "accepted_force_result_owner": "force/airgap-211",
        "force_result_sha256": "4" * 64, "accepted_force_result_sha256": "4" * 64,
    }
    return identity


def test_v34_public_positive_kelvin_and_force_closure():
    assert _gate(_identity_v34())["status"] == "ok"


def test_v34_public_kelvin_transform_open_boundary_radius_energy_flux_far_field_mismatch():
    identity = _identity_v34()
    identity["kelvin_transform_radius_permeability_jacobian_interface_energy_flux_far_field_mesh_owner_result_identity"].update({
        "radius_generation": "kelvin-open-boundary-210", "flux_generation": "kelvin-open-boundary-209", "result_generation": "kelvin-open-boundary-208",
        "result_kelvin_radius_m": 2.0, "result_mapped_permeability_relative": [1.0, 1.0, 1.0],
        "result_mapping_jacobian_determinants": [1.0, -0.25, 0.0],
        "result_interface_potential_jump": 0.1, "result_interface_normal_flux_jump_wb": 0.02,
        "result_magnetic_energy_j": -0.012, "result_outer_flux_wb": 0.01,
        "result_far_field_potential": [0.5, 0.6, 0.7], "result_kelvin_mesh_sha256": "9" * 64,
        "accepted_kelvin_result_owner": "stale/owner", "accepted_kelvin_result_sha256": "a" * 64})
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "kelvin_open_boundary_closes_mapping_interface_energy_flux_far_field_mesh_owner_and_result"
    ]


def test_v34_public_force_contour_path_airgap_stress_fourier_virtual_work_symmetry_mismatch():
    identity = _identity_v34()
    identity["force_airgap_contour_stress_fourier_virtual_work_symmetry_torque_origin_field_result_identity"].update({
        "contour_generation": "airgap-force-210", "virtual_work_generation": "airgap-force-209", "result_generation": "airgap-force-208",
        "result_airgap_contour_ids": [2, 1], "result_contour_forces_n": [[10.0, 0.0], [-5.0, 2.0]],
        "result_fourier_stress_harmonics_n": [[0, 10.0], [1, 4.0], [2, -3.0]],
        "result_virtual_work_force_n": -10.0, "result_virtual_work_displacement_m": -1.0e-5,
        "result_torque_origin_m": [1.0, 0.0], "result_torque_nm": 5.0,
        "result_force_symmetry": "none", "result_force_field_sha256": "b" * 64,
        "accepted_force_result_owner": "stale/force", "accepted_force_result_sha256": "c" * 64})
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "airgap_force_closes_contours_fourier_virtual_work_symmetry_torque_field_owner_and_result"
    ]


def test_v34_rejects_self_consistent_nondecaying_kelvin_far_field():
    identity = _identity_v34()
    row = identity["kelvin_transform_radius_permeability_jacobian_interface_energy_flux_far_field_mesh_owner_result_identity"]
    row["far_field_potential"] = [0.5, 0.5, 0.5]
    row["result_far_field_potential"] = [0.5, 0.5, 0.5]
    assert _gate(identity)["status"] == "needs_attention"


def test_v34_rejects_self_consistent_force_contour_disagreement():
    identity = _identity_v34()
    row = identity["force_airgap_contour_stress_fourier_virtual_work_symmetry_torque_origin_field_result_identity"]
    row["contour_forces_n"] = [[10.0, 0.0], [5.0, 0.0]]
    row["result_contour_forces_n"] = [[10.0, 0.0], [5.0, 0.0]]
    assert _gate(identity)["status"] == "needs_attention"
