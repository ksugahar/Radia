from __future__ import annotations

import math

from test_femm_generalization_v19 import _gate
from test_femm_generalization_v25 import _identity_v25
from test_force_coenergy_gate import _quadratic_case


def _identity_v26(sample_count):
    identity = _identity_v25(sample_count)
    identity["weighted_stress_tensor_mask_mesh_region_force_torque_energy_generation_identity"] = {
        "force_generation": "force-131", "mask_force_generation": "force-131",
        "mesh_force_generation": "force-131", "region_force_generation": "force-131",
        "energy_force_generation": "force-131", "torque_force_generation": "force-131",
        "result_force_generation": "force-131", "body_group_ids": [3, 4],
        "result_body_group_ids": [3, 4], "weighted_mask_node_ids": [101, 102, 103, 104],
        "result_weighted_mask_node_ids": [101, 102, 103, 104],
        "integration_region": "surrounding-air", "result_integration_region": "surrounding-air",
        "mesh_sha256": "1" * 64, "result_mesh_sha256": "1" * 64,
        "mask_sha256": "2" * 64, "result_mask_sha256": "2" * 64,
        "integration_region_sha256": "3" * 64, "result_integration_region_sha256": "3" * 64,
        "weighted_force_n": [12.5, -0.2, 0.0], "result_weighted_force_n": [12.5, -0.2, 0.0],
        "weighted_torque_nm": [0.0, 0.0, 1.8], "result_weighted_torque_nm": [0.0, 0.0, 1.8],
        "magnetic_energy_j": 0.44, "result_magnetic_energy_j": 0.44,
        "magnetic_coenergy_j": 0.61, "result_magnetic_coenergy_j": 0.61,
        "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64,
    }
    total_force = 1.2
    identity["axisymmetric_planar_depth_two_pi_r_force_normalization_coordinate_unit_generation_identity"] = {
        "normalization_generation": "normalization-131",
        "planar_depth_normalization_generation": "normalization-131",
        "radius_normalization_generation": "normalization-131",
        "coordinate_normalization_generation": "normalization-131",
        "unit_normalization_generation": "normalization-131",
        "mesh_normalization_generation": "normalization-131",
        "result_normalization_generation": "normalization-131",
        "planar_depth_m": 0.05, "result_planar_depth_m": 0.05,
        "planar_force_n_per_m": 24.0, "result_planar_force_n_per_m": 24.0,
        "planar_total_force_n": total_force, "result_planar_total_force_n": total_force,
        "axisymmetric_radius_m": 0.03, "result_axisymmetric_radius_m": 0.03,
        "axisymmetric_meridian_force_n_per_rad": total_force / (2.0 * math.pi),
        "result_axisymmetric_meridian_force_n_per_rad": total_force / (2.0 * math.pi),
        "axisymmetric_total_force_n": total_force, "result_axisymmetric_total_force_n": total_force,
        "radius_measure_convention": "2*pi*r", "result_radius_measure_convention": "2*pi*r",
        "coordinate_convention": "r_z_right_handed", "result_coordinate_convention": "r_z_right_handed",
        "force_unit": "N_total_3d", "result_force_unit": "N_total_3d",
        "mesh_sha256": "5" * 64, "result_mesh_sha256": "5" * 64,
    }
    return identity


def test_v26_public_positive_force_closure_and_normalization_identity():
    positions, _, _ = _quadratic_case()
    assert _gate(_identity_v26(len(positions)))["status"] == "ok"


def test_v26_public_weighted_stress_tensor_mask_mesh_region_force_torque_energy_generation_mismatch():
    positions, _, _ = _quadratic_case()
    identity = _identity_v26(len(positions))
    identity["weighted_stress_tensor_mask_mesh_region_force_torque_energy_generation_identity"].update({
        "mask_force_generation": "force-130", "mesh_force_generation": "force-129",
        "region_force_generation": "force-128", "result_body_group_ids": [4],
        "result_weighted_mask_node_ids": [101, 104], "result_integration_region": "rotor-steel",
        "result_mesh_sha256": "a" * 64, "result_weighted_force_n": [9.0, 0.4, 0.0],
        "result_weighted_torque_nm": [0.0, 0.0, -0.8], "result_magnetic_energy_j": 0.31,
    })
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["weighted_stress_tensor_uses_current_mask_mesh_region_force_torque_and_energy"]


def test_v26_public_axisymmetric_planar_depth_two_pi_r_force_normalization_coordinate_unit_mismatch():
    positions, _, _ = _quadratic_case()
    identity = _identity_v26(len(positions))
    identity["axisymmetric_planar_depth_two_pi_r_force_normalization_coordinate_unit_generation_identity"].update({
        "planar_depth_normalization_generation": "normalization-130",
        "radius_normalization_generation": "normalization-129",
        "result_planar_depth_m": 50.0, "result_planar_total_force_n": 1200.0,
        "result_axisymmetric_radius_m": 30.0,
        "result_axisymmetric_total_force_n": 1.2 / (2.0 * math.pi),
        "result_radius_measure_convention": "meridian_only",
        "result_coordinate_convention": "x_y_planar", "result_force_unit": "N_per_mm",
    })
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["axisymmetric_and_planar_force_share_depth_two_pi_r_coordinates_units_and_mesh"]
