"""Magnetostatic 2-D artifact-identity gates (v19..v28 chain, v45..v56).

One positive closure per gate; every negative test keeps its own failure
signal.  Payload builders live in ``_magnetostatic_2d_generalization_payloads``.
"""

from __future__ import annotations

import math
from copy import deepcopy

from _magnetostatic_2d_generalization_payloads import (
    BH_V46,
    BH_V47,
    CAPACITANCE_V55,
    CIRCUIT_V50,
    ELECTROSTATIC_V48,
    FORCE_V46,
    FORCE_V47,
    FORCE_V54,
    FORCE_V56,
    FROZEN_INDUCTANCE_V52,
    HARMONIC_V49,
    HEAT_V56,
    HYSTERESIS_V53,
    INCREMENTAL_V48,
    INCREMENTAL_V51,
    INDUCTION_V55,
    IRON_V50,
    MAGNETIC_PRESSURE_V52,
    NONLINEAR_V49,
    POWER_V54,
    VIRTUAL_WORK_V53,
    WEIGHTED_FORCE_V51,
    _gate,
    _identity_v19,
    _identity_v20,
    _identity_v21,
    _identity_v22,
    _identity_v23,
    _identity_v24,
    _identity_v25,
    _identity_v26,
    _identity_v27,
    _identity_v28,
    _identity_v45,
    _identity_v46,
    _identity_v47,
    _identity_v48,
    _identity_v49,
    _identity_v50,
    _identity_v51,
    _identity_v52,
    _identity_v53,
    _identity_v56,
    _payload_v54,
    _payload_v55,
    _quadratic_case,
)
from radia_mcp.radia_ngsolve.axisymmetric_v44_identity import (
    validate_public_identity as validate_public_identity_v45,
)
from radia_mcp.radia_ngsolve.axisymmetric_v46_identity import (
    validate_public_identity as validate_public_identity_v46,
)
from radia_mcp.radia_ngsolve.electromagnetic_artifact_identity_v49 import (
    validate_public_identity as validate_public_identity_v49,
)
from radia_mcp.radia_ngsolve.electromagnetic_artifact_identity_v50 import (
    validate_public_identity as validate_public_identity_v50,
)
from radia_mcp.radia_ngsolve.electromagnetic_artifact_identity_v51 import (
    validate_public_identity as validate_public_identity_v51,
)
from radia_mcp.radia_ngsolve.electromagnetic_artifact_identity_v52 import (
    validate_public_identity as validate_public_identity_v52,
)
from radia_mcp.radia_ngsolve.electromagnetic_artifact_identity_v53 import (
    validate_public_identity as validate_public_identity_v53,
)
from radia_mcp.radia_ngsolve.electromagnetic_artifact_identity_v54 import (
    validate_public_identity as validate_public_identity_v54,
)
from radia_mcp.radia_ngsolve.electromagnetic_artifact_identity_v55 import (
    validate_public_identity as validate_public_identity_v55,
)
from radia_mcp.radia_ngsolve.electromagnetic_force_heat_identity_v56 import (
    validate_public_identity as validate_public_identity_v56,
)
from radia_mcp.radia_ngsolve.electromagnetic_semantic_identity_v48 import (
    validate_public_identity as validate_public_identity_v48,
)
from radia_mcp.radia_ngsolve.magnetic_artifact_lineage_v47 import (
    validate_public_identity as validate_public_identity_v47,
)


# --- force_coenergy_displacement_gate (cumulative v19 .. v28 chain) --------


def test_v19_public_weighted_stress_tensor_air_mask_mesh_generation_mismatch():
    positions, _, _ = _quadratic_case()
    identity = _identity_v19(len(positions))
    identity[
        "weighted_stress_air_mask_nodal_weight_mesh_generation_identity"
    ].update(
        {
            "air_mask_mesh_generation": "mesh-20",
            "nodal_weight_mesh_generation": "mesh-20",
            "mask_air_region_ids": [2, 1],
            "force_weight_node_ids": [101, 103, 102],
            "force_air_mask_sha256": "f" * 64,
            "force_nodal_weight_sha256": "f" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "weighted_stress_air_mask_and_nodal_weights_use_current_mesh"
    ] is False


def test_v19_public_sliding_band_periodic_angle_rotor_position_generation_mismatch():
    positions, _, _ = _quadratic_case()
    identity = _identity_v19(len(positions))
    identity[
        "sliding_band_periodic_angle_rotor_position_generation_identity"
    ].update(
        {
            "sliding_band_rotor_position_generation": "rotor-position-20",
            "torque_periodic_angle_generation": "periodic-angle-20",
            "sliding_band_rotor_angle_deg": 10.0,
            "torque_periodic_angle_pairs_deg": [[30.0, 60.0], [0.0, 30.0]],
            "torque_sliding_band_map_sha256": "f" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "sliding_band_angles_and_torque_use_current_rotor_position"
    ] is False


def test_v20_public_coenergy_torque_angle_difference_remesh_state_generation_mismatch():
    positions, _, _ = _quadratic_case()
    identity = _identity_v20(len(positions))
    identity[
        "coenergy_torque_angle_difference_remesh_state_generation_identity"
    ].update(
        {
            "derivative_state_table_generation": "state-table-21",
            "mesh_remap_solve_generations": [
                "solve-21-a",
                "solve-21-b",
                "solve-21-c",
            ],
            "excitation_solve_generations": [
                "solve-22-a",
                "solve-21-b",
                "solve-22-c",
            ],
            "derivative_angles_deg": [14.0, 15.0, 17.0],
            "derivative_angle_spacing_generation": "angle-spacing-21",
            "derivative_state_table_sha256": "f" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "coenergy_torque_states_share_remesh_excitation_and_angle_generations"
    ] is False


def test_v20_public_axisymmetric_henrotte_hodge_radius_weight_coordinate_generation_mismatch():
    positions, _, _ = _quadratic_case()
    identity = _identity_v20(len(positions))
    identity[
        "axisymmetric_henrotte_hodge_radius_weight_coordinate_generation_identity"
    ].update(
        {
            "radius_weight_mesh_geometry_generation": "axi-mesh-21",
            "cylindrical_coordinate_mesh_geometry_generation": "axi-mesh-21",
            "hodge_radius_weight_m": [0.02, 0.03, 0.04],
            "cylindrical_r_coordinate_m": [0.03, 0.02, 0.01],
            "force_radius_weight_table_sha256": "f" * 64,
            "field_coordinate_table_sha256": "f" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "axisymmetric_hodge_radius_weights_use_current_mesh_coordinates"
    ] is False


def test_v21_public_weighted_stress_force_mask_material_mesh_generation_mismatch():
    positions, _, _ = _quadratic_case()
    identity = _identity_v21(len(positions))
    identity["weighted_stress_force_mask_material_mesh_generation_identity"].update(
        {
            "weighted_stress_mask_solve_generation": "magnetostatic-solve-30",
            "material_label_solve_generation": "magnetostatic-solve-29",
            "mesh_field_solve_generation": "magnetostatic-solve-28",
            "mask_body_group_ids": [4, 6],
            "resolved_material_labels": ["air", "magnet"],
            "reported_force_n": [10.1, -1.4],
            "force_mask_table_sha256": "a" * 64,
            "force_mesh_field_sha256": "b" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "weighted_stress_force_uses_current_mask_materials_and_mesh_field"
    ]


def test_v21_public_harmonic_loss_phase_frequency_lamination_generation_mismatch():
    positions, _, _ = _quadratic_case()
    identity = _identity_v21(len(positions))
    identity["harmonic_loss_phase_frequency_lamination_generation_identity"].update(
        {
            "phase_convention_analysis_generation": "harmonic-loss-30",
            "frequency_analysis_generation": "harmonic-loss-29",
            "lamination_analysis_generation": "harmonic-loss-28",
            "material_loss_analysis_generation": "harmonic-loss-27",
            "loss_phase_convention": "exp(-jwt)",
            "loss_frequency_hz": 60.0,
            "loss_lamination_orientations": ["stacking_z", "in_plane"],
            "loss_material_coefficients": [[0.01, 2.0], [0.02, 1.6]],
            "loss_material_table_sha256": "c" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "harmonic_loss_uses_current_phase_frequency_lamination_and_material_data"
    ]


def test_v22_public_axisymmetric_force_energy_two_pi_r_depth_normalization_generation_mismatch():
    positions, _, _ = _quadratic_case()
    identity = _identity_v22(len(positions))
    identity[
        "axisymmetric_force_energy_measure_depth_coordinate_generation_identity"
    ].update(
        {
            "measure_solve_generation": "axisym-solve-40",
            "depth_solve_generation": "axisym-solve-39",
            "coordinate_solve_generation": "axisym-solve-38",
            "result_problem_type": "planar",
            "result_measure_convention": "planar_depth",
            "result_planar_depth_m": 0.05,
            "result_coordinate_convention": "x_y",
            "reported_force_energy_values": [0.625, 0.00155],
            "result_normalization_table_sha256": "9" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "axisymmetric_force_energy_uses_current_measure_depth_and_coordinates"
    ]


def test_v22_public_nonlinear_incremental_permeability_force_perturbation_branch_generation_mismatch():
    positions, _, _ = _quadratic_case()
    identity = _identity_v22(len(positions))
    identity[
        "nonlinear_incremental_mu_force_branch_perturbation_generation_identity"
    ].update(
        {
            "branch_operating_point_generation": "nonlinear-op-40",
            "differential_mu_operating_point_generation": "nonlinear-op-39",
            "perturbation_operating_point_generation": "nonlinear-op-38",
            "force_branch_id": "down-sweep:17",
            "force_perturbation_current_a": 0.1,
            "force_differential_mu_sha256": "a" * 64,
            "reported_incremental_force_n": [-0.08, 0.04],
            "force_incremental_state_sha256": "b" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "nonlinear_incremental_force_uses_current_branch_mu_and_perturbation"
    ]


def test_v23_public_nonlinear_bh_branch_operating_point_force_generation_mismatch():
    positions, _, _ = _quadratic_case()
    identity = _identity_v23(len(positions))
    identity[
        "nonlinear_bh_branch_operating_point_force_mesh_generation_identity"
    ].update(
        {
            "branch_solve_generation": "nonlinear-force-50",
            "operating_point_solve_generation": "nonlinear-force-49",
            "permeability_solve_generation": "nonlinear-force-48",
            "force_mesh_solve_generation": "nonlinear-force-47",
            "force_result_solve_generation": "nonlinear-force-46",
            "force_branch_id": "descending:23",
            "force_operating_point_current_a": 6.0,
            "force_operating_point_flux_density_t": [1.1, 1.2],
            "force_permeability_state_sha256": "b" * 64,
            "integrated_force_mesh_sha256": "c" * 64,
            "reported_force_n": [15.0, 0.7],
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "nonlinear_force_uses_current_bh_branch_operating_point_mu_and_mesh"
    ]


def test_v23_public_sliding_band_angle_mesh_harmonic_torque_generation_mismatch():
    positions, _, _ = _quadratic_case()
    identity = _identity_v23(len(positions))
    identity["sliding_band_angle_mesh_harmonic_torque_generation_identity"].update(
        {
            "angle_sweep_generation": "sliding-sweep-50",
            "airgap_mesh_sweep_generation": "sliding-sweep-49",
            "phase_current_sweep_generation": "sliding-sweep-48",
            "torque_sample_sweep_generation": "sliding-sweep-47",
            "harmonic_sweep_generation": "sliding-sweep-46",
            "torque_rotor_angles_deg": [0.0, 10.0, 5.0, 15.0],
            "torque_airgap_mesh_sha256": "d" * 64,
            "torque_phase_current_table_sha256": "e" * 64,
            "harmonic_torque_samples_nm": [1.0, 0.9, 1.2, 1.1],
            "reported_harmonic_orders": [0, 2, 1],
            "reported_harmonic_amplitudes_nm": [1.05, 0.04, 0.12],
            "harmonic_sample_table_sha256": "f" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "sliding_band_harmonics_use_current_angles_mesh_currents_and_samples"
    ]


def test_v24_public_weighted_stress_force_energy_derivative_mesh_frame_unit_mismatch():
    positions, _, _ = _quadratic_case()
    identity = _identity_v24(len(positions))
    identity[
        "weighted_stress_energy_derivative_force_mesh_frame_unit_generation_identity"
    ].update(
        {
            "energy_derivative_force_generation": "force-100",
            "mesh_force_generation": "force-99",
            "displacement_frame_force_generation": "force-98",
            "unit_force_generation": "force-97",
            "energy_derivative_mesh_sha256": "b" * 64,
            "energy_derivative_displacement_frame_id": "rotor:x",
            "energy_derivative_displacement_unit": "mm",
            "energy_derivative_force_unit": "N/mm",
            "energy_derivative_force_n": [0.0124, -0.0002],
            "energy_derivative_result_sha256": "c" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "weighted_stress_and_energy_derivative_share_mesh_frame_units_and_generation"
    ]


def test_v24_public_axisymmetric_revolved_energy_force_2pir_jacobian_derham_generation_mismatch():
    positions, _, _ = _quadratic_case()
    identity = _identity_v24(len(positions))
    identity[
        "axisymmetric_revolved_energy_force_2pir_jacobian_derham_generation_identity"
    ].update(
        {
            "jacobian_axisymmetric_generation": "axisym-100",
            "field_axisymmetric_generation": "axisym-99",
            "material_axisymmetric_generation": "axisym-98",
            "mesh_axisymmetric_generation": "axisym-97",
            "revolved_result_axisymmetric_generation": "axisym-96",
            "revolved_jacobian_measure": "r",
            "revolved_field_state_sha256": "d" * 64,
            "revolved_material_map_sha256": "e" * 64,
            "revolved_source_mesh_sha256": "f" * 64,
            "revolution_angle_deg": 180.0,
            "revolved_energy_j": 0.0625,
            "revolved_force_n": [2.1, 0.0],
            "revolved_derham_sequence_id": "plain-H1",
            "revolved_result_sha256": "0" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "axisymmetric_revolved_energy_force_share_2pir_hodge_field_material_and_mesh"
    ]


def test_v25_public_nonlinear_bh_incremental_energy_coenergy_force_branch_mesh_mismatch():
    positions, _, _ = _quadratic_case()
    identity = _identity_v25(len(positions))
    identity[
        "nonlinear_bh_incremental_energy_coenergy_force_branch_mesh_generation_identity"
    ].update(
        {
            "bh_curve_nonlinear_generation": "nonlinear-110",
            "material_nonlinear_generation": "nonlinear-109",
            "branch_nonlinear_generation": "nonlinear-108",
            "incremental_state_nonlinear_generation": "nonlinear-107",
            "mesh_nonlinear_generation": "nonlinear-106",
            "energy_nonlinear_generation": "nonlinear-105",
            "coenergy_nonlinear_generation": "nonlinear-104",
            "force_nonlinear_generation": "nonlinear-103",
            "result_nonlinear_material_ids": ["core-old"],
            "result_bh_curve_sha256": "d" * 64,
            "result_material_map_sha256": "e" * 64,
            "result_branch_id": "descending:step-5",
            "result_load_current_a": 6.0,
            "result_incremental_state_sha256": "f" * 64,
            "result_mesh_sha256": "0" * 64,
            "result_magnetic_energy_j": 0.31,
            "result_magnetic_coenergy_j": 0.45,
            "result_incremental_force_n": [9.8, 0.2],
            "accepted_result_sha256": "1" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "nonlinear_bh_incremental_force_uses_current_branch_material_mesh_energy_and_coenergy"
    ]


def test_v25_public_open_boundary_domain_decay_multipole_moment_material_generation_mismatch():
    positions, _, _ = _quadratic_case()
    identity = _identity_v25(len(positions))
    identity[
        "open_boundary_domain_decay_multipole_moment_material_generation_identity"
    ].update(
        {
            "domain_boundary_generation": "open-boundary-110",
            "mesh_boundary_generation": "open-boundary-109",
            "material_boundary_generation": "open-boundary-108",
            "multipole_boundary_generation": "open-boundary-107",
            "decay_boundary_generation": "open-boundary-106",
            "result_boundary_type": "dirichlet_zero",
            "result_source_radius_m": 0.08,
            "result_outer_radius_m": 0.15,
            "result_multipole_order": 1,
            "result_material_map_sha256": "2" * 64,
            "result_mesh_sha256": "3" * 64,
            "result_multipole_moment_sha256": "4" * 64,
            "result_decay_sample_radii_m": [0.1, 0.12, 0.15],
            "result_decay_flux_density_t": [1.0e-3, 1.2e-3, 1.1e-3],
            "accepted_result_sha256": "5" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "open_boundary_uses_current_domain_decay_multipole_material_and_mesh"
    ]


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


def test_v27_public_nonlinear_bh_minor_loop_branch_interpolation_state_coenergy_force_mismatch():
    positions, _, _ = _quadratic_case()
    identity = _identity_v27(len(positions))
    identity[
        "nonlinear_bh_minor_loop_branch_interpolation_state_coenergy_force_generation_identity"
    ].update({
        "branch_nonlinear_generation": "bh-minor-loop-140",
        "state_nonlinear_generation": "bh-minor-loop-139",
        "result_bh_branch": "descending-major-loop",
        "result_interpolation_rule": "linear-b",
        "result_state_point_am": [-120.0, 700.0],
        "result_magnetic_coenergy_j": 0.28,
        "result_force_n": [12.0, 1.0],
        "result_bh_table_sha256": "a" * 64,
        "result_mesh_sha256": "b" * 64,
        "accepted_solution_sha256": "c" * 64,
    })
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "nonlinear_force_uses_current_minor_loop_branch_interpolation_state_coenergy_and_mesh"
    ]


def test_v27_public_harmonic_eddy_phasor_convention_conductivity_skin_depth_loss_mesh_mismatch():
    positions, _, _ = _quadratic_case()
    identity = _identity_v27(len(positions))
    identity[
        "harmonic_eddy_phasor_conductivity_skin_depth_frequency_loss_mesh_generation_identity"
    ].update({
        "phasor_eddy_generation": "harmonic-eddy-140",
        "conductivity_eddy_generation": "harmonic-eddy-139",
        "mesh_eddy_generation": "harmonic-eddy-138",
        "result_phasor_convention": "exp(-jwt)",
        "result_conductivity_s_m": 5.8e4,
        "result_frequency_hz": 50.0,
        "result_skin_depth_m": 0.02,
        "result_minimum_elements_per_skin_depth": 0.5,
        "result_joule_loss_w": 4.2,
        "result_mesh_sha256": "d" * 64,
        "accepted_loss_result_sha256": "e" * 64,
    })
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "harmonic_eddy_loss_uses_current_phasor_conductivity_skin_depth_frequency_and_mesh"
    ]


def test_v28_public_positive_axisymmetric_and_pm_identity():
    positions, _, _ = _quadratic_case()
    assert _gate(_identity_v28(len(positions)))["status"] == "ok"


def test_v28_public_axisymmetric_aphi_radial_weighting_region_force_energy_mesh_solution_mismatch():
    positions, _, _ = _quadratic_case()
    identity = _identity_v28(len(positions))
    identity["axisymmetric_aphi_radial_weight_region_energy_force_mesh_solution_generation_identity"].update(
        {"aphi_axisymmetric_generation": "axisym-force-150", "result_formulation": "Az",
         "result_radial_weighting": "1", "result_region_ids": [5, 6],
         "result_magnetic_coenergy_j": [0.49, 0.5, 0.52], "result_force_axis": "r",
         "result_force_from_energy_n": 55.0, "accepted_solution_sha256": "a" * 64}
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["axisymmetric_force_uses_current_aphi_weight_regions_energy_axis_mesh_and_solution"]


def test_v28_public_permanent_magnet_recoil_line_temperature_operating_point_demag_margin_force_mismatch():
    positions, _, _ = _quadratic_case()
    identity = _identity_v28(len(positions))
    identity["permanent_magnet_recoil_temperature_operating_point_demag_force_generation_identity"].update(
        {"recoil_magnet_generation": "pm-operating-point-150", "result_recoil_relative_permeability": 1.2,
         "result_remanence_t": 1.0, "result_magnet_temperature_c": 20.0,
         "result_operating_point_bh": [0.45, -500000.0], "result_demag_margin_a_m": -10000.0,
         "result_force_n": [8.0, -1.0], "accepted_result_sha256": "d" * 64}
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["permanent_magnet_force_uses_current_recoil_temperature_operating_point_frame_demag_and_mesh"]


# --- axisymmetric_v44_identity (v45) ----------------------------------------


def test_v45_public_magnetostatic_identity_accepts_closed_artifacts():
    checks = validate_public_identity_v45(_identity_v45())
    assert checks and all(checks.values())


def test_v45_public_magnetostatic_identity_rejects_force_owner_mutation():
    identity = _identity_v45()
    identity["v45_public_axisymmetric_force_torque_coenergy_stress_contour_owner_mismatch"]["result_mesh_owner"] = "stale:mesh"
    checks = validate_public_identity_v45(identity)
    assert checks and not all(checks.values())


def test_v45_public_magnetostatic_identity_rejects_fringe_energy_mutation():
    identity = _identity_v45()
    identity["v45_public_electrostatic_fringe_charge_energy_capacitance_interface_flux_axisfactor_mismatch"]["result_stored_energy_j"] = -1.0
    checks = validate_public_identity_v45(identity)
    assert checks and not all(checks.values())


# --- axisymmetric_v46_identity (v46) ----------------------------------------


def test_v46_public_magnetostatic_identity_accepts_closed_artifacts():
    checks = validate_public_identity_v46(_identity_v46())
    assert checks and all(checks.values())


def test_v46_public_magnetostatic_identity_rejects_force_partial_nan_mutation():
    identity = _identity_v46()
    identity[FORCE_V46]["result_partial_solve_status"] = "partial"
    identity[FORCE_V46]["result_force_n"] = [float("nan"), -3.0]
    checks = validate_public_identity_v46(identity)
    assert checks and not all(checks.values())


def test_v46_public_magnetostatic_identity_rejects_bh_unit_convergence_mutation():
    identity = _identity_v46()
    identity[BH_V46]["result_bh_unit"] = "gauss_oersted"
    identity[BH_V46]["result_convergence_status"] = "not_converged"
    checks = validate_public_identity_v46(identity)
    assert checks and not all(checks.values())


# --- magnetic_artifact_lineage_v47 (v47) ------------------------------------


def test_v47_positive_replays_are_accepted() -> None:
    assert all(validate_public_identity_v47(_identity_v47()).values())


def test_v47_force_owner_sign_pair_mutation_is_rejected() -> None:
    identity = _identity_v47()
    identity[FORCE_V47]["result_force_method"] = "contour_stress"
    identity[FORCE_V47]["result_body_owner"] = "group:fixed"
    identity[FORCE_V47]["result_force_sign_convention"] = "negative_displacement_direction"
    identity[FORCE_V47]["result_displacement_pair_m"] = [0.001, 0.0]
    assert not all(validate_public_identity_v47(identity).values())


def test_v47_bh_row_branch_history_mutation_is_rejected() -> None:
    identity = _identity_v47()
    identity[BH_V47]["result_operating_point_row_keys"] = ["current=3", "current=1", "current=2"]
    identity[BH_V47]["result_hysteresis_branches"] = ["descending", "ascending", "ascending"]
    identity[BH_V47]["result_excitation_history_sha256"] = "a" * 64
    assert not all(validate_public_identity_v47(identity).values())


# --- electromagnetic_semantic_identity_v48 (v48) ----------------------------


def test_v48_positive_incremental_and_electrostatic_artifacts_are_accepted() -> None:
    assert all(validate_public_identity_v48(_identity_v48()).values())


def test_v48_incremental_owner_and_phasor_mutations_are_rejected() -> None:
    identity = deepcopy(_identity_v48())
    identity[INCREMENTAL_V48]["result_harmonic_phasor_a"] = [1.5, 0.25]
    identity[INCREMENTAL_V48]["result_operating_point_owner"] = "operating-point:old"
    checks = validate_public_identity_v48(identity)
    assert checks["v48_incremental_bias_phasor_tangent_owner"] is False


def test_v48_electrostatic_row_and_owner_mutations_are_rejected() -> None:
    identity = deepcopy(_identity_v48())
    identity[ELECTROSTATIC_V48]["result_charge_c"] = [0.0, 2.0e-9, 6.0e-9, 4.0e-9]
    identity[ELECTROSTATIC_V48]["result_conductor_owner_rows"][2] = "conductor:electrode-b"
    checks = validate_public_identity_v48(identity)
    assert checks["v48_electrostatic_charge_energy_voltage_owner_closure"] is False


# --- electromagnetic_artifact_identity_v49 (v49) ----------------------------


def test_v49_positive_nonlinear_and_harmonic_artifacts_are_accepted() -> None:
    assert all(validate_public_identity_v49(_identity_v49()).values())


def test_v49_nonlinear_branch_temperature_and_owner_mutations_are_rejected() -> None:
    identity = deepcopy(_identity_v49())
    identity[NONLINEAR_V49]["result_bh_branch"] = "descending:minor-loop"
    identity[NONLINEAR_V49]["result_temperature_c"] = 20.0
    identity[NONLINEAR_V49]["result_material_owner"] = "material:old"
    assert validate_public_identity_v49(identity)["v49_nonlinear_bh_incremental_temperature_lamination_owner"] is False


def test_v49_harmonic_geometry_phase_loss_and_owner_mutations_are_rejected() -> None:
    identity = deepcopy(_identity_v49())
    identity[HARMONIC_V49]["result_geometry_type"] = "axisymmetric"
    identity[HARMONIC_V49]["result_phase_convention"] = "exp(-jwt)"
    identity[HARMONIC_V49]["result_losses_w"] = {"joule_w": 1.0, "core_w": 8.0}
    identity[HARMONIC_V49]["accepted_result_owner"] = "harmonic-result:old"
    assert validate_public_identity_v49(identity)["v49_harmonic_geometry_depth_frequency_phase_circuit_loss_owner"] is False


# --- electromagnetic_artifact_identity_v50 (v50) ----------------------------


def test_v50_positive_harmonic_artifacts_are_accepted() -> None:
    assert all(validate_public_identity_v50(_identity_v50()).values())


def test_v50_circuit_and_iron_mutations_are_rejected() -> None:
    identity = deepcopy(_identity_v50())
    identity[CIRCUIT_V50]["result_parallel_paths"] = 2
    identity[IRON_V50]["result_frequency_hz"] = 50.0
    assert not all(validate_public_identity_v50(identity).values())


def test_v50_self_consistent_invalid_waveform_and_loss_closure_are_rejected() -> None:
    identity = deepcopy(_identity_v50())
    identity[IRON_V50]["flux_density_waveform_t"] = identity[IRON_V50]["result_flux_density_waveform_t"] = [0.0, 1.0, 0.0]
    identity[IRON_V50]["loss_components"] = identity[IRON_V50]["result_loss_components"] = {
        "hysteresis_w": 1.8, "eddy_w": 0.9, "excess_w": 0.3, "total_w": 8.0
    }
    assert validate_public_identity_v50(identity)["v50_steinmetz_frequency_waveform_components_material_owner"] is False


# --- electromagnetic_artifact_identity_v51 (v51) ----------------------------


def test_v51_positive_public_artifacts_are_accepted() -> None:
    assert all(validate_public_identity_v51(_identity_v51()).values())


def test_v51_frozen_counterfactuals_are_rejected() -> None:
    identity = deepcopy(_identity_v51())
    identity[INCREMENTAL_V51].update({"result_analysis_mode": "frozen_permeability", "result_branch_state": "recoil_branch"})
    identity[WEIGHTED_FORCE_V51].update({"result_air_element_ids": [101, 105], "result_force_frame": "local_xy"})
    assert not all(validate_public_identity_v51(identity).values())


def test_v51_self_consistent_wrong_physics_are_rejected() -> None:
    identity = deepcopy(_identity_v51())
    identity[INCREMENTAL_V51]["analysis_mode"] = identity[INCREMENTAL_V51]["result_analysis_mode"] = "frozen_permeability"
    identity[WEIGHTED_FORCE_V51]["axisymmetric_factor_m"] = identity[WEIGHTED_FORCE_V51]["result_axisymmetric_factor_m"] = 1.0
    assert not all(validate_public_identity_v51(identity).values())


# --- electromagnetic_artifact_identity_v52 (v52) ----------------------------


def test_v52_positive_public_artifacts_are_accepted() -> None:
    assert all(validate_public_identity_v52(_identity_v52()).values())


def test_v52_frozen_counterfactuals_are_rejected() -> None:
    identity = deepcopy(_identity_v52())
    identity[MAGNETIC_PRESSURE_V52]["result_boundary_normal"] = [-0.6, -0.8, 0.0]
    identity[FROZEN_INDUCTANCE_V52]["result_bias_point"] = {"current_a": 0.0, "solution_sha256": "a" * 64}
    assert not all(validate_public_identity_v52(identity).values())


def test_v52_self_consistent_wrong_physics_are_rejected() -> None:
    identity = deepcopy(_identity_v52())
    identity[MAGNETIC_PRESSURE_V52]["boundary_normal"] = identity[MAGNETIC_PRESSURE_V52]["result_boundary_normal"] = [2.0, 0.0, 0.0]
    identity[FROZEN_INDUCTANCE_V52]["perturbation"] = identity[FROZEN_INDUCTANCE_V52]["result_perturbation"] = {"delta_current_a": 2.0, "frequency_hz": 400.0}
    assert not all(validate_public_identity_v52(identity).values())


# --- electromagnetic_artifact_identity_v53 (v53) ----------------------------


def test_v53_positive_public_artifacts_are_accepted():
    assert all(validate_public_identity_v53(_identity_v53()).values())


def test_v53_frozen_counterfactuals_are_rejected():
    identity = deepcopy(_identity_v53())
    identity[HYSTERESIS_V53]["result_phasor_convention"] = "exp(-j_omega_t)"
    identity[VIRTUAL_WORK_V53]["result_constraint_mode"] = "fixed_charge"
    assert not all(validate_public_identity_v53(identity).values())


def test_v53_self_consistent_wrong_physics_is_rejected():
    identity = deepcopy(_identity_v53())
    identity[HYSTERESIS_V53]["relative_permeability"]["imag"] = identity[HYSTERESIS_V53]["result_relative_permeability"]["imag"] = 18.0
    identity[VIRTUAL_WORK_V53]["force_n"] = identity[VIRTUAL_WORK_V53]["result_force_n"] = -4.0
    assert not all(validate_public_identity_v53(identity).values())


# --- electromagnetic_artifact_identity_v54 (v54) ----------------------------


def test_v54_positive_identities_are_accepted():
    assert all(validate_public_identity_v54(_payload_v54()).values())


def test_v54_frozen_mutations_are_rejected():
    payload = deepcopy(_payload_v54())
    payload[POWER_V54]["result_active_power_w"] = 80.0
    payload[FORCE_V54]["result_radius_weighting"] = "planar"
    assert not all(validate_public_identity_v54(payload).values())


def test_v54_self_consistent_nonphysical_records_are_rejected():
    payload = deepcopy(_payload_v54())
    payload[POWER_V54]["loss_components_w"] = payload[POWER_V54]["result_loss_components_w"] = {"copper": 10.0}
    payload[FORCE_V54]["force_direction_rz"] = payload[FORCE_V54]["result_force_direction_rz"] = [2.0, 0.0]
    assert not all(validate_public_identity_v54(payload).values())


def test_v54_malformed_values_reject_without_raising():
    payload = deepcopy(_payload_v54())
    payload[POWER_V54]["loss_components_w"] = {"copper": [80.0]}
    payload[FORCE_V54]["integration_selection"] = [["block:armature"]]
    assert not all(validate_public_identity_v54(payload).values())


# --- electromagnetic_artifact_identity_v55 (v55) ----------------------------


def test_v55_positive_identities_are_accepted():
    assert all(validate_public_identity_v55(_payload_v55()).values())


def test_v55_frozen_mutations_are_rejected():
    payload = deepcopy(_payload_v55())
    payload[INDUCTION_V55]["result_skin_depth_m"] = 0.1
    payload[CAPACITANCE_V55]["result_conductor_order"] = ["conductor:2", "conductor:1"]
    assert not all(validate_public_identity_v55(payload).values())


def test_v55_self_consistent_wrong_skin_depth_or_negative_loss_is_rejected():
    payload = deepcopy(_payload_v55())
    payload[INDUCTION_V55]["skin_depth_m"] = payload[INDUCTION_V55]["result_skin_depth_m"] = 0.1
    payload[INDUCTION_V55]["joule_loss_w"] = payload[INDUCTION_V55]["result_joule_loss_w"] = -1.0
    assert not all(validate_public_identity_v55(payload).values())


def test_v55_self_consistent_nonsymmetric_capacitance_or_bad_charge_is_rejected():
    payload = deepcopy(_payload_v55())
    bad_matrix = [[2.0e-12, 1.0e-12], [-1.0e-12, 2.0e-12]]
    payload[CAPACITANCE_V55]["capacitance_matrix_f"] = payload[CAPACITANCE_V55]["result_capacitance_matrix_f"] = bad_matrix
    payload[CAPACITANCE_V55]["charge_c"] = payload[CAPACITANCE_V55]["result_charge_c"] = [9.0, 9.0]
    assert not all(validate_public_identity_v55(payload).values())


def test_v55_numeric_sha256_values_are_rejected():
    payload = _payload_v55()
    numeric_digest = int("9" * 64)
    for row in payload.values():
        row["result_sha256"] = numeric_digest
        row["accepted_result_sha256"] = numeric_digest
    assert not all(validate_public_identity_v55(payload).values())


# --- electromagnetic_force_heat_identity_v56 (v56) --------------------------


def test_v56_positive_identity_is_accepted() -> None:
    assert all(validate_public_identity_v56(_identity_v56()).values())


def test_v56_frozen_result_mutations_are_rejected() -> None:
    identity = deepcopy(_identity_v56())
    identity[FORCE_V56]["result_force_n"] = -5.0
    identity[HEAT_V56]["result_outward_boundary_flux_w"] = -10.0
    assert not all(validate_public_identity_v56(identity).values())


def test_v56_self_consistent_physics_contradictions_are_rejected() -> None:
    identity = deepcopy(_identity_v56())
    identity[FORCE_V56]["force_n"] = identity[FORCE_V56]["result_force_n"] = 3.0
    identity[HEAT_V56]["outward_boundary_flux_w"] = identity[HEAT_V56]["result_outward_boundary_flux_w"] = 40.0
    identity[HEAT_V56]["balance_residual_w"] = identity[HEAT_V56]["result_balance_residual_w"] = 10.0
    assert not all(validate_public_identity_v56(identity).values())


def test_v56_malformed_samples_reject_without_raising() -> None:
    identity = deepcopy(_identity_v56())
    identity[FORCE_V56]["coenergy_samples"] = [[0.0]]
    assert not all(validate_public_identity_v56(identity).values())


def test_v56_numeric_digests_are_rejected() -> None:
    identity = deepcopy(_identity_v56())
    numeric_digest = int("1" * 64)
    for contract_name in (FORCE_V56, HEAT_V56):
        identity[contract_name]["result_sha256"] = numeric_digest
        identity[contract_name]["accepted_result_sha256"] = numeric_digest
    assert not all(validate_public_identity_v56(identity).values())
