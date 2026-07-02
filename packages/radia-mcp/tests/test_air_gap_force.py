"""Air-gap Maxwell pressure / holding-force helpers."""

import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.force import (  # noqa: E402
    MU0,
    air_gap_force_summary,
    air_gap_holding_force,
    air_gap_maxwell_pressure,
    air_gap_shear_stress,
    air_gap_shear_torque,
    air_gap_shear_torque_from_angle_samples,
    air_gap_shear_torque_summary,
    force_moment_resultant_summary,
    magnetic_field_probe_result_package_gate,
    maxwell_contour_force_2d,
    maxwell_contour_segment_balance_summary_2d,
    maxwell_line_segment_force_2d,
    maxwell_stress_tensor_air,
    maxwell_traction_air,
    maxwell_traction_summary,
    parallel_wire_force_result_package_gate,
    parallel_wire_lorentz_force_summary,
    parallel_wire_virtual_work_force_summary,
)


def test_air_gap_pressure_matches_maxwell_stress_at_one_tesla():
    expected = 1.0 / (2.0 * MU0)
    assert air_gap_maxwell_pressure(1.0) == pytest.approx(expected)
    assert air_gap_maxwell_pressure(-1.0) == pytest.approx(expected)


def test_maxwell_tensor_normal_field_reduces_to_air_gap_pressure():
    pressure = air_gap_maxwell_pressure(1.0)
    tensor = maxwell_stress_tensor_air((0.0, 0.0, 1.0))
    traction = maxwell_traction_air((0.0, 0.0, 1.0), (0.0, 0.0, 1.0))
    summary = maxwell_traction_summary((0.0, 0.0, 1.0), (0.0, 0.0, 2.0), area_m2=2.0e-4)

    assert tensor[0][0] == pytest.approx(-pressure)
    assert tensor[1][1] == pytest.approx(-pressure)
    assert tensor[2][2] == pytest.approx(pressure)
    assert traction == pytest.approx([0.0, 0.0, pressure])
    assert summary["normal_traction_Pa"] == pytest.approx(pressure)
    assert summary["normal_traction_identity_Pa"] == pytest.approx(pressure)
    assert summary["force_N"] == pytest.approx([0.0, 0.0, pressure * 2.0e-4])


def test_maxwell_tensor_tangential_field_is_magnetic_tension():
    pressure = air_gap_maxwell_pressure(1.0)
    traction = maxwell_traction_air((1.0, 0.0, 0.0), (0.0, 0.0, 1.0))
    summary = maxwell_traction_summary((1.0, 0.0, 0.0), (0.0, 0.0, 1.0))

    assert traction == pytest.approx([0.0, 0.0, -pressure])
    assert summary["B_normal_T"] == pytest.approx(0.0)
    assert summary["B_tangent_T"] == pytest.approx(1.0)
    assert summary["normal_traction_Pa"] == pytest.approx(-pressure)
    assert summary["tangential_traction_magnitude_Pa"] == pytest.approx(0.0)


def test_maxwell_traction_oblique_field_decomposes_into_normal_and_tangent():
    # B=(3,4,0), n=x gives Bn=3, Bt=4:
    # traction.n=(9-16)/(2mu), tangent traction=(Bn*Bt)/mu in y.
    summary = maxwell_traction_summary((3.0, 4.0, 0.0), (1.0, 0.0, 0.0))

    assert summary["B_normal_T"] == pytest.approx(3.0)
    assert summary["B_tangent_T"] == pytest.approx(4.0)
    assert summary["normal_traction_Pa"] == pytest.approx(-3.5 / MU0)
    assert summary["normal_traction_identity_Pa"] == pytest.approx(-3.5 / MU0)
    assert summary["tangential_traction_Pa"] == pytest.approx([0.0, 12.0 / MU0, 0.0])
    assert summary["tangential_traction_magnitude_Pa"] == pytest.approx(12.0 / MU0)


def test_air_gap_shear_stress_matches_maxwell_tangential_traction():
    Br = 0.8
    Bt = 0.1
    shear = air_gap_shear_stress(Br, Bt)
    traction = maxwell_traction_summary((Br, Bt, 0.0), (1.0, 0.0, 0.0))

    assert shear == pytest.approx(Br * Bt / MU0)
    assert traction["tangential_traction_Pa"] == pytest.approx([0.0, shear, 0.0])
    assert traction["tangential_traction_magnitude_Pa"] == pytest.approx(abs(shear))


def test_air_gap_shear_torque_scales_with_radius_length_angle_and_sign():
    Br = 0.8
    Bt = 0.1
    radius = 0.05
    length = 0.1
    full = air_gap_shear_torque(Br, Bt, radius, axial_length_m=length)
    half = air_gap_shear_torque(Br, Bt, radius, axial_length_m=length, angle_rad=math.pi)
    reverse = air_gap_shear_torque(Br, -Bt, radius, axial_length_m=length)
    summary = air_gap_shear_torque_summary(Br, Bt, radius, axial_length_m=length)

    assert full == pytest.approx(100.0)
    assert half == pytest.approx(50.0)
    assert reverse == pytest.approx(-100.0)
    assert summary["surface_area_m2"] == pytest.approx(2.0 * math.pi * radius * length)
    assert summary["tangential_force_N"] == pytest.approx(full / radius)
    assert summary["torque_Nm"] == pytest.approx(full)
    assert summary["torque_per_axial_length_N"] == pytest.approx(full / length)


def test_air_gap_sampled_shear_torque_uniform_matches_closed_form():
    Br = 0.8
    Bt = 0.1
    radius = 0.05
    length = 0.1
    angles = [0.0, 0.5 * math.pi, math.pi, 1.5 * math.pi]
    summary = air_gap_shear_torque_from_angle_samples(
        angles,
        [Br] * len(angles),
        [Bt] * len(angles),
        radius,
        axial_length_m=length,
    )

    expected = air_gap_shear_torque(Br, Bt, radius, axial_length_m=length)
    assert summary["n_segments"] == len(angles)
    assert summary["integrated_angle_rad"] == pytest.approx(2.0 * math.pi)
    assert summary["average_shear_stress_Pa"] == pytest.approx(Br * Bt / MU0)
    assert summary["torque_Nm"] == pytest.approx(expected)
    assert summary["tangential_force_N"] == pytest.approx(expected / radius)


def test_air_gap_sampled_shear_torque_sinusoidal_average():
    samples = 360
    harmonic = 2
    phase = math.pi / 3.0
    Br0 = 0.9
    Bt0 = 0.2
    radius = 0.04
    length = 0.12
    angles = [2.0 * math.pi * index / samples for index in range(samples)]
    br = [Br0 * math.cos(harmonic * angle) for angle in angles]
    bt = [Bt0 * math.cos(harmonic * angle + phase) for angle in angles]

    summary = air_gap_shear_torque_from_angle_samples(
        angles,
        br,
        bt,
        radius,
        axial_length_m=length,
    )
    expected_average_shear = 0.5 * Br0 * Bt0 * math.cos(phase) / MU0
    expected_torque = expected_average_shear * radius * radius * length * 2.0 * math.pi

    assert summary["average_shear_stress_Pa"] == pytest.approx(expected_average_shear, rel=1.0e-12)
    assert summary["torque_Nm"] == pytest.approx(expected_torque, rel=1.0e-12)
    assert summary["torque_per_axial_length_N"] == pytest.approx(expected_torque / length, rel=1.0e-12)


def test_motor_air_gap_harmonic_torque_phase_gate():
    samples = 720
    harmonic = 4
    Br0 = 0.82
    Bt0 = 0.18
    radius = 0.045
    length = 0.08
    angles = [2.0 * math.pi * index / samples for index in range(samples)]

    def torque_for_phase(phase):
        br = [Br0 * math.cos(harmonic * angle) for angle in angles]
        bt = [Bt0 * math.cos(harmonic * angle + phase) for angle in angles]
        return air_gap_shear_torque_from_angle_samples(
            angles,
            br,
            bt,
            radius,
            axial_length_m=length,
        )["torque_Nm"]

    expected = 0.5 * Br0 * Bt0 / MU0 * radius * radius * length * 2.0 * math.pi
    positive = torque_for_phase(0.0)
    quadrature = torque_for_phase(0.5 * math.pi)
    negative = torque_for_phase(math.pi)

    assert positive == pytest.approx(expected, rel=1.0e-12)
    assert abs(quadrature) < 1.0e-10
    assert negative == pytest.approx(-expected, rel=1.0e-12)


def test_force_moment_resultant_summary_handles_2d_force_couple():
    radius = 0.05
    force = 10.0
    summary = force_moment_resultant_summary(
        [(radius, 0.0), (-radius, 0.0)],
        [(0.0, force), (0.0, -force)],
    )
    shifted = force_moment_resultant_summary(
        [(radius, 0.0), (-radius, 0.0)],
        [(0.0, force), (0.0, -force)],
        pivot_m=(0.2, -0.1),
    )

    assert summary["dimension"] == 2
    assert summary["total_force"] == pytest.approx([0.0, 0.0])
    assert summary["total_moment"] == pytest.approx(2.0 * radius * force)
    assert summary["total_moment_magnitude"] == pytest.approx(2.0 * radius * force)
    assert shifted["total_moment"] == pytest.approx(summary["total_moment"])


def test_force_moment_resultant_summary_handles_3d_single_force():
    summary = force_moment_resultant_summary(
        [(0.0, 0.2, 0.0)],
        [(3.0, 0.0, 0.0)],
    )

    assert summary["dimension"] == 3
    assert summary["total_force"] == pytest.approx([3.0, 0.0, 0.0])
    assert summary["total_moment"] == pytest.approx([0.0, 0.0, -0.6])
    assert summary["total_moment_magnitude"] == pytest.approx(0.6)


def test_maxwell_line_segment_force_2d_matches_air_gap_pressure():
    pressure = air_gap_maxwell_pressure(1.0)
    row = maxwell_line_segment_force_2d(
        (0.0, -0.5),
        (0.0, 0.5),
        (1.0, 0.0),
        normal_side="right",
    )

    assert row["length_m"] == pytest.approx(1.0)
    assert row["unit_normal"] == pytest.approx([1.0, 0.0])
    assert row["traction_N_per_m2"] == pytest.approx([pressure, 0.0])
    assert row["force_per_depth_N_per_m"] == pytest.approx([pressure, 0.0])
    assert row["normal_force_per_depth_N_per_m"] == pytest.approx(pressure)


def test_maxwell_contour_force_2d_closed_uniform_field_cancels():
    contour = [(-1.0, -0.5), (1.0, -0.5), (1.0, 0.5), (-1.0, 0.5)]
    summary = maxwell_contour_force_2d(contour, (1.0, 0.0), orientation="ccw")

    assert summary["n_segments"] == 4
    assert summary["polygon_signed_area_m2"] == pytest.approx(2.0)
    assert summary["total_force_per_depth_N_per_m"] == pytest.approx([0.0, 0.0], abs=1.0e-9)
    assert summary["total_force_magnitude_per_depth_N_per_m"] == pytest.approx(0.0, abs=1.0e-9)
    assert summary["sum_abs_normal_force_per_depth_N_per_m"] > 0.0


def test_maxwell_contour_segment_balance_summary_reports_cancellation():
    pressure = air_gap_maxwell_pressure(1.0)
    contour = [(-1.0, -0.5), (1.0, -0.5), (1.0, 0.5), (-1.0, 0.5)]
    summary = maxwell_contour_segment_balance_summary_2d(
        contour,
        (1.0, 0.0),
        expected_force_per_depth_N_per_m=(0.0, 0.0),
    )

    assert summary["status"] == "ok"
    assert summary["reference_pass"] is True
    assert summary["orientation_consistent"] is True
    assert summary["total_force_per_depth_N_per_m"] == pytest.approx([0.0, 0.0], abs=1.0e-9)
    assert summary["sum_abs_normal_force_per_depth_N_per_m"] == pytest.approx(6.0 * pressure)
    assert summary["sum_abs_tangential_force_per_depth_N_per_m"] == pytest.approx(0.0)
    assert summary["cancellation_ratio"] == pytest.approx(0.0, abs=1.0e-14)
    assert summary["dominant_segment_index"] == 1
    assert [row["dominant_contribution"] for row in summary["segment_rows"]] == ["normal"] * 4


def test_parallel_wire_force_result_package_gate_tracks_ampere_sign_and_units():
    pair = parallel_wire_lorentz_force_summary(10.0, 5.0, 0.02)
    virtual = parallel_wire_virtual_work_force_summary(10.0, 5.0, 0.02)
    row = {
        "model_id": "slot188_parallel_wire_pair_v1",
        "operating_point_id": "I1_10A_I2_5A_d_20mm",
        "artifact_id": "femm_slot228_parallel_wire_recovered_table_A",
        "result_set_id": "rs_I1_10A_I2_5A_d_20mm_v1",
        "parameter_set_artifact_id": "femm_slot392_parallel_wire_parameter_set_v1.json",
        "parameter_set_digest": "sha256:slot392-parallel-wire-parameter-set-v1",
        "parameter_set_path": "artifacts/femm/slot392_parallel_wire_parameter_set.json",
        "model_input_artifact_id": "femm_slot378_parallel_wire_model_v1.fem",
        "model_input_digest": "sha256:slot378-parallel-wire-model-v1",
        "model_input_path": "artifacts/femm/femm_slot378_parallel_wire_model_v1.fem",
        "solution_artifact_id": "femm_slot252_parallel_wire_pair_v1.ans",
        "block_label_artifact_id": "femm_slot276_parallel_wire_block_labels_v1.json",
        "solution_loaded": True,
        "source_tool": "FEMM",
        "source_function": "block_integral_force_normalized_to_wire2_radial_axis",
        "source_group_id": 11,
        "target_group_id": 12,
        "source_center_xy_m": [0.0, 0.0],
        "target_center_xy_m": [0.02, 0.0],
        "source_region": "wire_1_positive_current",
        "target_region": "wire_2_positive_current",
        "source_material": "Copper",
        "target_material": "Copper",
        "selection_function": "mo_groupselectblock(12); mo_blockintegral(18); mo_blockintegral(19)",
        "force_component_frame": "global_cartesian_xy",
        "radial_projection_axis": "wire1_to_wire2_separation_axis_positive_away_from_wire1",
        "force_sign_convention": "positive_radial_force_points_away_from_wire1_attraction_is_negative",
        "force_extraction_method": "weighted_stress_block_integral_xy",
        "block_integral_types": [18, 19],
        "force_observable_id": "femm_slot300_weighted_stress_block_force_xy_v1",
        "force_observable_family": "femm_weighted_stress_block_force_xy",
        "force_convention_schema_id": "femm_parallel_wire_force_convention_v1",
        "force_component_basis_schema_id": "femm_global_xy_radial_projection_basis_v1",
        "postprocess_row_convention_schema_id": "femm_weighted_stress_force_row_convention_v1",
        "objective_observable_id": "femm_slot392_force_objective_radial_N_per_m_v1",
        "objective_observable_family": "force_minimize_abs_radial_N_per_m",
        "current1_A": 10.0,
        "current2_A": 5.0,
        "current_source_artifact_id": "femm_slot332_parallel_wire_current_definition_v1.json",
        "current_definition_method": "femm_circuit_current_snapshot",
        "separation_m": 0.02,
        "force_units": "N/m",
        "force_unit_basis": "per_length",
        "force_unit_basis_schema_id": "femm_planar_force_per_length_depth_basis_v1",
        "problem_depth_m": 1.0,
        "problem_type": "planar",
        "length_unit": "meters",
        "frequency_hz": 0.0,
        "solver_precision": 1.0e-8,
        "min_angle_deg": 30.0,
        "radial_force_on_wire2_N_per_m": -pair["force_magnitude_per_length_N_per_m"],
        "force_on_wire2_N_per_m": pair["force_on_wire2_N_per_m"],
        "interaction": "attraction",
    }

    gate = parallel_wire_force_result_package_gate(
        row,
        expected_model_id="slot188_parallel_wire_pair_v1",
        expected_operating_point_id="I1_10A_I2_5A_d_20mm",
        expected_artifact_id="femm_slot228_parallel_wire_recovered_table_A",
        expected_result_set_id="rs_I1_10A_I2_5A_d_20mm_v1",
        expected_parameter_set_artifact_id="femm_slot392_parallel_wire_parameter_set_v1.json",
        expected_parameter_set_digest="sha256:slot392-parallel-wire-parameter-set-v1",
        expected_parameter_set_path="artifacts/femm/slot392_parallel_wire_parameter_set.json",
        expected_model_input_artifact_id="femm_slot378_parallel_wire_model_v1.fem",
        expected_model_input_digest="sha256:slot378-parallel-wire-model-v1",
        expected_model_input_path="artifacts/femm/femm_slot378_parallel_wire_model_v1.fem",
        expected_solution_artifact_id="femm_slot252_parallel_wire_pair_v1.ans",
        expected_block_label_artifact_id="femm_slot276_parallel_wire_block_labels_v1.json",
        expected_source_tool="FEMM",
        expected_source_group_id=11,
        expected_target_group_id=12,
        expected_source_center_xy_m=(0.0, 0.0),
        expected_target_center_xy_m=(0.02, 0.0),
        expected_source_region="wire_1_positive_current",
        expected_target_region="wire_2_positive_current",
        expected_source_material="Copper",
        expected_target_material="Copper",
        expected_force_observable_id="femm_slot300_weighted_stress_block_force_xy_v1",
        expected_force_observable_family="femm_weighted_stress_block_force_xy",
        expected_force_convention_schema_id="femm_parallel_wire_force_convention_v1",
        expected_force_component_basis_schema_id="femm_global_xy_radial_projection_basis_v1",
        expected_force_unit_basis_schema_id="femm_planar_force_per_length_depth_basis_v1",
        expected_postprocess_row_convention_schema_id="femm_weighted_stress_force_row_convention_v1",
        expected_objective_observable_id="femm_slot392_force_objective_radial_N_per_m_v1",
        expected_objective_observable_family="force_minimize_abs_radial_N_per_m",
        expected_force_component_frame="global_cartesian_xy",
        expected_radial_projection_axis="wire1_to_wire2_separation_axis_positive_away_from_wire1",
        expected_force_sign_convention="positive_radial_force_points_away_from_wire1_attraction_is_negative",
        expected_force_extraction_method="weighted_stress_block_integral_xy",
        expected_block_integral_types=(18, 19),
        expected_current_source_artifact_id="femm_slot332_parallel_wire_current_definition_v1.json",
        expected_current_definition_method="femm_circuit_current_snapshot",
        expected_problem_type="planar",
        expected_length_unit="meters",
        expected_frequency_hz=0.0,
        expected_solver_precision=1.0e-8,
        max_solver_precision=1.0e-7,
        expected_min_angle_deg=30.0,
        require_solution_loaded=True,
        require_parameter_set_artifact=True,
        require_force_convention_schema=True,
        require_force_component_basis_schema=True,
        require_force_unit_basis_schema=True,
        require_postprocess_row_convention_schema=True,
    )

    assert pair["interaction"] == "attraction"
    assert pair["force_on_wire2_N_per_m"] == pytest.approx([-0.0005, 0.0])
    assert virtual["force_rel_error"] < 1.0e-8
    assert gate["policy"] == "parallel_wire_force_result_package_gate"
    assert gate["status"] == "ok"
    assert gate["artifact_id"] == "femm_slot228_parallel_wire_recovered_table_A"
    assert gate["result_set_id"] == "rs_I1_10A_I2_5A_d_20mm_v1"
    assert gate["checks"]["artifact_id_recorded"] is True
    assert gate["checks"]["result_set_id_recorded"] is True
    assert gate["checks"]["parameter_set_artifact_id_recorded"] is True
    assert gate["checks"]["parameter_set_digest_recorded"] is True
    assert gate["checks"]["parameter_set_path_recorded"] is True
    assert gate["checks"]["model_input_artifact_id_recorded"] is True
    assert gate["checks"]["model_input_digest_recorded"] is True
    assert gate["checks"]["model_input_path_recorded"] is True
    assert gate["checks"]["solution_artifact_id_recorded"] is True
    assert gate["checks"]["block_label_artifact_id_recorded"] is True
    assert gate["checks"]["expected_artifact_id_matches"] is True
    assert gate["checks"]["expected_result_set_id_matches"] is True
    assert gate["checks"]["expected_parameter_set_artifact_id_matches"] is True
    assert gate["checks"]["expected_parameter_set_digest_matches"] is True
    assert gate["checks"]["expected_parameter_set_path_matches"] is True
    assert gate["checks"]["expected_model_input_artifact_id_matches"] is True
    assert gate["checks"]["expected_model_input_digest_matches"] is True
    assert gate["checks"]["expected_model_input_path_matches"] is True
    assert gate["checks"]["expected_solution_artifact_id_matches"] is True
    assert gate["checks"]["expected_block_label_artifact_id_matches"] is True
    assert gate["checks"]["force_observable_id_recorded"] is True
    assert gate["checks"]["expected_force_observable_id_matches"] is True
    assert gate["checks"]["force_observable_family_recorded"] is True
    assert gate["checks"]["expected_force_observable_family_matches"] is True
    assert gate["checks"]["force_convention_schema_id_recorded"] is True
    assert gate["checks"]["expected_force_convention_schema_id_matches"] is True
    assert gate["checks"]["force_component_basis_schema_id_recorded"] is True
    assert gate["checks"]["expected_force_component_basis_schema_id_matches"] is True
    assert gate["checks"]["force_unit_basis_schema_id_recorded"] is True
    assert gate["checks"]["expected_force_unit_basis_schema_id_matches"] is True
    assert gate["checks"]["postprocess_row_convention_schema_id_recorded"] is True
    assert gate["checks"]["expected_postprocess_row_convention_schema_id_matches"] is True
    assert gate["checks"]["objective_observable_id_recorded"] is True
    assert gate["checks"]["expected_objective_observable_id_matches"] is True
    assert gate["checks"]["objective_observable_family_recorded"] is True
    assert gate["checks"]["expected_objective_observable_family_matches"] is True
    assert gate["parameter_set_artifact_id"] == "femm_slot392_parallel_wire_parameter_set_v1.json"
    assert gate["parameter_set_digest"] == "sha256:slot392-parallel-wire-parameter-set-v1"
    assert gate["parameter_set_path"].endswith("slot392_parallel_wire_parameter_set.json")
    assert gate["expected_parameter_set_artifact_id"] == "femm_slot392_parallel_wire_parameter_set_v1.json"
    assert gate["expected_parameter_set_digest"] == "sha256:slot392-parallel-wire-parameter-set-v1"
    assert gate["expected_parameter_set_path"].endswith("slot392_parallel_wire_parameter_set.json")
    assert gate["objective_observable_id"] == "femm_slot392_force_objective_radial_N_per_m_v1"
    assert gate["objective_observable_family"] == "force_minimize_abs_radial_N_per_m"
    assert gate["expected_objective_observable_id"] == "femm_slot392_force_objective_radial_N_per_m_v1"
    assert gate["expected_objective_observable_family"] == "force_minimize_abs_radial_N_per_m"
    assert gate["parameter_set_artifact_required"] is True
    assert gate["force_observable_id"] == "femm_slot300_weighted_stress_block_force_xy_v1"
    assert gate["force_observable_family"] == "femm_weighted_stress_block_force_xy"
    assert gate["force_convention_schema_id"] == "femm_parallel_wire_force_convention_v1"
    assert gate["expected_force_convention_schema_id"] == "femm_parallel_wire_force_convention_v1"
    assert gate["force_convention_schema_required"] is True
    assert gate["force_component_basis_schema_id"] == "femm_global_xy_radial_projection_basis_v1"
    assert gate["expected_force_component_basis_schema_id"] == "femm_global_xy_radial_projection_basis_v1"
    assert gate["force_component_basis_schema_required"] is True
    assert gate["force_unit_basis_schema_id"] == "femm_planar_force_per_length_depth_basis_v1"
    assert gate["expected_force_unit_basis_schema_id"] == "femm_planar_force_per_length_depth_basis_v1"
    assert gate["force_unit_basis_schema_required"] is True
    assert gate["postprocess_row_convention_schema_id"] == "femm_weighted_stress_force_row_convention_v1"
    assert gate["expected_postprocess_row_convention_schema_id"] == "femm_weighted_stress_force_row_convention_v1"
    assert gate["expected_force_component_frame"] == "global_cartesian_xy"
    assert gate["expected_radial_projection_axis"] == "wire1_to_wire2_separation_axis_positive_away_from_wire1"
    assert gate["expected_force_sign_convention"] == "positive_radial_force_points_away_from_wire1_attraction_is_negative"
    assert gate["checks"]["solution_loaded_recorded"] is True
    assert gate["checks"]["solution_loaded_before_postprocess"] is True
    assert gate["solution_artifact_id"] == "femm_slot252_parallel_wire_pair_v1.ans"
    assert gate["block_label_artifact_id"] == "femm_slot276_parallel_wire_block_labels_v1.json"
    assert gate["solution_loaded"] is True
    assert gate["solution_loaded_required"] is True
    assert gate["checks"]["expected_source_tool_matches"] is True
    assert gate["checks"]["force_units_are_per_length"] is True
    assert gate["force_unit_basis"] == "per_length"
    assert gate["problem_depth_m"] == pytest.approx(1.0)
    assert gate["problem_type"] == "planar"
    assert gate["length_unit"] == "meters"
    assert gate["frequency_hz"] == pytest.approx(0.0)
    assert gate["solver_precision"] == pytest.approx(1.0e-8)
    assert gate["min_angle_deg"] == pytest.approx(30.0)
    assert gate["expected_problem_type"] == "planar"
    assert gate["expected_length_unit"] == "meters"
    assert gate["expected_frequency_hz"] == pytest.approx(0.0)
    assert gate["expected_solver_precision"] == pytest.approx(1.0e-8)
    assert gate["max_solver_precision"] == pytest.approx(1.0e-7)
    assert gate["expected_min_angle_deg"] == pytest.approx(30.0)
    assert gate["checks"]["force_unit_basis_is_per_length"] is True
    assert gate["checks"]["depth_integrated_force_not_used_for_per_length_gate"] is True
    assert gate["checks"]["femm_planar_depth_recorded"] is True
    assert gate["checks"]["femm_planar_depth_positive"] is True
    assert gate["checks"]["problem_type_recorded"] is True
    assert gate["checks"]["expected_problem_type_matches"] is True
    assert gate["checks"]["length_unit_recorded"] is True
    assert gate["checks"]["expected_length_unit_matches"] is True
    assert gate["checks"]["frequency_hz_recorded"] is True
    assert gate["checks"]["expected_frequency_hz_matches"] is True
    assert gate["checks"]["solver_precision_recorded"] is True
    assert gate["checks"]["expected_solver_precision_matches"] is True
    assert gate["checks"]["solver_precision_within_max"] is True
    assert gate["checks"]["min_angle_deg_recorded"] is True
    assert gate["checks"]["expected_min_angle_deg_matches"] is True
    assert gate["checks"]["min_angle_deg_positive"] is True
    assert gate["source_group_id"] == "11"
    assert gate["target_group_id"] == "12"
    assert gate["source_center_xy_m"] == [0.0, 0.0]
    assert gate["target_center_xy_m"] == [0.02, 0.0]
    assert gate["center_distance_m"] == pytest.approx(0.02)
    assert gate["center_distance_error_m"] == pytest.approx(0.0)
    assert gate["checks"]["expected_source_group_id_matches"] is True
    assert gate["checks"]["expected_target_group_id_matches"] is True
    assert gate["checks"]["source_center_xy_recorded_when_expected"] is True
    assert gate["checks"]["target_center_xy_recorded_when_expected"] is True
    assert gate["checks"]["expected_source_center_xy_matches"] is True
    assert gate["checks"]["expected_target_center_xy_matches"] is True
    assert gate["checks"]["wire_center_separation_matches_separation_m"] is True
    assert gate["checks"]["expected_source_region_matches"] is True
    assert gate["checks"]["expected_target_region_matches"] is True
    assert gate["checks"]["expected_source_material_matches"] is True
    assert gate["checks"]["expected_target_material_matches"] is True
    assert gate["checks"]["source_target_groups_distinct"] is True
    assert gate["checks"]["selection_mentions_target_group"] is True
    assert gate["checks"]["force_component_frame_recorded"] is True
    assert gate["checks"]["radial_projection_axis_recorded"] is True
    assert gate["checks"]["radial_projection_axis_names_wire_pair"] is True
    assert gate["checks"]["force_component_frame_recorded_when_expected"] is True
    assert gate["checks"]["expected_force_component_frame_matches"] is True
    assert gate["checks"]["radial_projection_axis_recorded_when_expected"] is True
    assert gate["checks"]["expected_radial_projection_axis_matches"] is True
    assert gate["checks"]["force_sign_convention_recorded_when_expected"] is True
    assert gate["checks"]["expected_force_sign_convention_matches"] is True
    assert gate["checks"]["force_extraction_method_recorded_when_expected"] is True
    assert gate["checks"]["expected_force_extraction_method_matches"] is True
    assert gate["checks"]["weighted_stress_extraction_uses_force_xy_integrals"] is True
    assert gate["force_component_frame"] == "global_cartesian_xy"
    assert gate["radial_projection_axis"] == "wire1_to_wire2_separation_axis_positive_away_from_wire1"
    assert gate["force_sign_convention"] == "positive_radial_force_points_away_from_wire1_attraction_is_negative"
    assert gate["force_extraction_method"] == "weighted_stress_block_integral_xy"
    assert gate["expected_force_extraction_method"] == "weighted_stress_block_integral_xy"
    assert gate["current_source_artifact_id"] == "femm_slot332_parallel_wire_current_definition_v1.json"
    assert gate["current_definition_method"] == "femm_circuit_current_snapshot"
    assert gate["expected_current_source_artifact_id"] == "femm_slot332_parallel_wire_current_definition_v1.json"
    assert gate["expected_current_definition_method"] == "femm_circuit_current_snapshot"
    assert gate["checks"]["current_source_artifact_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_current_source_artifact_id_matches"] is True
    assert gate["checks"]["current_definition_method_recorded_when_expected"] is True
    assert gate["checks"]["expected_current_definition_method_matches"] is True
    assert gate["checks"]["block_integral_types_recorded"] is True
    assert gate["checks"]["block_integral_types_match_expected"] is True
    assert gate["checks"]["block_integral_types_are_force_xy"] is True
    assert gate["checks"]["block_integral_types_exclude_torque"] is True
    assert gate["checks"]["vector_force_matches_ampere"] is True
    assert gate["checks"]["radial_force_matches_ampere"] is True

    stale_component_basis = parallel_wire_force_result_package_gate(
        {
            **row,
            "force_component_basis_schema_id": "femm_local_rt_force_basis_v0",
        },
        expected_force_component_basis_schema_id="femm_global_xy_radial_projection_basis_v1",
        require_force_component_basis_schema=True,
    )
    assert stale_component_basis["status"] == "needs_attention"
    assert stale_component_basis["checks"]["expected_force_component_basis_schema_id_matches"] is False
    assert stale_component_basis["checks"]["vector_force_matches_ampere"] is True

    missing_component_basis = parallel_wire_force_result_package_gate(
        {
            key: value
            for key, value in row.items()
            if key != "force_component_basis_schema_id"
        },
        expected_force_component_basis_schema_id="femm_global_xy_radial_projection_basis_v1",
        require_force_component_basis_schema=True,
    )
    assert missing_component_basis["status"] == "needs_attention"
    assert missing_component_basis["checks"]["force_component_basis_schema_id_recorded"] is False
    assert missing_component_basis["checks"]["expected_force_component_basis_schema_id_matches"] is False

    stale_force_unit_basis = parallel_wire_force_result_package_gate(
        {
            **row,
            "force_unit_basis_schema_id": "femm_depth_integrated_force_unit_basis_v0",
        },
        expected_force_unit_basis_schema_id="femm_planar_force_per_length_depth_basis_v1",
        require_force_unit_basis_schema=True,
    )
    assert stale_force_unit_basis["status"] == "needs_attention"
    assert stale_force_unit_basis["checks"]["expected_force_unit_basis_schema_id_matches"] is False
    assert stale_force_unit_basis["checks"]["force_unit_basis_is_per_length"] is True
    assert stale_force_unit_basis["checks"]["vector_force_matches_ampere"] is True

    missing_force_unit_basis_row = dict(row)
    missing_force_unit_basis_row.pop("force_unit_basis_schema_id")
    missing_force_unit_basis = parallel_wire_force_result_package_gate(
        missing_force_unit_basis_row,
        expected_force_unit_basis_schema_id="femm_planar_force_per_length_depth_basis_v1",
        require_force_unit_basis_schema=True,
    )
    assert missing_force_unit_basis["status"] == "needs_attention"
    assert missing_force_unit_basis["checks"]["force_unit_basis_schema_id_recorded"] is False
    assert missing_force_unit_basis["checks"]["expected_force_unit_basis_schema_id_matches"] is False

    clear_selection = parallel_wire_force_result_package_gate(
        {
            **row,
            "artifact_id": "femm_slot244_parallel_wire_clear_selection_A",
            "result_set_id": "rs_I1_10A_I2_5A_d_20mm_clear_selection_v1",
            "selection_function": "mo_clearblock(); mo_groupselectblock(12); mo_blockintegral(18); mo_blockintegral(19)",
        },
        expected_artifact_id="femm_slot244_parallel_wire_clear_selection_A",
        expected_result_set_id="rs_I1_10A_I2_5A_d_20mm_clear_selection_v1",
        expected_target_group_id=12,
        require_selection_clear=True,
    )
    assert clear_selection["status"] == "ok"
    assert clear_selection["selection_clear_required"] is True
    assert clear_selection["checks"]["selection_clear_before_groupselect"] is True

    traced_postprocess = parallel_wire_force_result_package_gate(
        {
            **row,
            "artifact_id": "femm_slot284_parallel_wire_postprocess_trace_A",
            "result_set_id": "rs_I1_10A_I2_5A_d_20mm_trace_v1",
            "selection_function": "mo_clearblock(); mo_groupselectblock(12); mo_blockintegral(18); mo_blockintegral(19)",
            "postprocess_trace_id": "femm_slot284_postprocess_trace_v1",
            "postprocess_command_digest": "sha256:slot284_clear_select_force_xy",
            "postprocess_output_artifact_id": "femm_slot292_parallel_wire_force_table_v1.json",
            "postprocess_output_digest": "sha256:slot292_force_table_json",
            "postprocess_output_path": "artifacts/femm/slot292_force_table.json",
            "postprocess_output_schema_id": "femm_force_table_xy_v1",
            "postprocess_output_columns": [
                "target_group_id",
                "Fx_N_per_m",
                "Fy_N_per_m",
                "radial_force_on_wire2_N_per_m",
            ],
            "postprocess_output_units": {
                "target_group_id": "1",
                "Fx_N_per_m": "N/m",
                "Fy_N_per_m": "N/m",
                "radial_force_on_wire2_N_per_m": "N/m",
            },
            "postprocess_script_artifact_id": "femm_slot385_parallel_wire_postprocess_script_v1.py",
            "postprocess_script_digest": "sha256:slot385_parallel_wire_postprocess_script_v1",
            "postprocess_script_path": "artifacts/femm/femm_slot385_parallel_wire_postprocess.py",
            "postprocess_commands": [
                "mo_clearblock()",
                "mo_groupselectblock(12)",
                "mo_blockintegral(18)",
                "mo_blockintegral(19)",
            ],
        },
        expected_artifact_id="femm_slot284_parallel_wire_postprocess_trace_A",
        expected_result_set_id="rs_I1_10A_I2_5A_d_20mm_trace_v1",
        expected_target_group_id=12,
        expected_postprocess_trace_id="femm_slot284_postprocess_trace_v1",
        expected_postprocess_command_digest="sha256:slot284_clear_select_force_xy",
        expected_postprocess_output_artifact_id="femm_slot292_parallel_wire_force_table_v1.json",
        expected_postprocess_output_digest="sha256:slot292_force_table_json",
        expected_postprocess_output_schema_id="femm_force_table_xy_v1",
        expected_postprocess_output_columns=[
            "target_group_id",
            "Fx_N_per_m",
            "Fy_N_per_m",
            "radial_force_on_wire2_N_per_m",
        ],
        expected_postprocess_output_units={
            "target_group_id": "1",
            "Fx_N_per_m": "N/m",
            "Fy_N_per_m": "N/m",
            "radial_force_on_wire2_N_per_m": "N/m",
        },
        expected_postprocess_script_artifact_id="femm_slot385_parallel_wire_postprocess_script_v1.py",
        expected_postprocess_script_digest="sha256:slot385_parallel_wire_postprocess_script_v1",
        expected_postprocess_script_path="artifacts/femm/femm_slot385_parallel_wire_postprocess.py",
        expected_force_observable_id="femm_slot300_weighted_stress_block_force_xy_v1",
        expected_force_observable_family="femm_weighted_stress_block_force_xy",
        require_selection_clear=True,
        require_postprocess_command_trace=True,
        require_postprocess_output_artifact=True,
        require_postprocess_output_schema=True,
        require_postprocess_script_artifact=True,
    )
    assert traced_postprocess["status"] == "ok"
    assert traced_postprocess["postprocess_command_trace_required"] is True
    assert traced_postprocess["postprocess_output_artifact_required"] is True
    assert traced_postprocess["postprocess_script_artifact_required"] is True
    assert traced_postprocess["checks"]["postprocess_trace_id_recorded"] is True
    assert traced_postprocess["checks"]["expected_postprocess_trace_id_matches"] is True
    assert traced_postprocess["checks"]["postprocess_command_digest_recorded"] is True
    assert traced_postprocess["checks"]["expected_postprocess_command_digest_matches"] is True
    assert traced_postprocess["checks"]["postprocess_output_artifact_id_recorded"] is True
    assert traced_postprocess["checks"]["expected_postprocess_output_artifact_id_matches"] is True
    assert traced_postprocess["checks"]["postprocess_output_digest_recorded"] is True
    assert traced_postprocess["checks"]["expected_postprocess_output_digest_matches"] is True
    assert traced_postprocess["checks"]["postprocess_output_schema_id_recorded"] is True
    assert traced_postprocess["checks"]["expected_postprocess_output_schema_id_matches"] is True
    assert traced_postprocess["checks"]["postprocess_output_columns_recorded"] is True
    assert traced_postprocess["checks"]["expected_postprocess_output_columns_match"] is True
    assert traced_postprocess["checks"]["postprocess_output_units_recorded"] is True
    assert traced_postprocess["checks"]["expected_postprocess_output_units_match"] is True
    assert traced_postprocess["postprocess_output_schema_required"] is True
    assert traced_postprocess["postprocess_output_schema_id"] == "femm_force_table_xy_v1"
    assert traced_postprocess["checks"]["postprocess_script_artifact_id_recorded"] is True
    assert traced_postprocess["checks"]["postprocess_script_digest_recorded"] is True
    assert traced_postprocess["checks"]["postprocess_script_path_recorded"] is True
    assert traced_postprocess["checks"]["expected_postprocess_script_artifact_id_matches"] is True
    assert traced_postprocess["checks"]["expected_postprocess_script_digest_matches"] is True
    assert traced_postprocess["checks"]["expected_postprocess_script_path_matches"] is True
    assert traced_postprocess["checks"]["expected_force_observable_id_matches"] is True
    assert traced_postprocess["checks"]["expected_force_observable_family_matches"] is True
    assert traced_postprocess["postprocess_output_path"].endswith("slot292_force_table.json")
    assert traced_postprocess["postprocess_script_path"].endswith("femm_slot385_parallel_wire_postprocess.py")
    assert traced_postprocess["checks"]["postprocess_commands_recorded"] is True
    assert traced_postprocess["checks"]["postprocess_commands_clear_select_force_xy"] is True
    assert traced_postprocess["checks"]["postprocess_commands_mention_target_group"] is True

    execution_packaged = parallel_wire_force_result_package_gate(
        {
            **row,
            "artifact_id": "femm_slot371_parallel_wire_execution_package_A",
            "result_set_id": "rs_I1_10A_I2_5A_d_20mm_execution_v1",
            "created_at_utc": "2026-07-01T10:37:20Z",
            "run_timestamp_utc": "2026-07-01T10:37:00Z",
            "solver_version": "FEMM 4.2",
            "radia_mcp_version": "1.4.3",
            "run_duration_s": 2.0,
            "timing_breakdown_s": {
                "mesh_s": 0.8,
                "solve_s": 0.7,
                "postprocess_s": 0.3,
                "write_json_s": 0.1,
            },
        },
        expected_artifact_id="femm_slot371_parallel_wire_execution_package_A",
        expected_result_set_id="rs_I1_10A_I2_5A_d_20mm_execution_v1",
        expected_created_at_utc="2026-07-01T10:37:20Z",
        expected_run_timestamp_utc="2026-07-01T10:37:00Z",
        expected_solver_version="FEMM 4.2",
        expected_radia_mcp_version="1.4.3",
        max_created_run_skew_s=60,
        require_execution_metadata=True,
        require_timing_breakdown=True,
        min_timing_sections=4,
    )
    assert execution_packaged["status"] == "ok"
    assert execution_packaged["execution_metadata_required"] is True
    assert execution_packaged["timing_breakdown_required"] is True
    assert execution_packaged["checks"]["created_at_utc_recorded"] is True
    assert execution_packaged["checks"]["created_at_utc_parseable"] is True
    assert execution_packaged["checks"]["expected_created_at_utc_matches"] is True
    assert execution_packaged["checks"]["run_timestamp_utc_recorded"] is True
    assert execution_packaged["checks"]["run_timestamp_utc_parseable"] is True
    assert execution_packaged["checks"]["expected_run_timestamp_utc_matches"] is True
    assert execution_packaged["checks"]["created_run_timestamp_skew_within_limit"] is True
    assert execution_packaged["created_run_timestamp_skew_s"] == pytest.approx(20.0)
    assert execution_packaged["checks"]["solver_version_recorded"] is True
    assert execution_packaged["checks"]["expected_solver_version_matches"] is True
    assert execution_packaged["checks"]["radia_mcp_version_recorded"] is True
    assert execution_packaged["checks"]["expected_radia_mcp_version_matches"] is True
    assert execution_packaged["checks"]["run_duration_s_recorded"] is True
    assert execution_packaged["checks"]["run_duration_s_positive"] is True
    assert execution_packaged["checks"]["timing_breakdown_recorded"] is True
    assert execution_packaged["checks"]["timing_breakdown_has_required_sections"] is True
    assert execution_packaged["checks"]["timing_breakdown_top_sections_descending"] is True
    assert execution_packaged["timing_total_s"] == pytest.approx(1.9)

    model_input_packaged = parallel_wire_force_result_package_gate(
        {
            **row,
            "artifact_id": "femm_slot378_parallel_wire_model_input_package_A",
            "result_set_id": "rs_I1_10A_I2_5A_d_20mm_model_input_v1",
        },
        expected_artifact_id="femm_slot378_parallel_wire_model_input_package_A",
        expected_result_set_id="rs_I1_10A_I2_5A_d_20mm_model_input_v1",
        expected_model_input_artifact_id="femm_slot378_parallel_wire_model_v1.fem",
        expected_model_input_digest="sha256:slot378-parallel-wire-model-v1",
        expected_model_input_path="artifacts/femm/femm_slot378_parallel_wire_model_v1.fem",
        require_model_input_artifact=True,
    )
    assert model_input_packaged["status"] == "ok"
    assert model_input_packaged["model_input_artifact_required"] is True
    assert model_input_packaged["model_input_artifact_id"] == "femm_slot378_parallel_wire_model_v1.fem"
    assert model_input_packaged["checks"]["model_input_artifact_id_recorded"] is True
    assert model_input_packaged["checks"]["model_input_digest_recorded"] is True
    assert model_input_packaged["checks"]["model_input_path_recorded"] is True
    assert model_input_packaged["checks"]["expected_model_input_artifact_id_matches"] is True
    assert model_input_packaged["checks"]["expected_model_input_digest_matches"] is True
    assert model_input_packaged["checks"]["expected_model_input_path_matches"] is True

    stale_model_input_digest = parallel_wire_force_result_package_gate(
        {
            **row,
            "model_input_digest": "sha256:slot188-old-parallel-wire-model",
        },
        expected_model_input_artifact_id="femm_slot378_parallel_wire_model_v1.fem",
        expected_model_input_digest="sha256:slot378-parallel-wire-model-v1",
        expected_model_input_path="artifacts/femm/femm_slot378_parallel_wire_model_v1.fem",
        require_model_input_artifact=True,
    )
    assert stale_model_input_digest["status"] == "needs_attention"
    assert stale_model_input_digest["checks"]["expected_model_input_artifact_id_matches"] is True
    assert stale_model_input_digest["checks"]["expected_model_input_digest_matches"] is False
    assert stale_model_input_digest["checks"]["vector_force_matches_ampere"] is True

    missing_model_input_path_row = dict(row)
    missing_model_input_path_row.pop("model_input_path")
    missing_model_input_path = parallel_wire_force_result_package_gate(
        missing_model_input_path_row,
        require_model_input_artifact=True,
    )
    assert missing_model_input_path["status"] == "needs_attention"
    assert missing_model_input_path["checks"]["model_input_artifact_id_recorded"] is True
    assert missing_model_input_path["checks"]["model_input_digest_recorded"] is True
    assert missing_model_input_path["checks"]["model_input_path_recorded"] is False

    stale_execution_version = parallel_wire_force_result_package_gate(
        {
            **row,
            "run_timestamp_utc": "2026-07-01T10:37:00Z",
            "solver_version": "FEMM 4.1",
            "radia_mcp_version": "1.4.3",
            "run_duration_s": 2.0,
            "timing_breakdown_s": {
                "mesh_s": 0.8,
                "solve_s": 0.7,
                "postprocess_s": 0.3,
                "write_json_s": 0.1,
            },
        },
        expected_solver_version="FEMM 4.2",
        expected_radia_mcp_version="1.4.3",
        require_execution_metadata=True,
        require_timing_breakdown=True,
    )
    assert stale_execution_version["status"] == "needs_attention"
    assert stale_execution_version["checks"]["expected_solver_version_matches"] is False
    assert stale_execution_version["checks"]["vector_force_matches_ampere"] is True

    stale_created_run_skew = parallel_wire_force_result_package_gate(
        {
            **row,
            "created_at_utc": "2026-07-01T12:37:20Z",
            "run_timestamp_utc": "2026-07-01T10:37:00Z",
            "solver_version": "FEMM 4.2",
            "radia_mcp_version": "1.4.3",
            "run_duration_s": 2.0,
            "timing_breakdown_s": {
                "mesh_s": 0.8,
                "solve_s": 0.7,
                "postprocess_s": 0.3,
                "write_json_s": 0.1,
            },
        },
        expected_created_at_utc="2026-07-01T12:37:20Z",
        expected_run_timestamp_utc="2026-07-01T10:37:00Z",
        max_created_run_skew_s=60,
        require_execution_metadata=True,
        require_timing_breakdown=True,
    )
    assert stale_created_run_skew["status"] == "needs_attention"
    assert stale_created_run_skew["checks"]["created_at_utc_parseable"] is True
    assert stale_created_run_skew["checks"]["run_timestamp_utc_parseable"] is True
    assert stale_created_run_skew["checks"]["created_run_timestamp_skew_within_limit"] is False
    assert stale_created_run_skew["checks"]["vector_force_matches_ampere"] is True

    thin_timing_breakdown = parallel_wire_force_result_package_gate(
        {
            **row,
            "run_timestamp_utc": "2026-07-01T10:37:00Z",
            "solver_version": "FEMM 4.2",
            "radia_mcp_version": "1.4.3",
            "run_duration_s": 2.0,
            "timing_breakdown_s": {"solve_s": 1.5, "write_json_s": 0.1},
        },
        require_execution_metadata=True,
        require_timing_breakdown=True,
        min_timing_sections=4,
    )
    assert thin_timing_breakdown["status"] == "needs_attention"
    assert thin_timing_breakdown["checks"]["timing_breakdown_has_required_sections"] is False
    assert thin_timing_breakdown["checks"]["vector_force_matches_ampere"] is True

    impossible_timing_total = parallel_wire_force_result_package_gate(
        {
            **row,
            "run_timestamp_utc": "2026-07-01T10:37:00Z",
            "solver_version": "FEMM 4.2",
            "radia_mcp_version": "1.4.3",
            "run_duration_s": 1.0,
            "timing_breakdown_s": {
                "mesh_s": 0.8,
                "solve_s": 0.7,
                "postprocess_s": 0.3,
                "write_json_s": 0.1,
            },
        },
        require_execution_metadata=True,
        require_timing_breakdown=True,
    )
    assert impossible_timing_total["status"] == "needs_attention"
    assert impossible_timing_total["checks"]["timing_breakdown_total_within_run_duration"] is False
    assert impossible_timing_total["checks"]["vector_force_matches_ampere"] is True

    stale_trace = parallel_wire_force_result_package_gate(
        {
            **row,
            "postprocess_trace_id": "femm_slot244_postprocess_trace_old",
            "postprocess_command_digest": "sha256:slot284_clear_select_force_xy",
            "postprocess_commands": [
                "mo_clearblock()",
                "mo_groupselectblock(12)",
                "mo_blockintegral(18)",
                "mo_blockintegral(19)",
            ],
        },
        expected_postprocess_trace_id="femm_slot284_postprocess_trace_v1",
        expected_postprocess_command_digest="sha256:slot284_clear_select_force_xy",
        require_postprocess_command_trace=True,
    )
    assert stale_trace["status"] == "needs_attention"
    assert stale_trace["checks"]["expected_postprocess_trace_id_matches"] is False
    assert stale_trace["checks"]["vector_force_matches_ampere"] is True

    stale_postprocess_output = parallel_wire_force_result_package_gate(
        {
            **row,
            "postprocess_trace_id": "femm_slot284_postprocess_trace_v1",
            "postprocess_command_digest": "sha256:slot284_clear_select_force_xy",
            "postprocess_output_artifact_id": "femm_slot244_old_force_table.json",
            "postprocess_output_digest": "sha256:slot292_force_table_json",
            "postprocess_commands": [
                "mo_clearblock()",
                "mo_groupselectblock(12)",
                "mo_blockintegral(18)",
                "mo_blockintegral(19)",
            ],
        },
        expected_postprocess_trace_id="femm_slot284_postprocess_trace_v1",
        expected_postprocess_command_digest="sha256:slot284_clear_select_force_xy",
        expected_postprocess_output_artifact_id="femm_slot292_parallel_wire_force_table_v1.json",
        expected_postprocess_output_digest="sha256:slot292_force_table_json",
        require_postprocess_command_trace=True,
        require_postprocess_output_artifact=True,
    )
    assert stale_postprocess_output["status"] == "needs_attention"
    assert stale_postprocess_output["checks"]["expected_postprocess_output_artifact_id_matches"] is False
    assert stale_postprocess_output["checks"]["expected_postprocess_trace_id_matches"] is True
    assert stale_postprocess_output["checks"]["vector_force_matches_ampere"] is True

    stale_postprocess_output_schema = parallel_wire_force_result_package_gate(
        {
            **row,
            "postprocess_trace_id": "femm_slot284_postprocess_trace_v1",
            "postprocess_command_digest": "sha256:slot284_clear_select_force_xy",
            "postprocess_output_artifact_id": "femm_slot292_parallel_wire_force_table_v1.json",
            "postprocess_output_digest": "sha256:slot292_force_table_json",
            "postprocess_output_schema_id": "femm_force_table_scalar_v0",
            "postprocess_output_columns": [
                "target_group_id",
                "force_N",
            ],
            "postprocess_output_units": {
                "target_group_id": "1",
                "force_N": "N",
            },
            "postprocess_commands": [
                "mo_clearblock()",
                "mo_groupselectblock(12)",
                "mo_blockintegral(18)",
                "mo_blockintegral(19)",
            ],
        },
        expected_postprocess_trace_id="femm_slot284_postprocess_trace_v1",
        expected_postprocess_command_digest="sha256:slot284_clear_select_force_xy",
        expected_postprocess_output_artifact_id="femm_slot292_parallel_wire_force_table_v1.json",
        expected_postprocess_output_digest="sha256:slot292_force_table_json",
        expected_postprocess_output_schema_id="femm_force_table_xy_v1",
        expected_postprocess_output_columns=[
            "target_group_id",
            "Fx_N_per_m",
            "Fy_N_per_m",
            "radial_force_on_wire2_N_per_m",
        ],
        expected_postprocess_output_units={
            "target_group_id": "1",
            "Fx_N_per_m": "N/m",
            "Fy_N_per_m": "N/m",
            "radial_force_on_wire2_N_per_m": "N/m",
        },
        require_postprocess_command_trace=True,
        require_postprocess_output_artifact=True,
        require_postprocess_output_schema=True,
    )
    assert stale_postprocess_output_schema["status"] == "needs_attention"
    assert stale_postprocess_output_schema["checks"]["expected_postprocess_output_artifact_id_matches"] is True
    assert stale_postprocess_output_schema["checks"]["expected_postprocess_output_schema_id_matches"] is False
    assert stale_postprocess_output_schema["checks"]["expected_postprocess_output_columns_match"] is False
    assert stale_postprocess_output_schema["checks"]["expected_postprocess_output_units_match"] is False
    assert stale_postprocess_output_schema["checks"]["vector_force_matches_ampere"] is True

    stale_force_convention_schema = parallel_wire_force_result_package_gate(
        {
            **row,
            "force_convention_schema_id": "femm_force_value_only_convention_v0",
        },
        expected_force_convention_schema_id="femm_parallel_wire_force_convention_v1",
        require_force_convention_schema=True,
    )
    assert stale_force_convention_schema["status"] == "needs_attention"
    assert stale_force_convention_schema["checks"]["expected_force_convention_schema_id_matches"] is False
    assert stale_force_convention_schema["checks"]["force_convention_schema_id_recorded"] is True
    assert stale_force_convention_schema["checks"]["vector_force_matches_ampere"] is True

    missing_force_convention_schema_row = dict(row)
    missing_force_convention_schema_row.pop("force_convention_schema_id")
    missing_force_convention_schema = parallel_wire_force_result_package_gate(
        missing_force_convention_schema_row,
        expected_force_convention_schema_id="femm_parallel_wire_force_convention_v1",
        require_force_convention_schema=True,
    )
    assert missing_force_convention_schema["status"] == "needs_attention"
    assert missing_force_convention_schema["checks"]["force_convention_schema_id_recorded"] is False
    assert missing_force_convention_schema["checks"]["expected_force_convention_schema_id_matches"] is False
    assert missing_force_convention_schema["checks"]["vector_force_matches_ampere"] is True

    stale_postprocess_row_convention_schema = parallel_wire_force_result_package_gate(
        {
            **row,
            "postprocess_row_convention_schema_id": "femm_force_scalar_row_convention_v0",
        },
        expected_force_convention_schema_id="femm_parallel_wire_force_convention_v1",
        expected_postprocess_row_convention_schema_id="femm_weighted_stress_force_row_convention_v1",
        require_force_convention_schema=True,
        require_postprocess_row_convention_schema=True,
    )
    assert stale_postprocess_row_convention_schema["status"] == "needs_attention"
    assert stale_postprocess_row_convention_schema["checks"]["expected_force_convention_schema_id_matches"] is True
    assert (
        stale_postprocess_row_convention_schema["checks"][
            "expected_postprocess_row_convention_schema_id_matches"
        ]
        is False
    )
    assert stale_postprocess_row_convention_schema["checks"]["vector_force_matches_ampere"] is True

    missing_postprocess_row_convention_schema_row = dict(row)
    missing_postprocess_row_convention_schema_row.pop("postprocess_row_convention_schema_id")
    missing_postprocess_row_convention_schema = parallel_wire_force_result_package_gate(
        missing_postprocess_row_convention_schema_row,
        expected_postprocess_row_convention_schema_id="femm_weighted_stress_force_row_convention_v1",
        require_postprocess_row_convention_schema=True,
    )
    assert missing_postprocess_row_convention_schema["status"] == "needs_attention"
    assert (
        missing_postprocess_row_convention_schema["checks"][
            "postprocess_row_convention_schema_id_recorded"
        ]
        is False
    )
    assert (
        missing_postprocess_row_convention_schema["checks"][
            "expected_postprocess_row_convention_schema_id_matches"
        ]
        is False
    )
    assert missing_postprocess_row_convention_schema["checks"]["vector_force_matches_ampere"] is True

    stale_postprocess_script = parallel_wire_force_result_package_gate(
        {
            **row,
            "postprocess_trace_id": "femm_slot284_postprocess_trace_v1",
            "postprocess_command_digest": "sha256:slot284_clear_select_force_xy",
            "postprocess_output_artifact_id": "femm_slot292_parallel_wire_force_table_v1.json",
            "postprocess_output_digest": "sha256:slot292_force_table_json",
            "postprocess_script_artifact_id": "femm_slot385_parallel_wire_postprocess_script_v1.py",
            "postprocess_script_digest": "sha256:slot244_old_postprocess_script",
            "postprocess_script_path": "artifacts/femm/femm_slot385_parallel_wire_postprocess.py",
            "postprocess_commands": [
                "mo_clearblock()",
                "mo_groupselectblock(12)",
                "mo_blockintegral(18)",
                "mo_blockintegral(19)",
            ],
        },
        expected_postprocess_trace_id="femm_slot284_postprocess_trace_v1",
        expected_postprocess_command_digest="sha256:slot284_clear_select_force_xy",
        expected_postprocess_output_artifact_id="femm_slot292_parallel_wire_force_table_v1.json",
        expected_postprocess_output_digest="sha256:slot292_force_table_json",
        expected_postprocess_script_artifact_id="femm_slot385_parallel_wire_postprocess_script_v1.py",
        expected_postprocess_script_digest="sha256:slot385_parallel_wire_postprocess_script_v1",
        expected_postprocess_script_path="artifacts/femm/femm_slot385_parallel_wire_postprocess.py",
        require_postprocess_command_trace=True,
        require_postprocess_output_artifact=True,
        require_postprocess_script_artifact=True,
    )
    assert stale_postprocess_script["status"] == "needs_attention"
    assert stale_postprocess_script["checks"]["expected_postprocess_script_artifact_id_matches"] is True
    assert stale_postprocess_script["checks"]["expected_postprocess_script_digest_matches"] is False
    assert stale_postprocess_script["checks"]["expected_postprocess_output_artifact_id_matches"] is True
    assert stale_postprocess_script["checks"]["vector_force_matches_ampere"] is True

    missing_postprocess_script_path_row = {
        **row,
        "postprocess_trace_id": "femm_slot284_postprocess_trace_v1",
        "postprocess_command_digest": "sha256:slot284_clear_select_force_xy",
        "postprocess_script_artifact_id": "femm_slot385_parallel_wire_postprocess_script_v1.py",
        "postprocess_script_digest": "sha256:slot385_parallel_wire_postprocess_script_v1",
    }
    missing_postprocess_script_path = parallel_wire_force_result_package_gate(
        missing_postprocess_script_path_row,
        require_postprocess_script_artifact=True,
    )
    assert missing_postprocess_script_path["status"] == "needs_attention"
    assert missing_postprocess_script_path["checks"]["postprocess_script_artifact_id_recorded"] is True
    assert missing_postprocess_script_path["checks"]["postprocess_script_digest_recorded"] is True
    assert missing_postprocess_script_path["checks"]["postprocess_script_path_recorded"] is False

    stale_force_observable = parallel_wire_force_result_package_gate(
        {
            **row,
            "force_observable_id": "femm_slot212_torque_observable_old",
            "force_observable_family": "femm_weighted_stress_block_force_xy",
        },
        expected_force_observable_id="femm_slot300_weighted_stress_block_force_xy_v1",
        expected_force_observable_family="femm_weighted_stress_block_force_xy",
    )
    assert stale_force_observable["status"] == "needs_attention"
    assert stale_force_observable["checks"]["expected_force_observable_id_matches"] is False
    assert stale_force_observable["checks"]["expected_force_observable_family_matches"] is True
    assert stale_force_observable["checks"]["vector_force_matches_ampere"] is True

    wrong_force_observable_family = parallel_wire_force_result_package_gate(
        {
            **row,
            "force_observable_id": "femm_slot300_weighted_stress_block_force_xy_v1",
            "force_observable_family": "femm_weighted_stress_torque",
        },
        expected_force_observable_id="femm_slot300_weighted_stress_block_force_xy_v1",
        expected_force_observable_family="femm_weighted_stress_block_force_xy",
        expected_block_integral_types=(18, 19),
    )
    assert wrong_force_observable_family["status"] == "needs_attention"
    assert wrong_force_observable_family["checks"]["expected_force_observable_id_matches"] is True
    assert wrong_force_observable_family["checks"]["expected_force_observable_family_matches"] is False
    assert wrong_force_observable_family["checks"]["block_integral_types_match_expected"] is True

    stale_parameter_set_digest = parallel_wire_force_result_package_gate(
        {
            **row,
            "parameter_set_digest": "sha256:slot188-old-design-parameter-set",
        },
        expected_parameter_set_artifact_id="femm_slot392_parallel_wire_parameter_set_v1.json",
        expected_parameter_set_digest="sha256:slot392-parallel-wire-parameter-set-v1",
        expected_parameter_set_path="artifacts/femm/slot392_parallel_wire_parameter_set.json",
        require_parameter_set_artifact=True,
    )
    assert stale_parameter_set_digest["status"] == "needs_attention"
    assert stale_parameter_set_digest["checks"]["expected_parameter_set_artifact_id_matches"] is True
    assert stale_parameter_set_digest["checks"]["expected_parameter_set_digest_matches"] is False
    assert stale_parameter_set_digest["checks"]["expected_parameter_set_path_matches"] is True
    assert stale_parameter_set_digest["checks"]["vector_force_matches_ampere"] is True

    missing_parameter_set_path = parallel_wire_force_result_package_gate(
        {
            **row,
            "parameter_set_path": "",
        },
        require_parameter_set_artifact=True,
    )
    assert missing_parameter_set_path["status"] == "needs_attention"
    assert missing_parameter_set_path["checks"]["parameter_set_artifact_id_recorded"] is True
    assert missing_parameter_set_path["checks"]["parameter_set_digest_recorded"] is True
    assert missing_parameter_set_path["checks"]["parameter_set_path_recorded"] is False

    wrong_objective_family = parallel_wire_force_result_package_gate(
        {
            **row,
            "objective_observable_id": "femm_slot392_force_objective_radial_N_per_m_v1",
            "objective_observable_family": "torque_ripple_objective",
        },
        expected_objective_observable_id="femm_slot392_force_objective_radial_N_per_m_v1",
        expected_objective_observable_family="force_minimize_abs_radial_N_per_m",
    )
    assert wrong_objective_family["status"] == "needs_attention"
    assert wrong_objective_family["checks"]["expected_objective_observable_id_matches"] is True
    assert wrong_objective_family["checks"]["expected_objective_observable_family_matches"] is False
    assert wrong_objective_family["checks"]["radial_force_matches_ampere"] is True

    missing_force_y_command = parallel_wire_force_result_package_gate(
        {
            **row,
            "postprocess_trace_id": "femm_slot284_postprocess_trace_v1",
            "postprocess_command_digest": "sha256:slot284_clear_select_force_xy",
            "postprocess_commands": [
                "mo_clearblock()",
                "mo_groupselectblock(12)",
                "mo_blockintegral(18)",
            ],
        },
        expected_target_group_id=12,
        expected_postprocess_trace_id="femm_slot284_postprocess_trace_v1",
        expected_postprocess_command_digest="sha256:slot284_clear_select_force_xy",
        require_postprocess_command_trace=True,
    )
    assert missing_force_y_command["status"] == "needs_attention"
    assert missing_force_y_command["checks"]["postprocess_commands_clear_select_force_xy"] is False

    stale_selection_scope = parallel_wire_force_result_package_gate(
        row,
        expected_target_group_id=12,
        require_selection_clear=True,
    )
    assert stale_selection_scope["status"] == "needs_attention"
    assert stale_selection_scope["checks"]["selection_mentions_target_group"] is True
    assert stale_selection_scope["checks"]["selection_clear_before_groupselect"] is False

    stale_result = parallel_wire_force_result_package_gate(
        row,
        expected_model_id="slot188_parallel_wire_pair_v1",
        expected_operating_point_id="I1_10A_I2_5A_d_20mm",
        expected_artifact_id="femm_slot228_parallel_wire_recovered_table_A",
        expected_result_set_id="rs_I1_10A_I2_5A_d_20mm_v2",
        expected_source_tool="FEMM",
        expected_source_group_id=11,
        expected_target_group_id=12,
        expected_block_integral_types=(18, 19),
    )
    assert stale_result["status"] == "needs_attention"
    assert stale_result["checks"]["expected_result_set_id_matches"] is False
    assert stale_result["checks"]["vector_force_matches_ampere"] is True

    stale_solution = parallel_wire_force_result_package_gate(
        {
            **row,
            "solution_artifact_id": "femm_slot244_parallel_wire_pair_old.ans",
            "solution_loaded": True,
        },
        expected_solution_artifact_id="femm_slot252_parallel_wire_pair_v1.ans",
        require_solution_loaded=True,
    )
    assert stale_solution["status"] == "needs_attention"
    assert stale_solution["checks"]["expected_solution_artifact_id_matches"] is False
    assert stale_solution["checks"]["solution_loaded_before_postprocess"] is True
    assert stale_solution["checks"]["vector_force_matches_ampere"] is True

    missing_loaded_state = parallel_wire_force_result_package_gate(
        {**row, "solution_loaded": False},
        expected_solution_artifact_id="femm_slot252_parallel_wire_pair_v1.ans",
        require_solution_loaded=True,
    )
    assert missing_loaded_state["status"] == "needs_attention"
    assert missing_loaded_state["checks"]["expected_solution_artifact_id_matches"] is True
    assert missing_loaded_state["checks"]["solution_loaded_recorded"] is True
    assert missing_loaded_state["checks"]["solution_loaded_before_postprocess"] is False

    wrong_problem_type = parallel_wire_force_result_package_gate(
        {**row, "problem_type": "axisymmetric"},
        expected_problem_type="planar",
    )
    assert wrong_problem_type["status"] == "needs_attention"
    assert wrong_problem_type["checks"]["expected_problem_type_matches"] is False
    assert wrong_problem_type["checks"]["vector_force_matches_ampere"] is True

    wrong_length_unit = parallel_wire_force_result_package_gate(
        {**row, "length_unit": "millimeters"},
        expected_length_unit="meters",
    )
    assert wrong_length_unit["status"] == "needs_attention"
    assert wrong_length_unit["checks"]["expected_length_unit_matches"] is False
    assert wrong_length_unit["checks"]["vector_force_matches_ampere"] is True

    wrong_frequency = parallel_wire_force_result_package_gate(
        {**row, "frequency_hz": 50.0},
        expected_frequency_hz=0.0,
    )
    assert wrong_frequency["status"] == "needs_attention"
    assert wrong_frequency["checks"]["expected_frequency_hz_matches"] is False
    assert wrong_frequency["checks"]["vector_force_matches_ampere"] is True

    coarse_precision = parallel_wire_force_result_package_gate(
        {**row, "solver_precision": 1.0e-3},
        max_solver_precision=1.0e-7,
    )
    assert coarse_precision["status"] == "needs_attention"
    assert coarse_precision["checks"]["solver_precision_within_max"] is False
    assert coarse_precision["checks"]["vector_force_matches_ampere"] is True

    wrong_min_angle = parallel_wire_force_result_package_gate(
        {**row, "min_angle_deg": 5.0},
        expected_min_angle_deg=30.0,
    )
    assert wrong_min_angle["status"] == "needs_attention"
    assert wrong_min_angle["checks"]["expected_min_angle_deg_matches"] is False
    assert wrong_min_angle["checks"]["vector_force_matches_ampere"] is True

    stale_block_labels = parallel_wire_force_result_package_gate(
        {
            **row,
            "block_label_artifact_id": "femm_slot204_parallel_wire_block_labels_old.json",
        },
        expected_block_label_artifact_id="femm_slot276_parallel_wire_block_labels_v1.json",
        expected_source_region="wire_1_positive_current",
        expected_target_region="wire_2_positive_current",
        expected_source_material="Copper",
        expected_target_material="Copper",
    )
    assert stale_block_labels["status"] == "needs_attention"
    assert stale_block_labels["checks"]["expected_block_label_artifact_id_matches"] is False
    assert stale_block_labels["checks"]["expected_source_region_matches"] is True
    assert stale_block_labels["checks"]["vector_force_matches_ampere"] is True

    stale_target_material = parallel_wire_force_result_package_gate(
        {**row, "target_material": "Air"},
        expected_block_label_artifact_id="femm_slot276_parallel_wire_block_labels_v1.json",
        expected_target_material="Copper",
    )
    assert stale_target_material["status"] == "needs_attention"
    assert stale_target_material["checks"]["expected_block_label_artifact_id_matches"] is True
    assert stale_target_material["checks"]["expected_target_material_matches"] is False
    assert stale_target_material["checks"]["vector_force_matches_ampere"] is True

    repulsion = parallel_wire_force_result_package_gate(
        {
            **row,
            "current2_A": -5.0,
            "radial_force_on_wire2_N_per_m": pair["force_magnitude_per_length_N_per_m"],
            "force_on_wire2_N_per_m": [0.0005, 0.0],
            "interaction": "repulsion",
        }
    )
    assert repulsion["status"] == "ok"
    assert repulsion["analytic"]["interaction"] == "repulsion"

    wrong_sign = parallel_wire_force_result_package_gate(
        {
            **row,
            "radial_force_on_wire2_N_per_m": pair["force_magnitude_per_length_N_per_m"],
            "force_on_wire2_N_per_m": [0.0005, 0.0],
        }
    )
    assert wrong_sign["status"] == "needs_attention"
    assert wrong_sign["checks"]["vector_force_matches_ampere"] is False
    assert wrong_sign["checks"]["radial_force_matches_ampere"] is False

    missing_units = parallel_wire_force_result_package_gate({**row, "force_units": "N"})
    assert missing_units["status"] == "needs_attention"
    assert missing_units["checks"]["force_units_are_per_length"] is False

    total_force_row = dict(row)
    total_force_row.pop("force_unit_basis")
    total_force = parallel_wire_force_result_package_gate({**total_force_row, "force_units": "N"})
    assert total_force["status"] == "needs_attention"
    assert total_force["force_unit_basis"] == "depth_integrated"
    assert total_force["checks"]["force_unit_basis_is_per_length"] is False
    assert total_force["checks"]["depth_integrated_force_not_used_for_per_length_gate"] is False

    missing_depth_row = dict(row)
    missing_depth_row.pop("problem_depth_m")
    missing_depth = parallel_wire_force_result_package_gate(missing_depth_row)
    assert missing_depth["status"] == "needs_attention"
    assert missing_depth["checks"]["femm_planar_depth_recorded"] is False
    assert missing_depth["checks"]["force_unit_basis_is_per_length"] is True

    missing_projection_axis = parallel_wire_force_result_package_gate(
        {key: value for key, value in row.items() if key != "radial_projection_axis"}
    )
    assert missing_projection_axis["status"] == "needs_attention"
    assert missing_projection_axis["checks"]["radial_projection_axis_recorded"] is False
    assert missing_projection_axis["checks"]["vector_force_matches_ampere"] is True

    vague_projection_axis = parallel_wire_force_result_package_gate(
        {**row, "radial_projection_axis": "local_x"}
    )
    assert vague_projection_axis["status"] == "needs_attention"
    assert vague_projection_axis["checks"]["radial_projection_axis_recorded"] is True
    assert vague_projection_axis["checks"]["radial_projection_axis_names_wire_pair"] is False

    wrong_component_frame = parallel_wire_force_result_package_gate(
        {**row, "force_component_frame": "local_rt"},
        expected_force_component_frame="global_cartesian_xy",
    )
    assert wrong_component_frame["status"] == "needs_attention"
    assert wrong_component_frame["checks"]["force_component_frame_recorded_when_expected"] is True
    assert wrong_component_frame["checks"]["expected_force_component_frame_matches"] is False
    assert wrong_component_frame["checks"]["vector_force_matches_ampere"] is True

    wrong_projection_convention = parallel_wire_force_result_package_gate(
        {
            **row,
            "radial_projection_axis": "wire2_to_wire1_separation_axis_positive_away_from_wire2",
        },
        expected_radial_projection_axis="wire1_to_wire2_separation_axis_positive_away_from_wire1",
    )
    assert wrong_projection_convention["status"] == "needs_attention"
    assert wrong_projection_convention["checks"]["radial_projection_axis_recorded_when_expected"] is True
    assert wrong_projection_convention["checks"]["expected_radial_projection_axis_matches"] is False
    assert wrong_projection_convention["checks"]["radial_projection_axis_names_wire_pair"] is True

    wrong_sign_convention = parallel_wire_force_result_package_gate(
        {
            **row,
            "force_sign_convention": "positive_radial_force_points_toward_wire1_attraction_is_positive",
        },
        expected_force_sign_convention="positive_radial_force_points_away_from_wire1_attraction_is_negative",
    )
    assert wrong_sign_convention["status"] == "needs_attention"
    assert wrong_sign_convention["checks"]["force_sign_convention_recorded_when_expected"] is True
    assert wrong_sign_convention["checks"]["expected_force_sign_convention_matches"] is False
    assert wrong_sign_convention["checks"]["vector_force_matches_ampere"] is True

    missing_sign_convention_row = dict(row)
    missing_sign_convention_row.pop("force_sign_convention")
    missing_sign_convention = parallel_wire_force_result_package_gate(
        missing_sign_convention_row,
        expected_force_sign_convention="positive_radial_force_points_away_from_wire1_attraction_is_negative",
    )
    assert missing_sign_convention["status"] == "needs_attention"
    assert missing_sign_convention["checks"]["force_sign_convention_recorded_when_expected"] is False
    assert missing_sign_convention["checks"]["expected_force_sign_convention_matches"] is False
    assert missing_sign_convention["checks"]["vector_force_matches_ampere"] is True

    wrong_extraction_method = parallel_wire_force_result_package_gate(
        {**row, "force_extraction_method": "lorentz_body_force"},
        expected_force_extraction_method="weighted_stress_block_integral_xy",
        expected_block_integral_types=(18, 19),
    )
    assert wrong_extraction_method["status"] == "needs_attention"
    assert wrong_extraction_method["checks"]["force_extraction_method_recorded_when_expected"] is True
    assert wrong_extraction_method["checks"]["expected_force_extraction_method_matches"] is False
    assert wrong_extraction_method["checks"]["block_integral_types_match_expected"] is True
    assert wrong_extraction_method["checks"]["vector_force_matches_ampere"] is True

    weighted_stress_without_y_integral = parallel_wire_force_result_package_gate(
        {**row, "block_integral_types": [18]},
        expected_force_extraction_method="weighted_stress_block_integral_xy",
    )
    assert weighted_stress_without_y_integral["status"] == "needs_attention"
    assert weighted_stress_without_y_integral["checks"]["expected_force_extraction_method_matches"] is True
    assert weighted_stress_without_y_integral["checks"]["weighted_stress_extraction_uses_force_xy_integrals"] is False

    wrong_group = parallel_wire_force_result_package_gate(
        row,
        expected_source_group_id=11,
        expected_target_group_id=13,
    )
    assert wrong_group["status"] == "needs_attention"
    assert wrong_group["checks"]["expected_source_group_id_matches"] is True
    assert wrong_group["checks"]["expected_target_group_id_matches"] is False
    assert wrong_group["checks"]["selection_mentions_target_group"] is False

    stale_target_center = parallel_wire_force_result_package_gate(
        {**row, "target_center_xy_m": [0.03, 0.0]},
        expected_source_center_xy_m=(0.0, 0.0),
        expected_target_center_xy_m=(0.02, 0.0),
    )
    assert stale_target_center["status"] == "needs_attention"
    assert stale_target_center["checks"]["expected_source_center_xy_matches"] is True
    assert stale_target_center["checks"]["expected_target_center_xy_matches"] is False
    assert stale_target_center["checks"]["wire_center_separation_matches_separation_m"] is False
    assert stale_target_center["checks"]["vector_force_matches_ampere"] is True

    missing_source_center_row = dict(row)
    missing_source_center_row.pop("source_center_xy_m")
    missing_source_center = parallel_wire_force_result_package_gate(
        missing_source_center_row,
        expected_source_center_xy_m=(0.0, 0.0),
        expected_target_center_xy_m=(0.02, 0.0),
    )
    assert missing_source_center["status"] == "needs_attention"
    assert missing_source_center["checks"]["source_center_xy_recorded_when_expected"] is False
    assert missing_source_center["checks"]["target_center_xy_recorded_when_expected"] is True
    assert missing_source_center["checks"]["expected_source_center_xy_matches"] is False

    stale_center_distance = parallel_wire_force_result_package_gate(
        {**row, "target_center_xy_m": [0.021, 0.0]}
    )
    assert stale_center_distance["status"] == "needs_attention"
    assert stale_center_distance["checks"]["wire_center_separation_matches_separation_m"] is False
    assert stale_center_distance["checks"]["vector_force_matches_ampere"] is True

    stale_current_artifact = parallel_wire_force_result_package_gate(
        {**row, "current_source_artifact_id": "femm_slot188_old_current_table.json"},
        expected_current_source_artifact_id="femm_slot332_parallel_wire_current_definition_v1.json",
        expected_current_definition_method="femm_circuit_current_snapshot",
    )
    assert stale_current_artifact["status"] == "needs_attention"
    assert stale_current_artifact["checks"]["expected_current_source_artifact_id_matches"] is False
    assert stale_current_artifact["checks"]["expected_current_definition_method_matches"] is True
    assert stale_current_artifact["checks"]["vector_force_matches_ampere"] is True

    missing_current_method_row = dict(row)
    missing_current_method_row.pop("current_definition_method")
    missing_current_method = parallel_wire_force_result_package_gate(
        missing_current_method_row,
        expected_current_definition_method="femm_circuit_current_snapshot",
    )
    assert missing_current_method["status"] == "needs_attention"
    assert missing_current_method["checks"]["current_definition_method_recorded_when_expected"] is False
    assert missing_current_method["checks"]["expected_current_definition_method_matches"] is False
    assert missing_current_method["checks"]["vector_force_matches_ampere"] is True

    wrong_current_method = parallel_wire_force_result_package_gate(
        {**row, "current_definition_method": "rms_current_table"},
        expected_current_source_artifact_id="femm_slot332_parallel_wire_current_definition_v1.json",
        expected_current_definition_method="femm_circuit_current_snapshot",
    )
    assert wrong_current_method["status"] == "needs_attention"
    assert wrong_current_method["checks"]["expected_current_source_artifact_id_matches"] is True
    assert wrong_current_method["checks"]["expected_current_definition_method_matches"] is False
    assert wrong_current_method["checks"]["vector_force_matches_ampere"] is True

    same_group = parallel_wire_force_result_package_gate(
        {**row, "target_group_id": 11, "selection_function": "mo_groupselectblock(11); mo_blockintegral(18)"},
        expected_source_group_id=11,
        expected_target_group_id=11,
    )
    assert same_group["status"] == "needs_attention"
    assert same_group["checks"]["expected_source_group_id_matches"] is True
    assert same_group["checks"]["expected_target_group_id_matches"] is True
    assert same_group["checks"]["selection_mentions_target_group"] is True
    assert same_group["checks"]["source_target_groups_distinct"] is False

    torque_integral = parallel_wire_force_result_package_gate(
        {**row, "block_integral_types": [22], "selection_function": "mo_groupselectblock(12); mo_blockintegral(22)"},
        expected_source_group_id=11,
        expected_target_group_id=12,
        expected_block_integral_types=(18, 19),
    )
    assert torque_integral["status"] == "needs_attention"
    assert torque_integral["checks"]["block_integral_types_match_expected"] is False
    assert torque_integral["checks"]["block_integral_types_are_force_xy"] is False
    assert torque_integral["checks"]["block_integral_types_exclude_torque"] is False

    elf_pair = parallel_wire_lorentz_force_summary(5.0, 8.0, 0.02)
    elf_helper_row = {
        "model_id": "slot190_parallel_wire_analytic_reference",
        "operating_point_id": "I1_5A_I2_8A_d_20mm",
        "source_tool": "ELF/MAGIC",
        "source_function": "two_line_current_force_per_length",
        "current1_A": 5.0,
        "current2_A": 8.0,
        "separation_m": 0.02,
        "force_units": "N/m",
        "radial_force_on_wire2_N_per_m": -elf_pair["force_magnitude_per_length_N_per_m"],
        "interaction": "attraction",
    }
    elf_gate = parallel_wire_force_result_package_gate(elf_helper_row)
    assert elf_gate["status"] == "ok"
    assert elf_gate["source_tool"] == "ELF/MAGIC"
    assert elf_gate["source_function"] == "two_line_current_force_per_length"


def test_magnetic_field_probe_result_package_gate_binds_solution_point_and_output():
    row = {
        "model_id": "slot348_parallel_wire_probe_model",
        "operating_point_id": "I10A_I5A_static",
        "artifact_id": "femm_slot348_probe_case_v1",
        "solution_artifact_id": "femm_slot348_probe_case_v1.ans",
        "solution_digest": "sha256:femm-slot348-probe-case-ans-v1",
        "solution_path": r"artifacts/femm/femm_slot348_probe_case_v1.ans",
        "solution_loaded": True,
        "source_tool": "FEMM",
        "source_function": "mo_getb(0.010, 0.0)",
        "field_probe_id": "slot348_midgap_B_probe_v1",
        "probe_point_xy_m": [0.010, 0.0],
        "problem_length_unit": "meters",
        "probe_point_input_xy": [0.010, 0.0],
        "probe_point_input_unit": "meters",
        "coordinate_scale_to_m": 1.0,
        "B_T": [1.2e-4, -3.0e-5],
        "field_units": "T",
        "field_component_frame": "global_cartesian_xy",
        "field_probe_method": "femm_mo_getb_point_sample",
        "postprocess_trace_id": "slot348_mo_getb_trace_v1",
        "postprocess_command_digest": "sha256:slot348-mo-getb-trace-v1",
        "postprocess_commands": [
            "mi_loadsolution()",
            "mo_getb(0.010, 0.0)",
        ],
        "field_probe_output_artifact_id": "slot348_probe_table_v1",
        "field_probe_output_digest": "sha256:slot348-probe-table-v1",
        "field_probe_output_path": r"artifacts/femm/slot348_probe_table.json",
    }

    gate = magnetic_field_probe_result_package_gate(
        row,
        expected_model_id="slot348_parallel_wire_probe_model",
        expected_operating_point_id="I10A_I5A_static",
        expected_artifact_id="femm_slot348_probe_case_v1",
        expected_solution_artifact_id="femm_slot348_probe_case_v1.ans",
        expected_solution_digest="sha256:femm-slot348-probe-case-ans-v1",
        expected_solution_path=r"artifacts/femm/femm_slot348_probe_case_v1.ans",
        expected_source_tool="FEMM",
        expected_probe_id="slot348_midgap_B_probe_v1",
        expected_probe_point_xy_m=(0.010, 0.0),
        expected_problem_length_unit="meters",
        expected_probe_point_input_unit="meters",
        expected_coordinate_scale_to_m=1.0,
        expected_field_component_frame="global_cartesian_xy",
        expected_field_units="T",
        expected_field_probe_method="femm_mo_getb_point_sample",
        expected_postprocess_trace_id="slot348_mo_getb_trace_v1",
        expected_postprocess_command_digest="sha256:slot348-mo-getb-trace-v1",
        expected_probe_output_artifact_id="slot348_probe_table_v1",
        expected_probe_output_digest="sha256:slot348-probe-table-v1",
        require_solution_artifact=True,
        require_solution_loaded=True,
        require_postprocess_command_trace=True,
        require_probe_output_artifact=True,
        require_probe_coordinate_scale=True,
    )

    assert gate["policy"] == "magnetic_field_probe_result_package_gate"
    assert gate["status"] == "ok"
    assert gate["checks"]["expected_probe_point_xy_matches"] is True
    assert gate["checks"]["probe_point_input_scale_matches_probe_point_xy_m"] is True
    assert gate["checks"]["expected_coordinate_scale_to_m_matches"] is True
    assert gate["checks"]["field_vector_finite"] is True
    assert gate["checks"]["postprocess_commands_include_mo_getb"] is True
    assert gate["checks"]["expected_solution_digest_matches"] is True
    assert gate["checks"]["expected_solution_path_matches"] is True
    assert gate["checks"]["expected_probe_output_artifact_id_matches"] is True

    stale_solution = magnetic_field_probe_result_package_gate(
        {**row, "solution_artifact_id": "femm_slot340_old_force_case.ans"},
        expected_solution_artifact_id="femm_slot348_probe_case_v1.ans",
        require_solution_loaded=True,
    )
    assert stale_solution["status"] == "needs_attention"
    assert stale_solution["checks"]["expected_solution_artifact_id_matches"] is False

    stale_solution_digest = magnetic_field_probe_result_package_gate(
        {**row, "solution_digest": "sha256:old-slot348-ans"},
        expected_solution_artifact_id="femm_slot348_probe_case_v1.ans",
        expected_solution_digest="sha256:femm-slot348-probe-case-ans-v1",
        expected_solution_path=r"artifacts/femm/femm_slot348_probe_case_v1.ans",
        require_solution_artifact=True,
        require_solution_loaded=True,
    )
    assert stale_solution_digest["status"] == "needs_attention"
    assert stale_solution_digest["checks"]["expected_solution_artifact_id_matches"] is True
    assert stale_solution_digest["checks"]["expected_solution_digest_matches"] is False

    missing_solution_path = magnetic_field_probe_result_package_gate(
        {key: value for key, value in row.items() if key != "solution_path"},
        expected_solution_artifact_id="femm_slot348_probe_case_v1.ans",
        expected_solution_digest="sha256:femm-slot348-probe-case-ans-v1",
        expected_solution_path=r"artifacts/femm/femm_slot348_probe_case_v1.ans",
        require_solution_artifact=True,
        require_solution_loaded=True,
    )
    assert missing_solution_path["status"] == "needs_attention"
    assert missing_solution_path["checks"]["solution_path_recorded_when_required"] is False

    unloaded = magnetic_field_probe_result_package_gate(
        {**row, "solution_loaded": False},
        expected_solution_artifact_id="femm_slot348_probe_case_v1.ans",
        require_solution_loaded=True,
    )
    assert unloaded["status"] == "needs_attention"
    assert unloaded["checks"]["solution_loaded_true_when_required"] is False

    stale_point = magnetic_field_probe_result_package_gate(
        {**row, "probe_point_xy_m": [0.011, 0.0]},
        expected_probe_point_xy_m=(0.010, 0.0),
    )
    assert stale_point["status"] == "needs_attention"
    assert stale_point["checks"]["expected_probe_point_xy_matches"] is False

    wrong_coordinate_scale = magnetic_field_probe_result_package_gate(
        {**row, "coordinate_scale_to_m": 1.0e-3},
        expected_problem_length_unit="meters",
        expected_probe_point_input_unit="meters",
        expected_coordinate_scale_to_m=1.0,
        require_probe_coordinate_scale=True,
    )
    assert wrong_coordinate_scale["status"] == "needs_attention"
    assert wrong_coordinate_scale["checks"]["expected_coordinate_scale_to_m_matches"] is False
    assert wrong_coordinate_scale["checks"]["probe_point_input_scale_matches_probe_point_xy_m"] is False

    wrong_method = magnetic_field_probe_result_package_gate(
        {**row, "field_probe_method": "line_average_b_field"},
        expected_field_probe_method="femm_mo_getb_point_sample",
    )
    assert wrong_method["status"] == "needs_attention"
    assert wrong_method["checks"]["expected_field_probe_method_matches"] is False

    stale_output = magnetic_field_probe_result_package_gate(
        {**row, "field_probe_output_artifact_id": "slot340_old_probe_table"},
        expected_probe_output_artifact_id="slot348_probe_table_v1",
        expected_probe_output_digest="sha256:slot348-probe-table-v1",
        require_probe_output_artifact=True,
    )
    assert stale_output["status"] == "needs_attention"
    assert stale_output["checks"]["expected_probe_output_artifact_id_matches"] is False


def test_air_gap_force_scales_with_b_squared_area_and_faces():
    base = air_gap_holding_force(0.5, area_m2=2.0e-4)
    assert air_gap_holding_force(1.0, area_m2=2.0e-4) == pytest.approx(4.0 * base)
    assert air_gap_holding_force(0.5, area_m2=4.0e-4) == pytest.approx(2.0 * base)
    assert air_gap_holding_force(0.5, area_m2=2.0e-4, faces=2) == pytest.approx(2.0 * base)


def test_air_gap_force_summary_is_json_friendly_and_self_consistent():
    row = air_gap_force_summary(0.8, area_m2=1.5e-4, faces=2)
    pressure = 0.8 * 0.8 / (2.0 * MU0)
    assert row["B_T"] == pytest.approx(0.8)
    assert row["pressure_Pa"] == pytest.approx(pressure)
    assert row["energy_density_J_per_m3"] == pytest.approx(pressure)
    assert row["force_N"] == pytest.approx(pressure * 1.5e-4 * 2)
    assert row["force_per_area_N_per_m2"] == pytest.approx(pressure)


def test_air_gap_force_rejects_invalid_inputs():
    with pytest.raises(ValueError):
        air_gap_maxwell_pressure(1.0, mu=0.0)
    with pytest.raises(ValueError):
        air_gap_holding_force(1.0, area_m2=-1.0)
    with pytest.raises(ValueError):
        air_gap_holding_force(1.0, area_m2=1.0, faces=0)
    with pytest.raises(ValueError):
        maxwell_stress_tensor_air((1.0,), mu=MU0)
    with pytest.raises(ValueError):
        maxwell_traction_air((1.0, 0.0, 0.0), (0.0, 0.0, 0.0))
    with pytest.raises(ValueError):
        maxwell_traction_summary((1.0, 0.0), (1.0, 0.0, 0.0))
    with pytest.raises(ValueError):
        air_gap_shear_stress(1.0, 0.1, mu=0.0)
    with pytest.raises(ValueError):
        air_gap_shear_torque(1.0, 0.1, radius_m=-1.0)
    with pytest.raises(ValueError):
        air_gap_shear_torque_summary(1.0, 0.1, radius_m=1.0, axial_length_m=-1.0)
    with pytest.raises(ValueError):
        air_gap_shear_torque_from_angle_samples([0.0], [1.0], [0.1], radius_m=1.0)
    with pytest.raises(ValueError):
        air_gap_shear_torque_from_angle_samples([0.0, 0.1], [1.0], [0.1, 0.1], radius_m=1.0)
    with pytest.raises(ValueError):
        air_gap_shear_torque_from_angle_samples([0.0, 0.0], [1.0, 1.0], [0.1, 0.1], radius_m=1.0)
    with pytest.raises(ValueError):
        air_gap_shear_torque_from_angle_samples([0.0, 0.1], [1.0, 1.0], [0.1, 0.1], radius_m=-1.0)
    with pytest.raises(ValueError):
        air_gap_shear_torque_from_angle_samples([0.0, 0.1], [1.0, 1.0], [0.1, 0.1], radius_m=1.0, period_rad=0.0)
    with pytest.raises(ValueError):
        maxwell_line_segment_force_2d((0.0, 0.0), (0.0, 0.0), (1.0, 0.0))
    with pytest.raises(ValueError):
        maxwell_contour_force_2d([(0.0, 0.0), (1.0, 0.0)], (1.0, 0.0))
    with pytest.raises(ValueError):
        maxwell_contour_segment_balance_summary_2d(
            [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)],
            (1.0, 0.0),
            expected_force_per_depth_N_per_m=(0.0, 0.0, 0.0),
        )
    with pytest.raises(ValueError):
        force_moment_resultant_summary([], [])
    with pytest.raises(ValueError):
        force_moment_resultant_summary([(0.0, 0.0)], [(1.0, 0.0, 0.0)])


if __name__ == "__main__":
    test_air_gap_pressure_matches_maxwell_stress_at_one_tesla()
    test_maxwell_tensor_normal_field_reduces_to_air_gap_pressure()
    test_maxwell_tensor_tangential_field_is_magnetic_tension()
    test_maxwell_traction_oblique_field_decomposes_into_normal_and_tangent()
    test_air_gap_shear_stress_matches_maxwell_tangential_traction()
    test_air_gap_shear_torque_scales_with_radius_length_angle_and_sign()
    test_air_gap_sampled_shear_torque_uniform_matches_closed_form()
    test_air_gap_sampled_shear_torque_sinusoidal_average()
    test_motor_air_gap_harmonic_torque_phase_gate()
    test_force_moment_resultant_summary_handles_2d_force_couple()
    test_force_moment_resultant_summary_handles_3d_single_force()
    test_maxwell_line_segment_force_2d_matches_air_gap_pressure()
    test_maxwell_contour_force_2d_closed_uniform_field_cancels()
    test_maxwell_contour_segment_balance_summary_reports_cancellation()
    test_air_gap_force_scales_with_b_squared_area_and_faces()
    test_air_gap_force_summary_is_json_friendly_and_self_consistent()
    test_air_gap_force_rejects_invalid_inputs()
    print("[OK] air-gap Maxwell pressure and holding-force helpers validated.")
