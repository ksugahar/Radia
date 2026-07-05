import math

import pytest

from radia_mcp.radia_ngsolve.slot_gates import (
    acoustic_impedance_power_result_package_gate,
    acoustic_interface_result_package_gate,
    acoustic_plane_wave_intensity_convention_gate,
    acoustic_normal_incidence_interface_gate,
    balanced_back_emf_line_voltage_handoff_gate,
    box_projected_gradient_least_squares_gate,
    branch_line_hybrid_gate,
    carter_slot_opening_sweep_gate,
    coenergy_torque_periodic_summary,
    coaxial_rc_duality_gate,
    coaxial_pm_force_gap_sweep_gate,
    computed_reference_rows_gate,
    cross_validation_artifact_to_mcp_feedback_gate,
    source_native_seed_queue_gate,
    cst_abcd_cascade_solver_ready_manifest_gate,
    cst_export_manifest_solver_ready_gate,
    cst_result_export_package_gate,
    cst_touchstone_solver_ready_manifest_gate,
    dq_current_from_gamma_deg,
    dq_torque_table_health,
    dq_to_three_phase_currents,
    drive_cycle_weighted_efficiency_gate,
    double_layer_winding_pitch_harmonic_gate,
    femm_air_gap_sample_solver_ready_manifest_gate,
    femm_block_label_source_contract_gate,
    femm_group_motion_selection_gate,
    femm_motor_model_artifact_package_gate,
    femm_pm_magnetization_convention_gate,
    femm_source_current_solver_ready_manifest_gate,
    femm_static_current_circuit_rows_gate,
    femm_winding_current_package_gate,
    farfield_lobe_notebook_handoff_gate,
    farfield_pattern_metadata_gate,
    flux_linkage_back_emf_derivative_gate,
    geometric_integrator_energy_drift_gate,
    jmag_angle_alignment_contract_gate,
    jmag_current_torque_solver_ready_manifest_gate,
    jmag_export_case_package_gate,
    jmag_efficiency_operating_point_package_gate,
    jmag_force_table_metadata_gate,
    jmag_airgap_flux_sample_metadata_gate,
    jmag_airgap_torque_integration_package_gate,
    jmag_pm_short_circuit_fault_table_gate,
    ipm_saliency_torque_component_gate,
    inverter_dc_bus_voltage_limit_gate,
    jmag_motor_table_column_metadata_gate,
    jmag_symmetry_sweep_coverage_gate,
    lcurve_corner_choice,
    lumped_pm_dq_torque,
    maxwell_stress_surface_package_gate,
    mesh_import_quality_manifest_gate,
    motor_current_snapshot_table_contract_gate,
    morozov_discrepancy_choice,
    mqs_coulomb_gauge_efield_postprocess_gate,
    netgen_vol_boundary_orientation_trace_package_gate,
    netgen_vol_fem_bem_normal_flux_sign_package_gate,
    netgen_vol_first_order_fem_bem_trace_package_handoff,
    owned_solver_model_tag_lifecycle_gate,
    one_port_match_quality_gate,
    parallel_wire_force_per_length,
    pm_bem_surface_normal_metadata_gate,
    pm_demag_margin_screening_package_gate,
    pm_demag_package_identity_gate,
    pm_drive_loss_bucket_efficiency_gate,
    pm_drive_terminal_table_health,
    pm_loadline_metadata_gate,
    pm_recoil_demag_step_summary,
    quarter_wave_directional_coupler_gate,
    shared_solver_session_health_gate,
    solver_submodel_boundary_handoff_gate,
    solver_result_artifact_provenance_timing_gate,
    solver_result_table_metadata_gate,
    spwm_snapshot_current_handoff_summary,
    spherical_dirichlet_laplacian_eigen_gate,
    three_phase_currents_to_dq_summary,
    trace_surface_mass_energy_gate,
    thermal_annulus_conductance_gate,
    thermal_conduction_convection_robin_gate,
    thermal_layer_stack_conductance_gate,
    touchstone_frequency_grid_interpolation_gate,
    touchstone_frequency_unit_normalization_gate,
    touchstone_port_metadata_gate,
    touchstone_power_wave_balance_gate,
    touchstone_row_solver_ready_preflight_gate,
    touchstone_sparameter_to_complex,
    two_port_abcd_cascade_gate,
    two_port_s_to_yz_equivalent_gate,
    two_port_sparameter_health,
)


def test_parallel_wire_force_gate_signed_and_scaled():
    f = parallel_wire_force_per_length(10.0, 20.0, 0.05)
    assert f == pytest.approx(8.0e-4)
    assert parallel_wire_force_per_length(10.0, -20.0, 0.05) == pytest.approx(-f)
    assert parallel_wire_force_per_length(10.0, 20.0, 0.10) == pytest.approx(0.5 * f)
    with pytest.raises(ValueError):
        parallel_wire_force_per_length(1.0, 1.0, 0.0)


def test_coaxial_rc_duality_gate_keeps_capacitance_and_resistance_consistent():
    gate = coaxial_rc_duality_gate(
        inner_radius=0.01,
        outer_radius=0.05,
        eps_r=1.0,
        sigma=2.5,
        length=1.0,
        measured_capacitance=3.4566418309218263e-11,
        measured_resistance=math.log(5.0) / (2.0 * math.pi * 2.5),
        rtol=1.0e-6,
        atol=1.0e-18,
    )

    assert gate["status"] == "ok"
    assert gate["policy"] == "coaxial_rc_duality_gate"
    assert gate["capacitance_F"] == pytest.approx(3.456641746947256e-11)
    assert gate["resistance_ohm"] == pytest.approx(math.log(5.0) / (2.0 * math.pi * 2.5))
    assert gate["rc_reference_s"] == pytest.approx(8.8541878128e-12 / 2.5)
    assert gate["measured_capacitance_rel_error"] < 3.0e-8
    assert gate["checks"]["geometry_free_rc_duality_ok"] is True

    bad = coaxial_rc_duality_gate(
        inner_radius=0.01,
        outer_radius=0.05,
        measured_capacitance=4.0e-11,
        rtol=1.0e-6,
    )
    assert bad["status"] == "needs_attention"
    assert bad["checks"]["measured_capacitance_ok"] is False


def test_coenergy_torque_gate_uses_absolute_tolerance_at_zero_crossings():
    n = 64
    amp = 0.25
    theta = [2.0 * math.pi * i / n for i in range(n)]
    # W' = -A cos(theta); T = dW'/dtheta = A sin(theta).
    coenergy = [-amp * math.cos(t) for t in theta]
    torque = [amp * math.sin(t) for t in theta]

    summary = coenergy_torque_periodic_summary(
        theta,
        coenergy,
        torque,
        rtol=2.0e-3,
        atol=1.0e-12,
    )

    assert summary["status"] == "ok"
    assert summary["max_abs_error"] < 5.0e-4
    assert summary["near_zero_abs_tolerance_schema_id"] == "coenergy_torque_near_zero_abs_tolerance_v1"
    assert summary["near_zero_row_count"] >= 2
    assert summary["near_zero_rows_pass_absolute_tolerance"] is True
    assert summary["checks"]["near_zero_rows_use_absolute_tolerance"] is True
    zero_rows = [row for row in summary["rows"] if abs(row["reference_torque_nm"]) < 1.0e-12]
    assert zero_rows
    assert max(row["abs_error"] for row in zero_rows) < 1.0e-12


def test_coenergy_torque_gate_improves_with_angle_resolution():
    amp = 0.25
    summaries = []
    for n in (64, 256):
        theta = [2.0 * math.pi * i / n for i in range(n)]
        coenergy = [-amp * math.cos(t) for t in theta]
        torque = [amp * math.sin(t) for t in theta]
        summaries.append(coenergy_torque_periodic_summary(
            theta,
            coenergy,
            torque,
            rtol=2.0e-3,
            atol=1.0e-12,
        ))

    coarse, fine = summaries
    assert coarse["status"] == "ok"
    assert fine["status"] == "ok"
    assert fine["max_abs_error"] < 0.07 * coarse["max_abs_error"]


def test_carter_slot_opening_gate_tracks_femm_slotting_preflight():
    summary = carter_slot_opening_sweep_gate(0.012, 0.001, [0.0, 0.002, 0.004, 0.006])

    assert summary["status"] == "ok"
    assert summary["checks"]["zero_opening_identity"] is True
    assert summary["checks"]["kc_nondecreasing_with_opening"] is True
    assert summary["checks"]["permeance_nonincreasing_with_opening"] is True
    assert summary["rows"][0]["carter_coefficient"] == pytest.approx(1.0)
    assert summary["rows"][2]["carter_coefficient"] == pytest.approx(1.1739, abs=2.0e-3)
    assert summary["rows"][2]["permeance_factor"] == pytest.approx(1.0 / summary["rows"][2]["carter_coefficient"])
    assert summary["min_permeance_factor"] < 0.8


def test_double_layer_winding_pitch_harmonic_gate_tracks_femm_pitch_signs():
    gate = double_layer_winding_pitch_harmonic_gate(
        q=2,
        coil_pitch=5.0,
        pole_pitch=6.0,
        phases=3,
        harmonics=(1, 3, 5, 7),
        expected_kw_signs={1: 1, 5: 1, 7: -1},
        expected_kp_signs={5: 1, 7: 1},
    )

    rows = {row["harmonic"]: row for row in gate["rows"]}
    assert gate["status"] == "ok"
    assert gate["policy"] == "double_layer_winding_pitch_harmonic_gate"
    assert gate["slot_angle_electrical_deg"] == pytest.approx(30.0)
    assert rows[1]["distribution_factor_kd"] == pytest.approx(0.9659258262890683)
    assert rows[1]["pitch_factor_kp"] == pytest.approx(0.9659258262890683)
    assert rows[1]["winding_factor_kw"] == pytest.approx(0.9330127018922194)
    assert rows[5]["pitch_factor_kp"] == pytest.approx(0.25881904510252074)
    assert rows[5]["winding_factor_kw"] == pytest.approx(0.06698729810778076)
    assert rows[7]["pitch_factor_kp"] > 0.0
    assert rows[7]["distribution_factor_kd"] < 0.0
    assert rows[7]["winding_factor_kw"] == pytest.approx(-0.06698729810778058)

    bad = double_layer_winding_pitch_harmonic_gate(
        q=2,
        coil_pitch=5.0,
        pole_pitch=6.0,
        expected_kw_signs={7: 1},
    )
    assert bad["status"] == "needs_attention"
    assert bad["checks"]["kw_sign_h7_ok"] is False


def test_two_port_sparameter_health_checks_passivity_and_reciprocity():
    health = two_port_sparameter_health(0.1, 0.7, s12=0.7, s22=0.1)
    assert health["status"] == "ok"
    assert health["reciprocal"] is True
    assert health["passive"] is True
    assert health["passive_margin"] > 0.0

    nonreciprocal = two_port_sparameter_health(0.1, 0.7, s12=0.6, s22=0.1)
    assert nonreciprocal["reciprocal"] is False
    assert nonreciprocal["status"] == "needs_attention"

    active = two_port_sparameter_health(0.2, 1.1, s12=1.1, s22=0.2)
    assert active["passive"] is False
    assert active["status"] == "needs_attention"


def test_touchstone_power_wave_balance_gate_records_basis_and_absorbed_power():
    gate = touchstone_power_wave_balance_gate(
        [0.1, 0.0],
        [0.8, 0.0],
        s12=[0.8, 0.0],
        s22=[0.1, 0.0],
        data_format="MA",
        power_balance_basis="power_waves_unit_incident_port",
        export_artifact_id="cst_touchstone_export_A",
        expected_export_artifact_id="cst_touchstone_export_A",
        result_set_id="cst_resultset_rf_001",
        expected_result_set_id="cst_resultset_rf_001",
    )

    assert gate["policy"] == "touchstone_power_wave_balance_gate"
    assert gate["status"] == "ok"
    assert gate["rows"][0]["reflected_power_fraction"] == pytest.approx(0.01)
    assert gate["rows"][0]["transmitted_power_fraction"] == pytest.approx(0.64)
    assert gate["rows"][0]["absorbed_power_fraction"] == pytest.approx(0.35)
    assert gate["checks"]["power_balance_basis_matches_expected"] is True
    assert gate["checks"]["column_power_not_active"] is True
    assert gate["checks"]["sparameter_passivity_ok"] is True
    assert gate["checks"]["export_artifact_id_recorded"] is True
    assert gate["checks"]["expected_export_artifact_id_matches"] is True
    assert gate["checks"]["result_set_id_recorded"] is True
    assert gate["checks"]["expected_result_set_id_matches"] is True

    stale_export = touchstone_power_wave_balance_gate(
        [0.1, 0.0],
        [0.8, 0.0],
        s12=[0.8, 0.0],
        s22=[0.1, 0.0],
        data_format="MA",
        power_balance_basis="power_waves_unit_incident_port",
        export_artifact_id="cst_touchstone_export_A",
        expected_export_artifact_id="cst_touchstone_export_B",
        result_set_id="cst_resultset_rf_001",
        expected_result_set_id="cst_resultset_rf_001",
    )
    assert stale_export["status"] == "needs_attention"
    assert stale_export["checks"]["expected_export_artifact_id_matches"] is False
    assert stale_export["checks"]["expected_result_set_id_matches"] is True

    missing_basis = touchstone_power_wave_balance_gate(
        [0.1, 0.0],
        [0.8, 0.0],
        data_format="MA",
        power_balance_basis="",
    )
    assert missing_basis["status"] == "needs_attention"
    assert missing_basis["checks"]["power_balance_basis_recorded"] is False

    wrong_basis = touchstone_power_wave_balance_gate(
        [0.1, 0.0],
        [0.8, 0.0],
        data_format="MA",
        power_balance_basis="accepted_power_farfield",
    )
    assert wrong_basis["status"] == "needs_attention"
    assert wrong_basis["checks"]["power_balance_basis_matches_expected"] is False

    active = touchstone_power_wave_balance_gate(
        [0.0, 0.0],
        [1.1, 0.0],
        s12=[1.1, 0.0],
        s22=[0.0, 0.0],
        data_format="MA",
    )
    assert active["status"] == "needs_attention"
    assert active["checks"]["column_power_not_active"] is False
    assert active["checks"]["sparameter_passivity_ok"] is False


def test_one_port_match_quality_gate_round_trips_vswr_return_loss_and_mismatch_loss():
    gate = one_port_match_quality_gate(s11=1.0 / 3.0)

    assert gate["policy"] == "one_port_match_quality_gate"
    assert gate["status"] == "ok"
    assert gate["reflection_coefficient"] == pytest.approx(1.0 / 3.0)
    assert gate["vswr"] == pytest.approx(2.0)
    assert gate["return_loss_db"] == pytest.approx(9.542425094393248)
    assert gate["transmitted_power_fraction"] == pytest.approx(8.0 / 9.0)
    assert gate["mismatch_loss_db"] == pytest.approx(-10.0 * math.log10(8.0 / 9.0))
    assert gate["checks"]["vswr_round_trip_ok"] is True
    assert gate["checks"]["return_loss_round_trip_ok"] is True

    from_vswr = one_port_match_quality_gate(vswr=2.0)
    from_return_loss = one_port_match_quality_gate(return_loss_db=9.542425094393248)
    assert from_vswr["reflection_coefficient"] == pytest.approx(gate["reflection_coefficient"])
    assert from_return_loss["reflection_coefficient"] == pytest.approx(gate["reflection_coefficient"])

    matched = one_port_match_quality_gate(s11=0.0)
    assert matched["vswr"] == pytest.approx(1.0)
    assert math.isinf(matched["return_loss_db"])
    assert matched["mismatch_loss_db"] == pytest.approx(0.0)

    active = one_port_match_quality_gate(s11=1.02)
    assert active["status"] == "needs_attention"
    assert active["checks"]["passive_reflection_ok"] is False
    assert active["transmitted_power_fraction"] < 0.0


def test_touchstone_sparameter_to_complex_requires_explicit_format_before_health_gate():
    s21_db = touchstone_sparameter_to_complex([-3.0, 90.0], fmt="DB")
    assert abs(s21_db) == pytest.approx(10.0 ** (-3.0 / 20.0))
    assert s21_db.real == pytest.approx(0.0, abs=1.0e-15)
    assert s21_db.imag == pytest.approx(10.0 ** (-3.0 / 20.0))

    s21_ma = touchstone_sparameter_to_complex([10.0 ** (-3.0 / 20.0), 90.0], fmt="MA")
    assert s21_ma == pytest.approx(s21_db)

    health = two_port_sparameter_health(0.0, s21_db, s12=s21_db, s22=0.0)
    assert health["status"] == "ok"
    assert health["passive"] is True

    misread_db_as_ma = two_port_sparameter_health(
        0.0,
        touchstone_sparameter_to_complex([-3.0, 90.0], fmt="MA"),
        s12=touchstone_sparameter_to_complex([-3.0, 90.0], fmt="MA"),
        s22=0.0,
    )
    assert misread_db_as_ma["status"] == "needs_attention"
    assert misread_db_as_ma["passive"] is False


def test_touchstone_port_metadata_gate_freezes_ports_units_and_z0_before_rows():
    metadata = {
        "ports": ["P1", "P2"],
        "port_count": 2,
        "network_parameter": "S",
        "data_format": "MA",
        "frequency_unit": "GHz",
        "reference_impedance_ohm": 50.0,
        "reference_plane": "cst_port1_port2_deembedded_to_connector_faces",
        "port_mode_basis": "single_ended_power_wave_modes",
    }
    gate = touchstone_port_metadata_gate(
        metadata,
        required_ports=("P1", "P2"),
        data_format="MA",
        frequency_unit="GHz",
        reference_impedance_ohm=50.0,
        port_order=("P1", "P2"),
        reference_plane="cst_port1_port2_deembedded_to_connector_faces",
        port_mode_basis="single_ended_power_wave_modes",
    )

    assert gate["policy"] == "touchstone_port_metadata_gate"
    assert gate["status"] == "ok"
    assert gate["ports"] == ["P1", "P2"]
    assert gate["checks"]["port_order_matches_expected"] is True
    assert gate["checks"]["touchstone_format_matches_expected"] is True
    assert gate["checks"]["reference_impedance_matches_expected"] is True
    assert gate["checks"]["reference_plane_recorded"] is True
    assert gate["checks"]["reference_plane_matches_expected"] is True
    assert gate["checks"]["port_mode_basis_recorded"] is True
    assert gate["checks"]["port_mode_basis_matches_expected"] is True

    swapped = touchstone_port_metadata_gate(
        {**metadata, "ports": ["P2", "P1"]},
        required_ports=("P1", "P2"),
        data_format="MA",
        frequency_unit="GHz",
        reference_impedance_ohm=50.0,
        port_order=("P1", "P2"),
    )
    assert swapped["status"] == "needs_attention"
    assert swapped["checks"]["port_order_matches_expected"] is False

    missing_z0 = touchstone_port_metadata_gate(
        {key: value for key, value in metadata.items() if key != "reference_impedance_ohm"},
        data_format="MA",
        frequency_unit="GHz",
        reference_impedance_ohm=50.0,
    )
    assert missing_z0["status"] == "needs_attention"
    assert missing_z0["checks"]["reference_impedance_recorded"] is False

    bad_format = touchstone_port_metadata_gate({**metadata, "data_format": "DB"}, data_format="MA")
    assert bad_format["status"] == "needs_attention"
    assert bad_format["checks"]["touchstone_format_matches_expected"] is False

    missing_reference_plane = touchstone_port_metadata_gate(
        {key: value for key, value in metadata.items() if key != "reference_plane"},
        data_format="MA",
        frequency_unit="GHz",
        reference_impedance_ohm=50.0,
        reference_plane="cst_port1_port2_deembedded_to_connector_faces",
    )
    assert missing_reference_plane["status"] == "needs_attention"
    assert missing_reference_plane["checks"]["reference_plane_recorded"] is False

    wrong_mode_basis = touchstone_port_metadata_gate(
        {**metadata, "port_mode_basis": "mixed_mode_differential_common"},
        data_format="MA",
        frequency_unit="GHz",
        reference_impedance_ohm=50.0,
        port_mode_basis="single_ended_power_wave_modes",
    )
    assert wrong_mode_basis["status"] == "needs_attention"
    assert wrong_mode_basis["checks"]["port_mode_basis_recorded"] is True
    assert wrong_mode_basis["checks"]["port_mode_basis_matches_expected"] is False


def test_farfield_pattern_metadata_gate_freezes_units_cuts_and_polarization_before_rows():
    metadata = {
        "frequency_hz": 2.45e9,
        "angle_unit": "deg",
        "theta_values_deg": [0.0, 90.0, 180.0],
        "phi_values_deg": [0.0, 90.0],
        "coordinate_system": "spherical",
        "polarization_basis": "theta_phi",
        "quantity": "gain",
        "quantity_unit": "dBi",
        "normalization": "accepted_power",
        "field_components": ["Etheta", "Ephi"],
        "row_count": 6,
    }
    gate = farfield_pattern_metadata_gate(
        metadata,
        frequency_hz=2.45e9,
        required_phi_values_deg=(0.0, 90.0),
    )

    assert gate["policy"] == "farfield_pattern_metadata_gate"
    assert gate["status"] == "ok"
    assert gate["theta_span_deg"] == pytest.approx(180.0)
    assert gate["expected_grid_rows"] == 6
    assert gate["checks"]["required_components_present"] is True
    assert gate["checks"]["normalization_matches_expected"] is True
    assert gate["checks"]["required_phi_values_present"] is True

    bad_unit = farfield_pattern_metadata_gate({**metadata, "angle_unit": "rad"})
    assert bad_unit["status"] == "needs_attention"
    assert bad_unit["checks"]["angle_unit_matches_expected"] is False

    missing_component = farfield_pattern_metadata_gate({**metadata, "field_components": ["Etheta"]})
    assert missing_component["status"] == "needs_attention"
    assert missing_component["checks"]["required_components_present"] is False

    narrow_cut = farfield_pattern_metadata_gate({**metadata, "theta_values_deg": [0.0, 60.0, 120.0], "row_count": 6})
    assert narrow_cut["status"] == "needs_attention"
    assert narrow_cut["checks"]["theta_span_ok"] is False

    missing_norm = farfield_pattern_metadata_gate({key: value for key, value in metadata.items() if key != "normalization"})
    assert missing_norm["status"] == "needs_attention"
    assert missing_norm["checks"]["normalization_matches_expected"] is False


def test_farfield_lobe_notebook_handoff_gate_bundles_metadata_and_gain_row():
    metadata = {
        "frequency_hz": 2.45e9,
        "angle_unit": "deg",
        "theta_values_deg": [0.0, 90.0, 180.0],
        "phi_values_deg": [0.0, 90.0],
        "coordinate_system": "spherical",
        "polarization_basis": "theta_phi",
        "quantity": "gain",
        "quantity_unit": "dBi",
        "normalization": "accepted_power",
        "field_components": ["Etheta", "Ephi"],
        "row_count": 6,
    }
    directivity_dbi = 7.0
    eta = 0.65
    gain_dbi = 10.0 * math.log10((10.0 ** (directivity_dbi / 10.0)) * eta)
    row = {
        "lobe_id": "main",
        "frequency_hz": 2.45e9,
        "theta_deg": 90.0,
        "phi_deg": 0.0,
        "polarization_basis": "theta_phi",
        "normalization": "accepted_power",
        "gain_unit": "dBi",
        "directivity_unit": "dBi",
        "gain_dbi": gain_dbi,
        "directivity_dbi": directivity_dbi,
        "radiated_power_w": 6.5,
        "accepted_power_w": 10.0,
    }
    gate = farfield_lobe_notebook_handoff_gate(metadata, row)

    assert gate["policy"] == "farfield_lobe_notebook_handoff_gate"
    assert gate["status"] == "ok"
    assert gate["metadata_gate"]["status"] == "ok"
    assert gate["lobe_id"] == "main"
    assert gate["radiation_efficiency"] == pytest.approx(eta)
    assert gate["gain_relative_error"] < 1.0e-12
    assert all(gate["checks"].values())

    missing_lobe = farfield_lobe_notebook_handoff_gate(
        metadata,
        {key: value for key, value in row.items() if key != "lobe_id"},
    )
    assert missing_lobe["status"] == "needs_attention"
    assert missing_lobe["checks"]["lobe_id_recorded"] is False

    wrong_cut = farfield_lobe_notebook_handoff_gate(metadata, {**row, "phi_deg": 45.0})
    assert wrong_cut["status"] == "needs_attention"
    assert wrong_cut["checks"]["phi_on_export_grid"] is False

    too_high_gain = farfield_lobe_notebook_handoff_gate(
        metadata,
        {**row, "gain_dbi": directivity_dbi + 0.2},
    )
    assert too_high_gain["status"] == "needs_attention"
    assert too_high_gain["checks"]["gain_not_above_directivity"] is False


def test_cst_export_manifest_solver_ready_gate_freezes_identity_grid_and_files():
    manifest = {
        "project_id": "rf_widget_v1",
        "run_id": "run_2p45g_001",
        "export_id": "export_A",
        "source_tool": "CST Studio Suite",
        "solver_kind": "frequency_domain",
        "frequency_unit": "Hz",
        "frequency_grid_Hz": [2.40e9, 2.45e9, 2.50e9],
        "design_frequency_Hz": 2.45e9,
        "files": [
            {
                "kind": "touchstone",
                "path": "slot159_ports.s2p",
                "status": "exported",
                "data_format": "MA",
                "z0_ohm": 50.0,
                "port_order": ["in", "out"],
            },
            {
                "kind": "farfield",
                "path": "slot159_farfield.ffm",
                "status": "exported",
                "normalization": "accepted_power",
                "angle_unit": "deg",
                "polarization_basis": "theta_phi",
            },
        ],
    }

    gate = cst_export_manifest_solver_ready_gate(
        manifest,
        expected_project_id="rf_widget_v1",
        expected_run_id="run_2p45g_001",
        expected_export_id="export_A",
        expected_solver_kind="frequency_domain",
        expected_design_frequency_Hz=2.45e9,
    )

    assert gate["status"] == "ok"
    assert gate["policy"] == "cst_export_manifest_solver_ready_gate"
    assert gate["checks"]["frequency_grid_strictly_increasing"] is True
    assert gate["checks"]["design_frequency_bracketed"] is True
    assert gate["checks"]["touchstone_metadata_present"] is True
    assert gate["checks"]["farfield_metadata_present"] is True

    unbracketed = {**manifest, "frequency_grid_Hz": [2.30e9, 2.35e9, 2.40e9]}
    unbracketed_gate = cst_export_manifest_solver_ready_gate(unbracketed)
    assert unbracketed_gate["status"] == "needs_attention"
    assert unbracketed_gate["checks"]["design_frequency_bracketed"] is False

    missing_touchstone_metadata = {**manifest, "files": [dict(row) for row in manifest["files"]]}
    missing_touchstone_metadata["files"][0].pop("data_format")
    missing_touchstone_metadata_gate = cst_export_manifest_solver_ready_gate(missing_touchstone_metadata)
    assert missing_touchstone_metadata_gate["status"] == "needs_attention"
    assert missing_touchstone_metadata_gate["checks"]["touchstone_metadata_present"] is False

    wrong_source = {**manifest, "source_tool": "HFSS"}
    wrong_source_gate = cst_export_manifest_solver_ready_gate(wrong_source)
    assert wrong_source_gate["status"] == "needs_attention"
    assert wrong_source_gate["checks"]["source_tool_is_cst"] is False


def test_cst_result_export_package_gate_keeps_touchstone_and_farfield_rows_together():
    artifacts = [
        {
            "kind": "touchstone_metadata",
            "project_id": "rf_widget_v1",
            "run_id": "run_2p45g_001",
            "export_id": "export_A",
            "frequency_Hz": 2.45e9,
            "source_tool": "CST Studio Suite",
            "path": "slot151_ports.s2p",
            "gate_policy": "touchstone_port_metadata_gate",
            "status": "ok",
        },
        {
            "kind": "touchstone_row",
            "project_id": "rf_widget_v1",
            "run_id": "run_2p45g_001",
            "export_id": "export_A",
            "frequency_Hz": 2.45e9,
            "source_tool": "CST",
            "path": "slot151_ports_row.json",
            "gate_policy": "touchstone_row_solver_ready_preflight_gate",
            "status": "ok",
        },
        {
            "kind": "farfield_metadata",
            "project_id": "rf_widget_v1",
            "run_id": "run_2p45g_001",
            "export_id": "export_A",
            "frequency_Hz": 2.45e9,
            "source_tool": "CST",
            "path": "slot151_farfield_metadata.json",
            "gate_policy": "farfield_pattern_metadata_gate",
            "status": "ok",
        },
        {
            "kind": "farfield_lobe",
            "project_id": "rf_widget_v1",
            "run_id": "run_2p45g_001",
            "export_id": "export_A",
            "frequency_Hz": 2.45e9,
            "source_tool": "CST",
            "path": "slot151_farfield_lobe.json",
            "gate_policy": "farfield_lobe_notebook_handoff_gate",
            "status": "ok",
        },
    ]

    gate = cst_result_export_package_gate(
        artifacts,
        expected_project_id="rf_widget_v1",
        expected_run_id="run_2p45g_001",
        expected_export_id="export_A",
        expected_frequency_Hz=2.45e9,
    )

    assert gate["status"] == "ok"
    assert gate["policy"] == "cst_result_export_package_gate"
    assert gate["project_ids"] == ["rf_widget_v1"]
    assert gate["run_ids"] == ["run_2p45g_001"]
    assert gate["export_ids"] == ["export_A"]
    assert gate["checks"]["frequencies_match"] is True

    stale_export = [dict(row) for row in artifacts]
    stale_export[3]["export_id"] = "export_old"
    stale_export_gate = cst_result_export_package_gate(stale_export)
    assert stale_export_gate["status"] == "needs_attention"
    assert stale_export_gate["checks"]["export_ids_unique"] is False

    wrong_frequency = [dict(row) for row in artifacts]
    wrong_frequency[2]["frequency_Hz"] = 2.40e9
    wrong_frequency_gate = cst_result_export_package_gate(wrong_frequency)
    assert wrong_frequency_gate["status"] == "needs_attention"
    assert wrong_frequency_gate["checks"]["frequencies_match"] is False

    wrong_source = [dict(row) for row in artifacts]
    wrong_source[1]["source_tool"] = "JMAG"
    wrong_source_gate = cst_result_export_package_gate(wrong_source)
    assert wrong_source_gate["status"] == "needs_attention"
    assert wrong_source_gate["checks"]["source_tool_is_cst"] is False


def test_cst_touchstone_solver_ready_manifest_gate_keeps_grid_metadata_and_row_together():
    artifacts = [
        {
            "kind": "port_metadata",
            "project_id": "rf_filter_v2",
            "run_id": "run_s2p_001",
            "export_id": "touchstone_export_B",
            "touchstone_export_method": "cst_result_tree_touchstone_s2p_export",
            "source_tool": "CST Studio Suite",
            "path": "slot167_port_metadata.json",
            "gate_policy": "touchstone_port_metadata_gate",
            "status": "ok",
            "pass": True,
            "design_frequency_hz": 2.45e9,
            "network_kind": "S",
            "data_format": "MA",
            "reference_impedance_ohm": 50.0,
            "port_order": ["P1", "P2"],
            "reference_plane": "cst_port1_port2_deembedded_to_connector_faces",
            "reference_plane_geometry_digest": "sha256:cst_slot327_connector_face_plane_geometry_v1",
            "port_face_centers_xyz_m": [[0.0, 0.0, 0.0], [0.05, 0.0, 0.0]],
            "port_mode_basis": "single_ended_power_wave_modes",
            "incident_wave_convention": "unit_incident_power_wave_per_excited_port",
            "power_balance_basis": "power_waves_unit_incident_port",
            "frequency_grid_id": "grid_sweep_2p45g_5pts_v2",
            "interpolation_policy": "exact_design_row_no_interpolation",
            "touchstone_file_id": "touchstone_export_B_s2p_sha256_abc123",
            "touchstone_observable_id": "cst_slot303_s2p_single_ended_sparameter_v1",
            "touchstone_observable_family": "single_ended_sparameter",
            "touchstone_output_artifact_id": "touchstone_export_B_row_002_postprocess_v1",
            "touchstone_output_digest": "sha256:touchstone_export_B_row_002_postprocess_v1",
            "touchstone_output_path": "slot295_touchstone_row_002_postprocess.json",
            "renormalized_reference_impedance_ohm": 50.0,
            "renormalization_method": "not_renormalized_option_line_R50",
            "renormalization_artifact_id": "touchstone_export_B_renorm_none_R50_v1",
            "deembedding_method": "port_extension_reference_plane_shift",
            "deembedding_artifact_id": "touchstone_export_B_deembed_connector_faces_v1",
            "deembedding_length_m": 0.0125,
        },
        {
            "kind": "frequency_grid",
            "project_id": "rf_filter_v2",
            "run_id": "run_s2p_001",
            "export_id": "touchstone_export_B",
            "touchstone_export_method": "cst_result_tree_touchstone_s2p_export",
            "source_tool": "CST",
            "path": "slot167_frequency_grid.json",
            "gate_policy": "touchstone_frequency_grid_interpolation_gate",
            "status": "ok",
            "pass": True,
            "design_frequency_hz": 2.45e9,
            "network_kind": "S",
            "data_format": "MA",
            "reference_impedance_ohm": 50.0,
            "port_order": ["P1", "P2"],
            "reference_plane": "cst_port1_port2_deembedded_to_connector_faces",
            "reference_plane_geometry_digest": "sha256:cst_slot327_connector_face_plane_geometry_v1",
            "port_face_centers_xyz_m": [[0.0, 0.0, 0.0], [0.05, 0.0, 0.0]],
            "port_mode_basis": "single_ended_power_wave_modes",
            "incident_wave_convention": "unit_incident_power_wave_per_excited_port",
            "power_balance_basis": "power_waves_unit_incident_port",
            "frequency_grid_id": "grid_sweep_2p45g_5pts_v2",
            "interpolation_policy": "exact_design_row_no_interpolation",
            "touchstone_file_id": "touchstone_export_B_s2p_sha256_abc123",
            "touchstone_observable_id": "cst_slot303_s2p_single_ended_sparameter_v1",
            "touchstone_observable_family": "single_ended_sparameter",
            "touchstone_output_artifact_id": "touchstone_export_B_row_002_postprocess_v1",
            "touchstone_output_digest": "sha256:touchstone_export_B_row_002_postprocess_v1",
            "touchstone_output_path": "slot295_touchstone_row_002_postprocess.json",
            "renormalized_reference_impedance_ohm": 50.0,
            "renormalization_method": "not_renormalized_option_line_R50",
            "renormalization_artifact_id": "touchstone_export_B_renorm_none_R50_v1",
            "deembedding_method": "port_extension_reference_plane_shift",
            "deembedding_artifact_id": "touchstone_export_B_deembed_connector_faces_v1",
            "deembedding_length_m": 0.0125,
            "design_frequency_bracketed": True,
            "row_count": 5,
        },
        {
            "kind": "design_row",
            "project_id": "rf_filter_v2",
            "run_id": "run_s2p_001",
            "export_id": "touchstone_export_B",
            "touchstone_export_method": "cst_result_tree_touchstone_s2p_export",
            "source_tool": "CST",
            "path": "slot167_design_row.json",
            "gate_policy": "touchstone_row_solver_ready_preflight_gate",
            "status": "ok",
            "pass": True,
            "design_frequency_hz": 2.45e9,
            "network_kind": "S",
            "data_format": "MA",
            "reference_impedance_ohm": 50.0,
            "port_order": ["P1", "P2"],
            "reference_plane": "cst_port1_port2_deembedded_to_connector_faces",
            "reference_plane_geometry_digest": "sha256:cst_slot327_connector_face_plane_geometry_v1",
            "port_face_centers_xyz_m": [[0.0, 0.0, 0.0], [0.05, 0.0, 0.0]],
            "port_mode_basis": "single_ended_power_wave_modes",
            "incident_wave_convention": "unit_incident_power_wave_per_excited_port",
            "power_balance_basis": "power_waves_unit_incident_port",
            "frequency_grid_id": "grid_sweep_2p45g_5pts_v2",
            "interpolation_policy": "exact_design_row_no_interpolation",
            "touchstone_file_id": "touchstone_export_B_s2p_sha256_abc123",
            "touchstone_observable_id": "cst_slot303_s2p_single_ended_sparameter_v1",
            "touchstone_observable_family": "single_ended_sparameter",
            "touchstone_output_artifact_id": "touchstone_export_B_row_002_postprocess_v1",
            "touchstone_output_digest": "sha256:touchstone_export_B_row_002_postprocess_v1",
            "touchstone_output_path": "slot295_touchstone_row_002_postprocess.json",
            "renormalized_reference_impedance_ohm": 50.0,
            "renormalization_method": "not_renormalized_option_line_R50",
            "renormalization_artifact_id": "touchstone_export_B_renorm_none_R50_v1",
            "deembedding_method": "port_extension_reference_plane_shift",
            "deembedding_artifact_id": "touchstone_export_B_deembed_connector_faces_v1",
            "deembedding_length_m": 0.0125,
            "sparameter_passivity_ok": True,
            "sparameter_reciprocity_ok": True,
            "selected_row_index": 2,
            "selected_frequency_hz": 2.45e9,
        },
        {
            "kind": "touchstone_row",
            "project_id": "rf_filter_v2",
            "run_id": "run_s2p_001",
            "export_id": "touchstone_export_B",
            "touchstone_export_method": "cst_result_tree_touchstone_s2p_export",
            "source_tool": "CST",
            "path": "slot167_touchstone_row_002.json",
            "gate_policy": "touchstone_row_solver_ready_preflight_gate",
            "status": "ok",
            "pass": True,
            "design_frequency_hz": 2.45e9,
            "network_kind": "S",
            "data_format": "MA",
            "reference_impedance_ohm": 50.0,
            "port_order": ["P1", "P2"],
            "reference_plane": "cst_port1_port2_deembedded_to_connector_faces",
            "reference_plane_geometry_digest": "sha256:cst_slot327_connector_face_plane_geometry_v1",
            "port_face_centers_xyz_m": [[0.0, 0.0, 0.0], [0.05, 0.0, 0.0]],
            "port_mode_basis": "single_ended_power_wave_modes",
            "incident_wave_convention": "unit_incident_power_wave_per_excited_port",
            "power_balance_basis": "power_waves_unit_incident_port",
            "frequency_grid_id": "grid_sweep_2p45g_5pts_v2",
            "interpolation_policy": "exact_design_row_no_interpolation",
            "touchstone_file_id": "touchstone_export_B_s2p_sha256_abc123",
            "touchstone_observable_id": "cst_slot303_s2p_single_ended_sparameter_v1",
            "touchstone_observable_family": "single_ended_sparameter",
            "touchstone_output_artifact_id": "touchstone_export_B_row_002_postprocess_v1",
            "touchstone_output_digest": "sha256:touchstone_export_B_row_002_postprocess_v1",
            "touchstone_output_path": "slot295_touchstone_row_002_postprocess.json",
            "renormalized_reference_impedance_ohm": 50.0,
            "renormalization_method": "not_renormalized_option_line_R50",
            "renormalization_artifact_id": "touchstone_export_B_renorm_none_R50_v1",
            "deembedding_method": "port_extension_reference_plane_shift",
            "deembedding_artifact_id": "touchstone_export_B_deembed_connector_faces_v1",
            "deembedding_length_m": 0.0125,
            "row_index": 2,
            "row_frequency_hz": 2.45e9,
            "sparameter_passivity_ok": True,
            "sparameter_reciprocity_ok": True,
        },
    ]
    for row in artifacts:
        row["solver_setup_artifact_id"] = "cst_slot335_fd_solver_setup_v1"
        row["mesh_setup_artifact_id"] = "cst_slot335_adaptive_mesh_setup_v1"
        row["port_definition_artifact_id"] = "cst_slot351_waveguide_port_definition_v1"
        row["excitation_setup_artifact_id"] = "cst_slot351_unit_power_port_excitation_v1"
        row["result_tree_path"] = "1D Results\\S-Parameters\\Touchstone\\S2P Export"
        row["result_item_id"] = "cst_slot343_result_tree_sparams_s2p_v1"
        row["frequency_grid_digest"] = "sha256:grid-sweep-2p45g-5pts-v2"
        row["touchstone_option_line_artifact_id"] = "touchstone_export_B_option_line_v1"
        row["touchstone_option_line_digest"] = (
            "sha256:touchstone-export-B-option-line-ghz-s-ma-r50-v1"
        )
        row["touchstone_port_mode_basis_schema_id"] = (
            "cst_single_ended_power_wave_port_mode_basis_v1"
        )
        row["model_input_artifact_id"] = "cst_slot380_rf_filter_project_v1.cst"
        row["model_input_digest"] = "sha256:cst-slot380-rf-filter-project-v1"
        row["model_input_path"] = "artifacts/rf/cst_slot380_rf_filter_project_v1.cst"
        row["export_recipe_artifact_id"] = "cst_slot387_touchstone_export_macro_v1.bas"
        row["export_recipe_digest"] = "sha256:cst-slot387-touchstone-export-macro-v1"
        row["export_recipe_path"] = (
            "artifacts/rf/cst_slot387_touchstone_export_macro_v1.bas"
        )
        row["touchstone_output_schema_id"] = "cst_touchstone_s2p_row_table_v1"
        row["touchstone_convention_schema_id"] = "cst_touchstone_network_convention_v1"
        row["touchstone_postprocess_row_convention_schema_id"] = (
            "cst_touchstone_s2p_row_convention_v1"
        )
        row["touchstone_output_columns"] = [
            "frequency_hz",
            "S11",
            "S21",
            "S12",
            "S22",
        ]
        row["touchstone_output_units"] = {
            "frequency_hz": "Hz",
            "S11": "1",
            "S21": "1",
            "S12": "1",
            "S22": "1",
        }
        row["parameter_set_artifact_id"] = "cst_slot394_rf_filter_parameter_set_v1.json"
        row["parameter_set_digest"] = "sha256:cst-slot394-rf-filter-parameter-set-v1"
        row["parameter_set_path"] = (
            "artifacts/rf/cst_slot394_rf_filter_parameter_set.json"
        )
        row["objective_observable_id"] = "cst_slot394_s21_insertion_loss_objective_v1"
        row["objective_observable_family"] = "touchstone_s21_insertion_loss_objective"
        row["created_at_utc"] = "2026-07-01T14:05:20Z"
        row["run_timestamp_utc"] = "2026-07-01T14:05:00Z"

    gate = cst_touchstone_solver_ready_manifest_gate(
        artifacts,
        expected_project_id="rf_filter_v2",
        expected_run_id="run_s2p_001",
        expected_export_id="touchstone_export_B",
        expected_model_input_artifact_id="cst_slot380_rf_filter_project_v1.cst",
        expected_model_input_digest="sha256:cst-slot380-rf-filter-project-v1",
        expected_model_input_path="artifacts/rf/cst_slot380_rf_filter_project_v1.cst",
        expected_parameter_set_artifact_id="cst_slot394_rf_filter_parameter_set_v1.json",
        expected_parameter_set_digest="sha256:cst-slot394-rf-filter-parameter-set-v1",
        expected_parameter_set_path=(
            "artifacts/rf/cst_slot394_rf_filter_parameter_set.json"
        ),
        expected_objective_observable_id="cst_slot394_s21_insertion_loss_objective_v1",
        expected_objective_observable_family="touchstone_s21_insertion_loss_objective",
        expected_design_frequency_hz=2.45e9,
        expected_network_kind="S",
        expected_port_order=("P1", "P2"),
        expected_data_format="MA",
        expected_reference_impedance_ohm=50.0,
        expected_touchstone_option_line_artifact_id="touchstone_export_B_option_line_v1",
        expected_touchstone_option_line_digest=(
            "sha256:touchstone-export-B-option-line-ghz-s-ma-r50-v1"
        ),
        expected_reference_plane="cst_port1_port2_deembedded_to_connector_faces",
        expected_reference_plane_geometry_digest="sha256:cst_slot327_connector_face_plane_geometry_v1",
        expected_port_face_centers_xyz_m=((0.0, 0.0, 0.0), (0.05, 0.0, 0.0)),
        expected_port_mode_basis="single_ended_power_wave_modes",
        expected_touchstone_port_mode_basis_schema_id=(
            "cst_single_ended_power_wave_port_mode_basis_v1"
        ),
        expected_incident_wave_convention="unit_incident_power_wave_per_excited_port",
        expected_power_balance_basis="power_waves_unit_incident_port",
        expected_touchstone_export_method="cst_result_tree_touchstone_s2p_export",
        expected_export_recipe_artifact_id="cst_slot387_touchstone_export_macro_v1.bas",
        expected_export_recipe_digest="sha256:cst-slot387-touchstone-export-macro-v1",
        expected_export_recipe_path=(
            "artifacts/rf/cst_slot387_touchstone_export_macro_v1.bas"
        ),
        expected_result_tree_path="1D Results\\S-Parameters\\Touchstone\\S2P Export",
        expected_result_item_id="cst_slot343_result_tree_sparams_s2p_v1",
        expected_solver_setup_artifact_id="cst_slot335_fd_solver_setup_v1",
        expected_mesh_setup_artifact_id="cst_slot335_adaptive_mesh_setup_v1",
        expected_port_definition_artifact_id="cst_slot351_waveguide_port_definition_v1",
        expected_excitation_setup_artifact_id="cst_slot351_unit_power_port_excitation_v1",
        expected_frequency_grid_id="grid_sweep_2p45g_5pts_v2",
        expected_frequency_grid_digest="sha256:grid-sweep-2p45g-5pts-v2",
        expected_interpolation_policy="exact_design_row_no_interpolation",
        expected_touchstone_file_id="touchstone_export_B_s2p_sha256_abc123",
        expected_touchstone_observable_id="cst_slot303_s2p_single_ended_sparameter_v1",
        expected_touchstone_observable_family="single_ended_sparameter",
        expected_touchstone_output_artifact_id="touchstone_export_B_row_002_postprocess_v1",
        expected_touchstone_output_digest="sha256:touchstone_export_B_row_002_postprocess_v1",
        expected_touchstone_output_schema_id="cst_touchstone_s2p_row_table_v1",
        expected_touchstone_output_columns=[
            "frequency_hz",
            "S11",
            "S21",
            "S12",
            "S22",
        ],
        expected_touchstone_output_units={
            "frequency_hz": "Hz",
            "S11": "1",
            "S21": "1",
            "S12": "1",
            "S22": "1",
        },
        expected_touchstone_convention_schema_id="cst_touchstone_network_convention_v1",
        expected_touchstone_postprocess_row_convention_schema_id=(
            "cst_touchstone_s2p_row_convention_v1"
        ),
        require_touchstone_output_artifact=True,
        require_touchstone_output_schema=True,
        require_touchstone_convention_schema=True,
        require_touchstone_port_mode_basis_schema=True,
        require_touchstone_postprocess_row_convention_schema=True,
        require_export_recipe_artifact=True,
        expected_created_at_utc="2026-07-01T14:05:20Z",
        expected_run_timestamp_utc="2026-07-01T14:05:00Z",
        max_created_run_skew_s=60.0,
        require_execution_metadata=True,
        expected_renormalized_reference_impedance_ohm=50.0,
        expected_renormalization_method="not_renormalized_option_line_R50",
        expected_renormalization_artifact_id="touchstone_export_B_renorm_none_R50_v1",
        expected_deembedding_method="port_extension_reference_plane_shift",
        expected_deembedding_artifact_id="touchstone_export_B_deembed_connector_faces_v1",
        expected_deembedding_length_m=0.0125,
        require_model_input_artifact=True,
        require_parameter_set_artifact=True,
    )

    assert gate["status"] == "ok"
    assert gate["policy"] == "cst_touchstone_solver_ready_manifest_gate"
    assert gate["checks"]["port_metadata_complete"] is True
    assert gate["checks"]["frequency_grid_brackets_design"] is True
    assert gate["checks"]["design_row_passive"] is True
    assert gate["checks"]["touchstone_row_passive_when_present"] is True
    assert gate["checks"]["frequency_grid_row_count_recorded"] is True
    assert gate["checks"]["design_row_index_recorded"] is True
    assert gate["checks"]["design_row_index_within_grid"] is True
    assert gate["checks"]["touchstone_row_index_recorded_when_present"] is True
    assert gate["checks"]["touchstone_row_frequency_recorded_when_present"] is True
    assert gate["checks"]["touchstone_row_index_matches_design_row"] is True
    assert gate["checks"]["touchstone_row_frequency_matches_selected_frequency"] is True
    assert gate["checks"]["network_kind_unique"] is True
    assert gate["checks"]["expected_network_kind_matches"] is True
    assert gate["checks"]["port_order_unique"] is True
    assert gate["checks"]["expected_port_order_matches"] is True
    assert gate["checks"]["data_format_unique"] is True
    assert gate["checks"]["expected_data_format_matches"] is True
    assert gate["checks"]["reference_impedance_unique"] is True
    assert gate["checks"]["expected_reference_impedance_matches"] is True
    assert gate["checks"]["model_input_artifact_id_consistent_when_present"] is True
    assert gate["checks"]["model_input_digest_consistent_when_present"] is True
    assert gate["checks"]["model_input_path_consistent_when_present"] is True
    assert gate["checks"]["model_input_artifact_id_recorded_when_required"] is True
    assert gate["checks"]["model_input_digest_recorded_when_required"] is True
    assert gate["checks"]["model_input_path_recorded_when_required"] is True
    assert gate["checks"]["expected_model_input_artifact_id_matches"] is True
    assert gate["checks"]["expected_model_input_digest_matches"] is True
    assert gate["checks"]["expected_model_input_path_matches"] is True
    assert gate["model_input_artifact_ids"] == ["cst_slot380_rf_filter_project_v1.cst"]
    assert gate["model_input_digests"] == ["sha256:cst-slot380-rf-filter-project-v1"]
    assert gate["checks"]["parameter_set_artifact_id_consistent_when_present"] is True
    assert gate["checks"]["parameter_set_digest_consistent_when_present"] is True
    assert gate["checks"]["parameter_set_path_consistent_when_present"] is True
    assert gate["checks"]["objective_observable_id_consistent_when_present"] is True
    assert gate["checks"]["objective_observable_family_consistent_when_present"] is True
    assert gate["checks"]["parameter_set_artifact_id_recorded_when_required"] is True
    assert gate["checks"]["parameter_set_digest_recorded_when_required"] is True
    assert gate["checks"]["parameter_set_path_recorded_when_required"] is True
    assert gate["checks"]["expected_parameter_set_artifact_id_matches"] is True
    assert gate["checks"]["expected_parameter_set_digest_matches"] is True
    assert gate["checks"]["expected_parameter_set_path_matches"] is True
    assert gate["checks"]["expected_objective_observable_id_matches"] is True
    assert gate["checks"]["expected_objective_observable_family_matches"] is True
    assert gate["parameter_set_artifact_ids"] == [
        "cst_slot394_rf_filter_parameter_set_v1.json"
    ]
    assert gate["objective_observable_families"] == [
        "touchstone_s21_insertion_loss_objective"
    ]
    assert gate["checks"]["touchstone_option_line_artifact_id_consistent_when_present"] is True
    assert gate["checks"]["touchstone_option_line_digest_consistent_when_present"] is True
    assert gate["checks"]["touchstone_option_line_artifact_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_touchstone_option_line_artifact_id_matches"] is True
    assert gate["checks"]["touchstone_option_line_digest_recorded_when_expected"] is True
    assert gate["checks"]["expected_touchstone_option_line_digest_matches"] is True
    assert gate["touchstone_option_line_artifact_ids"] == [
        "touchstone_export_B_option_line_v1"
    ]
    assert gate["touchstone_option_line_digests"] == [
        "sha256:touchstone-export-B-option-line-ghz-s-ma-r50-v1"
    ]
    assert gate["checks"]["reference_plane_consistent_when_present"] is True
    assert gate["checks"]["reference_plane_recorded_when_expected"] is True
    assert gate["checks"]["expected_reference_plane_matches"] is True
    assert gate["checks"]["reference_plane_geometry_digest_consistent_when_present"] is True
    assert gate["checks"]["reference_plane_geometry_digest_recorded_when_expected"] is True
    assert gate["checks"]["expected_reference_plane_geometry_digest_matches"] is True
    assert gate["checks"]["port_face_centers_xyz_consistent_when_present"] is True
    assert gate["checks"]["port_face_centers_xyz_recorded_when_expected"] is True
    assert gate["checks"]["expected_port_face_centers_xyz_matches"] is True
    assert gate["reference_plane_geometry_digests"] == [
        "sha256:cst_slot327_connector_face_plane_geometry_v1"
    ]
    assert gate["port_face_centers_xyz_m"] == [
        [[0.0, 0.0, 0.0], [0.05, 0.0, 0.0]]
    ]
    assert gate["checks"]["port_mode_basis_consistent_when_present"] is True
    assert gate["checks"]["port_mode_basis_recorded_when_expected"] is True
    assert gate["checks"]["expected_port_mode_basis_matches"] is True
    assert (
        gate["checks"]["touchstone_port_mode_basis_schema_id_consistent_when_present"]
        is True
    )
    assert (
        gate["checks"]["touchstone_port_mode_basis_schema_id_recorded_when_required"]
        is True
    )
    assert (
        gate["checks"]["touchstone_port_mode_basis_schema_id_recorded_when_expected"]
        is True
    )
    assert (
        gate["checks"]["expected_touchstone_port_mode_basis_schema_id_matches"]
        is True
    )
    assert gate["touchstone_port_mode_basis_schema_ids"] == [
        "cst_single_ended_power_wave_port_mode_basis_v1"
    ]
    assert (
        gate["expected_touchstone_port_mode_basis_schema_id"]
        == "cst_single_ended_power_wave_port_mode_basis_v1"
    )
    assert gate["require_touchstone_port_mode_basis_schema"] is True
    assert gate["checks"]["incident_wave_convention_consistent_when_present"] is True
    assert gate["checks"]["incident_wave_convention_recorded_when_expected"] is True
    assert gate["checks"]["expected_incident_wave_convention_matches"] is True
    assert gate["checks"]["power_balance_basis_consistent_when_present"] is True
    assert gate["checks"]["power_balance_basis_recorded_when_expected"] is True
    assert gate["checks"]["expected_power_balance_basis_matches"] is True
    assert gate["checks"]["touchstone_export_method_consistent_when_present"] is True
    assert gate["checks"]["touchstone_export_method_recorded_when_expected"] is True
    assert gate["checks"]["expected_touchstone_export_method_matches"] is True
    assert gate["export_recipe_artifact_ids"] == [
        "cst_slot387_touchstone_export_macro_v1.bas"
    ]
    assert gate["export_recipe_digests"] == [
        "sha256:cst-slot387-touchstone-export-macro-v1"
    ]
    assert gate["export_recipe_paths"] == [
        "artifacts/rf/cst_slot387_touchstone_export_macro_v1.bas"
    ]
    assert gate["checks"]["export_recipe_artifact_id_consistent_when_present"] is True
    assert gate["checks"]["export_recipe_digest_consistent_when_present"] is True
    assert gate["checks"]["export_recipe_path_consistent_when_present"] is True
    assert gate["checks"]["export_recipe_artifact_id_recorded_when_required"] is True
    assert gate["checks"]["export_recipe_digest_recorded_when_required"] is True
    assert gate["checks"]["export_recipe_path_recorded_when_required"] is True
    assert gate["checks"]["export_recipe_artifact_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_export_recipe_artifact_id_matches"] is True
    assert gate["checks"]["export_recipe_digest_recorded_when_expected"] is True
    assert gate["checks"]["expected_export_recipe_digest_matches"] is True
    assert gate["checks"]["export_recipe_path_recorded_when_expected"] is True
    assert gate["checks"]["expected_export_recipe_path_matches"] is True
    assert gate["result_tree_paths"] == ["1D Results\\S-Parameters\\Touchstone\\S2P Export"]
    assert gate["result_item_ids"] == ["cst_slot343_result_tree_sparams_s2p_v1"]
    assert gate["checks"]["result_tree_path_consistent_when_present"] is True
    assert gate["checks"]["result_item_id_consistent_when_present"] is True
    assert gate["checks"]["result_tree_path_recorded_when_expected"] is True
    assert gate["checks"]["result_item_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_result_tree_path_matches"] is True
    assert gate["checks"]["expected_result_item_id_matches"] is True
    assert gate["solver_setup_artifact_ids"] == ["cst_slot335_fd_solver_setup_v1"]
    assert gate["mesh_setup_artifact_ids"] == ["cst_slot335_adaptive_mesh_setup_v1"]
    assert gate["port_definition_artifact_ids"] == ["cst_slot351_waveguide_port_definition_v1"]
    assert gate["excitation_setup_artifact_ids"] == ["cst_slot351_unit_power_port_excitation_v1"]
    assert gate["checks"]["solver_setup_artifact_id_consistent_when_present"] is True
    assert gate["checks"]["solver_setup_artifact_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_solver_setup_artifact_id_matches"] is True
    assert gate["checks"]["mesh_setup_artifact_id_consistent_when_present"] is True
    assert gate["checks"]["mesh_setup_artifact_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_mesh_setup_artifact_id_matches"] is True
    assert gate["checks"]["port_definition_artifact_id_consistent_when_present"] is True
    assert gate["checks"]["port_definition_artifact_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_port_definition_artifact_id_matches"] is True
    assert gate["checks"]["excitation_setup_artifact_id_consistent_when_present"] is True
    assert gate["checks"]["excitation_setup_artifact_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_excitation_setup_artifact_id_matches"] is True
    assert gate["checks"]["frequency_grid_id_consistent_when_present"] is True
    assert gate["checks"]["frequency_grid_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_frequency_grid_id_matches"] is True
    assert gate["checks"]["frequency_grid_digest_consistent_when_present"] is True
    assert gate["checks"]["frequency_grid_digest_recorded_when_expected"] is True
    assert gate["checks"]["expected_frequency_grid_digest_matches"] is True
    assert gate["checks"]["interpolation_policy_consistent_when_present"] is True
    assert gate["checks"]["interpolation_policy_recorded_when_expected"] is True
    assert gate["checks"]["expected_interpolation_policy_matches"] is True
    assert gate["checks"]["touchstone_file_id_consistent_when_present"] is True
    assert gate["checks"]["touchstone_file_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_touchstone_file_id_matches"] is True
    assert gate["checks"]["touchstone_observable_id_consistent_when_present"] is True
    assert gate["checks"]["touchstone_observable_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_touchstone_observable_id_matches"] is True
    assert gate["checks"]["touchstone_observable_family_consistent_when_present"] is True
    assert gate["checks"]["touchstone_observable_family_recorded_when_expected"] is True
    assert gate["checks"]["expected_touchstone_observable_family_matches"] is True
    assert gate["checks"]["touchstone_output_artifact_id_consistent_when_present"] is True
    assert gate["checks"]["touchstone_output_digest_consistent_when_present"] is True
    assert gate["checks"]["touchstone_output_path_consistent_when_present"] is True
    assert gate["checks"]["touchstone_output_artifact_id_recorded_when_required"] is True
    assert gate["checks"]["touchstone_output_digest_recorded_when_required"] is True
    assert gate["checks"]["touchstone_output_path_recorded_when_required"] is True
    assert gate["checks"]["touchstone_output_artifact_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_touchstone_output_artifact_id_matches"] is True
    assert gate["checks"]["touchstone_output_digest_recorded_when_expected"] is True
    assert gate["checks"]["expected_touchstone_output_digest_matches"] is True
    assert gate["checks"]["touchstone_output_path_recorded_when_expected"] is True
    assert gate["checks"]["touchstone_output_schema_id_consistent_when_present"] is True
    assert gate["checks"]["touchstone_output_columns_consistent_when_present"] is True
    assert gate["checks"]["touchstone_output_units_consistent_when_present"] is True
    assert gate["checks"]["touchstone_output_schema_id_recorded_when_required"] is True
    assert gate["checks"]["touchstone_output_columns_recorded_when_required"] is True
    assert gate["checks"]["touchstone_output_units_recorded_when_required"] is True
    assert gate["checks"]["touchstone_output_schema_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_touchstone_output_schema_id_matches"] is True
    assert gate["checks"]["touchstone_output_columns_recorded_when_expected"] is True
    assert gate["checks"]["expected_touchstone_output_columns_match"] is True
    assert gate["checks"]["touchstone_output_units_recorded_when_expected"] is True
    assert gate["checks"]["expected_touchstone_output_units_match"] is True
    assert gate["checks"]["touchstone_convention_schema_id_consistent_when_present"] is True
    assert gate["checks"]["touchstone_convention_schema_id_recorded_when_required"] is True
    assert gate["checks"]["touchstone_convention_schema_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_touchstone_convention_schema_id_matches"] is True
    assert gate["checks"]["touchstone_postprocess_row_convention_schema_id_consistent_when_present"] is True
    assert gate["checks"]["touchstone_postprocess_row_convention_schema_id_recorded_when_required"] is True
    assert gate["checks"]["touchstone_postprocess_row_convention_schema_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_touchstone_postprocess_row_convention_schema_id_matches"] is True
    assert gate["checks"]["renormalized_reference_impedance_consistent_when_present"] is True
    assert gate["checks"]["renormalization_method_consistent_when_present"] is True
    assert gate["checks"]["renormalization_artifact_id_consistent_when_present"] is True
    assert gate["checks"]["renormalized_reference_impedance_recorded_when_expected"] is True
    assert gate["checks"]["expected_renormalized_reference_impedance_matches"] is True
    assert gate["checks"]["renormalization_method_recorded_when_expected"] is True
    assert gate["checks"]["expected_renormalization_method_matches"] is True
    assert gate["checks"]["renormalization_artifact_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_renormalization_artifact_id_matches"] is True
    assert gate["checks"]["deembedding_method_consistent_when_present"] is True
    assert gate["checks"]["deembedding_artifact_id_consistent_when_present"] is True
    assert gate["checks"]["deembedding_length_consistent_when_present"] is True
    assert gate["checks"]["deembedding_method_recorded_when_expected"] is True
    assert gate["checks"]["expected_deembedding_method_matches"] is True
    assert gate["checks"]["deembedding_artifact_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_deembedding_artifact_id_matches"] is True
    assert gate["checks"]["deembedding_length_recorded_when_expected"] is True
    assert gate["checks"]["expected_deembedding_length_matches"] is True
    assert gate["checks"]["selected_frequency_recorded"] is True
    assert gate["checks"]["selected_frequency_matches_design"] is True
    assert gate["checks"]["selected_frequency_matches_expected_design"] is True
    assert gate["present_kinds"]["touchstone_row"] == 1
    assert gate["frequency_grid_ids"] == ["grid_sweep_2p45g_5pts_v2"]
    assert gate["frequency_grid_digests"] == ["sha256:grid-sweep-2p45g-5pts-v2"]
    assert gate["selected_frequencies_hz"] == [2.45e9]
    assert gate["touchstone_row_indices"] == [2]
    assert gate["touchstone_row_frequencies_hz"] == [2.45e9]
    assert gate["touchstone_file_ids"] == ["touchstone_export_B_s2p_sha256_abc123"]
    assert gate["touchstone_observable_ids"] == ["cst_slot303_s2p_single_ended_sparameter_v1"]
    assert gate["touchstone_observable_families"] == ["single_ended_sparameter"]
    assert gate["incident_wave_conventions"] == ["unit_incident_power_wave_per_excited_port"]
    assert gate["power_balance_bases"] == ["power_waves_unit_incident_port"]
    assert gate["touchstone_export_methods"] == ["cst_result_tree_touchstone_s2p_export"]
    assert gate["touchstone_output_artifact_ids"] == ["touchstone_export_B_row_002_postprocess_v1"]
    assert gate["touchstone_output_digests"] == ["sha256:touchstone_export_B_row_002_postprocess_v1"]
    assert gate["touchstone_output_paths"] == ["slot295_touchstone_row_002_postprocess.json"]
    assert gate["touchstone_output_schema_ids"] == ["cst_touchstone_s2p_row_table_v1"]
    assert gate["touchstone_convention_schema_ids"] == [
        "cst_touchstone_network_convention_v1"
    ]
    assert gate["touchstone_postprocess_row_convention_schema_ids"] == [
        "cst_touchstone_s2p_row_convention_v1"
    ]
    assert gate["touchstone_output_columns"] == [
        ["frequency_hz", "S11", "S21", "S12", "S22"]
    ]
    assert gate["touchstone_output_units"] == [
        {
            "frequency_hz": "Hz",
            "S11": "1",
            "S21": "1",
            "S12": "1",
            "S22": "1",
        }
    ]
    assert gate["checks"]["created_at_utc_consistent_when_present"] is True
    assert gate["checks"]["run_timestamp_utc_consistent_when_present"] is True
    assert gate["checks"]["created_at_utc_parseable_when_present"] is True
    assert gate["checks"]["run_timestamp_utc_parseable_when_present"] is True
    assert gate["checks"]["created_at_utc_recorded_when_required"] is True
    assert gate["checks"]["run_timestamp_utc_recorded_when_required"] is True
    assert gate["checks"]["created_at_utc_recorded_when_expected"] is True
    assert gate["checks"]["expected_created_at_utc_matches"] is True
    assert gate["checks"]["run_timestamp_utc_recorded_when_expected"] is True
    assert gate["checks"]["expected_run_timestamp_utc_matches"] is True
    assert gate["checks"]["created_run_timestamp_skew_recorded"] is True
    assert gate["checks"]["created_run_timestamp_skew_within_limit"] is True
    assert gate["created_at_utc_values"] == ["2026-07-01T14:05:20Z"]
    assert gate["run_timestamp_utc_values"] == ["2026-07-01T14:05:00Z"]
    assert gate["created_run_timestamp_skews_s"] == pytest.approx([20.0] * 4)
    assert gate["max_created_run_skew_s"] == pytest.approx(60.0)
    assert gate["renormalized_reference_impedances_ohm"] == [50.0]
    assert gate["renormalization_methods"] == ["not_renormalized_option_line_R50"]
    assert gate["renormalization_artifact_ids"] == ["touchstone_export_B_renorm_none_R50_v1"]
    assert gate["deembedding_methods"] == ["port_extension_reference_plane_shift"]
    assert gate["deembedding_artifact_ids"] == ["touchstone_export_B_deembed_connector_faces_v1"]
    assert gate["deembedding_lengths_m"] == [0.0125]
    assert "S-parameter evidence cannot mix" in gate["version_note"]
    assert "renormalized Z0 identity" in gate["version_note"]
    assert "de-embedding identity" in gate["version_note"]
    assert "output artifact identity" in gate["version_note"]
    assert "observable id/family" in gate["version_note"]
    assert "incident-wave convention" in gate["version_note"]
    assert "power-balance basis" in gate["version_note"]
    assert "Touchstone export method" in gate["version_note"]
    assert "export recipe/macro/script artifact id/digest/path" in gate["version_note"]
    assert "reference-plane geometry/port-face centers" in gate["version_note"]
    assert "Touchstone option-line artifact/digest" in gate["version_note"]
    assert "model-input artifact id/digest/path" in gate["version_note"]
    assert "execution created/run timestamp identity" in gate["version_note"]
    assert "output schema/column/unit identity" in gate["version_note"]

    active = [dict(row) for row in artifacts]
    active[2]["sparameter_passivity_ok"] = False
    active_gate = cst_touchstone_solver_ready_manifest_gate(active)
    assert active_gate["status"] == "needs_attention"
    assert active_gate["checks"]["design_row_passive"] is False

    unbracketed = [dict(row) for row in artifacts]
    unbracketed[1]["design_frequency_bracketed"] = False
    unbracketed_gate = cst_touchstone_solver_ready_manifest_gate(unbracketed)
    assert unbracketed_gate["status"] == "needs_attention"
    assert unbracketed_gate["checks"]["frequency_grid_brackets_design"] is False

    missing_format = [dict(row) for row in artifacts]
    missing_format[0].pop("data_format")
    metadata_gate = cst_touchstone_solver_ready_manifest_gate(missing_format)
    assert metadata_gate["status"] == "needs_attention"
    assert metadata_gate["checks"]["port_metadata_complete"] is False

    wrong_source = [dict(row) for row in artifacts]
    wrong_source[2]["source_tool"] = "HFSS"
    wrong_source_gate = cst_touchstone_solver_ready_manifest_gate(wrong_source)
    assert wrong_source_gate["status"] == "needs_attention"
    assert wrong_source_gate["checks"]["source_tool_is_cst"] is False

    stale_model_input_digest = [dict(row) for row in artifacts]
    stale_model_input_digest[2]["model_input_digest"] = "sha256:cst-slot380-rf-filter-project-old"
    stale_model_input_digest_gate = cst_touchstone_solver_ready_manifest_gate(
        stale_model_input_digest,
        expected_model_input_artifact_id="cst_slot380_rf_filter_project_v1.cst",
        expected_model_input_digest="sha256:cst-slot380-rf-filter-project-v1",
        expected_model_input_path="artifacts/rf/cst_slot380_rf_filter_project_v1.cst",
        require_model_input_artifact=True,
    )
    assert stale_model_input_digest_gate["status"] == "needs_attention"
    assert (
        stale_model_input_digest_gate["checks"]["model_input_digest_consistent_when_present"]
        is False
    )
    assert stale_model_input_digest_gate["checks"]["expected_model_input_digest_matches"] is False

    missing_model_input_path = [dict(row) for row in artifacts]
    missing_model_input_path[1].pop("model_input_path")
    missing_model_input_path_gate = cst_touchstone_solver_ready_manifest_gate(
        missing_model_input_path,
        expected_model_input_path="artifacts/rf/cst_slot380_rf_filter_project_v1.cst",
        require_model_input_artifact=True,
    )
    assert missing_model_input_path_gate["status"] == "needs_attention"
    assert missing_model_input_path_gate["checks"]["model_input_path_recorded_when_required"] is False
    assert missing_model_input_path_gate["checks"]["expected_model_input_path_matches"] is False

    stale_export_recipe_digest = [dict(row) for row in artifacts]
    stale_export_recipe_digest[2]["export_recipe_digest"] = (
        "sha256:cst-slot319-old-export-macro"
    )
    stale_export_recipe_digest_gate = cst_touchstone_solver_ready_manifest_gate(
        stale_export_recipe_digest,
        expected_export_recipe_artifact_id="cst_slot387_touchstone_export_macro_v1.bas",
        expected_export_recipe_digest="sha256:cst-slot387-touchstone-export-macro-v1",
        expected_export_recipe_path=(
            "artifacts/rf/cst_slot387_touchstone_export_macro_v1.bas"
        ),
        require_export_recipe_artifact=True,
    )
    assert stale_export_recipe_digest_gate["status"] == "needs_attention"
    assert (
        stale_export_recipe_digest_gate["checks"][
            "export_recipe_digest_consistent_when_present"
        ]
        is False
    )
    assert (
        stale_export_recipe_digest_gate["checks"]["expected_export_recipe_digest_matches"]
        is False
    )
    assert (
        stale_export_recipe_digest_gate["checks"][
            "expected_export_recipe_artifact_id_matches"
        ]
        is True
    )

    missing_export_recipe_path = [dict(row) for row in artifacts]
    missing_export_recipe_path[1].pop("export_recipe_path")
    missing_export_recipe_path_gate = cst_touchstone_solver_ready_manifest_gate(
        missing_export_recipe_path,
        expected_export_recipe_path=(
            "artifacts/rf/cst_slot387_touchstone_export_macro_v1.bas"
        ),
        require_export_recipe_artifact=True,
    )
    assert missing_export_recipe_path_gate["status"] == "needs_attention"
    assert (
        missing_export_recipe_path_gate["checks"][
            "export_recipe_path_recorded_when_required"
        ]
        is False
    )
    assert (
        missing_export_recipe_path_gate["checks"]["expected_export_recipe_path_matches"]
        is False
    )

    stale_parameter_set_digest = [dict(row) for row in artifacts]
    stale_parameter_set_digest[2]["parameter_set_digest"] = (
        "sha256:cst-slot380-old-rf-filter-parameter-set"
    )
    stale_parameter_set_digest_gate = cst_touchstone_solver_ready_manifest_gate(
        stale_parameter_set_digest,
        expected_parameter_set_artifact_id="cst_slot394_rf_filter_parameter_set_v1.json",
        expected_parameter_set_digest="sha256:cst-slot394-rf-filter-parameter-set-v1",
        expected_parameter_set_path=(
            "artifacts/rf/cst_slot394_rf_filter_parameter_set.json"
        ),
        require_parameter_set_artifact=True,
    )
    assert stale_parameter_set_digest_gate["status"] == "needs_attention"
    assert (
        stale_parameter_set_digest_gate["checks"][
            "parameter_set_digest_consistent_when_present"
        ]
        is False
    )
    assert (
        stale_parameter_set_digest_gate["checks"][
            "expected_parameter_set_digest_matches"
        ]
        is False
    )
    assert (
        stale_parameter_set_digest_gate["checks"][
            "expected_parameter_set_artifact_id_matches"
        ]
        is True
    )

    missing_parameter_set_path = [dict(row) for row in artifacts]
    missing_parameter_set_path[1].pop("parameter_set_path")
    missing_parameter_set_path_gate = cst_touchstone_solver_ready_manifest_gate(
        missing_parameter_set_path,
        expected_parameter_set_path=(
            "artifacts/rf/cst_slot394_rf_filter_parameter_set.json"
        ),
        require_parameter_set_artifact=True,
    )
    assert missing_parameter_set_path_gate["status"] == "needs_attention"
    assert (
        missing_parameter_set_path_gate["checks"][
            "parameter_set_path_recorded_when_required"
        ]
        is False
    )
    assert (
        missing_parameter_set_path_gate["checks"]["expected_parameter_set_path_matches"]
        is False
    )

    wrong_objective_family = [dict(row) for row in artifacts]
    wrong_objective_family[2]["objective_observable_family"] = (
        "touchstone_return_loss_objective"
    )
    wrong_objective_family_gate = cst_touchstone_solver_ready_manifest_gate(
        wrong_objective_family,
        expected_objective_observable_id="cst_slot394_s21_insertion_loss_objective_v1",
        expected_objective_observable_family="touchstone_s21_insertion_loss_objective",
    )
    assert wrong_objective_family_gate["status"] == "needs_attention"
    assert (
        wrong_objective_family_gate["checks"][
            "expected_objective_observable_id_matches"
        ]
        is True
    )
    assert (
        wrong_objective_family_gate["checks"][
            "objective_observable_family_consistent_when_present"
        ]
        is False
    )
    assert (
        wrong_objective_family_gate["checks"][
            "expected_objective_observable_family_matches"
        ]
        is False
    )

    missing_row_index = [dict(row) for row in artifacts]
    missing_row_index[2].pop("selected_row_index")
    missing_row_index_gate = cst_touchstone_solver_ready_manifest_gate(missing_row_index)
    assert missing_row_index_gate["status"] == "needs_attention"
    assert missing_row_index_gate["checks"]["design_row_index_recorded"] is False

    out_of_range = [dict(row) for row in artifacts]
    out_of_range[2]["selected_row_index"] = 99
    out_of_range_gate = cst_touchstone_solver_ready_manifest_gate(out_of_range)
    assert out_of_range_gate["status"] == "needs_attention"
    assert out_of_range_gate["checks"]["design_row_index_within_grid"] is False

    wrong_network_kind = [dict(row) for row in artifacts]
    wrong_network_kind[2]["network_kind"] = "Y"
    wrong_network_gate = cst_touchstone_solver_ready_manifest_gate(
        wrong_network_kind,
        expected_network_kind="S",
    )
    assert wrong_network_gate["status"] == "needs_attention"
    assert wrong_network_gate["checks"]["network_kind_unique"] is False
    assert wrong_network_gate["checks"]["expected_network_kind_matches"] is False

    swapped_port_order = [dict(row) for row in artifacts]
    swapped_port_order[2]["port_order"] = ["P2", "P1"]
    swapped_port_gate = cst_touchstone_solver_ready_manifest_gate(
        swapped_port_order,
        expected_port_order=("P1", "P2"),
    )
    assert swapped_port_gate["status"] == "needs_attention"
    assert swapped_port_gate["checks"]["port_order_unique"] is False
    assert swapped_port_gate["checks"]["expected_port_order_matches"] is False

    wrong_format = [dict(row) for row in artifacts]
    wrong_format[2]["data_format"] = "DB"
    wrong_format_gate = cst_touchstone_solver_ready_manifest_gate(
        wrong_format,
        expected_data_format="MA",
    )
    assert wrong_format_gate["status"] == "needs_attention"
    assert wrong_format_gate["checks"]["data_format_unique"] is False
    assert wrong_format_gate["checks"]["expected_data_format_matches"] is False

    wrong_z0 = [dict(row) for row in artifacts]
    wrong_z0[1]["reference_impedance_ohm"] = 75.0
    wrong_z0_gate = cst_touchstone_solver_ready_manifest_gate(
        wrong_z0,
        expected_reference_impedance_ohm=50.0,
    )
    assert wrong_z0_gate["status"] == "needs_attention"
    assert wrong_z0_gate["checks"]["reference_impedance_unique"] is False
    assert wrong_z0_gate["checks"]["expected_reference_impedance_matches"] is False

    stale_option_line = [dict(row) for row in artifacts]
    stale_option_line[2]["touchstone_option_line_artifact_id"] = (
        "touchstone_export_A_option_line_old"
    )
    stale_option_line_gate = cst_touchstone_solver_ready_manifest_gate(
        stale_option_line,
        expected_touchstone_option_line_artifact_id="touchstone_export_B_option_line_v1",
        expected_touchstone_option_line_digest=(
            "sha256:touchstone-export-B-option-line-ghz-s-ma-r50-v1"
        ),
    )
    assert stale_option_line_gate["status"] == "needs_attention"
    assert (
        stale_option_line_gate["checks"][
            "touchstone_option_line_artifact_id_consistent_when_present"
        ]
        is False
    )
    assert (
        stale_option_line_gate["checks"][
            "expected_touchstone_option_line_artifact_id_matches"
        ]
        is False
    )
    assert (
        stale_option_line_gate["checks"]["expected_touchstone_option_line_digest_matches"]
        is True
    )

    missing_option_line_digest = [dict(row) for row in artifacts]
    missing_option_line_digest[1].pop("touchstone_option_line_digest")
    missing_option_line_digest_gate = cst_touchstone_solver_ready_manifest_gate(
        missing_option_line_digest,
        expected_touchstone_option_line_digest=(
            "sha256:touchstone-export-B-option-line-ghz-s-ma-r50-v1"
        ),
    )
    assert missing_option_line_digest_gate["status"] == "needs_attention"
    assert (
        missing_option_line_digest_gate["checks"][
            "touchstone_option_line_digest_recorded_when_expected"
        ]
        is False
    )
    assert (
        missing_option_line_digest_gate["checks"][
            "expected_touchstone_option_line_digest_matches"
        ]
        is False
    )

    missing_reference_plane = [dict(row) for row in artifacts]
    missing_reference_plane[1].pop("reference_plane")
    missing_reference_plane_gate = cst_touchstone_solver_ready_manifest_gate(
        missing_reference_plane,
        expected_reference_plane="cst_port1_port2_deembedded_to_connector_faces",
    )
    assert missing_reference_plane_gate["status"] == "needs_attention"
    assert missing_reference_plane_gate["checks"]["reference_plane_recorded_when_expected"] is False
    assert missing_reference_plane_gate["checks"]["expected_reference_plane_matches"] is False

    stale_reference_plane_geometry = [dict(row) for row in artifacts]
    stale_reference_plane_geometry[2]["reference_plane_geometry_digest"] = (
        "sha256:cst_slot327_connector_face_plane_geometry_old"
    )
    stale_reference_plane_geometry_gate = cst_touchstone_solver_ready_manifest_gate(
        stale_reference_plane_geometry,
        expected_reference_plane_geometry_digest="sha256:cst_slot327_connector_face_plane_geometry_v1",
        expected_port_face_centers_xyz_m=((0.0, 0.0, 0.0), (0.05, 0.0, 0.0)),
    )
    assert stale_reference_plane_geometry_gate["status"] == "needs_attention"
    assert (
        stale_reference_plane_geometry_gate["checks"][
            "expected_reference_plane_geometry_digest_matches"
        ]
        is False
    )
    assert stale_reference_plane_geometry_gate["checks"]["expected_port_face_centers_xyz_matches"] is True

    shifted_port_face_center = [dict(row) for row in artifacts]
    shifted_port_face_center[1]["port_face_centers_xyz_m"] = [
        [0.0, 0.0, 0.0],
        [0.051, 0.0, 0.0],
    ]
    shifted_port_face_center_gate = cst_touchstone_solver_ready_manifest_gate(
        shifted_port_face_center,
        expected_reference_plane_geometry_digest="sha256:cst_slot327_connector_face_plane_geometry_v1",
        expected_port_face_centers_xyz_m=((0.0, 0.0, 0.0), (0.05, 0.0, 0.0)),
    )
    assert shifted_port_face_center_gate["status"] == "needs_attention"
    assert (
        shifted_port_face_center_gate["checks"][
            "expected_reference_plane_geometry_digest_matches"
        ]
        is True
    )
    assert shifted_port_face_center_gate["checks"]["expected_port_face_centers_xyz_matches"] is False

    missing_reference_plane_geometry = [dict(row) for row in artifacts]
    missing_reference_plane_geometry[1].pop("reference_plane_geometry_digest")
    missing_reference_plane_geometry[1].pop("port_face_centers_xyz_m")
    missing_reference_plane_geometry_gate = cst_touchstone_solver_ready_manifest_gate(
        missing_reference_plane_geometry,
        expected_reference_plane_geometry_digest="sha256:cst_slot327_connector_face_plane_geometry_v1",
        expected_port_face_centers_xyz_m=((0.0, 0.0, 0.0), (0.05, 0.0, 0.0)),
    )
    assert missing_reference_plane_geometry_gate["status"] == "needs_attention"
    assert (
        missing_reference_plane_geometry_gate["checks"][
            "reference_plane_geometry_digest_recorded_when_expected"
        ]
        is False
    )
    assert (
        missing_reference_plane_geometry_gate["checks"]["port_face_centers_xyz_recorded_when_expected"]
        is False
    )

    wrong_mode_basis = [dict(row) for row in artifacts]
    wrong_mode_basis[2]["port_mode_basis"] = "mixed_mode_differential_common"
    wrong_mode_basis_gate = cst_touchstone_solver_ready_manifest_gate(
        wrong_mode_basis,
        expected_port_mode_basis="single_ended_power_wave_modes",
    )
    assert wrong_mode_basis_gate["status"] == "needs_attention"
    assert wrong_mode_basis_gate["checks"]["port_mode_basis_consistent_when_present"] is False
    assert wrong_mode_basis_gate["checks"]["expected_port_mode_basis_matches"] is False

    wrong_incident_wave = [dict(row) for row in artifacts]
    wrong_incident_wave[2]["incident_wave_convention"] = "unit_voltage_wave_1v_peak"
    wrong_incident_wave_gate = cst_touchstone_solver_ready_manifest_gate(
        wrong_incident_wave,
        expected_incident_wave_convention="unit_incident_power_wave_per_excited_port",
    )
    assert wrong_incident_wave_gate["status"] == "needs_attention"
    assert wrong_incident_wave_gate["checks"]["incident_wave_convention_consistent_when_present"] is False
    assert wrong_incident_wave_gate["checks"]["expected_incident_wave_convention_matches"] is False

    wrong_power_basis = [dict(row) for row in artifacts]
    wrong_power_basis[1]["power_balance_basis"] = "voltage_wave_amplitude_not_power_normalized"
    wrong_power_basis_gate = cst_touchstone_solver_ready_manifest_gate(
        wrong_power_basis,
        expected_power_balance_basis="power_waves_unit_incident_port",
    )
    assert wrong_power_basis_gate["status"] == "needs_attention"
    assert wrong_power_basis_gate["checks"]["power_balance_basis_consistent_when_present"] is False
    assert wrong_power_basis_gate["checks"]["expected_power_balance_basis_matches"] is False

    wrong_touchstone_export_method = [dict(row) for row in artifacts]
    wrong_touchstone_export_method[2]["touchstone_export_method"] = "cst_ascii_table_copy_paste_export"
    wrong_touchstone_export_method_gate = cst_touchstone_solver_ready_manifest_gate(
        wrong_touchstone_export_method,
        expected_touchstone_export_method="cst_result_tree_touchstone_s2p_export",
    )
    assert wrong_touchstone_export_method_gate["status"] == "needs_attention"
    assert (
        wrong_touchstone_export_method_gate["checks"]["touchstone_export_method_consistent_when_present"]
        is False
    )
    assert (
        wrong_touchstone_export_method_gate["checks"]["expected_touchstone_export_method_matches"]
        is False
    )

    missing_touchstone_export_method = [dict(row) for row in artifacts]
    missing_touchstone_export_method[1].pop("touchstone_export_method")
    missing_touchstone_export_method_gate = cst_touchstone_solver_ready_manifest_gate(
        missing_touchstone_export_method,
        expected_touchstone_export_method="cst_result_tree_touchstone_s2p_export",
    )
    assert missing_touchstone_export_method_gate["status"] == "needs_attention"
    assert (
        missing_touchstone_export_method_gate["checks"]["touchstone_export_method_recorded_when_expected"]
        is False
    )
    assert (
        missing_touchstone_export_method_gate["checks"]["expected_touchstone_export_method_matches"]
        is False
    )

    stale_result_tree_path = [dict(row) for row in artifacts]
    stale_result_tree_path[2]["result_tree_path"] = "1D Results\\Tables\\Copied S-Parameters"
    stale_result_tree_path_gate = cst_touchstone_solver_ready_manifest_gate(
        stale_result_tree_path,
        expected_result_tree_path="1D Results\\S-Parameters\\Touchstone\\S2P Export",
        expected_result_item_id="cst_slot343_result_tree_sparams_s2p_v1",
    )
    assert stale_result_tree_path_gate["status"] == "needs_attention"
    assert (
        stale_result_tree_path_gate["checks"]["result_tree_path_consistent_when_present"]
        is False
    )
    assert stale_result_tree_path_gate["checks"]["expected_result_tree_path_matches"] is False
    assert stale_result_tree_path_gate["checks"]["expected_result_item_id_matches"] is True

    missing_result_item_id = [dict(row) for row in artifacts]
    missing_result_item_id[1].pop("result_item_id")
    missing_result_item_id_gate = cst_touchstone_solver_ready_manifest_gate(
        missing_result_item_id,
        expected_result_item_id="cst_slot343_result_tree_sparams_s2p_v1",
    )
    assert missing_result_item_id_gate["status"] == "needs_attention"
    assert missing_result_item_id_gate["checks"]["result_item_id_recorded_when_expected"] is False
    assert missing_result_item_id_gate["checks"]["expected_result_item_id_matches"] is False

    stale_grid_id = [dict(row) for row in artifacts]
    stale_grid_id[2]["frequency_grid_id"] = "grid_sweep_2p40g_5pts_old"
    stale_grid_gate = cst_touchstone_solver_ready_manifest_gate(
        stale_grid_id,
        expected_frequency_grid_id="grid_sweep_2p45g_5pts_v2",
    )
    assert stale_grid_gate["status"] == "needs_attention"
    assert stale_grid_gate["checks"]["frequency_grid_id_consistent_when_present"] is False
    assert stale_grid_gate["checks"]["expected_frequency_grid_id_matches"] is False

    stale_grid_digest = [dict(row) for row in artifacts]
    stale_grid_digest[2]["frequency_grid_digest"] = "sha256:grid-sweep-2p40g-5pts-old"
    stale_grid_digest_gate = cst_touchstone_solver_ready_manifest_gate(
        stale_grid_digest,
        expected_frequency_grid_digest="sha256:grid-sweep-2p45g-5pts-v2",
    )
    assert stale_grid_digest_gate["status"] == "needs_attention"
    assert (
        stale_grid_digest_gate["checks"]["frequency_grid_digest_consistent_when_present"]
        is False
    )
    assert stale_grid_digest_gate["checks"]["expected_frequency_grid_digest_matches"] is False

    missing_grid_digest = [dict(row) for row in artifacts]
    missing_grid_digest[1].pop("frequency_grid_digest")
    missing_grid_digest_gate = cst_touchstone_solver_ready_manifest_gate(
        missing_grid_digest,
        expected_frequency_grid_digest="sha256:grid-sweep-2p45g-5pts-v2",
    )
    assert missing_grid_digest_gate["status"] == "needs_attention"
    assert (
        missing_grid_digest_gate["checks"]["frequency_grid_digest_recorded_when_expected"]
        is False
    )
    assert missing_grid_digest_gate["checks"]["expected_frequency_grid_digest_matches"] is False

    stale_file_id = [dict(row) for row in artifacts]
    stale_file_id[3]["touchstone_file_id"] = "touchstone_export_A_s2p_sha256_old"
    stale_file_gate = cst_touchstone_solver_ready_manifest_gate(
        stale_file_id,
        expected_touchstone_file_id="touchstone_export_B_s2p_sha256_abc123",
    )
    assert stale_file_gate["status"] == "needs_attention"
    assert stale_file_gate["checks"]["touchstone_file_id_consistent_when_present"] is False
    assert stale_file_gate["checks"]["expected_touchstone_file_id_matches"] is False

    stale_observable_id = [dict(row) for row in artifacts]
    stale_observable_id[3]["touchstone_observable_id"] = "cst_slot287_old_mixed_mode_sparameter"
    stale_observable_id_gate = cst_touchstone_solver_ready_manifest_gate(
        stale_observable_id,
        expected_touchstone_observable_id="cst_slot303_s2p_single_ended_sparameter_v1",
        expected_touchstone_observable_family="single_ended_sparameter",
        expected_touchstone_output_artifact_id="touchstone_export_B_row_002_postprocess_v1",
        expected_touchstone_output_digest="sha256:touchstone_export_B_row_002_postprocess_v1",
        require_touchstone_output_artifact=True,
    )
    assert stale_observable_id_gate["status"] == "needs_attention"
    assert stale_observable_id_gate["checks"]["touchstone_observable_id_consistent_when_present"] is False
    assert stale_observable_id_gate["checks"]["expected_touchstone_observable_id_matches"] is False
    assert stale_observable_id_gate["checks"]["expected_touchstone_observable_family_matches"] is True
    assert stale_observable_id_gate["checks"]["expected_touchstone_output_artifact_id_matches"] is True

    wrong_observable_family = [dict(row) for row in artifacts]
    wrong_observable_family[2]["touchstone_observable_family"] = "mixed_mode_sparameter"
    wrong_observable_family_gate = cst_touchstone_solver_ready_manifest_gate(
        wrong_observable_family,
        expected_touchstone_observable_id="cst_slot303_s2p_single_ended_sparameter_v1",
        expected_touchstone_observable_family="single_ended_sparameter",
        expected_touchstone_output_artifact_id="touchstone_export_B_row_002_postprocess_v1",
        expected_touchstone_output_digest="sha256:touchstone_export_B_row_002_postprocess_v1",
        require_touchstone_output_artifact=True,
    )
    assert wrong_observable_family_gate["status"] == "needs_attention"
    assert wrong_observable_family_gate["checks"]["expected_touchstone_observable_id_matches"] is True
    assert wrong_observable_family_gate["checks"]["touchstone_observable_family_consistent_when_present"] is False
    assert wrong_observable_family_gate["checks"]["expected_touchstone_observable_family_matches"] is False
    assert wrong_observable_family_gate["checks"]["expected_touchstone_output_artifact_id_matches"] is True

    stale_output_artifact = [dict(row) for row in artifacts]
    stale_output_artifact[3]["touchstone_output_artifact_id"] = "touchstone_export_A_row_002_old"
    stale_output_artifact_gate = cst_touchstone_solver_ready_manifest_gate(
        stale_output_artifact,
        expected_touchstone_output_artifact_id="touchstone_export_B_row_002_postprocess_v1",
        expected_touchstone_output_digest="sha256:touchstone_export_B_row_002_postprocess_v1",
        require_touchstone_output_artifact=True,
    )
    assert stale_output_artifact_gate["status"] == "needs_attention"
    assert stale_output_artifact_gate["checks"]["touchstone_output_artifact_id_consistent_when_present"] is False
    assert stale_output_artifact_gate["checks"]["expected_touchstone_output_artifact_id_matches"] is False

    stale_output_digest = [dict(row) for row in artifacts]
    stale_output_digest[2]["touchstone_output_digest"] = "sha256:touchstone_export_A_row_002_old"
    stale_output_digest_gate = cst_touchstone_solver_ready_manifest_gate(
        stale_output_digest,
        expected_touchstone_output_artifact_id="touchstone_export_B_row_002_postprocess_v1",
        expected_touchstone_output_digest="sha256:touchstone_export_B_row_002_postprocess_v1",
        require_touchstone_output_artifact=True,
    )
    assert stale_output_digest_gate["status"] == "needs_attention"
    assert stale_output_digest_gate["checks"]["touchstone_output_digest_consistent_when_present"] is False
    assert stale_output_digest_gate["checks"]["expected_touchstone_output_digest_matches"] is False

    stale_output_schema = [dict(row) for row in artifacts]
    stale_output_schema[3]["touchstone_output_schema_id"] = (
        "cst_touchstone_scalar_s21_v0"
    )
    stale_output_schema[3]["touchstone_output_columns"] = ["frequency_hz", "S21"]
    stale_output_schema[3]["touchstone_output_units"] = {
        "frequency_hz": "GHz",
        "S21": "dB",
    }
    stale_output_schema_gate = cst_touchstone_solver_ready_manifest_gate(
        stale_output_schema,
        expected_touchstone_output_artifact_id="touchstone_export_B_row_002_postprocess_v1",
        expected_touchstone_output_digest="sha256:touchstone_export_B_row_002_postprocess_v1",
        expected_touchstone_output_schema_id="cst_touchstone_s2p_row_table_v1",
        expected_touchstone_output_columns=[
            "frequency_hz",
            "S11",
            "S21",
            "S12",
            "S22",
        ],
        expected_touchstone_output_units={
            "frequency_hz": "Hz",
            "S11": "1",
            "S21": "1",
            "S12": "1",
            "S22": "1",
        },
        require_touchstone_output_artifact=True,
        require_touchstone_output_schema=True,
    )
    assert stale_output_schema_gate["status"] == "needs_attention"
    assert (
        stale_output_schema_gate["checks"][
            "touchstone_output_artifact_id_consistent_when_present"
        ]
        is True
    )
    assert (
        stale_output_schema_gate["checks"][
            "touchstone_output_digest_consistent_when_present"
        ]
        is True
    )
    assert (
        stale_output_schema_gate["checks"][
            "touchstone_output_schema_id_consistent_when_present"
        ]
        is False
    )
    assert (
        stale_output_schema_gate["checks"]["expected_touchstone_output_schema_id_matches"]
        is False
    )
    assert (
        stale_output_schema_gate["checks"]["expected_touchstone_output_columns_match"]
        is False
    )
    assert (
        stale_output_schema_gate["checks"]["expected_touchstone_output_units_match"]
        is False
    )

    missing_output_schema = [dict(row) for row in artifacts]
    missing_output_schema[1].pop("touchstone_output_schema_id")
    missing_output_schema[1].pop("touchstone_output_columns")
    missing_output_schema[1].pop("touchstone_output_units")
    missing_output_schema_gate = cst_touchstone_solver_ready_manifest_gate(
        missing_output_schema,
        expected_touchstone_output_schema_id="cst_touchstone_s2p_row_table_v1",
        expected_touchstone_output_columns=[
            "frequency_hz",
            "S11",
            "S21",
            "S12",
            "S22",
        ],
        expected_touchstone_output_units={
            "frequency_hz": "Hz",
            "S11": "1",
            "S21": "1",
            "S12": "1",
            "S22": "1",
        },
        require_touchstone_output_schema=True,
    )
    assert missing_output_schema_gate["status"] == "needs_attention"
    assert (
        missing_output_schema_gate["checks"][
            "touchstone_output_schema_id_recorded_when_required"
        ]
        is False
    )
    assert (
        missing_output_schema_gate["checks"][
            "touchstone_output_columns_recorded_when_required"
        ]
        is False
    )
    assert (
        missing_output_schema_gate["checks"][
            "touchstone_output_units_recorded_when_required"
        ]
        is False
    )
    assert missing_output_schema_gate["missing_touchstone_output_schema_id_rows"] == [2]
    assert missing_output_schema_gate["missing_touchstone_output_columns_rows"] == [2]
    assert missing_output_schema_gate["missing_touchstone_output_units_rows"] == [2]

    stale_touchstone_convention_schema = [dict(row) for row in artifacts]
    stale_touchstone_convention_schema[3]["touchstone_convention_schema_id"] = (
        "cst_touchstone_value_only_convention_v0"
    )
    stale_touchstone_convention_schema_gate = cst_touchstone_solver_ready_manifest_gate(
        stale_touchstone_convention_schema,
        expected_touchstone_output_schema_id="cst_touchstone_s2p_row_table_v1",
        expected_touchstone_convention_schema_id="cst_touchstone_network_convention_v1",
        require_touchstone_output_schema=True,
        require_touchstone_convention_schema=True,
    )
    assert stale_touchstone_convention_schema_gate["status"] == "needs_attention"
    assert (
        stale_touchstone_convention_schema_gate["checks"][
            "expected_touchstone_output_schema_id_matches"
        ]
        is True
    )
    assert (
        stale_touchstone_convention_schema_gate["checks"][
            "touchstone_convention_schema_id_consistent_when_present"
        ]
        is False
    )
    assert (
        stale_touchstone_convention_schema_gate["checks"][
            "expected_touchstone_convention_schema_id_matches"
        ]
        is False
    )

    missing_touchstone_convention_schema = [dict(row) for row in artifacts]
    missing_touchstone_convention_schema[1].pop("touchstone_convention_schema_id")
    missing_touchstone_convention_schema_gate = cst_touchstone_solver_ready_manifest_gate(
        missing_touchstone_convention_schema,
        expected_touchstone_convention_schema_id="cst_touchstone_network_convention_v1",
        require_touchstone_convention_schema=True,
    )
    assert missing_touchstone_convention_schema_gate["status"] == "needs_attention"
    assert (
        missing_touchstone_convention_schema_gate["checks"][
            "touchstone_convention_schema_id_recorded_when_required"
        ]
        is False
    )
    assert (
        missing_touchstone_convention_schema_gate["checks"][
            "expected_touchstone_convention_schema_id_matches"
        ]
        is False
    )
    assert missing_touchstone_convention_schema_gate[
        "missing_touchstone_convention_schema_id_rows"
    ] == [2]

    stale_touchstone_postprocess_row_convention_schema = [dict(row) for row in artifacts]
    stale_touchstone_postprocess_row_convention_schema[3][
        "touchstone_postprocess_row_convention_schema_id"
    ] = "cst_touchstone_scalar_s21_row_v0"
    stale_touchstone_postprocess_row_convention_schema_gate = (
        cst_touchstone_solver_ready_manifest_gate(
            stale_touchstone_postprocess_row_convention_schema,
            expected_touchstone_output_schema_id="cst_touchstone_s2p_row_table_v1",
            expected_touchstone_convention_schema_id="cst_touchstone_network_convention_v1",
            expected_touchstone_postprocess_row_convention_schema_id=(
                "cst_touchstone_s2p_row_convention_v1"
            ),
            require_touchstone_output_schema=True,
            require_touchstone_convention_schema=True,
            require_touchstone_postprocess_row_convention_schema=True,
        )
    )
    assert stale_touchstone_postprocess_row_convention_schema_gate["status"] == "needs_attention"
    assert (
        stale_touchstone_postprocess_row_convention_schema_gate["checks"][
            "expected_touchstone_output_schema_id_matches"
        ]
        is True
    )
    assert (
        stale_touchstone_postprocess_row_convention_schema_gate["checks"][
            "expected_touchstone_convention_schema_id_matches"
        ]
        is True
    )
    assert (
        stale_touchstone_postprocess_row_convention_schema_gate["checks"][
            "touchstone_postprocess_row_convention_schema_id_consistent_when_present"
        ]
        is False
    )
    assert (
        stale_touchstone_postprocess_row_convention_schema_gate["checks"][
            "expected_touchstone_postprocess_row_convention_schema_id_matches"
        ]
        is False
    )

    missing_touchstone_postprocess_row_convention_schema = [dict(row) for row in artifacts]
    missing_touchstone_postprocess_row_convention_schema[1].pop(
        "touchstone_postprocess_row_convention_schema_id"
    )
    missing_touchstone_postprocess_row_convention_schema_gate = (
        cst_touchstone_solver_ready_manifest_gate(
            missing_touchstone_postprocess_row_convention_schema,
            expected_touchstone_postprocess_row_convention_schema_id=(
                "cst_touchstone_s2p_row_convention_v1"
            ),
            require_touchstone_postprocess_row_convention_schema=True,
        )
    )
    assert missing_touchstone_postprocess_row_convention_schema_gate["status"] == "needs_attention"
    assert (
        missing_touchstone_postprocess_row_convention_schema_gate["checks"][
            "touchstone_postprocess_row_convention_schema_id_recorded_when_required"
        ]
        is False
    )
    assert (
        missing_touchstone_postprocess_row_convention_schema_gate["checks"][
            "expected_touchstone_postprocess_row_convention_schema_id_matches"
        ]
        is False
    )
    assert missing_touchstone_postprocess_row_convention_schema_gate[
        "missing_touchstone_postprocess_row_convention_schema_id_rows"
    ] == [2]

    stale_touchstone_port_mode_basis_schema = [dict(row) for row in artifacts]
    stale_touchstone_port_mode_basis_schema[3][
        "touchstone_port_mode_basis_schema_id"
    ] = "cst_mixed_mode_port_mode_basis_v0"
    stale_touchstone_port_mode_basis_schema_gate = (
        cst_touchstone_solver_ready_manifest_gate(
            stale_touchstone_port_mode_basis_schema,
            expected_port_mode_basis="single_ended_power_wave_modes",
            expected_touchstone_port_mode_basis_schema_id=(
                "cst_single_ended_power_wave_port_mode_basis_v1"
            ),
            expected_touchstone_output_schema_id="cst_touchstone_s2p_row_table_v1",
            expected_touchstone_convention_schema_id="cst_touchstone_network_convention_v1",
            expected_touchstone_postprocess_row_convention_schema_id=(
                "cst_touchstone_s2p_row_convention_v1"
            ),
            require_touchstone_output_schema=True,
            require_touchstone_convention_schema=True,
            require_touchstone_port_mode_basis_schema=True,
            require_touchstone_postprocess_row_convention_schema=True,
        )
    )
    assert stale_touchstone_port_mode_basis_schema_gate["status"] == "needs_attention"
    assert (
        stale_touchstone_port_mode_basis_schema_gate["checks"][
            "expected_port_mode_basis_matches"
        ]
        is True
    )
    assert (
        stale_touchstone_port_mode_basis_schema_gate["checks"][
            "touchstone_port_mode_basis_schema_id_consistent_when_present"
        ]
        is False
    )
    assert (
        stale_touchstone_port_mode_basis_schema_gate["checks"][
            "expected_touchstone_port_mode_basis_schema_id_matches"
        ]
        is False
    )
    assert (
        stale_touchstone_port_mode_basis_schema_gate["checks"][
            "expected_touchstone_output_schema_id_matches"
        ]
        is True
    )
    assert (
        stale_touchstone_port_mode_basis_schema_gate["checks"][
            "expected_touchstone_convention_schema_id_matches"
        ]
        is True
    )
    assert (
        stale_touchstone_port_mode_basis_schema_gate["checks"][
            "expected_touchstone_postprocess_row_convention_schema_id_matches"
        ]
        is True
    )

    missing_touchstone_port_mode_basis_schema = [dict(row) for row in artifacts]
    missing_touchstone_port_mode_basis_schema[1].pop(
        "touchstone_port_mode_basis_schema_id"
    )
    missing_touchstone_port_mode_basis_schema_gate = (
        cst_touchstone_solver_ready_manifest_gate(
            missing_touchstone_port_mode_basis_schema,
            expected_touchstone_port_mode_basis_schema_id=(
                "cst_single_ended_power_wave_port_mode_basis_v1"
            ),
            require_touchstone_port_mode_basis_schema=True,
        )
    )
    assert missing_touchstone_port_mode_basis_schema_gate["status"] == "needs_attention"
    assert (
        missing_touchstone_port_mode_basis_schema_gate["checks"][
            "touchstone_port_mode_basis_schema_id_recorded_when_required"
        ]
        is False
    )
    assert (
        missing_touchstone_port_mode_basis_schema_gate["checks"][
            "expected_touchstone_port_mode_basis_schema_id_matches"
        ]
        is False
    )
    assert missing_touchstone_port_mode_basis_schema_gate[
        "missing_touchstone_port_mode_basis_schema_id_rows"
    ] == [2]

    missing_output_path = [dict(row) for row in artifacts]
    missing_output_path[1].pop("touchstone_output_path")
    missing_output_path_gate = cst_touchstone_solver_ready_manifest_gate(
        missing_output_path,
        expected_touchstone_output_artifact_id="touchstone_export_B_row_002_postprocess_v1",
        expected_touchstone_output_digest="sha256:touchstone_export_B_row_002_postprocess_v1",
        require_touchstone_output_artifact=True,
    )
    assert missing_output_path_gate["status"] == "needs_attention"
    assert missing_output_path_gate["checks"]["touchstone_output_path_recorded_when_required"] is False
    assert missing_output_path_gate["checks"]["touchstone_output_path_recorded_when_expected"] is False

    stale_execution_created = [dict(row) for row in artifacts]
    stale_execution_created[2]["created_at_utc"] = "2026-07-01T16:05:20Z"
    stale_execution_created_gate = cst_touchstone_solver_ready_manifest_gate(
        stale_execution_created,
        expected_created_at_utc="2026-07-01T14:05:20Z",
        expected_run_timestamp_utc="2026-07-01T14:05:00Z",
        max_created_run_skew_s=60.0,
        require_execution_metadata=True,
    )
    assert stale_execution_created_gate["status"] == "needs_attention"
    assert (
        stale_execution_created_gate["checks"]["created_at_utc_consistent_when_present"]
        is False
    )
    assert stale_execution_created_gate["checks"]["expected_created_at_utc_matches"] is False
    assert (
        stale_execution_created_gate["checks"]["created_run_timestamp_skew_within_limit"]
        is False
    )

    bad_run_timestamp = [dict(row) for row in artifacts]
    bad_run_timestamp[1]["run_timestamp_utc"] = "not-a-date"
    bad_run_timestamp_gate = cst_touchstone_solver_ready_manifest_gate(
        bad_run_timestamp,
        max_created_run_skew_s=60.0,
        require_execution_metadata=True,
    )
    assert bad_run_timestamp_gate["status"] == "needs_attention"
    assert bad_run_timestamp_gate["checks"]["run_timestamp_utc_parseable_when_present"] is False
    assert bad_run_timestamp_gate["checks"]["created_run_timestamp_skew_recorded"] is False

    stale_selected_frequency = [dict(row) for row in artifacts]
    stale_selected_frequency[2]["selected_frequency_hz"] = 2.40e9
    stale_frequency_gate = cst_touchstone_solver_ready_manifest_gate(
        stale_selected_frequency,
        expected_design_frequency_hz=2.45e9,
    )
    assert stale_frequency_gate["status"] == "needs_attention"
    assert stale_frequency_gate["checks"]["selected_frequency_matches_design"] is False
    assert stale_frequency_gate["checks"]["selected_frequency_matches_expected_design"] is False

    stale_touchstone_row_index = [dict(row) for row in artifacts]
    stale_touchstone_row_index[3]["row_index"] = 1
    stale_touchstone_row_index_gate = cst_touchstone_solver_ready_manifest_gate(stale_touchstone_row_index)
    assert stale_touchstone_row_index_gate["status"] == "needs_attention"
    assert stale_touchstone_row_index_gate["checks"]["touchstone_row_index_matches_design_row"] is False
    assert stale_touchstone_row_index_gate["touchstone_row_index_mismatch_rows"] == [4]

    stale_touchstone_row_frequency = [dict(row) for row in artifacts]
    stale_touchstone_row_frequency[3]["row_frequency_hz"] = 2.40e9
    stale_touchstone_row_frequency_gate = cst_touchstone_solver_ready_manifest_gate(stale_touchstone_row_frequency)
    assert stale_touchstone_row_frequency_gate["status"] == "needs_attention"
    assert stale_touchstone_row_frequency_gate["checks"]["touchstone_row_frequency_matches_selected_frequency"] is False
    assert stale_touchstone_row_frequency_gate["touchstone_row_frequency_mismatch_rows"] == [4]

    stale_renormalized_z0 = [dict(row) for row in artifacts]
    stale_renormalized_z0[2]["renormalized_reference_impedance_ohm"] = 75.0
    stale_renormalized_z0_gate = cst_touchstone_solver_ready_manifest_gate(
        stale_renormalized_z0,
        expected_renormalized_reference_impedance_ohm=50.0,
        expected_renormalization_method="not_renormalized_option_line_R50",
        expected_renormalization_artifact_id="touchstone_export_B_renorm_none_R50_v1",
    )
    assert stale_renormalized_z0_gate["status"] == "needs_attention"
    assert stale_renormalized_z0_gate["checks"]["renormalized_reference_impedance_consistent_when_present"] is False
    assert stale_renormalized_z0_gate["checks"]["expected_renormalized_reference_impedance_matches"] is False

    stale_renormalization_artifact = [dict(row) for row in artifacts]
    stale_renormalization_artifact[3]["renormalization_artifact_id"] = "touchstone_export_A_renorm_old"
    stale_renormalization_artifact_gate = cst_touchstone_solver_ready_manifest_gate(
        stale_renormalization_artifact,
        expected_renormalization_method="not_renormalized_option_line_R50",
        expected_renormalization_artifact_id="touchstone_export_B_renorm_none_R50_v1",
    )
    assert stale_renormalization_artifact_gate["status"] == "needs_attention"
    assert stale_renormalization_artifact_gate["checks"]["renormalization_artifact_id_consistent_when_present"] is False
    assert stale_renormalization_artifact_gate["checks"]["expected_renormalization_artifact_id_matches"] is False

    stale_deembedding_artifact = [dict(row) for row in artifacts]
    stale_deembedding_artifact[3]["deembedding_artifact_id"] = "touchstone_export_A_deembed_old"
    stale_deembedding_artifact_gate = cst_touchstone_solver_ready_manifest_gate(
        stale_deembedding_artifact,
        expected_deembedding_method="port_extension_reference_plane_shift",
        expected_deembedding_artifact_id="touchstone_export_B_deembed_connector_faces_v1",
        expected_deembedding_length_m=0.0125,
    )
    assert stale_deembedding_artifact_gate["status"] == "needs_attention"
    assert stale_deembedding_artifact_gate["checks"]["deembedding_artifact_id_consistent_when_present"] is False
    assert stale_deembedding_artifact_gate["checks"]["expected_deembedding_artifact_id_matches"] is False

    stale_deembedding_length = [dict(row) for row in artifacts]
    stale_deembedding_length[2]["deembedding_length_m"] = 0.009
    stale_deembedding_length_gate = cst_touchstone_solver_ready_manifest_gate(
        stale_deembedding_length,
        expected_deembedding_method="port_extension_reference_plane_shift",
        expected_deembedding_artifact_id="touchstone_export_B_deembed_connector_faces_v1",
        expected_deembedding_length_m=0.0125,
    )
    assert stale_deembedding_length_gate["status"] == "needs_attention"
    assert stale_deembedding_length_gate["checks"]["deembedding_length_consistent_when_present"] is False
    assert stale_deembedding_length_gate["checks"]["expected_deembedding_length_matches"] is False

    missing_deembedding_method = [dict(row) for row in artifacts]
    missing_deembedding_method[1].pop("deembedding_method")
    missing_deembedding_method_gate = cst_touchstone_solver_ready_manifest_gate(
        missing_deembedding_method,
        expected_deembedding_method="port_extension_reference_plane_shift",
        expected_deembedding_artifact_id="touchstone_export_B_deembed_connector_faces_v1",
        expected_deembedding_length_m=0.0125,
    )
    assert missing_deembedding_method_gate["status"] == "needs_attention"
    assert missing_deembedding_method_gate["checks"]["deembedding_method_recorded_when_expected"] is False

    stale_solver_setup = [dict(row) for row in artifacts]
    stale_solver_setup[2]["solver_setup_artifact_id"] = "cst_slot335_fd_solver_setup_old"
    stale_solver_setup_gate = cst_touchstone_solver_ready_manifest_gate(
        stale_solver_setup,
        expected_solver_setup_artifact_id="cst_slot335_fd_solver_setup_v1",
        expected_mesh_setup_artifact_id="cst_slot335_adaptive_mesh_setup_v1",
    )
    assert stale_solver_setup_gate["status"] == "needs_attention"
    assert stale_solver_setup_gate["checks"]["solver_setup_artifact_id_consistent_when_present"] is False
    assert stale_solver_setup_gate["checks"]["expected_solver_setup_artifact_id_matches"] is False
    assert stale_solver_setup_gate["checks"]["expected_mesh_setup_artifact_id_matches"] is True

    stale_mesh_setup = [dict(row) for row in artifacts]
    stale_mesh_setup[1]["mesh_setup_artifact_id"] = "cst_slot335_adaptive_mesh_setup_old"
    stale_mesh_setup_gate = cst_touchstone_solver_ready_manifest_gate(
        stale_mesh_setup,
        expected_solver_setup_artifact_id="cst_slot335_fd_solver_setup_v1",
        expected_mesh_setup_artifact_id="cst_slot335_adaptive_mesh_setup_v1",
    )
    assert stale_mesh_setup_gate["status"] == "needs_attention"
    assert stale_mesh_setup_gate["checks"]["mesh_setup_artifact_id_consistent_when_present"] is False
    assert stale_mesh_setup_gate["checks"]["expected_solver_setup_artifact_id_matches"] is True
    assert stale_mesh_setup_gate["checks"]["expected_mesh_setup_artifact_id_matches"] is False

    missing_mesh_setup = [dict(row) for row in artifacts]
    missing_mesh_setup[1].pop("mesh_setup_artifact_id")
    missing_mesh_setup_gate = cst_touchstone_solver_ready_manifest_gate(
        missing_mesh_setup,
        expected_mesh_setup_artifact_id="cst_slot335_adaptive_mesh_setup_v1",
    )
    assert missing_mesh_setup_gate["status"] == "needs_attention"
    assert missing_mesh_setup_gate["checks"]["mesh_setup_artifact_id_recorded_when_expected"] is False

    stale_port_definition = [dict(row) for row in artifacts]
    stale_port_definition[0]["port_definition_artifact_id"] = "cst_slot351_old_port_definition"
    stale_port_definition_gate = cst_touchstone_solver_ready_manifest_gate(
        stale_port_definition,
        expected_port_definition_artifact_id="cst_slot351_waveguide_port_definition_v1",
        expected_excitation_setup_artifact_id="cst_slot351_unit_power_port_excitation_v1",
    )
    assert stale_port_definition_gate["status"] == "needs_attention"
    assert stale_port_definition_gate["checks"]["port_definition_artifact_id_consistent_when_present"] is False
    assert stale_port_definition_gate["checks"]["expected_port_definition_artifact_id_matches"] is False
    assert stale_port_definition_gate["checks"]["expected_excitation_setup_artifact_id_matches"] is True

    stale_excitation_setup = [dict(row) for row in artifacts]
    stale_excitation_setup[3]["excitation_setup_artifact_id"] = "cst_slot351_voltage_wave_excitation_old"
    stale_excitation_setup_gate = cst_touchstone_solver_ready_manifest_gate(
        stale_excitation_setup,
        expected_port_definition_artifact_id="cst_slot351_waveguide_port_definition_v1",
        expected_excitation_setup_artifact_id="cst_slot351_unit_power_port_excitation_v1",
    )
    assert stale_excitation_setup_gate["status"] == "needs_attention"
    assert stale_excitation_setup_gate["checks"]["excitation_setup_artifact_id_consistent_when_present"] is False
    assert stale_excitation_setup_gate["checks"]["expected_port_definition_artifact_id_matches"] is True
    assert stale_excitation_setup_gate["checks"]["expected_excitation_setup_artifact_id_matches"] is False

    missing_excitation_setup = [dict(row) for row in artifacts]
    missing_excitation_setup[1].pop("excitation_setup_artifact_id")
    missing_excitation_setup_gate = cst_touchstone_solver_ready_manifest_gate(
        missing_excitation_setup,
        expected_excitation_setup_artifact_id="cst_slot351_unit_power_port_excitation_v1",
    )
    assert missing_excitation_setup_gate["status"] == "needs_attention"
    assert missing_excitation_setup_gate["checks"]["excitation_setup_artifact_id_recorded_when_expected"] is False


def test_touchstone_row_solver_ready_preflight_bundles_match_and_passivity():
    row = {
        "frequency": 1.0,
        "s11": [0.05, 0.0],
        "s21": [0.80, -10.0],
        "s12": [0.80, -10.0],
        "s22": [0.05, 0.0],
    }
    gate = touchstone_row_solver_ready_preflight_gate(
        row,
        data_format="MA",
        z0=50.0,
        return_loss_min_db=20.0,
        vswr_max=1.2,
    )

    assert gate["policy"] == "touchstone_row_solver_ready_preflight_gate"
    assert gate["status"] == "ok"
    assert gate["sparameter_health"]["status"] == "ok"
    assert gate["port_match"]["return_loss_db"] == pytest.approx(26.020599913279625)
    assert gate["checks"]["touchstone_format_recorded"] is True
    assert gate["checks"]["reference_impedance_recorded"] is True
    assert gate["checks"]["return_loss_limit_ok"] is True
    assert gate["checks"]["vswr_limit_ok"] is True

    active = touchstone_row_solver_ready_preflight_gate(
        {"frequency": 1.0, "s11": [0.05, 0.0], "s21": [1.2, 0.0], "s12": [1.2, 0.0], "s22": [0.05, 0.0]},
        data_format="MA",
    )
    assert active["status"] == "needs_attention"
    assert active["checks"]["sparameter_passivity_ok"] is False


def test_touchstone_frequency_grid_interpolation_gate_requires_design_bracket():
    gate = touchstone_frequency_grid_interpolation_gate(
        [0.95e9, 0.99e9, 1.01e9, 1.05e9],
        1.0e9,
        max_relative_spacing=0.03,
    )

    assert gate["policy"] == "touchstone_frequency_grid_interpolation_gate"
    assert gate["status"] == "ok"
    assert gate["lower_index"] == 1
    assert gate["upper_index"] == 2
    assert gate["bracket_gap_rel"] == pytest.approx(0.02)
    assert gate["checks"]["design_frequency_bracketed"] is True
    assert gate["checks"]["design_spacing_ok"] is True

    exact = touchstone_frequency_grid_interpolation_gate([0.9e9, 1.0e9, 1.1e9], 1.0e9)
    assert exact["status"] == "ok"
    assert exact["lower_index"] == 1
    assert exact["upper_index"] == 1
    assert exact["bracket_gap_hz"] == pytest.approx(0.0)

    coarse = touchstone_frequency_grid_interpolation_gate(
        [0.8e9, 1.2e9],
        1.0e9,
        max_relative_spacing=0.05,
    )
    assert coarse["status"] == "needs_attention"
    assert coarse["checks"]["design_frequency_bracketed"] is True
    assert coarse["checks"]["design_spacing_ok"] is False

    outside = touchstone_frequency_grid_interpolation_gate([0.8e9, 0.9e9], 1.0e9)
    assert outside["status"] == "needs_attention"
    assert outside["checks"]["design_frequency_bracketed"] is False

    with pytest.raises(ValueError):
        touchstone_frequency_grid_interpolation_gate([0.9e9, 0.9e9, 1.0e9], 1.0e9)


def test_touchstone_frequency_unit_normalization_gate_keeps_raw_unit_and_selected_row_together():
    gate = touchstone_frequency_unit_normalization_gate(
        [0.95, 1.0, 1.05],
        frequency_unit="GHz",
        design_frequency=1.0,
        expected_frequency_unit="GHz",
        expected_design_frequency_hz=1.0e9,
        selected_row_index=1,
        max_relative_spacing=0.06,
    )

    assert gate["policy"] == "touchstone_frequency_unit_normalization_gate"
    assert gate["status"] == "ok"
    assert gate["frequency_hz"] == [0.95e9, 1.0e9, 1.05e9]
    assert gate["selected_frequency_hz"] == pytest.approx(1.0e9)
    assert gate["grid_contract"]["checks"]["design_frequency_bracketed"] is True
    assert gate["checks"]["frequency_unit_matches_expected"] is True
    assert gate["checks"]["selected_row_matches_design_frequency"] is True

    design_in_mhz = touchstone_frequency_unit_normalization_gate(
        [0.95, 1.0, 1.05],
        frequency_unit="GHz",
        design_frequency=1000.0,
        design_frequency_unit="MHz",
        expected_design_frequency_hz=1.0e9,
        selected_row_index=1,
        max_relative_spacing=0.06,
    )
    assert design_in_mhz["status"] == "ok"
    assert design_in_mhz["design_frequency_hz"] == pytest.approx(1.0e9)

    wrong_unit = touchstone_frequency_unit_normalization_gate(
        [0.95, 1.0, 1.05],
        frequency_unit="GHz",
        design_frequency=1.0,
        expected_frequency_unit="Hz",
        selected_row_index=1,
        max_relative_spacing=0.06,
    )
    assert wrong_unit["status"] == "needs_attention"
    assert wrong_unit["checks"]["frequency_unit_matches_expected"] is False

    stale_index = touchstone_frequency_unit_normalization_gate(
        [0.95, 1.0, 1.05],
        frequency_unit="GHz",
        design_frequency=1.0,
        selected_row_index=0,
        max_relative_spacing=0.06,
    )
    assert stale_index["status"] == "needs_attention"
    assert stale_index["checks"]["selected_row_matches_design_frequency"] is False

    sparse = touchstone_frequency_unit_normalization_gate(
        [0.8, 1.2],
        frequency_unit="GHz",
        design_frequency=1.0,
        max_relative_spacing=0.05,
    )
    assert sparse["status"] == "needs_attention"
    assert sparse["grid_contract"]["checks"]["design_spacing_ok"] is False


def test_shared_solver_session_health_gate_separates_reuse_from_physics():
    gate = shared_solver_session_health_gate(
        connected=True,
        api_visible=True,
        discovered_engines=["MATLAB_10416"],
        matlab_engine_find_matlab=["MATLAB_10416"],
        shared_engine_name="MATLAB_10416",
        livelink_matlab_pid=10416,
        status="already-connected",
        started_new_process=False,
        killed_process=False,
        direct_discovery_status="no MATLAB session discovered",
        passive_diagnostic_verdict="livelink-matlab-present",
        version="R2026a Update 3",
        model_tags=["cc_ht_stack_probe_73"],
        passive_server_pid=12284,
        passive_matlab_pid=9748,
        passive_worker_pid=10416,
        passive_matlab_parent_pid=12284,
        passive_worker_parent_pid=9748,
        target_port=2036,
        established_connection_count=1,
        shared_engine_eval="artifacts/livelink/eval_shared_matlab.py",
        livelink_out_fields=["connected", "port", "status", "reason"],
        matlab_version_source="version and version('-release')",
        solver_version_source="ModelUtil.getComsolVersion()",
        passive_diagnostic_timestamp="2026-06-30T16:41:41+09:00",
        passive_machine_policy="Passive inspection only: no COMSOL TCP probes, no MATLAB/COMSOL start, no kill.",
        passive_port_owner_pid=12284,
        shared_engine_eval_status="ok",
        previous_shared_engine_eval_status="timeout",
        matlab_process_count=5,
        matlab_mcp_server_count=20,
        livelink_candidate_count=1,
        session_selection_basis="parent process chain plus established port connection",
    )

    assert gate["policy"] == "shared_solver_session_health_gate"
    assert gate["status_label"] == "ok"
    assert gate["matlab_engine_find_matlab"] == ["MATLAB_10416"]
    assert gate["shared_engine_name"] == "MATLAB_10416"
    assert gate["livelink_matlab_pid"] == 10416
    assert gate["checks"]["session_connected"] is True
    assert gate["checks"]["api_visible"] is True
    assert gate["checks"]["engine_discovered"] is True
    assert gate["checks"]["matlab_engine_find_matlab_recorded"] is True
    assert gate["checks"]["selected_shared_engine_visible_in_find_matlab"] is True
    assert gate["checks"]["status_allows_reuse"] is True
    assert gate["checks"]["started_no_new_process"] is True
    assert gate["checks"]["direct_discovery_false_negative_reconciled"] is True
    assert gate["checks"]["direct_discovery_false_negative_has_selected_engine"] is True
    assert gate["checks"]["direct_discovery_false_negative_has_find_matlab_engine"] is True
    assert gate["checks"]["direct_discovery_false_negative_has_ok_shared_eval"] is True
    assert gate["checks"]["shared_engine_name_recorded"] is True
    assert gate["checks"]["shared_engine_name_discovered"] is True
    assert gate["checks"]["shared_engine_name_matches_pid"] is True
    assert gate["checks"]["livelink_matlab_pid_matches_worker_pid"] is True
    assert gate["passive_diagnostic_verdict"] == "livelink-matlab-present"
    assert gate["version"] == "R2026a Update 3"
    assert gate["model_tags"] == ["cc_ht_stack_probe_73"]
    assert gate["passive_server_pid"] == 12284
    assert gate["passive_matlab_pid"] == 9748
    assert gate["passive_worker_pid"] == 10416
    assert gate["passive_matlab_parent_pid"] == 12284
    assert gate["passive_worker_parent_pid"] == 9748
    assert gate["target_port"] == 2036
    assert gate["established_connection_count"] == 1
    assert gate["checks"]["passive_session_evidence_complete"] is True
    assert gate["checks"]["shared_engine_eval_recorded"] is True
    assert gate["checks"]["livelink_core_fields_present"] is True
    assert gate["checks"]["livelink_version_field_not_required"] is True
    assert gate["checks"]["matlab_version_source_recorded"] is True
    assert gate["checks"]["solver_version_source_recorded"] is True
    assert gate["checks"]["passive_no_tcp_probe_policy_recorded"] is True
    assert gate["checks"]["passive_matlab_parent_pid_recorded"] is True
    assert gate["checks"]["livelink_matlab_parent_is_server"] is True
    assert gate["checks"]["passive_worker_parent_pid_recorded"] is True
    assert gate["checks"]["worker_parent_is_livelink_matlab"] is True
    assert gate["checks"]["target_port_owned_by_server_pid"] is True
    assert gate["checks"]["shared_engine_eval_status_recorded"] is True
    assert gate["shared_engine_eval_status"] == "ok"
    assert gate["previous_shared_engine_eval_status"] == "timeout"
    assert gate["checks"]["previous_shared_engine_eval_status_recorded"] is True
    assert gate["checks"]["shared_engine_eval_recovered_after_timeout"] is True
    assert gate["checks"]["multiple_matlab_processes_do_not_create_ambiguity"] is True
    assert gate["checks"]["selection_uses_parent_or_port_chain"] is True
    assert gate["checks"]["model_tags_are_introspection_only"] is True
    assert "preflight" in " ".join(gate["notes"])
    assert "model tags" in " ".join(gate["notes"])

    bad = shared_solver_session_health_gate(
        connected=True,
        api_visible=False,
        discovered_engines=[],
        status="already-connected",
        started_new_process=True,
    )
    assert bad["status_label"] == "needs_attention"
    assert bad["checks"]["api_visible"] is False
    assert bad["checks"]["engine_discovered"] is False
    assert bad["checks"]["started_no_new_process"] is False

    incomplete = shared_solver_session_health_gate(
        connected=True,
        api_visible=True,
        discovered_engines=["MATLAB_10416"],
        matlab_engine_find_matlab=["MATLAB_10416"],
        status="already-connected",
        passive_diagnostic_verdict="livelink-matlab-present",
        passive_server_pid=12284,
        target_port=2036,
    )
    assert incomplete["status_label"] == "needs_attention"
    assert incomplete["checks"]["passive_session_evidence_complete"] is False
    assert incomplete["checks"]["passive_server_pid_recorded"] is True
    assert incomplete["checks"]["passive_worker_pid_recorded"] is False

    false_negative_without_eval = shared_solver_session_health_gate(
        connected=True,
        api_visible=True,
        discovered_engines=["MATLAB_10416"],
        shared_engine_name="MATLAB_10416",
        livelink_matlab_pid=10416,
        status="already-connected",
        direct_discovery_status="no MATLAB session discovered",
        passive_diagnostic_verdict="livelink-matlab-present",
        passive_worker_pid=10416,
    )
    assert false_negative_without_eval["status_label"] == "needs_attention"
    assert false_negative_without_eval["checks"]["direct_discovery_false_negative_reconciled"] is True
    assert false_negative_without_eval["checks"]["direct_discovery_false_negative_has_selected_engine"] is True
    assert false_negative_without_eval["checks"]["direct_discovery_false_negative_has_ok_shared_eval"] is False

    false_negative_without_find_matlab_identity = shared_solver_session_health_gate(
        connected=True,
        api_visible=True,
        discovered_engines=["MATLAB_10416"],
        matlab_engine_find_matlab=["MATLAB_77777"],
        shared_engine_name="MATLAB_10416",
        livelink_matlab_pid=10416,
        status="already-connected",
        direct_discovery_status="no MATLAB session discovered",
        passive_diagnostic_verdict="livelink-matlab-present",
        passive_worker_pid=10416,
        shared_engine_eval="artifacts/livelink/eval_shared_matlab.py",
        shared_engine_eval_status="ok",
    )
    assert false_negative_without_find_matlab_identity["status_label"] == "needs_attention"
    assert false_negative_without_find_matlab_identity["checks"]["matlab_engine_find_matlab_recorded"] is True
    assert false_negative_without_find_matlab_identity["checks"]["selected_shared_engine_visible_in_find_matlab"] is False
    assert false_negative_without_find_matlab_identity["checks"]["direct_discovery_false_negative_reconciled"] is False
    assert false_negative_without_find_matlab_identity["checks"]["direct_discovery_false_negative_has_find_matlab_engine"] is False

    wrong_port_owner = shared_solver_session_health_gate(
        connected=True,
        api_visible=True,
        discovered_engines=["MATLAB_10416"],
        status="already-connected",
        passive_server_pid=12284,
        passive_matlab_pid=9748,
        passive_worker_pid=10416,
        target_port=2036,
        established_connection_count=1,
        shared_engine_eval="artifacts/livelink/eval_shared_matlab.py",
        passive_port_owner_pid=99999,
    )
    assert wrong_port_owner["status_label"] == "needs_attention"
    assert wrong_port_owner["checks"]["target_port_owned_by_server_pid"] is False

    wrong_worker_parent = shared_solver_session_health_gate(
        connected=True,
        api_visible=True,
        discovered_engines=["MATLAB_10416"],
        status="already-connected",
        passive_server_pid=12284,
        passive_matlab_pid=9748,
        passive_worker_pid=10416,
        passive_matlab_parent_pid=12284,
        passive_worker_parent_pid=12284,
        target_port=2036,
        established_connection_count=1,
        shared_engine_eval="artifacts/livelink/eval_shared_matlab.py",
    )
    assert wrong_worker_parent["status_label"] == "needs_attention"
    assert wrong_worker_parent["checks"]["livelink_matlab_parent_is_server"] is True
    assert wrong_worker_parent["checks"]["worker_parent_is_livelink_matlab"] is False

    bad_version_assumption = shared_solver_session_health_gate(
        connected=True,
        api_visible=True,
        discovered_engines=["MATLAB_10416"],
        status="already-connected",
        livelink_out_fields=["connected", "port", "status", "reason"],
    )
    assert bad_version_assumption["status_label"] == "needs_attention"
    assert bad_version_assumption["checks"]["livelink_version_field_not_required"] is False

    timed_out_eval = shared_solver_session_health_gate(
        connected=True,
        api_visible=False,
        discovered_engines=["MATLAB_10416"],
        status="already-connected",
        started_new_process=False,
        killed_process=False,
        passive_diagnostic_verdict="livelink-matlab-present",
        passive_server_pid=12284,
        passive_matlab_pid=9748,
        passive_worker_pid=10416,
        target_port=2036,
        established_connection_count=1,
        shared_engine_eval="artifacts/livelink/eval_shared_matlab.py",
        shared_engine_eval_status="timeout",
        shared_engine_eval_timeout_s=20,
        shared_engine_eval_timeout_mode="child-process",
    )
    assert timed_out_eval["status_label"] == "needs_attention"
    assert timed_out_eval["checks"]["api_visible"] is False
    assert timed_out_eval["checks"]["shared_engine_eval_status_recorded"] is True
    assert timed_out_eval["checks"]["shared_engine_eval_timeout_recorded"] is True
    assert timed_out_eval["checks"]["shared_engine_eval_timeout_mode_recorded"] is True
    assert timed_out_eval["checks"]["shared_engine_eval_timeout_is_diagnostic"] is True
    assert timed_out_eval["shared_engine_eval_timeout_s"] == 20.0

    unrecovered_timeout = shared_solver_session_health_gate(
        connected=True,
        api_visible=False,
        discovered_engines=["MATLAB_10416"],
        status="already-connected",
        shared_engine_eval_status="timeout",
        previous_shared_engine_eval_status="timeout",
        shared_engine_eval_timeout_s=20,
        shared_engine_eval_timeout_mode="child-process",
    )
    assert unrecovered_timeout["status_label"] == "needs_attention"
    assert unrecovered_timeout["checks"]["previous_shared_engine_eval_status_recorded"] is True
    assert unrecovered_timeout["checks"]["shared_engine_eval_recovered_after_timeout"] is False

    ambiguous_process_selection = shared_solver_session_health_gate(
        connected=True,
        api_visible=True,
        discovered_engines=["MATLAB_10416"],
        status="already-connected",
        matlab_process_count=5,
        matlab_mcp_server_count=20,
        livelink_candidate_count=1,
        session_selection_basis="first MATLAB process in task list",
    )
    assert ambiguous_process_selection["status_label"] == "needs_attention"
    assert ambiguous_process_selection["checks"]["selection_uses_parent_or_port_chain"] is False
    assert ambiguous_process_selection["checks"]["multiple_matlab_processes_do_not_create_ambiguity"] is False

    wrong_engine_pid = shared_solver_session_health_gate(
        connected=True,
        api_visible=True,
        discovered_engines=["MATLAB_10416"],
        shared_engine_name="MATLAB_99999",
        livelink_matlab_pid=10416,
        status="already-connected",
        passive_worker_pid=10416,
    )
    assert wrong_engine_pid["status_label"] == "needs_attention"
    assert wrong_engine_pid["checks"]["shared_engine_name_discovered"] is False
    assert wrong_engine_pid["checks"]["shared_engine_name_matches_pid"] is False
    assert wrong_engine_pid["checks"]["livelink_matlab_pid_matches_worker_pid"] is True


def test_solver_submodel_boundary_handoff_gate_tracks_boundary_transfer_error():
    gate = solver_submodel_boundary_handoff_gate(
        parent_model_id="global_plate_bending_coarse_v1",
        parent_mesh_id="mesh_global_h2",
        submodel_region_id="zoom_region_tip_01",
        local_mesh_id="mesh_zoom_h8_adaptive",
        zoom_boundary_id="zoom_boundary_tip_01",
        boundary_trace_id="trace_global_to_zoom_01",
        boundary_condition_source="parent displacement and slope trace",
        boundary_transfer_quantity="displacement+slope",
        expected_boundary_transfer_quantity="displacement+slope",
        boundary_transfer_error_estimate=0.018,
        boundary_transfer_error_unit="relative",
        max_boundary_transfer_error=0.02,
        local_refinement_rule="single-pass adaptive refinement from parent error indicator",
        target_observable_id="tip_bending_moment",
    )

    assert gate["policy"] == "solver_submodel_boundary_handoff_gate"
    assert gate["status"] == "ok"
    assert gate["checks"]["parent_model_id_recorded"] is True
    assert gate["checks"]["submodel_region_id_recorded"] is True
    assert gate["checks"]["zoom_boundary_id_recorded"] is True
    assert gate["checks"]["boundary_transfer_error_nonnegative"] is True
    assert gate["checks"]["boundary_transfer_error_unit_recorded"] is True
    assert gate["checks"]["boundary_transfer_error_within_limit"] is True
    assert gate["checks"]["boundary_transfer_quantity_matches_expected"] is True
    assert gate["checks"]["parent_local_mesh_identity_separated"] is True
    assert gate["checks"]["boundary_handoff_not_value_only"] is True
    assert "boundary transfer" in " ".join(gate["notes"])

    missing_handoff = solver_submodel_boundary_handoff_gate(
        parent_model_id="global_plate_bending_coarse_v1",
        submodel_region_id="zoom_region_tip_01",
        zoom_boundary_id="",
        boundary_transfer_quantity="",
        boundary_transfer_error_estimate=None,
        local_refinement_rule="refined local mesh",
        target_observable_id="tip_bending_moment",
    )
    assert missing_handoff["status"] == "needs_attention"
    assert missing_handoff["checks"]["zoom_boundary_id_recorded"] is False
    assert missing_handoff["checks"]["boundary_transfer_quantity_recorded"] is False
    assert missing_handoff["checks"]["boundary_transfer_error_estimate_recorded"] is False
    assert missing_handoff["checks"]["boundary_handoff_not_value_only"] is False

    over_budget = solver_submodel_boundary_handoff_gate(
        parent_model_id="global_plate_bending_coarse_v1",
        submodel_region_id="zoom_region_tip_01",
        zoom_boundary_id="zoom_boundary_tip_01",
        boundary_transfer_quantity="displacement+slope",
        expected_boundary_transfer_quantity="magnetic_vector_potential_trace",
        boundary_transfer_error_estimate=0.031,
        boundary_transfer_error_unit="relative",
        max_boundary_transfer_error=0.02,
        local_refinement_rule="single-pass adaptive refinement from parent error indicator",
        target_observable_id="tip_bending_moment",
    )
    assert over_budget["status"] == "needs_attention"
    assert over_budget["checks"]["boundary_transfer_error_within_limit"] is False
    assert over_budget["checks"]["boundary_transfer_quantity_matches_expected"] is False


def test_solver_result_artifact_provenance_timing_gate_records_versions_dates_and_heavy_stages():
    artifact = {
        "schema": "cae-ai-lab.crossval.v1",
        "created_at_utc": "2026-06-30T04:00:00+00:00",
        "tool_slot": "COMSOL",
        "versions": {
            "solver": "COMSOL 6.4.0.378",
            "matlab": "R2026a",
            "mcp": "comsol-mcp 0.10.0",
        },
        "execution": {
            "run_date_utc": "2026-06-30T04:00:01Z",
            "shared_engine": "MATLAB_10416",
        },
        "result_output_schema_id": "matlab_fem_bem_result_table_v1",
        "result_output_columns": [
            "alpha",
            "trace_residual_norm",
            "solution_norm",
            "objective_value",
        ],
        "result_output_units": {
            "alpha": "1",
            "trace_residual_norm": "1",
            "solution_norm": "1",
            "objective_value": "1",
        },
        "timing_breakdown_s": {
            "attach_livelink": 0.12,
            "model_build": 0.30,
            "solve": 1.50,
            "postprocess": 0.20,
        },
    }

    gate = solver_result_artifact_provenance_timing_gate(
        artifact,
        required_versions=("solver", "matlab"),
        required_timing_stages=("attach_livelink", "solve", "postprocess"),
        min_timing_stages=4,
        expected_created_at_utc="2026-06-30T04:00:00+00:00",
        expected_run_date_utc="2026-06-30T04:00:01Z",
        max_created_run_skew_s=5.0,
        expected_execution_session_id="MATLAB_10416",
        require_execution_session_id=True,
        expected_result_output_schema_id="matlab_fem_bem_result_table_v1",
        expected_result_output_columns=(
            "alpha",
            "trace_residual_norm",
            "solution_norm",
            "objective_value",
        ),
        expected_result_output_units={
            "alpha": "1",
            "trace_residual_norm": "1",
            "solution_norm": "1",
            "objective_value": "1",
        },
        require_result_output_schema=True,
    )

    assert gate["status"] == "ok"
    assert gate["policy"] == "solver_result_artifact_provenance_timing_gate"
    assert gate["missing_versions"] == []
    assert gate["missing_timing_stages"] == []
    assert gate["dominant_timing_stages"][0]["stage"] == "solve"
    assert gate["checks"]["timing_stage_count_reasonable"] is True
    assert gate["total_recorded_timing_s"] == pytest.approx(2.12)
    assert gate["checks"]["created_at_utc_recorded"] is True
    assert gate["checks"]["created_at_utc_parseable"] is True
    assert gate["checks"]["run_date_utc_recorded_when_required"] is True
    assert gate["checks"]["run_date_utc_parseable"] is True
    assert gate["checks"]["expected_created_at_utc_matches"] is True
    assert gate["checks"]["expected_run_date_utc_matches"] is True
    assert gate["checks"]["execution_session_id_recorded_when_required"] is True
    assert gate["checks"]["expected_execution_session_id_matches"] is True
    assert gate["execution_session_id"] == "MATLAB_10416"
    assert gate["checks"]["result_output_schema_id_recorded_when_required"] is True
    assert gate["checks"]["result_output_columns_recorded_when_required"] is True
    assert gate["checks"]["result_output_units_recorded_when_required"] is True
    assert gate["checks"]["expected_result_output_schema_id_matches"] is True
    assert gate["checks"]["expected_result_output_columns_match"] is True
    assert gate["checks"]["expected_result_output_units_match"] is True
    assert gate["result_output_schema_id"] == "matlab_fem_bem_result_table_v1"
    assert gate["result_output_columns"] == [
        "alpha",
        "trace_residual_norm",
        "solution_norm",
        "objective_value",
    ]
    assert gate["result_output_units"] == {
        "alpha": "1",
        "trace_residual_norm": "1",
        "solution_norm": "1",
        "objective_value": "1",
    }
    assert gate["checks"]["created_run_timestamp_skew_within_limit"] is True
    assert gate["created_run_skew_s"] == pytest.approx(1.0)

    missing_solver = dict(artifact)
    missing_solver["versions"] = {"matlab": "R2026a"}
    missing_gate = solver_result_artifact_provenance_timing_gate(
        missing_solver,
        required_versions=("solver", "matlab"),
        required_timing_stages=("attach_livelink", "solve", "postprocess"),
        min_timing_stages=4,
    )
    assert missing_gate["status"] == "needs_attention"
    assert missing_gate["checks"]["required_versions_recorded"] is False

    sparse_timing = dict(artifact)
    sparse_timing["timing_breakdown_s"] = {"solve": 1.5}
    sparse_gate = solver_result_artifact_provenance_timing_gate(
        sparse_timing,
        required_versions=("solver", "matlab"),
        required_timing_stages=("attach_livelink", "solve", "postprocess"),
        min_timing_stages=4,
    )
    assert sparse_gate["status"] == "needs_attention"
    assert sparse_gate["checks"]["timing_breakdown_recorded"] is False
    assert sparse_gate["checks"]["required_timing_stages_recorded"] is False

    bad_date = dict(artifact)
    bad_date["created_at_utc"] = "not-a-date"
    bad_date_gate = solver_result_artifact_provenance_timing_gate(bad_date, required_versions=("solver",))
    assert bad_date_gate["status"] == "needs_attention"
    assert bad_date_gate["checks"]["created_at_utc_parseable"] is False

    missing_run_date = dict(artifact)
    missing_run_date.pop("execution")
    missing_run_date_gate = solver_result_artifact_provenance_timing_gate(
        missing_run_date,
        required_versions=("solver", "matlab"),
        required_timing_stages=("attach_livelink", "solve", "postprocess"),
        min_timing_stages=4,
    )
    assert missing_run_date_gate["status"] == "needs_attention"
    assert missing_run_date_gate["checks"]["run_date_utc_recorded_when_required"] is False

    missing_session = dict(artifact)
    missing_session["execution"] = {"run_date_utc": "2026-06-30T04:00:01Z"}
    missing_session_gate = solver_result_artifact_provenance_timing_gate(
        missing_session,
        required_versions=("solver", "matlab"),
        required_timing_stages=("attach_livelink", "solve", "postprocess"),
        min_timing_stages=4,
        require_execution_session_id=True,
    )
    assert missing_session_gate["status"] == "needs_attention"
    assert missing_session_gate["checks"]["execution_session_id_recorded_when_required"] is False

    stale_session = dict(artifact)
    stale_session["execution"] = dict(artifact["execution"])
    stale_session["execution"]["shared_engine"] = "MATLAB_99999"
    stale_session_gate = solver_result_artifact_provenance_timing_gate(
        stale_session,
        expected_execution_session_id="MATLAB_10416",
    )
    assert stale_session_gate["status"] == "needs_attention"
    assert stale_session_gate["checks"]["expected_execution_session_id_matches"] is False

    stale_result_schema = dict(artifact)
    stale_result_schema["result_output_schema_id"] = "matlab_scalar_trace_residual_v0"
    stale_result_schema["result_output_columns"] = ["alpha", "trace_residual_norm"]
    stale_result_schema["result_output_units"] = {
        "alpha": "1",
        "trace_residual_norm": "dB",
    }
    stale_result_schema_gate = solver_result_artifact_provenance_timing_gate(
        stale_result_schema,
        required_versions=("solver", "matlab"),
        required_timing_stages=("attach_livelink", "solve", "postprocess"),
        min_timing_stages=4,
        expected_result_output_schema_id="matlab_fem_bem_result_table_v1",
        expected_result_output_columns=(
            "alpha",
            "trace_residual_norm",
            "solution_norm",
            "objective_value",
        ),
        expected_result_output_units={
            "alpha": "1",
            "trace_residual_norm": "1",
            "solution_norm": "1",
            "objective_value": "1",
        },
        require_result_output_schema=True,
    )
    assert stale_result_schema_gate["status"] == "needs_attention"
    assert stale_result_schema_gate["checks"]["expected_result_output_schema_id_matches"] is False
    assert stale_result_schema_gate["checks"]["expected_result_output_columns_match"] is False
    assert stale_result_schema_gate["checks"]["expected_result_output_units_match"] is False

    missing_result_schema = dict(artifact)
    missing_result_schema.pop("result_output_schema_id")
    missing_result_schema.pop("result_output_columns")
    missing_result_schema.pop("result_output_units")
    missing_result_schema_gate = solver_result_artifact_provenance_timing_gate(
        missing_result_schema,
        required_versions=("solver", "matlab"),
        required_timing_stages=("attach_livelink", "solve", "postprocess"),
        min_timing_stages=4,
        require_result_output_schema=True,
    )
    assert missing_result_schema_gate["status"] == "needs_attention"
    assert (
        missing_result_schema_gate["checks"][
            "result_output_schema_id_recorded_when_required"
        ]
        is False
    )
    assert missing_result_schema_gate["checks"]["result_output_columns_recorded_when_required"] is False
    assert missing_result_schema_gate["checks"]["result_output_units_recorded_when_required"] is False

    stale_run_date = dict(artifact)
    stale_run_date["execution"] = dict(artifact["execution"])
    stale_run_date["execution"]["run_date_utc"] = "2026-06-30T04:10:00Z"
    stale_run_date_gate = solver_result_artifact_provenance_timing_gate(
        stale_run_date,
        required_versions=("solver", "matlab"),
        required_timing_stages=("attach_livelink", "solve", "postprocess"),
        min_timing_stages=4,
        expected_run_date_utc="2026-06-30T04:00:01Z",
        max_created_run_skew_s=5.0,
    )
    assert stale_run_date_gate["status"] == "needs_attention"
    assert stale_run_date_gate["checks"]["expected_run_date_utc_matches"] is False
    assert stale_run_date_gate["checks"]["created_run_timestamp_skew_within_limit"] is False

    noisy_timing = dict(artifact)
    noisy_timing["timing_breakdown_s"] = {
        "attach_livelink": 0.12,
        "model_build": 0.30,
        "mesh": 0.40,
        "assembly": 0.50,
        "solve": 1.50,
        "postprocess": 0.20,
    }
    noisy_gate = solver_result_artifact_provenance_timing_gate(
        noisy_timing,
        required_versions=("solver", "matlab"),
        required_timing_stages=("attach_livelink", "solve", "postprocess"),
        min_timing_stages=4,
    )
    assert noisy_gate["status"] == "needs_attention"
    assert noisy_gate["checks"]["timing_stage_count_reasonable"] is False
    assert noisy_gate["timing_stage_count"] == 6
    assert len(noisy_gate["dominant_timing_stages"]) == 4


def test_source_native_seed_queue_gate_separates_preflight_from_crossval_learning():
    queue = {
        "created_at": "2026-07-03T00:00:00Z",
        "rounds": 2,
        "total_slots": 4,
        "rotation": ["COMSOL", "FEMM"],
        "slots": [
            {
                "tool": "COMSOL",
                "source_native_example": "public-example:acdc/coaxial-cable",
                "source_type": "public_doc",
                "lesson_axis": "coaxial EC/ES duality",
                "intended_validation": "radia-ngsolve public analogue",
                "lap": 1,
                "slot_id": "seed-1",
                "status": "queued_source_native_preflight",
            },
            {
                "tool": "FEMM",
                "source_native_example": "local-fixture:femm/examples/magnetostatic",
                "source_type": "local_path",
                "local_exists": True,
                "lesson_axis": "magnetization sign convention",
                "intended_validation": "public sign gate",
                "lap": 1,
                "slot_id": "seed-2",
                "status": "queued_source_native_preflight",
            },
            {
                "tool": "COMSOL",
                "source_native_example": "public-example:acdc/current-carrying-wire",
                "source_type": "public_doc",
                "lesson_axis": "parallel-wire force sign",
                "intended_validation": "analytic wire force gate",
                "lap": 2,
                "slot_id": "seed-3",
                "status": "queued_source_native_preflight",
            },
            {
                "tool": "FEMM",
                "source_native_example": "upstream-example:femm/coilgun",
                "source_type": "upstream_example",
                "lesson_axis": "force/coenergy handoff",
                "intended_validation": "public coenergy derivative gate",
                "lap": 2,
                "slot_id": "seed-4",
                "status": "queued_source_native_preflight",
            },
        ],
    }

    gate = source_native_seed_queue_gate(
        queue,
        expected_tools=("COMSOL", "FEMM"),
        expected_rounds=2,
        expected_total_slots=4,
        require_public_safe_sources=True,
    )
    assert gate["status"] == "ok"
    assert gate["learning_stage"] == "queued_not_learned"
    assert gate["checks"]["no_solver_or_learning_overclaim"] is True
    assert gate["tool_counts"] == {"COMSOL": 2, "FEMM": 2}

    overclaim = {**queue, "slots": [dict(slot) for slot in queue["slots"]]}
    overclaim["slots"][0]["status"] = "verified_crossval_passed"
    overclaim_gate = source_native_seed_queue_gate(overclaim)
    assert overclaim_gate["status"] == "needs_attention"
    assert overclaim_gate["checks"]["no_solver_or_learning_overclaim"] is False
    assert overclaim_gate["solver_claim_slots"][0]["slot_id"] == "seed-1"

    missing_field = {**queue, "slots": [dict(slot) for slot in queue["slots"]]}
    missing_field["slots"][1]["lesson_axis"] = ""
    missing_gate = source_native_seed_queue_gate(missing_field)
    assert missing_gate["status"] == "needs_attention"
    assert missing_gate["checks"]["required_slot_fields_present"] is False
    assert missing_gate["missing_fields"][0]["missing"] == ["lesson_axis"]

    private_source = {**queue, "slots": [dict(slot) for slot in queue["slots"]]}
    private_source["slots"][0]["source_native_example"] = "internal://licensed/source"
    private_gate = source_native_seed_queue_gate(
        private_source,
        require_public_safe_sources=True,
    )
    assert private_gate["status"] == "needs_attention"
    assert private_gate["checks"]["public_safe_sources_when_required"] is False

    feedback_artifact = {
        "schema": "radia.crossval.v1",
        "tool_slot": "radia-mcp",
        "pass": True,
        "created_at_utc": "2026-07-03T00:01:00Z",
        "versions": {
            "solver": "source-native seed queue gate v1",
            "radia_mcp": "1.4.3",
        },
        "execution": {"run_date_utc": "2026-07-03T00:01:01Z"},
        "result_artifact_id": "source_native_seed_queue_gate_20260703",
        "result_output_schema_id": "source_native_seed_queue_gate_v1",
        "result_output_columns": ["tool", "queued_slots", "gate_status"],
        "result_output_units": {
            "tool": "1",
            "queued_slots": "1",
            "gate_status": "1",
        },
        "timing_breakdown_s": {"queue_gate": 0.01, "feedback_gate": 0.01},
        "learning_lanes": {"public": "verified", "source_tool": "verified"},
        "public_lesson": (
            "Source-native seed queues are replay material; radia-mcp may call "
            "them learned only after a promoted result artifact passes feedback."
        ),
        "learning_targets": [
            "radia-mcp: source_native_seed_queue_gate",
            "radia_ngsolve.slot_gates",
        ],
        "verification": {
            "public": "python -m pytest packages/radia-mcp/tests/test_loop_slot_gates.py -q -k source_native_seed_queue",
            "commands": [
                {
                    "command": "python -m pytest packages/radia-mcp/tests/test_loop_slot_gates.py -q -k source_native_seed_queue",
                    "result": "passed",
                }
            ],
        },
        "mcp_feedback": {
            "public_summary": (
                "Added a gate that accepts source-native loop seeds while "
                "blocking solver-learning overclaims."
            ),
            "encoded_targets": ["radia-mcp: source_native_seed_queue_gate"],
        },
        "next_slot_allowed": True,
    }
    feedback_gate = cross_validation_artifact_to_mcp_feedback_gate(
        feedback_artifact,
        require_replayable_verification_commands=True,
    )
    assert feedback_gate["status"] == "ok"
    assert feedback_gate["learning_stage"] == "learned"
    assert feedback_gate["provenance_gate_status"] == "ok"


def test_computed_reference_rows_gate_checks_real_result_rows():
    artifact = {
        "results": [
            {
                "name": "coax_line_inductance",
                "checks": [
                    {
                        "quantity": "external_inductance_per_m",
                        "computed": 2.7599084905172834e-7,
                        "reference": 2.7725887222397814e-7,
                        "unit": "H/m",
                        "rel_error": 0.004573426855842696,
                        "tolerance": 0.01,
                        "pass": True,
                    },
                    {
                        "quantity": "characteristic_impedance",
                        "computed": 82.92982911875659,
                        "reference": 83.17766166719343,
                        "unit": "ohm",
                        "rel_error": 0.0029795565716726068,
                        "tolerance": 0.01,
                        "pass": True,
                    },
                ],
            },
            {
                "name": "parallel_plate_pressure",
                "checks": [
                    {
                        "quantity": "pressure",
                        "computed": 1.1067734765999921,
                        "reference": 1.1067734766000001,
                        "unit": "Pa",
                        "rel_error": 7.222440676711394e-15,
                        "tolerance": 1.0e-3,
                        "pass": True,
                    }
                ],
            },
        ]
    }

    gate = computed_reference_rows_gate(
        artifact,
        max_global_rel_error=0.01,
    )
    assert gate["status"] == "ok"
    assert gate["row_count"] == 3
    assert gate["valid_row_count"] == 3
    assert gate["checks"]["row_errors_within_tolerance"] is True
    assert gate["checks"]["rel_error_matches_computed_reference"] is True
    assert gate["max_rel_error"] == pytest.approx(0.004573426855842696)

    too_loose_claim = {
        "rows": [
            {
                "case": "coax_line_inductance",
                "quantity": "external_inductance_per_m",
                "computed": 2.7599084905172834e-7,
                "reference": 2.7725887222397814e-7,
                "rel_error": 0.0,
                "tolerance": 0.01,
                "pass": True,
            }
        ]
    }
    mismatch_gate = computed_reference_rows_gate(too_loose_claim)
    assert mismatch_gate["status"] == "needs_attention"
    assert mismatch_gate["checks"]["rel_error_matches_computed_reference"] is False

    failing_row = {
        "rows": [
            {
                "case": "coax_line_inductance",
                "quantity": "external_inductance_per_m",
                "computed": 2.7599084905172834e-7,
                "reference": 2.7725887222397814e-7,
                "tolerance": 1.0e-4,
                "pass": False,
            }
        ]
    }
    fail_gate = computed_reference_rows_gate(failing_row)
    assert fail_gate["status"] == "needs_attention"
    assert fail_gate["checks"]["row_errors_within_tolerance"] is False
    assert fail_gate["checks"]["pass_flags_true_when_required"] is False

    missing = {"rows": [{"quantity": "pressure", "computed": 1.0, "reference": 1.0}]}
    missing_gate = computed_reference_rows_gate(missing)
    assert missing_gate["status"] == "needs_attention"
    assert missing_gate["checks"]["required_fields_present"] is False


def test_cross_validation_artifact_to_mcp_feedback_gate_requires_lesson_target_and_verification():
    artifact = {
        "schema": "radia.crossval.v1",
        "pass": True,
        "created_at_utc": "2026-07-02T12:00:00Z",
        "versions": {
            "solver": "radia 4.95.2",
            "radia_mcp": "0.82.0",
        },
        "execution": {"run_date_utc": "2026-07-02T12:00:01Z"},
        "result_artifact_id": "vim_hdiv_mmmm_20260702",
        "notebook": {
            "notebook_source_artifact_id": "vim_hdiv_mmmm_panel_ipynb_v1",
            "notebook_source_digest": "sha256:notebook-source-v1",
            "notebook_source_path": "docs/vim_hdiv_mmmm_panel.ipynb",
        },
        "result_output_schema_id": "vim_demag_table_v1",
        "result_output_columns": ["method", "dof", "wall_s", "rel_error"],
        "result_output_units": {
            "method": "1",
            "dof": "1",
            "wall_s": "s",
            "rel_error": "1",
        },
        "timing_breakdown_s": {"solve": 2.5, "postprocess": 0.3},
        "learning_lanes": {"public": "verified", "source_tool": "candidate"},
        "public_lesson": "HDiv/HDiv-VIM comparison artifacts must record timing and the public solver convention.",
        "learning_targets": ["radia-mcp", "radia_ngsolve.loop_learning"],
        "verification": {"public": "pytest packages/radia-mcp/tests/test_loop_slot_gates.py"},
    }

    gate = cross_validation_artifact_to_mcp_feedback_gate(
        artifact,
        require_notebook_source=True,
    )
    assert gate["status"] == "ok"
    assert gate["learning_stage"] == "learned"
    assert gate["source_tool_lane_status"] == "candidate"
    assert gate["result_artifact_id"] == "vim_hdiv_mmmm_20260702"
    assert gate["notebook_source_artifact_id"] == "vim_hdiv_mmmm_panel_ipynb_v1"
    assert gate["checks"]["public_lesson_recorded_when_required"] is True
    assert gate["checks"]["public_learning_target_recorded_when_required"] is True
    assert gate["checks"]["public_verification_recorded_when_required"] is True
    assert gate["checks"]["notebook_source_digest_recorded_when_required"] is True
    assert gate["provenance_gate_status"] == "ok"

    missing_lesson = dict(artifact)
    missing_lesson.pop("public_lesson")
    missing_lesson["notes"] = []
    missing_gate = cross_validation_artifact_to_mcp_feedback_gate(missing_lesson)
    assert missing_gate["status"] == "needs_attention"
    assert missing_gate["checks"]["public_lesson_recorded_when_required"] is False

    missing_target = dict(artifact)
    missing_target["learning_targets"] = ["private-tool-memory"]
    target_gate = cross_validation_artifact_to_mcp_feedback_gate(missing_target)
    assert target_gate["status"] == "needs_attention"
    assert target_gate["checks"]["public_learning_target_recorded_when_required"] is False

    missing_notebook_digest = {
        **artifact,
        "notebook": {
            **artifact["notebook"],
            "notebook_source_digest": "",
        },
    }
    notebook_gate = cross_validation_artifact_to_mcp_feedback_gate(
        missing_notebook_digest,
        require_notebook_source=True,
    )
    assert notebook_gate["status"] == "needs_attention"
    assert notebook_gate["checks"]["notebook_source_digest_recorded_when_required"] is False

    replayable = {
        **artifact,
        "verification": {
            "public": "pytest packages/radia-mcp/tests/test_loop_slot_gates.py",
            "commands": [
                {
                    "command": "python -m pytest packages/radia-mcp/tests/test_loop_slot_gates.py -q",
                    "result": "passed",
                }
            ],
        },
    }
    replayable_gate = cross_validation_artifact_to_mcp_feedback_gate(
        replayable,
        require_replayable_verification_commands=True,
    )
    assert replayable_gate["status"] == "ok"
    assert replayable_gate["checks"]["replay_command_recorded_when_required"] is True
    assert replayable_gate["checks"]["replay_commands_normalized_when_required"] is True
    assert replayable_gate["normalized_replay_commands"] == [
        "python -m pytest packages/radia-mcp/tests/test_loop_slot_gates.py -q"
    ]

    annotated_command = {
        **artifact,
        "verification": {
            "commands": [
                {
                    "command": "python -m pytest packages/radia-mcp/tests/test_loop_slot_gates.py -q -> passed"
                }
            ],
        },
    }
    annotated_gate = cross_validation_artifact_to_mcp_feedback_gate(
        annotated_command,
        require_replayable_verification_commands=True,
    )
    assert annotated_gate["status"] == "needs_attention"
    assert annotated_gate["checks"]["replay_command_recorded_when_required"] is True
    assert annotated_gate["checks"]["replay_commands_normalized_when_required"] is False


def test_owned_solver_model_tag_lifecycle_gate_requires_owned_cleanup_before_reuse():
    artifact = {
        "connected": True,
        "status": "already-connected",
        "artifact_id": "comsol_slot233_model_tag_identity",
        "model_tag": "cc_slot233_tag_identity_20260630_203333",
        "owned_model_tag_prefix": "cc_slot233_",
        "model_tag_owner": "codex_slot233",
        "created_present": True,
        "removed_after_probe": True,
        "preexisting_tags_preserved": True,
        "tags_before_count": 1,
        "tags_created_count": 2,
        "tags_after_count": 1,
        "tags_before": ["cc_ht_stack_probe_73"],
        "tags_after": ["cc_ht_stack_probe_73"],
        "started_new_matlab": False,
        "started_new_comsol": False,
        "killed_process": False,
    }

    gate = owned_solver_model_tag_lifecycle_gate(
        artifact,
        expected_artifact_id="comsol_slot233_model_tag_identity",
        expected_model_tag_prefix="cc_slot233_",
        expected_model_tag_owner="codex_slot233",
    )

    assert gate["status"] == "ok"
    assert gate["policy"] == "owned_solver_model_tag_lifecycle_gate"
    assert gate["checks"]["tag_was_created_and_visible"] is True
    assert gate["checks"]["owned_tag_removed_after_probe"] is True
    assert gate["checks"]["preexisting_tags_preserved"] is True
    assert gate["checks"]["tag_absent_after_cleanup"] is True
    assert gate["checks"]["created_count_increased_by_one"] is True
    assert gate["checks"]["after_count_restored"] is True
    assert gate["checks"]["expected_model_tag_owner_matches"] is True
    assert gate["checks"]["owned_prefix_absent_after_cleanup"] is True
    assert gate["owned_prefix_tags_after"] == []

    stale_tag = dict(artifact)
    stale_tag["tags_after"] = ["cc_ht_stack_probe_73", artifact["model_tag"]]
    stale_gate = owned_solver_model_tag_lifecycle_gate(stale_tag)
    assert stale_gate["status"] == "needs_attention"
    assert stale_gate["checks"]["tag_absent_after_cleanup"] is False

    stale_prefix = dict(artifact)
    stale_prefix["tags_after"] = ["cc_ht_stack_probe_73", "cc_slot233_abandoned_probe"]
    stale_prefix_gate = owned_solver_model_tag_lifecycle_gate(
        stale_prefix,
        expected_model_tag_prefix="cc_slot233_",
    )
    assert stale_prefix_gate["status"] == "needs_attention"
    assert stale_prefix_gate["checks"]["tag_absent_after_cleanup"] is True
    assert stale_prefix_gate["checks"]["owned_prefix_absent_after_cleanup"] is False
    assert stale_prefix_gate["owned_prefix_tags_after"] == ["cc_slot233_abandoned_probe"]

    missing_owner = dict(artifact)
    missing_owner.pop("model_tag_owner")
    missing_owner_gate = owned_solver_model_tag_lifecycle_gate(missing_owner)
    assert missing_owner_gate["status"] == "needs_attention"
    assert missing_owner_gate["checks"]["model_tag_owner_recorded"] is False

    wrong_counts = dict(artifact)
    wrong_counts["tags_created_count"] = 1
    wrong_counts_gate = owned_solver_model_tag_lifecycle_gate(wrong_counts)
    assert wrong_counts_gate["status"] == "needs_attention"
    assert wrong_counts_gate["checks"]["created_count_increased_by_one"] is False


def test_solver_result_table_metadata_gate_requires_columns_units_axis_and_rows():
    metadata = {
        "source": "COMSOL",
        "dataset": "dset1",
        "solution_tag": "sol1",
        "solution_artifact_id": "comsol_slot353_solution_data_sol1_v1",
        "solution_digest": "sha256:comsol_slot353_solution_data_sol1_v1",
        "sweep_axis_id": "freq_grid_slot361_v1",
        "sweep_axis_digest": "sha256:freq_grid_slot361_v1",
        "sweep_axis_row_count": 4,
        "parameter_set_artifact_id": "comsol_slot389_parameter_set_v1",
        "parameter_set_digest": "sha256:comsol-slot389-parameter-set-v1",
        "parameter_set_path": "validation/comsol/slot389_parameter_set.json",
        "objective_observable_id": "comsol_slot389_impedance_absorption_objective_v1",
        "objective_observable_family": "impedance_absorption_objective",
        "solver_configuration_artifact_id": "comsol_slot368_solver_config_sol1_v1",
        "solver_configuration_digest": "sha256:comsol_slot368_solver_config_sol1_v1",
        "solver_sequence_tag": "sol1",
        "linear_solver": "direct_pardiso",
        "relative_tolerance": 1.0e-6,
        "study_tag": "std1",
        "study_step_tag": "stat",
        "table_id": "tbl_probe_impedance",
        "selection_tags": ["bnd_probe"],
        "entity_dimensions": ["boundary"],
        "expressions": ["freq", "intop_bnd(acpr.p_t)", "acpr.Q_tot"],
        "operator_tags": ["intop_bnd"],
        "result_table_schema_id": "comsol_mphtable_probe_impedance_sweep_v1",
        "result_output_artifact_id": "tbl_probe_impedance_output_v1",
        "result_output_digest": "sha256:tbl_probe_impedance_output_v1",
        "result_output_path": "tbl_probe_impedance_output.json",
        "result_observable_id": "tbl_probe_impedance_boundary_absorption_v1",
        "result_observable_family": "boundary_absorption_table",
        "result_row_convention": "frequency_sweep_one_boundary_selection_per_row",
        "result_normalization_basis": "absorption_normalized_by_incident_power",
        "result_evaluation_method": "mphtable_after_dataset_bound_derived_value",
        "physics_convention_schema_id": "comsol_acoustic_impedance_derived_value_convention_v1",
        "result_postprocess_row_convention_schema_id": (
            "comsol_mphtable_boundary_absorption_postprocess_row_v1"
        ),
        "result_component_basis_schema_id": (
            "comsol_mphtable_probe_impedance_component_basis_v1"
        ),
        "result_artifact_id": "comsol_slot345_result_table_package_v1",
        "run_started_at": "2026-07-01T14:00:00+09:00",
        "comsol_version": "COMSOL 6.4.0.378",
        "timing_breakdown_s": {
            "attach_livelink_s": 0.05,
            "model_build_s": 0.20,
            "solve_s": 0.85,
            "table_eval_s": 0.04,
        },
        "execution": {
            "resultArtifactId": "comsol_slot345_result_table_package_v1",
            "solutionArtifactId": "comsol_slot353_solution_data_sol1_v1",
            "solutionDigest": "sha256:comsol_slot353_solution_data_sol1_v1",
            "sweepAxisId": "freq_grid_slot361_v1",
            "sweepAxisDigest": "sha256:freq_grid_slot361_v1",
            "sweepAxisRowCount": 4,
            "parameterSetArtifactId": "comsol_slot389_parameter_set_v1",
            "parameterSetDigest": "sha256:comsol-slot389-parameter-set-v1",
            "parameterSetPath": "validation/comsol/slot389_parameter_set.json",
            "objectiveObservableId": "comsol_slot389_impedance_absorption_objective_v1",
            "objectiveObservableFamily": "impedance_absorption_objective",
            "solverConfigurationArtifactId": "comsol_slot368_solver_config_sol1_v1",
            "solverConfigurationDigest": "sha256:comsol_slot368_solver_config_sol1_v1",
            "solverSequenceTag": "sol1",
            "linearSolver": "direct_pardiso",
            "relativeTolerance": 1.0e-6,
            "resultTableSchemaId": "comsol_mphtable_probe_impedance_sweep_v1",
            "physicsConventionSchemaId": "comsol_acoustic_impedance_derived_value_convention_v1",
            "resultPostprocessRowConventionSchemaId": (
                "comsol_mphtable_boundary_absorption_postprocess_row_v1"
            ),
            "resultComponentBasisSchemaId": (
                "comsol_mphtable_probe_impedance_component_basis_v1"
            ),
            "runStartedAt": "2026-07-01T14:00:00+09:00",
            "comsolVersion": "COMSOL 6.4.0.378",
            "timingBreakdown": {
                "attach_livelink_s": 0.05,
                "model_build_s": 0.20,
                "solve_s": 0.85,
                "table_eval_s": 0.04,
            },
        },
        "optimization": {
            "parameterSetArtifactId": "comsol_slot389_parameter_set_v1",
            "parameterSetDigest": "sha256:comsol-slot389-parameter-set-v1",
            "parameterSetPath": "validation/comsol/slot389_parameter_set.json",
            "objectiveObservableId": "comsol_slot389_impedance_absorption_objective_v1",
            "objectiveObservableFamily": "impedance_absorption_objective",
        },
        "columns": ["freq_Hz", "Z11_ohm", "absorption"],
        "units": {"freq_Hz": "Hz", "Z11_ohm": "ohm", "absorption": "1"},
        "independent_axis": "freq_Hz",
        "row_count": 4,
    }
    gate = solver_result_table_metadata_gate(
        metadata,
        required_columns=("freq_Hz", "Z11_ohm", "absorption"),
        required_units={"freq_Hz": "Hz", "Z11_ohm": "ohm"},
        independent_axis="freq_Hz",
        expected_source="COMSOL",
        expected_dataset_id="dset1",
        expected_solution_tag="sol1",
        expected_solution_artifact_id="comsol_slot353_solution_data_sol1_v1",
        expected_solution_digest="sha256:comsol_slot353_solution_data_sol1_v1",
        expected_sweep_axis_id="freq_grid_slot361_v1",
        expected_sweep_axis_digest="sha256:freq_grid_slot361_v1",
        expected_sweep_axis_row_count=4,
        expected_parameter_set_artifact_id="comsol_slot389_parameter_set_v1",
        expected_parameter_set_digest="sha256:comsol-slot389-parameter-set-v1",
        expected_parameter_set_path="validation/comsol/slot389_parameter_set.json",
        expected_objective_observable_id="comsol_slot389_impedance_absorption_objective_v1",
        expected_objective_observable_family="impedance_absorption_objective",
        expected_solver_configuration_artifact_id="comsol_slot368_solver_config_sol1_v1",
        expected_solver_configuration_digest="sha256:comsol_slot368_solver_config_sol1_v1",
        expected_solver_sequence_tag="sol1",
        expected_linear_solver="direct_pardiso",
        expected_relative_tolerance=1.0e-6,
        expected_study_tag="std1",
        expected_study_step_tag="stat",
        expected_table_id="tbl_probe_impedance",
        expected_selection_tags=("bnd_probe",),
        expected_entity_dimensions=("boundary",),
        expected_expressions=("intop_bnd(acpr.p_t)", "acpr.Q_tot"),
        expected_operator_tags=("intop_bnd",),
        expected_result_table_schema_id="comsol_mphtable_probe_impedance_sweep_v1",
        expected_result_output_artifact_id="tbl_probe_impedance_output_v1",
        expected_result_output_digest="sha256:tbl_probe_impedance_output_v1",
        expected_result_observable_id="tbl_probe_impedance_boundary_absorption_v1",
        expected_result_observable_family="boundary_absorption_table",
        expected_result_row_convention="frequency_sweep_one_boundary_selection_per_row",
        expected_result_normalization_basis="absorption_normalized_by_incident_power",
        expected_result_evaluation_method="mphtable_after_dataset_bound_derived_value",
        expected_physics_convention_schema_id="comsol_acoustic_impedance_derived_value_convention_v1",
        expected_result_postprocess_row_convention_schema_id=(
            "comsol_mphtable_boundary_absorption_postprocess_row_v1"
        ),
        expected_result_component_basis_schema_id=(
            "comsol_mphtable_probe_impedance_component_basis_v1"
        ),
        expected_result_artifact_id="comsol_slot345_result_table_package_v1",
        expected_comsol_version="COMSOL 6.4.0.378",
        require_solver_configuration=True,
        require_parameter_set_artifact=True,
        require_result_provenance=True,
        require_result_table_schema=True,
        require_result_output_artifact=True,
        require_physics_convention_schema=True,
        require_result_postprocess_row_convention_schema=True,
        require_result_component_basis_schema=True,
        min_rows=3,
    )

    assert gate["policy"] == "solver_result_table_metadata_gate"
    assert gate["status"] == "ok"
    assert gate["checks"]["required_columns_present"] is True
    assert gate["checks"]["expected_units_match"] is True
    assert gate["checks"]["independent_axis_is_column"] is True
    assert gate["checks"]["expected_dataset_id_matches"] is True
    assert gate["checks"]["expected_solution_tag_matches"] is True
    assert gate["checks"]["solution_artifact_id_consistent_when_present"] is True
    assert gate["checks"]["solution_digest_consistent_when_present"] is True
    assert gate["checks"]["solution_artifact_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_solution_artifact_id_matches"] is True
    assert gate["checks"]["solution_digest_recorded_when_expected"] is True
    assert gate["checks"]["expected_solution_digest_matches"] is True
    assert gate["checks"]["sweep_axis_id_consistent_when_present"] is True
    assert gate["checks"]["sweep_axis_digest_consistent_when_present"] is True
    assert gate["checks"]["sweep_axis_row_count_consistent_when_present"] is True
    assert gate["checks"]["sweep_axis_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_sweep_axis_id_matches"] is True
    assert gate["checks"]["sweep_axis_digest_recorded_when_expected"] is True
    assert gate["checks"]["expected_sweep_axis_digest_matches"] is True
    assert gate["checks"]["sweep_axis_row_count_recorded_when_expected"] is True
    assert gate["checks"]["expected_sweep_axis_row_count_matches"] is True
    assert gate["checks"]["sweep_axis_row_count_matches_table_rows_when_present"] is True
    assert gate["checks"]["parameter_set_artifact_id_consistent_when_present"] is True
    assert gate["checks"]["parameter_set_digest_consistent_when_present"] is True
    assert gate["checks"]["parameter_set_path_consistent_when_present"] is True
    assert gate["checks"]["objective_observable_id_consistent_when_present"] is True
    assert gate["checks"]["objective_observable_family_consistent_when_present"] is True
    assert gate["checks"]["parameter_set_artifact_id_recorded_when_required"] is True
    assert gate["checks"]["parameter_set_digest_recorded_when_required"] is True
    assert gate["checks"]["parameter_set_path_recorded_when_required"] is True
    assert gate["checks"]["parameter_set_artifact_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_parameter_set_artifact_id_matches"] is True
    assert gate["checks"]["parameter_set_digest_recorded_when_expected"] is True
    assert gate["checks"]["expected_parameter_set_digest_matches"] is True
    assert gate["checks"]["parameter_set_path_recorded_when_expected"] is True
    assert gate["checks"]["expected_parameter_set_path_matches"] is True
    assert gate["checks"]["objective_observable_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_objective_observable_id_matches"] is True
    assert gate["checks"]["objective_observable_family_recorded_when_expected"] is True
    assert gate["checks"]["expected_objective_observable_family_matches"] is True
    assert gate["checks"]["solver_configuration_artifact_id_consistent_when_present"] is True
    assert gate["checks"]["solver_configuration_digest_consistent_when_present"] is True
    assert gate["checks"]["solver_sequence_tag_consistent_when_present"] is True
    assert gate["checks"]["linear_solver_consistent_when_present"] is True
    assert gate["checks"]["relative_tolerance_consistent_when_present"] is True
    assert gate["checks"]["solver_configuration_artifact_id_recorded_when_required"] is True
    assert gate["checks"]["solver_configuration_digest_recorded_when_required"] is True
    assert gate["checks"]["solver_sequence_tag_recorded_when_required"] is True
    assert gate["checks"]["relative_tolerance_recorded_when_required"] is True
    assert gate["checks"]["expected_solver_configuration_artifact_id_matches"] is True
    assert gate["checks"]["expected_solver_configuration_digest_matches"] is True
    assert gate["checks"]["expected_solver_sequence_tag_matches"] is True
    assert gate["checks"]["expected_linear_solver_matches"] is True
    assert gate["checks"]["expected_relative_tolerance_matches"] is True
    assert gate["checks"]["relative_tolerance_finite_positive_when_present"] is True
    assert gate["checks"]["expected_study_tag_matches"] is True
    assert gate["checks"]["expected_study_step_tag_matches"] is True
    assert gate["checks"]["expected_table_id_matches"] is True
    assert gate["checks"]["expected_selection_tags_match"] is True
    assert gate["checks"]["expected_entity_dimensions_match"] is True
    assert gate["checks"]["selection_entity_scope_consistent"] is True
    assert gate["checks"]["expected_expressions_present"] is True
    assert gate["checks"]["expected_operator_tags_match"] is True
    assert gate["checks"]["result_table_schema_id_recorded_when_required"] is True
    assert gate["checks"]["result_table_schema_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_result_table_schema_id_matches"] is True
    assert gate["checks"]["physics_convention_schema_id_consistent_when_present"] is True
    assert gate["checks"]["physics_convention_schema_id_recorded_when_required"] is True
    assert gate["checks"]["physics_convention_schema_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_physics_convention_schema_id_matches"] is True
    assert gate["checks"]["result_postprocess_row_convention_schema_id_consistent_when_present"] is True
    assert gate["checks"]["result_postprocess_row_convention_schema_id_recorded_when_required"] is True
    assert gate["checks"]["result_postprocess_row_convention_schema_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_result_postprocess_row_convention_schema_id_matches"] is True
    assert gate["checks"]["result_component_basis_schema_id_consistent_when_present"] is True
    assert gate["checks"]["result_component_basis_schema_id_recorded_when_required"] is True
    assert gate["checks"]["result_component_basis_schema_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_result_component_basis_schema_id_matches"] is True
    assert gate["checks"]["result_output_artifact_id_recorded_when_required"] is True
    assert gate["checks"]["result_output_digest_recorded_when_required"] is True
    assert gate["checks"]["result_output_path_recorded_when_required"] is True
    assert gate["checks"]["expected_result_output_artifact_id_matches"] is True
    assert gate["checks"]["expected_result_output_digest_matches"] is True
    assert gate["checks"]["result_observable_id_consistent_when_present"] is True
    assert gate["checks"]["result_observable_family_consistent_when_present"] is True
    assert gate["checks"]["result_observable_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_result_observable_id_matches"] is True
    assert gate["checks"]["result_observable_family_recorded_when_expected"] is True
    assert gate["checks"]["expected_result_observable_family_matches"] is True
    assert gate["checks"]["result_row_convention_recorded_when_expected"] is True
    assert gate["checks"]["expected_result_row_convention_matches"] is True
    assert gate["checks"]["result_normalization_basis_recorded_when_expected"] is True
    assert gate["checks"]["expected_result_normalization_basis_matches"] is True
    assert gate["checks"]["result_evaluation_method_recorded_when_expected"] is True
    assert gate["checks"]["expected_result_evaluation_method_matches"] is True
    assert gate["checks"]["result_artifact_id_recorded_when_required"] is True
    assert gate["checks"]["expected_result_artifact_id_matches"] is True
    assert gate["checks"]["run_started_at_recorded_when_required"] is True
    assert gate["checks"]["run_started_at_parseable_when_present"] is True
    assert gate["checks"]["comsol_version_recorded_when_required"] is True
    assert gate["checks"]["expected_comsol_version_matches"] is True
    assert gate["checks"]["timing_breakdown_recorded_when_required"] is True
    assert gate["checks"]["timing_breakdown_has_at_least_four_items"] is True
    assert gate["checks"]["timing_breakdown_values_finite_nonnegative"] is True
    assert gate["dataset_id"] == "dset1"
    assert gate["solution_tag"] == "sol1"
    assert gate["solution_artifact_id"] == "comsol_slot353_solution_data_sol1_v1"
    assert gate["solution_digest"] == "sha256:comsol_slot353_solution_data_sol1_v1"
    assert gate["sweep_axis_id"] == "freq_grid_slot361_v1"
    assert gate["sweep_axis_digest"] == "sha256:freq_grid_slot361_v1"
    assert gate["sweep_axis_row_count"] == 4
    assert gate["parameter_set_artifact_id"] == "comsol_slot389_parameter_set_v1"
    assert gate["parameter_set_digest"] == "sha256:comsol-slot389-parameter-set-v1"
    assert gate["parameter_set_path"] == "validation/comsol/slot389_parameter_set.json"
    assert gate["objective_observable_id"] == "comsol_slot389_impedance_absorption_objective_v1"
    assert gate["objective_observable_family"] == "impedance_absorption_objective"
    assert gate["solver_configuration_artifact_id"] == "comsol_slot368_solver_config_sol1_v1"
    assert gate["solver_configuration_digest"] == "sha256:comsol_slot368_solver_config_sol1_v1"
    assert gate["solver_sequence_tag"] == "sol1"
    assert gate["linear_solver"] == "direct_pardiso"
    assert gate["relative_tolerance"] == pytest.approx(1.0e-6)
    assert gate["study_tag"] == "std1"
    assert gate["study_step_tag"] == "stat"
    assert gate["table_id"] == "tbl_probe_impedance"
    assert gate["selection_tags"] == ["bnd_probe"]
    assert gate["entity_dimensions"] == ["boundary"]
    assert gate["expressions"] == ["freq", "intop_bnd(acpr.p_t)", "acpr.Q_tot"]
    assert gate["operator_tags"] == ["intop_bnd"]
    assert gate["result_table_schema_id"] == "comsol_mphtable_probe_impedance_sweep_v1"
    assert gate["physics_convention_schema_id"] == "comsol_acoustic_impedance_derived_value_convention_v1"
    assert gate["require_physics_convention_schema"] is True
    assert gate["result_postprocess_row_convention_schema_id"] == (
        "comsol_mphtable_boundary_absorption_postprocess_row_v1"
    )
    assert gate["result_postprocess_row_convention_schema_ids"] == [
        "comsol_mphtable_boundary_absorption_postprocess_row_v1"
    ]
    assert gate["require_result_postprocess_row_convention_schema"] is True
    assert gate["result_component_basis_schema_id"] == (
        "comsol_mphtable_probe_impedance_component_basis_v1"
    )
    assert gate["result_component_basis_schema_ids"] == [
        "comsol_mphtable_probe_impedance_component_basis_v1"
    ]
    assert gate["require_result_component_basis_schema"] is True
    assert gate["result_output_artifact_id"] == "tbl_probe_impedance_output_v1"
    assert gate["result_output_digest"] == "sha256:tbl_probe_impedance_output_v1"
    assert gate["result_output_path"] == "tbl_probe_impedance_output.json"
    assert gate["result_observable_id"] == "tbl_probe_impedance_boundary_absorption_v1"
    assert gate["result_observable_family"] == "boundary_absorption_table"
    assert gate["result_row_convention"] == "frequency_sweep_one_boundary_selection_per_row"
    assert gate["result_normalization_basis"] == "absorption_normalized_by_incident_power"
    assert gate["result_evaluation_method"] == "mphtable_after_dataset_bound_derived_value"
    assert gate["result_artifact_id"] == "comsol_slot345_result_table_package_v1"
    assert gate["run_started_at"] == "2026-07-01T14:00:00+09:00"
    assert gate["comsol_version"] == "COMSOL 6.4.0.378"
    assert gate["timing_breakdown_names"] == [
        "attach_livelink_s",
        "model_build_s",
        "solve_s",
        "table_eval_s",
    ]
    assert gate["timing_breakdown_seconds"]["solve_s"] == pytest.approx(0.85)

    stale_table_schema = dict(metadata)
    stale_table_schema["result_table_schema_id"] = "comsol_scalar_probe_impedance_v0"
    stale_table_schema["execution"] = {
        **metadata["execution"],
        "resultTableSchemaId": "comsol_scalar_probe_impedance_v0",
    }
    stale_table_schema_gate = solver_result_table_metadata_gate(
        stale_table_schema,
        required_columns=("freq_Hz", "Z11_ohm", "absorption"),
        expected_result_table_schema_id="comsol_mphtable_probe_impedance_sweep_v1",
        expected_result_output_artifact_id="tbl_probe_impedance_output_v1",
        expected_result_output_digest="sha256:tbl_probe_impedance_output_v1",
        require_result_table_schema=True,
    )
    assert stale_table_schema_gate["status"] == "needs_attention"
    assert (
        stale_table_schema_gate["checks"]["expected_result_table_schema_id_matches"]
        is False
    )
    assert (
        stale_table_schema_gate["checks"]["expected_result_output_artifact_id_matches"]
        is True
    )
    assert stale_table_schema_gate["checks"]["expected_result_output_digest_matches"] is True

    missing_table_schema = {
        key: value for key, value in metadata.items() if key != "result_table_schema_id"
    }
    missing_table_schema["execution"] = dict(metadata["execution"])
    missing_table_schema["execution"].pop("resultTableSchemaId")
    missing_table_schema_gate = solver_result_table_metadata_gate(
        missing_table_schema,
        required_columns=("freq_Hz", "Z11_ohm", "absorption"),
        require_result_table_schema=True,
    )
    assert missing_table_schema_gate["status"] == "needs_attention"
    assert (
        missing_table_schema_gate["checks"]["result_table_schema_id_recorded_when_required"]
        is False
    )

    stale_physics_convention_schema = dict(metadata)
    stale_physics_convention_schema["physics_convention_schema_id"] = (
        "comsol_value_only_convention_v0"
    )
    stale_physics_convention_schema["execution"] = {
        **metadata["execution"],
        "physicsConventionSchemaId": "comsol_value_only_convention_v0",
    }
    stale_physics_convention_schema_gate = solver_result_table_metadata_gate(
        stale_physics_convention_schema,
        required_columns=("freq_Hz", "Z11_ohm", "absorption"),
        expected_result_table_schema_id="comsol_mphtable_probe_impedance_sweep_v1",
        expected_physics_convention_schema_id="comsol_acoustic_impedance_derived_value_convention_v1",
        require_result_table_schema=True,
        require_physics_convention_schema=True,
    )
    assert stale_physics_convention_schema_gate["status"] == "needs_attention"
    assert stale_physics_convention_schema_gate["checks"]["expected_result_table_schema_id_matches"] is True
    assert stale_physics_convention_schema_gate["checks"]["physics_convention_schema_id_consistent_when_present"] is True
    assert stale_physics_convention_schema_gate["checks"]["physics_convention_schema_id_recorded_when_required"] is True
    assert stale_physics_convention_schema_gate["checks"]["expected_physics_convention_schema_id_matches"] is False

    missing_physics_convention_schema = {
        key: value for key, value in metadata.items() if key != "physics_convention_schema_id"
    }
    missing_physics_convention_schema["execution"] = dict(metadata["execution"])
    missing_physics_convention_schema["execution"].pop("physicsConventionSchemaId")
    missing_physics_convention_schema_gate = solver_result_table_metadata_gate(
        missing_physics_convention_schema,
        required_columns=("freq_Hz", "Z11_ohm", "absorption"),
        expected_physics_convention_schema_id="comsol_acoustic_impedance_derived_value_convention_v1",
        require_physics_convention_schema=True,
    )
    assert missing_physics_convention_schema_gate["status"] == "needs_attention"
    assert (
        missing_physics_convention_schema_gate["checks"][
            "physics_convention_schema_id_recorded_when_required"
        ]
        is False
    )
    assert (
        missing_physics_convention_schema_gate["checks"][
            "physics_convention_schema_id_recorded_when_expected"
        ]
        is False
    )

    stale_postprocess_row_convention_schema = dict(metadata)
    stale_postprocess_row_convention_schema["result_postprocess_row_convention_schema_id"] = (
        "comsol_scalar_residual_row_v0"
    )
    stale_postprocess_row_convention_schema["execution"] = {
        **metadata["execution"],
        "resultPostprocessRowConventionSchemaId": "comsol_scalar_residual_row_v0",
    }
    stale_postprocess_row_convention_schema_gate = solver_result_table_metadata_gate(
        stale_postprocess_row_convention_schema,
        required_columns=("freq_Hz", "Z11_ohm", "absorption"),
        expected_result_table_schema_id="comsol_mphtable_probe_impedance_sweep_v1",
        expected_physics_convention_schema_id="comsol_acoustic_impedance_derived_value_convention_v1",
        expected_result_postprocess_row_convention_schema_id=(
            "comsol_mphtable_boundary_absorption_postprocess_row_v1"
        ),
        require_result_table_schema=True,
        require_physics_convention_schema=True,
        require_result_postprocess_row_convention_schema=True,
    )
    assert stale_postprocess_row_convention_schema_gate["status"] == "needs_attention"
    assert stale_postprocess_row_convention_schema_gate["checks"]["expected_result_table_schema_id_matches"] is True
    assert stale_postprocess_row_convention_schema_gate["checks"]["expected_physics_convention_schema_id_matches"] is True
    assert stale_postprocess_row_convention_schema_gate["checks"]["result_postprocess_row_convention_schema_id_recorded_when_required"] is True
    assert (
        stale_postprocess_row_convention_schema_gate["checks"][
            "expected_result_postprocess_row_convention_schema_id_matches"
        ]
        is False
    )

    missing_postprocess_row_convention_schema = {
        key: value
        for key, value in metadata.items()
        if key != "result_postprocess_row_convention_schema_id"
    }
    missing_postprocess_row_convention_schema["execution"] = dict(metadata["execution"])
    missing_postprocess_row_convention_schema["execution"].pop(
        "resultPostprocessRowConventionSchemaId"
    )
    missing_postprocess_row_convention_schema_gate = solver_result_table_metadata_gate(
        missing_postprocess_row_convention_schema,
        required_columns=("freq_Hz", "Z11_ohm", "absorption"),
        expected_result_postprocess_row_convention_schema_id=(
            "comsol_mphtable_boundary_absorption_postprocess_row_v1"
        ),
        require_result_postprocess_row_convention_schema=True,
    )
    assert missing_postprocess_row_convention_schema_gate["status"] == "needs_attention"
    assert (
        missing_postprocess_row_convention_schema_gate["checks"][
            "result_postprocess_row_convention_schema_id_recorded_when_required"
        ]
        is False
    )
    assert (
        missing_postprocess_row_convention_schema_gate["checks"][
            "result_postprocess_row_convention_schema_id_recorded_when_expected"
        ]
        is False
    )

    stale_component_basis_schema = dict(metadata)
    stale_component_basis_schema["result_component_basis_schema_id"] = (
        "comsol_mphtable_probe_impedance_abs_phase_v0"
    )
    stale_component_basis_schema["execution"] = {
        **metadata["execution"],
        "resultComponentBasisSchemaId": "comsol_mphtable_probe_impedance_abs_phase_v0",
    }
    stale_component_basis_schema_gate = solver_result_table_metadata_gate(
        stale_component_basis_schema,
        required_columns=("freq_Hz", "Z11_ohm", "absorption"),
        expected_result_table_schema_id="comsol_mphtable_probe_impedance_sweep_v1",
        expected_physics_convention_schema_id="comsol_acoustic_impedance_derived_value_convention_v1",
        expected_result_postprocess_row_convention_schema_id=(
            "comsol_mphtable_boundary_absorption_postprocess_row_v1"
        ),
        expected_result_component_basis_schema_id=(
            "comsol_mphtable_probe_impedance_component_basis_v1"
        ),
        require_result_component_basis_schema=True,
    )
    assert stale_component_basis_schema_gate["status"] == "needs_attention"
    assert stale_component_basis_schema_gate["checks"]["expected_result_table_schema_id_matches"] is True
    assert stale_component_basis_schema_gate["checks"]["expected_physics_convention_schema_id_matches"] is True
    assert (
        stale_component_basis_schema_gate["checks"][
            "expected_result_postprocess_row_convention_schema_id_matches"
        ]
        is True
    )
    assert stale_component_basis_schema_gate["checks"]["result_component_basis_schema_id_recorded_when_required"] is True
    assert (
        stale_component_basis_schema_gate["checks"][
            "expected_result_component_basis_schema_id_matches"
        ]
        is False
    )

    missing_component_basis_schema = {
        key: value
        for key, value in metadata.items()
        if key != "result_component_basis_schema_id"
    }
    missing_component_basis_schema["execution"] = dict(metadata["execution"])
    missing_component_basis_schema["execution"].pop("resultComponentBasisSchemaId")
    missing_component_basis_schema_gate = solver_result_table_metadata_gate(
        missing_component_basis_schema,
        required_columns=("freq_Hz", "Z11_ohm", "absorption"),
        expected_result_component_basis_schema_id=(
            "comsol_mphtable_probe_impedance_component_basis_v1"
        ),
        require_result_component_basis_schema=True,
    )
    assert missing_component_basis_schema_gate["status"] == "needs_attention"
    assert (
        missing_component_basis_schema_gate["checks"][
            "result_component_basis_schema_id_recorded_when_required"
        ]
        is False
    )
    assert (
        missing_component_basis_schema_gate["checks"][
            "result_component_basis_schema_id_recorded_when_expected"
        ]
        is False
    )

    missing_unit = solver_result_table_metadata_gate(
        {**metadata, "units": {"freq_Hz": "Hz"}},
        required_columns=("freq_Hz", "Z11_ohm"),
        required_units={"freq_Hz": "Hz", "Z11_ohm": "ohm"},
        independent_axis="freq_Hz",
    )
    assert missing_unit["status"] == "needs_attention"
    assert missing_unit["checks"]["units_recorded_for_required_columns"] is False
    assert missing_unit["checks"]["expected_units_match"] is False

    duplicate = solver_result_table_metadata_gate(
        {**metadata, "columns": ["freq_Hz", "freq_Hz", "Z11_ohm"], "row_count": 1},
        required_columns=("freq_Hz", "Z11_ohm"),
        min_rows=3,
    )
    assert duplicate["status"] == "needs_attention"
    assert duplicate["checks"]["columns_unique"] is False
    assert duplicate["checks"]["row_count_at_least_minimum"] is False

    missing_dataset = solver_result_table_metadata_gate(
        {key: value for key, value in metadata.items() if key != "dataset"},
        required_columns=("freq_Hz", "Z11_ohm"),
        expected_dataset_id="dset1",
    )
    assert missing_dataset["status"] == "needs_attention"
    assert missing_dataset["checks"]["dataset_id_recorded"] is False

    stale_solution = solver_result_table_metadata_gate(
        {**metadata, "solution_tag": "sol_old"},
        required_columns=("freq_Hz", "Z11_ohm"),
        expected_dataset_id="dset1",
        expected_solution_tag="sol1",
    )
    assert stale_solution["status"] == "needs_attention"
    assert stale_solution["checks"]["expected_dataset_id_matches"] is True
    assert stale_solution["checks"]["expected_solution_tag_matches"] is False

    stale_solution_artifact = solver_result_table_metadata_gate(
        {**metadata, "solution_artifact_id": "comsol_slot345_old_solution_data"},
        required_columns=("freq_Hz", "Z11_ohm"),
        expected_dataset_id="dset1",
        expected_solution_tag="sol1",
        expected_solution_artifact_id="comsol_slot353_solution_data_sol1_v1",
        expected_solution_digest="sha256:comsol_slot353_solution_data_sol1_v1",
    )
    assert stale_solution_artifact["status"] == "needs_attention"
    assert stale_solution_artifact["checks"]["expected_solution_tag_matches"] is True
    assert stale_solution_artifact["checks"]["solution_artifact_id_consistent_when_present"] is False
    assert stale_solution_artifact["checks"]["expected_solution_artifact_id_matches"] is False
    assert stale_solution_artifact["checks"]["expected_solution_digest_matches"] is True

    missing_solution_digest = {
        key: value for key, value in metadata.items() if key != "solution_digest"
    }
    missing_solution_digest["execution"] = dict(metadata["execution"])
    missing_solution_digest["execution"].pop("solutionDigest")
    missing_solution_digest_gate = solver_result_table_metadata_gate(
        missing_solution_digest,
        required_columns=("freq_Hz", "Z11_ohm"),
        expected_solution_digest="sha256:comsol_slot353_solution_data_sol1_v1",
    )
    assert missing_solution_digest_gate["status"] == "needs_attention"
    assert missing_solution_digest_gate["checks"]["solution_digest_recorded_when_expected"] is False

    stale_sweep_axis_digest = dict(metadata)
    stale_sweep_axis_digest["sweep_axis_digest"] = "sha256:freq_grid_old"
    stale_sweep_axis_digest["execution"] = {
        **metadata["execution"],
        "sweepAxisDigest": "sha256:freq_grid_old",
    }
    stale_sweep_axis_digest_gate = solver_result_table_metadata_gate(
        stale_sweep_axis_digest,
        required_columns=("freq_Hz", "Z11_ohm"),
        expected_sweep_axis_id="freq_grid_slot361_v1",
        expected_sweep_axis_digest="sha256:freq_grid_slot361_v1",
        expected_sweep_axis_row_count=4,
    )
    assert stale_sweep_axis_digest_gate["status"] == "needs_attention"
    assert stale_sweep_axis_digest_gate["checks"]["expected_sweep_axis_id_matches"] is True
    assert stale_sweep_axis_digest_gate["checks"]["expected_sweep_axis_digest_matches"] is False
    assert stale_sweep_axis_digest_gate["checks"]["expected_sweep_axis_row_count_matches"] is True

    missing_sweep_axis_digest = {
        key: value for key, value in metadata.items() if key != "sweep_axis_digest"
    }
    missing_sweep_axis_digest["execution"] = dict(metadata["execution"])
    missing_sweep_axis_digest["execution"].pop("sweepAxisDigest")
    missing_sweep_axis_digest_gate = solver_result_table_metadata_gate(
        missing_sweep_axis_digest,
        required_columns=("freq_Hz", "Z11_ohm"),
        expected_sweep_axis_digest="sha256:freq_grid_slot361_v1",
    )
    assert missing_sweep_axis_digest_gate["status"] == "needs_attention"
    assert missing_sweep_axis_digest_gate["checks"]["sweep_axis_digest_recorded_when_expected"] is False

    wrong_sweep_axis_row_count = dict(metadata)
    wrong_sweep_axis_row_count["sweep_axis_row_count"] = 3
    wrong_sweep_axis_row_count["execution"] = {
        **metadata["execution"],
        "sweepAxisRowCount": 3,
    }
    wrong_sweep_axis_row_count_gate = solver_result_table_metadata_gate(
        wrong_sweep_axis_row_count,
        required_columns=("freq_Hz", "Z11_ohm"),
        expected_sweep_axis_row_count=4,
    )
    assert wrong_sweep_axis_row_count_gate["status"] == "needs_attention"
    assert wrong_sweep_axis_row_count_gate["checks"]["expected_sweep_axis_row_count_matches"] is False
    assert wrong_sweep_axis_row_count_gate["checks"]["sweep_axis_row_count_matches_table_rows_when_present"] is False

    stale_parameter_set_digest = dict(metadata)
    stale_parameter_set_digest["parameter_set_digest"] = "sha256:old-parameter-set"
    stale_parameter_set_digest["execution"] = {
        **metadata["execution"],
        "parameterSetDigest": "sha256:old-parameter-set",
    }
    stale_parameter_set_digest["optimization"] = {
        **metadata["optimization"],
        "parameterSetDigest": "sha256:old-parameter-set",
    }
    stale_parameter_set_digest_gate = solver_result_table_metadata_gate(
        stale_parameter_set_digest,
        required_columns=("freq_Hz", "Z11_ohm"),
        expected_parameter_set_digest="sha256:comsol-slot389-parameter-set-v1",
        require_parameter_set_artifact=True,
    )
    assert stale_parameter_set_digest_gate["status"] == "needs_attention"
    assert stale_parameter_set_digest_gate["checks"]["parameter_set_digest_consistent_when_present"] is True
    assert stale_parameter_set_digest_gate["checks"]["expected_parameter_set_digest_matches"] is False

    missing_parameter_set_path = {
        key: value for key, value in metadata.items() if key != "parameter_set_path"
    }
    missing_parameter_set_path["execution"] = dict(metadata["execution"])
    missing_parameter_set_path["execution"].pop("parameterSetPath")
    missing_parameter_set_path["optimization"] = dict(metadata["optimization"])
    missing_parameter_set_path["optimization"].pop("parameterSetPath")
    missing_parameter_set_path_gate = solver_result_table_metadata_gate(
        missing_parameter_set_path,
        required_columns=("freq_Hz", "Z11_ohm"),
        require_parameter_set_artifact=True,
    )
    assert missing_parameter_set_path_gate["status"] == "needs_attention"
    assert missing_parameter_set_path_gate["checks"]["parameter_set_artifact_id_recorded_when_required"] is True
    assert missing_parameter_set_path_gate["checks"]["parameter_set_digest_recorded_when_required"] is True
    assert missing_parameter_set_path_gate["checks"]["parameter_set_path_recorded_when_required"] is False

    wrong_objective_family = dict(metadata)
    wrong_objective_family["objective_observable_family"] = "remote_field_map"
    wrong_objective_family["execution"] = {
        **metadata["execution"],
        "objectiveObservableFamily": "remote_field_map",
    }
    wrong_objective_family["optimization"] = {
        **metadata["optimization"],
        "objectiveObservableFamily": "remote_field_map",
    }
    wrong_objective_family_gate = solver_result_table_metadata_gate(
        wrong_objective_family,
        required_columns=("freq_Hz", "Z11_ohm"),
        expected_objective_observable_family="impedance_absorption_objective",
    )
    assert wrong_objective_family_gate["status"] == "needs_attention"
    assert wrong_objective_family_gate["checks"]["objective_observable_family_consistent_when_present"] is True
    assert wrong_objective_family_gate["checks"]["expected_objective_observable_family_matches"] is False

    stale_solver_configuration_digest = dict(metadata)
    stale_solver_configuration_digest["solver_configuration_digest"] = "sha256:solver_config_old"
    stale_solver_configuration_digest["execution"] = {
        **metadata["execution"],
        "solverConfigurationDigest": "sha256:solver_config_old",
    }
    stale_solver_configuration_digest_gate = solver_result_table_metadata_gate(
        stale_solver_configuration_digest,
        required_columns=("freq_Hz", "Z11_ohm"),
        expected_solver_configuration_artifact_id="comsol_slot368_solver_config_sol1_v1",
        expected_solver_configuration_digest="sha256:comsol_slot368_solver_config_sol1_v1",
        expected_solver_sequence_tag="sol1",
        expected_relative_tolerance=1.0e-6,
        require_solver_configuration=True,
    )
    assert stale_solver_configuration_digest_gate["status"] == "needs_attention"
    assert stale_solver_configuration_digest_gate["checks"]["expected_solver_configuration_artifact_id_matches"] is True
    assert stale_solver_configuration_digest_gate["checks"]["expected_solver_configuration_digest_matches"] is False
    assert stale_solver_configuration_digest_gate["checks"]["expected_solver_sequence_tag_matches"] is True
    assert stale_solver_configuration_digest_gate["checks"]["expected_relative_tolerance_matches"] is True

    missing_solver_configuration_digest = {
        key: value for key, value in metadata.items() if key != "solver_configuration_digest"
    }
    missing_solver_configuration_digest["execution"] = dict(metadata["execution"])
    missing_solver_configuration_digest["execution"].pop("solverConfigurationDigest")
    missing_solver_configuration_digest_gate = solver_result_table_metadata_gate(
        missing_solver_configuration_digest,
        required_columns=("freq_Hz", "Z11_ohm"),
        expected_solver_configuration_digest="sha256:comsol_slot368_solver_config_sol1_v1",
        require_solver_configuration=True,
    )
    assert missing_solver_configuration_digest_gate["status"] == "needs_attention"
    assert missing_solver_configuration_digest_gate["checks"]["solver_configuration_digest_recorded_when_required"] is False
    assert missing_solver_configuration_digest_gate["checks"]["expected_solver_configuration_digest_matches"] is False

    wrong_relative_tolerance = {
        **metadata,
        "relative_tolerance": 1.0e-3,
        "execution": {
            **metadata["execution"],
            "relativeTolerance": 1.0e-3,
        },
    }
    wrong_relative_tolerance_gate = solver_result_table_metadata_gate(
        wrong_relative_tolerance,
        required_columns=("freq_Hz", "Z11_ohm"),
        expected_solver_configuration_digest="sha256:comsol_slot368_solver_config_sol1_v1",
        expected_relative_tolerance=1.0e-6,
    )
    assert wrong_relative_tolerance_gate["status"] == "needs_attention"
    assert wrong_relative_tolerance_gate["checks"]["expected_solver_configuration_digest_matches"] is True
    assert wrong_relative_tolerance_gate["checks"]["expected_relative_tolerance_matches"] is False

    stale_solver_sequence = {
        **metadata,
        "solver_sequence_tag": "sol_old",
        "execution": {
            **metadata["execution"],
            "solverSequenceTag": "sol_old",
        },
    }
    stale_solver_sequence_gate = solver_result_table_metadata_gate(
        stale_solver_sequence,
        required_columns=("freq_Hz", "Z11_ohm"),
        expected_solver_sequence_tag="sol1",
        expected_linear_solver="direct_pardiso",
    )
    assert stale_solver_sequence_gate["status"] == "needs_attention"
    assert stale_solver_sequence_gate["checks"]["expected_solver_sequence_tag_matches"] is False
    assert stale_solver_sequence_gate["checks"]["expected_linear_solver_matches"] is True

    stale_study = solver_result_table_metadata_gate(
        {**metadata, "study_tag": "std_old"},
        required_columns=("freq_Hz", "Z11_ohm"),
        expected_dataset_id="dset1",
        expected_solution_tag="sol1",
        expected_study_tag="std1",
        expected_study_step_tag="stat",
    )
    assert stale_study["status"] == "needs_attention"
    assert stale_study["checks"]["expected_solution_tag_matches"] is True
    assert stale_study["checks"]["expected_study_tag_matches"] is False
    assert stale_study["checks"]["expected_study_step_tag_matches"] is True

    missing_study_step = solver_result_table_metadata_gate(
        {key: value for key, value in metadata.items() if key != "study_step_tag"},
        required_columns=("freq_Hz", "Z11_ohm"),
        expected_study_tag="std1",
        expected_study_step_tag="stat",
    )
    assert missing_study_step["status"] == "needs_attention"
    assert missing_study_step["checks"]["expected_study_tag_matches"] is True
    assert missing_study_step["checks"]["study_step_tag_recorded"] is False

    stale_selection = solver_result_table_metadata_gate(
        {**metadata, "selection_tags": ["bnd_old"]},
        required_columns=("freq_Hz", "Z11_ohm"),
        expected_selection_tags=("bnd_probe",),
        expected_entity_dimensions=("boundary",),
    )
    assert stale_selection["status"] == "needs_attention"
    assert stale_selection["checks"]["expected_selection_tags_match"] is False
    assert stale_selection["checks"]["expected_entity_dimensions_match"] is True

    missing_entity_dimension = solver_result_table_metadata_gate(
        {key: value for key, value in metadata.items() if key != "entity_dimensions"},
        required_columns=("freq_Hz", "Z11_ohm"),
        expected_selection_tags=("bnd_probe",),
        expected_entity_dimensions=("boundary",),
    )
    assert missing_entity_dimension["status"] == "needs_attention"
    assert missing_entity_dimension["checks"]["expected_selection_tags_match"] is True
    assert missing_entity_dimension["checks"]["entity_dimensions_recorded"] is False

    inconsistent_scope = solver_result_table_metadata_gate(
        {**metadata, "selection_tags": ["bnd_probe", "bnd_far"], "entity_dimensions": ["boundary", "domain", "edge"]},
        required_columns=("freq_Hz", "Z11_ohm"),
    )
    assert inconsistent_scope["status"] == "needs_attention"
    assert inconsistent_scope["checks"]["selection_entity_scope_consistent"] is False

    missing_expressions = solver_result_table_metadata_gate(
        {key: value for key, value in metadata.items() if key != "expressions"},
        required_columns=("freq_Hz", "Z11_ohm"),
        expected_expressions=("intop_bnd(acpr.p_t)",),
    )
    assert missing_expressions["status"] == "needs_attention"
    assert missing_expressions["checks"]["expressions_recorded_when_expected"] is False

    stale_expression = solver_result_table_metadata_gate(
        {**metadata, "expressions": ["freq", "intop_old(acpr.p_t)", "acpr.Q_tot"]},
        required_columns=("freq_Hz", "Z11_ohm"),
        expected_expressions=("intop_bnd(acpr.p_t)", "acpr.Q_tot"),
    )
    assert stale_expression["status"] == "needs_attention"
    assert stale_expression["checks"]["expected_expressions_present"] is False

    stale_operator = solver_result_table_metadata_gate(
        {**metadata, "operator_tags": ["intop_old"]},
        required_columns=("freq_Hz", "Z11_ohm"),
        expected_operator_tags=("intop_bnd",),
    )
    assert stale_operator["status"] == "needs_attention"
    assert stale_operator["checks"]["expected_operator_tags_match"] is False

    stale_output_artifact = solver_result_table_metadata_gate(
        {**metadata, "result_output_artifact_id": "stale_tbl_output"},
        required_columns=("freq_Hz", "Z11_ohm"),
        expected_result_output_artifact_id="tbl_probe_impedance_output_v1",
    )
    assert stale_output_artifact["status"] == "needs_attention"
    assert stale_output_artifact["checks"]["expected_result_output_artifact_id_matches"] is False

    stale_output_digest = solver_result_table_metadata_gate(
        {**metadata, "result_output_digest": "sha256:stale_tbl_output"},
        required_columns=("freq_Hz", "Z11_ohm"),
        expected_result_output_digest="sha256:tbl_probe_impedance_output_v1",
    )
    assert stale_output_digest["status"] == "needs_attention"
    assert stale_output_digest["checks"]["expected_result_output_digest_matches"] is False

    missing_output_path_metadata = {
        key: value for key, value in metadata.items() if key != "result_output_path"
    }
    missing_output_path = solver_result_table_metadata_gate(
        missing_output_path_metadata,
        required_columns=("freq_Hz", "Z11_ohm"),
        require_result_output_artifact=True,
    )
    assert missing_output_path["status"] == "needs_attention"
    assert missing_output_path["checks"]["result_output_path_recorded_when_required"] is False

    stale_observable_id = solver_result_table_metadata_gate(
        {**metadata, "result_observable_id": "stale_remote_field_map_observable"},
        required_columns=("freq_Hz", "Z11_ohm"),
        expected_result_observable_id="tbl_probe_impedance_boundary_absorption_v1",
        expected_result_observable_family="boundary_absorption_table",
        expected_result_output_artifact_id="tbl_probe_impedance_output_v1",
    )
    assert stale_observable_id["status"] == "needs_attention"
    assert stale_observable_id["checks"]["expected_result_observable_id_matches"] is False
    assert stale_observable_id["checks"]["expected_result_observable_family_matches"] is True
    assert stale_observable_id["checks"]["expected_result_output_artifact_id_matches"] is True

    wrong_observable_family = solver_result_table_metadata_gate(
        {**metadata, "result_observable_family": "remote_field_map"},
        required_columns=("freq_Hz", "Z11_ohm"),
        expected_result_observable_id="tbl_probe_impedance_boundary_absorption_v1",
        expected_result_observable_family="boundary_absorption_table",
        expected_result_output_artifact_id="tbl_probe_impedance_output_v1",
    )
    assert wrong_observable_family["status"] == "needs_attention"
    assert wrong_observable_family["checks"]["expected_result_observable_id_matches"] is True
    assert wrong_observable_family["checks"]["expected_result_observable_family_matches"] is False
    assert wrong_observable_family["checks"]["expected_result_output_artifact_id_matches"] is True

    stale_row_convention = solver_result_table_metadata_gate(
        {**metadata, "result_row_convention": "domain_average_rows"},
        required_columns=("freq_Hz", "Z11_ohm"),
        expected_result_observable_id="tbl_probe_impedance_boundary_absorption_v1",
        expected_result_observable_family="boundary_absorption_table",
        expected_result_row_convention="frequency_sweep_one_boundary_selection_per_row",
        expected_result_normalization_basis="absorption_normalized_by_incident_power",
    )
    assert stale_row_convention["status"] == "needs_attention"
    assert stale_row_convention["checks"]["expected_result_observable_id_matches"] is True
    assert stale_row_convention["checks"]["expected_result_observable_family_matches"] is True
    assert stale_row_convention["checks"]["expected_result_row_convention_matches"] is False
    assert stale_row_convention["checks"]["expected_result_normalization_basis_matches"] is True

    wrong_normalization_basis = solver_result_table_metadata_gate(
        {**metadata, "result_normalization_basis": "raw_boundary_integral_no_power_normalization"},
        required_columns=("freq_Hz", "Z11_ohm"),
        expected_result_row_convention="frequency_sweep_one_boundary_selection_per_row",
        expected_result_normalization_basis="absorption_normalized_by_incident_power",
    )
    assert wrong_normalization_basis["status"] == "needs_attention"
    assert wrong_normalization_basis["checks"]["expected_result_row_convention_matches"] is True
    assert wrong_normalization_basis["checks"]["expected_result_normalization_basis_matches"] is False

    wrong_evaluation_method = solver_result_table_metadata_gate(
        {**metadata, "result_evaluation_method": "mphglobal_current_dataset_default"},
        required_columns=("freq_Hz", "Z11_ohm"),
        expected_result_row_convention="frequency_sweep_one_boundary_selection_per_row",
        expected_result_normalization_basis="absorption_normalized_by_incident_power",
        expected_result_evaluation_method="mphtable_after_dataset_bound_derived_value",
    )
    assert wrong_evaluation_method["status"] == "needs_attention"
    assert wrong_evaluation_method["checks"]["expected_result_row_convention_matches"] is True
    assert wrong_evaluation_method["checks"]["expected_result_normalization_basis_matches"] is True
    assert wrong_evaluation_method["checks"]["expected_result_evaluation_method_matches"] is False

    missing_evaluation_method_metadata = {
        key: value for key, value in metadata.items() if key != "result_evaluation_method"
    }
    missing_evaluation_method = solver_result_table_metadata_gate(
        missing_evaluation_method_metadata,
        required_columns=("freq_Hz", "Z11_ohm"),
        expected_result_evaluation_method="mphtable_after_dataset_bound_derived_value",
    )
    assert missing_evaluation_method["status"] == "needs_attention"
    assert missing_evaluation_method["checks"]["result_evaluation_method_recorded_when_expected"] is False
    assert missing_evaluation_method["checks"]["expected_result_evaluation_method_matches"] is False

    stale_result_artifact = solver_result_table_metadata_gate(
        {**metadata, "result_artifact_id": "comsol_slot345_old_result_package_v0"},
        required_columns=("freq_Hz", "Z11_ohm"),
        expected_result_artifact_id="comsol_slot345_result_table_package_v1",
        expected_comsol_version="COMSOL 6.4.0.378",
        require_result_provenance=True,
    )
    assert stale_result_artifact["status"] == "needs_attention"
    assert stale_result_artifact["checks"]["result_artifact_id_recorded_when_required"] is True
    assert stale_result_artifact["checks"]["expected_result_artifact_id_matches"] is False

    sparse_timing = solver_result_table_metadata_gate(
        {**metadata, "timing_breakdown_s": {"solve_s": 0.85}, "execution": {}},
        required_columns=("freq_Hz", "Z11_ohm"),
        expected_result_artifact_id="comsol_slot345_result_table_package_v1",
        expected_comsol_version="COMSOL 6.4.0.378",
        require_result_provenance=True,
    )
    assert sparse_timing["status"] == "needs_attention"
    assert sparse_timing["checks"]["timing_breakdown_recorded_when_required"] is True
    assert sparse_timing["checks"]["timing_breakdown_has_at_least_four_items"] is False

    missing_run_started = solver_result_table_metadata_gate(
        {
            key: value
            for key, value in metadata.items()
            if key not in {"run_started_at", "execution"}
        },
        required_columns=("freq_Hz", "Z11_ohm"),
        expected_result_artifact_id="comsol_slot345_result_table_package_v1",
        expected_comsol_version="COMSOL 6.4.0.378",
        require_result_provenance=True,
    )
    assert missing_run_started["status"] == "needs_attention"
    assert missing_run_started["checks"]["run_started_at_recorded_when_required"] is False


def test_femm_static_current_circuit_rows_gate_requires_solver_ready_rows():
    theta = math.radians(17.5)
    currents = dq_to_three_phase_currents(-2.5, 11.0, theta)
    rows = {
        "U": {
            "circuit_name": "phase_U",
            "current_A": currents["U"],
            "turns": 18,
            "current_kind": "instantaneous",
        },
        "V": {
            "circuit_name": "phase_V",
            "current_A": currents["V"],
            "turns": -18,
            "current_kind": "instantaneous",
        },
        "W": {
            "circuit_name": "phase_W",
            "current_A": currents["W"],
            "turns": 18,
            "current_kind": "instantaneous",
        },
    }

    gate = femm_static_current_circuit_rows_gate(
        currents,
        theta,
        rows,
        expected_id=-2.5,
        expected_iq=11.0,
    )

    assert gate["status"] == "ok"
    assert gate["policy"] == "femm_static_current_circuit_rows_gate"
    assert gate["dq"]["id"] == pytest.approx(-2.5)
    assert gate["dq"]["iq"] == pytest.approx(11.0)
    assert gate["checks"]["turns_nonzero"] is True
    assert gate["max_circuit_current_abs_error_A"] == pytest.approx(0.0)

    rms_rows = {phase: dict(row, current_kind="rms") for phase, row in rows.items()}
    rms_gate = femm_static_current_circuit_rows_gate(
        currents,
        theta,
        rms_rows,
        expected_id=-2.5,
        expected_iq=11.0,
    )
    assert rms_gate["status"] == "needs_attention"
    assert rms_gate["dq"]["status"] == "ok"
    assert rms_gate["checks"]["current_kind_matches"] is False

    missing_rows = dict(rows)
    del missing_rows["W"]
    missing_gate = femm_static_current_circuit_rows_gate(currents, theta, missing_rows)
    assert missing_gate["status"] == "needs_attention"
    assert missing_gate["checks"]["phase_set_ok"] is False

    mismatch_rows = {phase: dict(row) for phase, row in rows.items()}
    mismatch_rows["V"]["current_A"] += 1.0e-3
    mismatch_gate = femm_static_current_circuit_rows_gate(currents, theta, mismatch_rows)
    assert mismatch_gate["status"] == "needs_attention"
    assert mismatch_gate["checks"]["circuit_currents_match"] is False


def test_femm_block_label_source_contract_gate_keeps_sources_explicit():
    rows = [
        {"region": "air", "material": "Air", "group": 0, "source_kind": "air"},
        {"region": "stator", "material": "M-19 Steel", "group": 1, "source_kind": "passive"},
        {"region": "pm_N", "material": "NdFeB 40", "group": 2, "source_kind": "pm", "magnetization_deg": 0.0},
        {"region": "phase_U_plus", "material": "Copper", "group": 10, "source_kind": "coil", "circuit_name": "phase_U", "turns": 18},
        {"region": "phase_U_minus", "material": "Copper", "group": 11, "source_kind": "coil", "circuit_name": "phase_U", "turns": -18},
    ]

    gate = femm_block_label_source_contract_gate(
        rows,
        required_regions=["air", "stator", "pm_N", "phase_U_plus", "phase_U_minus"],
    )

    assert gate["status"] == "ok"
    assert gate["policy"] == "femm_block_label_source_contract_gate"
    assert gate["source_counts"] == {"air": 1, "coil": 2, "passive": 1, "pm": 1}
    assert gate["checks"]["coil_rows_have_nonzero_turns"] is True
    assert gate["checks"]["pm_rows_have_magnetization_direction"] is True
    assert "coils and PMs stay explicit sources" in gate["version_note"]

    bad = femm_block_label_source_contract_gate(
        [
            {"region": "air", "material": "Air", "group": 0, "source_kind": "air", "circuit_name": "phase_U"},
            {"region": "pm_N", "material": "NdFeB 40", "group": 2, "source_kind": "pm"},
            {"region": "phase_U_plus", "material": "Copper", "group": 10, "source_kind": "coil", "turns": 0},
            {"region": "phase_U_plus", "material": "", "source_kind": "coil", "circuit_name": "phase_U", "turns": 18},
            {"region": "mystery", "material": "Steel", "group": 99, "source_kind": "current_sheet"},
        ],
        required_regions=["air", "stator", "pm_N", "phase_U_plus"],
    )

    assert bad["status"] == "needs_attention"
    assert bad["missing_required_regions"] == ["stator"]
    assert bad["duplicate_regions"] == ["phase_U_plus"]
    assert bad["unknown_source_kinds"] == ["current_sheet"]
    assert bad["checks"]["materials_present"] is False
    assert bad["checks"]["groups_present"] is False
    assert bad["checks"]["coil_rows_have_circuit"] is False
    assert bad["checks"]["coil_rows_have_nonzero_turns"] is False
    assert bad["checks"]["pm_rows_have_magnetization_direction"] is False
    assert bad["checks"]["air_passive_rows_have_no_source_metadata"] is False


def test_femm_group_motion_selection_gate_keeps_rotor_entities_moving_together():
    rows = [
        {
            "entity_kind": "block_label",
            "name": "rotor_pm_N",
            "group_id": 7,
            "selected_for_motion": True,
            "motion_command": "mi_moverotate",
        },
        {
            "entity_kind": "segment",
            "entity_id": 12,
            "group_id": 7,
            "selected_for_motion": True,
            "motion_command": "mi_moverotate",
        },
        {
            "entity_kind": "arc_segment",
            "entity_id": 4,
            "group_id": 7,
            "selected_for_motion": True,
            "motion_command": "mi_moverotate",
        },
    ]

    gate = femm_group_motion_selection_gate(rows, expected_group_id=7)

    assert gate["policy"] == "femm_group_motion_selection_gate"
    assert gate["status"] == "ok"
    assert gate["checks"]["required_entity_kinds_present"] is True
    assert gate["checks"]["all_rows_use_expected_group"] is True
    assert gate["checks"]["no_group_zero_rows"] is True

    bad_rows = [dict(row) for row in rows]
    bad_rows[1]["group_id"] = 0
    bad_rows[2]["selected_for_motion"] = False
    bad = femm_group_motion_selection_gate(bad_rows, expected_group_id=7)
    assert bad["status"] == "needs_attention"
    assert bad["checks"]["all_rows_use_expected_group"] is False
    assert bad["checks"]["no_group_zero_rows"] is False
    assert bad["checks"]["all_rows_selected_for_motion"] is False


def test_femm_pm_magnetization_convention_gate_requires_degree_frame_and_strength():
    rows = [
        {
            "region": "pm_N",
            "magdir_deg": 0.0,
            "angle_unit": "deg",
            "frame": "global_xy",
            "Br_T": 1.2,
        },
        {
            "region": "pm_S",
            "magnetization_deg": 180.0,
            "angle_unit": "degrees",
            "coordinate_frame": "global_xy",
            "Hc_A_per_m": 900000.0,
        },
    ]

    gate = femm_pm_magnetization_convention_gate(
        rows,
        required_regions=["pm_N", "pm_S"],
    )

    assert gate["status"] == "ok"
    assert gate["policy"] == "femm_pm_magnetization_convention_gate"
    assert gate["rows"][0]["unit_vector_xy"] == pytest.approx([1.0, 0.0], abs=1.0e-15)
    assert gate["rows"][1]["unit_vector_xy"] == pytest.approx([-1.0, 0.0], abs=1.0e-15)
    assert gate["checks"]["angles_are_degrees_and_finite"] is True
    assert gate["checks"]["strength_present"] is True
    assert "FEMM magdir must be degrees" in gate["version_note"]

    bad = femm_pm_magnetization_convention_gate(
        [
            {"region": "pm_N", "magdir_deg": 0.0, "angle_unit": "rad", "frame": "global_xy", "Br_T": 1.2},
            {"region": "pm_N", "magdir_deg": "east", "frame": "rotor_xy", "Hc_A_per_m": 900000.0},
            {"region": "pm_S", "magdir_deg": 180.0, "frame": "screen_xy"},
        ],
        required_regions=["pm_N", "pm_S", "pm_aux"],
    )

    assert bad["status"] == "needs_attention"
    assert bad["duplicate_regions"] == ["pm_N"]
    assert bad["missing_required_regions"] == ["pm_aux"]
    assert bad["bad_angle_regions"] == ["pm_N", "pm_N"]
    assert bad["bad_frame_regions"] == ["pm_S"]
    assert bad["missing_strength_regions"] == ["pm_S"]


def test_jmag_motor_table_column_metadata_gate_closes_units_before_values():
    metadata = {
        "columns": ["RotorAngle_deg", "Torque_Nm", "Id_A", "Iq_A", "Speed_rpm"],
        "angle_column": "RotorAngle_deg",
        "torque_column": "Torque_Nm",
        "angle_unit": "deg",
        "angle_basis": "mechanical",
        "pole_pairs": 4,
        "symmetry_factor": 6,
        "current_basis": "peak",
        "torque_sign_convention": "positive_motoring",
    }

    gate = jmag_motor_table_column_metadata_gate(
        metadata,
        required_columns=["RotorAngle_deg", "Torque_Nm", "Id_A", "Iq_A"],
    )

    assert gate["status"] == "ok"
    assert gate["policy"] == "jmag_motor_table_column_metadata_gate"
    assert gate["angle_unit"] == "deg"
    assert gate["angle_basis"] == "mechanical"
    assert gate["checks"]["required_columns_present"] is True
    assert gate["checks"]["current_basis_valid"] is True
    assert "before JMAG torque/current/efficiency table parsing" in gate["version_note"]

    bad = jmag_motor_table_column_metadata_gate(
        {
            "columns": ["RotorAngle", "Id_A"],
            "angle_column": "RotorAngle_deg",
            "torque_column": "Torque_Nm",
            "angle_unit": "degreee",
            "angle_basis": "rotor",
            "pole_pairs": 0,
            "symmetry_factor": 0,
            "current_basis": "line_rms",
            "torque_sign_convention": "unknown",
        },
        required_columns=["RotorAngle_deg", "Torque_Nm", "Id_A", "Iq_A"],
    )

    assert bad["status"] == "needs_attention"
    assert bad["missing_required_columns"] == ["RotorAngle_deg", "Torque_Nm", "Iq_A"]
    assert bad["checks"]["angle_column_present"] is False
    assert bad["checks"]["angle_unit_valid"] is False
    assert bad["checks"]["angle_basis_valid"] is False
    assert bad["checks"]["pole_pairs_positive"] is False
    assert bad["checks"]["symmetry_factor_positive"] is False
    assert bad["checks"]["current_basis_valid"] is False
    assert bad["checks"]["torque_sign_convention_valid"] is False


def test_jmag_force_table_metadata_gate_closes_units_frame_and_sign_before_values():
    metadata = {
        "columns": [
            "Gap_mm",
            "Fx_N_per_m",
            "Fy_N_per_m",
            "CaseId",
            "OperatingPointId",
            "TargetBodyId",
            "TargetBodyName",
            "TargetMaterialName",
            "TargetRegionArtifactId",
            "AnalysisType",
            "FrequencyHz",
        ],
        "position_columns": ["Gap_mm"],
        "force_columns": ["Fx_N_per_m", "Fy_N_per_m"],
        "identity_columns": ["TargetBodyId"],
        "case_id": "case_gap_force_sweep_A",
        "study_id": "magnetostatic_force",
        "analysis_type": "magnetostatic",
        "frequency_hz": 0.0,
        "operating_point_id": "op_gap_sweep_A",
        "force_unit": "N/m",
        "position_unit": "mm",
        "component_frame": "global_xy",
        "projection_axis": "moving_core_to_fixed_core_gap_normal_positive_attraction",
        "source_tool": "JMAG-Designer",
        "export_artifact_id": "jmag_force_export_artifact_A",
        "result_set_id": "resultset_force_20260630_A",
        "mesh_id": "jmag_mesh_gap_force_A_v3",
        "solver_run_id": "jmag_solver_run_20260630_A",
        "result_revision_id": "jmag_result_revision_force_A_001",
        "solver_setup_artifact_id": "jmag_solver_setup_force_A_strict_iccg_nr_v1.json",
        "material_state_artifact_id": "jmag_material_state_force_A_bh_curve_v2.json",
        "material_state_digest": "sha256:jmag_slot333_material_state_digest_v1",
        "excitation_source_artifact_id": "jmag_excitation_source_force_A_peak_current_v1.json",
        "current_definition_method": "peak_phase_current_table",
        "export_trace_id": "jmag_slot285_force_export_macro_trace_v1",
        "export_command_digest": "sha256:jmag_slot285_force_write_all_case_table",
        "export_commands": [
            "select result set resultset_force_20260630_A",
            "open force report moving_core",
            "WriteAllCaseTable force_table.csv",
        ],
        "export_output_artifact_id": "jmag_slot293_force_table_csv_v1",
        "export_output_digest": "sha256:jmag_slot293_force_table_csv",
        "export_output_path": "artifacts/motor/jmag_slot293_force_table.csv",
        "force_observable_id": "jmag_slot301_maxwell_stress_force_report_xy_v1",
        "force_observable_family": "jmag_force_report_maxwell_stress_xy",
        "force_report_method": "maxwell_stress_force_report_xy",
        "target_region_id": "moving_core",
        "target_region_name": "moving_core_block_label",
        "target_material": "magnetic_steel_core",
        "target_region_artifact_id": "jmag_region_labels_force_A_v2.json",
        "target_region_geometry_digest": "sha256:jmag_slot325_moving_core_geometry_digest_v1",
        "target_region_centroid_xyz_m": [0.012, 0.0, 0.0],
        "symmetry_factor": 8,
        "force_sign_convention": "positive_attraction",
        "force_kind": "maxwell_stress",
        "quantity_dimension": "2d_per_length",
    }

    gate = jmag_force_table_metadata_gate(
        metadata,
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_case_id="case_gap_force_sweep_A",
        expected_study_id="magnetostatic_force",
        expected_operating_point_id="op_gap_sweep_A",
        expected_analysis_type="magnetostatic",
        expected_frequency_hz=0.0,
        expected_export_artifact_id="jmag_force_export_artifact_A",
        expected_result_set_id="resultset_force_20260630_A",
        expected_mesh_id="jmag_mesh_gap_force_A_v3",
        expected_solver_run_id="jmag_solver_run_20260630_A",
        expected_result_revision_id="jmag_result_revision_force_A_001",
        expected_solver_setup_artifact_id="jmag_solver_setup_force_A_strict_iccg_nr_v1.json",
        expected_material_state_artifact_id="jmag_material_state_force_A_bh_curve_v2.json",
        expected_material_state_digest="sha256:jmag_slot333_material_state_digest_v1",
        expected_excitation_source_artifact_id="jmag_excitation_source_force_A_peak_current_v1.json",
        expected_current_definition_method="peak_phase_current_table",
        expected_export_trace_id="jmag_slot285_force_export_macro_trace_v1",
        expected_export_command_digest="sha256:jmag_slot285_force_write_all_case_table",
        expected_export_output_artifact_id="jmag_slot293_force_table_csv_v1",
        expected_export_output_digest="sha256:jmag_slot293_force_table_csv",
        expected_force_observable_id="jmag_slot301_maxwell_stress_force_report_xy_v1",
        expected_force_observable_family="jmag_force_report_maxwell_stress_xy",
        expected_component_frame="global_xy",
        expected_projection_axis="moving_core_to_fixed_core_gap_normal_positive_attraction",
        expected_force_sign_convention="positive_attraction",
        expected_force_report_method="maxwell_stress_force_report_xy",
        require_export_command_trace=True,
        require_export_output_artifact=True,
        expected_target_region_id="moving_core",
        expected_target_region_name="moving_core_block_label",
        expected_target_material="magnetic_steel_core",
        expected_target_region_artifact_id="jmag_region_labels_force_A_v2.json",
        expected_target_region_geometry_digest="sha256:jmag_slot325_moving_core_geometry_digest_v1",
        expected_target_region_centroid_xyz_m=(0.012, 0.0, 0.0),
    )

    assert gate["status"] == "ok"
    assert gate["policy"] == "jmag_force_table_metadata_gate"
    assert gate["checks"]["force_columns_present"] is True
    assert gate["checks"]["force_unit_valid"] is True
    assert gate["checks"]["quantity_dimension_valid"] is True
    assert gate["checks"]["force_unit_matches_quantity_dimension_when_present"] is True
    assert gate["quantity_dimension"] == "2d_per_length"
    assert gate["checks"]["component_frame_valid"] is True
    assert gate["checks"]["component_frame_recorded_when_expected"] is True
    assert gate["checks"]["expected_component_frame_matches"] is True
    assert gate["checks"]["force_columns_match_component_frame"] is True
    assert gate["checks"]["force_projection_axis_recorded"] is True
    assert gate["checks"]["force_projection_axis_descriptive"] is True
    assert gate["checks"]["projection_axis_recorded_when_expected"] is True
    assert gate["checks"]["expected_projection_axis_matches"] is True
    assert gate["projection_axis"] == "moving_core_to_fixed_core_gap_normal_positive_attraction"
    assert gate["force_column_axes_present"]["x"] is True
    assert gate["force_column_axes_present"]["y"] is True
    assert gate["expected_force_axes_for_frame"] == ["x", "y"]
    assert gate["checks"]["source_tool_is_jmag"] is True
    assert gate["checks"]["force_sign_convention_recorded_when_expected"] is True
    assert gate["checks"]["expected_force_sign_convention_matches"] is True
    assert gate["expected_component_frame"] == "global_xy"
    assert gate["expected_projection_axis"] == "moving_core_to_fixed_core_gap_normal_positive_attraction"
    assert gate["expected_force_sign_convention"] == "positive_attraction"
    assert gate["case_id"] == "case_gap_force_sweep_A"
    assert gate["study_id"] == "magnetostatic_force"
    assert gate["analysis_type"] == "magnetostatic"
    assert gate["frequency_hz"] == pytest.approx(0.0)
    assert gate["operating_point_id"] == "op_gap_sweep_A"
    assert gate["checks"]["case_id_recorded"] is True
    assert gate["checks"]["study_id_recorded"] is True
    assert gate["checks"]["operating_point_id_recorded"] is True
    assert gate["checks"]["expected_case_id_matches"] is True
    assert gate["checks"]["expected_study_id_matches"] is True
    assert gate["checks"]["expected_operating_point_id_matches"] is True
    assert gate["checks"]["analysis_type_recorded"] is True
    assert gate["checks"]["expected_analysis_type_matches"] is True
    assert gate["checks"]["frequency_hz_recorded"] is True
    assert gate["checks"]["expected_frequency_hz_matches"] is True
    assert gate["export_artifact_id"] == "jmag_force_export_artifact_A"
    assert gate["checks"]["export_artifact_id_recorded"] is True
    assert gate["checks"]["expected_export_artifact_id_matches"] is True
    assert gate["checks"]["expected_result_set_id_matches"] is True
    assert gate["mesh_id"] == "jmag_mesh_gap_force_A_v3"
    assert gate["solver_run_id"] == "jmag_solver_run_20260630_A"
    assert gate["result_revision_id"] == "jmag_result_revision_force_A_001"
    assert gate["solver_setup_artifact_id"] == "jmag_solver_setup_force_A_strict_iccg_nr_v1.json"
    assert gate["material_state_artifact_id"] == "jmag_material_state_force_A_bh_curve_v2.json"
    assert gate["material_state_digest"] == "sha256:jmag_slot333_material_state_digest_v1"
    assert gate["excitation_source_artifact_id"] == "jmag_excitation_source_force_A_peak_current_v1.json"
    assert gate["current_definition_method"] == "peak_phase_current_table"
    assert gate["checks"]["mesh_id_recorded"] is True
    assert gate["checks"]["solver_run_id_recorded"] is True
    assert gate["checks"]["result_revision_id_recorded"] is True
    assert gate["checks"]["solver_setup_artifact_id_recorded"] is True
    assert gate["checks"]["material_state_artifact_id_recorded"] is True
    assert gate["checks"]["material_state_digest_recorded"] is True
    assert gate["checks"]["excitation_source_artifact_id_recorded"] is True
    assert gate["checks"]["current_definition_method_recorded"] is True
    assert gate["checks"]["expected_mesh_id_matches"] is True
    assert gate["checks"]["expected_solver_run_id_matches"] is True
    assert gate["checks"]["expected_result_revision_id_matches"] is True
    assert gate["checks"]["expected_solver_setup_artifact_id_matches"] is True
    assert gate["checks"]["expected_material_state_artifact_id_matches"] is True
    assert gate["checks"]["expected_material_state_digest_matches"] is True
    assert gate["checks"]["expected_excitation_source_artifact_id_matches"] is True
    assert gate["checks"]["expected_current_definition_method_matches"] is True
    assert gate["export_trace_id"] == "jmag_slot285_force_export_macro_trace_v1"
    assert gate["export_command_digest"] == "sha256:jmag_slot285_force_write_all_case_table"
    assert gate["checks"]["export_trace_id_recorded"] is True
    assert gate["checks"]["expected_export_trace_id_matches"] is True
    assert gate["checks"]["export_command_digest_recorded"] is True
    assert gate["checks"]["expected_export_command_digest_matches"] is True
    assert gate["checks"]["export_commands_recorded"] is True
    assert gate["checks"]["export_commands_include_table_export"] is True
    assert gate["checks"]["export_commands_reference_force_report"] is True
    assert gate["export_output_artifact_id"] == "jmag_slot293_force_table_csv_v1"
    assert gate["export_output_digest"] == "sha256:jmag_slot293_force_table_csv"
    assert gate["export_output_path"] == "artifacts/motor/jmag_slot293_force_table.csv"
    assert gate["checks"]["export_output_artifact_id_recorded"] is True
    assert gate["checks"]["expected_export_output_artifact_id_matches"] is True
    assert gate["checks"]["export_output_digest_recorded"] is True
    assert gate["checks"]["expected_export_output_digest_matches"] is True
    assert gate["checks"]["export_output_path_recorded"] is True
    assert gate["force_observable_id"] == "jmag_slot301_maxwell_stress_force_report_xy_v1"
    assert gate["force_observable_family"] == "jmag_force_report_maxwell_stress_xy"
    assert gate["checks"]["force_observable_id_recorded"] is True
    assert gate["checks"]["expected_force_observable_id_matches"] is True
    assert gate["checks"]["force_observable_family_recorded"] is True
    assert gate["checks"]["expected_force_observable_family_matches"] is True
    assert gate["force_report_method"] == "maxwell_stress_force_report_xy"
    assert gate["expected_force_report_method"] == "maxwell_stress_force_report_xy"
    assert gate["checks"]["force_report_method_recorded_when_expected"] is True
    assert gate["checks"]["expected_force_report_method_matches"] is True
    assert gate["checks"]["expected_target_region_id_matches"] is True
    assert gate["target_region_name"] == "moving_core_block_label"
    assert gate["target_material"] == "magnetic_steel_core"
    assert gate["target_region_artifact_id"] == "jmag_region_labels_force_A_v2.json"
    assert gate["target_region_geometry_digest"] == "sha256:jmag_slot325_moving_core_geometry_digest_v1"
    assert gate["target_region_centroid_xyz_m"] == [0.012, 0.0, 0.0]
    assert gate["checks"]["target_region_name_recorded"] is True
    assert gate["checks"]["expected_target_region_name_matches"] is True
    assert gate["checks"]["target_material_recorded"] is True
    assert gate["checks"]["expected_target_material_matches"] is True
    assert gate["checks"]["target_region_artifact_id_recorded"] is True
    assert gate["checks"]["expected_target_region_artifact_id_matches"] is True
    assert gate["checks"]["target_region_geometry_digest_recorded"] is True
    assert gate["checks"]["expected_target_region_geometry_digest_matches"] is True
    assert gate["checks"]["target_region_centroid_xyz_recorded_when_expected"] is True
    assert gate["checks"]["expected_target_region_centroid_xyz_matches"] is True
    assert gate["checks"]["identity_columns_present"] is True
    assert gate["target_region_id"] == "moving_core"
    assert "before JMAG force/NVH/contact-force table parsing" in gate["version_note"]

    stale_mesh_metadata = jmag_force_table_metadata_gate(
        {**metadata, "mesh_id": "jmag_mesh_gap_force_old"},
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_mesh_id="jmag_mesh_gap_force_A_v3",
    )
    assert stale_mesh_metadata["status"] == "needs_attention"
    assert stale_mesh_metadata["checks"]["expected_mesh_id_matches"] is False

    stale_solver_run_metadata = jmag_force_table_metadata_gate(
        {**metadata, "solver_run_id": "jmag_solver_run_old"},
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_solver_run_id="jmag_solver_run_20260630_A",
    )
    assert stale_solver_run_metadata["status"] == "needs_attention"
    assert stale_solver_run_metadata["checks"]["expected_solver_run_id_matches"] is False

    stale_result_revision_metadata = jmag_force_table_metadata_gate(
        {**metadata, "result_revision_id": "jmag_result_revision_old"},
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_result_revision_id="jmag_result_revision_force_A_001",
    )
    assert stale_result_revision_metadata["status"] == "needs_attention"
    assert stale_result_revision_metadata["checks"]["expected_result_revision_id_matches"] is False

    stale_solver_setup_metadata = jmag_force_table_metadata_gate(
        {**metadata, "solver_setup_artifact_id": "jmag_solver_setup_force_A_loose_auto_v1.json"},
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_solver_setup_artifact_id="jmag_solver_setup_force_A_strict_iccg_nr_v1.json",
        expected_material_state_artifact_id="jmag_material_state_force_A_bh_curve_v2.json",
        expected_material_state_digest="sha256:jmag_slot333_material_state_digest_v1",
    )
    assert stale_solver_setup_metadata["status"] == "needs_attention"
    assert stale_solver_setup_metadata["checks"]["expected_solver_setup_artifact_id_matches"] is False
    assert stale_solver_setup_metadata["checks"]["expected_material_state_artifact_id_matches"] is True
    assert stale_solver_setup_metadata["checks"]["expected_material_state_digest_matches"] is True

    stale_material_state_metadata = jmag_force_table_metadata_gate(
        {**metadata, "material_state_artifact_id": "jmag_material_state_force_A_old.json"},
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_solver_setup_artifact_id="jmag_solver_setup_force_A_strict_iccg_nr_v1.json",
        expected_material_state_artifact_id="jmag_material_state_force_A_bh_curve_v2.json",
        expected_material_state_digest="sha256:jmag_slot333_material_state_digest_v1",
    )
    assert stale_material_state_metadata["status"] == "needs_attention"
    assert stale_material_state_metadata["checks"]["expected_solver_setup_artifact_id_matches"] is True
    assert stale_material_state_metadata["checks"]["expected_material_state_artifact_id_matches"] is False
    assert stale_material_state_metadata["checks"]["expected_material_state_digest_matches"] is True

    stale_material_digest_metadata = jmag_force_table_metadata_gate(
        {**metadata, "material_state_digest": "sha256:jmag_slot333_material_state_digest_old"},
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_solver_setup_artifact_id="jmag_solver_setup_force_A_strict_iccg_nr_v1.json",
        expected_material_state_artifact_id="jmag_material_state_force_A_bh_curve_v2.json",
        expected_material_state_digest="sha256:jmag_slot333_material_state_digest_v1",
    )
    assert stale_material_digest_metadata["status"] == "needs_attention"
    assert stale_material_digest_metadata["checks"]["expected_solver_setup_artifact_id_matches"] is True
    assert stale_material_digest_metadata["checks"]["expected_material_state_artifact_id_matches"] is True
    assert stale_material_digest_metadata["checks"]["expected_material_state_digest_matches"] is False

    stale_excitation_source_metadata = jmag_force_table_metadata_gate(
        {**metadata, "excitation_source_artifact_id": "jmag_excitation_source_force_A_old.json"},
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_excitation_source_artifact_id="jmag_excitation_source_force_A_peak_current_v1.json",
        expected_current_definition_method="peak_phase_current_table",
    )
    assert stale_excitation_source_metadata["status"] == "needs_attention"
    assert stale_excitation_source_metadata["checks"]["expected_excitation_source_artifact_id_matches"] is False
    assert stale_excitation_source_metadata["checks"]["expected_current_definition_method_matches"] is True

    wrong_current_definition_metadata = jmag_force_table_metadata_gate(
        {**metadata, "current_definition_method": "rms_phase_current_table"},
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_excitation_source_artifact_id="jmag_excitation_source_force_A_peak_current_v1.json",
        expected_current_definition_method="peak_phase_current_table",
    )
    assert wrong_current_definition_metadata["status"] == "needs_attention"
    assert wrong_current_definition_metadata["checks"]["expected_excitation_source_artifact_id_matches"] is True
    assert wrong_current_definition_metadata["checks"]["expected_current_definition_method_matches"] is False

    missing_current_definition_metadata = jmag_force_table_metadata_gate(
        {key: value for key, value in metadata.items() if key != "current_definition_method"},
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_current_definition_method="peak_phase_current_table",
    )
    assert missing_current_definition_metadata["status"] == "needs_attention"
    assert missing_current_definition_metadata["checks"]["current_definition_method_recorded"] is False
    assert missing_current_definition_metadata["checks"]["expected_current_definition_method_matches"] is False

    stale_export_trace_metadata = jmag_force_table_metadata_gate(
        metadata,
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_export_trace_id="jmag_slot285_force_export_macro_trace_old",
        expected_export_command_digest="sha256:jmag_slot285_force_write_all_case_table",
        require_export_command_trace=True,
    )
    assert stale_export_trace_metadata["status"] == "needs_attention"
    assert stale_export_trace_metadata["checks"]["expected_export_trace_id_matches"] is False
    assert stale_export_trace_metadata["checks"]["expected_export_command_digest_matches"] is True
    assert stale_export_trace_metadata["checks"]["export_commands_include_table_export"] is True

    missing_table_export_command = jmag_force_table_metadata_gate(
        {
            **metadata,
            "export_commands": [
                "select result set resultset_force_20260630_A",
                "open force report moving_core",
            ],
        },
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_export_trace_id="jmag_slot285_force_export_macro_trace_v1",
        expected_export_command_digest="sha256:jmag_slot285_force_write_all_case_table",
        require_export_command_trace=True,
    )
    assert missing_table_export_command["status"] == "needs_attention"
    assert missing_table_export_command["checks"]["expected_export_trace_id_matches"] is True
    assert missing_table_export_command["checks"]["export_commands_recorded"] is True
    assert missing_table_export_command["checks"]["export_commands_include_table_export"] is False
    assert missing_table_export_command["checks"]["export_commands_reference_force_report"] is True

    stale_export_output_artifact = jmag_force_table_metadata_gate(
        metadata,
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_export_output_artifact_id="jmag_slot293_force_table_csv_old",
        expected_export_output_digest="sha256:jmag_slot293_force_table_csv",
        require_export_output_artifact=True,
    )
    assert stale_export_output_artifact["status"] == "needs_attention"
    assert stale_export_output_artifact["checks"]["expected_export_output_artifact_id_matches"] is False
    assert stale_export_output_artifact["checks"]["expected_export_output_digest_matches"] is True
    assert stale_export_output_artifact["checks"]["export_output_path_recorded"] is True

    stale_export_output_digest = jmag_force_table_metadata_gate(
        metadata,
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_export_output_artifact_id="jmag_slot293_force_table_csv_v1",
        expected_export_output_digest="sha256:jmag_slot293_force_table_csv_old",
        require_export_output_artifact=True,
    )
    assert stale_export_output_digest["status"] == "needs_attention"
    assert stale_export_output_digest["checks"]["expected_export_output_artifact_id_matches"] is True
    assert stale_export_output_digest["checks"]["expected_export_output_digest_matches"] is False

    stale_force_observable = jmag_force_table_metadata_gate(
        metadata,
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_force_observable_id="jmag_slot212_old_force_report_xy",
        expected_force_observable_family="jmag_force_report_maxwell_stress_xy",
        expected_export_output_artifact_id="jmag_slot293_force_table_csv_v1",
        expected_export_output_digest="sha256:jmag_slot293_force_table_csv",
    )
    assert stale_force_observable["status"] == "needs_attention"
    assert stale_force_observable["checks"]["expected_force_observable_id_matches"] is False
    assert stale_force_observable["checks"]["expected_force_observable_family_matches"] is True
    assert stale_force_observable["checks"]["expected_export_output_artifact_id_matches"] is True
    assert stale_force_observable["checks"]["expected_export_output_digest_matches"] is True

    wrong_force_observable_family = jmag_force_table_metadata_gate(
        {**metadata, "force_observable_family": "jmag_torque_report"},
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_force_observable_id="jmag_slot301_maxwell_stress_force_report_xy_v1",
        expected_force_observable_family="jmag_force_report_maxwell_stress_xy",
    )
    assert wrong_force_observable_family["status"] == "needs_attention"
    assert wrong_force_observable_family["checks"]["expected_force_observable_id_matches"] is True
    assert wrong_force_observable_family["checks"]["expected_force_observable_family_matches"] is False
    assert wrong_force_observable_family["checks"]["force_kind_valid"] is True

    wrong_force_report_method = jmag_force_table_metadata_gate(
        {**metadata, "force_report_method": "nodal_force_report_xy"},
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_force_observable_id="jmag_slot301_maxwell_stress_force_report_xy_v1",
        expected_force_observable_family="jmag_force_report_maxwell_stress_xy",
        expected_force_report_method="maxwell_stress_force_report_xy",
    )
    assert wrong_force_report_method["status"] == "needs_attention"
    assert wrong_force_report_method["checks"]["expected_force_observable_id_matches"] is True
    assert wrong_force_report_method["checks"]["expected_force_observable_family_matches"] is True
    assert wrong_force_report_method["checks"]["force_report_method_recorded_when_expected"] is True
    assert wrong_force_report_method["checks"]["expected_force_report_method_matches"] is False

    missing_force_report_method = jmag_force_table_metadata_gate(
        {key: value for key, value in metadata.items() if key != "force_report_method"},
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_force_report_method="maxwell_stress_force_report_xy",
    )
    assert missing_force_report_method["status"] == "needs_attention"
    assert missing_force_report_method["checks"]["force_report_method_recorded_when_expected"] is False
    assert missing_force_report_method["checks"]["expected_force_report_method_matches"] is False

    wrong_expected_frame = jmag_force_table_metadata_gate(
        {**metadata, "component_frame": "as_exported"},
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_component_frame="global_xy",
        expected_projection_axis="moving_core_to_fixed_core_gap_normal_positive_attraction",
        expected_force_sign_convention="positive_attraction",
    )
    assert wrong_expected_frame["status"] == "needs_attention"
    assert wrong_expected_frame["checks"]["component_frame_valid"] is True
    assert wrong_expected_frame["checks"]["force_columns_match_component_frame"] is True
    assert wrong_expected_frame["checks"]["expected_component_frame_matches"] is False
    assert wrong_expected_frame["checks"]["expected_projection_axis_matches"] is True
    assert wrong_expected_frame["checks"]["expected_force_sign_convention_matches"] is True

    wrong_expected_projection = jmag_force_table_metadata_gate(
        {**metadata, "projection_axis": "fixed_core_to_moving_core_gap_normal_positive_repulsion"},
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_component_frame="global_xy",
        expected_projection_axis="moving_core_to_fixed_core_gap_normal_positive_attraction",
        expected_force_sign_convention="positive_attraction",
    )
    assert wrong_expected_projection["status"] == "needs_attention"
    assert wrong_expected_projection["checks"]["force_projection_axis_descriptive"] is True
    assert wrong_expected_projection["checks"]["expected_component_frame_matches"] is True
    assert wrong_expected_projection["checks"]["expected_projection_axis_matches"] is False
    assert wrong_expected_projection["checks"]["expected_force_sign_convention_matches"] is True

    wrong_expected_sign = jmag_force_table_metadata_gate(
        {**metadata, "force_sign_convention": "positive_repulsion"},
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_component_frame="global_xy",
        expected_projection_axis="moving_core_to_fixed_core_gap_normal_positive_attraction",
        expected_force_sign_convention="positive_attraction",
    )
    assert wrong_expected_sign["status"] == "needs_attention"
    assert wrong_expected_sign["checks"]["force_sign_convention_valid"] is True
    assert wrong_expected_sign["checks"]["expected_component_frame_matches"] is True
    assert wrong_expected_sign["checks"]["expected_projection_axis_matches"] is True
    assert wrong_expected_sign["checks"]["expected_force_sign_convention_matches"] is False

    stale_target_name_metadata = jmag_force_table_metadata_gate(
        {**metadata, "target_region_name": "fixed_core_block_label"},
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_target_region_name="moving_core_block_label",
    )
    assert stale_target_name_metadata["status"] == "needs_attention"
    assert stale_target_name_metadata["checks"]["expected_target_region_name_matches"] is False

    stale_target_material_metadata = jmag_force_table_metadata_gate(
        {**metadata, "target_material": "air"},
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_target_material="magnetic_steel_core",
    )
    assert stale_target_material_metadata["status"] == "needs_attention"
    assert stale_target_material_metadata["checks"]["expected_target_material_matches"] is False

    stale_target_artifact_metadata = jmag_force_table_metadata_gate(
        {**metadata, "target_region_artifact_id": "jmag_region_labels_old.json"},
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_target_region_artifact_id="jmag_region_labels_force_A_v2.json",
    )
    assert stale_target_artifact_metadata["status"] == "needs_attention"
    assert stale_target_artifact_metadata["checks"]["expected_target_region_artifact_id_matches"] is False

    stale_target_geometry_digest = jmag_force_table_metadata_gate(
        {**metadata, "target_region_geometry_digest": "sha256:jmag_slot325_moving_core_geometry_old"},
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_target_region_geometry_digest="sha256:jmag_slot325_moving_core_geometry_digest_v1",
        expected_target_region_centroid_xyz_m=(0.012, 0.0, 0.0),
    )
    assert stale_target_geometry_digest["status"] == "needs_attention"
    assert stale_target_geometry_digest["checks"]["expected_target_region_geometry_digest_matches"] is False
    assert stale_target_geometry_digest["checks"]["expected_target_region_centroid_xyz_matches"] is True

    stale_target_centroid = jmag_force_table_metadata_gate(
        {**metadata, "target_region_centroid_xyz_m": [0.013, 0.0, 0.0]},
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_target_region_geometry_digest="sha256:jmag_slot325_moving_core_geometry_digest_v1",
        expected_target_region_centroid_xyz_m=(0.012, 0.0, 0.0),
    )
    assert stale_target_centroid["status"] == "needs_attention"
    assert stale_target_centroid["checks"]["expected_target_region_geometry_digest_matches"] is True
    assert stale_target_centroid["checks"]["expected_target_region_centroid_xyz_matches"] is False

    missing_target_centroid_metadata = dict(metadata)
    missing_target_centroid_metadata.pop("target_region_centroid_xyz_m")
    missing_target_centroid = jmag_force_table_metadata_gate(
        missing_target_centroid_metadata,
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_target_region_centroid_xyz_m=(0.012, 0.0, 0.0),
    )
    assert missing_target_centroid["status"] == "needs_attention"
    assert missing_target_centroid["checks"]["target_region_centroid_xyz_recorded_when_expected"] is False
    assert missing_target_centroid["checks"]["expected_target_region_centroid_xyz_matches"] is False

    stale_export_artifact = jmag_force_table_metadata_gate(
        metadata,
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_export_artifact_id="jmag_force_export_artifact_B",
        expected_result_set_id="resultset_force_20260630_A",
        expected_target_region_id="moving_core",
    )
    assert stale_export_artifact["status"] == "needs_attention"
    assert stale_export_artifact["checks"]["expected_export_artifact_id_matches"] is False
    assert stale_export_artifact["checks"]["expected_result_set_id_matches"] is True

    row_identity = jmag_force_table_metadata_gate(
        {**metadata, "identity_columns": ["TargetBodyId", "Gap_mm"]},
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_export_artifact_id="jmag_force_export_artifact_A",
        expected_result_set_id="resultset_force_20260630_A",
        expected_target_region_id="moving_core",
        min_table_rows=2,
        require_unique_row_identity=True,
        table_rows=[
            {
                "Gap_mm": 1.0,
                "Fx_N_per_m": -4.0,
                "Fy_N_per_m": 0.0,
                "CaseId": "case_gap_force_sweep_A",
                "OperatingPointId": "op_gap_sweep_A",
                "TargetBodyId": "moving_core",
                "TargetBodyName": "moving_core_block_label",
                "TargetMaterialName": "magnetic_steel_core",
                "TargetRegionArtifactId": "jmag_region_labels_force_A_v2.json",
                "AnalysisType": "magnetostatic",
                "FrequencyHz": 0.0,
            },
            {
                "Gap_mm": 2.0,
                "Fx_N_per_m": -2.0,
                "Fy_N_per_m": 0.0,
                "CaseId": "case_gap_force_sweep_A",
                "OperatingPointId": "op_gap_sweep_A",
                "TargetBodyId": "moving_core",
                "TargetBodyName": "moving_core_block_label",
                "TargetMaterialName": "magnetic_steel_core",
                "TargetRegionArtifactId": "jmag_region_labels_force_A_v2.json",
                "AnalysisType": "magnetostatic",
                "FrequencyHz": 0.0,
            },
        ],
    )
    assert row_identity["status"] == "ok"
    assert row_identity["table_row_count"] == 2
    assert row_identity["min_table_rows"] == 2
    assert row_identity["require_unique_row_identity"] is True
    assert row_identity["checks"]["table_rows_meet_minimum"] is True
    assert row_identity["checks"]["row_identity_columns_recorded_for_uniqueness"] is True
    assert row_identity["checks"]["row_identity_columns_populated"] is True
    assert row_identity["checks"]["row_identity_unique_when_required"] is True
    assert row_identity["checks"]["row_identity_matches_target_region"] is True
    assert row_identity["checks"]["row_case_id_matches_package"] is True
    assert row_identity["checks"]["row_operating_point_id_matches_package"] is True
    assert row_identity["checks"]["row_analysis_type_matches_package"] is True
    assert row_identity["checks"]["row_frequency_hz_matches_package"] is True
    assert row_identity["checks"]["row_target_region_name_matches_package"] is True
    assert row_identity["checks"]["row_target_material_matches_package"] is True
    assert row_identity["checks"]["row_target_region_artifact_id_matches_package"] is True

    row_run_identity = jmag_force_table_metadata_gate(
        {
            **metadata,
            "columns": metadata["columns"] + [
                "MeshId",
                "SolverRunId",
                "ResultRevisionId",
                "SolverSetupArtifactId",
                "MaterialStateArtifactId",
                "MaterialStateDigest",
                "ExcitationSourceArtifactId",
                "CurrentDefinitionMethod",
            ],
            "identity_columns": ["TargetBodyId", "Gap_mm"],
        },
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_mesh_id="jmag_mesh_gap_force_A_v3",
        expected_solver_run_id="jmag_solver_run_20260630_A",
        expected_result_revision_id="jmag_result_revision_force_A_001",
        expected_solver_setup_artifact_id="jmag_solver_setup_force_A_strict_iccg_nr_v1.json",
        expected_material_state_artifact_id="jmag_material_state_force_A_bh_curve_v2.json",
        expected_material_state_digest="sha256:jmag_slot333_material_state_digest_v1",
        expected_excitation_source_artifact_id="jmag_excitation_source_force_A_peak_current_v1.json",
        expected_current_definition_method="peak_phase_current_table",
        table_rows=[
            {
                "Gap_mm": 1.0,
                "Fx_N_per_m": -4.0,
                "Fy_N_per_m": 0.0,
                "CaseId": "case_gap_force_sweep_A",
                "OperatingPointId": "op_gap_sweep_A",
                "TargetBodyId": "moving_core",
                "TargetBodyName": "moving_core_block_label",
                "TargetMaterialName": "magnetic_steel_core",
                "TargetRegionArtifactId": "jmag_region_labels_force_A_v2.json",
                "AnalysisType": "magnetostatic",
                "FrequencyHz": 0.0,
                "MeshId": "jmag_mesh_gap_force_A_v3",
                "SolverRunId": "jmag_solver_run_20260630_A",
                "ResultRevisionId": "jmag_result_revision_force_A_001",
                "SolverSetupArtifactId": "jmag_solver_setup_force_A_strict_iccg_nr_v1.json",
                "MaterialStateArtifactId": "jmag_material_state_force_A_bh_curve_v2.json",
                "MaterialStateDigest": "sha256:jmag_slot333_material_state_digest_v1",
                "ExcitationSourceArtifactId": "jmag_excitation_source_force_A_peak_current_v1.json",
                "CurrentDefinitionMethod": "peak_phase_current_table",
            },
        ],
    )
    assert row_run_identity["status"] == "ok"
    assert row_run_identity["checks"]["row_mesh_id_matches_package"] is True
    assert row_run_identity["checks"]["row_solver_run_id_matches_package"] is True
    assert row_run_identity["checks"]["row_result_revision_id_matches_package"] is True
    assert row_run_identity["checks"]["row_solver_setup_artifact_id_matches_package"] is True
    assert row_run_identity["checks"]["row_material_state_artifact_id_matches_package"] is True
    assert row_run_identity["checks"]["row_material_state_digest_matches_package"] is True
    assert row_run_identity["checks"]["row_excitation_source_artifact_id_matches_package"] is True
    assert row_run_identity["checks"]["row_current_definition_method_matches_package"] is True
    assert row_run_identity["checks"]["row_target_region_name_matches_package"] is True
    assert row_run_identity["checks"]["row_target_material_matches_package"] is True
    assert row_run_identity["checks"]["row_target_region_artifact_id_matches_package"] is True

    stale_solver_material_row = jmag_force_table_metadata_gate(
        {
            **metadata,
            "columns": metadata["columns"] + [
                "MeshId",
                "SolverRunId",
                "ResultRevisionId",
                "SolverSetupArtifactId",
                "MaterialStateArtifactId",
                "MaterialStateDigest",
                "ExcitationSourceArtifactId",
                "CurrentDefinitionMethod",
            ],
            "identity_columns": ["TargetBodyId", "Gap_mm"],
        },
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_solver_setup_artifact_id="jmag_solver_setup_force_A_strict_iccg_nr_v1.json",
        expected_material_state_artifact_id="jmag_material_state_force_A_bh_curve_v2.json",
        expected_material_state_digest="sha256:jmag_slot333_material_state_digest_v1",
        expected_excitation_source_artifact_id="jmag_excitation_source_force_A_peak_current_v1.json",
        expected_current_definition_method="peak_phase_current_table",
        table_rows=[
            {
                "Gap_mm": 1.0,
                "Fx_N_per_m": -4.0,
                "Fy_N_per_m": 0.0,
                "CaseId": "case_gap_force_sweep_A",
                "OperatingPointId": "op_gap_sweep_A",
                "TargetBodyId": "moving_core",
                "AnalysisType": "magnetostatic",
                "FrequencyHz": 0.0,
                "MeshId": "jmag_mesh_gap_force_A_v3",
                "SolverRunId": "jmag_solver_run_20260630_A",
                "ResultRevisionId": "jmag_result_revision_force_A_001",
                "SolverSetupArtifactId": "jmag_solver_setup_force_A_loose_auto_v1.json",
                "MaterialStateArtifactId": "jmag_material_state_force_A_old.json",
                "MaterialStateDigest": "sha256:jmag_slot333_material_state_digest_old",
                "ExcitationSourceArtifactId": "jmag_excitation_source_force_A_old.json",
                "CurrentDefinitionMethod": "rms_phase_current_table",
            },
        ],
    )
    assert stale_solver_material_row["status"] == "needs_attention"
    assert stale_solver_material_row["checks"]["expected_solver_setup_artifact_id_matches"] is True
    assert stale_solver_material_row["checks"]["expected_material_state_artifact_id_matches"] is True
    assert stale_solver_material_row["checks"]["expected_material_state_digest_matches"] is True
    assert stale_solver_material_row["checks"]["row_solver_setup_artifact_id_matches_package"] is False
    assert stale_solver_material_row["checks"]["row_material_state_artifact_id_matches_package"] is False
    assert stale_solver_material_row["checks"]["row_material_state_digest_matches_package"] is False
    assert stale_solver_material_row["checks"]["row_excitation_source_artifact_id_matches_package"] is False
    assert stale_solver_material_row["checks"]["row_current_definition_method_matches_package"] is False
    assert stale_solver_material_row["row_solver_setup_artifact_id_mismatch_rows"] == [0]
    assert stale_solver_material_row["row_material_state_artifact_id_mismatch_rows"] == [0]
    assert stale_solver_material_row["row_material_state_digest_mismatch_rows"] == [0]
    assert stale_solver_material_row["row_excitation_source_artifact_id_mismatch_rows"] == [0]
    assert stale_solver_material_row["row_current_definition_method_mismatch_rows"] == [0]

    stale_target_material_row = jmag_force_table_metadata_gate(
        {
            **metadata,
            "columns": metadata["columns"] + ["MeshId", "SolverRunId", "ResultRevisionId"],
            "identity_columns": ["TargetBodyId", "Gap_mm"],
        },
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_target_region_name="moving_core_block_label",
        expected_target_material="magnetic_steel_core",
        expected_target_region_artifact_id="jmag_region_labels_force_A_v2.json",
        table_rows=[
            {
                "Gap_mm": 1.0,
                "Fx_N_per_m": -4.0,
                "Fy_N_per_m": 0.0,
                "CaseId": "case_gap_force_sweep_A",
                "OperatingPointId": "op_gap_sweep_A",
                "TargetBodyId": "moving_core",
                "TargetBodyName": "moving_core_block_label",
                "TargetMaterialName": "air",
                "TargetRegionArtifactId": "jmag_region_labels_force_A_v2.json",
                "AnalysisType": "magnetostatic",
                "FrequencyHz": 0.0,
                "MeshId": "jmag_mesh_gap_force_A_v3",
                "SolverRunId": "jmag_solver_run_20260630_A",
                "ResultRevisionId": "jmag_result_revision_force_A_001",
            },
        ],
    )
    assert stale_target_material_row["status"] == "needs_attention"
    assert stale_target_material_row["checks"]["expected_target_material_matches"] is True
    assert stale_target_material_row["checks"]["row_target_region_name_matches_package"] is True
    assert stale_target_material_row["checks"]["row_target_material_matches_package"] is False
    assert stale_target_material_row["checks"]["row_target_region_artifact_id_matches_package"] is True
    assert stale_target_material_row["row_target_material_mismatch_rows"] == [0]

    stale_mesh_row = jmag_force_table_metadata_gate(
        {
            **metadata,
            "columns": metadata["columns"] + ["MeshId", "SolverRunId", "ResultRevisionId"],
            "identity_columns": ["TargetBodyId", "Gap_mm"],
        },
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_mesh_id="jmag_mesh_gap_force_A_v3",
        expected_solver_run_id="jmag_solver_run_20260630_A",
        expected_result_revision_id="jmag_result_revision_force_A_001",
        table_rows=[
            {
                "Gap_mm": 1.0,
                "Fx_N_per_m": -4.0,
                "Fy_N_per_m": 0.0,
                "CaseId": "case_gap_force_sweep_A",
                "OperatingPointId": "op_gap_sweep_A",
                "TargetBodyId": "moving_core",
                "AnalysisType": "magnetostatic",
                "FrequencyHz": 0.0,
                "MeshId": "jmag_mesh_gap_force_old",
                "SolverRunId": "jmag_solver_run_20260630_A",
                "ResultRevisionId": "jmag_result_revision_force_A_001",
            },
        ],
    )
    assert stale_mesh_row["status"] == "needs_attention"
    assert stale_mesh_row["checks"]["expected_mesh_id_matches"] is True
    assert stale_mesh_row["checks"]["row_mesh_id_matches_package"] is False
    assert stale_mesh_row["checks"]["row_solver_run_id_matches_package"] is True
    assert stale_mesh_row["checks"]["row_result_revision_id_matches_package"] is True
    assert stale_mesh_row["row_mesh_id_mismatch_rows"] == [0]

    stale_case_row = jmag_force_table_metadata_gate(
        {**metadata, "identity_columns": ["TargetBodyId", "Gap_mm"]},
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_case_id="case_gap_force_sweep_A",
        expected_operating_point_id="op_gap_sweep_A",
        expected_export_artifact_id="jmag_force_export_artifact_A",
        expected_result_set_id="resultset_force_20260630_A",
        expected_target_region_id="moving_core",
        min_table_rows=1,
        require_unique_row_identity=True,
        table_rows=[
            {
                "Gap_mm": 1.0,
                "Fx_N_per_m": -4.0,
                "Fy_N_per_m": 0.0,
                "CaseId": "case_gap_force_sweep_old",
                "OperatingPointId": "op_gap_sweep_A",
                "TargetBodyId": "moving_core",
                "AnalysisType": "magnetostatic",
                "FrequencyHz": 0.0,
            },
        ],
    )
    assert stale_case_row["status"] == "needs_attention"
    assert stale_case_row["checks"]["expected_case_id_matches"] is True
    assert stale_case_row["checks"]["row_case_id_matches_package"] is False
    assert stale_case_row["checks"]["row_operating_point_id_matches_package"] is True
    assert stale_case_row["row_case_id_mismatch_rows"] == [0]

    stale_analysis_metadata = jmag_force_table_metadata_gate(
        {**metadata, "analysis_type": "transient"},
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_analysis_type="magnetostatic",
    )
    assert stale_analysis_metadata["status"] == "needs_attention"
    assert stale_analysis_metadata["checks"]["expected_analysis_type_matches"] is False

    stale_frequency_metadata = jmag_force_table_metadata_gate(
        {**metadata, "frequency_hz": 50.0},
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_frequency_hz=0.0,
    )
    assert stale_frequency_metadata["status"] == "needs_attention"
    assert stale_frequency_metadata["checks"]["expected_frequency_hz_matches"] is False

    stale_analysis_row = jmag_force_table_metadata_gate(
        {**metadata, "identity_columns": ["TargetBodyId", "Gap_mm"]},
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_analysis_type="magnetostatic",
        expected_frequency_hz=0.0,
        table_rows=[
            {
                "Gap_mm": 1.0,
                "Fx_N_per_m": -4.0,
                "Fy_N_per_m": 0.0,
                "CaseId": "case_gap_force_sweep_A",
                "OperatingPointId": "op_gap_sweep_A",
                "TargetBodyId": "moving_core",
                "AnalysisType": "frequency_response",
                "FrequencyHz": 50.0,
            },
        ],
    )
    assert stale_analysis_row["status"] == "needs_attention"
    assert stale_analysis_row["checks"]["expected_analysis_type_matches"] is True
    assert stale_analysis_row["checks"]["expected_frequency_hz_matches"] is True
    assert stale_analysis_row["checks"]["row_analysis_type_matches_package"] is False
    assert stale_analysis_row["checks"]["row_frequency_hz_matches_package"] is False
    assert stale_analysis_row["row_analysis_type_mismatch_rows"] == [0]
    assert stale_analysis_row["row_frequency_hz_mismatch_rows"] == [0]

    empty_force_rows = jmag_force_table_metadata_gate(
        metadata,
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_export_artifact_id="jmag_force_export_artifact_A",
        expected_result_set_id="resultset_force_20260630_A",
        expected_target_region_id="moving_core",
        min_table_rows=1,
        table_rows=[],
    )
    assert empty_force_rows["status"] == "needs_attention"
    assert empty_force_rows["table_row_count"] == 0
    assert empty_force_rows["checks"]["table_rows_meet_minimum"] is False

    duplicate_row_identity = jmag_force_table_metadata_gate(
        {**metadata, "identity_columns": ["TargetBodyId", "Gap_mm"]},
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_export_artifact_id="jmag_force_export_artifact_A",
        expected_result_set_id="resultset_force_20260630_A",
        expected_target_region_id="moving_core",
        min_table_rows=2,
        require_unique_row_identity=True,
        table_rows=[
            {
                "Gap_mm": 1.0,
                "Fx_N_per_m": -4.0,
                "Fy_N_per_m": 0.0,
                "CaseId": "case_gap_force_sweep_A",
                "OperatingPointId": "op_gap_sweep_A",
                "TargetBodyId": "moving_core",
            },
            {
                "Gap_mm": 1.0,
                "Fx_N_per_m": -4.1,
                "Fy_N_per_m": 0.0,
                "CaseId": "case_gap_force_sweep_A",
                "OperatingPointId": "op_gap_sweep_A",
                "TargetBodyId": "moving_core",
            },
        ],
    )
    assert duplicate_row_identity["status"] == "needs_attention"
    assert duplicate_row_identity["checks"]["table_rows_meet_minimum"] is True
    assert duplicate_row_identity["checks"]["row_identity_unique_when_required"] is False
    assert duplicate_row_identity["duplicate_row_identity_values"] == [
        {"first_row": 0, "row": 1, "values": ["moving_core", "1.0"]}
    ]

    wrong_row_identity = jmag_force_table_metadata_gate(
        metadata,
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_export_artifact_id="jmag_force_export_artifact_A",
        expected_result_set_id="resultset_force_20260630_A",
        expected_target_region_id="moving_core",
        table_rows=[
            {
                "Gap_mm": 1.0,
                "Fx_N_per_m": -4.0,
                "Fy_N_per_m": 0.0,
                "CaseId": "case_gap_force_sweep_A",
                "OperatingPointId": "op_gap_sweep_A",
                "TargetBodyId": "fixed_core",
            },
        ],
    )
    assert wrong_row_identity["status"] == "needs_attention"
    assert wrong_row_identity["checks"]["row_identity_columns_populated"] is True
    assert wrong_row_identity["checks"]["row_identity_matches_target_region"] is False
    assert wrong_row_identity["row_identity_mismatch_rows"] == [0]

    wrong_force_dimension = jmag_force_table_metadata_gate(
        {**metadata, "force_unit": "N", "quantity_dimension": "2d_per_length"},
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_export_artifact_id="jmag_force_export_artifact_A",
        expected_result_set_id="resultset_force_20260630_A",
        expected_target_region_id="moving_core",
    )
    assert wrong_force_dimension["status"] == "needs_attention"
    assert wrong_force_dimension["checks"]["force_unit_valid"] is True
    assert wrong_force_dimension["checks"]["quantity_dimension_valid"] is True
    assert wrong_force_dimension["checks"]["force_unit_matches_quantity_dimension_when_present"] is False

    missing_projection_axis = dict(metadata)
    missing_projection_axis.pop("projection_axis")
    missing_projection = jmag_force_table_metadata_gate(
        missing_projection_axis,
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_export_artifact_id="jmag_force_export_artifact_A",
        expected_result_set_id="resultset_force_20260630_A",
        expected_target_region_id="moving_core",
    )
    assert missing_projection["status"] == "needs_attention"
    assert missing_projection["checks"]["force_projection_axis_recorded"] is False
    assert missing_projection["checks"]["force_columns_match_component_frame"] is True

    vague_projection_axis = jmag_force_table_metadata_gate(
        {**metadata, "projection_axis": "screen"},
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_export_artifact_id="jmag_force_export_artifact_A",
        expected_result_set_id="resultset_force_20260630_A",
        expected_target_region_id="moving_core",
    )
    assert vague_projection_axis["status"] == "needs_attention"
    assert vague_projection_axis["checks"]["force_projection_axis_recorded"] is True
    assert vague_projection_axis["checks"]["force_projection_axis_descriptive"] is False

    wrong_frame_columns = jmag_force_table_metadata_gate(
        {**metadata, "component_frame": "cylindrical_rt"},
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_export_artifact_id="jmag_force_export_artifact_A",
        expected_result_set_id="resultset_force_20260630_A",
        expected_target_region_id="moving_core",
    )
    assert wrong_frame_columns["status"] == "needs_attention"
    assert wrong_frame_columns["checks"]["component_frame_valid"] is True
    assert wrong_frame_columns["checks"]["force_columns_match_component_frame"] is False
    assert wrong_frame_columns["expected_force_axes_for_frame"] == ["radial", "tangential"]

    wrong_identity = jmag_force_table_metadata_gate(
        metadata,
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_export_artifact_id="jmag_force_export_artifact_A",
        expected_result_set_id="resultset_force_20260630_A",
        expected_target_region_id="fixed_core",
    )
    assert wrong_identity["status"] == "needs_attention"
    assert wrong_identity["checks"]["expected_result_set_id_matches"] is True
    assert wrong_identity["checks"]["expected_target_region_id_matches"] is False

    bad = jmag_force_table_metadata_gate(
        {
            "columns": ["Gap", "Fx"],
            "position_columns": ["Gap_mm"],
            "force_columns": ["Fx_N_per_m", "Fy_N_per_m"],
            "identity_columns": ["TargetBodyId"],
            "force_unit": "Nmm",
            "position_unit": "inch",
            "component_frame": "screen",
            "source_tool": "FEMM",
            "export_artifact_id": "old_export",
            "result_set_id": "old",
            "symmetry_factor": 0,
            "force_sign_convention": "unknown",
            "force_kind": "screenshot",
            "quantity_dimension": "screen_force",
        },
        required_columns=["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"],
        expected_export_artifact_id="jmag_force_export_artifact_A",
        expected_result_set_id="resultset_force_20260630_A",
        expected_target_region_id="moving_core",
    )

    assert bad["status"] == "needs_attention"
    assert bad["missing_required_columns"] == ["Gap_mm", "Fx_N_per_m", "Fy_N_per_m"]
    assert bad["missing_force_columns"] == ["Fx_N_per_m", "Fy_N_per_m"]
    assert bad["missing_position_columns"] == ["Gap_mm"]
    assert bad["checks"]["force_unit_valid"] is False
    assert bad["checks"]["position_unit_valid"] is False
    assert bad["checks"]["component_frame_valid"] is False
    assert bad["checks"]["force_columns_match_component_frame"] is True
    assert bad["checks"]["source_tool_is_jmag"] is False
    assert bad["checks"]["expected_export_artifact_id_matches"] is False
    assert bad["checks"]["symmetry_factor_positive"] is False
    assert bad["checks"]["force_sign_convention_valid"] is False
    assert bad["checks"]["force_kind_valid"] is False
    assert bad["checks"]["quantity_dimension_valid"] is False
    assert bad["checks"]["expected_result_set_id_matches"] is False
    assert bad["checks"]["target_region_id_recorded"] is False
    assert bad["checks"]["identity_columns_present"] is False


def test_jmag_airgap_flux_sample_metadata_gate_binds_probe_output_identity():
    rows = []
    for angle, br, bt in [(0.0, 0.72, 0.11), (5.0, 0.71, 0.12), (10.0, 0.70, 0.13)]:
        rows.append({
            "source_tool": "JMAG-Designer",
            "result_set_id": "jmag_slot349_resultset_airgap_A",
            "export_artifact_id": "jmag_slot349_airgap_line_probe_export_v1",
            "field_probe_id": "jmag_slot349_mid_airgap_probe_v1",
            "field_probe_method": "jmag_airgap_line_probe_export",
            "field_probe_output_artifact_id": "jmag_slot349_airgap_flux_table_v1.csv",
            "field_probe_output_digest": "sha256:jmag-slot349-airgap-flux-table-v1",
            "field_probe_output_path": "artifacts/motor/jmag_slot349_airgap_flux_table_v1.csv",
            "sample_grid_id": "jmag_slot357_airgap_angle_grid_v1",
            "sample_grid_digest": "sha256:jmag-slot357-airgap-angle-grid-v1",
            "sample_count": 3,
            "RotorAngle_deg": angle,
            "angle_unit": "deg",
            "angle_basis": "mechanical",
            "component_frame": "cylindrical_rt",
            "Br_T": br,
            "Bt_T": bt,
            "radius_m": 0.0505,
            "axial_length_m": 0.1,
            "symmetry_factor": 6,
            "torque_sign_convention": "positive_motoring",
        })

    gate = jmag_airgap_flux_sample_metadata_gate(
        rows,
        expected_result_set_id="jmag_slot349_resultset_airgap_A",
        expected_export_artifact_id="jmag_slot349_airgap_line_probe_export_v1",
        expected_field_probe_id="jmag_slot349_mid_airgap_probe_v1",
        expected_field_probe_method="jmag_airgap_line_probe_export",
        expected_field_probe_output_artifact_id="jmag_slot349_airgap_flux_table_v1.csv",
        expected_field_probe_output_digest="sha256:jmag-slot349-airgap-flux-table-v1",
        expected_sample_grid_id="jmag_slot357_airgap_angle_grid_v1",
        expected_sample_grid_digest="sha256:jmag-slot357-airgap-angle-grid-v1",
        expected_sample_count=3,
        expected_angle_unit="deg",
        expected_angle_basis="mechanical",
        expected_component_frame="cylindrical_rt",
        expected_torque_sign_convention="positive_motoring",
        require_field_probe_output_artifact=True,
    )

    assert gate["policy"] == "jmag_airgap_flux_sample_metadata_gate"
    assert gate["status"] == "ok"
    assert gate["checks"]["source_tool_is_jmag"] is True
    assert gate["checks"]["br_t_recorded_and_finite"] is True
    assert gate["checks"]["bt_t_recorded_and_finite"] is True
    assert gate["checks"]["expected_result_set_id_matches"] is True
    assert gate["checks"]["expected_export_artifact_id_matches"] is True
    assert gate["checks"]["expected_field_probe_id_matches"] is True
    assert gate["checks"]["expected_field_probe_method_matches"] is True
    assert gate["checks"]["expected_field_probe_output_artifact_id_matches"] is True
    assert gate["checks"]["expected_field_probe_output_digest_matches"] is True
    assert gate["checks"]["field_probe_output_path_recorded"] is True
    assert gate["checks"]["expected_sample_grid_id_matches"] is True
    assert gate["checks"]["expected_sample_grid_digest_matches"] is True
    assert gate["checks"]["expected_sample_count_matches"] is True
    assert gate["sample_grid_ids"] == ["jmag_slot357_airgap_angle_grid_v1"]
    assert gate["field_probe_ids"] == ["jmag_slot349_mid_airgap_probe_v1"]
    assert gate["field_probe_methods"] == ["jmag_airgap_line_probe_export"]

    stale_result = jmag_airgap_flux_sample_metadata_gate(
        [{**row, "result_set_id": "jmag_slot341_old_resultset"} for row in rows],
        expected_result_set_id="jmag_slot349_resultset_airgap_A",
    )
    assert stale_result["status"] == "needs_attention"
    assert stale_result["checks"]["expected_result_set_id_matches"] is False

    stale_export = jmag_airgap_flux_sample_metadata_gate(
        [{**row, "export_artifact_id": "jmag_slot333_old_airgap_export"} for row in rows],
        expected_export_artifact_id="jmag_slot349_airgap_line_probe_export_v1",
    )
    assert stale_export["status"] == "needs_attention"
    assert stale_export["checks"]["expected_export_artifact_id_matches"] is False

    stale_probe = jmag_airgap_flux_sample_metadata_gate(
        [{**row, "field_probe_id": "jmag_slot326_old_probe"} for row in rows],
        expected_field_probe_id="jmag_slot349_mid_airgap_probe_v1",
    )
    assert stale_probe["status"] == "needs_attention"
    assert stale_probe["checks"]["expected_field_probe_id_matches"] is False

    wrong_method = jmag_airgap_flux_sample_metadata_gate(
        [{**row, "field_probe_method": "point_b_probe"} for row in rows],
        expected_field_probe_method="jmag_airgap_line_probe_export",
    )
    assert wrong_method["status"] == "needs_attention"
    assert wrong_method["checks"]["expected_field_probe_method_matches"] is False

    stale_output = jmag_airgap_flux_sample_metadata_gate(
        [{**row, "field_probe_output_artifact_id": "jmag_slot293_old_airgap.csv"} for row in rows],
        expected_field_probe_output_artifact_id="jmag_slot349_airgap_flux_table_v1.csv",
        expected_field_probe_output_digest="sha256:jmag-slot349-airgap-flux-table-v1",
        require_field_probe_output_artifact=True,
    )
    assert stale_output["status"] == "needs_attention"
    assert stale_output["checks"]["expected_field_probe_output_artifact_id_matches"] is False
    assert stale_output["checks"]["expected_field_probe_output_digest_matches"] is True

    stale_grid = jmag_airgap_flux_sample_metadata_gate(
        [{**row, "sample_grid_id": "jmag_slot349_old_angle_grid"} for row in rows],
        expected_sample_grid_id="jmag_slot357_airgap_angle_grid_v1",
        expected_sample_grid_digest="sha256:jmag-slot357-airgap-angle-grid-v1",
        expected_sample_count=3,
    )
    assert stale_grid["status"] == "needs_attention"
    assert stale_grid["checks"]["expected_sample_grid_id_matches"] is False
    assert stale_grid["checks"]["expected_sample_grid_digest_matches"] is True

    wrong_sample_count = jmag_airgap_flux_sample_metadata_gate(
        [{**row, "sample_count": 4} for row in rows],
        expected_sample_count=3,
    )
    assert wrong_sample_count["status"] == "needs_attention"
    assert wrong_sample_count["checks"]["expected_sample_count_matches"] is False

    wrong_frame = jmag_airgap_flux_sample_metadata_gate(
        [{**row, "component_frame": "cartesian_xy"} for row in rows],
        expected_component_frame="cylindrical_rt",
    )
    assert wrong_frame["status"] == "needs_attention"
    assert wrong_frame["checks"]["expected_component_frame_matches"] is False

    missing_bt = jmag_airgap_flux_sample_metadata_gate(
        [{k: v for k, v in row.items() if k != "Bt_T"} for row in rows],
    )
    assert missing_bt["status"] == "needs_attention"
    assert missing_bt["checks"]["bt_t_recorded_and_finite"] is False


def test_jmag_airgap_torque_integration_package_gate_binds_field_and_grid_identity():
    rows = []
    for angle, br, bt in [(0.0, 0.72, 0.11), (5.0, 0.71, 0.12), (10.0, 0.70, 0.13)]:
        rows.append({
            "source_tool": "JMAG-Designer",
            "result_set_id": "jmag_slot365_resultset_airgap_A",
            "export_artifact_id": "jmag_slot365_airgap_line_probe_export_v1",
            "field_probe_id": "jmag_slot365_mid_airgap_probe_v1",
            "field_probe_method": "jmag_airgap_line_probe_export",
            "field_probe_output_artifact_id": "jmag_slot365_airgap_flux_table_v1.csv",
            "field_probe_output_digest": "sha256:jmag-slot365-airgap-flux-table-v1",
            "field_probe_output_path": "artifacts/motor/jmag_slot365_airgap_flux_table_v1.csv",
            "sample_grid_id": "jmag_slot365_airgap_angle_grid_v1",
            "sample_grid_digest": "sha256:jmag-slot365-airgap-angle-grid-v1",
            "sample_count": 3,
            "RotorAngle_deg": angle,
            "angle_unit": "deg",
            "angle_basis": "mechanical",
            "component_frame": "cylindrical_rt",
            "Br_T": br,
            "Bt_T": bt,
            "radius_m": 0.0505,
            "axial_length_m": 0.1,
            "symmetry_factor": 6,
            "torque_sign_convention": "positive_motoring",
        })
    sample_gate = jmag_airgap_flux_sample_metadata_gate(
        rows,
        expected_result_set_id="jmag_slot365_resultset_airgap_A",
        expected_export_artifact_id="jmag_slot365_airgap_line_probe_export_v1",
        expected_field_probe_id="jmag_slot365_mid_airgap_probe_v1",
        expected_field_probe_method="jmag_airgap_line_probe_export",
        expected_field_probe_output_artifact_id="jmag_slot365_airgap_flux_table_v1.csv",
        expected_field_probe_output_digest="sha256:jmag-slot365-airgap-flux-table-v1",
        expected_sample_grid_id="jmag_slot365_airgap_angle_grid_v1",
        expected_sample_grid_digest="sha256:jmag-slot365-airgap-angle-grid-v1",
        expected_sample_count=3,
        expected_angle_basis="mechanical",
        expected_torque_sign_convention="positive_motoring",
        require_field_probe_output_artifact=True,
    )
    assert sample_gate["status"] == "ok"

    package = {
        "input_field_table_artifact_id": "jmag_slot365_airgap_flux_table_v1.csv",
        "input_field_table_digest": "sha256:jmag-slot365-airgap-flux-table-v1",
        "input_field_table_path": "artifacts/motor/jmag_slot365_airgap_flux_table_v1.csv",
        "model_input_artifact_id": "jmag_slot379_airgap_motor_project_v1.jproj",
        "model_input_digest": "sha256:jmag-slot379-airgap-motor-project-v1",
        "model_input_path": "artifacts/motor/jmag_slot379_airgap_motor_project_v1.jproj",
        "export_recipe_artifact_id": "jmag_slot386_airgap_export_recipe_v1.py",
        "export_recipe_digest": "sha256:jmag-slot386-airgap-export-recipe-v1",
        "export_recipe_path": "artifacts/motor/jmag_slot386_airgap_export_recipe.py",
        "parameter_set_artifact_id": "jmag_slot393_airgap_torque_parameter_set_v1.json",
        "parameter_set_digest": "sha256:jmag-slot393-airgap-torque-parameter-set-v1",
        "parameter_set_path": "artifacts/motor/jmag_slot393_airgap_torque_parameter_set.json",
        "objective_observable_id": "jmag_slot393_airgap_torque_objective_v1",
        "objective_observable_family": "airgap_torque_maxwell_shear_objective",
        "sample_grid_id": "jmag_slot365_airgap_angle_grid_v1",
        "sample_grid_digest": "sha256:jmag-slot365-airgap-angle-grid-v1",
        "sample_count": 3,
        "integration_method": "maxwell_shear_from_br_bt_samples",
        "integration_policy": "trapezoid_periodic_full_sector_expand_symmetry",
        "component_frame": "cylindrical_rt",
        "torque_sign_convention": "positive_motoring",
        "torque_output_artifact_id": "jmag_slot365_airgap_torque_integration_v1.json",
        "torque_output_digest": "sha256:jmag-slot365-airgap-torque-integration-v1",
        "torque_output_path": "artifacts/motor/jmag_slot365_airgap_torque_integration_v1.json",
        "torque_output_schema_id": "jmag_airgap_torque_table_v1",
        "torque_output_columns": [
            "RotorAngle_deg",
            "Br_T",
            "Bt_T",
            "torque_density_N_per_m2",
            "torque_Nm",
        ],
        "torque_output_units": {
            "RotorAngle_deg": "deg",
            "Br_T": "T",
            "Bt_T": "T",
            "torque_density_N_per_m2": "N/m^2",
            "torque_Nm": "N*m",
        },
        "torque_convention_schema_id": "jmag_airgap_torque_convention_v1",
        "torque_component_basis_schema_id": "jmag_airgap_cylindrical_rt_component_basis_v1",
        "torque_postprocess_row_convention_schema_id": "jmag_airgap_torque_row_convention_v1",
        "torque_Nm": 0.04125,
    }
    gate = jmag_airgap_torque_integration_package_gate(
        sample_gate,
        package,
        expected_input_field_table_artifact_id="jmag_slot365_airgap_flux_table_v1.csv",
        expected_input_field_table_digest="sha256:jmag-slot365-airgap-flux-table-v1",
        expected_model_input_artifact_id="jmag_slot379_airgap_motor_project_v1.jproj",
        expected_model_input_digest="sha256:jmag-slot379-airgap-motor-project-v1",
        expected_model_input_path="artifacts/motor/jmag_slot379_airgap_motor_project_v1.jproj",
        expected_export_recipe_artifact_id="jmag_slot386_airgap_export_recipe_v1.py",
        expected_export_recipe_digest="sha256:jmag-slot386-airgap-export-recipe-v1",
        expected_export_recipe_path="artifacts/motor/jmag_slot386_airgap_export_recipe.py",
        expected_parameter_set_artifact_id="jmag_slot393_airgap_torque_parameter_set_v1.json",
        expected_parameter_set_digest="sha256:jmag-slot393-airgap-torque-parameter-set-v1",
        expected_parameter_set_path="artifacts/motor/jmag_slot393_airgap_torque_parameter_set.json",
        expected_objective_observable_id="jmag_slot393_airgap_torque_objective_v1",
        expected_objective_observable_family="airgap_torque_maxwell_shear_objective",
        expected_sample_grid_id="jmag_slot365_airgap_angle_grid_v1",
        expected_sample_grid_digest="sha256:jmag-slot365-airgap-angle-grid-v1",
        expected_torque_output_artifact_id="jmag_slot365_airgap_torque_integration_v1.json",
        expected_torque_output_digest="sha256:jmag-slot365-airgap-torque-integration-v1",
        expected_torque_output_schema_id="jmag_airgap_torque_table_v1",
        expected_torque_output_columns=[
            "RotorAngle_deg",
            "Br_T",
            "Bt_T",
            "torque_density_N_per_m2",
            "torque_Nm",
        ],
        expected_torque_output_units={
            "RotorAngle_deg": "deg",
            "Br_T": "T",
            "Bt_T": "T",
            "torque_density_N_per_m2": "N/m^2",
            "torque_Nm": "N*m",
        },
        expected_torque_convention_schema_id="jmag_airgap_torque_convention_v1",
        expected_torque_component_basis_schema_id="jmag_airgap_cylindrical_rt_component_basis_v1",
        expected_torque_postprocess_row_convention_schema_id="jmag_airgap_torque_row_convention_v1",
        expected_torque_nm=0.04125,
        torque_abs_tol=1.0e-14,
        require_export_recipe_artifact=True,
        require_parameter_set_artifact=True,
        require_torque_output_artifact=True,
        require_torque_output_schema=True,
        require_torque_convention_schema=True,
        require_torque_component_basis_schema=True,
        require_torque_postprocess_row_convention_schema=True,
    )

    assert gate["policy"] == "jmag_airgap_torque_integration_package_gate"
    assert gate["status"] == "ok"
    assert gate["checks"]["sample_metadata_gate_ok"] is True
    assert gate["checks"]["input_field_table_artifact_matches_sample_gate"] is True
    assert gate["checks"]["input_field_table_digest_matches_sample_gate"] is True
    assert gate["checks"]["model_input_artifact_id_recorded"] is True
    assert gate["checks"]["model_input_digest_recorded"] is True
    assert gate["checks"]["model_input_path_recorded"] is True
    assert gate["checks"]["expected_model_input_artifact_id_matches"] is True
    assert gate["checks"]["expected_model_input_digest_matches"] is True
    assert gate["checks"]["expected_model_input_path_matches"] is True
    assert gate["checks"]["export_recipe_artifact_id_recorded"] is True
    assert gate["checks"]["export_recipe_digest_recorded"] is True
    assert gate["checks"]["export_recipe_path_recorded"] is True
    assert gate["checks"]["expected_export_recipe_artifact_id_matches"] is True
    assert gate["checks"]["expected_export_recipe_digest_matches"] is True
    assert gate["checks"]["expected_export_recipe_path_matches"] is True
    assert gate["export_recipe_artifact_required"] is True
    assert gate["export_recipe_artifact_id"] == "jmag_slot386_airgap_export_recipe_v1.py"
    assert gate["parameter_set_artifact_required"] is True
    assert gate["checks"]["parameter_set_artifact_id_recorded"] is True
    assert gate["checks"]["parameter_set_digest_recorded"] is True
    assert gate["checks"]["parameter_set_path_recorded"] is True
    assert gate["checks"]["expected_parameter_set_artifact_id_matches"] is True
    assert gate["checks"]["expected_parameter_set_digest_matches"] is True
    assert gate["checks"]["expected_parameter_set_path_matches"] is True
    assert gate["checks"]["objective_observable_id_recorded"] is True
    assert gate["checks"]["expected_objective_observable_id_matches"] is True
    assert gate["checks"]["objective_observable_family_recorded"] is True
    assert gate["checks"]["expected_objective_observable_family_matches"] is True
    assert gate["parameter_set_artifact_id"] == "jmag_slot393_airgap_torque_parameter_set_v1.json"
    assert gate["objective_observable_family"] == "airgap_torque_maxwell_shear_objective"
    assert gate["checks"]["sample_grid_id_matches_sample_gate"] is True
    assert gate["checks"]["sample_grid_digest_matches_sample_gate"] is True
    assert gate["checks"]["sample_count_matches_sample_gate"] is True
    assert gate["checks"]["expected_integration_method_matches"] is True
    assert gate["torque_output_schema_required"] is True
    assert gate["checks"]["torque_output_schema_id_recorded"] is True
    assert gate["checks"]["expected_torque_output_schema_id_matches"] is True
    assert gate["checks"]["expected_torque_output_columns_match"] is True
    assert gate["checks"]["expected_torque_output_units_match"] is True
    assert gate["torque_output_schema_id"] == "jmag_airgap_torque_table_v1"
    assert gate["torque_convention_schema_required"] is True
    assert gate["checks"]["torque_convention_schema_id_recorded"] is True
    assert gate["checks"]["expected_torque_convention_schema_id_matches"] is True
    assert gate["torque_convention_schema_id"] == "jmag_airgap_torque_convention_v1"
    assert gate["torque_component_basis_schema_required"] is True
    assert gate["checks"]["torque_component_basis_schema_id_recorded"] is True
    assert gate["checks"]["expected_torque_component_basis_schema_id_matches"] is True
    assert gate["torque_component_basis_schema_id"] == "jmag_airgap_cylindrical_rt_component_basis_v1"
    assert gate["expected_torque_component_basis_schema_id"] == "jmag_airgap_cylindrical_rt_component_basis_v1"
    assert gate["torque_postprocess_row_convention_schema_required"] is True
    assert gate["checks"]["torque_postprocess_row_convention_schema_id_recorded"] is True
    assert gate["checks"]["expected_torque_postprocess_row_convention_schema_id_matches"] is True
    assert gate["torque_postprocess_row_convention_schema_id"] == "jmag_airgap_torque_row_convention_v1"
    assert gate["checks"]["expected_torque_nm_matches"] is True

    execution_package = {
        **package,
        "torque_output_artifact_id": "jmag_slot372_airgap_torque_integration_v1.json",
        "torque_output_digest": "sha256:jmag-slot372-airgap-torque-integration-v1",
        "torque_output_path": "artifacts/motor/jmag_slot372_airgap_torque_integration_v1.json",
        "created_at_utc": "2026-07-01T11:12:20Z",
        "run_timestamp_utc": "2026-07-01T11:12:00Z",
        "solver_version": "JMAG-Designer 24.2.1",
        "radia_mcp_version": "1.4.3",
        "run_duration_s": 9.4,
        "timing_breakdown_s": {
            "solve_s": 5.2,
            "export_s": 2.1,
            "integrate_s": 1.4,
            "write_json_s": 0.3,
        },
    }
    execution_gate = jmag_airgap_torque_integration_package_gate(
        sample_gate,
        execution_package,
        expected_input_field_table_artifact_id="jmag_slot365_airgap_flux_table_v1.csv",
        expected_input_field_table_digest="sha256:jmag-slot365-airgap-flux-table-v1",
        expected_sample_grid_id="jmag_slot365_airgap_angle_grid_v1",
        expected_sample_grid_digest="sha256:jmag-slot365-airgap-angle-grid-v1",
        expected_torque_output_artifact_id="jmag_slot372_airgap_torque_integration_v1.json",
        expected_torque_output_digest="sha256:jmag-slot372-airgap-torque-integration-v1",
        expected_torque_nm=0.04125,
        expected_created_at_utc="2026-07-01T11:12:20Z",
        expected_run_timestamp_utc="2026-07-01T11:12:00Z",
        expected_solver_version="JMAG-Designer 24.2.1",
        expected_radia_mcp_version="1.4.3",
        max_created_run_skew_s=60,
        require_torque_output_artifact=True,
        require_execution_metadata=True,
        require_timing_breakdown=True,
        min_timing_sections=4,
    )
    assert execution_gate["status"] == "ok"
    assert execution_gate["execution_metadata_required"] is True
    assert execution_gate["timing_breakdown_required"] is True
    assert execution_gate["checks"]["created_at_utc_recorded"] is True
    assert execution_gate["checks"]["created_at_utc_parseable"] is True
    assert execution_gate["checks"]["expected_created_at_utc_matches"] is True
    assert execution_gate["checks"]["run_timestamp_utc_recorded"] is True
    assert execution_gate["checks"]["run_timestamp_utc_parseable"] is True
    assert execution_gate["checks"]["created_run_timestamp_skew_within_limit"] is True
    assert execution_gate["created_run_timestamp_skew_s"] == pytest.approx(20.0)
    assert execution_gate["checks"]["expected_solver_version_matches"] is True
    assert execution_gate["checks"]["expected_radia_mcp_version_matches"] is True
    assert execution_gate["checks"]["timing_breakdown_has_required_sections"] is True
    assert execution_gate["checks"]["timing_breakdown_top_sections_descending"] is True
    assert execution_gate["timing_total_s"] == pytest.approx(9.0)

    model_input_package = {
        **execution_package,
        "torque_output_artifact_id": "jmag_slot379_airgap_torque_integration_v1.json",
        "torque_output_digest": "sha256:jmag-slot379-airgap-torque-integration-v1",
        "torque_output_path": "artifacts/motor/jmag_slot379_airgap_torque_integration_v1.json",
    }
    model_input_gate = jmag_airgap_torque_integration_package_gate(
        sample_gate,
        model_input_package,
        expected_input_field_table_artifact_id="jmag_slot365_airgap_flux_table_v1.csv",
        expected_input_field_table_digest="sha256:jmag-slot365-airgap-flux-table-v1",
        expected_model_input_artifact_id="jmag_slot379_airgap_motor_project_v1.jproj",
        expected_model_input_digest="sha256:jmag-slot379-airgap-motor-project-v1",
        expected_model_input_path="artifacts/motor/jmag_slot379_airgap_motor_project_v1.jproj",
        expected_sample_grid_id="jmag_slot365_airgap_angle_grid_v1",
        expected_sample_grid_digest="sha256:jmag-slot365-airgap-angle-grid-v1",
        expected_torque_output_artifact_id="jmag_slot379_airgap_torque_integration_v1.json",
        expected_torque_output_digest="sha256:jmag-slot379-airgap-torque-integration-v1",
        expected_torque_nm=0.04125,
        require_model_input_artifact=True,
        require_torque_output_artifact=True,
        require_execution_metadata=True,
        require_timing_breakdown=True,
    )
    assert model_input_gate["status"] == "ok"
    assert model_input_gate["model_input_artifact_required"] is True
    assert model_input_gate["model_input_artifact_id"] == "jmag_slot379_airgap_motor_project_v1.jproj"
    assert model_input_gate["checks"]["model_input_artifact_id_recorded"] is True
    assert model_input_gate["checks"]["model_input_digest_recorded"] is True
    assert model_input_gate["checks"]["model_input_path_recorded"] is True
    assert model_input_gate["checks"]["expected_model_input_digest_matches"] is True

    stale_project_digest = jmag_airgap_torque_integration_package_gate(
        sample_gate,
        {**model_input_package, "model_input_digest": "sha256:jmag-slot365-old-project"},
        expected_model_input_artifact_id="jmag_slot379_airgap_motor_project_v1.jproj",
        expected_model_input_digest="sha256:jmag-slot379-airgap-motor-project-v1",
        expected_model_input_path="artifacts/motor/jmag_slot379_airgap_motor_project_v1.jproj",
        require_model_input_artifact=True,
    )
    assert stale_project_digest["status"] == "needs_attention"
    assert stale_project_digest["checks"]["expected_model_input_artifact_id_matches"] is True
    assert stale_project_digest["checks"]["expected_model_input_digest_matches"] is False

    missing_project_path = jmag_airgap_torque_integration_package_gate(
        sample_gate,
        {key: value for key, value in model_input_package.items() if key != "model_input_path"},
        require_model_input_artifact=True,
    )
    assert missing_project_path["status"] == "needs_attention"
    assert missing_project_path["checks"]["model_input_artifact_id_recorded"] is True
    assert missing_project_path["checks"]["model_input_digest_recorded"] is True
    assert missing_project_path["checks"]["model_input_path_recorded"] is False

    stale_export_recipe_digest = jmag_airgap_torque_integration_package_gate(
        sample_gate,
        {**model_input_package, "export_recipe_digest": "sha256:jmag-slot365-old-export-recipe"},
        expected_export_recipe_artifact_id="jmag_slot386_airgap_export_recipe_v1.py",
        expected_export_recipe_digest="sha256:jmag-slot386-airgap-export-recipe-v1",
        expected_export_recipe_path="artifacts/motor/jmag_slot386_airgap_export_recipe.py",
        require_export_recipe_artifact=True,
    )
    assert stale_export_recipe_digest["status"] == "needs_attention"
    assert stale_export_recipe_digest["checks"]["expected_export_recipe_artifact_id_matches"] is True
    assert stale_export_recipe_digest["checks"]["expected_export_recipe_digest_matches"] is False

    missing_export_recipe_path = jmag_airgap_torque_integration_package_gate(
        sample_gate,
        {key: value for key, value in model_input_package.items() if key != "export_recipe_path"},
        require_export_recipe_artifact=True,
    )
    assert missing_export_recipe_path["status"] == "needs_attention"
    assert missing_export_recipe_path["checks"]["export_recipe_artifact_id_recorded"] is True
    assert missing_export_recipe_path["checks"]["export_recipe_digest_recorded"] is True
    assert missing_export_recipe_path["checks"]["export_recipe_path_recorded"] is False

    stale_parameter_set_digest = jmag_airgap_torque_integration_package_gate(
        sample_gate,
        {**model_input_package, "parameter_set_digest": "sha256:jmag-slot365-old-parameter-set"},
        expected_parameter_set_artifact_id="jmag_slot393_airgap_torque_parameter_set_v1.json",
        expected_parameter_set_digest="sha256:jmag-slot393-airgap-torque-parameter-set-v1",
        expected_parameter_set_path="artifacts/motor/jmag_slot393_airgap_torque_parameter_set.json",
        require_parameter_set_artifact=True,
    )
    assert stale_parameter_set_digest["status"] == "needs_attention"
    assert stale_parameter_set_digest["checks"]["expected_parameter_set_artifact_id_matches"] is True
    assert stale_parameter_set_digest["checks"]["expected_parameter_set_digest_matches"] is False
    assert stale_parameter_set_digest["checks"]["expected_parameter_set_path_matches"] is True
    assert stale_parameter_set_digest["checks"]["expected_torque_nm_matches"] is True

    missing_parameter_set_path = jmag_airgap_torque_integration_package_gate(
        sample_gate,
        {key: value for key, value in model_input_package.items() if key != "parameter_set_path"},
        require_parameter_set_artifact=True,
    )
    assert missing_parameter_set_path["status"] == "needs_attention"
    assert missing_parameter_set_path["checks"]["parameter_set_artifact_id_recorded"] is True
    assert missing_parameter_set_path["checks"]["parameter_set_digest_recorded"] is True
    assert missing_parameter_set_path["checks"]["parameter_set_path_recorded"] is False

    wrong_objective_family = jmag_airgap_torque_integration_package_gate(
        sample_gate,
        {**model_input_package, "objective_observable_family": "flux_density_ripple_objective"},
        expected_objective_observable_id="jmag_slot393_airgap_torque_objective_v1",
        expected_objective_observable_family="airgap_torque_maxwell_shear_objective",
    )
    assert wrong_objective_family["status"] == "needs_attention"
    assert wrong_objective_family["checks"]["expected_objective_observable_id_matches"] is True
    assert wrong_objective_family["checks"]["expected_objective_observable_family_matches"] is False
    assert wrong_objective_family["checks"]["expected_torque_nm_matches"] is True

    stale_execution_version = jmag_airgap_torque_integration_package_gate(
        sample_gate,
        {**execution_package, "solver_version": "JMAG-Designer 23.0 stale"},
        expected_solver_version="JMAG-Designer 24.2.1",
        require_execution_metadata=True,
        require_timing_breakdown=True,
    )
    assert stale_execution_version["status"] == "needs_attention"
    assert stale_execution_version["checks"]["expected_solver_version_matches"] is False

    stale_created_run_skew = jmag_airgap_torque_integration_package_gate(
        sample_gate,
        {**execution_package, "created_at_utc": "2026-07-01T13:12:20Z"},
        expected_created_at_utc="2026-07-01T13:12:20Z",
        expected_run_timestamp_utc="2026-07-01T11:12:00Z",
        max_created_run_skew_s=60,
        require_execution_metadata=True,
        require_timing_breakdown=True,
    )
    assert stale_created_run_skew["status"] == "needs_attention"
    assert stale_created_run_skew["checks"]["created_at_utc_parseable"] is True
    assert stale_created_run_skew["checks"]["run_timestamp_utc_parseable"] is True
    assert stale_created_run_skew["checks"]["created_run_timestamp_skew_within_limit"] is False
    assert stale_created_run_skew["checks"]["expected_torque_nm_matches"] is True

    thin_timing_breakdown = jmag_airgap_torque_integration_package_gate(
        sample_gate,
        {**execution_package, "timing_breakdown_s": {"solve_s": 8.0, "write_json_s": 0.2}},
        require_execution_metadata=True,
        require_timing_breakdown=True,
        min_timing_sections=4,
    )
    assert thin_timing_breakdown["status"] == "needs_attention"
    assert thin_timing_breakdown["checks"]["timing_breakdown_has_required_sections"] is False

    impossible_timing_total = jmag_airgap_torque_integration_package_gate(
        sample_gate,
        {**execution_package, "run_duration_s": 4.0},
        require_execution_metadata=True,
        require_timing_breakdown=True,
    )
    assert impossible_timing_total["status"] == "needs_attention"
    assert impossible_timing_total["checks"]["timing_breakdown_total_within_run_duration"] is False

    stale_field = jmag_airgap_torque_integration_package_gate(
        sample_gate,
        {**package, "input_field_table_artifact_id": "jmag_slot357_airgap_flux_table_old.csv"},
        expected_input_field_table_artifact_id="jmag_slot365_airgap_flux_table_v1.csv",
    )
    assert stale_field["status"] == "needs_attention"
    assert stale_field["checks"]["input_field_table_artifact_matches_sample_gate"] is False
    assert stale_field["checks"]["expected_input_field_table_artifact_matches"] is False

    stale_grid = jmag_airgap_torque_integration_package_gate(
        sample_gate,
        {**package, "sample_grid_id": "jmag_slot357_airgap_angle_grid_v1"},
        expected_sample_grid_id="jmag_slot365_airgap_angle_grid_v1",
    )
    assert stale_grid["status"] == "needs_attention"
    assert stale_grid["checks"]["sample_grid_id_matches_sample_gate"] is False
    assert stale_grid["checks"]["expected_sample_grid_id_matches"] is False

    wrong_method = jmag_airgap_torque_integration_package_gate(
        sample_gate,
        {**package, "integration_method": "plain_mean_of_magnitude_only"},
    )
    assert wrong_method["status"] == "needs_attention"
    assert wrong_method["checks"]["expected_integration_method_matches"] is False

    missing_output_digest = jmag_airgap_torque_integration_package_gate(
        sample_gate,
        {key: value for key, value in package.items() if key != "torque_output_digest"},
        expected_torque_output_artifact_id="jmag_slot365_airgap_torque_integration_v1.json",
        expected_torque_output_digest="sha256:jmag-slot365-airgap-torque-integration-v1",
        require_torque_output_artifact=True,
    )
    assert missing_output_digest["status"] == "needs_attention"
    assert missing_output_digest["checks"]["torque_output_digest_recorded"] is False

    stale_torque_output_schema = jmag_airgap_torque_integration_package_gate(
        sample_gate,
        {
            **package,
            "torque_output_schema_id": "jmag_airgap_scalar_torque_v0",
            "torque_output_columns": ["torque_Nm"],
            "torque_output_units": {"torque_Nm": "N*m"},
        },
        expected_torque_output_artifact_id="jmag_slot365_airgap_torque_integration_v1.json",
        expected_torque_output_digest="sha256:jmag-slot365-airgap-torque-integration-v1",
        expected_torque_output_schema_id="jmag_airgap_torque_table_v1",
        expected_torque_output_columns=[
            "RotorAngle_deg",
            "Br_T",
            "Bt_T",
            "torque_density_N_per_m2",
            "torque_Nm",
        ],
        expected_torque_output_units={
            "RotorAngle_deg": "deg",
            "Br_T": "T",
            "Bt_T": "T",
            "torque_density_N_per_m2": "N/m^2",
            "torque_Nm": "N*m",
        },
        expected_torque_nm=0.04125,
        require_torque_output_artifact=True,
        require_torque_output_schema=True,
    )
    assert stale_torque_output_schema["status"] == "needs_attention"
    assert stale_torque_output_schema["checks"]["expected_torque_output_artifact_matches"] is True
    assert stale_torque_output_schema["checks"]["expected_torque_output_digest_matches"] is True
    assert stale_torque_output_schema["checks"]["expected_torque_output_schema_id_matches"] is False
    assert stale_torque_output_schema["checks"]["expected_torque_output_columns_match"] is False
    assert stale_torque_output_schema["checks"]["expected_torque_output_units_match"] is False
    assert stale_torque_output_schema["checks"]["expected_torque_nm_matches"] is True

    stale_torque_convention_schema = jmag_airgap_torque_integration_package_gate(
        sample_gate,
        {
            **package,
            "torque_convention_schema_id": "jmag_airgap_value_only_convention_v0",
        },
        expected_torque_output_schema_id="jmag_airgap_torque_table_v1",
        expected_torque_convention_schema_id="jmag_airgap_torque_convention_v1",
        expected_torque_nm=0.04125,
        require_torque_output_schema=True,
        require_torque_convention_schema=True,
    )
    assert stale_torque_convention_schema["status"] == "needs_attention"
    assert stale_torque_convention_schema["checks"]["expected_torque_output_schema_id_matches"] is True
    assert stale_torque_convention_schema["checks"]["expected_torque_convention_schema_id_matches"] is False
    assert stale_torque_convention_schema["checks"]["expected_torque_nm_matches"] is True

    missing_torque_convention_schema = jmag_airgap_torque_integration_package_gate(
        sample_gate,
        {key: value for key, value in package.items() if key != "torque_convention_schema_id"},
        expected_torque_convention_schema_id="jmag_airgap_torque_convention_v1",
        expected_torque_nm=0.04125,
        require_torque_convention_schema=True,
    )
    assert missing_torque_convention_schema["status"] == "needs_attention"
    assert missing_torque_convention_schema["checks"]["torque_convention_schema_id_recorded"] is False
    assert missing_torque_convention_schema["checks"]["expected_torque_nm_matches"] is True

    stale_torque_component_basis_schema = jmag_airgap_torque_integration_package_gate(
        sample_gate,
        {
            **package,
            "torque_component_basis_schema_id": "jmag_airgap_cartesian_xy_component_basis_v0",
        },
        expected_torque_output_schema_id="jmag_airgap_torque_table_v1",
        expected_torque_convention_schema_id="jmag_airgap_torque_convention_v1",
        expected_torque_component_basis_schema_id="jmag_airgap_cylindrical_rt_component_basis_v1",
        expected_torque_postprocess_row_convention_schema_id="jmag_airgap_torque_row_convention_v1",
        expected_torque_nm=0.04125,
        require_torque_output_schema=True,
        require_torque_convention_schema=True,
        require_torque_component_basis_schema=True,
        require_torque_postprocess_row_convention_schema=True,
    )
    assert stale_torque_component_basis_schema["status"] == "needs_attention"
    assert stale_torque_component_basis_schema["checks"]["expected_torque_output_schema_id_matches"] is True
    assert stale_torque_component_basis_schema["checks"]["expected_torque_convention_schema_id_matches"] is True
    assert (
        stale_torque_component_basis_schema["checks"][
            "expected_torque_component_basis_schema_id_matches"
        ]
        is False
    )
    assert (
        stale_torque_component_basis_schema["checks"][
            "expected_torque_postprocess_row_convention_schema_id_matches"
        ]
        is True
    )
    assert stale_torque_component_basis_schema["checks"]["expected_torque_nm_matches"] is True

    missing_torque_component_basis_schema = jmag_airgap_torque_integration_package_gate(
        sample_gate,
        {
            key: value
            for key, value in package.items()
            if key != "torque_component_basis_schema_id"
        },
        expected_torque_component_basis_schema_id="jmag_airgap_cylindrical_rt_component_basis_v1",
        expected_torque_nm=0.04125,
        require_torque_component_basis_schema=True,
    )
    assert missing_torque_component_basis_schema["status"] == "needs_attention"
    assert (
        missing_torque_component_basis_schema["checks"][
            "torque_component_basis_schema_id_recorded"
        ]
        is False
    )
    assert missing_torque_component_basis_schema["checks"]["expected_torque_nm_matches"] is True

    stale_torque_row_convention_schema = jmag_airgap_torque_integration_package_gate(
        sample_gate,
        {
            **package,
            "torque_postprocess_row_convention_schema_id": "jmag_airgap_scalar_torque_row_v0",
        },
        expected_torque_output_schema_id="jmag_airgap_torque_table_v1",
        expected_torque_convention_schema_id="jmag_airgap_torque_convention_v1",
        expected_torque_postprocess_row_convention_schema_id="jmag_airgap_torque_row_convention_v1",
        expected_torque_nm=0.04125,
        require_torque_output_schema=True,
        require_torque_convention_schema=True,
        require_torque_postprocess_row_convention_schema=True,
    )
    assert stale_torque_row_convention_schema["status"] == "needs_attention"
    assert stale_torque_row_convention_schema["checks"]["expected_torque_output_schema_id_matches"] is True
    assert stale_torque_row_convention_schema["checks"]["expected_torque_convention_schema_id_matches"] is True
    assert (
        stale_torque_row_convention_schema["checks"][
            "expected_torque_postprocess_row_convention_schema_id_matches"
        ]
        is False
    )
    assert stale_torque_row_convention_schema["checks"]["expected_torque_nm_matches"] is True

    missing_torque_row_convention_schema = jmag_airgap_torque_integration_package_gate(
        sample_gate,
        {
            key: value
            for key, value in package.items()
            if key != "torque_postprocess_row_convention_schema_id"
        },
        expected_torque_postprocess_row_convention_schema_id="jmag_airgap_torque_row_convention_v1",
        expected_torque_nm=0.04125,
        require_torque_postprocess_row_convention_schema=True,
    )
    assert missing_torque_row_convention_schema["status"] == "needs_attention"
    assert (
        missing_torque_row_convention_schema["checks"][
            "torque_postprocess_row_convention_schema_id_recorded"
        ]
        is False
    )
    assert missing_torque_row_convention_schema["checks"]["expected_torque_nm_matches"] is True


def test_jmag_symmetry_sweep_coverage_gate_checks_sector_span_before_values():
    rows = [{"RotorAngle_deg": float(theta)} for theta in range(0, 61, 5)]

    gate = jmag_symmetry_sweep_coverage_gate(
        rows,
        pole_pairs=4,
        symmetry_factor=6,
        angle_column="RotorAngle_deg",
        angle_unit="deg",
        angle_basis="mechanical",
        endpoint_policy="included",
    )

    assert gate["status"] == "ok"
    assert gate["policy"] == "jmag_symmetry_sweep_coverage_gate"
    assert gate["expected_mechanical_span_deg"] == pytest.approx(60.0)
    assert gate["covered_mechanical_span_deg"] == pytest.approx(60.0)
    assert gate["covered_electrical_span_deg"] == pytest.approx(240.0)
    assert gate["mean_step_deg_mechanical"] == pytest.approx(5.0)
    assert "before torque/current/harmonic parsing" in gate["version_note"]

    electrical_rows = [
        {"Theta_e_rad": math.radians(float(theta) * 4.0)}
        for theta in range(0, 60, 5)
    ]
    excluded = jmag_symmetry_sweep_coverage_gate(
        electrical_rows,
        pole_pairs=4,
        symmetry_factor=6,
        angle_column="Theta_e_rad",
        angle_unit="rad",
        angle_basis="electrical",
        endpoint_policy="excluded",
    )
    assert excluded["status"] == "ok"
    assert excluded["covered_mechanical_span_deg"] == pytest.approx(60.0)

    nonuniform = [
        {"RotorAngle_deg": value}
        for value in (0.0, 5.0, 10.25, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0)
    ]
    bad_step = jmag_symmetry_sweep_coverage_gate(nonuniform, pole_pairs=4, symmetry_factor=6)
    assert bad_step["status"] == "needs_attention"
    assert bad_step["checks"]["angle_step_uniform"] is False

    wrong_span = [{"RotorAngle_deg": float(theta)} for theta in range(0, 56, 5)]
    bad_span = jmag_symmetry_sweep_coverage_gate(wrong_span, pole_pairs=4, symmetry_factor=6)
    assert bad_span["status"] == "needs_attention"
    assert bad_span["checks"]["mechanical_sector_span_matches_symmetry"] is False


def test_jmag_angle_alignment_contract_gate_closes_gamma_and_rotor_offsets():
    rows = [
        {
            "theta_mech_deg": 0.0,
            "theta_e_deg": 7.5,
            "gamma_jmag_deg": 35.0,
            "gamma_reference_deg": 30.0,
            "symmetry_factor": 8,
        },
        {
            "theta_mech_deg": 2.5,
            "theta_e_deg": 17.5,
            "gamma_jmag_deg": 35.0,
            "gamma_reference_deg": 30.0,
            "symmetry_factor": 8,
        },
        {
            "theta_mech_deg": 5.0,
            "theta_e_deg": 27.5,
            "gamma_jmag_deg": 35.0,
            "gamma_reference_deg": 30.0,
            "symmetry_factor": 8,
        },
    ]
    gate = jmag_angle_alignment_contract_gate(
        rows,
        pole_pairs=4,
        expected_gamma_offset_deg=5.0,
        rotor_electrical_offset_deg=7.5,
        expected_symmetry_factor=8,
    )

    assert gate["policy"] == "jmag_angle_alignment_contract_gate"
    assert gate["status"] == "ok"
    assert gate["max_theta_e_error_deg"] == pytest.approx(0.0)
    assert gate["max_gamma_offset_error_deg"] == pytest.approx(0.0)
    assert gate["mean_mechanical_step_deg"] == pytest.approx(2.5)
    assert all(gate["checks"].values())

    wrong_gamma = [dict(row) for row in rows]
    wrong_gamma[1]["gamma_jmag_deg"] = 30.0
    gamma_gate = jmag_angle_alignment_contract_gate(
        wrong_gamma,
        pole_pairs=4,
        expected_gamma_offset_deg=5.0,
        rotor_electrical_offset_deg=7.5,
        expected_symmetry_factor=8,
    )
    assert gamma_gate["status"] == "needs_attention"
    assert gamma_gate["checks"]["gamma_offset_matches_expected"] is False

    wrong_theta = [dict(row) for row in rows]
    wrong_theta[2]["theta_e_deg"] = 20.0
    theta_gate = jmag_angle_alignment_contract_gate(
        wrong_theta,
        pole_pairs=4,
        expected_gamma_offset_deg=5.0,
        rotor_electrical_offset_deg=7.5,
        expected_symmetry_factor=8,
    )
    assert theta_gate["status"] == "needs_attention"
    assert theta_gate["checks"]["theta_e_matches_pole_pairs_and_offset"] is False

    wrong_symmetry = [dict(row) for row in rows]
    wrong_symmetry[0]["symmetry_factor"] = 1
    symmetry_gate = jmag_angle_alignment_contract_gate(
        wrong_symmetry,
        pole_pairs=4,
        expected_gamma_offset_deg=5.0,
        rotor_electrical_offset_deg=7.5,
        expected_symmetry_factor=8,
    )
    assert symmetry_gate["status"] == "needs_attention"
    assert symmetry_gate["checks"]["symmetry_factor_matches_expected"] is False


def test_jmag_export_case_package_gate_keeps_case_study_and_result_set_ids():
    artifacts = [
        {
            "kind": "column_metadata",
            "case_id": "case_fw_004",
            "study_id": "pm_drive_map",
            "result_set_id": "resultset_20260629_A",
            "source_tool": "JMAG-Designer",
            "path": "slot149_columns.json",
            "gate_policy": "jmag_motor_table_column_metadata_gate",
            "status": "ok",
        },
        {
            "kind": "symmetry_coverage",
            "case_id": "case_fw_004",
            "study_id": "pm_drive_map",
            "result_set_id": "resultset_20260629_A",
            "source_tool": "JMAG",
            "path": "slot149_symmetry.json",
            "gate_policy": "jmag_symmetry_sweep_coverage_gate",
            "status": "ok",
        },
        {
            "kind": "value_table",
            "case_id": "case_fw_004",
            "study_id": "pm_drive_map",
            "result_set_id": "resultset_20260629_A",
            "source_tool": "JMAG",
            "path": "slot149_pm_drive_table.csv",
            "gate_policy": "pm_drive_terminal_table_health_gate",
            "status": "ok",
            "operating_point_ids": ["MTPA", "FW"],
        },
        {
            "kind": "notebook_row",
            "case_id": "case_fw_004",
            "study_id": "pm_drive_map",
            "result_set_id": "resultset_20260629_A",
            "source_tool": "JMAG",
            "path": "slot149_notebook_row.json",
            "gate_policy": "pm_drive_operating_point_notebook_handoff_gate",
            "status": "ok",
            "operating_point_id": "FW",
        },
    ]

    gate = jmag_export_case_package_gate(
        artifacts,
        expected_case_id="case_fw_004",
        expected_study_id="pm_drive_map",
        expected_result_set_id="resultset_20260629_A",
    )

    assert gate["status"] == "ok"
    assert gate["policy"] == "jmag_export_case_package_gate"
    assert gate["case_ids"] == ["case_fw_004"]
    assert gate["study_ids"] == ["pm_drive_map"]
    assert gate["result_set_ids"] == ["resultset_20260629_A"]
    assert gate["notebook_operating_point_ids"] == ["FW"]
    assert gate["checks"]["notebook_operating_point_in_value_table"] is True

    stale_result = [dict(row) for row in artifacts]
    stale_result[3]["result_set_id"] = "resultset_old"
    stale_result_gate = jmag_export_case_package_gate(stale_result)
    assert stale_result_gate["status"] == "needs_attention"
    assert stale_result_gate["checks"]["result_set_ids_unique"] is False

    missing_case = [dict(row) for row in artifacts]
    missing_case[0].pop("case_id")
    missing_case_gate = jmag_export_case_package_gate(missing_case)
    assert missing_case_gate["status"] == "needs_attention"
    assert missing_case_gate["checks"]["case_ids_present"] is False

    missing_op_in_table = [dict(row) for row in artifacts]
    missing_op_in_table[2]["operating_point_ids"] = ["MTPA"]
    missing_op_gate = jmag_export_case_package_gate(missing_op_in_table)
    assert missing_op_gate["status"] == "needs_attention"
    assert missing_op_gate["checks"]["notebook_operating_point_in_value_table"] is False

    wrong_source = [dict(row) for row in artifacts]
    wrong_source[1]["source_tool"] = "FEMM"
    wrong_source_gate = jmag_export_case_package_gate(wrong_source)
    assert wrong_source_gate["status"] == "needs_attention"
    assert wrong_source_gate["checks"]["source_tool_is_jmag"] is False


def test_jmag_current_torque_solver_ready_manifest_keeps_current_and_torque_locked():
    artifacts = [
        {
            "kind": "column_metadata",
            "case_id": "case_torque_006",
            "result_set_id": "resultset_20260630_B",
            "source_tool": "JMAG-Designer",
            "path": "slot165_columns.json",
            "gate_policy": "jmag_motor_table_column_metadata_gate",
            "status": "ok",
        },
        {
            "kind": "symmetry_coverage",
            "case_id": "case_torque_006",
            "result_set_id": "resultset_20260630_B",
            "source_tool": "JMAG",
            "path": "slot165_symmetry.json",
            "gate_policy": "jmag_symmetry_sweep_coverage_gate",
            "status": "ok",
        },
        {
            "kind": "current_snapshot",
            "case_id": "case_torque_006",
            "result_set_id": "resultset_20260630_B",
            "operating_point_id": "id-4_iq18_theta22p5",
            "source_tool": "JMAG",
            "path": "slot165_current_snapshot.json",
            "gate_policy": "motor_current_snapshot_table_contract_gate",
            "status": "ok",
            "current_kind": "instantaneous",
            "phase_set": ["U", "V", "W"],
        },
        {
            "kind": "torque_table",
            "case_id": "case_torque_006",
            "result_set_id": "resultset_20260630_B",
            "operating_point_id": "id-4_iq18_theta22p5",
            "source_tool": "JMAG",
            "path": "slot165_torque_table.json",
            "gate_policy": "torque_angle_table_export_health",
            "status": "ok",
            "phase_set": ["U", "V", "W"],
            "rotor_current_phase_locked": True,
        },
    ]

    gate = jmag_current_torque_solver_ready_manifest_gate(
        artifacts,
        expected_case_id="case_torque_006",
        expected_result_set_id="resultset_20260630_B",
        expected_operating_point_id="id-4_iq18_theta22p5",
    )

    assert gate["policy"] == "jmag_current_torque_solver_ready_manifest_gate"
    assert gate["status"] == "ok"
    assert gate["present_kinds"] == {
        "column_metadata": 1,
        "current_snapshot": 1,
        "symmetry_coverage": 1,
        "torque_table": 1,
    }
    assert gate["checks"]["operating_point_ids_unique"] is True
    assert gate["checks"]["torque_table_locked_to_current_phase"] is True
    assert "before table values are joined" in gate["version_note"]

    stale = [dict(row) for row in artifacts]
    stale[3]["result_set_id"] = "resultset_old"
    stale_gate = jmag_current_torque_solver_ready_manifest_gate(stale)
    assert stale_gate["status"] == "needs_attention"
    assert stale_gate["checks"]["result_set_ids_unique"] is False

    rms = [dict(row) for row in artifacts]
    rms[2]["current_kind"] = "rms"
    rms_gate = jmag_current_torque_solver_ready_manifest_gate(rms)
    assert rms_gate["status"] == "needs_attention"
    assert rms_gate["checks"]["current_snapshot_is_instantaneous"] is False

    unlocked = [dict(row) for row in artifacts]
    unlocked[3]["rotor_current_phase_locked"] = False
    unlocked_gate = jmag_current_torque_solver_ready_manifest_gate(unlocked)
    assert unlocked_gate["status"] == "needs_attention"
    assert unlocked_gate["checks"]["torque_table_locked_to_current_phase"] is False

    missing_torque_op = [dict(row) for row in artifacts]
    missing_torque_op[3].pop("operating_point_id")
    missing_gate = jmag_current_torque_solver_ready_manifest_gate(missing_torque_op)
    assert missing_gate["status"] == "needs_attention"
    assert missing_gate["checks"]["operating_point_ids_present_for_current_and_torque"] is False


def test_jmag_efficiency_operating_point_package_keeps_map_rows_aligned():
    artifacts = [
        {
            "kind": "terminal_table",
            "case_id": "case_fw_004",
            "result_set_id": "resultset_20260629_A",
            "source_tool": "JMAG",
            "path": "slot157_terminal.json",
            "gate_policy": "pm_drive_terminal_table_health",
            "status": "ok",
            "operating_point_ids": ["MTPA", "FW", "high_current"],
        },
        {
            "kind": "loss_bucket_table",
            "case_id": "case_fw_004",
            "result_set_id": "resultset_20260629_A",
            "source_tool": "JMAG-Designer",
            "path": "slot157_losses.json",
            "gate_policy": "pm_drive_loss_bucket_efficiency_gate",
            "status": "ok",
            "operating_point_ids": ["MTPA", "FW", "high_current"],
        },
        {
            "kind": "drive_cycle_weights",
            "case_id": "case_fw_004",
            "result_set_id": "resultset_20260629_A",
            "source_tool": "JMAG",
            "path": "slot157_drive_cycle.json",
            "gate_policy": "drive_cycle_weighted_efficiency_gate",
            "status": "ok",
            "operating_point_ids": ["MTPA", "FW", "high_current"],
        },
        {
            "kind": "notebook_row",
            "case_id": "case_fw_004",
            "result_set_id": "resultset_20260629_A",
            "source_tool": "JMAG",
            "path": "slot157_notebook_row.json",
            "gate_policy": "pm_drive_operating_point_notebook_handoff_gate",
            "status": "ok",
            "operating_point_id": "FW",
        },
    ]

    gate = jmag_efficiency_operating_point_package_gate(
        artifacts,
        expected_case_id="case_fw_004",
        expected_result_set_id="resultset_20260629_A",
    )

    assert gate["policy"] == "jmag_efficiency_operating_point_package_gate"
    assert gate["status"] == "ok"
    assert gate["reference_operating_point_ids"] == ["FW", "MTPA", "high_current"]
    assert gate["checks"]["table_operating_point_sets_match"] is True
    assert gate["checks"]["notebook_operating_points_in_tables"] is True

    stale_result = [dict(row) for row in artifacts]
    stale_result[1]["result_set_id"] = "resultset_old"
    stale_gate = jmag_efficiency_operating_point_package_gate(stale_result)
    assert stale_gate["status"] == "needs_attention"
    assert stale_gate["checks"]["result_set_ids_unique"] is False

    missing_loss_op = [dict(row) for row in artifacts]
    missing_loss_op[1]["operating_point_ids"] = ["MTPA", "FW"]
    mismatch_gate = jmag_efficiency_operating_point_package_gate(missing_loss_op)
    assert mismatch_gate["status"] == "needs_attention"
    assert mismatch_gate["checks"]["table_operating_point_sets_match"] is False

    notebook_outside = [dict(row) for row in artifacts]
    notebook_outside[3]["operating_point_id"] = "burst"
    outside_gate = jmag_efficiency_operating_point_package_gate(notebook_outside)
    assert outside_gate["status"] == "needs_attention"
    assert outside_gate["checks"]["notebook_operating_points_in_tables"] is False


def test_two_port_s_to_yz_equivalent_gate_keeps_reference_impedance_contract():
    gate_50 = two_port_s_to_yz_equivalent_gate(0.0, 0.0, s12=0.0, s22=0.0, z0=50.0)
    gate_75 = two_port_s_to_yz_equivalent_gate(0.0, 0.0, s12=0.0, s22=0.0, z0=75.0)

    assert gate_50["policy"] == "two_port_s_to_yz_equivalent_gate"
    assert gate_50["status"] == "ok"
    assert gate_50["checks"]["sparameter_passivity_ok"] is True
    assert gate_50["checks"]["y_matrix_defined"] is True
    assert gate_50["checks"]["z_matrix_defined"] is True
    assert gate_50["y_shunt1"]["real"] == pytest.approx(1.0 / 50.0)
    assert gate_50["y_shunt2"]["real"] == pytest.approx(1.0 / 50.0)
    assert gate_50["y_series"]["abs"] == pytest.approx(0.0)
    assert gate_50["z11"]["real"] == pytest.approx(50.0)
    assert gate_50["z_series1"]["real"] == pytest.approx(50.0)
    assert gate_50["z_series2"]["real"] == pytest.approx(50.0)
    assert gate_50["z_shunt"]["abs"] == pytest.approx(0.0)
    assert gate_75["status"] == "ok"
    assert gate_75["y_shunt1"]["real"] == pytest.approx(1.0 / 75.0)
    assert gate_75["z11"]["real"] == pytest.approx(75.0)

    with pytest.raises(ValueError):
        two_port_s_to_yz_equivalent_gate(0.0, 0.0, z0=0.0)


def test_rc_snubber_one_pole_rows_stay_passive_reciprocal_and_lagging():
    cutoff = 1.0
    freqs = [0.25, 0.5, 1.0, 2.0, 4.0]
    rows = []
    for freq in freqs:
        h = 1.0 / complex(1.0, freq / cutoff)
        rows.append({
            "frequency": freq,
            "h": h,
            "health": two_port_sparameter_health(0.0, h, s12=h, s22=0.0),
        })

    magnitudes = [abs(row["h"]) for row in rows]
    phases = [math.degrees(math.atan2(row["h"].imag, row["h"].real)) for row in rows]
    assert all(row["health"]["status"] == "ok" for row in rows)
    assert all(a >= b for a, b in zip(magnitudes, magnitudes[1:]))
    assert all(phase <= 0.0 for phase in phases)
    assert magnitudes[2] == pytest.approx(1.0 / math.sqrt(2.0))
    assert phases[2] == pytest.approx(-45.0)

    active = two_port_sparameter_health(0.0, 1.02, s12=1.02, s22=0.0)
    assert active["status"] == "needs_attention"
    assert active["passive"] is False


def test_quarter_wave_directional_coupler_gate_for_cst_slot_learning():
    gate = quarter_wave_directional_coupler_gate(1.0 / math.sqrt(2.0), z0=50.0)
    assert gate["status"] == "ok"
    assert gate["impedance_product"] == pytest.approx(50.0 ** 2)
    assert gate["power_sum"] == pytest.approx(1.0)
    assert gate["matched"] is True
    assert gate["isolated"] is True
    assert gate["lossless"] is True
    assert gate["z0_even"] == pytest.approx(120.71067811865474)
    assert gate["z0_odd"] == pytest.approx(20.710678118654755)
    assert gate["s21"]["abs"] == pytest.approx(1.0 / math.sqrt(2.0))
    assert gate["s31"]["abs"] == pytest.approx(1.0 / math.sqrt(2.0))


def test_branch_line_hybrid_gate_for_cst_port_order_and_phase():
    gate = branch_line_hybrid_gate(z0=50.0)

    assert gate["status"] == "ok"
    assert gate["through_branch_impedance"] == pytest.approx(50.0)
    assert gate["shunt_branch_impedance"] == pytest.approx(50.0 / math.sqrt(2.0))
    assert gate["s11"]["abs"] == pytest.approx(0.0)
    assert gate["s41"]["abs"] == pytest.approx(0.0)
    assert gate["s21"]["abs"] == pytest.approx(1.0 / math.sqrt(2.0))
    assert gate["s31"]["abs"] == pytest.approx(1.0 / math.sqrt(2.0))
    assert gate["split_power_sum"] == pytest.approx(1.0)
    assert gate["through_phase_deg"] == pytest.approx(-90.0)
    assert gate["coupled_phase_deg"] == pytest.approx(180.0)
    assert gate["phase_difference_deg"] == pytest.approx(-90.0)
    assert gate["max_column_power_error"] == pytest.approx(0.0)
    assert gate["max_orthogonality_error"] == pytest.approx(0.0)


def test_three_phase_dq_current_handoff_roundtrip_and_phase_sequence_guard():
    id_ref = -3.0
    iq_ref = 12.0
    theta = 0.0
    abc = dq_to_three_phase_currents(id_ref, iq_ref, theta)

    summary = three_phase_currents_to_dq_summary(
        abc,
        theta,
        expected_id=id_ref,
        expected_iq=iq_ref,
    )

    assert abc == pytest.approx({
        "U": -3.0,
        "V": 11.892304845413264,
        "W": -8.892304845413264,
    })
    assert summary["status"] == "ok"
    assert summary["id"] == pytest.approx(id_ref)
    assert summary["iq"] == pytest.approx(iq_ref)
    assert summary["zero_sequence_abs"] == pytest.approx(0.0)
    assert summary["abc_square_sum"] == pytest.approx(1.5 * (id_ref ** 2 + iq_ref ** 2))

    swapped = {"U": abc["U"], "V": abc["W"], "W": abc["V"]}
    wrong = three_phase_currents_to_dq_summary(
        swapped,
        theta,
        expected_id=id_ref,
        expected_iq=iq_ref,
    )

    assert wrong["checks"]["zero_sequence_ok"] is True
    assert wrong["checks"]["iq_ok"] is False
    assert wrong["iq"] == pytest.approx(-iq_ref)
    assert wrong["iq_abs_error"] == pytest.approx(24.0)


def test_spwm_snapshot_current_handoff_samples_one_electrical_period():
    summary = spwm_snapshot_current_handoff_summary(
        id_current=-2.5,
        iq_current=11.0,
        sample_count=24,
        sample_offset_fraction=0.5,
        carrier_ratio=12,
    )

    assert summary["status"] == "ok"
    assert summary["policy"] == "spwm_snapshot_current_handoff_gate"
    assert summary["sample_count"] == 24
    assert summary["carrier_ratio"] == pytest.approx(12.0)
    assert summary["current_amplitude"] == pytest.approx(math.hypot(-2.5, 11.0))
    assert summary["expected_phase_rms"] == pytest.approx(math.hypot(-2.5, 11.0) / math.sqrt(2.0))
    assert summary["max_id_abs_error"] < 1.0e-14
    assert summary["max_iq_abs_error"] < 1.0e-14
    assert summary["max_zero_sequence_abs"] < 1.0e-14
    assert summary["max_abc_square_sum_error"] < 1.0e-13
    assert summary["max_phase_rms_abs_error"] < 1.0e-14
    assert summary["rows"][0]["theta_e_deg"] == pytest.approx(7.5)
    assert summary["rows"][0]["dq"]["status"] == "ok"


def test_spwm_snapshot_sampling_mode_metadata_checks_timer_offset():
    symmetric = spwm_snapshot_current_handoff_summary(
        id_current=-2.5,
        iq_current=11.0,
        sample_count=24,
        sample_offset_fraction=0.5,
        sampling_mode="symmetrical",
        timer_alignment="center-aligned",
        carrier_ratio=12,
    )
    asymmetric = spwm_snapshot_current_handoff_summary(
        id_current=-2.5,
        iq_current=11.0,
        sample_count=24,
        sample_offset_fraction=0.0,
        sampling_mode="asymmetrical",
        timer_alignment="edge-aligned",
        carrier_ratio=12,
    )
    wrong = spwm_snapshot_current_handoff_summary(
        id_current=-2.5,
        iq_current=11.0,
        sample_count=24,
        sample_offset_fraction=0.0,
        sampling_mode="symmetrical",
        timer_alignment="center-aligned",
        carrier_ratio=12,
    )

    assert symmetric["status"] == "ok"
    assert symmetric["sampling_mode"] == "symmetrical"
    assert symmetric["timer_alignment"] == "center_aligned"
    assert symmetric["checks"]["sampling_offset_matches_mode"] is True
    assert asymmetric["status"] == "ok"
    assert asymmetric["sampling_mode"] == "asymmetrical"
    assert asymmetric["checks"]["sampling_offset_matches_mode"] is True
    assert wrong["status"] == "needs_attention"
    assert wrong["checks"]["sampling_offset_matches_mode"] is False


def test_spwm_snapshot_rows_feed_jmag_style_dq_torque_table_contract():
    summary = spwm_snapshot_current_handoff_summary(
        id_current=-2.5,
        iq_current=11.0,
        sample_count=24,
        sample_offset_fraction=0.5,
        carrier_ratio=12,
    )
    current = summary["current_amplitude"]
    gamma_deg = math.degrees(math.atan2(-summary["id"], summary["iq"]))
    lambda_m = 0.10
    Ld = 4.0e-3
    Lq = 9.0e-3
    pole_pairs = 3
    torque = lumped_pm_dq_torque(lambda_m, Ld, Lq, summary["id"], summary["iq"], pole_pairs)

    row = {
        "sample_index": summary["rows"][0]["sample"],
        "theta_e_deg": summary["rows"][0]["theta_e_deg"],
        "theta_mech_deg": summary["rows"][0]["theta_e_deg"] / pole_pairs,
        "current_U_A": summary["rows"][0]["currents"]["U"],
        "current_V_A": summary["rows"][0]["currents"]["V"],
        "current_W_A": summary["rows"][0]["currents"]["W"],
        "id_A": summary["id"],
        "iq_A": summary["iq"],
        "gamma_deg": gamma_deg,
        "torque_Nm": torque,
        "carrier_ratio": summary["carrier_ratio"],
        "sample_offset_fraction": summary["sample_offset_fraction"],
        "phase_order": "U,V,W",
    }
    health = dq_torque_table_health(
        [
            {"gamma_deg": -30.0, "id_A": 0.5 * current, "iq_A": math.sqrt(3.0) * current / 2.0,
             "torque_Nm": lumped_pm_dq_torque(lambda_m, Ld, Lq, 0.5 * current, math.sqrt(3.0) * current / 2.0, pole_pairs)},
            {"gamma_deg": 0.0, "id_A": 0.0, "iq_A": current,
             "torque_Nm": lumped_pm_dq_torque(lambda_m, Ld, Lq, 0.0, current, pole_pairs)},
            {"gamma_deg": row["gamma_deg"], "id_A": row["id_A"], "iq_A": row["iq_A"], "torque_Nm": row["torque_Nm"]},
            {"gamma_deg": 30.0, "id_A": -0.5 * current, "iq_A": math.sqrt(3.0) * current / 2.0,
             "torque_Nm": lumped_pm_dq_torque(lambda_m, Ld, Lq, -0.5 * current, math.sqrt(3.0) * current / 2.0, pole_pairs)},
        ],
        lambda_m=lambda_m,
        Ld=Ld,
        Lq=Lq,
        current=current,
        pole_pairs=pole_pairs,
    )

    assert summary["status"] == "ok"
    assert row["theta_mech_deg"] == pytest.approx(2.5)
    assert row["gamma_deg"] == pytest.approx(12.80426606528675)
    assert row["torque_Nm"] == pytest.approx(5.56875)
    assert row["carrier_ratio"] == pytest.approx(12.0)
    assert health["status"] == "ok"
    assert health["peak_row"]["gamma_deg"] == pytest.approx(30.0)


def test_motor_current_snapshot_table_contract_closes_angle_and_metadata_rows():
    summary = spwm_snapshot_current_handoff_summary(
        id_current=-2.5,
        iq_current=11.0,
        sample_count=24,
        sample_offset_fraction=0.5,
        sampling_mode="symmetrical",
        timer_alignment="center-aligned",
        carrier_ratio=12,
    )
    pole_pairs = 3
    gamma_deg = math.degrees(math.atan2(-summary["id"], summary["iq"]))
    current = math.hypot(summary["id"], summary["iq"])
    rows = []
    for item in summary["rows"][:4]:
        rows.append({
            "sample_index": item["sample"],
            "theta_e_deg": item["theta_e_deg"],
            "theta_mech_deg": item["theta_e_deg"] / pole_pairs,
            "current_U_A": item["currents"]["U"],
            "current_V_A": item["currents"]["V"],
            "current_W_A": item["currents"]["W"],
            "id_A": summary["id"],
            "iq_A": summary["iq"],
            "gamma_deg": gamma_deg,
            "current_A": current,
            "current_kind": "instantaneous",
            "sampling_mode": summary["sampling_mode"],
            "timer_alignment": summary["timer_alignment"],
            "carrier_ratio": summary["carrier_ratio"],
            "sample_offset_fraction": summary["sample_offset_fraction"],
        })

    gate = motor_current_snapshot_table_contract_gate(rows, pole_pairs=pole_pairs)

    assert gate["status"] == "ok"
    assert gate["policy"] == "motor_current_snapshot_table_contract_gate"
    assert gate["n_rows"] == 4
    assert gate["max_angle_abs_error_rad"] == pytest.approx(0.0, abs=1.0e-15)
    assert gate["max_id_abs_error_A"] < 1.0e-14
    assert gate["max_iq_abs_error_A"] < 1.0e-14
    assert gate["max_gamma_abs_error_deg"] == pytest.approx(0.0)
    assert gate["max_current_abs_error_A"] == pytest.approx(0.0)

    wrong_mech = [dict(row) for row in rows]
    wrong_mech[0]["theta_mech_deg"] = wrong_mech[0]["theta_e_deg"]
    wrong_mech_gate = motor_current_snapshot_table_contract_gate(wrong_mech, pole_pairs=pole_pairs)
    assert wrong_mech_gate["status"] == "needs_attention"
    assert wrong_mech_gate["checks"]["mechanical_electrical_angle_ok"] is False

    rms_rows = [dict(row, current_kind="rms") for row in rows]
    rms_gate = motor_current_snapshot_table_contract_gate(rms_rows, pole_pairs=pole_pairs)
    assert rms_gate["status"] == "needs_attention"
    assert rms_gate["checks"]["current_kind_matches"] is False

    swapped = [dict(row) for row in rows]
    swapped[0]["current_V_A"], swapped[0]["current_W_A"] = swapped[0]["current_W_A"], swapped[0]["current_V_A"]
    swapped_gate = motor_current_snapshot_table_contract_gate(swapped, pole_pairs=pole_pairs)
    assert swapped_gate["status"] == "needs_attention"
    assert swapped_gate["checks"]["dq_recovery_ok"] is False


def test_femm_motor_model_artifact_package_keeps_model_and_operating_point_ids():
    artifacts = [
        {
            "kind": "block_labels",
            "model_id": "ipm_teaching_v1",
            "source_tool": "FEMM",
            "path": "slot148_block_labels.json",
            "gate_policy": "femm_block_label_source_contract_gate",
            "status": "ok",
        },
        {
            "kind": "current_snapshot",
            "model_id": "ipm_teaching_v1",
            "operating_point_id": "id-3_iq12_theta0",
            "source_tool": "FEMM",
            "path": "slot148_current_snapshot.json",
            "gate_policy": "motor_current_snapshot_table_contract_gate",
            "status": "ok",
            "current_kind": "instantaneous",
        },
        {
            "kind": "torque_table",
            "model_id": "ipm_teaching_v1",
            "operating_point_id": "id-3_iq12_theta0",
            "source_tool": "FEMM",
            "path": "slot148_torque_angle.csv",
            "gate_policy": "torque_angle_table_export_health",
            "status": "ok",
            "angle_basis": "mechanical",
            "source_function": "mo_blockintegral(22)",
            "rotor_current_phase_locked": True,
        },
    ]

    gate = femm_motor_model_artifact_package_gate(
        artifacts,
        expected_model_id="ipm_teaching_v1",
        expected_operating_point_id="id-3_iq12_theta0",
    )

    assert gate["status"] == "ok"
    assert gate["policy"] == "femm_motor_model_artifact_package_gate"
    assert gate["model_ids"] == ["ipm_teaching_v1"]
    assert gate["operating_point_ids"] == ["id-3_iq12_theta0"]
    assert gate["checks"]["required_kinds_present"] is True
    assert gate["checks"]["upstream_gate_policy_known"] is True
    assert "different motor models or operating points" in gate["version_note"]

    wrong_model = [dict(row) for row in artifacts]
    wrong_model[2]["model_id"] = "stale_motor_v0"
    wrong_model_gate = femm_motor_model_artifact_package_gate(wrong_model)
    assert wrong_model_gate["status"] == "needs_attention"
    assert wrong_model_gate["checks"]["model_ids_unique"] is False

    missing_op = [dict(row) for row in artifacts]
    missing_op[1].pop("operating_point_id")
    missing_op_gate = femm_motor_model_artifact_package_gate(missing_op)
    assert missing_op_gate["status"] == "needs_attention"
    assert missing_op_gate["checks"]["operating_point_ids_present_for_current_and_torque"] is False

    bad_torque = [dict(row) for row in artifacts]
    bad_torque[2]["angle_basis"] = "electrical"
    bad_torque[2]["rotor_current_phase_locked"] = False
    bad_torque_gate = femm_motor_model_artifact_package_gate(bad_torque)
    assert bad_torque_gate["status"] == "needs_attention"
    assert bad_torque_gate["checks"]["torque_table_metadata_solver_ready"] is False


def test_femm_winding_current_package_keeps_model_and_phase_ids():
    artifacts = [
        {
            "kind": "winding_table",
            "model_id": "ipm_teaching_v1",
            "source_tool": "analytic",
            "path": "slot156_winding_factor.json",
            "gate_policy": "double_layer_winding_pitch_harmonic_gate",
            "status": "ok",
            "slots": 24,
            "poles": 4,
            "phase_count": 3,
            "phase_set": ["U", "V", "W"],
        },
        {
            "kind": "block_labels",
            "model_id": "ipm_teaching_v1",
            "source_tool": "FEMM",
            "path": "slot156_block_labels.json",
            "gate_policy": "femm_block_label_source_contract_gate",
            "status": "ok",
            "phase_set": ["U", "V", "W"],
        },
        {
            "kind": "current_snapshot",
            "model_id": "ipm_teaching_v1",
            "source_tool": "FEMM",
            "path": "slot156_current_snapshot.json",
            "gate_policy": "femm_static_current_circuit_rows_gate",
            "status": "ok",
            "current_kind": "instantaneous",
            "phase_set": ["U", "V", "W"],
        },
    ]

    gate = femm_winding_current_package_gate(
        artifacts,
        expected_model_id="ipm_teaching_v1",
    )

    assert gate["policy"] == "femm_winding_current_package_gate"
    assert gate["status"] == "ok"
    assert gate["checks"]["winding_geometry_metadata_present"] is True
    assert gate["checks"]["phase_sets_match_expected"] is True
    assert "winding-factor table cannot be mixed" in gate["version_note"]

    stale = [dict(row) for row in artifacts]
    stale[0]["model_id"] = "old_winding_v0"
    stale_gate = femm_winding_current_package_gate(stale)
    assert stale_gate["status"] == "needs_attention"
    assert stale_gate["checks"]["model_ids_unique"] is False

    wrong_phase = [dict(row) for row in artifacts]
    wrong_phase[1]["phase_set"] = ["U", "W", "V"]
    wrong_phase_gate = femm_winding_current_package_gate(wrong_phase)
    assert wrong_phase_gate["status"] == "needs_attention"
    assert wrong_phase_gate["checks"]["phase_sets_match_expected"] is False

    rms_current = [dict(row) for row in artifacts]
    rms_current[2]["current_kind"] = "rms"
    rms_gate = femm_winding_current_package_gate(rms_current)
    assert rms_gate["status"] == "needs_attention"
    assert rms_gate["checks"]["current_snapshot_is_instantaneous"] is False


def test_femm_source_current_solver_ready_manifest_keeps_presolve_sources_together():
    artifacts = [
        {
            "kind": "block_labels",
            "model_id": "ipm_teaching_v1",
            "source_tool": "FEMM",
            "path": "slot164_block_labels.json",
            "gate_policy": "femm_block_label_source_contract_gate",
            "status": "ok",
            "phase_set": ["U", "V", "W"],
        },
        {
            "kind": "pm_magnetization",
            "model_id": "ipm_teaching_v1",
            "source_tool": "FEMM",
            "path": "slot164_pm_magnetization.json",
            "gate_policy": "femm_pm_magnetization_convention_gate",
            "status": "ok",
        },
        {
            "kind": "current_snapshot",
            "model_id": "ipm_teaching_v1",
            "operating_point_id": "id-3_iq12_theta17p5",
            "source_tool": "pyFEMM",
            "path": "slot164_current_snapshot.json",
            "gate_policy": "femm_static_current_circuit_rows_gate",
            "status": "ok",
            "current_kind": "instantaneous",
            "phase_set": ["U", "V", "W"],
        },
    ]

    gate = femm_source_current_solver_ready_manifest_gate(
        artifacts,
        expected_model_id="ipm_teaching_v1",
        expected_operating_point_id="id-3_iq12_theta17p5",
    )

    assert gate["policy"] == "femm_source_current_solver_ready_manifest_gate"
    assert gate["status"] == "ok"
    assert gate["present_kinds"] == {
        "block_labels": 1,
        "current_snapshot": 1,
        "pm_magnetization": 1,
    }
    assert gate["checks"]["required_kinds_present"] is True
    assert gate["checks"]["current_snapshot_is_instantaneous"] is True
    assert "before FEMM solve/result tables" in gate["version_note"]

    missing_pm = femm_source_current_solver_ready_manifest_gate(
        [artifacts[0], artifacts[2]],
        expected_model_id="ipm_teaching_v1",
    )
    assert missing_pm["status"] == "needs_attention"
    assert missing_pm["checks"]["required_kinds_present"] is False

    rms_current = [dict(row) for row in artifacts]
    rms_current[2]["current_kind"] = "rms"
    rms_gate = femm_source_current_solver_ready_manifest_gate(rms_current)
    assert rms_gate["status"] == "needs_attention"
    assert rms_gate["checks"]["current_snapshot_is_instantaneous"] is False

    stale_pm = [dict(row) for row in artifacts]
    stale_pm[1]["model_id"] = "stale_pm_model"
    stale_gate = femm_source_current_solver_ready_manifest_gate(stale_pm)
    assert stale_gate["status"] == "needs_attention"
    assert stale_gate["checks"]["model_ids_unique"] is False


def test_femm_air_gap_sample_manifest_keeps_source_samples_and_torque_locked():
    artifacts = [
        {
            "kind": "source_current_manifest",
            "model_id": "ipm_teaching_v1",
            "operating_point_id": "id-3_iq12_theta17p5",
            "source_tool": "FEMM",
            "path": "slot172_source_current_manifest.json",
            "gate_policy": "femm_source_current_solver_ready_manifest_gate",
            "status": "ok",
        },
        {
            "kind": "air_gap_sample_table",
            "model_id": "ipm_teaching_v1",
            "operating_point_id": "id-3_iq12_theta17p5",
            "source_tool": "pyFEMM",
            "path": "slot172_gap_samples.json",
            "gate_policy": "femm_air_gap_sample_metadata_contract",
            "status": "ok",
            "angle_unit": "deg",
            "component_frame": "cylindrical_rt",
            "radius_m": 0.045,
            "axial_length_m": 0.080,
            "torque_sign_convention": "positive_ccw",
        },
        {
            "kind": "torque_summary",
            "model_id": "ipm_teaching_v1",
            "operating_point_id": "id-3_iq12_theta17p5",
            "source_tool": "radia-ngsolve",
            "path": "slot172_air_gap_shear_torque.json",
            "gate_policy": "air_gap_shear_torque_from_angle_samples",
            "status": "ok",
            "torque_sign_convention": "positive_ccw",
        },
    ]

    gate = femm_air_gap_sample_solver_ready_manifest_gate(
        artifacts,
        expected_model_id="ipm_teaching_v1",
        expected_operating_point_id="id-3_iq12_theta17p5",
    )

    assert gate["policy"] == "femm_air_gap_sample_solver_ready_manifest_gate"
    assert gate["status"] == "ok"
    assert gate["present_kinds"] == {
        "air_gap_sample_table": 1,
        "source_current_manifest": 1,
        "torque_summary": 1,
    }
    assert gate["checks"]["sample_angle_unit_is_deg"] is True
    assert gate["checks"]["sample_component_frame_is_cylindrical_rt"] is True
    assert gate["checks"]["sample_radius_and_length_positive"] is True
    assert gate["checks"]["torque_sign_convention_consistent"] is True
    assert "before promoting Br/Bt rows" in gate["version_note"]

    missing_sample = femm_air_gap_sample_solver_ready_manifest_gate(
        [artifacts[0], artifacts[2]],
        expected_model_id="ipm_teaching_v1",
    )
    assert missing_sample["status"] == "needs_attention"
    assert missing_sample["checks"]["required_kinds_present"] is False

    stale_op = [dict(row) for row in artifacts]
    stale_op[1]["operating_point_id"] = "theta_old"
    stale_gate = femm_air_gap_sample_solver_ready_manifest_gate(stale_op)
    assert stale_gate["status"] == "needs_attention"
    assert stale_gate["checks"]["operating_point_ids_unique"] is False

    rad_angle = [dict(row) for row in artifacts]
    rad_angle[1]["angle_unit"] = "rad"
    rad_gate = femm_air_gap_sample_solver_ready_manifest_gate(rad_angle)
    assert rad_gate["status"] == "needs_attention"
    assert rad_gate["checks"]["sample_angle_unit_is_deg"] is False

    cartesian_frame = [dict(row) for row in artifacts]
    cartesian_frame[1]["component_frame"] = "cartesian_xy"
    frame_gate = femm_air_gap_sample_solver_ready_manifest_gate(cartesian_frame)
    assert frame_gate["status"] == "needs_attention"
    assert frame_gate["checks"]["sample_component_frame_is_cylindrical_rt"] is False

    stale_status = [dict(row) for row in artifacts]
    stale_status[0]["status"] = "needs_attention"
    status_gate = femm_air_gap_sample_solver_ready_manifest_gate(stale_status)
    assert status_gate["status"] == "needs_attention"
    assert status_gate["checks"]["upstream_gate_status_ok"] is False

    sign_flip = [dict(row) for row in artifacts]
    sign_flip[2]["torque_sign_convention"] = "positive_generator"
    sign_gate = femm_air_gap_sample_solver_ready_manifest_gate(sign_flip)
    assert sign_gate["status"] == "needs_attention"
    assert sign_gate["checks"]["torque_sign_convention_consistent"] is False


def test_balanced_back_emf_line_voltage_gate_cancels_triplen_from_line_line():
    harmonics = {1: 20.0, 3: 7.0, 5: 1.6, 7: 0.8}
    gate = balanced_back_emf_line_voltage_handoff_gate(harmonics)

    phase_rms = math.sqrt((20.0 ** 2 + 7.0 ** 2 + 1.6 ** 2 + 0.8 ** 2) / 2.0)
    line_rms = math.sqrt((
        (math.sqrt(3.0) * 20.0) ** 2
        + (math.sqrt(3.0) * 1.6) ** 2
        + (math.sqrt(3.0) * 0.8) ** 2
    ) / 2.0)

    assert gate["policy"] == "balanced_back_emf_line_voltage_handoff_gate"
    assert gate["phase_rms_total"] == pytest.approx(phase_rms)
    assert gate["line_line_rms_total"] == pytest.approx(line_rms)
    assert gate["max_triplen_line_line_peak"] == pytest.approx(0.0, abs=1.0e-12)
    assert gate["checks"]["triplen_cancel_from_line_line"] is True

    wrong_line = balanced_back_emf_line_voltage_handoff_gate(
        harmonics,
        measured_line_line_rms=math.sqrt(3.0) * phase_rms,
    )
    assert wrong_line["status"] == "needs_attention"
    assert wrong_line["checks"]["measured_line_line_rms_ok"] is False


def test_flux_linkage_back_emf_derivative_gate_tracks_motor_pickup_rows():
    n = 96
    lambda_peak = 0.031
    omega = 220.0
    theta = [2.0 * math.pi * i / n for i in range(n)]
    flux = [lambda_peak * math.cos(angle) for angle in theta]
    back_emf = [omega * lambda_peak * math.sin(angle) for angle in theta]

    gate = flux_linkage_back_emf_derivative_gate(
        theta,
        flux,
        back_emf,
        omega,
        rtol=1.0e-3,
        atol=1.0e-10,
    )

    assert gate["status"] == "ok"
    assert gate["policy"] == "flux_linkage_back_emf_derivative_gate"
    assert gate["sign_convention"] == "negative_derivative"
    assert gate["n_samples"] == n
    assert gate["covered_period_rad"] == pytest.approx(2.0 * math.pi)
    assert gate["max_abs_error_v"] / (omega * lambda_peak) < 1.0e-3
    assert gate["checks"]["one_period_without_duplicate_endpoint"] is True

    wrong_sign = flux_linkage_back_emf_derivative_gate(
        theta,
        flux,
        [-value for value in back_emf],
        omega,
        rtol=1.0e-3,
    )
    assert wrong_sign["status"] == "needs_attention"
    assert wrong_sign["checks"]["all_rows_within_tolerance"] is False

    with pytest.raises(ValueError, match="duplicate endpoint"):
        flux_linkage_back_emf_derivative_gate(
            theta + [2.0 * math.pi],
            flux + [lambda_peak],
            back_emf + [0.0],
            omega,
        )


def test_inverter_dc_bus_voltage_limit_gate_separates_spwm_and_svpwm_margin():
    harmonics = {1: 20.0, 3: 7.0, 5: 1.6, 7: 0.8}
    emf = balanced_back_emf_line_voltage_handoff_gate(harmonics)
    line_rms = emf["line_line_rms_total"]

    spwm = inverter_dc_bus_voltage_limit_gate(
        dc_bus_v=48.0,
        modulation_index=0.92,
        method="spwm",
        measured_line_line_rms=line_rms,
    )
    low_spwm = inverter_dc_bus_voltage_limit_gate(
        dc_bus_v=48.0,
        modulation_index=0.80,
        method="spwm",
        measured_line_line_rms=line_rms,
    )
    svpwm = inverter_dc_bus_voltage_limit_gate(
        dc_bus_v=48.0,
        modulation_index=0.80,
        method="svpwm",
        measured_line_line_rms=line_rms,
    )

    assert spwm["policy"] == "inverter_dc_bus_voltage_limit_gate"
    assert spwm["method"] == "spwm"
    assert spwm["line_line_rms_factor_at_m1"] == pytest.approx(math.sqrt(3.0) / (2.0 * math.sqrt(2.0)))
    assert spwm["line_line_rms_limit"] == pytest.approx(27.042366760326285)
    assert spwm["line_line_rms_margin"] == pytest.approx(2.449684922023252)
    assert spwm["status"] == "ok"
    assert low_spwm["status"] == "needs_attention"
    assert low_spwm["checks"]["measured_line_line_rms_within_limit"] is False
    assert svpwm["status"] == "ok"
    assert svpwm["line_line_rms_limit"] == pytest.approx(27.152900397563425)


def test_dq_torque_table_health_closes_jmag_map_column_contract():
    lambda_m = 0.10
    Ld = 4.0e-3
    Lq = 9.0e-3
    current = 20.0
    pole_pairs = 3
    rows = []
    for gamma_deg in range(-60, 65, 5):
        id_current, iq_current = dq_current_from_gamma_deg(current, gamma_deg)
        rows.append({
            "gamma_deg": float(gamma_deg),
            "id_A": id_current,
            "iq_A": iq_current,
            "torque_Nm": lumped_pm_dq_torque(lambda_m, Ld, Lq, id_current, iq_current, pole_pairs),
        })

    health = dq_torque_table_health(
        rows,
        lambda_m=lambda_m,
        Ld=Ld,
        Lq=Lq,
        current=current,
        pole_pairs=pole_pairs,
    )

    assert health["status"] == "ok"
    assert health["row_count"] == 25
    assert health["max_current_abs_error_A"] < 1.0e-12
    assert health["max_torque_abs_error_Nm"] < 1.0e-12
    assert health["pure_q_row"]["gamma_deg"] == pytest.approx(0.0)
    assert health["pure_q_row"]["torque_Nm"] == pytest.approx(9.0)
    assert health["peak_row"]["gamma_deg"] == pytest.approx(30.0)
    assert health["peak_row"]["id_A"] == pytest.approx(-10.0)
    assert health["peak_row"]["iq_A"] == pytest.approx(17.320508075688775)
    assert health["peak_row"]["torque_Nm"] == pytest.approx(11.691342951089922)
    assert health["peak_row"]["torque_Nm"] / health["pure_q_row"]["torque_Nm"] - 1.0 == pytest.approx(
        0.299038105676658
    )

    wrong_peak = [dict(row) for row in rows]
    for row in wrong_peak:
        if row["gamma_deg"] == 30.0:
            row["torque_Nm"] = row["torque_Nm"] - 4.0
    bad = dq_torque_table_health(wrong_peak, lambda_m, Ld, Lq, current, pole_pairs)
    assert bad["status"] == "needs_attention"
    assert bad["checks"]["torque_column_ok"] is False
    assert bad["checks"]["peak_row_matches_closed_form"] is False


def test_jmag_pm_short_circuit_fault_table_gate_tracks_dq_fault_physics():
    R = 0.05
    Ld = 8.0e-3
    Lq = 16.0e-3
    lambda_m = 0.10
    pole_pairs = 4
    ich = lambda_m / Ld

    def row(omega_e):
        den = R * R + omega_e * omega_e * Ld * Lq
        id_current = -omega_e * omega_e * Lq * lambda_m / den
        iq_current = -omega_e * R * lambda_m / den
        torque = lumped_pm_dq_torque(lambda_m, Ld, Lq, id_current, iq_current, pole_pairs)
        vd = R * id_current - omega_e * Lq * iq_current
        vq = R * iq_current + omega_e * (Ld * id_current + lambda_m)
        return {
            "omega_e": omega_e,
            "omega_mech": omega_e / pole_pairs,
            "id_A": id_current,
            "iq_A": iq_current,
            "torque_Nm": torque,
            "vd_residual": vd,
            "vq_residual": vq,
            "current_ratio_to_characteristic": math.hypot(id_current, iq_current) / ich,
            "d_axis_demag_fraction": -id_current / ich,
        }

    rows = [row(1.0), row(6.25), row(50.0), row(5000.0)]
    gate = jmag_pm_short_circuit_fault_table_gate(rows, R, Ld, Lq, lambda_m, pole_pairs)

    assert gate["policy"] == "jmag_pm_short_circuit_fault_table_gate"
    assert gate["status"] == "ok"
    assert gate["characteristic_current_A"] == pytest.approx(12.5)
    assert gate["peak_braking_row"]["omega_e"] == pytest.approx(6.25)
    assert gate["high_speed_row"]["d_axis_demag_fraction"] == pytest.approx(0.9999992187506104)
    assert gate["checks"]["short_terminal_residuals_ok"] is True
    assert gate["checks"]["braking_torque_peaks_at_intermediate_speed"] is True
    assert gate["checks"]["high_speed_demag_fraction_near_one"] is True

    bad_residual = [dict(item) for item in rows]
    bad_residual[1]["vd_residual"] = 1.0e-3
    residual_gate = jmag_pm_short_circuit_fault_table_gate(bad_residual, R, Ld, Lq, lambda_m, pole_pairs)
    assert residual_gate["status"] == "needs_attention"
    assert residual_gate["checks"]["short_terminal_residuals_ok"] is False

    stale_speed = [dict(item) for item in rows]
    stale_speed[2]["omega_mech"] = stale_speed[2]["omega_e"]
    speed_gate = jmag_pm_short_circuit_fault_table_gate(stale_speed, R, Ld, Lq, lambda_m, pole_pairs)
    assert speed_gate["status"] == "needs_attention"
    assert speed_gate["checks"]["omega_mech_matches_pole_pairs"] is False

    wrong_demag = [dict(item) for item in rows]
    wrong_demag[-1]["d_axis_demag_fraction"] = 0.20
    demag_gate = jmag_pm_short_circuit_fault_table_gate(wrong_demag, R, Ld, Lq, lambda_m, pole_pairs)
    assert demag_gate["status"] == "needs_attention"
    assert demag_gate["checks"]["d_axis_demag_fraction_matches_characteristic_current"] is False
    assert demag_gate["checks"]["high_speed_demag_fraction_near_one"] is False

    wrong_torque = [dict(item) for item in rows]
    wrong_torque[1]["torque_Nm"] = abs(wrong_torque[1]["torque_Nm"])
    torque_gate = jmag_pm_short_circuit_fault_table_gate(wrong_torque, R, Ld, Lq, lambda_m, pole_pairs)
    assert torque_gate["status"] == "needs_attention"
    assert torque_gate["checks"]["torque_column_matches_closed_form"] is False


def test_pm_drive_terminal_table_health_allows_r_included_voltage_report_over_one():
    pole_pairs = 2
    vmax = 6.0
    rows = []
    for omega_e, id_current, iq_current, region in (
        (50.0, 0.0, 10.0, "MTPA"),
        (100.0, 0.0, 10.0, "FW"),
    ):
        omega_mech = omega_e / pole_pairs
        vd = -omega_e * 0.002 * iq_current
        vq = 0.1 * iq_current + omega_e * 0.05
        vmag = math.hypot(vd, vq)
        imag = math.hypot(id_current, iq_current)
        torque = 1.5
        p_em = torque * omega_mech
        p_cu = 1.5 * 0.1 * imag * imag
        p_in = p_em + p_cu
        rows.append({
            "omega_e": omega_e,
            "omega_mech": omega_mech,
            "speed_multiple": omega_e / 100.0,
            "region": region,
            "id_A": id_current,
            "iq_A": iq_current,
            "vd_V": vd,
            "vq_V": vq,
            "Vmag_V": vmag,
            "Imag_A": imag,
            "torque_Nm": torque,
            "P_in_W": p_in,
            "P_em_W": p_em,
            "P_cu_W": p_cu,
            "power_factor": p_in / (1.5 * vmag * imag),
            "efficiency": p_em / p_in,
            "Vmax_V": vmax,
            "voltage_utilization": vmag / vmax,
            "voltage_utilization_lossless": omega_e * math.hypot(0.05, 0.002 * iq_current) / vmax,
        })

    health = pm_drive_terminal_table_health(rows, pole_pairs=pole_pairs)

    assert health["status"] == "ok"
    assert health["policy"] == "pm_drive_terminal_table_health_gate"
    assert health["max_speed_contract_rel_error"] == pytest.approx(0.0)
    assert health["max_power_balance_rel_error"] == pytest.approx(0.0)
    assert health["max_torque_speed_rel_error"] == pytest.approx(0.0)
    assert health["max_power_factor_abs_error"] == pytest.approx(0.0)
    assert health["max_lossless_voltage_utilization"] < 1.0
    assert health["max_terminal_voltage_utilization"] > 1.0
    assert health["terminal_voltage_over_limit_row_count"] == 1
    assert health["checks"]["lossless_voltage_constraint_ok"] is True


def test_jmag_terminal_table_pairs_with_dc_bus_voltage_limit_handoff():
    pole_pairs = 2
    rows = []
    for omega_e, id_current, iq_current, region in (
        (50.0, 0.0, 10.0, "MTPA"),
        (100.0, 0.0, 10.0, "FW"),
    ):
        omega_mech = omega_e / pole_pairs
        vd = -omega_e * 0.002 * iq_current
        vq = 0.1 * iq_current + omega_e * 0.05
        vmag = math.hypot(vd, vq)
        imag = math.hypot(id_current, iq_current)
        torque = 1.5
        p_em = torque * omega_mech
        p_cu = 1.5 * 0.1 * imag * imag
        p_in = p_em + p_cu
        rows.append({
            "omega_e": omega_e,
            "omega_mech": omega_mech,
            "speed_multiple": omega_e / 100.0,
            "region": region,
            "id_A": id_current,
            "iq_A": iq_current,
            "vd_V": vd,
            "vq_V": vq,
            "Vmag_V": vmag,
            "Imag_A": imag,
            "torque_Nm": torque,
            "P_in_W": p_in,
            "P_em_W": p_em,
            "P_cu_W": p_cu,
            "power_factor": p_in / (1.5 * vmag * imag),
            "efficiency": p_em / p_in,
            "Vmax_V": 6.0,
            "voltage_utilization": vmag / 6.0,
            "voltage_utilization_lossless": omega_e * math.hypot(0.05, 0.002 * iq_current) / 6.0,
        })
    terminal = pm_drive_terminal_table_health(rows, pole_pairs=pole_pairs)
    emf = balanced_back_emf_line_voltage_handoff_gate({1: 20.0, 3: 7.0, 5: 1.6, 7: 0.8})
    dc_bus = inverter_dc_bus_voltage_limit_gate(
        dc_bus_v=48.0,
        modulation_index=0.92,
        method="spwm",
        measured_line_line_rms=emf["line_line_rms_total"],
    )

    assert terminal["status"] == "ok"
    assert terminal["terminal_voltage_over_limit_policy"] == (
        "report_only_use_lossless_voltage_utilization_for_selector_feasibility"
    )
    assert terminal["max_terminal_voltage_utilization"] > 1.0
    assert terminal["max_lossless_voltage_utilization"] < 1.0
    assert dc_bus["status"] == "ok"
    assert dc_bus["line_line_rms_utilization"] == pytest.approx(0.9094130506342708)
    assert dc_bus["line_line_rms_margin"] == pytest.approx(2.449684922023252)


def test_efficiency_map_cell_voltage_margin_matches_dc_bus_gate():
    emf = balanced_back_emf_line_voltage_handoff_gate({1: 20.0, 3: 7.0, 5: 1.6, 7: 0.8})
    dc_bus = inverter_dc_bus_voltage_limit_gate(
        dc_bus_v=48.0,
        modulation_index=0.92,
        method="spwm",
        measured_line_line_rms=emf["line_line_rms_total"],
    )
    map_cell = {
        "point_id": "s00_t00",
        "speed_rpm": 1000.0,
        "torque_Nm": 0.2,
        "efficiency": 0.91,
        "total_loss_W": 2.0,
        "voltage_margin_v": dc_bus["line_line_rms_margin"],
    }

    assert dc_bus["status"] == "ok"
    assert map_cell["voltage_margin_v"] == pytest.approx(2.449684922023252)
    assert map_cell["voltage_margin_v"] == pytest.approx(dc_bus["line_line_rms_margin"])
    assert map_cell["voltage_margin_v"] > 0.0


def test_pm_drive_loss_bucket_efficiency_gate_checks_power_bookkeeping():
    rows = []
    for operating_point, omega_mech, torque, p_cu, p_iron, p_magnet, p_mech in (
        ("MTPA", 100.0, 12.0, 60.0, 35.0, 8.0, 12.0),
        ("FW", 200.0, 12.0, 120.0, 55.0, 15.0, 18.0),
        ("high_current", 150.0, 12.0, 260.0, 80.0, 22.0, 20.0),
    ):
        p_out = omega_mech * torque
        p_in = p_out + p_cu + p_iron + p_magnet + p_mech
        rows.append({
            "operating_point": operating_point,
            "omega_mech": omega_mech,
            "torque_Nm": torque,
            "P_out_W": p_out,
            "P_in_W": p_in,
            "P_cu_W": p_cu,
            "P_iron_W": p_iron,
            "P_magnet_W": p_magnet,
            "P_mechanical_loss_W": p_mech,
            "efficiency": p_out / p_in,
        })

    gate = pm_drive_loss_bucket_efficiency_gate(rows)

    assert gate["status"] == "ok"
    assert gate["policy"] == "pm_drive_loss_bucket_efficiency_gate"
    assert gate["row_count"] == 3
    assert gate["max_power_balance_rel_error"] == pytest.approx(0.0)
    assert gate["max_efficiency_abs_error"] == pytest.approx(0.0)
    assert gate["max_torque_speed_rel_error"] == pytest.approx(0.0)
    assert gate["max_efficiency_row"]["operating_point"] == "FW"
    assert gate["max_efficiency"] == pytest.approx(2400.0 / 2608.0)
    assert gate["max_loss_fraction_row"]["operating_point"] == "high_current"
    assert gate["max_loss_fraction_row"]["dominant_loss_bucket"] == "P_cu_W"

    omitted_magnet_loss = [dict(row) for row in rows]
    omitted_magnet_loss[1]["P_in_W"] -= omitted_magnet_loss[1]["P_magnet_W"]
    bad = pm_drive_loss_bucket_efficiency_gate(omitted_magnet_loss)
    assert bad["status"] == "needs_attention"
    assert bad["checks"]["power_balance_ok"] is False
    assert bad["checks"]["efficiency_formula_ok"] is False

    elf_style_rows = []
    for row in rows:
        elf_style_rows.append({
            "operating_point": row["operating_point"],
            "omega_mech": row["omega_mech"],
            "torque_Nm": row["torque_Nm"],
            "P_out_W": row["P_out_W"],
            "P_in_W": row["P_in_W"],
            "copper_loss_w": row["P_cu_W"],
            "iron_loss_w": row["P_iron_W"],
            "magnet_loss_w": row["P_magnet_W"],
            "mechanical_loss_w": row["P_mechanical_loss_W"],
            "efficiency": row["efficiency"],
        })
    elf_gate = pm_drive_loss_bucket_efficiency_gate(elf_style_rows)
    assert elf_gate["status"] == "ok"
    assert elf_gate["loss_columns"] == [
        "copper_loss_w",
        "iron_loss_w",
        "magnet_loss_w",
        "mechanical_loss_w",
    ]
    assert elf_gate["max_efficiency"] == pytest.approx(gate["max_efficiency"])


def test_drive_cycle_weighted_efficiency_gate_scores_elf_style_rows():
    rows = [
        {
            "point_id": "hover",
            "weight": 50.0,
            "P_out_W": 100.0,
            "P_in_W": 112.0,
            "copper_loss_w": 7.0,
            "iron_loss_w": 3.0,
            "magnet_loss_w": 1.0,
            "mechanical_loss_w": 1.0,
            "efficiency": 100.0 / 112.0,
        },
        {
            "point_id": "climb",
            "weight": 25.0,
            "P_out_W": 220.0,
            "P_in_W": 250.0,
            "copper_loss_w": 18.0,
            "iron_loss_w": 7.0,
            "magnet_loss_w": 3.0,
            "mechanical_loss_w": 2.0,
            "efficiency": 220.0 / 250.0,
        },
        {
            "point_id": "burst",
            "weight": 25.0,
            "P_out_W": 360.0,
            "P_in_W": 420.0,
            "copper_loss_w": 35.0,
            "iron_loss_w": 15.0,
            "magnet_loss_w": 6.0,
            "mechanical_loss_w": 4.0,
            "efficiency": 360.0 / 420.0,
        },
    ]

    gate = drive_cycle_weighted_efficiency_gate(rows)

    assert gate["status"] == "ok"
    assert gate["policy"] == "drive_cycle_weighted_efficiency_gate"
    assert gate["input_weight_sum"] == pytest.approx(100.0)
    assert gate["weighted_output_W"] == pytest.approx(195.0)
    assert gate["weighted_input_W"] == pytest.approx(223.5)
    assert gate["weighted_total_loss_W"] == pytest.approx(28.5)
    assert gate["cycle_efficiency"] == pytest.approx(195.0 / 223.5)
    assert gate["dominant_weighted_loss_bucket"] == "P_cu_W"
    assert gate["weighted_losses_W"]["P_cu_W"] == pytest.approx(16.75)
    assert gate["worst_efficiency_row"]["point_id"] == "burst"

    bad = [dict(row) for row in rows]
    bad[1]["P_in_W"] -= 3.0
    bad_gate = drive_cycle_weighted_efficiency_gate(bad)
    assert bad_gate["status"] == "needs_attention"
    assert bad_gate["checks"]["power_balance_ok"] is False


def test_pm_recoil_demag_step_summary_tracks_irreversible_crossing():
    safe = pm_recoil_demag_step_summary(
        {0: -2.0e5, 1: -4.0e5, 2: -2.2e5},
        H_knee_A_per_m=-5.0e5,
    )
    assert safe["status"] == "ok"
    assert safe["irreversible_demag"] is False
    assert safe["recoil_remanence_ratio_proxy"] == pytest.approx(1.0)

    crossed = pm_recoil_demag_step_summary(
        {0: -2.0e5, 1: -6.0e5, 2: -3.0e5},
        H_knee_A_per_m=-5.0e5,
    )
    assert crossed["status"] == "ok"
    assert crossed["irreversible_demag"] is True
    assert crossed["margins_A_per_m"]["1"] == pytest.approx(-1.0e5)
    assert crossed["recoil_remanence_ratio_proxy"] == pytest.approx(0.8)

    overclaimed_recovery = pm_recoil_demag_step_summary(
        {0: -2.0e5, 1: -6.0e5, 2: -1.0e5},
        H_knee_A_per_m=-5.0e5,
    )
    assert overclaimed_recovery["status"] == "needs_attention"
    assert overclaimed_recovery["checks"]["step2_not_stronger_than_nominal_after_crossing"] is False


def test_pm_recoil_demag_step_summary_matches_elf_step_contract_terms():
    gate = pm_recoil_demag_step_summary(
        {0: -2.0e5, 1: -6.0e5, 2: -3.0e5},
        H_knee_A_per_m=-5.0e5,
    )

    assert gate["policy"] == "pm_recoil_demag_three_step_gate"
    assert list(gate["step_fields_A_per_m"].keys()) == ["0", "1", "2"]
    assert gate["checks"]["three_steps_present"] is True
    assert gate["checks"]["step1_is_worst_field"] is True
    assert gate["irreversible_demag"] is True
    assert gate["recoil_remanence_ratio_proxy"] == pytest.approx(0.8)


def test_pm_loadline_metadata_gate_closes_units_and_sign_before_values():
    metadata = {
        "columns": ["temperature_C", "B_gap_T", "H_pm_A_per_m", "H_knee_A_per_m"],
        "h_field_unit": "A/m",
        "b_flux_density_unit": "T",
        "temperature_unit": "C",
        "field_sign_convention": "negative_is_demag",
        "magnetization_axis": "radial",
        "knee_reference": "intrinsic_bh_knee",
        "recoil_mu_r": 1.05,
    }

    gate = pm_loadline_metadata_gate(metadata)

    assert gate["status"] == "ok"
    assert gate["policy"] == "pm_loadline_metadata_gate"
    assert gate["h_field_unit"] == "a_per_m"
    assert gate["field_sign_convention"] == "negative_is_demag"
    assert all(gate["checks"].values())

    bad = pm_loadline_metadata_gate({
        "columns": ["temperature", "B", "H"],
        "h_field_unit": "Oe",
        "b_flux_density_unit": "Gauss",
        "temperature_unit": "F",
        "field_sign_convention": "toward_magnet",
        "magnetization_axis": "north",
        "knee_reference": "unknown",
        "recoil_mu_r": 0.0,
    })

    assert bad["status"] == "needs_attention"
    assert bad["missing_required_columns"] == [
        "temperature_C",
        "B_gap_T",
        "H_pm_A_per_m",
        "H_knee_A_per_m",
    ]
    assert bad["checks"] == {
        "required_columns_present": False,
        "h_field_unit_valid": False,
        "b_flux_density_unit_valid": False,
        "temperature_unit_valid": False,
        "field_sign_convention_valid": False,
        "magnetization_axis_valid": False,
        "knee_reference_valid": False,
        "recoil_mu_r_positive": False,
    }


def test_pm_bem_surface_normal_metadata_gate_balances_closed_pm_surface_charge():
    identity = {
        "surface_source_artifact_id": "elf_slot278_pm_A_surface_sources_v1.json",
        "magnetization_source_id": "elf_slot278_pm_A_magnetization_map_v1.json",
        "material_id": "pm_A",
        "material_name": "NdFeB_recoil_A",
    }
    rows = [
        {"surface": "top", "area_m2": 2.0, "normal": [0.0, 0.0, 1.0], "magnetization": [0.0, 0.0, 1.0], "normal_convention": "outward_from_magnet", **identity},
        {"surface": "bottom", "area_m2": 2.0, "normal": [0.0, 0.0, -1.0], "magnetization": [0.0, 0.0, 1.0], "normal_convention": "outward_from_magnet", **identity},
        {"surface": "side_xp", "area_m2": 1.0, "normal": [1.0, 0.0, 0.0], "magnetization": [0.0, 0.0, 1.0], "normal_convention": "outward_from_magnet", **identity},
        {"surface": "side_xm", "area_m2": 1.0, "normal": [-1.0, 0.0, 0.0], "magnetization": [0.0, 0.0, 1.0], "normal_convention": "outward_from_magnet", **identity},
    ]

    gate = pm_bem_surface_normal_metadata_gate(
        rows,
        expected_surface_source_artifact_id="elf_slot278_pm_A_surface_sources_v1.json",
        expected_magnetization_source_id="elf_slot278_pm_A_magnetization_map_v1.json",
        expected_material_id="pm_A",
        expected_material_name="NdFeB_recoil_A",
    )

    assert gate["status"] == "ok"
    assert gate["policy"] == "pm_bem_surface_normal_metadata_gate"
    assert gate["signed_charge_proxy_sum"] == pytest.approx(0.0, abs=1.0e-15)
    by_surface = {row["surface"]: row for row in gate["rows"]}
    assert by_surface["top"]["m_dot_n"] == pytest.approx(1.0)
    assert by_surface["bottom"]["m_dot_n"] == pytest.approx(-1.0)
    assert by_surface["side_xp"]["m_dot_n"] == pytest.approx(0.0)
    assert by_surface["top"]["surface_source_artifact_id"] == "elf_slot278_pm_A_surface_sources_v1.json"
    assert by_surface["top"]["magnetization_source_id"] == "elf_slot278_pm_A_magnetization_map_v1.json"
    assert by_surface["top"]["material_id"] == "pm_A"
    assert by_surface["top"]["material_name"] == "NdFeB_recoil_A"
    assert gate["observed_surface_source_artifact_ids"] == ["elf_slot278_pm_A_surface_sources_v1.json"]
    assert gate["checks"]["surface_source_artifact_id_matches_expected"] is True
    assert gate["checks"]["magnetization_source_id_matches_expected"] is True
    assert gate["checks"]["material_id_matches_expected"] is True
    assert gate["checks"]["material_name_matches_expected"] is True
    assert "magnetic-charge BEM assembly" in gate["version_note"]
    assert "material identity" in gate["version_note"]

    bad = pm_bem_surface_normal_metadata_gate([
        {"surface": "top", "area_m2": 2.0, "normal": [0.0, 0.0, 2.0], "magnetization": [0.0, 0.0, 1.0], "normal_convention": "outward_from_magnet"},
        {"surface": "top", "area_m2": 2.0, "normal": [0.0, 0.0, 1.0], "normal_orientation": "inward", "magnetization": [0.0, 0.0, 1.0]},
        {"surface": "side", "area_m2": 0.0, "normal": [1.0, 0.0, 0.0], "normal_convention": "outward_from_magnet"},
    ])

    assert bad["status"] == "needs_attention"
    assert bad["duplicate_surfaces"] == ["top"]
    assert bad["bad_normal_surfaces"] == ["top"]
    assert bad["missing_or_wrong_convention_surfaces"] == ["top"]
    assert bad["nonpositive_area_surfaces"] == ["side"]
    assert bad["missing_magnetization_surfaces"] == ["side"]
    assert bad["checks"]["closed_surface_charge_balances"] is False

    stale_source = [dict(row) for row in rows]
    stale_source[1]["surface_source_artifact_id"] = "elf_slot278_old_surface_sources.json"
    stale_source_gate = pm_bem_surface_normal_metadata_gate(
        stale_source,
        expected_surface_source_artifact_id="elf_slot278_pm_A_surface_sources_v1.json",
        expected_magnetization_source_id="elf_slot278_pm_A_magnetization_map_v1.json",
        expected_material_id="pm_A",
        expected_material_name="NdFeB_recoil_A",
    )
    assert stale_source_gate["status"] == "needs_attention"
    assert stale_source_gate["checks"]["closed_surface_charge_balances"] is True
    assert stale_source_gate["checks"]["surface_source_artifact_id_matches_expected"] is False
    assert stale_source_gate["mismatched_identity_surfaces"]["surface_source_artifact_id"] == ["bottom"]

    stale_material = [dict(row) for row in rows]
    stale_material[2]["material_name"] = "Air"
    stale_material_gate = pm_bem_surface_normal_metadata_gate(
        stale_material,
        expected_surface_source_artifact_id="elf_slot278_pm_A_surface_sources_v1.json",
        expected_magnetization_source_id="elf_slot278_pm_A_magnetization_map_v1.json",
        expected_material_id="pm_A",
        expected_material_name="NdFeB_recoil_A",
    )
    assert stale_material_gate["status"] == "needs_attention"
    assert stale_material_gate["checks"]["closed_surface_charge_balances"] is True
    assert stale_material_gate["checks"]["material_name_matches_expected"] is False
    assert stale_material_gate["mismatched_identity_surfaces"]["material_name"] == ["side_xp"]


def test_pm_demag_package_identity_gate_bundles_run_result_loadline_bem_and_recoil():
    artifacts = [
        {
            "kind": "run_result",
            "case_id": "demag_case_007",
            "magnet_id": "pm_A",
            "path": "slot150_run_result.json",
            "gate_policy": "elf_python_run_result_parse_path",
            "status": "ok",
            "normalized_columns": ["case_id", "H_pm_A_per_m", "H_knee_A_per_m", "safe_against_knee"],
        },
        {
            "kind": "loadline_metadata",
            "case_id": "demag_case_007",
            "magnet_id": "pm_A",
            "path": "slot150_loadline_metadata.json",
            "gate_policy": "pm_loadline_metadata_gate",
            "status": "ok",
        },
        {
            "kind": "bem_surface",
            "case_id": "demag_case_007",
            "magnet_id": "pm_A",
            "path": "slot150_bem_surface.json",
            "gate_policy": "pm_bem_surface_normal_metadata_gate",
            "status": "ok",
        },
        {
            "kind": "recoil_steps",
            "case_id": "demag_case_007",
            "magnet_id": "pm_A",
            "path": "slot150_recoil_steps.json",
            "gate_policy": "pm_recoil_demag_three_step_gate",
            "status": "ok",
            "steps": [0, 1, 2],
        },
    ]

    gate = pm_demag_package_identity_gate(
        artifacts,
        expected_case_id="demag_case_007",
        expected_magnet_id="pm_A",
    )

    assert gate["status"] == "ok"
    assert gate["policy"] == "pm_demag_package_identity_gate"
    assert gate["case_ids"] == ["demag_case_007"]
    assert gate["magnet_ids"] == ["pm_A"]
    assert gate["checks"]["run_result_has_loadline_columns"] is True
    assert gate["checks"]["recoil_steps_are_three_step"] is True

    stale_magnet = [dict(row) for row in artifacts]
    stale_magnet[2]["magnet_id"] = "pm_B"
    stale_magnet_gate = pm_demag_package_identity_gate(stale_magnet)
    assert stale_magnet_gate["status"] == "needs_attention"
    assert stale_magnet_gate["checks"]["magnet_ids_unique"] is False

    missing_case = [dict(row) for row in artifacts]
    missing_case[0].pop("case_id")
    missing_case_gate = pm_demag_package_identity_gate(missing_case)
    assert missing_case_gate["status"] == "needs_attention"
    assert missing_case_gate["checks"]["case_ids_present"] is False

    missing_h = [dict(row) for row in artifacts]
    missing_h[0]["normalized_columns"] = ["case_id", "H_pm_A_per_m"]
    missing_h_gate = pm_demag_package_identity_gate(missing_h)
    assert missing_h_gate["status"] == "needs_attention"
    assert missing_h_gate["checks"]["run_result_has_loadline_columns"] is False

    bad_steps = [dict(row) for row in artifacts]
    bad_steps[3]["steps"] = [0, 1]
    bad_steps_gate = pm_demag_package_identity_gate(bad_steps)
    assert bad_steps_gate["status"] == "needs_attention"
    assert bad_steps_gate["checks"]["recoil_steps_are_three_step"] is False


def test_pm_demag_margin_screening_package_gate_bundles_sweep_fault_bem_and_package():
    artifacts = [
        {
            "kind": "loadline_sweep",
            "case_id": "demag_case_007",
            "magnet_id": "pm_A",
            "temperature_C": 120,
            "path": "slot158_loadline_sweep.json",
            "gate_policy": "pm_temperature_demag_sweep_summary",
            "status": "ok",
            "first_unsafe_gap_m": 0.008,
            "minimum_demag_margin_A_per_m": -12500.0,
        },
        {
            "kind": "fault_current_screening",
            "case_id": "demag_case_007",
            "magnet_id": "pm_A",
            "temperature_C": 120,
            "path": "slot158_fault_current_screening.json",
            "gate_policy": "fault_current_demag_screening",
            "status": "ok",
            "negative_id_is_demag_direction": True,
            "recommended_observable_keys": [
                "field_probe",
                "demag_margin_A_per_m",
                "recoil_remanence_ratio_proxy",
            ],
            "field_probe_id": "fault_probe_pm_A_demag_axis_v1",
            "field_probe_family": "demag_margin_field_probe",
            "observation_region_id": "pm_A_volume",
            "observation_component": "H_parallel_demag_axis",
            "field_axis_convention": "magnetization_axis_positive",
            "field_sign_convention": "negative_h_parallel_is_demag",
            "averaging_rule": "volume_average",
            "field_probe_output_artifact_id": "fault_probe_pm_A_table_v1",
            "field_probe_output_digest": "sha256:fault_probe_pm_A_table_v1",
            "field_probe_output_path": "slot158_fault_probe_table.json",
        },
        {
            "kind": "bem_source_balance",
            "case_id": "demag_case_007",
            "magnet_id": "pm_A",
            "temperature_C": 120,
            "path": "slot158_bem_source_balance.json",
            "gate_policy": "pm_bem_surface_source_balance_gate",
            "status": "ok",
            "surface_mesh_id": "bem_surface_mesh_pm_A_v1",
            "surface_mesh_digest": "sha256:bem-surface-mesh-pm-A-v1",
            "surface_row_count": 6,
            "source_balance_artifact_id": "bem_source_balance_pm_A_v1",
            "signed_charge_balance_rel_tol": 1.0e-9,
            "signed_charge_balance_abs": 0.0,
            "total_area_m2": 0.024,
            "normal_convention": "outward_from_magnet",
            "source_balance_unit": "area_weighted_m_dot_n",
        },
        {
            "kind": "demag_package",
            "case_id": "demag_case_007",
            "magnet_id": "pm_A",
            "temperature_C": 120,
            "path": "slot158_demag_package.json",
            "gate_policy": "pm_demag_package_identity_gate",
            "status": "ok",
            "required_artifacts": [
                "run_result",
                "loadline_metadata",
                "bem_surface",
                "recoil_steps",
            ],
        },
    ]

    gate = pm_demag_margin_screening_package_gate(
        artifacts,
        expected_case_id="demag_case_007",
        expected_magnet_id="pm_A",
        expected_temperature_c=120,
        expected_bem_surface_mesh_id="bem_surface_mesh_pm_A_v1",
        expected_bem_surface_mesh_digest="sha256:bem-surface-mesh-pm-A-v1",
        expected_bem_surface_row_count=6,
        expected_field_probe_id="fault_probe_pm_A_demag_axis_v1",
        expected_field_probe_family="demag_margin_field_probe",
        expected_observation_region_id="pm_A_volume",
        expected_observation_component="h_parallel_demag_axis",
        expected_field_axis_convention="magnetization_axis_positive",
        expected_field_sign_convention="negative_h_parallel_is_demag",
        expected_averaging_rule="volume_average",
        expected_field_probe_output_artifact_id="fault_probe_pm_A_table_v1",
        expected_field_probe_output_digest="sha256:fault_probe_pm_A_table_v1",
        require_field_probe_identity=True,
        require_field_probe_output_artifact=True,
    )

    assert gate["status"] == "ok"
    assert gate["policy"] == "pm_demag_margin_screening_package_gate"
    assert gate["case_ids"] == ["demag_case_007"]
    assert gate["magnet_ids"] == ["pm_A"]
    assert gate["temperatures_C"] == [120.0]
    assert gate["checks"]["fault_current_observables_complete"] is True
    assert gate["checks"]["bem_surface_area_recorded"] is True
    assert gate["checks"]["bem_normal_convention_recorded"] is True
    assert gate["checks"]["bem_normal_convention_matches_expected"] is True
    assert gate["expected_bem_normal_convention"] == "outward_from_magnet"
    assert gate["checks"]["bem_source_unit_recorded"] is True
    assert gate["checks"]["bem_balance_value_recorded"] is True
    assert gate["checks"]["bem_surface_mesh_id_recorded"] is True
    assert gate["checks"]["bem_surface_mesh_digest_recorded"] is True
    assert gate["checks"]["bem_surface_row_count_recorded"] is True
    assert gate["checks"]["expected_bem_surface_mesh_id_matches"] is True
    assert gate["checks"]["expected_bem_surface_mesh_digest_matches"] is True
    assert gate["checks"]["expected_bem_surface_row_count_matches"] is True
    assert gate["bem_surface_mesh_ids"] == ["bem_surface_mesh_pm_A_v1"]
    assert gate["bem_surface_mesh_digests"] == ["sha256:bem-surface-mesh-pm-A-v1"]
    assert gate["bem_surface_row_counts"] == [6]
    assert gate["checks"]["bem_source_balance_artifact_id_recorded"] is True
    assert gate["checks"]["demag_package_artifacts_complete"] is True
    assert gate["field_probe_ids"] == ["fault_probe_pm_A_demag_axis_v1"]
    assert gate["field_probe_families"] == ["demag_margin_field_probe"]
    assert gate["observation_region_ids"] == ["pm_A_volume"]
    assert gate["observation_components"] == ["h_parallel_demag_axis"]
    assert gate["field_axis_conventions"] == ["magnetization_axis_positive"]
    assert gate["field_sign_conventions"] == ["negative_h_parallel_is_demag"]
    assert gate["averaging_rules"] == ["volume_average"]
    assert gate["field_probe_output_artifact_ids"] == ["fault_probe_pm_A_table_v1"]
    assert gate["field_probe_output_digests"] == ["sha256:fault_probe_pm_A_table_v1"]
    assert gate["field_probe_output_paths"] == ["slot158_fault_probe_table.json"]
    assert gate["checks"]["field_probe_id_recorded"] is True
    assert gate["checks"]["field_probe_family_recorded"] is True
    assert gate["checks"]["observation_region_id_recorded"] is True
    assert gate["checks"]["observation_component_recorded"] is True
    assert gate["checks"]["field_axis_convention_recorded"] is True
    assert gate["checks"]["field_sign_convention_recorded"] is True
    assert gate["checks"]["averaging_rule_recorded"] is True
    assert gate["checks"]["field_probe_output_artifact_id_recorded"] is True
    assert gate["checks"]["field_probe_output_digest_recorded"] is True
    assert gate["checks"]["field_probe_output_path_recorded"] is True
    assert gate["checks"]["expected_field_probe_id_matches"] is True
    assert gate["checks"]["expected_field_probe_family_matches"] is True
    assert gate["checks"]["expected_observation_region_id_matches"] is True
    assert gate["checks"]["expected_observation_component_matches"] is True
    assert gate["checks"]["expected_field_axis_convention_matches"] is True
    assert gate["checks"]["expected_field_sign_convention_matches"] is True
    assert gate["checks"]["expected_averaging_rule_matches"] is True
    assert gate["checks"]["expected_field_probe_output_artifact_id_matches"] is True
    assert gate["checks"]["expected_field_probe_output_digest_matches"] is True

    stale_temperature = [dict(row) for row in artifacts]
    stale_temperature[1]["temperature_C"] = 80
    stale_temperature_gate = pm_demag_margin_screening_package_gate(stale_temperature)
    assert stale_temperature_gate["status"] == "needs_attention"
    assert stale_temperature_gate["checks"]["temperatures_unique"] is False

    missing_fault_direction = [dict(row) for row in artifacts]
    missing_fault_direction[1]["negative_id_is_demag_direction"] = False
    missing_fault_gate = pm_demag_margin_screening_package_gate(missing_fault_direction)
    assert missing_fault_gate["status"] == "needs_attention"
    assert missing_fault_gate["checks"]["fault_current_demag_direction_recorded"] is False

    missing_bem_tolerance = [dict(row) for row in artifacts]
    missing_bem_tolerance[2].pop("signed_charge_balance_rel_tol")
    missing_bem_gate = pm_demag_margin_screening_package_gate(missing_bem_tolerance)
    assert missing_bem_gate["status"] == "needs_attention"
    assert missing_bem_gate["checks"]["bem_balance_tolerance_recorded"] is False

    missing_bem_metadata = [dict(row) for row in artifacts]
    missing_bem_metadata[2].pop("total_area_m2")
    missing_bem_metadata[2].pop("source_balance_unit")
    missing_bem_metadata_gate = pm_demag_margin_screening_package_gate(missing_bem_metadata)
    assert missing_bem_metadata_gate["status"] == "needs_attention"
    assert missing_bem_metadata_gate["checks"]["bem_surface_area_recorded"] is False
    assert missing_bem_metadata_gate["checks"]["bem_source_unit_recorded"] is False

    wrong_bem_normal = [dict(row) for row in artifacts]
    wrong_bem_normal[2]["normal_convention"] = "inward_to_magnet"
    wrong_bem_normal_gate = pm_demag_margin_screening_package_gate(wrong_bem_normal)
    assert wrong_bem_normal_gate["status"] == "needs_attention"
    assert wrong_bem_normal_gate["checks"]["bem_normal_convention_recorded"] is True
    assert wrong_bem_normal_gate["checks"]["bem_normal_convention_matches_expected"] is False
    assert wrong_bem_normal_gate["wrong_bem_normal_convention_rows"][0]["normal_convention"] == "inward_to_magnet"

    missing_bem_identity = [dict(row) for row in artifacts]
    missing_bem_identity[2].pop("surface_mesh_id")
    missing_bem_identity[2].pop("source_balance_artifact_id")
    missing_bem_identity_gate = pm_demag_margin_screening_package_gate(missing_bem_identity)
    assert missing_bem_identity_gate["status"] == "needs_attention"
    assert missing_bem_identity_gate["checks"]["bem_surface_mesh_id_recorded"] is False
    assert missing_bem_identity_gate["checks"]["bem_source_balance_artifact_id_recorded"] is False

    stale_bem_mesh = [dict(row) for row in artifacts]
    stale_bem_mesh[2]["surface_mesh_digest"] = "sha256:old-bem-surface-mesh"
    stale_bem_mesh_gate = pm_demag_margin_screening_package_gate(
        stale_bem_mesh,
        expected_bem_surface_mesh_id="bem_surface_mesh_pm_A_v1",
        expected_bem_surface_mesh_digest="sha256:bem-surface-mesh-pm-A-v1",
        expected_bem_surface_row_count=6,
    )
    assert stale_bem_mesh_gate["status"] == "needs_attention"
    assert stale_bem_mesh_gate["checks"]["expected_bem_surface_mesh_digest_matches"] is False

    wrong_bem_row_count = [dict(row) for row in artifacts]
    wrong_bem_row_count[2]["surface_row_count"] = 5
    wrong_bem_row_count_gate = pm_demag_margin_screening_package_gate(
        wrong_bem_row_count,
        expected_bem_surface_row_count=6,
    )
    assert wrong_bem_row_count_gate["status"] == "needs_attention"
    assert wrong_bem_row_count_gate["checks"]["expected_bem_surface_row_count_matches"] is False

    incomplete_package = [dict(row) for row in artifacts]
    incomplete_package[3]["required_artifacts"] = ["run_result", "loadline_metadata"]
    incomplete_package_gate = pm_demag_margin_screening_package_gate(incomplete_package)
    assert incomplete_package_gate["status"] == "needs_attention"
    assert incomplete_package_gate["checks"]["demag_package_artifacts_complete"] is False

    stale_probe = [dict(row) for row in artifacts]
    stale_probe[1]["field_probe_id"] = "fault_probe_pm_A_old"
    stale_probe_gate = pm_demag_margin_screening_package_gate(
        stale_probe,
        expected_field_probe_id="fault_probe_pm_A_demag_axis_v1",
        require_field_probe_identity=True,
    )
    assert stale_probe_gate["status"] == "needs_attention"
    assert stale_probe_gate["checks"]["expected_field_probe_id_matches"] is False

    wrong_probe_family = [dict(row) for row in artifacts]
    wrong_probe_family[1]["field_probe_family"] = "thermal_field_probe"
    wrong_probe_family_gate = pm_demag_margin_screening_package_gate(
        wrong_probe_family,
        expected_field_probe_id="fault_probe_pm_A_demag_axis_v1",
        expected_field_probe_family="demag_margin_field_probe",
        require_field_probe_identity=True,
    )
    assert wrong_probe_family_gate["status"] == "needs_attention"
    assert wrong_probe_family_gate["checks"]["expected_field_probe_id_matches"] is True
    assert wrong_probe_family_gate["checks"]["expected_field_probe_family_matches"] is False

    wrong_axis = [dict(row) for row in artifacts]
    wrong_axis[1]["field_axis_convention"] = "global_z_positive"
    wrong_axis_gate = pm_demag_margin_screening_package_gate(
        wrong_axis,
        expected_field_axis_convention="magnetization_axis_positive",
        require_field_probe_identity=True,
    )
    assert wrong_axis_gate["status"] == "needs_attention"
    assert wrong_axis_gate["checks"]["field_axis_convention_recorded"] is True
    assert wrong_axis_gate["checks"]["expected_field_axis_convention_matches"] is False

    wrong_sign = [dict(row) for row in artifacts]
    wrong_sign[1]["field_sign_convention"] = "positive_h_parallel_is_demag"
    wrong_sign_gate = pm_demag_margin_screening_package_gate(
        wrong_sign,
        expected_field_sign_convention="negative_h_parallel_is_demag",
        require_field_probe_identity=True,
    )
    assert wrong_sign_gate["status"] == "needs_attention"
    assert wrong_sign_gate["checks"]["field_sign_convention_recorded"] is True
    assert wrong_sign_gate["checks"]["expected_field_sign_convention_matches"] is False

    missing_average = [dict(row) for row in artifacts]
    missing_average[1].pop("averaging_rule")
    missing_average_gate = pm_demag_margin_screening_package_gate(
        missing_average,
        require_field_probe_identity=True,
    )
    assert missing_average_gate["status"] == "needs_attention"
    assert missing_average_gate["checks"]["averaging_rule_recorded"] is False


def test_pm_demag_margin_screening_package_gate_records_slot166_elf_manifest_identity():
    artifacts = [
        {
            "kind": "loadline_sweep",
            "case_id": "demag_case_166",
            "magnet_id": "pm_spm_A",
            "temperature_C": 120.0,
            "material_state": {"Br_T": 1.068, "H_knee_A_per_m": -495000.0, "recoil_mu_r": 1.05},
            "material_state_artifact_id": "elf_slot334_pm_spm_A_hbrm_hbcn_state_v1.json",
            "material_state_digest": "sha256:elf_slot334_pm_spm_A_hbrm_hbcn_state_v1",
            "load_step_id": "elf_slot342_loadline_hot_120c_v1",
            "path": "slot166_loadline_sweep.json",
            "gate_policy": "pm_temperature_demag_sweep_summary",
            "status": "ok",
            "minimum_demag_margin_A_per_m": -18000.0,
            "risk_label": "red",
        },
        {
            "kind": "fault_current_screening",
            "case_id": "demag_case_166",
            "magnet_id": "pm_spm_A",
            "temperature_C": 120.0,
            "material_state": {"Br_T": 1.068, "H_knee_A_per_m": -495000.0, "recoil_mu_r": 1.05},
            "material_state_artifact_id": "elf_slot334_pm_spm_A_hbrm_hbcn_state_v1.json",
            "material_state_digest": "sha256:elf_slot334_pm_spm_A_hbrm_hbcn_state_v1",
            "fault_step_id": "elf_slot342_negative_id_fault_step_v1",
            "path": "slot166_fault_current_screening.json",
            "gate_policy": "fault_current_demag_screening",
            "status": "ok",
            "negative_id_is_demag_direction": True,
            "recommended_observable_keys": [
                "field_probe",
                "demag_margin_A_per_m",
                "recoil_remanence_ratio_proxy",
            ],
            "field_probe_id": "elf_slot286_field_probe_pm_spm_A_demag_axis_v1",
            "field_probe_family": "elf_demag_margin_field_probe",
            "observation_region_id": "pm_spm_A_volume",
            "observation_component": "H_parallel_demag_axis",
            "field_axis_convention": "magnetization_axis_positive",
            "field_sign_convention": "negative_h_parallel_is_demag",
            "field_probe_method": "elf_volume_average_h_parallel_probe",
            "averaging_rule": "volume_average",
            "field_probe_geometry_digest": "sha256:elf_slot326_field_probe_volume_pm_spm_A",
            "field_probe_point_xyz_m": [0.028, 0.0, 0.0],
            "field_probe_output_artifact_id": "elf_slot294_field_probe_table_pm_spm_A_v1",
            "field_probe_output_digest": "sha256:elf_slot294_field_probe_table_pm_spm_A",
            "field_probe_output_path": "artifacts/field/elf_slot294_field_probe_table.json",
        },
        {
            "kind": "bem_source_balance",
            "case_id": "demag_case_166",
            "magnet_id": "pm_spm_A",
            "temperature_C": 120.0,
            "material_state": {"Br_T": 1.068, "H_knee_A_per_m": -495000.0, "recoil_mu_r": 1.05},
            "material_state_artifact_id": "elf_slot334_pm_spm_A_hbrm_hbcn_state_v1.json",
            "material_state_digest": "sha256:elf_slot334_pm_spm_A_hbrm_hbcn_state_v1",
            "path": "slot166_bem_source_balance.json",
            "gate_policy": "pm_bem_surface_source_balance_gate",
            "status": "ok",
            "surface_mesh_id": "bem_surface_mesh_pm_spm_A_v1",
            "source_balance_artifact_id": "bem_source_balance_pm_spm_A_v1",
            "source_balance_digest": "sha256:bem_source_balance_pm_spm_A_v1",
            "source_convention": "sigma_m_equals_m_dot_n",
            "signed_charge_balance_rel_tol": 1.0e-9,
            "signed_charge_balance_abs": 0.0,
            "total_area_m2": 0.031,
            "normal_convention": "outward_from_magnet",
            "source_balance_unit": "area_weighted_m_dot_n",
        },
        {
            "kind": "demag_package",
            "case_id": "demag_case_166",
            "magnet_id": "pm_spm_A",
            "temperature_C": 120.0,
            "material_state": {"Br_T": 1.068, "H_knee_A_per_m": -495000.0, "recoil_mu_r": 1.05},
            "material_state_artifact_id": "elf_slot334_pm_spm_A_hbrm_hbcn_state_v1.json",
            "material_state_digest": "sha256:elf_slot334_pm_spm_A_hbrm_hbcn_state_v1",
            "demag_step_id": "elf_slot342_hbcn_demag_step2_v1",
            "path": "slot166_demag_package.json",
            "gate_policy": "pm_demag_package_identity_gate",
            "status": "ok",
            "required_artifacts": [
                "run_result",
                "loadline_metadata",
                "bem_surface",
                "recoil_steps",
            ],
        },
    ]

    gate = pm_demag_margin_screening_package_gate(
        artifacts,
        expected_case_id="demag_case_166",
        expected_magnet_id="pm_spm_A",
        expected_temperature_c=120.0,
        expected_bem_source_balance_artifact_id="bem_source_balance_pm_spm_A_v1",
        expected_bem_source_balance_digest="sha256:bem_source_balance_pm_spm_A_v1",
        expected_bem_source_convention="sigma_m_equals_m_dot_n",
        expected_field_probe_id="elf_slot286_field_probe_pm_spm_A_demag_axis_v1",
        expected_field_probe_family="elf_demag_margin_field_probe",
        expected_observation_region_id="pm_spm_A_volume",
        expected_observation_component="H_parallel_demag_axis",
        expected_field_axis_convention="magnetization_axis_positive",
        expected_field_sign_convention="negative_h_parallel_is_demag",
        expected_field_probe_method="elf_volume_average_h_parallel_probe",
        expected_averaging_rule="volume_average",
        expected_field_probe_geometry_digest="sha256:elf_slot326_field_probe_volume_pm_spm_A",
        expected_field_probe_point_xyz_m=(0.028, 0.0, 0.0),
        expected_field_probe_output_artifact_id="elf_slot294_field_probe_table_pm_spm_A_v1",
        expected_field_probe_output_digest="sha256:elf_slot294_field_probe_table_pm_spm_A",
        expected_material_state_artifact_id="elf_slot334_pm_spm_A_hbrm_hbcn_state_v1.json",
        expected_material_state_digest="sha256:elf_slot334_pm_spm_A_hbrm_hbcn_state_v1",
        expected_load_step_id="elf_slot342_loadline_hot_120c_v1",
        expected_fault_step_id="elf_slot342_negative_id_fault_step_v1",
        expected_demag_step_id="elf_slot342_hbcn_demag_step2_v1",
        require_field_probe_identity=True,
        require_field_probe_output_artifact=True,
    )

    assert gate["status"] == "ok"
    assert gate["case_ids"] == ["demag_case_166"]
    assert gate["magnet_ids"] == ["pm_spm_A"]
    assert gate["material_states"] == [
        {"Br_T": 1.068, "H_knee_A_per_m": -495000.0, "recoil_mu_r": 1.05}
    ]
    assert gate["material_state_artifact_ids"] == [
        "elf_slot334_pm_spm_A_hbrm_hbcn_state_v1.json"
    ]
    assert gate["material_state_digests"] == [
        "sha256:elf_slot334_pm_spm_A_hbrm_hbcn_state_v1"
    ]
    assert gate["load_step_ids"] == ["elf_slot342_loadline_hot_120c_v1"]
    assert gate["fault_step_ids"] == ["elf_slot342_negative_id_fault_step_v1"]
    assert gate["demag_step_ids"] == ["elf_slot342_hbcn_demag_step2_v1"]
    assert gate["checks"]["expected_magnet_id_matches"] is True
    assert gate["checks"]["bem_surface_area_recorded"] is True
    assert gate["checks"]["bem_normal_convention_recorded"] is True
    assert gate["checks"]["bem_source_unit_recorded"] is True
    assert gate["checks"]["bem_balance_value_recorded"] is True
    assert gate["checks"]["bem_surface_mesh_id_recorded"] is True
    assert gate["checks"]["bem_source_balance_artifact_id_recorded"] is True
    assert gate["checks"]["bem_source_balance_artifact_id_unique_when_present"] is True
    assert gate["checks"]["bem_source_balance_digest_recorded"] is True
    assert gate["checks"]["bem_source_balance_digest_unique_when_present"] is True
    assert gate["checks"]["bem_source_convention_recorded"] is True
    assert gate["checks"]["bem_source_convention_unique_when_present"] is True
    assert gate["bem_source_balance_artifact_ids"] == ["bem_source_balance_pm_spm_A_v1"]
    assert gate["bem_source_balance_digests"] == ["sha256:bem_source_balance_pm_spm_A_v1"]
    assert gate["bem_source_conventions"] == ["sigma_m_equals_m_dot_n"]
    assert gate["checks"]["material_state_complete_when_present"] is True
    assert gate["checks"]["material_state_unique_when_present"] is True
    assert gate["checks"]["material_state_artifact_id_recorded"] is True
    assert gate["checks"]["material_state_artifact_id_unique_when_present"] is True
    assert gate["checks"]["material_state_digest_recorded"] is True
    assert gate["checks"]["material_state_digest_unique_when_present"] is True
    assert gate["field_probe_ids"] == ["elf_slot286_field_probe_pm_spm_A_demag_axis_v1"]
    assert gate["field_probe_families"] == ["elf_demag_margin_field_probe"]
    assert gate["observation_region_ids"] == ["pm_spm_A_volume"]
    assert gate["observation_components"] == ["h_parallel_demag_axis"]
    assert gate["field_axis_conventions"] == ["magnetization_axis_positive"]
    assert gate["field_sign_conventions"] == ["negative_h_parallel_is_demag"]
    assert gate["field_probe_methods"] == ["elf_volume_average_h_parallel_probe"]
    assert gate["averaging_rules"] == ["volume_average"]
    assert gate["field_probe_geometry_digests"] == [
        "sha256:elf_slot326_field_probe_volume_pm_spm_A"
    ]
    assert gate["field_probe_points_xyz_m"] == [[0.028, 0.0, 0.0]]
    assert gate["field_probe_output_artifact_ids"] == [
        "elf_slot294_field_probe_table_pm_spm_A_v1"
    ]
    assert gate["field_probe_output_digests"] == [
        "sha256:elf_slot294_field_probe_table_pm_spm_A"
    ]
    assert gate["field_probe_output_paths"] == [
        "artifacts/field/elf_slot294_field_probe_table.json"
    ]
    assert gate["checks"]["field_probe_id_recorded"] is True
    assert gate["checks"]["field_probe_family_recorded"] is True
    assert gate["checks"]["observation_region_id_recorded"] is True
    assert gate["checks"]["observation_component_recorded"] is True
    assert gate["checks"]["field_axis_convention_recorded"] is True
    assert gate["checks"]["field_sign_convention_recorded"] is True
    assert gate["checks"]["field_probe_method_recorded_when_expected"] is True
    assert gate["checks"]["averaging_rule_recorded"] is True
    assert gate["checks"]["field_probe_geometry_digest_recorded"] is True
    assert gate["checks"]["field_probe_point_xyz_recorded_when_expected"] is True
    assert gate["checks"]["field_probe_output_artifact_id_recorded"] is True
    assert gate["checks"]["field_probe_output_digest_recorded"] is True
    assert gate["checks"]["field_probe_output_path_recorded"] is True
    assert gate["checks"]["expected_field_probe_id_matches"] is True
    assert gate["checks"]["expected_field_probe_family_matches"] is True
    assert gate["checks"]["expected_observation_region_id_matches"] is True
    assert gate["checks"]["expected_observation_component_matches"] is True
    assert gate["checks"]["expected_field_axis_convention_matches"] is True
    assert gate["checks"]["expected_field_sign_convention_matches"] is True
    assert gate["checks"]["expected_field_probe_method_matches"] is True
    assert gate["checks"]["expected_averaging_rule_matches"] is True
    assert gate["checks"]["expected_field_probe_geometry_digest_matches"] is True
    assert gate["checks"]["expected_field_probe_point_xyz_matches"] is True
    assert gate["checks"]["expected_field_probe_output_artifact_id_matches"] is True
    assert gate["checks"]["expected_field_probe_output_digest_matches"] is True
    assert gate["checks"]["expected_bem_source_balance_artifact_id_matches"] is True
    assert gate["checks"]["expected_bem_source_balance_digest_matches"] is True
    assert gate["checks"]["expected_bem_source_convention_matches"] is True
    assert gate["checks"]["expected_material_state_artifact_id_matches"] is True
    assert gate["checks"]["expected_material_state_digest_matches"] is True
    assert gate["checks"]["load_step_id_recorded"] is True
    assert gate["checks"]["fault_step_id_recorded"] is True
    assert gate["checks"]["demag_step_id_recorded"] is True
    assert gate["checks"]["expected_load_step_id_matches"] is True
    assert gate["checks"]["expected_fault_step_id_matches"] is True
    assert gate["checks"]["expected_demag_step_id_matches"] is True
    assert "demag-margin panels cannot mix" in gate["version_note"]
    assert "load/fault/demag step identities" in gate["version_note"]
    assert "material-state artifact id and digest" in gate["version_note"]

    stale_magnet = [dict(row) for row in artifacts]
    stale_magnet[3]["magnet_id"] = "pm_old"
    stale_gate = pm_demag_margin_screening_package_gate(stale_magnet)
    assert stale_gate["status"] == "needs_attention"
    assert stale_gate["checks"]["magnet_ids_unique"] is False

    positive_id = [dict(row) for row in artifacts]
    positive_id[1]["negative_id_is_demag_direction"] = False
    positive_gate = pm_demag_margin_screening_package_gate(positive_id)
    assert positive_gate["status"] == "needs_attention"
    assert positive_gate["checks"]["fault_current_demag_direction_recorded"] is False

    missing_bem_tol = [dict(row) for row in artifacts]
    missing_bem_tol[2].pop("signed_charge_balance_rel_tol")
    bem_gate = pm_demag_margin_screening_package_gate(missing_bem_tol)
    assert bem_gate["status"] == "needs_attention"
    assert bem_gate["checks"]["bem_balance_tolerance_recorded"] is False

    missing_bem_normal = [dict(row) for row in artifacts]
    missing_bem_normal[2].pop("normal_convention")
    missing_bem_normal_gate = pm_demag_margin_screening_package_gate(missing_bem_normal)
    assert missing_bem_normal_gate["status"] == "needs_attention"
    assert missing_bem_normal_gate["checks"]["bem_normal_convention_recorded"] is False

    incomplete_package = [dict(row) for row in artifacts]
    incomplete_package[3]["required_artifacts"] = ["run_result", "loadline_metadata", "bem_surface"]
    incomplete_gate = pm_demag_margin_screening_package_gate(incomplete_package)
    assert incomplete_gate["status"] == "needs_attention"
    assert incomplete_gate["checks"]["demag_package_artifacts_complete"] is False

    stale_material = [dict(row) for row in artifacts]
    stale_material[0] = {
        **stale_material[0],
        "material_state": {"Br_T": 1.02, "H_knee_A_per_m": -495000.0, "recoil_mu_r": 1.05},
    }
    stale_material_gate = pm_demag_margin_screening_package_gate(stale_material)
    assert stale_material_gate["status"] == "needs_attention"
    assert stale_material_gate["checks"]["material_state_unique_when_present"] is False

    stale_material_artifact = [dict(row) for row in artifacts]
    stale_material_artifact[0]["material_state_artifact_id"] = (
        "elf_slot334_pm_spm_A_hbrm_hbcn_state_old.json"
    )
    stale_material_artifact_gate = pm_demag_margin_screening_package_gate(
        stale_material_artifact,
        expected_material_state_artifact_id="elf_slot334_pm_spm_A_hbrm_hbcn_state_v1.json",
        expected_material_state_digest="sha256:elf_slot334_pm_spm_A_hbrm_hbcn_state_v1",
    )
    assert stale_material_artifact_gate["status"] == "needs_attention"
    assert (
        stale_material_artifact_gate["checks"]["expected_material_state_artifact_id_matches"]
        is False
    )
    assert stale_material_artifact_gate["checks"]["expected_material_state_digest_matches"] is True

    stale_material_digest = [dict(row) for row in artifacts]
    stale_material_digest[0]["material_state_digest"] = (
        "sha256:elf_slot334_pm_spm_A_hbrm_hbcn_state_old"
    )
    stale_material_digest_gate = pm_demag_margin_screening_package_gate(
        stale_material_digest,
        expected_material_state_artifact_id="elf_slot334_pm_spm_A_hbrm_hbcn_state_v1.json",
        expected_material_state_digest="sha256:elf_slot334_pm_spm_A_hbrm_hbcn_state_v1",
    )
    assert stale_material_digest_gate["status"] == "needs_attention"
    assert stale_material_digest_gate["checks"]["expected_material_state_artifact_id_matches"] is True
    assert stale_material_digest_gate["checks"]["expected_material_state_digest_matches"] is False

    stale_bem_source_artifact = [dict(row) for row in artifacts]
    stale_bem_source_artifact[2]["source_balance_artifact_id"] = "bem_source_balance_pm_spm_A_old"
    stale_bem_source_artifact_gate = pm_demag_margin_screening_package_gate(
        stale_bem_source_artifact,
        expected_bem_source_balance_artifact_id="bem_source_balance_pm_spm_A_v1",
        expected_bem_source_balance_digest="sha256:bem_source_balance_pm_spm_A_v1",
        expected_bem_source_convention="sigma_m_equals_m_dot_n",
    )
    assert stale_bem_source_artifact_gate["status"] == "needs_attention"
    assert (
        stale_bem_source_artifact_gate["checks"]["expected_bem_source_balance_artifact_id_matches"]
        is False
    )
    assert (
        stale_bem_source_artifact_gate["checks"]["expected_bem_source_balance_digest_matches"]
        is True
    )
    assert stale_bem_source_artifact_gate["checks"]["expected_bem_source_convention_matches"] is True

    stale_bem_source_digest = [dict(row) for row in artifacts]
    stale_bem_source_digest[2]["source_balance_digest"] = "sha256:bem_source_balance_pm_spm_A_old"
    stale_bem_source_digest_gate = pm_demag_margin_screening_package_gate(
        stale_bem_source_digest,
        expected_bem_source_balance_artifact_id="bem_source_balance_pm_spm_A_v1",
        expected_bem_source_balance_digest="sha256:bem_source_balance_pm_spm_A_v1",
        expected_bem_source_convention="sigma_m_equals_m_dot_n",
    )
    assert stale_bem_source_digest_gate["status"] == "needs_attention"
    assert (
        stale_bem_source_digest_gate["checks"]["expected_bem_source_balance_artifact_id_matches"]
        is True
    )
    assert (
        stale_bem_source_digest_gate["checks"]["expected_bem_source_balance_digest_matches"]
        is False
    )
    assert stale_bem_source_digest_gate["checks"]["expected_bem_source_convention_matches"] is True

    wrong_bem_source_convention = [dict(row) for row in artifacts]
    wrong_bem_source_convention[2]["source_convention"] = "sigma_m_equals_negative_m_dot_n"
    wrong_bem_source_convention_gate = pm_demag_margin_screening_package_gate(
        wrong_bem_source_convention,
        expected_bem_source_balance_artifact_id="bem_source_balance_pm_spm_A_v1",
        expected_bem_source_balance_digest="sha256:bem_source_balance_pm_spm_A_v1",
        expected_bem_source_convention="sigma_m_equals_m_dot_n",
    )
    assert wrong_bem_source_convention_gate["status"] == "needs_attention"
    assert (
        wrong_bem_source_convention_gate["checks"]["expected_bem_source_balance_artifact_id_matches"]
        is True
    )
    assert (
        wrong_bem_source_convention_gate["checks"]["expected_bem_source_balance_digest_matches"]
        is True
    )
    assert wrong_bem_source_convention_gate["checks"]["expected_bem_source_convention_matches"] is False

    missing_material = [dict(row) for row in artifacts]
    missing_material[2].pop("material_state")
    missing_material_gate = pm_demag_margin_screening_package_gate(missing_material)
    assert missing_material_gate["status"] == "needs_attention"
    assert missing_material_gate["checks"]["material_state_complete_when_present"] is False

    missing_material_identity = [dict(row) for row in artifacts]
    missing_material_identity[2].pop("material_state_artifact_id")
    missing_material_identity[2].pop("material_state_digest")
    missing_material_identity_gate = pm_demag_margin_screening_package_gate(
        missing_material_identity,
        expected_material_state_artifact_id="elf_slot334_pm_spm_A_hbrm_hbcn_state_v1.json",
        expected_material_state_digest="sha256:elf_slot334_pm_spm_A_hbrm_hbcn_state_v1",
    )
    assert missing_material_identity_gate["status"] == "needs_attention"
    assert (
        missing_material_identity_gate["checks"]["material_state_artifact_id_recorded"]
        is False
    )
    assert missing_material_identity_gate["checks"]["material_state_digest_recorded"] is False

    stale_fault_step = [dict(row) for row in artifacts]
    stale_fault_step[1]["fault_step_id"] = "elf_slot342_fault_step_from_old_manifest"
    stale_fault_step_gate = pm_demag_margin_screening_package_gate(
        stale_fault_step,
        expected_fault_step_id="elf_slot342_negative_id_fault_step_v1",
    )
    assert stale_fault_step_gate["status"] == "needs_attention"
    assert stale_fault_step_gate["checks"]["fault_step_id_recorded"] is True
    assert stale_fault_step_gate["checks"]["expected_fault_step_id_matches"] is False

    missing_demag_step = [dict(row) for row in artifacts]
    missing_demag_step[3].pop("demag_step_id")
    missing_demag_step_gate = pm_demag_margin_screening_package_gate(
        missing_demag_step,
        expected_demag_step_id="elf_slot342_hbcn_demag_step2_v1",
    )
    assert missing_demag_step_gate["status"] == "needs_attention"
    assert missing_demag_step_gate["checks"]["demag_step_id_recorded"] is False
    assert missing_demag_step_gate["checks"]["expected_demag_step_id_matches"] is False

    wrong_component = [dict(row) for row in artifacts]
    wrong_component[1]["observation_component"] = "B_magnitude"
    wrong_component_gate = pm_demag_margin_screening_package_gate(
        wrong_component,
        expected_observation_component="H_parallel_demag_axis",
        require_field_probe_identity=True,
    )
    assert wrong_component_gate["status"] == "needs_attention"
    assert wrong_component_gate["checks"]["expected_observation_component_matches"] is False

    missing_region = [dict(row) for row in artifacts]
    missing_region[1].pop("observation_region_id")
    missing_region_gate = pm_demag_margin_screening_package_gate(
        missing_region,
        require_field_probe_identity=True,
    )
    assert missing_region_gate["status"] == "needs_attention"
    assert missing_region_gate["checks"]["observation_region_id_recorded"] is False

    stale_probe_geometry = [dict(row) for row in artifacts]
    stale_probe_geometry[1]["field_probe_geometry_digest"] = (
        "sha256:elf_slot326_field_probe_volume_pm_spm_A_old"
    )
    stale_probe_geometry_gate = pm_demag_margin_screening_package_gate(
        stale_probe_geometry,
        expected_field_probe_geometry_digest="sha256:elf_slot326_field_probe_volume_pm_spm_A",
        expected_field_probe_point_xyz_m=(0.028, 0.0, 0.0),
        require_field_probe_identity=True,
    )
    assert stale_probe_geometry_gate["status"] == "needs_attention"
    assert stale_probe_geometry_gate["checks"]["expected_field_probe_geometry_digest_matches"] is False
    assert stale_probe_geometry_gate["checks"]["expected_field_probe_point_xyz_matches"] is True

    stale_probe_point = [dict(row) for row in artifacts]
    stale_probe_point[1]["field_probe_point_xyz_m"] = [0.029, 0.0, 0.0]
    stale_probe_point_gate = pm_demag_margin_screening_package_gate(
        stale_probe_point,
        expected_field_probe_geometry_digest="sha256:elf_slot326_field_probe_volume_pm_spm_A",
        expected_field_probe_point_xyz_m=(0.028, 0.0, 0.0),
        require_field_probe_identity=True,
    )
    assert stale_probe_point_gate["status"] == "needs_attention"
    assert stale_probe_point_gate["checks"]["expected_field_probe_geometry_digest_matches"] is True
    assert stale_probe_point_gate["checks"]["expected_field_probe_point_xyz_matches"] is False

    missing_probe_geometry = [dict(row) for row in artifacts]
    missing_probe_geometry[1].pop("field_probe_geometry_digest")
    missing_probe_geometry[1].pop("field_probe_point_xyz_m")
    missing_probe_geometry_gate = pm_demag_margin_screening_package_gate(
        missing_probe_geometry,
        expected_field_probe_geometry_digest="sha256:elf_slot326_field_probe_volume_pm_spm_A",
        expected_field_probe_point_xyz_m=(0.028, 0.0, 0.0),
        require_field_probe_identity=True,
    )
    assert missing_probe_geometry_gate["status"] == "needs_attention"
    assert missing_probe_geometry_gate["checks"]["field_probe_geometry_digest_recorded"] is False
    assert missing_probe_geometry_gate["checks"]["field_probe_point_xyz_recorded_when_expected"] is False

    stale_output_artifact = [dict(row) for row in artifacts]
    stale_output_artifact[1]["field_probe_output_artifact_id"] = "elf_slot286_old_probe_table"
    stale_output_artifact_gate = pm_demag_margin_screening_package_gate(
        stale_output_artifact,
        expected_field_probe_output_artifact_id="elf_slot294_field_probe_table_pm_spm_A_v1",
        expected_field_probe_output_digest="sha256:elf_slot294_field_probe_table_pm_spm_A",
        require_field_probe_output_artifact=True,
    )
    assert stale_output_artifact_gate["status"] == "needs_attention"
    assert (
        stale_output_artifact_gate["checks"]["expected_field_probe_output_artifact_id_matches"]
        is False
    )
    assert (
        stale_output_artifact_gate["checks"]["expected_field_probe_output_digest_matches"]
        is True
    )

    stale_output_digest = [dict(row) for row in artifacts]
    stale_output_digest[1]["field_probe_output_digest"] = "sha256:old_probe_table"
    stale_output_digest_gate = pm_demag_margin_screening_package_gate(
        stale_output_digest,
        expected_field_probe_output_artifact_id="elf_slot294_field_probe_table_pm_spm_A_v1",
        expected_field_probe_output_digest="sha256:elf_slot294_field_probe_table_pm_spm_A",
        require_field_probe_output_artifact=True,
    )
    assert stale_output_digest_gate["status"] == "needs_attention"
    assert (
        stale_output_digest_gate["checks"]["expected_field_probe_output_artifact_id_matches"]
        is True
    )
    assert (
        stale_output_digest_gate["checks"]["expected_field_probe_output_digest_matches"]
        is False
    )

    wrong_probe_family = [dict(row) for row in artifacts]
    wrong_probe_family[1]["field_probe_family"] = "elf_force_probe"
    wrong_probe_family_gate = pm_demag_margin_screening_package_gate(
        wrong_probe_family,
        expected_field_probe_id="elf_slot286_field_probe_pm_spm_A_demag_axis_v1",
        expected_field_probe_family="elf_demag_margin_field_probe",
        expected_field_probe_output_artifact_id="elf_slot294_field_probe_table_pm_spm_A_v1",
        expected_field_probe_output_digest="sha256:elf_slot294_field_probe_table_pm_spm_A",
        require_field_probe_identity=True,
        require_field_probe_output_artifact=True,
    )
    assert wrong_probe_family_gate["status"] == "needs_attention"
    assert wrong_probe_family_gate["checks"]["expected_field_probe_id_matches"] is True
    assert wrong_probe_family_gate["checks"]["expected_field_probe_family_matches"] is False
    assert (
        wrong_probe_family_gate["checks"]["expected_field_probe_output_artifact_id_matches"]
        is True
    )
    assert wrong_probe_family_gate["checks"]["expected_field_probe_output_digest_matches"] is True

    wrong_probe_method = [dict(row) for row in artifacts]
    wrong_probe_method[1]["field_probe_method"] = "elf_point_sample_h_parallel_probe"
    wrong_probe_method_gate = pm_demag_margin_screening_package_gate(
        wrong_probe_method,
        expected_field_probe_id="elf_slot286_field_probe_pm_spm_A_demag_axis_v1",
        expected_field_probe_family="elf_demag_margin_field_probe",
        expected_field_probe_method="elf_volume_average_h_parallel_probe",
        require_field_probe_identity=True,
    )
    assert wrong_probe_method_gate["status"] == "needs_attention"
    assert wrong_probe_method_gate["checks"]["expected_field_probe_id_matches"] is True
    assert wrong_probe_method_gate["checks"]["expected_field_probe_family_matches"] is True
    assert wrong_probe_method_gate["checks"]["field_probe_method_recorded_when_expected"] is True
    assert wrong_probe_method_gate["checks"]["expected_field_probe_method_matches"] is False

    missing_probe_method = [dict(row) for row in artifacts]
    missing_probe_method[1].pop("field_probe_method")
    missing_probe_method_gate = pm_demag_margin_screening_package_gate(
        missing_probe_method,
        expected_field_probe_method="elf_volume_average_h_parallel_probe",
        require_field_probe_identity=True,
    )
    assert missing_probe_method_gate["status"] == "needs_attention"
    assert missing_probe_method_gate["checks"]["field_probe_method_recorded_when_expected"] is False
    assert missing_probe_method_gate["checks"]["expected_field_probe_method_matches"] is False

    missing_output_path = [dict(row) for row in artifacts]
    missing_output_path[1].pop("field_probe_output_path")
    missing_output_path_gate = pm_demag_margin_screening_package_gate(
        missing_output_path,
        require_field_probe_output_artifact=True,
    )
    assert missing_output_path_gate["status"] == "needs_attention"
    assert missing_output_path_gate["checks"]["field_probe_output_path_recorded"] is False


def test_box_projected_gradient_gate_closes_matlab_optimization_contract():
    gate = box_projected_gradient_least_squares_gate(
        matrix=[[1.0, 0.0], [0.0, 1.0]],
        rhs=[2.0, -1.0],
        lower=[0.0, 0.0],
        upper=[1.0, 3.0],
        initial=[0.0, 0.0],
        step_size=1.0,
        max_iterations=5,
    )

    assert gate["status"] == "ok"
    assert gate["x"] == pytest.approx([1.0, 0.0])
    assert gate["gradient"] == pytest.approx([-1.0, 1.0])
    assert gate["active_upper"] == [True, False]
    assert gate["active_lower"] == [False, True]
    assert gate["objective_history"][0] == pytest.approx(2.5)
    assert gate["objective"] == pytest.approx(1.0)
    assert gate["max_kkt_residual"] == pytest.approx(0.0)
    assert gate["projected_gradient_residual"] == pytest.approx(0.0)

    bad_bounds = dict(
        matrix=[[1.0]],
        rhs=[1.0],
        lower=[2.0],
        upper=[1.0],
    )
    with pytest.raises(ValueError):
        box_projected_gradient_least_squares_gate(**bad_bounds)


def test_geometric_integrator_energy_drift_gate_closes_matlab_teaching_contract():
    rows = [
        {
            "method": "explicit_euler",
            "energy_initial": 0.5,
            "energy_final": 0.746,
            "max_rel_energy_drift": 0.491,
            "steps": 1000,
            "step_size_s": 0.02,
            "omega_rad_per_s": 1.0,
        },
        {
            "method": "symplectic_euler",
            "energy_initial": 0.5,
            "energy_final": 0.502,
            "max_rel_energy_drift": 0.0102,
            "steps": 1000,
            "step_size_s": 0.02,
            "omega_rad_per_s": 1.0,
        },
        {
            "method": "implicit_midpoint",
            "energy_initial": 0.5,
            "energy_final": 0.5,
            "max_rel_energy_drift": 2.0e-14,
            "steps": 1000,
            "step_size_s": 0.02,
            "omega_rad_per_s": 1.0,
        },
    ]

    gate = geometric_integrator_energy_drift_gate(
        rows,
        max_geometric_rel_drift=2.0e-2,
        min_explicit_to_geometric_drift_ratio=20.0,
    )

    assert gate["status"] == "ok"
    assert gate["policy"] == "geometric_integrator_energy_drift_gate"
    assert gate["common_steps"] == 1000
    assert gate["common_step_size_s"] == pytest.approx(0.02)
    assert gate["common_omega_rad_per_s"] == pytest.approx(1.0)
    assert gate["best_geometric_method"] == "implicit_midpoint"
    assert gate["max_geometric_rel_energy_drift"] == pytest.approx(0.0102)
    assert gate["explicit_to_worst_geometric_drift_ratio"] > 20.0
    assert gate["checks"]["geometric_energy_drift_bounded"] is True

    bad_symplectic = [dict(row) for row in rows]
    bad_symplectic[1]["max_rel_energy_drift"] = 0.20
    bad_gate = geometric_integrator_energy_drift_gate(
        bad_symplectic,
        max_geometric_rel_drift=2.0e-2,
        min_explicit_to_geometric_drift_ratio=20.0,
    )
    assert bad_gate["status"] == "needs_attention"
    assert bad_gate["checks"]["geometric_energy_drift_bounded"] is False

    weak_negative_control = [dict(row) for row in rows]
    weak_negative_control[0]["max_rel_energy_drift"] = 0.011
    weak_gate = geometric_integrator_energy_drift_gate(
        weak_negative_control,
        max_geometric_rel_drift=2.0e-2,
        min_explicit_to_geometric_drift_ratio=20.0,
    )
    assert weak_gate["status"] == "needs_attention"
    assert weak_gate["checks"]["explicit_drift_larger_than_geometric"] is False


def test_morozov_discrepancy_choice_closes_matlab_regularization_contract():
    gate = morozov_discrepancy_choice(
        alphas=[0.0, 1.0e-3, 1.0e-2, 1.0e-1, 1.0],
        residual_norms=[0.02, 0.08, 0.24, 0.71, 1.80],
        noise_norm=0.70,
    )

    assert gate["status"] == "ok"
    assert gate["policy"] == "morozov_discrepancy_regularization_choice"
    assert gate["selected_index"] == 4
    assert gate["selected_alpha"] == pytest.approx(1.0e-1)
    assert gate["selected_residual_norm"] == pytest.approx(0.71)
    assert gate["lower_bracket_index"] == 3
    assert gate["upper_bracket_index"] == 4
    assert gate["checks"]["noise_bracketed"] is True
    assert gate["checks"]["residuals_nondecreasing"] is True

    unbracketed = morozov_discrepancy_choice([0.0, 1.0], [0.02, 0.08], noise_norm=0.70)
    assert unbracketed["status"] == "needs_attention"
    assert unbracketed["checks"]["noise_bracketed"] is False

    nonmonotone = morozov_discrepancy_choice(
        [0.0, 1.0e-3, 1.0e-2],
        [0.02, 0.20, 0.10],
        noise_norm=0.11,
    )
    assert nonmonotone["status"] == "needs_attention"
    assert nonmonotone["checks"]["residuals_nondecreasing"] is False
    assert nonmonotone["selected_index0"] == 2


def test_lcurve_corner_choice_closes_matlab_regularization_path_contract():
    gate = lcurve_corner_choice(
        alphas=[1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1, 1.0],
        residual_norms=[
            0.010120936526662955,
            0.10019082610814475,
            0.9810590139741418,
            8.018627954622158,
            28.166442261676118,
        ],
        solution_norms=[
            54.76057558223518,
            54.75908791041406,
            54.61101899497311,
            45.92556510691625,
            14.405110742701566,
        ],
    )

    assert gate["status"] == "ok"
    assert gate["policy"] == "lcurve_corner_regularization_choice"
    assert gate["selected_index"] == 4
    assert gate["selected_alpha"] == pytest.approx(1.0e-1)
    assert gate["selected_residual_norm"] == pytest.approx(8.018627954622158)
    assert gate["selected_solution_norm"] == pytest.approx(45.92556510691625)
    assert gate["curvature"][0] == pytest.approx(0.0)
    assert gate["curvature"][-1] == pytest.approx(0.0)
    assert gate["checks"]["interior_selected"] is True
    assert gate["checks"]["residuals_nondecreasing"] is True
    assert gate["checks"]["solution_norms_nonincreasing"] is True

    overfit_endpoint = lcurve_corner_choice(
        alphas=[1.0e-4, 1.0e-3, 1.0e-2],
        residual_norms=[0.30, 0.20, 0.10],
        solution_norms=[1.0, 2.0, 3.0],
    )
    assert overfit_endpoint["status"] == "needs_attention"
    assert overfit_endpoint["checks"]["residuals_nondecreasing"] is False
    assert overfit_endpoint["checks"]["solution_norms_nonincreasing"] is False


def test_trace_surface_mass_energy_gate_closes_matlab_fem_bem_trace_contract():
    trace = [
        [1.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 1.0, 0.0],
    ]
    surface_mass = [[0.0 for _ in range(4)] for _ in range(4)]
    for tri, area in (
        ((0, 1, 2), 0.5),
        ((0, 3, 1), 0.5),
        ((0, 2, 3), 0.5),
        ((1, 3, 2), math.sqrt(3.0) / 2.0),
    ):
        local = area / 12.0
        for a, ia in enumerate(tri):
            for b, ib in enumerate(tri):
                surface_mass[ia][ib] += local * (2.0 if a == b else 1.0)

    gate = trace_surface_mass_energy_gate(
        trace,
        surface_mass,
        fem_values=[1.0, 2.0, -1.0, 0.5, 9.0],
    )

    assert gate["status"] == "ok"
    assert gate["trace_shape"] == [4, 5]
    assert gate["interior_fem_node_ids"] == [5]
    assert gate["surface_area_from_mass"] == pytest.approx(1.5 + math.sqrt(3.0) / 2.0)
    assert gate["energy_abs_error"] == pytest.approx(0.0)
    assert gate["max_interior_boundary_action"] == pytest.approx(0.0)

    bad_trace = [row[:] for row in trace]
    bad_trace[0][4] = 0.1
    bad = trace_surface_mass_energy_gate(bad_trace, surface_mass, [1.0, 2.0, -1.0, 0.5, 9.0])
    assert bad["status"] == "needs_attention"
    assert bad["checks"]["constant_trace_ok"] is False


def test_netgen_vol_first_order_fem_bem_trace_package_handoff_keeps_node_identity():
    package = {
        "mesh_id": "unit_tet_mesh",
        "export_id": "coreform_netgen_unit_tet",
        "trace_artifact_id": "unit_tet_h1_to_scalar_bem_trace_p1",
        "trace_operator_artifact_id": "unit_tet_h1_to_scalar_bem_trace_operator_p1",
        "trace_operator_policy": "one_hot_boundary_node_injection_from_vol_node_ids",
        "trace_output_artifact_id": "unit_tet_h1_to_scalar_bem_trace_output_p1",
        "trace_output_digest": "sha256:unit_tet_h1_to_scalar_bem_trace_output_p1",
        "trace_output_path": "unit_tet_h1_to_scalar_bem_trace_output.json",
        "trace_observable_id": "unit_tet_h1_to_scalar_bem_boundary_trace_v1",
        "trace_observable_family": "fem_bem_boundary_trace",
        "coupled_system_artifact_id": "unit_tet_laplace_fem_bem_schur_system_v1",
        "coupled_system_digest": "sha256:unit_tet_laplace_fem_bem_schur_system_v1",
        "result_artifact_id": "matlab_slot344_fem_bem_manifest_result_v1",
        "linear_solver_report_artifact_id": "matlab_slot367_fem_bem_linear_solver_report_v1",
        "linear_solver_report_digest": "sha256:matlab_slot367_fem_bem_linear_solver_report_v1",
        "linear_solver_name": "minimum_norm_pinv_rank_deficient",
        "linear_solver_tolerance": 1.0e-10,
        "linear_solver_residual_norm": 2.0e-13,
        "linear_solver_iteration_count": 1,
        "parameter_set_artifact_id": "matlab_slot388_fem_bem_parameter_set_v1",
        "parameter_set_digest": "sha256:matlab-slot388-fem-bem-parameter-set-v1",
        "parameter_set_path": "docs/fem_bem/first_order_fem_bem_parameter_set.json",
        "objective_observable_id": "matlab_slot388_trace_lsq_residual_objective_v1",
        "objective_observable_family": "fem_bem_trace_least_squares_objective",
        "run_started_at": "2026-07-01T13:50:00+09:00",
        "matlab_version": "R2026a",
        "timing_breakdown": {
            "mesh_read_s": 0.001,
            "trace_assembly_s": 0.002,
            "manifest_build_s": 0.003,
            "json_write_s": 0.004,
        },
        "surface_mesh_id": "unit_tet_boundary_tri_p1",
        "source_file_id": "sha256:unit_tet_vol_source_abc123",
        "source_path": "unit_tet.vol",
        "source_format": ".vol",
        "policy": "netgen_vol_tri_tet_only_shared_one_based_nodes",
        "polynomial_order": 1,
        "curved_element_count": 0,
        "coupling_kind": "h1_p1_to_scalar_bem_p1_trace",
        "formulation_id": "laplace_single_layer_teaching",
        "bem_kernel_family": "laplace_single_layer",
        "coupling_convention_schema_id": "matlab_first_order_fem_bem_coupling_convention_v1",
        "fem_bem_postprocess_row_convention_schema_id": "matlab_fem_bem_trace_lsq_row_convention_v1",
        "trace_basis_schema_id": "matlab_h1_p1_to_scalar_bem_p1_trace_basis_v1",
        "assembly_rule_id": "first_order_tet_h1_trace_tri_p1_bem_teaching_v1",
        "quadrature_rule_id": "tri_p1_exact_mass_regular_kernel_teaching_v1",
        "volume_space": "H1_P1_tet",
        "surface_space": "scalar_P1_tri",
        "geo": {
            "N": 5,
            "conn_matrix": [[1, 2, 3, 4], [1, 2, 3, 5]],
        },
        "gypsilab": {
            "elt": [[1, 2, 3], [1, 4, 2], [1, 3, 4], [2, 4, 3]],
            "col": [1, 1, 1, 1],
            "boundary_names": ["outer", "outer", "outer", "outer"],
            "boundary_row_identity": [
                {"surface_triangle_index": 1, "surface_triangle_nodes": [1, 2, 3], "boundary_number": 1, "boundary_name": "outer"},
                {"surface_triangle_index": 2, "surface_triangle_nodes": [1, 4, 2], "boundary_number": 1, "boundary_name": "outer"},
                {"surface_triangle_index": 3, "surface_triangle_nodes": [1, 3, 4], "boundary_number": 1, "boundary_name": "outer"},
                {"surface_triangle_index": 4, "surface_triangle_nodes": [2, 4, 3], "boundary_number": 1, "boundary_name": "outer"},
            ],
        },
        "trace": {
            "fem_node_ids": [1, 2, 3, 4],
            "bem_node_ids": [1, 2, 3, 4],
            "surface_triangles": [[1, 2, 3], [1, 4, 2], [1, 3, 4], [2, 4, 3]],
            "boundary_numbers": [1, 1, 1, 1],
            "boundary_names": ["outer", "outer", "outer", "outer"],
            "boundary_row_identity": [
                {"surface_triangle_index": 1, "surface_triangle_nodes": [1, 2, 3], "boundary_number": 1, "boundary_name": "outer"},
                {"surface_triangle_index": 2, "surface_triangle_nodes": [1, 4, 2], "boundary_number": 1, "boundary_name": "outer"},
                {"surface_triangle_index": 3, "surface_triangle_nodes": [1, 3, 4], "boundary_number": 1, "boundary_name": "outer"},
                {"surface_triangle_index": 4, "surface_triangle_nodes": [2, 4, 3], "boundary_number": 1, "boundary_name": "outer"},
            ],
            "source_file_id": "sha256:unit_tet_vol_source_abc123",
            "trace_output_artifact_id": "unit_tet_h1_to_scalar_bem_trace_output_p1",
            "trace_output_digest": "sha256:unit_tet_h1_to_scalar_bem_trace_output_p1",
            "trace_output_path": "unit_tet_h1_to_scalar_bem_trace_output.json",
            "trace_observable_id": "unit_tet_h1_to_scalar_bem_boundary_trace_v1",
            "trace_observable_family": "fem_bem_boundary_trace",
            "coupling_convention_schema_id": "matlab_first_order_fem_bem_coupling_convention_v1",
            "fem_bem_postprocess_row_convention_schema_id": "matlab_fem_bem_trace_lsq_row_convention_v1",
            "trace_basis_schema_id": "matlab_h1_p1_to_scalar_bem_p1_trace_basis_v1",
            "assembly_rule_id": "first_order_tet_h1_trace_tri_p1_bem_teaching_v1",
            "quadrature_rule_id": "tri_p1_exact_mass_regular_kernel_teaching_v1",
            "trace_row_identity": [
                {"trace_row_index": 1, "fem_node_id": 1, "bem_node_id": 1, "surface_node_index": 1},
                {"trace_row_index": 2, "fem_node_id": 2, "bem_node_id": 2, "surface_node_index": 2},
                {"trace_row_index": 3, "fem_node_id": 3, "bem_node_id": 3, "surface_node_index": 3},
                {"trace_row_index": 4, "fem_node_id": 4, "bem_node_id": 4, "surface_node_index": 4},
            ],
            "trace_matrix": [
                [1.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 1.0, 0.0],
            ],
        },
        "operators": {
            "coupled_system_artifact_id": "unit_tet_laplace_fem_bem_schur_system_v1",
            "coupled_system_digest": "sha256:unit_tet_laplace_fem_bem_schur_system_v1",
            "coupling_convention_schema_id": "matlab_first_order_fem_bem_coupling_convention_v1",
            "fem_bem_postprocess_row_convention_schema_id": "matlab_fem_bem_trace_lsq_row_convention_v1",
            "trace_basis_schema_id": "matlab_h1_p1_to_scalar_bem_p1_trace_basis_v1",
            "trace": {
                "boundary_row_identity": [
                    {"surface_triangle_index": 1, "surface_triangle_nodes": [1, 2, 3], "boundary_number": 1, "boundary_name": "outer"},
                    {"surface_triangle_index": 2, "surface_triangle_nodes": [1, 4, 2], "boundary_number": 1, "boundary_name": "outer"},
                    {"surface_triangle_index": 3, "surface_triangle_nodes": [1, 3, 4], "boundary_number": 1, "boundary_name": "outer"},
                    {"surface_triangle_index": 4, "surface_triangle_nodes": [2, 4, 3], "boundary_number": 1, "boundary_name": "outer"},
                ],
                "trace_row_identity": [
                    {"trace_row_index": 1, "fem_node_id": 1, "bem_node_id": 1, "surface_node_index": 1},
                    {"trace_row_index": 2, "fem_node_id": 2, "bem_node_id": 2, "surface_node_index": 2},
                    {"trace_row_index": 3, "fem_node_id": 3, "bem_node_id": 3, "surface_node_index": 3},
                    {"trace_row_index": 4, "fem_node_id": 4, "bem_node_id": 4, "surface_node_index": 4},
                ],
                "assembly_rule_id": "first_order_tet_h1_trace_tri_p1_bem_teaching_v1",
                "quadrature_rule_id": "tri_p1_exact_mass_regular_kernel_teaching_v1",
                "fem_bem_postprocess_row_convention_schema_id": "matlab_fem_bem_trace_lsq_row_convention_v1",
                "trace_basis_schema_id": "matlab_h1_p1_to_scalar_bem_p1_trace_basis_v1",
            },
            "bem": {
                "coupling_convention_schema_id": "matlab_first_order_fem_bem_coupling_convention_v1",
                "assembly_rule_id": "first_order_tet_h1_trace_tri_p1_bem_teaching_v1",
                "quadrature_rule_id": "tri_p1_exact_mass_regular_kernel_teaching_v1",
            },
        },
        "execution": {
            "resultArtifactId": "matlab_slot344_fem_bem_manifest_result_v1",
            "linearSolverReportArtifactId": "matlab_slot367_fem_bem_linear_solver_report_v1",
            "linearSolverReportDigest": "sha256:matlab_slot367_fem_bem_linear_solver_report_v1",
            "linearSolverName": "minimum_norm_pinv_rank_deficient",
            "linearSolverTolerance": 1.0e-10,
            "linearSolverResidualNorm": 2.0e-13,
            "linearSolverIterationCount": 1,
            "parameterSetArtifactId": "matlab_slot388_fem_bem_parameter_set_v1",
            "parameterSetDigest": "sha256:matlab-slot388-fem-bem-parameter-set-v1",
            "parameterSetPath": "docs/fem_bem/first_order_fem_bem_parameter_set.json",
            "objectiveObservableId": "matlab_slot388_trace_lsq_residual_objective_v1",
            "objectiveObservableFamily": "fem_bem_trace_least_squares_objective",
            "runStartedAt": "2026-07-01T13:50:00+09:00",
            "matlabVersion": "R2026a",
            "timingBreakdown": {
                "mesh_read_s": 0.001,
                "trace_assembly_s": 0.002,
                "manifest_build_s": 0.003,
                "json_write_s": 0.004,
            },
        },
        "optimization": {
            "parameterSetArtifactId": "matlab_slot388_fem_bem_parameter_set_v1",
            "parameterSetDigest": "sha256:matlab-slot388-fem-bem-parameter-set-v1",
            "parameterSetPath": "docs/fem_bem/first_order_fem_bem_parameter_set.json",
            "objectiveObservableId": "matlab_slot388_trace_lsq_residual_objective_v1",
            "objectiveObservableFamily": "fem_bem_trace_least_squares_objective",
        },
    }

    gate = netgen_vol_first_order_fem_bem_trace_package_handoff(
        package,
        expected_mesh_id="unit_tet_mesh",
        expected_export_id="coreform_netgen_unit_tet",
        expected_trace_artifact_id="unit_tet_h1_to_scalar_bem_trace_p1",
        expected_surface_mesh_id="unit_tet_boundary_tri_p1",
        expected_source_file_id="sha256:unit_tet_vol_source_abc123",
        expected_coupling_kind="h1_p1_to_scalar_bem_p1_trace",
        expected_formulation_id="laplace_single_layer_teaching",
        expected_bem_kernel_family="laplace_single_layer",
        expected_coupling_convention_schema_id="matlab_first_order_fem_bem_coupling_convention_v1",
        expected_fem_bem_postprocess_row_convention_schema_id="matlab_fem_bem_trace_lsq_row_convention_v1",
        expected_trace_basis_schema_id="matlab_h1_p1_to_scalar_bem_p1_trace_basis_v1",
        expected_assembly_rule_id="first_order_tet_h1_trace_tri_p1_bem_teaching_v1",
        expected_quadrature_rule_id="tri_p1_exact_mass_regular_kernel_teaching_v1",
        expected_volume_space="H1_P1_tet",
        expected_surface_space="scalar_P1_tri",
        expected_boundary_numbers=[1],
        expected_boundary_names=["outer"],
        expected_boundary_row_identity=package["trace"]["boundary_row_identity"],
        expected_trace_operator_artifact_id="unit_tet_h1_to_scalar_bem_trace_operator_p1",
        expected_trace_operator_policy="one_hot_boundary_node_injection_from_vol_node_ids",
        expected_trace_output_artifact_id="unit_tet_h1_to_scalar_bem_trace_output_p1",
        expected_trace_output_digest="sha256:unit_tet_h1_to_scalar_bem_trace_output_p1",
        expected_trace_observable_id="unit_tet_h1_to_scalar_bem_boundary_trace_v1",
        expected_trace_observable_family="fem_bem_boundary_trace",
        expected_coupled_system_artifact_id="unit_tet_laplace_fem_bem_schur_system_v1",
        expected_coupled_system_digest="sha256:unit_tet_laplace_fem_bem_schur_system_v1",
        expected_result_artifact_id="matlab_slot344_fem_bem_manifest_result_v1",
        expected_matlab_version="R2026a",
        expected_linear_solver_report_artifact_id="matlab_slot367_fem_bem_linear_solver_report_v1",
        expected_linear_solver_report_digest="sha256:matlab_slot367_fem_bem_linear_solver_report_v1",
        expected_linear_solver_name="minimum_norm_pinv_rank_deficient",
        expected_linear_solver_tolerance=1.0e-10,
        expected_linear_solver_residual_norm_max=1.0e-12,
        expected_parameter_set_artifact_id="matlab_slot388_fem_bem_parameter_set_v1",
        expected_parameter_set_digest="sha256:matlab-slot388-fem-bem-parameter-set-v1",
        expected_parameter_set_path="docs/fem_bem/first_order_fem_bem_parameter_set.json",
        expected_objective_observable_id="matlab_slot388_trace_lsq_residual_objective_v1",
        expected_objective_observable_family="fem_bem_trace_least_squares_objective",
        require_result_provenance=True,
        require_trace_output_artifact=True,
        require_linear_solver_report=True,
        require_parameter_set_artifact=True,
        require_coupling_convention_schema=True,
        require_fem_bem_postprocess_row_convention_schema=True,
        require_trace_basis_schema=True,
    )

    assert gate["status"] == "ok"
    assert gate["policy"] == "netgen_vol_first_order_fem_bem_trace_package_handoff"
    assert gate["trace_shape"] == [4, 5]
    assert gate["boundary_nodes_from_triangles"] == [1, 2, 3, 4]
    assert gate["checks"]["polynomial_order_first_order"] is True
    assert gate["checks"]["curvedelements_absent"] is True
    assert gate["checks"]["fem_bem_node_ids_identical"] is True
    assert gate["checks"]["trace_matrix_matches_fem_node_ids"] is True
    assert gate["checks"]["expected_export_id_matches"] is True
    assert gate["checks"]["trace_artifact_id_recorded"] is True
    assert gate["checks"]["trace_operator_artifact_id_recorded"] is True
    assert gate["checks"]["trace_operator_policy_recorded"] is True
    assert gate["checks"]["trace_output_artifact_id_recorded_when_required"] is True
    assert gate["checks"]["trace_output_digest_recorded_when_required"] is True
    assert gate["checks"]["trace_output_path_recorded_when_required"] is True
    assert gate["checks"]["trace_output_artifact_id_consistent_when_present"] is True
    assert gate["checks"]["trace_output_digest_consistent_when_present"] is True
    assert gate["checks"]["trace_observable_id_consistent_when_present"] is True
    assert gate["checks"]["trace_observable_family_consistent_when_present"] is True
    assert gate["checks"]["trace_output_path_consistent_when_present"] is True
    assert gate["checks"]["surface_mesh_id_recorded"] is True
    assert gate["checks"]["expected_trace_artifact_id_matches"] is True
    assert gate["checks"]["expected_trace_operator_artifact_id_matches"] is True
    assert gate["checks"]["expected_trace_operator_policy_matches"] is True
    assert gate["checks"]["expected_trace_output_artifact_id_matches"] is True
    assert gate["checks"]["expected_trace_output_digest_matches"] is True
    assert gate["checks"]["trace_observable_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_trace_observable_id_matches"] is True
    assert gate["checks"]["trace_observable_family_recorded_when_expected"] is True
    assert gate["checks"]["expected_trace_observable_family_matches"] is True
    assert gate["checks"]["coupled_system_artifact_id_consistent_when_present"] is True
    assert gate["checks"]["coupled_system_digest_consistent_when_present"] is True
    assert gate["checks"]["coupled_system_artifact_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_coupled_system_artifact_id_matches"] is True
    assert gate["checks"]["coupled_system_digest_recorded_when_expected"] is True
    assert gate["checks"]["expected_coupled_system_digest_matches"] is True
    assert gate["checks"]["coupling_convention_schema_id_consistent_when_present"] is True
    assert gate["checks"]["coupling_convention_schema_id_recorded_when_required"] is True
    assert gate["checks"]["coupling_convention_schema_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_coupling_convention_schema_id_matches"] is True
    assert gate["coupling_convention_schema_id"] == "matlab_first_order_fem_bem_coupling_convention_v1"
    assert gate["coupling_convention_schema_ids"] == ["matlab_first_order_fem_bem_coupling_convention_v1"]
    assert gate["require_coupling_convention_schema"] is True
    assert gate["checks"]["fem_bem_postprocess_row_convention_schema_id_consistent_when_present"] is True
    assert gate["checks"]["fem_bem_postprocess_row_convention_schema_id_recorded_when_required"] is True
    assert gate["checks"]["fem_bem_postprocess_row_convention_schema_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_fem_bem_postprocess_row_convention_schema_id_matches"] is True
    assert gate["fem_bem_postprocess_row_convention_schema_id"] == "matlab_fem_bem_trace_lsq_row_convention_v1"
    assert gate["fem_bem_postprocess_row_convention_schema_ids"] == [
        "matlab_fem_bem_trace_lsq_row_convention_v1"
    ]
    assert gate["require_fem_bem_postprocess_row_convention_schema"] is True
    assert gate["checks"]["trace_basis_schema_id_consistent_when_present"] is True
    assert gate["checks"]["trace_basis_schema_id_recorded_when_required"] is True
    assert gate["checks"]["trace_basis_schema_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_trace_basis_schema_id_matches"] is True
    assert gate["trace_basis_schema_id"] == "matlab_h1_p1_to_scalar_bem_p1_trace_basis_v1"
    assert gate["trace_basis_schema_ids"] == [
        "matlab_h1_p1_to_scalar_bem_p1_trace_basis_v1"
    ]
    assert gate["require_trace_basis_schema"] is True
    assert gate["checks"]["linear_solver_report_artifact_id_consistent_when_present"] is True
    assert gate["checks"]["linear_solver_report_digest_consistent_when_present"] is True
    assert gate["checks"]["linear_solver_name_consistent_when_present"] is True
    assert gate["checks"]["linear_solver_tolerance_consistent_when_present"] is True
    assert gate["checks"]["linear_solver_residual_norm_consistent_when_present"] is True
    assert gate["checks"]["linear_solver_iteration_count_consistent_when_present"] is True
    assert gate["checks"]["linear_solver_report_artifact_id_recorded_when_required"] is True
    assert gate["checks"]["linear_solver_report_digest_recorded_when_required"] is True
    assert gate["checks"]["linear_solver_name_recorded_when_required"] is True
    assert gate["checks"]["linear_solver_tolerance_recorded_when_required"] is True
    assert gate["checks"]["linear_solver_residual_norm_recorded_when_required"] is True
    assert gate["checks"]["linear_solver_report_artifact_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_linear_solver_report_artifact_id_matches"] is True
    assert gate["checks"]["linear_solver_report_digest_recorded_when_expected"] is True
    assert gate["checks"]["expected_linear_solver_report_digest_matches"] is True
    assert gate["checks"]["linear_solver_name_recorded_when_expected"] is True
    assert gate["checks"]["expected_linear_solver_name_matches"] is True
    assert gate["checks"]["linear_solver_tolerance_recorded_when_expected"] is True
    assert gate["checks"]["expected_linear_solver_tolerance_matches"] is True
    assert gate["checks"]["linear_solver_residual_norm_recorded_when_expected"] is True
    assert gate["checks"]["linear_solver_residual_norm_below_expected_max"] is True
    assert gate["checks"]["linear_solver_iteration_count_nonnegative_when_present"] is True
    assert gate["checks"]["parameter_set_artifact_id_consistent_when_present"] is True
    assert gate["checks"]["parameter_set_digest_consistent_when_present"] is True
    assert gate["checks"]["parameter_set_path_consistent_when_present"] is True
    assert gate["checks"]["objective_observable_id_consistent_when_present"] is True
    assert gate["checks"]["objective_observable_family_consistent_when_present"] is True
    assert gate["checks"]["parameter_set_artifact_id_recorded_when_required"] is True
    assert gate["checks"]["parameter_set_digest_recorded_when_required"] is True
    assert gate["checks"]["parameter_set_path_recorded_when_required"] is True
    assert gate["checks"]["parameter_set_artifact_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_parameter_set_artifact_id_matches"] is True
    assert gate["checks"]["parameter_set_digest_recorded_when_expected"] is True
    assert gate["checks"]["expected_parameter_set_digest_matches"] is True
    assert gate["checks"]["parameter_set_path_recorded_when_expected"] is True
    assert gate["checks"]["expected_parameter_set_path_matches"] is True
    assert gate["checks"]["objective_observable_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_objective_observable_id_matches"] is True
    assert gate["checks"]["objective_observable_family_recorded_when_expected"] is True
    assert gate["checks"]["expected_objective_observable_family_matches"] is True
    assert gate["checks"]["result_artifact_id_recorded_when_required"] is True
    assert gate["checks"]["expected_result_artifact_id_matches"] is True
    assert gate["checks"]["run_started_at_recorded_when_required"] is True
    assert gate["checks"]["matlab_version_recorded_when_required"] is True
    assert gate["checks"]["expected_matlab_version_matches"] is True
    assert gate["checks"]["timing_breakdown_recorded_when_required"] is True
    assert gate["checks"]["timing_breakdown_has_at_least_four_items"] is True
    assert gate["checks"]["timing_breakdown_values_finite_nonnegative"] is True
    assert gate["checks"]["expected_surface_mesh_id_matches"] is True
    assert gate["checks"]["source_file_id_consistent_when_present"] is True
    assert gate["checks"]["source_file_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_source_file_id_matches"] is True
    assert gate["checks"]["expected_trace_source_file_id_matches"] is True
    assert gate["checks"]["expected_coupling_kind_matches"] is True
    assert gate["checks"]["expected_formulation_id_matches"] is True
    assert gate["checks"]["expected_bem_kernel_family_matches"] is True
    assert gate["checks"]["assembly_rule_id_consistent_when_present"] is True
    assert gate["checks"]["quadrature_rule_id_consistent_when_present"] is True
    assert gate["checks"]["assembly_rule_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_assembly_rule_id_matches"] is True
    assert gate["checks"]["quadrature_rule_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_quadrature_rule_id_matches"] is True
    assert gate["checks"]["expected_volume_space_matches"] is True
    assert gate["checks"]["expected_surface_space_matches"] is True
    assert gate["checks"]["boundary_numbers_match_surface_triangles"] is True
    assert gate["checks"]["boundary_numbers_positive"] is True
    assert gate["checks"]["boundary_names_match_surface_triangles"] is True
    assert gate["checks"]["boundary_names_recorded"] is True
    assert gate["checks"]["boundary_number_name_pairs_consistent"] is True
    assert gate["checks"]["boundary_numbers_recorded_when_expected"] is True
    assert gate["checks"]["expected_boundary_numbers_match"] is True
    assert gate["checks"]["boundary_names_recorded_when_expected"] is True
    assert gate["checks"]["expected_boundary_names_match"] is True
    assert gate["checks"]["boundary_row_identity_rows_match_surface_triangles"] is True
    assert gate["checks"]["boundary_row_identity_matches_surface_triangles"] is True
    assert gate["checks"]["boundary_row_identity_boundary_numbers_match"] is True
    assert gate["checks"]["boundary_row_identity_boundary_names_match"] is True
    assert gate["checks"]["operator_boundary_row_identity_matches_trace_identity"] is True
    assert gate["checks"]["boundary_row_identity_recorded_when_expected"] is True
    assert gate["checks"]["expected_boundary_row_identity_matches"] is True
    assert gate["trace_row_identity_present"] is True
    assert gate["checks"]["trace_row_identity_rows_match_trace_rows"] is True
    assert gate["checks"]["trace_row_identity_row_indices_match"] is True
    assert gate["checks"]["trace_row_identity_fem_nodes_match"] is True
    assert gate["checks"]["trace_row_identity_bem_nodes_match"] is True
    assert gate["checks"]["trace_row_identity_matches_trace_matrix"] is True
    assert gate["operator_trace_row_identity_present"] is True
    assert gate["boundary_row_identity_present"] is True
    assert gate["operator_boundary_row_identity_present"] is True
    assert gate["checks"]["operator_trace_row_identity_rows_match_trace_rows"] is True
    assert gate["checks"]["operator_trace_row_identity_row_indices_match"] is True
    assert gate["checks"]["operator_trace_row_identity_fem_nodes_match"] is True
    assert gate["checks"]["operator_trace_row_identity_bem_nodes_match"] is True
    assert gate["checks"]["operator_trace_row_identity_matches_trace_matrix"] is True
    assert gate["checks"]["operator_trace_row_identity_matches_trace_identity"] is True
    assert gate["trace_row_identity_mismatch_rows"] == []
    assert gate["operator_trace_row_identity_mismatch_rows"] == []
    assert gate["boundary_row_identity_mismatch_rows"] == []
    assert gate["operator_boundary_row_identity_mismatch_rows"] == []
    assert gate["source_file_id"] == "sha256:unit_tet_vol_source_abc123"
    assert gate["trace_source_file_id"] == "sha256:unit_tet_vol_source_abc123"
    assert gate["trace_operator_artifact_id"] == "unit_tet_h1_to_scalar_bem_trace_operator_p1"
    assert gate["trace_operator_policy"] == "one_hot_boundary_node_injection_from_vol_node_ids"
    assert gate["trace_output_artifact_id"] == "unit_tet_h1_to_scalar_bem_trace_output_p1"
    assert gate["trace_output_digest"] == "sha256:unit_tet_h1_to_scalar_bem_trace_output_p1"
    assert gate["trace_output_path"] == "unit_tet_h1_to_scalar_bem_trace_output.json"
    assert gate["trace_observable_id"] == "unit_tet_h1_to_scalar_bem_boundary_trace_v1"
    assert gate["trace_observable_family"] == "fem_bem_boundary_trace"
    assert gate["coupled_system_artifact_id"] == "unit_tet_laplace_fem_bem_schur_system_v1"
    assert gate["coupled_system_digest"] == "sha256:unit_tet_laplace_fem_bem_schur_system_v1"
    assert gate["linear_solver_report_artifact_id"] == "matlab_slot367_fem_bem_linear_solver_report_v1"
    assert gate["linear_solver_report_digest"] == "sha256:matlab_slot367_fem_bem_linear_solver_report_v1"
    assert gate["linear_solver_name"] == "minimum_norm_pinv_rank_deficient"
    assert gate["linear_solver_tolerance"] == 1.0e-10
    assert gate["linear_solver_residual_norm"] == 2.0e-13
    assert gate["linear_solver_iteration_count"] == 1
    assert gate["linear_solver_report_artifact_ids"] == ["matlab_slot367_fem_bem_linear_solver_report_v1"]
    assert gate["linear_solver_report_digests"] == ["sha256:matlab_slot367_fem_bem_linear_solver_report_v1"]
    assert gate["parameter_set_artifact_id"] == "matlab_slot388_fem_bem_parameter_set_v1"
    assert gate["parameter_set_digest"] == "sha256:matlab-slot388-fem-bem-parameter-set-v1"
    assert gate["parameter_set_path"] == "docs/fem_bem/first_order_fem_bem_parameter_set.json"
    assert gate["objective_observable_id"] == "matlab_slot388_trace_lsq_residual_objective_v1"
    assert gate["objective_observable_family"] == "fem_bem_trace_least_squares_objective"
    assert gate["result_artifact_id"] == "matlab_slot344_fem_bem_manifest_result_v1"
    assert gate["run_started_at"] == "2026-07-01T13:50:00+09:00"
    assert gate["matlab_version"] == "R2026a"
    assert gate["timing_breakdown_names"] == [
        "json_write_s",
        "manifest_build_s",
        "mesh_read_s",
        "trace_assembly_s",
    ]
    assert gate["timing_breakdown_seconds"]["trace_assembly_s"] == 0.002
    assert gate["assembly_rule_id"] == "first_order_tet_h1_trace_tri_p1_bem_teaching_v1"
    assert gate["quadrature_rule_id"] == "tri_p1_exact_mass_regular_kernel_teaching_v1"
    assert gate["boundary_numbers"] == [1, 1, 1, 1]
    assert gate["boundary_names"] == ["outer", "outer", "outer", "outer"]
    assert gate["boundary_row_identity"][1]["surface_triangle_nodes"] == [1, 4, 2]

    matlab_trace = dict(package["trace"], boundary_names="outer")
    matlab_gypsilab = dict(package["gypsilab"], boundary_names="outer")
    matlab_json_style = dict(package)
    matlab_json_style["trace"] = [matlab_trace, dict(matlab_trace)]
    matlab_json_style["gypsilab"] = [matlab_gypsilab, dict(matlab_gypsilab)]
    matlab_json_style["operators"] = {"trace": [dict(package["operators"]["trace"])]}
    matlab_json_gate = netgen_vol_first_order_fem_bem_trace_package_handoff(
        matlab_json_style,
        expected_trace_operator_artifact_id="unit_tet_h1_to_scalar_bem_trace_operator_p1",
        expected_trace_operator_policy="one_hot_boundary_node_injection_from_vol_node_ids",
    )
    assert matlab_json_gate["status"] == "ok"
    assert matlab_json_gate["operator_trace_row_identity_present"] is True
    assert matlab_json_gate["boundary_names"] == ["outer", "outer", "outer", "outer"]

    source_less = dict(package)
    source_less.pop("source_path")
    missing_source = netgen_vol_first_order_fem_bem_trace_package_handoff(source_less)
    assert missing_source["status"] == "needs_attention"
    assert missing_source["checks"]["source_path_recorded"] is False

    swapped_bem = dict(package)
    swapped_bem["trace"] = dict(package["trace"], bem_node_ids=[1, 2, 4, 3])
    swapped = netgen_vol_first_order_fem_bem_trace_package_handoff(swapped_bem)
    assert swapped["status"] == "needs_attention"
    assert swapped["checks"]["fem_bem_node_ids_identical"] is False

    trace_id_less = dict(package)
    trace_id_less.pop("trace_artifact_id")
    trace_id_less_gate = netgen_vol_first_order_fem_bem_trace_package_handoff(trace_id_less)
    assert trace_id_less_gate["status"] == "needs_attention"
    assert trace_id_less_gate["checks"]["trace_artifact_id_recorded"] is False

    trace_operator_id_less = dict(package)
    trace_operator_id_less.pop("trace_operator_artifact_id")
    trace_operator_id_less_gate = netgen_vol_first_order_fem_bem_trace_package_handoff(trace_operator_id_less)
    assert trace_operator_id_less_gate["status"] == "needs_attention"
    assert trace_operator_id_less_gate["checks"]["trace_operator_artifact_id_recorded"] is False

    stale_trace_operator = dict(package, trace_operator_artifact_id="stale_trace_operator")
    stale_trace_operator_gate = netgen_vol_first_order_fem_bem_trace_package_handoff(
        stale_trace_operator,
        expected_trace_operator_artifact_id="unit_tet_h1_to_scalar_bem_trace_operator_p1",
    )
    assert stale_trace_operator_gate["status"] == "needs_attention"
    assert stale_trace_operator_gate["checks"]["expected_trace_operator_artifact_id_matches"] is False

    wrong_trace_operator_policy = dict(package, trace_operator_policy="remote_field_observation_map")
    wrong_trace_operator_policy_gate = netgen_vol_first_order_fem_bem_trace_package_handoff(
        wrong_trace_operator_policy,
        expected_trace_operator_policy="one_hot_boundary_node_injection_from_vol_node_ids",
    )
    assert wrong_trace_operator_policy_gate["status"] == "needs_attention"
    assert wrong_trace_operator_policy_gate["checks"]["expected_trace_operator_policy_matches"] is False

    stale_trace_output = dict(package, trace_output_artifact_id="stale_trace_output")
    stale_trace_output_gate = netgen_vol_first_order_fem_bem_trace_package_handoff(
        stale_trace_output,
        expected_trace_output_artifact_id="unit_tet_h1_to_scalar_bem_trace_output_p1",
    )
    assert stale_trace_output_gate["status"] == "needs_attention"
    assert stale_trace_output_gate["checks"]["trace_output_artifact_id_consistent_when_present"] is False
    assert stale_trace_output_gate["checks"]["expected_trace_output_artifact_id_matches"] is False

    stale_trace_output_digest = dict(package, trace_output_digest="sha256:old_trace_output")
    stale_trace_output_digest_gate = netgen_vol_first_order_fem_bem_trace_package_handoff(
        stale_trace_output_digest,
        expected_trace_output_digest="sha256:unit_tet_h1_to_scalar_bem_trace_output_p1",
    )
    assert stale_trace_output_digest_gate["status"] == "needs_attention"
    assert stale_trace_output_digest_gate["checks"]["trace_output_digest_consistent_when_present"] is False
    assert stale_trace_output_digest_gate["checks"]["expected_trace_output_digest_matches"] is False

    missing_trace_output_path = dict(package)
    missing_trace_output_path.pop("trace_output_path")
    missing_trace_output_path["trace"] = dict(package["trace"])
    missing_trace_output_path["trace"].pop("trace_output_path")
    missing_trace_output_path_gate = netgen_vol_first_order_fem_bem_trace_package_handoff(
        missing_trace_output_path,
        require_trace_output_artifact=True,
    )
    assert missing_trace_output_path_gate["status"] == "needs_attention"
    assert missing_trace_output_path_gate["checks"]["trace_output_path_recorded_when_required"] is False

    stale_trace_observable = dict(package, trace_observable_id="stale_remote_field_map_observable")
    stale_trace_observable_gate = netgen_vol_first_order_fem_bem_trace_package_handoff(
        stale_trace_observable,
        expected_trace_observable_id="unit_tet_h1_to_scalar_bem_boundary_trace_v1",
        expected_trace_observable_family="fem_bem_boundary_trace",
        expected_trace_output_artifact_id="unit_tet_h1_to_scalar_bem_trace_output_p1",
        expected_trace_output_digest="sha256:unit_tet_h1_to_scalar_bem_trace_output_p1",
    )
    assert stale_trace_observable_gate["status"] == "needs_attention"
    assert stale_trace_observable_gate["checks"]["trace_observable_id_consistent_when_present"] is False
    assert stale_trace_observable_gate["checks"]["expected_trace_observable_id_matches"] is False
    assert stale_trace_observable_gate["checks"]["expected_trace_observable_family_matches"] is True
    assert stale_trace_observable_gate["checks"]["expected_trace_output_artifact_id_matches"] is True

    wrong_trace_observable_family = dict(package)
    wrong_trace_observable_family["trace"] = dict(package["trace"], trace_observable_family="remote_field_map")
    wrong_trace_observable_family_gate = netgen_vol_first_order_fem_bem_trace_package_handoff(
        wrong_trace_observable_family,
        expected_trace_observable_id="unit_tet_h1_to_scalar_bem_boundary_trace_v1",
        expected_trace_observable_family="fem_bem_boundary_trace",
        expected_trace_output_artifact_id="unit_tet_h1_to_scalar_bem_trace_output_p1",
        expected_trace_output_digest="sha256:unit_tet_h1_to_scalar_bem_trace_output_p1",
    )
    assert wrong_trace_observable_family_gate["status"] == "needs_attention"
    assert wrong_trace_observable_family_gate["checks"]["expected_trace_observable_id_matches"] is True
    assert wrong_trace_observable_family_gate["checks"]["trace_observable_family_consistent_when_present"] is False
    assert wrong_trace_observable_family_gate["checks"]["expected_trace_observable_family_matches"] is False
    assert wrong_trace_observable_family_gate["checks"]["expected_trace_output_artifact_id_matches"] is True

    stale_coupled_system = dict(package, coupled_system_artifact_id="stale_schur_system_v0")
    stale_coupled_system["operators"] = dict(
        package["operators"],
        coupled_system_artifact_id="stale_schur_system_v0",
    )
    stale_coupled_system_gate = netgen_vol_first_order_fem_bem_trace_package_handoff(
        stale_coupled_system,
        expected_coupled_system_artifact_id="unit_tet_laplace_fem_bem_schur_system_v1",
        expected_coupled_system_digest="sha256:unit_tet_laplace_fem_bem_schur_system_v1",
    )
    assert stale_coupled_system_gate["status"] == "needs_attention"
    assert stale_coupled_system_gate["checks"]["coupled_system_artifact_id_consistent_when_present"] is True
    assert stale_coupled_system_gate["checks"]["expected_coupled_system_artifact_id_matches"] is False
    assert stale_coupled_system_gate["checks"]["expected_coupled_system_digest_matches"] is True

    missing_coupled_system_digest = dict(package)
    missing_coupled_system_digest.pop("coupled_system_digest")
    missing_coupled_system_digest["operators"] = dict(package["operators"])
    missing_coupled_system_digest["operators"].pop("coupled_system_digest")
    missing_coupled_system_digest_gate = netgen_vol_first_order_fem_bem_trace_package_handoff(
        missing_coupled_system_digest,
        expected_coupled_system_digest="sha256:unit_tet_laplace_fem_bem_schur_system_v1",
    )
    assert missing_coupled_system_digest_gate["status"] == "needs_attention"
    assert missing_coupled_system_digest_gate["checks"]["coupled_system_digest_recorded_when_expected"] is False

    stale_solver_report_digest = dict(package, linear_solver_report_digest="sha256:old_solver_report")
    stale_solver_report_digest["execution"] = dict(
        package["execution"],
        linearSolverReportDigest="sha256:old_solver_report",
    )
    stale_solver_report_digest_gate = netgen_vol_first_order_fem_bem_trace_package_handoff(
        stale_solver_report_digest,
        expected_linear_solver_report_digest="sha256:matlab_slot367_fem_bem_linear_solver_report_v1",
        require_linear_solver_report=True,
    )
    assert stale_solver_report_digest_gate["status"] == "needs_attention"
    assert stale_solver_report_digest_gate["checks"]["linear_solver_report_digest_consistent_when_present"] is True
    assert stale_solver_report_digest_gate["checks"]["expected_linear_solver_report_digest_matches"] is False

    missing_solver_report = dict(package)
    for key in (
        "linear_solver_report_artifact_id",
        "linear_solver_report_digest",
        "linear_solver_name",
        "linear_solver_tolerance",
        "linear_solver_residual_norm",
        "linear_solver_iteration_count",
    ):
        missing_solver_report.pop(key)
    missing_solver_report["execution"] = dict(package["execution"])
    for key in (
        "linearSolverReportArtifactId",
        "linearSolverReportDigest",
        "linearSolverName",
        "linearSolverTolerance",
        "linearSolverResidualNorm",
        "linearSolverIterationCount",
    ):
        missing_solver_report["execution"].pop(key)
    missing_solver_report_gate = netgen_vol_first_order_fem_bem_trace_package_handoff(
        missing_solver_report,
        require_linear_solver_report=True,
    )
    assert missing_solver_report_gate["status"] == "needs_attention"
    assert missing_solver_report_gate["checks"]["linear_solver_report_artifact_id_recorded_when_required"] is False
    assert missing_solver_report_gate["checks"]["linear_solver_report_digest_recorded_when_required"] is False
    assert missing_solver_report_gate["checks"]["linear_solver_residual_norm_recorded_when_required"] is False

    high_residual_solver_report = dict(package, linear_solver_residual_norm=2.0e-6)
    high_residual_solver_report["execution"] = dict(
        package["execution"],
        linearSolverResidualNorm=2.0e-6,
    )
    high_residual_solver_report_gate = netgen_vol_first_order_fem_bem_trace_package_handoff(
        high_residual_solver_report,
        expected_linear_solver_residual_norm_max=1.0e-12,
        require_linear_solver_report=True,
    )
    assert high_residual_solver_report_gate["status"] == "needs_attention"
    assert high_residual_solver_report_gate["checks"]["linear_solver_residual_norm_below_expected_max"] is False

    wrong_solver_name = dict(package, linear_solver_name="direct_lu")
    wrong_solver_name["execution"] = dict(package["execution"], linearSolverName="direct_lu")
    wrong_solver_name_gate = netgen_vol_first_order_fem_bem_trace_package_handoff(
        wrong_solver_name,
        expected_linear_solver_name="minimum_norm_pinv_rank_deficient",
        require_linear_solver_report=True,
    )
    assert wrong_solver_name_gate["status"] == "needs_attention"
    assert wrong_solver_name_gate["checks"]["expected_linear_solver_name_matches"] is False

    stale_parameter_digest = dict(package, parameter_set_digest="sha256:old-parameter-set")
    stale_parameter_digest["execution"] = dict(
        package["execution"],
        parameterSetDigest="sha256:old-parameter-set",
    )
    stale_parameter_digest["optimization"] = dict(
        package["optimization"],
        parameterSetDigest="sha256:old-parameter-set",
    )
    stale_parameter_digest_gate = netgen_vol_first_order_fem_bem_trace_package_handoff(
        stale_parameter_digest,
        expected_parameter_set_digest="sha256:matlab-slot388-fem-bem-parameter-set-v1",
        require_parameter_set_artifact=True,
    )
    assert stale_parameter_digest_gate["status"] == "needs_attention"
    assert stale_parameter_digest_gate["checks"]["parameter_set_digest_consistent_when_present"] is True
    assert stale_parameter_digest_gate["checks"]["expected_parameter_set_digest_matches"] is False

    missing_parameter_path = dict(package)
    missing_parameter_path.pop("parameter_set_path")
    missing_parameter_path["execution"] = dict(package["execution"])
    missing_parameter_path["execution"].pop("parameterSetPath")
    missing_parameter_path["optimization"] = dict(package["optimization"])
    missing_parameter_path["optimization"].pop("parameterSetPath")
    missing_parameter_path_gate = netgen_vol_first_order_fem_bem_trace_package_handoff(
        missing_parameter_path,
        require_parameter_set_artifact=True,
    )
    assert missing_parameter_path_gate["status"] == "needs_attention"
    assert missing_parameter_path_gate["checks"]["parameter_set_artifact_id_recorded_when_required"] is True
    assert missing_parameter_path_gate["checks"]["parameter_set_digest_recorded_when_required"] is True
    assert missing_parameter_path_gate["checks"]["parameter_set_path_recorded_when_required"] is False

    wrong_objective_family = dict(package, objective_observable_family="remote_field_map")
    wrong_objective_family["execution"] = dict(
        package["execution"],
        objectiveObservableFamily="remote_field_map",
    )
    wrong_objective_family["optimization"] = dict(
        package["optimization"],
        objectiveObservableFamily="remote_field_map",
    )
    wrong_objective_family_gate = netgen_vol_first_order_fem_bem_trace_package_handoff(
        wrong_objective_family,
        expected_objective_observable_family="fem_bem_trace_least_squares_objective",
    )
    assert wrong_objective_family_gate["status"] == "needs_attention"
    assert wrong_objective_family_gate["checks"]["objective_observable_family_consistent_when_present"] is True
    assert wrong_objective_family_gate["checks"]["expected_objective_observable_family_matches"] is False

    stale_coupling_convention_schema = dict(
        package,
        coupling_convention_schema_id="matlab_fem_bem_value_only_convention_v0",
    )
    stale_coupling_convention_schema["trace"] = dict(
        package["trace"],
        coupling_convention_schema_id="matlab_fem_bem_value_only_convention_v0",
    )
    stale_coupling_convention_schema["operators"] = dict(
        package["operators"],
        coupling_convention_schema_id="matlab_fem_bem_value_only_convention_v0",
    )
    stale_coupling_convention_schema["operators"]["bem"] = dict(
        package["operators"]["bem"],
        coupling_convention_schema_id="matlab_fem_bem_value_only_convention_v0",
    )
    stale_coupling_convention_schema_gate = netgen_vol_first_order_fem_bem_trace_package_handoff(
        stale_coupling_convention_schema,
        expected_coupling_kind="h1_p1_to_scalar_bem_p1_trace",
        expected_coupling_convention_schema_id="matlab_first_order_fem_bem_coupling_convention_v1",
        require_coupling_convention_schema=True,
    )
    assert stale_coupling_convention_schema_gate["status"] == "needs_attention"
    assert stale_coupling_convention_schema_gate["checks"]["coupling_convention_schema_id_consistent_when_present"] is True
    assert stale_coupling_convention_schema_gate["checks"]["coupling_convention_schema_id_recorded_when_required"] is True
    assert stale_coupling_convention_schema_gate["checks"]["expected_coupling_convention_schema_id_matches"] is False
    assert stale_coupling_convention_schema_gate["checks"]["expected_coupling_kind_matches"] is True

    missing_coupling_convention_schema = dict(package)
    missing_coupling_convention_schema.pop("coupling_convention_schema_id")
    missing_coupling_convention_schema["trace"] = dict(package["trace"])
    missing_coupling_convention_schema["trace"].pop("coupling_convention_schema_id")
    missing_coupling_convention_schema["operators"] = dict(package["operators"])
    missing_coupling_convention_schema["operators"].pop("coupling_convention_schema_id")
    missing_coupling_convention_schema["operators"]["bem"] = dict(package["operators"]["bem"])
    missing_coupling_convention_schema["operators"]["bem"].pop("coupling_convention_schema_id")
    missing_coupling_convention_schema_gate = netgen_vol_first_order_fem_bem_trace_package_handoff(
        missing_coupling_convention_schema,
        expected_coupling_convention_schema_id="matlab_first_order_fem_bem_coupling_convention_v1",
        require_coupling_convention_schema=True,
    )
    assert missing_coupling_convention_schema_gate["status"] == "needs_attention"
    assert missing_coupling_convention_schema_gate["checks"]["coupling_convention_schema_id_recorded_when_required"] is False
    assert missing_coupling_convention_schema_gate["checks"]["coupling_convention_schema_id_recorded_when_expected"] is False

    stale_postprocess_row_convention_schema = dict(
        package,
        fem_bem_postprocess_row_convention_schema_id="matlab_fem_bem_scalar_residual_row_v0",
    )
    stale_postprocess_row_convention_schema["trace"] = dict(
        package["trace"],
        fem_bem_postprocess_row_convention_schema_id="matlab_fem_bem_scalar_residual_row_v0",
    )
    stale_postprocess_row_convention_schema["operators"] = dict(
        package["operators"],
        fem_bem_postprocess_row_convention_schema_id="matlab_fem_bem_scalar_residual_row_v0",
    )
    stale_postprocess_row_convention_schema["operators"]["trace"] = dict(
        package["operators"]["trace"],
        fem_bem_postprocess_row_convention_schema_id="matlab_fem_bem_scalar_residual_row_v0",
    )
    stale_postprocess_row_convention_schema_gate = netgen_vol_first_order_fem_bem_trace_package_handoff(
        stale_postprocess_row_convention_schema,
        expected_coupling_convention_schema_id="matlab_first_order_fem_bem_coupling_convention_v1",
        expected_fem_bem_postprocess_row_convention_schema_id="matlab_fem_bem_trace_lsq_row_convention_v1",
        require_coupling_convention_schema=True,
        require_fem_bem_postprocess_row_convention_schema=True,
    )
    assert stale_postprocess_row_convention_schema_gate["status"] == "needs_attention"
    assert stale_postprocess_row_convention_schema_gate["checks"]["expected_coupling_convention_schema_id_matches"] is True
    assert stale_postprocess_row_convention_schema_gate["checks"]["fem_bem_postprocess_row_convention_schema_id_recorded_when_required"] is True
    assert (
        stale_postprocess_row_convention_schema_gate["checks"][
            "expected_fem_bem_postprocess_row_convention_schema_id_matches"
        ]
        is False
    )
    assert stale_postprocess_row_convention_schema_gate["checks"]["trace_matrix_matches_fem_node_ids"] is True

    missing_postprocess_row_convention_schema = dict(package)
    missing_postprocess_row_convention_schema.pop("fem_bem_postprocess_row_convention_schema_id")
    missing_postprocess_row_convention_schema["trace"] = dict(package["trace"])
    missing_postprocess_row_convention_schema["trace"].pop("fem_bem_postprocess_row_convention_schema_id")
    missing_postprocess_row_convention_schema["operators"] = dict(package["operators"])
    missing_postprocess_row_convention_schema["operators"].pop("fem_bem_postprocess_row_convention_schema_id")
    missing_postprocess_row_convention_schema["operators"]["trace"] = dict(package["operators"]["trace"])
    missing_postprocess_row_convention_schema["operators"]["trace"].pop("fem_bem_postprocess_row_convention_schema_id")
    missing_postprocess_row_convention_schema_gate = netgen_vol_first_order_fem_bem_trace_package_handoff(
        missing_postprocess_row_convention_schema,
        expected_fem_bem_postprocess_row_convention_schema_id="matlab_fem_bem_trace_lsq_row_convention_v1",
        require_fem_bem_postprocess_row_convention_schema=True,
    )
    assert missing_postprocess_row_convention_schema_gate["status"] == "needs_attention"
    assert (
        missing_postprocess_row_convention_schema_gate["checks"][
            "fem_bem_postprocess_row_convention_schema_id_recorded_when_required"
        ]
        is False
    )
    assert (
        missing_postprocess_row_convention_schema_gate["checks"][
            "fem_bem_postprocess_row_convention_schema_id_recorded_when_expected"
        ]
        is False
    )

    stale_trace_basis_schema = dict(
        package,
        trace_basis_schema_id="matlab_h1_trace_basis_value_only_v0",
    )
    stale_trace_basis_schema["trace"] = dict(
        package["trace"],
        trace_basis_schema_id="matlab_h1_trace_basis_value_only_v0",
    )
    stale_trace_basis_schema["operators"] = dict(
        package["operators"],
        trace_basis_schema_id="matlab_h1_trace_basis_value_only_v0",
    )
    stale_trace_basis_schema["operators"]["trace"] = dict(
        package["operators"]["trace"],
        trace_basis_schema_id="matlab_h1_trace_basis_value_only_v0",
    )
    stale_trace_basis_schema_gate = netgen_vol_first_order_fem_bem_trace_package_handoff(
        stale_trace_basis_schema,
        expected_trace_basis_schema_id="matlab_h1_p1_to_scalar_bem_p1_trace_basis_v1",
        require_trace_basis_schema=True,
    )
    assert stale_trace_basis_schema_gate["status"] == "needs_attention"
    assert stale_trace_basis_schema_gate["checks"]["trace_basis_schema_id_consistent_when_present"] is True
    assert stale_trace_basis_schema_gate["checks"]["trace_basis_schema_id_recorded_when_required"] is True
    assert stale_trace_basis_schema_gate["checks"]["expected_trace_basis_schema_id_matches"] is False
    assert stale_trace_basis_schema_gate["checks"]["trace_matrix_matches_fem_node_ids"] is True

    missing_trace_basis_schema = dict(package)
    missing_trace_basis_schema.pop("trace_basis_schema_id")
    missing_trace_basis_schema["trace"] = dict(package["trace"])
    missing_trace_basis_schema["trace"].pop("trace_basis_schema_id")
    missing_trace_basis_schema["operators"] = dict(package["operators"])
    missing_trace_basis_schema["operators"].pop("trace_basis_schema_id")
    missing_trace_basis_schema["operators"]["trace"] = dict(package["operators"]["trace"])
    missing_trace_basis_schema["operators"]["trace"].pop("trace_basis_schema_id")
    missing_trace_basis_schema_gate = netgen_vol_first_order_fem_bem_trace_package_handoff(
        missing_trace_basis_schema,
        expected_trace_basis_schema_id="matlab_h1_p1_to_scalar_bem_p1_trace_basis_v1",
        require_trace_basis_schema=True,
    )
    assert missing_trace_basis_schema_gate["status"] == "needs_attention"
    assert missing_trace_basis_schema_gate["checks"]["trace_basis_schema_id_recorded_when_required"] is False
    assert missing_trace_basis_schema_gate["checks"]["trace_basis_schema_id_recorded_when_expected"] is False

    wrong_surface_id = netgen_vol_first_order_fem_bem_trace_package_handoff(
        package,
        expected_surface_mesh_id="stale_boundary_tri_p1",
    )
    assert wrong_surface_id["status"] == "needs_attention"
    assert wrong_surface_id["checks"]["expected_surface_mesh_id_matches"] is False

    stale_source = dict(package)
    stale_source["trace"] = dict(package["trace"], source_file_id="sha256:old_unit_tet_vol_source")
    stale_source_gate = netgen_vol_first_order_fem_bem_trace_package_handoff(
        stale_source,
        expected_source_file_id="sha256:unit_tet_vol_source_abc123",
    )
    assert stale_source_gate["status"] == "needs_attention"
    assert stale_source_gate["checks"]["source_file_id_consistent_when_present"] is False
    assert stale_source_gate["checks"]["expected_source_file_id_matches"] is True
    assert stale_source_gate["checks"]["expected_trace_source_file_id_matches"] is False

    wrong_kernel = netgen_vol_first_order_fem_bem_trace_package_handoff(
        package,
        expected_bem_kernel_family="helmholtz_single_layer",
    )
    assert wrong_kernel["status"] == "needs_attention"
    assert wrong_kernel["checks"]["expected_bem_kernel_family_matches"] is False

    stale_assembly = dict(package, assembly_rule_id="remote_field_assembly_v0")
    stale_assembly["trace"] = dict(package["trace"], assembly_rule_id="remote_field_assembly_v0")
    stale_assembly["operators"] = {
        "trace": dict(package["operators"]["trace"], assembly_rule_id="remote_field_assembly_v0"),
        "bem": dict(package["operators"]["bem"], assembly_rule_id="remote_field_assembly_v0"),
    }
    stale_assembly_gate = netgen_vol_first_order_fem_bem_trace_package_handoff(
        stale_assembly,
        expected_assembly_rule_id="first_order_tet_h1_trace_tri_p1_bem_teaching_v1",
    )
    assert stale_assembly_gate["status"] == "needs_attention"
    assert stale_assembly_gate["checks"]["assembly_rule_id_consistent_when_present"] is True
    assert stale_assembly_gate["checks"]["expected_assembly_rule_id_matches"] is False

    missing_quadrature = dict(package)
    missing_quadrature.pop("quadrature_rule_id")
    missing_quadrature["trace"] = dict(package["trace"])
    missing_quadrature["trace"].pop("quadrature_rule_id")
    missing_quadrature["operators"] = {
        "trace": dict(package["operators"]["trace"]),
        "bem": dict(package["operators"]["bem"]),
    }
    missing_quadrature["operators"]["trace"].pop("quadrature_rule_id")
    missing_quadrature["operators"]["bem"].pop("quadrature_rule_id")
    missing_quadrature_gate = netgen_vol_first_order_fem_bem_trace_package_handoff(
        missing_quadrature,
        expected_quadrature_rule_id="tri_p1_exact_mass_regular_kernel_teaching_v1",
    )
    assert missing_quadrature_gate["status"] == "needs_attention"
    assert missing_quadrature_gate["checks"]["quadrature_rule_id_recorded_when_expected"] is False

    stale_result_artifact = dict(package)
    stale_result_artifact["result_artifact_id"] = "matlab_slot344_old_result_v0"
    stale_result_artifact_gate = netgen_vol_first_order_fem_bem_trace_package_handoff(
        stale_result_artifact,
        expected_result_artifact_id="matlab_slot344_fem_bem_manifest_result_v1",
        require_result_provenance=True,
    )
    assert stale_result_artifact_gate["status"] == "needs_attention"
    assert stale_result_artifact_gate["checks"]["result_artifact_id_recorded_when_required"] is True
    assert stale_result_artifact_gate["checks"]["expected_result_artifact_id_matches"] is False

    missing_timing_breakdown = dict(package)
    missing_timing_breakdown["timing_breakdown"] = {"mesh_read_s": 0.001}
    missing_timing_breakdown["execution"] = dict(package["execution"])
    missing_timing_breakdown["execution"]["timingBreakdown"] = {"mesh_read_s": 0.001}
    missing_timing_breakdown_gate = netgen_vol_first_order_fem_bem_trace_package_handoff(
        missing_timing_breakdown,
        require_result_provenance=True,
    )
    assert missing_timing_breakdown_gate["status"] == "needs_attention"
    assert missing_timing_breakdown_gate["checks"]["timing_breakdown_recorded_when_required"] is True
    assert missing_timing_breakdown_gate["checks"]["timing_breakdown_has_at_least_four_items"] is False

    stale_boundary_name = dict(package)
    stale_boundary_name["trace"] = dict(package["trace"], boundary_names=["outer", "coil", "outer", "outer"])
    stale_boundary_name_gate = netgen_vol_first_order_fem_bem_trace_package_handoff(
        stale_boundary_name,
        expected_boundary_names=["outer"],
    )
    assert stale_boundary_name_gate["status"] == "needs_attention"
    assert stale_boundary_name_gate["checks"]["boundary_number_name_pairs_consistent"] is False
    assert stale_boundary_name_gate["checks"]["expected_boundary_names_match"] is False

    stale_boundary_number = dict(package)
    stale_boundary_number["trace"] = dict(package["trace"], boundary_numbers=[1, 2, 1, 1])
    stale_boundary_number_gate = netgen_vol_first_order_fem_bem_trace_package_handoff(
        stale_boundary_number,
        expected_boundary_numbers=[1],
    )
    assert stale_boundary_number_gate["status"] == "needs_attention"
    assert stale_boundary_number_gate["checks"]["expected_boundary_numbers_match"] is False
    assert stale_boundary_number_gate["checks"]["boundary_row_identity_boundary_numbers_match"] is False

    stale_boundary_row_identity = dict(package)
    stale_boundary_row_identity["trace"] = dict(
        package["trace"],
        boundary_row_identity=[
            {"surface_triangle_index": 1, "surface_triangle_nodes": [1, 2, 3], "boundary_number": 1, "boundary_name": "outer"},
            {"surface_triangle_index": 2, "surface_triangle_nodes": [1, 4, 3], "boundary_number": 1, "boundary_name": "outer"},
            {"surface_triangle_index": 3, "surface_triangle_nodes": [1, 3, 4], "boundary_number": 1, "boundary_name": "outer"},
            {"surface_triangle_index": 4, "surface_triangle_nodes": [2, 4, 3], "boundary_number": 1, "boundary_name": "outer"},
        ],
    )
    stale_boundary_row_identity_gate = netgen_vol_first_order_fem_bem_trace_package_handoff(
        stale_boundary_row_identity,
        expected_boundary_row_identity=package["trace"]["boundary_row_identity"],
    )
    assert stale_boundary_row_identity_gate["status"] == "needs_attention"
    assert stale_boundary_row_identity_gate["checks"]["boundary_row_identity_matches_surface_triangles"] is False
    assert stale_boundary_row_identity_gate["checks"]["expected_boundary_row_identity_matches"] is False
    assert stale_boundary_row_identity_gate["boundary_row_identity_mismatch_rows"] == [2]

    wrong_trace = dict(package)
    wrong_trace["trace"] = dict(
        package["trace"],
        trace_matrix=[
            [1.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0],
        ],
    )
    wrong = netgen_vol_first_order_fem_bem_trace_package_handoff(wrong_trace)
    assert wrong["status"] == "needs_attention"
    assert wrong["checks"]["trace_matrix_matches_fem_node_ids"] is False

    stale_row_identity = dict(package)
    stale_row_identity["trace"] = dict(
        package["trace"],
        trace_row_identity=[
            {"trace_row_index": 1, "fem_node_id": 1, "bem_node_id": 1, "surface_node_index": 1},
            {"trace_row_index": 2, "fem_node_id": 2, "bem_node_id": 2, "surface_node_index": 2},
            {"trace_row_index": 3, "fem_node_id": 4, "bem_node_id": 3, "surface_node_index": 3},
            {"trace_row_index": 4, "fem_node_id": 4, "bem_node_id": 4, "surface_node_index": 4},
        ],
    )
    stale = netgen_vol_first_order_fem_bem_trace_package_handoff(stale_row_identity)
    assert stale["status"] == "needs_attention"
    assert stale["checks"]["trace_row_identity_fem_nodes_match"] is False
    assert stale["checks"]["trace_row_identity_matches_trace_matrix"] is False
    assert stale["trace_row_identity_mismatch_rows"] == [3]

    stale_operator_row_identity = dict(package)
    stale_operator_row_identity["operators"] = {
        "trace": {
            "trace_row_identity": [
                {"trace_row_index": 1, "fem_node_id": 1, "bem_node_id": 1, "surface_node_index": 1},
                {"trace_row_index": 2, "fem_node_id": 2, "bem_node_id": 2, "surface_node_index": 2},
                {"trace_row_index": 3, "fem_node_id": 4, "bem_node_id": 3, "surface_node_index": 3},
                {"trace_row_index": 4, "fem_node_id": 4, "bem_node_id": 4, "surface_node_index": 4},
            ],
        },
    }
    stale_operator = netgen_vol_first_order_fem_bem_trace_package_handoff(stale_operator_row_identity)
    assert stale_operator["status"] == "needs_attention"
    assert stale_operator["checks"]["operator_trace_row_identity_fem_nodes_match"] is False
    assert stale_operator["checks"]["operator_trace_row_identity_matches_trace_matrix"] is False
    assert stale_operator["checks"]["operator_trace_row_identity_matches_trace_identity"] is False
    assert stale_operator["operator_trace_row_identity_mismatch_rows"] == [3]

    curved = dict(package, polynomial_order=2, curved_element_count=1)
    curved_gate = netgen_vol_first_order_fem_bem_trace_package_handoff(curved)
    assert curved_gate["status"] == "needs_attention"
    assert curved_gate["checks"]["polynomial_order_first_order"] is False
    assert curved_gate["checks"]["curvedelements_absent"] is False


def test_netgen_vol_boundary_orientation_trace_package_gate_requires_explicit_normal_convention():
    package = {
        "mesh_id": "unit_tet_mesh",
        "export_id": "coreform_netgen_unit_tet",
        "gypsilab": {
            "elt": [[1, 2, 3], [1, 4, 2], [2, 4, 3], [3, 4, 1]],
        },
        "trace": {
            "surface_triangles": [[1, 2, 3], [1, 4, 2], [2, 4, 3], [3, 4, 1]],
            "boundary_orientation": "inward",
            "triangle_orientation_signs_to_outward": [-1, -1, -1, -1],
            "adjacent_tet_indices": [1, 1, 1, 1],
        },
    }

    gate = netgen_vol_boundary_orientation_trace_package_gate(
        package,
        expected_boundary_orientation="inward",
    )

    assert gate["status"] == "ok"
    assert gate["policy"] == "netgen_vol_boundary_orientation_trace_package_gate"
    assert gate["checks"]["boundary_orientation_recorded"] is True
    assert gate["checks"]["uniform_orientation_matches_signs"] is True
    assert gate["checks"]["expected_boundary_orientation_matches"] is True
    assert gate["triangle_orientation_signs_to_outward"] == [-1, -1, -1, -1]

    missing_orientation = {"trace": {"surface_triangles": package["trace"]["surface_triangles"]}}
    missing = netgen_vol_boundary_orientation_trace_package_gate(missing_orientation)
    assert missing["status"] == "needs_attention"
    assert missing["checks"]["boundary_orientation_recorded"] is False

    wrong_expected = netgen_vol_boundary_orientation_trace_package_gate(
        package,
        expected_boundary_orientation="outward",
    )
    assert wrong_expected["status"] == "needs_attention"
    assert wrong_expected["checks"]["expected_boundary_orientation_matches"] is False

    orphan = dict(package)
    orphan["trace"] = dict(package["trace"], adjacent_tet_indices=[1, 1, 0, 1])
    orphan_gate = netgen_vol_boundary_orientation_trace_package_gate(orphan)
    assert orphan_gate["status"] == "needs_attention"
    assert orphan_gate["checks"]["no_orphan_boundary_triangles"] is False


def test_netgen_vol_fem_bem_normal_flux_sign_package_gate_uses_outward_signs():
    package = {
        "trace": {
            "triangle_orientation_signs_to_outward": [-1, -1, -1, -1],
        },
        "normal_flux": {
            "normal_convention": "outward_from_volume",
            "stored_normal_flux": [1.5, 1.0, -3.0, 0.5],
            "orientation_corrected_normal_flux": [-1.5, -1.0, 3.0, -0.5],
            "outward_normal_flux_reference": [-1.5, -1.0, 3.0, -0.5],
            "closed_surface_flux_sum": 0.0,
            "expected_closed_surface_flux_sum": 0.0,
        },
    }

    gate = netgen_vol_fem_bem_normal_flux_sign_package_gate(
        package,
        expected_normal_convention="outward_from_volume",
    )

    assert gate["policy"] == "netgen_vol_fem_bem_normal_flux_sign_package_gate"
    assert gate["status"] == "ok"
    assert gate["sign_corrected_normal_flux"] == pytest.approx([-1.5, -1.0, 3.0, -0.5])
    assert gate["checks"]["sign_corrected_flux_matches_rows"] is True
    assert gate["checks"]["corrected_flux_matches_outward_reference"] is True
    assert gate["checks"]["closed_surface_flux_balance_ok"] is True
    assert gate["checks"]["normal_convention_recorded"] is True
    assert gate["checks"]["normal_convention_matches_expected"] is True

    stored_normals_reused = {
        "trace": package["trace"],
        "normal_flux": dict(
            package["normal_flux"],
            orientation_corrected_normal_flux=[1.5, 1.0, -3.0, 0.5],
        ),
    }
    bad = netgen_vol_fem_bem_normal_flux_sign_package_gate(stored_normals_reused)
    assert bad["status"] == "needs_attention"
    assert bad["checks"]["sign_corrected_flux_matches_rows"] is False
    assert bad["checks"]["corrected_flux_matches_outward_reference"] is False

    missing_signs = {"normal_flux": package["normal_flux"]}
    missing = netgen_vol_fem_bem_normal_flux_sign_package_gate(missing_signs)
    assert missing["status"] == "needs_attention"
    assert missing["checks"]["orientation_signs_recorded"] is False

    wrong_normal = {
        "trace": package["trace"],
        "normal_flux": dict(package["normal_flux"], normal_convention="stored_from_vol"),
    }
    wrong_normal_gate = netgen_vol_fem_bem_normal_flux_sign_package_gate(
        wrong_normal,
        expected_normal_convention="outward_from_volume",
    )
    assert wrong_normal_gate["status"] == "needs_attention"
    assert wrong_normal_gate["checks"]["normal_convention_recorded"] is True
    assert wrong_normal_gate["checks"]["normal_convention_matches_expected"] is False
    assert wrong_normal_gate["checks"]["corrected_flux_matches_outward_reference"] is True


def test_mqs_coulomb_gauge_efield_postprocess_gate_bundles_mqs_gauge_and_validity():
    artifacts = [
        {
            "kind": "mqs_solution",
            "case_id": "coil_mqs_001",
            "mesh_id": "unit_tet_mesh",
            "frequency_Hz": 1.0e6,
            "path": "slot160_mqs_solution.json",
            "gate_policy": "mqs_a_phi_solution",
            "status": "ok",
            "formulation": "A_phi",
        },
        {
            "kind": "coulomb_gauge",
            "case_id": "coil_mqs_001",
            "mesh_id": "unit_tet_mesh",
            "frequency_Hz": 1.0e6,
            "path": "slot160_coulomb_gauge.json",
            "gate_policy": "coulomb_gauge_postprocess",
            "status": "ok",
            "gauge_condition": "div_A_zero",
        },
        {
            "kind": "spatial_potential",
            "case_id": "coil_mqs_001",
            "mesh_id": "unit_tet_mesh",
            "frequency_Hz": 1.0e6,
            "path": "slot160_spatial_potential.json",
            "gate_policy": "electrostatic_potential_postprocess",
            "status": "ok",
            "boundary_condition_source": "conductor_surface_potential_from_coulomb_gauge",
        },
        {
            "kind": "electric_field",
            "case_id": "coil_mqs_001",
            "mesh_id": "unit_tet_mesh",
            "frequency_Hz": 1.0e6,
            "path": "slot160_efield.json",
            "gate_policy": "efield_gradient_recovery",
            "status": "ok",
            "E_unit": "V_per_m",
        },
        {
            "kind": "validity_envelope",
            "case_id": "coil_mqs_001",
            "mesh_id": "unit_tet_mesh",
            "frequency_Hz": 1.0e6,
            "path": "slot160_validity.json",
            "gate_policy": "mqs_darwin_fullwave_validity_envelope",
            "status": "ok",
            "frequency_ratio_to_fullwave_limit": 0.02,
            "dominant_inductive": True,
            "comparison_reference": "Darwin_full_wave",
        },
    ]

    gate = mqs_coulomb_gauge_efield_postprocess_gate(
        artifacts,
        expected_case_id="coil_mqs_001",
        expected_mesh_id="unit_tet_mesh",
        expected_frequency_Hz=1.0e6,
    )

    assert gate["status"] == "ok"
    assert gate["policy"] == "mqs_coulomb_gauge_efield_postprocess_gate"
    assert gate["checks"]["mqs_a_phi_formulation_recorded"] is True
    assert gate["checks"]["validity_envelope_ok"] is True

    stale_frequency = [dict(row) for row in artifacts]
    stale_frequency[3]["frequency_Hz"] = 1.2e6
    stale_frequency_gate = mqs_coulomb_gauge_efield_postprocess_gate(stale_frequency)
    assert stale_frequency_gate["status"] == "needs_attention"
    assert stale_frequency_gate["checks"]["frequencies_match"] is False

    missing_bc = [dict(row) for row in artifacts]
    missing_bc[2].pop("boundary_condition_source")
    missing_bc_gate = mqs_coulomb_gauge_efield_postprocess_gate(missing_bc)
    assert missing_bc_gate["status"] == "needs_attention"
    assert missing_bc_gate["checks"]["spatial_potential_bc_recorded"] is False

    invalid_validity = [dict(row) for row in artifacts]
    invalid_validity[4]["frequency_ratio_to_fullwave_limit"] = 0.5
    invalid_validity_gate = mqs_coulomb_gauge_efield_postprocess_gate(invalid_validity)
    assert invalid_validity_gate["status"] == "needs_attention"
    assert invalid_validity_gate["checks"]["validity_envelope_ok"] is False


def test_acoustic_plane_wave_intensity_gate_closes_comsol_amplitude_contract():
    gate = acoustic_plane_wave_intensity_convention_gate(
        pressure_peak=2.0,
        rho=1.2,
        c=343.0,
        area=0.5,
    )

    z0 = 1.2 * 343.0
    expected_intensity = 0.5 * 2.0 * (2.0 / z0)
    assert gate["status"] == "ok"
    assert gate["specific_impedance_Pa_s_per_m"] == pytest.approx(z0)
    assert gate["pressure_rms_Pa"] == pytest.approx(math.sqrt(2.0))
    assert gate["intensity_from_peak_W_per_m2"] == pytest.approx(expected_intensity)
    assert gate["intensity_from_rms_W_per_m2"] == pytest.approx(expected_intensity)
    assert gate["power_from_peak_W"] == pytest.approx(0.5 * expected_intensity)
    assert gate["peak_rms_power_residual_W"] == pytest.approx(0.0)

    with pytest.raises(ValueError):
        acoustic_plane_wave_intensity_convention_gate(1.0, rho=0.0)


def test_acoustic_normal_incidence_interface_gate_closes_comsol_fem_bem_contract():
    gate = acoustic_normal_incidence_interface_gate(
        rho_left=1.2,
        c_left=343.0,
        rho_right=1000.0,
        c_right=1480.0,
        pressure_incident_peak=2.0,
        area=0.25,
    )

    z1 = 1.2 * 343.0
    z2 = 1000.0 * 1480.0
    expected_r = (z2 - z1) / (z2 + z1)
    expected_t = 2.0 * z2 / (z2 + z1)
    assert gate["status"] == "ok"
    assert gate["pressure_reflection_coefficient"] == pytest.approx(expected_r)
    assert gate["pressure_transmission_coefficient"] == pytest.approx(expected_t)
    assert gate["pressure_jump_Pa"] == pytest.approx(0.0)
    assert gate["velocity_jump_m_per_s"] == pytest.approx(0.0)
    assert gate["reflected_power_ratio"] + gate["transmitted_power_ratio"] == pytest.approx(1.0)
    assert gate["transmitted_power_ratio"] == pytest.approx(4.0 * z1 * z2 / (z1 + z2) ** 2)

    matched = acoustic_normal_incidence_interface_gate(1.2, 343.0, 1.2, 343.0)
    assert matched["pressure_reflection_coefficient"] == pytest.approx(0.0)
    assert matched["transmitted_power_ratio"] == pytest.approx(1.0)


def test_acoustic_interface_result_package_gate_keeps_rows_solver_ready():
    interface = acoustic_normal_incidence_interface_gate(
        rho_left=1.2,
        c_left=343.0,
        rho_right=1000.0,
        c_right=1480.0,
        pressure_incident_peak=2.0,
        area=0.25,
    )
    shared = {
        "case_id": "air_water_normal_incidence",
        "run_id": "slot177",
        "export_id": "interface_package_v1",
        "frequency_hz": 1500.0,
        "source_tool": "analytic_reference",
        "status": "ok",
        "pass": True,
    }
    rows = [
        {
            **shared,
            "kind": "material_impedance",
            "path": "memory://slot177/material_impedance.json",
            "gate_policy": "acoustic_material_impedance_contract",
            "rho_left_kg_per_m3": interface["rho_left_kg_per_m3"],
            "c_left_m_per_s": interface["c_left_m_per_s"],
            "rho_right_kg_per_m3": interface["rho_right_kg_per_m3"],
            "c_right_m_per_s": interface["c_right_m_per_s"],
            "z_left_Pa_s_per_m": interface["z_left_Pa_s_per_m"],
            "z_right_Pa_s_per_m": interface["z_right_Pa_s_per_m"],
        },
        {
            **shared,
            "kind": "interface_continuity",
            "path": "memory://slot177/interface_continuity.json",
            "gate_policy": "acoustic_interface_continuity_contract",
            "pressure_jump_Pa": interface["pressure_jump_Pa"],
            "normal_velocity_continuity_residual_m_per_s": interface["velocity_jump_m_per_s"],
        },
        {
            **shared,
            "kind": "power_split",
            "path": "memory://slot177/power_split.json",
            "gate_policy": "acoustic_power_split_contract",
            "reflected_power_ratio": interface["reflected_power_ratio"],
            "transmitted_power_ratio": interface["transmitted_power_ratio"],
            "power_balance_error": interface["power_balance_error"],
        },
    ]

    gate = acoustic_interface_result_package_gate(
        rows,
        expected_case_id="air_water_normal_incidence",
        expected_run_id="slot177",
        expected_export_id="interface_package_v1",
        expected_frequency_hz=1500.0,
        residual_tol=1.0e-12,
    )

    assert gate["status"] == "ok"
    assert gate["present_kinds"] == {
        "interface_continuity": 1,
        "material_impedance": 1,
        "power_split": 1,
    }
    assert gate["max_pressure_residual_Pa"] == pytest.approx(0.0)
    assert gate["max_normal_velocity_residual_m_per_s"] == pytest.approx(0.0)
    assert gate["max_abs_power_balance_error"] == pytest.approx(0.0)
    assert gate["checks"]["passive_power_split_ok"] is True

    bad_rows = [dict(row) for row in rows]
    bad_rows[2]["transmitted_power_ratio"] += 0.05
    bad_rows[2].pop("power_balance_error")
    bad = acoustic_interface_result_package_gate(bad_rows, residual_tol=1.0e-12)
    assert bad["status"] == "needs_attention"
    assert bad["checks"]["passive_power_split_ok"] is False


def test_acoustic_impedance_power_result_package_gate_keeps_one_port_rows_solver_ready():
    z_load = 2.0 + 0.5j
    z0 = 1.0
    gamma = (z_load - z0) / (z_load + z0)
    incident = 0.5
    reflected = incident * abs(gamma) ** 2
    active = incident - reflected
    row = {
        "kind": "impedance_power",
        "case_id": "duct_one_port",
        "run_id": "run_acoustic_001",
        "export_id": "power_export_A",
        "frequency_hz": 1000.0,
        "source_tool": "open_acoustic_solver",
        "path": "memory://duct_one_port/power_export_A",
        "gate_policy": "acoustic_impedance_power_contract",
        "status": "ok",
        "specific_impedance": {"real": z_load.real, "imag": z_load.imag},
        "characteristic_normal_impedance": z0,
        "pressure_reflection_coefficient": {"real": gamma.real, "imag": gamma.imag},
        "absorption_coefficient": 1.0 - abs(gamma) ** 2,
        "incident_intensity": incident,
        "reflected_intensity": reflected,
        "boundary_active_intensity_into_load": active,
        "power_balance_residual": 0.0,
    }

    gate = acoustic_impedance_power_result_package_gate(
        [row],
        expected_case_id="duct_one_port",
        expected_run_id="run_acoustic_001",
        expected_export_id="power_export_A",
        expected_frequency_hz=1000.0,
        residual_tol=1.0e-12,
    )

    assert gate["policy"] == "acoustic_impedance_power_result_package_gate"
    assert gate["status"] == "ok"
    assert gate["checks"]["impedance_passive_or_reactive"] is True
    assert gate["checks"]["absorption_matches_reflection"] is True
    assert gate["checks"]["power_balance_ok"] is True
    assert gate["min_boundary_active_power"] == pytest.approx(active)

    active_impedance = dict(row)
    active_impedance["specific_impedance"] = {"real": -2.0, "imag": 0.0}
    active_impedance["pressure_reflection_coefficient"] = {"real": 3.0, "imag": 0.0}
    active_impedance["absorption_coefficient"] = -8.0
    active_impedance["boundary_active_intensity_into_load"] = -4.0
    active_impedance["reflected_intensity"] = 4.5
    active_gate = acoustic_impedance_power_result_package_gate([active_impedance])
    assert active_gate["status"] == "needs_attention"
    assert active_gate["checks"]["impedance_passive_or_reactive"] is False
    assert active_gate["checks"]["boundary_active_power_nonnegative"] is False

    stale = dict(row)
    stale["export_id"] = "power_export_old"
    mixed = acoustic_impedance_power_result_package_gate([row, stale])
    assert mixed["status"] == "needs_attention"
    assert mixed["checks"]["single_export_id"] is False


def test_thermal_annulus_conductance_gate_closes_comsol_heat_transfer_contract():
    gate = thermal_annulus_conductance_gate(
        inner_radius=0.01,
        outer_radius=0.05,
        conductivity=50.0,
        delta_temperature=10.0,
        length=1.0,
        measured_temperature=5.0,
        measured_conductance=2.0 * math.pi * 50.0 / math.log(5.0),
        rtol=1.0e-12,
    )

    assert gate["status"] == "ok"
    assert gate["policy"] == "thermal_annulus_radial_conductance_gate"
    assert gate["probe_radius_m"] == pytest.approx(math.sqrt(0.01 * 0.05))
    assert gate["temperature_at_probe_K"] == pytest.approx(5.0)
    assert gate["geometric_mean_temperature_K"] == pytest.approx(5.0)
    assert gate["thermal_conductance_W_per_K"] == pytest.approx(2.0 * math.pi * 50.0 / math.log(5.0))
    assert gate["heat_rate_W"] == pytest.approx(10.0 * 2.0 * math.pi * 50.0 / math.log(5.0))
    assert gate["temperature_rel_error"] == pytest.approx(0.0)
    assert gate["conductance_rel_error"] == pytest.approx(0.0)

    bad = thermal_annulus_conductance_gate(
        0.01,
        0.05,
        50.0,
        delta_temperature=10.0,
        measured_temperature=4.0,
        measured_conductance=190.0,
        rtol=1.0e-12,
    )
    assert bad["status"] == "needs_attention"
    assert bad["checks"]["measured_temperature_ok"] is False
    assert bad["checks"]["measured_conductance_ok"] is False


def test_thermal_layer_stack_conductance_gate_closes_comsol_heat_stack_contract():
    gate = thermal_layer_stack_conductance_gate(
        area=0.02,
        thicknesses=[0.0008, 0.0012],
        conductivities=[2.0, 6.0],
        delta_temperature=10.0,
        measured_conductance=0.02 / (0.0008 / 2.0 + 0.0012 / 6.0),
        measured_interface_temperatures=[10.0 - 10.0 * (0.0008 / 2.0) / (0.0008 / 2.0 + 0.0012 / 6.0)],
        measured_heat_rate=0.02 * 10.0 / (0.0008 / 2.0 + 0.0012 / 6.0),
    )

    assert gate["status"] == "ok"
    assert gate["policy"] == "thermal_layer_stack_series_resistance_gate"
    assert gate["total_resistance_K_per_W"] == pytest.approx(0.03)
    assert gate["thermal_conductance_W_per_K"] == pytest.approx(33.333333333333336)
    assert gate["heat_rate_W"] == pytest.approx(333.33333333333337)
    assert gate["heat_flux_W_per_m2"] == pytest.approx(16666.666666666668)
    assert gate["layer_temperature_drops_K"] == pytest.approx([6.666666666666667, 3.3333333333333335])
    assert gate["interface_temperatures_K"] == pytest.approx([3.333333333333333])
    assert gate["effective_conductivity_W_per_m_K"] == pytest.approx(3.3333333333333335)

    bad = thermal_layer_stack_conductance_gate(
        0.02,
        [0.0008, 0.0012],
        [2.0, 6.0],
        10.0,
        measured_interface_temperatures=[6.666666666666667],
    )
    assert bad["status"] == "needs_attention"
    assert bad["checks"]["measured_interface_temperatures_ok"] is False


def test_thermal_conduction_convection_robin_gate_closes_comsol_robin_contract():
    gate = thermal_conduction_convection_robin_gate(
        area=1.0,
        thicknesses=[0.20],
        conductivities=[10.0],
        heat_transfer_coefficient=15.0,
        hot_temperature=100.0,
        ambient_temperature=0.0,
        measured_conductance=11.538461541399737,
        measured_heat_flux=1153.8461541399738,
        measured_surface_temperature=76.923076917294878,
        measured_heat_transfer_coefficient=15.00000000046324,
        rtol=2.0e-9,
        atol=1.0e-8,
    )

    assert gate["status"] == "ok"
    assert gate["policy"] == "thermal_conduction_convection_robin_gate"
    assert gate["conduction_resistances_K_per_W"] == pytest.approx([0.02])
    assert gate["convection_resistance_K_per_W"] == pytest.approx(1.0 / 15.0)
    assert gate["total_resistance_K_per_W"] == pytest.approx(0.08666666666666667)
    assert gate["thermal_conductance_W_per_K"] == pytest.approx(11.538461538461538)
    assert gate["heat_flux_W_per_m2"] == pytest.approx(1153.8461538461538)
    assert gate["robin_surface_temperature_K"] == pytest.approx(76.92307692307692)
    assert gate["measured_heat_transfer_coefficient_rel_error"] < 1.0e-10

    bad = thermal_conduction_convection_robin_gate(
        1.0,
        [0.20],
        [10.0],
        15.0,
        100.0,
        measured_surface_temperature=90.0,
        measured_heat_transfer_coefficient=11.0,
        rtol=1.0e-9,
    )
    assert bad["status"] == "needs_attention"
    assert bad["checks"]["measured_robin_surface_temperature_ok"] is False
    assert bad["checks"]["measured_heat_transfer_coefficient_ok"] is False


def test_spherical_dirichlet_laplacian_eigen_gate_closes_cubit_hex_sphere_contract():
    gate = spherical_dirichlet_laplacian_eigen_gate(
        radius=1.0,
        eigenvalues=[9.86268248871533, 20.178703158644527],
        rtol_first=2.0e-3,
        rtol_l1=2.0e-3,
    )

    assert gate["status"] == "ok"
    assert gate["policy"] == "spherical_dirichlet_laplacian_eigen_gate"
    assert gate["reference"]["first_radial_lambda"] == pytest.approx(math.pi**2)
    assert gate["reference"]["first_l1_lambda"] == pytest.approx(4.493409457909064**2)
    assert gate["first_rel_error"] < 1.0e-3
    assert gate["l1_rel_error"] < 1.0e-3

    bad = spherical_dirichlet_laplacian_eigen_gate(
        radius=1.0,
        eigenvalues=[9.0, 22.0],
        rtol_first=2.0e-3,
        rtol_l1=2.0e-3,
    )
    assert bad["status"] == "needs_attention"
    assert bad["checks"]["first_mode_matches_ball"] is False
    assert bad["checks"]["l1_mode_matches_ball"] is False


def test_ipm_saliency_torque_component_gate_closes_femm_ld_lq_contract():
    gate = ipm_saliency_torque_component_gate(
        lambda_m=0.08,
        Ld=0.004,
        Lq=0.0085,
        id_current=-6.0,
        iq_current=14.0,
        pole_pairs=4,
    )

    assert gate["status"] == "ok"
    assert gate["policy"] == "ipm_saliency_torque_component_gate"
    assert gate["saliency_ratio_Lq_over_Ld"] == pytest.approx(2.125)
    assert gate["magnet_torque_Nm"] == pytest.approx(6.72)
    assert gate["reluctance_torque_Nm"] == pytest.approx(2.268)
    assert gate["total_torque_Nm"] == pytest.approx(8.988)
    assert gate["reluctance_fraction_of_total"] == pytest.approx(2.268 / 8.988)

    wrong_id = ipm_saliency_torque_component_gate(
        lambda_m=0.08,
        Ld=0.004,
        Lq=0.0085,
        id_current=6.0,
        iq_current=14.0,
        pole_pairs=4,
    )
    assert wrong_id["status"] == "needs_attention"
    assert wrong_id["checks"]["field_weakening_negative_id"] is False
    assert wrong_id["checks"]["reluctance_torque_adds_to_magnet"] is False


def test_ipm_saliency_component_row_can_seed_jmag_torque_table_contract():
    lambda_m = 0.08
    Ld = 0.004
    Lq = 0.0085
    id_current = -6.0
    iq_current = 14.0
    pole_pairs = 4
    current = math.hypot(id_current, iq_current)
    gamma_deg = math.degrees(math.asin(-id_current / current))

    component = ipm_saliency_torque_component_gate(
        lambda_m=lambda_m,
        Ld=Ld,
        Lq=Lq,
        id_current=id_current,
        iq_current=iq_current,
        pole_pairs=pole_pairs,
    )
    rows = []
    for gamma in (-30.0, -15.0, 0.0, gamma_deg, 30.0, 45.0):
        row_id, row_iq = dq_current_from_gamma_deg(current, gamma)
        rows.append({
            "gamma_deg": gamma,
            "id_A": row_id,
            "iq_A": row_iq,
            "torque_Nm": lumped_pm_dq_torque(lambda_m, Ld, Lq, row_id, row_iq, pole_pairs),
        })
    table = dq_torque_table_health(rows, lambda_m, Ld, Lq, current, pole_pairs, tol=1.0e-11)
    operating = next(row for row in table["rows"] if abs(row["gamma_deg"] - gamma_deg) < 1.0e-12)

    assert component["status"] == "ok"
    assert gamma_deg == pytest.approx(23.198590513648185)
    assert operating["id_A"] == pytest.approx(-6.0)
    assert operating["iq_A"] == pytest.approx(14.0)
    assert operating["torque_Nm"] == pytest.approx(component["total_torque_Nm"])
    assert table["status"] == "ok"


def test_coaxial_pm_force_gap_sweep_gate_closes_elf_force_observable_contract():
    constant = 1.6e-12
    rows = [
        {"gap_m": 0.002, "force_N": -constant / 0.002**4},
        {"gap_m": 0.004, "force_N": -constant / 0.004**4},
        {"gap_m": 0.008, "force_N": -constant / 0.008**4},
    ]
    gate = coaxial_pm_force_gap_sweep_gate(rows, rtol_invariant=1.0e-12)

    assert gate["status"] == "ok"
    assert gate["force_ratio_first_last"] == pytest.approx(256.0)
    assert gate["expected_force_ratio_first_last"] == pytest.approx(256.0)
    assert gate["max_force_gap4_invariant_rel_error"] < 1.0e-12

    wrong_sign = [dict(row) for row in rows]
    wrong_sign[1]["force_N"] = abs(wrong_sign[1]["force_N"])
    bad = coaxial_pm_force_gap_sweep_gate(wrong_sign, expected_sign="attractive_negative")
    assert bad["status"] == "needs_attention"
    assert bad["checks"]["force_sign_ok"] is False


def test_maxwell_stress_surface_package_gate_keeps_surface_and_method_together():
    artifacts = [
        {
            "kind": "stress_surface",
            "case_id": "elf_force_case_182",
            "stress_surface_id": "MCM_rotor_001",
            "result_set_id": "elf_force_result_198",
            "run_artifact_id": "elf_run_manifest_270_A",
            "result_revision_id": "elf_result_revision_270_A",
            "observable_id": "fort_z_force",
            "closed_surface": True,
            "normal_orientation": "outward",
            "formulation_id": "elf_magic_moment_mcm_fort_v1",
            "kernel_family": "laplace_magnetostatic_mom",
            "singular_treatment": "product_default_surface_quadrature",
            "symmetry_factor": 4,
            "status": "ok",
        },
        {
            "kind": "solve_command",
            "case_id": "elf_force_case_182",
            "stress_surface_id": "MCM_rotor_001",
            "result_set_id": "elf_force_result_198",
            "run_artifact_id": "elf_run_manifest_270_A",
            "result_revision_id": "elf_result_revision_270_A",
            "observable_id": "fort_z_force",
            "sol_command": "SOL FORT",
            "axis": "z",
            "formulation_id": "elf_magic_moment_mcm_fort_v1",
            "kernel_family": "laplace_magnetostatic_mom",
            "singular_treatment": "product_default_surface_quadrature",
            "symmetry_factor": 4,
            "status": "ok",
        },
        {
            "kind": "observable_table",
            "case_id": "elf_force_case_182",
            "stress_surface_id": "MCM_rotor_001",
            "result_set_id": "elf_force_result_198",
            "run_artifact_id": "elf_run_manifest_270_A",
            "result_revision_id": "elf_result_revision_270_A",
            "observable_id": "fort_z_force",
            "force_method": "FORT",
            "axis": "z",
            "sign_convention": "positive_z_force",
            "quantity_dimension": "3d_total",
            "force_unit": "N",
            "formulation_id": "elf_magic_moment_mcm_fort_v1",
            "kernel_family": "laplace_magnetostatic_mom",
            "singular_treatment": "product_default_surface_quadrature",
            "symmetry_factor": 4,
            "status": "ok",
        },
    ]
    gate = maxwell_stress_surface_package_gate(
        artifacts,
        expected_case_id="elf_force_case_182",
        expected_surface_id="MCM_rotor_001",
        expected_result_set_id="elf_force_result_198",
        expected_observable_id="fort_z_force",
        required_axis="z",
        expected_normal_orientation="outward",
        expected_formulation_id="elf_magic_moment_mcm_fort_v1",
        expected_kernel_family="laplace_magnetostatic_mom",
        expected_run_artifact_id="elf_run_manifest_270_A",
        expected_result_revision_id="elf_result_revision_270_A",
    )

    assert gate["policy"] == "maxwell_stress_surface_package_gate"
    assert gate["status"] == "ok"
    assert gate["methods"] == ["fort", "sol_fort"]
    assert gate["axes"] == ["z"]
    assert gate["result_set_ids"] == ["elf_force_result_198"]
    assert gate["run_artifact_ids"] == ["elf_run_manifest_270_A"]
    assert gate["result_revision_ids"] == ["elf_result_revision_270_A"]
    assert gate["observable_ids"] == ["fort_z_force"]
    assert gate["quantity_dimensions"] == ["3d_total"]
    assert gate["force_units"] == ["n"]
    assert gate["expected_normal_orientation"] == "outward"
    assert gate["formulation_ids"] == ["elf_magic_moment_mcm_fort_v1"]
    assert gate["kernel_families"] == ["laplace_magnetostatic_mom"]
    assert gate["singular_treatments"] == ["product_default_surface_quadrature"]
    assert gate["checks"]["result_set_id_matches_expected"] is True
    assert gate["checks"]["run_artifact_id_matches_expected"] is True
    assert gate["checks"]["result_revision_id_matches_expected"] is True
    assert gate["checks"]["observable_id_matches_expected"] is True
    assert gate["checks"]["observable_id_axis_matches_axis_when_named"] is True
    assert gate["checks"]["normal_orientation_matches_expected"] is True
    assert gate["checks"]["formulation_id_matches_expected"] is True
    assert gate["checks"]["kernel_family_matches_expected"] is True
    assert gate["checks"]["singular_treatment_consistent_when_present"] is True
    assert gate["checks"]["force_unit_matches_quantity_dimension_when_present"] is True
    assert all(gate["checks"].values())

    open_surface = [dict(row) for row in artifacts]
    open_surface[0]["closed_surface"] = False
    open_gate = maxwell_stress_surface_package_gate(open_surface)
    assert open_gate["status"] == "needs_attention"
    assert open_gate["checks"]["closed_surface_confirmed"] is False

    stale_surface = [dict(row) for row in artifacts]
    stale_surface[2]["stress_surface_id"] = "MCM_stale"
    stale_gate = maxwell_stress_surface_package_gate(stale_surface)
    assert stale_gate["status"] == "needs_attention"
    assert stale_gate["checks"]["stress_surface_id_consistent"] is False

    wrong_normal = [dict(row) for row in artifacts]
    wrong_normal[0]["normal_orientation"] = "inward"
    normal_gate = maxwell_stress_surface_package_gate(
        wrong_normal,
        expected_normal_orientation="outward",
    )
    assert normal_gate["status"] == "needs_attention"
    assert normal_gate["checks"]["normal_orientation_recorded"] is True
    assert normal_gate["checks"]["normal_orientation_matches_expected"] is False

    stale_result = [dict(row) for row in artifacts]
    stale_result[2]["result_set_id"] = "elf_force_old"
    stale_result_gate = maxwell_stress_surface_package_gate(
        stale_result,
        expected_result_set_id="elf_force_result_198",
        expected_observable_id="fort_z_force",
    )
    assert stale_result_gate["status"] == "needs_attention"
    assert stale_result_gate["checks"]["result_set_id_consistent_when_present"] is False
    assert stale_result_gate["checks"]["result_set_id_matches_expected"] is False

    stale_run = [dict(row) for row in artifacts]
    stale_run[1]["run_artifact_id"] = "elf_run_manifest_old"
    stale_run_gate = maxwell_stress_surface_package_gate(
        stale_run,
        expected_run_artifact_id="elf_run_manifest_270_A",
        expected_result_revision_id="elf_result_revision_270_A",
    )
    assert stale_run_gate["status"] == "needs_attention"
    assert stale_run_gate["checks"]["run_artifact_id_consistent_when_present"] is False
    assert stale_run_gate["checks"]["run_artifact_id_matches_expected"] is False
    assert stale_run_gate["checks"]["result_revision_id_matches_expected"] is True

    wrong_observable_axis = [dict(row) for row in artifacts]
    wrong_observable_axis[1]["axis"] = "y"
    wrong_observable_axis[2]["axis"] = "y"
    observable_axis_gate = maxwell_stress_surface_package_gate(wrong_observable_axis)
    assert observable_axis_gate["status"] == "needs_attention"
    assert observable_axis_gate["checks"]["axis_consistent"] is True
    assert observable_axis_gate["checks"]["observable_id_consistent_when_present"] is True
    assert observable_axis_gate["checks"]["observable_id_axis_matches_axis_when_named"] is False
    assert observable_axis_gate["observable_axis_mismatches"] == [
        {"observable_id": "fort_z_force", "axis": "y", "observable_axis_tokens": ["z"]}
    ]

    wrong_kernel = [dict(row) for row in artifacts]
    wrong_kernel[2]["kernel_family"] = "helmholtz_time_harmonic_mom"
    kernel_gate = maxwell_stress_surface_package_gate(
        wrong_kernel,
        expected_formulation_id="elf_magic_moment_mcm_fort_v1",
        expected_kernel_family="laplace_magnetostatic_mom",
    )
    assert kernel_gate["status"] == "needs_attention"
    assert kernel_gate["checks"]["formulation_id_matches_expected"] is True
    assert kernel_gate["checks"]["kernel_family_consistent_when_present"] is False
    assert kernel_gate["checks"]["kernel_family_matches_expected"] is False

    wrong_method = [dict(row) for row in artifacts]
    wrong_method[1]["sol_command"] = "SOL FORC"
    method_gate = maxwell_stress_surface_package_gate(wrong_method)
    assert method_gate["status"] == "needs_attention"
    assert method_gate["checks"]["maxwell_stress_method_declared"] is False

    unit_mismatch = [dict(row) for row in artifacts]
    unit_mismatch[2]["quantity_dimension"] = "2d_per_length"
    unit_gate = maxwell_stress_surface_package_gate(unit_mismatch)
    assert unit_gate["status"] == "needs_attention"
    assert unit_gate["checks"]["force_unit_matches_quantity_dimension_when_present"] is False
    assert unit_gate["bad_force_unit_dimension_pairs"][0]["expected_force_units"] == ["n/m", "n_per_m"]


def test_two_port_abcd_cascade_gate_closes_cst_network_export_contract():
    z0 = 50.0

    def quarter_wave(zline):
        return (0.0, 1j * zline, 1j / zline, 0.0)

    gate = two_port_abcd_cascade_gate(
        [quarter_wave(75.0), quarter_wave(75.0)],
        z0=z0,
        expect_lossless=True,
    )

    assert gate["status"] == "ok"
    assert gate["n_sections"] == 2
    assert gate["A"]["real"] == pytest.approx(-1.0)
    assert gate["D"]["real"] == pytest.approx(-1.0)
    assert gate["B"]["abs"] < 1.0e-12
    assert gate["C"]["abs"] < 1.0e-12
    assert gate["s11"]["abs"] < 1.0e-12
    assert gate["s21"]["abs"] == pytest.approx(1.0)
    assert gate["lossless_power_sum"] == pytest.approx(1.0)

    active = two_port_abcd_cascade_gate(
        [{"A": 1.0, "B": -25.0, "C": 0.0, "D": 1.0}],
        z0=z0,
    )
    assert active["status"] == "needs_attention"
    assert active["checks"]["sparameter_passivity_ok"] is False


def test_cst_abcd_cascade_solver_ready_manifest_gate_keeps_export_identity():
    z0 = 50.0

    def quarter_wave(zline):
        return (0.0, 1j * zline, 1j / zline, 0.0)

    cascade = two_port_abcd_cascade_gate(
        [quarter_wave(75.0), quarter_wave(75.0)],
        z0=z0,
        expect_lossless=True,
    )
    artifacts = [
        {
            "kind": "port_metadata",
            "project_id": "rf_chain_v1",
            "run_id": "run_abcd_001",
            "export_id": "abcd_export_C",
            "source_tool": "CST",
            "path": "slot175_ports.s2p",
            "gate_policy": "touchstone_port_metadata_gate",
            "status": "ok",
            "pass": True,
            "design_frequency_hz": 1.0e9,
            "data_format": "MA",
            "reference_impedance_ohm": z0,
            "port_order": ["P1", "P2"],
        },
        {
            "kind": "abcd_cascade",
            "project_id": "rf_chain_v1",
            "run_id": "run_abcd_001",
            "export_id": "abcd_export_C",
            "source_tool": "CST",
            "path": "slot175_abcd_chain.json",
            "gate_policy": cascade["policy"],
            "status": cascade["status"],
            "pass": cascade["status"] == "ok",
            "design_frequency_hz": 1.0e9,
            "n_sections": cascade["n_sections"],
            "z0_ohm": cascade["z0_ohm"],
            "expect_lossless": cascade["expect_lossless"],
            "checks": cascade["checks"],
        },
    ]

    gate = cst_abcd_cascade_solver_ready_manifest_gate(
        artifacts,
        expected_project_id="rf_chain_v1",
        expected_run_id="run_abcd_001",
        expected_export_id="abcd_export_C",
        expected_design_frequency_hz=1.0e9,
    )

    assert gate["status"] == "ok"
    assert gate["policy"] == "cst_abcd_cascade_solver_ready_manifest_gate"
    assert gate["checks"]["required_kinds_present"] is True
    assert gate["checks"]["z0_consistent"] is True
    assert gate["checks"]["abcd_sparameters_passive"] is True
    assert gate["checks"]["abcd_lossless_when_requested"] is True

    stale_export = [dict(row) for row in artifacts]
    stale_export[1]["export_id"] = "old_export"
    stale_gate = cst_abcd_cascade_solver_ready_manifest_gate(stale_export)
    assert stale_gate["status"] == "needs_attention"
    assert stale_gate["checks"]["export_ids_unique"] is False

    wrong_source = [dict(row) for row in artifacts]
    wrong_source[0]["source_tool"] = "HFSS"
    wrong_source_gate = cst_abcd_cascade_solver_ready_manifest_gate(wrong_source)
    assert wrong_source_gate["status"] == "needs_attention"
    assert wrong_source_gate["checks"]["source_tool_is_cst"] is False

    wrong_z0 = [dict(row) for row in artifacts]
    wrong_z0[0]["reference_impedance_ohm"] = 75.0
    wrong_z0_gate = cst_abcd_cascade_solver_ready_manifest_gate(wrong_z0)
    assert wrong_z0_gate["status"] == "needs_attention"
    assert wrong_z0_gate["checks"]["z0_consistent"] is False

    active = two_port_abcd_cascade_gate([{"A": 1.0, "B": -25.0, "C": 0.0, "D": 1.0}], z0=z0)
    active_artifacts = [dict(artifacts[0]), dict(artifacts[1])]
    active_artifacts[1].update({
        "status": active["status"],
        "pass": False,
        "n_sections": active["n_sections"],
        "z0_ohm": active["z0_ohm"],
        "expect_lossless": active["expect_lossless"],
        "checks": active["checks"],
    })
    active_gate = cst_abcd_cascade_solver_ready_manifest_gate(active_artifacts)
    assert active_gate["status"] == "needs_attention"
    assert active_gate["checks"]["abcd_sparameters_passive"] is False


def test_mesh_import_quality_manifest_gate_keeps_tri_tet_policy_visible():
    gate = mesh_import_quality_manifest_gate(
        surface_element_types=("triangle", "tri3"),
        volume_element_types=("tetrahedron", "tet4"),
        order=1,
        min_scaled_jacobian_before=0.36,
        min_scaled_jacobian_after=0.75,
        min_scaled_jacobian_threshold=0.2,
        negative_jacobian_count_before=2,
        negative_jacobian_count_after=0,
        cad_connectivity_recorded=True,
        cad_compliance_recorded=True,
        boundary_conformity_tolerance=1.0e-6,
        max_boundary_distance=4.0e-7,
    )

    assert gate["status"] == "ok"
    assert gate["policy"] == "mesh_import_quality_manifest_gate"
    assert gate["surface_element_types"] == ["tri", "tri"]
    assert gate["volume_element_types"] == ["tet", "tet"]
    assert gate["min_scaled_jacobian_improvement"] == pytest.approx(0.39)
    assert gate["checks"]["first_order_tri_tet_policy_honored"] is True
    assert gate["checks"]["boundary_conformity_within_tolerance"] is True

    high_order = mesh_import_quality_manifest_gate(order=4)
    assert high_order["status"] == "needs_attention"
    assert high_order["checks"]["first_order_only"] is False
    assert "first_order_only" in high_order["issues"]

    wrong_topology = mesh_import_quality_manifest_gate(
        surface_element_types=("tri", "quad"),
        volume_element_types=("tet", "hex"),
        min_scaled_jacobian_after=0.05,
        min_scaled_jacobian_threshold=0.1,
        negative_jacobian_count_after=1,
        cad_connectivity_recorded=False,
        cad_compliance_recorded=False,
        max_boundary_distance=1.0e-5,
    )
    assert wrong_topology["status"] == "needs_attention"
    assert wrong_topology["checks"]["surface_triangles_only"] is False
    assert wrong_topology["checks"]["volume_tetrahedra_only"] is False
    assert wrong_topology["checks"]["no_negative_jacobian_after"] is False
    assert wrong_topology["checks"]["boundary_conformity_within_tolerance"] is False
