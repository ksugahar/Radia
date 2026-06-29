import math

import pytest

from radia_mcp.radia_ngsolve.slot_gates import (
    acoustic_plane_wave_intensity_convention_gate,
    acoustic_normal_incidence_interface_gate,
    balanced_back_emf_line_voltage_handoff_gate,
    box_projected_gradient_least_squares_gate,
    branch_line_hybrid_gate,
    carter_slot_opening_sweep_gate,
    coenergy_torque_periodic_summary,
    coaxial_rc_duality_gate,
    coaxial_pm_force_gap_sweep_gate,
    dq_current_from_gamma_deg,
    dq_torque_table_health,
    dq_to_three_phase_currents,
    drive_cycle_weighted_efficiency_gate,
    double_layer_winding_pitch_harmonic_gate,
    femm_block_label_source_contract_gate,
    femm_pm_magnetization_convention_gate,
    femm_static_current_circuit_rows_gate,
    farfield_pattern_metadata_gate,
    flux_linkage_back_emf_derivative_gate,
    ipm_saliency_torque_component_gate,
    inverter_dc_bus_voltage_limit_gate,
    jmag_motor_table_column_metadata_gate,
    jmag_symmetry_sweep_coverage_gate,
    lcurve_corner_choice,
    lumped_pm_dq_torque,
    motor_current_snapshot_table_contract_gate,
    morozov_discrepancy_choice,
    one_port_match_quality_gate,
    parallel_wire_force_per_length,
    pm_bem_surface_normal_metadata_gate,
    pm_drive_loss_bucket_efficiency_gate,
    pm_drive_terminal_table_health,
    pm_loadline_metadata_gate,
    pm_recoil_demag_step_summary,
    quarter_wave_directional_coupler_gate,
    shared_solver_session_health_gate,
    solver_result_table_metadata_gate,
    spwm_snapshot_current_handoff_summary,
    spherical_dirichlet_laplacian_eigen_gate,
    three_phase_currents_to_dq_summary,
    trace_surface_mass_energy_gate,
    thermal_annulus_conductance_gate,
    thermal_conduction_convection_robin_gate,
    thermal_layer_stack_conductance_gate,
    touchstone_frequency_grid_interpolation_gate,
    touchstone_port_metadata_gate,
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
    }
    gate = touchstone_port_metadata_gate(
        metadata,
        required_ports=("P1", "P2"),
        data_format="MA",
        frequency_unit="GHz",
        reference_impedance_ohm=50.0,
        port_order=("P1", "P2"),
    )

    assert gate["policy"] == "touchstone_port_metadata_gate"
    assert gate["status"] == "ok"
    assert gate["ports"] == ["P1", "P2"]
    assert gate["checks"]["port_order_matches_expected"] is True
    assert gate["checks"]["touchstone_format_matches_expected"] is True
    assert gate["checks"]["reference_impedance_matches_expected"] is True

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


def test_shared_solver_session_health_gate_separates_reuse_from_physics():
    gate = shared_solver_session_health_gate(
        connected=True,
        api_visible=True,
        discovered_engines=["MATLAB_10416"],
        status="already-connected",
        started_new_process=False,
        killed_process=False,
    )

    assert gate["policy"] == "shared_solver_session_health_gate"
    assert gate["status_label"] == "ok"
    assert gate["checks"]["session_connected"] is True
    assert gate["checks"]["api_visible"] is True
    assert gate["checks"]["engine_discovered"] is True
    assert gate["checks"]["status_allows_reuse"] is True
    assert gate["checks"]["started_no_new_process"] is True
    assert "preflight" in " ".join(gate["notes"])

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


def test_solver_result_table_metadata_gate_requires_columns_units_axis_and_rows():
    metadata = {
        "source": "COMSOL",
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
        min_rows=3,
    )

    assert gate["policy"] == "solver_result_table_metadata_gate"
    assert gate["status"] == "ok"
    assert gate["checks"]["required_columns_present"] is True
    assert gate["checks"]["expected_units_match"] is True
    assert gate["checks"]["independent_axis_is_column"] is True

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
    rows = [
        {"surface": "top", "area_m2": 2.0, "normal": [0.0, 0.0, 1.0], "magnetization": [0.0, 0.0, 1.0], "normal_convention": "outward_from_magnet"},
        {"surface": "bottom", "area_m2": 2.0, "normal": [0.0, 0.0, -1.0], "magnetization": [0.0, 0.0, 1.0], "normal_convention": "outward_from_magnet"},
        {"surface": "side_xp", "area_m2": 1.0, "normal": [1.0, 0.0, 0.0], "magnetization": [0.0, 0.0, 1.0], "normal_convention": "outward_from_magnet"},
        {"surface": "side_xm", "area_m2": 1.0, "normal": [-1.0, 0.0, 0.0], "magnetization": [0.0, 0.0, 1.0], "normal_convention": "outward_from_magnet"},
    ]

    gate = pm_bem_surface_normal_metadata_gate(rows)

    assert gate["status"] == "ok"
    assert gate["policy"] == "pm_bem_surface_normal_metadata_gate"
    assert gate["signed_charge_proxy_sum"] == pytest.approx(0.0, abs=1.0e-15)
    by_surface = {row["surface"]: row for row in gate["rows"]}
    assert by_surface["top"]["m_dot_n"] == pytest.approx(1.0)
    assert by_surface["bottom"]["m_dot_n"] == pytest.approx(-1.0)
    assert by_surface["side_xp"]["m_dot_n"] == pytest.approx(0.0)
    assert "magnetic-charge BEM assembly" in gate["version_note"]

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
