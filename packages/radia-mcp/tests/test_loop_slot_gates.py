import math

import pytest

from radia_mcp.radia_ngsolve.slot_gates import (
    acoustic_plane_wave_intensity_convention_gate,
    box_projected_gradient_least_squares_gate,
    branch_line_hybrid_gate,
    coenergy_torque_periodic_summary,
    dq_current_from_gamma_deg,
    dq_torque_table_health,
    dq_to_three_phase_currents,
    lumped_pm_dq_torque,
    parallel_wire_force_per_length,
    pm_recoil_demag_step_summary,
    quarter_wave_directional_coupler_gate,
    spwm_snapshot_current_handoff_summary,
    three_phase_currents_to_dq_summary,
    two_port_sparameter_health,
)


def test_parallel_wire_force_gate_signed_and_scaled():
    f = parallel_wire_force_per_length(10.0, 20.0, 0.05)
    assert f == pytest.approx(8.0e-4)
    assert parallel_wire_force_per_length(10.0, -20.0, 0.05) == pytest.approx(-f)
    assert parallel_wire_force_per_length(10.0, 20.0, 0.10) == pytest.approx(0.5 * f)
    with pytest.raises(ValueError):
        parallel_wire_force_per_length(1.0, 1.0, 0.0)


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
