from __future__ import annotations

import copy
import json
import math

import pytest

from radia_mcp.radia_ngsolve.rotational_eddy_brake_energy_gate import (
    rotational_eddy_brake_energy_gate as gate,
)
from radia_mcp.radia_ngsolve.server import rotational_eddy_brake_energy_gate as mcp_gate


def _summary() -> dict:
    density, radius, thickness = 7800.0, 0.1, 0.01
    inertia = 0.5 * density * math.pi * radius**4 * thickness
    times = [0.05 * index for index in range(201)]
    rate = 0.2
    omega = [100.0 * math.exp(-rate * value) for value in times]
    torque = [inertia * rate * value for value in omega]
    joule = [force * speed for force, speed in zip(torque, omega, strict=True)]
    generation = "solve-generation-42"
    frame = "laboratory-frame-generation-7"

    def replay(label: str) -> dict:
        return {
            "label": label,
            "solve_seconds": 2.0,
            "time_s": list(times),
            "angular_velocity_rad_s": list(omega),
            "braking_torque_nm": list(torque),
            "joule_loss_w": list(joule),
            "artifact_generations": {
                "time_s": generation,
                "angular_velocity_rad_s": generation,
                "braking_torque_nm": generation,
                "joule_loss_w": generation,
            },
            "artifact_coordinate_frames": {
                "time_s": frame,
                "angular_velocity_rad_s": frame,
                "braking_torque_nm": frame,
                "joule_loss_w": frame,
            },
        }

    boundary = 100
    accumulated_joule = sum(
        0.5 * (joule[index] + joule[index + 1])
        * (times[index + 1] - times[index])
        for index in range(boundary)
    )
    stored_energy = 0.5 * inertia * omega[boundary] ** 2 + 0.3

    return {
        "contract": {
            "body": "uniform_solid_conducting_disc",
            "inertia_reference": "analytic_uniform_solid_disc",
            "angular_momentum_balance": "inertia_delta_angular_velocity_plus_integrated_braking_torque_equals_zero",
            "instantaneous_power_comparison": "diagnostic_only_when_field_energy_rate_is_not_sampled_on_the_probe_grid",
            "energy_balance": "initial_kinetic_plus_magnetic_equals_final_kinetic_plus_magnetic_plus_joule",
        },
        "units": {
            "time": "s",
            "angular_velocity": "rad/s",
            "torque": "N*m",
            "power": "W",
            "inertia": "kg*m^2",
            "energy": "J",
            "density": "kg/m^3",
            "length": "m",
        },
        "disc": {
            "density_kg_m3": density,
            "radius_m": radius,
            "thickness_m": thickness,
        },
        "reported_inertia_kg_m2": inertia,
        "replays": [replay("one"), replay("two")],
        "energy_replay": {
            **replay("energy"),
            "field_energy_time_s": list(times),
            "magnetic_energy_j": [0.3 for _ in times],
            "artifact_generations": {
                "time_s": generation,
                "angular_velocity_rad_s": generation,
                "braking_torque_nm": generation,
                "joule_loss_w": generation,
                "field_energy_time_s": generation,
                "magnetic_energy_j": generation,
            },
            "artifact_coordinate_frames": {
                "time_s": frame,
                "angular_velocity_rad_s": frame,
                "braking_torque_nm": frame,
                "joule_loss_w": frame,
                "field_energy_time_s": frame,
                "magnetic_energy_j": frame,
            },
            "restart_boundaries": [
                {
                    "left_index": boundary,
                    "right_index": boundary + 1,
                    "generation_before": generation,
                    "generation_after": "solve-generation-42-restart-1",
                    "stored_energy_before_j": stored_energy,
                    "stored_energy_after_j": stored_energy,
                    "accumulated_joule_before_j": accumulated_joule,
                    "accumulated_joule_offset_after_j": accumulated_joule,
                }
            ],
        },
        "convergence_provenance": {
            "solution_generation": "solution-generation-42",
            "result_iteration_generation": "nonlinear-iteration-8",
            "convergence_table_iteration_generation": "nonlinear-iteration-8",
            "terminal_state": "converged",
            "terminal_relative_residual": 2.0e-9,
        },
        "timing_breakdown_s": {
            "attach": 0.1,
            "solve": 6.0,
            "extract": 0.1,
            "verify": 0.1,
        },
    }


def test_accepts_free_brake_with_field_energy_and_mcp_dispatch() -> None:
    result = gate(_summary())
    assert result["status"] == "ok"
    assert result["checks"]["kinetic_magnetic_joule_energy_closes"] is True
    assert json.loads(mcp_gate(json.dumps(_summary())))["status"] == "ok"


def test_rejects_wrong_inertia_and_angular_impulse() -> None:
    summary = copy.deepcopy(_summary())
    summary["reported_inertia_kg_m2"] *= 1.25
    summary["contract"]["inertia_reference"] = "unchecked_reported_value"
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert result["checks"]["analytic_disc_inertia_matches_reported"] is False
    assert result["checks"]["angular_impulse_balance_closes"] is False


def test_rejects_missing_field_storage_and_false_power_claim() -> None:
    summary = copy.deepcopy(_summary())
    summary["contract"]["instantaneous_power_comparison"] = "always_equal_to_joule"
    summary["energy_replay"]["magnetic_energy_j"] = []
    try:
        result = gate(summary)
    except ValueError as exc:
        assert "magnetic_energy_j" in str(exc)
    else:
        assert result["status"] == "needs_attention"


def test_rejects_nonreplaying_torque_history() -> None:
    summary = copy.deepcopy(_summary())
    summary["replays"][1]["braking_torque_nm"] = list(
        summary["replays"][1]["braking_torque_nm"]
    )
    summary["replays"][1]["braking_torque_nm"][50] *= 1.2
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert result["checks"]["fresh_replay_fields_match"] is False


def test_rejects_joule_history_inconsistent_with_motion_energy() -> None:
    summary = copy.deepcopy(_summary())
    summary["energy_replay"]["joule_loss_w"] = [
        1.4 * value for value in summary["energy_replay"]["joule_loss_w"]
    ]
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert result["checks"]["kinetic_magnetic_joule_energy_closes"] is False


def test_rejects_isolated_magnetic_energy_spike() -> None:
    summary = copy.deepcopy(_summary())
    summary["energy_replay"]["magnetic_energy_j"][10] *= 4.0
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "field_energy_history_is_nonnegative_and_has_no_isolated_jump"
        ]
        is False
    )


@pytest.mark.parametrize(
    "case_id",
    ["disc_radius", "momentum_contract", "torque_replay", "torque_unit", "field_energy"],
)
def test_counterfactual_curriculum90_public(case_id: str) -> None:
    summary = copy.deepcopy(_summary())
    if case_id == "disc_radius":
        summary["disc"]["radius_m"] *= 1.2
    elif case_id == "momentum_contract":
        summary["contract"]["angular_momentum_balance"] = "unchecked"
    elif case_id == "torque_replay":
        summary["replays"][1]["braking_torque_nm"] = list(
            summary["replays"][1]["braking_torque_nm"]
        )
        summary["replays"][1]["braking_torque_nm"][25] *= 1.25
    elif case_id == "torque_unit":
        summary["units"]["torque"] = "kg"
    else:
        summary["energy_replay"]["magnetic_energy_j"][10] *= 4.0
    assert gate(summary)["status"] == "needs_attention"


def test_generalization_v3s_rejects_negative_field_energy() -> None:
    summary = copy.deepcopy(_summary())
    summary["energy_replay"]["magnetic_energy_j"][17] = -0.1
    assert gate(summary)["status"] == "needs_attention"


@pytest.mark.parametrize(
    "case_id",
    [
        "v4_time_duplicate",
        "v4_velocity_local_rise",
        "v4_positive_braking_torque",
        "v4_negative_joule_loss",
        "v4_field_time_shift",
    ],
)
def test_counterfactual_curriculum90_v4_public(case_id: str) -> None:
    summary = copy.deepcopy(_summary())
    if case_id == "v4_time_duplicate":
        summary["replays"][0]["time_s"][20] = summary["replays"][0]["time_s"][19]
    elif case_id == "v4_velocity_local_rise":
        summary["replays"][0]["angular_velocity_rad_s"][25] *= 1.5
    elif case_id == "v4_positive_braking_torque":
        summary["replays"][0]["braking_torque_nm"][20] = 1.0
    elif case_id == "v4_negative_joule_loss":
        summary["energy_replay"]["joule_loss_w"][12] = -1.0
    else:
        summary["energy_replay"]["field_energy_time_s"][30] *= 1.1
    assert gate(summary)["status"] == "needs_attention"


def test_generalization_v5_rejects_non_si_energy_unit() -> None:
    summary = copy.deepcopy(_summary())
    summary["units"]["energy"] = "mJ"
    assert gate(summary)["status"] == "needs_attention"


@pytest.mark.parametrize(
    "case_id",
    ["v6_public_energy_state_drift", "v6_public_replay_joule_drift"],
)
def test_generalization_v6_public(case_id: str) -> None:
    summary = copy.deepcopy(_summary())
    if case_id == "v6_public_energy_state_drift":
        summary["energy_replay"]["magnetic_energy_j"][20] *= 1.10
    else:
        summary["replays"][1]["joule_loss_w"] = list(
            summary["replays"][1]["joule_loss_w"]
        )
        summary["replays"][1]["joule_loss_w"][20] *= 1.05
    assert gate(summary)["status"] == "needs_attention"


@pytest.mark.parametrize(
    "case_id",
    [
        "v7_public_segment_overlap_after_restart",
        "v7_public_field_loss_timebase_skew",
    ],
)
def test_generalization_v7_public(case_id: str) -> None:
    summary = copy.deepcopy(_summary())
    if case_id == "v7_public_segment_overlap_after_restart":
        for replay in [*summary["replays"], summary["energy_replay"]]:
            replay["time_s"][101] = replay["time_s"][99]
        summary["energy_replay"]["field_energy_time_s"][101] = summary[
            "energy_replay"
        ]["field_energy_time_s"][99]
    else:
        for index in range(60, 140):
            summary["energy_replay"]["field_energy_time_s"][index] += 0.0125
    assert gate(summary)["status"] == "needs_attention"


@pytest.mark.parametrize(
    "case_id",
    [
        "v8_public_restart_energy_baseline_reset",
        "v8_public_loss_history_prior_solve_generation",
    ],
)
def test_generalization_v8_public(case_id: str) -> None:
    summary = copy.deepcopy(_summary())
    if case_id == "v8_public_restart_energy_baseline_reset":
        summary["energy_replay"]["restart_boundaries"][0][
            "accumulated_joule_offset_after_j"
        ] = 0.0
        expected_check = "restart_energy_offsets_are_continuous"
    else:
        summary["energy_replay"]["artifact_generations"][
            "joule_loss_w"
        ] = "solve-generation-41"
        expected_check = "artifact_series_share_their_solve_generation"
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert result["checks"][expected_check] is False


@pytest.mark.parametrize(
    "case_id",
    [
        "v9_public_torque_energy_coordinate_frame_mismatch",
        "v9_public_convergence_table_prior_nonlinear_iteration",
    ],
)
def test_generalization_v9_public(case_id: str) -> None:
    summary = copy.deepcopy(_summary())
    if case_id == "v9_public_torque_energy_coordinate_frame_mismatch":
        summary["energy_replay"]["artifact_coordinate_frames"][
            "braking_torque_nm"
        ] = "rotor-local-frame-generation-7"
        expected_check = "artifact_series_share_one_coordinate_frame"
    else:
        summary["convergence_provenance"][
            "convergence_table_iteration_generation"
        ] = "nonlinear-iteration-7"
        expected_check = "convergence_table_matches_result_iteration"
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert result["checks"][expected_check] is False


def _with_v10_ownership(summary: dict) -> dict:
    summary["force_selection_identity"] = {
        "geometry_generation": "geometry-generation-12",
        "solution_geometry_generation": "geometry-generation-12",
        "integration_selection_generation": "geometry-generation-12",
        "selection_entity_digest": "d" * 64,
        "force_result_selection_digest": "d" * 64,
    }
    summary["excitation_basis_identity"] = {
        "sweep_generation": "sweep-generation-42",
        "solve_amplitude_basis": "rms",
        "extract_amplitude_basis": "rms",
        "solve_scale_to_rms": 1.0,
        "extract_scale_to_rms": 1.0,
        "torque_loss_normalization_basis": "rms_excitation",
    }
    return summary


def test_v10_public_force_selection_generation_changed() -> None:
    summary = _with_v10_ownership(copy.deepcopy(_summary()))
    summary["force_selection_identity"].update(
        {
            "integration_selection_generation": "geometry-generation-13",
            "selection_entity_digest": "c" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"]["force_integral_uses_current_geometry_selection"]
        is False
    )


def test_v10_public_excitation_peak_rms_sweep_basis_mismatch() -> None:
    summary = _with_v10_ownership(copy.deepcopy(_summary()))
    summary["excitation_basis_identity"].update(
        {
            "extract_amplitude_basis": "peak",
            "extract_scale_to_rms": 2.0**-0.5,
            "torque_loss_normalization_basis": "compensated_mixed_basis",
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert result["checks"]["sweep_excitation_uses_one_rms_basis"] is False
