from __future__ import annotations

import math

from radia_mcp.radia_ngsolve.magnetic_force_method_profile_gate import (
    magnetic_force_method_profile_gate,
)
from test_magnetic_force_generalization_v37_profile import _summary_v37


_PROMOTED_CASE_IDS = (
    "v38_public_eddy_current_maglev_plate_velocity_skin_depth_lift_drag_loss_power_mismatch",
    "v38_public_pm_coupling_torque_angle_periodicity_energy_derivative_action_reaction_mismatch",
)


def _summary_v38():
    summary = _summary_v37()
    identity = summary["artifact_identity"]
    generation = "eddy-maglev-258"
    velocity = 20.0
    pole_pitch = 0.1
    frequency = velocity / pole_pitch
    conductivity = 3.5e7
    relative_permeability = 1.0
    skin_depth = math.sqrt(
        2.0
        / (
            2.0
            * math.pi
            * frequency
            * (4.0e-7 * math.pi)
            * relative_permeability
            * conductivity
        )
    )
    lift = 100.0
    drag = 20.0
    drag_power = drag * velocity
    identity[
        "eddy_current_maglev_plate_velocity_frequency_conductivity_skin_depth_lift_drag_loss_power_mesh_owner_result_identity"
    ] = {
        "maglev_generation": generation,
        **{
            key: generation
            for key in (
                "velocity_generation", "frequency_generation",
                "conductivity_generation", "skin_generation", "force_generation",
                "loss_generation", "power_generation", "mesh_generation",
                "owner_generation", "result_generation",
            )
        },
        "plate_velocity_m_s": velocity,
        "result_plate_velocity_m_s": velocity,
        "pole_pitch_m": pole_pitch,
        "result_pole_pitch_m": pole_pitch,
        "excitation_frequency_hz": frequency,
        "result_excitation_frequency_hz": frequency,
        "plate_conductivity_s_m": conductivity,
        "result_plate_conductivity_s_m": conductivity,
        "relative_permeability": relative_permeability,
        "result_relative_permeability": relative_permeability,
        "skin_depth_m": skin_depth,
        "result_skin_depth_m": skin_depth,
        "lift_force_n": lift,
        "result_lift_force_n": lift,
        "drag_force_n": drag,
        "result_drag_force_n": drag,
        "joule_loss_w": drag_power,
        "result_joule_loss_w": drag_power,
        "mechanical_drag_power_w": drag_power,
        "result_mechanical_drag_power_w": drag_power,
        "power_balance_residual_w": 0.0,
        "result_power_balance_residual_w": 0.0,
        "power_tolerance_w": 1.0e-9,
        "result_power_tolerance_w": 1.0e-9,
        "mesh_owner": "mesh:eddy-maglev-258",
        "accepted_mesh_owner": "mesh:eddy-maglev-258",
        "maglev_result_sha256": "1" * 64,
        "accepted_maglev_result_sha256": "1" * 64,
    }

    generation = "pm-coupling-258"
    pole_pairs = 4
    period = 2.0 * math.pi / pole_pairs
    angle = math.pi / 16.0
    delta = 1.0e-4

    def energy(theta: float) -> float:
        return -0.5 * math.cos(pole_pairs * theta)

    energy_minus = energy(angle - delta)
    energy_center = energy(angle)
    energy_plus = energy(angle + delta)
    derivative_torque = -(energy_plus - energy_minus) / (2.0 * delta)
    identity[
        "pm_coupling_angle_pole_periodicity_energy_derivative_driver_driven_torque_action_reaction_frame_mesh_owner_result_identity"
    ] = {
        "coupling_generation": generation,
        **{
            key: generation
            for key in (
                "angle_generation", "periodicity_generation", "energy_generation",
                "derivative_generation", "torque_generation", "reaction_generation",
                "frame_generation", "mesh_generation", "owner_generation",
                "result_generation",
            )
        },
        "relative_angle_rad": angle,
        "result_relative_angle_rad": angle,
        "pole_pairs": pole_pairs,
        "result_pole_pairs": pole_pairs,
        "pole_period_rad": period,
        "result_pole_period_rad": period,
        "angle_perturbation_rad": delta,
        "result_angle_perturbation_rad": delta,
        "energy_minus_j": energy_minus,
        "result_energy_minus_j": energy_minus,
        "energy_center_j": energy_center,
        "result_energy_center_j": energy_center,
        "energy_plus_j": energy_plus,
        "result_energy_plus_j": energy_plus,
        "periodic_energy_j": energy(angle + period),
        "result_periodic_energy_j": energy(angle + period),
        "energy_derivative_torque_nm": derivative_torque,
        "result_energy_derivative_torque_nm": derivative_torque,
        "driver_torque_nm": derivative_torque,
        "result_driver_torque_nm": derivative_torque,
        "driven_torque_nm": -derivative_torque,
        "result_driven_torque_nm": -derivative_torque,
        "torque_frame": "relative_angle_driver_positive",
        "result_torque_frame": "relative_angle_driver_positive",
        "mesh_owner": "mesh:pm-coupling-258",
        "accepted_mesh_owner": "mesh:pm-coupling-258",
        "coupling_result_sha256": "2" * 64,
        "accepted_coupling_result_sha256": "2" * 64,
    }
    return summary


def test_v38_public_positive_eddy_maglev_and_pm_coupling_closure():
    assert magnetic_force_method_profile_gate(_summary_v38())["status"] == "ok"


def test_v38_public_eddy_current_maglev_plate_velocity_skin_depth_lift_drag_loss_power_mismatch():
    summary = _summary_v38()
    row = summary["artifact_identity"][
        "eddy_current_maglev_plate_velocity_frequency_conductivity_skin_depth_lift_drag_loss_power_mesh_owner_result_identity"
    ]
    row.update(
        {
            "skin_generation": "eddy-maglev-257",
            "power_generation": "eddy-maglev-256",
            "result_generation": "eddy-maglev-255",
            "result_plate_velocity_m_s": -20.0,
            "result_excitation_frequency_hz": -200.0,
            "result_plate_conductivity_s_m": -3.5e7,
            "result_skin_depth_m": -1.0,
            "result_lift_force_n": -100.0,
            "result_drag_force_n": -20.0,
            "result_joule_loss_w": -400.0,
            "result_mechanical_drag_power_w": 40.0,
            "result_power_balance_residual_w": 440.0,
            "accepted_mesh_owner": "stale:maglev",
            "accepted_maglev_result_sha256": "a" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "eddy_maglev_closes_velocity_frequency_skin_depth_lift_drag_joule_power_mesh_owner_and_result"
    ]


def test_v38_public_pm_coupling_torque_angle_periodicity_energy_derivative_action_reaction_mismatch():
    summary = _summary_v38()
    row = summary["artifact_identity"][
        "pm_coupling_angle_pole_periodicity_energy_derivative_driver_driven_torque_action_reaction_frame_mesh_owner_result_identity"
    ]
    row.update(
        {
            "periodicity_generation": "pm-coupling-257",
            "reaction_generation": "pm-coupling-256",
            "result_generation": "pm-coupling-255",
            "result_relative_angle_rad": -1.0,
            "result_pole_pairs": 0,
            "result_pole_period_rad": -1.0,
            "result_periodic_energy_j": 9.0,
            "result_energy_derivative_torque_nm": -9.0,
            "result_driver_torque_nm": 5.0,
            "result_driven_torque_nm": 5.0,
            "result_torque_frame": "left_handed_local",
            "accepted_mesh_owner": "stale:coupling",
            "accepted_coupling_result_sha256": "b" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "pm_coupling_closes_pole_periodic_energy_derivative_action_reaction_frame_mesh_owner_and_result"
    ]


def test_v38_public_rejects_self_consistent_wrong_maglev_skin_depth():
    summary = _summary_v38()
    row = summary["artifact_identity"][
        "eddy_current_maglev_plate_velocity_frequency_conductivity_skin_depth_lift_drag_loss_power_mesh_owner_result_identity"
    ]
    row["skin_depth_m"] *= 2.0
    row["result_skin_depth_m"] = row["skin_depth_m"]
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v38_public_rejects_self_consistent_same_sign_coupling_torque():
    summary = _summary_v38()
    row = summary["artifact_identity"][
        "pm_coupling_angle_pole_periodicity_energy_derivative_driver_driven_torque_action_reaction_frame_mesh_owner_result_identity"
    ]
    row["driven_torque_nm"] = row["driver_torque_nm"]
    row["result_driven_torque_nm"] = row["driven_torque_nm"]
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"
