from __future__ import annotations

import math

from radia_mcp.radia_ngsolve.pwm_controlled_motor_loss_gate import pwm_controlled_motor_loss_gate
from test_motor_generalization_v40 import _payload_v40


_SRM = "srm_inductance_position_current_coenergy_torque_ripple_power_model_result_generation_identity"
_AXIAL = "axial_flux_sector_periodicity_skew_end_effect_flux_torque_loss_power_mesh_result_generation_identity"
_PROMOTED_CASE_IDS = (
    "v41_public_srm_inductance_position_current_coenergy_torque_ripple_power_mismatch",
    "v41_public_axialflux_sector_periodicity_skew_endeffect_flux_torque_loss_power_mismatch",
)


def _derivative(values: list[float], coordinates: list[float]) -> list[float]:
    result = []
    for index in range(len(values)):
        if index == 0:
            left, right = 0, 1
        elif index == len(values) - 1:
            left, right = index - 1, index
        else:
            left, right = index - 1, index + 1
        result.append(
            (values[right] - values[left])
            / (coordinates[right] - coordinates[left])
        )
    return result


def _payload_v41():
    payload = _payload_v40()
    identity = payload["artifact_identity"]

    generation = "srm-map-724"
    positions = [0.0, 0.1, 0.2, 0.3, 0.4]
    currents = [5.0, 10.0]
    inductance = [
        [0.0105, 0.0125, 0.0150, 0.0170, 0.0185],
        [0.0100, 0.0120, 0.0145, 0.0165, 0.0180],
    ]
    selected_current = 10.0
    coenergy = [0.5 * item * selected_current**2 for item in inductance[1]]
    torque = _derivative(coenergy, positions)
    average_torque = sum(torque) / len(torque)
    phase_count = 3
    phase_resistance = 0.05
    copper_loss = phase_count * phase_resistance * selected_current**2
    mechanical_speed = 100.0
    values = {
        "rotor_position_rad": positions,
        "current_samples_a": currents,
        "inductance_h_by_current": inductance,
        "selected_current_a": selected_current,
        "coenergy_j": coenergy,
        "torque_nm": torque,
        "average_torque_nm": average_torque,
        "torque_ripple_nm": max(torque) - min(torque),
        "phase_count": phase_count,
        "phase_resistance_ohm": phase_resistance,
        "copper_loss_w": copper_loss,
        "mechanical_speed_rad_s": mechanical_speed,
        "mechanical_power_w": average_torque * mechanical_speed,
        "electrical_power_w": average_torque * mechanical_speed + copper_loss,
    }
    identity[_SRM] = {
        "srm_generation": generation,
        **{
            key: generation
            for key in (
                "position_generation", "current_generation",
                "inductance_generation", "coenergy_generation",
                "torque_generation", "ripple_generation", "loss_generation",
                "power_generation", "model_generation", "result_generation",
            )
        },
        **values,
        **{f"result_{key}": value for key, value in values.items()},
        "model_owner": "motor:srm-map-724",
        "accepted_model_owner": "motor:srm-map-724",
        "srm_result_sha256": "5" * 64,
        "accepted_srm_result_sha256": "5" * 64,
    }

    generation = "axial-flux-724"
    sector_count = 8
    pole_pairs = 4
    skew_angle = 2.0
    skew_argument = pole_pairs * math.radians(skew_angle) / 2.0
    skew_factor = math.sin(skew_argument) / skew_argument
    end_effect = 0.95
    uncorrected_flux = 0.5
    corrected_flux = uncorrected_flux * skew_factor * end_effect
    full_torque = 26.2
    current_q = full_torque / (1.5 * pole_pairs * corrected_flux)
    phase_resistance = 0.1
    copper_loss = 3.0 * phase_resistance * current_q**2
    iron_loss = 50.0
    mechanical_speed = 100.0
    values = {
        "sector_count": sector_count,
        "sector_angle_deg": 360.0 / sector_count,
        "periodicity": "periodic",
        "pole_pairs": pole_pairs,
        "skew_angle_deg": skew_angle,
        "skew_factor": skew_factor,
        "end_effect_factor": end_effect,
        "uncorrected_airgap_flux_wb": uncorrected_flux,
        "corrected_airgap_flux_wb": corrected_flux,
        "current_q_a": current_q,
        "sector_torque_nm": full_torque / sector_count,
        "full_torque_nm": full_torque,
        "phase_resistance_ohm": phase_resistance,
        "copper_loss_w": copper_loss,
        "iron_loss_w": iron_loss,
        "mechanical_speed_rad_s": mechanical_speed,
        "mechanical_power_w": full_torque * mechanical_speed,
        "electrical_power_w": full_torque * mechanical_speed + copper_loss + iron_loss,
    }
    identity[_AXIAL] = {
        "axial_flux_generation": generation,
        **{
            key: generation
            for key in (
                "sector_generation", "periodicity_generation", "skew_generation",
                "end_effect_generation", "flux_generation", "torque_generation",
                "loss_generation", "power_generation", "mesh_generation",
                "result_generation",
            )
        },
        **values,
        **{f"result_{key}": value for key, value in values.items()},
        "mesh_owner": "mesh:axial-flux-724",
        "accepted_mesh_owner": "mesh:axial-flux-724",
        "axial_flux_result_sha256": "6" * 64,
        "accepted_axial_flux_result_sha256": "6" * 64,
    }
    return payload


def test_v41_public_positive_srm_and_axial_flux_closure():
    assert pwm_controlled_motor_loss_gate(_payload_v41())["status"] == "ok"


def test_v41_public_srm_inductance_position_current_coenergy_torque_ripple_power_mismatch():
    payload = _payload_v41()
    payload["artifact_identity"][_SRM].update(
        {
            "inductance_generation": "srm-map-723",
            "result_rotor_position_rad": [0.4, 0.3, 0.2, 0.1, 0.0],
            "result_torque_nm": [-1.0],
            "accepted_model_owner": "stale:model",
        }
    )
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v41_public_axialflux_sector_periodicity_skew_endeffect_flux_torque_loss_power_mismatch():
    payload = _payload_v41()
    payload["artifact_identity"][_AXIAL].update(
        {
            "periodicity_generation": "axial-flux-723",
            "result_periodicity": "antiperiodic",
            "result_corrected_airgap_flux_wb": -1.0,
            "accepted_mesh_owner": "stale:mesh",
        }
    )
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v41_public_rejects_self_consistent_wrong_srm_coenergy():
    payload = _payload_v41()
    row = payload["artifact_identity"][_SRM]
    row["coenergy_j"] = row["result_coenergy_j"] = [1.0] * 5
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v41_public_rejects_self_consistent_wrong_axial_flux_power():
    payload = _payload_v41()
    row = payload["artifact_identity"][_AXIAL]
    row["electrical_power_w"] = row["result_electrical_power_w"] = row["mechanical_power_w"]
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"
