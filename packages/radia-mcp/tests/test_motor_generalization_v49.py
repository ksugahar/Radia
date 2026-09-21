from __future__ import annotations

from copy import deepcopy

from radia_mcp.radia_ngsolve.motor_artifact_identity_v49 import DEMAG, IRON, validate_public_identity


PROMOTED_CASE_IDS = {
    "v49_public_demag_temperature_current_angle_irreversible_magnet_state_owner_mismatch",
    "v49_public_iron_loss_harmonic_time_frequency_hysteresis_eddy_excess_coefficient_owner_mismatch",
}


def _identity() -> dict[str, object]:
    demag_generation = "demag-operating-point-v49-901"
    iron_generation = "iron-loss-spectrum-v49-901"
    state = {"magnet:north": 0.998, "magnet:south": 0.997}
    harmonics = [
        {"order": 1, "frequency_hz": 400.0, "hysteresis_w": 8.0, "eddy_w": 3.0, "excess_w": 1.0},
        {"order": 3, "frequency_hz": 1200.0, "hysteresis_w": 1.5, "eddy_w": 1.2, "excess_w": 0.3},
    ]
    coefficients = {"hysteresis": 1.0, "eddy": 0.020, "excess": 0.004}
    return {
        DEMAG: {
            "generation": demag_generation,
            "temperature_generation": demag_generation,
            "current_generation": demag_generation,
            "angle_generation": demag_generation,
            "state_generation": demag_generation,
            "result_generation": demag_generation,
            "temperature_c": 120.0,
            "result_temperature_c": 120.0,
            "phase_current_a": [80.0, -40.0, -40.0],
            "result_phase_current_a": [80.0, -40.0, -40.0],
            "rotor_angle_electrical_deg": 90.0,
            "result_rotor_angle_electrical_deg": 90.0,
            "irreversible_magnet_state": state,
            "result_irreversible_magnet_state": state,
            "operating_point_owner": "operating-point:demag-v49-901",
            "result_operating_point_owner": "operating-point:demag-v49-901",
            "result_sha256": "a" * 64,
            "accepted_result_sha256": "a" * 64,
        },
        IRON: {
            "generation": iron_generation,
            "time_generation": iron_generation,
            "frequency_generation": iron_generation,
            "harmonic_generation": iron_generation,
            "coefficient_generation": iron_generation,
            "result_generation": iron_generation,
            "time_window_s": [0.0, 0.01],
            "result_time_window_s": [0.0, 0.01],
            "harmonic_rows": harmonics,
            "result_harmonic_rows": harmonics,
            "loss_coefficients": coefficients,
            "result_loss_coefficients": coefficients,
            "loss_owner": "loss-table:iron-v49-901",
            "result_loss_owner": "loss-table:iron-v49-901",
            "result_sha256": "b" * 64,
            "accepted_result_sha256": "b" * 64,
        },
    }


def test_v49_positive_demag_and_iron_loss_artifacts_are_accepted() -> None:
    assert all(validate_public_identity(_identity()).values())


def test_v49_demag_operating_point_and_state_mutations_are_rejected() -> None:
    identity = deepcopy(_identity())
    identity[DEMAG]["result_temperature_c"] = 20.0
    identity[DEMAG]["result_irreversible_magnet_state"] = {"magnet:north": 1.0, "magnet:south": 1.0}
    identity[DEMAG]["result_operating_point_owner"] = "operating-point:old"
    assert validate_public_identity(identity)["motor_v49_demag_temperature_current_angle_state_owner"] is False


def test_v49_iron_loss_window_harmonic_coefficient_and_owner_mutations_are_rejected() -> None:
    identity = deepcopy(_identity())
    identity[IRON]["result_time_window_s"] = [0.0, 0.005]
    identity[IRON]["result_harmonic_rows"] = list(reversed(identity[IRON]["harmonic_rows"]))
    identity[IRON]["result_loss_coefficients"] = {"hysteresis": 0.5, "eddy": 0.04, "excess": 0.0}
    identity[IRON]["result_loss_owner"] = "loss-table:old"
    assert validate_public_identity(identity)["motor_v49_iron_loss_harmonic_time_frequency_coefficient_owner"] is False
