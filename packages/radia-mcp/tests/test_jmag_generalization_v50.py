from __future__ import annotations

from copy import deepcopy

from radia_mcp.radia_ngsolve.motor_artifact_identity_v50 import DQ, THERMAL, validate_public_identity


PROMOTED_CASE_IDS = {
    "v50_public_dq_parameter_park_angle_current_phase_saliency_operating_point_owner_mismatch",
    "v50_public_thermal_loss_map_boundary_convection_temperature_material_owner_mismatch",
}


def _identity() -> dict[str, object]:
    dq_generation = "dq-operating-point-v50-901"
    thermal_generation = "electrothermal-v50-901"
    currents = {"id_a": -35.0, "iq_a": 82.0, "phase_order": "uvw"}
    saliency = {"ld_h": 0.0018, "lq_h": 0.0032, "psi_pm_wb": 0.092}
    loss = {"stator_iron_w": 42.0, "rotor_iron_w": 11.0, "copper_w": 68.0}
    convection = [{"boundary": "housing", "h_w_m2k": 18.0, "ambient_c": 25.0}]
    temperatures = {"winding_c": 96.0, "magnet_c": 78.0, "housing_c": 54.0}
    materials = {"copper": "cu:v3", "steel": "steel:v7", "magnet": "pm:v5"}
    return {
        DQ: {
            "generation": dq_generation, "angle_generation": dq_generation, "current_generation": dq_generation,
            "saliency_generation": dq_generation, "operating_point_generation": dq_generation, "result_generation": dq_generation,
            "park_angle_electrical_deg": 37.5, "result_park_angle_electrical_deg": 37.5,
            "dq_currents": currents, "result_dq_currents": currents,
            "saliency_parameters": saliency, "result_saliency_parameters": saliency,
            "operating_point_id": "operating-point:dq-v50-901", "result_operating_point_id": "operating-point:dq-v50-901",
            "result_owner": "dq-result:motor-v50-901", "accepted_result_owner": "dq-result:motor-v50-901",
            "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64,
        },
        THERMAL: {
            "generation": thermal_generation, "loss_generation": thermal_generation, "boundary_generation": thermal_generation,
            "temperature_generation": thermal_generation, "material_generation": thermal_generation, "result_generation": thermal_generation,
            "loss_map": loss, "replayed_loss_map": loss,
            "convection_boundaries": convection, "replayed_convection_boundaries": convection,
            "temperature_map": temperatures, "replayed_temperature_map": temperatures,
            "material_revisions": materials, "replayed_material_revisions": materials,
            "thermal_owner": "thermal-result:motor-v50-901", "replayed_thermal_owner": "thermal-result:motor-v50-901",
            "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
        },
    }


def test_v50_positive_dq_and_electrothermal_artifacts_are_accepted() -> None:
    assert all(validate_public_identity(_identity()).values())


def test_v50_dq_angle_current_saliency_operating_point_and_owner_drift_is_rejected() -> None:
    identity = deepcopy(_identity())
    identity[DQ]["result_park_angle_electrical_deg"] = 7.5
    identity[DQ]["result_dq_currents"] = {"id_a": 82.0, "iq_a": -35.0, "phase_order": "uwv"}
    identity[DQ]["result_operating_point_id"] = "operating-point:old"
    identity[DQ]["accepted_result_owner"] = "dq-result:foreign"
    assert validate_public_identity(identity)["motor_v50_dq_park_current_saliency_operating_point_owner"] is False


def test_v50_thermal_loss_boundary_temperature_material_and_owner_drift_is_rejected() -> None:
    identity = deepcopy(_identity())
    identity[THERMAL]["replayed_loss_map"] = {"copper_w": 34.0}
    identity[THERMAL]["replayed_convection_boundaries"] = [{"boundary": "shaft", "h_w_m2k": 5.0, "ambient_c": 40.0}]
    identity[THERMAL]["replayed_material_revisions"] = {"copper": "cu:old"}
    identity[THERMAL]["replayed_thermal_owner"] = "thermal-result:foreign"
    assert validate_public_identity(identity)["motor_v50_thermal_loss_convection_temperature_material_owner"] is False
