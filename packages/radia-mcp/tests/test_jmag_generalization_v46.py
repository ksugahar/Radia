from __future__ import annotations

import math

from radia_mcp.radia_ngsolve.jmag_v46_identity import validate_public_identity


TORQUE = "v46_public_transient_torque_sampling_periodic_window_partial_solution_mismatch"
THERMAL = "v46_public_thermal_coupling_unit_scale_temperature_coordinate_frame_mismatch"


def _identity():
    generation = "test-torque-v46"
    torque = {
        "generation": generation,
        **{key: generation for key in ("torque_generation", "sampling_generation", "periodic_window_generation", "partial_solution_generation", "result_generation")},
        "sample_times_s": [0.0, 1.0e-4, 2.0e-4], "result_sample_times_s": [0.0, 1.0e-4, 2.0e-4],
        "periodic_window_deg": 360.0, "result_periodic_window_deg": 360.0,
        "partial_transient_status": "complete", "result_partial_transient_status": "complete",
        "mesh_owner": "mesh:test-torque-v46", "result_mesh_owner": "mesh:test-torque-v46",
        "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64,
    }
    generation = "test-thermal-v46"
    thermal = {
        "generation": generation,
        **{key: generation for key in ("thermal_generation", "unit_generation", "temperature_generation", "coordinate_frame_generation", "result_generation")},
        "thermal_unit": "kelvin", "result_thermal_unit": "kelvin", "temperature_frame": "absolute_kelvin", "result_temperature_frame": "absolute_kelvin",
        "coordinate_frame": "global_cartesian", "result_coordinate_frame": "global_cartesian", "temperature_k": 353.15, "result_temperature_k": 353.15,
        "study_owner": "study:test-thermal-v46", "result_study_owner": "study:test-thermal-v46", "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
    }
    return {TORQUE: torque, THERMAL: thermal}


def test_v46_public_jmag_identity_accepts_closed_artifacts():
    checks = validate_public_identity(_identity())
    assert checks and all(checks.values())


def test_v46_public_jmag_identity_rejects_sampling_mutation():
    identity = _identity()
    identity[TORQUE]["result_periodic_window_deg"] = 180.0
    identity[TORQUE]["result_sample_times_s"] = [0.0]
    assert not all(validate_public_identity(identity).values())


def test_v46_public_jmag_identity_rejects_thermal_frame_mutation():
    identity = _identity()
    identity[THERMAL]["result_thermal_unit"] = "celsius"
    identity[THERMAL]["result_coordinate_frame"] = "rotor_local"
    assert not all(validate_public_identity(identity).values())


def test_v46_case_ids_are_frozen():
    assert TORQUE.startswith("v46_public_")
    assert THERMAL.startswith("v46_public_")
