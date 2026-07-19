from __future__ import annotations

from radia_mcp.radia_ngsolve.motor_artifact_lineage_v47 import DQ, WINDOW, validate_public_identity


PROMOTED_CASE_IDS = {
    "v47_public_dq_phase_order_electrical_angle_pole_pair_mapping_mismatch",
    "v47_public_torque_loss_integration_window_parameter_row_key_mismatch",
}


def _identity() -> dict[str, object]:
    dq_generation = "dq-v47"
    window_generation = "window-v47"
    return {
        DQ: {
            "generation": dq_generation,
            **{
                key: dq_generation
                for key in (
                    "dq_generation",
                    "phase_order_generation",
                    "electrical_angle_generation",
                    "pole_pair_generation",
                    "transform_generation",
                    "result_generation",
                )
            },
            "phase_order": ["A", "B", "C"],
            "result_phase_order": ["A", "B", "C"],
            "electrical_angle_origin_deg": 0.0,
            "result_electrical_angle_origin_deg": 0.0,
            "pole_pairs": 4,
            "result_pole_pairs": 4,
            "transform_identity": "power_invariant_park",
            "result_transform_identity": "power_invariant_park",
            "angle_direction": "electrical_ccw",
            "result_angle_direction": "electrical_ccw",
            "result_sha256": "1" * 64,
            "accepted_result_sha256": "1" * 64,
        },
        WINDOW: {
            "generation": window_generation,
            **{
                key: window_generation
                for key in (
                    "integration_window_generation",
                    "parameter_row_generation",
                    "torque_generation",
                    "loss_generation",
                    "result_generation",
                )
            },
            "integration_window_s": [0.01, 0.02],
            "result_integration_window_s": [0.01, 0.02],
            "parameter_row_key": "speed=3000rpm,current=100A",
            "result_parameter_row_key": "speed=3000rpm,current=100A",
            "torque_mean_nm": 12.5,
            "result_torque_mean_nm": 12.5,
            "loss_total_w": 345.0,
            "result_loss_total_w": 345.0,
            "study_owner": "study:motor-v47",
            "result_study_owner": "study:motor-v47",
            "result_sha256": "2" * 64,
            "accepted_result_sha256": "2" * 64,
        },
    }


def test_v47_positive_motor_artifacts_are_accepted() -> None:
    assert all(validate_public_identity(_identity()).values())


def test_v47_dq_mapping_mutation_is_rejected() -> None:
    identity = _identity()
    identity[DQ]["result_phase_order"] = ["A", "C", "B"]
    identity[DQ]["result_electrical_angle_origin_deg"] = 30.0
    identity[DQ]["result_pole_pairs"] = 3
    identity[DQ]["result_transform_identity"] = "amplitude_invariant_park"
    assert not all(validate_public_identity(identity).values())


def test_v47_window_parameter_row_mutation_is_rejected() -> None:
    identity = _identity()
    identity[WINDOW]["result_integration_window_s"] = [0.02, 0.03]
    identity[WINDOW]["result_parameter_row_key"] = "speed=6000rpm,current=50A"
    assert not all(validate_public_identity(identity).values())
