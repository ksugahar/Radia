from __future__ import annotations

from radia_mcp.radia_ngsolve.magnetic_force_artifact_lineage_v47 import (
    FORCE,
    MOTOR,
    MOTOR_LANES,
    validate_public_v47_identity,
)


PROMOTED_CASE_IDS = {
    "v47_public_motor_dual_lane_geometry_material_excitation_operating_point_identity_mismatch",
    "v47_public_force_coenergy_displacement_pair_body_owner_aggregation_mismatch",
}


def _identity() -> dict[str, object]:
    motor_generation = "motor-v47"
    force_generation = "force-v47"
    excitation = {"phase_order": ["A", "B", "C"], "current_a": [10.0, -5.0, -5.0]}
    components = {"core": 80.0, "magnet": 20.0}
    return {
        MOTOR: {
            "generation": motor_generation,
            **{
                key: motor_generation
                for key in (
                    "geometry_generation",
                    "material_generation",
                    "excitation_generation",
                    "operating_point_generation",
                    "lane_a_generation",
                    "lane_b_generation",
                    "result_generation",
                )
            },
            "lane_ids": MOTOR_LANES,
            "result_lane_ids": MOTOR_LANES,
            "geometry_identity_sha256": "1" * 64,
            "lane_a_geometry_identity_sha256": "1" * 64,
            "lane_b_geometry_identity_sha256": "1" * 64,
            "material_identity": "material:motor-v47",
            "lane_a_material_identity": "material:motor-v47",
            "lane_b_material_identity": "material:motor-v47",
            "excitation_identity": excitation,
            "lane_a_excitation_identity": excitation,
            "lane_b_excitation_identity": excitation,
            "operating_point_key": "speed=3000rpm,current=10A",
            "lane_a_operating_point_key": "speed=3000rpm,current=10A",
            "lane_b_operating_point_key": "speed=3000rpm,current=10A",
            "result_sha256": "2" * 64,
            "accepted_result_sha256": "2" * 64,
        },
        FORCE: {
            "generation": force_generation,
            **{
                key: force_generation
                for key in (
                    "coenergy_generation",
                    "displacement_generation",
                    "body_owner_generation",
                    "aggregation_generation",
                    "result_generation",
                )
            },
            "displacement_pair_m": [0.0, 0.001],
            "result_displacement_pair_m": [0.0, 0.001],
            "coenergy_pair_j": [1.0, 1.1],
            "result_coenergy_pair_j": [1.0, 1.1],
            "body_owner": "body:moving-assembly",
            "result_body_owner": "body:moving-assembly",
            "component_force_n": components,
            "result_component_force_n": components,
            "aggregated_force_n": 100.0,
            "result_aggregated_force_n": 100.0,
            "result_sha256": "3" * 64,
            "accepted_result_sha256": "3" * 64,
        },
    }


def test_v47_positive_dual_lane_and_force_artifacts_are_accepted() -> None:
    assert all(validate_public_v47_identity(_identity()).values())


def test_v47_dual_lane_shared_physics_mutation_is_rejected() -> None:
    identity = _identity()
    identity[MOTOR]["lane_b_geometry_identity_sha256"] = "a" * 64
    identity[MOTOR]["lane_b_material_identity"] = "material:other"
    identity[MOTOR]["lane_b_excitation_identity"] = {
        "phase_order": ["A", "C", "B"],
        "current_a": [10.0, -5.0, -5.0],
    }
    identity[MOTOR]["lane_b_operating_point_key"] = "speed=6000rpm,current=5A"
    assert not all(validate_public_v47_identity(identity).values())


def test_v47_force_pair_body_aggregation_mutation_is_rejected() -> None:
    identity = _identity()
    identity[FORCE]["result_displacement_pair_m"] = [0.001, 0.0]
    identity[FORCE]["result_body_owner"] = "body:fixed"
    identity[FORCE]["result_component_force_n"] = {"core": 80.0}
    identity[FORCE]["result_aggregated_force_n"] = 80.0
    assert not all(validate_public_v47_identity(identity).values())
