from __future__ import annotations

from radia_mcp.radia_ngsolve.force_coenergy_gate import force_coenergy_displacement_gate
from test_force_coenergy_gate import _artifact_identity, _quadratic_case


def _identity_v19(sample_count):
    identity = _artifact_identity(sample_count)
    identity["weighted_stress_air_mask_nodal_weight_mesh_generation_identity"] = {
        "field_mesh_generation": "mesh-21",
        "air_mask_mesh_generation": "mesh-21",
        "nodal_weight_mesh_generation": "mesh-21",
        "force_integral_mesh_generation": "mesh-21",
        "air_region_ids": [1, 2],
        "mask_air_region_ids": [1, 2],
        "nodal_weight_node_ids": [101, 102, 103],
        "force_weight_node_ids": [101, 102, 103],
        "air_mask_sha256": "a" * 64,
        "force_air_mask_sha256": "a" * 64,
        "nodal_weight_sha256": "b" * 64,
        "force_nodal_weight_sha256": "b" * 64,
    }
    identity["sliding_band_periodic_angle_rotor_position_generation_identity"] = {
        "rotor_position_generation": "rotor-position-21",
        "sliding_band_rotor_position_generation": "rotor-position-21",
        "torque_rotor_position_generation": "rotor-position-21",
        "periodic_angle_generation": "periodic-angle-21",
        "sliding_band_periodic_angle_generation": "periodic-angle-21",
        "torque_periodic_angle_generation": "periodic-angle-21",
        "rotor_angle_deg": 15.0,
        "sliding_band_rotor_angle_deg": 15.0,
        "periodic_angle_pairs_deg": [[0.0, 30.0], [30.0, 60.0]],
        "torque_periodic_angle_pairs_deg": [[0.0, 30.0], [30.0, 60.0]],
        "sliding_band_map_sha256": "c" * 64,
        "torque_sliding_band_map_sha256": "c" * 64,
    }
    return identity


def _gate(identity):
    positions, coenergy, forces = _quadratic_case()
    return force_coenergy_displacement_gate(
        positions, coenergy, forces, artifact_identity=identity
    )


def test_v19_public_positive_force_mask_and_sliding_band_identity():
    positions, _, _ = _quadratic_case()
    result = _gate(_identity_v19(len(positions)))
    assert result["status"] == "ok"
    assert result["checks"][
        "weighted_stress_air_mask_and_nodal_weights_use_current_mesh"
    ]
    assert result["checks"][
        "sliding_band_angles_and_torque_use_current_rotor_position"
    ]


def test_v19_public_weighted_stress_tensor_air_mask_mesh_generation_mismatch():
    positions, _, _ = _quadratic_case()
    identity = _identity_v19(len(positions))
    identity[
        "weighted_stress_air_mask_nodal_weight_mesh_generation_identity"
    ].update(
        {
            "air_mask_mesh_generation": "mesh-20",
            "nodal_weight_mesh_generation": "mesh-20",
            "mask_air_region_ids": [2, 1],
            "force_weight_node_ids": [101, 103, 102],
            "force_air_mask_sha256": "f" * 64,
            "force_nodal_weight_sha256": "f" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "weighted_stress_air_mask_and_nodal_weights_use_current_mesh"
    ] is False


def test_v19_public_sliding_band_periodic_angle_rotor_position_generation_mismatch():
    positions, _, _ = _quadratic_case()
    identity = _identity_v19(len(positions))
    identity[
        "sliding_band_periodic_angle_rotor_position_generation_identity"
    ].update(
        {
            "sliding_band_rotor_position_generation": "rotor-position-20",
            "torque_periodic_angle_generation": "periodic-angle-20",
            "sliding_band_rotor_angle_deg": 10.0,
            "torque_periodic_angle_pairs_deg": [[30.0, 60.0], [0.0, 30.0]],
            "torque_sliding_band_map_sha256": "f" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "sliding_band_angles_and_torque_use_current_rotor_position"
    ] is False
