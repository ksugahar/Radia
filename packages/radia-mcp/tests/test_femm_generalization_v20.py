from __future__ import annotations

from test_femm_generalization_v19 import _gate, _identity_v19
from test_force_coenergy_gate import _quadratic_case


def _identity_v20(sample_count):
    identity = _identity_v19(sample_count)
    solve_generations = ["solve-22-a", "solve-22-b", "solve-22-c"]
    identity["coenergy_torque_angle_difference_remesh_state_generation_identity"] = {
        "torque_generation": "torque-22",
        "derivative_torque_generation": "torque-22",
        "state_table_generation": "state-table-22",
        "derivative_state_table_generation": "state-table-22",
        "coenergy_solve_generations": solve_generations,
        "mesh_remap_solve_generations": solve_generations,
        "excitation_solve_generations": solve_generations,
        "angle_state_solve_generations": solve_generations,
        "angles_deg": [14.0, 15.0, 16.0],
        "derivative_angles_deg": [14.0, 15.0, 16.0],
        "angle_spacing_generation": "angle-spacing-22",
        "derivative_angle_spacing_generation": "angle-spacing-22",
        "coenergy_state_table_sha256": "1" * 64,
        "derivative_state_table_sha256": "1" * 64,
    }
    identity[
        "axisymmetric_henrotte_hodge_radius_weight_coordinate_generation_identity"
    ] = {
        "mesh_geometry_generation": "axi-mesh-22",
        "field_mesh_geometry_generation": "axi-mesh-22",
        "radius_weight_mesh_geometry_generation": "axi-mesh-22",
        "cylindrical_coordinate_mesh_geometry_generation": "axi-mesh-22",
        "node_ids": [101, 102, 103],
        "field_node_ids": [101, 102, 103],
        "radius_m": [0.01, 0.02, 0.03],
        "hodge_radius_weight_m": [0.01, 0.02, 0.03],
        "cylindrical_r_coordinate_m": [0.01, 0.02, 0.03],
        "radius_weight_table_sha256": "2" * 64,
        "force_radius_weight_table_sha256": "2" * 64,
        "coordinate_table_sha256": "3" * 64,
        "field_coordinate_table_sha256": "3" * 64,
    }
    return identity


def test_v20_public_positive_torque_state_and_axisymmetric_hodge_identity():
    positions, _, _ = _quadratic_case()
    result = _gate(_identity_v20(len(positions)))
    assert result["status"] == "ok"
    assert result["checks"][
        "coenergy_torque_states_share_remesh_excitation_and_angle_generations"
    ]
    assert result["checks"][
        "axisymmetric_hodge_radius_weights_use_current_mesh_coordinates"
    ]


def test_v20_public_coenergy_torque_angle_difference_remesh_state_generation_mismatch():
    positions, _, _ = _quadratic_case()
    identity = _identity_v20(len(positions))
    identity[
        "coenergy_torque_angle_difference_remesh_state_generation_identity"
    ].update(
        {
            "derivative_state_table_generation": "state-table-21",
            "mesh_remap_solve_generations": [
                "solve-21-a",
                "solve-21-b",
                "solve-21-c",
            ],
            "excitation_solve_generations": [
                "solve-22-a",
                "solve-21-b",
                "solve-22-c",
            ],
            "derivative_angles_deg": [14.0, 15.0, 17.0],
            "derivative_angle_spacing_generation": "angle-spacing-21",
            "derivative_state_table_sha256": "f" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "coenergy_torque_states_share_remesh_excitation_and_angle_generations"
    ] is False


def test_v20_public_axisymmetric_henrotte_hodge_radius_weight_coordinate_generation_mismatch():
    positions, _, _ = _quadratic_case()
    identity = _identity_v20(len(positions))
    identity[
        "axisymmetric_henrotte_hodge_radius_weight_coordinate_generation_identity"
    ].update(
        {
            "radius_weight_mesh_geometry_generation": "axi-mesh-21",
            "cylindrical_coordinate_mesh_geometry_generation": "axi-mesh-21",
            "hodge_radius_weight_m": [0.02, 0.03, 0.04],
            "cylindrical_r_coordinate_m": [0.03, 0.02, 0.01],
            "force_radius_weight_table_sha256": "f" * 64,
            "field_coordinate_table_sha256": "f" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "axisymmetric_hodge_radius_weights_use_current_mesh_coordinates"
    ] is False
