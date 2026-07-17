from __future__ import annotations

from test_femm_generalization_v19 import _gate
from test_femm_generalization_v23 import _identity_v23
from test_force_coenergy_gate import _quadratic_case


def _identity_v24(sample_count):
    identity = _identity_v23(sample_count)
    identity[
        "weighted_stress_energy_derivative_force_mesh_frame_unit_generation_identity"
    ] = {
        "force_generation": "force-101",
        "weighted_stress_force_generation": "force-101",
        "energy_derivative_force_generation": "force-101",
        "mesh_force_generation": "force-101",
        "displacement_frame_force_generation": "force-101",
        "unit_force_generation": "force-101",
        "result_force_generation": "force-101",
        "weighted_stress_mesh_sha256": "1" * 64,
        "energy_derivative_mesh_sha256": "1" * 64,
        "displacement_frame_id": "world:x",
        "energy_derivative_displacement_frame_id": "world:x",
        "displacement_unit": "m",
        "energy_derivative_displacement_unit": "m",
        "force_unit": "N",
        "energy_derivative_force_unit": "N",
        "weighted_stress_force_n": [12.4, -0.2],
        "energy_derivative_force_n": [12.4, -0.2],
        "force_result_sha256": "2" * 64,
        "energy_derivative_result_sha256": "2" * 64,
    }
    identity[
        "axisymmetric_revolved_energy_force_2pir_jacobian_derham_generation_identity"
    ] = {
        "axisymmetric_generation": "axisym-101",
        "jacobian_axisymmetric_generation": "axisym-101",
        "field_axisymmetric_generation": "axisym-101",
        "material_axisymmetric_generation": "axisym-101",
        "mesh_axisymmetric_generation": "axisym-101",
        "revolved_result_axisymmetric_generation": "axisym-101",
        "jacobian_measure": "2*pi*r",
        "revolved_jacobian_measure": "2*pi*r",
        "field_state_sha256": "3" * 64,
        "revolved_field_state_sha256": "3" * 64,
        "material_map_sha256": "4" * 64,
        "revolved_material_map_sha256": "4" * 64,
        "axisymmetric_mesh_sha256": "5" * 64,
        "revolved_source_mesh_sha256": "5" * 64,
        "revolution_angle_deg": 360.0,
        "axisymmetric_energy_j": 0.125,
        "revolved_energy_j": 0.125,
        "axisymmetric_force_n": [4.2, 0.0],
        "revolved_force_n": [4.2, 0.0],
        "derham_sequence_id": "H1-axisym:Hodge-2pir",
        "revolved_derham_sequence_id": "H1-axisym:Hodge-2pir",
        "result_sha256": "6" * 64,
        "revolved_result_sha256": "6" * 64,
    }
    return identity


def test_v24_public_positive_force_and_axisymmetric_revolution_identity():
    positions, _, _ = _quadratic_case()
    assert _gate(_identity_v24(len(positions)))["status"] == "ok"


def test_v24_public_weighted_stress_force_energy_derivative_mesh_frame_unit_mismatch():
    positions, _, _ = _quadratic_case()
    identity = _identity_v24(len(positions))
    identity[
        "weighted_stress_energy_derivative_force_mesh_frame_unit_generation_identity"
    ].update(
        {
            "energy_derivative_force_generation": "force-100",
            "mesh_force_generation": "force-99",
            "displacement_frame_force_generation": "force-98",
            "unit_force_generation": "force-97",
            "energy_derivative_mesh_sha256": "b" * 64,
            "energy_derivative_displacement_frame_id": "rotor:x",
            "energy_derivative_displacement_unit": "mm",
            "energy_derivative_force_unit": "N/mm",
            "energy_derivative_force_n": [0.0124, -0.0002],
            "energy_derivative_result_sha256": "c" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "weighted_stress_and_energy_derivative_share_mesh_frame_units_and_generation"
    ]


def test_v24_public_axisymmetric_revolved_energy_force_2pir_jacobian_derham_generation_mismatch():
    positions, _, _ = _quadratic_case()
    identity = _identity_v24(len(positions))
    identity[
        "axisymmetric_revolved_energy_force_2pir_jacobian_derham_generation_identity"
    ].update(
        {
            "jacobian_axisymmetric_generation": "axisym-100",
            "field_axisymmetric_generation": "axisym-99",
            "material_axisymmetric_generation": "axisym-98",
            "mesh_axisymmetric_generation": "axisym-97",
            "revolved_result_axisymmetric_generation": "axisym-96",
            "revolved_jacobian_measure": "r",
            "revolved_field_state_sha256": "d" * 64,
            "revolved_material_map_sha256": "e" * 64,
            "revolved_source_mesh_sha256": "f" * 64,
            "revolution_angle_deg": 180.0,
            "revolved_energy_j": 0.0625,
            "revolved_force_n": [2.1, 0.0],
            "revolved_derham_sequence_id": "plain-H1",
            "revolved_result_sha256": "0" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "axisymmetric_revolved_energy_force_share_2pir_hodge_field_material_and_mesh"
    ]
