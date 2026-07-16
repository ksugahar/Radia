import json
import math

import pytest

from radia_mcp.radia_ngsolve.force_coenergy_gate import force_coenergy_displacement_gate
from radia_mcp.radia_ngsolve.server import force_coenergy_displacement_gate as mcp_gate


def _quadratic_case():
    positions = [0.002 * index for index in range(7)]
    coenergy = [2.0 - 40.0 * x + 100.0 * x * x for x in positions]
    forces = [-40.0 + 200.0 * x for x in positions]
    return positions, coenergy, forces


def _artifact_identity(sample_count):
    return {
        "direct_force_snapshot": {
            "load_step_id": "load-step-42",
            "time_s": 0.025,
        },
        "coenergy_derivative_snapshot": {
            "load_step_id": "load-step-42",
            "time_s": 0.025,
        },
        "coenergy_mesh_family_generations": ["mesh-family-7"] * sample_count,
        "displacement_axis": {
            "numeric_unit": "m",
            "derivative_unit": "m",
            "scale_to_si": 1.0,
        },
        "force_frame": {
            "direct_frame_id": "model-frame",
            "derivative_frame_id": "model-frame",
            "direct_axis": [0.0, -1.0, 0.0],
            "derivative_axis": [0.0, -1.0, 0.0],
            "reflection_applied": True,
        },
        "force_normalization": {
            "formulation": "axisymmetric",
            "solver_result_scope": "total_3d_force",
            "reported_result_scope": "total_3d_force",
            "revolution_factor_application_count": 0,
        },
        "force_body_selection": {
            "target_group_ids": [1],
            "weighted_stress_selected_group_ids": [1],
            "material_roles": {"0": "air", "1": "magnetic_body"},
            "excluded_air_group_ids": [0],
            "selection_generation": "selection-12",
        },
    }


def test_force_coenergy_gate_accepts_constant_current_virtual_work_identity():
    positions, coenergy, forces = _quadratic_case()
    result = force_coenergy_displacement_gate(positions, coenergy, forces)
    assert result["status"] == "ok"
    assert result["max_central_relative_error"] < 1.0e-12
    assert result["endpoint_errors_are_diagnostic_only"] is True
    assert result["rows"][0]["stencil"] == "forward"


def test_force_coenergy_gate_accepts_bound_snapshot_and_mesh_family_identity():
    positions, coenergy, forces = _quadratic_case()
    result = force_coenergy_displacement_gate(
        positions,
        coenergy,
        forces,
        artifact_identity=_artifact_identity(len(positions)),
    )
    assert result["status"] == "ok"
    assert result["warnings"] == []


def test_force_coenergy_gate_rejects_force_with_wrong_projection_sign():
    positions, coenergy, forces = _quadratic_case()
    result = force_coenergy_displacement_gate(positions, coenergy, [-f for f in forces])
    assert result["status"] == "needs_attention"
    assert result["checks"]["central_virtual_work_matches_direct_force"] is False


def test_force_coenergy_mcp_tool_dispatches_json_and_handles_bad_shape():
    positions, coenergy, forces = _quadratic_case()
    result = json.loads(mcp_gate(positions, coenergy, forces))
    assert result["status"] == "ok"
    bad = json.loads(mcp_gate(positions, coenergy[:-1], forces))
    assert bad["status"] == "invalid_input"


def test_force_coenergy_gate_requires_constant_current_semantics():
    positions, coenergy, forces = _quadratic_case()
    result = force_coenergy_displacement_gate(
        positions, coenergy, forces, energy_kind="stored_energy_at_fixed_flux"
    )
    assert result["status"] == "needs_attention"
    assert result["checks"]["constant_current_coenergy_recorded"] is False


@pytest.mark.parametrize(
    "case_id",
    [
        "v7_public_force_energy_derivative_sign_conflict",
        "v7_public_axisymmetric_two_pi_double_count",
    ],
)
def test_generalization_v7_public(case_id):
    positions, coenergy, forces = _quadratic_case()
    if case_id == "v7_public_force_energy_derivative_sign_conflict":
        forces = [-force for force in forces]
    else:
        forces = [2.0 * math.pi * force for force in forces]
    result = force_coenergy_displacement_gate(positions, coenergy, forces)
    assert result["status"] == "needs_attention"
    assert result["checks"]["central_virtual_work_matches_direct_force"] is False


@pytest.mark.parametrize(
    ("case_id", "failed_check"),
    [
        (
            "v8_public_force_snapshot_time_skew",
            "force_and_coenergy_share_load_step_snapshot",
        ),
        (
            "v8_public_coenergy_derivative_mesh_generation_mix",
            "coenergy_stencil_uses_one_mesh_family_generation",
        ),
    ],
)
def test_generalization_v8_public(case_id, failed_check):
    positions, coenergy, forces = _quadratic_case()
    identity = _artifact_identity(len(positions))
    if case_id == "v8_public_force_snapshot_time_skew":
        identity["direct_force_snapshot"].update(
            {"load_step_id": "load-step-43", "time_s": 0.026}
        )
    else:
        identity["coenergy_mesh_family_generations"][4] = "mesh-family-8"
    result = force_coenergy_displacement_gate(
        positions,
        coenergy,
        forces,
        artifact_identity=identity,
    )
    assert result["status"] == "needs_attention"
    assert result["checks"][failed_check] is False


@pytest.mark.parametrize(
    ("case_id", "failed_check"),
    [
        (
            "v9_public_virtual_work_displacement_unit_mismatch",
            "displacement_axis_uses_one_si_unit",
        ),
        (
            "v9_public_force_direction_frame_reflected",
            "force_vectors_share_transformed_frame",
        ),
    ],
)
def test_generalization_v9_public(case_id, failed_check):
    positions, coenergy, forces = _quadratic_case()
    identity = _artifact_identity(len(positions))
    if case_id == "v9_public_virtual_work_displacement_unit_mismatch":
        identity["displacement_axis"].update(
            {"numeric_unit": "mm", "derivative_unit": "m", "scale_to_si": 1.0}
        )
    else:
        identity["force_frame"].update(
            {
                "direct_frame_id": "reflected-frame",
                "direct_axis": [0.0, 1.0, 0.0],
                "reflection_applied": False,
            }
        )
    result = force_coenergy_displacement_gate(
        positions,
        coenergy,
        forces,
        artifact_identity=identity,
    )
    assert result["status"] == "needs_attention"
    assert result["checks"][failed_check] is False


@pytest.mark.parametrize(
    ("case_id", "failed_check"),
    [
        (
            "v10_public_axisymmetric_force_revolution_factor_twice",
            "axisymmetric_force_is_already_total_3d",
        ),
        (
            "v10_public_force_selection_includes_air_component",
            "weighted_stress_selects_only_target_magnetic_body",
        ),
    ],
)
def test_generalization_v10_public(case_id, failed_check):
    positions, coenergy, forces = _quadratic_case()
    identity = _artifact_identity(len(positions))
    if case_id == "v10_public_axisymmetric_force_revolution_factor_twice":
        identity["force_normalization"]["revolution_factor_application_count"] = 1
    else:
        identity["force_body_selection"][
            "weighted_stress_selected_group_ids"
        ] = [0, 1]
    result = force_coenergy_displacement_gate(
        positions, coenergy, forces, artifact_identity=identity
    )
    assert result["status"] == "needs_attention"
    assert result["checks"][failed_check] is False
