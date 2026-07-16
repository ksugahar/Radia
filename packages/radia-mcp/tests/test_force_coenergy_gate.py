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
        "virtual_work_constraint_basis": {
            "direct_force_constraint": "fixed_current",
            "coenergy_derivative_constraint": "fixed_current",
            "current_control_generation": "current-control-13",
            "derivative_control_generation": "current-control-13",
            "flux_constraint_active": False,
        },
        "eddy_loss_harmonic_basis": {
            "harmonic_frequency_hz": [50.0, 150.0, 250.0],
            "skin_depth_state_frequency_hz": [50.0, 150.0, 250.0],
            "amplitude_basis_id": "harmonic-basis-13",
            "material_state_basis_id": "harmonic-basis-13",
            "solve_generation": "solve-13",
            "material_state_solve_generation": "solve-13",
        },
        "axisymmetric_force_measure_identity": {
            "formulation": "axisymmetric",
            "integration_measure": "2*pi*r*dr*dz",
            "reference_integration_measure": "2*pi*r*dr*dz",
            "radius_weighting_basis_id": "axisymmetric-radius-weighted-v1",
            "force_result_basis_id": "axisymmetric-radius-weighted-v1",
            "radius_coordinate_frame": "cylindrical-rz",
            "force_component_frame": "global-z",
            "solve_generation": "solve-14",
            "integration_solve_generation": "solve-14",
        },
        "eddy_loss_material_frequency_identity": {
            "field_solution_frequency_hz": 1000.0,
            "loss_evaluation_frequency_hz": 1000.0,
            "material_state_frequency_hz": 1000.0,
            "material_assignment_generation": "materials-14",
            "conductivity_material_generation": "materials-14",
            "lamination_material_generation": "materials-14",
            "solve_generation": "solve-14",
            "material_state_solve_generation": "solve-14",
            "conductivity_sha256": "1" * 64,
            "lamination_state_sha256": "2" * 64,
        },
        "weighted_stress_mask_mesh_identity": {
            "active_air_mesh_generation": "air-mesh-15",
            "field_solution_mesh_generation": "air-mesh-15",
            "weighted_mask_mesh_generation": "air-mesh-15",
            "force_integration_mesh_generation": "air-mesh-15",
            "weighted_mask_sha256": "5" * 64,
            "force_mask_sha256": "5" * 64,
            "mask_basis": "nodal_weighting_function",
            "force_method": "weighted_stress_tensor",
        },
        "complex_current_phasor_basis_identity": {
            "source_current_basis": "rms_phasor",
            "field_current_basis": "rms_phasor",
            "force_loss_current_basis": "rms_phasor",
            "source_scale_to_rms": 1.0,
            "field_scale_to_rms": 1.0,
            "force_loss_scale_to_rms": 1.0,
            "complex_time_convention": "exp(+jwt)",
            "result_time_convention": "exp(+jwt)",
            "solve_generation": "solve-15",
            "result_generation": "solve-15",
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


@pytest.mark.parametrize(
    ("case_id", "failed_check"),
    [
        (
            "v11_public_force_coenergy_constraint_basis_mismatch",
            "force_and_coenergy_use_same_fixed_current_constraint",
        ),
        (
            "v11_public_eddy_loss_frequency_harmonic_basis_mismatch",
            "eddy_loss_harmonics_share_frequency_and_material_basis",
        ),
    ],
)
def test_generalization_v11_public(case_id, failed_check):
    positions, coenergy, forces = _quadratic_case()
    identity = _artifact_identity(len(positions))
    if case_id == "v11_public_force_coenergy_constraint_basis_mismatch":
        identity["virtual_work_constraint_basis"].update(
            {
                "coenergy_derivative_constraint": "fixed_flux",
                "derivative_control_generation": "flux-control-13",
                "flux_constraint_active": True,
            }
        )
    else:
        identity["eddy_loss_harmonic_basis"].update(
            {
                "skin_depth_state_frequency_hz": [60.0, 180.0, 300.0],
                "material_state_basis_id": "harmonic-basis-12",
                "material_state_solve_generation": "solve-12",
            }
        )
    result = force_coenergy_displacement_gate(
        positions, coenergy, forces, artifact_identity=identity
    )
    assert result["status"] == "needs_attention"
    assert result["checks"][failed_check] is False


def test_v12_public_axisymmetric_force_radius_weighting_basis_mismatch() -> None:
    positions, coenergy, forces = _quadratic_case()
    identity = _artifact_identity(len(positions))
    identity["axisymmetric_force_measure_identity"].update(
        {
            "integration_measure": "dx*dy",
            "force_result_basis_id": "planar-per-depth-v1",
        }
    )
    result = force_coenergy_displacement_gate(
        positions, coenergy, forces, artifact_identity=identity
    )
    assert result["status"] == "needs_attention"
    assert result["checks"]["axisymmetric_force_uses_radius_weighted_measure"] is False


def test_v12_public_eddy_loss_frequency_material_generation_mismatch() -> None:
    positions, coenergy, forces = _quadratic_case()
    identity = _artifact_identity(len(positions))
    identity["eddy_loss_material_frequency_identity"].update(
        {
            "conductivity_material_generation": "materials-13",
            "lamination_material_generation": "materials-13",
            "material_state_solve_generation": "solve-13",
        }
    )
    result = force_coenergy_displacement_gate(
        positions, coenergy, forces, artifact_identity=identity
    )
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "eddy_loss_uses_current_frequency_and_material_generation"
        ]
        is False
    )


def test_v13_public_weighted_stress_mask_mesh_generation_mismatch() -> None:
    positions, coenergy, forces = _quadratic_case()
    identity = _artifact_identity(len(positions))
    identity["weighted_stress_mask_mesh_identity"].update(
        {
            "weighted_mask_mesh_generation": "air-mesh-14",
            "weighted_mask_sha256": "7" * 64,
        }
    )
    result = force_coenergy_displacement_gate(
        positions, coenergy, forces, artifact_identity=identity
    )
    assert result["status"] == "needs_attention"
    assert (
        result["checks"]["weighted_stress_mask_matches_current_air_mesh_generation"]
        is False
    )


def test_v13_public_complex_current_peak_rms_phasor_basis_mismatch() -> None:
    positions, coenergy, forces = _quadratic_case()
    identity = _artifact_identity(len(positions))
    identity["complex_current_phasor_basis_identity"].update(
        {
            "field_current_basis": "peak_phasor",
            "field_scale_to_rms": 1.0 / math.sqrt(2.0),
        }
    )
    result = force_coenergy_displacement_gate(
        positions, coenergy, forces, artifact_identity=identity
    )
    assert result["status"] == "needs_attention"
    assert (
        result["checks"]["complex_current_force_and_loss_share_phasor_basis"]
        is False
    )
