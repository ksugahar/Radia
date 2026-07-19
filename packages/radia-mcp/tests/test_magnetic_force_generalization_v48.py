from __future__ import annotations

from copy import deepcopy
import math

from radia_mcp.radia_ngsolve.bem_hysteresis_identity_v48 import BEM, HYSTERESIS, validate_public_identity


PROMOTED_CASE_IDS = {
    "v48_public_bem_near_far_panel_quadrature_normal_solid_angle_mesh_revision_mismatch",
    "v48_public_hysteresis_minor_loop_return_point_state_temperature_frequency_owner_mismatch",
}


def _identity() -> dict[str, object]:
    bem_generation = "bem-panel-v48-901"
    hysteresis_generation = "minor-loop-v48-901"
    panels = ["panel:1", "panel:2", "panel:3", "panel:4"]
    normals = [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, -1.0, 0.0]]
    angles = [math.pi, math.pi, math.pi, math.pi]
    points = [[0.2, 120.0], [0.8, 480.0], [0.3, 160.0]]
    state = {"branch": "descending", "last_reversal": 0.8, "memory_depth": 2}
    return {
        BEM: {
            "generation": bem_generation,
            "quadrature_generation": bem_generation,
            "normal_generation": bem_generation,
            "solid_angle_generation": bem_generation,
            "mesh_generation": bem_generation,
            "result_generation": bem_generation,
            "panel_ids": panels,
            "result_panel_ids": panels,
            "near_quadrature_order": [8, 8, 8, 8],
            "result_near_quadrature_order": [8, 8, 8, 8],
            "far_quadrature_order": [4, 4, 4, 4],
            "result_far_quadrature_order": [4, 4, 4, 4],
            "panel_normals": normals,
            "result_panel_normals": normals,
            "panel_solid_angles_sr": angles,
            "result_panel_solid_angles_sr": angles,
            "mesh_revision": "mesh:bem-v48-901",
            "result_mesh_revision": "mesh:bem-v48-901",
            "result_sha256": "6" * 64,
            "accepted_result_sha256": "6" * 64,
        },
        HYSTERESIS: {
            "generation": hysteresis_generation,
            "return_point_generation": hysteresis_generation,
            "state_generation": hysteresis_generation,
            "environment_generation": hysteresis_generation,
            "material_generation": hysteresis_generation,
            "result_generation": hysteresis_generation,
            "return_points": points,
            "result_return_points": points,
            "internal_state": state,
            "result_internal_state": state,
            "temperature_k": 353.15,
            "result_temperature_k": 353.15,
            "frequency_hz": 50.0,
            "result_frequency_hz": 50.0,
            "material_owner": "material:hysteresis-v48-901",
            "result_material_owner": "material:hysteresis-v48-901",
            "result_sha256": "7" * 64,
            "accepted_result_sha256": "7" * 64,
        },
    }


def test_v48_positive_bem_and_hysteresis_artifacts_are_accepted() -> None:
    assert all(validate_public_identity(_identity()).values())


def test_v48_bem_mesh_and_quadrature_mutations_are_rejected() -> None:
    identity = deepcopy(_identity())
    identity[BEM]["result_near_quadrature_order"] = [8, 4, 8, 8]
    identity[BEM]["result_mesh_revision"] = "mesh:old"
    assert validate_public_identity(identity)["bem_v48_panel_quadrature_normal_solid_angle_mesh"] is False


def test_v48_hysteresis_history_and_owner_mutations_are_rejected() -> None:
    identity = deepcopy(_identity())
    identity[HYSTERESIS]["result_internal_state"] = {"branch": "ascending", "last_reversal": 0.3, "memory_depth": 1}
    identity[HYSTERESIS]["result_material_owner"] = "material:old"
    assert validate_public_identity(identity)["hysteresis_v48_return_state_environment_material_owner"] is False
