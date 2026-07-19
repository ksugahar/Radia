from copy import deepcopy

from radia_mcp.build123d.feature_replay_identity_v50 import (
    FEATURE,
    HEALING,
    MATE,
    STL,
    validate_public_identity,
    validate_source_identity,
)


PROMOTED_CASE_IDS = {
    "v50_public_assembly_mate_constraint_dof_frame_occurrence_transform_owner_mismatch",
    "v50_public_fillet_chamfer_edge_selector_radius_topology_history_owner_mismatch",
    "v50_source_tool_occt_tolerance_healing_sewing_shell_solid_orientation_owner_mismatch",
    "v50_source_tool_stl_tessellation_linear_angular_deflection_triangle_normal_unit_owner_mismatch",
}


def _payloads() -> tuple[dict[str, object], dict[str, object]]:
    generation = "build123d-v50"
    result = {"result_sha256": "f" * 64, "accepted_result_sha256": "f" * 64}
    frames = {
        "base": [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        "shaft": [0.0, 0.0, 0.02, 1.0, 0.0, 0.0, 0.0],
    }
    transforms = {
        "occurrence:base": frames["base"],
        "occurrence:shaft": frames["shaft"],
    }
    selectors = {"fillet": ["edge:11", "edge:12"], "chamfer": ["edge:21"]}
    history = {
        "edge:11": ["edge:31", "face:41"],
        "edge:12": ["edge:32", "face:42"],
        "edge:21": ["edge:33", "face:43"],
    }
    row = {
        MATE: {
            "generation": generation,
            "constraint_generation": generation,
            "dof_generation": generation,
            "frame_generation": generation,
            "occurrence_generation": generation,
            "transform_generation": generation,
            "owner_generation": generation,
            "result_generation": generation,
            "mate_constraints": ["fixed:base", "revolute:shaft-base"],
            "result_mate_constraints": ["fixed:base", "revolute:shaft-base"],
            "remaining_dof": 1,
            "result_remaining_dof": 1,
            "mate_frames": frames,
            "result_mate_frames": frames,
            "occurrence_transforms": transforms,
            "result_occurrence_transforms": transforms,
            "assembly_owner": "assembly:mates-v50",
            "result_assembly_owner": "assembly:mates-v50",
            **result,
        },
        FEATURE: {
            "generation": generation,
            "selector_generation": generation,
            "radius_generation": generation,
            "topology_generation": generation,
            "history_generation": generation,
            "owner_generation": generation,
            "result_generation": generation,
            "edge_selectors": selectors,
            "result_edge_selectors": selectors,
            "fillet_radius_m": 0.002,
            "result_fillet_radius_m": 0.002,
            "chamfer_distance_m": 0.001,
            "result_chamfer_distance_m": 0.001,
            "topology_history": history,
            "result_topology_history": history,
            "shape_owner": "shape:features-v50",
            "result_shape_owner": "shape:features-v50",
            **result,
        },
    }
    orientations = {"shell:outer": 1, "shell:inner": -1}
    healing = {
        "generation": generation,
        "tolerance_generation": generation,
        "healing_generation": generation,
        "sewing_generation": generation,
        "shell_generation": generation,
        "solid_generation": generation,
        "orientation_generation": generation,
        "owner_generation": generation,
        "result_generation": generation,
        "input_tolerance_m": 1e-7,
        "result_input_tolerance_m": 1e-7,
        "healing_applied": True,
        "result_healing_applied": True,
        "sewing_tolerance_m": 5e-7,
        "result_sewing_tolerance_m": 5e-7,
        "shell_count": 2,
        "result_shell_count": 2,
        "solid_count": 1,
        "result_solid_count": 1,
        "shell_orientation_signs": orientations,
        "result_shell_orientation_signs": orientations,
        "shape_owner": "shape:healed-v50",
        "result_shape_owner": "shape:healed-v50",
        **result,
    }
    normals = {"outward": 12480, "inward": 0, "degenerate": 0}
    stl = {
        "generation": generation,
        "linear_generation": generation,
        "angular_generation": generation,
        "triangle_generation": generation,
        "normal_generation": generation,
        "unit_generation": generation,
        "owner_generation": generation,
        "result_generation": generation,
        "linear_deflection_m": 1e-4,
        "result_linear_deflection_m": 1e-4,
        "angular_deflection_rad": 0.1,
        "result_angular_deflection_rad": 0.1,
        "triangle_count": 12480,
        "result_triangle_count": 12480,
        "normal_counts": normals,
        "result_normal_counts": normals,
        "length_unit": "m",
        "result_length_unit": "m",
        "mesh_owner": "mesh:stl-v50",
        "result_mesh_owner": "mesh:stl-v50",
        **result,
    }
    public = {"reference": [row], "measured": {"cad": [deepcopy(row)]}}
    source = {"replay_identity": {HEALING: healing, STL: stl}}
    return public, source


def test_v50_positive_public_and_source_replays_are_accepted() -> None:
    public, source = _payloads()
    assert validate_public_identity(public)["status"] == "ok"
    assert validate_source_identity(source)["status"] == "ok"


def test_v50_public_mutations_are_rejected() -> None:
    public, _ = _payloads()
    public["reference"][0][MATE]["result_remaining_dof"] = 5
    public["reference"][0][FEATURE]["result_fillet_radius_m"] = 0.02
    result = validate_public_identity(public)
    assert result["status"] == "needs_attention"
    assert len(result["issues"]) == 2


def test_v50_source_mutations_are_rejected() -> None:
    _, source = _payloads()
    source["replay_identity"][HEALING]["result_solid_count"] = 0
    source["replay_identity"][STL]["result_length_unit"] = "mm"
    result = validate_source_identity(source)
    assert result["status"] == "needs_attention"
    assert len(result["issues"]) == 2


def test_v50_self_consistent_invalid_frames_and_feature_history_are_rejected() -> None:
    public, _ = _payloads()
    for row in [public["reference"][0], public["measured"]["cad"][0]]:
        row[MATE]["mate_frames"]["shaft"] = row[MATE]["result_mate_frames"]["shaft"] = [0, 0, 0.02, 2, 0, 0, 0]
        row[FEATURE]["topology_history"] = row[FEATURE]["result_topology_history"] = {"edge:11": ["edge:31"]}
    assert validate_public_identity(public)["status"] == "needs_attention"


def test_v50_self_consistent_open_shell_and_bad_stl_normals_are_rejected() -> None:
    _, source = _payloads()
    healing = source["replay_identity"][HEALING]
    healing["solid_count"] = healing["result_solid_count"] = 0
    stl = source["replay_identity"][STL]
    stl["normal_counts"] = stl["result_normal_counts"] = {"outward": 12000, "inward": 480, "degenerate": 0}
    assert validate_source_identity(source)["status"] == "needs_attention"
