from copy import deepcopy

from radia_mcp.radia_ngsolve.solver_state_identity_v49 import validate_public_v49_identity


PROMOTED_CASE_IDS = {
    "v49_public_nonlinear_material_interpolation_branch_unit_temperature_extrapolation_owner_mismatch",
    "v49_public_moving_mesh_sliding_interface_frame_time_remesh_solution_owner_mismatch",
}


def _records() -> dict[str, object]:
    material_generation = "nonlinear-material-v49"
    sliding_generation = "sliding-interface-v49"
    return {
        "nonlinear_material_interpolation_branch_unit_temperature_extrapolation_owner_identity": {
            "generation": material_generation,
            "material_generation": material_generation,
            "branch_generation": material_generation,
            "temperature_generation": material_generation,
            "interpolation_generation": material_generation,
            "result_generation": material_generation,
            "interpolation_branch": "ascending-major-loop",
            "result_interpolation_branch": "ascending-major-loop",
            "input_units": {"magnetic_flux_density": "T", "magnetic_field": "A/m", "temperature": "K"},
            "result_input_units": {"magnetic_flux_density": "T", "magnetic_field": "A/m", "temperature": "K"},
            "temperature_value": 353.15,
            "result_temperature_value": 353.15,
            "extrapolation_policy": "reject",
            "result_extrapolation_policy": "reject",
            "material_owner": "material:nonlinear-steel-v49",
            "result_material_owner": "material:nonlinear-steel-v49",
            "result_sha256": "1" * 64,
            "accepted_result_sha256": "1" * 64,
        },
        "moving_mesh_sliding_interface_frame_time_remesh_solution_owner_identity": {
            "generation": sliding_generation,
            "mesh_generation": sliding_generation,
            "interface_generation": sliding_generation,
            "frame_generation": sliding_generation,
            "time_generation": sliding_generation,
            "remesh_generation": sliding_generation,
            "solution_generation": sliding_generation,
            "result_generation": sliding_generation,
            "sliding_interface_map_sha256": "2" * 64,
            "result_sliding_interface_map_sha256": "2" * 64,
            "coordinate_frame": "spatial",
            "result_coordinate_frame": "spatial",
            "time_value_s": 0.0125,
            "result_time_value_s": 0.0125,
            "remesh_revision": "remesh-v49-r3",
            "result_remesh_revision": "remesh-v49-r3",
            "solution_owner": "solution:moving-mesh-v49",
            "result_solution_owner": "solution:moving-mesh-v49",
            "result_sha256": "3" * 64,
            "accepted_result_sha256": "3" * 64,
        },
    }


def test_v49_public_positive_replay_is_accepted() -> None:
    result = validate_public_v49_identity(_records())
    assert result["status"] == "ok"
    assert all(result["checks"].values())


def test_v49_public_mixed_material_state_is_rejected() -> None:
    records = deepcopy(_records())
    row = records["nonlinear_material_interpolation_branch_unit_temperature_extrapolation_owner_identity"]
    row["result_interpolation_branch"] = "descending-recoil-loop"
    row["result_input_units"] = {"temperature": "degC"}
    row["result_temperature_value"] = 80.0
    row["result_extrapolation_policy"] = "linear"
    row["result_material_owner"] = "material:old"
    assert validate_public_v49_identity(records)["status"] == "needs_attention"


def test_v49_public_mixed_sliding_mesh_state_is_rejected() -> None:
    records = deepcopy(_records())
    row = records["moving_mesh_sliding_interface_frame_time_remesh_solution_owner_identity"]
    row["result_sliding_interface_map_sha256"] = "a" * 64
    row["result_coordinate_frame"] = "material"
    row["result_time_value_s"] = 0.01
    row["result_remesh_revision"] = "remesh-v49-r2"
    row["result_solution_owner"] = "solution:old"
    assert validate_public_v49_identity(records)["status"] == "needs_attention"

