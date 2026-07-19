from __future__ import annotations

from test_build123d_generalization_v44 import _public_v44, _source_v44
from test_build123d_generalization_v29 import _public_result, _source_result


_PROMOTED_CASE_IDS = (
    "v45_public_boolean_fillet_shell_massproperties_volume_area_brep_owner_mismatch",
    "v45_public_loft_sweep_section_frame_tangent_continuity_inertia_export_digest_mismatch",
    "v45_source_sketch_constraint_order_plane_frame_solver_cache_shape_generation_owner_mismatch",
    "v45_source_step_units_tessellation_tolerance_face_topology_brep_export_owner_mismatch",
)


def _public_v45():
    reference, measured = _public_v44()
    for rows in (reference, *measured.values()):
        for row in rows:
            row["boolean_fillet_shell_massproperties_volume_area_brep_owner_identity"] = {
                "generation": "boolean-shell-fillet-v45-812", "operation": "cut", "result_operation": "cut",
                "shell_thickness_m": 0.001, "result_shell_thickness_m": 0.001, "fillet_radius_m": 0.0005, "result_fillet_radius_m": 0.0005,
                "center_of_mass_m": [0.01, 0.02, 0.03], "result_center_of_mass_m": [0.01, 0.02, 0.03], "volume_m3": 9.9e-4, "result_volume_m3": 9.9e-4, "surface_area_m2": 6.1e-2, "result_surface_area_m2": 6.1e-2,
                "topology_signature": {"solid": 1, "shell": 1, "face": 10}, "result_topology_signature": {"solid": 1, "shell": 1, "face": 10}, "shape_owner": "part:boolean-shell-fillet-v45-812", "result_shape_owner": "part:boolean-shell-fillet-v45-812", "boolean_brep_sha256": "1" * 64, "accepted_boolean_brep_sha256": "1" * 64,
            }
            row["loft_sweep_section_frame_tangent_continuity_inertia_export_digest_identity"] = {
                "generation": "loft-sweep-v45-812", "section_count": 3, "result_section_count": 3, "section_orientation": "consistent_ccw", "result_section_orientation": "consistent_ccw", "tangent_continuity": "C1", "result_tangent_continuity": "C1", "section_area_m2": [1e-4, 1.2e-4, 1e-4], "result_section_area_m2": [1e-4, 1.2e-4, 1e-4], "volume_m3": 2.5e-5, "result_volume_m3": 2.5e-5, "inertia_tensor_kg_m2": [[1e-8, 0.0, 0.0], [0.0, 2e-8, 0.0], [0.0, 0.0, 3e-8]], "result_inertia_tensor_kg_m2": [[1e-8, 0.0, 0.0], [0.0, 2e-8, 0.0], [0.0, 0.0, 3e-8]], "shape_owner": "part:loft-sweep-v45-812", "result_shape_owner": "part:loft-sweep-v45-812", "loft_brep_sha256": "2" * 64, "accepted_loft_brep_sha256": "2" * 64,
            }
    return reference, measured


def _source_v45():
    row = _source_v44()
    replay = row["replay_identity"]
    replay["sketch_constraint_order_plane_frame_solver_cache_shape_generation_owner_identity"] = {"generation": "sketch-v45-812", "constraint_order": ["coincident:1", "distance:2"], "replayed_constraint_order": ["coincident:1", "distance:2"], "plane_frame": "XY", "replayed_plane_frame": "XY", "parameter_cache_key": "width=0.1", "replayed_parameter_cache_key": "width=0.1", "solver_status": "solved", "replayed_solver_status": "solved", "shape_generation_id": 812, "replayed_shape_generation_id": 812, "shape_owner": "headless:sketch-v45-812", "replayed_shape_owner": "headless:sketch-v45-812", "result_sha256": "3" * 64, "accepted_result_sha256": "3" * 64}
    replay["step_units_tessellation_tolerance_face_topology_brep_export_owner_identity"] = {"generation": "step-v45-812", "length_unit": "mm", "replayed_length_unit": "mm", "tessellation_tolerance_m": 1e-5, "replayed_tessellation_tolerance_m": 1e-5, "face_topology": {"solid": 1, "shell": 1, "face": 42}, "replayed_face_topology": {"solid": 1, "shell": 1, "face": 42}, "brep_generation_id": 812, "replayed_brep_generation_id": 812, "export_owner": "headless:step-v45-812", "replayed_export_owner": "headless:step-v45-812", "step_digest_sha256": "4" * 64, "replayed_step_digest_sha256": "4" * 64, "accepted_result_sha256": "4" * 64}
    return row


def test_v45_positive_identity_contracts():
    reference, measured = _public_v45()
    assert _public_result(reference, measured)["status"] == "ok"
    assert _source_result(_source_v45())["status"] == "ok"
    assert len(_PROMOTED_CASE_IDS) == 4


def test_v45_rejects_boolean_brep_mismatch():
    reference, measured = _public_v45()
    measured["external_cad"][0]["boolean_fillet_shell_massproperties_volume_area_brep_owner_identity"]["result_volume_m3"] = -1.0
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v45_rejects_loft_frame_mismatch():
    reference, measured = _public_v45()
    measured["external_cad"][0]["loft_sweep_section_frame_tangent_continuity_inertia_export_digest_identity"]["result_section_orientation"] = "inconsistent"
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v45_rejects_sketch_cache_mismatch():
    row = _source_v45()
    row["replay_identity"]["sketch_constraint_order_plane_frame_solver_cache_shape_generation_owner_identity"]["replayed_parameter_cache_key"] = "old"
    assert _source_result(row)["status"] == "needs_attention"


def test_v45_rejects_step_topology_mismatch():
    row = _source_v45()
    row["replay_identity"]["step_units_tessellation_tolerance_face_topology_brep_export_owner_identity"]["replayed_length_unit"] = "m"
    assert _source_result(row)["status"] == "needs_attention"
