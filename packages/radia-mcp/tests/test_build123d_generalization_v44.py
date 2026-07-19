from __future__ import annotations

from test_build123d_generalization_v29 import _public_result, _source_result
from test_build123d_generalization_v43 import _public_v43, _source_v43


_BOOLEAN = "boolean_shell_fillet_massproperties_centerofmass_volume_brep_owner_generation_identity"
_LOFT = "loft_sweep_section_orientation_tangent_area_volume_inertia_brep_generation_identity"
_SKETCH = "sketch_constraint_solver_order_plane_frame_parameter_cache_shape_owner_generation_identity"
_STEP = "step_export_units_tessellation_tolerance_facecount_brep_digest_owner_generation_identity"

_CASE_IDS = (
    "v44_public_boolean_shell_fillet_massproperties_centerofmass_volume_brep_owner_mismatch",
    "v44_public_loft_sweep_section_orientation_tangent_area_volume_inertia_brep_mismatch",
    "v44_source_sketch_constraint_solver_order_plane_frame_parameter_cache_shape_owner_mismatch",
    "v44_source_step_export_units_tessellation_tolerance_facecount_brep_digest_owner_mismatch",
)


def _public_v44():
    reference, measured = _public_v43()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            digest = str(index + 1) * 64
            generation = "boolean-shell-fillet-v44-731"
            row[_BOOLEAN] = {
                "boolean_generation": generation,
                **{name: generation for name in ("history_generation", "shell_generation", "fillet_generation", "mass_generation", "center_generation", "volume_generation", "owner_generation", "brep_generation", "result_generation")},
                "operation": "cut", "result_operation": "cut",
                "shell_thickness_m": 0.001, "result_shell_thickness_m": 0.001,
                "fillet_radius_m": 0.0005, "result_fillet_radius_m": 0.0005,
                "center_of_mass_m": [0.01, 0.02, 0.03], "result_center_of_mass_m": [0.01, 0.02, 0.03],
                "volume_m3": 9.9e-4, "result_volume_m3": 9.9e-4,
                "surface_area_m2": 6.1e-2, "result_surface_area_m2": 6.1e-2,
                "topology_signature": {"solid": 1, "shell": 1, "face": 10}, "result_topology_signature": {"solid": 1, "shell": 1, "face": 10},
                "shape_owner": "part:boolean-shell-fillet-731", "result_shape_owner": "part:boolean-shell-fillet-731",
                "boolean_brep_sha256": digest, "accepted_boolean_brep_sha256": digest,
            }
            generation = "loft-sweep-v44-731"
            row[_LOFT] = {
                "loft_generation": generation,
                **{name: generation for name in ("section_generation", "orientation_generation", "tangent_generation", "area_generation", "volume_generation", "inertia_generation", "owner_generation", "brep_generation", "result_generation")},
                "section_count": 3, "result_section_count": 3,
                "section_orientation": "consistent_ccw", "result_section_orientation": "consistent_ccw",
                "tangent_continuity": "C1", "result_tangent_continuity": "C1",
                "section_area_m2": [1e-4, 1.2e-4, 1e-4], "result_section_area_m2": [1e-4, 1.2e-4, 1e-4],
                "volume_m3": 2.5e-5, "result_volume_m3": 2.5e-5,
                "inertia_tensor_kg_m2": [[1e-8, 0.0, 0.0], [0.0, 2e-8, 0.0], [0.0, 0.0, 3e-8]],
                "result_inertia_tensor_kg_m2": [[1e-8, 0.0, 0.0], [0.0, 2e-8, 0.0], [0.0, 0.0, 3e-8]],
                "shape_owner": "part:loft-sweep-731", "result_shape_owner": "part:loft-sweep-731",
                "loft_brep_sha256": digest, "accepted_loft_brep_sha256": digest,
            }
    return reference, measured


def _source_v44():
    row = _source_v43()
    identity = row["replay_identity"]
    generation = "sketch-replay-v44-731"
    identity[_SKETCH] = {
        "sketch_generation": generation,
        **{name: generation for name in ("constraint_generation", "solver_order_generation", "plane_frame_generation", "parameter_cache_generation", "shape_generation", "owner_generation", "result_generation")},
        "constraint_order": ["coincident:1", "distance:2", "horizontal:3"], "replayed_constraint_order": ["coincident:1", "distance:2", "horizontal:3"],
        "plane_frame": "XY", "replayed_plane_frame": "XY", "parameter_cache_key": "width=0.1;radius=0.02", "replayed_parameter_cache_key": "width=0.1;radius=0.02",
        "solver_status": "solved", "replayed_solver_status": "solved", "shape_generation_id": 731, "replayed_shape_generation_id": 731,
        "shape_owner": "headless:sketch-731", "replayed_shape_owner": "headless:sketch-731", "sketch_result_sha256": "a" * 64, "accepted_sketch_result_sha256": "a" * 64,
    }
    generation = "step-export-v44-731"
    identity[_STEP] = {
        "step_generation": generation,
        **{name: generation for name in ("unit_generation", "tessellation_generation", "tolerance_generation", "facecount_generation", "brep_generation", "digest_generation", "owner_generation", "result_generation")},
        "length_unit": "mm", "replayed_length_unit": "mm", "tessellation_tolerance_m": 1e-5, "replayed_tessellation_tolerance_m": 1e-5,
        "face_count": 42, "replayed_face_count": 42, "topology_signature": {"solid": 1, "shell": 1, "face": 42}, "replayed_topology_signature": {"solid": 1, "shell": 1, "face": 42},
        "brep_generation_id": 731, "replayed_brep_generation_id": 731, "export_owner": "headless:step-export-731", "replayed_export_owner": "headless:step-export-731",
        "step_digest_sha256": "b" * 64, "replayed_step_digest_sha256": "b" * 64, "accepted_step_result_sha256": "c" * 64,
    }
    return row


def test_v44_positive_public_and_source_contracts():
    reference, measured = _public_v44()
    assert _public_result(reference, measured)["status"] == "ok"
    assert _source_result(_source_v44())["status"] == "ok"
    assert len(_CASE_IDS) == 4


def test_v44_rejects_boolean_shell_fillet_mismatch():
    reference, measured = _public_v44()
    measured["external_cad"][0][_BOOLEAN]["result_volume_m3"] = -1.0
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v44_rejects_loft_orientation_mismatch():
    reference, measured = _public_v44()
    measured["external_cad"][0][_LOFT]["result_section_orientation"] = "inconsistent"
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v44_rejects_sketch_replay_mismatch():
    row = _source_v44()
    row["replay_identity"][_SKETCH]["replayed_constraint_order"] = ["distance:2", "coincident:1", "horizontal:3"]
    assert _source_result(row)["status"] == "needs_attention"


def test_v44_rejects_step_export_mismatch():
    row = _source_v44()
    row["replay_identity"][_STEP]["replayed_length_unit"] = "m"
    assert _source_result(row)["status"] == "needs_attention"


def test_v44_rejects_malformed_numeric_identity_without_raising():
    reference, measured = _public_v44()
    measured["external_cad"][0][_BOOLEAN]["volume_m3"] = {"bad": "value"}
    assert _public_result(reference, measured)["status"] == "needs_attention"

    row = _source_v44()
    row["replay_identity"][_STEP]["face_count"] = {"bad": "value"}
    assert _source_result(row)["status"] == "needs_attention"


def test_v44_rejects_self_consistent_invalid_mass_properties_and_owner():
    reference, measured = _public_v44()
    identity = measured["external_cad"][0][_BOOLEAN]
    identity["center_of_mass_m"] = [float("nan"), 0.0, 0.0]
    identity["result_center_of_mass_m"] = identity["center_of_mass_m"]
    identity["shape_owner"] = identity["result_shape_owner"] = ""
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v44_rejects_self_consistent_invalid_loft_geometry():
    reference, measured = _public_v44()
    identity = measured["external_cad"][0][_LOFT]
    identity["section_area_m2"] = [1.0e-4, -1.0, 1.0e-4]
    identity["result_section_area_m2"] = identity["section_area_m2"]
    identity["inertia_tensor_kg_m2"] = [[-1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]]
    identity["result_inertia_tensor_kg_m2"] = identity["inertia_tensor_kg_m2"]
    assert _public_result(reference, measured)["status"] == "needs_attention"
