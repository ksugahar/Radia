import json

from radia_mcp.ih.server import ih_team36_contract, ih_team36_validate
from radia_mcp.ih.team36_gate import (
    EXPECTED_EXCITATION,
    EXPECTED_GEOMETRY_M,
    SCHEMA,
    SOURCE_URL,
    evaluate_team36_artifact,
)


EM_HASH = "a" * 64
THERMAL_HASH = "b" * 64
IDENTITY = {
    "geometry_sha256": "c" * 64,
    "material_tables_sha256": "d" * 64,
    "excitation_sha256": "e" * 64,
    "coordinate_system": "axisymmetric_r_z",
}


def artifact():
    return {
        "radia_version": "4.0.0",
        "executed_at_utc": "2026-07-24T00:00:00Z",
        "host": "validation-worker",
        "artifact_schema": SCHEMA,
        "benchmark_source": SOURCE_URL,
        "coordinate_system": "axisymmetric_r_z",
        "identity": dict(IDENTITY),
        "geometry_m": dict(EXPECTED_GEOMETRY_M),
        "excitation": dict(EXPECTED_EXCITATION),
        "material_model": {
            "resistivity_point_count": 16,
            "mu20_point_count": 20,
            "conductivity_point_count": 17,
            "heat_capacity_point_count": 36,
            "curie_temperature_c": 770.0,
            "transition_width_c": 20.0,
        },
        "meshes": {
            "electromagnetic": {
                "topology_sha256": EM_HASH,
                "vertex_count": 2503,
                "element_count": 4916,
                "element_kinds": ["triangle"],
                "element_order": 1,
                "billet_skin_layer_count": 5,
            },
            "thermal": {
                "topology_sha256": THERMAL_HASH,
                "vertex_count": 981,
                "element_count": 1654,
                "element_kinds": ["triangle"],
                "element_order": 1,
            },
        },
        "coupling": {
            "temperature_to_em": {
                "source_mesh_sha256": THERMAL_HASH,
                "target_mesh_sha256": EM_HASH,
                "sample_count": 800,
                "outside_count": 0,
            },
            "joule_power_to_thermal": {
                "source_mesh_sha256": EM_HASH,
                "target_mesh_sha256": THERMAL_HASH,
                "sample_count": 1000,
                "maximum_relative_error": 1.0e-12,
                "maximum_scale_deviation": 0.04,
            },
        },
        "history": [
            {
                "time_s": 0.0,
                "axis_temperature_c": 20.0,
                "surface_temperature_c": 20.0,
                "maximum_temperature_c": 20.0,
                "induced_power_w": 0.0,
            },
            {
                "time_s": 250.0,
                "axis_temperature_c": 835.0,
                "surface_temperature_c": 850.0,
                "maximum_temperature_c": 1200.0,
                "induced_power_w": 26700.0,
                "em_converged": True,
                "thermal_converged": True,
            },
        ],
        "provenance": {
            "executed_at_utc": "2026-07-24T00:00:00Z",
            "host": "validation-worker",
            "radia_version": "4.0.0",
            "ngsolve_version": "6.2.0",
            "git_commit": "f" * 64,
        },
        "timing_s": {
            "mesh": 1.0,
            "electromagnetic": 10.0,
            "mapping": 2.0,
            "thermal": 8.0,
        },
    }


def reference():
    return {
        "identity": dict(IDENTITY),
        "comparisons": [
            {
                "observable": "surface_temperature_c_at_250_s",
                "radia_value": 850.0,
                "reference_value": 854.0,
                "relative_tolerance": 0.01,
            }
        ],
    }


def test_complete_run_is_solver_ready_but_not_cross_validated_without_reference():
    result = evaluate_team36_artifact(artifact())
    assert result["accepted_for_solver_execution"]
    assert not result["accepted_for_cross_validation"]
    assert not result["accepted_for_mcp_learning"]
    assert result["failed_checks"] == [
        "cross_reference_supplied",
        "cross_reference_identity_matches",
        "cross_reference_observables_match",
    ]


def test_identity_matched_reference_promotes_cross_validation():
    result = evaluate_team36_artifact(artifact(), reference=reference())
    assert result["accepted_for_solver_execution"]
    assert result["accepted_for_cross_validation"]
    assert result["accepted_for_mcp_learning"]
    assert result["failed_checks"] == []


def test_reference_values_are_bound_to_the_artifact_history():
    value = artifact()
    value["history"][-1]["surface_temperature_c"] = 700.0
    result = evaluate_team36_artifact(value, reference=reference())
    assert result["accepted_for_solver_execution"]
    assert not result["accepted_for_cross_validation"]
    assert not result["checks"]["cross_reference_observables_match"]


def test_unknown_reference_observable_is_rejected():
    value = reference()
    value["comparisons"][0]["observable"] = "unregistered_value"
    result = evaluate_team36_artifact(artifact(), reference=value)
    assert not result["accepted_for_cross_validation"]
    assert not result["checks"]["cross_reference_observables_match"]


def test_coincident_mesh_and_nonconservative_map_are_rejected():
    value = artifact()
    value["meshes"]["thermal"]["topology_sha256"] = EM_HASH
    value["meshes"]["thermal"]["vertex_count"] = 2503
    value["meshes"]["thermal"]["element_count"] = 4916
    value["coupling"]["temperature_to_em"]["source_mesh_sha256"] = EM_HASH
    value["coupling"]["joule_power_to_thermal"]["target_mesh_sha256"] = EM_HASH
    value["coupling"]["joule_power_to_thermal"]["maximum_relative_error"] = 0.03
    result = evaluate_team36_artifact(value, reference=reference())
    assert not result["accepted_for_solver_execution"]
    assert not result["checks"]["meshes_are_noncoincident"]
    assert not result["checks"]["joule_mapping_is_conservative"]


def test_mcp_boundary_accepts_json_not_paths():
    contract = ih_team36_contract()
    assert contract["artifact_schema"] == SCHEMA
    result = ih_team36_validate(json.dumps(artifact()), json.dumps(reference()))
    assert result["accepted_for_cross_validation"]
