import copy
import json

import pytest

from radia_mcp.build123d.server import (
    build123d_loft_example_source_replay_gate,
    build123d_lofted_shell_handoff_gate,
)


def summary() -> dict:
    native = {
        "solid_count": 1,
        "shell_count": 1,
        "face_count": 4,
        "edge_count": 6,
        "vertex_count": 4,
        "volume_mm3": 1306.3405290340206,
        "surface_area_mm2": 5285.090189964201,
        "is_valid": True,
    }
    source_replay = {
        "metrics": native,
        "slice_count": 10,
        "top_bottom_face_count": 2,
        "shell_offset_mm": 0.5,
        "official_expected_volume_mm3": 1306.3405290344635,
        "official_assertion_delta_mm3": 4.43e-10,
        "official_assertion_tolerance_mm3": 0.013063405290344635,
    }
    external_replays = []
    for index in (1, 2):
        for mode in ("noheal", "heal"):
            external_replays.append(
                {
                    "index": index,
                    "mode": mode,
                    "command": {"returned": True, "exception": None},
                    "snapshot": {
                        "volume_count": 1,
                        "positive_volume_count": 1,
                        "volume_sum_mm3": 1310.4631569386177,
                        "surface_count": 4,
                        "positive_surface_count": 4,
                        "surface_area_sum_mm2": 5285.489161360994,
                        "curve_count": 6,
                        "vertex_count": 4,
                    },
                }
            )
    source_checks = {
        "source_replays_are_deterministic": True,
        "official_volume_assertion_reproduced": True,
        "eleven_profiles_and_two_openings_recorded": True,
        "native_is_valid_single_solid_shell": True,
        "step_topology_matches_native": True,
        "brep_topology_matches_native": True,
        "step_mass_properties_match_native": True,
        "brep_mass_properties_match_native": True,
    }
    external_checks = {
        "all_import_commands_succeeded": True,
        "both_modes_replay_deterministically": True,
        "all_imports_are_one_positive_volume": True,
        "all_imports_preserve_topology": True,
        "external_volume_drift_is_nonzero_but_bounded": True,
        "external_area_drift_is_bounded": True,
        "heal_and_noheal_mass_properties_match": True,
        "mesh_intentionally_not_attempted": True,
        "headless_process_exit_classified": True,
        "fresh_result_artifact": True,
        "no_unexpected_error_lines": True,
    }
    return {
        "source_kind": (
            "upstream-tagged-example-exact-replay-plus-headless-external-cad"
        ),
        "build": {
            "pass": True,
            "source": {
                "repository": "gumyr/build123d",
                "tag": "v0.10.0",
                "commit": "a" * 40,
                "path": "examples/loft.py",
                "sha256": "b" * 64,
                "license": "Apache-2.0",
                "source_preserved": True,
                "display_stubbed_only": True,
            },
            "replays": [source_replay, copy.deepcopy(source_replay)],
            "native": native,
            "same_kernel_roundtrips": {
                "step": {
                    **native,
                    "volume_mm3": 1306.3404838415981,
                    "surface_area_mm2": 5285.0901805904705,
                },
                "brep": {
                    **native,
                    "volume_mm3": 1306.3405290340188,
                },
            },
            "checks": source_checks,
        },
        "external": {
            "pass": True,
            "execution_mode": "headless_python_api_synchronous_commands",
            "headless_flags": ["-nographics", "-batch"],
            "gui_daemon_enabled": False,
            "replays": external_replays,
            "cad_handoff_ready": True,
            "mesh_attempted": False,
            "solver_ready": False,
            "checks": external_checks,
            "process": {
                "acceptable": True,
                "result_artifact_fresh": True,
                "unexpected_error_lines": [],
                "owned_processes_remaining": 0,
            },
        },
        "timing_breakdown_s": {
            "two_source_replays_and_roundtrips": 8.1,
            "four_headless_external_imports": 10.7,
            "cross_kernel_diagnosis": 0.01,
            "evidence_finalization": 0.01,
        },
    }


def test_public_gate_accepts_bounded_repeatable_cad_handoff():
    result = json.loads(build123d_lofted_shell_handoff_gate(json.dumps(summary())))
    assert result["status"] == "ok"
    assert result["cad_handoff_ready"] is True
    assert result["mesh_ready"] is False
    assert result["solver_ready"] is False
    assert result["checks"]["external_volume_drift_is_nonzero_but_bounded"]


def test_public_gate_rejects_large_cross_kernel_volume_drift():
    row = summary()
    for replay in row["external"]["replays"]:
        replay["snapshot"]["volume_sum_mm3"] = 1350.0
    result = json.loads(build123d_lofted_shell_handoff_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["cad_handoff_ready"] is False
    assert result["checks"]["external_volume_drift_is_nonzero_but_bounded"] is False


def test_public_gate_rejects_solver_ready_overclaim():
    row = summary()
    row["external"]["solver_ready"] = True
    result = json.loads(build123d_lofted_shell_handoff_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["cad_handoff_without_mesh_or_solver_overclaim"] is False


def test_public_gate_rejects_nonpositive_external_entities():
    row = summary()
    row["external"]["replays"][2]["snapshot"]["positive_volume_count"] = 0
    result = json.loads(build123d_lofted_shell_handoff_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert (
        result["checks"]["external_positive_entity_counts_match_topology"]
        is False
    )


def test_source_gate_accepts_immutable_official_loft_replay():
    result = json.loads(build123d_loft_example_source_replay_gate(json.dumps(summary())))
    assert result["status"] == "ok"
    assert result["source"] == "examples/loft.py"
    assert result["solver_ready"] is False
    assert result["checks"]["eleven_profiles_two_openings_and_half_mm_shell_recorded"]


def test_source_gate_rejects_hidden_source_geometry_change():
    row = summary()
    row["build"]["replays"][0]["slice_count"] = 9
    result = json.loads(build123d_loft_example_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["two_exact_source_replays_match"] is False
    assert (
        result["checks"]["eleven_profiles_two_openings_and_half_mm_shell_recorded"]
        is False
    )


def test_source_gate_rejects_stale_external_process():
    row = summary()
    row["external"]["process"]["result_artifact_fresh"] = False
    result = json.loads(build123d_loft_example_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["external_process_is_fresh_classified_and_clean"] is False


def test_source_gate_rejects_modified_source_digest_contract():
    row = summary()
    row["build"]["source"]["source_preserved"] = False
    row["build"]["source"]["sha256"] = "0" * 64
    result = json.loads(build123d_loft_example_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["viewer_stub_only_and_source_preserved"] is False


def test_server_rejects_missing_external_replays():
    row = summary()
    del row["external"]["replays"]
    result = json.loads(build123d_lofted_shell_handoff_gate(json.dumps(row)))
    assert result["status"] == "invalid_input"
    assert "external.replays" in result["error"]


@pytest.mark.parametrize(
    "case_id",
    ["volume_budget", "topology", "command_failure", "solver_overclaim", "external_pass"],
)
def test_counterfactual_curriculum90_public(case_id):
    row = summary()
    if case_id == "volume_budget":
        row["external"]["replays"][0]["snapshot"]["volume_sum_mm3"] *= 1.1
    elif case_id == "topology":
        row["external"]["replays"][0]["snapshot"]["surface_count"] *= 2
    elif case_id == "command_failure":
        row["external"]["replays"][0]["command"]["returned"] = False
    elif case_id == "solver_overclaim":
        row["external"]["solver_ready"] = True
    else:
        row["external"]["pass"] = False
    result = json.loads(build123d_lofted_shell_handoff_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"


@pytest.mark.parametrize(
    "case_id",
    ["repository", "tag", "path", "license", "process"],
)
def test_counterfactual_curriculum90_source(case_id):
    row = summary()
    if case_id == "repository":
        row["build"]["source"]["repository"] = "example/other"
    elif case_id == "tag":
        row["build"]["source"]["tag"] = "v0.9.0"
    elif case_id == "path":
        row["build"]["source"]["path"] = "examples/other.py"
    elif case_id == "license":
        row["build"]["source"]["license"] = "unknown"
    else:
        row["external"]["process"]["acceptable"] = False
    result = json.loads(build123d_loft_example_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"


def test_generalization_v3s_rejects_nonpositive_external_surfaces():
    row = summary()
    row["external"]["replays"][1]["snapshot"]["positive_surface_count"] = 0
    result = json.loads(build123d_lofted_shell_handoff_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"


def test_generalization_v3s_rejects_nonstubbed_viewer_claim():
    row = summary()
    row["build"]["source"]["display_stubbed_only"] = False
    result = json.loads(build123d_loft_example_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
