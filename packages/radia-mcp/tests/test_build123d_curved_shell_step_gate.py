import copy
import json

from radia_mcp.build123d.server import (
    build123d_curved_shell_step_semantics_gate,
    build123d_tea_cup_source_contract_gate,
)


def summary() -> dict:
    topology = {
        "solid_count": 1,
        "face_count": 28,
        "edge_count": 68,
        "vertex_count": 42,
        "is_valid": True,
    }
    imports = [
        {
            "mode": mode,
            "volume_count": 1,
            "surface_count": 28,
            "curve_count": 68,
            "vertex_count": 42,
            "native_volume_relative_error": 0.1844044769569793,
            "native_area_relative_error": 0.014877834492298582,
        }
        for mode in ("noheal", "heal")
    ]
    return {
        "source_kind": "upstream_native_build123d_example",
        "build": {
            "pass": True,
            "versions": {
                "build123d": "0.10.0",
                "upstream_commit": "a" * 40,
            },
            "source_native_example": {
                "source": "examples/tea_cup.py",
                "source_sha256": "b" * 64,
                "viewer_change": "ocp_vscode.show stubbed; geometry executed unchanged",
                "features": [
                    "spline_profile_revolve",
                    "open_face_shell_offset",
                    "bottom_fusion",
                    "edge_fillet",
                    "intersection_derived_handle_contacts",
                    "nonplanar_handle_sweep",
                ],
            },
            "native": {**topology, "volume_mm3": 130326.769, "surface_area_mm2": 87984.572},
            "step_reimport": {
                **topology,
                "volume_mm3": 130327.080,
                "surface_area_mm2": 87984.588,
            },
            "errors": {
                "volume_relative": 2.3824373892369937e-6,
                "surface_area_relative": 1.7742000314622792e-7,
            },
            "checks": {
                "source_immutable": True,
                "official_volume_assertion_reproduced": True,
            },
            "disposition": "same_kernel_topology_preserved_but_mass_not_exact",
        },
        "external": {
            "pass": True,
            "execution_mode": "python_api_headless_synchronous_commands",
            "headless_flags": ["-nographics", "-batch"],
            "persistent_gui_started": False,
            "imports": imports,
            "mode_volume_relative_spread": 0.0,
            "mode_area_relative_spread": 0.0,
            "disposition": "external_kernel_mass_loss_with_topology_preserved_not_solver_ready",
            "process": {
                "exit_code": 3,
                "acceptable": True,
                "result_artifact_fresh": True,
                "known_headless_diagnostics_only": True,
                "unexpected_error_lines": [],
                "owned_processes_remaining": 0,
            },
        },
        "timing_breakdown_s": {
            "source_build_and_roundtrip": 1.9,
            "external_noheal": 0.34,
            "external_heal": 0.33,
            "diagnosis": 0.01,
        },
    }


def test_public_gate_accepts_validated_negative_with_topology_preserved():
    result = json.loads(build123d_curved_shell_step_semantics_gate(json.dumps(summary())))
    assert result["status"] == "ok"
    assert result["solver_ready"] is False
    assert result["diagnosis"] == "topology_preserved_external_curved_surface_mass_loss"
    assert result["checks"]["external_topology_still_matches_native"] is True


def test_public_gate_rejects_topology_count_as_mass_equivalence():
    row = summary()
    for import_row in row["external"]["imports"]:
        import_row["native_volume_relative_error"] = 1.0e-8
        import_row["native_area_relative_error"] = 1.0e-8
    result = json.loads(build123d_curved_shell_step_semantics_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["external_topology_still_matches_native"] is True
    assert result["checks"]["external_mass_loss_is_material"] is False


def test_public_gate_rejects_heal_noheal_spread_masquerading_as_semantic_loss():
    row = summary()
    row["external"]["mode_volume_relative_spread"] = 0.02
    result = json.loads(build123d_curved_shell_step_semantics_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["healing_does_not_change_external_result"] is False


def test_source_gate_accepts_immutable_upstream_replay_and_negative_disposition():
    result = json.loads(build123d_tea_cup_source_contract_gate(json.dumps(summary())))
    assert result["status"] == "ok"
    assert result["solver_ready"] is False
    assert result["checks"]["viewer_only_stub_preserves_geometry"] is True
    assert result["checks"]["public_semantics_gate_passed"] is True


def test_source_gate_rejects_geometry_edit_hidden_as_viewer_stub():
    row = summary()
    row["build"]["source_native_example"]["viewer_change"] = "geometry simplified"
    result = json.loads(build123d_tea_cup_source_contract_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["viewer_only_stub_preserves_geometry"] is False


def test_source_gate_rejects_stale_external_process():
    row = copy.deepcopy(summary())
    row["external"]["process"]["result_artifact_fresh"] = False
    result = json.loads(build123d_tea_cup_source_contract_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["external_process_is_fresh_and_classified"] is False


def test_server_rejects_missing_external_imports():
    row = summary()
    del row["external"]["imports"]
    result = json.loads(build123d_curved_shell_step_semantics_gate(json.dumps(row)))
    assert result["status"] == "invalid_input"
    assert "external.imports" in result["error"]
