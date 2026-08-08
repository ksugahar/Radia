from __future__ import annotations

from copy import deepcopy
import json

from radia_mcp.build123d.face_first_perforation_gate import (
    face_first_perforation_handoff_gate,
    face_first_perforation_source_replay_gate,
)
from radia_mcp.build123d.server import (
    build123d_face_first_perforation_handoff_gate,
    build123d_face_first_perforation_source_replay_gate,
)


VOLUME = 168472.3918002237


def _public_summary() -> dict[str, object]:
    return {
        "construction": {
            "mode": "face_from_outer_wire_and_hole_wires_then_single_extrude",
            "requested_hole_count": 625,
            "generated_location_count": 625,
            "hole_wire_count": 625,
            "face_wire_count": 626,
            "hole_side_count": 6,
            "outer_side_count": 4,
        },
        "native": {"body_count": 1, "surface_count": 3756, "volume": VOLUME},
        "self_roundtrips": [
            {
                "format": "step",
                "body_count": 1,
                "surface_count": 3756,
                "volume": 168472.3918012877,
            },
            {
                "format": "brep",
                "body_count": 1,
                "surface_count": 3756,
                "volume": 168472.39180022335,
            },
        ],
        "external_imports": [
            {
                "mode": mode,
                "body_count": 1,
                "surface_count": 3756,
                "volume": 168472.39180031797,
            }
            for mode in ("noheal", "heal")
        ],
    }


def _source_summary() -> dict[str, object]:
    return {
        "source": {
            "kind": "upstream_native_example",
            "tag": "v0.10.0",
            "commit": "1" * 40,
            "expected_sha256": "2" * 64,
            "observed_sha256": "2" * 64,
            "preserved": True,
        },
        "execution": {
            "mode": "exact_source_with_display_stub",
            "native_run_count": 2,
            "native_runs_deterministic": True,
            "source_counts_observed": True,
        },
        "exports": {"step_sha256": "3" * 64, "brep_sha256": "4" * 64},
        "external": {
            "execution_mode": "python_api_headless",
            "headless_flags": ["-nographics", "-batch"],
            "import_modes": ["noheal", "heal"],
            "gui_daemon_enabled": False,
            "owned_processes_remaining": 0,
            "artifact_fresh": True,
            "process_exit_code": 1,
            "process_exit_policy": (
                "fresh_pass_artifact_plus_allowlisted_startup_diagnostics"
            ),
            "allowlisted_startup_diagnostics_only": True,
        },
        "public_gate_status": "ok",
    }


def test_live_shaped_public_and_source_evidence_are_accepted() -> None:
    public = face_first_perforation_handoff_gate(_public_summary())
    source = face_first_perforation_source_replay_gate(_source_summary())
    assert public["status"] == "ok"
    assert public["expected_surface_count"] == 3756
    assert public["metrics"]["max_volume_relative_error"] < 1.0e-10
    assert source["status"] == "ok"
    assert source["process_exit"]["classified_nonzero"] is True


def test_actual_hole_count_drift_is_rejected_even_with_plausible_import() -> None:
    payload = _public_summary()
    payload["construction"]["generated_location_count"] = 624
    payload["construction"]["hole_wire_count"] = 624
    result = face_first_perforation_handoff_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"]["requested_locations_match"] is False
    assert result["checks"]["generated_hole_wires_match"] is False
    assert result["checks"]["external_topology_matches"] is True


def test_self_consistent_but_wrong_downstream_topology_is_rejected() -> None:
    payload = _public_summary()
    for row in payload["self_roundtrips"] + payload["external_imports"]:
        row["surface_count"] = 3750
    result = face_first_perforation_handoff_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"]["self_roundtrip_topology_matches"] is False
    assert result["checks"]["external_topology_matches"] is False


def test_missing_external_import_mode_is_rejected() -> None:
    payload = _public_summary()
    payload["external_imports"] = payload["external_imports"][:1]
    result = face_first_perforation_handoff_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"]["heal_and_noheal_imported"] is False


def test_source_digest_gui_and_process_contradictions_are_rejected() -> None:
    payload = _source_summary()
    payload["source"]["observed_sha256"] = "5" * 64
    payload["external"]["gui_daemon_enabled"] = True
    payload["external"]["owned_processes_remaining"] = 1
    result = face_first_perforation_source_replay_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"]["source_digest_bound"] is False
    assert result["checks"]["no_gui_daemon"] is False
    assert result["checks"]["owned_processes_cleaned"] is False


def test_unclassified_nonzero_exit_is_rejected() -> None:
    payload = _source_summary()
    payload["external"]["allowlisted_startup_diagnostics_only"] = False
    result = face_first_perforation_source_replay_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"]["process_outcome_classified"] is False


def test_mcp_tools_dispatch_positive_and_negative_json() -> None:
    positive_public = json.loads(
        build123d_face_first_perforation_handoff_gate(
            json.dumps(_public_summary())
        )
    )
    negative = deepcopy(_public_summary())
    negative["construction"]["hole_wire_count"] = 624
    negative_public = json.loads(
        build123d_face_first_perforation_handoff_gate(json.dumps(negative))
    )
    positive_source = json.loads(
        build123d_face_first_perforation_source_replay_gate(
            json.dumps(_source_summary())
        )
    )
    assert positive_public["status"] == "ok"
    assert negative_public["status"] == "needs_attention"
    assert positive_source["status"] == "ok"


def test_mcp_tools_report_invalid_json() -> None:
    public = json.loads(build123d_face_first_perforation_handoff_gate("{"))
    source = json.loads(build123d_face_first_perforation_source_replay_gate("{"))
    assert public["status"] == "invalid_input"
    assert source["status"] == "invalid_input"
