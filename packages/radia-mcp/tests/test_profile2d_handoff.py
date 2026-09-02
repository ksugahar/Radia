from __future__ import annotations

import hashlib
import json

import pytest

from radia_mcp.radia_ngsolve.profile2d_handoff import profile2d_handoff_gate


def _canonical(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def _artifact(content: str, fmt: str) -> dict[str, str]:
    return {
        "content": content,
        "format": fmt,
        "sha256": hashlib.sha256(content.encode()).hexdigest(),
    }


def _square() -> dict:
    semantics = {
        "material_ids": ["air"],
        "boundary_ids": ["outer"],
        "conductor_ids": [],
    }
    return {
        "schema": "radia.profile2d-handoff.v1",
        "length_unit": "m",
        "semantics": semantics,
        "profile": {
            "nodes": [
                {"id": 0, "x_m": 0.0, "y_m": 0.0},
                {"id": 1, "x_m": 2.0, "y_m": 0.0},
                {"id": 2, "x_m": 2.0, "y_m": 1.0},
                {"id": 3, "x_m": 0.0, "y_m": 1.0},
            ],
            "segments": [
                {"id": 0, "start": 0, "end": 1, "boundary_id": "outer"},
                {"id": 1, "start": 1, "end": 2, "boundary_id": "outer"},
                {"id": 2, "start": 2, "end": 3, "boundary_id": "outer"},
                {"id": 3, "start": 3, "end": 0, "boundary_id": "outer"},
            ],
            "arcs": [],
            "regions": [{"id": 0, "x_m": 1.0, "y_m": 0.5, "material_id": "air"}],
        },
    }


def _abi(transport: str = "simulink_s_function") -> dict:
    result = {
        "schema": "radia.fixed-scalar-io.v1",
        "transport": transport,
        "inputs": [{"name": "drive_v", "unit": "V"}],
        "outputs": [
            {"name": "reaction_a", "unit": "A"},
            {"name": "loss_w", "unit": "W"},
        ],
        "input_width": 1,
        "output_width": 2,
        "fixed_width": True,
        "dynamic_paths": False,
        "request_contract_sha256": "a" * 64,
    }
    if transport == "simulink_s_function":
        result["sample_period_s"] = 1.0e-4
    return result


def test_line_loop_has_exact_green_area_and_perimeter() -> None:
    result = profile2d_handoff_gate(_square())
    assert result["status"] == "verified"
    assert result["measurements"]["exact_closed_loop_area_m2"] == pytest.approx(2.0)
    assert result["measurements"]["exact_closed_loop_perimeter_m"] == pytest.approx(6.0)


def test_signed_arc_loop_has_exact_circle_measurements() -> None:
    packet = _square()
    packet["profile"] = {
        "nodes": [
            {"id": 0, "x_m": -1.0, "y_m": 0.0},
            {"id": 1, "x_m": 1.0, "y_m": 0.0},
        ],
        "segments": [],
        "arcs": [
            {"id": 0, "start": 0, "end": 1, "sweep_deg": 180.0, "max_segment_deg": 5.0, "boundary_id": "outer"},
            {"id": 1, "start": 1, "end": 0, "sweep_deg": 180.0, "max_segment_deg": 5.0, "boundary_id": "outer"},
        ],
        "regions": [{"id": 0, "x_m": 0.0, "y_m": 0.0, "material_id": "air"}],
    }
    result = profile2d_handoff_gate(packet)
    assert result["measurements"]["exact_closed_loop_area_m2"] == pytest.approx(3.141592653589793)
    assert result["measurements"]["exact_closed_loop_perimeter_m"] == pytest.approx(2.0 * 3.141592653589793)


def test_step_and_semantics_are_separate_verified_artifacts() -> None:
    packet = _square()
    sidecar = _canonical(packet["semantics"])
    step = "ISO-10303-21;\nHEADER;ENDSEC;\nDATA;ENDSEC;\nEND-ISO-10303-21;"
    profile_sha = profile2d_handoff_gate(packet)["profile_sha256"]
    step_sha = hashlib.sha256(step.encode()).hexdigest()
    cad_measurement = _canonical({
        "schema": "radia.cad2d-measurement.v1",
        "step_sha256": step_sha,
        "profile_sha256": profile_sha,
        "area_m2": 2.0,
        "perimeter_m": 6.0,
    })
    packet["artifacts"] = {
        "step": _artifact(step, "STEP AP242"),
        "semantic_sidecar": _artifact(sidecar, "radia semantic sidecar v1"),
        "cad_measurement": _artifact(cad_measurement, "radia CAD measurement v1"),
    }
    result = profile2d_handoff_gate(packet)
    assert result["cad_ready"] is True
    assert result["policy"]["step_scope"] == "geometry_only"


def test_fixed_width_simulink_abi_is_realtime_ready() -> None:
    packet = _square()
    packet["execution_abi"] = _abi()
    result = profile2d_handoff_gate(packet)
    assert result["realtime_ready"] is True
    assert result["execution_abi"]["input_width"] == 1
    assert result["execution_abi"]["output_width"] == 2


def test_open_profile_is_rejected() -> None:
    packet = _square()
    packet["profile"]["segments"].pop()
    with pytest.raises(ValueError, match="dangling edge"):
        profile2d_handoff_gate(packet)


def test_isolated_point_property_is_allowed_but_untyped_unused_node_is_rejected() -> None:
    packet = _square()
    packet["semantics"]["point_property_ids"] = ["ring_current"]
    packet["profile"]["nodes"].append(
        {
            "id": 4,
            "x_m": 1.0,
            "y_m": 0.5,
            "point_property_id": "ring_current",
        }
    )
    result = profile2d_handoff_gate(packet)
    assert result["status"] == "verified"
    assert result["checks"]["isolated_nodes_have_point_semantics"] is True

    packet["profile"]["nodes"][-1].pop("point_property_id")
    with pytest.raises(ValueError, match="unused node"):
        profile2d_handoff_gate(packet)


def test_invalid_arc_is_rejected() -> None:
    packet = _square()
    packet["profile"] = {
        "nodes": [
            {"id": 0, "x_m": 0.0, "y_m": 0.0},
            {"id": 1, "x_m": 1.0, "y_m": 0.0},
            {"id": 2, "x_m": 0.0, "y_m": 1.0},
        ],
        "segments": [
            {"id": 1, "start": 1, "end": 2},
            {"id": 2, "start": 2, "end": 0},
        ],
        "arcs": [{"id": 0, "start": 0, "end": 1, "sweep_deg": 270.0, "max_segment_deg": 5.0}],
        "regions": [{"id": 0, "x_m": 0.2, "y_m": 0.2, "material_id": "air"}],
    }
    with pytest.raises(ValueError, match="sweep_deg"):
        profile2d_handoff_gate(packet)


def test_step_without_semantic_sidecar_is_rejected() -> None:
    packet = _square()
    step = "ISO-10303-21;\nHEADER;ENDSEC;\nDATA;ENDSEC;\nEND-ISO-10303-21;"
    packet["artifacts"] = {"step": _artifact(step, "STEP AP242")}
    with pytest.raises(ValueError, match="semantic sidecar"):
        profile2d_handoff_gate(packet)


def test_dynamic_or_variable_width_abi_is_rejected() -> None:
    packet = _square()
    packet["execution_abi"] = _abi()
    packet["execution_abi"]["dynamic_paths"] = True
    packet["execution_abi"]["output_width"] = 3
    with pytest.raises(ValueError, match="fixed-runtime"):
        profile2d_handoff_gate(packet)


def test_fem_mcp_registers_profile_handoff_gate() -> None:
    from radia_mcp.fem.server import mcp

    names = {tool.name for tool in mcp._tool_manager.list_tools()}
    assert "fem_validation_catalog" in names
    catalog = mcp._tool_manager._tools["fem_validation_catalog"].fn()
    operations = {item["name"] for item in catalog["operations"]}
    assert "fem_profile2d_handoff_gate" in operations
