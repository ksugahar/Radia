from __future__ import annotations

import asyncio
import copy
import hashlib
import json
import os
import sys
from pathlib import Path

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

from radia_mcp.fem.legacy_signature_migration import (
    validate_legacy_signature_migration_evidence,
)
from radia_mcp.fem.server import fem_legacy_signature_migration_gate, mcp


def _sha(packet: dict) -> str:
    payload = dict(packet)
    payload.pop("evidence_payload_sha256", None)
    raw = json.dumps(
        payload, ensure_ascii=True, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _packet() -> dict:
    packet = {
        "schema": "radia.legacy-signature-migration-evidence.v1",
        "executed_at_utc": "2026-07-29T01:02:03+00:00",
        "execution_version": {
            "producer": "cae-ai-lab-signature-migration",
            "producer_version": "1",
            "solver": "not-launched",
        },
        "solver_launched": False,
        "live_dependency_required": False,
        "signature_count": 72,
        "assigned_signature_count": 72,
        "lane_counts": {
            "axifem_eddy": 5,
            "axifem_henrotte": 27,
            "planar_hcurl_eddy_bubble": 12,
            "planar_magnetostatic_age": 24,
            "scalar_h1_electrostatic": 4,
        },
        "requirement_counts": {
            "age_motion_interface": 1,
            "boundary_operator_identity": 35,
            "circuit_excitation_identity": 19,
            "curved_profile_geometry": 56,
            "kelvin_open_boundary": 12,
            "nonlinear_material_curve": 14,
            "point_source_or_probe_identity": 38,
        },
        "native_emitter_candidate_count": 0,
        "false_solver_ready_count": 0,
        "identity_mismatch_count": 0,
        "conversion_error_count": 0,
        "script_compile_failure_count": 0,
        "signature_catalog_sha256": "a" * 64,
    }
    packet["evidence_payload_sha256"] = _sha(packet)
    return packet


def _resign(packet: dict) -> dict:
    packet["evidence_payload_sha256"] = _sha(packet)
    return packet


def test_solver_neutral_signature_packet_is_accepted() -> None:
    result = validate_legacy_signature_migration_evidence(_packet())
    assert result["status"] == "accepted"
    assert result["signature_count"] == 72
    assert result["native_emitter_candidate_count"] == 0


def test_false_ready_or_count_drift_is_rejected() -> None:
    false_ready = _packet()
    false_ready["false_solver_ready_count"] = 1
    result = validate_legacy_signature_migration_evidence(_resign(false_ready))
    assert result["status"] == "rejected"
    assert not result["checks"]["no_false_ready_or_conversion_failure"]

    drift = _packet()
    drift["lane_counts"]["axifem_eddy"] += 1
    result = validate_legacy_signature_migration_evidence(_resign(drift))
    assert result["status"] == "rejected"
    assert not result["checks"]["all_signatures_assigned"]


def test_closed_world_packet_rejects_paths_or_unknown_lanes() -> None:
    path_leak = _packet()
    path_leak["corpus_root"] = "X:/private/source"
    result = validate_legacy_signature_migration_evidence(_resign(path_leak))
    assert not result["checks"]["closed_world_top_level"]

    lane = _packet()
    lane["lane_counts"]["unknown_solver"] = 0
    result = validate_legacy_signature_migration_evidence(_resign(lane))
    assert not result["checks"]["lane_vocabulary_closed"]


def test_digest_mismatch_is_rejected() -> None:
    packet = _packet()
    packet["signature_count"] += 1
    result = validate_legacy_signature_migration_evidence(packet)
    assert not result["checks"]["evidence_payload_sha256"]


def test_mcp_tool_returns_gate_result_and_is_read_only() -> None:
    result = json.loads(fem_legacy_signature_migration_gate(json.dumps(_packet())))
    assert result["status"] == "accepted"
    tool = mcp._tool_manager._tools["fem_validation_run"]
    catalog = mcp._tool_manager._tools["fem_validation_catalog"].fn(
        query="legacy_signature"
    )
    assert catalog["operations"][0]["name"] == "fem_legacy_signature_migration_gate"
    assert tool.annotations.readOnlyHint is True
    assert tool.annotations.destructiveHint is False
    assert tool.annotations.idempotentHint is True


def test_mcp_tool_rejects_non_object_input() -> None:
    result = json.loads(fem_legacy_signature_migration_gate("[]"))
    assert result["status"] == "invalid_input"


async def _probe_stdio(packet: dict) -> dict:
    package_root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    source_root = package_root / "src"
    env["PYTHONPATH"] = os.pathsep.join(
        [str(source_root), env.get("PYTHONPATH", "")]
    ).rstrip(os.pathsep)
    env["RADIA_MCP_TOOL_PROFILE"] = "core"
    params = StdioServerParameters(
        command=sys.executable,
        args=["-m", "radia_mcp.fem.server"],
        cwd=str(package_root),
        env=env,
    )
    async with stdio_client(params) as (read_stream, write_stream):
        async with ClientSession(read_stream, write_stream) as session:
            initialized = await session.initialize()
            listed = await session.list_tools()
            called = await session.call_tool(
                "fem_validation_run",
                {
                    "name": "fem_legacy_signature_migration_gate",
                    "arguments": {"evidence_json": json.dumps(packet)},
                },
            )
            return {
                "server_name": initialized.serverInfo.name,
                "listed": any(
                    tool.name == "fem_validation_run"
                    for tool in listed.tools
                ),
                "is_error": bool(called.isError),
                "has_content": bool(called.content),
            }


def test_changed_tool_passes_real_stdio_protocol() -> None:
    result = asyncio.run(asyncio.wait_for(_probe_stdio(_packet()), timeout=45))
    assert result == {
        "server_name": "mcp-server-fem",
        "listed": True,
        "is_error": False,
        "has_content": True,
    }
