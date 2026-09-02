from __future__ import annotations

import copy
import asyncio
import hashlib
import json
import os
import sys
from pathlib import Path

import pytest
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

from radia_mcp.fem.legacy_corpus_absorption import (
    validate_legacy_corpus_evidence,
)
from radia_mcp.fem.server import fem_legacy_corpus_absorption_gate, mcp


def _sha(packet: dict) -> str:
    payload = dict(packet)
    payload.pop("evidence_payload_sha256", None)
    raw = json.dumps(
        payload, ensure_ascii=True, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _packet() -> dict:
    packet = {
        "schema": "radia.legacy-solver-corpus-evidence.v1",
        "executed_at_utc": "2026-07-29T01:02:03+00:00",
        "execution_version": {
            "producer": "cae-ai-lab-corpus-absorption",
            "producer_version": "1",
            "solver": "not-launched",
        },
        "solver_launched": False,
        "live_dependency_required": False,
        "model_surface": {
            "input_count": 100,
            "unique_content_count": 90,
            "semantic_signature_count": 12,
            "parse_error_count": 0,
            "catalog_sha256": "1" * 64,
        },
        "automation_surface": {
            "script_count": 20,
            "unique_command_count": 8,
            "classified_command_count": 8,
            "unknown_command_count": 0,
            "catalog_sha256": "2" * 64,
        },
        "document_surface": {
            "document_count": 30,
            "unique_content_count": 25,
            "covered_document_count": 30,
            "failure_count": 0,
            "catalog_sha256": "3" * 64,
        },
        "topic_surface": {
            "discovered_topic_count": 9,
            "catalogued_topic_count": 9,
            "missing_topic_count": 0,
            "catalog_sha256": "4" * 64,
        },
    }
    packet["evidence_payload_sha256"] = _sha(packet)
    return packet


def _resign(packet: dict) -> dict:
    packet["evidence_payload_sha256"] = _sha(packet)
    return packet


def test_solver_neutral_corpus_packet_is_accepted() -> None:
    result = validate_legacy_corpus_evidence(_packet())
    assert result["status"] == "accepted"
    assert result["pass"] is True
    assert result["coverage"]["model_inputs"] == 100
    assert result["live_dependency_required"] is False


@pytest.mark.parametrize(
    ("surface", "field", "value"),
    [
        ("model_surface", "parse_error_count", 1),
        ("automation_surface", "unknown_command_count", 1),
        ("document_surface", "covered_document_count", 29),
        ("topic_surface", "missing_topic_count", 1),
    ],
)
def test_incomplete_surface_is_rejected(surface: str, field: str, value: int) -> None:
    packet = _packet()
    packet[surface][field] = value
    result = validate_legacy_corpus_evidence(_resign(packet))
    assert result["status"] == "rejected"
    assert result["pass"] is False


def test_digest_mismatch_is_rejected() -> None:
    packet = _packet()
    packet["model_surface"]["input_count"] += 1
    result = validate_legacy_corpus_evidence(packet)
    assert not result["checks"]["evidence_payload_sha256"]


def test_closed_world_packet_rejects_paths_or_provenance_fields() -> None:
    packet = _packet()
    packet["corpus_root"] = "X:/private/source"
    result = validate_legacy_corpus_evidence(_resign(packet))
    assert not result["checks"]["closed_world_top_level"]


def test_boolean_counters_are_rejected() -> None:
    packet = _packet()
    packet["model_surface"]["input_count"] = True
    result = validate_legacy_corpus_evidence(_resign(packet))
    assert not result["checks"]["model_surface_complete"]


def test_mcp_tool_returns_gate_result_and_is_read_only() -> None:
    result = json.loads(fem_legacy_corpus_absorption_gate(json.dumps(_packet())))
    assert result["status"] == "accepted"
    tool = mcp._tool_manager._tools["fem_validation_run"]
    catalog = mcp._tool_manager._tools["fem_validation_catalog"].fn(
        query="legacy_corpus"
    )
    assert catalog["operations"][0]["name"] == "fem_legacy_corpus_absorption_gate"
    assert tool.annotations.readOnlyHint is True
    assert tool.annotations.destructiveHint is False
    assert tool.annotations.idempotentHint is True


def test_mcp_tool_rejects_non_object_input() -> None:
    result = json.loads(fem_legacy_corpus_absorption_gate("[]"))
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
                    "name": "fem_legacy_corpus_absorption_gate",
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
