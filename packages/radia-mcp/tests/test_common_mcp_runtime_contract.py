"""Fleet-wide runtime-contract tests for radia-mcp servers."""

from __future__ import annotations

import asyncio
import importlib.util
import json
from pathlib import Path

from mcp.server.fastmcp import FastMCP
from mcp.types import ToolAnnotations

from radia_mcp.common.mcp_contract import (
    ARTIFACT_IDENTITY_SCHEMA,
    CAPABILITY_MATURITY_LEVELS,
    SCHEMA,
    validate_solver_artifact_identity,
)
from radia_mcp.common.server_hardening import install_call_log
from radia_mcp.common.status import register_status_tool


def test_status_registration_completes_and_exposes_runtime_contract(
    monkeypatch,
):
    monkeypatch.setenv("RADIA_MCP_CALL_LOG", "0")
    mcp = FastMCP("contract-test")

    @mcp.tool()
    def inspect_value() -> dict[str, object]:
        """Inspect and return local values."""
        return {"status": "ok"}

    @mcp.tool()
    def mystery_operation():
        """An intentionally ambiguous operation."""
        return "ok"

    explicit = ToolAnnotations(
        readOnlyHint=True,
        destructiveHint=False,
        idempotentHint=False,
        openWorldHint=True,
    )

    @mcp.tool(annotations=explicit)
    def upstream_reference() -> dict[str, object]:
        """Read one remote reference."""
        return {"status": "ok"}

    register_status_tool(
        mcp,
        server_name="mcp-server-contract-test",
        description="runtime contract fixture",
        subpackage="radia_mcp.meta",
    )

    tools = mcp._tool_manager._tools
    assert tools["inspect_value"].annotations.readOnlyHint is True
    assert tools["mystery_operation"].annotations.destructiveHint is True
    assert (
        tools["mystery_operation"].meta["caeai.annotation_source"]
        == "conservative-default"
    )
    assert tools["upstream_reference"].annotations == explicit
    assert tools["upstream_reference"].meta["caeai.annotation_source"] \
        == "explicit-decorator"

    status_tool = tools["contract_test_status"]
    assert status_tool.fn_metadata.output_schema is not None
    assert status_tool.meta["caeai.control_plane"] == "status"
    assert status_tool.description.startswith("Status / introspection")
    payload = status_tool.fn()
    assert payload["schema"] == "radia-mcp.server-status.v2"
    assert payload["status"] == "ready"
    assert payload["runtime_contract"]["schema"] == SCHEMA
    assert payload["runtime_contract"]["complete"] is True
    assert payload["runtime_contract"]["n_tools"] == len(tools)
    assert payload["client_connection"] == {
        "connection_count": 0,
        "latest": None,
        "recent": [],
    }
    assert tools["inspect_value"].meta["caeai.maturity"] == "artifact_gate"
    assert (
        tools["contract_test_status"].meta["caeai.maturity"]
        == "knowledge_only"
    )
    assert set(payload["runtime_contract"]["capability_maturity"]) <= set(
        CAPABILITY_MATURITY_LEVELS
    )
    assert payload["runtime_contract"]["artifact_identity_contract"][
        "canonical_container"
    ] == ".hdf5"
    assert payload["runtime_provenance"]["module_file"].endswith("server.py")
    assert len(
        payload["runtime_provenance"]["module_sha256_at_registration"]
    ) == 64
    assert len(payload["runtime_provenance"]["module_file_sha256"]) == 64
    assert payload["runtime_provenance"][
        "source_changed_since_registration"
    ] is False
    assert payload["runtime_provenance"]["mcp_sdk_version"]


def test_call_log_is_lazy_idempotent_and_does_not_record_values(
    monkeypatch, tmp_path
):
    monkeypatch.setenv("RADIA_MCP_STATE_DIR", str(tmp_path))
    monkeypatch.setenv("RADIA_MCP_CALL_LOG", "1")
    mcp = FastMCP("log-test")

    @mcp.tool()
    def echo(secret: str, values: list[int]) -> dict[str, object]:
        return {"length": len(secret), "count": len(values)}

    assert install_call_log(mcp, "calls.jsonl") is True
    assert install_call_log(mcp, "second.jsonl") is False
    assert not (tmp_path / "logs").exists()

    asyncio.run(
        mcp.call_tool(
            "echo",
            {"secret": "super-secret-value", "values": [1, 2, 3]},
        )
    )
    log = tmp_path / "logs" / "calls.jsonl"
    text = log.read_text(encoding="utf-8")
    assert "super-secret-value" not in text
    assert "[1, 2, 3]" not in text
    record = json.loads(text)
    assert record["schema"] == "radia-mcp.tool-call.v1"
    assert record["args"] == {
        "secret": {"type": "str", "length": 18},
        "values": {"type": "list", "length": 3},
    }
    assert record["ok"] is True


def test_meta_server_passes_real_stdio_runtime_contract():
    script = Path(__file__).resolve().parents[1] / "tools" / "smoke_mcp_stdio.py"
    spec = importlib.util.spec_from_file_location("radia_mcp_stdio_probe", script)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    result = module.probe_server("meta", timeout=45)
    assert result["server_name"] == "mcp-server-radia-meta"
    assert result["n_tools"] >= 5
    assert result["structured_status"] is True
    assert result["module_sha256_at_registration"]
    assert result["module_file_sha256"]
    assert result["source_changed_since_registration"] is False
    assert result["server_version"] == result["distribution"]["version"]
    assert result["client_connection"]["connection_count"] == 1
    assert result["client_connection"]["latest"]["name"] \
        == "cae-lab-radia-contract-probe"
    assert result["artifact_gate"] == {
        "positive": "accepted",
        "negative": "rejected",
    }


def _valid_solver_artifact() -> dict:
    return {
        "schema": ARTIFACT_IDENTITY_SCHEMA,
        "artifact_id": "field-run-001",
        "created_at_utc": "2026-09-03T00:00:00Z",
        "producer": {"name": "example-runner", "version": "1.2.3"},
        "solver": {"name": "example-solver", "version": "4.5.6"},
        "status": "complete",
        "coordinate_system": "cartesian-right-handed",
        "unit_system": "SI",
        "input_sha256": "a" * 64,
        "result_sha256": "b" * 64,
        "total_compute_seconds": 12.5,
        "timing_breakdown_s": {"solve": 10.0, "verify": 2.5},
        "mesh_based": True,
        "mesh_sha256": "c" * 64,
        "study_type": "transient",
        "time_axis": {"unit": "s", "count": 11, "start": 0.0, "end": 1.0},
    }


def test_solver_artifact_identity_accepts_complete_hdf5_contract():
    result = validate_solver_artifact_identity(_valid_solver_artifact())
    assert result["accepted"] is True
    assert result["status"] == "accepted"
    assert result["errors"] == []


def test_solver_artifact_identity_rejects_plausible_but_ambiguous_result():
    artifact = _valid_solver_artifact()
    artifact["unit_system"] = "unknown"
    artifact["result_sha256"] = "stale-result"
    artifact.pop("mesh_sha256")
    artifact["time_axis"] = {"unit": "s", "count": 11}
    result = validate_solver_artifact_identity(artifact)
    assert result["accepted"] is False
    assert result["status"] == "rejected"
    assert "unit_system:SI_required" in result["errors"]
    assert "result_sha256:sha256_required" in result["errors"]
    assert "mesh_sha256:required_for_mesh_based_result" in result["errors"]
    assert "time_axis:unit_count_start_end_required" in result["errors"]


def test_fleet_audit_writes_portable_hdf5_report(tmp_path):
    h5py = __import__("h5py")
    from tools.audit_mcp_fleet import write_report

    report = {
        "schema": "cae-ai-lab.mcp-fleet-audit.v1",
        "artifact_id": "fleet-fixture",
        "created_at_utc": "2026-09-03T00:00:00+00:00",
        "producer": {"name": "audit_mcp_fleet", "version": "1.0"},
        "status": "passed",
        "all_passed": True,
        "n_servers": 1,
        "total_compute_seconds": 0.1,
        "servers": [{"id": "fixture", "ok": True}],
    }
    target = tmp_path / "fleet.hdf5"
    write_report(report, target)
    with h5py.File(target, "r") as handle:
        assert handle.attrs["schema"] == report["schema"]
        assert handle.attrs["artifact_id"] == "fleet-fixture"
        assert "report_json" in handle
        assert bool(handle["servers/fixture"].attrs["ok"]) is True
