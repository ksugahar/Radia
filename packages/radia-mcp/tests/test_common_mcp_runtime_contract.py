"""Fleet-wide runtime-contract tests for radia-mcp servers."""

from __future__ import annotations

import asyncio
import json

from mcp.server.fastmcp import FastMCP
from mcp.types import ToolAnnotations

from radia_mcp.common.mcp_contract import SCHEMA
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
    from tools.smoke_mcp_stdio import probe_server

    result = probe_server("meta", timeout=45)
    assert result["server_name"] == "mcp-server-radia-meta"
    assert result["n_tools"] >= 5
    assert result["structured_status"] is True
    assert result["module_sha256_at_registration"]
    assert result["module_file_sha256"]
    assert result["source_changed_since_registration"] is False
    assert result["server_version"] == result["distribution"]["version"]
