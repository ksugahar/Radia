"""Keep package metadata and minimal CI on the supported FastMCP SDK line."""

import asyncio
import importlib.metadata
import json
import os
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.10
    from importlib import import_module

    tomllib = import_module("tomli")


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
SDK_REQUIREMENT = "mcp>=1.20.0,<2"


def test_mcp_sdk_dependency_declares_supported_floor():
    project = tomllib.loads(
        (PACKAGE_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    )
    assert SDK_REQUIREMENT in project["project"]["dependencies"]



def test_supported_sdk_registers_metadata_lists_schema_and_calls_tool():
    from mcp.server.fastmcp import FastMCP
    from cae_mcp_core.common.server_hardening import ANN_READONLY

    expected = os.environ.get("RADIA_MCP_EXPECTED_SDK")
    if expected:
        assert importlib.metadata.version("mcp") == expected
    server = FastMCP("sdk-contract")

    def increment(value: int, step: int = 1) -> dict:
        """Increment a number without external state."""
        return {"value": value + step}

    # Same metadata arguments used when Radia refreshes registered tools.
    server.add_tool(increment, annotations=ANN_READONLY, icons=None,
                    meta={"radia": {"contract": "sdk-floor"}})

    async def exercise():
        tools = await server.list_tools()
        tool = next(item for item in tools if item.name == "increment")
        assert tool.meta == {"radia": {"contract": "sdk-floor"}}
        assert tool.annotations.readOnlyHint is True
        assert tool.inputSchema["required"] == ["value"]
        assert tool.inputSchema["properties"]["step"]["default"] == 1
        result = await server.call_tool("increment", {"value": "4"})
        # FastMCP may return content alone or (content, structured output).
        content = result[0] if isinstance(result, tuple) else result
        assert json.loads(content[0].text) == {"value": 5}

    asyncio.run(exercise())
