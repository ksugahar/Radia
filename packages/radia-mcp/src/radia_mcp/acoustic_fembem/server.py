"""MCP server for the readable MATLAB acoustic FEM-BEM education solver."""

from __future__ import annotations

import json
import sys

from mcp.server.fastmcp import FastMCP

from ..common import register_status_tool
from . import (
    acoustic_fembem_agent_guide as _agent_guide,
    acoustic_fembem_extension_contract as _extension_contract,
    acoustic_fembem_server_config as _server_config,
)


mcp = FastMCP("mcp-server-acoustic-fembem")


@mcp.tool()
def acoustic_fembem_agent_guide() -> str:
    """Guide for the readable P1 MATLAB acoustic FEM-BEM education solver."""
    return _agent_guide()


@mcp.tool()
def acoustic_fembem_extension_contract() -> str:
    """Inspect the official MATLAB MCP extension shipped by Radia."""
    return json.dumps(_extension_contract(), ensure_ascii=False, indent=2)


@mcp.tool()
def acoustic_fembem_server_config(
    project_root: str = "",
    profile: str = "existing",
) -> str:
    """Compose official MATLAB MCP arguments for the education solver."""
    try:
        payload = _server_config(project_root, profile)
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        payload = {
            "schema": "radia-mcp.acoustic-fembem-server-config/v1",
            "status": "error",
            "error": str(exc),
        }
    return json.dumps(payload, ensure_ascii=False, indent=2)


register_status_tool(
    mcp,
    server_name="mcp-server-acoustic-fembem",
    description="Readable MATLAB P1 FEM/BEM, CQ, and validation education solver",
    subpackage="radia_mcp.acoustic_fembem",
    related_servers=["radia-ngsolve", "bem", "fem"],
)


def main() -> None:
    if "--selftest" in sys.argv:
        contract = _extension_contract()
        assert contract["ok"] and contract["tool_count"] == 10
        print("acoustic FEM-BEM MCP server self-test: OK")
        return
    mcp.run()


if __name__ == "__main__":
    main()
