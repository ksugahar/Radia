"""MCP server for NGSolve core FEM spaces, solvers, and taskmanager (ngsolve-mcp.core)."""

from __future__ import annotations

import argparse
import asyncio

from mcp.server.fastmcp import FastMCP

from radia_mcp.fem.fem_reference_topics import CORE_TOPICS as _TOPICS

from radia_mcp.common import register_status_tool
from radia_mcp.common.windows_stdio import run_stdio
from radia_mcp.meta.catalog import CATALOG

mcp = FastMCP("ngsolve-mcp.core")


@mcp.tool()
def ngsolve_usage(topic: str = "index") -> str:
    """Return NGSolve core usage guidance for the given topic."""
    key = topic.strip().lower()
    if key in _TOPICS:
        return _TOPICS[key]
    return f"Unknown topic '{topic}'. Available: {', '.join(sorted(_TOPICS))}"


register_status_tool(
    mcp,
    "ngsolve-mcp.core",
    "NGSolve finite element spaces, solvers and TaskManager.",
    "ngsolve_mcp.core",
    related_servers=CATALOG["ngsolve-mcp.core"]["related"],
)


def _selftest() -> None:
    tools = asyncio.run(mcp.list_tools())
    if "ngsolve_usage" not in {t.name for t in tools}:
        raise SystemExit("ngsolve_usage tool missing")
    print("ngsolve-mcp.core selftest passed")


def main() -> None:
    parser = argparse.ArgumentParser(description="ngsolve-mcp.core MCP server")
    parser.add_argument("--selftest", action="store_true", help="Verify tool registration")
    if parser.parse_args().selftest:
        _selftest()
        return
    run_stdio(mcp)


if __name__ == "__main__":
    main()
