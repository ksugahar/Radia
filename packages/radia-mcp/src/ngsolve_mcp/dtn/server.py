"""MCP server for Dirichlet-to-Neumann (DtN) open-boundary methods (ngsolve-mcp.dtn)."""

from __future__ import annotations

import argparse
import asyncio

from mcp.server.fastmcp import FastMCP

from radia_mcp.fem.fem_reference_topics import DTN_TOPICS as _TOPICS

from radia_mcp.common import register_status_tool

mcp = FastMCP("ngsolve-mcp.dtn")


@mcp.tool()
def dtn_usage(topic: str = "overview") -> str:
    """Return Dirichlet-to-Neumann (DtN) guidance for the given topic."""
    key = topic.strip().lower()
    if key in _TOPICS:
        return _TOPICS[key]
    return f"Unknown topic '{topic}'. Available: {', '.join(sorted(_TOPICS))}"


# Registered after every tool: the status registration applies the
# annotation contract to the tools that exist at this point.
register_status_tool(
    mcp,
    "ngsolve-mcp.dtn",
    "Dirichlet-to-Neumann (DtN) open-boundary operator and coarse-mesh coupling for NGSolve.",
    "ngsolve_mcp.dtn",
    related_servers=["radia-mcp.ngsolve", "fem", "ngsolve-mcp.kelvin"],
)


def _selftest() -> None:
    tools = asyncio.run(mcp.list_tools())
    if "dtn_usage" not in {t.name for t in tools}:
        raise SystemExit("dtn_usage tool missing")
    print("ngsolve-mcp.dtn selftest passed")


def main() -> None:
    parser = argparse.ArgumentParser(description="ngsolve-mcp.dtn MCP server")
    parser.add_argument("--selftest", action="store_true", help="Verify tool registration")
    if parser.parse_args().selftest:
        from radia_mcp.common.utf8_stdout import use_utf8_stdout

        use_utf8_stdout()
        _selftest()
        return
    from radia_mcp.common.windows_stdio import run_stdio

    run_stdio(mcp)


if __name__ == "__main__":
    main()
