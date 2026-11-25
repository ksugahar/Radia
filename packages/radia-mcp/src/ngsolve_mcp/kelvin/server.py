"""MCP server for Kelvin transformation open-boundary methods (ngsolve-mcp.kelvin)."""

from __future__ import annotations

import argparse
import asyncio

from mcp.server.fastmcp import FastMCP

from radia_mcp.fem.fem_reference_topics import KELVIN_TOPICS as _TOPICS

from radia_mcp.common import register_status_tool
from radia_mcp.meta.catalog import CATALOG

mcp = FastMCP("ngsolve-mcp.kelvin")


@mcp.tool()
def kelvin_usage(topic: str = "overview") -> str:
    """Return Kelvin transformation guidance for the given topic."""
    key = topic.strip().lower()
    if key in _TOPICS:
        return _TOPICS[key]
    return f"Unknown topic '{topic}'. Available: {', '.join(sorted(_TOPICS))}"


# Registered after every tool: the status registration applies the
# annotation contract to the tools that exist at this point.
register_status_tool(
    mcp,
    "ngsolve-mcp.kelvin",
    "Kelvin transformation open-boundary exterior inversion for NGSolve/Radia.",
    "ngsolve_mcp.kelvin",
    related_servers=CATALOG["ngsolve-mcp.kelvin"]["related"],
)


def _selftest() -> None:
    tools = asyncio.run(mcp.list_tools())
    if "kelvin_usage" not in {t.name for t in tools}:
        raise SystemExit("kelvin_usage tool missing")
    print("ngsolve-mcp.kelvin selftest passed")


def main() -> None:
    parser = argparse.ArgumentParser(description="ngsolve-mcp.kelvin MCP server")
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
