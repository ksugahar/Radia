"""MCP server for electromagnetic gauging (ngsolve-mcp.gauging)."""

from __future__ import annotations

import argparse
import asyncio

from mcp.server.fastmcp import FastMCP

from radia_mcp.fem.fem_reference_topics import GAUGING_TOPICS as _TOPICS

from radia_mcp.common import register_status_tool
from radia_mcp.meta.catalog import CATALOG

mcp = FastMCP("ngsolve-mcp.gauging")


@mcp.tool()
def gauging_usage(topic: str = "overview") -> str:
    """Return electromagnetic gauging guidance for the given topic."""
    key = topic.strip().lower()
    if key in _TOPICS:
        return _TOPICS[key]
    return f"Unknown topic '{topic}'. Available: {', '.join(sorted(_TOPICS))}"


# Registered after every tool: the status registration applies the
# annotation contract to the tools that exist at this point.
register_status_tool(
    mcp,
    "ngsolve-mcp.gauging",
    "Electromagnetic gauging (Tree-Cotree, Coulomb/Lorenz, open-boundary, gauge-invariance gate).",
    "ngsolve_mcp.gauging",
    related_servers=CATALOG["ngsolve-mcp.gauging"]["related"],
)


def _selftest() -> None:
    tools = asyncio.run(mcp.list_tools())
    names = {t.name for t in tools}
    if "gauging_usage" not in names:
        raise SystemExit("gauging_usage tool missing")
    print("ngsolve-mcp.gauging selftest passed")


def main() -> None:
    parser = argparse.ArgumentParser(description="ngsolve-mcp.gauging MCP server")
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
