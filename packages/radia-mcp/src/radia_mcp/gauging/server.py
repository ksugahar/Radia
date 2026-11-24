"""MCP server for electromagnetic gauging (ngsolve-mcp.gauging)."""

from __future__ import annotations

import argparse
import asyncio

from mcp.server.fastmcp import FastMCP

from radia_mcp import __version__
from radia_mcp.common import register_status_tool
from radia_mcp.common.mcp_contract import apply_tool_contract

mcp = FastMCP("ngsolve-mcp.gauging")

register_status_tool(
    mcp,
    "ngsolve-mcp.gauging",
    "Electromagnetic gauging (Tree-Cotree, Coulomb/Lorenz, open-boundary, gauge-invariance gate).",
    "radia_mcp.gauging",
    related_servers=["radia-ngsolve", "fem", "differential-forms"],
)

_TOPICS = {
    "overview": (
        "ngsolve-mcp.gauging: Electromagnetic gauging reference and gates.\n"
        "Topics: tree_cotree, coulomb, open_boundary."
    ),
    "tree_cotree": (
        "Tree-Cotree gauging for HCurl magnetostatics/eddy-current formulations:\n"
        "Remove gradient-space nullspace (ker curl) by eliminating spanning-tree\n"
        "edge DoFs on the mesh graph before solving."
    ),
    "coulomb": (
        "Coulomb gauge (div A = 0) and grad-div stabilization:\n"
        "Enforce div A = 0 via Lagrange multiplier or augmented grad-div term."
    ),
    "open_boundary": (
        "Open-boundary gauging with Kelvin transformation or BEM-FEM coupling:\n"
        "Ensure exterior potential decay and compatible gauge across interfaces."
    ),
}


@mcp.tool()
def gauging_usage(topic: str = "overview") -> str:
    """Return electromagnetic gauging guidance for the given topic."""
    key = topic.strip().lower()
    if key in _TOPICS:
        return _TOPICS[key]
    return f"Unknown topic '{topic}'. Available: {', '.join(sorted(_TOPICS))}"


apply_tool_contract(mcp, server_name="ngsolve-mcp.gauging", version=__version__)


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
