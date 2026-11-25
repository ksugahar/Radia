"""MCP server for electromagnetic gauging (ngsolve-mcp.gauging)."""

from __future__ import annotations

import argparse
import asyncio

from mcp.server.fastmcp import FastMCP

from radia_mcp.common import register_status_tool

mcp = FastMCP("ngsolve-mcp.gauging")

_TOPICS = {
    "overview": (
        "ngsolve-mcp.gauging: Electromagnetic gauging reference and gates.\n"
        "Topics: tree_cotree, coil_av_tree, coulomb, open_boundary."
    ),
    "tree_cotree": (
        "Tree-Cotree gauging for HCurl magnetostatics/eddy-current formulations:\n"
        "Remove gradient-space nullspace (ker curl) by eliminating spanning-tree\n"
        "edge DoFs on the mesh graph before solving."
    ),
    "coil_av_tree": """# Tree gauge for a voltage-driven coil

Keep source and sink ports separate in the A-V space. Build the tree from
active free edges, excluding Dirichlet edges and periodic slaves. Contract
periodic copies consistently. Preserve tangential A on each connected SIBC
trace by contracting each trace component separately in the tree graph.

Choose a fixed admissible mesh-edge path for the terminal voltage. Keep it off
the SIBC trace and reject the setup if no admissible path exists. Do not use a
tree that short-circuits the ports or eliminates the prescribed terminal drop.
For disconnected free-edge graphs use one spanning tree per component.

Check terminal impedance and field energy under a second admissible tree.
Compare against an independent air-coil reference and check loss convergence
separately. The IH voltage scaling, production recipe and comparison criteria
remain in the IH application guidance.""",
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


# Registered after every tool: the status registration applies the
# annotation contract to the tools that exist at this point.
register_status_tool(
    mcp,
    "ngsolve-mcp.gauging",
    "Electromagnetic gauging (Tree-Cotree, Coulomb/Lorenz, open-boundary, gauge-invariance gate).",
    "ngsolve_mcp.gauging",
    related_servers=["radia-mcp.ngsolve", "fem", "differential-forms"],
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
