"""MCP server for Dirichlet-to-Neumann (DtN) open-boundary methods (ngsolve-mcp.dtn)."""

from __future__ import annotations

import argparse
import asyncio

from mcp.server.fastmcp import FastMCP

from radia_mcp.common import register_status_tool

mcp = FastMCP("ngsolve-mcp.dtn")

register_status_tool(
    mcp,
    "ngsolve-mcp.dtn",
    "Dirichlet-to-Neumann (DtN) open-boundary operator and coarse-mesh coupling for NGSolve.",
    "radia_mcp.dtn",
    related_servers=["radia-ngsolve", "fem", "ngsolve-mcp.kelvin"],
)

_TOPICS = {
    "overview": (
        "ngsolve-mcp.dtn: Dirichlet-to-Neumann (DtN) open-boundary condition.\n"
        "Topics: spherical_harmonics, coarse_mesh, bem_fem_comparison."
    ),
    "spherical_harmonics": (
        "Exact DtN map on a spherical artificial boundary r = R0 expands exterior\n"
        "harmonic modes in spherical harmonics Y_l^m with radial decay (R0/r)^(l+1)."
    ),
    "coarse_mesh": (
        "Use curved high-order boundary elements on r = R0 so DtN quadrature\n"
        "remains accurate even on coarse interior meshes."
    ),
    "bem_fem_comparison": (
        "DtN is diagonal in spherical modes (ideal for spherical truncation);\n"
        "BEM-FEM Schur or Kelvin transformation handles arbitrary boundary shapes."
    ),
}


@mcp.tool()
def dtn_usage(topic: str = "overview") -> str:
    """Return Dirichlet-to-Neumann (DtN) guidance for the given topic."""
    key = topic.strip().lower()
    if key in _TOPICS:
        return _TOPICS[key]
    return f"Unknown topic '{topic}'. Available: {', '.join(sorted(_TOPICS))}"


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
