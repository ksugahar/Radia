"""MCP server for Air-Gap Element (AGE) rotating-machine coupling (ngsolve-mcp.age)."""

from __future__ import annotations

import argparse
import asyncio

from mcp.server.fastmcp import FastMCP

from radia_mcp import __version__
from radia_mcp.common import register_status_tool
from radia_mcp.common.mcp_contract import apply_tool_contract

mcp = FastMCP("ngsolve-mcp.age")

register_status_tool(
    mcp,
    "ngsolve-mcp.age",
    "Air-Gap Element (AGE) analytical/harmonic rotor-stator coupling for NGSolve.",
    "radia_mcp.age",
    related_servers=["radia-ngsolve", "fem", "motor"],
)

_TOPICS = {
    "overview": (
        "ngsolve-mcp.age: Air-Gap Element (AGE) for rotating electrical machines.\n"
        "Topics: harmonic_coupling, torque, periodic_motion."
    ),
    "harmonic_coupling": (
        "AGE replaces the physical air-gap mesh with an analytical Fourier-series\n"
        "harmonic coupling between stator and rotor cylindrical interfaces."
    ),
    "torque": (
        "Electromagnetic torque is computed directly from the AGE Fourier coefficients\n"
        "without remeshing across arbitrary rotor angles."
    ),
    "periodic_motion": (
        "Supports fractional-slot / anti-periodic sector symmetry and smooth\n"
        "angular sweeps for cogging torque and back-EMF."
    ),
}


@mcp.tool()
def age_usage(topic: str = "overview") -> str:
    """Return Air-Gap Element (AGE) guidance for the given topic."""
    key = topic.strip().lower()
    if key in _TOPICS:
        return _TOPICS[key]
    return f"Unknown topic '{topic}'. Available: {', '.join(sorted(_TOPICS))}"


apply_tool_contract(mcp, server_name="ngsolve-mcp.age", version=__version__)


def _selftest() -> None:
    tools = asyncio.run(mcp.list_tools())
    if "age_usage" not in {t.name for t in tools}:
        raise SystemExit("age_usage tool missing")
    print("ngsolve-mcp.age selftest passed")


def main() -> None:
    parser = argparse.ArgumentParser(description="ngsolve-mcp.age MCP server")
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
