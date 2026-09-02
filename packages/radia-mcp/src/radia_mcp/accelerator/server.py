"""
Accelerator MCP Server (radia_mcp.accelerator) — accelerator magnet
design, end-pole analytical chamfers, Radia validation case studies.
"""

import sys

from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool, register_topics_tool
from ..common.tool_group import CoarseToolRegistry

from .knowledge import get_accelerator_documentation, TOPICS
from ..common.lazy_call import lazy_callable
build_magnetic_trajectory_pair_gate = lazy_callable(".magnetic_trajectory_gate", "magnetic_trajectory_pair_gate", __package__)

mcp = FastMCP("mcp-server-accelerator")
_validation = CoarseToolRegistry(mcp, namespace="accelerator")


@mcp.tool()
def accelerator(topic: str = "all") -> str:
    """Accelerator magnet design with Radia + radia-mcp.

    Topics:
      "all"
      "end_pole"      - Analytical chamfer r(z) = ∆(1/2-z/L_f)^(1/n)
                        (Delferriere-de Menezes-Duperrier, SOLEIL)
      "kolkata"       - Kolkata SC Cyclotron case study (Pradhan 2007)
                        — Radia + TOSCA + Mathematica integration
      "rotating_coil" - Multipole measurement + 3D field reconstruction
      "two_plane_design" / "endpack_two_plane" / "endpack_cobake" /
      "sector_saturation" - Clebsch-hodograph accelerator pole-face design
                        examples used by the radia-em Clebsch hodograph line
    """
    return get_accelerator_documentation(topic)


@_validation.tool()
def accelerator_magnetic_trajectory_pair_gate(summary_json: str) -> str:
    """Gate paired charged-particle trajectories with magnetic field off/on.

    The gate requires measurable transverse deflection while speed, kinetic
    energy, transported current, boundary hits, and collision power remain
    closed, expressing that the magnetic Lorentz force does no work.
    """
    return build_magnetic_trajectory_pair_gate(summary_json)


_validation.install()


register_status_tool(
    mcp,
    server_name='mcp-server-accelerator',
    description='Accelerator physics: beam optics, dipole/quad/sext magnets, undulator/wiggler',
    subpackage='radia_mcp.accelerator',
    related_servers=["electromagnet", "fusion-reactor"],
    optional_deps=["radia"],
)


register_topics_tool(
    mcp,
    server_name='mcp-server-accelerator',
    topics=TOPICS,
)


def main():
    if "--selftest" in sys.argv:
        print("Accelerator MCP server self-test:")
        print(f"  accelerator('all'): {len(accelerator('all'))} chars")
        print("  PASSED")
        return
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
