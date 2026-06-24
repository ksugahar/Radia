"""
Accelerator MCP Server (radia_mcp.accelerator) — accelerator magnet
design, end-pole analytical chamfers, Radia validation case studies.
"""

import sys

from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool, register_topics_tool

from .knowledge import get_accelerator_documentation, TOPICS

mcp = FastMCP("mcp-server-accelerator")


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
