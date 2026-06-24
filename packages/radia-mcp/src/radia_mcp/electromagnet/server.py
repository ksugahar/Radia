"""
Accelerator Electromagnet MCP Server (radia_mcp.electromagnet)

Provides tools for:
- CoilBuilder workflow (racetrack, saddle, racetrack coils)
- Cubit hex mesh -> NGSolve curved -> Kelvin transform pipeline
- Omega-reduced scalar potential formulation
- Hantila polarization method (LU once, back-substitution iteration)
- B-input Play/Energy hysteresis models
- IMA (Image Method of Analysis) sign selection
- Field harmonics / multipole analysis
- radia-em Clebsch hodograph panel mode and accelerator design references

Usage:
    mcp-server-electromagnet              # Start MCP server (stdio transport)
    mcp-server-electromagnet --selftest   # Run self-test

Promoted from s:/mcp-server/mcp-server-electromagnet/ to public
radia-mcp on 2026-04-24 (single Radia monorepo).  General Kelvin
and NGSolve/BEM knowledge live in radia_mcp.radia_ngsolve.*; this
subpackage holds electromagnet-specific topics (CoilBuilder,
Hantila polarization, B-input hysteresis, IMA sign selection,
field-harmonics / multipole analysis).
"""

import sys

from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool, register_topics_tool

from .em_knowledge import get_electromagnet_documentation, TOPICS

mcp = FastMCP("mcp-server-electromagnet")


# ============================================================
# MCP Tools
# ============================================================

@mcp.tool()
def electromagnet_usage(topic: str = "overview") -> str:
    """
    Get accelerator electromagnet analysis documentation.

    Complete pipeline for accelerator magnet design and analysis:
    Radia CoilBuilder (Biot-Savart source) -> Cubit hex mesh ->
    NGSolve FEM (Omega-reduced + Kelvin) -> GMSH visualization.

    Covers dipole, quadrupole, and sextupole magnets with nonlinear
    iron yoke, hysteresis, and open boundary (Kelvin transformation).

    Args:
        topic: Documentation topic. Options:
            "overview"         - Pipeline architecture, unique capabilities
            "coilbuilder"      - CoilBuilder fluent API, examples
            "kelvin_workflow"  - Cubit -> NGSolve -> Kelvin pipeline
            "hantila"          - Hantila polarization (LU once)
            "hysteresis"       - B-input Play/Energy models
            "ima"              - Image Method sign selection
            "harmonics"        - Multipole analysis, FFT extraction
            "clebsch_hodograph" - radia-em Clebsch hodograph mode and links
                                  to accelerator pole-face design topics
            "all"              - Complete documentation
    """
    return get_electromagnet_documentation(topic)


# ============================================================
# MCP Prompts
# ============================================================

@mcp.prompt()
def new_electromagnet_simulation(magnet_type: str = "dipole") -> str:
    """Set up a new accelerator electromagnet simulation."""
    type_lower = magnet_type.strip().lower()

    base = (
        f"Set up an accelerator {magnet_type} magnet simulation.\n\n"
        "Use the electromagnet_usage tool for detailed documentation.\n\n"
        "General workflow:\n"
        "1. Define coil with CoilBuilder (no coil mesh needed)\n"
        "   - electromagnet_usage('coilbuilder') for API reference\n"
        "2. Create yoke + air + Kelvin domain in Cubit\n"
        "   - export netgen 'model.vol' order 2 overwrite\n"
        "3. NGSolve FEM solve (Omega-reduced scalar potential)\n"
        "   - electromagnet_usage('kelvin_workflow') for full code\n"
        "4. Visualize with GmshPostExport + coil STEP overlay\n"
        "5. Extract field harmonics at reference radius\n"
        "   - electromagnet_usage('harmonics') for FFT method\n\n"
    )

    if type_lower == "dipole":
        base += (
            "Dipole-specific notes:\n"
            "- 2-fold symmetry: allowed harmonics b_1, b_3, b_5, ...\n"
            "- Quarter model with IMA '+x-z' (Bz parallel to X, perp to Z)\n"
            "- Racetrack coil: add_straight + add_arc(180 deg)\n"
            "- Gap field estimate: B ~ mu_0 * NI / gap\n"
        )
    elif type_lower == "quadrupole":
        base += (
            "Quadrupole-specific notes:\n"
            "- 4-fold symmetry: allowed harmonics b_2, b_6, b_10, ...\n"
            "- 1/8 model with IMA (octant symmetry)\n"
            "- Saddle coils: 4 coils at +/-45 deg\n"
            "- Gradient: G = dBy/dx [T/m]\n"
        )
    elif type_lower == "sextupole":
        base += (
            "Sextupole-specific notes:\n"
            "- 6-fold symmetry: allowed harmonics b_3, b_9, b_15, ...\n"
            "- 1/12 model with IMA\n"
            "- 6 coils at 30 deg intervals\n"
        )
    else:
        base += (
            f"Custom magnet type '{magnet_type}':\n"
            "- Determine symmetry order and allowed harmonics\n"
            "- Choose IMA signs based on dominant field component\n"
            "  electromagnet_usage('ima') for sign selection rules\n"
            "- For nonlinear iron: electromagnet_usage('hantila')\n"
            "- For hysteresis: electromagnet_usage('hysteresis')\n"
        )

    return base


# ============================================================
# Entry point
# ============================================================



register_status_tool(
    mcp,
    server_name='mcp-server-electromagnet',
    description='Accelerator electromagnet: CoilBuilder, Hantila, Play/Energy hysteresis',
    subpackage='radia_mcp.electromagnet',
    related_servers=["motor", "accelerator", "magnetic-materials"],
    optional_deps=["radia", "ngsolve"],
)


register_topics_tool(
    mcp,
    server_name='mcp-server-electromagnet',
    topics=TOPICS,
)


def main():
    if "--selftest" in sys.argv:
        print("Electromagnet MCP server self-test:")
        topics = [
            "overview", "coilbuilder", "kelvin_workflow",
            "hantila", "hysteresis", "ima", "harmonics",
            "clebsch_hodograph",
        ]
        for t in topics:
            result = electromagnet_usage(t)
            print(f"  electromagnet_usage('{t}'): {len(result)} chars")
            assert len(result) > 100, f"Topic '{t}' too short"
        # Test prompt
        prompt = new_electromagnet_simulation("dipole")
        print(f"  new_electromagnet_simulation('dipole'): {len(prompt)} chars")
        assert "CoilBuilder" in prompt
        print("  PASSED")
        return

    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
