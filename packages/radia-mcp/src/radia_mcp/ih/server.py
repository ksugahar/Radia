"""
Induction Heating MCP Server (radia_mcp.ih)

Provides tools for:
- IH-specific linting (eddy current, SIBC, ESIM, Bessel functions)
- SIBC method selection (BEM Scalar BIE vs FEM scattered-field)
- ESIM nonlinear surface impedance documentation
- Biot-Savart coil field computation
- Screening physics and Karl iteration
- IH simulation workflow (EM -> thermal)

Usage:
    mcp-server-ih              # Start MCP server (stdio transport)
    mcp-server-ih --selftest   # Run self-test

Promoted from s:/mcp-server/mcp-server-ih/ to public radia-mcp on
2026-04-24 (single Radia monorepo).  General knowledge lives in
radia_mcp.radia_ngsolve.*; this subpackage holds IH-only topics
(induction heating workflow, ESIM cell problem, workpiece SIBC,
Karl iteration, screening physics).
"""

import os
import sys
from pathlib import Path

from mcp.server.fastmcp import FastMCP

from .ih_knowledge import get_induction_heating_documentation
from .sibc_knowledge import get_ih_sibc_documentation

mcp = FastMCP("mcp-server-ih")


# ============================================================
# MCP Tools
# ============================================================

@mcp.tool()
def induction_heating(topic: str = "all") -> str:
    """
    Get induction heating simulation documentation.

    Complete workflow for electromagnetic induction heating analysis:
    EM eddy current solve -> Joule heat computation -> transient
    thermal analysis, including mesh loading, rotating workpiece,
    VTK output, and post-processing patterns.

    Args:
        topic: Documentation topic. Options:
            "all"           - Complete documentation
            "overview"      - Physics overview, parameters, skin depth
            "gmsh_mesh"     - Mesh loading (.vol), physical groups
            "eddy_current"  - A-Phi formulation
            "thermal"       - Transient heat equation
            "rotating"      - Rotating workpiece
            "postprocess"   - VTK output, field evaluation
            "pitfalls"      - Common mistakes
            "esim_kelvin"   - ESIM + Kelvin for IH
    """
    return get_induction_heating_documentation(topic)


@mcp.tool()
def ih_sibc(topic: str = "all") -> str:
    """
    Get IH solver architecture and SIBC documentation.

    Panel pipeline (2026-04-19):
    - PEEC+BEM (1-way)        -> calc_peec_bem.py       (P_wp focus, fast)
    - FEM A-V + wp SIBC + Kelvin -> calc_fem_coilmesh.py (L+P_wp+P_coil, exact)

    Surface Impedance Boundary Condition approaches:
    - SIBC: linear surface impedance (Cu, Al)
    - ESIM: nonlinear Z_s(H) for steel/ferrite (BH curve)

    Args:
        topic: Options:
            "all"         - Complete documentation
            "peec_fem"    - (historical, see ih_knowledge AV_COIL_SIGMA)
            "overview"    - SIBC method selection table
            "esim"        - ESIM cell problem, Karl iteration
            "biot_savart" - Coil field computation (phi_inc, A_inc, H_inc)
            "screening"   - Screening physics, dimensionless parameter
    """
    return get_ih_sibc_documentation(topic)


# ============================================================
# MCP Prompts
# ============================================================

@mcp.prompt()
def new_ih_simulation(geometry: str, material: str = "steel") -> str:
    """Set up a new induction heating simulation."""
    return (
        f"Set up an induction heating simulation for: {geometry}\n"
        f"Material: {material}\n\n"
        "Use the ih_sibc and ih_knowledge tools.\n"
        "Key points (panel pipeline 2026-04-19):\n"
        "1. P_wp fast + rotating WP: PEEC+BEM (1-way) / calc_peec_bem.py\n"
        "2. L + P_wp + P_coil exact: FEM A-V / calc_fem_coilmesh.py\n"
        "3. Both need GAPPED torus (real port terminations)\n"
        "4. For nonlinear steel: ESIM + Karl iteration (FEM path)\n"
        "5. Sample .jou: ih_peec_bem_coarse.jou / ih_fem_kelvin_skin_fine.jou\n"
    )


# ============================================================
# Entry point
# ============================================================

def main():
    if "--selftest" in sys.argv:
        print("IH MCP server self-test:")
        print(f"  ih_sibc('overview'): {len(ih_sibc('overview'))} chars")
        print(f"  induction_heating('overview'): "
              f"{len(induction_heating('overview'))} chars")
        print("  PASSED")
        return

    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
