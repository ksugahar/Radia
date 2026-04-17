"""
Induction Heating MCP Server

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

    Two solver paths (v4.6.0):
    - PEEC+FEM (default): coil PEEC filaments + FEM-SIBC+Kelvin workpiece
    - ALL FEM (reference): full volume FEM with Kelvin

    Surface Impedance Boundary Condition approaches:
    - SIBC: linear surface impedance (Cu, Al)
    - ESIM: nonlinear Z_s(H) for steel/ferrite (BH curve)

    Args:
        topic: Options:
            "all"         - Complete documentation
            "peec_fem"    - PEEC+FEM vs ALL FEM architecture (START HERE)
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
        "Use the ih_sibc tool for SIBC method selection.\n"
        "Key points:\n"
        "1. For P_total: use BEM (ScalarBIESIBCSolver) - calc_heating_bem.py\n"
        "2. For L/B: use FEM (HCurl + Kelvin) - calc_fem_kelvin.py\n"
        "3. For nonlinear steel: use ESIM + Karl iteration\n"
        "4. Coil = Biot-Savart filament (no coil mesh needed for BEM)\n"
        "5. Workpiece surface mesh from OCC or Cubit .vol export\n"
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
