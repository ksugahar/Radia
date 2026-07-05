"""
Induction Heating MCP Server (radia_mcp.ih)

Provides tools for:
- IH-specific linting (eddy current, SIBC, ESIM, Bessel functions)
- SIBC method selection (BEM Scalar BIE vs FEM coilmesh)
- ESIM nonlinear surface impedance documentation
- Biot-Savart coil field computation
- Screening physics and Karl iteration
- IH simulation workflow (EM -> thermal)

Usage:
    mcp-server-ih              # Start MCP server (stdio transport)
    mcp-server-ih --selftest   # Run self-test

Promoted from legacy private source tree to public radia-mcp on
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
from .esim_knowledge import get_ih_esim_documentation
from ..common import register_status_tool

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

    Panel pipeline (v4.25.0 unified Layer 4):
    - PEEC+BEM (weak coupling)   -> calc_inductance.py --coil-solver peec --vol <wp>  (Telegen ΔL + P_wp)
    - BEM-A+BEM (weak coupling)  -> calc_inductance.py --coil-solver bem-a --vol <wp> (n_peri-free coil)
    - PEEC + FEM wp + Kelvin     -> calc_fem_kelvin.py --formulation total --peec-step ...
    - Full FEM A-V + wp SIBC + Kelvin -> calc_fem_coilmesh.py    (L+P_wp+P_coil, exact)

    Surface Impedance Boundary Condition approaches:
    - SIBC: linear surface impedance (Cu, Al)
    - ESIM: nonlinear Z_s(H) for steel/ferrite (BH curve)

    Args:
        topic: Options:
            "all"                       - Complete documentation
            "peec_fem"                  - (historical, see ih_knowledge AV_COIL_SIGMA)
            "esim"                      - ESIM cell problem, Karl iteration
            "screening"                 - Screening physics, dimensionless parameter
            "mathematica_verification"  - Wolfram Language recipes
                                          (Leontovich Z_s derivation, Smythe
                                          sphere reference, per-panel Mitzner
                                          check) for symbolic side-by-side
                                          validation of BEM/FEM SIBC results.
            "ngsolve_recipes"           - Production NGSolve code patterns:
                                          Scalar BIE+SIBC with COCR+ComplexCompactAMS,
                                          per-node Z_s from per-panel curvature,
                                          FEM+Kelvin+Robin hole approach,
                                          Verify-First Policy FES checks.
    """
    return get_ih_sibc_documentation(topic)


@mcp.tool()
def ih_esim(topic: str = "all") -> str:
    """
    Induction-heating ESIM (Effective Surface Impedance Method) usage.

    Practical how-to-use guide for the nonlinear ESIM cell-problem
    solver in the IH workflow.  Sibling to `ih_sibc` (which covers
    linear SIBC + general IH-solver architecture).

    Covers the three production CLIs (calc_inductance.py,
    calc_fem_kelvin.py, calc_fem_coilmesh.py): CLI flag reference,
    BH-curve file format, scalar-vs-per-element decision, Karl-
    iteration tuning, JSON output schema, and troubleshooting.

    For the underlying physics + formulation, see the canonical docs
    at repo:/docs/esim/ (especially MATHEMATICAL_ANALYSIS.md
    and IMPLEMENTATION.md) or the sibling tool `esim` in
    radia_ngsolve server (physics-side overview).

    Args:
        topic: Section to return.  Options:
            "all"                  - All sections concatenated
            "overview"             - When to use ESIM vs linear SIBC
                                     (decision table, cost vs benefit)
            "em_table_coupling"    - sigma(T)/mu(T) thermal-EM coupling
                                     via a precomputed Z_s(|H_t|,T) table
                                     (calc_em_table.py + calc_heat_with_
                                     em_table.py; kelvin vs biot ht-source)
            "bh_file"              - BH-curve file format spec
            "inductance_cli"       - calc_inductance.py CLI flags + example
            "fem_kelvin_cli"       - calc_fem_kelvin.py CLI flags + example
                                     (note inconsistent flag names!)
            "fem_coilmesh_cli"     - calc_fem_coilmesh.py CLI flags + example
            "per_element"          - --esim-per-panel: when it pays off,
                                     convergence caveats, backend restrictions
            "convergence"          - --esim-relax tuning, reading esim_history
            "json_output"          - Output JSON schema (esim_history, etc.)
            "troubleshooting"      - Common errors + sanity checks
            "scalar_vs_vector_bem" - Why scalar BIE-SIBC + curved Tri6 is
                                     the right combination (vs vector BEM-A
                                     / FEM-Kelvin / FEM-coilmesh).  Decision
                                     matrix, error-rate match, Karl per-iter
                                     cost comparison, reviewer Q&A.
                                     (IGTE 2026 paper marketing line.)
            "headline_numbers"     - Locked-in numerical results for paper
                                     citation: dense 108-case sweep,
                                     per-DOF vs uniform headline,
                                     three-path consistency table,
                                     Bessel cell validation table.
    """
    return get_ih_esim_documentation(topic)


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
        "Key points (panel pipeline, v4.25.0+ unified Layer 4):\n"
        "1. L + P_wp via weak coupling (PEEC coil, Telegen ΔL):\n"
        "     calc_inductance.py --coil-solver peec --coil-step <coil.step>\n"
        "                         --vol <wp.vol>\n"
        "2. L + P_wp + P_coil exact: calc_fem_coilmesh.py (FEM A-V)\n"
        "3. All paths assume GAPPED torus (real port terminations)\n"
        "4. For nonlinear steel: ESIM + Karl iteration\n"
        "   (--impedance-model esim --bh-file ...)\n"
        "5. Sample .jou: ih_peec_bem_coarse.jou / ih_fem_kelvin_skin_fine.jou\n"
    )


# ============================================================
# Entry point
# ============================================================

register_status_tool(
    mcp,
    server_name="mcp-server-ih",
    description=("Induction heating: SIBC + ESIM + Karl iteration "
                  "+ workpiece coupling. Production calc scripts: "
                  "calc_inductance.py / calc_fem_kelvin.py / calc_fem_coilmesh.py."),
    subpackage="radia_mcp.ih",
    related_servers=["radia-peec", "radia-magnetic-materials",
                      "radia-radia-ngsolve", "radia-litz-transmission"],
    optional_deps=["radia", "ngsolve"],
)


def main():
    if "--selftest" in sys.argv:
        print("IH MCP server self-test:")
        print(f"  ih_sibc('overview'): {len(ih_sibc('overview'))} chars")
        print(f"  induction_heating('overview'): "
              f"{len(induction_heating('overview'))} chars")
        print(f"  ih_esim('overview'): "
              f"{len(ih_esim('overview'))} chars")
        print("  PASSED")
        return

    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
