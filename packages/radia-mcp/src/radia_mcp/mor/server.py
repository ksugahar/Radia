"""
MOR (Model Order Reduction) MCP Server (radia_mcp.mor)

Knowledge layer for eddy-current FEM model order reduction, centered
on the **Cauer Ladder Network (CLN)** method — a LAB SPECIALTY of
菅原研 (Sugahara Lab) where Sugahara is co-author on the canonical
foundational papers (2018-2020).

Usage:
    mcp-server-mor              # Start MCP server (stdio)
    mcp-server-mor --selftest   # Self-test
"""

import sys

from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool

from .cln_knowledge import get_cln_documentation
from .systematic_knowledge import get_systematic_mor_knowledge

try:
    from .bibliography_index_knowledge import get_bibliography_index
except ImportError:
    def get_bibliography_index(query: str = "") -> str:
        return "MOR bibliography index not yet generated."


mcp = FastMCP("mcp-server-mor")


@mcp.tool()
def mor_cln(topic: str = "all") -> str:
    """Cauer Ladder Network (CLN) -- the lab-specialty MOR for eddy-current FEM.

    Topics:
      "all"
      "overview"     - CLN big picture, Sugahara Lab's role
      "recursion"    - Basic A-formulation recursion + Cauer circuit
      "multiple"     - Multiple expansion points (Kuriyama 2019)
      "nonlinear"    - Nonlinear ferromagnetic extension
      "applications" - Industrial inductors, WPT, hybrid twin
    """
    return get_cln_documentation(topic)


@mcp.tool()
def mor_systematic(topic: str = "mor_taxonomy") -> str:
    """
    Systematic MOR knowledge -- distilled from the deGruyter 3-volume
    MOR Handbook (1221 pages) + Kiss-Orosz 2024 Energies rotating-
    machines review + Ioan Vol3 Ch5 EM-specific MOR.

    Args:
        topic: One of:
            "mor_taxonomy"          - Three families: projection / system / data
            "projection_pod_rb"     - POD + Reduced Basis (Vol 2 Ch 2, 4)
            "projection_krylov"     - Krylov / Arnoldi / Lanczos (Vol 1 Ch 3)
            "pgd"                   - Proper Generalized Decomposition (Vol 2 Ch 3)
            "hyperreduction"        - DEIM / EIM for nonlinear (Vol 2 Ch 5)
            "system_theoretic_bt"   - Balanced Truncation, H2/H_inf (Vol 1 Ch 2)
            "data_driven_dmd_oi"    - DMD, OI, Loewner (Vol 1 Ch 6 + Vol 2 Ch 7)
            "parametric_pmor"       - Parametric MOR (Vol 2 Ch 1)
            "cln_sugahara"          - CLN positioned in MOR taxonomy
            "em_specific_ioan"      - Ioan Vol 3 Ch 5 (56p EM-specific)
            "rotating_machines_kiss_orosz_2024" - 2024 review for motors
            "software_lab"          - pyMOR / MORLab / lab tools
            "lab_recommendation"    - Decision guide for radia + NGSolve
            "all"                   - Everything
    """
    return get_systematic_mor_knowledge(topic)


@mcp.tool()
def mor_bibliography(query: str = "") -> str:
    """Search the MOR bibliography catalog (87 papers in lab library)."""
    return get_bibliography_index(query)




register_status_tool(
    mcp,
    server_name='mcp-server-mor',
    description='Model Order Reduction: PRIMA, Cauer Ladder Network, hyperreduction (DEIM)',
    subpackage='radia_mcp.mor',
    related_servers=["radia-ngsolve", "rna-mec"],
)


def main():
    if "--selftest" in sys.argv:
        print("MOR MCP server self-test:")
        print(f"  cln('all'): {len(mor_cln('all'))} chars")
        from .systematic_knowledge import SECTIONS as S
        for k in S:
            r = mor_systematic(k)
            assert len(r) > 500, f"{k} short ({len(r)})"
            print(f"  systematic({k!r}): {len(r)} chars")
        print("  PASSED")
        return
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
