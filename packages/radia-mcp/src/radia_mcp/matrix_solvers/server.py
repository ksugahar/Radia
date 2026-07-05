"""MCP Server: radia_mcp.matrix_solvers

Sparse linear solver theory + decision layer for EM analysis.

Covers direct solvers (PARDISO/MUMPS/LU), Krylov methods (CG/BiCGSTAB/
GMRES/COCG/COCR/IDR(s)), preconditioners (Jacobi/IC/MIC/ILU/AMG/AMS
Hiptmair-Xu), and EM-specific topics (tree-cotree gauging, Biro-Preis
A-V, shifted preconditioner for air+conductor).

This is the **theory/genealogy** layer.  For concrete code usage of the
lab stack (CompactAMS, COCR, etc.), see the `sparsesolv` tool in
`radia_mcp.radia_ngsolve`.

Usage:
    mcp-server-matrix-solvers              # stdio
    mcp-server-matrix-solvers --selftest   # self-test
"""

import sys

from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool

from .overview_knowledge import get_overview_knowledge
from .direct_solvers_knowledge import get_direct_solvers_knowledge
from .krylov_knowledge import get_krylov_knowledge
from .preconditioners_knowledge import get_preconditioners_knowledge
from .em_specific_knowledge import get_em_specific_knowledge


mcp = FastMCP("mcp-server-matrix-solvers")


@mcp.tool()
def matrix_solvers_overview(topic: str = "decision_tree") -> str:
    """
    Solver landscape: lab stack, decision tree, genealogy.

    ★ DEFAULT topic = 'decision_tree' (which solver+preconditioner for
    which problem class).

    Args:
        topic: One of:
            "lab_stack"       - radia.sparsesolv_ngsolve components +
                                 which to pick (★ CompactAMS / COCR)
            "decision_tree"   - Solver+preconditioner decision (DEFAULT)
            "history"         - Genealogy 1952 CG → 2007 Hiptmair-Xu
            "all"             - Everything
    """
    return get_overview_knowledge(topic)


@mcp.tool()
def matrix_solvers_direct(topic: str = "overview") -> str:
    """
    Direct sparse solvers: LU, PARDISO, MUMPS, SuperLU.

    Args:
        topic: One of:
            "overview"      - When direct beats iterative (DEFAULT)
            "pardiso"       - Intel MKL PARDISO (NGSolve default)
            "mumps"         - Distributed MUMPS for >1 node
            "lu_radia"      - Built-in LU for Radia HDiv-VIM (N<500)
            "all"           - Everything
    """
    return get_direct_solvers_knowledge(topic)


@mcp.tool()
def matrix_solvers_krylov(topic: str = "catalog") -> str:
    """
    Krylov subspace methods: CG, BiCGSTAB, GMRES, COCG, COCR, IDR(s).

    ★ Lab core: COCR (Sogabe-Zhang 2007) for complex symmetric eddy
    current — see `cocg_cocr` topic.

    Args:
        topic: One of:
            "catalog"           - All Krylov methods at a glance (DEFAULT)
            "cg"                - Hestenes-Stiefel 1952 CG (★ lab CGSolver)
            "gmres"             - Saad-Schultz 1986 GMRES
            "bicgstab"          - vdVorst 1992 BiCGSTAB
            "cocg_cocr"         - Complex symmetric: COCG (vdVorst 1990)
                                   and COCR (Sogabe-Zhang 2007) ★ lab CORE
            "idr"               - Sonneveld-vGijzen 2009 IDR(s)
            "all"               - Everything
    """
    return get_krylov_knowledge(topic)


@mcp.tool()
def matrix_solvers_preconditioners(topic: str = "catalog") -> str:
    """
    Preconditioner catalog: classical, AMG, AMS (Hiptmair-Xu).

    ★ Lab core: CompactAMS (Hiptmair-Xu 2007) for HCurl problems —
    see `ams_hiptmair_xu` topic.

    Args:
        topic: One of:
            "catalog"            - Family overview (DEFAULT)
            "classical"          - Jacobi/SOR/IC/MIC/ILU
                                   (Meijerink-vdVorst 1977 family)
            "amg"                - Algebraic multigrid (Ruge-Stuben,
                                   BoomerAMG, CompactAMG)
            "ams_hiptmair_xu"    - ★ AMS (Hiptmair-Xu 2007) = lab CompactAMS
            "shifted"            - Shifted preconditioner for air+conductor
                                   HCurl (lab production config)
            "all"                - Everything
    """
    return get_preconditioners_knowledge(topic)


@mcp.tool()
def matrix_solvers_em_specific(topic: str = "gauging") -> str:
    """
    EM-specific solver / formulation choices.

    Tree-cotree gauging (Albanese-Rubinacci, Ostrowski-Hiptmair),
    Biro-Preis A-V formulation, eddy-current low-frequency
    stabilization.

    Args:
        topic: One of:
            "gauging"            - Gauge overview + when to use which (DEFAULT)
            "tree_cotree"        - Albanese-Rubinacci + Two-Step Maxwell
                                    (Ostrowski-Hiptmair 2021)
            "biro_preis_av"      - A-V reduced potential (multiply connected)
            "eddy_current_stab"  - Low-freq + air-region stabilization
                                    (lab production config)
            "all"                - Everything
    """
    return get_em_specific_knowledge(topic)


# ============================================================
# MCP Prompts
# ============================================================

@mcp.prompt()
def pick_a_solver(problem_class: str) -> str:
    """Recommend a solver+preconditioner combination for a class of EM problem."""
    guidance = {
        "magnetostatic_h1": (
            "Magnetostatic (H1 scalar potential):\n"
            "1. SPD matrix → CG\n"
            "   → matrix_solvers_krylov('cg')\n"
            "2. Preconditioner: classical IC for N<100k, AMG for larger\n"
            "   → matrix_solvers_preconditioners('classical')\n"
            "   → matrix_solvers_preconditioners('amg')\n"
            "3. For Radia HDiv-VIM (N<500): LU direct (method=0)\n"
            "   → matrix_solvers_direct('lu_radia')\n"
        ),
        "magnetostatic_hcurl": (
            "Magnetostatic (HCurl A-formulation):\n"
            "1. Gauge choice — default: ungauged + CompactAMS\n"
            "   → matrix_solvers_em_specific('gauging')\n"
            "2. Solver: CG (real SPD)\n"
            "   → matrix_solvers_krylov('cg')\n"
            "3. Preconditioner: ★ CompactAMS (Hiptmair-Xu)\n"
            "   → matrix_solvers_preconditioners('ams_hiptmair_xu')\n"
            "4. Code recipe via radia.sparsesolv_ngsolve\n"
            "   → radia_ngsolve MCP, tool sparsesolv('compact_ams')\n"
        ),
        "eddy_current_mqs": (
            "Eddy current MQS (complex symmetric HCurl):\n"
            "1. Solver: ★ COCR (Sogabe-Zhang 2007)\n"
            "   → matrix_solvers_krylov('cocg_cocr')\n"
            "2. Preconditioner: ★ ComplexCompactAMS\n"
            "   → matrix_solvers_preconditioners('ams_hiptmair_xu')\n"
            "3. Air-region stabilization: shifted preconditioner\n"
            "   → matrix_solvers_preconditioners('shifted')\n"
            "   → matrix_solvers_em_specific('eddy_current_stab')\n"
            "4. Code recipe via radia.sparsesolv_ngsolve\n"
            "   → radia_ngsolve MCP, tool sparsesolv('example_compact_ams')\n"
        ),
        "frequency_sweep": (
            "Frequency sweep (same A factored once, many RHS):\n"
            "1. Direct solver PARDISO — amortizes factor cost\n"
            "   → matrix_solvers_direct('pardiso')\n"
            "2. For very large N (>500k): iterative with recycling Krylov\n"
            "   (not yet in lab stack — defer to NGSolve `MGmm` examples)\n"
            "3. For S-parameter extraction in PEEC: PRIMA Lanczos MOR\n"
            "   → radia_mcp.mor MCP\n"
        ),
        "non_symmetric": (
            "Non-symmetric real (advection-diffusion, etc.):\n"
            "1. Mild non-symmetry → BiCGSTAB + ILU\n"
            "   → matrix_solvers_krylov('bicgstab')\n"
            "   → matrix_solvers_preconditioners('classical')\n"
            "2. Severe + restart-OK → GMRES(m) + AMG\n"
            "   → matrix_solvers_krylov('gmres')\n"
            "3. Hard cases → IDR(s) (not shipped in lab stack)\n"
            "   → matrix_solvers_krylov('idr')\n"
        ),
        "indefinite_saddle_point": (
            "Indefinite saddle-point (Stokes-like, mixed FE):\n"
            "1. MINRES + block-diagonal preconditioner\n"
            "2. Direct: MUMPS (handles symmetric indefinite cleanly)\n"
            "   → matrix_solvers_direct('mumps')\n"
            "3. NOT in lab core (defer to NGSolve examples)\n"
        ),
    }
    return guidance.get(problem_class, (
        f"Unknown problem_class '{problem_class}'. Available:\n"
        "  magnetostatic_h1, magnetostatic_hcurl, eddy_current_mqs,\n"
        "  frequency_sweep, non_symmetric, indefinite_saddle_point.\n"
        "For overview, see matrix_solvers_overview('decision_tree').\n"
    ))




register_status_tool(
    mcp,
    server_name='mcp-server-matrix-solvers',
    description='Sparse solver theory + decision tree: Krylov (CG/BiCGSTAB/GMRES/COCG/COCR/IDR), preconditioners (AMG, Hiptmair-Xu AMS), Biro-Preis A-V,...',
    subpackage='radia_mcp.matrix_solvers',
    related_servers=["radia-ngsolve", "fem"],
)


def main():
    """Entry point for mcp-server-matrix-solvers."""
    if "--selftest" in sys.argv:
        print("matrix_solvers MCP server self-test:")
        print(f"  overview: {len(get_overview_knowledge('all'))} chars")
        print(f"  direct_solvers: {len(get_direct_solvers_knowledge('all'))} chars")
        print(f"  krylov: {len(get_krylov_knowledge('all'))} chars")
        print(f"  preconditioners: {len(get_preconditioners_knowledge('all'))} chars")
        print(f"  em_specific: {len(get_em_specific_knowledge('all'))} chars")
        print("OK")
        return
    mcp.run()


if __name__ == "__main__":
    main()
