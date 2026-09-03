"""MCP Server: radia_mcp.bem

Moment Method / Boundary Element Method theory + decision layer.

Covers: MoM foundations (RWG/Galerkin), surface IE (EFIE/MFIE/CFIE/PMCHWT),
low-frequency stabilization (Loop-Star/Calderon/Weggler ST), H-matrix
acceleration (HACApK), and FEM-BEM hybrid.

This is the **theory/genealogy** layer.  For concrete code usage:
- Radia HDiv-VIM: `radia_mcp.radia_ngsolve.hdiv_vim`
- PEEC: `radia_mcp.peec`
- ngsolve.bem: `radia_mcp.radia_ngsolve.ngsbem_inductance`
- HACApK: `radia_mcp.matrix_solvers` (preconditioners_ams)

Usage:
    mcp-server-bem               # stdio
    mcp-server-bem --selftest    # self-test
"""

import sys

from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool

from .overview_knowledge import get_overview_knowledge
from .mom_foundations_knowledge import get_mom_foundations_knowledge
from .surface_ie_knowledge import get_surface_ie_knowledge
from .low_freq_knowledge import get_low_freq_knowledge
from .h_matrix_knowledge import get_h_matrix_knowledge
from .fem_bem_hybrid_knowledge import get_fem_bem_hybrid_knowledge


mcp = FastMCP("mcp-server-bem")


@mcp.tool()
def bem_overview(topic: str = "decision_tree") -> str:
    """
    BEM/MoM landscape: lab stack, decision tree, genealogy.

    ★ DEFAULT topic = 'decision_tree' (which BEM/MoM technique for which problem).

    Args:
        topic: One of:
            "lab_stack"       - Sugahara lab BEM/MoM stack
                                 (HDiv-VIM, ngsolve.bem, HACApK)
            "decision_tree"   - Which technique for which problem (DEFAULT)
            "history"         - Genealogy 1968 Harrington -> 2024
            "all"             - Everything
    """
    return get_overview_knowledge(topic)


@mcp.tool()
def bem_mom_foundations(topic: str = "introduction") -> str:
    """
    MoM foundations: Harrington 1968, RWG 1982, wire-grid (NEC).

    Args:
        topic: One of:
            "introduction"   - The Method in 3 lines (DEFAULT)
            "rwg_basis"      - Rao-Wilton-Glisson 1982 (RWG)
            "wire_grid"      - Legacy NEC-style wire-grid MoM
            "all"            - Everything
    """
    return get_mom_foundations_knowledge(topic)


@mcp.tool()
def bem_surface_ie(topic: str = "efie") -> str:
    """
    Surface Integral Equations: EFIE, MFIE, CFIE, PMCHWT.

    Args:
        topic: One of:
            "efie"      - Electric Field IE (DEFAULT)
            "mfie"      - Magnetic Field IE
            "cfie"      - Combined Field IE (interior resonance free)
            "pmchwt"    - Penetrable dielectric (Poggio-Miller-Chang-...)
            "all"       - Everything
    """
    return get_surface_ie_knowledge(topic)


@mcp.tool()
def bem_low_freq(topic: str = "problem") -> str:
    """
    Low-frequency BEM stabilization.

    ★ Critical for lab IH (kHz), PEEC, ngsolve.bem at DC.

    Args:
        topic: One of:
            "problem"                - Why LF breakdown happens (DEFAULT)
            "loop_star"              - Loop-Star decomposition (PEEC native)
            "calderon"               - Calderon preconditioner (Andriulli 2008)
            "weggler"                - ★ Lucy Weggler's product-space EFIE
                                        stabilization (ngsolve.bem canonical
                                        LF fix; HDivSurface × SurfaceL2 with
                                        κ² V_κ re-scaling).
            "single_trace"           - Alias for "weggler" (legacy name).
            "equivalence_source_lf"  - LF for the equivalence-theorem one-
                                        shot evaluator (NearFieldSource):
                                        Weggler does NOT directly apply
                                        (no matrix to invert); Laplace ML
                                        routing does.  Two-breakdown
                                        distinction.
            "all"                    - Everything
    """
    return get_low_freq_knowledge(topic)


@mcp.tool()
def bem_h_matrix(topic: str = "overview") -> str:
    """
    H-matrix / ACA acceleration for BEM.

    HACApK is a production H-matrix library used by owning HDiv/PEEC/BEM
    operators; it is not selected by a rad.Solve method integer.

    Args:
        topic: One of:
            "overview"   - H-matrix vs FMM family overview (DEFAULT)
            "aca"        - Adaptive Cross Approximation (Bebendorf)
            "hacapk"     - ★ HACApK library (lab production)
            "all"        - Everything
    """
    return get_h_matrix_knowledge(topic)


@mcp.tool()
def bem_fem_bem_hybrid(topic: str = "overview") -> str:
    """
    FEM-BEM hybrid methods for open-boundary EM.

    Note: lab typically prefers Kelvin transform (pure FEM) over FEM-BEM.

    Args:
        topic: One of:
            "overview"   - When to use hybrid vs alternatives (DEFAULT)
            "coupling"   - Iterative vs monolithic coupling
            "all"        - Everything
    """
    return get_fem_bem_hybrid_knowledge(topic)


# ============================================================
# MCP Prompts
# ============================================================

@mcp.prompt()
def pick_a_bem_method(problem_class: str) -> str:
    """Recommend a BEM/MoM technique for a class of problem."""
    guidance = {
        "magnetostatic_pm_iron": (
            "Magnetostatic with PM + nonlinear iron (★ Radia core use case):\n"
            "1. Use HDiv-VIM as the Radia magnetic-material path\n"
            "   -> radia_mcp.radia_ngsolve.hdiv_vim\n"
            "2. Keep reduced-FEM / NGSolve coupling in the same workflow\n"
            "   -> radia_mcp.radia_ngsolve\n"
            "3. Large charge-Gram builds -> HACApK acceleration\n"
            "   -> bem_h_matrix('hacapk')\n"
            "4. Code recipe via radia\n"
            "   -> radia_ngsolve MCP, radia_usage\n"
        ),
        "ih_workpiece_eddy": (
            "Induction heating workpiece (★ lab IH path):\n"
            "1. FEM-SIBC inside workpiece + Kelvin transform for open BC\n"
            "   -> calc_fem_kelvin.py (NOT pure BEM)\n"
            "2. SIBC Robin BC theory\n"
            "   -> radia_mcp.ih (sibc topic)\n"
            "3. If switching to pure surface BEM: use Weggler ST (ngsolve.bem)\n"
            "   -> bem_low_freq('single_trace')\n"
        ),
        "peec_coil_inductance": (
            "PEEC coil / filament inductance:\n"
            "1. Loop-Star is natural in PEEC by construction\n"
            "   -> bem_low_freq('loop_star')\n"
            "2. Coil filaments + panels via PEECBuilder\n"
            "   -> radia_mcp.peec\n"
            "3. Magnetic material coupling should route through HDiv-VIM / reduced FEM\n"
            "   -> radia_mcp.radia_ngsolve.hdiv_vim\n"
        ),
        "high_freq_scattering": (
            "High-frequency PEC/dielectric scattering (>100 MHz):\n"
            "   NOTE: Lab is Laplace-kernel only (MQS/Darwin) per CLAUDE.md\n"
            "1. Use ngsolve.bem with Helmholtz kernel\n"
            "2. EFIE for open surfaces, CFIE for closed\n"
            "   -> bem_surface_ie('efie')\n"
            "   -> bem_surface_ie('cfie')\n"
            "3. Calderon preconditioner for ill-conditioning\n"
            "   -> bem_low_freq('calderon')  (also helps at HF)\n"
            "4. FMM or H-matrix acceleration for large N\n"
            "   -> bem_h_matrix('aca')\n"
        ),
        "transformer_open_boundary": (
            "Transformer with open boundary:\n"
            "1. LAB approach: FEM + Kelvin transform (NOT FEM-BEM hybrid)\n"
            "   -> radia_mcp.radia_ngsolve.kelvin_transformation\n"
            "2. Use HDiv-VIM / reduced-FEM coupling for magnetic materials\n"
            "   -> radia_mcp.radia_ngsolve.hdiv_vim\n"
            "3. Hybrid FEM-BEM only as last resort\n"
            "   -> bem_fem_bem_hybrid('coupling')\n"
        ),
    }
    return guidance.get(problem_class, (
        f"Unknown problem_class '{problem_class}'. Available:\n"
        "  magnetostatic_pm_iron, ih_workpiece_eddy, peec_coil_inductance,\n"
        "  high_freq_scattering, transformer_open_boundary.\n"
        "For overview, see bem_overview('decision_tree').\n"
    ))




register_status_tool(
    mcp,
    server_name='mcp-server-bem',
    description='MoM/BEM theory: RWG, EFIE/MFIE/CFIE/PMCHWT, Loop-Star, Calderon, HACApK, FEM-BEM',
    subpackage='radia_mcp.bem',
    related_servers=["radia-ngsolve", "peec"],
)


def main():
    """Entry point for mcp-server-bem."""
    if "--selftest" in sys.argv:
        print("bem MCP server self-test:")
        print(f"  overview: {len(get_overview_knowledge('all'))} chars")
        print(f"  mom_foundations: {len(get_mom_foundations_knowledge('all'))} chars")
        print(f"  surface_ie: {len(get_surface_ie_knowledge('all'))} chars")
        print(f"  low_freq: {len(get_low_freq_knowledge('all'))} chars")
        print(f"  h_matrix: {len(get_h_matrix_knowledge('all'))} chars")
        print(f"  fem_bem_hybrid: {len(get_fem_bem_hybrid_knowledge('all'))} chars")
        print("OK")
        return
    mcp.run()


if __name__ == "__main__":
    main()
