"""MCP Server: radia_mcp.fem

FEM formulations theory + decision layer for EM analysis.

Covers all major potential formulations (A-Omega, T-Omega, H, Reduced,
Darwin), element technology (edge, high-order, XFEM, isogeometric, DG),
gauging (tree-cotree, ungauged + AMS), open boundary (Kelvin, ABC, PML),
domain-specific (axisym Henrotte, TD-FEM, harmonic balance, HF),
circuit coupling, and large-scale + multi-scale (Hollaus MSFEM).

This is the **theory/genealogy** layer.  For code usage:
- NGSolve API: see `radia_mcp.radia_ngsolve`
- Axisymmetric: see `radia_mcp.radia_ngsolve.axifemm_documentation`
- Hollaus lamination: see `radia_mcp.motor.hollaus_eddy`
- Solver+preconditioner: see `radia_mcp.matrix_solvers`
- BEM/MoM: see `radia_mcp.bem`

Usage:
    mcp-server-fem              # stdio
    mcp-server-fem --selftest   # self-test
"""

import sys

from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool

from .overview_knowledge import get_overview_knowledge
from .potential_formulations_knowledge import get_potential_formulations_knowledge
from .elements_knowledge import get_elements_knowledge
from .gauge_open_boundary_knowledge import get_gauge_open_boundary_knowledge
from .time_domain_axisym_knowledge import get_time_domain_axisym_knowledge
from .large_scale_special_knowledge import get_large_scale_special_knowledge
from .ngsolve_hierarchy_knowledge import get_ngsolve_hierarchy_knowledge


mcp = FastMCP("mcp-server-fem")


@mcp.tool()
def fem_overview(topic: str = "decision_tree") -> str:
    """
    FEM landscape: lab stack, decision tree, genealogy.

    Args:
        topic: One of:
            "lab_stack"      - Production NGSolve + sparsesolv + Radia stack
            "decision_tree"  - Which formulation for which problem (DEFAULT)
            "history"        - Genealogy 1980 Nedelec -> 2025
            "all"            - Everything
    """
    return get_overview_knowledge(topic)


@mcp.tool()
def fem_potential_formulations(topic: str = "catalog") -> str:
    """
    Potential formulations: A-Omega, T-Omega, H, Reduced, Darwin.

    Args:
        topic: One of:
            "catalog"        - All formulations + decision tree (DEFAULT)
            "a_omega"        - A-Omega (mixed vector-scalar)
            "t_omega"        - T-Omega (electric T + magnetic Omega)
            "h_formulation"  - H-formulation (HCurl direct, SC / nonlinear)
            "reduced"        - ★ Reduced potential (accelerator magnet)
            "darwin"         - Darwin model (DC to MHz transition)
            "all"            - Everything
    """
    return get_potential_formulations_knowledge(topic)


@mcp.tool()
def fem_elements(topic: str = "catalog") -> str:
    """
    Element technology: edge (Nedelec), high-order, XFEM, isogeometric, DG.

    Args:
        topic: One of:
            "catalog"        - H1/HCurl/HDiv/L2 + hierarchical (DEFAULT)
            "edge"           - Edge elements (Nedelec) - foundation of HCurl
            "high_order"     - p-version (hierarchical basis)
            "xfem"           - eXtended FEM (cracks, interfaces)
            "isogeometric"   - Isogeometric Analysis (NURBS)
            "dg"             - Discontinuous Galerkin
            "all"            - Everything
    """
    return get_elements_knowledge(topic)


@mcp.tool()
def fem_gauge_open_boundary(topic: str = "gauge") -> str:
    """
    Gauging + open boundary techniques.

    Args:
        topic: One of:
            "gauge"             - Gauge choices for HCurl A (DEFAULT)
            "kelvin_transform"  - ★ Kelvin transformation (lab default)
            "abc_pml"           - Asymptotic BC and PML alternatives
            "all"               - Everything
    """
    return get_gauge_open_boundary_knowledge(topic)


@mcp.tool()
def fem_time_domain_axisym(topic: str = "henrotte_axisym") -> str:
    """
    Time-domain, axisymmetric (Henrotte), harmonic balance, HF, circuit coupling.

    Args:
        topic: One of:
            "time_domain_fem"    - TD-FEM time-stepping
            "henrotte_axisym"    - ★ Henrotte basis for axisym (lab CORE)
            "harmonic_balance"   - HB method for periodic steady-state
            "high_frequency"     - HF Helmholtz (reference, out of lab scope)
            "circuit_coupling"   - FE-circuit coupling (Hanser 2025 Schur)
            "all"                - Everything
    """
    return get_time_domain_axisym_knowledge(topic)


@mcp.tool()
def fem_ngsolve_hierarchy(topic: str = "decision_tree") -> str:
    """
    NGSolve hierarchical H(curl) bases - Zaglmayr / nograds / tree-cotree.

    Migrated 2026-05-23 from mcp_server_document/mathematica/notes_fem_hcurl.md
    per "all FEM knowledge lives in radia-mcp" policy.

    Args:
        topic: One of:
            "decision_tree"            - Which strategy for which question (DEFAULT)
            "zaglmayr_split"           - Hierarchical gradient + rotational split
            "nograds_option"           - NGSolve `nograds=True` behavior
            "tree_cotree"              - Tree-cotree gauge + generalizations
            "mathematica_verification" - Symbolic verification recipes (Wolfram)
            "all"                      - Everything
    """
    return get_ngsolve_hierarchy_knowledge(topic)


@mcp.tool()
def fem_large_scale_special(topic: str = "large_scale") -> str:
    """
    Large-scale, error theory, multi-scale (Hollaus MSFEM), misc techniques.

    Args:
        topic: One of:
            "large_scale"       - Parallel scaling, AMS/AMG (DEFAULT)
            "error_theory"      - A-posteriori estimators, adaptivity
            "hollaus_msfem"     - ★ Multi-scale FEM for laminated cores
            "misc_techniques"   - Fixed-point, polarization, Helmholtz-Weyl
            "all"               - Everything
    """
    return get_large_scale_special_knowledge(topic)


# ============================================================
# MCP Prompts
# ============================================================

@mcp.prompt()
def pick_a_fem_formulation(problem_class: str) -> str:
    """Recommend FEM formulation for a class of EM problem."""
    guidance = {
        "accelerator_magnet": (
            "Accelerator electromagnet (★ lab core use case):\n"
            "1. Coil source: Biot-Savart from CoilBuilder (no coil mesh!)\n"
            "   → fem_potential_formulations('reduced')\n"
            "2. Iron yoke: Reduced Omega + nonlinear iron\n"
            "   → fem_potential_formulations('reduced')\n"
            "3. Open boundary: Kelvin transform\n"
            "   → fem_gauge_open_boundary('kelvin_transform')\n"
            "4. Solver: CG + CompactAMS (HCurl) or direct LU (small)\n"
            "   → matrix_solvers MCP\n"
            "5. Code recipe via radia\n"
            "   → radia_ngsolve.radia_usage, electromagnet MCP\n"
        ),
        "induction_heating_workpiece": (
            "IH workpiece (★ lab IH path):\n"
            "1. Workpiece: bulk eddy current with SIBC\n"
            "   → fem_potential_formulations('a_omega')\n"
            "   → radia_mcp.ih.sibc\n"
            "2. Coil: PEEC filaments (not FEM)\n"
            "   → radia_mcp.peec\n"
            "3. Open boundary: Kelvin\n"
            "   → fem_gauge_open_boundary('kelvin_transform')\n"
            "4. Frequency-domain: complex symmetric → COCR + ComplexCompactAMS\n"
            "   → matrix_solvers_krylov('cocg_cocr')\n"
            "   → matrix_solvers_preconditioners('ams_hiptmair_xu')\n"
            "5. Code: calc_fem_kelvin.py\n"
        ),
        "motor_transient": (
            "Motor transient (PMSM/IM/SRM):\n"
            "1. Vector A with HCurl + Picard nonlinear (iron BH)\n"
            "   → fem_potential_formulations('a_omega')\n"
            "2. Time-stepping: implicit Euler / Crank-Nicolson\n"
            "   → fem_time_domain_axisym('time_domain_fem')\n"
            "3. Rotor rotation: rotor PM via Br(theta) CoefficientFunction\n"
            "   → radia_mcp.motor\n"
            "4. Circuit coupling: Hanser 2025 Schur voltage drive\n"
            "   → fem_time_domain_axisym('circuit_coupling')\n"
            "5. Laminated stator: Hollaus MSFEM effective material\n"
            "   → fem_large_scale_special('hollaus_msfem')\n"
            "6. Code: calc_motor_transient.py, calc_motor_lamination.py\n"
        ),
        "axisymmetric_magnetic": (
            "Axisymmetric magnetic A_phi (★ Henrotte CORE):\n"
            "1. Use Henrotte basis {1, r^2, z}\n"
            "   → fem_time_domain_axisym('henrotte_axisym')\n"
            "2. radia.radia_axifemm.H1Henrotte\n"
            "3. Standard H1 for scalar T, V, etc. (NOT Henrotte)\n"
            "4. Code: radia_ngsolve.axifemm_documentation\n"
        ),
        "open_boundary_static": (
            "Open boundary magnetostatic (★ lab default):\n"
            "1. Kelvin transform (radia + NGSolve, pure FEM)\n"
            "   → fem_gauge_open_boundary('kelvin_transform')\n"
            "2. NOT FEM-BEM hybrid (extra complexity)\n"
            "3. NOT PML (fails at DC)\n"
            "4. Code: calc_fem_kelvin.py, kelvin_transformation MCP\n"
        ),
        "eddy_current_mqs": (
            "Eddy current MQS (frequency domain):\n"
            "1. Bulk conductor, multiply-connected → A-V Biro-Preis\n"
            "   → fem_potential_formulations('catalog')\n"
            "2. Surface skin only → A + SIBC Robin BC (lab IH)\n"
            "   → radia_mcp.ih.sibc\n"
            "3. Solver: COCR + ComplexCompactAMS\n"
            "   → matrix_solvers_krylov('cocg_cocr')\n"
            "4. Low-freq stabilization: shifted preconditioner\n"
            "   → matrix_solvers_preconditioners('shifted')\n"
        ),
    }
    return guidance.get(problem_class, (
        f"Unknown problem_class '{problem_class}'. Available:\n"
        "  accelerator_magnet, induction_heating_workpiece, motor_transient,\n"
        "  axisymmetric_magnetic, open_boundary_static, eddy_current_mqs.\n"
        "For overview, see fem_overview('decision_tree').\n"
    ))




register_status_tool(
    mcp,
    server_name='mcp-server-fem',
    description='FEM formulations theory layer (A-Omega / T-Omega / H / Reduced / Darwin, edge / HO / XFEM / IGA / DG, gauging + Kelvin, MSFEM, Schur...',
    subpackage='radia_mcp.fem',
    related_servers=["radia-ngsolve", "bem", "differential-forms"],
)


def main():
    """Entry point for mcp-server-fem."""
    if "--selftest" in sys.argv:
        print("fem MCP server self-test:")
        print(f"  overview: {len(get_overview_knowledge('all'))} chars")
        print(f"  potential_formulations: {len(get_potential_formulations_knowledge('all'))} chars")
        print(f"  elements: {len(get_elements_knowledge('all'))} chars")
        print(f"  gauge_open_boundary: {len(get_gauge_open_boundary_knowledge('all'))} chars")
        print(f"  time_domain_axisym: {len(get_time_domain_axisym_knowledge('all'))} chars")
        print(f"  large_scale_special: {len(get_large_scale_special_knowledge('all'))} chars")
        print(f"  ngsolve_hierarchy: {len(get_ngsolve_hierarchy_knowledge('all'))} chars")
        print("OK")
        return
    mcp.run()


if __name__ == "__main__":
    main()
