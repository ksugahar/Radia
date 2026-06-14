"""Motor MCP Server (radia_mcp.motor)

Rotating electric machine analysis: 2D magnetodynamic A-formulation,
moving-band air-gap coupling, motor-type vocabulary (PMSM/IM/SRM/WFSM),
SynRM topology optimization (Wakao 2025 autoencoder + LS), and the
Darwin model time-domain solver (Kaimori-Mifune-Kameari-Wakao 2024).

Distilled from:
- S:/ONELAB/ElectricMachines/   (Sabariego-Gyselinck-Geuzaine)
- W:/04_卒論論文関係/2025年度/136_劉馨遙/   (Liu Xinyao thesis)

Usage:
    mcp-server-motor              # Start MCP server (stdio transport)
    mcp-server-motor --selftest   # Run self-test
"""

import sys

from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool

from .onelab_knowledge import get_onelab_knowledge
from .topology_opt_knowledge import get_topology_opt_knowledge
from .darwin_model_knowledge import get_darwin_knowledge
from .femm_transient_knowledge import get_femm_transient_knowledge
from .henrotte_lineage_knowledge import get_henrotte_lineage_knowledge
from .hollaus_eddy_knowledge import get_hollaus_eddy_knowledge
from .hollaus_genealogy_knowledge import get_hollaus_genealogy
from .tritool_cross_reference_knowledge import get_tritool_cross_reference

try:
    from .bibliography_index_knowledge import get_bibliography_index
except ImportError:
    def get_bibliography_index(query: str = "") -> str:
        return "Motor bibliography index not yet generated."


mcp = FastMCP("mcp-server-motor")


# ============================================================
# MCP Tools
# ============================================================

@mcp.tool()
def motor_onelab(topic: str = "overview") -> str:
    """
    ONELAB/GetDP electric-machine reference template knowledge.

    Args:
        topic: One of:
            "overview"        - Bundle layout + reference machine table
            "groups"          - Region / FunctionSpace / Constraint vocab
            "formulation"     - 2D A-formulation weak form + non-linear ν
            "magstadyn_source"- VERIFIED weak form/torque/slip/Park from .pro source
            "twod_corrections"- 2D→3D corrections: end-effects/skew/lamination/etc
            "motor_types"     - PMSM, IM, SRM, WFSM, shaded-pole table
            "ngsolve_xlate"   - GetDP → NGSolve translation patterns
            "analysis_modes"  - Static / Time-Domain / Frequency-Domain
            "circuit"         - External-circuit coupling (MNA)
            "post"            - Torque (Arkkio), iron-loss, FFT recipes
            "liu_thesis"      - 2025 SynRM thesis: ONELAB vs an independent ref
            "all"             - Everything
    """
    return get_onelab_knowledge(topic)


@mcp.tool()
def motor_topology_optimization(topic: str = "wakao_ae_ls") -> str:
    """
    SynRM topology optimization (Wakao 2025 autoencoder + level-set).

    Args:
        topic: One of:
            "wakao_ae_ls"            - Autoencoder + LS hybrid method
            "liu_thesis_application" - Liu Xinyao end-to-end recipe
            "pareto_navigation"      - 2D latent-space Pareto walk
            "all"                    - Everything
    """
    return get_topology_opt_knowledge(topic)


@mcp.tool()
def motor_darwin_model(topic: str = "overview") -> str:
    """
    Darwin-model time-domain formulation (capacitive + inductive coupling).

    Args:
        topic: One of:
            "overview"     - When to use Darwin vs MQS / EQS / full Maxwell
            "formulation"  - Coulomb-gauge A-φ-χ symmetric system
            "vs_peec"      - Darwin (FE) vs PEEC (IE) practical comparison
            "all"          - Everything
    """
    return get_darwin_knowledge(topic)


@mcp.tool()
def motor_femm_transient(topic: str = "lab_recommendation") -> str:
    """
    FEMM newbuild transient solver — Lange-Henrotte-Hameyer 2009
    incremental-permeability linearization with sliding-band air-gap BC.

    Source: David Meeker, https://www.femm.info/doku/doku.php?id=newbuild
    + E. Lange, F. Henrotte, K. Hameyer, IEEE Trans. Mag. 45(3):1258-1261,
    2009 (DOI: 10.1109/TMAG.2009.2012585).

    Distinguishing feature: factor the nonlinear FE solve from time
    integration.  At a chosen rotor angle, 1 nonlinear FE solve +
    N_phase linear incremental-permeability solves (shared factorization)
    → extract `(L_inc(i, theta), e_bemf, T_em)` → circuit simulator
    integrates `L_inc di/dt + R i = v - e_bemf` at a much smaller PWM
    time step without re-running the FE.  Rotor angle advances
    analytically between FE refreshes.  Much faster than fully-coupled
    nonlinear FE+circuit (ONELAB) for controller-in-the-loop simulations.

    Args:
        topic: One of:
            "lab_recommendation"     - ★ Lab-canonical path + Hameyer/Henrotte refs
            "lange_2009"             - Underlying 2009 paper formulation
            "henrotte_no_rotation"   - Why the rotor never rotates in the FE
            "femm_newbuild"          - Meeker's FEMM implementation + antunes.fem
            "sliding_vs_moving_band" - Sliding-band BC vs moving-band remesh
            "ngsolve_recipe"         - Port to NGSolve with shared factorization
            "simulink_coupling"      - Simulink S-function + HIL/RT notes
            "all"                    - Everything
    """
    return get_femm_transient_knowledge(topic)


@mcp.tool()
def motor_henrotte_lineage(topic: str = "research_arc") -> str:
    """
    The Henrotte–Hameyer–RWTH research arc (energy-consistent E&M FE).

    Consolidates 7 key papers from the 1993-2018 ULg/UCL/RWTH school
    that culminated in the Lange-Henrotte-Hameyer 2009 field-circuit
    coupling.  All available in the Sugahara lab library.

    Args:
        topic: One of:
            "axisym_1993"       - Henrotte basis s=r²/2 (foundation of axifem)
            "source_field_1997" - Dular-Henrotte cut-aware source field
            "energy_hys_2006"   - Energy-based vector hysteresis (friction)
            "variational_2013"  - Variational FE-friendly extension
            "jacques_2018"      - 254p PhD monograph (canonical reference)
            "carstensen_2007"   - 187p SRM winding eddy currents PhD
            "research_arc"      - Synthesis: how the papers connect + roadmap
            "all"               - Everything
    """
    return get_henrotte_lineage_knowledge(topic)


@mcp.tool()
def motor_hollaus_eddy(topic: str = "intro") -> str:
    """
    Karl Hollaus / TU Wien MSFEM for laminated-iron eddy currents.

    THE answer to "Lange-Henrotte-Hameyer doesn't handle eddy currents".
    MSFEM homogenizes the through-the-thickness lamination structure
    so a 3D-laminated machine reduces to a 2D-1D or single-cell
    problem with ~100x fewer DOFs at <1% loss error.

    Args:
        topic: One of:
            "intro"             - Why standard FE fails; MSFEM idea
            "msfem_2d1d"        - 2D-1D MSFEM T-formulation (rotating machines)
            "t_phi_open"        - T-Phi-Phi MSFEM for OPEN circuits (motors)
            "effective_material"- EM homogenization (practical workhorse)
            "error_estimator"   - Equilibrated error estimator (Prager-Synge)
            "msfem_3d"          - Full 3D MSFEM for end-region asymmetric cases
            "nonlinear_esi"     - Nonlinear ESI in magnetic scalar potential
                                  (cross-references radia_mcp.ih.esim_cell_problem)
            "msfem_plus_mor"    - Cascade with POD/snapshots for transient PWM
            "nonasymptotic_homogenization" - Schobinger-Hollaus-Tsukerman 2020:
                                             two mu_eff definitions, resonance
            "msfem_mor_deim_detail" - Hollaus-Schobinger 2024: MSFEM+MOR+DEIM
                                       with structure-preserving DEIM
            "hanser_2025_circuit_coupling" - Hanser 2025 COMPEL: A-formulation
                                              + Schur for circuit coupling
                                              (300x-5000x speedup, IMPLEMENTATION
                                              ROADMAP for calc_motor_lamination v3)
            "hollaus_2014_nonlinear_two_scale" - Foundational two-scale paper
                                                   (Hollaus-Hannukainen-Schoberl 2014)
            "hierarchical_error_estimator" - Schobinger-Hollaus 2021 hierarchical
                                              error estimator for 3D MSFEM
            "frljic_2026_perpendicular_flux" - Frljic et al. 2026: anisotropic
                                                 mu_eff for open-type cores with
                                                 perpendicular flux
            "motor_strategy"    - SYNTHESIS: hybrid stack for eddy-current-capable
                                  motor analysis combining Lange-HH + MSFEM + ESIM
            "all"               - Everything
    """
    return get_hollaus_eddy_knowledge(topic)


@mcp.tool()
def motor_em_force_recipe(topic: str = "method_choice") -> str:
    """
    Practical NGSolve EM-force recipe for motor analysis.

    Forwards to `differential_forms_em_force_recipe` (radia_mcp.
    differential_forms.em_force_ngsolve_recipe_knowledge).  Tells you
    which method to use for motor torque / Maxwell stress on stator
    teeth / Lorentz on coil ends / etc.

    Achieved 0.5% Newton-3 violation on iron+coil 2D benchmark via
    Mertens-Hameyer 2007 v_order matching + far-field outer Dirichlet
    box (10x better than Pile-2018 "production" 5%).

    Args:
        topic: One of:
            "method_choice"         - Decision tree
            "high_order_recipe"     - Mertens-Hameyer 2007 prescription
            "common_pitfalls"       - 6 traps that destroy accuracy
            "validation_protocol"   - Newton-3 + path-independence
            "implementation_status" - calc_em_force.py status
            "full_example"          - Working code template
            "all"                   - Everything
    """
    # Lazy import to keep cross-package dependency soft
    from radia_mcp.differential_forms.em_force_ngsolve_recipe_knowledge \
        import get_em_force_ngsolve_recipe
    return get_em_force_ngsolve_recipe(topic)


@mcp.tool()
def motor_em_force_extras(topic: str = "all") -> str:
    """
    Forward to `differential_forms_em_force_extras` -- advanced EM force
    topics beyond the 7-method catalog.

    Args:
        topic: One of:
            "permanent_magnet_force"  - PM force formulations (bound
                                          current K = M x n, Kelvin,
                                          scalar charge)
            "energy_method_classical" - constant flux vs current,
                                          Legendre W <-> W'
            "sensitivity_vwp"         - shape derivative (Pironneau-
                                          Allaire), proper VWP fix
            "lorentz_canonical"       - Lorentz as first-class method
            "meissner_force"          - Type-I superconductor mu_r -> 0
            "all"                     - Everything
    """
    from radia_mcp.differential_forms.em_force_extras_knowledge \
        import get_em_force_extras
    return get_em_force_extras(topic)


@mcp.tool()
def motor_hollaus_genealogy(view: str = "by_topic") -> str:
    """
    Visualize the Karl Hollaus / TU Wien MSFEM research genealogy
    (1992-2026 paper lineage with citation DAG).

    The MSFEM stack now spans 13+ papers (see motor_hollaus_eddy for the
    physics content). This tool maps the dependency graph: which paper
    builds on which, grouped by research thread.

    Args:
        view: One of:
            "by_topic"  - Grouped by research thread (default; best for
                          orienting a new reader)
            "chrono"    - Chronological flat list with DOIs and parents
            "ascii"     - ASCII art tree (roots -> descendants)
            "dot"       - Graphviz DOT format (paste into
                          https://dreampuf.github.io/GraphvizOnline/)
            "all"       - All four views
    """
    return get_hollaus_genealogy(view)


@mcp.tool()
def motor_bibliography(query: str = "") -> str:
    """Search the motor analysis bibliography catalog."""
    return get_bibliography_index(query)


@mcp.tool()
def motor_tritool_cross_reference(topic: str = "overview") -> str:
    """
    Tri-tool cross-reference: FEMM / JMAG / radia-ngsolve (相互学習).

    Cross-learning that ties the lab's three motor-FEA tools together so each
    strengthens the others. Mirrors the shared cross_ref.json carried by the
    two dedicated lab-internal knowledge servers (mcp-server-femm and
    mcp-server-jmag), whose commercial-tool specifics stay lab-private.

    Args:
        topic: One of:
            "overview"          - The ecosystem + the two new MCP servers
            "capability_matrix" - Per-capability FEMM/JMAG/radia table + strongest
            "radia_can_exceed"  - Where radia-ngsolve matches or beats the others
            "jmag_only"         - Genuine JMAG-only capabilities (gaps)
            "femm_role"         - FEMM as the open-source 2D yardstick + .fem strategy
            "roadmap"           - radia-motor strengthening roadmap (ranked)
            "all"               - Everything
    """
    return get_tritool_cross_reference(topic)


# ============================================================
# MCP Prompts
# ============================================================

@mcp.prompt()
def new_motor_simulation(motor_type: str = "pmsm") -> str:
    """Set up a new electric-machine simulation."""
    t = motor_type.strip().lower()
    base = (
        f"Set up a {motor_type} simulation using the radia-mcp motor toolchain.\n\n"
        "Recommended path:\n"
        "1. motor_femm_transient('lab_recommendation') for the **lab-canonical**\n"
        "   transient approach (Lange-Henrotte-Hameyer 2009).\n"
        "2. motor_onelab('overview') for ONELAB ElectricMachines layout\n"
        "   (reference geometry templates).\n"
        "3. motor_onelab('motor_types') to confirm Region group names.\n"
        "4. motor_onelab('formulation') for the 2D weak form.\n"
        "5. motor_onelab('ngsolve_xlate') for the NGSolve port.\n"
        "6. motor_onelab('analysis_modes') to pick Static / TD / FD.\n\n"
    )
    if t == "pmsm":
        base += (
            "PMSM-specific:\n"
            "- Rotor: surface or interior permanent magnets (Br vector / region).\n"
            "- Anti-periodicity factor: NbrPolePairs (1 sector = 1 pole-pair).\n"
            "- Cogging torque map: motor_onelab('analysis_modes') Static, sweep θ.\n"
            "- d/q inductance: Static sweep I_d, I_q.\n"
        )
    elif t == "synrm":
        base += (
            "SynRM-specific:\n"
            "- No PM, no rotor coil — pure reluctance machine.\n"
            "- Dominant design objectives: T_ave, T_rip.\n"
            "- Use motor_topology_optimization('wakao_ae_ls') for shape design.\n"
            "- motor_topology_optimization('liu_thesis_application') for end-to-end recipe.\n"
        )
    elif t in ("im", "induction"):
        base += (
            "Induction-motor-specific:\n"
            "- Use Frequency-Domain analysis with slip-frequency for steady state.\n"
            "- Use Time-Domain for transient run-up.\n"
            "- Cage bars are RotorC region; end-rings via external circuit.\n"
            "- For broken-bar fault, see ONELAB `im` bundle.\n"
        )
    elif t == "srm":
        base += (
            "SRM-specific:\n"
            "- Doubly salient, no PM, no rotor coil, no cage.\n"
            "- Phase-pulsed excitation (square-wave or current-controlled).\n"
            "- Strongly non-linear ν(B) — use Newton with damping.\n"
        )
    elif t == "wfsm":
        base += (
            "WFSM-specific:\n"
            "- Rotor field winding (DC) + damper cage.\n"
            "- For sub-transient analysis, treat damper bars as RotorC.\n"
        )
    else:
        base += (
            f"Custom motor type '{motor_type}':\n"
            "- Map your rotor source to one of: Rotor_Magnets, RotorC, Rotor_Inds.\n"
            "- Choose anti-periodicity factor from pole count.\n"
        )
    return base


# ============================================================
# Entry point
# ============================================================



register_status_tool(
    mcp,
    server_name='mcp-server-motor',
    description='Motor analysis: ONELAB transient, Hollaus effective material (lamination), Wakao autoencoder topology, Kaimori-Mifune Darwin TD',
    subpackage='radia_mcp.motor',
    related_servers=["electromagnet", "topology-optimization", "magnetic-materials"],
    optional_deps=["radia", "ngsolve"],
)


def main():
    if "--selftest" in sys.argv:
        print("Motor MCP server self-test:")
        from .onelab_knowledge import SECTIONS as O_SEC
        from .topology_opt_knowledge import SECTIONS as T_SEC
        from .darwin_model_knowledge import SECTIONS as D_SEC
        from .femm_transient_knowledge import SECTIONS as F_SEC
        from .henrotte_lineage_knowledge import SECTIONS as H_SEC
        from .hollaus_eddy_knowledge import SECTIONS as E_SEC
        from .tritool_cross_reference_knowledge import SECTIONS as X_SEC
        for k in O_SEC:
            r = motor_onelab(k)
            print(f"  motor_onelab({k!r}): {len(r)} chars")
            assert len(r) > 100, f"ONELAB topic {k} too short"
        for k in T_SEC:
            r = motor_topology_optimization(k)
            print(f"  motor_topology_optimization({k!r}): {len(r)} chars")
            assert len(r) > 100, f"Topology-opt topic {k} too short"
        for k in D_SEC:
            r = motor_darwin_model(k)
            print(f"  motor_darwin_model({k!r}): {len(r)} chars")
            assert len(r) > 100, f"Darwin topic {k} too short"
        for k in F_SEC:
            r = motor_femm_transient(k)
            print(f"  motor_femm_transient({k!r}): {len(r)} chars")
            assert len(r) > 100, f"FEMM transient topic {k} too short"
        for k in H_SEC:
            r = motor_henrotte_lineage(k)
            print(f"  motor_henrotte_lineage({k!r}): {len(r)} chars")
            assert len(r) > 100, f"Henrotte lineage topic {k} too short"
        for k in E_SEC:
            r = motor_hollaus_eddy(k)
            print(f"  motor_hollaus_eddy({k!r}): {len(r)} chars")
            assert len(r) > 100, f"Hollaus eddy topic {k} too short"
        for k in X_SEC:
            r = motor_tritool_cross_reference(k)
            print(f"  motor_tritool_cross_reference({k!r}): {len(r)} chars")
            assert len(r) > 100, f"Tri-tool topic {k} too short"
        p = new_motor_simulation("synrm")
        print(f"  new_motor_simulation('synrm'): {len(p)} chars")
        assert "wakao_ae_ls" in p
        print("  PASSED")
        return

    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
