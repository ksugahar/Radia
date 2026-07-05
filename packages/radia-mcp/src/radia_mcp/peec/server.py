"""
PEEC (Partial Element Equivalent Circuit) MCP Server (radia_mcp.peec)

Provides tools for:
- PEEC architecture: Loop-Star, filament-panel, node-segment topology
- PyPEECBuilder and PEECCircuitSolver API guidance
- Multi-filament subdivision (nwinc/nhinc) for skin/proximity effect
- Surface impedance methods (Bessel, Dowell, ESIM)
- Magnetic-material coupling is handled by HDiv-VIM / reduced FEM, not PEEC
- FastHenry .inp file parsing
- PRIMA model order reduction and SPICE extraction
- ngsolve.bem integration notes

Usage:
    mcp-server-peec              # Start MCP server (stdio transport)
    mcp-server-peec --selftest   # Run self-test

Promoted from legacy private source tree to public radia-mcp on
2026-04-24 (single Radia monorepo).  General PEEC inductance formulae
live in radia_mcp.radia_ngsolve.peec_inductance_knowledge; this
subpackage holds the LAB PEEC workflow (Loop-Star, FastHenry parsing,
PyPEECBuilder/PEECCircuitSolver API, Bessel/Dowell/ESIM surface
impedance, PRIMA MOR, SPICE extraction).
"""

import sys

from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool

from .peec_knowledge import get_peec_documentation
from .carstensen_ac_copper_loss_knowledge import get_carstensen_knowledge
from .hoibc_knowledge import get_hoibc_documentation

mcp = FastMCP("mcp-server-peec")


# ============================================================
# MCP Tools
# ============================================================

@mcp.tool()
def peec_carstensen_ac_loss(topic: str = "motivation") -> str:
    """
    Carstensen-Dowell analytical AC copper-loss formulas for stranded
    slot windings.  Avoids meshing individual strands by using
    Kelvin-function-based skin + proximity factors per layer.

    Reference: C. Carstensen, "Eddy Currents in Windings of Switched
    Reluctance Machines", PhD diss., RWTH Aachen ISEA, 2007 (187p).
    Lab library: public-safe curated corpus

    Args:
        topic: One of:
            "motivation"        - Why analytic AC loss matters for slot windings
            "skin_round"        - Bessel-function skin effect in round wire
            "proximity_round"   - Proximity effect from external H
            "dowell_layered"    - Classical Dowell formula for layered slot
            "combined_skin_prox"- Carstensen per-layer breakdown
            "field_mapping"     - FE -> per-layer H_slot extraction
            "python_recipe"     - Drop-in Python post-processor
            "validation_cases"  - Reference cases for verification
            "motor_application" - End-to-end recipe with calc_motor_transient.py
            "all"               - Everything
    """
    return get_carstensen_knowledge(topic)


@mcp.tool()
def peec_usage(topic: str = "all") -> str:
    """
    Get PEEC (Partial Element Equivalent Circuit) documentation.

    Complete workflow for extracting circuit parameters (L, R, C, M)
    from conductor geometry using integral equations. Covers node-segment
    topology, MNA port impedance, multi-filament subdivision, surface
    impedance, PRIMA model order reduction, and
    SPICE netlist extraction.

    Args:
        topic: Documentation topic. Options:
            "all"            - Complete documentation
            "overview"       - PEEC architecture, system equation, Loop-Star
            "builder"        - PyPEECBuilder API (add_node_at, add_segment, add_port)
            "solver"         - PEECCircuitSolver (MNA, port impedance, sweep)
            "multi_filament" - nwinc/nhinc subdivision for skin/proximity
            "fasthenry"      - FastHenry .inp parser usage
            "magnetic"       - Magnetic-material coupling policy
            "sibc"           - Surface impedance (Bessel, Dowell, ESIM)
            "prima"          - PRIMA model order reduction, SPICE export
            "ngsolve_bem"    - ngsolve.bem integration, frequency range guide
            "pitfalls"       - Common mistakes and gotchas
    """
    return get_peec_documentation(topic)


@mcp.tool()
def peec_hoibc(topic: str = "vs_sibc_decision") -> str:
    """
    HOIBC (Higher Order Impedance Boundary Conditions) — extension of
    the 1st-order Leontovich SIBC to capture surface curvature and
    remain accurate when skin depth becomes comparable to the radius
    of curvature.

    The standard SIBC in `peec_usage(topic='sibc')` is the 1st-order
    Leontovich approximation.  For:
      - Thin wires at low freq (delta ~ wire radius)
      - High-mu_r magnetic conductors (effective delta = mu_r * delta)
      - Coil end-regions with sharp curvature
    you need HIGHER ORDER (Mitzner correction or beyond).

    Reference papers (public-safe curated corpus):
      - Course G2ELab 2018 lesson 2 (derivation)
      - Dong-Di Rienzo 2020 IEEE Access 8:186496 (3D FEM/BEM)
      - Bilicz-Badics-Pavo 2023 ISEM (wide-band nonlocal)

    Args:
        topic: One of:
            "derivation"            - Perturbation theory (Rytov-Senior-
                                       Mitzner, small parameter p̃ =
                                       mu_r * delta / D)
            "dong_di_rienzo_2020"   - 3-Laplace-BVP algorithm (PEC →
                                       Leontovich → Mitzner); multi-freq
                                       advantage; numerical validation
                                       table
            "bilicz_2023_nonlocal"  - Wide-band alternative: 2D cross-
                                       section eddy-current BVP per
                                       frequency (works at any delta)
            "vs_sibc_decision"      - Decision tree: when SIBC suffices,
                                       when HOIBC needed, when nonlocal
                                       IBC needed
            "radia_application"     - How to extend calc_fem_kelvin.py
                                       with --ibc-order N flag
            "sibc_family"           - Variant landscape: XFEM SIBC,
                                       FDTD SIBC, Nonlinear SIBC,
                                       Multi-layer / coated conductor,
                                       WPT wide-band (lab public-safe curated corpus)
            "mathematica_derivation"- Wolfram Language recipes: small-
                                       parameter p_tilde, BVP_0/1/2
                                       Neumann RHS symbolic, Mitzner
                                       sphere closed form, prolate-
                                       ellipsoid principal radii, freq-
                                       domain switch (s -> j*omega).
            "ngsolve_recipes"       - Production NGSolve code patterns:
                                       BVP_0 (PEC), BVP_1 (Leontovich),
                                       BVP_2 (Mitzner) Laplace cascade;
                                       multi-frequency zero-cost sweep;
                                       per-panel curvature via _compute_
                                       panel_local_radii; prolate-
                                       ellipsoid golden test.
            "all"                   - Everything
    """
    return get_hoibc_documentation(topic)


# ============================================================
# MCP Prompts
# ============================================================

@mcp.prompt()
def new_peec_simulation(conductor_type: str) -> str:
    """Set up a new PEEC simulation for circuit parameter extraction."""
    guidance = {
        "wire": (
            "Round wire conductor:\n"
            "- Use circular_sibc (Bessel I0/I1) for skin effect\n"
            "- nwinc/nhinc not needed for round cross-section\n"
        ),
        "pcb": (
            "PCB trace / rectangular conductor:\n"
            "- Use Dowell SIBC or nwinc/nhinc multi-filament\n"
            "- Typical: nwinc=3, nhinc=1 for thin traces\n"
        ),
        "busbar": (
            "Busbar / power conductor:\n"
            "- Use nwinc=5, nhinc=5 for thick conductors\n"
            "- Consider Dowell SIBC for d << w geometry\n"
        ),
        "coil_on_core": (
            "Coil on magnetic core (transformer/inductor):\n"
            "- Use PEEC for conductors only\n"
            "- Use HDiv-VIM / reduced FEM for the magnetic core\n"
            "- Exchange fields through the application workflow\n"
        ),
        "litz": (
            "Litz wire (parallel strands):\n"
            "- Model each strand as separate segment (same nodes)\n"
            "- Proximity effect captured by mutual inductance\n"
        ),
    }

    conductor_lower = conductor_type.lower().strip()
    specific = guidance.get(conductor_lower, "")
    if not specific:
        types_list = ", ".join(guidance.keys())
        specific = (
            f"Unknown conductor type: '{conductor_type}'.\n"
            f"Known types: {types_list}\n"
            "Proceeding with general PEEC setup guidance.\n"
        )

    return (
        f"Set up a PEEC simulation for: {conductor_type}\n\n"
        f"{specific}\n"
        "General workflow:\n"
        "1. Define geometry: PyPEECBuilder (add_node_at, add_connected_segment)\n"
        "2. Define ports: builder.add_port(n_pos, n_neg)\n"
        "3. Build topology: topo = builder.build_topology()\n"
        "4. Create solver: PEECCircuitSolver(topo)\n"
        "5. Compute: Z = solver.compute_port_impedance(freq)\n"
        "6. Sweep: Z_sweep = solver.frequency_sweep(freqs, Zs_func)\n"
        "7. Optional: PRIMA reduction for SPICE netlist\n\n"
        "Use the peec_usage tool for detailed API documentation.\n"
        "Key topics: 'builder', 'solver', 'sibc', 'magnetic', 'prima'\n"
    )


# ============================================================
# Entry point
# ============================================================



register_status_tool(
    mcp,
    server_name='mcp-server-peec',
    description='PEEC filament/panel, FastHenry parser, HOIBC, Carstensen AC copper loss',
    subpackage='radia_mcp.peec',
    related_servers=["radia-ngsolve", "ih", "litz-transmission"],
    optional_deps=["radia", "scipy"],
)


def main():
    if "--selftest" in sys.argv:
        print("PEEC MCP server self-test:")
        for topic in ["overview", "builder", "solver", "multi_filament",
                       "fasthenry", "magnetic", "sibc", "prima",
                       "ngsolve_bem", "pitfalls"]:
            doc = peec_usage(topic)
            print(f"  peec_usage('{topic}'): {len(doc)} chars")
        all_doc = peec_usage("all")
        print(f"  peec_usage('all'): {len(all_doc)} chars")
        print("  PASSED")
        return

    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
