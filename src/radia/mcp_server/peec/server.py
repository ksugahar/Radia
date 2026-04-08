"""
PEEC (Partial Element Equivalent Circuit) MCP Server

Provides tools for:
- PEEC architecture: Loop-Star, filament-panel, node-segment topology
- PyPEECBuilder and PEECCircuitSolver API guidance
- Multi-filament subdivision (nwinc/nhinc) for skin/proximity effect
- Surface impedance methods (Bessel, Dowell, ESIM)
- Coupled PEEC + MMM for conductors near magnetic materials
- FastHenry .inp file parsing
- PRIMA model order reduction and SPICE extraction
- ngsolve.bem integration notes

Usage:
    mcp-server-peec              # Start MCP server (stdio transport)
    mcp-server-peec --selftest   # Run self-test
"""

import sys

from mcp.server.fastmcp import FastMCP

from .peec_knowledge import get_peec_documentation

mcp = FastMCP("mcp-server-peec")


# ============================================================
# MCP Tools
# ============================================================

@mcp.tool()
def peec_usage(topic: str = "all") -> str:
    """
    Get PEEC (Partial Element Equivalent Circuit) documentation.

    Complete workflow for extracting circuit parameters (L, R, C, M)
    from conductor geometry using integral equations. Covers node-segment
    topology, MNA port impedance, multi-filament subdivision, surface
    impedance, coupled PEEC+MMM, PRIMA model order reduction, and
    SPICE netlist extraction.

    Args:
        topic: Documentation topic. Options:
            "all"            - Complete documentation
            "overview"       - PEEC architecture, system equation, Loop-Star
            "builder"        - PyPEECBuilder API (add_node_at, add_segment, add_port)
            "solver"         - PEECCircuitSolver (MNA, port impedance, sweep)
            "multi_filament" - nwinc/nhinc subdivision for skin/proximity
            "fasthenry"      - FastHenry .inp parser usage
            "coupled"        - CoupledPEECSolver (PEEC + MMM, Delta_L)
            "sibc"           - Surface impedance (Bessel, Dowell, ESIM)
            "prima"          - PRIMA model order reduction, SPICE export
            "ngsolve_bem"    - ngsolve.bem integration, frequency range guide
            "pitfalls"       - Common mistakes and gotchas
    """
    return get_peec_documentation(topic)


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
            "- Use CoupledPEECSolver with Radia magnetic objects\n"
            "- compute_coupling_matrix() for Delta_L\n"
            "- mu_r_imag for core loss\n"
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
        "4. Create solver: PEECCircuitSolver(topo) or CoupledPEECSolver(topo, [core])\n"
        "5. Compute: Z = solver.compute_port_impedance(freq)\n"
        "6. Sweep: Z_sweep = solver.frequency_sweep(freqs, Zs_func)\n"
        "7. Optional: PRIMA reduction for SPICE netlist\n\n"
        "Use the peec_usage tool for detailed API documentation.\n"
        "Key topics: 'builder', 'solver', 'sibc', 'coupled', 'prima'\n"
    )


# ============================================================
# Entry point
# ============================================================

def main():
    if "--selftest" in sys.argv:
        print("PEEC MCP server self-test:")
        for topic in ["overview", "builder", "solver", "multi_filament",
                       "fasthenry", "coupled", "sibc", "prima",
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
