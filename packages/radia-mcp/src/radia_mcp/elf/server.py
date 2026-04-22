"""
ELF MCP Server

Provides documentation and knowledge for the ELF600 electromagnetic field
analysis suite (MAGIC / ELFIN / BEAM solvers) — file formats, commands,
element types, solver options, and mesh scripting.

Topics:
- .mai file format (analysis/solver control input)
- .mei file format (mesh generation script)
- .meg file format (compiled mesh data)
- MAGIC solver (magnetostatic BEM)
- ELFIN solver (eddy current BEM)
- BEAM solver (charged particle tracking)
- Element type naming convention
- B-H curve specification
- IPM motor analysis
- Inductance computation (LscLl)
- Magnetization / demagnetization

Usage:
    mcp-server-elf              # Start MCP server (stdio transport)
    mcp-server-elf --selftest   # Run self-test
"""

import sys

from mcp.server.fastmcp import FastMCP

from .elf_knowledge import get_elf_documentation

mcp = FastMCP("mcp-server-elf")


# ============================================================
# MCP Tools
# ============================================================

@mcp.tool()
def elf_usage(topic: str = "all") -> str:
    """
    Get ELF600 electromagnetic field analysis documentation.

    ELF600 is a BEM-based electromagnetic analysis suite with three solvers:
    MAGIC (magnetostatic), ELFIN (electrostatic), BEAM (particle tracking).
    Input consists of .mai (analysis control), .mei (mesh script), and
    .meg (compiled mesh) files.

    Args:
        topic: Documentation topic. Options:
            "all"              - Complete documentation (~60k chars)
            "overview"         - ELF600 suite overview, solver modules, tools
            "mai_format"       - .mai file format (PRE keywords, all solvers)
            "mei_format"       - .mei file format (mesh scripting structure)
            "meg_format"       - .meg file format (symmetry, nodes, elements)
            "magic"            - MAGIC magnetostatic solver details
            "elfin"            - ELFIN electrostatic solver (D-E curves)
            "beam"             - BEAM particle tracking (RELA, BINT, GRAV)
            "element_types"    - All element types (MAGIC/ELFIN/BEAM, DOF)
            "bh_curves"        - B-H curves, recoil, extrapolation
            "sol_commands"     - SOL blocks, NONL strategy, PASS optimization
            "mei_commands"     - IEmesh commands (AA/AN/AE/ME/BE/TE/NB/...)
            "ipm_motor"        - IPM motor Ld/Lq calculation workflow
            "inductance"       - Lsc (JIS) and Ll (IEEJ) with 6 samples
            "magnetization"    - Magnetization (MAGNE2) and demagnetization
            "examples"         - Example file catalog
            "meg_export"       - MEG export from Cubit (labels, DIM)
            "treasure_box"     - Quick-reference tables (element-PRE map, etc.)
            "sinusoidal"       - SOL MOMC, complex permeability, AC force
            "anisotropy"       - HBA1/HBA2, VEC1/VEC3, lamination
            "sted"             - Steady-state eddy current motion
            "meshing"          - Element quality, gaps, connectivity rules
            "convergence"      - Nonlinear convergence troubleshooting
            "force_methods"    - FORC vs FORT vs FIXB comparison
            "errors"           - Error messages (ELF-Q/E/W, 160+ codes)
            "iemesh"           - IEmesh tool overview and command families
            "tools"            - Wmap3, MagFilter2, MaiEditor3, ELF/Bench
    """
    return get_elf_documentation(topic)


# ============================================================
# MCP Prompts
# ============================================================

@mcp.prompt()
def new_elf_analysis(geometry: str, solver: str = "MAGIC") -> str:
    """Set up a new ELF analysis."""
    solver = solver.upper()
    if solver == "MAGIC":
        return (
            f"Set up a MAGIC magnetostatic analysis for: {geometry}\n\n"
            "Workflow:\n"
            "1. Create .mei mesh script (MODEL MAGIC T/K/R)\n"
            "2. Define coils (MCL), iron (MWL/MMB), magnets (HBRM)\n"
            "3. Add evaluation contour (MCO) after SPACE\n"
            "4. Run IEmesh to generate .meg\n"
            "5. Create .mai: USE MAGIC 3.50, PRE block (B-H, coils), SOL MOME\n"
            "6. Add SOL FIEL, SOL FORC/FORT as needed\n"
        )
    elif solver == "ELFIN":
        return (
            f"Set up an ELFIN eddy current analysis for: {geometry}\n\n"
            "Workflow:\n"
            "1. Create .mei mesh script (MODEL ELFIN T/K/R)\n"
            "2. Define conductors (ESC), magnetic bodies (EMB)\n"
            "3. Add evaluation contour (ECO) after SPACE\n"
            "4. Run IEmesh to generate .meg\n"
            "5. Create .mai: USE ELFIN 3.50, PRE block (VOL1, EDSC), SOL MOME\n"
            "6. Add SOL FIEL for field computation\n"
        )
    else:
        return (
            f"Set up a BEAM particle tracking analysis for: {geometry}\n\n"
            "Workflow:\n"
            "1. Compute EM fields with MAGIC or ELFIN first\n"
            "2. Create .mei mesh script (MODEL BEAM T/K/R)\n"
            "3. Define beam source nodes\n"
            "4. Create .mai: USE BEAM 3.50, PRE (CHAR, MASS, VOLB, FILE)\n"
            "5. SOL BEAM for tracking\n"
        )


# ============================================================
# Entry point
# ============================================================

def main():
    if "--selftest" in sys.argv:
        from radia_mcp.common.utf8_stdout import use_utf8_stdout
        use_utf8_stdout()
        print("ELF MCP server self-test:")
        topics = [
            "overview", "mai_format", "mei_format", "meg_format",
            "magic", "elfin", "beam", "element_types", "bh_curves",
            "sol_commands", "mei_commands", "ipm_motor", "inductance",
            "magnetization", "examples", "meg_export",
            "treasure_box", "sinusoidal", "anisotropy", "sted",
            "meshing", "convergence", "force_methods", "errors",
            "iemesh", "tools",
        ]
        for t in topics:
            result = elf_usage(t)
            print(f"  elf_usage('{t}'): {len(result)} chars")
            assert len(result) > 50, f"Topic '{t}' too short"
        print("  PASSED")
        return

    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
