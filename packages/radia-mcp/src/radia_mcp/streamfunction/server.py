"""
Stream-Function (SF) Coil Design MCP Server (radia_mcp.streamfunction)

SF-focused front door for the Radia stream-function coil-design framework:

- Kernel-agnostic (ACA+)+TSVD least-norm solver (coils Biot-Savart OR
  magnets via the Radia MMM material kernel -- the matrix entry is a callback)
- FE-direct psi (continuous H1 GridFunction) on ANY surface: plane / cylinder
  / sphere / conformal / 3D-printed former
- Regularisation menu (L2 / H1 / sigma-weighted / inductance / L-inf) folded
  onto ONE ACA factorisation via RegularizedTSVD; folded TIKHONOV ("+ alpha I"
  core) -> the (field homogeneity, PEAK current density) PARETO FRONT
- Four stackable levers to push the front: Tikhonov alpha (L-curve), L-inf
  minimax, geometry (former size / cylinder length), SHEET-METAL (板金) forming
- Single-stroke ("一筆書き") chain + sheet-metal wire distortion (one feed)

The detailed knowledge lives in this server's own
``radia_mcp.streamfunction.knowledge.aca_tsvd`` (moved from radia-ngsolve in
2026-06: SF coil design is not a general NGSolve usage).  The ``streamfunction``
tool adds a dedicated SF overview + topic map over that knowledge.

Usage:
    mcp-server-radia-streamfunction              # Start MCP server (stdio)
    mcp-server-radia-streamfunction --selftest   # Run self-test
"""

import sys

from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool, register_topics_tool

from .streamfunction_knowledge import get_streamfunction_documentation, TOPICS

mcp = FastMCP("mcp-server-radia-streamfunction")


# ============================================================
# MCP Tools
# ============================================================

@mcp.tool()
def streamfunction(topic: str = "overview") -> str:
    """
    Get Stream-Function (SF) coil-design documentation.

    Design a surface coil that produces a target field over a design volume:
    a surface current K = n_hat x grad(psi) is written from a scalar stream
    function psi; its equal-increment iso-contours ARE the wire paths.  The
    least-norm system A psi = B_target is solved with the kernel-agnostic
    (ACA+)+TSVD solver and a chosen regularisation, then connected into one
    continuous wire (single-stroke).

    Args:
        topic: Documentation topic. Options:
            "overview"        - framework, pipeline, topic map (front door)
            "theory"/"method" - SFM + (ACA+)+TSVD least-norm math
            "api"             - radia.stream_function API
            "kernel_agnostic" - callback matrix entry (coils OR magnets)
            "regularized"     - regularisation menu + folded Tikhonov closed form
            "pareto"          - (homogeneity, peak current density) Pareto
                                front + 4 levers + sheet-metal (板金) forming
            "material_aware"  - SF design WITH IRON (--iron-vol): Kelvin-FEM
                                reaction folded into the whole pipeline,
                                obs-adjoint scalability, exact-quad source
            "low_turn"        - low/many-turn discrete refinement: greedy
                                constructive (--greedy-turns, monotone) +
                                connector-aware, --pin-tiling driven array,
                                --optimize-levels, --distort, single-pass
            "single_stroke"   - one-wire chain, FE-direct on arbitrary formers,
                                sheet-metal wire distortion
            "fe_direct"       - FE-direct H1 psi on plane/cylinder/sphere
            "deformation"     - surface-deformation outer loop (accuracy)
            "cmaes"           - SA-25-020 CMA-ES outer loop
            "performance"     - (ACA+)+TSVD amortisation numbers
            "validation"      - analytic-benchmark checks
            "literature"      - SFM lineage (Turner / Peeren / current potential)
            "workflow"        - end-to-end demo recipes
            "fusion"          - stellarator Stage-2 (REGCOIL/NESCOIL/FOCUS)
            "clebsch_3d"      - [IN PROGRESS] 3D stream function + cohomology:
                                de Rham unification, ACA+TSVD on current-Clebsch,
                                Tampere bidirectional map, honest 3D frontier

    Returns:
        Markdown/plain-text documentation for the requested topic.
    """
    return get_streamfunction_documentation(topic)


# ============================================================
# MCP Prompts
# ============================================================

@mcp.prompt()
def new_sf_coil_design(target: str = "uniform Bz") -> str:
    """Set up a new stream-function coil design for a target field."""
    return (
        f"Design a stream-function surface coil for a '{target}' target.\n\n"
        "Use the streamfunction(topic) tool for documentation.\n\n"
        "Workflow:\n"
        "1. Choose the former surface (plane / cylinder / sphere) and the DSV\n"
        "   target points.  streamfunction('kernel_agnostic') for the A(i,j)\n"
        "   matrix-entry callback (coils Biot-Savart OR magnets MMM).\n"
        "2. Solve A psi = B with (ACA+)+TSVD + a regularisation.\n"
        "   streamfunction('regularized') -- pick H1 (min surface-current\n"
        "   energy) as the physical default; sigma-weighted for ohmic;\n"
        "   L-inf for a peak cap.\n"
        "3. If you care about (field homogeneity vs PEAK current density),\n"
        "   trace and PUSH the Pareto front: streamfunction('pareto')\n"
        "   -- Tikhonov alpha (L-curve) + L-inf + geometry + sheet-metal (板金)\n"
        "   forming.  Folded Tikhonov sweeps the front at ~50 us/point.\n"
        "4. Extract equal-increment iso-contours and connect them into ONE\n"
        "   wire: streamfunction('single_stroke') (field_aware / Kuijpers).\n"
        "5. Optional: single-current sheet-metal WIRE distortion to cancel the\n"
        "   single-stroke residual without extra feeds.\n\n"
        "FE-direct psi (streamfunction('fe_direct')) is required for non-grid\n"
        "formers (sphere / conformal / 3D-printed).\n"
    )


# ============================================================
# Status + topics + entry point
# ============================================================

register_status_tool(
    mcp,
    server_name='mcp-server-radia-streamfunction',
    description='Stream-function coil design: (ACA+)+TSVD least-norm, FE-direct '
                'psi, regularisation / folded-Tikhonov Pareto front, '
                'single-stroke chain, sheet-metal (板金) levers',
    subpackage='radia_mcp.streamfunction',
    related_servers=["radia-ngsolve", "bayesian-opt", "nmr-mri",
                     "electromagnet", "peec"],
    optional_deps=["radia", "ngsolve"],
)

register_topics_tool(
    mcp,
    server_name='mcp-server-radia-streamfunction',
    topics=TOPICS,
)


def main():
    if "--selftest" in sys.argv:
        print("Stream-function MCP server self-test:")
        for t in ["overview", "theory", "api", "kernel_agnostic",
                  "regularized", "pareto", "material_aware", "low_turn",
                  "single_stroke", "fe_direct",
                  "deformation", "cmaes", "performance", "validation",
                  "literature", "workflow", "fusion", "clebsch_3d"]:
            result = streamfunction(t)
            print(f"  streamfunction('{t}'): {len(result)} chars")
            assert len(result) > 100, f"Topic '{t}' too short"
        prompt = new_sf_coil_design("uniform Bz")
        print(f"  new_sf_coil_design('uniform Bz'): {len(prompt)} chars")
        assert "single_stroke" in prompt
        # the Pareto / sheet-metal section must be reachable
        assert "板金" in streamfunction("pareto") or \
               "sheet-metal" in streamfunction("pareto")
        print("  PASSED")
        return

    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
