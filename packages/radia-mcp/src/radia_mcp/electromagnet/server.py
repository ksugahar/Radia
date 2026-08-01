"""
Accelerator Electromagnet MCP Server (radia_mcp.electromagnet)

Provides tools for:
- CoilBuilder workflow (racetrack, saddle, racetrack coils)
- Cubit hex mesh -> NGSolve curved -> Kelvin transform pipeline
- Omega-reduced scalar potential formulation
- Hantila polarization method (LU once, back-substitution iteration)
- B-input Play/Energy hysteresis models
- IMA (Image Method of Analysis) sign selection
- Field harmonics / multipole analysis
- radia-em Clebsch hodograph panel mode and accelerator design references
- Accelerator fundamentals, beam-optics handoff, ramped/superconducting
  magnet engineering, measurement, and a curated textbook source guide

Usage:
    mcp-server-electromagnet              # Start MCP server (stdio transport)
    mcp-server-electromagnet --selftest   # Run self-test

Promoted from legacy private source tree to public
radia-mcp on 2026-04-24 (single Radia monorepo).  General Kelvin
and NGSolve/BEM knowledge live in radia_mcp.radia_ngsolve.*; this
subpackage holds electromagnet-specific topics (CoilBuilder,
Hantila polarization, B-input hysteresis, IMA sign selection,
field-harmonics / multipole analysis).
"""

import importlib.util
import sys
from pathlib import Path

from mcp.server.fastmcp import FastMCP

from ..common import register_status_tool, register_topics_tool
from .accelerator_fundamentals_knowledge import get_accelerator_source_guide
from .em_knowledge import TOPICS, get_electromagnet_documentation

mcp = FastMCP("mcp-server-electromagnet")


# ============================================================
# MCP Tools
# ============================================================

@mcp.tool()
def electromagnet_usage(topic: str = "overview") -> str:
    """
    Get accelerator electromagnet analysis documentation.

    Complete pipeline for accelerator magnet design and analysis:
    Radia CoilBuilder (Biot-Savart source) -> Cubit hex mesh ->
    NGSolve FEM (Omega-reduced + Kelvin) -> GMSH visualization.

    Covers dipole, quadrupole, and sextupole magnets with nonlinear
    iron yoke, hysteresis, and open boundary (Kelvin transformation).

    Args:
        topic: Documentation topic. Options:
            "overview"         - Pipeline architecture, unique capabilities
            "coilbuilder"      - CoilBuilder fluent API, examples
            "kelvin_workflow"  - Cubit -> NGSolve -> Kelvin pipeline
            "hantila"          - Hantila polarization (LU once)
            "hysteresis"       - B-input Play/Energy models
            "ima"              - Image Method sign selection
            "harmonics"        - Multipole analysis, FFT extraction
            "clebsch_hodograph" - radia-em Clebsch hodograph mode and links
                                  to accelerator pole-face design topics
            "accelerator_fundamentals" - Magnetic rigidity and requirements
            "beam_optics_contract" - Twiss/dispersion/tune field handoff
            "accelerator_magnet_types" - Magnet roles and technology choice
            "accelerator_magnet_design" - End-to-end engineering workflow
            "rapid_cycling_magnets" - Eddy current and power-supply design
            "superconducting_accelerator_magnets" - Conductor/quench design
            "accelerator_magnet_measurement" - Field QA and commissioning
            "accelerator_model_boundaries" - Coupled accelerator analyses
            "accelerator_sources" - Curated 12-textbook source guide
            "all"              - Complete documentation
    """
    return get_electromagnet_documentation(topic)


@mcp.tool()
def electromagnet_accelerator_sources(query: str = "") -> str:
    """Search the curated accelerator textbook source guide.

    The guide covers 12 textbooks and lecture notes used to build the
    accelerator fundamentals in this server. It returns bibliographic and
    topical locators only; source PDF text is not distributed.

    Args:
        query: Optional keyword filter such as ``beam optics``,
               ``rapid cycling``, ``space charge``, ``measurement``, or
               ``superconducting``. Empty returns the complete guide.
    """
    return get_accelerator_source_guide(query)


def _load_coils_for_audit(coil_script: str):
    """Load the existing panel ``build_coil`` contract for audit tools."""
    path = Path(coil_script).expanduser().resolve(strict=True)
    if path.suffix.lower() != ".py":
        raise ValueError("coil_script must be a Python .py file")
    spec = importlib.util.spec_from_file_location(
        f"radia_mcp_coil_audit_{abs(hash(path))}", path
    )
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load coil script: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    builder = getattr(module, "build_coils", None)
    if builder is None:
        builder = getattr(module, "build_coil", None)
    if not callable(builder):
        raise ValueError(
            "coil_script must define build_coil() or build_coils()"
        )
    return builder()


@mcp.tool()
def electromagnet_coil_yoke_clearance_audit(
    coil_script: str,
    yoke_step: str,
    minimum_clearance: float = 0.0,
    intersection_volume_tolerance: float = 1.0e-15,
    fail_on_error: bool = False,
) -> dict:
    """Reject coil/yoke overlap and insufficient manufacturing clearance.

    The coil script uses the same trusted local ``build_coil() ->
    CoilBuilder`` contract as the Electromagnet application. It may instead
    define ``build_coils()`` and return multiple CoilBuilder objects. The yoke
    is read from STEP, intersected with the swept copper solid, and checked
    before any field solve is started.
    """
    from radia.coil_builder import audit_coil_yoke_clearance

    coils = _load_coils_for_audit(coil_script)
    yoke_path = Path(yoke_step).expanduser().resolve(strict=True)
    report = audit_coil_yoke_clearance(
        coils,
        yoke_path,
        minimum_clearance=minimum_clearance,
        intersection_volume_tolerance=intersection_volume_tolerance,
    )
    if fail_on_error and not report["passed"]:
        raise RuntimeError(
            "Coil/yoke clearance audit failed: "
            f"intersection_volume={report['intersection_volume']:.9g}, "
            f"measured_clearance={report['measured_clearance']:.9g}, "
            f"required_clearance={report['minimum_clearance']:.9g}"
        )
    return report


@mcp.tool()
def electromagnet_coil_field_audit(
    coil_script: str,
    sample_points: list[list[float]],
    n_arc: int = 200,
    arc_max_segment_length: float | None = None,
    relative_tolerance: float = 0.02,
    absolute_tolerance_T: float = 1.0e-9,
    closure_tolerance: float = 1.0e-9,
    fail_on_error: bool = False,
) -> dict:
    """Cross-check CoilBuilder solid-current and FE filament field sources.

    Sample points must lie in the beam/field observation region, outside the
    conductor. This also rejects open current paths. It is intended as a
    pre-solve source-model gate, not as a replacement for the FE solve.
    """
    from radia.coil_builder import audit_coil_field_consistency

    coils = _load_coils_for_audit(coil_script)
    report = audit_coil_field_consistency(
        coils,
        sample_points,
        n_arc=n_arc,
        arc_max_segment_length=arc_max_segment_length,
        relative_tolerance=relative_tolerance,
        absolute_tolerance_T=absolute_tolerance_T,
        closure_tolerance=closure_tolerance,
    )
    if fail_on_error and not report["passed"]:
        raise RuntimeError(
            "Coil field audit failed: "
            f"closed={report['closed']}, "
            f"max_relative_error={report['max_relative_error']:.9g}"
        )
    return report


# ============================================================
# MCP Prompts
# ============================================================

@mcp.prompt()
def new_electromagnet_simulation(magnet_type: str = "dipole") -> str:
    """Set up a new accelerator electromagnet simulation."""
    type_lower = magnet_type.strip().lower()

    base = (
        f"Set up an accelerator {magnet_type} magnet simulation.\n\n"
        "Use the electromagnet_usage tool for detailed documentation.\n\n"
        "General workflow:\n"
        "0. Translate lattice requirements through magnetic rigidity\n"
        "   - electromagnet_usage('accelerator_fundamentals')\n"
        "   - electromagnet_usage('beam_optics_contract')\n"
        "1. Define coil with CoilBuilder (no coil mesh needed)\n"
        "   - electromagnet_usage('coilbuilder') for API reference\n"
        "2. Export the yoke STEP and run the coil/yoke clearance audit\n"
        "   - fail on overlap or insufficient manufacturing clearance\n"
        "3. Create yoke + air + Kelvin domain in Cubit\n"
        "   - export netgen 'model.vol' order 2 overwrite\n"
        "4. NGSolve FEM solve (Omega-reduced scalar potential)\n"
        "   - electromagnet_usage('kelvin_workflow') for full code\n"
        "5. Visualize with GmshPostExport + coil STEP overlay\n"
        "6. Extract field harmonics at reference radius\n"
        "   - electromagnet_usage('harmonics') for FFT method\n\n"
    )

    if type_lower == "dipole":
        base += (
            "Dipole-specific notes:\n"
            "- 2-fold symmetry: allowed harmonics b_1, b_3, b_5, ...\n"
            "- Quarter model with IMA '+x-z' (Bz parallel to X, perp to Z)\n"
            "- Racetrack coil: add_straight + add_arc(180 deg)\n"
            "- Gap field estimate: B ~ mu_0 * NI / gap\n"
        )
    elif type_lower == "quadrupole":
        base += (
            "Quadrupole-specific notes:\n"
            "- 4-fold symmetry: allowed harmonics b_2, b_6, b_10, ...\n"
            "- 1/8 model with IMA (octant symmetry)\n"
            "- Saddle coils: 4 coils at +/-45 deg\n"
            "- Gradient: G = dBy/dx [T/m]\n"
        )
    elif type_lower == "sextupole":
        base += (
            "Sextupole-specific notes:\n"
            "- 6-fold symmetry: allowed harmonics b_3, b_9, b_15, ...\n"
            "- 1/12 model with IMA\n"
            "- 6 coils at 30 deg intervals\n"
        )
    else:
        base += (
            f"Custom magnet type '{magnet_type}':\n"
            "- Determine symmetry order and allowed harmonics\n"
            "- Choose IMA signs based on dominant field component\n"
            "  electromagnet_usage('ima') for sign selection rules\n"
            "- For nonlinear iron: electromagnet_usage('hantila')\n"
            "- For hysteresis: electromagnet_usage('hysteresis')\n"
        )

    return base


# ============================================================
# Entry point
# ============================================================



register_status_tool(
    mcp,
    server_name='mcp-server-electromagnet',
    description='Accelerator electromagnet: CoilBuilder, Hantila, Play/Energy hysteresis',
    subpackage='radia_mcp.electromagnet',
    related_servers=["motor", "accelerator", "magnetic-materials"],
    optional_deps=["radia", "ngsolve"],
)


register_topics_tool(
    mcp,
    server_name='mcp-server-electromagnet',
    topics=TOPICS,
)


def main():
    if "--selftest" in sys.argv:
        print("Electromagnet MCP server self-test:")
        topics = [
            "overview", "coilbuilder", "kelvin_workflow",
            "hantila", "hysteresis", "ima", "harmonics",
            "clebsch_hodograph",
            "accelerator_fundamentals", "beam_optics_contract",
            "accelerator_magnet_types", "accelerator_magnet_design",
            "rapid_cycling_magnets", "superconducting_accelerator_magnets",
            "accelerator_magnet_measurement", "accelerator_model_boundaries",
            "accelerator_sources",
        ]
        for t in topics:
            result = electromagnet_usage(t)
            print(f"  electromagnet_usage('{t}'): {len(result)} chars")
            assert len(result) > 100, f"Topic '{t}' too short"
        # Test prompt
        prompt = new_electromagnet_simulation("dipole")
        print(f"  new_electromagnet_simulation('dipole'): {len(prompt)} chars")
        assert "CoilBuilder" in prompt
        sources = electromagnet_accelerator_sources("rapid cycling")
        print(f"  electromagnet_accelerator_sources('rapid cycling'): {len(sources)} chars")
        assert "High-Intensity Proton Synchrotrons" in sources
        print("  PASSED")
        return

    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
