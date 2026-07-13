"""
build123d MCP Server — CAE-Oriented CAD Modeling

Provides tools for:
- build123d API reference (primitives, operations, export/import)
- CAE-specific guidelines (clean geometry, mesh quality)
- Script execution with geometry validation
- STEP/BREP export for Netgen and Cubit pipelines

Usage:
    mcp-server-build123d              # Start MCP server (stdio transport)
    mcp-server-build123d --selftest   # Run self-test
"""

import sys
import json
import logging
import traceback
from pathlib import Path
from io import StringIO

# Silence build123d INFO logging before it floods stderr.
# FastMCP uses stdio transport (stdin/stdout = MCP protocol frames,
# stderr = logs). build123d 0.10.x emits hundreds of INFO lines per
# Text/extrude/Sketch call, overflowing the OS pipe buffer (~64 KB)
# and stalling the MCP tool-call round-trip. Users (LAB + 100号機)
# observed execute_build123d timing out without a visible error.
# Setting the level to WARNING restores normal throughput.
logging.getLogger("build123d").setLevel(logging.WARNING)

# Eager-import build123d in the main thread. OCCT (the underlying C++
# kernel) is not safe to import for the first time from an arbitrary
# worker thread — doing so can deadlock. Loading it here means all
# subsequent `from build123d import *` inside execute_build123d hits the
# sys.modules cache without triggering OCCT init on the worker thread.
try:
    import build123d as _preload_build123d  # noqa: F401
except ImportError:
    pass

from mcp.server.fastmcp import FastMCP

from .build123d_knowledge import get_build123d_documentation
from .mass_topology_diagnosis_gate import cross_kernel_mass_topology_diagnosis_gate as _cross_kernel_mass_topology_diagnosis_gate, upstream_source_external_cad_contract_gate as _upstream_source_external_cad_contract_gate
from .dual_api_prismatic_gate import dual_api_prismatic_pattern_gate as _dual_api_prismatic_pattern_gate, dual_api_source_replay_gate as _dual_api_source_replay_gate
from .drafted_housing_gate import drafted_housing_cross_kernel_gate as _drafted_housing_cross_kernel_gate, drafted_housing_source_replay_gate as _drafted_housing_source_replay_gate
from .jointed_assembly_gate import jointed_assembly_heal_invariance_gate as _jointed_assembly_heal_invariance_gate, jointed_assembly_step_closure_gate as _jointed_assembly_step_closure_gate, jointed_assembly_source_replay_gate as _jointed_assembly_source_replay_gate
from .nested_assembly_volume_gate import nested_assembly_volume_gate as _nested_assembly_volume_gate, stud_wall_source_replay_gate as _stud_wall_source_replay_gate
from .patterned_compound_gate import patterned_compound_translation_gate as _patterned_compound_translation_gate, wrap_faces_rotational_source_replay_gate as _wrap_faces_rotational_source_replay_gate
from .reflection_handoff_gate import (
    build123d_heat_exchanger_source_recovery_gate as _build123d_heat_exchanger_source_recovery_gate,
    build123d_reflection_rotation_handoff_gate as _build123d_reflection_rotation_handoff_gate,
)
from .curved_shell_step_gate import (
    build123d_curved_shell_step_semantics_gate as _build123d_curved_shell_step_semantics_gate,
    build123d_tea_cup_source_contract_gate as _build123d_tea_cup_source_contract_gate,
    build123d_vase_external_solid_contract_gate as _build123d_vase_external_solid_contract_gate,
)
from .repeated_cavity_gate import (
    build123d_repeated_cavity_dual_api_gate as _build123d_repeated_cavity_dual_api_gate,
    build123d_repeated_cavity_source_replay_gate as _build123d_repeated_cavity_source_replay_gate,
)
from .faceted_edit_gate import (
    build123d_faceted_edit_portability_gate as _build123d_faceted_edit_portability_gate,
    build123d_faceted_source_replay_gate as _build123d_faceted_source_replay_gate,
)
from .lofted_shell_handoff_gate import (
    build123d_loft_example_source_replay_gate as _build123d_loft_example_source_replay_gate,
    build123d_lofted_shell_handoff_gate as _build123d_lofted_shell_handoff_gate,
)
from .rules import ALL_RULES as _B3D_LINT_RULES
from ..common import failure_log as _fl, register_status_tool
from ..common import web_docs as _wd
from ..common import examples as _ex

mcp = FastMCP("mcp-server-build123d")


@mcp.tool()
def build123d_lofted_shell_handoff_gate(summary_json: str) -> str:
    """Gate a bounded lofted-shell CAD handoff without solver overclaim."""
    try:
        result = _build123d_lofted_shell_handoff_gate(json.loads(summary_json))
    except (json.JSONDecodeError, TypeError, ValueError, KeyError) as exc:
        result = {
            "policy": "build123d_lofted_shell_handoff_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def build123d_loft_example_source_replay_gate(summary_json: str) -> str:
    """Gate the immutable upstream loft source and headless CAD replay."""
    try:
        result = _build123d_loft_example_source_replay_gate(
            json.loads(summary_json)
        )
    except (json.JSONDecodeError, TypeError, ValueError, KeyError) as exc:
        result = {
            "policy": "build123d_loft_example_source_replay_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)

@mcp.tool()
def build123d_faceted_edit_portability_gate(summary_json: str) -> str:
    """Separate faceted CAD portability from downstream mesh readiness."""
    try: result=_build123d_faceted_edit_portability_gate(json.loads(summary_json))
    except (json.JSONDecodeError,TypeError,ValueError,KeyError) as exc: result={"policy":"build123d_faceted_edit_portability_gate_v1","status":"invalid_input","error":str(exc)}
    return json.dumps(result,indent=2,sort_keys=True)

@mcp.tool()
def build123d_faceted_source_replay_gate(summary_json: str) -> str:
    """Gate tagged source, dependent STL, viewer stub, and external replay."""
    try: result=_build123d_faceted_source_replay_gate(json.loads(summary_json))
    except (json.JSONDecodeError,TypeError,ValueError,KeyError) as exc: result={"policy":"build123d_faceted_source_replay_gate_v1","status":"invalid_input","error":str(exc)}
    return json.dumps(result,indent=2,sort_keys=True)

@mcp.tool()
def build123d_cross_kernel_mass_topology_diagnosis_gate(summary_json: str) -> str:
    """Diagnose STEP portability while separating evidence quality from acceptance."""
    try: result=_cross_kernel_mass_topology_diagnosis_gate(json.loads(summary_json))
    except (json.JSONDecodeError,TypeError,ValueError,KeyError) as exc: result={"policy":"build123d_cross_kernel_mass_topology_diagnosis_gate_v1","status":"invalid_input","error":str(exc)}
    return json.dumps(result,indent=2,sort_keys=True)

@mcp.tool()
def build123d_upstream_source_external_cad_contract_gate(summary_json: str) -> str:
    """Gate immutable upstream execution and an explicit external-CAD decision."""
    try: result=_upstream_source_external_cad_contract_gate(json.loads(summary_json))
    except (json.JSONDecodeError,TypeError,ValueError) as exc: result={"policy":"build123d_upstream_source_external_cad_contract_gate_v1","status":"invalid_input","error":str(exc)}
    return json.dumps(result,indent=2,sort_keys=True)

@mcp.tool()
def build123d_dual_api_prismatic_pattern_gate(summary_json: str) -> str:
    """Gate native dual-API parity separately from external STEP-kernel bias."""
    try: result=_dual_api_prismatic_pattern_gate(json.loads(summary_json))
    except (json.JSONDecodeError,TypeError,ValueError,KeyError) as exc: result={"policy":"build123d_dual_api_prismatic_pattern_gate_v1","status":"invalid_input","error":str(exc)}
    return json.dumps(result,indent=2,sort_keys=True)

@mcp.tool()
def build123d_dual_api_source_replay_gate(summary_json: str) -> str:
    """Gate immutable upstream dual-API execution and headless CAD replay."""
    try: result=_dual_api_source_replay_gate(json.loads(summary_json))
    except (json.JSONDecodeError,TypeError,ValueError,KeyError) as exc: result={"policy":"build123d_dual_api_source_replay_gate_v1","status":"invalid_input","error":str(exc)}
    return json.dumps(result,indent=2,sort_keys=True)

@mcp.tool()
def build123d_drafted_housing_cross_kernel_gate(summary_json: str) -> str:
    """Gate drafted housing mass/topology across B-rep, STEP, Cubit, and Gmsh."""
    try: result=_drafted_housing_cross_kernel_gate(json.loads(summary_json))
    except (json.JSONDecodeError,TypeError,ValueError,KeyError) as exc: result={"policy":"build123d_drafted_housing_cross_kernel_gate_v1","status":"invalid_input","error":str(exc)}
    return json.dumps(result,indent=2,sort_keys=True)

@mcp.tool()
def build123d_drafted_housing_source_replay_gate(summary_json: str) -> str:
    """Gate tagged draft/fillet/hole source and headless mesh-companion replay."""
    try: result=_drafted_housing_source_replay_gate(json.loads(summary_json))
    except (json.JSONDecodeError,TypeError,ValueError,KeyError) as exc: result={"policy":"build123d_drafted_housing_source_replay_gate_v1","status":"invalid_input","error":str(exc)}
    return json.dumps(result,indent=2,sort_keys=True)

@mcp.tool()
def build123d_jointed_assembly_step_closure_gate(summary_json: str) -> str:
    """Diagnose a component-level solid closure loss in a jointed STEP assembly."""
    try: result=_jointed_assembly_step_closure_gate(json.loads(summary_json))
    except (json.JSONDecodeError,TypeError,ValueError,KeyError) as exc: result={"policy":"build123d_jointed_assembly_step_closure_gate_v1","status":"invalid_input","error":str(exc)}
    return json.dumps(result,indent=2,sort_keys=True)

@mcp.tool()
def build123d_jointed_assembly_heal_invariance_gate(summary_json: str) -> str:
    """Verify that STEP solid closure loss persists across heal/noheal imports."""
    try: result=_jointed_assembly_heal_invariance_gate(json.loads(summary_json))
    except (json.JSONDecodeError,TypeError,ValueError,KeyError) as exc: result={"policy":"build123d_jointed_assembly_heal_invariance_gate_v1","status":"invalid_input","error":str(exc)}
    return json.dumps(result,indent=2,sort_keys=True)

@mcp.tool()
def build123d_jointed_assembly_source_replay_gate(summary_json: str) -> str:
    """Gate immutable source, joint graph, and headless external-CAD evidence."""
    try: result=_jointed_assembly_source_replay_gate(json.loads(summary_json))
    except (json.JSONDecodeError,TypeError,ValueError,KeyError) as exc: result={"policy":"build123d_jointed_assembly_source_replay_gate_v1","status":"invalid_input","error":str(exc)}
    return json.dumps(result,indent=2,sort_keys=True)

@mcp.tool()
def build123d_nested_assembly_volume_gate(summary_json: str) -> str:
    """Distinguish a zero parent Compound from an empty CAD handoff."""
    try: result=_nested_assembly_volume_gate(json.loads(summary_json))
    except (json.JSONDecodeError,TypeError,ValueError,KeyError) as exc: result={"policy":"build123d_nested_assembly_volume_gate_v1","status":"invalid_input","error":str(exc)}
    return json.dumps(result,indent=2,sort_keys=True)

@mcp.tool()
def build123d_stud_wall_source_replay_gate(summary_json: str) -> str:
    """Gate exact stud-wall source, RigidJoints, and headless CAD replay."""
    try: result=_stud_wall_source_replay_gate(json.loads(summary_json))
    except (json.JSONDecodeError,TypeError,ValueError,KeyError) as exc: result={"policy":"build123d_stud_wall_source_replay_gate_v1","status":"invalid_input","error":str(exc)}
    return json.dumps(result,indent=2,sort_keys=True)

@mcp.tool()
def build123d_patterned_compound_translation_gate(summary_json: str) -> str:
    """Diagnose dominant curved-body STEP bias without solver-ready overclaim."""
    try: result=_patterned_compound_translation_gate(json.loads(summary_json))
    except (json.JSONDecodeError,TypeError,ValueError,KeyError) as exc: result={"policy":"build123d_patterned_compound_translation_gate_v1","status":"invalid_input","error":str(exc)}
    return json.dumps(result,indent=2,sort_keys=True)

@mcp.tool()
def build123d_wrap_faces_rotational_source_replay_gate(summary_json: str) -> str:
    """Gate immutable wrap_faces, thicken, and rotational-pattern replay."""
    try: result=_wrap_faces_rotational_source_replay_gate(json.loads(summary_json))
    except (json.JSONDecodeError,TypeError,ValueError,KeyError) as exc: result={"policy":"build123d_wrap_faces_rotational_source_replay_gate_v1","status":"invalid_input","error":str(exc)}
    return json.dumps(result,indent=2,sort_keys=True)


@mcp.tool()
def build123d_reflection_rotation_handoff_gate(summary_json: str) -> str:
    """Gate reflection failures and a proper-rotation two-body STEP handoff."""
    try:
        result = _build123d_reflection_rotation_handoff_gate(json.loads(summary_json))
    except (json.JSONDecodeError, TypeError, ValueError, KeyError) as exc:
        result = {
            "policy": "build123d_reflection_rotation_handoff_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def build123d_heat_exchanger_source_recovery_gate(summary_json: str) -> str:
    """Gate the upstream heat-exchanger replay and rotation recovery."""
    try:
        result = _build123d_heat_exchanger_source_recovery_gate(
            json.loads(summary_json)
        )
    except (json.JSONDecodeError, TypeError, ValueError, KeyError) as exc:
        result = {
            "policy": "build123d_heat_exchanger_source_recovery_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def build123d_curved_shell_step_semantics_gate(summary_json: str) -> str:
    """Diagnose topology-preserving curved STEP mass loss across CAD kernels."""
    try:
        result = _build123d_curved_shell_step_semantics_gate(json.loads(summary_json))
    except (json.JSONDecodeError, TypeError, ValueError, KeyError) as exc:
        result = {
            "policy": "build123d_curved_shell_step_semantics_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def build123d_tea_cup_source_contract_gate(summary_json: str) -> str:
    """Gate the upstream tea-cup source and headless portability diagnosis."""
    try:
        result = _build123d_tea_cup_source_contract_gate(json.loads(summary_json))
    except (json.JSONDecodeError, TypeError, ValueError, KeyError) as exc:
        result = {
            "policy": "build123d_tea_cup_source_contract_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def build123d_vase_external_solid_contract_gate(summary_json: str) -> str:
    """Gate an exact vase replay and reject zero-volume external solids."""
    try:
        result = _build123d_vase_external_solid_contract_gate(
            json.loads(summary_json)
        )
    except (json.JSONDecodeError, TypeError, ValueError, KeyError) as exc:
        result = {
            "policy": "build123d_vase_external_solid_contract_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def build123d_repeated_cavity_dual_api_gate(summary_json: str) -> str:
    """Gate dual APIs and four STEP imports for a repeated-feature cavity solid."""
    try:
        result = _build123d_repeated_cavity_dual_api_gate(json.loads(summary_json))
    except (json.JSONDecodeError, TypeError, ValueError, KeyError) as exc:
        result = {
            "policy": "build123d_repeated_cavity_dual_api_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def build123d_repeated_cavity_source_replay_gate(summary_json: str) -> str:
    """Gate immutable dual sources, STEP identities, and headless CAD replay."""
    try:
        result = _build123d_repeated_cavity_source_replay_gate(json.loads(summary_json))
    except (json.JSONDecodeError, TypeError, ValueError, KeyError) as exc:
        result = {
            "policy": "build123d_repeated_cavity_source_replay_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


# ============================================================
# Failure heuristics
# ============================================================

_B3D_ERROR_HEURISTICS: list[tuple[str, str]] = [
    ("nameerror",                "Undefined name — check imports. build123d exposes everything via `from build123d import *` (already prepended)."),
    ("attributeerror",           "Wrong object type. Use `build123d_lookup(query=...)` to find the right method."),
    ("typeerror",                "Argument type mismatch. Consult `build123d_usage(topic='api')` or `build123d_web_docs(query='<symbol>')`."),
    ("brep_api_makeshape failed", "OCCT builder failed — likely degenerate geometry (zero-length edge, self-intersection). Simplify inputs."),
    ("standard_constructionerror","Degenerate construction (radius 0, collinear points, ...). Check numeric inputs."),
    ("no such mode",             "Mode builder wrong (BuildPart vs BuildSketch vs BuildLine). See `build123d_usage(topic='builders')`."),
    ("plane",                    "Plane/workplane mismatch. Use `Plane.XY` / `Plane.XZ` explicitly or pass `Location`."),
    ("fillet",                   "Fillet radius too large vs adjacent edge. Reduce radius or pick fewer edges."),
    ("chamfer",                  "Chamfer length too large. Reduce or pick non-adjacent edges."),
    ("import_step",              "STEP import failed. Check the path exists and the file is not corrupt."),
]


def _build123d_failure_hint(tb: str) -> str:
    """Pick a short hint from common build123d traceback patterns."""
    blob = tb.lower()
    for needle, advice in _B3D_ERROR_HEURISTICS:
        if needle in blob:
            return advice
    return ("No specific heuristic matched. Try "
            "`build123d_lookup(query=...)` (bundled docs + failure log) or "
            "`build123d_web_docs(query=...)` (live readthedocs).")


# ============================================================
# MCP Tools
# ============================================================

@mcp.tool()
def build123d_usage(topic: str = "overview") -> str:
    """
    Get build123d CAD modeling documentation for CAE workflows.

    build123d is a Python-native parametric CAD library on the OCCT kernel.
    This server focuses on CAE (FEM/BEM) geometry creation, not 3D printing.

    Args:
        topic: Documentation topic. Options:
            "overview"        - What build123d is, CAE pipeline, safe subset
            "parametric_library" - Tested ops (wedge/tube/racetrack/arrays/mirror/assembly) + EM
                                archetypes (Halbach ring, C-core, magnets, solenoid): modeling.py + archetypes.py
            "lab_policy"      - Role split vs Cubit (tet vs hex), translation guidance
            "cubit_rosetta"   - Cubit `.jou` verb ↔ build123d mapping table
            "primitives_3d"   - Box, Cylinder, Cone, Sphere, Torus, Wedge + boolean
            "primitives_2d"   - Circle, Rectangle, Polygon, etc. (sketch objects)
            "operations"      - extrude, revolve, sweep, loft, fillet, mirror, split
            "curves"          - Line, Polyline, Spline, Arc types, Helix, Bezier
            "export_import"   - STEP, BREP, STL, glTF, SVG + CAE pipeline examples
            "topology"        - faces/edges/vertices queries, labels, selectors
            "cae_guidelines"  - Clean geometry rules, quality checks, label conventions
            "examples"        - IH coil, E-core transformer, dipole yoke, parametric
            "examples_intro"  - 36 introductory examples (Box/Cylinder -> loft/revolve/sweep)
            "examples_gallery"- 21 gallery examples (Benchy, Heat Exchanger, Vase, ...)
            "examples_lab_patterns" - Lab-specific: IH/Halbach/dipole/PEEC/maglev archetypes
            "coil_modeling"   - PEEC filament extraction from CAD
            "all"             - Complete documentation
    """
    return get_build123d_documentation(topic)


@mcp.tool()
def build123d_motor_housing_thermal_reference(
    inner_radius: float = 42.0,
    outer_radius: float = 50.0,
    length: float = 90.0,
    fin_count: int = 8,
    fin_height: float = 12.0,
    fin_thickness: float = 3.0,
) -> str:
    """Return analytic volume/area/body-count data for a finned motor housing.

    The result is a multi-body sleeve + radial-fin reference row that can feed
    ``build123d_volume_crosscheck`` and a Cubit STEP sidecar. Dimensions use
    millimetres; volume and area are returned in mm^3 and mm^2.
    """
    from .modeling import motor_housing_radial_fin_reference_row

    return json.dumps(
        motor_housing_radial_fin_reference_row(
            inner_radius,
            outer_radius,
            length,
            fin_count,
            fin_height,
            fin_thickness,
        ),
        indent=2,
        sort_keys=True,
    )


@mcp.tool()
def build123d_volume_crosscheck(
    reference_rows_json: str,
    measured_sets_json: str,
    rtol: float = 1.0e-5,
) -> str:
    """Compare build123d reference volumes with Cubit or external-CAD volumes.

    This is the lowest-friction CAD round-trip gate: both sides only need rows
    shaped like ``{"name": "...", "volume": 123.0}``.  ``measured_sets_json``
    may be a mapping such as ``{"cubit": [...], "external_cad": [...]}`` or a
    list of ``{"source": label, "rows": [...]}`` objects.  Use this before
    spending solver time on STEP files that may have imported at the wrong unit
    scale or lost a boolean feature.
    """

    try:
        reference_rows = json.loads(reference_rows_json)
        measured_sets = json.loads(measured_sets_json)
    except json.JSONDecodeError as e:
        return json.dumps({
            "status": "error",
            "stage": "parse_json",
            "error": f"{type(e).__name__}: {e}",
        }, indent=2)

    try:
        from .modeling import shape_volume_crosscheck_summary
        summary = shape_volume_crosscheck_summary(
            reference_rows,
            measured_sets,
            rtol=rtol,
        )
    except Exception as e:  # noqa: BLE001
        return json.dumps({
            "status": "error",
            "stage": "volume_crosscheck",
            "error": f"{type(e).__name__}: {e}",
        }, indent=2)

    return json.dumps(summary, indent=2)


@mcp.tool()
def build123d_brep_mass_topology_roundtrip_gate(
    reference_json: str,
    measured_rows_json: str,
    expected_volume: float | None = None,
    expected_volume_rtol: float = 1.0e-6,
    volume_rtol: float = 1.0e-4,
    area_rtol: float = 2.0e-3,
    bbox_atol: float = 0.11,
    centroid_atol: float = 2.0e-4,
) -> str:
    """Gate CAD roundtrips by mass properties, centroid semantics and B-rep Euler topology."""

    try:
        from .brep_roundtrip_gate import brep_mass_topology_roundtrip_gate

        result = brep_mass_topology_roundtrip_gate(
            json.loads(reference_json),
            json.loads(measured_rows_json),
            expected_volume=expected_volume,
            expected_volume_rtol=expected_volume_rtol,
            volume_rtol=volume_rtol,
            area_rtol=area_rtol,
            bbox_atol=bbox_atol,
            centroid_atol=centroid_atol,
        )
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "policy": "build123d_brep_mass_topology_roundtrip_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def build123d_platonic_solid_family_gate(summary_json: str) -> str:
    """Gate all five Platonic solids by topology, analytic volume and CAD replay."""
    try:
        from .platonic_family_gate import platonic_solid_family_gate

        result = platonic_solid_family_gate(json.loads(summary_json))
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "policy": "build123d_platonic_solid_family_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def build123d_dual_api_perforated_board_gate(summary_json: str) -> str:
    """Gate equivalent Builder/Algebra perforated boards through two CAD imports."""
    try:
        from .dual_api_board_gate import dual_api_perforated_board_gate

        result = dual_api_perforated_board_gate(json.loads(summary_json))
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "policy": "build123d_dual_api_perforated_board_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def build123d_upstream_example_roundtrip_gate(
    result_json: str,
    mass_property_rtol: float = 1.0e-12,
    centroid_atol: float = 1.0e-12,
) -> str:
    """Gate source identity and STEP self-roundtrip for an upstream build123d example."""

    try:
        from .external_cad_gate import build123d_upstream_example_roundtrip_gate as gate

        result = gate(
            json.loads(result_json),
            mass_property_rtol=mass_property_rtol,
            centroid_atol=centroid_atol,
        )
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "policy": "build123d_upstream_example_roundtrip_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def build123d_external_cad_mass_topology_gate(
    reference_json: str,
    external_json: str,
    volume_rtol: float = 2.0e-6,
    area_rtol: float = 1.0e-10,
    bbox_atol: float = 1.0e-10,
    centroid_atol: float = 1.0e-6,
) -> str:
    """Crosscheck two CAD kernels without confusing entity centers with mass centroids."""

    try:
        from .external_cad_gate import external_cad_mass_topology_crosscheck_gate

        result = external_cad_mass_topology_crosscheck_gate(
            json.loads(reference_json),
            json.loads(external_json),
            volume_rtol=volume_rtol,
            area_rtol=area_rtol,
            bbox_atol=bbox_atol,
            centroid_atol=centroid_atol,
        )
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "policy": "external_cad_mass_topology_crosscheck_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def build123d_step_portability_diagnosis_gate(
    result_json: str,
    self_roundtrip_rtol: float = 1.0e-12,
    external_volume_rtol: float = 2.0e-6,
    import_mode_rtol: float = 1.0e-12,
) -> str:
    """Diagnose whether STEP mass loss occurs in export or external import.

    The structured input records native and same-kernel STEP volumes plus
    independent ``heal`` and ``noheal`` import rows.
    """

    try:
        from .external_cad_gate import step_portability_diagnosis_gate

        result = step_portability_diagnosis_gate(
            json.loads(result_json),
            self_roundtrip_rtol=self_roundtrip_rtol,
            external_volume_rtol=external_volume_rtol,
            import_mode_rtol=import_mode_rtol,
        )
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "policy": "build123d_step_portability_diagnosis_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def build123d_curved_step_topology_crosscheck_gate(
    result_json: str,
    self_mass_rtol: float = 1.0e-7,
    external_volume_rtol: float = 1.0e-4,
    import_mode_rtol: float = 1.0e-10,
    bbox_atol: float = 1.0e-9,
) -> str:
    """Gate curved STEP mass and exact topology across independent imports."""
    try:
        from .external_cad_gate import curved_step_topology_crosscheck_gate

        result = curved_step_topology_crosscheck_gate(
            json.loads(result_json),
            self_mass_rtol=self_mass_rtol,
            external_volume_rtol=external_volume_rtol,
            import_mode_rtol=import_mode_rtol,
            bbox_atol=bbox_atol,
        )
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "policy": "build123d_curved_step_topology_crosscheck_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def build123d_path_sweep_handoff_gate(
    result_json: str,
    path_length_rtol: float = 1.0e-12,
    external_volume_rtol: float = 5.0e-5,
    external_area_rtol: float = 1.0e-10,
    bbox_atol: float = 1.0e-10,
) -> str:
    """Gate a curved build123d sweep through analytic path and STEP/CAD checks.

    The analytic ``section area * path length`` value is reported only as a
    diagnostic. The finished trimmed solid must agree with an independent CAD
    kernel in volume, area, bounding box, STEP digest, and closed topology.
    """

    try:
        from .path_sweep_gate import build123d_path_sweep_handoff_gate as gate

        result = gate(
            json.loads(result_json),
            path_length_rtol=path_length_rtol,
            external_volume_rtol=external_volume_rtol,
            external_area_rtol=external_area_rtol,
            bbox_atol=bbox_atol,
        )
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "policy": "build123d_path_sweep_handoff_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def build123d_path_sweep_source_contract_gate(result_json: str) -> str:
    """Gate the source-native build123d sweep idiom and ``is_valid`` API form."""

    try:
        from .path_sweep_gate import build123d_path_sweep_source_contract_gate as gate

        result = gate(json.loads(result_json))
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "policy": "build123d_path_sweep_source_contract_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def build123d_perforated_prism_roundtrip_gate(
    reference_volume: float,
    imported_volume: float,
    hole_count: int,
    hole_side_count: int,
    imported_surface_count: int,
    imported_body_count: int = 1,
    outer_side_count: int = 4,
    volume_rtol: float = 1.0e-9,
) -> str:
    """Check STEP roundtrip volume and through-hole boundary topology."""

    try:
        from .modeling import shape_perforated_prism_roundtrip_gate

        result = shape_perforated_prism_roundtrip_gate(
            reference_volume=reference_volume,
            imported_volume=imported_volume,
            hole_count=hole_count,
            hole_side_count=hole_side_count,
            imported_surface_count=imported_surface_count,
            imported_body_count=imported_body_count,
            outer_side_count=outer_side_count,
            volume_rtol=volume_rtol,
        )
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "build123d_perforated_prism_roundtrip_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2)


@mcp.tool()
def build123d_volume_crosscheck_with_units(
    reference_rows_json: str,
    measured_sets_json: str,
    target_volume_unit: str = "mm^3",
    rtol: float = 1.0e-5,
) -> str:
    """Normalize explicit cubic units before comparing CAD volumes.

    Every row must carry ``volume_unit``. This prevents a numerically plausible
    comparison from silently mixing build123d's usual millimetre model units
    with external CAD exports reported in cubic centimetres or cubic metres.
    """

    factors_to_mm3 = {
        "mm^3": 1.0,
        "mm3": 1.0,
        "cm^3": 1.0e3,
        "cm3": 1.0e3,
        "m^3": 1.0e9,
        "m3": 1.0e9,
        "in^3": 25.4 ** 3,
        "in3": 25.4 ** 3,
    }

    try:
        reference_rows = json.loads(reference_rows_json)
        measured_sets = json.loads(measured_sets_json)
    except json.JSONDecodeError as exc:
        return json.dumps({
            "policy": "build123d_volume_crosscheck_with_units",
            "status": "invalid_input",
            "error": f"{type(exc).__name__}: {exc}",
        }, indent=2)

    target = str(target_volume_unit or "").strip().lower()
    if target not in factors_to_mm3:
        return json.dumps({
            "policy": "build123d_volume_crosscheck_with_units",
            "status": "invalid_input",
            "error": f"unsupported target_volume_unit: {target_volume_unit}",
        }, indent=2)

    issues: list[str] = []

    def normalize_rows(rows, source: str):
        normalized = []
        for index, row in enumerate(rows or []):
            item = dict(row)
            unit = str(item.get("volume_unit", "")).strip().lower()
            if unit not in factors_to_mm3:
                issues.append(f"{source}[{index}] missing or unsupported volume_unit")
                continue
            try:
                volume = float(item["volume"])
            except (KeyError, TypeError, ValueError):
                issues.append(f"{source}[{index}] missing numeric volume")
                continue
            item["source_volume"] = volume
            item["source_volume_unit"] = unit
            item["volume"] = volume * factors_to_mm3[unit] / factors_to_mm3[target]
            item["volume_unit"] = target
            normalized.append(item)
        return normalized

    normalized_reference = normalize_rows(reference_rows, "reference")
    if isinstance(measured_sets, dict):
        normalized_measured = {
            str(source): normalize_rows(rows, str(source))
            for source, rows in measured_sets.items()
        }
    elif isinstance(measured_sets, list):
        normalized_measured = []
        for index, record in enumerate(measured_sets):
            item = dict(record)
            source = str(item.get("source", f"source_{index}"))
            item["rows"] = normalize_rows(item.get("rows", []), source)
            normalized_measured.append(item)
    else:
        issues.append("measured_sets must be a mapping or list")
        normalized_measured = {}

    if issues:
        return json.dumps({
            "policy": "build123d_volume_crosscheck_with_units",
            "status": "needs_attention",
            "target_volume_unit": target,
            "issues": issues,
        }, indent=2)

    try:
        from .modeling import shape_volume_crosscheck_summary
        summary = shape_volume_crosscheck_summary(
            normalized_reference,
            normalized_measured,
            rtol=rtol,
        )
    except Exception as exc:  # noqa: BLE001
        return json.dumps({
            "policy": "build123d_volume_crosscheck_with_units",
            "status": "error",
            "error": f"{type(exc).__name__}: {exc}",
        }, indent=2)

    return json.dumps({
        "policy": "build123d_volume_crosscheck_with_units",
        "status": summary.get("status", "needs_attention"),
        "target_volume_unit": target,
        "normalized_reference_rows": normalized_reference,
        "normalized_measured_sets": normalized_measured,
        "crosscheck": summary,
        "issues": [],
    }, indent=2)


@mcp.tool()
def build123d_volume_crosscheck_source_coverage_gate(
    volume_summary_json: str,
    required_sources_json: str = "",
    max_allowed_volume_rel_error: float | None = None,
) -> str:
    """Require Cubit/external CAD source coverage after a volume crosscheck."""

    try:
        volume_summary = json.loads(volume_summary_json)
        required_sources = (
            None if not required_sources_json.strip()
            else json.loads(required_sources_json)
        )
    except json.JSONDecodeError as e:
        return json.dumps({
            "status": "error",
            "stage": "parse_json",
            "error": f"{type(e).__name__}: {e}",
        }, indent=2)

    try:
        from .modeling import shape_volume_crosscheck_source_coverage_gate
        kwargs = {}
        if required_sources is not None:
            kwargs["required_sources"] = required_sources
        if max_allowed_volume_rel_error is not None:
            kwargs["max_allowed_volume_rel_error"] = max_allowed_volume_rel_error
        summary = shape_volume_crosscheck_source_coverage_gate(
            volume_summary,
            **kwargs,
        )
    except Exception as e:  # noqa: BLE001
        return json.dumps({
            "status": "error",
            "stage": "volume_crosscheck_source_coverage_gate",
            "error": f"{type(e).__name__}: {e}",
        }, indent=2)

    return json.dumps(summary, indent=2)


@mcp.tool()
def build123d_volume_crosscheck_source_identity_gate(
    volume_summary_json: str,
    expected_measurement_methods_json: str = "",
    expected_body_identity_keys_json: str = "",
    expected_source_artifact_ids_json: str = "",
    expected_parameter_set_artifact_ids_json: str = "",
    expected_parameter_set_digests_json: str = "",
    expected_parameter_set_paths_json: str = "",
    expected_objective_observable_ids_json: str = "",
    expected_objective_observable_families_json: str = "",
) -> str:
    """Require source identity metadata after an external CAD volume crosscheck.

    Use this after ``build123d_volume_crosscheck`` when external rows came from
    Cubit, CST, or another CAD importer and the validation claim depends on the
    measurement method, body identity key, or source artifact id.
    """

    try:
        volume_summary = json.loads(volume_summary_json)
        expected_measurement_methods = (
            None if not expected_measurement_methods_json.strip()
            else json.loads(expected_measurement_methods_json)
        )
        expected_body_identity_keys = (
            None if not expected_body_identity_keys_json.strip()
            else json.loads(expected_body_identity_keys_json)
        )
        expected_source_artifact_ids = (
            None if not expected_source_artifact_ids_json.strip()
            else json.loads(expected_source_artifact_ids_json)
        )
        expected_parameter_set_artifact_ids = (
            None if not expected_parameter_set_artifact_ids_json.strip()
            else json.loads(expected_parameter_set_artifact_ids_json)
        )
        expected_parameter_set_digests = (
            None if not expected_parameter_set_digests_json.strip()
            else json.loads(expected_parameter_set_digests_json)
        )
        expected_parameter_set_paths = (
            None if not expected_parameter_set_paths_json.strip()
            else json.loads(expected_parameter_set_paths_json)
        )
        expected_objective_observable_ids = (
            None if not expected_objective_observable_ids_json.strip()
            else json.loads(expected_objective_observable_ids_json)
        )
        expected_objective_observable_families = (
            None if not expected_objective_observable_families_json.strip()
            else json.loads(expected_objective_observable_families_json)
        )
    except json.JSONDecodeError as e:
        return json.dumps({
            "status": "error",
            "stage": "parse_json",
            "error": f"{type(e).__name__}: {e}",
        }, indent=2)

    try:
        from .modeling import shape_volume_crosscheck_source_identity_gate
        summary = shape_volume_crosscheck_source_identity_gate(
            volume_summary,
            expected_measurement_methods=expected_measurement_methods,
            expected_body_identity_keys=expected_body_identity_keys,
            expected_source_artifact_ids=expected_source_artifact_ids,
            expected_parameter_set_artifact_ids=expected_parameter_set_artifact_ids,
            expected_parameter_set_digests=expected_parameter_set_digests,
            expected_parameter_set_paths=expected_parameter_set_paths,
            expected_objective_observable_ids=expected_objective_observable_ids,
            expected_objective_observable_families=expected_objective_observable_families,
        )
    except Exception as e:  # noqa: BLE001
        return json.dumps({
            "status": "error",
            "stage": "volume_crosscheck_source_identity_gate",
            "error": f"{type(e).__name__}: {e}",
        }, indent=2)

    return json.dumps(summary, indent=2)


@mcp.tool()
def build123d_mass_property_crosscheck(
    reference_rows_json: str,
    measured_sets_json: str,
    rtol: float = 1.0e-5,
    bbox_atol: float = 1.0e-6,
) -> str:
    """Compare build123d volume/area/bbox rows with one or more CAD sources.

    Use this stronger gate when the external CAD kernel can report rows shaped
    like ``{"name": "...", "volume": 1.0, "area": 6.0, "bounding_box": ...}``.
    It keeps volume, surface area, and optional bounding-box checks together for
    build123d -> STEP -> Cubit/CST/external-CAD round trips.
    """

    try:
        reference_rows = json.loads(reference_rows_json)
        measured_sets = json.loads(measured_sets_json)
    except json.JSONDecodeError as e:
        return json.dumps({
            "status": "error",
            "stage": "parse_json",
            "error": f"{type(e).__name__}: {e}",
        }, indent=2)

    try:
        from .modeling import shape_mass_property_crosscheck_summary
        summary = shape_mass_property_crosscheck_summary(
            reference_rows,
            measured_sets,
            rtol=rtol,
            bbox_atol=bbox_atol,
        )
    except Exception as e:  # noqa: BLE001
        return json.dumps({
            "status": "error",
            "stage": "mass_property_crosscheck",
            "error": f"{type(e).__name__}: {e}",
        }, indent=2)

    return json.dumps(summary, indent=2)


@mcp.tool()
def build123d_cad_route_source_contract(
    shape_rows_json: str,
    external_crosscheck_summary_json: str,
    required_source_groups_json: str = "",
    expected_route: str = "cubit_hex_or_mixed_path",
    expected_length_unit: str = "",
    expected_area_unit: str = "",
    expected_volume_unit: str = "",
    required_metadata_fields_json: str = "",
) -> str:
    """Gate a build123d CAD package before Cubit hex/mixed route promotion.

    ``shape_rows_json`` should contain build123d measurement rows enriched with
    ``geometry_id``, ``authoring_source`` or ``source_kind``, and
    ``mesh_route``.  ``external_crosscheck_summary_json`` is the JSON output of
    ``build123d_volume_crosscheck`` or ``build123d_mass_property_crosscheck``.
    Optionally pass ``required_source_groups_json`` as e.g.
    ``[["cubit", "coreform_cubit"], ["external_cad", "cst_import"]]``.
    Pass ``required_metadata_fields_json`` as e.g.
    ``["recipe_id", "cad_kernel", "cad_kernel_version", "script_path", "export_id"]``
    when solver-ready promotion must bind CAD values to authoring provenance.
    """

    try:
        shape_rows = json.loads(shape_rows_json)
        external_crosscheck_summary = json.loads(external_crosscheck_summary_json)
        required_source_groups = (
            None if not required_source_groups_json.strip()
            else json.loads(required_source_groups_json)
        )
        required_metadata_fields = (
            None if not required_metadata_fields_json.strip()
            else json.loads(required_metadata_fields_json)
        )
    except json.JSONDecodeError as e:
        return json.dumps({
            "status": "error",
            "stage": "parse_json",
            "error": f"{type(e).__name__}: {e}",
        }, indent=2)

    try:
        from .modeling import shape_cad_route_source_contract_gate
        kwargs = {"expected_route": expected_route}
        if expected_length_unit.strip():
            kwargs["expected_length_unit"] = expected_length_unit
        if expected_area_unit.strip():
            kwargs["expected_area_unit"] = expected_area_unit
        if expected_volume_unit.strip():
            kwargs["expected_volume_unit"] = expected_volume_unit
        if required_source_groups is not None:
            kwargs["required_source_groups"] = required_source_groups
        if required_metadata_fields is not None:
            kwargs["required_metadata_fields"] = required_metadata_fields
        summary = shape_cad_route_source_contract_gate(
            shape_rows,
            external_crosscheck_summary,
            **kwargs,
        )
    except Exception as e:  # noqa: BLE001
        return json.dumps({
            "status": "error",
            "stage": "cad_route_source_contract",
            "error": f"{type(e).__name__}: {e}",
        }, indent=2)

    return json.dumps(summary, indent=2)


@mcp.tool()
def build123d_external_cad_volume_evidence_package(
    shape_rows_json: str,
    volume_summary_json: str,
    required_sources_json: str = "",
    expected_measurement_methods_json: str = "",
    expected_body_identity_keys_json: str = "",
    expected_source_artifact_ids_json: str = "",
    expected_parameter_set_artifact_ids_json: str = "",
    expected_parameter_set_digests_json: str = "",
    expected_parameter_set_paths_json: str = "",
    expected_objective_observable_ids_json: str = "",
    expected_objective_observable_families_json: str = "",
    expected_route: str = "cubit_hex_or_mixed_path",
    expected_length_unit: str = "",
    expected_area_unit: str = "",
    expected_volume_unit: str = "",
    required_metadata_fields_json: str = "",
    max_allowed_volume_rel_error: float | None = None,
) -> str:
    """Bundle dual-source CAD volume evidence before reuse.

    ``volume_summary_json`` is the output of ``build123d_volume_crosscheck``.
    This package gate then requires source coverage, source identity metadata,
    and a build123d route contract in one reusable JSON object.
    """

    try:
        shape_rows = json.loads(shape_rows_json)
        volume_summary = json.loads(volume_summary_json)
        required_sources = (
            None if not required_sources_json.strip()
            else json.loads(required_sources_json)
        )
        expected_measurement_methods = (
            None if not expected_measurement_methods_json.strip()
            else json.loads(expected_measurement_methods_json)
        )
        expected_body_identity_keys = (
            None if not expected_body_identity_keys_json.strip()
            else json.loads(expected_body_identity_keys_json)
        )
        expected_source_artifact_ids = (
            None if not expected_source_artifact_ids_json.strip()
            else json.loads(expected_source_artifact_ids_json)
        )
        expected_parameter_set_artifact_ids = (
            None if not expected_parameter_set_artifact_ids_json.strip()
            else json.loads(expected_parameter_set_artifact_ids_json)
        )
        expected_parameter_set_digests = (
            None if not expected_parameter_set_digests_json.strip()
            else json.loads(expected_parameter_set_digests_json)
        )
        expected_parameter_set_paths = (
            None if not expected_parameter_set_paths_json.strip()
            else json.loads(expected_parameter_set_paths_json)
        )
        expected_objective_observable_ids = (
            None if not expected_objective_observable_ids_json.strip()
            else json.loads(expected_objective_observable_ids_json)
        )
        expected_objective_observable_families = (
            None if not expected_objective_observable_families_json.strip()
            else json.loads(expected_objective_observable_families_json)
        )
        required_metadata_fields = (
            None if not required_metadata_fields_json.strip()
            else json.loads(required_metadata_fields_json)
        )
    except json.JSONDecodeError as e:
        return json.dumps({
            "status": "error",
            "stage": "parse_json",
            "error": f"{type(e).__name__}: {e}",
        }, indent=2)

    try:
        from .modeling import shape_external_cad_volume_evidence_package_gate
        kwargs = {"expected_route": expected_route}
        if required_sources is not None:
            kwargs["required_sources"] = required_sources
        if max_allowed_volume_rel_error is not None:
            kwargs["max_allowed_volume_rel_error"] = max_allowed_volume_rel_error
        if expected_measurement_methods is not None:
            kwargs["expected_measurement_methods"] = expected_measurement_methods
        if expected_body_identity_keys is not None:
            kwargs["expected_body_identity_keys"] = expected_body_identity_keys
        if expected_source_artifact_ids is not None:
            kwargs["expected_source_artifact_ids"] = expected_source_artifact_ids
        if expected_parameter_set_artifact_ids is not None:
            kwargs["expected_parameter_set_artifact_ids"] = expected_parameter_set_artifact_ids
        if expected_parameter_set_digests is not None:
            kwargs["expected_parameter_set_digests"] = expected_parameter_set_digests
        if expected_parameter_set_paths is not None:
            kwargs["expected_parameter_set_paths"] = expected_parameter_set_paths
        if expected_objective_observable_ids is not None:
            kwargs["expected_objective_observable_ids"] = expected_objective_observable_ids
        if expected_objective_observable_families is not None:
            kwargs["expected_objective_observable_families"] = expected_objective_observable_families
        if expected_length_unit.strip():
            kwargs["expected_length_unit"] = expected_length_unit
        if expected_area_unit.strip():
            kwargs["expected_area_unit"] = expected_area_unit
        if expected_volume_unit.strip():
            kwargs["expected_volume_unit"] = expected_volume_unit
        if required_metadata_fields is not None:
            kwargs["required_metadata_fields"] = required_metadata_fields
        summary = shape_external_cad_volume_evidence_package_gate(
            shape_rows,
            volume_summary,
            **kwargs,
        )
    except Exception as e:  # noqa: BLE001
        return json.dumps({
            "status": "error",
            "stage": "external_cad_volume_evidence_package",
            "error": f"{type(e).__name__}: {e}",
        }, indent=2)

    return json.dumps(summary, indent=2)


@mcp.tool()
def build123d_cubit_solver_route_handoff(
    shape_rows_json: str,
    solver_route_gate_json: str,
    geometry_id_key: str = "geometry_id",
    expected_route: str = "cubit_hex_or_mixed_path",
    expected_solver_route_package_id: str = "",
    expected_solver_route_convention_schema_id: str = "",
    require_solver_route_convention_schema: bool = False,
    require_no_implicit_tetization: bool = True,
) -> str:
    """Bind build123d CAD rows to a Cubit mixed solver-route manifest gate."""

    try:
        shape_rows = json.loads(shape_rows_json)
        solver_route_gate = json.loads(solver_route_gate_json)
    except json.JSONDecodeError as e:
        return json.dumps({
            "status": "error",
            "stage": "parse_json",
            "error": f"{type(e).__name__}: {e}",
        }, indent=2)

    try:
        from .modeling import shape_cubit_solver_route_handoff_gate
        kwargs = {
            "geometry_id_key": geometry_id_key,
            "expected_route": expected_route,
            "require_no_implicit_tetization": require_no_implicit_tetization,
        }
        if expected_solver_route_package_id.strip():
            kwargs["expected_solver_route_package_id"] = expected_solver_route_package_id
        if expected_solver_route_convention_schema_id.strip():
            kwargs["expected_solver_route_convention_schema_id"] = expected_solver_route_convention_schema_id
        if require_solver_route_convention_schema:
            kwargs["require_solver_route_convention_schema"] = True
        summary = shape_cubit_solver_route_handoff_gate(
            shape_rows,
            solver_route_gate,
            **kwargs,
        )
    except Exception as e:  # noqa: BLE001
        return json.dumps({
            "status": "error",
            "stage": "cubit_solver_route_handoff",
            "error": f"{type(e).__name__}: {e}",
        }, indent=2)

    return json.dumps(summary, indent=2)


@mcp.tool()
def build123d_cubit_quality_ledger_handoff(
    shape_rows_json: str,
    quality_ledger_gate_json: str,
    geometry_id_key: str = "geometry_id",
    expected_quality_artifact_id: str = "",
    expected_quality_digest: str = "",
    expected_metric_set_id: str = "",
    expected_export_id: str = "",
    expected_mesh_artifact_id: str = "",
    expected_mesh_digest: str = "",
    expected_route: str = "cubit_hex_or_mixed_path",
    min_scaled_jacobian_threshold: float = 0.2,
    require_hex_or_mixed_route: bool = True,
) -> str:
    """Bind build123d CAD rows to a Cubit mesh-quality ledger identity gate."""

    try:
        shape_rows = json.loads(shape_rows_json)
        quality_ledger_gate = json.loads(quality_ledger_gate_json)
    except json.JSONDecodeError as e:
        return json.dumps({
            "status": "error",
            "stage": "parse_json",
            "error": f"{type(e).__name__}: {e}",
        }, indent=2)

    try:
        from .modeling import shape_cubit_quality_ledger_handoff_gate
        kwargs = {
            "geometry_id_key": geometry_id_key,
            "expected_route": expected_route,
            "min_scaled_jacobian_threshold": min_scaled_jacobian_threshold,
            "require_hex_or_mixed_route": require_hex_or_mixed_route,
        }
        if expected_quality_artifact_id.strip():
            kwargs["expected_quality_artifact_id"] = expected_quality_artifact_id
        if expected_quality_digest.strip():
            kwargs["expected_quality_digest"] = expected_quality_digest
        if expected_metric_set_id.strip():
            kwargs["expected_metric_set_id"] = expected_metric_set_id
        if expected_export_id.strip():
            kwargs["expected_export_id"] = expected_export_id
        if expected_mesh_artifact_id.strip():
            kwargs["expected_mesh_artifact_id"] = expected_mesh_artifact_id
        if expected_mesh_digest.strip():
            kwargs["expected_mesh_digest"] = expected_mesh_digest
        summary = shape_cubit_quality_ledger_handoff_gate(
            shape_rows,
            quality_ledger_gate,
            **kwargs,
        )
    except Exception as e:  # noqa: BLE001
        return json.dumps({
            "status": "error",
            "stage": "cubit_quality_ledger_handoff",
            "error": f"{type(e).__name__}: {e}",
        }, indent=2)

    return json.dumps(summary, indent=2)


@mcp.tool()
def build123d_cad_handoff_manifest(
    shape_rows_json: str,
    file_manifest_json: str,
    external_volume_summary_json: str = "",
    cubit_export_handoff_json: str = "",
    cubit_quality_handoff_json: str = "",
    cubit_quality_ledger_handoff_json: str = "",
    cubit_solver_route_handoff_json: str = "",
    required_file_kinds_json: str = "",
    expected_geometry_ids_json: str = "",
    expected_cad_output_artifact_id: str = "",
    expected_cad_output_digest: str = "",
    require_cad_output_artifact: bool = False,
    expected_cad_observable_id: str = "",
    expected_cad_observable_family: str = "",
    require_cad_observable: bool = False,
    expected_length_unit: str = "",
    expected_area_unit: str = "",
    expected_volume_unit: str = "",
    expected_measurement_convention: str = "",
    expected_measurement_postprocess_row_convention_schema_id: str = "",
    expected_measurement_component_basis_schema_id: str = "",
    require_measurement_postprocess_row_convention_schema: bool = False,
    require_measurement_component_basis_schema: bool = False,
) -> str:
    """Run the final build123d CAD handoff manifest gate from JSON inputs."""

    try:
        shape_rows = json.loads(shape_rows_json)
        file_manifest = json.loads(file_manifest_json)
        external_volume_summary = (
            None if not external_volume_summary_json.strip()
            else json.loads(external_volume_summary_json)
        )
        cubit_export_handoff = (
            None if not cubit_export_handoff_json.strip()
            else json.loads(cubit_export_handoff_json)
        )
        cubit_quality_handoff = (
            None if not cubit_quality_handoff_json.strip()
            else json.loads(cubit_quality_handoff_json)
        )
        cubit_quality_ledger_handoff = (
            None if not cubit_quality_ledger_handoff_json.strip()
            else json.loads(cubit_quality_ledger_handoff_json)
        )
        cubit_solver_route_handoff = (
            None if not cubit_solver_route_handoff_json.strip()
            else json.loads(cubit_solver_route_handoff_json)
        )
        required_file_kinds = (
            None if not required_file_kinds_json.strip()
            else json.loads(required_file_kinds_json)
        )
        expected_geometry_ids = (
            None if not expected_geometry_ids_json.strip()
            else json.loads(expected_geometry_ids_json)
        )
    except json.JSONDecodeError as e:
        return json.dumps({
            "status": "error",
            "stage": "parse_json",
            "error": f"{type(e).__name__}: {e}",
        }, indent=2)

    try:
        from .modeling import shape_cad_handoff_manifest_gate
        kwargs = {
            "file_manifest": file_manifest,
            "require_cad_output_artifact": require_cad_output_artifact,
            "require_cad_observable": require_cad_observable,
            "require_measurement_postprocess_row_convention_schema": (
                require_measurement_postprocess_row_convention_schema
            ),
            "require_measurement_component_basis_schema": (
                require_measurement_component_basis_schema
            ),
        }
        if external_volume_summary is not None:
            kwargs["external_volume_summary"] = external_volume_summary
        if cubit_export_handoff is not None:
            kwargs["cubit_export_handoff"] = cubit_export_handoff
        if cubit_quality_handoff is not None:
            kwargs["cubit_quality_handoff"] = cubit_quality_handoff
        if cubit_quality_ledger_handoff is not None:
            kwargs["cubit_quality_ledger_handoff"] = cubit_quality_ledger_handoff
        if cubit_solver_route_handoff is not None:
            kwargs["cubit_solver_route_handoff"] = cubit_solver_route_handoff
        if required_file_kinds is not None:
            kwargs["required_file_kinds"] = required_file_kinds
        if expected_geometry_ids is not None:
            kwargs["expected_geometry_ids"] = expected_geometry_ids
        if expected_cad_output_artifact_id.strip():
            kwargs["expected_cad_output_artifact_id"] = expected_cad_output_artifact_id
        if expected_cad_output_digest.strip():
            kwargs["expected_cad_output_digest"] = expected_cad_output_digest
        if expected_cad_observable_id.strip():
            kwargs["expected_cad_observable_id"] = expected_cad_observable_id
        if expected_cad_observable_family.strip():
            kwargs["expected_cad_observable_family"] = expected_cad_observable_family
        if expected_length_unit.strip():
            kwargs["expected_length_unit"] = expected_length_unit
        if expected_area_unit.strip():
            kwargs["expected_area_unit"] = expected_area_unit
        if expected_volume_unit.strip():
            kwargs["expected_volume_unit"] = expected_volume_unit
        if expected_measurement_convention.strip():
            kwargs["expected_measurement_convention"] = expected_measurement_convention
        if expected_measurement_postprocess_row_convention_schema_id.strip():
            kwargs["expected_measurement_postprocess_row_convention_schema_id"] = (
                expected_measurement_postprocess_row_convention_schema_id
            )
        if expected_measurement_component_basis_schema_id.strip():
            kwargs["expected_measurement_component_basis_schema_id"] = (
                expected_measurement_component_basis_schema_id
            )
        summary = shape_cad_handoff_manifest_gate(shape_rows, **kwargs)
    except Exception as e:  # noqa: BLE001
        return json.dumps({
            "status": "error",
            "stage": "cad_handoff_manifest",
            "error": f"{type(e).__name__}: {e}",
        }, indent=2)

    return json.dumps(summary, indent=2)


@mcp.tool()
async def execute_build123d(
    script: str,
    export_dir: str = "",
    export_format: str = "step",
) -> str:
    """
    Execute a build123d Python script and return geometry information.

    The script should create build123d objects. The last Part, Sketch, or
    Compound assigned to a variable will be inspected and optionally exported.

    Args:
        script: Python code using build123d. Must be self-contained.
        export_dir: Directory for exported files. Empty string = no export.
        export_format: Export format: "step", "brep", or "stl".

    Returns:
        JSON with geometry info: validity, volume, area, face/edge counts,
        min edge length, bounding box, labels, and export path (if requested).
    """
    # Run the blocking work in a worker thread so we don't stall the MCP
    # stdio event loop. mcp 1.27.x FastMCP still calls sync @mcp.tool()
    # functions directly inside the asyncio loop (see
    # https://github.com/modelcontextprotocol/python-sdk/issues/1839 ),
    # which makes any blocking I/O — including cold `from build123d import *` —
    # freeze the whole server for every request.
    import anyio
    return await anyio.to_thread.run_sync(
        lambda: _execute_build123d_sync(script, export_dir, export_format)
    )


def _execute_build123d_sync(
    script: str,
    export_dir: str = "",
    export_format: str = "step",
) -> str:
    # Capture stdout
    old_stdout = sys.stdout
    sys.stdout = captured = StringIO()

    # Pre-flight: warn if a very similar script recently failed
    preflight = _b3d_preflight_warn(script)

    namespace = {}
    try:
        exec(  # noqa: S102
            "from build123d import *\n" + script,
            namespace,
        )
    except Exception as e:
        sys.stdout = old_stdout
        tb = traceback.format_exc()
        _fl.record_failure("build123d", {
            "input": script,
            "error": f"{type(e).__name__}: {e}",
            "traceback": tb[-1200:],
            "stdout": captured.getvalue()[-400:],
        })
        err_payload: dict = {
            "status": "error",
            "error": tb,
            "stdout": captured.getvalue(),
            "hint": _build123d_failure_hint(tb),
        }
        if preflight:
            err_payload["preflight"] = preflight
        return json.dumps(err_payload, indent=2)
    finally:
        sys.stdout = old_stdout

    stdout_text = captured.getvalue()

    # Find the last Shape-like object in namespace
    from build123d import Shape, Compound, Part, Sketch, Curve

    target = None
    target_name = None
    for name, obj in namespace.items():
        if name.startswith("_"):
            continue
        if isinstance(obj, Shape):
            target = obj
            target_name = name

    if target is None:
        return json.dumps({
            "status": "ok",
            "message": "Script executed but no Shape object found in namespace.",
            "stdout": stdout_text,
        }, indent=2)

    # Inspect geometry
    info = {
        "status": "ok",
        "variable": target_name,
        "type": type(target).__name__,
        "is_valid": target.is_valid,
        "label": target.label or "(none)",
    }

    try:
        info["volume"] = round(target.volume, 8)
    except Exception:
        info["volume"] = None

    try:
        info["area"] = round(target.area, 8)
    except Exception:
        info["area"] = None

    edges = target.edges()
    faces = target.faces()
    info["face_count"] = len(faces)
    info["edge_count"] = len(edges)
    info["vertex_count"] = len(target.vertices())

    if edges:
        edge_lengths = [e.length for e in edges]
        info["min_edge_length"] = round(min(edge_lengths), 8)
        info["max_edge_length"] = round(max(edge_lengths), 8)

    try:
        bb = target.bounding_box()
        info["bounding_box"] = {
            "min": [round(bb.min.X, 6), round(bb.min.Y, 6), round(bb.min.Z, 6)],
            "max": [round(bb.max.X, 6), round(bb.max.Y, 6), round(bb.max.Z, 6)],
            "size": [round(s, 6) for s in bb.size],
        }
    except Exception:
        pass

    # CAE quality warnings
    warnings = []
    if edges:
        min_e = info.get("min_edge_length", 0)
        bb_size = info.get("bounding_box", {}).get("size", [1, 1, 1])
        char_len = max(bb_size) if bb_size else 1
        if char_len > 0 and min_e / char_len < 0.001:
            warnings.append(
                f"Micro-edge detected: {min_e:.2e} "
                f"(ratio {min_e/char_len:.2e} of characteristic length {char_len:.4f})"
            )
    if not target.is_valid:
        warnings.append("Shape is not valid — may cause meshing issues")
    if warnings:
        info["cae_warnings"] = warnings

    # Export if requested
    if export_dir:
        from build123d import export_step, export_brep, export_stl

        export_path = Path(export_dir)
        export_path.mkdir(parents=True, exist_ok=True)
        base_name = target.label if target.label else target_name
        fmt = export_format.lower().strip()

        if fmt == "step":
            fpath = export_path / f"{base_name}.step"
            export_step(target, str(fpath))
            info["exported"] = str(fpath)
        elif fmt == "brep":
            fpath = export_path / f"{base_name}.brep"
            export_brep(target, str(fpath))
            info["exported"] = str(fpath)
        elif fmt == "stl":
            fpath = export_path / f"{base_name}.stl"
            export_stl(target, str(fpath))
            info["exported"] = str(fpath)
        else:
            info["export_error"] = f"Unknown format: {fmt}"

    if stdout_text:
        info["stdout"] = stdout_text

    if preflight:
        info["preflight"] = preflight

    return json.dumps(info, indent=2)


@mcp.tool()
def preview_shape_in_cubit(script: str, label: str = "preview") -> str:
    """
    Run a build123d script and show the resulting Shape in the
    **persistent Cubit Viewer** (GUI window), replacing OCP CAD Viewer.

    Lab direction (2026-04-19): Cubit Viewer is the single CAD preview
    target — industrial-grade visualization / measurement. OCP CAD
    Viewer has been retired; Plan A's persistent Cubit GUI handles
    both AI-driven and user-interactive workflows in one window.

    Pipeline:
      1. exec the build123d script (like execute_build123d).
      2. Extract the last Shape found in the namespace.
      3. Export to a temp STEP file.
      4. Hand to the persistent Cubit session (`cubit_show`-equivalent);
         the daemon imports STEP into its live Cubit GUI window.
      5. Temp STEP is deleted by default (Cubit now owns the geometry).

    Performance:
      - First call: ~0.7-1.0 s (Cubit daemon cold-start) + STEP roundtrip
      - Subsequent calls: <1 s (daemon warm + STEP import ~0.5 s)

    Args:
        script: self-contained build123d Python code.
        label: label attached to the shape (used for STEP filename).

    Returns:
        JSON with shape stats (volume / bbox / faces / edges) plus the
        Cubit session response (per-command ok flags).
    """
    import sys as _sys
    import tempfile as _tf
    import traceback as _tb
    from io import StringIO
    from pathlib import Path as _P

    namespace = {}
    _buf = StringIO()
    _orig_stdout = sys.stdout
    sys.stdout = _buf
    try:
        exec(  # noqa: S102
            "from build123d import *\n" + script,
            namespace,
        )
    except Exception:
        sys.stdout = _orig_stdout
        return json.dumps({
            "status": "error",
            "stage": "script_exec",
            "error": _tb.format_exc(),
            "stdout": _buf.getvalue(),
        }, indent=2)
    finally:
        sys.stdout = _orig_stdout
    stdout_text = _buf.getvalue()

    from build123d import Shape, export_step

    target = None
    target_name = None
    for n, obj in namespace.items():
        if n.startswith("_"):
            continue
        if isinstance(obj, Shape):
            target, target_name = obj, n
    if target is None:
        return json.dumps({
            "status": "error",
            "stage": "extract",
            "error": "Script executed but no build123d Shape was found.",
            "stdout": stdout_text,
        }, indent=2)

    # Write to a temp STEP; Cubit daemon imports, then we clean up.
    tmp = _P(_tf.mkdtemp(prefix="b123d_to_cubit_"))
    step_path = tmp / f"{label}.step"
    try:
        export_step(target, str(step_path))
    except Exception:
        return json.dumps({
            "status": "error",
            "stage": "step_export",
            "error": _tb.format_exc(),
        }, indent=2)

    # Hand to Cubit viewer (lazy-start daemon).
    try:
        from radia_mcp.cubit import session as _cs
    except ImportError as e:
        return json.dumps({
            "status": "error",
            "stage": "cubit_import",
            "error": f"radia_mcp.cubit not available: {e}",
        }, indent=2)

    try:
        sess = _cs.CubitSession.get()
        sess.ensure_started()
    except _cs.CubitSessionError as e:
        return json.dumps({
            "status": "error",
            "stage": "cubit_session",
            "error": str(e),
            "hint": ("Is Coreform Cubit installed? Override path via "
                     "`CUBIT_BIN_DIR` env var or `set_cubit_bin_dir()`."),
        }, indent=2)

    step_fwd = str(step_path).replace("\\", "/")
    try:
        r = sess.call("cmd", [f'import step "{step_fwd}"'])
    except _cs.CubitSessionError as e:
        return json.dumps({
            "status": "error",
            "stage": "cubit_rpc",
            "error": str(e),
        }, indent=2)

    info = {
        "status": "ok",
        "stage": "shown_in_cubit",
        "viewer": "cubit_persistent_session",
        "label": label,
        "variable": target_name,
        "type": type(target).__name__,
        "is_valid": target.is_valid,
        "cubit_result": r,
    }
    try:
        info["volume"] = round(target.volume, 6)
    except Exception:
        info["volume"] = None
    try:
        info["face_count"] = len(target.faces())
        info["edge_count"] = len(target.edges())
    except Exception:
        pass
    try:
        bb = target.bounding_box()
        info["bounding_box"] = {
            "min": [round(bb.min.X, 4), round(bb.min.Y, 4), round(bb.min.Z, 4)],
            "max": [round(bb.max.X, 4), round(bb.max.Y, 4), round(bb.max.Z, 4)],
        }
    except Exception:
        pass
    if stdout_text:
        info["stdout"] = stdout_text

    # Temp STEP can be cleaned now — Cubit has imported it.
    try:
        step_path.unlink(missing_ok=True)
        tmp.rmdir()
    except OSError:
        pass

    return json.dumps(info, indent=2)


@mcp.tool()
def inspect_geometry(file_path: str) -> str:
    """
    Inspect a STEP or BREP file for CAE quality.

    Loads the geometry and reports: validity, volume, area, face/edge/vertex
    counts, minimum edge length, bounding box, and CAE quality warnings.

    Args:
        file_path: Path to .step or .brep file.

    Returns:
        JSON with geometry inspection results and CAE quality warnings.
    """
    p = Path(file_path)
    if not p.exists():
        return json.dumps({"status": "error", "error": f"File not found: {file_path}"})

    try:
        suffix = p.suffix.lower()
        if suffix in (".step", ".stp"):
            from build123d import import_step
            shape = import_step(str(p))
        elif suffix in (".brep",):
            from build123d import import_brep
            shape = import_brep(str(p))
        else:
            return json.dumps({
                "status": "error",
                "error": f"Unsupported format: {suffix}. Use .step or .brep",
            })
    except Exception:
        return json.dumps({
            "status": "error",
            "error": traceback.format_exc(),
        })

    info = {
        "status": "ok",
        "file": str(p),
        "type": type(shape).__name__,
        "is_valid": shape.is_valid,
        "label": shape.label or "(none)",
    }

    try:
        info["volume"] = round(shape.volume, 8)
    except Exception:
        info["volume"] = None

    try:
        info["area"] = round(shape.area, 8)
    except Exception:
        info["area"] = None

    edges = shape.edges()
    faces = shape.faces()
    info["face_count"] = len(faces)
    info["edge_count"] = len(edges)
    info["vertex_count"] = len(shape.vertices())

    if edges:
        edge_lengths = sorted([e.length for e in edges])
        info["min_edge_length"] = round(edge_lengths[0], 8)
        info["max_edge_length"] = round(edge_lengths[-1], 8)
        info["edge_length_histogram"] = {
            "p5": round(edge_lengths[max(0, len(edge_lengths) // 20)], 8),
            "p25": round(edge_lengths[len(edge_lengths) // 4], 8),
            "median": round(edge_lengths[len(edge_lengths) // 2], 8),
            "p75": round(edge_lengths[3 * len(edge_lengths) // 4], 8),
        }

    try:
        bb = shape.bounding_box()
        info["bounding_box"] = {
            "min": [round(bb.min.X, 6), round(bb.min.Y, 6), round(bb.min.Z, 6)],
            "max": [round(bb.max.X, 6), round(bb.max.Y, 6), round(bb.max.Z, 6)],
            "size": [round(s, 6) for s in bb.size],
        }
    except Exception:
        pass

    # CAE quality analysis
    warnings = []
    if edges:
        min_e = info["min_edge_length"]
        bb_size = info.get("bounding_box", {}).get("size", [1, 1, 1])
        char_len = max(bb_size) if bb_size else 1
        if char_len > 0 and min_e / char_len < 0.001:
            warnings.append(
                f"Micro-edge: {min_e:.2e} "
                f"(ratio {min_e/char_len:.2e})"
            )
        if min_e / char_len < 0.01:
            warnings.append(
                f"Short edge: {min_e:.2e} may require local mesh refinement"
            )
    if not shape.is_valid:
        warnings.append("Shape is not valid")

    # Per-face info (labels, area, center, normal)
    if faces:
        face_info = []
        for idx, face in enumerate(faces):
            fi = {"index": idx, "area": round(face.area, 6)}
            if face.label:
                fi["label"] = face.label
            try:
                c = face.center()
                fi["center"] = [round(c.X, 6), round(c.Y, 6), round(c.Z, 6)]
            except Exception:
                pass
            try:
                n = face.normal_at(face.center())
                fi["normal"] = [round(n.X, 6), round(n.Y, 6), round(n.Z, 6)]
            except Exception:
                pass
            face_info.append(fi)
        info["faces"] = face_info

    # Check for children/assembly structure
    try:
        children = shape.children
        if children:
            info["children"] = []
            for child in children:
                child_info = {
                    "label": child.label or "(none)",
                    "type": type(child).__name__,
                }
                try:
                    child_info["volume"] = round(child.volume, 8)
                except Exception:
                    pass
                try:
                    child_info["face_count"] = len(child.faces())
                except Exception:
                    pass
                info["children"].append(child_info)
    except Exception:
        pass

    if warnings:
        info["cae_warnings"] = warnings

    return json.dumps(info, indent=2)


@mcp.tool()
def section_along_path(
    file_path: str,
    path_json: str,
    n_sections: int = 0,
) -> str:
    """
    Section a STEP/BREP coil solid along a discrete path and extract
    per-segment cross-section dimensions.

    Given a solid and a centerline path (list of [x, y, z] points in
    the same units as the CAD file, typically mm), the tool:
    1. Computes segment midpoints and tangent vectors
    2. Sections the solid perpendicular to the tangent at each midpoint
    3. Picks the face closest to the path center (multi-turn safe)
    4. Returns area, estimated width/height per segment

    This is the CAD-side half of the PEEC filament extraction pipeline.
    Feed the output into PEECBuilder.add_connected_segment() with
    the per-segment (w, h) values (convert to meters if CAD is in mm).

    Args:
        file_path: Path to .step or .brep file.
        path_json: JSON string encoding a list of [x, y, z] path points
            (N+1 points for N segments), in CAD units (typically mm).
            Example: "[[50,0,0],[0,50,4],[-50,0,8],[0,-50,12],[50,0,16]]"
        n_sections: If > 0, resample path to this many equidistant
            segments before sectioning.  0 = use path points as-is.

    Returns:
        JSON with per-segment: index, midpoint, area, w_est, h_est,
        plus a summary of the taper range.
    """
    import math
    import numpy as np

    p = Path(file_path)
    if not p.exists():
        return json.dumps({"status": "error", "error": f"File not found: {file_path}"})

    try:
        pts = json.loads(path_json)
    except Exception:
        return json.dumps({"status": "error", "error": "path_json is not valid JSON"})

    pts = np.array(pts, dtype=np.float64)
    if pts.ndim != 2 or pts.shape[1] != 3 or pts.shape[0] < 2:
        return json.dumps({
            "status": "error",
            "error": f"path must be (N, 3) with N >= 2, got shape {pts.shape}",
        })

    # Optional resample
    if n_sections > 0 and n_sections != pts.shape[0] - 1:
        from scipy.interpolate import interp1d
        cum = np.zeros(len(pts))
        for i in range(1, len(pts)):
            cum[i] = cum[i - 1] + np.linalg.norm(pts[i] - pts[i - 1])
        t_orig = cum / cum[-1]
        t_new = np.linspace(0, 1, n_sections + 1)
        pts = np.column_stack([
            interp1d(t_orig, pts[:, k], kind="cubic")(t_new)
            for k in range(3)
        ])

    n_seg = pts.shape[0] - 1
    midpoints = 0.5 * (pts[:-1] + pts[1:])
    diffs = np.diff(pts, axis=0)
    norms = np.linalg.norm(diffs, axis=1, keepdims=True)
    norms = np.maximum(norms, 1e-30)
    tangents = diffs / norms

    # Load solid
    suffix = p.suffix.lower()
    try:
        if suffix in (".step", ".stp"):
            from build123d import import_step
            solid = import_step(str(p))
        else:
            from build123d import import_brep
            solid = import_brep(str(p))
    except Exception:
        return json.dumps({"status": "error", "error": traceback.format_exc()})

    from build123d import section, Plane, Vector

    segments = []
    for i in range(n_seg):
        c = midpoints[i]
        t = tangents[i]
        origin = Vector(float(c[0]), float(c[1]), float(c[2]))
        z_dir = Vector(float(t[0]), float(t[1]), float(t[2]))
        sec_plane = Plane(origin=origin, z_dir=z_dir)
        try:
            cross = section(solid, section_by=sec_plane)
        except Exception:
            segments.append({"index": i, "error": "section failed"})
            continue

        if not cross or len(cross.faces()) == 0:
            segments.append({"index": i, "error": "no faces"})
            continue

        best_face = min(cross.faces(),
                        key=lambda f: (f.center() - origin).length)
        area = best_face.area
        side = math.sqrt(area)
        seg_len = float(norms[i, 0])

        segments.append({
            "index": i,
            "midpoint": [round(c[0], 6), round(c[1], 6), round(c[2], 6)],
            "length": round(seg_len, 6),
            "area": round(area, 6),
            "w_est": round(side, 6),
            "h_est": round(side, 6),
        })

    areas = [s["area"] for s in segments if "area" in s]
    ws = [s["w_est"] for s in segments if "w_est" in s]
    summary = {
        "n_segments": n_seg,
        "n_ok": len(areas),
        "area_range": [round(min(areas), 4), round(max(areas), 4)] if areas else None,
        "w_range": [round(min(ws), 4), round(max(ws), 4)] if ws else None,
    }

    return json.dumps({"status": "ok", "summary": summary, "segments": segments}, indent=2)


@mcp.tool()
def generate_helix_coil(
    radius: float = 50.0,
    pitch: float = 10.0,
    n_turns: float = 5.0,
    w_start: float = 4.0,
    h_start: float = 4.0,
    w_end: float = 4.0,
    h_end: float = 4.0,
    sections_per_turn: int = 12,
    export_dir: str = "",
    label: str = "helix_coil",
) -> str:
    """
    Generate a helical coil with optional cross-section taper.

    Creates a solid coil by lofting rectangular cross-sections along a
    helix path.  Cross-section dimensions interpolate linearly from
    (w_start, h_start) at the bottom to (w_end, h_end) at the top.

    All dimensions are in the CAD working units (typically mm).

    For PEEC filament extraction, use section_along_path() on the
    exported file with the helix centerline path.

    Args:
        radius: Helix radius (center of wire to axis).
        pitch: Axial advance per turn.
        n_turns: Number of turns (can be fractional).
        w_start: Cross-section width at bottom.
        h_start: Cross-section height at bottom.
        w_end: Cross-section width at top.
        h_end: Cross-section height at top.
        sections_per_turn: Loft sections per turn (12 = 30 deg each).
        export_dir: Directory for STEP export (empty = no export).
        label: Name for the solid and export filename.

    Returns:
        JSON with geometry info (volume, area, bbox) and the helix
        path as a list of [x, y, z] points (for section_along_path).
    """
    import math
    from build123d import (BuildSketch, Rectangle, Plane, Vector, loft,
                           export_step)

    n_sec = int(n_turns * sections_per_turn) + 1
    height = pitch * n_turns
    sketches = []
    path_pts = []

    for i in range(n_sec):
        t = i / (n_sec - 1) if n_sec > 1 else 0
        angle = 2 * math.pi * n_turns * t
        z = height * t
        cx = radius * math.cos(angle)
        cy = radius * math.sin(angle)
        w = w_start + (w_end - w_start) * t
        h = h_start + (h_end - h_start) * t

        dx = -radius * math.sin(angle) * 2 * math.pi * n_turns
        dy = radius * math.cos(angle) * 2 * math.pi * n_turns
        dz = height
        norm = math.sqrt(dx**2 + dy**2 + dz**2)
        tangent = (dx / norm, dy / norm, dz / norm)

        # Cross-section perpendicular to helix tangent
        x_dir_raw = (tangent[1], -tangent[0], 0)
        x_norm = math.sqrt(x_dir_raw[0]**2 + x_dir_raw[1]**2)
        if x_norm < 1e-10:
            x_dir_raw = (1, 0, 0)
            x_norm = 1.0
        x_dir = (x_dir_raw[0] / x_norm, x_dir_raw[1] / x_norm, 0)

        plane = Plane(origin=(cx, cy, z), x_dir=x_dir, z_dir=tangent)
        with BuildSketch(plane) as sk:
            Rectangle(w, h)
        sketches.append(sk.sketch)
        path_pts.append([round(cx, 6), round(cy, 6), round(z, 6)])

    coil = loft(sketches)
    coil.label = label

    info = {
        "status": "ok",
        "label": label,
        "n_turns": n_turns,
        "radius": radius,
        "pitch": pitch,
        "cross_section_start": [w_start, h_start],
        "cross_section_end": [w_end, h_end],
        "is_valid": coil.is_valid,
        "volume": round(coil.volume, 4),
        "face_count": len(coil.faces()),
        "edge_count": len(coil.edges()),
    }

    try:
        bb = coil.bounding_box()
        info["bounding_box"] = {
            "min": [round(bb.min.X, 4), round(bb.min.Y, 4), round(bb.min.Z, 4)],
            "max": [round(bb.max.X, 4), round(bb.max.Y, 4), round(bb.max.Z, 4)],
        }
    except Exception:
        pass

    info["path_points"] = path_pts

    if export_dir:
        export_path = Path(export_dir)
        export_path.mkdir(parents=True, exist_ok=True)
        fpath = export_path / f"{label}.step"
        export_step(coil, str(fpath))
        info["exported"] = str(fpath)

    return json.dumps(info, indent=2)


# ============================================================
# MCP Prompts
# ============================================================

@mcp.prompt()
def new_cae_geometry(geometry_type: str = "ih_coil") -> str:
    """Create a new CAE geometry with build123d."""
    base = (
        "Create a CAE-ready geometry using build123d.\n\n"
        "FIRST CHECK THE TESTED LIBRARY: build123d_usage('parametric_library') -- ready-made,\n"
        "Netgen-meshable, region-labelled generators in radia_mcp.build123d.modeling (annular_segment,\n"
        "tube, racetrack_coil, polar_array, linear_array, mirrored, assembly) and .archetypes\n"
        "(cylindrical_magnet/block_magnet, halbach_ring, pole_tip, c_core, multipole_yoke, h_dipole,\n"
        "solenoid, helmholtz_pair, cos_theta_dipole, magnetization_map).  Prefer these over hand-rolled\n"
        "geometry; fall back to primitives + boolean only for shapes they don't cover.\n\n"
        "Use build123d_usage tool for API reference.\n\n"
        "Workflow:\n"
        "1. Compose from the parametric_library, or primitives + boolean (build123d_usage('primitives_3d'))\n"
        "2. Assign labels for material regions on the Compound's CHILDREN (build123d_usage('topology'))\n"
        "3. Validate: check is_valid, min edge length (build123d_usage('cae_guidelines'))\n"
        "4. Export BREP (for Netgen) or STEP (for Cubit)\n"
        "5. Use execute_build123d tool to run and validate\n\n"
    )

    geo = geometry_type.strip().lower()

    if geo in ("ih", "ih_coil", "induction_heating"):
        base += (
            "Induction Heating model:\n"
            "- Cylindrical workpiece (steel) + surrounding coil (copper)\n"
            "- Air domain enclosing both\n"
            "- Label: 'workpiece', 'coil', 'air'\n"
            "- See build123d_usage('examples') for complete code\n"
        )
    elif geo in ("transformer", "e_core", "ecore"):
        base += (
            "E-core Transformer:\n"
            "- archetypes.e_core(width, height, depth, leg_width, back_thickness) -- E-shaped iron\n"
            "- add winding regions (archetypes.solenoid / modeling.racetrack_coil) around the legs\n"
            "- Label: 'iron_core', 'winding_primary', 'winding_secondary'\n"
        )
    elif geo in ("motor", "pmsm", "stator", "rotor"):
        base += (
            "PM synchronous machine (PMSM) cross-section:\n"
            "- archetypes.slotted_stator(r_bore, r_yoke, n_slots, slot_depth, slot_span_deg, h)\n"
            "- archetypes.spm_rotor(r_shaft, r_rotor, n_poles, magnet_thickness, magnet_span_deg, h)\n"
            "  (radial alternating-N/S magnets; easy axes in the labels -> magnetization_map)\n"
            "- drops into the AGE solver: radia_ngsolve.airgap_motor_workflow (incl. hysteresis sweep)\n"
        )
    elif geo in ("dipole", "magnet", "accelerator"):
        base += (
            "Accelerator Dipole Magnet:\n"
            "- archetypes.h_dipole(width, height, depth, leg, pole_width, gap) -- H-frame yoke, OR\n"
            "  archetypes.c_core(...) for a C-yoke, OR archetypes.multipole_yoke(n_poles=2, ...)\n"
            "- add a coil (archetypes.solenoid / cos_theta_dipole) and a 'bore' region\n"
            "- Label: 'yoke', 'coil', 'bore', 'air'\n"
        )
    elif geo in ("halbach", "halbach_ring", "pm_array", "undulator"):
        base += (
            "Halbach permanent-magnet ring:\n"
            "- archetypes.halbach_ring(r_in, r_out, h, n_segments, pole_pairs=1) -- dipole (uniform\n"
            "  bore field), pole_pairs=2 quadrupole, etc.\n"
            "- archetypes.magnetization_map(ring, Br) -> {label: (Mx,My)} drives the solver per segment\n"
            "- each segment label encodes the Mallinson easy axis; see build123d_usage('parametric_library')\n"
        )
    elif geo in ("solenoid", "coil", "helmholtz", "litz"):
        base += (
            "Coil / winding region:\n"
            "- archetypes.solenoid(r_in, r_out, h) bundle, archetypes.helmholtz_pair(..., separation),\n"
            "  modeling.racetrack_coil(...), or archetypes.cos_theta_dipole(...) for a transverse field\n"
            "- archetypes.litz_wire(n_strands, strand_radius, bundle_radius, length, pitch) -- twisted\n"
            "  strands as SEPARATE conductor regions for AC-loss (skin/proximity) analysis (PEEC / FE)\n"
            "- Label as 'coil'; the solver sets the current density\n"
        )
    else:
        base += (
            f"Custom geometry type: {geometry_type}\n"
            "- Start with build123d_usage('overview') for CAE-safe subset\n"
            "- Use primitives + boolean for clean geometry\n"
            "- Validate with execute_build123d before meshing\n"
        )

    return base


# ============================================================
# Convenience tools — common shapes routed straight to Cubit
# ============================================================
#
# 2026-04-19: `preview_shape_in_cubit` already handles the build123d →
# STEP → Cubit pipeline, but writing the build123d script every time
# for simple shapes (text, extrude, boolean) is noisy for the LLM /
# user.  These helpers produce the script internally and delegate.

def _hand_to_cubit(part_script: str, label: str) -> str:
    """Internal: run a build123d snippet, export STEP, send to Cubit.

    Reuses the full `preview_shape_in_cubit` path. Returns its JSON
    verbatim so callers can forward directly.
    """
    return preview_shape_in_cubit(script=part_script, label=label)


@mcp.tool()
def preview_text(text: str, font_size: float = 30.0,
                 thickness: float = 5.0,
                 font: str = "") -> str:
    """Generate 3D extruded text and send it to the live Cubit viewer.

    Builds a build123d `Text` sketch at `font_size`, extrudes by
    `thickness`, exports STEP, and delegates to Cubit. Each letter
    becomes its own prismatic volume — ideal for demonstrating hex
    meshing via `scheme auto` / `scheme sweep`.

    Args:
        text: the string to render. ASCII works reliably; Unicode depends
            on the font.
        font_size: text height (build123d units, typically mm).
        thickness: extrusion depth along Z.
        font: optional font family name (e.g., "Arial", "Times New Roman").
            Empty uses build123d's platform default.

    Returns the same JSON shape as `preview_shape_in_cubit`.
    """
    safe = text.replace("\\", "\\\\").replace('"', '\\"')
    font_line = f', font="{font}"' if font else ""
    # Expose the resulting Shape as a top-level name so preview_shape_in_cubit's
    # namespace scan ("isinstance(obj, Shape)") picks it up.
    script = (
        "with BuildPart() as _bp:\n"
        "    with BuildSketch():\n"
        f'        Text("{safe}", font_size={float(font_size)}{font_line}, '
        "align=(Align.CENTER, Align.CENTER))\n"
        f"    extrude(amount={float(thickness)})\n"
        "part = _bp.part\n"
    )
    return _hand_to_cubit(script, label=f"text_{safe[:16].replace(' ', '_')}")


@mcp.tool()
def preview_extrude(sketch_code: str, amount: float = 10.0,
                    label: str = "extrude") -> str:
    """Extrude a build123d sketch by `amount` and show in Cubit.

    You provide just the sketch body (the code that goes inside a
    `with BuildSketch():` block), this tool wraps it with BuildPart /
    BuildSketch / extrude and sends the result through.

    Args:
        sketch_code: lines placed inside `with BuildSketch():`. E.g.,
            "Circle(radius=5); Rectangle(10,5)".
        amount: extrusion depth (build123d units).
        label: name attached to the shape.

    Example:
        preview_extrude("Rectangle(20, 10); fillet=Rectangle(2, 2)",
                        amount=5.0)
    """
    # Indent user body inside the with-block.
    body_lines = sketch_code.splitlines() or [sketch_code]
    indented = "\n".join("        " + line.strip() if line.strip() else ""
                         for line in body_lines)
    script = (
        "with BuildPart() as _bp:\n"
        "    with BuildSketch():\n"
        f"{indented}\n"
        f"    extrude(amount={float(amount)})\n"
        "part = _bp.part\n"
    )
    return _hand_to_cubit(script, label=label)


@mcp.tool()
def preview_boolean(op: str, a_script: str, b_script: str,
                    label: str = "boolean") -> str:
    """Perform a boolean (union/subtract/intersect) of two build123d
    snippets and send the result to Cubit.

    Args:
        op: "union" / "subtract" / "intersect".
        a_script: build123d code producing the first part (must assign
            to `a`).  E.g., "a = Box(10, 10, 10)".
        b_script: same for second part, assigning to `b`.
        label: name attached to the shape.

    Example:
        preview_boolean(
            "subtract",
            "a = Box(20, 20, 20)",
            "b = Cylinder(radius=5, height=25)",
        )
    """
    op_map = {
        "union": "+",
        "subtract": "-",
        "minus": "-",
        "intersect": "&",
        "intersection": "&",
    }
    sym = op_map.get(op.strip().lower())
    if sym is None:
        return json.dumps({
            "status": "error",
            "stage": "op",
            "error": f"unknown op {op!r}; use union / subtract / intersect",
        })
    script = (
        f"{a_script.strip()}\n"
        f"{b_script.strip()}\n"
        f"part = a {sym} b\n"
    )
    return _hand_to_cubit(script, label=label)


# ============================================================
# Retrieval + web docs + failure log
# ============================================================

import math as _math
import re as _re

_STOP_WORDS = frozenset({
    "the", "a", "an", "of", "to", "in", "and", "or", "is", "are", "be",
    "on", "at", "for", "with", "as", "by", "from", "that", "this",
    "it", "its", "if", "so", "how", "what", "which", "when", "where",
})
_WORD_RE = _re.compile(r"[a-zA-Z_][a-zA-Z0-9_]*")


def _tokenize(s: str) -> list[str]:
    return [m.group(0).lower() for m in _WORD_RE.finditer(s)]


def _chunk_by_heading(text: str, target_lines: int = 40) -> list[tuple[int, str]]:
    """Same heading-anchored chunker as the cubit server."""
    lines = text.splitlines()
    chunks: list[tuple[int, str]] = []
    cur_start = 0
    cur: list[str] = []
    for i, line in enumerate(lines):
        is_heading = (line.startswith("# ") or line.startswith("## ")
                      or line.startswith("### "))
        if is_heading and cur and (i - cur_start) >= max(5, target_lines // 4):
            chunks.append((cur_start + 1, "\n".join(cur)))
            cur = []
            cur_start = i
        cur.append(line)
        if len(cur) >= target_lines * 2:
            chunks.append((cur_start + 1, "\n".join(cur)))
            cur = []
            cur_start = i + 1
    if cur:
        chunks.append((cur_start + 1, "\n".join(cur)))
    return chunks


def _heading_of(chunk: str) -> str:
    for line in chunk.splitlines():
        s = line.strip()
        if s:
            return s
    return ""


def _b3d_knowledge_topics() -> list[str]:
    """Dynamically enumerate every bundled topic so new ones (added to
    `build123d_knowledge._TOPICS`) show up in `build123d_lookup`
    without a code edit here.
    """
    from .build123d_knowledge import _TOPICS as _T
    return list(_T.keys())


_B3D_KNOWLEDGE_TOPICS = _b3d_knowledge_topics()


def _b3d_build_corpus() -> list[dict]:
    if getattr(_b3d_build_corpus, "_cache", None) is not None:
        return _b3d_build_corpus._cache  # type: ignore[attr-defined]
    corpus: list[dict] = []
    for topic in _B3D_KNOWLEDGE_TOPICS:
        try:
            text = get_build123d_documentation(topic)
        except Exception:
            continue
        if not text:
            continue
        for start, chunk in _chunk_by_heading(text):
            tokens = _tokenize(chunk)
            if not tokens:
                continue
            corpus.append({
                "source": f"kb:{topic}",
                "start_line": start,
                "chunk": chunk,
                "heading": _heading_of(chunk),
                "tokens": tokens,
                "token_set": set(tokens),
            })
    # Failure log as extra source
    failures_text = _fl.failures_as_text("build123d", limit=50)
    if failures_text:
        for start, chunk in _chunk_by_heading(failures_text):
            tokens = _tokenize(chunk)
            if not tokens:
                continue
            corpus.append({
                "source": "failures",
                "start_line": start,
                "chunk": chunk,
                "heading": _heading_of(chunk),
                "tokens": tokens,
                "token_set": set(tokens),
            })
    _b3d_build_corpus._cache = corpus  # type: ignore[attr-defined]
    return corpus


def _score_tfidf(query_terms: list[str], doc: dict, idf: dict[str, float]) -> float:
    tokens = doc["tokens"]
    if not tokens:
        return 0.0
    tf: dict[str, int] = {}
    for tok in tokens:
        tf[tok] = tf.get(tok, 0) + 1
    inv_len = 1.0 / (1.0 + _math.log(1 + len(tokens)))
    score = 0.0
    hits = 0
    for term in query_terms:
        f = tf.get(term, 0)
        if f == 0:
            continue
        hits += 1
        score += (1.0 + _math.log(f)) * idf.get(term, 1.0) * inv_len
    if hits == 0:
        return 0.0
    head_lc = doc["heading"].lower()
    head_hits = sum(1 for t in query_terms if t in head_lc)
    if head_hits:
        score *= 1.0 + 0.75 * head_hits
    score *= 1.0 + 0.25 * (hits / max(1, len(query_terms)))
    return score


@mcp.tool()
def build123d_lookup(query: str, max_results: int = 5) -> str:
    """Retrieve relevant sections from the bundled build123d knowledge +
    failure log, scored by tf-idf with heading boost.

    Complements `build123d_usage(topic=...)` (full topic dump) by letting
    the LLM query with free-form keywords (e.g., "extrude until face",
    "loft profiles", "Plane rotation").

    Args:
        query: keywords (AND-scored).
        max_results: cap on returned chunks.
    """
    q = (query or "").strip()
    if not q:
        return json.dumps({"status": "error", "error": "empty query"})
    terms = [t for t in _tokenize(q) if t not in _STOP_WORDS]
    if not terms:
        return json.dumps({"status": "error",
                            "error": "no searchable terms (all stop words?)"})
    corpus = _b3d_build_corpus()
    n = max(1, len(corpus))
    idf = {
        t: _math.log((n + 1) / (sum(1 for d in corpus if t in d["token_set"]) + 1)) + 1.0
        for t in terms
    }
    scored: list[tuple[float, dict]] = []
    for doc in corpus:
        s = _score_tfidf(terms, doc, idf)
        if s > 0:
            scored.append((s, doc))
    scored.sort(key=lambda x: (-x[0], x[1]["source"], x[1]["start_line"]))
    top = scored[: max(1, int(max_results))]
    return json.dumps({
        "status": "ok",
        "query": q,
        "terms": terms,
        "candidates_total": len(scored),
        "returned": len(top),
        "results": [
            {
                "source": d["source"],
                "start_line": d["start_line"],
                "heading": d["heading"][:120],
                "score": round(s, 3),
                "chunk": d["chunk"],
            }
            for s, d in top
        ],
    }, indent=2)


@mcp.tool()
def build123d_recent_failures(limit: int = 10) -> str:
    """Return the last N failed `execute_build123d` invocations (from log)."""
    recs = _fl.recent_failures("build123d", limit=limit)
    return json.dumps({
        "status": "ok",
        "count": len(recs),
        "log_path": str(_fl._log_path("build123d")),
        "entries": recs,
    }, indent=2)


@mcp.tool()
def build123d_web_docs(query: str, max_hits: int = 5,
                       force_refresh: bool = False) -> str:
    """Fetch live build123d documentation (readthedocs) and grep for `query`.

    Complements `build123d_lookup` when a symbol is missing from the
    bundled kb or you need the authoritative API signature. Pages are
    cached for 24 h.

    Args:
        query: keywords (AND-matched, case-insensitive).
        max_hits: cap on matching lines across all indexed pages.
        force_refresh: bypass cache and re-fetch.
    """
    q = (query or "").strip()
    if not q:
        return json.dumps({"status": "error", "error": "empty query"})
    urls = _wd.DOCS_INDEX.get("build123d", [])
    all_hits: list[dict] = []
    errors: list[dict] = []
    remaining = max(1, int(max_hits))
    for url in urls:
        if remaining <= 0:
            break
        data = _wd.fetch(url, force_refresh=force_refresh)
        if data.get("status") != "ok":
            errors.append({"url": url, "error": data.get("error")})
            continue
        hits = _wd.search_text(data["text"], q, max_hits=remaining)
        if hits:
            for h in hits:
                h["url"] = data["url"]
                h["from_cache"] = data.get("from_cache", False)
                all_hits.append(h)
            remaining -= len(hits)
    return json.dumps({
        "status": "ok",
        "query": q,
        "hits": all_hits,
        "errors": errors,
    }, indent=2)


@mcp.tool()
def build123d_examples(query: str, limit: int = 3,
                       force_refresh: bool = False) -> str:
    """Search build123d + **bd_warehouse** + **GitHub Issues** (unioned).

    Three sources, one query:
      1. **`gumyr/build123d/examples/`** (~65 curated Python scripts:
         bicycle tire, Benchy, fillets, lofts, joints, roller coaster,
         vase, …).
      2. **`gumyr/bd_warehouse`** (industrial parts lib: `bearing`,
         `fastener`, `flange`, `gear`, `open_builds`, `pipe`,
         `sprocket`) — fills gaps build123d itself doesn't cover
         (e.g., `gear involute` lives here, not in core examples).
      3. **`gumyr/build123d` Issues** (GitHub issues + comment threads)
         — the de-facto forum: Q&A, bug reports with reproducers,
         feature discussion, multi-turn threads. Scraped via REST
         anonymously (60 req/h — default 60 issues indexed).

    All fetched via GitHub API; 7-day cache. Hits carry `source` so
    you can tell examples apart from issue threads.

    Args:
        query: keywords (e.g., "helix sweep", "gear involute",
            "loft profile", "Plane.__mul__", "None return type").
        limit: cap on returned examples.
        force_refresh: re-scrape all three GitHub sources.
    """
    r = _ex.search_examples("build123d", query, limit=limit,
                             force_refresh=force_refresh)
    return json.dumps(r, indent=2, ensure_ascii=False)


@mcp.tool()
def build123d_examples_refresh(limit: int = 0,
                                max_issues: int = 60,
                                include_youtube: bool = True) -> str:
    """Force-refresh all build123d sources from GitHub + YouTube.

    Args:
        limit: cap for `examples/` dir fetch (0 = all ~65).
        max_issues: cap for issue-thread scrape.
        include_youtube: also scrape YouTube transcripts (build123d
            tutorial videos via `youtube-transcript-api`; needs
            `pip install radia-mcp[youtube]`).
    """
    out = {
        "build123d": _ex.refresh_build123d_examples(limit=limit),
        "bd_warehouse": _ex.refresh_bd_warehouse_examples(),
        "build123d_discussions": _ex.refresh_build123d_discussions(
            max_issues=max_issues,
        ),
    }
    if include_youtube:
        out["build123d_youtube"] = _ex.refresh_build123d_youtube()
    return json.dumps(out, indent=2, ensure_ascii=False)


def _b3d_preflight_warn(script: str) -> dict | None:
    """Return a compact warning if the incoming script resembles a
    recent build123d failure (token Jaccard ≥ 0.6). Non-blocking."""
    try:
        recs = _fl.recent_failures("build123d", limit=30)
    except Exception:
        return None
    if not recs:
        return None
    incoming = set(_tokenize(script))
    if not incoming:
        return None
    best, best_j = None, 0.0
    for r in recs:
        rin = r.get("input")
        toks = set(_tokenize(str(rin) if rin else ""))
        if not toks:
            continue
        j = len(incoming & toks) / max(1, len(incoming | toks))
        if j > best_j:
            best_j, best = j, r
    if best is None or best_j < 0.6:
        return None
    return {
        "similar_past_failure": True,
        "jaccard": round(best_j, 2),
        "past_error": str(best.get("error") or "")[:300],
        "past_iso": best.get("iso"),
        "advice": ("Similar script failed recently. "
                   "Check `build123d_recent_failures(limit=10)` before "
                   "re-running."),
    }


@mcp.tool()
def build123d_ask(query: str, limit: int = 6,
                  include_web: bool = False) -> str:
    """One-shot search across every build123d knowledge surface.

    Unions the result sets of:
      - `build123d_lookup`   (bundled kb + failure log)
      - `build123d_examples` (GitHub examples/ + bd_warehouse +
        issues + GraphQL Discussions — unioned via the `build123d`
        family in radia_mcp.common.examples)
      - optionally `build123d_web_docs` (live readthedocs)

    Hits carry `layer` (`kb` / `examples` / `web`) so the caller
    sees provenance. Same tokenizer / idf across layers so a single
    score-sorted list is meaningful.

    Args:
        query: keywords.
        limit: cap on returned hits across all layers.
        include_web: if True, also hit `build123d_web_docs` (live HTTP).
    """
    q = (query or "").strip()
    if not q:
        return json.dumps({"status": "error", "error": "empty query"})

    results: list[dict] = []

    # Layer 1: bundled KB
    try:
        r = json.loads(build123d_lookup(query=q, max_results=limit))
        for h in r.get("results", []):
            h["layer"] = "kb"
            results.append(h)
    except Exception as e:  # noqa: BLE001
        results.append({"layer": "kb", "error": str(e)})

    # Layer 2: scraped examples (family union)
    try:
        ex_r = _ex.search_examples("build123d", q, limit=limit,
                                    auto_refresh_if_empty=False)
        for h in ex_r.get("results", []):
            h["layer"] = "examples"
            results.append(h)
    except Exception as e:  # noqa: BLE001
        results.append({"layer": "examples", "error": str(e)})

    # Layer 3: optional live web
    if include_web:
        try:
            web = json.loads(build123d_web_docs(query=q, max_hits=3))
            for h in web.get("hits", []):
                h["layer"] = "web"
                results.append(h)
        except Exception as e:  # noqa: BLE001
            results.append({"layer": "web", "error": str(e)})

    def _score_of(h: dict) -> float:
        try:
            return float(h.get("score") or 0)
        except (TypeError, ValueError):
            return 0.0
    results.sort(key=lambda h: -_score_of(h))
    return json.dumps({
        "status": "ok",
        "query": q,
        "total": len(results),
        "results": results[: max(1, int(limit))],
    }, indent=2, ensure_ascii=False)


@mcp.tool()
def build123d_discussions(query: str, limit: int = 5,
                          force_refresh: bool = False) -> str:
    """Search the build123d GitHub Issues archive (de-facto forum).

    Every hit is an issue thread: title + body + every comment,
    linked back to its `github.com/gumyr/build123d/issues/N` URL.
    Use this when you need community context (workarounds, edge-case
    reports, feature discussion) that isn't in the bundled docs or
    the curated examples dir.

    Args:
        query: keywords (e.g., "Plane Location mul", "3.13 vtk pin",
            "extrude until face").
        limit: cap on returned threads.
        force_refresh: re-scrape (anonymous 60 req/h — use sparingly).
    """
    r = _ex.search_examples("build123d_discussions", query, limit=limit,
                             force_refresh=force_refresh)
    return json.dumps(r, indent=2, ensure_ascii=False)


# ============================================================
# CadQuery interop (OCCT sibling — shared STEP boundary)
# ============================================================

def _find_cadquery_shape(namespace: dict):
    """Extract the last cadquery Shape-like object from an exec'd namespace.

    Priority:
      1. a `result` variable (cadquery convention for Workplane.val())
      2. the last `Workplane` assigned
      3. the last `Shape`/`Solid`/`Compound`
    """
    try:
        import cadquery as cq
    except ImportError:
        return None, "cadquery not installed"
    target = None
    target_name = None
    # 1. `result` wins (cadquery-cli convention)
    if "result" in namespace:
        target = namespace["result"]
        target_name = "result"
    else:
        # 2 + 3. scan for last Workplane / Shape
        for name, obj in namespace.items():
            if name.startswith("_"):
                continue
            if isinstance(obj, (cq.Workplane,)):
                target = obj
                target_name = name
            elif hasattr(cq, "Shape") and isinstance(obj, cq.Shape):
                target = obj
                target_name = name
    if target is None:
        return None, "no cadquery Workplane or Shape found in script namespace"
    # Normalize Workplane → underlying Compound / Solid
    if hasattr(target, "val") and hasattr(target, "vals"):
        vals = target.vals()
        if not vals:
            return None, "Workplane is empty"
        if len(vals) == 1:
            return (vals[0], target_name), None
        # Multiple shapes — combine into a Compound
        return (cq.Compound.makeCompound(vals), target_name), None
    return (target, target_name), None


@mcp.tool()
def execute_cadquery(script: str,
                     export_dir: str = "",
                     export_format: str = "step") -> str:
    """Execute a **CadQuery** Python script (sibling OCCT CAD lib).

    Enables interop with the broader CadQuery community and with
    `cadquery-mcp` (rishigundakaram/cadquery-mcp-server): AI can
    author or fetch CadQuery scripts and route them through our
    Radia pipeline without hand-porting to build123d.

    Scripts run with `from cadquery import *` prepended. By convention
    they assign the final geometry to `result` (plain Workplane or
    Shape). Falls back to "last Workplane/Shape in namespace" if
    `result` isn't set.

    Returns the same geometry-info shape as `execute_build123d`
    (validity, volume, area, face/edge/vertex counts, bounding box,
    min edge length, CAE warnings) so downstream tools don't care
    which backend authored the shape.

    Args:
        script: Python code using cadquery. Self-contained.
        export_dir: optional directory for STEP/BREP/STL export
            (empty = no export).
        export_format: "step" / "brep" / "stl".
    """
    try:
        import cadquery as cq
    except ImportError:
        return json.dumps({
            "status": "error",
            "error": ("cadquery not installed. Install with "
                      "`pip install radia-mcp[cadquery]` or "
                      "`pip install cadquery`."),
        }, indent=2)

    old_stdout = sys.stdout
    sys.stdout = captured = StringIO()

    namespace = {}
    try:
        exec(  # noqa: S102
            "from cadquery import *\n" + script,
            namespace,
        )
    except Exception as e:
        sys.stdout = old_stdout
        tb = traceback.format_exc()
        _fl.record_failure("cadquery", {
            "input": script,
            "error": f"{type(e).__name__}: {e}",
            "traceback": tb[-1200:],
            "stdout": captured.getvalue()[-400:],
        })
        return json.dumps({
            "status": "error",
            "error": tb,
            "stdout": captured.getvalue(),
            "hint": "Check imports (cadquery exposes `Workplane`, `cq` module, etc.)",
        }, indent=2)
    finally:
        sys.stdout = old_stdout

    stdout_text = captured.getvalue()

    found, err = _find_cadquery_shape(namespace)
    if found is None:
        return json.dumps({
            "status": "ok" if err is None else "error",
            "message": err or "no shape found",
            "stdout": stdout_text,
        }, indent=2)
    target, target_name = found

    info = {
        "status": "ok",
        "variable": target_name,
        "type": type(target).__name__,
    }
    try:
        info["is_valid"] = bool(target.isValid())
    except Exception:
        info["is_valid"] = None
    try:
        info["volume"] = round(target.Volume(), 8)
    except Exception:
        info["volume"] = None
    try:
        info["area"] = round(target.Area(), 8)
    except Exception:
        info["area"] = None
    try:
        faces = target.Faces()
        edges = target.Edges()
        info["face_count"] = len(faces)
        info["edge_count"] = len(edges)
        info["vertex_count"] = len(target.Vertices())
        if edges:
            lens = [e.Length() for e in edges]
            info["min_edge_length"] = round(min(lens), 8)
            info["max_edge_length"] = round(max(lens), 8)
    except Exception:
        pass
    try:
        bb = target.BoundingBox()
        info["bounding_box"] = {
            "min": [round(bb.xmin, 6), round(bb.ymin, 6), round(bb.zmin, 6)],
            "max": [round(bb.xmax, 6), round(bb.ymax, 6), round(bb.zmax, 6)],
            "size": [round(bb.xlen, 6), round(bb.ylen, 6), round(bb.zlen, 6)],
        }
    except Exception:
        pass

    # CAE quality warnings — same logic as execute_build123d
    warnings = []
    if info.get("min_edge_length"):
        bb_size = info.get("bounding_box", {}).get("size", [1, 1, 1])
        char_len = max(bb_size) if bb_size else 1
        min_e = info["min_edge_length"]
        if char_len > 0 and min_e / char_len < 0.001:
            warnings.append(
                f"Micro-edge detected: {min_e:.2e} "
                f"(ratio {min_e/char_len:.2e} of characteristic length {char_len:.4f})"
            )
    if info.get("is_valid") is False:
        warnings.append("Shape is not valid — may cause meshing issues")
    if warnings:
        info["cae_warnings"] = warnings

    # Export
    if export_dir:
        export_path = Path(export_dir)
        export_path.mkdir(parents=True, exist_ok=True)
        base_name = target_name
        fmt = export_format.lower().strip()
        try:
            if fmt == "step":
                fpath = export_path / f"{base_name}.step"
                target.exportStep(str(fpath))
                info["exported"] = str(fpath)
            elif fmt == "brep":
                fpath = export_path / f"{base_name}.brep"
                target.exportBrep(str(fpath))
                info["exported"] = str(fpath)
            elif fmt == "stl":
                fpath = export_path / f"{base_name}.stl"
                target.exportStl(str(fpath))
                info["exported"] = str(fpath)
            else:
                info["export_error"] = f"Unknown format: {fmt}"
        except Exception as e:
            info["export_error"] = f"{type(e).__name__}: {e}"

    if stdout_text:
        info["stdout"] = stdout_text

    return json.dumps(info, indent=2)


@mcp.tool()
def cadquery_to_cubit_hex(script: str,
                          target_size: float = 1.0,
                          prefer: str = "hex",
                          commit_to_gui: bool = True,
                          timeout_s: int = 600) -> str:
    """End-to-end: cadquery script → STEP → Cubit `cubit_mesh_auto`
    (batch-validated scheme ladder → winning recipe replayed in GUI).

    One call, three libraries, one output. Intended for the flow:
    "I have a cadquery geometry (from cadquery-mcp, a community
    gist, or my own script) — hex-mesh it and show me in Cubit."

    Args:
        script: cadquery script that assigns the final shape to
            `result` (or leaves one Workplane/Shape in the namespace).
        target_size: mesh element size (passed to `volume all size`).
        prefer: "hex" (require hex>0 to accept) or "any".
        commit_to_gui: replay winning recipe in the live Cubit GUI.
        timeout_s: per-rung timeout for the mesh-auto ladder.

    Returns:
        JSON with cadquery info (validity/volume/bbox) + STEP path +
        cubit_mesh_auto result (attempts ladder + winner + GUI replay).
    """
    try:
        import cadquery as cq  # noqa: F401
    except ImportError:
        return json.dumps({"status": "error",
                           "error": "cadquery not installed"}, indent=2)

    import tempfile
    tmp_dir = Path(tempfile.mkdtemp(prefix="cq_to_cubit_"))
    # Step 1: run cadquery → STEP
    cq_info_raw = execute_cadquery(script=script,
                                    export_dir=str(tmp_dir),
                                    export_format="step")
    cq_info = json.loads(cq_info_raw)
    if cq_info.get("status") != "ok" or "exported" not in cq_info:
        return json.dumps({
            "status": "error",
            "stage": "cadquery",
            "cadquery": cq_info,
        }, indent=2)

    step_path = cq_info["exported"].replace("\\", "/")

    # Step 2: hand off to cubit_mesh_auto (imported lazily — keeps
    # build123d server standalone if user installs without cubit)
    try:
        from ..cubit.server import cubit_mesh_auto
    except ImportError as e:
        return json.dumps({
            "status": "error",
            "stage": "cubit_import",
            "error": str(e),
            "cadquery": cq_info,
            "step_path": step_path,
        }, indent=2)

    mesh_raw = cubit_mesh_auto(step_path=step_path,
                                target_size=target_size,
                                prefer=prefer,
                                commit_to_gui=commit_to_gui,
                                timeout_s=timeout_s)
    mesh_info = json.loads(mesh_raw)
    return json.dumps({
        "status": mesh_info.get("status", "ok"),
        "cadquery": {
            k: cq_info[k] for k in
            ("is_valid", "volume", "face_count", "edge_count",
             "bounding_box", "cae_warnings")
            if k in cq_info
        },
        "step_path": step_path,
        "cubit_mesh_auto": mesh_info,
    }, indent=2)


# ============================================================
# Lint: static analysis of build123d scripts
# ============================================================

def _b3d_lint_file(filepath: str) -> list[dict]:
    try:
        with open(filepath, "r", encoding="utf-8", errors="replace") as f:
            lines = f.readlines()
    except OSError as e:
        return [{"line": 0, "severity": "ERROR", "rule": "read-error",
                 "message": f"cannot read: {e}"}]
    findings: list[dict] = []
    for rule in _B3D_LINT_RULES:
        try:
            findings.extend(rule(filepath, lines))
        except Exception as e:  # noqa: BLE001
            findings.append({"line": 0, "severity": "ERROR",
                             "rule": rule.__name__,
                             "message": f"rule crashed: {e}"})
    order = {"CRITICAL": 0, "HIGH": 1, "MODERATE": 2, "LOW": 3, "ERROR": -1}
    findings.sort(key=lambda f: (order.get(f["severity"], 9), f["line"]))
    return findings


@mcp.tool()
def lint_build123d_script(filepath: str) -> str:
    """Static-analysis lint for a single build123d `.py` script.

    Rules (see `radia_mcp.build123d.rules`):
      * `missing-buildpart-context` — `extrude/sweep/revolve` outside
        `with BuildPart() as _bp:` (Builder API) without an assigned
        result (Algebra API).
      * `sweep-no-path` — positional `sweep(sketch)` missing `path=`.
      * `polyline-not-closed` — Polyline in BuildSketch never `close()`-d.
      * `buildsketch-ambiguous-arg` — positional plane arg that can
        collide with `mode` kwarg.
      * `missing-export` — script has geometry but no `export_*` call.
      * `cadquery-in-build123d` — `cq.Workplane(...)` in a b3d file.
      * `micro-fillet-radius` — fillet/chamfer < 1e-3 mm (CAE quality).

    Args:
        filepath: absolute or cwd-relative path to a `.py` file.
    """
    import os as _os
    if not _os.path.isfile(filepath):
        return json.dumps({"status": "error",
                           "error": f"not a file: {filepath}"})
    findings = _b3d_lint_file(filepath)
    return json.dumps({
        "status": "ok",
        "path": filepath,
        "count": len(findings),
        "findings": findings,
    }, indent=2, ensure_ascii=False)


@mcp.tool()
def lint_build123d_directory(directory: str = "examples") -> str:
    """Lint every `.py` file under `directory` (non-recursive default).

    Returns a roll-up of findings per file. Handy before committing
    a batch of build123d scripts to a repo.
    """
    import os as _os
    if not _os.path.isdir(directory):
        return json.dumps({"status": "error",
                           "error": f"not a directory: {directory}"})
    by_file: list[dict] = []
    total = 0
    for name in sorted(_os.listdir(directory)):
        if not name.endswith(".py"):
            continue
        p = _os.path.join(directory, name)
        findings = _b3d_lint_file(p)
        by_file.append({"path": p, "count": len(findings),
                        "findings": findings})
        total += len(findings)
    return json.dumps({
        "status": "ok",
        "directory": directory,
        "total_findings": total,
        "files": by_file,
    }, indent=2, ensure_ascii=False)


# ============================================================
# Script generator: boilerplate templates for common shapes
# ============================================================

_B3D_TEMPLATES: dict[str, str] = {
    "helix_coil": '''"""Helical coil as a single swept prism — hex-meshable in Cubit.

3-body composition (lead_in + helix + lead_out) emits a clean STEP
with 3 topologically-simple prisms. `scheme auto` in Cubit sweeps
each trivially → pure hex mesh.
"""
from build123d import (
    BuildLine, BuildSketch, BuildPart, Line, Helix, Rectangle,
    Plane, Compound, sweep, export_step,
)

R, PITCH, N_TURNS = 15.0, 8.0, 4
H = PITCH * N_TURNS
LEAD_LEN, XS = 20.0, 2.0

p_in_start = (R + LEAD_LEN, 0.0, 0.0)
p_helix_s  = (R, 0.0, 0.0)
p_helix_e  = (R, 0.0, H)
p_out_end  = (R + LEAD_LEN, 0.0, H)

def _one(path_fn, plane):
    with BuildLine() as pth:
        path_fn()
    with BuildSketch(plane) as xs:
        Rectangle(XS, XS)
    with BuildPart() as bp:
        sweep(xs.sketch, path=pth.line)
    return bp.part

lead_in   = _one(lambda: Line(p_in_start, p_helix_s),
                 Plane(origin=p_in_start, x_dir=(0, 1, 0), z_dir=(-1, 0, 0)))
helix_body = _one(lambda: Helix(pitch=PITCH, height=H, radius=R),
                  Plane(origin=p_helix_s, x_dir=(1, 0, 0), z_dir=(0, 1, 0)))
lead_out  = _one(lambda: Line(p_helix_e, p_out_end),
                 Plane(origin=p_helix_e, x_dir=(0, 1, 0), z_dir=(1, 0, 0)))

coil = Compound(children=[lead_in, helix_body, lead_out])
export_step(coil, "coil.step")
''',

    "l_bracket": '''"""L-shaped bracket — 2D polyline extruded, clean prism for hex."""
from build123d import (
    BuildPart, BuildSketch, BuildLine, Polyline, make_face, extrude,
    export_step,
)

pts = [(0, 0), (60, 0), (60, 10), (10, 10), (10, 40), (0, 40), (0, 0)]
with BuildPart() as bp:
    with BuildSketch() as sk:
        with BuildLine() as ln:
            Polyline(*pts)
        make_face()
    extrude(amount=15)
part = bp.part
export_step(part, "bracket.step")
''',

    "cae_block": '''"""CAE-oriented uniform block — named region + clean hex candidate."""
from build123d import BuildPart, Box, export_step

with BuildPart() as bp:
    Box(30, 20, 10)
part = bp.part
part.label = "magnet"          # carries into STEP as the block name
export_step(part, "magnet.step")
''',

    "gear_bd_warehouse": '''"""Spur gear via bd_warehouse (install `pip install bd-warehouse`)."""
from build123d import export_step
from bd_warehouse.gear import SpurGear

gear = SpurGear(
    module=1.0, tooth_count=20, thickness=5.0,
    pressure_angle=20.0, root_fillet=0.3,
)
export_step(gear, "gear.step")
''',

    "fastener_assembly": '''"""Bolted plate — plate with hole + bd_warehouse SocketHeadCapScrew."""
from build123d import BuildPart, Box, Cylinder, Mode, Location, Plane
from bd_warehouse.fastener import SocketHeadCapScrew
from build123d import export_step

with BuildPart() as plate:
    Box(40, 40, 8)
    with BuildPart(Plane.XY.offset(4), mode=Mode.SUBTRACT) as _hole:
        Cylinder(radius=3.2, height=20)

screw = SocketHeadCapScrew(size="M6-1", length=20)
screw.locate(Location((0, 0, 12)))

export_step(plate.part, "plate.step")
export_step(screw, "screw.step")
''',

    "sweep_square_path": '''"""Sweep a square cross-section along a custom path — hex-friendly."""
from build123d import (
    BuildLine, BuildSketch, BuildPart, Line, Spline, Rectangle,
    Plane, sweep, export_step,
)

with BuildLine() as path:
    Line((0, 0, 0), (20, 0, 0))
    Spline([(20, 0, 0), (30, 10, 0), (40, 10, 10)])

with BuildSketch(Plane.YZ) as xs:
    Rectangle(5, 5)

with BuildPart() as bp:
    sweep(xs.sketch, path=path.line)

part = bp.part
export_step(part, "swept.step")
''',

    # ------------------------------------------------------------------
    # Radia-specific templates — magnet arrays, iron cores, coils
    # ------------------------------------------------------------------

    "magnet_ring": '''"""Ring of cuboid magnets arranged around the Z axis.

Each cuboid is labelled so STEP → Cubit preserves it as a block name.
Typical use: rotor / stator magnets in a motor, or seed geometry for
Halbach variants. All volumes are prismatic → hex-friendly.
"""
from build123d import BuildPart, Box, Compound, Location, Plane, export_step
import math

N_MAGNETS = 12        # number of segments
R_INNER = 20.0        # mm — ring inner radius
MAG_W, MAG_H, MAG_D = 8.0, 5.0, 10.0   # width (tangential), height (radial), depth (axial)

magnets = []
for i in range(N_MAGNETS):
    angle_deg = 360.0 * i / N_MAGNETS
    angle = math.radians(angle_deg)
    cx = (R_INNER + MAG_H / 2) * math.cos(angle)
    cy = (R_INNER + MAG_H / 2) * math.sin(angle)
    with BuildPart() as m:
        Box(MAG_W, MAG_H, MAG_D)
    m.part.label = f"magnet_{i:02d}"
    # Rotate around Z by angle_deg, then translate to ring position
    m.part.locate(Location((cx, cy, 0), (0, 0, angle_deg)))
    magnets.append(m.part)

ring = Compound(children=magnets)
export_step(ring, "magnet_ring.step")
''',

    "halbach_array": '''"""Linear Halbach array — N segments with M rotating 360°/N per segment.

The magnetization direction is carried only as a label (build123d is
pure CAD); downstream magnet solver (Radia) reads back by label.
Rotating M in 90° steps (N_PERIODS × 4 segs/period) gives the classic
strong-side / null-side pattern.
"""
from build123d import BuildPart, Box, Compound, Location, export_step

SEG_L, SEG_W, SEG_H = 10.0, 20.0, 10.0   # length (x) × width (y) × height (z)
N_PERIODS = 3
SEGS_PER_PERIOD = 4
N = N_PERIODS * SEGS_PER_PERIOD

segments = []
for i in range(N):
    with BuildPart() as seg:
        Box(SEG_L, SEG_W, SEG_H)
    # Magnetization direction cycles: +Z, +X, -Z, -X, ...
    m_dir = (90 * i) % 360                # degrees
    seg.part.label = f"halbach_{i:02d}_mdir{m_dir}"
    seg.part.locate(Location((SEG_L * (i + 0.5), 0, 0)))
    segments.append(seg.part)

array = Compound(children=segments)
export_step(array, "halbach_array.step")
''',

    "c_core": '''"""C-shaped iron core with an air gap. Classic dipole magnet yoke.

Cross-section is a C, extruded along Y.  The gap sits on +X side.
Labels: 'core' for the C-body, 'air_gap' for the subtracted gap box
(kept as a separate volume for named-region meshing in Cubit).
"""
from build123d import (
    BuildPart, BuildLine, BuildSketch, Polyline, make_face, extrude,
    Box, Compound, Location, Mode, export_step,
)

# Outer C envelope (in XZ plane)
# Coordinates (x, z): trace the C clockwise
outer = [(0, 0), (60, 0), (60, 50), (0, 50)]
inner = [(10, 10), (40, 10), (40, 25), (60, 25), (60, 35), (40, 35), (40, 40), (10, 40)]
THICKNESS_Y = 20.0

with BuildPart() as core:
    with BuildSketch() as sk:
        with BuildLine() as ln_out:
            Polyline(*outer, close=True)
        make_face()
        with BuildLine() as ln_in:
            Polyline(*inner, close=True)
        make_face(mode=Mode.SUBTRACT)
    extrude(amount=THICKNESS_Y, both=True)   # ±Y half-thickness each
core.part.label = "core"

# Air gap volume (between pole faces at z=25..35 inside core).
# Box default is Align.CENTER on all axes — no need to pass align.
with BuildPart() as gap:
    Box(20, THICKNESS_Y * 2, 10)
gap.part.locate(Location((50, 0, 30)))
gap.part.label = "air_gap"

asm = Compound(children=[core.part, gap.part])
export_step(asm, "c_core.step")
''',

    "e_core": '''"""E-shaped transformer core (central leg + two outer legs).

Mirror pair of C-cores sharing a back iron — common in single-phase
transformers. Three legs labelled 'leg_center', 'leg_left', 'leg_right'
for named-winding-region meshing.
"""
from build123d import (
    BuildPart, Box, Compound, Location, export_step,
)

BACK_W, BACK_H, DEPTH = 80.0, 12.0, 30.0
LEG_W, LEG_H = 12.0, 40.0
CENTER_LEG_W = 20.0                       # wider central leg

with BuildPart() as back:
    Box(BACK_W, BACK_H, DEPTH)
back.part.locate(Location((0, LEG_H / 2 + BACK_H / 2, 0)))
back.part.label = "back_iron"

with BuildPart() as center:
    Box(CENTER_LEG_W, LEG_H, DEPTH)
center.part.label = "leg_center"

with BuildPart() as left:
    Box(LEG_W, LEG_H, DEPTH)
left.part.locate(Location((-(BACK_W / 2 - LEG_W / 2), 0, 0)))
left.part.label = "leg_left"

with BuildPart() as right:
    Box(LEG_W, LEG_H, DEPTH)
right.part.locate(Location((BACK_W / 2 - LEG_W / 2, 0, 0)))
right.part.label = "leg_right"

asm = Compound(children=[back.part, center.part, left.part, right.part])
export_step(asm, "e_core.step")
''',

    "pole_piece": '''"""Tapered pole shoe (wedge-on-box) — shapes field at the gap.

Rectangular base + a frustum cap that narrows toward the pole face.
Useful as the working face of a dipole C-core or MRI pole.
"""
from build123d import (
    BuildPart, BuildSketch, Rectangle, extrude, loft, Plane, Compound,
    Location, export_step,
)

BASE_W, BASE_D, BASE_H = 50.0, 50.0, 10.0
TIP_W, TIP_D = 30.0, 30.0
FRUSTUM_H = 15.0

with BuildPart() as pole:
    # Base slab
    with BuildSketch() as base_sk:
        Rectangle(BASE_W, BASE_D)
    extrude(amount=BASE_H)
    # Frustum via loft between two rectangles
    with BuildSketch(Plane.XY.offset(BASE_H)) as s0:
        Rectangle(BASE_W, BASE_D)
    with BuildSketch(Plane.XY.offset(BASE_H + FRUSTUM_H)) as s1:
        Rectangle(TIP_W, TIP_D)
    loft()

pole.part.label = "pole_piece"
export_step(pole.part, "pole_piece.step")
''',

    "stator_lamination": '''"""Slotted annular stator lamination — donut with N rectangular
slots cut inward from the bore. Single lamination thickness; stack
for a real stator core.

Approach: build the annulus as a Part, then loop-subtract each slot
as a separately-placed Box. Keeps each operation simple and keeps
the BuildSketch context uncluttered.
"""
from build123d import (
    BuildPart, BuildSketch, Circle, Mode, Box, Location, extrude,
    export_step,
)
import math

R_OUT, R_IN = 60.0, 30.0
N_SLOTS = 24
SLOT_TANG = 4.0             # tangential slot width at the bore
SLOT_RADIAL = 10.0          # radial depth of the slot (into the iron)
LAM_THICK = 0.5             # single lamination thickness

with BuildPart() as lam:
    # 1. Annular disk
    with BuildSketch() as sk:
        Circle(R_OUT)
        Circle(R_IN, mode=Mode.SUBTRACT)
    extrude(amount=LAM_THICK)

    # 2. Carve slots — one SUBTRACT Box per slot, rotated + placed
    for i in range(N_SLOTS):
        angle_deg = 360.0 * i / N_SLOTS
        angle = math.radians(angle_deg)
        cx = (R_IN + SLOT_RADIAL / 2) * math.cos(angle)
        cy = (R_IN + SLOT_RADIAL / 2) * math.sin(angle)
        with BuildPart(Location((cx, cy, LAM_THICK / 2),
                                 (0, 0, angle_deg)),
                       mode=Mode.SUBTRACT):
            Box(SLOT_TANG, SLOT_RADIAL, LAM_THICK * 1.1)

lam.part.label = "stator_lamination"
export_step(lam.part, "stator.step")
''',

    "racetrack_coil": '''"""Racetrack coil (obround sweep) — two straight segments joined by
semicircles at each end, swept with a rectangular cross-section.
Four-body decomposition keeps each segment a clean prism for hex
meshing (matches the helix_coil pattern).
"""
from build123d import (
    BuildLine, BuildSketch, BuildPart, Line, RadiusArc, Rectangle,
    Plane, Compound, sweep, export_step,
)

L_STRAIGHT = 40.0          # mm, straight section length
R_BEND = 10.0              # mm, semicircle radius
XS_W, XS_H = 4.0, 3.0      # cross-section (width radial, height axial)

# Four path segments in the XY plane at z=0
p0 = (-L_STRAIGHT / 2, -R_BEND, 0)
p1 = ( L_STRAIGHT / 2, -R_BEND, 0)
p2 = ( L_STRAIGHT / 2,  R_BEND, 0)
p3 = (-L_STRAIGHT / 2,  R_BEND, 0)

def _body(path_points_or_fn, plane):
    with BuildLine() as pth:
        if callable(path_points_or_fn):
            path_points_or_fn()
        else:
            a, b = path_points_or_fn
            Line(a, b)
    with BuildSketch(plane) as xs:
        Rectangle(XS_W, XS_H)
    with BuildPart() as bp:
        sweep(xs.sketch, path=pth.line)
    return bp.part

# Two straight sides
bottom = _body((p0, p1),
               Plane(origin=p0, x_dir=(0, 0, 1), z_dir=(1, 0, 0)))
top    = _body((p2, p3),
               Plane(origin=p2, x_dir=(0, 0, 1), z_dir=(-1, 0, 0)))
# Two semicircular arcs (need RadiusArc inside BuildLine)
right = _body(lambda: RadiusArc(p1, p2, R_BEND),
              Plane(origin=p1, x_dir=(0, 0, 1), z_dir=(0, -1, 0)))
left  = _body(lambda: RadiusArc(p3, p0, R_BEND),
              Plane(origin=p3, x_dir=(0, 0, 1), z_dir=(0, 1, 0)))

coil = Compound(children=[bottom, right, top, left])
export_step(coil, "racetrack_coil.step")
''',
}


@mcp.tool()
def generate_build123d_script(pattern: str = "helix_coil") -> str:
    """Return a ready-to-run build123d boilerplate for a common pattern.

    **General patterns**:
      * `helix_coil` — 4-turn helix with radial leads (hex-friendly
        3-body STEP; the same structure used in `spiral_coil_v3.py`).
      * `l_bracket` — 2D polyline + extrude, clean prism.
      * `cae_block` — simple Box with `.label` set (survives STEP
        round-trip as a block name in Cubit).
      * `gear_bd_warehouse` — SpurGear from bd_warehouse.
      * `fastener_assembly` — plate with hole + M6 SocketHeadCapScrew.
      * `sweep_square_path` — sweep rectangle along Line+Spline path.

    **Radia-specific magnet / iron / coil patterns**:
      * `magnet_ring` — N cuboid magnets around a Z-axis ring
        (rotor seed; each magnet labelled per-index).
      * `halbach_array` — linear N-period Halbach with M-direction
        carried as a label (downstream solver reads by name).
      * `c_core` — iron dipole yoke with explicit air-gap volume.
      * `e_core` — single-phase transformer E-core (3 legs +
        back iron, each labelled).
      * `pole_piece` — tapered pole shoe (base slab + loft frustum).
      * `stator_lamination` — slotted annular lamination (24 slots
        by default, single-thickness; stack for a real stator).
      * `racetrack_coil` — obround coil (2 straight + 2 arc segments,
        4-body decomposition for hex meshing).

    Returns the script text so the caller can feed it directly into
    `execute_build123d` / `build123d_try` or save to disk.
    """
    key = (pattern or "helix_coil").strip().lower()
    if key not in _B3D_TEMPLATES:
        return json.dumps({"status": "error",
                           "error": f"unknown pattern {pattern!r}",
                           "known": list(_B3D_TEMPLATES.keys())})
    return json.dumps({
        "status": "ok",
        "pattern": key,
        "script": _B3D_TEMPLATES[key],
    }, indent=2, ensure_ascii=False)


# ============================================================
# Subprocess-isolated dry-run (build123d_try)
# ============================================================

@mcp.tool()
def build123d_try(script: str, timeout_s: int = 90) -> str:
    """Dry-run a build123d script in a **fresh Python subprocess**.

    Unlike `execute_build123d` (same process, mutates namespace),
    this spawns a clean interpreter so:
      - Import state cannot leak between calls.
      - An OCCT segfault / abort on bad geometry is contained.
      - `sys.modules` pollution during iterative prompt engineering
        doesn't accumulate.

    Captures stdout/stderr and reports subprocess return code. The
    script runs in a temp dir so any export paths are sandboxed.

    Args:
        script: Python using build123d (self-contained).
        timeout_s: subprocess timeout (default 90).
    """
    import subprocess as _sp
    import tempfile as _tmp
    preflight = _b3d_preflight_warn(script)
    tmp = Path(_tmp.mkdtemp(prefix="b3d_try_"))
    script_path = tmp / "_try.py"
    wrapped = (
        "import sys, json, traceback\n"
        "try:\n"
        "    from build123d import *\n"
        + "\n".join("    " + ln for ln in script.splitlines())
        + "\n"
        "    print('__B3D_TRY_OK__')\n"
        "except Exception as e:\n"
        "    print('__B3D_TRY_ERR__')\n"
        "    traceback.print_exc()\n"
    )
    script_path.write_text(wrapped, encoding="utf-8")
    try:
        r = _sp.run(
            [sys.executable, str(script_path)],
            cwd=str(tmp),
            capture_output=True, text=True,
            timeout=timeout_s,
        )
    except _sp.TimeoutExpired as e:
        return json.dumps({"status": "error", "stage": "timeout",
                           "error": f"exceeded {timeout_s}s",
                           "stdout": e.stdout or "",
                           "stderr": e.stderr or ""}, indent=2)
    except OSError as e:
        return json.dumps({"status": "error", "stage": "spawn",
                           "error": str(e)}, indent=2)

    out = r.stdout or ""
    err = r.stderr or ""
    passed = (r.returncode == 0) and ("__B3D_TRY_OK__" in out)
    payload = {
        "status": "ok" if passed else "error",
        "returncode": r.returncode,
        "stdout": out[-4000:],
        "stderr": err[-4000:],
        "artefacts": [p.name for p in tmp.iterdir() if p.name != "_try.py"],
        "tmpdir": str(tmp),
    }
    if preflight:
        payload["preflight"] = preflight
    if not passed:
        _fl.record_failure("build123d", {
            "input": script,
            "error": err[-1200:] if err else out[-1200:],
            "via": "build123d_try",
        })
    return json.dumps(payload, indent=2, ensure_ascii=False)


# ============================================================
# Parallel race: try N build123d variants concurrently, best/first wins
# ============================================================

@mcp.tool()
def build123d_try_race(scripts: list,
                        max_concurrent: int = 4,
                        timeout_s: int = 90,
                        prefer: str = "valid_largest_volume") -> str:
    """Race N build123d script variants in parallel subprocesses,
    return the per-variant report + chosen winner.

    Each variant runs in its own clean Python (`build123d_try` style)
    so a degenerate-geometry crash on one variant doesn't kill the
    others. Use this for parameter sweeps — "try 3 fillet radii in
    parallel and keep whichever is valid + has the largest volume."

    Args:
        scripts: list of {"name": str, "script": str} dicts (or plain
            strings — auto-named).
        max_concurrent: subprocess pool size.
        timeout_s: per-variant timeout.
        prefer: winner-selection rule:
          * "valid_largest_volume" — must be valid; pick max volume.
          * "first_valid" — first one that's valid (race semantics).
          * "first_ok" — first that finishes without error.

    Returns JSON: list of per-variant info + chosen winner index.
    """
    import concurrent.futures as _cf

    if not scripts:
        return json.dumps({"status": "error",
                           "error": "scripts list cannot be empty"})

    # Normalize entries
    norm: list[dict] = []
    for i, s in enumerate(scripts):
        if isinstance(s, str):
            norm.append({"name": f"variant_{i}", "script": s})
        elif isinstance(s, dict) and "script" in s:
            norm.append({"name": s.get("name", f"variant_{i}"),
                          "script": s["script"]})
        else:
            return json.dumps({"status": "error",
                               "error": f"scripts[{i}] must be str or {{name, script}}"})

    def _one(rec: dict) -> dict:
        raw = build123d_try(script=rec["script"], timeout_s=timeout_s)
        info = json.loads(raw)
        # build123d_try returns artefacts + status; we also want
        # validity / volume which come from execute_build123d. Run
        # that too for accepted variants — only inside the same
        # subprocess pool (cheap, in-process, but we already paid the
        # subprocess cost above; for simplicity reuse the in-process
        # execute_build123d which gives us validity/volume directly).
        if info.get("status") != "ok":
            return {**rec, "status": "error", "info": info,
                    "is_valid": None, "volume": None}
        # Now in-process inspect for validity / volume
        inspect_raw = _execute_build123d_sync(script=rec["script"], export_dir="")
        inspect = json.loads(inspect_raw)
        return {
            "name": rec["name"],
            "status": inspect.get("status"),
            "is_valid": inspect.get("is_valid"),
            "volume": inspect.get("volume"),
            "face_count": inspect.get("face_count"),
            "bounding_box": inspect.get("bounding_box"),
            "subprocess_artefacts": info.get("artefacts", []),
        }

    results: list[dict] = []
    with _cf.ThreadPoolExecutor(max_workers=max(1, int(max_concurrent))) as pool:
        futures = {pool.submit(_one, r): r for r in norm}
        for fut in _cf.as_completed(futures):
            results.append(fut.result())

    # Stable order by submission index, but we lost it across threads;
    # re-key by name → keep results as-is, callers can sort.

    # Pick winner per `prefer` rule
    winner_idx = None
    rule = (prefer or "valid_largest_volume").strip().lower()
    if rule == "first_ok":
        for i, r in enumerate(results):
            if r.get("status") == "ok":
                winner_idx = i
                break
    elif rule == "first_valid":
        for i, r in enumerate(results):
            if r.get("is_valid") is True:
                winner_idx = i
                break
    else:  # valid_largest_volume
        valid = [(i, r) for i, r in enumerate(results)
                 if r.get("is_valid") is True
                 and isinstance(r.get("volume"), (int, float))]
        if valid:
            winner_idx = max(valid, key=lambda iv: iv[1]["volume"])[0]

    return json.dumps({
        "status": "ok",
        "rule": rule,
        "results": results,
        "winner_index": winner_idx,
        "winner": results[winner_idx] if winner_idx is not None else None,
    }, indent=2, ensure_ascii=False)


# ============================================================
# One-shot pipeline: build123d script → STEP → cubit_mesh_auto
# ============================================================

@mcp.tool()
def build123d_to_cubit_hex(script: str,
                           target_size: float = 1.0,
                           prefer: str = "hex",
                           commit_to_gui: bool = True,
                           timeout_s: int = 600) -> str:
    """End-to-end: build123d script → STEP → `cubit_mesh_auto` (batch-
    validated scheme ladder → winner replayed in live Cubit GUI).

    Mirror of `cadquery_to_cubit_hex` for the build123d side. Use this
    when the user says "make a hex mesh of <build123d shape>" and you
    want one-call coverage of the whole pipeline.

    Args:
        script: build123d script; must leave a Shape-like object in
            the namespace (or assign to `part`) so it can be exported.
        target_size: mesh size (`volume all size`).
        prefer: "hex" (need hex>0 to accept) or "any".
        commit_to_gui: replay the winning recipe in the live Cubit GUI.
        timeout_s: per-rung timeout.

    Returns JSON with: build123d info (validity/volume/bbox) + STEP
    path + cubit_mesh_auto ladder report + GUI replay state.
    """
    import tempfile
    tmp_dir = Path(tempfile.mkdtemp(prefix="b3d_to_cubit_"))
    exec_raw = _execute_build123d_sync(script=script,
                                        export_dir=str(tmp_dir),
                                        export_format="step")
    exec_info = json.loads(exec_raw)
    if exec_info.get("status") != "ok" or "exported" not in exec_info:
        return json.dumps({"status": "error",
                           "stage": "build123d",
                           "build123d": exec_info}, indent=2)
    step_path = exec_info["exported"].replace("\\", "/")

    try:
        from ..cubit.server import cubit_mesh_auto
    except ImportError as e:
        return json.dumps({"status": "error",
                           "stage": "cubit_import",
                           "error": str(e),
                           "build123d": exec_info,
                           "step_path": step_path}, indent=2)

    mesh_raw = cubit_mesh_auto(step_path=step_path,
                                target_size=target_size,
                                prefer=prefer,
                                commit_to_gui=commit_to_gui,
                                timeout_s=timeout_s)
    mesh_info = json.loads(mesh_raw)
    return json.dumps({
        "status": mesh_info.get("status", "ok"),
        "build123d": {
            k: exec_info[k] for k in
            ("is_valid", "volume", "face_count", "edge_count",
             "bounding_box", "cae_warnings")
            if k in exec_info
        },
        "step_path": step_path,
        "cubit_mesh_auto": mesh_info,
    }, indent=2, ensure_ascii=False)


# ============================================================
# State-aware next-step suggestions
# ============================================================

@mcp.tool()
def build123d_suggest_next(goal: str = "cae_pipeline",
                           script: str = "") -> str:
    """Suggest concrete next build123d steps toward a goal.

    Rule-based, fast, no LLM. Inspects the provided `script` (empty
    string = generic advice) and proposes 2-5 next edits with rationale.

    Goals:
      * `cae_pipeline` — get to a clean STEP ready for Cubit hex mesh
      * `clean` — fix common geometry issues (micro-edges, open polys)
      * `export` — add exports and naming
      * `refine` — improve mesh quality at the CAD boundary
      * `assemble` — combine parts / add joints / fasteners

    Args:
        goal: one of the keys above.
        script: optional script text; if provided, suggestions adapt
            to what is / isn't already present.
    """
    g = (goal or "cae_pipeline").strip().lower()
    s = script or ""
    has_bp = "with BuildPart" in s or "BuildPart(" in s
    has_export = any(fn in s for fn in ("export_step", "export_brep",
                                         "export_stl"))
    has_label = ".label" in s
    has_fillet = "fillet(" in s or "chamfer(" in s

    suggestions: list[dict] = []
    def add(cmd: str, why: str) -> None:
        suggestions.append({"step": cmd, "why": why})

    if g == "cae_pipeline":
        if not has_bp:
            add("with BuildPart() as bp: ...",
                "Wrap geometry in BuildPart so `.part` / `.sketch` are "
                "well-defined and sweep/extrude see a pending curve.")
        if not has_export:
            add('export_step(bp.part, "out.step")',
                "Pipeline needs a STEP handoff. Cubit imports STEP via "
                "`import step \"...\"` downstream.")
        if not has_label:
            add('part.label = "magnet"',
                "Labels survive STEP round-trip as Cubit block names — "
                "critical for named physical groups.")
        add("build123d_try(script)  # before hitting Cubit",
            "Run in a subprocess first to catch OCCT errors without "
            "polluting the live MCP session.")
        add("build123d_to_cubit_hex(script, target_size=1.0)",
            "One-shot pipeline: build123d → STEP → Cubit scheme ladder "
            "→ live GUI replay.")
    elif g == "clean":
        add("part.is_valid",
            "Check OCCT validity before export; invalid parts break mesh.")
        if has_fillet:
            add("lint_build123d_script(path)",
                "Lint checks fillet/chamfer radius vs characteristic "
                "length; micro-edges ruin mesh quality.")
        add("part.edges().filter_by(lambda e: e.length < 1e-3)",
            "Enumerate micro-edges; remove or redesign.")
    elif g == "export":
        add('export_step(part, "out.step")',
            "Default exchange format for Cubit / ngsolve.")
        add('export_brep(part, "out.brep")',
            "Lossless OCCT native; use when STEP round-trip loses "
            "parametric data.")
        add('export_stl(part, "out.stl", tolerance=0.01)',
            "Visualization / printing only — loses curves.")
    elif g == "refine":
        add("fillet(part.edges(), radius=0.5)",
            "Soft small edges so mesh has no stress concentrations.")
        add("part.faces().filter_by(Axis.Z)",
            "Select faces for named groups before export (better block "
            "names downstream).")
    elif g == "assemble":
        add("from bd_warehouse.fastener import SocketHeadCapScrew",
            "Use bd_warehouse for industrial fasteners/bearings/gears.")
        add("part.locate(Location((x, y, z), (rx, ry, rz)))",
            "Absolute placement via Location for assembly poses.")
        add("RigidJoint(label='a_to_b', to_part=b, joint_location=...)",
            "For kinematic/assembly relationships across parts.")
    else:
        add("# unknown goal; try: cae_pipeline / clean / export / refine / assemble",
            "")

    return json.dumps({
        "status": "ok",
        "goal": g,
        "script_analysis": {
            "has_BuildPart_context": has_bp,
            "has_export": has_export,
            "has_label": has_label,
            "has_fillet_or_chamfer": has_fillet,
        },
        "suggestions": suggestions,
    }, indent=2, ensure_ascii=False)


# ============================================================
# API reference — focused lookup over the auto-generated reference
# ============================================================

@mcp.tool()
def build123d_api(query: str, max_results: int = 5) -> str:
    """Search the bundled build123d API reference (auto-generated from
    `inspect.getmembers(build123d)` — 142 classes + 65 functions).

    Unlike `build123d_usage` (returns a whole topic) or
    `build123d_lookup` (tf-idf over all knowledge), this tool is
    **API-first**: results are per-class / per-function entries
    with ctor signature + docstring + method list.

    Use when the user mentions a specific class or function name
    and you want signature + behaviour info without pulling full
    topics.

    Args:
        query: class or function name, or free-form keywords.
        max_results: cap on chunks returned.
    """
    q = (query or "").strip()
    if not q:
        return json.dumps({"status": "error", "error": "empty query"})
    # Reuse the existing build123d_lookup mechanism but restrict to
    # the api_reference topic by directly chunking + scoring it.
    from .build123d_knowledge import _TOPICS
    api_text = _TOPICS.get("api_reference", "")
    if not api_text:
        return json.dumps({"status": "error",
                           "error": "api_reference topic missing"})
    chunks = _chunk_by_heading(api_text, target_lines=20)
    terms = [t for t in _tokenize(q) if t not in _STOP_WORDS]
    if not terms:
        return json.dumps({"status": "error",
                            "error": "no searchable terms (all stop words?)"})
    docs = []
    for start, text in chunks:
        tokens = _tokenize(text)
        if not tokens:
            continue
        docs.append({
            "start_line": start,
            "chunk": text,
            "heading": _heading_of(text),
            "tokens": tokens,
            "token_set": set(tokens),
        })
    n = max(1, len(docs))
    idf = {
        t: _math.log((n + 1) / (sum(1 for d in docs if t in d["token_set"]) + 1)) + 1.0
        for t in terms
    }
    scored = []
    for doc in docs:
        s = _score_tfidf(terms, doc, idf)
        if s > 0:
            scored.append((s, doc))
    scored.sort(key=lambda x: (-x[0], x[1]["start_line"]))
    top = scored[: max(1, int(max_results))]
    return json.dumps({
        "status": "ok",
        "query": q,
        "terms": terms,
        "total_scored": len(scored),
        "returned": len(top),
        "results": [
            {
                "start_line": d["start_line"],
                "heading": d["heading"][:120],
                "score": round(s, 3),
                "chunk": d["chunk"],
            }
            for s, d in top
        ],
    }, indent=2, ensure_ascii=False)


# ============================================================
# STEP inspect + heal (Cubit-side sanity gate)
# ============================================================

@mcp.tool()
def build123d_inspect_step(path: str) -> str:
    """Inspect an external STEP file via the build123d / OCCT importer
    and report a CAE-readiness summary.

    Reports, per root shape (the STEP may contain a Compound of many):
      * label (if the source preserved it)
      * type (Solid / Compound / Shell / …)
      * is_valid (OCCT geometric validity check)
      * volume / area / face_count / edge_count / vertex_count
      * bounding_box + characteristic length
      * min_edge_length + micro-edge ratio (flags if < 1e-3)
      * cae_warnings (invalid shape, micro-edges, suspect topology)

    Use this to gate a STEP file before handing to Cubit — catches
    noisy externals (cadquery_to_cubit_hex output, vendor STEPs)
    before meshing fails on it.

    Args:
        path: absolute or cwd-relative STEP file.
    """
    p = Path(path)
    if not p.exists():
        return json.dumps({"status": "error",
                           "error": f"not found: {path}"})
    if not p.is_file():
        return json.dumps({"status": "error",
                           "error": f"not a file: {path}"})
    try:
        from build123d import import_step, Compound, Shape
    except ImportError as e:
        return json.dumps({"status": "error",
                           "error": f"build123d unavailable: {e}"})
    try:
        imported = import_step(str(p))
    except Exception as e:
        return json.dumps({"status": "error",
                           "stage": "import",
                           "error": f"{type(e).__name__}: {e}"})

    # Flatten: a Compound yields children; a plain Shape yields itself
    shapes: list = []
    if isinstance(imported, Compound) and imported.children:
        shapes = list(imported.children)
    else:
        shapes = [imported]

    summaries: list[dict] = []
    warnings: list[str] = []
    for i, sh in enumerate(shapes):
        info: dict = {"index": i, "type": type(sh).__name__,
                      "label": (sh.label if getattr(sh, "label", None)
                                else "(none)")}
        try:
            info["is_valid"] = bool(sh.is_valid)
        except Exception:
            info["is_valid"] = None
        try:
            info["volume"] = round(float(sh.volume), 6)
        except Exception:
            info["volume"] = None
        try:
            info["area"] = round(float(sh.area), 6)
        except Exception:
            info["area"] = None
        try:
            faces = sh.faces()
            edges = sh.edges()
            info["face_count"] = len(faces)
            info["edge_count"] = len(edges)
            info["vertex_count"] = len(sh.vertices())
            if edges:
                lens = [float(e.length) for e in edges]
                info["min_edge_length"] = round(min(lens), 8)
                info["max_edge_length"] = round(max(lens), 8)
        except Exception:
            pass
        try:
            bb = sh.bounding_box()
            info["bounding_box"] = {
                "min": [round(bb.min.X, 6), round(bb.min.Y, 6),
                        round(bb.min.Z, 6)],
                "max": [round(bb.max.X, 6), round(bb.max.Y, 6),
                        round(bb.max.Z, 6)],
                "size": [round(v, 6) for v in bb.size],
            }
            char_len = max(bb.size) if bb.size else 1
            info["characteristic_length"] = round(float(char_len), 6)
            min_e = info.get("min_edge_length")
            if min_e and char_len > 0:
                ratio = min_e / char_len
                info["micro_edge_ratio"] = round(ratio, 8)
                if ratio < 1e-3:
                    warnings.append(
                        f"shape[{i}] micro-edge: ratio {ratio:.2e} "
                        f"(min_edge {min_e:.2e} / char_len {char_len:.4f})"
                    )
        except Exception:
            pass
        if info.get("is_valid") is False:
            warnings.append(f"shape[{i}] OCCT is_valid=False")
        if info.get("volume") is not None and info["volume"] <= 0:
            warnings.append(f"shape[{i}] zero / negative volume "
                            f"({info['volume']})")
        summaries.append(info)

    return json.dumps({
        "status": "ok",
        "path": str(p.resolve()),
        "size_bytes": p.stat().st_size,
        "root_type": type(imported).__name__,
        "child_count": len(shapes),
        "shapes": summaries,
        "cae_warnings": warnings,
        "recommendation": (
            "ready for Cubit" if not warnings else
            "run `build123d_heal(step_in, step_out)` and re-inspect "
            "before meshing"),
    }, indent=2, ensure_ascii=False)


@mcp.tool()
def build123d_heal(step_in: str, step_out: str = "",
                   fix_small_edges: bool = True,
                   fix_faces: bool = True) -> str:
    """Run OCCT `ShapeFix_Shape` on a STEP file, write a healed copy.

    Typical fixes: zero-length edges, small gaps between adjacent
    faces, wrong face orientation, redundant vertices, degenerate
    wires. Skips anything ShapeFix can't address; the output is
    always a valid STEP (possibly unchanged from the input).

    Re-inspect after healing via `build123d_inspect_step(step_out)`
    to confirm warnings are cleared.

    Args:
        step_in: source STEP file.
        step_out: destination. Defaults to `<step_in>__healed.step`.
        fix_small_edges: run the small-edge elimination pass.
        fix_faces: run face-orientation + face-gap fixes.
    """
    p_in = Path(step_in)
    if not p_in.exists():
        return json.dumps({"status": "error",
                           "error": f"not found: {step_in}"})
    p_out = Path(step_out) if step_out else \
        p_in.with_name(p_in.stem + "__healed.step")
    p_out.parent.mkdir(parents=True, exist_ok=True)

    try:
        from build123d import import_step, export_step, Compound
    except ImportError as e:
        return json.dumps({"status": "error",
                           "error": f"build123d unavailable: {e}"})

    # ShapeFix lives in OCP (build123d's OCCT binding).
    try:
        from OCP.ShapeFix import ShapeFix_Shape
        from OCP.TopoDS import TopoDS_Shape
    except ImportError as e:
        return json.dumps({"status": "error",
                           "error": f"OCP ShapeFix unavailable: {e}"})

    try:
        imported = import_step(str(p_in))
    except Exception as e:
        return json.dumps({"status": "error", "stage": "import",
                           "error": f"{type(e).__name__}: {e}"})

    # Accept Compound-of-shapes or a single shape
    children: list = (list(imported.children)
                      if isinstance(imported, Compound) and imported.children
                      else [imported])

    healed_children = []
    per_shape_reports: list[dict] = []
    for i, sh in enumerate(children):
        try:
            raw_shape: TopoDS_Shape = sh.wrapped
        except Exception as e:
            per_shape_reports.append({"index": i, "status": "skip",
                                      "reason": f"no .wrapped: {e}"})
            healed_children.append(sh)
            continue
        fixer = ShapeFix_Shape(raw_shape)
        try:
            if fix_small_edges:
                # ShapeFix tolerance — kept conservative
                fixer.SetPrecision(1e-6)
                fixer.SetMaxTolerance(1e-3)
            fixer.Perform()
            fixed_raw = fixer.Shape()
        except Exception as e:
            per_shape_reports.append({"index": i, "status": "fixer_error",
                                      "error": str(e)[:200]})
            healed_children.append(sh)
            continue
        # Re-wrap as a build123d Shape so export_step accepts it
        try:
            new_shape = type(sh)(fixed_raw)
        except Exception:
            # Fallback: subclass construction varies — use base class
            try:
                from build123d import Shape as _Shape
                new_shape = _Shape(fixed_raw)
            except Exception:
                new_shape = sh
        if getattr(sh, "label", None):
            try:
                new_shape.label = sh.label
            except Exception:
                pass
        healed_children.append(new_shape)
        per_shape_reports.append({
            "index": i,
            "status": "healed",
            "label": getattr(sh, "label", "") or "(none)",
        })

    # Re-package into a Compound when we started with one
    if isinstance(imported, Compound) and imported.children:
        try:
            out_shape = Compound(children=healed_children)
        except Exception:
            out_shape = healed_children[0] if healed_children else imported
    else:
        out_shape = healed_children[0] if healed_children else imported

    try:
        export_step(out_shape, str(p_out))
    except Exception as e:
        return json.dumps({"status": "error", "stage": "export",
                           "error": f"{type(e).__name__}: {e}"})

    return json.dumps({
        "status": "ok",
        "input": str(p_in.resolve()),
        "output": str(p_out.resolve()),
        "input_size": p_in.stat().st_size,
        "output_size": p_out.stat().st_size,
        "shapes_processed": len(children),
        "shape_reports": per_shape_reports,
        "next": f'build123d_inspect_step("{p_out}") — confirm warnings cleared',
    }, indent=2, ensure_ascii=False)


# ============================================================
# Entry point
# ============================================================



register_status_tool(
    mcp,
    server_name='mcp-server-build123d',
    description='build123d STEP authoring (CAD-as-code) + Cubit interop',
    subpackage='radia_mcp.build123d',
    related_servers=["cubit", "interop"],
    optional_deps=["build123d"],
)


def main():
    if "--selftest" in sys.argv:
        from radia_mcp.common.utf8_stdout import use_utf8_stdout
        use_utf8_stdout()
        print("build123d MCP server self-test:")

        # Test knowledge base topics (always available — pure text)
        topics = [
            "overview", "primitives_3d", "primitives_2d", "operations",
            "curves", "export_import", "topology", "cae_guidelines",
            "examples", "coil_modeling",
        ]
        for t in topics:
            result = build123d_usage(t)
            print(f"  build123d_usage('{t}'): {len(result)} chars")
            assert len(result) > 100, f"Topic '{t}' too short"

        # Tests that require the build123d library itself (optional
        # extra). Skip gracefully in minimal installs — the MCP server
        # CLI is still considered healthy without build123d present;
        # the execute/generate tools report a structured "not installed"
        # error at runtime when the user actually calls them.
        try:
            import build123d  # noqa: F401
            _have_b3d = True
        except ImportError:
            _have_b3d = False

        if _have_b3d:
            import json as _json

            # Test execute (use the sync impl directly — selftest runs
            # outside the FastMCP event loop and cannot `await` the tool)
            result = _execute_build123d_sync("box = Box(10, 20, 30)")
            print(f"  execute_build123d('Box(10,20,30)'): {len(result)} chars")
            data = _json.loads(result)
            assert data["status"] == "ok", f"Execute failed: {data}"
            assert data["is_valid"], "Box should be valid"
            assert abs(data["volume"] - 6000.0) < 1, (
                f"Volume wrong: {data['volume']}")

            # Test generate_helix_coil
            result = generate_helix_coil(radius=30, pitch=8, n_turns=2,
                                         w_start=3, h_start=3,
                                         w_end=2, h_end=2,
                                         sections_per_turn=8)
            data = _json.loads(result)
            print(f"  generate_helix_coil: volume={data.get('volume')}, "
                  f"path_pts={len(data.get('path_points', []))}")
            assert data["status"] == "ok", (
                f"generate_helix_coil failed: {data}")
            assert data["is_valid"], "Helix coil should be valid"
            assert data["volume"] > 0, "Volume should be positive"
            assert len(data["path_points"]) > 10, "Path should have points"
        else:
            print("  build123d library not installed — "
                  "skipping execute/generate tests (install "
                  "`radia-mcp[build123d]` for full coverage)")

        # Test prompt (no build123d dep)
        prompt = new_cae_geometry("ih_coil")
        print(f"  new_cae_geometry('ih_coil'): {len(prompt)} chars")
        assert "workpiece" in prompt

        print("  PASSED")
        return

    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
