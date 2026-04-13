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
import traceback
from pathlib import Path
from io import StringIO

from mcp.server.fastmcp import FastMCP

from .build123d_knowledge import get_build123d_documentation

mcp = FastMCP("mcp-server-build123d")


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
            "primitives_3d"   - Box, Cylinder, Cone, Sphere, Torus, Wedge + boolean
            "primitives_2d"   - Circle, Rectangle, Polygon, etc. (sketch objects)
            "operations"      - extrude, revolve, sweep, loft, fillet, mirror, split
            "curves"          - Line, Polyline, Spline, Arc types, Helix, Bezier
            "export_import"   - STEP, BREP, STL, glTF, SVG + CAE pipeline examples
            "topology"        - faces/edges/vertices queries, labels, selectors
            "cae_guidelines"  - Clean geometry rules, quality checks, label conventions
            "examples"        - IH coil, E-core transformer, dipole yoke, parametric
            "all"             - Complete documentation
    """
    return get_build123d_documentation(topic)


@mcp.tool()
def execute_build123d(
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
    # Capture stdout
    old_stdout = sys.stdout
    sys.stdout = captured = StringIO()

    namespace = {}
    try:
        exec(  # noqa: S102
            "from build123d import *\n" + script,
            namespace,
        )
    except Exception:
        sys.stdout = old_stdout
        return json.dumps({
            "status": "error",
            "error": traceback.format_exc(),
            "stdout": captured.getvalue(),
        }, indent=2)
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


# ============================================================
# MCP Prompts
# ============================================================

@mcp.prompt()
def new_cae_geometry(geometry_type: str = "ih_coil") -> str:
    """Create a new CAE geometry with build123d."""
    base = (
        "Create a CAE-ready geometry using build123d.\n\n"
        "Use build123d_usage tool for API reference.\n\n"
        "Workflow:\n"
        "1. Create geometry with primitives + boolean (build123d_usage('primitives_3d'))\n"
        "2. Assign labels for material regions (build123d_usage('topology'))\n"
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
            "- E-shaped iron core (parametric: width, height, leg thickness)\n"
            "- Winding regions around center and/or outer legs\n"
            "- Label: 'iron_core', 'winding_primary', 'winding_secondary'\n"
            "- See build123d_usage('examples') for e_core function\n"
        )
    elif geo in ("dipole", "magnet", "accelerator"):
        base += (
            "Accelerator Dipole Magnet (quarter model):\n"
            "- Quarter-annular yoke (symmetry exploitation)\n"
            "- Bore region for beam pipe\n"
            "- Label: 'iron_yoke', 'bore', 'air'\n"
            "- See build123d_usage('examples') for yoke code\n"
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
# Entry point
# ============================================================

def main():
    if "--selftest" in sys.argv:
        print("build123d MCP server self-test:")

        # Test knowledge base topics
        topics = [
            "overview", "primitives_3d", "primitives_2d", "operations",
            "curves", "export_import", "topology", "cae_guidelines", "examples",
        ]
        for t in topics:
            result = build123d_usage(t)
            print(f"  build123d_usage('{t}'): {len(result)} chars")
            assert len(result) > 100, f"Topic '{t}' too short"

        # Test execute
        result = execute_build123d("box = Box(10, 20, 30)")
        print(f"  execute_build123d('Box(10,20,30)'): {len(result)} chars")
        import json as _json
        data = _json.loads(result)
        assert data["status"] == "ok", f"Execute failed: {data}"
        assert data["is_valid"], "Box should be valid"
        assert abs(data["volume"] - 6000.0) < 1, f"Volume wrong: {data['volume']}"

        # Test prompt
        prompt = new_cae_geometry("ih_coil")
        print(f"  new_cae_geometry('ih_coil'): {len(prompt)} chars")
        assert "workpiece" in prompt

        print("  PASSED")
        return

    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
