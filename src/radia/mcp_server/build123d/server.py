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
def preview_shape(script: str, name: str = "preview") -> str:
    """
    Run a build123d script and show the resulting Shape in the OCP CAD
    Viewer (VSCode panel).

    Lab policy (2026-04-19): build123d previews go **direct** —
    `ocp_vscode.show(part)` on the in-memory Part. No STEP roundtrip.
    This is mandatory; the Cubit path uses STEP because Cubit is
    out-of-process, but build123d has no such constraint.

    Prerequisite: the "OCP CAD Viewer" extension must be running in
    VSCode (default port 3939). If it's not, `show()` raises and this
    tool returns an error JSON — fix by opening the viewer panel.

    Args:
        script: Python code using build123d. Must be self-contained.
            The last Shape-like object assigned in the namespace is
            sent to the viewer.
        name: Display name for the shape in the viewer sidebar.

    Returns:
        JSON with preview status, shape summary (volume / faces /
        bbox), and any error traceback.
    """
    import sys as _sys
    import traceback as _tb
    from io import StringIO

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

    from build123d import Shape

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
            "error": "Script executed but no build123d Shape was found "
                     "in the namespace.",
            "stdout": stdout_text,
        }, indent=2)

    # Direct Part -> OCP viewer. No STEP.
    try:
        from ocp_vscode import show
    except ImportError as e:
        return json.dumps({
            "status": "error",
            "stage": "ocp_import",
            "error": f"ocp_vscode not installed: {e}. "
                     f"Install via: pip install ocp-vscode",
        }, indent=2)

    try:
        show(target, names=[name])
    except Exception:
        return json.dumps({
            "status": "error",
            "stage": "show",
            "error": _tb.format_exc(),
            "hint": "Is the OCP CAD Viewer panel open in VSCode? "
                    "(Command Palette -> OCP CAD Viewer: Open)",
        }, indent=2)

    info = {
        "status": "ok",
        "stage": "shown",
        "viewer": "ocp_vscode",
        "name": name,
        "variable": target_name,
        "type": type(target).__name__,
        "is_valid": target.is_valid,
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
            "min": [round(bb.min.X, 4), round(bb.min.Y, 4),
                    round(bb.min.Z, 4)],
            "max": [round(bb.max.X, 4), round(bb.max.Y, 4),
                    round(bb.max.Z, 4)],
        }
    except Exception:
        pass
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
            "curves", "export_import", "topology", "cae_guidelines",
            "examples", "coil_modeling",
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

        # Test generate_helix_coil
        result = generate_helix_coil(radius=30, pitch=8, n_turns=2,
                                     w_start=3, h_start=3, w_end=2, h_end=2,
                                     sections_per_turn=8)
        data = _json.loads(result)
        print(f"  generate_helix_coil: volume={data.get('volume')}, "
              f"path_pts={len(data.get('path_points', []))}")
        assert data["status"] == "ok", f"generate_helix_coil failed: {data}"
        assert data["is_valid"], "Helix coil should be valid"
        assert data["volume"] > 0, "Volume should be positive"
        assert len(data["path_points"]) > 10, "Path should have points"

        # Test prompt
        prompt = new_cae_geometry("ih_coil")
        print(f"  new_cae_geometry('ih_coil'): {len(prompt)} chars")
        assert "workpiece" in prompt

        print("  PASSED")
        return

    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
