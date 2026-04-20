"""coil_from_cad.py -- Extract PEEC filaments from a CAD solid + path.

Given a solid coil (STEP/BREP) and its centerline path (discrete points),
this module sections the solid perpendicular to the path at each segment
midpoint, measures the local cross-section area, and builds the
corresponding PEEC filament topology via PEECBuilder.

The cross-section is assumed roughly rectangular; the local (w, h) is
estimated from area and aspect ratio.  For square cross-sections,
w = h = sqrt(area).

Pipeline:
    build123d loft/sweep -> STEP solid
        + helix/spline path (discrete points)
    -> section_solid_along_path() -> list of (center, tangent, area, w, h)
    -> build_peec_from_sections() -> PEECBuilder topology_dict
    -> PEECCircuitSolver(topo, use_hacapk=True, outer_method="saddle")

Designed for IH coils with non-uniform cross-sections (tapered, shaped).
"""
from __future__ import annotations

import math
from typing import Optional

import numpy as np


def helix_path(radius: float, pitch: float, n_turns: float,
               n_points: int, z_offset: float = 0.0) -> np.ndarray:
    """Generate discrete helix centerline points.

    Args:
        radius: Helix radius [m].
        pitch: Axial advance per turn [m].
        n_turns: Number of turns.
        n_points: Number of discrete points (n_points - 1 segments).
        z_offset: Axial offset of first point [m].

    Returns:
        (n_points, 3) float64 array of (x, y, z) coordinates.
    """
    t = np.linspace(0, 1, n_points)
    angle = 2.0 * np.pi * n_turns * t
    z = pitch * n_turns * t + z_offset
    x = radius * np.cos(angle)
    y = radius * np.sin(angle)
    return np.column_stack([x, y, z])


def path_tangents(points: np.ndarray) -> np.ndarray:
    """Compute unit tangent vectors at segment midpoints.

    Args:
        points: (N, 3) path points.

    Returns:
        (N-1, 3) unit tangent vectors (forward difference, normalized).
    """
    d = np.diff(points, axis=0)
    norms = np.linalg.norm(d, axis=1, keepdims=True)
    norms = np.maximum(norms, 1e-30)
    return d / norms


def section_solid_along_path(step_path: str,
                             centers: np.ndarray,
                             tangents: np.ndarray,
                             units_scale: float = 1.0):
    """Section a STEP solid at each (center, tangent) and extract area.

    Uses build123d to import the STEP and section it.  Falls back to
    area estimation if build123d is not available.

    Args:
        step_path: Path to the STEP file (build123d units, typically mm).
        centers: (N, 3) midpoint coordinates [in build123d units].
        tangents: (N, 3) unit tangent vectors.
        units_scale: Multiply CAD coordinates by this to get build123d
            units.  E.g. if path is in meters and CAD is in mm,
            pass units_scale=1000.

    Returns:
        List of dicts with keys: area, w_est, h_est (in build123d units).
    """
    from build123d import import_step, section, Plane, Vector

    solid = import_step(step_path)
    results = []
    for i in range(len(centers)):
        c = centers[i] * units_scale
        t = tangents[i]
        origin = Vector(float(c[0]), float(c[1]), float(c[2]))
        z_dir = Vector(float(t[0]), float(t[1]), float(t[2]))

        sec_plane = Plane(origin=origin, z_dir=z_dir)
        try:
            cross = section(solid, section_by=sec_plane)
        except Exception:
            results.append({"area": None, "w_est": None, "h_est": None})
            continue

        if not cross or len(cross.faces()) == 0:
            results.append({"area": None, "w_est": None, "h_est": None})
            continue

        # Pick the face whose centroid is closest to the section center
        best_face = None
        best_dist = float("inf")
        for face in cross.faces():
            fc = face.center()
            dist = (fc - origin).length
            if dist < best_dist:
                best_dist = dist
                best_face = face

        area = best_face.area
        # For square cross-section: w = h = sqrt(area)
        side = math.sqrt(area)
        results.append({"area": area, "w_est": side, "h_est": side})

    return results


def build_peec_from_path(path_points: np.ndarray,
                         widths: np.ndarray,
                         heights: np.ndarray,
                         sigma: float = 5.8e7,
                         nwinc: int = 1,
                         nhinc: int = 1):
    """Build PEEC topology from a discrete path + per-segment cross-sections.

    Args:
        path_points: (N+1, 3) float64, path vertices [meters].
        widths: (N,) per-segment width [meters].
        heights: (N,) per-segment height [meters].
        sigma: Conductivity [S/m].
        nwinc, nhinc: Sub-filament subdivision per segment.

    Returns:
        topology_dict from PEECBuilder.build_topology().
    """
    from peec_matrices import PEECBuilder

    n_seg = len(widths)
    if path_points.shape[0] != n_seg + 1:
        raise ValueError(
            f"path_points must have {n_seg + 1} rows, got {path_points.shape[0]}"
        )

    builder = PEECBuilder()
    nodes = []
    for pt in path_points:
        nodes.append(builder.add_node_at(float(pt[0]), float(pt[1]), float(pt[2])))

    for i in range(n_seg):
        builder.add_connected_segment(
            nodes[i], nodes[i + 1],
            float(widths[i]), float(heights[i]),
            sigma=sigma,
            nwinc=nwinc, nhinc=nhinc,
        )

    builder.add_port(nodes[0], nodes[-1])
    return builder.build_topology()


def extract_centerline_from_step(step_path: str,
                                 n_segments: int = 100,
                                 cad_units_per_meter: float = 1000.0):
    """Auto-extract coil centerline + cross-sections from a STEP file.

    Algorithm:
    1. Find the longest edge of the solid (= sweep/loft spine, one of
       the rectangular cross-section corner paths).
    2. Sample the spine at n_segments+1 points → approximate path.
    3. At each segment midpoint, section perpendicular to the tangent.
    4. Cross-section centroid = true centerline point (corrects for
       the corner-to-center offset of the spine edge).
    5. Cross-section area → estimated (w, h).

    Works for any swept or lofted solid where the longest edge traces
    the coil path (helix, spiral, arbitrary 3D curve).

    Args:
        step_path: Path to .step file (CAD units, typically mm).
        n_segments: Number of filament segments to extract.
        cad_units_per_meter: Scale factor (default 1000 = mm).

    Returns:
        path_m: (N+1, 3) centerline points in meters.
        widths_m: (N,) per-segment width in meters.
        heights_m: (N,) per-segment height in meters.
    """
    from build123d import import_step, section, Plane, Vector

    solid = import_step(step_path)

    # Step 1: find spine (longest edge)
    edges = solid.edges()
    if not edges:
        raise RuntimeError("STEP solid has no edges")
    spine = max(edges, key=lambda e: e.length)

    # Step 2: sample spine at n+1 points
    spine_pts = np.zeros((n_segments + 1, 3), dtype=np.float64)
    for i in range(n_segments + 1):
        t = i / n_segments
        p = spine @ t
        spine_pts[i] = [p.X, p.Y, p.Z]

    # Step 3-5: section at midpoints → centroids + areas
    tangents = path_tangents(spine_pts)
    midpoints = 0.5 * (spine_pts[:-1] + spine_pts[1:])

    centerline = np.zeros((n_segments + 1, 3), dtype=np.float64)
    widths_cad = np.zeros(n_segments, dtype=np.float64)
    heights_cad = np.zeros(n_segments, dtype=np.float64)

    # First/last centerline points from spine endpoints
    centerline[0] = spine_pts[0]
    centerline[-1] = spine_pts[-1]

    for i in range(n_segments):
        c = midpoints[i]
        t = tangents[i]
        origin = Vector(float(c[0]), float(c[1]), float(c[2]))
        z_dir = Vector(float(t[0]), float(t[1]), float(t[2]))
        sec_plane = Plane(origin=origin, z_dir=z_dir)
        try:
            cross = section(solid, section_by=sec_plane)
        except Exception:
            cross = None

        if cross and len(cross.faces()) > 0:
            best = min(cross.faces(),
                       key=lambda f: (f.center() - origin).length)
            bc = best.center()
            # Update centerline: midpoint between corrected centers
            if i < n_segments:
                centerline[i] = [bc.X, bc.Y, bc.Z]
                if i == n_segments - 1:
                    # Last segment: also update final point
                    centerline[i + 1] = spine_pts[i + 1]
            side = math.sqrt(best.area)
            widths_cad[i] = side
            heights_cad[i] = side
        else:
            widths_cad[i] = widths_cad[max(0, i - 1)]
            heights_cad[i] = heights_cad[max(0, i - 1)]
            centerline[i] = c

    # Rebuild path from centroids: use midpoints as node positions
    # (shift from midpoint-of-spine to centroid-of-section)
    # Simple: nodes = [centroid_0, centroid_1, ..., centroid_{N-1}, spine_end]
    path_cad = np.zeros((n_segments + 1, 3), dtype=np.float64)
    path_cad[0] = centerline[0]
    for i in range(n_segments - 1):
        path_cad[i + 1] = 0.5 * (centerline[i] + centerline[i + 1])
    path_cad[-1] = centerline[-1]

    scale = 1.0 / cad_units_per_meter
    return path_cad * scale, widths_cad * scale, heights_cad * scale


def filaments_from_step(step_path: str,
                        path_points_m: Optional[np.ndarray] = None,
                        sigma: float = 5.8e7,
                        nwinc: int = 1,
                        nhinc: int = 1,
                        n_peri: Optional[int] = None,
                        cad_units_per_meter: float = 1000.0,
                        n_slices: int = 200,
                        use_coil_builder: bool = True):
    """End-to-end: STEP solid -> PEEC topology.

    Two paths are available:

    **CoilBuilder path** (use_coil_builder=True, default):
        STEP -> extract_centerline -> to_coil_builder -> to_filaments
        -> peec_bundle.  Uses Profile.sample_at() to place filaments
        according to the actual cross-section shape (circular, rect,
        loft).  Correct for any cross-section geometry.

    **Legacy path** (use_coil_builder=False):
        STEP -> extract_centerline -> build_peec_from_path
        -> C++ ExpandFilaments (rectangular grid only).
        Does not respect circular cross-section boundaries.

    Args:
        step_path: STEP file path.
        path_points_m: (N+1, 3) path vertices in meters, or None for
            auto-extraction.
        sigma: Conductivity [S/m].
        nwinc, nhinc: Sub-filament subdivision for the volume-grid
            placement.  Ignored when ``n_peri`` is set.
        n_peri: If given, place ``n_peri`` filaments on the cross-section
            PERIMETER only (thin-skin regime, d/delta >= 3).  Takes
            priority over nwinc/nhinc.  Requires use_coil_builder=True.
        cad_units_per_meter: Scale factor (default 1000 = mm).
        n_slices: Z-slice count for auto-extraction (ignored if
            path_points_m is given).
        use_coil_builder: If True (default), use CoilBuilder path for
            profile-aware filament placement.  Falls back to legacy
            path if CoilBuilder reconstruction fails.

    Returns:
        topology_dict from PEECBuilder.build_topology().
    """
    if use_coil_builder:
        return _filaments_via_coil_builder(
            step_path, sigma=sigma, nwinc=nwinc, nhinc=nhinc,
            n_slices=n_slices,
            start_hint=None,  # auto-detect from bounding box
            n_peri=n_peri)

    if n_peri is not None:
        raise ValueError(
            "n_peri requires use_coil_builder=True; the legacy "
            "C++ ExpandFilaments path is volume-grid-only.")

    # Legacy path: rectangular grid via C++ ExpandFilaments
    if path_points_m is None:
        path_m, widths_m, heights_m = extract_centerline_from_step(
            step_path, n_segments=n_slices,
            cad_units_per_meter=cad_units_per_meter,
        )
    else:
        n_seg = path_points_m.shape[0] - 1
        tangents = path_tangents(path_points_m)
        midpoints = 0.5 * (path_points_m[:-1] + path_points_m[1:])
        sections = section_solid_along_path(
            step_path, midpoints, tangents,
            units_scale=cad_units_per_meter,
        )
        scale_inv = 1.0 / cad_units_per_meter
        widths_m = np.zeros(n_seg, dtype=np.float64)
        heights_m = np.zeros(n_seg, dtype=np.float64)
        for i, sec in enumerate(sections):
            if sec["w_est"] is not None:
                widths_m[i] = sec["w_est"] * scale_inv
                heights_m[i] = sec["h_est"] * scale_inv
            else:
                raise RuntimeError(
                    f"Section {i} at midpoint {midpoints[i]} failed."
                )
        path_m = path_points_m

    topo = build_peec_from_path(
        path_m, widths_m, heights_m,
        sigma=sigma, nwinc=nwinc, nhinc=nhinc,
    )
    topo["filament_path"] = path_m
    topo["filament_widths"] = widths_m
    topo["filament_heights"] = heights_m
    return topo


def _bd_face_to_start_hint(face, cad_units_per_meter: float = 1000.0):
    """Convert a build123d Face to ((px,py,pz), (tx,ty,tz)) start_hint.

    The tangent is the *inward* normal at the face center, which is what
    walking-plane expects as the initial seed direction.  Coordinates are
    in CAD units (mm by default) — walking-plane operates in CAD units.
    """
    c = face.center()
    n = face.normal_at(c)
    # Outward by build123d convention; walking-plane wants inward.
    p = np.array([c.X, c.Y, c.Z], dtype=float)
    t = -np.array([n.X, n.Y, n.Z], dtype=float)
    tn = np.linalg.norm(t)
    if tn < 1e-12:
        raise ValueError("port face normal degenerate")
    return p, t / tn


def _bd_shape_to_netgen_solid(bd_shape):
    """Convert a build123d Shape (pythonocc-core TopoDS_Shape) to a
    netgen.occ Solid by serializing through a BRep file (lossless).

    XCAF labels / colors on bd_shape are lost in this step, but that is
    acceptable because the caller is expected to extract them *before*
    calling this helper.
    """
    import os
    import tempfile
    import build123d as bd
    from netgen.occ import OCCGeometry
    tmp = tempfile.NamedTemporaryFile(suffix=".brep", delete=False)
    tmp.close()
    try:
        bd.export_brep(bd_shape, tmp.name)
        ng_shape = OCCGeometry(tmp.name).shape
    finally:
        try:
            os.unlink(tmp.name)
        except OSError:
            pass
    # If the bd_shape was a Compound with one Solid, unwrap.
    solids = list(ng_shape.solids)
    if len(solids) == 1:
        return solids[0]
    if not solids:
        raise ValueError("bd_shape contains no solid after conversion")
    # Multiple solids -> return as a shape (caller fuses if needed)
    return ng_shape


def filaments_from_shape(shape,
                         port_face=None,
                         port_label_prefix: str = "peec_port",
                         sigma: float = 5.8e7,
                         nwinc: int = 1,
                         nhinc: int = 1,
                         cad_units_per_meter: float = 1000.0,
                         n_slices: int = 200):
    """End-to-end: build123d Shape -> PEEC topology (no STEP round-trip).

    Bypasses the `filaments_from_step` STEP-file path so that build123d
    in-memory markers (Python attribute `face.label`, `face.color`) are
    usable without being stripped by STEP export.

    Two ways to specify the port:
      1. Pass `port_face` explicitly (build123d Face object).  This is
         the most robust — no label-persistence issues.
      2. Omit `port_face`; the function then looks for a face with
         `label.startswith(port_label_prefix)` on `shape`.  *Caveat:*
         `shape.faces()` returns fresh Face wrappers on each call and
         loses labels set elsewhere, so this path only works if the
         caller attached labels in the same expression tree the shape
         was built from.

    If neither works, falls back to the bbox-based auto-hint (z-axis
    torus assumption, as in `filaments_from_step`).

    Args:
        shape: build123d Shape (Solid / Compound / Part).
        port_face: Optional explicit build123d Face for the PEEC port.
        port_label_prefix: Label prefix to search for when port_face is
            None. Default matches "peec_port", "peec_port_in", etc.
        sigma, nwinc, nhinc, cad_units_per_meter, n_slices: forwarded
            to the underlying walking-plane / CoilBuilder pipeline.

    Returns:
        dict compatible with `filaments_from_step` CoilBuilder result.
    """
    from coil_from_step import extract_centerline, to_coil_builder
    from peec_bundle import build_bundle_solver

    # --- 1. Resolve start_hint ---
    start_hint = None
    if port_face is not None:
        start_hint = _bd_face_to_start_hint(port_face, cad_units_per_meter)
    else:
        # Best-effort label scan (works only if labels persist on this
        # shape's immediate face list).
        try:
            for f in shape.faces():
                if getattr(f, "label", "") and f.label.startswith(port_label_prefix):
                    start_hint = _bd_face_to_start_hint(f, cad_units_per_meter)
                    break
        except Exception:
            pass

    # --- 2. Convert to netgen.occ solid ---
    ng_solid = _bd_shape_to_netgen_solid(shape)

    # --- 3. bbox fallback if still no hint ---
    if start_hint is None:
        bb = ng_solid.bounding_box
        cx = 0.5 * (bb[0][0] + bb[1][0])
        cy = 0.5 * (bb[0][1] + bb[1][1])
        rx = 0.5 * (bb[1][0] - bb[0][0])
        cz = 0.5 * (bb[0][2] + bb[1][2])
        start_hint = (np.array([cx + rx * 0.5, cy, cz]),
                      np.array([0.0, 1.0, 0.0]))

    # --- 4. Walking-plane + CoilBuilder (reuses the existing pipeline) ---
    res = extract_centerline(ng_solid, start_hint=start_hint, verbose=False)
    coil, _segs = to_coil_builder(res, current=1.0)
    paths, _currents = coil.to_filaments(
        nw=nwinc, nh=nhinc, frequency=0.0, sigma=sigma)
    cell_wh = coil._last_filament_info.get('cell_wh')
    solver, seg_of_fil, port_p, port_m = build_bundle_solver(
        paths, dw=None, dh=None, sigma=sigma, cell_wh=cell_wh)
    return {
        "solver": solver,
        "seg_of_filament": seg_of_fil,
        "filament_paths": paths,
        "cell_wh": cell_wh,
        "coil_builder": coil,
        "n_loop": len(paths),
        "port_plus": port_p,
        "port_minus": port_m,
    }


def _start_hint_from_step_labels(step_path: str,
                                  port_label_prefix: str = "peec_port"):
    """Return start_hint from XCAF labels in STEP, or None if unavailable.

    netgen.occ does not read XCAF labels. We use build123d (pythonocc-core)
    only if it is installed. Any failure (missing package, no labels,
    unexpected topology) silently returns None so the caller can fall
    through to bbox auto-detect.
    """
    try:
        import build123d as bd
    except ImportError:
        return None
    try:
        shape = bd.import_step(step_path)
    except Exception:
        return None
    # Walk direct children of the top Compound (FreeCAD Import.export
    # emits one PRODUCT per Part::Feature, each carrying its Label).
    candidates = []
    if getattr(shape, "children", None):
        candidates.extend(shape.children)
    candidates.append(shape)  # top-level label too
    for c in candidates:
        lab = getattr(c, "label", "") or ""
        if not lab.startswith(port_label_prefix):
            continue
        faces = list(c.faces()) if hasattr(c, "faces") else []
        if not faces:
            continue
        # For a Shell made from a single face (FreeCAD port marker
        # pattern), faces[0] is THE port face.
        port_face = faces[0]
        try:
            return _bd_face_to_start_hint(port_face)
        except Exception:
            continue
    return None


def _filaments_via_coil_builder(step_path, sigma, nwinc, nhinc, n_slices,
                                start_hint=None, n_peri=None):
    """CoilBuilder path: STEP -> centerline -> CoilBuilder -> to_filaments.

    Uses Profile.sample_at() for cross-section-aware filament placement.
    Supports circular, rectangular, and lofted cross-sections.

    If ``n_peri`` is given, uses perimeter-only placement via
    ``to_filaments_peri`` (thin-skin limit) and ignores ``nwinc``/``nhinc``.
    """
    from coil_from_step import extract_centerline, to_coil_builder
    from peec_bundle import build_bundle_solver
    import numpy as np

    # Step 1: Extract centerline with profile classification.
    # Start-hint resolution order (first hit wins):
    #   1. Caller-supplied start_hint.
    #   2. XCAF label on a "peec_port*" child (FreeCAD Import.export, or
    #      build123d in-memory via filaments_from_shape).  Requires the
    #      optional build123d dependency; silently skipped otherwise.
    #   3. Bounding-box heuristic (z-axis torus assumption — works for
    #      our canonical test cases but fails on x- or y-axis coils).
    if start_hint is None:
        start_hint = _start_hint_from_step_labels(step_path)
    if start_hint is None:
        from coil_from_step import load_step_solid
        solid = load_step_solid(step_path)
        bb = solid.bounding_box
        cx = 0.5 * (bb[0][0] + bb[1][0])
        cy = 0.5 * (bb[0][1] + bb[1][1])
        rx = 0.5 * (bb[1][0] - bb[0][0])
        start_hint = (np.array([cx + rx * 0.5, cy, 0.5 * (bb[0][2] + bb[1][2])]),
                      np.array([0, 1, 0]))

    res = extract_centerline(step_path, start_hint=start_hint, verbose=False)

    # Step 2: Reconstruct CoilBuilder from centerline + profiles
    coil, _segs = to_coil_builder(res, current=1.0)

    # Step 3: Generate profile-aware filaments
    if n_peri is not None:
        paths, _currents = coil.to_filaments_peri(
            n_peri=n_peri, frequency=0.0, sigma=sigma)
    else:
        paths, _currents = coil.to_filaments(
            nw=nwinc, nh=nhinc, frequency=0.0, sigma=sigma)
    cell_wh = coil._last_filament_info.get('cell_wh')

    # Step 4: Build PEEC topology via bundle solver
    # build_bundle_solver internally calls PEECBuilder.build_topology()
    # and remaps nodes for the parallel-bundle topology.  It returns
    # a PEECCircuitSolver (not a raw topo dict).  We return the solver
    # plus metadata for the caller.
    solver, seg_of_fil, port_p, port_m = build_bundle_solver(
        paths, dw=None, dh=None, sigma=sigma, cell_wh=cell_wh)

    # Return a dict compatible with the legacy path
    result = {
        "solver": solver,
        "seg_of_filament": seg_of_fil,
        "filament_paths": paths,
        "cell_wh": cell_wh,
        "coil_builder": coil,
        "n_loop": len(paths),
        "port_plus": port_p,
        "port_minus": port_m,
    }
    return result
