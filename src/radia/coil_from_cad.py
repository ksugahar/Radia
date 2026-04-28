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


def _collect_loft_cross_sections(solid,
                                  min_count: int = 5,
                                  area_variation_tol: float = 0.5):
    """Return planar cross-section faces of a loft-of-profiles solid.

    A loft solid (multi-turn coil, generic spline coil) has N planar
    cross-section end-cap faces (circles / rects / custom profiles) +
    N-1 spline lateral surfaces.  The cross-section faces share very
    similar area (same profile).  Returns the filtered planar faces,
    or [] if the solid does NOT look like a loft-of-profiles (too few
    planar faces, or planar face area varies wildly).
    """
    try:
        from build123d import GeomType
        planar = [f for f in solid.faces() if f.geom_type == GeomType.PLANE]
    except Exception:
        return []
    if len(planar) < min_count:
        return []
    areas = np.array([f.area for f in planar], dtype=np.float64)
    median_area = float(np.median(areas))
    if median_area <= 0:
        return []
    keep = np.abs(areas - median_area) / median_area < area_variation_tol
    cross = [f for f, k in zip(planar, keep) if k]
    if len(cross) < min_count:
        return []
    return cross


def _chain_centroids_nn(pts: np.ndarray) -> np.ndarray:
    """Order a set of 3D points along a smooth polyline via nearest-neighbor.

    Endpoint detection: interior points have 2 close neighbors at
    similar distance; endpoints have only 1.  Rank points by
    second-nearest / first-nearest distance ratio — endpoints have
    the largest ratio.

    Then greedy nearest-neighbor walk from one endpoint.  For tight
    spiral / pancake geometries where 2 adjacent turns can be nearly
    as close as consecutive cross-sections within a turn, add a
    tangent-continuity bias: once we have 2 points, prefer neighbors
    whose direction continues the current tangent (penalize backward
    / perpendicular jumps to a different turn).
    """
    n = len(pts)
    if n < 2:
        return pts.copy()
    # Pairwise distances (vectorised, O(N^2) — fine for N ≤ 1000 or so)
    diff = pts[:, None, :] - pts[None, :, :]
    D = np.linalg.norm(diff, axis=-1)
    D_self_masked = D.copy()
    np.fill_diagonal(D_self_masked, np.inf)
    sorted_d = np.sort(D_self_masked, axis=1)
    nearest = sorted_d[:, 0]
    second = sorted_d[:, 1]
    endpoint_score = second / np.maximum(nearest, 1e-12)
    start = int(np.argmax(endpoint_score))

    visited = [start]
    remaining = set(range(n)) - {start}
    while remaining:
        curr = visited[-1]
        if len(visited) >= 2:
            prev = visited[-2]
            tangent = pts[curr] - pts[prev]
            t_norm = np.linalg.norm(tangent)
            tangent = tangent / t_norm if t_norm > 1e-12 else np.zeros(3)
            # Score: distance + backward-motion penalty
            best_idx = None
            best_score = float("inf")
            for i in remaining:
                vec = pts[i] - pts[curr]
                d = np.linalg.norm(vec)
                if d < 1e-12:
                    continue
                # Project vec onto tangent; negative cos means backward
                cos_a = float(np.dot(vec, tangent) / d)
                # Score combines nearness + forward alignment
                # cos = 1 (forward) => score = d;  cos = 0 (perp) => d*2
                # cos = -1 (backward) => d*inf effectively
                penalty = max(1.0 - cos_a, 0.01)
                score = d * penalty
                if score < best_score:
                    best_score = score
                    best_idx = i
            nxt = best_idx
        else:
            nxt = min(remaining, key=lambda i: D[curr][i])
        visited.append(nxt)
        remaining.discard(nxt)
    return pts[visited]


def _centerline_from_cross_sections(solid,
                                      cad_units_per_meter: float = 1.0):
    """Loft-of-profiles centerline: chain of cross-section face centroids.

    For lofted coils (multi-turn, tight pancake, Cubit `create volume
    loft surface i j` chains), each input profile becomes a planar
    end-cap face in the STEP.  Their centroids define the centerline
    at each profile station.  This bypasses the "longest edge" spine
    heuristic, which picks a cross-section circle (not a spine) and
    gives wildly wrong cross-section area at section() on multi-turn
    lofts.

    Returns (path_m, widths_m, heights_m) with widths = heights =
    equivalent-square-side from mean cross-section area.
    """
    cross = _collect_loft_cross_sections(solid)
    if not cross:
        raise ValueError(
            "solid does not look like a loft-of-profiles "
            "(fewer than 5 consistent-area planar faces)")
    centroids_raw = np.array([[c.X, c.Y, c.Z] for c in
                              (f.center() for f in cross)], dtype=np.float64)
    # Dedupe near-duplicate centroids: Cubit lofts share their end-cap
    # circle between adjacent loft volumes, and the shared surface can
    # appear as two planar faces (one per owning volume) at nearly the
    # same centroid.  Merge any two faces whose centroids are closer
    # than 10% of the cross-section equivalent radius.
    mean_area = float(np.mean([f.area for f in cross]))
    eq_radius = math.sqrt(mean_area / math.pi)
    dedup_tol = 0.1 * eq_radius

    kept = []
    for c in centroids_raw:
        if not any(np.linalg.norm(c - k) < dedup_tol for k in kept):
            kept.append(c)
    centroids = np.array(kept, dtype=np.float64)
    if len(centroids) < 3:
        raise ValueError(
            f"cross-section centroid dedupe left {len(centroids)} points "
            f"(<3); solid may not be a loft-of-profiles")
    ordered = _chain_centroids_nn(centroids)
    side = math.sqrt(mean_area)

    n_seg = len(ordered) - 1
    widths_cad = np.full(n_seg, side, dtype=np.float64)
    heights_cad = np.full(n_seg, side, dtype=np.float64)

    scale = 1.0 / cad_units_per_meter
    return ordered * scale, widths_cad * scale, heights_cad * scale


def _centerline_from_open_spine(solid, n_segments: int,
                                 cad_units_per_meter: float):
    """Single-loop coil centerline: longest open edge (spine) sampling.

    For swept / bent coils (gapped torus, single-turn helix), the STEP
    solid has an open boundary edge that traces the coil spine.  We
    pick the longest OPEN (non-closed-loop) edge — closed circle
    edges are cross-section boundaries, not spines.

    Raises ValueError if no open spine edge exists (e.g. loft of
    circles: all edges are closed cross-section circles).  Caller
    should fall back to ``_centerline_from_cross_sections``.
    """
    from build123d import section, Plane, Vector
    edges = solid.edges()
    if not edges:
        raise RuntimeError("STEP solid has no edges")

    def _is_closed(e):
        try:
            return bool(getattr(e, "is_closed", False))
        except Exception:
            try:
                start, end = e.start_point(), e.end_point()
                return (start - end).length < 1e-9
            except Exception:
                return False

    open_edges = [e for e in edges if not _is_closed(e)]
    if not open_edges:
        raise ValueError("no open spine edges (all edges are closed loops)")
    spine = max(open_edges, key=lambda e: e.length)

    spine_pts = np.zeros((n_segments + 1, 3), dtype=np.float64)
    for i in range(n_segments + 1):
        t = i / n_segments
        p = spine @ t
        spine_pts[i] = [p.X, p.Y, p.Z]

    tangents = path_tangents(spine_pts)
    midpoints = 0.5 * (spine_pts[:-1] + spine_pts[1:])

    centerline = np.zeros((n_segments + 1, 3), dtype=np.float64)
    widths_cad = np.zeros(n_segments, dtype=np.float64)
    heights_cad = np.zeros(n_segments, dtype=np.float64)
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
            if i < n_segments:
                centerline[i] = [bc.X, bc.Y, bc.Z]
                if i == n_segments - 1:
                    centerline[i + 1] = spine_pts[i + 1]
            side = math.sqrt(best.area)
            widths_cad[i] = side
            heights_cad[i] = side
        else:
            widths_cad[i] = widths_cad[max(0, i - 1)]
            heights_cad[i] = heights_cad[max(0, i - 1)]
            centerline[i] = c

    path_cad = np.zeros((n_segments + 1, 3), dtype=np.float64)
    path_cad[0] = centerline[0]
    for i in range(n_segments - 1):
        path_cad[i + 1] = 0.5 * (centerline[i] + centerline[i + 1])
    path_cad[-1] = centerline[-1]

    scale = 1.0 / cad_units_per_meter
    return path_cad * scale, widths_cad * scale, heights_cad * scale


def _centerline_from_revolution_sweep(solid, n_segments: int,
                                        cad_units_per_meter: float):
    """Single-loop swept coil: analytical arc from revolution surfaces.

    Handles any coil built by sweeping a profile (circle / rect /
    polygon) around an axis by an arbitrary angle — typical Cubit
    ``sweep surface N axis ... angle A`` workflow.  The lateral
    surfaces are revolution-type:

    * **Circle profile**: ``GeomType.TORUS`` faces (MajorR + MinorR)
    * **Rect profile**  : ``GeomType.CYLINDER`` faces (one per rect
      edge at its own radius) + ``GeomType.PLANE`` top/bottom caps
    * **Polygon profile**: mix of CYLINDER + PLANE faces
    * **Generic profile**: ``GeomType.REVOLUTION`` spline faces

    Algorithm:
      1. Collect revolution surfaces (TORUS / CYLINDER / CONE /
         REVOLUTION) that share a common axis.
      2. Require at least one PLANE end-cap face (for cross-section
         area + spine-radius extraction).  Single-loop swept coils
         always have 2 caps (at sweep angle start and end).
      3. Sweep angle = union of U intervals (max − min) across all
         revolution surfaces.
      4. Spine radius = distance from axis to cap centroid.
      5. Cross-section area = cap face area → equivalent-square side
         for the filaments_from_polyline downstream.

    Raises ``ValueError`` if no revolution surfaces, no end caps, or
    axis disagreement between revolution surfaces (not a simple
    single-loop sweep).
    """
    from build123d import GeomType
    try:
        from OCP.BRepAdaptor import BRepAdaptor_Surface
        from OCP.GeomAbs import (GeomAbs_Torus, GeomAbs_Cylinder,
                                  GeomAbs_Cone, GeomAbs_SurfaceOfRevolution)
    except ImportError as exc:
        raise ValueError(
            f"OCP not available for revolution surface extraction: {exc}")

    axis_loc = None
    axis_dir = None
    u_intervals = []
    rev_types_seen = []
    for f in solid.faces():
        ad = BRepAdaptor_Surface(f.wrapped)
        t = ad.GetType()
        if t == GeomAbs_Torus:
            surf = ad.Torus()
        elif t == GeomAbs_Cylinder:
            surf = ad.Cylinder()
        elif t == GeomAbs_Cone:
            surf = ad.Cone()
        elif t == GeomAbs_SurfaceOfRevolution:
            surf = ad  # SurfaceOfRevolution exposes Axis() directly via
            # the adaptor in OCP; fall back gracefully
        else:
            continue

        try:
            ax = surf.Axis()
        except AttributeError:
            continue
        loc = ax.Location()
        dr = ax.Direction()
        c = np.array([loc.X(), loc.Y(), loc.Z()], dtype=np.float64)
        a = np.array([dr.X(), dr.Y(), dr.Z()], dtype=np.float64)

        if axis_loc is None:
            axis_loc = c
            axis_dir = a
        else:
            if (np.linalg.norm(c - axis_loc) > 1e-6 or
                    abs(abs(np.dot(a, axis_dir)) - 1.0) > 1e-6):
                continue  # skip surfaces on a different axis (imprint quirks)

        u_intervals.append((ad.FirstUParameter(), ad.LastUParameter()))
        rev_types_seen.append(t)

    if axis_loc is None or not u_intervals:
        raise ValueError("no revolution surfaces with a consistent axis")

    # Sweep angle = union of U intervals (max - min).
    u_min = min(a for a, _ in u_intervals)
    u_max = max(b for _, b in u_intervals)

    axis_dir = axis_dir / max(np.linalg.norm(axis_dir), 1e-30)

    # End-cap faces: planar faces whose normal is roughly parallel to
    # sweep U direction at the cap location.  Simpler proxy: planar
    # faces whose centroid is FARTHER from the axis than 0 and whose
    # area is much smaller than the lateral surfaces (cap area is the
    # cross-section area; lateral area is sweep-angle × circumference
    # × height, typically much larger).
    planar = [f for f in solid.faces() if f.geom_type == GeomType.PLANE]
    if not planar:
        raise ValueError("no PLANE end-cap faces")

    # Pick the smallest-area planar face as the cross-section cap.
    # For a rect-torus: end caps = rect (area small); top/bottom flat
    # "cap" of the lateral surface = much larger planar region.
    planar_sorted = sorted(planar, key=lambda f: f.area)
    cap = planar_sorted[0]
    cap_center = cap.center()
    cap_c = np.array([cap_center.X, cap_center.Y, cap_center.Z])

    # Project cap centroid to the axis-perpendicular plane at axis_loc
    offset = cap_c - axis_loc
    along_axis = np.dot(offset, axis_dir) * axis_dir
    radial = offset - along_axis
    R_spine = float(np.linalg.norm(radial))
    if R_spine < 1e-12:
        raise ValueError("cap centroid lies on the sweep axis — degenerate")

    # Axis-normal basis: u_hat points from axis to cap centroid
    u_hat = radial / R_spine
    v_hat = np.cross(axis_dir, u_hat)

    # Sweep angle U=0 corresponds to u_hat direction in our basis; shift
    # so sample at u_min lands near the cap centroid.
    us = np.linspace(u_min, u_max, n_segments + 1)
    path_raw = (axis_loc
                + R_spine * (np.cos(us)[:, None] * u_hat
                              + np.sin(us)[:, None] * v_hat)
                + along_axis)

    # Cross-section: equivalent-square side from cap area.
    side = math.sqrt(cap.area)
    widths_cad = np.full(n_segments, side, dtype=np.float64)
    heights_cad = np.full(n_segments, side, dtype=np.float64)

    scale = 1.0 / cad_units_per_meter
    return path_raw * scale, widths_cad * scale, heights_cad * scale


# Backward-compat alias for the torus-specific name used in tests.
_centerline_from_torus_sweep = _centerline_from_revolution_sweep


def extract_centerline_from_step(step_path: str,
                                 n_segments: int = 100,
                                 cad_units_per_meter: float = 1.0):
    """Auto-extract coil centerline + cross-sections from a STEP file.

    Dispatches on solid topology:

    * **Loft of profiles** (multi-turn coil, tight pancake spiral):
      solid has >= 5 planar end-cap faces of consistent area →
      chain their centroids via nearest-neighbor + tangent continuity.
      No section() calls, robust for tight geometries where adjacent
      turns are close enough to confuse the spine method.

    * **Single-loop swept coil** (gapped torus, simple bend):
      solid has an open boundary edge that traces the coil spine →
      sample the longest OPEN edge, section at midpoints, centroid
      of each section gives the true centerline.  Open-edge filter
      excludes closed circle edges (cross-section boundaries).

    Args:
        step_path: Path to .step file.  Coordinates MUST be in metres
            unless ``cad_units_per_meter`` is explicitly overridden.
        n_segments: Number of filament segments for the open-spine
            path.  Ignored for the loft path (uses N planar faces).
        cad_units_per_meter: Scale factor.  Default 1.0 = STEP
            coordinates are in metres (CLAUDE.md "Unit System Policy:
            Radia always uses meters").  Pass 1000.0 if the STEP is
            in millimetres.  No auto-detection -- caller is
            responsible for knowing the input unit (Fail Fast Loud).

    Returns:
        path_m: (N+1, 3) centerline points in meters.
        widths_m: (N,) per-segment width in meters.
        heights_m: (N,) per-segment height in meters.
    """
    from build123d import import_step

    solid = import_step(step_path)

    # Try loft-of-profiles first: if the solid clearly has many
    # consistent-area cross-sections, chain their centroids.
    # (This is the robust path for Kubota's 3turncoil.stp class.)
    cross_faces = _collect_loft_cross_sections(solid)
    if cross_faces:
        return _centerline_from_cross_sections(solid, cad_units_per_meter)

    # Torus-shaped single-loop coil (gapped torus, full torus): extract
    # major / minor radius + axis + sweep angle analytically from the
    # TORUS face parameters.  Handles the single-loop case that the
    # open-spine fallback gets wrong (picks a cross-section arc as
    # "longest open edge", path length comes out half the real arc).
    try:
        return _centerline_from_torus_sweep(solid, n_segments, cad_units_per_meter)
    except ValueError:
        pass

    # Fall back to open-spine method for generic single-loop swept solids
    # (non-torus: helical bend, spline coil, rectangular cross-section).
    return _centerline_from_open_spine(solid, n_segments, cad_units_per_meter)


def filaments_from_step(step_path: str,
                        path_points_m: Optional[np.ndarray] = None,
                        sigma: float = 5.8e7,
                        nwinc: int = 1,
                        nhinc: int = 1,
                        n_peri: Optional[int] = None,
                        cad_units_per_meter: float = 1.0,
                        n_slices: int = 200,
                        use_coil_builder: bool = True):
    """End-to-end: STEP solid -> PEEC topology.

    Two paths are available:

    **CoilBuilder path** (use_coil_builder=True, n_peri=None, default
    for nwinc/nhinc volume-grid mode):
        STEP -> extract_centerline -> to_coil_builder -> to_filaments
        -> peec_bundle.  Uses Profile.sample_at() to place filaments
        according to the actual cross-section shape (circular, rect,
        loft).  Correct for any cross-section geometry.

    **Longest-edge path** (n_peri given, perimeter-only placement):
        STEP -> extract_centerline_from_step (spine = longest edge,
        global sampling) -> filaments_from_polyline (parallel-transport
        frame, equivalent-circle profile from mean area).
        This is the PRIMARY path for ``n_peri``.  Chosen over the walker
        because the walker hangs or natively crashes on multi-turn loft
        STEPs (Kubota's 3turncoil.stp: walker hangs netgen.occ > 5 min;
        on 100号機 the subprocess exits with an unhandleable native
        error code).  Longest-edge is robust because it samples the
        whole spine in one pass rather than walking step-by-step.

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
        cad_units_per_meter: Scale factor.  Default 1.0 = STEP coordinates
            are in metres (CLAUDE.md "Unit System Policy: Radia always uses
            meters").  Pass 1000.0 if the STEP is in millimetres.
        n_slices: Z-slice count for auto-extraction (ignored if
            path_points_m is given).
        use_coil_builder: If True (default), use CoilBuilder path for
            profile-aware filament placement.  Falls back to legacy
            path if CoilBuilder reconstruction fails.

    Returns:
        topology_dict from PEECBuilder.build_topology().
    """
    if use_coil_builder and n_peri is None:
        # Volume-grid placement (nwinc/nhinc) needs the profile-aware
        # walker path.  Only the walker knows per-segment (w, h).
        return _filaments_via_coil_builder(
            step_path, sigma=sigma, nwinc=nwinc, nhinc=nhinc,
            n_slices=n_slices,
            start_hint=None,  # auto-detect from bounding box
            n_peri=None)

    if n_peri is not None and use_coil_builder:
        # Perimeter-only placement: use the longest-edge extractor as the
        # PRIMARY path.  The walker can hang or crash on multi-turn
        # loft STEPs (observed: Kubota's 3turncoil.stp hangs netgen.occ
        # indefinitely; on 100号機 the subprocess exits with a native
        # error code that Python cannot catch via try/except).
        # Longest-edge samples the spine globally and works for both
        # simple loops and multi-turn lofts.  For n_peri the equivalent
        # circle from mean cross-section area is accurate enough in the
        # thin-skin regime this mode targets.
        import numpy as np
        # filaments_from_polyline lives at module bottom (was in
        # coil_from_jou.py until 4.13.0 .jou-path retirement).
        path_m, w_m, h_m = extract_centerline_from_step(
            step_path, n_segments=n_slices,
            cad_units_per_meter=cad_units_per_meter)
        mean_area = float(np.mean(w_m * h_m))
        r_m = float(np.sqrt(mean_area / np.pi))
        topo = filaments_from_polyline(
            path_m, r_m,
            sigma=sigma, n_peri=n_peri,
            source_tag="step_longest_edge")
        return topo

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


def _bd_face_to_start_hint(face, cad_units_per_meter: float = 1.0):
    """Convert a build123d Face to ((px,py,pz), (tx,ty,tz)) start_hint.

    The tangent is the *inward* normal at the face center, which is what
    walking-plane expects as the initial seed direction.  Coordinates are
    returned in raw CAD units -- walking-plane operates in CAD units, the
    caller scales by 1/cad_units_per_meter at the boundary.
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
                         sigma: float = 5.8e7,
                         nwinc: int = 1,
                         nhinc: int = 1,
                         cad_units_per_meter: float = 1.0,
                         n_slices: int = 200):
    """End-to-end: build123d Shape -> PEEC topology (no STEP round-trip).

    Canonical entry point for build123d workflows.  The in-memory path
    bypasses STEP export entirely — important because
    `build123d.export_step` strips both face labels and child-solid
    labels.  Pass the port face directly instead of relying on
    through-STEP markers.

    Port specification:
      - `port_face` (build123d Face, optional): explicit seed.  Center
        + inward normal are used to construct the walking-plane start
        hint.  Most reliable option; recommended for all build123d
        callers.
      - If omitted, the function falls back to the bbox-based auto-hint
        (z-axis torus assumption; works for the classical single-axis
        helical coils but not for off-axis geometry).

    STEP-file workflows should use `filaments_from_step(step_path)`
    instead, which picks up XCAF labels written by external CAD tools
    (e.g. FreeCAD's Import.export).

    Args:
        shape: build123d Shape (Solid / Compound / Part).
        port_face: Optional explicit build123d Face for the PEEC port.
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


def export_step_with_labels(items, path: str, schema: str = "AP214IS"):
    """Write a STEP with per-shape XCAF labels preserved.

    build123d's own `export_step` drops Compound-child labels in the
    current release (only the top-level Compound label survives; see
    2026-04-20 tests).  This helper bypasses the build123d writer and
    uses pythonocc-core (OCP) XCAF directly, producing STEP files with
    one `PRODUCT('<label>', ...)` entity per labeled sub-shape — the
    same form FreeCAD emits via Import.export, and the form the rest of
    Radia's label-detection pipeline expects.

    Note: build123d `Compound([a, b])` returns `.children == ()` (the
    `children` attribute is for explicit Part/Assembly trees), so we
    cannot walk a Compound to recover its labeled members.  Callers must
    pass the list of labeled shapes directly.

    Args:
        items: iterable of build123d shapes (Solid / Shell / Part), each
            with its `label` attribute set.  A single Shape is also
            accepted and written as one XCAF entity.
        path: Output STEP file path.
        schema: STEP schema, default "AP214IS" (same as FreeCAD default).

    Side-effect: writes `path`. Raises `RuntimeError` if the OCP writer
    does not return success.
    """
    try:
        from OCP.TDocStd import TDocStd_Document
        from OCP.TCollection import TCollection_ExtendedString
        from OCP.XCAFApp import XCAFApp_Application
        from OCP.XCAFDoc import XCAFDoc_DocumentTool
        from OCP.TDataStd import TDataStd_Name
        from OCP.STEPCAFControl import STEPCAFControl_Writer
        from OCP.STEPControl import STEPControl_AsIs
        from OCP.Interface import Interface_Static
        from OCP.IFSelect import IFSelect_RetDone
    except ImportError as exc:
        raise RuntimeError(
            "export_step_with_labels requires build123d (pythonocc-core)"
        ) from exc

    doc = TDocStd_Document(TCollection_ExtendedString("radia-coil"))
    XCAFApp_Application.GetApplication_s().InitDocument(doc)
    tool = XCAFDoc_DocumentTool.ShapeTool_s(doc.Main())

    # Normalize: if a single shape was passed, wrap in a list.
    if hasattr(items, "wrapped"):
        items_iter = [items]
    else:
        items_iter = list(items)

    for item in items_iter:
        wrapped = getattr(item, "wrapped", None)
        if wrapped is None:
            continue
        lab = tool.AddShape(wrapped, False, True)
        name = getattr(item, "label", "") or ""
        if name:
            TDataStd_Name.Set_s(lab, TCollection_ExtendedString(name))

    Interface_Static.SetCVal_s("write.step.schema", schema)
    writer = STEPCAFControl_Writer()
    writer.Transfer(doc, STEPControl_AsIs)
    res = writer.Write(path)
    if res != IFSelect_RetDone:
        raise RuntimeError(f"STEP write failed: {res}")


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


# ----------------------------------------------------------------------
# Polyline -> PEEC topology (used by the longest-edge STEP path)
#
# Until 4.13.0 these helpers lived in coil_from_jou.py to support the
# legacy `.jou` explicit-centerline input.  That input path was retired
# (CLAUDE.md "Radia always uses meters" + No-Fallbacks: a single STEP
# input is the canonical PEEC source, and cross-section centroids are
# auto-extracted from B-Rep).  The polyline helpers themselves are
# still useful internally -- the longest-edge path samples the spine
# and feeds it through `filaments_from_polyline`.  Moving them here
# lets coil_from_jou.py be deleted entirely.
# ----------------------------------------------------------------------

def _parallel_transport_frame(pts: np.ndarray):
    """Compute (tangent, u_hat, v_hat) per vertex of a 3D polyline.

    Uses parallel transport: start with an arbitrary u perpendicular to
    the first tangent, then rotate u per segment by the bend angle using
    Rodrigues' formula.  Avoids the Frenet-Serret twist that appears
    when curvature vanishes (straight segments).

    Returns 3 arrays of shape (N, 3) where N = len(pts).
    """
    n = len(pts)
    if n < 2:
        raise ValueError(f"need at least 2 points, got {n}")

    seg_t = np.diff(pts, axis=0)
    seg_len = np.linalg.norm(seg_t, axis=1, keepdims=True)
    seg_len = np.maximum(seg_len, 1e-30)
    seg_t = seg_t / seg_len

    tangent = np.zeros((n, 3))
    tangent[0] = seg_t[0]
    tangent[-1] = seg_t[-1]
    if n > 2:
        tangent[1:-1] = seg_t[:-1] + seg_t[1:]
        mid_norm = np.linalg.norm(tangent[1:-1], axis=1, keepdims=True)
        mid_norm = np.maximum(mid_norm, 1e-30)
        tangent[1:-1] = tangent[1:-1] / mid_norm

    u_hat = np.zeros((n, 3))
    v_hat = np.zeros((n, 3))
    t0 = tangent[0]
    pick = int(np.argmin(np.abs(t0)))
    cand = np.zeros(3)
    cand[pick] = 1.0
    u0 = cand - np.dot(cand, t0) * t0
    u0_norm = np.linalg.norm(u0)
    if u0_norm < 1e-12:
        cand = np.array([0.0, 1.0, 0.0]) if pick == 0 else np.array([1.0, 0.0, 0.0])
        u0 = cand - np.dot(cand, t0) * t0
        u0_norm = np.linalg.norm(u0)
    u_hat[0] = u0 / u0_norm
    v_hat[0] = np.cross(tangent[0], u_hat[0])

    for i in range(1, n):
        t_prev = tangent[i - 1]
        t_curr = tangent[i]
        axis = np.cross(t_prev, t_curr)
        sin_a = np.linalg.norm(axis)
        cos_a = float(np.clip(np.dot(t_prev, t_curr), -1.0, 1.0))
        if sin_a < 1e-12:
            u_hat[i] = u_hat[i - 1]
        else:
            k = axis / sin_a
            u_prev = u_hat[i - 1]
            u_new = (u_prev * cos_a
                     + np.cross(k, u_prev) * sin_a
                     + k * np.dot(k, u_prev) * (1.0 - cos_a))
            u_hat[i] = u_new
        u_hat[i] -= np.dot(u_hat[i], t_curr) * t_curr
        u_hat[i] /= max(np.linalg.norm(u_hat[i]), 1e-30)
        v_hat[i] = np.cross(tangent[i], u_hat[i])

    return tangent, u_hat, v_hat


def filaments_from_polyline(pts_m: np.ndarray,
                             radius_m: float,
                             *,
                             sigma: float = 5.8e7,
                             n_peri: int = 16,
                             source_tag: str = "polyline"):
    """Build a PEEC topology from an explicit 3D polyline + circular profile.

    Internal helper.  Used by ``filaments_from_step`` longest-edge path
    after sampling the spine of a clean swept-loop STEP solid.  The
    n_peri filaments are placed at equal arc-length around the
    cross-section circle perimeter (thin-skin perimeter PEEC).

    Args:
        pts_m: (N, 3) float64 centerline points in METERS.
        radius_m: circular cross-section radius in METERS.
        sigma, n_peri: solver parameters.
        source_tag: free-form string copied into the result dict's
            ``source`` field (e.g. "step_longest_edge").

    Returns:
        topology_dict, same shape as ``filaments_from_step``.
    """
    pts = np.asarray(pts_m, dtype=float)
    n_pts = len(pts)
    if n_pts < 2:
        raise ValueError(f"need at least 2 centerline points, got {n_pts}")

    _, u_hat, v_hat = _parallel_transport_frame(pts)

    theta = 2.0 * np.pi * (np.arange(n_peri) + 0.5) / float(n_peri)
    u_offset = radius_m * np.cos(theta)
    v_offset = radius_m * np.sin(theta)

    A_cell = (math.pi * radius_m * radius_m) / float(n_peri)
    side = float(math.sqrt(A_cell))

    filament_paths = []
    cell_wh = []
    for k in range(n_peri):
        fil_pts = pts + u_offset[k] * u_hat + v_offset[k] * v_hat
        segs = [(tuple(fil_pts[i]), tuple(fil_pts[i + 1]))
                for i in range(n_pts - 1)]
        filament_paths.append(segs)
        cell_wh.append([(side, side)] * (n_pts - 1))

    import radia  # noqa: F401
    from peec_bundle import build_bundle_solver
    solver, seg_of_fil, port_p, port_m = build_bundle_solver(
        filament_paths, dw=None, dh=None, sigma=sigma, cell_wh=cell_wh)

    return {
        "solver": solver,
        "seg_of_filament": seg_of_fil,
        "filament_paths": filament_paths,
        "cell_wh": cell_wh,
        "n_loop": len(filament_paths),
        "port_plus": port_p,
        "port_minus": port_m,
        "source": source_tag,
        "n_path_pts": n_pts,
        "cross_section_radius_m": radius_m,
    }
