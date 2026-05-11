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
                                  cluster_jump_threshold: float = 1.5):
    """Return planar cross-section CAP faces of a loft-of-profiles solid.

    A loft solid has N planar cross-section end-caps + (for non-
    circular cross-sections) additional planar top/bottom faces from
    each loft segment.  For a 6 x 4 mm rect swept around a 30 mm arc
    with 9.7 mm segments, each cross-section cap has 24 mm^2 (4 LINE
    edges) while each top/bottom face has ~58 mm^2 (2 LINE + 2 BSPLINE
    edges).  A pure area-cluster filter would still mix them when the
    sizes overlap, so combine TWO filters:

      1. Smallest area cluster (cap is smallest planar face by area)
      2. Boundary made entirely of LINE / CIRCLE / ELLIPSE edges
         (cross-section caps are bounded by the original 2D profile
         curves -- straight lines or arcs.  Top / bottom strips of a
         swept solid acquire BSPLINE edges along the spine direction
         due to the curvature, so they fail this filter.)

    Returns [] if the solid does NOT look like a loft-of-profiles
    (too few planar faces of clean shape).
    """
    try:
        from build123d import GeomType
        planar_all = [f for f in solid.faces()
                      if f.geom_type == GeomType.PLANE]
    except Exception:
        return []
    if len(planar_all) < min_count:
        return []

    # Filter to faces whose outer wire has only "clean" cross-section
    # edges (no BSPLINE).  Cross-section caps were drawn as 2D shapes
    # (rect = 4 LINE, circle = 1 CIRCLE, polygon = N LINE, fillet rect
    # = LINE + ELLIPSE/CIRCLE arcs); they have NO BSPLINE on their
    # boundary.  Top/bottom strips ALWAYS pick up BSPLINE edges along
    # the swept direction unless the spine is straight.
    clean_types = {GeomType.LINE, GeomType.CIRCLE, GeomType.ELLIPSE}
    planar = []
    for f in planar_all:
        try:
            edges = f.outer_wire().edges()
        except Exception:
            continue
        if not edges:
            continue
        if all(e.geom_type in clean_types for e in edges):
            planar.append(f)

    if len(planar) < min_count:
        # No clean caps detected -- this solid is not a clean loft of
        # 2D profiles (or its cross-section uses freeform / spline
        # boundaries).  Caller falls through to other paths.
        return []

    areas = np.array([f.area for f in planar], dtype=np.float64)
    if not np.all(areas > 0):
        return []

    # Within the clean-bounded subset, keep the SMALLEST area cluster
    # (handles the very rare case where some clean-bounded faces are
    # still extras -- e.g. a flat tab welded onto the coil).
    sorted_idx = np.argsort(areas)
    sorted_areas = areas[sorted_idx]
    cluster_end = 1
    for i in range(1, len(sorted_areas)):
        if sorted_areas[i] / sorted_areas[i - 1] > cluster_jump_threshold:
            break
        cluster_end = i + 1
    if cluster_end < min_count:
        return []
    cross = [planar[i] for i in sorted_idx[:cluster_end]]
    return cross


def _find_lateral_surface(solid):
    """Pick the SINGLE dominant lateral surface from a coil solid.

    A clean swept / lofted coil with a single continuous lateral
    surface has exactly ONE lateral face that dominates the surface
    area.  This function returns it for direct (u, v) sampling.

    Returns ``None`` (caller falls through) when:
      - There are no lateral-type candidates
      - There are too many candidates (multi-fragment loft, united
        multi-turn pancake) -- a single face only covers part of
        the spine, sampling it would give a fragmented result
      - The largest candidate's area is not clearly dominant
        (gapped torus split into 4 TORUS quadrants by Cubit's
        webcut at the gap, etc.)

    The threshold "largest >= 80% of total" is a strong signal that
    we have a single-piece lateral surface.  Anything below falls
    through to the per-station-faces path or the legacy equivalent-
    circle path.
    """
    from build123d import GeomType
    candidates = []
    for f in solid.faces():
        gt = f.geom_type
        if gt in (GeomType.BSPLINE, GeomType.CYLINDER, GeomType.TORUS,
                   GeomType.REVOLUTION, GeomType.EXTRUSION):
            candidates.append((float(f.area), f))
    if not candidates:
        return None
    candidates.sort(key=lambda x: -x[0])
    largest_area = candidates[0][0]
    total_area = sum(a for a, _ in candidates)
    if total_area <= 0:
        return None
    dominance = largest_area / total_area
    # Multi-fragment threshold.  Single-face lateral: dominance ~ 1.0.
    # Anything below 0.8 indicates fragmentation -- fall through.
    if dominance < 0.8:
        return None
    # Also reject if there are too many candidates even if the largest
    # is technically dominant (defensive bound).
    if len(candidates) > 3:
        return None
    return candidates[0][1]


def _sample_lateral_surface_uv(face, n_stations: int, n_peri: int):
    """Sample a coil's lateral surface at a (u, v) parametric grid.

    OCC's parametric surface ``Geom_Surface.Value(u, v)`` evaluates
    the 3D position at any (u, v) within the surface's parameter
    range.  For a coil-style lateral surface, ONE axis is closed
    (the cross-section perimeter wraps around the spine) and the
    OTHER axis is open (the spine).  We auto-detect which is which
    via ``IsUClosed`` / ``IsVClosed``.

    Args:
        face: build123d Face whose underlying ``Geom_Surface`` is
            sampled (typically the largest BSPLINE / TORUS / CYLINDER
            / REVOLUTION face from ``_find_lateral_surface``).
        n_stations: number of samples along the spine direction.
        n_peri: number of samples around the perimeter direction.

    Returns:
        (n_stations, n_peri, 3) array of 3D points in CAD raw units.
        ``pts[i, k]`` is the position of filament k at station i.
    """
    from OCP.BRep import BRep_Tool
    from OCP.BRepTools import BRepTools

    surface = BRep_Tool.Surface_s(face.wrapped)
    u_min, u_max, v_min, v_max = BRepTools.UVBounds_s(face.wrapped)

    u_closed = bool(surface.IsUClosed())
    v_closed = bool(surface.IsVClosed())

    if u_closed and not v_closed:
        peri_lo, peri_hi = u_min, u_max
        spine_lo, spine_hi = v_min, v_max
        peri_axis = "u"
    elif v_closed and not u_closed:
        peri_lo, peri_hi = v_min, v_max
        spine_lo, spine_hi = u_min, u_max
        peri_axis = "v"
    else:
        # Closed both / neither: cannot identify perimeter unambiguously
        raise ValueError(
            f"lateral surface is unsuitable for UV sampling "
            f"(IsUClosed={u_closed}, IsVClosed={v_closed}); coil "
            f"shape may not be a clean swept loft")

    pts = np.zeros((n_stations, n_peri, 3), dtype=np.float64)
    spine_span = spine_hi - spine_lo
    peri_span = peri_hi - peri_lo
    for i in range(n_stations):
        s = i / (n_stations - 1) if n_stations > 1 else 0.0
        spine_p = spine_lo + spine_span * s
        for k in range(n_peri):
            # Cell-centred perimeter samples (avoid hitting the
            # closed-loop seam exactly which can be a degenerate point
            # for some surfaces).
            t = (k + 0.5) / n_peri
            peri_p = peri_lo + peri_span * t
            if peri_axis == "u":
                p = surface.Value(peri_p, spine_p)
            else:
                p = surface.Value(spine_p, peri_p)
            pts[i, k] = (p.X(), p.Y(), p.Z())
    return pts


def _filaments_from_lateral_surface_uv(face,
                                         cad_units_per_meter: float,
                                         sigma: float,
                                         n_stations: int,
                                         n_peri: int,
                                         source_tag: str = "step_uv"):
    """Build a PEEC topology from a lateral surface's UV grid.

    For each filament k (k in [0, n_peri-1]):
      polyline = [pts[0,k], pts[1,k], ..., pts[n_stations-1, k]]

    Per-segment cell area is computed from the local cross-section
    polygon (the n_peri samples at one station), so this handles
    variable cross-section automatically (rect, polygon, smooth
    transition between shapes).
    """
    pts_cad = _sample_lateral_surface_uv(face, n_stations, n_peri)
    pts_m = pts_cad / cad_units_per_meter

    # Per-station cross-section area: shoelace polygon area in 3D.
    # Project the n_peri points at each station onto a plane fit to
    # them (= station cross-section plane), then 2D shoelace.
    area_per_station = np.zeros(n_stations, dtype=np.float64)
    for i in range(n_stations):
        ring = pts_m[i]                         # (n_peri, 3)
        center = ring.mean(axis=0)              # (3,)
        # Best-fit plane normal via SVD on (ring - center)
        delta = ring - center
        _, _, vh = np.linalg.svd(delta, full_matrices=False)
        normal = vh[-1] / np.linalg.norm(vh[-1])
        # In-plane axes
        u_axis = vh[0] / np.linalg.norm(vh[0])
        v_axis = np.cross(normal, u_axis)
        u_2d = delta @ u_axis
        v_2d = delta @ v_axis
        # Shoelace
        area_per_station[i] = 0.5 * abs(float(np.sum(
            u_2d * np.roll(v_2d, -1) - np.roll(u_2d, -1) * v_2d)))

    # Build n_peri filaments
    filament_paths = []
    cell_wh = []
    for k in range(n_peri):
        fil = pts_m[:, k, :]
        segs = [(tuple(fil[i]), tuple(fil[i + 1]))
                for i in range(n_stations - 1)]
        filament_paths.append(segs)
        seg_areas = 0.5 * (area_per_station[:-1] + area_per_station[1:])
        cell_wh.append([(float(math.sqrt(a / n_peri)),
                          float(math.sqrt(a / n_peri)))
                         for a in seg_areas])

    import radia  # noqa: F401
    from peec_bundle import build_bundle_solver
    solver, seg_of_fil, port_p, port_m_idx = build_bundle_solver(
        filament_paths, dw=None, dh=None, sigma=sigma, cell_wh=cell_wh)

    return {
        "solver": solver,
        "seg_of_filament": seg_of_fil,
        "filament_paths": filament_paths,
        "cell_wh": cell_wh,
        "n_loop": len(filament_paths),
        "port_plus": port_p,
        "port_minus": port_m_idx,
        "source": source_tag,
        "n_path_pts": n_stations,
        "cross_section_area_m2_mean": float(np.mean(area_per_station)),
        "cross_section_area_m2_min":  float(np.min(area_per_station)),
        "cross_section_area_m2_max":  float(np.max(area_per_station)),
    }


def _filaments_from_circle_edges_per_station(solid,
                                                cad_units_per_meter: float,
                                                sigma: float,
                                                n_peri: int,
                                                source_tag: str = "step_circle_uv"):
    """Phase C-light: per-station variable-radius circle filaments
    from CIRCLE edges in a united multi-loft solid.

    Like ``_collect_circle_edge_centers`` + ``_centerline_from_circle_edge_centers``
    but instead of collapsing every cross-section to the median radius
    (= constant equivalent-circle, 4.14.0), place filaments around each
    station's OWN circle radius via parallel-transport (u_hat, v_hat).
    For constant-radius coils (Kubota 3turncoil) this is identical to
    4.14.0; for tapered / varying circular cross-sections (some IH coil
    designs) this respects the actual local radius.

    Spine policy (2026-05-02 unified):
      * PRIMARY: ``coil_topology.extract_coil_topology`` +
        ``generate_spine`` -- OPEN/CLOSED-aware.  CLOSED full torus
        (no caps) gets a full 360deg spine; OPEN gapped torus gets
        cap_a -> cap_b along the LONG arc.  R_spine is refined from
        the actual mean of the cross-section centroids around the
        rotation axis (more accurate than the bbox 0.85 fallback).
      * FALLBACK: legacy ``_chain_centroids_nn_index`` on the dedupd
        circle centroids.  Used only when topology extraction fails.
        NN-chain handles OPEN reliably (endpoint score detects caps)
        but NOT CLOSED full torus (no endpoints -> ~half-loop only).

    Returns ``None`` (caller falls through to the constant equivalent-
    circle path) when the solid does not have a clean population of
    consistent-radius circle edges.
    """
    from build123d import GeomType
    from OCP.BRepAdaptor import BRepAdaptor_Curve

    circles = [e for e in solid.edges() if e.geom_type == GeomType.CIRCLE]
    if len(circles) < 5:
        return None

    raw = []
    for e in circles:
        try:
            adapt = BRepAdaptor_Curve(e.wrapped)
            circ = adapt.Circle()
            c = circ.Location()
            r = float(circ.Radius())
            ax = circ.Axis().Direction()
            raw.append({
                "center": np.array([c.X(), c.Y(), c.Z()], dtype=np.float64),
                "radius": r,
                "normal": np.array([ax.X(), ax.Y(), ax.Z()],
                                     dtype=np.float64),
            })
        except Exception:
            continue
    if not raw:
        return None

    radii = np.array([d["radius"] for d in raw])
    median_r = float(np.median(radii))
    if median_r <= 0:
        return None
    # 30% radius spread allowed (tapered windings); rejects lead arcs
    # of unrelated radius (e.g. coil terminals) and other artifacts.
    consistent = [d for d, r in zip(raw, radii)
                  if abs(r - median_r) / median_r < 0.3]
    if len(consistent) < 5:
        return None

    # Dedupe semicircle pairs by center proximity (semicircles of one
    # cross-section share the same arc_center).
    dedup_tol = 0.1 * median_r
    kept = []
    for d in consistent:
        if not any(np.linalg.norm(d["center"] - k["center"]) < dedup_tol
                   for k in kept):
            kept.append(d)
    if len(kept) < 5:
        return None

    centers_cad = np.array([d["center"] for d in kept])
    radii_cad_kept = np.array([d["radius"] for d in kept])

    # Multi-turn helix guard: ``coil_topology`` builds a single planar
    # arc spine in the rotation-axis (= z) plane.  For a 3turnCoil-class
    # helix the actual conductor spans multiple turns at different
    # z-levels, so the planar topology spine is geometrically wrong --
    # nearest-neighbour mapping would still produce a single-turn loop
    # at z = 0 with the right cross-section radius (collapsing 3 turns
    # to 1 -- L drops by a factor ~6).  Detect via z-extent of the
    # consistent cross-section centroids: a single-turn coil has
    # z_extent <= cross-section radius; a multi-turn helix has
    # z_extent >> cross-section radius.
    z_extent = float(centers_cad[:, 2].max() - centers_cad[:, 2].min())
    is_multi_turn_helix = (z_extent > 2.0 * median_r)

    # PRIMARY: unified topology spine (OPEN/CLOSED-aware) for single-
    # turn coils.  Multi-turn helix uses the legacy NN-chain (tangent
    # continuity correctly walks the spiral).
    centroids_m = None
    radii_m = None
    if not is_multi_turn_helix:
        try:
            from radia.coil_topology import (
                extract_coil_topology as _extract_topo,
                generate_spine as _gen_spine,
                CoilTopology as _CoilTopology,
            )
            topo = _extract_topo(solid)
            # Refine R_spine from the actual cross-section centroid
            # distance to the rotation axis (axis = z).  The bbox-based
            # 0.85 * R_outer fallback in coil_topology is conservative;
            # using the mean of the cross-section centroid radii is exact
            # for a well-behaved swept geometry.
            R_spine_refined_cad = float(np.mean(
                np.linalg.norm(centers_cad[:, :2], axis=1)))
            topo_refined = _CoilTopology(
                is_open=topo.is_open,
                cap_a=topo.cap_a, cap_b=topo.cap_b,
                theta_a=topo.theta_a, theta_b=topo.theta_b,
                sweep_deg=topo.sweep_deg,
                axis=topo.axis,
                R_spine=R_spine_refined_cad,
                cross_section_kind=topo.cross_section_kind,
            )
            # Use at least as many spine stations as we have circle
            # cross-sections in the original geometry; cap to a sensible
            # upper bound so dense 3turncoil-class coils don't blow up.
            n_stations = max(20, len(kept))
            n_stations = min(n_stations, 200)
            spine_cad = _gen_spine(topo_refined, n_stations)
            # Map each spine station to its nearest circle-edge centroid
            # to recover the per-station radius.  For constant-radius
            # geometry this is identical to a uniform median; for tapered
            # geometry this respects the local radius.
            radii_at_stations_cad = np.zeros(n_stations, dtype=np.float64)
            for i in range(n_stations):
                d2 = np.sum((centers_cad - spine_cad[i]) ** 2, axis=1)
                radii_at_stations_cad[i] = radii_cad_kept[int(np.argmin(d2))]
            centroids_m = spine_cad / cad_units_per_meter
            radii_m = radii_at_stations_cad / cad_units_per_meter
        except Exception:
            centroids_m = None
            radii_m = None

    if centroids_m is None:
        # FALLBACK: legacy NN-chain on circle centroids.  Works for
        # OPEN coils (endpoint score finds caps) but on CLOSED full
        # torus (no endpoints) traces only ~half the loop.
        perm = _chain_centroids_nn_index(centers_cad)
        ordered = [kept[i] for i in perm]
        centroids_m = centers_cad[perm] / cad_units_per_meter
        radii_m = np.array([d["radius"] for d in ordered]) / cad_units_per_meter

    # Parallel-transport frame on the chained centerline.  By
    # construction of the loft, each cross-section's normal is along
    # the local spine tangent, so u_hat / v_hat are already in the
    # cross-section plane and we can use them directly to place
    # filaments around each station's circle.
    _, u_hat, v_hat = _parallel_transport_frame(centroids_m)

    n_path = len(centroids_m)
    theta = 2.0 * np.pi * (np.arange(n_peri) + 0.5) / float(n_peri)
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)

    filament_paths = []
    cell_wh = []
    for k in range(n_peri):
        fil_pts = np.zeros((n_path, 3), dtype=np.float64)
        for i in range(n_path):
            r = radii_m[i]
            fil_pts[i] = (centroids_m[i]
                           + r * cos_t[k] * u_hat[i]
                           + r * sin_t[k] * v_hat[i])
        segs = [(tuple(fil_pts[i]), tuple(fil_pts[i + 1]))
                for i in range(n_path - 1)]
        filament_paths.append(segs)
        # Per-segment cell side from local circle area.  For tapered
        # windings the area changes per segment.
        seg_areas = math.pi * 0.5 * (radii_m[:-1] ** 2 + radii_m[1:] ** 2)
        cell_wh.append([(float(math.sqrt(a / n_peri)),
                          float(math.sqrt(a / n_peri)))
                         for a in seg_areas])

    import radia  # noqa: F401
    from peec_bundle import build_bundle_solver
    solver, seg_of_fil, port_p, port_m_idx = build_bundle_solver(
        filament_paths, dw=None, dh=None, sigma=sigma, cell_wh=cell_wh)

    return {
        "solver": solver,
        "seg_of_filament": seg_of_fil,
        "filament_paths": filament_paths,
        "cell_wh": cell_wh,
        "n_loop": len(filament_paths),
        "port_plus": port_p,
        "port_minus": port_m_idx,
        "source": source_tag,
        "n_path_pts": n_path,
        "cross_section_radius_m_min": float(radii_m.min()),
        "cross_section_radius_m_max": float(radii_m.max()),
        "cross_section_radius_m_mean": float(radii_m.mean()),
    }


def _collect_circle_edge_centers(solid):
    """Cross-section circle centers from CIRCLE edges (united-loft fallback).

    For UNITED multi-turn pancake STEP files where boolean ``unite``
    has consumed the planar end-cap faces but the cross-section
    CIRCLE edges are still present (often split into 2 semicircles
    per cross-section by the unite operation), pull the circle centre
    + radius from each CIRCLE edge's underlying ``Geom_Circle``,
    filter by consistent radius, and dedupe semicircle pairs.

    Returns ``None`` if the solid does not look like a multi-turn
    coil with a clear consistent-radius circle population (lets the
    caller fall through to ``_centerline_from_torus_sweep`` or
    ``_centerline_from_open_spine``).

    Returns ``(centers_cad, median_radius_cad)`` on success:
        centers_cad: list of (3,) np arrays in raw CAD units
        median_radius_cad: float in raw CAD units
    """
    from build123d import GeomType
    from OCP.BRepAdaptor import BRepAdaptor_Curve

    circles = [e for e in solid.edges() if e.geom_type == GeomType.CIRCLE]
    if len(circles) < 5:
        return None

    raw_data = []
    for e in circles:
        try:
            adapt = BRepAdaptor_Curve(e.wrapped)
            circ = adapt.Circle()
            c = circ.Location()
            r = float(circ.Radius())
            raw_data.append((np.array([c.X(), c.Y(), c.Z()],
                                       dtype=np.float64), r))
        except Exception:
            continue

    if not raw_data:
        return None

    radii = np.array([r for _, r in raw_data])
    median_r = float(np.median(radii))
    if median_r <= 0:
        return None

    # Filter by consistent radius (cross-section circles all share the
    # nominal wire radius; any mavericks are e.g. lead arcs of a
    # different radius).
    mask = np.abs(radii - median_r) / median_r < 0.1
    consistent = [c for (c, _r), keep in zip(raw_data, mask) if keep]
    if len(consistent) < 5:
        return None

    # Dedupe near-duplicate centres (semicircles of the same circle
    # share the same arc_center).
    dedup_tol = 0.1 * median_r
    kept = []
    for c in consistent:
        if not any(np.linalg.norm(c - k) < dedup_tol for k in kept):
            kept.append(c)
    if len(kept) < 5:
        return None
    return kept, median_r


def _centerline_from_circle_edge_centers(centers_cad: list,
                                           median_radius_cad: float,
                                           cad_units_per_meter: float = 1.0):
    """Centerline from a list of cross-section circle centres.

    Companion of ``_collect_circle_edge_centers``: chain the centres
    via NN + tangent continuity, expose a constant-circle equivalent
    cross-section ``(width, height)`` in metres for downstream
    ``filaments_from_polyline`` (which uses an equivalent-square side
    derived from the cross-section area).

    Returns the same ``(path_m, widths_m, heights_m)`` shape as the
    other ``_centerline_from_*`` helpers.
    """
    centroids = np.array(centers_cad, dtype=np.float64)
    perm = _chain_centroids_nn_index(centroids)
    ordered = centroids[perm]

    # Cross-section is a circle of radius `median_radius_cad`.
    mean_area_cad2 = math.pi * median_radius_cad * median_radius_cad
    side_cad = math.sqrt(mean_area_cad2)
    n_seg = len(ordered) - 1
    widths_cad = np.full(n_seg, side_cad, dtype=np.float64)
    heights_cad = np.full(n_seg, side_cad, dtype=np.float64)

    scale = 1.0 / cad_units_per_meter
    return ordered * scale, widths_cad * scale, heights_cad * scale


def _chain_centroids_nn_index(pts: np.ndarray) -> list:
    """Return a permutation index that orders 3D points along a smooth
    polyline via nearest-neighbor + tangent continuity.

    Endpoint detection: interior points have 2 close neighbors at
    similar distance; endpoints have only 1.  Rank points by
    second-nearest / first-nearest distance ratio -- endpoints have
    the largest ratio.

    Then greedy nearest-neighbor walk from one endpoint.  For tight
    spiral / pancake geometries where 2 adjacent turns can be nearly
    as close as consecutive cross-sections within a turn, add a
    tangent-continuity bias: once we have 2 points, prefer neighbors
    whose direction continues the current tangent (penalize backward
    / perpendicular jumps to a different turn).

    Returns an index list ``visited`` such that ``pts[visited]`` is
    the ordered polyline.  Caller can apply the same index to any
    aligned per-point data (e.g. cross-section face list).
    """
    n = len(pts)
    if n < 2:
        return list(range(n))
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
            best_idx = None
            best_score = float("inf")
            for i in remaining:
                vec = pts[i] - pts[curr]
                d = np.linalg.norm(vec)
                if d < 1e-12:
                    continue
                cos_a = float(np.dot(vec, tangent) / d)
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
    return visited


def _chain_centroids_nn(pts: np.ndarray) -> np.ndarray:
    """Convenience wrapper: order 3D points using the index walker.

    Returns the reordered array directly.  For callers that need the
    permutation (e.g. to also order an aligned face list), use
    ``_chain_centroids_nn_index``.
    """
    return pts[_chain_centroids_nn_index(pts)]


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


def _centerline_from_topology_spine(solid, n_segments: int,
                                     cad_units_per_meter: float):
    """OPEN/CLOSED-aware spine via ``coil_topology`` + sectioning for w/h.

    Used when ``_centerline_from_revolution_sweep`` raises (typically
    "no PLANE end-cap faces" on a CLOSED full revolution): the
    unified ``coil_topology.extract_coil_topology`` returns a
    correct ``CoilTopology`` for both OPEN and CLOSED coils, and
    ``generate_spine`` produces the correct spine arc (full 360deg
    for CLOSED, cap-aware long arc for OPEN).  Cross-section
    width/height come from a single sectioning at the spine midpoint
    -- equivalent-square side from the section area.

    Raises ``ValueError`` if topology extraction fails or sectioning
    at the midpoint produces no face.  Caller should then fall
    through to ``_centerline_from_open_spine`` (longest-open-edge
    fallback, the historical Path 3 behaviour).
    """
    from build123d import section, Plane, Vector
    try:
        from radia.coil_topology import (
            extract_coil_topology as _extract_topo,
            generate_spine as _gen_spine,
        )
    except Exception as exc:
        raise ValueError(
            f"coil_topology import failed: {exc}") from exc
    try:
        topo = _extract_topo(solid)
    except Exception as exc:
        raise ValueError(
            f"extract_coil_topology failed: {exc}") from exc

    spine_cad = _gen_spine(topo, n_segments + 1)
    n_path = spine_cad.shape[0]

    # Section at the polyline midpoint to get a representative
    # cross-section area.  For a uniform sweep this is the only
    # information we need; tapered sweeps would degrade gracefully
    # to the median area, which is a sound default for the
    # equivalent-circle filaments_from_polyline.
    mid = n_path // 2
    if 0 < mid < n_path - 1:
        c = 0.5 * (spine_cad[mid] + spine_cad[mid + 1])
        t = spine_cad[mid + 1] - spine_cad[mid]
    else:
        c = spine_cad[0]
        t = spine_cad[1] - spine_cad[0] if n_path > 1 else np.array(
            [1.0, 0.0, 0.0])
    t_n = float(np.linalg.norm(t))
    if t_n < 1e-30:
        raise ValueError("topology spine has zero-length midpoint tangent")
    t = t / t_n
    origin = Vector(float(c[0]), float(c[1]), float(c[2]))
    z_dir = Vector(float(t[0]), float(t[1]), float(t[2]))
    sec_plane = Plane(origin=origin, z_dir=z_dir)
    try:
        cross = section(solid, section_by=sec_plane)
    except Exception as exc:
        raise ValueError(
            f"midpoint sectioning failed: {exc}") from exc
    faces = cross.faces() if cross is not None else []
    if not faces:
        raise ValueError("midpoint sectioning produced no face")

    best = min(faces, key=lambda f: (f.center() - origin).length)
    side = math.sqrt(float(best.area))

    widths_cad = np.full(n_segments, side, dtype=np.float64)
    heights_cad = np.full(n_segments, side, dtype=np.float64)

    scale = 1.0 / cad_units_per_meter
    return spine_cad * scale, widths_cad * scale, heights_cad * scale


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
    # (This is the robust path for non-united multi-loft coils.)
    cross_faces = _collect_loft_cross_sections(solid)
    if cross_faces:
        return _centerline_from_cross_sections(solid, cad_units_per_meter)

    # United multi-turn pancake fallback: when boolean unite has
    # consumed the planar end-caps but the cross-section CIRCLE edges
    # are still present (split into 2 semicircles per cross-section
    # by Cubit's unite operation).  Group circles by arc_center, dedupe
    # the semicircle pair into one per cross-section, chain via NN +
    # tangent continuity.  Handles Kubota's 3turncoil united.stp class.
    circle_centers_radius = _collect_circle_edge_centers(solid)
    if circle_centers_radius is not None:
        centers_cad, median_r_cad = circle_centers_radius
        return _centerline_from_circle_edge_centers(
            centers_cad, median_r_cad, cad_units_per_meter)

    # Torus-shaped single-loop coil (gapped torus, full torus): extract
    # major / minor radius + axis + sweep angle analytically from the
    # TORUS face parameters.  Handles the single-loop case that the
    # open-spine fallback gets wrong (picks a cross-section arc as
    # "longest open edge", path length comes out half the real arc).
    try:
        return _centerline_from_torus_sweep(solid, n_segments, cad_units_per_meter)
    except ValueError:
        pass

    # Unified-topology spine: catches CLOSED full-revolution coils that
    # ``_centerline_from_torus_sweep`` rejects for "no PLANE end-cap
    # faces".  ``coil_topology.extract_coil_topology`` returns
    # is_open=False for CLOSED and ``generate_spine`` produces the
    # correct full 360deg spine; sectioning at the spine midpoint
    # recovers the cross-section area (equivalent-circle radius).
    # Without this path, CLOSED coils fall through to
    # ``_centerline_from_open_spine`` which picks a half-arc seam edge
    # and produces a 178deg spine -- the regression-test failure on
    # ih_closed_torus_coil.step.
    try:
        return _centerline_from_topology_spine(
            solid, n_segments, cad_units_per_meter)
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
        # Perimeter-only placement.  Three paths, in priority order:
        #
        # 1. **Lateral surface UV sampling** (4.15.0+):
        #    OCC's parametric Geom_Surface.Value(u, v) sampled on a
        #    grid (n_stations x n_peri).  Handles ANY cross-section
        #    shape (circle, rect, polygon, arbitrary) and ANY
        #    variation along the spine (tapered, transition between
        #    shapes), AND united multi-loft solids whose lateral
        #    surfaces are merged into one big BSPLINE.  This is the
        #    most general path.
        #
        # 2. **Per-station from loft cross-section faces** (4.14.0):
        #    if the STEP is a NON-united loft-of-profiles solid,
        #    sample each cross-section face's outer wire directly.
        #    Triggered when path 1 cannot identify a single dominant
        #    closed-perimeter lateral surface (e.g. multi-loft with
        #    each lateral as a separate face).
        #
        # 3. **Constant equivalent-circle** (legacy):
        #    when neither (1) nor (2) apply (gapped torus via TORUS
        #    sweep, swept helix via open-spine), use mean cross-
        #    section AREA + equivalent-circle radius.  Pre-4.14.0
        #    behaviour.
        import numpy as np
        from build123d import import_step
        solid = import_step(step_path)

        # Path 1: BSPLINE / TORUS / CYLINDER lateral surface UV grid.
        lateral = _find_lateral_surface(solid)
        if lateral is not None:
            try:
                return _filaments_from_lateral_surface_uv(
                    lateral, cad_units_per_meter=cad_units_per_meter,
                    sigma=sigma,
                    n_stations=max(20, n_slices // 2),
                    n_peri=n_peri,
                    source_tag="step_uv")
            except ValueError:
                # Surface couldn't be sampled (no closed UV axis,
                # degenerate shape).  Fall through.
                pass

        # Path 2: per-station planar end-cap faces (NON-united loft).
        loft = _try_extract_loft_with_profile(
            solid, cad_units_per_meter=cad_units_per_meter)
        if loft is not None:
            path_m, faces_ordered = loft
            return _filaments_from_per_station_faces(
                path_m, faces_ordered,
                sigma=sigma, n_peri=n_peri,
                source_tag="step_per_station")

        # Path 2b: per-station VARIABLE-RADIUS CIRCLE from CIRCLE
        # edges (united multi-loft with circular cross-section, possibly
        # tapered).  Kubota 3turncoil class.  Stronger than the
        # equivalent-circle fallback because each station gets its own
        # radius (handles tapered windings); identical for constant-
        # radius coils.
        topo_circle = _filaments_from_circle_edges_per_station(
            solid, cad_units_per_meter=cad_units_per_meter,
            sigma=sigma, n_peri=n_peri,
            source_tag="step_circle_uv")
        if topo_circle is not None:
            return topo_circle

        # Path 2c (Phase C-heavy, v4.24.0+): united multi-loft with
        # NON-circular cross-section (rect / polygon / shape transition).
        # Sections the solid at n_stations along the spine to recover
        # cross-section faces, then reuses Tier 2's per-station-faces
        # filament placement.  Slow (BRepAlgoAPI_Section per station)
        # but works on cases where Tier 1 / Tier 2 / Tier 2b all trip.
        # GATED: only fire if the lateral has many PLANE faces (>= 4)
        # AND no revolution-type surfaces (TORUS / CYLINDER / REVOLUTION)
        # -- those are handled by Path 1 / Path 2b at higher accuracy.
        # Without this gate we'd intercept gapped-torus (circular
        # cross-section, TORUS lateral) and degrade its accuracy.
        from build123d import GeomType
        face_geom_counts = {}
        for f in solid.faces():
            face_geom_counts[f.geom_type] = face_geom_counts.get(f.geom_type, 0) + 1
        n_plane = face_geom_counts.get(GeomType.PLANE, 0)
        n_revolution_like = sum(face_geom_counts.get(g, 0) for g in (
            GeomType.TORUS, GeomType.CYLINDER, GeomType.CONE,
            GeomType.REVOLUTION,
        ))
        if n_plane >= 4 and n_revolution_like == 0:
            topo_section = _filaments_from_section_planes(
                solid, cad_units_per_meter=cad_units_per_meter,
                sigma=sigma, n_peri=n_peri,
                n_stations=20,
                source_tag="step_section_planes")
            if topo_section is not None:
                return topo_section

        # Path 3: constant equivalent-circle via existing dispatch.
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


# ----------------------------------------------------------------------
# Variable cross-section path (4.14.0+):
#   Sample each cross-section face's outer boundary directly in the
#   parallel-transport frame so filament placement tracks rect /
#   polygon / arbitrary cross-section shape AND any variation along
#   the spine (tapered conductor, varying-width pancake, etc.).
# ----------------------------------------------------------------------

def _sample_face_perimeter_in_pt_frame(face, centroid_3d: np.ndarray,
                                          u_hat: np.ndarray,
                                          v_hat: np.ndarray,
                                          n_peri: int) -> np.ndarray:
    """Return (n_peri, 2) UV samples of a planar face's outer boundary.

    The face's outer wire is sampled at n_peri arc-length-equispaced
    positions and the 3D points are projected onto the parallel-
    transport frame ``(u_hat, v_hat)`` rooted at ``centroid_3d``.

    Args:
        face: build123d Face whose outer boundary is sampled.
        centroid_3d: (3,) face centroid in the same coordinate system.
        u_hat: (3,) parallel-transport in-plane axis at this station.
        v_hat: (3,) parallel-transport in-plane axis (= n x u_hat).
        n_peri: number of samples along the perimeter.

    Returns:
        (n_peri, 2) array of (u, v) offsets, ordered by arc-length
        traversal of the boundary starting from the boundary point
        closest to +u_hat.
    """
    from build123d import PositionMode

    wire = face.outer_wire()
    edges = wire.edges()
    seg_len = np.array([float(e.length) for e in edges], dtype=np.float64)
    total_len = float(seg_len.sum())
    if total_len <= 0 or len(edges) == 0:
        raise ValueError("face has zero-perimeter outer wire")

    cumlen = np.cumsum(np.concatenate([[0.0], seg_len]))
    s_targets = (np.arange(n_peri) + 0.5) / n_peri * total_len

    samples_3d = np.zeros((n_peri, 3), dtype=np.float64)
    for k, s in enumerate(s_targets):
        # Find which edge holds s
        idx = int(np.searchsorted(cumlen[1:], s, side="right"))
        idx = min(idx, len(edges) - 1)
        s_local = float(s - cumlen[idx])
        edge = edges[idx]
        # build123d position_at uses parameter [0, 1] by default; pass
        # PositionMode.LENGTH to sample by arc length within the edge.
        pos = edge.position_at(s_local, position_mode=PositionMode.LENGTH)
        samples_3d[k] = (pos.X, pos.Y, pos.Z)

    # Project to parallel-transport UV
    delta = samples_3d - centroid_3d
    u_vals = delta @ u_hat
    v_vals = delta @ v_hat
    uv = np.column_stack([u_vals, v_vals])

    # Roll the array so index 0 lands at the boundary point closest
    # to +u_hat (smallest signed angle from +u in the (u,v) plane).
    # This stabilises filament k across stations: all stations agree
    # on which corner / side of the rect filament 0 is on.
    angles = np.arctan2(v_vals, u_vals)
    # Distance from theta=0, wrapped: |angle| in (-pi, pi)
    abs_ang = np.abs((angles + np.pi) % (2 * np.pi) - np.pi)
    start = int(np.argmin(abs_ang))
    uv = np.roll(uv, -start, axis=0)
    return uv


def _filaments_from_per_station_faces(centroids_m: np.ndarray,
                                        faces_ordered: list,
                                        sigma: float = 5.8e7,
                                        n_peri: int = 16,
                                        source_tag: str = "step_per_station"):
    """Build a PEEC topology from per-station planar faces.

    Each face's outer boundary is sampled at n_peri arc-length-equi-
    spaced points in the parallel-transport frame at that station.
    Filament k connects the k-th sample at station i to the k-th
    sample at station i+1.  This handles **variable cross-section**
    along the spine (rect / polygon / arbitrary shape, possibly
    tapered).  Also handles constant-circular cross-section (gives
    same result as ``filaments_from_polyline`` to within numerical
    noise).

    Args:
        centroids_m: (N, 3) cross-section centroids in metres,
            ordered along the spine.
        faces_ordered: list of N build123d Face objects (planar
            cross-section caps) aligned with centroids.  faces[i] is
            the cap whose centroid is centroids_m[i] (BEFORE the
            cad_units_per_meter scale; the caller is responsible for
            the centroid scale -- the face still lives in raw CAD
            units).
        sigma: conductivity.
        n_peri: number of perimeter filaments.
        source_tag: free-form string for the result dict's "source".

    Returns:
        topology_dict, same shape as filaments_from_step.
    """
    n_path = len(centroids_m)
    if n_path < 2:
        raise ValueError(
            f"need at least 2 stations, got {n_path}")
    if len(faces_ordered) != n_path:
        raise ValueError(
            f"faces_ordered length {len(faces_ordered)} does not match "
            f"centroids length {n_path}")

    _, u_hat, v_hat = _parallel_transport_frame(centroids_m)

    # Per-station UV (n_path, n_peri, 2) and per-station perimeter (n_path,)
    uv_per_station = np.zeros((n_path, n_peri, 2), dtype=np.float64)
    area_per_station = np.zeros(n_path, dtype=np.float64)
    for i in range(n_path):
        face = faces_ordered[i]
        c = face.center()
        c_np = np.array([c.X, c.Y, c.Z], dtype=np.float64)
        # The face is in raw CAD units; centroids_m is already scaled
        # to metres.  But the UV are RELATIVE to the face centroid, so
        # the relative offsets (samples - centroid) need to be scaled
        # by the same factor as the centerline.  Since the relative
        # offsets are dimensionless w.r.t. the centerline scale, we
        # pull the scale from the matching centerline -> face delta.
        # In practice for shipped samples both live in metres (Cubit
        # mks).  We compute UV in CAD units and then scale below.
        uv_cad = _sample_face_perimeter_in_pt_frame(
            face, c_np, u_hat[i], v_hat[i], n_peri)
        uv_per_station[i] = uv_cad
        area_per_station[i] = float(face.area)

    # Per-station UV is in CAD units of the face.  centroids_m is in
    # metres after the caller's cad_units_per_meter scale.  Scale UV
    # to metres via the inverse: detect by comparing one face centroid
    # in raw units vs the matching centerline point in metres.  This
    # is cleaner than threading cad_units_per_meter through here.
    f0 = faces_ordered[0]
    c0_cad = np.array([f0.center().X, f0.center().Y, f0.center().Z],
                       dtype=np.float64)
    c0_m = centroids_m[0]
    cad_to_m = np.linalg.norm(c0_m) / max(np.linalg.norm(c0_cad), 1e-30) \
        if np.linalg.norm(c0_cad) > 1e-30 else 1.0
    uv_per_station *= cad_to_m
    area_per_station *= cad_to_m * cad_to_m

    # Build n_peri filaments, each as a polyline through the matching
    # k-th sample at every station.
    filament_paths = []
    cell_wh = []
    for k in range(n_peri):
        fil_pts = np.zeros((n_path, 3), dtype=np.float64)
        for i in range(n_path):
            fil_pts[i] = (centroids_m[i]
                           + uv_per_station[i, k, 0] * u_hat[i]
                           + uv_per_station[i, k, 1] * v_hat[i])
        segs = [(tuple(fil_pts[i]), tuple(fil_pts[i + 1]))
                for i in range(n_path - 1)]
        filament_paths.append(segs)
        # Per-segment filament cell area = local face area / n_peri
        # (face provides total cross-section); side = sqrt(area / n_peri)
        seg_areas = 0.5 * (area_per_station[:-1] + area_per_station[1:])
        cell_wh.append([(float(math.sqrt(a / n_peri)),
                          float(math.sqrt(a / n_peri)))
                         for a in seg_areas])

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
        "n_path_pts": n_path,
        "cross_section_area_m2_mean": float(np.mean(area_per_station)),
        "cross_section_area_m2_min": float(np.min(area_per_station)),
        "cross_section_area_m2_max": float(np.max(area_per_station)),
    }


def _section_solid_at_plane(solid, point_xyz, normal_xyz):
    """Section a build123d Solid by a plane at point_xyz with given normal.

    Returns a build123d Face if the section produced a single closed
    cross-section, or None if the section failed / produced multiple
    or no closed wires.

    Uses OCP BRepAlgoAPI_Section + BRepBuilderAPI_MakeWire +
    BRepBuilderAPI_MakeFace.  Phase C-heavy helper (2026-05-02).
    """
    try:
        from OCP.BRepAlgoAPI import BRepAlgoAPI_Section
        from OCP.gp import gp_Pln, gp_Pnt, gp_Dir
        from OCP.BRepBuilderAPI import (
            BRepBuilderAPI_MakeWire, BRepBuilderAPI_MakeFace,
        )
        from OCP.TopExp import TopExp_Explorer
        from OCP.TopAbs import TopAbs_EDGE
        from OCP.TopoDS import TopoDS
        from build123d import Face
    except ImportError:
        return None

    # Build cutting plane
    gp_origin = gp_Pnt(float(point_xyz[0]),
                       float(point_xyz[1]),
                       float(point_xyz[2]))
    gp_normal = gp_Dir(float(normal_xyz[0]),
                       float(normal_xyz[1]),
                       float(normal_xyz[2]))
    plane = gp_Pln(gp_origin, gp_normal)

    # Section: returns a Compound of edges
    sec = BRepAlgoAPI_Section(solid.wrapped, plane, False)
    sec.ComputePCurveOn1(True)
    sec.Approximation(False)
    sec.Build()
    if not sec.IsDone():
        return None

    # Collect edges from the section result
    section_shape = sec.Shape()
    edges = []
    exp = TopExp_Explorer(section_shape, TopAbs_EDGE)
    while exp.More():
        edges.append(TopoDS.Edge_s(exp.Current()))
        exp.Next()
    if not edges:
        return None

    # Build wires from the edges. ConnectEdgesToWires does the ordering.
    # NOTE: a single cutting plane can intersect the solid at MULTIPLE
    # locations (e.g. the y=0 plane cuts a torus at theta=0 AND
    # theta=180).  ConnectEdgesToWires returns one wire per intersection
    # locus.  We pick the wire whose centroid is closest to the
    # query point, NOT the largest wire (which was the previous bug).
    try:
        from OCP.ShapeAnalysis import ShapeAnalysis_FreeBounds
        from OCP.TopTools import TopTools_HSequenceOfShape
        from OCP.GProp import GProp_GProps
        from OCP.BRepGProp import BRepGProp
        edge_seq = TopTools_HSequenceOfShape()
        for e in edges:
            edge_seq.Append(e)
        wires_seq = TopTools_HSequenceOfShape()
        ShapeAnalysis_FreeBounds.ConnectEdgesToWires_s(
            edge_seq, 1e-6, False, wires_seq)
        if wires_seq.Length() < 1:
            return None

        query = np.asarray(point_xyz, dtype=np.float64)
        best_wire = None
        best_dist = float("inf")
        for i in range(1, wires_seq.Length() + 1):
            w_ds = TopoDS.Wire_s(wires_seq.Value(i))
            # Wire centroid via LinearProperties center-of-mass
            props = GProp_GProps()
            BRepGProp.LinearProperties_s(w_ds, props)
            com = props.CentreOfMass()
            d = np.linalg.norm(np.array([com.X(), com.Y(), com.Z()]) - query)
            if d < best_dist:
                best_dist = d
                best_wire = w_ds
        if best_wire is None:
            return None

        # Make a planar face from the wire
        face_maker = BRepBuilderAPI_MakeFace(plane, best_wire, True)
        if not face_maker.IsDone():
            return None
        ds_face = face_maker.Face()
        return Face(ds_face)
    except Exception:
        return None


def _filaments_from_section_planes(solid,
                                    cad_units_per_meter: float,
                                    sigma: float,
                                    n_peri: int,
                                    n_stations: int,
                                    source_tag: str = "step_section_planes"):
    """Phase C-heavy fallback: united multi-loft (any cross-section).

    Used when:
      * Tier 1 lateral-surface UV trips (multi-fragment lateral)
      * Tier 2 per-station faces trips (no surviving end-cap faces, e.g.
        after ``unite()`` on a multi-loft)
      * Tier 2b circle-edge per-station trips (non-circular cross-section)

    Algorithm:
      1. Extract spine polyline + per-station tangents (existing
         ``extract_centerline_from_step`` falls back to longest-edge
         spine; we reuse its result and recompute tangents here).
      2. At each station, build a cutting plane perpendicular to the
         local tangent, section the solid, recover the cross-section
         Face from the resulting edges.
      3. Pass (centroids, faces) to ``_filaments_from_per_station_faces``
         (existing Tier 2 algorithm) for the actual filament placement.

    Returns None if sectioning fails on >50 % of stations (the spine
    is unreliable) or if any structural assumption breaks.

    Cost: one ``BRepAlgoAPI_Section`` call per station.  ~5-30 s for
    a 10 MB STEP at n_stations=20.  Acceptable for production; the
    much faster Tier 1 / Tier 2 paths are tried first.
    """
    if solid is None:
        return None
    # Spine polyline.  PRIMARY PATH (2026-05-02): the unified
    # ``coil_topology.extract_coil_topology`` + ``generate_spine``
    # which is OPEN/CLOSED-aware (cap detection + cap-aware arc
    # endpoints).  This replaces the prior call to
    # ``_centerline_from_solid_geometry`` whose
    # ``_spine_from_rotation_axis_z`` sub-path always sampled
    # ``np.linspace(0, 2*pi, n, endpoint=False)``, clipping ~14 deg
    # of conductor on a 355 deg gapped torus.  Fallback to the legacy
    # multi-path extractor if topology extraction fails so previously-
    # working geometries don't regress.
    #
    # ``spine_is_topology_ordered`` records whether the spine came from
    # the unified extractor (which already orders stations cap_a ->
    # ... -> cap_b along the LONG arc).  When True, we SKIP the
    # downstream NN-chain re-ordering on the recovered cross-section
    # centroids: NN-chain assumes Euclidean nearest neighbours, which
    # for OPEN coils with a small angular gap WRONGLY hops cap_a ->
    # cap_b across the gap (e.g. 5 deg vs 18 deg/station along the
    # conductor) and the resulting filament zigzags through the gap.
    # The unified spine already gives the correct order; only the
    # legacy longest-edge / rotation_axis_z paths needed NN-chain
    # because those samplers visited centroids in a non-spine order.
    spine_is_topology_ordered = False
    topo = None
    try:
        from radia.coil_topology import (
            extract_coil_topology as _extract_topo,
            generate_spine as _gen_spine,
        )
        topo = _extract_topo(solid)
        path_cad = _gen_spine(topo, n_stations)
        path_m = path_cad / cad_units_per_meter
        spine_is_topology_ordered = True
    except Exception:
        topo = None
        try:
            path_m_via_solid = _centerline_from_solid_geometry(
                solid, n_segments=n_stations,
                cad_units_per_meter=cad_units_per_meter)
            if path_m_via_solid is None:
                return None
            path_m = path_m_via_solid
        except Exception:
            return None

    n_path = len(path_m)
    if n_path < 3:
        return None

    # Convert path_m back to CAD units for sectioning (which works in
    # the solid's native frame).
    path_cad = path_m * cad_units_per_meter

    # Tangents via central differences
    tangents = np.zeros_like(path_cad)
    for i in range(n_path):
        if i == 0:
            t = path_cad[1] - path_cad[0]
        elif i == n_path - 1:
            t = path_cad[-1] - path_cad[-2]
        else:
            t = path_cad[i + 1] - path_cad[i - 1]
        nrm = np.linalg.norm(t)
        if nrm > 1e-30:
            tangents[i] = t / nrm
        else:
            tangents[i] = np.array([1.0, 0.0, 0.0])

    # OPEN coils: the spine endpoints sit on cap_a / cap_b, so the
    # cutting plane the section call would build is COINCIDENT with
    # the cap face.  BRepAlgoAPI_Section then returns either a
    # degenerate (~0 area) or multi-disjoint wire from grazing the
    # cap, the area-outlier filter below drops the station, and the
    # loft end cells get no filaments.  Skip the section call at the
    # endpoints and use the cap face directly -- it IS exactly the
    # cross-section we want there (centroid, area, normal all match
    # the conductor geometry by construction).
    cap_a_face = None
    cap_b_face = None
    if topo is not None and getattr(topo, "is_open", False):
        cap_a_face = getattr(topo, "cap_a", None)
        cap_b_face = getattr(topo, "cap_b", None)

    # Section at each station
    faces_attempted = []
    centroids_attempted = []
    areas_attempted = []
    for i in range(n_path):
        if i == 0 and cap_a_face is not None:
            face = cap_a_face
        elif i == n_path - 1 and cap_b_face is not None:
            face = cap_b_face
        else:
            face = _section_solid_at_plane(
                solid, path_cad[i], tangents[i])
        if face is None:
            continue
        faces_attempted.append(face)
        centroids_attempted.append(path_m[i])
        areas_attempted.append(float(face.area))

    if len(faces_attempted) < max(3, int(0.5 * n_path)):
        return None

    # Robust outlier filter: at INTERIOR spine stations the cutting
    # plane can catch the gap or multiple disjoint pieces, producing
    # an inflated cross-section area.  Drop stations whose area
    # deviates by >30 % from the median.  Cap-face endpoints
    # (replaced above) are excluded from the filter so they always
    # contribute a filament.
    cap_endpoint_idx = set()
    if cap_a_face is not None and len(faces_attempted) > 0 and \
            faces_attempted[0] is cap_a_face:
        cap_endpoint_idx.add(0)
    if cap_b_face is not None and len(faces_attempted) > 0 and \
            faces_attempted[-1] is cap_b_face:
        cap_endpoint_idx.add(len(faces_attempted) - 1)
    areas_arr = np.array(areas_attempted)
    interior_mask = np.array(
        [i not in cap_endpoint_idx for i in range(len(areas_arr))],
        dtype=bool)
    if int(np.sum(interior_mask)) > 0:
        median_area = float(np.median(areas_arr[interior_mask]))
    else:
        median_area = float(np.median(areas_arr))
    if median_area <= 0:
        return None
    keep_mask = np.abs(areas_arr - median_area) <= 0.3 * median_area
    # Always keep cap-face endpoints regardless of area (cap area can
    # legitimately differ slightly from the lofted-section median).
    for i in cap_endpoint_idx:
        keep_mask[i] = True
    if int(np.sum(keep_mask)) < max(3, int(0.5 * len(areas_arr))):
        # Outlier filter would drop too many; bail out so a downstream
        # tier (or the equivalent-circle Path 3) can take a shot.
        return None

    faces_kept = [faces_attempted[i]
                  for i in range(len(faces_attempted)) if keep_mask[i]]
    centroids_from_sections_cad = np.array(
        [[f.center().X, f.center().Y, f.center().Z]
         for f in faces_kept],
        dtype=np.float64,
    )
    if spine_is_topology_ordered:
        # The unified topology spine already orders stations cap_a ->
        # ... -> cap_b along the LONG arc.  Sectioning preserves this
        # order (sections returned in spine traversal order); NN-chain
        # would WRONGLY hop cap_a -> cap_b across the small angular
        # gap and break the filament path.  See the spine-construction
        # block above for the full rationale.
        centroids_ordered_cad = centroids_from_sections_cad
        faces_ordered = faces_kept
    else:
        # Legacy spine (longest-edge / rotation_axis_z fallback): the
        # sectioning order is unreliable, so re-chain by nearest-
        # neighbour to recover the spine traversal.
        perm = _chain_centroids_nn_index(centroids_from_sections_cad)
        centroids_ordered_cad = centroids_from_sections_cad[perm]
        faces_ordered = [faces_kept[i] for i in perm]
    centroids_kept_np = centroids_ordered_cad / cad_units_per_meter
    return _filaments_from_per_station_faces(
        centroids_kept_np, faces_ordered,
        sigma=sigma, n_peri=n_peri, source_tag=source_tag)


def _spine_from_rotation_axis_z(solid,
                                 n_segments: int = 30,
                                 cad_units_per_meter: float = 1.0):
    """Generate a spine arc around the z-axis from solid bbox + centroid.

    Used by Phase C-heavy when the solid is a sweep / loft around the
    z-axis (covers gapped torus, multi-turn pancake, rect / polygon
    cross-section, united / non-united multi-loft alike).

    Algorithm:
      1. Bbox in (x, y); set rotation axis = z-axis through origin
         (heuristic — assumes the user has the coil's revolution axis
         on the z-axis, which the Cubit and build123d helpers default
         to).
      2. Estimate spine radius R from the average of bbox extremes
         (typical (x_max - x_min) / 2 for a coil centered at origin).
      3. Sweep ``n_segments`` angular samples uniformly in [0, 2π).
         Sections that fail or produce outlier areas are dropped
         downstream (in `_filaments_from_section_planes`).

    Returns the path in metres (n_segments points, NOT closed).
    """
    bbox = solid.bounding_box()
    # Spine radius: use mean of horizontal half-widths
    # (works for coils centered on the z-axis).
    R_x = max(abs(bbox.max.X), abs(bbox.min.X))
    R_y = max(abs(bbox.max.Y), abs(bbox.min.Y))
    R_outer = max(R_x, R_y)
    # The spine sits INSIDE the cross-section's outer edge.  For a
    # rect cross-section of half-width w, R_spine = R_outer - w.
    # We don't know w yet; use 0.85 * R_outer as a first guess
    # (typical for coils where w is ~10-20 % of R).
    R_spine = 0.85 * R_outer
    if R_spine <= 0:
        return None
    thetas = np.linspace(0.0, 2.0 * math.pi, n_segments, endpoint=False)
    path = np.column_stack([
        R_spine * np.cos(thetas),
        R_spine * np.sin(thetas),
        np.zeros_like(thetas),
    ])
    # Convert to metres (path is already in CAD units)
    return path / cad_units_per_meter


def _centerline_from_solid_geometry(solid,
                                     n_segments: int = 30,
                                     cad_units_per_meter: float = 1.0):
    """Spine polyline from a build123d Solid (no step_path dependency).

    Tries (1) revolution-sweep, (2) z-axis rotation-symmetry, (3)
    longest-edge fallback, in that order.  Returns the path in metres
    (n_segments+1 points).
    """
    # Try revolution-sweep
    try:
        path_m, _, _ = _centerline_from_revolution_sweep(
            solid, n_segments, cad_units_per_meter)
        if path_m is not None and len(path_m) >= 3:
            return path_m
    except Exception:
        pass
    # Try z-axis rotation-symmetric spine (rect / polygon cross-section
    # sweeps around z-axis: NOT triggered by revolution-sweep because
    # the lateral surfaces are PLANE, not TORUS / CYLINDER / REVOLUTION).
    try:
        path_m = _spine_from_rotation_axis_z(
            solid, n_segments, cad_units_per_meter)
        if path_m is not None and len(path_m) >= 3:
            return path_m
    except Exception:
        pass
    # Try open-spine longest-edge
    try:
        path_m, _, _ = _centerline_from_open_spine(
            solid, n_segments, cad_units_per_meter)
        if path_m is not None and len(path_m) >= 3:
            return path_m
    except Exception:
        pass
    return None


def _try_extract_loft_with_profile(solid,
                                     cad_units_per_meter: float = 1.0):
    """Loft path that ALSO returns the ordered cross-section faces.

    Mirrors ``_centerline_from_cross_sections`` but returns
    ``(path_m, faces_ordered)`` where faces[i] aligns with path[i].

    Returns None if the solid does not look like a loft of profiles
    (caller falls back to the regular ``extract_centerline_from_step``).
    """
    cross = _collect_loft_cross_sections(solid)
    if not cross:
        return None
    # Centroids in raw CAD units
    centroids_raw = np.array(
        [[c.X, c.Y, c.Z] for c in (f.center() for f in cross)],
        dtype=np.float64)
    mean_area = float(np.mean([f.area for f in cross]))
    eq_radius = math.sqrt(mean_area / math.pi)
    dedup_tol = 0.1 * eq_radius

    kept_centroids = []
    kept_faces = []
    for c, face in zip(centroids_raw, cross):
        if not any(np.linalg.norm(c - k) < dedup_tol
                   for k in kept_centroids):
            kept_centroids.append(c)
            kept_faces.append(face)
    centroids = np.array(kept_centroids, dtype=np.float64)
    if len(centroids) < 3:
        return None

    perm = _chain_centroids_nn_index(centroids)
    ordered_centroids = centroids[perm]
    ordered_faces = [kept_faces[i] for i in perm]

    scale = 1.0 / cad_units_per_meter
    path_m = ordered_centroids * scale
    return path_m, ordered_faces
