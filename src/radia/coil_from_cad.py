"""coil_from_cad.py -- Extract PEEC filaments from a coil STEP solid.

Public entry points:

* ``filaments_from_step(step_path, ...)`` -- end-to-end STEP -> PEEC
  filament topology.  Caller-selected modes (No-Fallbacks: no
  automatic switchover): the default coil_builder walker path
  (volume-grid nwinc x nhinc), the ``n_peri`` perimeter-only UV tiers
  (thin-skin IH production), and the legacy C++ ExpandFilaments grid
  (``use_coil_builder=False``).
* ``filaments_from_shape(shape, ...)`` -- the same walker path for an
  in-memory build123d shape.
* ``extract_centerline_from_step(step_path, ...)`` -- centerline +
  per-segment cross-sections.  ONE geometric marching engine (the
  walking-plane march) plus exact CAD-feature fast paths; see its
  docstring for the predicate architecture.
* ``build_peec_from_path(path, widths, heights, ...)`` -- polyline +
  cross-sections -> PEECBuilder topology (consumed by the panels'
  parametric sources and the legacy grid mode).

Units: paths and cross-sections are returned in METRES (CLAUDE.md
"Radia always uses meters"); STEP-native CAD units enter through the
explicit ``cad_units_per_meter`` (no auto-detection).
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
    from radia.peec_matrices import PEECBuilder

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
        from radia._b3d_shim import GeomType
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
    area AND has at least one closed UV axis (so that the cross-
    section perimeter wraps around the spine).  This function returns
    it for direct (u, v) sampling by ``_filaments_from_lateral_surface_uv``.

    Returns ``None`` (caller dispatches to a different predicate) when:
      - There are no lateral-type candidates
      - There are too many candidates (multi-fragment loft, united
        multi-turn pancake) -- a single face only covers part of
        the spine, sampling it would give a fragmented result
      - The largest candidate's area is not clearly dominant
        (gapped torus split into 4 TORUS quadrants by Cubit's
        webcut at the gap, etc.)
      - The candidate's UV parameter range has neither U nor V closed
        (cannot identify a perimeter direction; e.g. keiko-class
        "arc + lead bars" loft where the lateral is split into two
        BSPLINE halves at the z=0 equator)

    The threshold "largest >= 80% of total" is a strong signal that
    we have a single-piece lateral surface.  Anything below dispatches
    to the per-station-faces / circle-edge / open-spine paths
    instead.  This predicate must be PRECISE -- once it returns a
    face, the downstream UV sampling must succeed.  Per CLAUDE.md
    "No Fallbacks -- Fail Fast, Fail Loud", we do NOT try-and-recover
    inside ``filaments_from_step``.
    """
    from radia._b3d_shim import GeomType
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
    # Anything below 0.8 indicates fragmentation -- different dispatch.
    if dominance < 0.8:
        return None
    # Also reject if there are too many candidates even if the largest
    # is technically dominant (defensive bound).
    if len(candidates) > 3:
        return None
    # UV-closure check: the lateral must have at least one closed UV
    # axis for ``_sample_lateral_surface_uv`` to identify a perimeter
    # direction.  Doing this here makes ``_find_lateral_surface`` a
    # complete precondition for Path 1 -- no try/except downstream.
    candidate = candidates[0][1]
    from OCP.BRep import BRep_Tool
    surface = BRep_Tool.Surface_s(candidate.wrapped)
    if not (bool(surface.IsUClosed()) or bool(surface.IsVClosed())):
        return None
    return candidate


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
    from radia.peec_bundle import build_bundle_solver
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


def _detect_lead_bars_cad(solid, median_radius_cad: float):
    """Detect straight lead bars in a coil STEP as CYLINDER faces.

    Each detected lead is a cylindrical wire segment whose radius
    matches the helix cross-section radius (within 10%) and whose
    axial length is at least 5x the radius (i.e. not a short cylindrical
    fillet).  Multiple CYLINDER faces may share an axis (e.g. when the
    boolean unite splits the lateral cylinder into 2 half-faces); they
    are merged into a single lead by axis identity.

    For 3turnCoil_work_coil.step this returns the 2 lead bars at
    y=+/-12.5 mm extending 60-61 mm in +/-X from the helix end to the
    lead terminal.

    Args:
        solid: build123d Solid (CAD units, NOT yet scaled to meters).
        median_radius_cad: median cross-section radius of the helix
            (CAD units).  CYLINDER faces with a different radius are
            rejected.

    Returns:
        list of lead dicts, each with:
            "loc":      axis Location() as np.array(3,) in CAD units
            "dir":      unit axis direction as np.array(3,)
            "length":   total axial length (CAD units)
            "radius":   cylinder radius (CAD units)
            "cap_a":    one cap center (CAD units, np.array(3,))
            "cap_b":    other cap center (CAD units)
        Empty list if no qualifying lead bars are present.
    """
    from radia._b3d_shim import GeomType
    from OCP.BRepAdaptor import BRepAdaptor_Surface

    cyl_faces = [f for f in solid.faces() if f.geom_type == GeomType.CYLINDER]
    leads = []
    for f in cyl_faces:
        try:
            adapt = BRepAdaptor_Surface(f.wrapped)
            cyl = adapt.Cylinder()
            ax = cyl.Axis()
            loc = np.array([ax.Location().X(),
                            ax.Location().Y(),
                            ax.Location().Z()], dtype=np.float64)
            direc = np.array([ax.Direction().X(),
                              ax.Direction().Y(),
                              ax.Direction().Z()], dtype=np.float64)
            r = float(cyl.Radius())
        except Exception:
            continue
        # Radius must match median cross-section radius (10% spread).
        if median_radius_cad > 0 and abs(r - median_radius_cad) / median_radius_cad > 0.1:
            continue
        # Compute axial extent from face bounding box (project corners).
        bb = f.bounding_box()
        bb_min = np.array([bb.min.X, bb.min.Y, bb.min.Z], dtype=np.float64)
        bb_max = np.array([bb.max.X, bb.max.Y, bb.max.Z], dtype=np.float64)
        ts = []
        for cx in (bb_min[0], bb_max[0]):
            for cy in (bb_min[1], bb_max[1]):
                for cz in (bb_min[2], bb_max[2]):
                    ts.append(float(np.dot(
                        np.array([cx, cy, cz], dtype=np.float64) - loc,
                        direc)))
        t_min, t_max = min(ts), max(ts)
        length = t_max - t_min
        # Reject short cylinders (fillets, cap radii).  A real lead bar
        # is at least ~5 wire-radii long; helix profile cylinders are
        # never that long.
        if length < 5.0 * r:
            continue
        cap_a = loc + direc * t_min
        cap_b = loc + direc * t_max
        # Merge with an existing lead that shares this axis (same axis
        # line, same cap positions up to direction sign).  Two CYLINDER
        # faces with cap_a/cap_b swapped (because the unite split the
        # cylinder into 2 lateral half-faces with opposite axis_dir
        # parameterisation) get merged here.
        merged = False
        for existing in leads:
            same_cap = (
                (np.linalg.norm(cap_a - existing["cap_a"]) < 1e-6 and
                 np.linalg.norm(cap_b - existing["cap_b"]) < 1e-6) or
                (np.linalg.norm(cap_a - existing["cap_b"]) < 1e-6 and
                 np.linalg.norm(cap_b - existing["cap_a"]) < 1e-6))
            if same_cap:
                merged = True
                break
        if not merged:
            # Store axis as (cap_a + axis_dir_normed * t) with t in
            # [0, length].  This way _point_on_lead_axis can use t in
            # [0, length] without caring about the BRepAdaptor's loc
            # parameterisation (which may put loc at EITHER end).
            length_ab = float(np.linalg.norm(cap_b - cap_a))
            if length_ab > 1e-12:
                axis_dir_ab = (cap_b - cap_a) / length_ab
            else:
                axis_dir_ab = direc
            leads.append({"loc": loc, "dir": direc,  # raw BRepAdaptor values
                          "length": length_ab,
                          "radius": r,
                          "cap_a": cap_a, "cap_b": cap_b,
                          "axis_dir_ab": axis_dir_ab})
    return leads


def _classify_lead_caps(leads: list, chain_centroids_cad: np.ndarray):
    """For each lead, decide which cap is INNER (attached to helix) vs OUTER (tip).

    Inner = closer to a helix-only centroid (chain vertex that is NOT
    itself at a lead cap position).  Without the helix-only filter, the
    chain vertices at lead cap positions create distance=0 ties on both
    caps and argmin picks the wrong cap (whichever appears first), which
    flips the orientation.  See the W:/kubota/3turncoil.stp case where
    chain[0] = inner_cap AND chain[-1] = outer_cap of the same lead --
    both produce d=0 and argmin returns the first match.

    Mutates each lead dict in place, adding "inner_cap" and "outer_cap".
    """
    if chain_centroids_cad.shape[0] == 0:
        for lead in leads:
            lead["inner_cap"] = lead["cap_a"]
            lead["outer_cap"] = lead["cap_b"]
        return
    # Build helix-only set: chain vertices NOT at any lead cap position.
    # Tolerance = 0.5 * smallest lead radius (catches exact cap matches
    # but not nearby helix points).
    if leads:
        cap_match_tol = 0.5 * min(lead["radius"] for lead in leads)
    else:
        cap_match_tol = 1e-6
    helix_mask = np.ones(chain_centroids_cad.shape[0], dtype=bool)
    for lead in leads:
        for cap in (lead["cap_a"], lead["cap_b"]):
            d = np.linalg.norm(chain_centroids_cad - cap, axis=1)
            helix_mask &= (d > cap_match_tol)
    if not np.any(helix_mask):
        # Fallback: every vertex is a lead cap (degenerate)
        helix_only = chain_centroids_cad
    else:
        helix_only = chain_centroids_cad[helix_mask]
    for lead in leads:
        d_a = float(np.min(np.linalg.norm(helix_only - lead["cap_a"], axis=1)))
        d_b = float(np.min(np.linalg.norm(helix_only - lead["cap_b"], axis=1)))
        if d_a <= d_b:
            lead["inner_cap"] = lead["cap_a"]
            lead["outer_cap"] = lead["cap_b"]
        else:
            lead["inner_cap"] = lead["cap_b"]
            lead["outer_cap"] = lead["cap_a"]


def _point_on_lead_axis(p: np.ndarray, lead: dict,
                        perp_tol_factor: float = 1.5) -> bool:
    """Return True if point ``p`` lies within ``perp_tol_factor*radius`` of
    the lead's axis line, AND its along-axis parameter is within [0, length].

    Uses ``cap_a`` as origin and ``axis_dir_ab`` as direction, NOT the raw
    BRepAdaptor ``loc``/``dir`` which may parameterise the axis with loc
    at EITHER cap and dir pointing in either direction along the axis line.
    """
    v = p - lead["cap_a"]
    t = float(np.dot(v, lead["axis_dir_ab"]))
    if t < -1e-9 or t > lead["length"] + 1e-9:
        return False
    perp = v - t * lead["axis_dir_ab"]
    return float(np.linalg.norm(perp)) < perp_tol_factor * lead["radius"]


def _segment_on_lead_body(p0: np.ndarray, p1: np.ndarray, lead: dict) -> bool:
    """Return True if both segment endpoints lie on the same lead axis."""
    return _point_on_lead_axis(p0, lead) and _point_on_lead_axis(p1, lead)


def _augment_chain_with_lead_bars(centroids_cad: np.ndarray,
                                  radii_cad: np.ndarray,
                                  leads: list,
                                  median_seg_len_cad: float):
    """Repair an NN-chain that suffers from lead-related artifacts.

    The NN-chain on dedup'd cross-section CIRCLE centers handles the
    helix densely (sub-mm samples) but the straight lead bars are
    represented by only their cap circles -- producing 1-segment leads
    and, when a lead's body is entirely skipped, a "fake" segment
    jumping through air between two unrelated cap positions.

    This helper post-processes the chain to:
      1. Subdivide each long segment whose endpoints both lie on the
         same lead axis -- this IS the lead body, just under-sampled.
      2. Drop orphan trailing vertices that produce a fake air-jump
         (e.g. v[-1] is dangling at a lead's outer tip with no path
         through wire reaching it).
      3. Prepend / append missing lead bodies for leads whose body is
         not present anywhere in the chain.  An end of the chain that
         coincides with a lead's inner cap becomes the attachment
         point.

    Args:
        centroids_cad: (N, 3) NN-chain centroids in CAD units.
        radii_cad: (N,) per-centroid radius in CAD units.
        leads: list of lead dicts as returned by _detect_lead_bars_cad
            and classified by _classify_lead_caps.
        median_seg_len_cad: median segment length of the input chain
            (CAD units); used to set the long-segment threshold and
            the subdivision density.

    Returns:
        (new_centroids_cad, new_radii_cad): repaired chain.
    """
    if not leads or len(centroids_cad) < 2:
        return centroids_cad, radii_cad

    long_thresh = 5.0 * median_seg_len_cad
    n_stations_per_lead_body = lambda L: max(2, min(50, int(round(
        L / max(median_seg_len_cad, 1e-12)))))
    # Cap inclusion tolerance: how close a chain vertex must be to a
    # lead's cap to count as a hit (in CAD units).
    cap_match_tol = max(2.0 * median_seg_len_cad,
                        1.5 * max(lead["radius"] for lead in leads))

    centroids = list(centroids_cad)
    radii = list(radii_cad)

    # --- Step 1: subdivide long segments that lie on a single lead's body.
    # Walk the chain backwards (so indices remain valid as we splice).
    i = len(centroids) - 1
    while i >= 1:
        p0 = centroids[i - 1]
        p1 = centroids[i]
        seg_len = float(np.linalg.norm(p1 - p0))
        if seg_len > long_thresh:
            host_lead = None
            for lead in leads:
                if _segment_on_lead_body(p0, p1, lead):
                    host_lead = lead
                    break
            if host_lead is not None:
                # Subdivide into n new stations between p0 and p1.
                n_new = n_stations_per_lead_body(seg_len)
                for k in range(n_new - 1, 0, -1):
                    alpha = k / n_new
                    p_new = (1.0 - alpha) * p0 + alpha * p1
                    centroids.insert(i, p_new)
                    radii.insert(i, host_lead["radius"])
        i -= 1

    # --- Step 2: drop a dangling LAST vertex that produces a fake jump.
    # After step 1 the chain has all valid lead-body segs subdivided,
    # so any remaining long segment whose endpoints are not co-axial on
    # a single lead is a fake jump.  We only handle the boundary case
    # where the LAST segment is fake (the only one observed in
    # 3turnCoil_work_coil.step's broken chain).  Interior fake jumps
    # would require a graph-split which is out of scope here.
    if len(centroids) >= 2:
        last_seg_len = float(np.linalg.norm(centroids[-1] - centroids[-2]))
        if last_seg_len > long_thresh:
            on_a_lead = any(
                _segment_on_lead_body(centroids[-2], centroids[-1], lead)
                for lead in leads)
            if not on_a_lead:
                # Drop the orphan tip.  The chain now ends at what used
                # to be v[-2], which (for 3turnCoil) is a lead outer
                # tip -- a valid port location.
                centroids.pop()
                radii.pop()
    # Symmetric check at the chain start.
    if len(centroids) >= 2:
        first_seg_len = float(np.linalg.norm(centroids[1] - centroids[0]))
        if first_seg_len > long_thresh:
            on_a_lead = any(
                _segment_on_lead_body(centroids[0], centroids[1], lead)
                for lead in leads)
            if not on_a_lead:
                centroids.pop(0)
                radii.pop(0)

    # --- Step 3: prepend/append a missing lead body to either chain end.
    # For each lead, decide whether its body is already in the chain.
    # We say "in the chain" if some segment has both endpoints on this
    # lead's axis.  If not, the lead is missing; we attach it at the
    # chain end whose vertex matches its inner_cap.
    centroids_arr = np.array(centroids)
    for lead in leads:
        # Already represented?
        already = False
        for j in range(len(centroids) - 1):
            if _segment_on_lead_body(centroids[j], centroids[j + 1], lead):
                already = True
                break
        if already:
            continue
        # Find which chain endpoint matches this lead's inner cap.
        d_start = float(np.linalg.norm(centroids[0] - lead["inner_cap"]))
        d_end = float(np.linalg.norm(centroids[-1] - lead["inner_cap"]))
        n_new = n_stations_per_lead_body(lead["length"])
        if d_start <= d_end and d_start <= cap_match_tol:
            # Prepend: outer_tip -> stations -> inner_cap (= existing chain start).
            new_stations = []
            new_radii = []
            for k in range(n_new + 1):  # include both endpoints
                alpha = k / n_new
                p_new = (1.0 - alpha) * lead["outer_cap"] + alpha * lead["inner_cap"]
                new_stations.append(p_new)
                new_radii.append(lead["radius"])
            # Drop duplicate inner_cap -- existing chain[0] is the inner end.
            new_stations = new_stations[:-1]
            new_radii = new_radii[:-1]
            centroids = new_stations + centroids
            radii = new_radii + radii
        elif d_end < d_start and d_end <= cap_match_tol:
            # Append: inner_cap (= existing chain end) -> stations -> outer_tip.
            new_stations = []
            new_radii = []
            for k in range(n_new + 1):
                alpha = k / n_new
                p_new = (1.0 - alpha) * lead["inner_cap"] + alpha * lead["outer_cap"]
                new_stations.append(p_new)
                new_radii.append(lead["radius"])
            new_stations = new_stations[1:]  # drop duplicate inner_cap
            new_radii = new_radii[1:]
            centroids = centroids + new_stations
            radii = radii + new_radii

    return np.array(centroids), np.array(radii)


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

    Spine policy (2026-05-02 unified; 2026-07-29 classified, no
    fallback):
      * SINGLE-TURN coil -> ``coil_topology.extract_coil_topology`` +
        ``generate_spine`` -- OPEN/CLOSED-aware.  CLOSED full torus
        (no caps) gets a full 360deg spine; OPEN gapped torus gets
        cap_a -> cap_b along the LONG arc.  R_spine is refined from
        the actual mean of the cross-section centroids around the
        rotation axis (more accurate than the bbox 0.85 default).
        NOTE: ``generate_spine`` supports only the z rotation axis. If
        non-z geometry reaches this tier, the caller's bbox, inside-
        solid, and near-surface checks reject the mis-seeded topology;
        there is no automatic alternative extraction after that
        positive classification.
      * MULTI-TURN HELIX (z-extent of the cross-section centroids
        > 2x the wire radius) -> ``_chain_centroids_nn_index`` with
        tangent continuity, which walks the spiral correctly; the
        planar topology arc would collapse N turns into one
        (L drops ~N^2).
    This is a CLASSIFIED dispatch on measured geometry, not a
    try/except chain -- extraction failures propagate.

    Returns ``None`` (caller falls through to the constant equivalent-
    circle path) when the solid does not have a clean population of
    consistent-radius circle edges.
    """
    from radia._b3d_shim import GeomType
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
    # Vectorised O(N*K) greedy dedup -- on the 3turncoil sample this drops
    # the Python-level any()/norm() count from ~199k to ~5k, recovering
    # ~0.9 s of the 9.4 s filaments_from_step cold run (2026-05-19 profile).
    dedup_tol = 0.1 * median_r
    centers_arr = np.asarray([d["center"] for d in consistent], dtype=float)
    kept_idx = []
    kept_arr = np.empty((0, 3), dtype=float)
    for i, c in enumerate(centers_arr):
        if kept_arr.shape[0] == 0:
            kept_idx.append(i)
            kept_arr = c[None, :].copy()
            continue
        # Vectorised: one numpy norm call per new candidate instead of
        # a Python-loop generator over K previously-kept centres.
        d_min = np.linalg.norm(kept_arr - c, axis=1).min()
        if d_min >= dedup_tol:
            kept_idx.append(i)
            kept_arr = np.vstack([kept_arr, c[None, :]])
    kept = [consistent[i] for i in kept_idx]
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

    # Spine: CLASSIFIED single dispatch on measured geometry
    # (No-Fallbacks, 2026-07-29 -- the old ``except Exception ->
    # NN-chain`` swallow routed a broken topology extraction into a
    # silently different algorithm; extraction failures now propagate).
    if is_multi_turn_helix:
        # Multi-turn helix: NN-chain with tangent continuity walks the
        # spiral correctly; the planar topology arc would collapse N
        # turns into one (L drops ~N^2).  Handles OPEN reliably
        # (endpoint score detects the caps).
        perm = _chain_centroids_nn_index(centers_cad)
        ordered = [kept[i] for i in perm]
        centroids_m = centers_cad[perm] / cad_units_per_meter
        radii_m = np.array([d["radius"] for d in ordered]) / cad_units_per_meter
    else:
        from radia.coil_topology import (
            extract_coil_topology as _extract_topo,
            generate_spine as _gen_spine,
            CoilTopology as _CoilTopology,
        )
        topo = _extract_topo(solid)
        # Refine R_spine from the actual cross-section centroid
        # distance to the rotation axis (axis = z).  The bbox-based
        # 0.85 * R_outer default in coil_topology is conservative;
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

    # Repair lead-related artifacts in the NN-chain:
    # - Single-segment lead bodies (no intermediate samples)
    # - Fake air-jumps between unrelated lead tips
    # - Lead bodies missing from the chain entirely
    # See _augment_chain_with_lead_bars docstring for details.
    leads_cad = _detect_lead_bars_cad(solid, median_r)
    if leads_cad:
        centroids_cad_chain = centroids_m * cad_units_per_meter
        radii_cad_chain = radii_m * cad_units_per_meter
        _classify_lead_caps(leads_cad, centroids_cad_chain)
        seg_lens_cad = np.linalg.norm(
            np.diff(centroids_cad_chain, axis=0), axis=1)
        # Use median over the SHORT segments only (lead jumps would
        # poison the median otherwise).  Filter by length < 5x mean.
        if len(seg_lens_cad) > 0:
            mean_seg = float(np.mean(seg_lens_cad))
            short_segs = seg_lens_cad[seg_lens_cad < 5.0 * mean_seg]
            median_seg_cad = float(np.median(short_segs)) if len(short_segs) > 0 \
                else float(np.median(seg_lens_cad))
        else:
            median_seg_cad = 1.0
        centroids_cad_chain, radii_cad_chain = _augment_chain_with_lead_bars(
            centroids_cad_chain, radii_cad_chain,
            leads_cad, median_seg_cad)
        centroids_m = centroids_cad_chain / cad_units_per_meter
        radii_m = radii_cad_chain / cad_units_per_meter

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
    from radia.peec_bundle import build_bundle_solver
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
    coil with a clear consistent-radius circle population.  This is
    a predicate-style negative ("the cap-circle pattern does not
    match"); the classification dispatch in
    ``extract_centerline_from_step`` then tries the next predicate.
    Returning None is NOT a fallback -- predicates are positive
    matches by design.

    Returns ``(centers_cad, median_radius_cad)`` on success:
        centers_cad: list of (3,) np arrays in raw CAD units
        median_radius_cad: float in raw CAD units
    """
    from radia._b3d_shim import GeomType
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


def _walked_stations_to_path_cad(res, ng_solid, predicate_name: str):
    """Walker ``CenterlineResult`` -> ``(path_cad, widths_cad)``.

    The walked stations are kept at their adaptive density: uniform
    steps on straight or smooth stretches and halved steps through
    bends.  Only the omitted wrap edge of a closed path is normalized
    to roughly one local step so the PEEC excitation port cannot become
    vanishingly small or span several cells.  Per-vertex bend angles
    stay small through bends as a result, which keeps the downstream
    ``_check_spine_no_singular_corner`` semantics sound (a bend cluster
    is short-but-shallow; a corner-turn jump is sharp-but-long; neither
    trips the sharp-AND-short condition).

    Widths: per-segment equivalent-square side.  The RELATIVE
    structure comes from the walker's per-station slab areas (clamped
    to [0.5, 2] x median so an oblique bend cut or a fused boss cannot
    leak a wildly wrong local value), but the ABSOLUTE scale is pinned
    to the exact global truth ``ng_solid.mass / walked_length``:
    OCC thin-slab booleans are numerically fragile on BSPLINE-loft
    laterals (measured 2026-07-29 on keiko_outsideline: a uniform
    ~0.65x area deficit at the working slab thickness, while the
    NEIGHBOURING thicknesses return EMPTY intersections and a
    build123d planar section at the same point gives the true
    26.6 mm^2).  ``.mass`` of the full solid is reliable, so
    mass / length recovers the exact arclength-mean section area.
    """
    pts = np.asarray(res.polyline, dtype=float).copy()
    if pts.ndim != 2 or pts.shape[1] != 3:
        raise ValueError(
            f"extract_centerline_from_step({predicate_name}): walker "
            f"polyline must have shape (N, 3), got {pts.shape}.")
    if not np.all(np.isfinite(pts)):
        raise ValueError(
            f"extract_centerline_from_step({predicate_name}): walker "
            "polyline contains non-finite coordinates.")
    if pts.shape[0] < 4:
        raise ValueError(
            f"extract_centerline_from_step({predicate_name}): the walk "
            f"produced only {pts.shape[0]} stations -- too coarse for "
            "a spine.")
    if res.areas is None:
        raise ValueError(
            f"extract_centerline_from_step({predicate_name}): walker "
            "result carries no per-station areas.")
    areas = np.asarray(res.areas, dtype=float).copy()
    if areas.shape != (pts.shape[0],):
        raise ValueError(
            f"extract_centerline_from_step({predicate_name}): walker "
            "must provide one cross-section area per station; got "
            f"{areas.shape} for {pts.shape[0]} stations.")
    good = np.isfinite(areas) & (areas > 0.0)
    if not np.any(good):
        raise ValueError(
            f"extract_centerline_from_step({predicate_name}): the walk "
            "recovered no positive cross-section areas; cannot derive "
            "the conductor cross-section.")
    med = float(np.median(areas[good]))
    areas = np.where(good, areas, med)
    seg_len = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    if not np.all(np.isfinite(seg_len)) or np.any(seg_len <= 0.0):
        raise ValueError(
            f"extract_centerline_from_step({predicate_name}): walker "
            "polyline contains a non-finite or zero-length segment.")

    # A closed PEEC path deliberately leaves the final wrap edge out of
    # the conductor segments and places the excitation port across it.
    # The walker stops anywhere inside close_tol, so that edge can range
    # from almost zero to about two regular steps.  Preserve every
    # interior adaptive station, but trim a near-duplicate closure point
    # or split an oversized wrap chord so the omitted port edge stays at
    # a stable local-cell scale.
    if res.closed:
        step_ref = float(np.median(seg_len))
        close_len = float(np.linalg.norm(pts[0] - pts[-1]))
        while (pts.shape[0] > 4
               and close_len < 0.5 * step_ref):
            pts = pts[:-1]
            areas = areas[:-1]
            seg_len = np.linalg.norm(np.diff(pts, axis=0), axis=1)
            step_ref = float(np.median(seg_len))
            close_len = float(np.linalg.norm(pts[0] - pts[-1]))
        if close_len <= 0.0 or not math.isfinite(close_len):
            raise ValueError(
                f"extract_centerline_from_step({predicate_name}): "
                "closed walker path has a zero or non-finite port gap.")
        if close_len < 0.5 * step_ref:
            raise ValueError(
                f"extract_centerline_from_step({predicate_name}): "
                "closed walker path has too few distinct stations to "
                "form a stable PEEC port gap.")
        if close_len > 1.5 * step_ref:
            n_parts = max(2, int(math.ceil(close_len / step_ref)))
            frac = np.arange(1, n_parts, dtype=float) / n_parts
            inserts = (pts[-1][None, :]
                       + frac[:, None] * (pts[0] - pts[-1])[None, :])
            area_inserts = areas[-1] + frac * (areas[0] - areas[-1])
            pts = np.vstack([pts, inserts])
            areas = np.concatenate([areas, area_inserts])
            seg_len = np.linalg.norm(np.diff(pts, axis=0), axis=1)

    areas_clamped = np.clip(areas, 0.5 * med, 2.0 * med)
    seg_areas = 0.5 * (areas_clamped[:-1] + areas_clamped[1:])
    full_len = seg_len
    full_areas = seg_areas
    if res.closed:
        close_len = float(np.linalg.norm(pts[0] - pts[-1]))
        close_area = 0.5 * (areas_clamped[-1] + areas_clamped[0])
        full_len = np.append(seg_len, close_len)
        full_areas = np.append(seg_areas, close_area)
    length = float(np.sum(full_len))
    if length <= 0.0:
        raise ValueError(
            f"extract_centerline_from_step({predicate_name}): walked "
            "centerline has zero length.")
    mass = float(ng_solid.mass)
    if not math.isfinite(mass) or mass <= 0.0:
        raise ValueError(
            f"extract_centerline_from_step({predicate_name}): conductor "
            f"solid has non-positive or non-finite volume {mass!r}.")
    mean_area_true = mass / length
    walker_mean = float(np.sum(full_areas * full_len) / length)
    if not math.isfinite(walker_mean) or walker_mean <= 0.0:
        raise ValueError(
            f"extract_centerline_from_step({predicate_name}): walker "
            "mean section area is non-positive.")
    k = mean_area_true / walker_mean
    if not (0.4 <= k <= 2.5):
        raise ValueError(
            f"extract_centerline_from_step({predicate_name}): slab-"
            f"area absolute-scale correction k = {k:.3g} outside "
            "[0.4, 2.5] -- the walked length and the solid volume are "
            "inconsistent, i.e. the walk did NOT cover the whole "
            "conductor (e.g. a closed sub-loop short-circuiting a "
            "lead-pair junction), the solid is multi-body, or a "
            "boolean broke.  (Empirical bounds: the OCC thin-slab "
            "boolean under-measures BSPLINE-loft sections by up to "
            "~0.65x -> k ~1.5; a modest non-wire feature such as a "
            "boss adds a little more.)  Regenerate the STEP or switch "
            "to --coil-solver bem-a --coil-vol <pre-meshed.vol>.")
    widths_cad = np.sqrt(seg_areas * k)
    return pts, widths_cad


def _centerline_from_open_spine(solid, n_segments: int,
                                 cad_units_per_meter: float):
    """Single-loop coil centerline: longest open edge (spine) sampling.

    For swept / bent coils (gapped torus, single-turn helix), the STEP
    solid has an open boundary edge that traces the coil spine.  We
    pick the longest OPEN (non-closed-loop) edge — closed circle
    edges are cross-section boundaries, not spines.

    Raises ValueError if no open spine edge exists (e.g. loft of
    circles: all edges are closed cross-section circles).  This is
    a hard error -- the classification dispatch in
    ``extract_centerline_from_step`` chooses this extractor based on
    a positive predicate (``topo.is_open``); a failure here means
    the dispatcher's predicate was wrong, not that "the next path
    should be tried" (CLAUDE.md "No Fallbacks").
    """
    from radia._b3d_shim import section, Plane, Vector
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
    # Deterministic edge selection (CLAUDE.md "Cubit Probe, Don't Guess"):
    # multiple equal-length BSPLINE rim halves on lofted-cap geometries
    # produce ties on `key=lambda e: e.length`.  OCC / Cubit listing
    # order is non-deterministic across versions, so we break the tie
    # on (length, centroid_x, centroid_y, centroid_z) lexicographic
    # sort.  This keeps the same edge across runs / machine versions.
    def _spine_sort_key(e):
        try:
            c = e.center()
            return (-e.length, float(c.X), float(c.Y), float(c.Z))
        except Exception:
            # Curve types without center() (analytic curves): fall back
            # to length-only.  Such curves rarely tie on length.
            return (-e.length, 0.0, 0.0, 0.0)
    spine = sorted(open_edges, key=_spine_sort_key)[0]

    # Adaptive sampling (v4.53.0, keiko's request 2026-05-16): the
    # caller's ``n_segments`` is treated as an UPPER bound only.  The
    # actual segment count is chosen so that the resulting spine
    # segment length is at least ``1.10 * wire_radius_estimate``,
    # which prevents the downstream ``_check_spine_no_singular_corner``
    # check from firing on benign smooth bends that happen to be
    # over-sampled.  Without this, a 208 mm spine with caller-default
    # ``n_segments = 100`` gives 2.08 mm segments -- with a 2.9 mm
    # wire radius the ratio is 0.72 (below the 1.0 threshold), and any
    # > 60 deg bend along the spine (e.g. the lead-arc junction on
    # keiko's 1turn_coil_loft) trips the singular-corner check even
    # though the geometry is physically fine.  Adaptive sampling makes
    # the check semantically correct: it fires on TRUE corners (sharp
    # AND densely-sampled by user request) but not on smooth bends
    # accidentally over-sampled.
    #
    # Wire radius is estimated by a CHEAP midpoint section before the
    # full sampling loop runs.
    spine_length_cad = float(spine.length)
    from radia._b3d_shim import section as _bd_section
    from radia._b3d_shim import Plane as _bd_Plane
    from radia._b3d_shim import Vector as _bd_Vector
    mid_p = spine @ 0.5
    mid_tangent = spine @ 0.51 - mid_p
    tn = math.sqrt(mid_tangent.X ** 2 + mid_tangent.Y ** 2
                    + mid_tangent.Z ** 2)
    if tn <= 1e-12:
        raise ValueError(
            "_centerline_from_open_spine: spine midpoint tangent is "
            "degenerate (|spine@0.51 - spine@0.50| < 1e-12).  The "
            "selected open edge is not a usable spine -- regenerate "
            "the STEP with a smooth single-piece BSPLINE lateral so "
            "Predicate 1 (UV-map sampling) handles it directly."
        )
    probe_plane = _bd_Plane(
        origin=_bd_Vector(float(mid_p.X), float(mid_p.Y),
                            float(mid_p.Z)),
        z_dir=_bd_Vector(float(mid_tangent.X / tn),
                           float(mid_tangent.Y / tn),
                           float(mid_tangent.Z / tn)))
    probe_section = _bd_section(solid, section_by=probe_plane)
    probe_faces = (probe_section.faces() if probe_section is not None
                    else [])
    if not probe_faces:
        raise ValueError(
            "_centerline_from_open_spine: adaptive-resampling probe "
            "at spine midpoint produced no cross-section face.  The "
            "open-spine assumption breaks here -- regenerate the "
            "STEP with a smooth single-piece BSPLINE lateral so "
            "Predicate 1 (UV-map sampling) handles it directly."
        )
    best_probe = min(
        probe_faces,
        key=lambda f: ((f.center() - _bd_Vector(
            float(mid_p.X), float(mid_p.Y), float(mid_p.Z))
        ).length))
    cross_area_cad = float(best_probe.area)
    wire_r_cad = math.sqrt(cross_area_cad / math.pi)
    min_seg_cad = 1.10 * wire_r_cad
    if min_seg_cad > 0:
        max_segments_by_density = max(
            3, int(spine_length_cad / min_seg_cad))
        n_segments = min(n_segments, max_segments_by_density)

    spine_pts = np.zeros((n_segments + 1, 3), dtype=np.float64)
    for i in range(n_segments + 1):
        t = i / n_segments
        p = spine @ t
        spine_pts[i] = [p.X, p.Y, p.Z]

    # Corner densification (v4.54.0): uniform sampling gives smooth
    # bend angles per segment, but at the lead-arc junctions of "arc
    # + leads" coils a single bend can absorb 60+ deg of curvature in
    # one station.  This makes the cross-section frame rotate ~60 deg
    # in one step, which (a) visually shows as filament "bunching"
    # when GMSH viewers project the tilted cross-section plane onto an
    # axis-aligned view, and (b) makes the parallel-bundle filament
    # path very piecewise-linear at the corner.  Insert intermediate
    # stations along the spine arc-length so any one segment-to-segment
    # bend is at most ``max_bend_per_step_deg`` (default 20 deg).
    # The insertion runs along the OCC spine curve (smoothly param-
    # eterized), so the inserted points lie on the real spine, not on
    # a polyline interpolation.
    max_bend_per_step_deg = 20.0
    spine_pts = _densify_at_corners(
        spine, spine_pts, max_bend_per_step_deg)
    n_segments = len(spine_pts) - 1

    tangents = path_tangents(spine_pts)
    midpoints = 0.5 * (spine_pts[:-1] + spine_pts[1:])

    centerline = np.zeros((n_segments + 1, 3), dtype=np.float64)
    widths_cad = np.zeros(n_segments, dtype=np.float64)
    heights_cad = np.zeros(n_segments, dtype=np.float64)
    centerline[0] = spine_pts[0]
    centerline[-1] = spine_pts[-1]

    # Per-station OCC sectioning.  Failures (degenerate plane, no face,
    # OCC native exception) are HARD -- the open-spine path is built on
    # the assumption that the longest open lateral edge approximates
    # the spine, and a failed section means that assumption is wrong
    # on this particular STEP (likely a tangent-discontinuity or a
    # self-intersecting region).  Per CLAUDE.md "No Fallbacks - Fail
    # Fast, Fail Loud", we raise instead of silently using the previous
    # station's width (the v4.48.x cascade-era residue that produced
    # widths_cad[0] = 0.0 when station 0 failed -> downstream NaN).
    for i in range(n_segments):
        c = midpoints[i]
        t = tangents[i]
        origin = Vector(float(c[0]), float(c[1]), float(c[2]))
        z_dir = Vector(float(t[0]), float(t[1]), float(t[2]))
        sec_plane = Plane(origin=origin, z_dir=z_dir)
        try:
            cross = section(solid, section_by=sec_plane)
        except Exception as exc:
            raise ValueError(
                f"_centerline_from_open_spine: section() failed at "
                f"station {i}/{n_segments} "
                f"(midpoint=({c[0]:.4g},{c[1]:.4g},{c[2]:.4g}), "
                f"tangent=({t[0]:.3f},{t[1]:.3f},{t[2]:.3f})): "
                f"{type(exc).__name__}: {exc}.  "
                f"This usually means the longest open lateral edge "
                f"is not a clean spine on this geometry (tangent "
                f"discontinuity, self-intersection, or numerically "
                f"degenerate plane).  Regenerate the STEP with a "
                f"smooth single-piece BSPLINE lateral so Predicate 1 "
                f"(UV-map sampling) handles it directly."
            ) from exc

        if cross is None or len(cross.faces()) == 0:
            raise ValueError(
                f"_centerline_from_open_spine: section() at station "
                f"{i}/{n_segments} produced no face "
                f"(midpoint=({c[0]:.4g},{c[1]:.4g},{c[2]:.4g})).  "
                f"The sectioning plane likely missed the solid -- "
                f"the open-spine assumption breaks here.  Regenerate "
                f"the STEP with a smooth single-piece BSPLINE lateral "
                f"so Predicate 1 (UV-map sampling) handles it directly."
            )

        best = min(cross.faces(),
                   key=lambda f: (f.center() - origin).length)
        bc = best.center()
        centerline[i] = [bc.X, bc.Y, bc.Z]
        if i == n_segments - 1:
            centerline[i + 1] = spine_pts[i + 1]
        side = math.sqrt(best.area)
        widths_cad[i] = side
        heights_cad[i] = side

    path_cad = np.zeros((n_segments + 1, 3), dtype=np.float64)
    path_cad[0] = centerline[0]
    for i in range(n_segments - 1):
        path_cad[i + 1] = 0.5 * (centerline[i] + centerline[i + 1])
    path_cad[-1] = centerline[-1]

    # Cap-centroid endpoint anchoring (v4.55.0, keiko viz report
    # 2026-05-16): the "longest open edge" spine traces the
    # conductor's LATERAL RIM (e.g. z = +-wire_radius on a flat
    # coil), not the centroid axis.  Interior path_cad points come
    # from midpoint sectioning -> face centroids = correct centroid
    # path, but the ENDPOINTS path_cad[0] and path_cad[-1] are pinned
    # to spine_pts[0]/[-1] (rim endpoints near the cap edge), causing
    # a 41 deg cap-direction kink + ~48% |I| spread.  Fix: replace
    # the rim endpoints with the cap-face centroids from
    # coil_topology.  Predicate 4 (OPEN) routes here precisely
    # because cap detection succeeded, so cap_a/cap_b MUST be
    # available.  If they are not, the caller's classification is
    # inconsistent -- raise loudly per CLAUDE.md "No Fallbacks".
    from radia.coil_topology import extract_coil_topology as _ext_topo
    _topo = _ext_topo(solid)
    if not (_topo.is_open and _topo.cap_a is not None
            and _topo.cap_b is not None):
        raise ValueError(
            "_centerline_from_open_spine: coil_topology.extract_coil_topology "
            "returned is_open={} cap_a={} cap_b={} -- inconsistent with "
            "Predicate 4 (OPEN) classification that routed here.  This "
            "indicates a regression in coil_topology cap detection; "
            "regenerate the STEP with planar end-caps or file a bug."
            .format(_topo.is_open, _topo.cap_a is not None,
                    _topo.cap_b is not None))
    ca = _topo.cap_a.center()
    cb = _topo.cap_b.center()
    cap_a_xyz = np.array([float(ca.X), float(ca.Y), float(ca.Z)])
    cap_b_xyz = np.array([float(cb.X), float(cb.Y), float(cb.Z)])
    d0_a = float(np.linalg.norm(path_cad[0] - cap_a_xyz))
    d0_b = float(np.linalg.norm(path_cad[0] - cap_b_xyz))
    if d0_a <= d0_b:
        path_cad[0] = cap_a_xyz
        path_cad[-1] = cap_b_xyz
    else:
        path_cad[0] = cap_b_xyz
        path_cad[-1] = cap_a_xyz

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
    from radia._b3d_shim import GeomType
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
    """CLOSED coil spine by walking the actual solid (slab marching).

    v4.95.x: replaces the retired ``coil_topology.generate_spine``
    bbox-circle map, which assumed BOTH a z rotation axis AND a
    circular planar spine at 0.85 * bbox outer radius.  Any closed
    NON-circular or non-z-axis coil received a geometrically-unrelated
    circle -- observed on B_rect_sweep.step / J_boss_fused.step (a
    127.5 mm circle for a 103 x 214 mm racetrack, caught only
    afterwards by the bbox checks) -- and even a TRUE circular torus
    got the wrong radius (0.85 * (R + r) instead of R).  The walker
    marches cross-sections of the true solid, so ANY orientation and
    ANY closed spine shape is traced.

    Restricted to CLOSED coils by the ``extract_centerline_from_step``
    dispatch.  If the walk does not close on itself this raises:
    either the coil is actually OPEN and cap-face detection missed the
    caps (regenerate the STEP with clean planar end caps so
    Predicate 4 routes), or the walker halted at a sharp spine corner.

    ``_walked_stations_to_path_cad`` preserves the walker's adaptive
    interior density, normalizes only the omitted PEEC port edge, and
    derives per-segment widths from the measured areas.  ``n_segments``
    is unused (parameter kept for dispatch signature parity).
    """
    from radia.coil_from_step import extract_centerline, _axis_agnostic_seed

    ng_solid = _bd_shape_to_netgen_solid(solid)
    seed = _axis_agnostic_seed(ng_solid)
    res = extract_centerline(ng_solid, start_hint=seed, verbose=False)

    if not res.closed:
        pts_arr = np.asarray(res.polyline, dtype=float)
        span = ((pts_arr.max(axis=0) - pts_arr.min(axis=0))
                if pts_arr.size else None)
        raise ValueError(
            "extract_centerline_from_step(topology_spine): the coil was "
            "classified CLOSED (no cap faces detected) but the "
            f"centerline walk did not close on itself ({pts_arr.shape[0]} "
            f"stations, walked span {span} CAD units).  Either the coil "
            "is actually OPEN and cap-face detection missed the caps "
            "(regenerate the STEP with clean planar end caps so "
            "Predicate 4 (open_spine) routes), or the walker halted at "
            "a sharp spine corner (fillet the spine corners, or switch "
            "to --coil-solver bem-a --coil-vol <pre-meshed.vol> which "
            "bypasses spine extraction entirely).")

    pts, widths_cad = _walked_stations_to_path_cad(
        res, ng_solid, "topology_spine")
    scale = 1.0 / cad_units_per_meter
    widths_m = widths_cad * scale
    return pts * scale, widths_m, widths_m.copy()


_MIN_PLAUSIBLE_COIL_EXTENT_M = 1e-5
_MAX_PLAUSIBLE_COIL_EXTENT_M = 5.0


def _normalize_cad_units_per_meter(cad_units_per_meter, source_tag):
    """Return a finite positive CAD-units-per-metre value as ``float``."""
    try:
        scale = float(cad_units_per_meter)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"{source_tag}: cad_units_per_meter must be finite and > 0, "
            f"got {cad_units_per_meter!r}.") from exc
    if not math.isfinite(scale) or scale <= 0:
        raise ValueError(
            f"{source_tag}: cad_units_per_meter must be finite and > 0, got "
            f"{cad_units_per_meter!r}.")
    return scale


def _check_solid_extent_plausible(solid, cad_units_per_meter, source_tag):
    """Fail-loud guard against a units mismatch (CLAUDE.md "Fail Fast,
    Fail Loud" -- no silent auto-correction).

    The extractor scales CAD coordinates to METRES via
    ``cad_units_per_meter``.  If the caller leaves the default 1.0
    (STEP assumed in metres) but the STEP is actually in millimetres,
    every downstream length is 1000x too large and the inductance is
    ~1000x wrong SILENTLY (observed on the I_wrong_units variant:
    L=114109 nH instead of ~114 nH).  This does NOT guess the unit; it
    raises when the implied physical coil size falls outside the
    plausible range [0.01 mm, 5 m] and names the cad_units_per_meter
    that WOULD make it plausible.
    """
    scale = _normalize_cad_units_per_meter(cad_units_per_meter, source_tag)
    bb = solid.bounding_box()
    extents_cad = (
        float(bb.max.X - bb.min.X),
        float(bb.max.Y - bb.min.Y),
        float(bb.max.Z - bb.min.Z),
    )
    ext_cad = max(extents_cad)
    if not all(math.isfinite(value) for value in extents_cad) or ext_cad <= 0:
        raise ValueError(
            f"{source_tag}: STEP bounding-box extents must be finite and "
            f"non-zero, got {extents_cad!r} CAD units.")
    ext_m = ext_cad / scale
    # Plausible coil size window.  The MAIN mistake this guards is
    # mm-read-as-metres, which inflates the coil ~1000x (upper bound).
    # The lower bound stays very permissive (10 um) so genuinely small
    # synthetic/test coils are not rejected.
    if _MIN_PLAUSIBLE_COIL_EXTENT_M <= ext_m <= _MAX_PLAUSIBLE_COIL_EXTENT_M:
        return scale  # plausible coil size -- pass and normalize
    suggest = None
    for cand in (1.0, 1e3, 1e-3, 1e2, 1e6):
        if (_MIN_PLAUSIBLE_COIL_EXTENT_M
                <= ext_cad / cand
                <= _MAX_PLAUSIBLE_COIL_EXTENT_M):
            suggest = cand
            break
    raise ValueError(
        f"{source_tag}: implied coil extent {ext_m:.3g} m is implausible "
        f"(expected 0.01 mm .. 5 m).  The STEP bbox spans {ext_cad:.3g} "
        f"CAD units and cad_units_per_meter={scale:g}.  "
        f"This is almost certainly a UNIT mismatch: a millimetre STEP "
        f"read as metres yields a 1000x-too-large coil and a ~1000x-wrong "
        f"inductance." +
        (f"  Pass cad_units_per_meter={suggest:g} for this STEP."
         if suggest else "  Verify the STEP's export units."))


def extract_centerline_from_step(step_path: str,
                                 n_segments: int = 100,
                                 cad_units_per_meter: float = 1.0):
    """Auto-extract coil centerline + cross-sections from a STEP file.

    ARCHITECTURE (v4.95.x): ONE geometric marching engine + exact
    CAD-feature fast paths.  The walking-plane march
    (``coil_from_step.extract_centerline``: axis-agnostic seed,
    bidirectional, corner-turning) is the geometry-marching engine --
    Predicate 5 (CLOSED) delegates to it, as does the default
    ``filaments_from_step`` path.  Predicates 1-4 are NOT alternative
    engines: they are positive matches on exact CAD features that
    specific authoring workflows leave in the STEP, read without
    marching:

    * **P1 loft-of-profiles** (>= 5 consistent-area planar station
      faces, e.g. Cubit station-loft): chain the exact station
      centroids.  Robust for tight pancakes where adjacent turns are
      closer than a march step.
    * **P2 united multi-turn pancake** (>= 5 consistent-radius CIRCLE
      edges surviving a boolean unite): chain the exact circle
      centers.  A marching engine risks snapping across turns when
      the inter-turn pitch is below ~4 march steps, so the exact
      reader stays canonical for this class.
    * **P3 revolution sweep** (TORUS/CYLINDER/CONE/REVOLUTION lateral
      + planar caps): axis + major radius + sweep angle extracted
      ANALYTICALLY -- exact, no sampling at all.
    * **P4 OPEN coil** (cap faces detected): sample the longest open
      LATERAL RIM EDGE (an exact CAD feature), section at midpoints,
      pin endpoints to the cap centroids.  The rim edge is IMMUNE to
      lead-pair junction proximity, where a marching engine
      short-circuits the junction into a closed sub-loop (measured
      2026-07-29 on D_spline_sharp: the walk 'closed' across leads
      ~2 wire-diameters apart while the rim extractor traces the true
      open path) -- so the exact reader stays canonical here too.

    Everything else marches:

    * **P5 CLOSED loop** (no caps): walk; the walk must close on
      itself, span the solid, and volume-reconcile (k-guard).

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
    from radia._b3d_shim import import_step, GeomType, Compound

    solid = import_step(step_path)

    # Multi-solid entry guard (v4.49.0): single-coil PEEC assumes one
    # solid in the STEP.  build123d's import_step silently flattens
    # a multi-solid Compound and `solid.faces()` enumerates ALL
    # solids' faces -- the predicates then mis-classify (e.g. cap
    # detection sees 2N planes and Predicate 1 dominance plummets).
    # CLAUDE.md "No Fallbacks - Fail Fast, Fail Loud": raise loud.
    if isinstance(solid, Compound):
        # Compound may wrap one or many solids; only raise when > 1
        sub_solids = list(solid.solids())
        if len(sub_solids) > 1:
            bb = [(s.bounding_box().min, s.bounding_box().max)
                   for s in sub_solids]
            summary = "; ".join(
                f"#{i}: x=[{m.X * 1e3:+.1f},{x.X * 1e3:+.1f}] "
                f"y=[{m.Y * 1e3:+.1f},{x.Y * 1e3:+.1f}] mm"
                for i, (m, x) in enumerate(bb[:5]))
            raise ValueError(
                f"extract_centerline_from_step: STEP contains "
                f"{len(sub_solids)} solids, expected exactly 1.  "
                f"Bbox summary: {summary}"
                f"{'; ...' if len(bb) > 5 else ''}.  "
                f"Single-coil PEEC handles one solid per file -- split "
                f"the STEP into per-solid files and call this function "
                f"on each, OR boolean-unite the conductor solids in "
                f"CAD before export so the result is a single watertight "
                f"solid.")
        if len(sub_solids) == 1:
            solid = sub_solids[0]
        # else: solid stays the Compound; downstream predicates will
        # likely raise on empty faces -- that is fine, hard error.

    # Units sanity: fail-loud on a mm-as-metres mismatch (silent
    # ~1000x-wrong L otherwise; No-Fallbacks / Fail-Loud).
    cad_units_per_meter = _check_solid_extent_plausible(
        solid, cad_units_per_meter, "extract_centerline_from_step")

    # Classification-based dispatch (No-Fallback policy, CLAUDE.md
    # "No Fallbacks - Fail Fast, Fail Loud"): inspect the solid's
    # features upfront and route to exactly ONE specialized centerline
    # extractor.  No try/except cascades; each predicate is positive-
    # match only.  If a predicate misclassifies, the failure surfaces
    # as a hard error from the specialised extractor rather than being
    # silently swallowed by "try the next path".

    # Each predicate's extractor returns (path_m, widths_m, heights_m).
    # Per CLAUDE.md "No Fallbacks", AFTER the chosen extractor returns
    # we run THREE orthogonal positive checks (NOT a fallback chain --
    # all must pass):
    #   1. _check_centerline_inside_solid: bbox-containment, fast O(N)
    #      sanity (catches gross wrong-location spines).
    #   2. _check_centerline_near_solid_surface (v4.51.0): STRONG
    #      per-point distance-to-solid check, catches spines that
    #      exit the wire tube envelope (Predicate 4 surface-rim,
    #      Predicate 5 wrong-radius); subsamples for performance.
    #   3. (bbox-cover check runs later in filaments_from_step on the
    #      filament paths, since filaments may extend slightly beyond
    #      the bare centerline)
    def _dispatch_and_verify(extractor, predicate_name, *args, **kwargs):
        result = extractor(*args, **kwargs)
        path_m, widths_m, heights_m = result
        _check_centerline_inside_solid(
            solid, path_m,
            f"extract_centerline_from_step({predicate_name})",
            cad_units_per_meter=cad_units_per_meter)
        # Strong distance check (v4.51.0): wire radius from mean
        # cross-section area (equivalent-circle assumption).
        mean_area_m2 = float(np.mean(widths_m * heights_m))
        wire_r_m = 0.0
        if mean_area_m2 > 0:
            wire_r_m = float(np.sqrt(mean_area_m2 / np.pi))
            _check_centerline_near_solid_surface(
                solid, path_m, wire_r_m,
                f"extract_centerline_from_step({predicate_name})",
                cad_units_per_meter=cad_units_per_meter)
        # Coverage check (v4.95.x): the STEP holds ONE solid = the
        # coil, so the centerline must SPAN the solid bbox to within a
        # wire-radius-scale slack on every axis.  An under-covering
        # path means the extractor traced only part of the conductor
        # -- e.g. a closed walk that short-circuited a lead-pair
        # junction into a sub-loop -- and would silently hand PEEC a
        # partial coil (with mass/length width pinning, an INFLATED
        # cross-section as well).
        bb = solid.bounding_box()
        s_min = np.array([float(bb.min.X), float(bb.min.Y),
                          float(bb.min.Z)]) / cad_units_per_meter
        s_max = np.array([float(bb.max.X), float(bb.max.Y),
                          float(bb.max.Z)]) / cad_units_per_meter
        slack = max(2.5 * wire_r_m, 5e-4)
        msgs = []
        for ax, name in enumerate(("x", "y", "z")):
            if (s_max[ax] - s_min[ax]) < 1e-9:
                continue
            gap_lo = float(path_m[:, ax].min() - s_min[ax])
            gap_hi = float(s_max[ax] - path_m[:, ax].max())
            for gap, side in ((gap_lo, "min"), (gap_hi, "max")):
                if gap > slack:
                    msgs.append(
                        f"{name}-{side}: centerline stops "
                        f"{gap * 1e3:.1f} mm short of the solid "
                        f"(slack {slack * 1e3:.1f} mm)")
        if msgs:
            raise ValueError(
                f"extract_centerline_from_step({predicate_name}): the "
                "extracted centerline does not span the conductor "
                "solid -- only part of the coil was traced.\n  "
                + "\n  ".join(msgs)
                + "\nCommon cause: a closed walk short-circuited a "
                "lead-pair junction into a sub-loop.  Regenerate the "
                "STEP (separate the leads by more than a wire "
                "diameter), or switch to --coil-solver bem-a "
                "--coil-vol <pre-meshed.vol> which bypasses spine "
                "extraction entirely.")
        return result

    # Predicate 1: multi-station loft of profiles (NON-united).  When
    # the solid has many consistent-area cross-section PLANE faces,
    # chain their centroids.  Robust for non-united multi-loft coils.
    if _collect_loft_cross_sections(solid):
        return _dispatch_and_verify(
            _centerline_from_cross_sections,
            "loft_cross_sections",
            solid, cad_units_per_meter)

    # Predicate 2: united multi-turn pancake.  Boolean unite consumes
    # the planar end-caps but cross-section CIRCLE edges remain (split
    # into 2 semicircles per cross-section by Cubit's unite).  Group
    # circles by arc_center, dedupe the semicircle pair into one per
    # cross-section, chain via NN + tangent continuity.  Handles
    # Kubota's 3turncoil united.stp class.
    circle_centers_radius = _collect_circle_edge_centers(solid)
    if circle_centers_radius is not None:
        centers_cad, median_r_cad = circle_centers_radius
        return _dispatch_and_verify(
            _centerline_from_circle_edge_centers,
            "circle_edge_centers",
            centers_cad, median_r_cad, cad_units_per_meter)

    # Predicate 3: single-loop revolution sweep.  Solid has at least
    # one TORUS / CYLINDER / CONE / REVOLUTION lateral surface AND at
    # least one PLANE end-cap face -- the gapped-torus / rect-torus
    # case.  ``_centerline_from_revolution_sweep`` extracts axis +
    # major-R + sweep angle analytically.  Note: a CLOSED full torus
    # has revolution surfaces but NO planar end-caps and therefore
    # falls through to predicate 5 below.
    has_revolution_face = any(
        f.geom_type in (GeomType.TORUS, GeomType.CYLINDER, GeomType.CONE,
                         GeomType.REVOLUTION)
        for f in solid.faces())
    has_plane_face = any(f.geom_type == GeomType.PLANE for f in solid.faces())
    if has_revolution_face and has_plane_face:
        return _dispatch_and_verify(
            _centerline_from_revolution_sweep,
            "revolution_sweep",
            solid, n_segments, cad_units_per_meter)

    # Predicates 4 / 5 require coil-topology classification (OPEN vs
    # CLOSED based on cap-face detection).
    from radia.coil_topology import extract_coil_topology as _extract_topo
    topo = _extract_topo(solid)

    # Predicate 4: OPEN coil (2 cap faces detected) without simple
    # revolution surfaces -- BSPLINE-lofted "arc + leads" geometries
    # such as keiko's 1turn_coil_loft_outsideline.step.
    # ``_centerline_from_open_spine`` samples the longest open lateral
    # rim edge as the spine, which correctly traces lead extensions.
    if topo.is_open:
        return _dispatch_and_verify(
            _centerline_from_open_spine,
            "open_spine",
            solid, n_segments, cad_units_per_meter)

    # Predicate 5: CLOSED loop (no caps, no revolution surface with
    # planar end-caps -- e.g. CLOSED bspline-lofted torus, racetrack
    # sweep).  ``_centerline_from_topology_spine`` walks the actual
    # solid (slab marching from the axis-agnostic seed), so any
    # orientation and any closed spine shape is traced; the median
    # station area recovers the cross-section.  This path is
    # CLOSED-only; OPEN coils route to predicate 4 above.
    return _dispatch_and_verify(
        _centerline_from_topology_spine,
        "topology_spine",
        solid, n_segments, cad_units_per_meter)


def _check_centerline_inside_solid(solid, path_m, source_tag,
                                     cad_units_per_meter=1.0,
                                     slack_factor=0.05):
    """Universal positive proof that the extracted centerline lies
    within the conductor's bounding box, plus a small slack
    (v4.50.0, CLAUDE.md "No Fallbacks - Fail Fast, Fail Loud").

    Background: the existing ``_check_filaments_cover_solid_bbox`` is
    an EXCLUSION proof on the FILAMENT extents -- it catches "spine
    doesn't extend to bbox extents" (under-coverage / fallback-radius
    overshoot).  It does NOT catch the orthogonal failure where the
    spine sits at the wrong location entirely (e.g. Predicate 5
    mapping a non-axisymmetric CLOSED racetrack loop to a planar
    circle of radius 0.85 * R_outer that misses the rectangular
    corners; the circle spine and the conductor bbox have similar
    extents, but the circle points outside the conductor's corner
    regions lie outside the conductor bbox).

    This check verifies each centerline point falls within the
    solid's bounding box + ``slack_factor`` of the bbox diagonal.
    It is a WEAK positive proof -- the bbox is convex and any
    non-convex conductor (e.g. a coil with a gap) has bbox points
    NOT inside the conductor.  A stronger check would use
    ``BRepClass3d_SolidClassifier`` to verify INSIDE-ness, but that
    classifier is unreliable on BSpline solids (tested 2026-05-16
    on a smooth sweep coil: 78%% of true-interior centerline points
    classified as OUT, including the wire axis along the lead).

    Catches:
    - Predicate 5 racetrack-as-circle (circle corners far outside
      bbox)
    - Predicate 4 picking an obviously-wrong edge whose extent
      exceeds the solid bbox
    - Predicate 1 UV sampling a wrong face (e.g. a stray helper
      surface) that bounces outside the bbox

    Does NOT catch (acceptable false-negatives on the weak bbox check;
    the bbox-cover check + corner detect cover most cases):
    - Predicate 4 picking a surface-rim edge that LIES inside the
      solid bbox but on the solid's surface (would need the
      stronger SolidClassifier, which is unreliable on BSpline)
    - Spines with wrong centerline within bbox (need per-point
      distance-to-surface check; deferred to v4.51.0+ pending a
      reliable OCC inside-test API)

    Args:
        solid: build123d Solid in CAD units.
        path_m: (N, 3) centerline polyline in METERS.
        source_tag: free-form string for the diagnostic (extractor name).
        cad_units_per_meter: scale to convert path_m back to CAD units.
        slack_factor: bbox padding as fraction of bbox diagonal.
            Default 0.05 (5%%) accommodates wire-radius extent at
            cap faces; spines that exceed 5%% are clearly wrong.
    """
    path = np.asarray(path_m, dtype=float)
    n_pts = len(path)
    if n_pts == 0:
        return

    pts_cad = path * cad_units_per_meter
    bb = solid.bounding_box()
    bb_min = np.array([float(bb.min.X), float(bb.min.Y), float(bb.min.Z)])
    bb_max = np.array([float(bb.max.X), float(bb.max.Y), float(bb.max.Z)])
    diag = float(np.linalg.norm(bb_max - bb_min))
    slack = slack_factor * diag

    over_min = np.maximum(bb_min - slack - pts_cad, 0.0)
    over_max = np.maximum(pts_cad - (bb_max + slack), 0.0)
    excursion = np.maximum(over_min, over_max).max(axis=1)
    bad = excursion > 0.0
    n_bad = int(bad.sum())
    if n_bad == 0:
        return  # PASS

    bad_idx = np.where(bad)[0]
    samples = []
    for idx in bad_idx[:5]:
        p_m = path[idx]
        e_mm = float(excursion[idx]) * 1e3 / cad_units_per_meter
        samples.append(
            f"  point {idx}/{n_pts - 1}: xyz=({p_m[0] * 1e3:+.2f}, "
            f"{p_m[1] * 1e3:+.2f}, {p_m[2] * 1e3:+.2f}) mm, "
            f"excursion beyond bbox+slack = {e_mm:.2f} mm")
    raise ValueError(
        f"{source_tag}: centerline extends beyond solid bbox + "
        f"{slack_factor:.0%} slack ({n_bad}/{n_pts} points outside).  "
        f"Solid bbox: x=[{bb_min[0] * 1e3 / cad_units_per_meter:+.1f},"
        f"{bb_max[0] * 1e3 / cad_units_per_meter:+.1f}] "
        f"y=[{bb_min[1] * 1e3 / cad_units_per_meter:+.1f},"
        f"{bb_max[1] * 1e3 / cad_units_per_meter:+.1f}] "
        f"z=[{bb_min[2] * 1e3 / cad_units_per_meter:+.1f},"
        f"{bb_max[2] * 1e3 / cad_units_per_meter:+.1f}] mm; "
        f"slack=({slack * 1e3 / cad_units_per_meter:.1f} mm = "
        f"{slack_factor:.0%} of diag).\n"
        f"Sample outside points:\n" + "\n".join(samples) +
        f"\nThis means the extracted spine has the wrong shape for "
        f"the conductor's topology.  Common cause: Predicate 5 maps "
        f"a non-axisymmetric CLOSED loop (e.g. racetrack) to a "
        f"planar circle whose corners lie outside the racetrack bbox.  "
        f"FIX: regenerate the STEP with a clean single-piece BSPLINE "
        f"lateral so Predicate 1 (UV-map sampling) handles it "
        f"directly, OR ensure the CAD topology matches one of the "
        f"supported predicate classes (gapped torus / loft-of-circles "
        f"/ united multi-turn pancake / sweep)."
    )


def _check_centerline_near_solid_surface(solid, path_m, wire_radius_m,
                                            source_tag,
                                            cad_units_per_meter=1.0,
                                            distance_tolerance_factor=1.10,
                                            max_violation_fraction=0.05,
                                            sample_count=20):
    """STRONG positive proof that the extracted centerline lies within
    the wire tube envelope (v4.51.0, CLAUDE.md "No Fallbacks - Fail
    Fast, Fail Loud").

    Background: ``_check_centerline_inside_solid`` (v4.50.0) only
    checks bbox-containment, which is convex and admits any spine
    that touches all axis extents.  This stronger check uses
    ``BRepExtrema_DistShapeShape`` to compute the actual distance
    from each centerline point to the solid boundary -- points
    INSIDE the solid return 0, points OUTSIDE return the distance to
    the nearest surface.  A correctly-placed centerline should be
    either INSIDE the wire tube (d=0) or within
    ``distance_tolerance_factor * wire_radius_m`` of the boundary
    (accounting for parallel-transport displacement on smooth
    spline sweeps -- verified empirically: 100% of points on a
    smooth build123d sweep coil's per-station-mean centerline fall
    within wire_radius of the lateral surface).

    Performance: BRepExtrema_DistShapeShape is O(face_count) per
    point.  For a typical 100-point centerline on a 700-face STEP
    this is ~10-100 ms/point.  We sub-sample ``sample_count`` points
    evenly along the centerline so the check stays bounded
    (default 20 points -> ~1-2 s on a 700-face STEP).  Sub-sampling
    is sufficient because failure modes (wrong-radius spine, surface-
    rim spine) affect contiguous regions, not isolated points.

    Catches that ``_check_centerline_inside_solid`` MISSES:
    - Predicate 4 picking a surface-rim edge whose centerline lies
      ON the solid surface but INSIDE the bbox
    - Predicate 5 racetrack-as-circle where corners are inside bbox
      but outside the actual conductor cross-section
    - Wrong-radius spine that stays within bbox but exits the wire

    Args:
        solid: build123d Solid in CAD units.
        path_m: (N, 3) centerline polyline in METERS.
        wire_radius_m: nominal wire radius in METERS (from
            ``cross_section_radius_m`` in the topo dict, or
            ``sqrt(mean(w*h)/pi)`` from Path 3 widths/heights).
        source_tag: free-form string for the diagnostic.
        cad_units_per_meter: scale to convert path_m back to CAD units.
        distance_tolerance_factor: max allowed distance as a multiple
            of wire_radius_m.  Default 1.10 (10% slack for numerical
            noise on the swept tube boundary).
        max_violation_fraction: max fraction of sub-sampled points
            allowed to exceed the tolerance.  Default 0.05 (5%).
        sample_count: number of points to sub-sample for the OCC
            distance query.  Default 20.
    """
    from OCP.BRepExtrema import BRepExtrema_DistShapeShape
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
    from OCP.gp import gp_Pnt

    path = np.asarray(path_m, dtype=float)
    n_pts = len(path)
    if n_pts == 0 or wire_radius_m <= 0:
        return

    # Sub-sample evenly.
    if n_pts <= sample_count:
        sample_idx = np.arange(n_pts)
    else:
        sample_idx = np.linspace(0, n_pts - 1, sample_count, dtype=int)
    pts_cad = path[sample_idx] * cad_units_per_meter

    tolerance_cad = (distance_tolerance_factor * wire_radius_m
                      * cad_units_per_meter)

    violations = []
    for k, i in enumerate(sample_idx):
        vtx = BRepBuilderAPI_MakeVertex(
            gp_Pnt(float(pts_cad[k][0]),
                    float(pts_cad[k][1]),
                    float(pts_cad[k][2]))).Vertex()
        ext = BRepExtrema_DistShapeShape(vtx, solid.wrapped)
        ext.Perform()
        d_cad = ext.Value()
        if d_cad > tolerance_cad:
            violations.append((int(i), d_cad / cad_units_per_meter))

    n_sampled = len(sample_idx)
    n_violated = len(violations)
    violation_fraction = n_violated / n_sampled
    if violation_fraction <= max_violation_fraction:
        return  # PASS

    samples = []
    for i, d_m in violations[:5]:
        p_m = path[i]
        samples.append(
            f"  point {i}/{n_pts - 1}: xyz=({p_m[0] * 1e3:+.2f}, "
            f"{p_m[1] * 1e3:+.2f}, {p_m[2] * 1e3:+.2f}) mm, "
            f"distance to solid = {d_m * 1e3:.2f} mm "
            f"(> {distance_tolerance_factor:.1f}x wire_r = "
            f"{distance_tolerance_factor * wire_radius_m * 1e3:.2f} mm)")
    raise ValueError(
        f"{source_tag}: centerline exits the wire tube envelope "
        f"({n_violated}/{n_sampled} sampled points > "
        f"{distance_tolerance_factor:.0%} of wire radius = "
        f"{distance_tolerance_factor * wire_radius_m * 1e3:.2f} mm "
        f"from the solid boundary; threshold "
        f"{max_violation_fraction:.0%}).\n"
        f"Wire radius: {wire_radius_m * 1e3:.2f} mm.\n"
        f"Sample violating points:\n" + "\n".join(samples) +
        f"\nThis is the STRONG positive proof that the spine actually "
        f"traces the conductor cross-section centroid path -- a spine "
        f"that exits the wire tube means the extractor mis-classified "
        f"the geometry.  Common causes: (a) Predicate 4 picked an edge "
        f"that lies on the solid's lateral surface rather than its "
        f"centroid; (b) Predicate 5 (topology_spine) used the bbox-"
        f"derived radius (0.85 * R_outer) which differs from the "
        f"conductor's actual centroid radius; (c) CAD has self-"
        f"intersecting or non-manifold topology.  "
        f"FIX: regenerate the STEP with a clean single-piece BSPLINE "
        f"lateral so Predicate 1 (UV-map sampling) handles it directly."
    )


def _centerline_from_filament_paths(filament_paths):
    """Derive a centerline from a filament_paths list by averaging the
    n_peri filaments at each station.

    For Paths 1/2/2b/2c the filaments are placed around the cross-
    section perimeter; their mean per station IS the centerline.
    Used to feed `_check_centerline_inside_solid` from extractors
    that build filaments directly without an explicit centerline.
    """
    import numpy as np
    paths_arr = np.asarray(filament_paths)
    # Shape: (n_fil, n_seg, 2, 3) [(start, end) per segment per fil]
    if paths_arr.ndim != 4 or paths_arr.shape[2] != 2 or paths_arr.shape[3] != 3:
        # Fall back to flat sample
        return paths_arr.reshape(-1, 3)
    n_fil, n_seg, _, _ = paths_arr.shape
    # Centerline = per-station mean of filament endpoints
    # Station 0 = mean of seg-0 starts; station 1..n_seg = mean of seg-i ends
    stations = np.empty((n_seg + 1, 3), dtype=np.float64)
    stations[0] = paths_arr[:, 0, 0, :].mean(axis=0)
    for i in range(n_seg):
        stations[i + 1] = paths_arr[:, i, 1, :].mean(axis=0)
    return stations


def _check_filaments_cover_solid_bbox(topo, solid_bbox_min, solid_bbox_max,
                                       tier: str, slack_factor: float = 1.5):
    """Sanity-check that the extracted filaments span the conductor solid.

    Two failure modes both detected:

    1. **Coverage gap** (lead skipped): filament max in some axis falls
       short of solid bbox max by more than ``slack`` -- typical when
       a planar bbox-radius spine arc bypasses a straight lead bar.
    2. **Overshoot** (wrong-radius fallback): filament min/max
       extends BEYOND the solid bbox by more than ``slack`` -- typical
       when the bbox-radius fallback produces a spine OUTSIDE the
       actual conductor centerline.

    Both cases mean PEEC silently produces a topologically-wrong path.
    Raise with a hint.

    ``slack_factor`` is in units of cross_section_radius.  Default 1.5
    accommodates loft chamfers / round-downs without false-positive on
    clean coils, while still catching keiko's lead bar (6 mm gap,
    wire_radius=3 mm, ratio 2.0).
    """
    import numpy as np
    paths = topo.get("filament_paths")
    if paths is None:
        return
    paths_arr = np.asarray(paths)
    if paths_arr.size == 0:
        return
    # Shape variants: (n_fil, n_seg, 2, 3) or (n_fil, n_pts, 3).
    if paths_arr.ndim >= 3:
        pts = paths_arr.reshape(-1, 3)
    else:
        raise ValueError(
            "_check_filaments_cover_solid_bbox: filament_paths has "
            f"unexpected ndim={paths_arr.ndim} shape={paths_arr.shape}; "
            "expected (n_fil, n_seg, 2, 3) or (n_fil, n_pts, 3)."
        )
    fil_min = pts.min(axis=0)
    fil_max = pts.max(axis=0)
    wire_r = float(topo.get("cross_section_radius_m") or 0.0)
    if wire_r > 0:
        slack = slack_factor * wire_r
    else:
        # Tiers that don't report cross_section_radius (e.g. Tier 2b
        # circle-uv).  Estimate from the smallest non-zero solid bbox
        # extent (typically the wire thickness for a planar coil).
        extents = np.asarray(solid_bbox_max) - np.asarray(solid_bbox_min)
        nonzero = extents[extents > 1e-6]
        if nonzero.size:
            slack = slack_factor * float(nonzero.min()) / 2.0
        else:
            slack = 1e-3
    slack = max(slack, 5e-4)  # absolute floor 0.5 mm

    msgs = []
    for axis, name in enumerate(("x", "y", "z")):
        s_min = float(solid_bbox_min[axis])
        s_max = float(solid_bbox_max[axis])
        s_ext = s_max - s_min
        if s_ext < 1e-9:
            continue
        f_min, f_max = float(fil_min[axis]), float(fil_max[axis])
        # Coverage gaps.
        gap_lo = max(f_min - s_min, 0.0)
        gap_hi = max(s_max - f_max, 0.0)
        # Overshoot beyond solid.
        over_lo = max(s_min - f_min, 0.0)
        over_hi = max(f_max - s_max, 0.0)
        if gap_lo > slack:
            msgs.append(
                f"{name}: filament min {f_min:+.4f} but solid min "
                f"{s_min:+.4f} (gap {gap_lo*1e3:.1f} mm > "
                f"slack {slack*1e3:.1f} mm) -- lead/extension not traced")
        if gap_hi > slack:
            msgs.append(
                f"{name}: filament max {f_max:+.4f} but solid max "
                f"{s_max:+.4f} (gap {gap_hi*1e3:.1f} mm > "
                f"slack {slack*1e3:.1f} mm) -- lead/extension not traced")
        if over_lo > slack:
            msgs.append(
                f"{name}: filament min {f_min:+.4f} BELOW solid min "
                f"{s_min:+.4f} (overshoot {over_lo*1e3:.1f} mm > "
                f"slack {slack*1e3:.1f} mm) -- bbox-radius fallback spine")
        if over_hi > slack:
            msgs.append(
                f"{name}: filament max {f_max:+.4f} ABOVE solid max "
                f"{s_max:+.4f} (overshoot {over_hi*1e3:.1f} mm > "
                f"slack {slack*1e3:.1f} mm) -- bbox-radius fallback spine")
    if msgs:
        raise ValueError(
            f"Filament path does not match the conductor solid bbox "
            f"(tier={tier!r}, wire_radius={wire_r*1e3:.1f} mm, "
            f"slack={slack*1e3:.1f} mm). "
            f"Spine extraction likely bypassed a lead or used the "
            f"bbox-radius fallback at a wrong radius -- PEEC would "
            f"silently produce a number that does not reflect the "
            f"intended coil.\nDiagnostics:\n  " + "\n  ".join(msgs)
            + "\nHINT: regenerate the STEP with a clean centerline "
              "(e.g. consistent loft vertex alignment so the lateral "
              "surface is a single dominant BSPLINE / TORUS instead of "
              "two equal half-faces split at the equator), OR switch "
              "to --coil-solver bem-a --coil-vol <pre-meshed.vol> "
              "which bypasses spine extraction entirely.")


_PEEC_CACHE_FORMAT_VERSION = 2


def _peec_cache_path(step_path: str, n_peri, sigma, nwinc, nhinc,
                     n_slices, cad_units_per_meter, use_coil_builder):
    """Disk-cache filename for filaments_from_step output.

    One cache file per (step_path, params) combination so different
    parameter sets (different n_peri / sigma / units) don't trample
    each other.
    """
    import os
    base, _ext = os.path.splitext(step_path)
    if use_coil_builder and n_peri is not None:
        flavour = f"peri{int(n_peri)}"
    elif use_coil_builder:
        flavour = f"grid{int(nwinc)}x{int(nhinc)}"
    else:
        flavour = f"legacy{int(nwinc)}x{int(nhinc)}"
    return (f"{base}.peec_v{_PEEC_CACHE_FORMAT_VERSION}_{flavour}"
            f"_n{int(n_slices)}_sigma{sigma:.2e}"
            f"_u{cad_units_per_meter:g}.json")


def _peec_cache_step_sha256(step_path: str) -> str:
    """SHA256 of the STEP file's bytes."""
    import hashlib
    h = hashlib.sha256()
    with open(step_path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _peec_cache_load(cache_path: str, expected_sha: str, sigma: float):
    """Return rebuilt topo dict from cache, or None on miss / mismatch.

    The cache JSON stores the lightweight filament topology data
    (paths + cell_wh + cross-section radii).  The PEECCircuitSolver
    is rebuilt via build_bundle_solver -- a fast (~0.5 s) call that
    re-computes mutual L from the cached filament paths.
    """
    import json, os
    if not os.path.isfile(cache_path):
        return None
    try:
        with open(cache_path, encoding="utf-8") as f:
            data = json.load(f)
    except Exception:
        return None
    # Format 1 predates the default-walker CAD-unit boundary and may
    # contain millimetre coordinates cached as metres.  Never rebuild a
    # solver from those silently-wrong paths; force one checked extraction
    # and rewrite the cache in the current format.
    if data.get("_format_version") != _PEEC_CACHE_FORMAT_VERSION:
        return None
    if data.get("step_sha256") != expected_sha:
        return None
    paths_serial = data.get("filament_paths")
    cell_wh_serial = data.get("cell_wh")
    if not paths_serial or not cell_wh_serial:
        return None
    # Deserialise filament_paths: list of polyline pieces, each piece is
    # ((p1, p2), (p2, p3), ...).  JSON lists -> tuples expected by
    # PEECBuilder.add_connected_segment.
    filament_paths = [
        [tuple(tuple(pt) for pt in pair) for pair in fil]
        for fil in paths_serial
    ]
    cell_wh = [
        [tuple(piece) for piece in fil]
        for fil in cell_wh_serial
    ]
    # Rebuild the solver from cached paths.
    from radia.peec_bundle import build_bundle_solver
    solver, seg_of_fil, port_p, port_m = build_bundle_solver(
        filament_paths, dw=None, dh=None, sigma=sigma, cell_wh=cell_wh)
    topo = dict(data.get("aux", {}))
    topo.update({
        "solver": solver,
        "seg_of_filament": seg_of_fil,
        "filament_paths": filament_paths,
        "cell_wh": cell_wh,
        "port_plus": port_p,
        "port_minus": port_m,
        "n_loop": len(filament_paths),
        "_peec_cache_hit": True,
    })
    return topo


def _peec_cache_save(cache_path: str, topo: dict, step_sha: str,
                       params: dict):
    """Persist topo to disk.  Solver and any non-serialisable fields
    are dropped; they are rebuilt on load via build_bundle_solver."""
    import json
    # JSON-safe copy of filament_paths and cell_wh.
    fil = []
    for filk in topo.get("filament_paths") or []:
        fil.append([[list(pt) for pt in pair] for pair in filk])
    cwh = []
    for filk in topo.get("cell_wh") or []:
        cwh.append([list(piece) for piece in filk])
    # Keep the small scalar/string aux fields that downstream code
    # reads (cross_section_radius_m_mean, source, n_path_pts, ...).
    aux = {}
    for k, v in topo.items():
        if k in ("solver", "filament_paths", "cell_wh", "seg_of_filament",
                 "port_plus", "port_minus", "n_loop", "_peec_cache_hit"):
            continue
        try:
            json.dumps(v)
        except TypeError:
            continue  # skip non-serialisable extras
        aux[k] = v
    data = {
        "_format_version": _PEEC_CACHE_FORMAT_VERSION,
        "step_sha256": step_sha,
        "params": params,
        "filament_paths": fil,
        "cell_wh": cwh,
        "aux": aux,
    }
    try:
        with open(cache_path, "w", encoding="utf-8") as f:
            json.dump(data, f)
    except Exception:
        # Cache-write failure is non-fatal -- the next run just
        # re-computes.  Do not raise here.
        pass


def filaments_from_step(step_path: str,
                        sigma: float = 5.8e7,
                        nwinc: int = 1,
                        nhinc: int = 1,
                        n_peri: Optional[int] = None,
                        cad_units_per_meter: float = 1.0,
                        n_slices: int = 200,
                        use_coil_builder: bool = True):
    """End-to-end STEP -> PEEC topology, with on-disk cache.

    Set RADIA_PEEC_CACHE_DISABLE=1 to bypass the cache entirely.
    The cache file is
    ``<step>.peec_v<format>_<flavour>_n<N>_sigma<S>_u<U>.json``
    next to the .step.  Keyed by the SHA256 of the .step file's
    bytes plus the parameter set, so any edit to the CAD invalidates
    the cache automatically.  Cache hit rebuilds the PEEC solver
    (~0.5 s) instead of re-running OCC topology analysis (~25 s).

    See ``_filaments_from_step_compute`` for the actual extraction
    logic and arguments documentation.
    """
    import os
    # Normalize before formatting the cache key.  Geometry plausibility is
    # checked on cache creation in _filaments_from_step_compute; a current-
    # format hit is tied to the exact STEP SHA and normalized unit scale.
    cad_units_per_meter = _normalize_cad_units_per_meter(
        cad_units_per_meter, "filaments_from_step")
    cache_disabled = os.environ.get("RADIA_PEEC_CACHE_DISABLE", "0") != "0"
    cache_path = None
    step_sha = None
    if not cache_disabled and os.path.isfile(step_path):
        cache_path = _peec_cache_path(
            step_path, n_peri, sigma, nwinc, nhinc,
            n_slices, cad_units_per_meter, use_coil_builder)
        step_sha = _peec_cache_step_sha256(step_path)
        cached = _peec_cache_load(cache_path, step_sha, sigma)
        if cached is not None:
            return cached

    topo = _filaments_from_step_compute(
        step_path,
        sigma=sigma,
        nwinc=nwinc,
        nhinc=nhinc,
        n_peri=n_peri,
        cad_units_per_meter=cad_units_per_meter,
        n_slices=n_slices,
        use_coil_builder=use_coil_builder,
    )

    if cache_path is not None and step_sha is not None:
        _peec_cache_save(cache_path, topo, step_sha, {
            "n_peri": n_peri,
            "sigma": sigma,
            "nwinc": nwinc,
            "nhinc": nhinc,
            "n_slices": n_slices,
            "cad_units_per_meter": cad_units_per_meter,
            "use_coil_builder": use_coil_builder,
        })
    return topo


def _filaments_from_step_compute(step_path: str,
                                  sigma: float = 5.8e7,
                                  nwinc: int = 1,
                                  nhinc: int = 1,
                                  n_peri: Optional[int] = None,
                                  cad_units_per_meter: float = 1.0,
                                  n_slices: int = 200,
                                  use_coil_builder: bool = True):
    """End-to-end: STEP solid -> PEEC topology (no cache).

    The centerline is ALWAYS auto-detected from the STEP solid via
    ``extract_centerline_from_step`` (dispatches across multiple
    topology-aware paths: loft cross-sections, circle-edge stations,
    torus sweep, topology spine, open-spine longest-edge).  There is
    no caller-provided ``path_points_m`` override -- if the auto-
    detection produces a spine that does not cover the conductor's
    bounding box, ``_check_filaments_cover_solid_bbox`` raises
    fail-fast so the caller fixes the CAD rather than papering over
    the bad geometry with a hand-crafted centerline JSON.

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
        sigma: Conductivity [S/m].
        nwinc, nhinc: Sub-filament subdivision for the volume-grid
            placement.  Ignored when ``n_peri`` is set.
        n_peri: If given, place ``n_peri`` filaments on the cross-section
            PERIMETER only (thin-skin regime, d/delta >= 3).  Takes
            priority over nwinc/nhinc.  Requires use_coil_builder=True.
        cad_units_per_meter: Scale factor.  Default 1.0 = STEP coordinates
            are in metres (CLAUDE.md "Unit System Policy: Radia always uses
            meters").  Pass 1000.0 if the STEP is in millimetres.
        n_slices: Z-slice count for auto-extraction.
        use_coil_builder: If True (default), use CoilBuilder path for
            profile-aware filament placement (volume-grid via
            nwinc/nhinc OR n_peri perimeter-only).  If False, use the
            legacy C++ ExpandFilaments path which only supports a
            rectangular sub-filament grid.  These are SEPARATE entry
            paths chosen by the caller -- per CLAUDE.md "No Fallbacks",
            there is no automatic switchover when one fails; failures
            propagate as ValueError.

    Returns:
        topology_dict from PEECBuilder.build_topology().
    """
    # Units sanity: fail-loud on a mm-as-metres mismatch (silent
    # ~1000x-wrong L otherwise; No-Fallbacks / Fail-Loud).  Covers BOTH
    # filament paths (coil-builder walker and n_peri UV tiers).
    from radia._b3d_shim import import_step as _import_step_units_chk
    cad_units_per_meter = _check_solid_extent_plausible(
        _import_step_units_chk(step_path), cad_units_per_meter,
        "filaments_from_step")

    if use_coil_builder and n_peri is None:
        # Volume-grid placement (nwinc/nhinc) needs the profile-aware
        # walker path.  Only the walker knows per-segment (w, h).
        return _filaments_via_coil_builder(
            step_path, sigma=sigma, nwinc=nwinc, nhinc=nhinc,
            n_slices=n_slices,
            cad_units_per_meter=cad_units_per_meter,
            start_hint=None,  # auto-detect (labels -> axis-agnostic seed)
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
        from radia._b3d_shim import import_step
        solid = import_step(step_path)

        # Compute the solid bounding box ONCE; reused by the
        # `_verify_topo` helper below to run the bbox-cover check on
        # every Path 1/2/2b/2c filament output (in addition to the
        # centerline-inside-bbox check from v4.50.0).  Both checks
        # must pass (orthogonal failure modes, not a fallback chain):
        # bbox-cover catches under-coverage / wrong-radius overshoot;
        # inside-bbox catches wrong-location spine.  Single-extractor
        # geometries that fail either get a hard ValueError pointing
        # at CAD regeneration -- this prevents the keiko-class
        # ``_outsideline`` regression (silent NaN L from an
        # under-covered spine, 2026-05-15) and the Predicate 5
        # racetrack-as-circle class (2026-05-16).
        bb = solid.bounding_box()
        solid_bbox_min = np.array([float(bb.min.X), float(bb.min.Y),
                                    float(bb.min.Z)]) / cad_units_per_meter
        solid_bbox_max = np.array([float(bb.max.X), float(bb.max.Y),
                                    float(bb.max.Z)]) / cad_units_per_meter

        # Helper: derive an effective centerline from the per-station
        # filament mean, then run THREE orthogonal positive checks
        # (bbox-cover for under-coverage, inside-solid for
        # wrong-location, near-surface for wrong-radius/surface-rim).
        # Per CLAUDE.md "No Fallbacks", all three must pass or this
        # raises.  Not a fallback chain -- they catch disjoint failure
        # modes.
        def _verify_topo(topo_dict, tier_name):
            _check_filaments_cover_solid_bbox(
                topo_dict, solid_bbox_min, solid_bbox_max, tier=tier_name)
            path_eff = _centerline_from_filament_paths(
                topo_dict.get("filament_paths", []))
            _check_centerline_inside_solid(
                solid, path_eff,
                f"filaments_from_step({tier_name})",
                cad_units_per_meter=cad_units_per_meter)
            # Strong distance check (v4.51.0): wire radius from topo
            # dict (cross_section_radius_m) or derived from cell_wh.
            wire_r_m = float(topo_dict.get("cross_section_radius_m") or 0)
            if wire_r_m <= 0:
                cell_wh = topo_dict.get("cell_wh") or []
                if cell_wh:
                    areas = []
                    for fil in cell_wh:
                        for (w, h) in fil:
                            areas.append(float(w) * float(h))
                    if areas:
                        wire_r_m = float(np.sqrt(np.mean(areas) / np.pi))
            if wire_r_m > 0:
                _check_centerline_near_solid_surface(
                    solid, path_eff, wire_r_m,
                    f"filaments_from_step({tier_name})",
                    cad_units_per_meter=cad_units_per_meter)

        # Path 1: BSPLINE / TORUS / CYLINDER lateral surface UV grid.
        # ``_find_lateral_surface`` is a strict predicate: when it
        # returns a face, UV sampling MUST succeed (UV-closure check
        # already passed).  Per CLAUDE.md "No Fallbacks", any failure
        # here propagates -- we do NOT silently try Path 2.
        lateral = _find_lateral_surface(solid)
        if lateral is not None:
            topo_uv = _filaments_from_lateral_surface_uv(
                lateral, cad_units_per_meter=cad_units_per_meter,
                sigma=sigma,
                n_stations=max(20, n_slices // 2),
                n_peri=n_peri,
                source_tag="step_uv")
            _verify_topo(topo_uv, "step_uv")
            return topo_uv

        # Path 2: per-station planar end-cap faces (NON-united loft).
        loft = _try_extract_loft_with_profile(
            solid, cad_units_per_meter=cad_units_per_meter)
        if loft is not None:
            path_m, faces_ordered = loft
            topo_pst = _filaments_from_per_station_faces(
                path_m, faces_ordered,
                sigma=sigma, n_peri=n_peri,
                source_tag="step_per_station")
            _verify_topo(topo_pst, "step_per_station")
            return topo_pst

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
            _verify_topo(topo_circle, "step_circle_uv")
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
        from radia._b3d_shim import GeomType
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
                _verify_topo(topo_section, "step_section_planes")
                return topo_section

        # Path 3: constant equivalent-circle via existing dispatch.
        # extract_centerline_from_step ALREADY runs
        # _check_centerline_inside_solid on its dispatched extractor's
        # output, so the centerline is verified before we get here.
        # We still run the bbox-cover check on the filament_paths
        # (orthogonal: catches "filaments extend beyond solid", which
        # the centerline check misses since the centerline itself is
        # at the conductor centroid, not its surface).
        path_m, w_m, h_m = extract_centerline_from_step(
            step_path, n_segments=n_slices,
            cad_units_per_meter=cad_units_per_meter)
        mean_area = float(np.mean(w_m * h_m))
        r_m = float(np.sqrt(mean_area / np.pi))
        topo = filaments_from_polyline(
            path_m, r_m,
            sigma=sigma, n_peri=n_peri,
            source_tag="step_longest_edge")
        _check_filaments_cover_solid_bbox(
            topo, solid_bbox_min, solid_bbox_max,
            tier="step_longest_edge")
        return topo

    if n_peri is not None:
        raise ValueError(
            "n_peri requires use_coil_builder=True; the legacy "
            "C++ ExpandFilaments path is volume-grid-only.")

    # Legacy path: rectangular grid via C++ ExpandFilaments.
    # Centerline always comes from STEP auto-detection (no caller-
    # provided override -- see module docstring "Fail Fast Loud").
    path_m, widths_m, heights_m = extract_centerline_from_step(
        step_path, n_segments=n_slices,
        cad_units_per_meter=cad_units_per_meter,
    )

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
      - If omitted, the orientation-agnostic seed probes the wire
        tangent by minimum cross-section area (any coil orientation,
        open or closed loops).

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
    from radia.coil_from_step import _axis_agnostic_seed

    # Validate and normalize before BRep conversion or walking.  The STEP
    # entry point has the same guard; the in-memory path must not turn an
    # invalid zero scale into a late ZeroDivisionError or accept a unit
    # mismatch that silently rescales PEEC geometry.
    cad_units_per_meter = _check_solid_extent_plausible(
        shape, cad_units_per_meter, "filaments_from_shape")

    # --- 1. Resolve start_hint (CAD units) ---
    start_hint = None
    if port_face is not None:
        start_hint = _bd_face_to_start_hint(port_face, cad_units_per_meter)

    # --- 2. Convert to netgen.occ solid ---
    ng_solid = _bd_shape_to_netgen_solid(shape)

    # --- 3. orientation-agnostic seed if still no hint ---
    if start_hint is None:
        start_hint = _axis_agnostic_seed(ng_solid)

    # --- 4. Shared walking-plane core (scaling + verification inside).
    # Tight AABB from the build123d input shape (netgen bbox unreliable,
    # see _filaments_from_walked_centerline docstring).
    bb = shape.bounding_box()
    solid_bbox_cad = (
        np.array([bb.min.X, bb.min.Y, bb.min.Z], dtype=float),
        np.array([bb.max.X, bb.max.Y, bb.max.Z], dtype=float))
    return _filaments_from_walked_centerline(
        ng_solid, start_hint, sigma=sigma, nwinc=nwinc, nhinc=nhinc,
        n_peri=None, cad_units_per_meter=cad_units_per_meter,
        solid_bbox_cad=solid_bbox_cad)


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


def _filaments_from_walked_centerline(ng_solid, start_hint, *, sigma,
                                      nwinc, nhinc, n_peri,
                                      cad_units_per_meter, solid_bbox_cad):
    """Shared walking-plane core: netgen solid -> centerline -> metres ->
    CoilBuilder -> filaments -> bundle solver -> coverage verification.

    Used by ``_filaments_via_coil_builder`` (STEP-path entry) and
    ``filaments_from_shape`` (in-memory build123d entry) so the two
    wrappers cannot drift.

    Unit boundary: the walker operates in the solid's native CAD units
    (see ``_bd_face_to_start_hint``); this core scales the walked
    centerline by ``1 / cad_units_per_meter`` BEFORE CoilBuilder, so
    everything downstream (filaments, PEEC, L) is in metres.  Before
    v4.95.27 this scaling was missing: a mm-authored STEP flowed into
    PEEC as if metres and returned a silently ~1000x-too-large L.

    Verification (No-Fallbacks / Fail-Loud): the walker can halt
    mid-solid (e.g. at a sharp spine corner) and return a PARTIAL
    centerline; ``_check_filaments_cover_solid_bbox`` turns that into a
    hard ValueError instead of a silently-wrong partial coil (observed
    on B_rect_sweep: a single ~97 mm leg of a ~600 mm racetrack).

    ``solid_bbox_cad`` is the TIGHT conductor AABB in CAD units,
    ``(min_xyz, max_xyz)``, computed by the caller from build123d.
    Do NOT derive it from ``ng_solid.bounding_box``: the netgen.occ
    bbox of a loaded shape is unreliable for verification -- measured
    2026-07-29 on a STEP-roundtripped torus, it inflates curved x/y
    extents by ~8% (bspline control-polygon hull) and collapses z to
    +-1e-7 (falsely degenerate) while ``.mass`` stays correct.
    """
    from radia.coil_from_step import (extract_centerline,
                                      scale_centerline_result,
                                      to_coil_builder)
    from radia.peec_bundle import build_bundle_solver
    import numpy as np

    res = extract_centerline(ng_solid, start_hint=start_hint, verbose=False)
    res = scale_centerline_result(res, 1.0 / cad_units_per_meter)

    # Reconstruct CoilBuilder from centerline + profiles (metres).
    coil, _segs = to_coil_builder(res, current=1.0)

    # Generate profile-aware filaments.
    if n_peri is not None:
        paths, _currents = coil.to_filaments_peri(
            n_peri=n_peri, frequency=0.0, sigma=sigma)
    else:
        paths, _currents = coil.to_filaments(
            nw=nwinc, nh=nhinc, frequency=0.0, sigma=sigma)
    cell_wh = coil._last_filament_info.get('cell_wh')

    # Build PEEC topology via bundle solver.  build_bundle_solver
    # internally calls PEECBuilder.build_topology() and remaps nodes for
    # the parallel-bundle topology.  It returns a PEECCircuitSolver (not
    # a raw topo dict).  We return the solver plus metadata.
    solver, seg_of_fil, port_p, port_m = build_bundle_solver(
        paths, dw=None, dh=None, sigma=sigma, cell_wh=cell_wh)

    # Equivalent wire radius for the coverage-check slack: half the
    # largest cross-section bounding dimension over all stations (covers
    # flat-ribbon profiles where the equivalent-circle radius would
    # under-estimate the legitimate centerline-to-surface offset).
    r_eq = 0.0
    for prof in res.profiles:
        w, h = prof.bounding_wh()
        r_eq = max(r_eq, 0.5 * float(max(w, h)))

    result = {
        "solver": solver,
        "seg_of_filament": seg_of_fil,
        "filament_paths": paths,
        "cell_wh": cell_wh,
        "coil_builder": coil,
        "n_loop": len(paths),
        "port_plus": port_p,
        "port_minus": port_m,
        "cross_section_radius_m": r_eq,
    }

    bbox_min = np.asarray(solid_bbox_cad[0], dtype=float) / cad_units_per_meter
    bbox_max = np.asarray(solid_bbox_cad[1], dtype=float) / cad_units_per_meter
    # slack_factor 2.5 (vs the UV tiers' 1.5): a closed-loop CENTERLINE
    # legitimately stops r_wire short of the solid extent on each side,
    # and the current 1x1 to_filaments sample on a CircleProfile sits at
    # r/sqrt(2) off the centerline (area-preserving (0.5, 0.5) is not
    # the geometric center), so the legitimate shortfall reaches
    # r * (1 + 1/sqrt(2)) ~= 1.71 * r_eq.  The tier's actual failure
    # mode -- the walker halting mid-solid -- leaves O(10 * r_eq) gaps,
    # so 2.5 keeps full detection margin.
    _check_filaments_cover_solid_bbox(
        result, bbox_min, bbox_max, tier="coil_builder", slack_factor=2.5)
    return result


def _filaments_via_coil_builder(step_path, sigma, nwinc, nhinc, n_slices,
                                cad_units_per_meter,
                                start_hint=None, n_peri=None):
    """CoilBuilder path: STEP -> centerline -> CoilBuilder -> to_filaments.

    Uses Profile.sample_at() for cross-section-aware filament placement.
    Supports circular, rectangular, and lofted cross-sections.

    If ``n_peri`` is given, uses perimeter-only placement via
    ``to_filaments_peri`` (thin-skin limit) and ignores ``nwinc``/``nhinc``.
    """
    from radia.coil_from_step import load_step_solid, _axis_agnostic_seed
    from radia._b3d_shim import import_step as _import_step_bbox

    # Load the conductor solid ONCE (seed + walk share it).
    ng_solid = load_step_solid(step_path)

    # Tight conductor AABB for coverage verification (build123d; the
    # netgen.occ bbox of a loaded shape is unreliable -- see
    # _filaments_from_walked_centerline docstring).
    bb = _import_step_bbox(step_path).bounding_box()
    solid_bbox_cad = (
        np.array([bb.min.X, bb.min.Y, bb.min.Z], dtype=float),
        np.array([bb.max.X, bb.max.Y, bb.max.Z], dtype=float))

    # Start-hint resolution order (first hit wins):
    #   1. Caller-supplied start_hint (CAD units).
    #   2. XCAF label on a "peec_port*" child (FreeCAD Import.export, or
    #      build123d in-memory via filaments_from_shape).  Requires the
    #      optional build123d dependency; silently skipped otherwise.
    #   3. Orientation-agnostic seed (v4.95.x): probes the wire tangent
    #      by minimum cross-section area, so any coil orientation (and
    #      open OR closed loops) seeds correctly.  Replaced the legacy
    #      z-axis-torus bbox heuristic that broke on x/y-axis coils.
    if start_hint is None:
        start_hint = _start_hint_from_step_labels(step_path)
    if start_hint is None:
        start_hint = _axis_agnostic_seed(ng_solid)

    return _filaments_from_walked_centerline(
        ng_solid, start_hint, sigma=sigma, nwinc=nwinc, nhinc=nhinc,
        n_peri=n_peri, cad_units_per_meter=cad_units_per_meter,
        solid_bbox_cad=solid_bbox_cad)


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

def _densify_at_corners(spine, init_pts, max_bend_per_step_deg=20.0,
                          max_total_points=500):
    """Insert intermediate spine points until adjacent segments bend
    by no more than ``max_bend_per_step_deg`` at any interior vertex.

    Resamples the spine ON the OCC curve (using ``spine @ t``) so
    inserted points lie exactly on the real spine, not on a polyline
    interpolation of ``init_pts``.

    Used by ``_centerline_from_open_spine`` to smooth visual filament
    bunching at sharp lead-arc corners on "arc + leads" 1-turn coil
    geometries (v4.54.0 response to keiko 2026-05-16 viz report).

    Args:
        spine: OCC edge supporting ``spine @ t`` for t in [0, 1].
        init_pts: (N+1, 3) initial uniform sampling of the spine.
        max_bend_per_step_deg: bend cap per segment transition.
            Default 20 deg keeps cross-section frame rotation gradual
            (3 stations to traverse a 60 deg corner).
        max_total_points: hard cap to prevent infinite blowup on
            spines with cusps or numerical jitter.  Default 500.

    Returns:
        (M+1, 3) densified spine points, M >= N.
    """
    n_init = len(init_pts) - 1
    # Map each sample to its OCC parameter t in [0, 1].  Uniform
    # init samples were taken at t_i = i / n_init.
    t_vals = [i / n_init for i in range(n_init + 1)]
    pts_list = [tuple(p) for p in init_pts]
    cos_max = math.cos(math.radians(max_bend_per_step_deg))

    # Iterative refinement: at each pass, find the worst interior
    # vertex (highest bend angle) and bisect the larger neighbouring
    # segment.  Stop when no vertex exceeds the threshold OR we hit
    # the total-point cap.
    pass_count = 0
    while len(pts_list) < max_total_points and pass_count < max_total_points:
        pass_count += 1
        arr = np.asarray(pts_list)
        seg_v = np.diff(arr, axis=0)
        seg_n = np.linalg.norm(seg_v, axis=1, keepdims=True)
        seg_n = np.maximum(seg_n, 1e-30)
        seg_u = seg_v / seg_n
        cos_bend = np.einsum('ij,ij->i', seg_u[:-1], seg_u[1:])
        cos_bend = np.clip(cos_bend, -1.0, 1.0)
        worst_i = int(np.argmin(cos_bend))  # interior vertex idx +1
        if cos_bend[worst_i] >= cos_max:
            break  # all bends are acceptable

        # Insert ONE intermediate sample on either side of vertex
        # (worst_i + 1) -- bisect the LONGER neighbour to balance
        # segment lengths.
        v_idx = worst_i + 1
        left_len = float(seg_n[worst_i, 0])
        right_len = float(seg_n[worst_i + 1, 0])
        if left_len >= right_len:
            t_a = t_vals[v_idx - 1]
            t_b = t_vals[v_idx]
            insert_pos = v_idx
        else:
            t_a = t_vals[v_idx]
            t_b = t_vals[v_idx + 1]
            insert_pos = v_idx + 1
        t_mid = 0.5 * (t_a + t_b)
        p_mid = spine @ t_mid
        pts_list.insert(insert_pos, (p_mid.X, p_mid.Y, p_mid.Z))
        t_vals.insert(insert_pos, t_mid)

    return np.array(pts_list, dtype=np.float64)


def _parallel_transport_frame(pts: np.ndarray):
    """Compute (tangent, u_hat, v_hat) per vertex of a 3D polyline
    using a Rotation-Minimizing Frame (Wang-Joe double-reflection,
    Wang et al. 2008, "Computation of rotation minimizing frames",
    ACM TOG 27(1):2).

    The double-reflection method minimizes the integral of squared
    angular velocity along the curve -- the resulting frame has the
    smallest possible accumulated twist between adjacent vertices.
    This matters for filament cross-section orientation at sharp
    spine bends: PT (Rodrigues) and RMF agree for smooth curves but
    RMF has better numerical stability + provably-minimum twist on
    polylines with kinks (the v4.53.0 keiko 1turn_coil case).

    Note: visual "bunching" of perimeter filaments at sharp corners
    when viewed from an axis-aligned angle is FORESHORTENING (the
    cross-section plane rotates with the bend; viewed perpendicular
    it projects to a narrower line) -- NOT a twist issue.  RMF and
    PT give the same visual at a 64 deg bend.  To address visual
    bunching, increase spine sampling density near corners or use
    a smooth-curve viz renderer; the frame choice is independent.

    Returns 3 arrays of shape (N, 3) where N = len(pts).

    Function name retained from the v4.13.0 PT implementation for
    backward compat; behaviour changed to RMF in v4.54.0.
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

    # Initial frame at p_0: pick an arbitrary u perpendicular to t_0
    # (avoid alignment with t_0's largest component direction).
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

    # Wang-Joe double-reflection: at each step (p_i, t_i) -> (p_{i+1}, t_{i+1}):
    #   v1 = p_{i+1} - p_i;  c1 = v1.v1
    #   r_L = r_i - (2/c1)*(v1.r_i)*v1     (reflection through plane normal to v1)
    #   t_L = t_i - (2/c1)*(v1.t_i)*v1
    #   v2 = t_{i+1} - t_L;  c2 = v2.v2
    #   r_{i+1} = r_L - (2/c2)*(v2.r_L)*v2 (reflection through plane normal to v2)
    #   s_{i+1} = t_{i+1} x r_{i+1}
    for i in range(n - 1):
        p_i = pts[i]
        p_next = pts[i + 1]
        t_i = tangent[i]
        t_next = tangent[i + 1]
        r_i = u_hat[i]

        v1 = p_next - p_i
        c1 = float(np.dot(v1, v1))
        if c1 < 1e-30:
            u_hat[i + 1] = r_i
        else:
            r_L = r_i - (2.0 / c1) * float(np.dot(v1, r_i)) * v1
            t_L = t_i - (2.0 / c1) * float(np.dot(v1, t_i)) * v1
            v2 = t_next - t_L
            c2 = float(np.dot(v2, v2))
            if c2 < 1e-30:
                u_hat[i + 1] = r_L
            else:
                u_hat[i + 1] = r_L - (2.0 / c2) * float(np.dot(v2, r_L)) * v2

        # Re-orthogonalise (numerical defensiveness; reflections are
        # exact in theory but float rounding can drift over long paths).
        u_hat[i + 1] -= np.dot(u_hat[i + 1], t_next) * t_next
        u_hat[i + 1] /= max(np.linalg.norm(u_hat[i + 1]), 1e-30)
        v_hat[i + 1] = np.cross(t_next, u_hat[i + 1])

    return tangent, u_hat, v_hat


def _check_spine_no_singular_corner(pts_m, radius_m, source_tag):
    """Fail-fast pre-flight: detect spine corners that will produce a
    singular Ruehli L matrix when filament cross-sections are placed by
    parallel transport.

    Background (v4.49.0, 2026-05-16): keiko's
    ``1turn_coil_loft_outsideline.step`` has a 64 deg corner at the
    lead-cap junction where adjacent station spacing (~0.5 mm) is much
    smaller than the wire radius (~2.9 mm).  Parallel transport places
    perimeter filament samples at offset ``radius_m`` perpendicular to
    the local tangent; across a sharp corner whose adjacent segments
    are shorter than the offset, the rotated samples cross over each
    other and adjacent filament segments become near-coincident.  The
    Ruehli mutual-inductance kernel is singular on coincident pairs
    and previously returned NaN/Inf.  v4.48.2 caught it post-assembly
    in ``peec_bundle._assert_solver_L_finite``; v4.49.0 catches it HERE
    BEFORE the O(N^2) Ruehli build, with a far more actionable
    diagnostic (the offending spine vertex coordinate, the bend angle,
    and the segment-length / wire-radius ratio that triggered it).

    Per CLAUDE.md "No Fallbacks - Fail Fast, Fail Loud", this is a
    hard ValueError -- no automatic chamfer insertion, no silent
    re-mesh.  The user fixes the CAD per the FIX hint.

    Detection logic: for each interior spine vertex i (1..N-1),
    compute the bend angle alpha between segments (i-1, i) and
    (i, i+1).  If ``alpha > 60 deg`` AND the minimum adjacent segment
    length is less than the wire radius, the perimeter filaments
    will cross.  Both conditions must hold; either alone is fine
    (a smooth bend of 90 deg with long segments, or a tight short
    spine with no bend, both produce well-conditioned L).

    Threshold rationale:
    - 60 deg: empirical, keiko's failure at 64 deg.  10 deg below the
      observed failure to leave a small safety margin while not
      flagging benign 30-45 deg bends.
    - adj_segment_length < radius_m: at offset ``radius_m * sin(alpha)``
      the perimeter samples shift sideways by O(radius_m); if the
      forward segment length is shorter than that shift, the next
      filament segment lands "behind" the previous one's end ->
      crossing.
    """
    n_pts = len(pts_m)
    if n_pts < 3:
        return
    seg_vec = np.diff(pts_m, axis=0)
    seg_len = np.linalg.norm(seg_vec, axis=1)
    if (seg_len < 1e-12).any():
        bad_i = int(np.argmin(seg_len))
        raise ValueError(
            f"{source_tag}: spine has zero-length segment at station "
            f"{bad_i} -- centerline extractor produced duplicate points.")
    seg_unit = seg_vec / seg_len[:, None]
    cos_bend = np.einsum('ij,ij->i', seg_unit[:-1], seg_unit[1:])
    cos_bend = np.clip(cos_bend, -1.0, 1.0)
    bend_deg = np.degrees(np.arccos(cos_bend))
    adj_min_len = np.minimum(seg_len[:-1], seg_len[1:])
    ratio = adj_min_len / radius_m
    bad = (bend_deg > 60.0) & (ratio < 1.0)
    if not bad.any():
        return
    bad_idx = np.where(bad)[0]
    first = int(bad_idx[0])
    v = pts_m[first + 1]
    raise ValueError(
        f"{source_tag}: spine has a sharp corner that will produce a "
        f"singular PEEC L matrix.  {len(bad_idx)} singular corner(s) "
        f"detected; worst at spine vertex {first + 1}/{n_pts - 1} "
        f"(xyz=({v[0] * 1e3:+.2f}, {v[1] * 1e3:+.2f}, {v[2] * 1e3:+.2f}) mm): "
        f"bend={bend_deg[first]:.1f} deg, adjacent segment length / wire "
        f"radius = {ratio[first]:.3f} (< 1.0 means perimeter filaments "
        f"would cross across the corner under parallel transport).  "
        f"FIX: regenerate the STEP with a smooth single-piece BSPLINE "
        f"lateral so Predicate 1 (UV-map sampling) handles it directly "
        f"without parallel transport.  Alternative: insert a chamfer "
        f"at the corner so adjacent segment length exceeds the wire "
        f"radius.  This check covers both dense-L and HACApK paths "
        f"because it runs BEFORE the kernel build.")


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

    # Pre-flight singular-corner check at the construction layer.
    # Catches keiko-class "sharp corner + short segments" geometries
    # before the O(N^2) Ruehli build; covers both dense-L and HACApK
    # paths because it runs BEFORE solver assembly.
    _check_spine_no_singular_corner(pts, radius_m, source_tag)

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
    from radia.peec_bundle import build_bundle_solver
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
    from radia._b3d_shim import PositionMode

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

    # CCW winding normalisation (v4.53.0, keiko's bug 2026-05-16):
    # per-segment Cubit lofts produce shared cross-section faces with
    # alternating outer-wire orientation -- one volume's "end cap" is
    # the same OCC face as its neighbour's "start cap", but the wire
    # parametrization is reversed.  Without normalisation
    # ``_sample_face_perimeter_in_pt_frame`` returns CW samples at
    # alternating stations, the downstream parallel-bundle solver
    # connects sample k at station i to sample k at station i+1
    # ASSUMING consistent orientation, so adjacent stations zigzag in
    # opposite directions around the cross-section -- the resulting
    # filament paths self-intersect and the Ruehli L matrix
    # degenerates to NaN.  Normalise to CCW (signed area > 0) here
    # so the caller never sees orientation flip.
    signed_area = float(np.sum(uv[:-1, 0] * uv[1:, 1]
                                - uv[:-1, 1] * uv[1:, 0]))
    if signed_area < 0:
        uv = uv[::-1, :].copy()
        u_vals = uv[:, 0]
        v_vals = uv[:, 1]

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

    # For planar (z = 0) circular-arc spines around the z axis (the
    # rect_torus / gapped-torus / loft-around-z geometry family),
    # parallel transport gives u_hat = +z (axial, correct) but
    # v_hat is the chord-perpendicular which is ROTATED relative to
    # the true radial direction at each station by half the
    # angular step (~9 deg at n_stations=20).  The rect cross-
    # section's natural edges align with (+z axial, radial), so PT
    # samples the rect along DIAGONAL axes, shrinking the radial
    # extent by cos(half_step) (~1.5 % at 18 deg/station).
    #
    # Detect the planar-arc case (all centroids at the same z within
    # tolerance) and override v_hat to the true radial direction at
    # each station.  Fixes task #30 rect_torus L = 192 nH -> ~ analytical.
    z_vals = centroids_m[:, 2]
    z_spread = float(np.max(z_vals) - np.min(z_vals))
    radial_xy_min = float(np.min(np.linalg.norm(centroids_m[:, :2], axis=1)))
    if z_spread < 1e-3 * max(radial_xy_min, 1e-30):
        # Planar arc in xy.  Use the geometry-aware frame.
        u_hat = np.tile(np.array([0.0, 0.0, 1.0]), (n_path, 1))
        radial = centroids_m.copy()
        radial[:, 2] = 0.0
        rn = np.linalg.norm(radial, axis=1, keepdims=True)
        rn = np.maximum(rn, 1e-30)
        v_hat = radial / rn

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
    # metres after the caller's cad_units_per_meter scale.  Recover
    # the scale by comparing the station-to-station SPACING in CAD vs
    # in metres -- this is robust to the origin location (the previous
    # `norm(c0_m) / norm(c0_cad)` form silently degenerated to 1.0
    # when c0 happened to lie near the origin, e.g. quarter-symmetry
    # coils with one cap at (0, +y, 0)).  CLAUDE.md "No Fallbacks":
    # we raise instead of silently using cad_to_m=1.0 if the spacing
    # is degenerate too (single-station path).
    if n_path < 2:
        raise ValueError(
            "_filaments_from_per_station_faces: need >= 2 stations to "
            "infer cad-to-m scale, got n_path={}".format(n_path))
    f0c = faces_ordered[0].center()
    f1c = faces_ordered[1].center()
    span_cad = math.sqrt((f1c.X - f0c.X) ** 2
                          + (f1c.Y - f0c.Y) ** 2
                          + (f1c.Z - f0c.Z) ** 2)
    span_m = float(np.linalg.norm(centroids_m[1] - centroids_m[0]))
    if span_cad < 1e-30 or span_m < 1e-30:
        raise ValueError(
            "_filaments_from_per_station_faces: stations 0 and 1 are "
            "coincident (span_cad={:.3e}, span_m={:.3e}); cannot "
            "recover cad-to-m scale".format(span_cad, span_m))
    cad_to_m = span_m / span_cad
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
    from radia.peec_bundle import build_bundle_solver
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
        from radia._b3d_shim import Face
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
    # Spine polyline: CLASSIFIED dispatch (2026-05-02 unified;
    # comment refreshed 2026-07-29 -- there is NO fallback here,
    # extraction failures propagate).  ``extract_coil_topology`` +
    # ``generate_spine`` is OPEN/CLOSED-aware (cap detection +
    # cap-aware arc endpoints); the open-no-clean-cap-pair branch
    # routes through ``_centerline_from_solid_geometry``, which
    # MIRRORS the extract_centerline_from_step classification
    # dispatch on the in-memory solid (it is not an alternative
    # engine).  NOTE: ``generate_spine`` assumes the z rotation axis
    # and only SEEDS the cutting planes -- the filament anchors come
    # from the ACTUAL section-face centroids below, and a non-z coil
    # fails the counted station gate so the dispatch declines to the
    # orientation-agnostic longest-edge tier.
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
    # No-Fallback policy: route OPEN coils through the topology spine
    # (which has cap-aware long-arc handling) and CLOSED coils through
    # the same classification path used by ``extract_centerline_from_step``.
    # Errors propagate; we do NOT try-and-catch from one spine method
    # into another.
    from radia.coil_topology import (
        extract_coil_topology as _extract_topo,
        generate_spine as _gen_spine,
    )
    topo = _extract_topo(solid)
    spine_is_topology_ordered = True
    if topo.is_open and topo.cap_a is not None and topo.cap_b is not None:
        # OPEN swept-cross-section coil with detected caps
        # (rect_torus_lofted_united, gapped torus, etc).  topo has
        # cap_a/cap_b; the spine is the planar long-arc between them
        # at R_spine.  Use the analytical arc spine from
        # ``_gen_spine`` -- ``_centerline_from_open_spine`` (the rim
        # tracer used for coils with straight LEADS extending
        # tangentially out of the plane) returns 7 unevenly-spaced
        # points with a 21.96 mm jump near cap_a on the rect_torus
        # fixture, tripping the spacing-vs-median check below.
        #
        # Note: R_spine from the bbox heuristic may differ from the
        # actual conductor R (e.g. 45.9 vs 50 mm on rect_torus).
        # That's fine -- ``centroids_attempted`` below uses the
        # ACTUAL face centroid (not the spine point), so the spine
        # only seeds the cutting planes, not the filament anchor.
        # Fixes task #30 (2026-05-20).
        path_cad = _gen_spine(topo, n_stations)
        path_m = path_cad / cad_units_per_meter
    elif topo.is_open:
        # OPEN coil with no clean cap pair (e.g. helix with leads):
        # use the rim tracer.
        path_m = _centerline_from_solid_geometry(
            solid, n_segments=n_stations,
            cad_units_per_meter=cad_units_per_meter)
        spine_is_topology_ordered = False  # NN-chain re-order downstream
    else:
        path_cad = _gen_spine(topo, n_stations)
        path_m = path_cad / cad_units_per_meter

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

    # Section at each station.  ``centroids_attempted`` MUST be the
    # face's actual centroid (in metres) per the contract of
    # ``_filaments_from_per_station_faces`` -- that downstream
    # function relies on centroids[i] and faces[i].center() pointing
    # at the SAME physical point (one in m, one in raw CAD units)
    # to recover cad_to_m via ``span_m / span_cad``.  Passing the
    # spine point ``path_m[i]`` here breaks the contract when the
    # spine R differs from the conductor R (e.g. R_spine=45.9 mm
    # from the bbox heuristic vs cap_a R=50 mm on the rect_torus
    # fixture); cad_to_m then degenerates to (R_spine/R_face) and
    # silently shrinks UV by ~8 %.  Always use ``face.center()``
    # for the centroid.  Fixes task #30 (2026-05-20).
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
        fc = face.center()
        face_center_m = np.array(
            [fc.X, fc.Y, fc.Z], dtype=np.float64) / cad_units_per_meter
        faces_attempted.append(face)
        centroids_attempted.append(face_center_m)
        areas_attempted.append(float(face.area))

    if len(faces_attempted) < max(3, int(0.5 * n_path)):
        return None

    # Spine-fits-conductor sanity check.  When the planar-arc spine
    # generated by ``coil_topology.generate_spine`` skips a chunk of
    # the conductor (typical failure: 1-turn coil with straight LEADS
    # extending tangentially out of the loop plane -- the spine arc
    # cuts across the loop INSIDE without sampling the leads), the
    # face centroid jump between cap_a (station 0) and the first
    # interior station is much larger than the typical interior
    # station-to-station spacing.  Detect this and surface a hint
    # instead of producing silently-wrong filaments (keiko / mdx
    # 2026-05-11 report: 'Gmsh にリードが表示されない / filament も
    # リードの形状に追随されてない').
    if (cap_a_face is not None and cap_b_face is not None
            and len(faces_attempted) >= 5):
        face_centroids = np.array(
            [[f.center().X, f.center().Y, f.center().Z]
             for f in faces_attempted], dtype=np.float64)
        spacings = np.linalg.norm(
            face_centroids[1:] - face_centroids[:-1], axis=1)
        # Interior spacings: skip the first and last (cap_a -> first
        # interior, last interior -> cap_b).
        if len(spacings) >= 4:
            interior_spacings = spacings[1:-1]
            interior_median = float(np.median(interior_spacings))
            cap_a_jump = float(spacings[0])
            cap_b_jump = float(spacings[-1])
            for label, jump in (("cap_a", cap_a_jump),
                                  ("cap_b", cap_b_jump)):
                if interior_median > 0 and jump > 2.0 * interior_median:
                    raise ValueError(
                        f"_filaments_from_section_planes: spacing "
                        f"between {label} and the adjacent interior "
                        f"station is {jump * 1e3:.2f} mm vs typical "
                        f"interior spacing "
                        f"{interior_median * 1e3:.2f} mm "
                        f"({jump / interior_median:.1f}x larger) -- "
                        f"the planar-arc spine skips a chunk of the "
                        f"conductor near {label}.  HINT: this "
                        f"typically means the coil has straight "
                        f"LEADS extending tangentially out of the "
                        f"loop plane (1-turn coil with Y-axis leads, "
                        f"etc.) which the current centerline tracer "
                        f"does not sample.  Workarounds: (a) "
                        f"regenerate the .step using a "
                        f"vertex='inner' or 'outer' loft orientation "
                        f"(see gen_*_loft.py -- the vertex=None "
                        f"default produces twisted lateral surfaces "
                        f"that confuse the tracer; keiko/mdx case), "
                        f"or (b) use a Cubit-meshed .vol with "
                        f"--coil-solver bem-a (BEM-A consumes the "
                        f"surface mesh directly and does not need "
                        f"PEEC filament centerline extraction).")

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

    Mirrors the classification dispatch in
    ``extract_centerline_from_step`` (No-Fallback policy, CLAUDE.md
    "No Fallbacks - Fail Fast, Fail Loud"): inspect the solid's
    features and route to ONE extractor.  No try/except cascade -- if
    the chosen extractor fails, the failure propagates.  Returns the
    path in metres (n_segments+1 points).

    Note: the rect / polygon swept-around-z case (lateral surfaces are
    PLANE, not TORUS / CYLINDER) is handled by predicate 3 below
    because such a solid has revolution-style topology even though
    its individual faces are planar -- the
    ``_spine_from_rotation_axis_z`` helper remains available as an
    explicit utility but the main dispatch does not chain it as a
    fallback.
    """
    from radia._b3d_shim import GeomType

    if _collect_loft_cross_sections(solid):
        # _centerline_from_cross_sections returns path only.
        return _centerline_from_cross_sections(
            solid, cad_units_per_meter)[0]

    circle_centers_radius = _collect_circle_edge_centers(solid)
    if circle_centers_radius is not None:
        centers_cad, median_r_cad = circle_centers_radius
        return _centerline_from_circle_edge_centers(
            centers_cad, median_r_cad, cad_units_per_meter)[0]

    has_revolution_face = any(
        f.geom_type in (GeomType.TORUS, GeomType.CYLINDER, GeomType.CONE,
                         GeomType.REVOLUTION)
        for f in solid.faces())
    has_plane_face = any(f.geom_type == GeomType.PLANE for f in solid.faces())
    if has_revolution_face and has_plane_face:
        return _centerline_from_revolution_sweep(
            solid, n_segments, cad_units_per_meter)[0]

    from radia.coil_topology import extract_coil_topology as _extract_topo
    topo = _extract_topo(solid)
    if topo.is_open:
        return _centerline_from_open_spine(
            solid, n_segments, cad_units_per_meter)[0]
    return _centerline_from_topology_spine(
        solid, n_segments, cad_units_per_meter)[0]


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
