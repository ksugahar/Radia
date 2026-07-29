"""coil_from_step.py

Extract coil centerline + per-station profile from an arbitrary STEP
swept-solid via a walking-plane sweep.  Used exclusively by the PEEC
path (``coil_from_cad.filaments_from_step``); the EM / analytical
Biot-Savart path is CoilBuilder-only (`.py` ``build_coil()``), not
STEP.

    STEP file (swept solid with cross-section)
        |
        +-> load_step_solid()
        +-> extract_centerline()
                |
                +-> to_coil_builder()  (intermediate -- NOT called
                                          by any panel directly)
                        then subdivided by coil_from_cad into
                        PEECBuilder filaments (nwinc x nhinc) for
                        skin / proximity effects.

Legacy note (2026-04-25): earlier revisions exposed
``coil_builder_from_step`` and ``coil_builder_from_wire_step`` as
one-shot wrappers for the EM panel.  Both were removed when the EM
panel adopted a CoilBuilder-only policy (user: "AnalyticalCoilは、
CoilBuilder経由で、PEECはstepファイル経由がよいよ").  If you need
CoilBuilder-from-STEP today, author the coil as a `.py`
``build_coil()`` module and call it directly; there is no longer a
STEP -> CoilBuilder path through this module.

Walking-plane algorithm
-----------------------
1. Pick a starting centerline point + tangent (from a small "seed" face,
   or the user-supplied start_pos / start_dir).
2. At each step k:
     a. Cut the solid by an infinite plane perpendicular to the current
        tangent through the current centerline point.
     b. Among the resulting cross-section faces, keep the one closest
        to the previous centerline point (the others are phantom
        sections through distant parts of a closed loop coil).
     c. The face's area-centroid becomes the next centerline point.
     d. Estimate the new tangent from the last few centerline points
        (finite difference + smoothing).
3. Stop when the centerline returns to the start (closed loop) or
   exits the solid (open chain).

Profile classification
----------------------
Each cross-section face is classified into one of:
  - RectProfile  (aspect ratio ~1 with rectangular bbox match)
  - CircleProfile (perimeter^2 / area ~ 4*pi)
  - InterpolatedProfile (general polygon outline as fallback)

The classifier uses the face's mass (=area), perimeter (=outer wire
length) and bounding-box dimensions.
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass

import numpy as np

from radia.coil_profile import (Profile, RectProfile, CircleProfile,
                                InterpolatedProfile)


# ---------- Data containers -------------------------------------------------


@dataclass
class Station:
    """One station along the extracted centerline."""
    s: float                        # arc length from start
    point: np.ndarray               # (3,) world coordinates of centroid
    tangent: np.ndarray             # (3,) unit tangent
    area: float                     # cross-section area
    profile: Profile                # classified Profile object
    polygon: np.ndarray             # (N,2) cross-section outline in local UV


def _pnt_to_np(p):
    """Convert a netgen.occ gp_Pnt (or sequence) to np.ndarray of length 3."""
    return np.array([p[0], p[1], p[2]], dtype=float)


@dataclass
class CenterlineResult:
    """Output of extract_centerline()."""
    polyline: np.ndarray            # (M,3) centerline points
    tangents: np.ndarray            # (M,3) unit tangents
    profiles: list                  # M Profile objects
    polygons: list                  # M (N_k,2) outlines in local UV
    arclen: np.ndarray              # (M,) cumulative arclength
    closed: bool                    # True if returned to start


def scale_centerline_result(res: "CenterlineResult",
                            scale: float) -> "CenterlineResult":
    """Return ``res`` with every geometric quantity multiplied by ``scale``.

    This is the UNIT BOUNDARY of the walking-plane pipeline: the walker
    operates in the STEP file's native CAD units (typically mm), while
    CoilBuilder / PEEC consume metres (CLAUDE.md "Radia always uses
    meters").  Callers pass ``scale = 1 / cad_units_per_meter``
    immediately after ``extract_centerline`` and BEFORE
    ``to_coil_builder``.

    Tangents are direction cosines (scale-free).  Profiles carry
    cross-section dimensions, so each is rebuilt with scaled dimensions;
    an unsupported Profile subclass raises so a new walker profile type
    cannot silently skip the unit conversion.
    """
    try:
        s = float(scale)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"scale_centerline_result: scale must be a finite positive "
            f"number, got {scale!r}") from exc
    if not math.isfinite(s) or s <= 0.0:
        raise ValueError(
            f"scale_centerline_result: scale must be finite and > 0, "
            f"got {scale!r}")
    if s == 1.0:
        return res

    def _scaled_profile(prof):
        if isinstance(prof, CircleProfile):
            return CircleProfile(prof.r * s)
        if isinstance(prof, RectProfile):
            return RectProfile(prof.w * s, prof.h * s)
        raise TypeError(
            f"scale_centerline_result: unsupported Profile type "
            f"{type(prof).__name__}; extend the scaling dispatch when "
            f"the walker learns to emit this profile class.")

    return CenterlineResult(
        polyline=np.asarray(res.polyline, dtype=float) * s,
        tangents=np.asarray(res.tangents, dtype=float).copy(),
        profiles=[_scaled_profile(p) for p in res.profiles],
        polygons=[np.asarray(p, dtype=float) * s for p in res.polygons],
        arclen=np.asarray(res.arclen, dtype=float) * s,
        closed=res.closed,
    )


# ---------- STEP loading ----------------------------------------------------


def _load_step_via_build123d(step_path):
    """Primary STEP reader: build123d (pythonocc-core / XCAF).

    Returns a netgen.occ Solid for the walking-plane pipeline.  The
    transit is `build123d.import_step -> export_brep -> OCCGeometry`:
    build123d reads XCAF (labels, colors, assembly structure) at the
    front, BRep carries lossless geometry to the back-end mesher.

    Multiple solids are fused so the downstream walking-plane sees one
    body.  Returns None if build123d is not installed, so the caller can
    fall through to the legacy netgen.occ path.
    """
    try:
        import build123d as bd
    except ImportError:
        return None
    import tempfile
    from netgen.occ import OCCGeometry
    shape = bd.import_step(step_path)
    # Collect all Solids recursively (Compound.children may be Parts/Solids)
    solids = list(shape.solids())
    if not solids:
        raise ValueError(f"STEP contains no solids: {step_path}")
    if len(solids) == 1:
        merged = solids[0]
    else:
        merged = solids[0]
        for s in solids[1:]:
            merged = merged.fuse(s)
    tmp = tempfile.NamedTemporaryFile(suffix=".brep", delete=False)
    tmp.close()
    try:
        bd.export_brep(merged, tmp.name)
        ng_shape = OCCGeometry(tmp.name).shape
    finally:
        try:
            os.unlink(tmp.name)
        except OSError:
            pass
    ng_solids = list(ng_shape.solids)
    if not ng_solids:
        raise ValueError(
            f"BRep transit produced no solid (empty geometry?) from {step_path}"
        )
    return ng_solids[0] if len(ng_solids) == 1 else ng_shape


def _load_step_via_netgen(step_path):
    """Fallback STEP reader: netgen.occ direct (legacy path, no XCAF).

    Used when build123d is not installed.  Cannot read labels / colors
    because netgen.occ does not expose the XCAF TDF attribute layer.
    """
    from netgen.occ import OCCGeometry
    geo = OCCGeometry(step_path)
    shape = geo.shape
    solids = list(shape.solids)
    if not solids:
        raise ValueError(f"STEP contains no solids: {step_path}")
    if len(solids) == 1:
        return solids[0]
    fused = solids[0]
    for s in solids[1:]:
        fused = fused + s
    return fused


def load_step_solid(step_path):
    """Load a STEP file and return a single Solid (netgen.occ).

    Reader selection:
      1. build123d (pythonocc-core / XCAF) — preferred, preserves label
         and color metadata on the build123d side before BRep transit.
      2. netgen.occ direct — fallback when build123d is not installed.

    Unifying on build123d as the primary STEP reader means that label
    / color auto-detection (see _start_hint_from_step_labels in
    coil_from_cad) and geometry loading use the same parsed document
    state, regardless of which CAD tool produced the STEP.

    Note: in principle the XCAF reading should live in netgen.occ
    itself — there is an open question to Joachim (2026-04-20) about
    exposing XCAF labels on OCCGeometry.  Until that lands, build123d
    is the pragmatic substitute.  See to_developers/ngsolve/2026_04_20
    _to_joachim_glue_step_roundtrip.ipynb for the related thread.

    Multiple solids in the file are fused into one.
    """
    if not os.path.isfile(step_path):
        raise FileNotFoundError(step_path)
    result = _load_step_via_build123d(step_path)
    if result is None:
        result = _load_step_via_netgen(step_path)
    return result


# ---------- Cross-section helpers -------------------------------------------


def _to_axes(point, normal, helper):
    """Build a netgen.occ Axes from numpy arrays."""
    from netgen.occ import Axes, Pnt, Dir
    return Axes(p=Pnt(*point), n=Dir(*normal), h=Dir(*helper))


def _pick_helper(normal):
    """Pick an in-plane heading vector orthogonal to normal."""
    n = np.asarray(normal, dtype=float)
    n /= np.linalg.norm(n)
    # Use whichever world axis is least parallel to n.
    pick = np.argmin(np.abs(n))
    cand = np.zeros(3); cand[pick] = 1.0
    h = cand - np.dot(cand, n) * n
    h /= np.linalg.norm(h)
    return h


def _slab_disk(point, normal, in_plane_radius, thickness):
    """Build a thin disk (Cylinder) centered at point, normal-aligned.

    Used as a boolean-intersection probe for cross-section extraction.
    Works for arbitrary normal orientation (Cylinder takes axis_dir).
    """
    from netgen.occ import Cylinder, Pnt, Dir
    n = np.asarray(normal, dtype=float)
    n /= np.linalg.norm(n)
    p = np.asarray(point, dtype=float)
    base = p - 0.5 * thickness * n
    return Cylinder(Pnt(*base), Dir(*n), in_plane_radius, thickness)


def cross_section_pieces(solid, point, normal, in_plane_radius=None,
                         thickness=None):
    """Slice the solid with a thin slab and return the connected pieces.

    Each piece is reported as ``(world_centroid, area_estimate, sub_solid)``.
    Area estimate = piece.mass / thickness.

    Args:
        in_plane_radius: half-extent of the slab in the cutting plane.
            Defaults to the solid's bbox diagonal.
        thickness: slab thickness along the normal. Defaults to
            ``in_plane_radius / 1000``.
    """
    if in_plane_radius is None:
        bb = solid.bounding_box
        diag = math.sqrt(sum((bb[1][i] - bb[0][i]) ** 2 for i in range(3)))
        in_plane_radius = max(diag, 1e-3)
    if thickness is None:
        thickness = max(in_plane_radius / 1000.0, 1e-7)

    slab = _slab_disk(point, normal, in_plane_radius, thickness)
    inter = solid * slab
    pieces = []
    for sub in inter.solids:
        c = _pnt_to_np(sub.center)
        m = sub.mass
        pieces.append((c, m / thickness, sub))
    return pieces


def pick_local_piece(pieces, prev_point, max_dist=None):
    """Among pieces, return the one whose centroid is closest to prev.

    Returns (centroid, area, sub_solid) or (None, None, None).
    """
    if not pieces:
        return None, None, None
    prev = np.asarray(prev_point, dtype=float)
    best = None
    best_d = math.inf
    for c, a, sub in pieces:
        d = np.linalg.norm(c - prev)
        if d < best_d:
            best, best_d = (c, a, sub), d
    if max_dist is not None and best_d > max_dist:
        return None, None, None
    return best


# ---------- Profile classification ------------------------------------------


def _piece_bbox_uv(sub_solid, point, normal):
    """Return (u_min, u_max, v_min, v_max) of sub_solid in plane local UV.

    The 8 corners of sub.bounding_box (world AABB) are projected onto
    the cutting plane's (u, v) axes; min/max in each direction is the
    plane-local bbox of the cross-section.
    """
    n = np.asarray(normal, dtype=float); n /= np.linalg.norm(n)
    u_axis = _pick_helper(n)
    v_axis = np.cross(n, u_axis)
    origin = np.asarray(point, dtype=float)
    bb = sub_solid.bounding_box
    lo = np.array([bb[0][0], bb[0][1], bb[0][2]])
    hi = np.array([bb[1][0], bb[1][1], bb[1][2]])
    corners = np.array([[x, y, z]
                        for x in (lo[0], hi[0])
                        for y in (lo[1], hi[1])
                        for z in (lo[2], hi[2])])
    rel = corners - origin
    us = rel @ u_axis
    vs = rel @ v_axis
    return float(us.min()), float(us.max()), float(vs.min()), float(vs.max())


def classify_profile(area, bbox_uv):
    """Classify a cross-section by its area + plane-local bbox.

    bbox_uv = (u_min, u_max, v_min, v_max).

    Heuristics:
      * area == bbox_area  -> RectProfile (exact fill)
      * area / bbox_area ≈ pi/4 (with bbox aspect ~1) -> CircleProfile
      * else -> RectProfile fallback (TODO: InterpolatedProfile when
        we can densely sample the boundary).
    """
    u_min, u_max, v_min, v_max = bbox_uv
    bw = u_max - u_min
    bh = v_max - v_min
    bbox_area = bw * bh
    if area <= 0 or bw <= 1e-9 or bh <= 1e-9 or bbox_area <= 0:
        return RectProfile(max(bw, 1e-6), max(bh, 1e-6))

    fill = area / bbox_area
    aspect = max(bw, bh) / min(bw, bh)

    # Exact rectangle fill
    if abs(fill - 1.0) < 0.05:
        return RectProfile(bw, bh)

    # Circle / ellipse: pi/4 = 0.785 fill, aspect ~1
    if abs(fill - math.pi / 4) < 0.05 and abs(aspect - 1.0) < 0.10:
        r = math.sqrt(area / math.pi)
        return CircleProfile(r)

    # Fallback: rect approximation
    return RectProfile(bw, bh)


# ---------- Centerline walking ----------------------------------------------


def _initial_seed(solid, start_hint=None):
    """Pick a starting point + tangent for the walk.

    If start_hint=(point, dir) is supplied, use it. Otherwise use the
    smallest-area face (typically a port end-cap on an open coil); for
    a closed loop, the user must supply a hint.
    """
    if start_hint is not None:
        p, t = start_hint
        return np.asarray(p, dtype=float), np.asarray(t, dtype=float)

    smallest = None
    smallest_area = math.inf
    for f in solid.faces:
        a = f.mass
        if a < smallest_area and a > 0:
            smallest_area = a
            smallest = f
    if smallest is None:
        raise RuntimeError("no usable seed face found")
    p = _pnt_to_np(smallest.center)
    # Tangent = face normal (approximate, since end-cap normal == axis)
    # We approximate by area weighted normal via integrating outward;
    # netgen.occ doesn't expose face normal directly, so fall back to
    # vector from solid centroid to face centroid (works for end caps).
    sc = _pnt_to_np(solid.center)
    t = p - sc
    n = np.linalg.norm(t)
    if n < 1e-12:
        raise RuntimeError("seed tangent ill-defined; supply start_hint")
    t /= n
    return p, t


def _axis_agnostic_seed(solid):
    """Orientation-agnostic (point, tangent) walk seed for a coil solid.

    Replaces the legacy z-axis-torus heuristic (start on +x, tangent +y)
    that the default filament path forced -- overriding ``_initial_seed``'s
    own axis-agnostic smallest-face seed -- and so broke on x/y-axis coils.

    Strategy = the minimum-cross-section-area principle.  A thin slab taken
    PERPENDICULAR to the local wire direction yields the smallest single
    connected cross-section, so the section normal that minimises the piece
    area at a point on the coil IS the wire tangent there, and the winning
    piece centroid is a point on the wire.  This is independent of the
    model's principal orientation and works for open or closed loops
    (no port-cap face required).
    """
    bb = solid.bounding_box
    lo = np.array([bb[0][0], bb[0][1], bb[0][2]], dtype=float)
    hi = np.array([bb[1][0], bb[1][1], bb[1][2]], dtype=float)
    center = 0.5 * (lo + hi)
    ext = hi - lo
    diag = float(np.linalg.norm(ext))
    slab_r = max(diag, 1e-3)
    thick = max(diag / 5000.0, 1e-7)
    axes = np.eye(3)

    # Step 1: a point ON the wire.  Section through the bbox centre with the
    # longest-extent axis as the plane normal (that plane crosses the loop
    # for any orientation); take the first connected piece centroid.
    p_on = None
    for ai in list(np.argsort(-ext)):
        if ext[ai] <= 0:
            continue
        pieces = cross_section_pieces(solid, center, axes[ai],
                                      in_plane_radius=slab_r, thickness=thick)
        if pieces:
            p_on = np.asarray(pieces[0][0], dtype=float)
            break
    if p_on is None:
        raise RuntimeError(
            "_axis_agnostic_seed: no cross-section through the bbox centre "
            "along any axis; the solid may not be a coil -- supply "
            "start_hint explicitly.")

    # Step 2: wire tangent at p_on = the axis normal that MINIMISES the
    # single-piece section area (smallest cut == perpendicular to the wire).
    best = None  # (area, centroid, normal)
    for ti in range(3):
        pieces = cross_section_pieces(solid, p_on, axes[ti],
                                      in_plane_radius=slab_r, thickness=thick)
        c, area, _sub = pick_local_piece(pieces, p_on, max_dist=slab_r)
        if c is None or area is None or area <= 0:
            continue
        if best is None or area < best[0]:
            best = (area, np.asarray(c, dtype=float), axes[ti].copy())
    if best is None:
        raise RuntimeError(
            "_axis_agnostic_seed: no valid tangent section at the wire "
            "point; supply start_hint explicitly.")
    _area, p, t = best
    return p, t


def extract_centerline(step_path_or_solid, *,
                       start_hint=None,
                       step_size=None,
                       max_stations=2000,
                       close_tol=None,
                       verbose=False):
    """Walk a coil solid and return its centerline + per-station profiles.

    Args:
        step_path_or_solid: STEP filepath, or a netgen.occ Solid.
        start_hint: optional ((px,py,pz), (tx,ty,tz)) seed.
        step_size: walking step [m]. Defaults to the cube root of
            (solid.mass / 100), giving ~100 stations for a typical coil.
        max_stations: hard cap to prevent runaway loops.
        close_tol: distance below which the walk is considered to have
            closed back to its start point. Defaults to 2*step_size.
        verbose: print per-step progress.

    Returns:
        CenterlineResult.
    """
    if isinstance(step_path_or_solid, str):
        solid = load_step_solid(step_path_or_solid)
    else:
        solid = step_path_or_solid

    if step_size is None:
        step_size = (solid.mass / 100.0) ** (1.0 / 3.0)
    if close_tol is None:
        close_tol = 2.0 * step_size

    p, t = _initial_seed(solid, start_hint)
    points, tangents, profiles, polygons, areas = [], [], [], [], []
    closed = False

    bb = solid.bounding_box
    diag = math.sqrt(sum((bb[1][i] - bb[0][i]) ** 2 for i in range(3)))
    slab_radius = max(diag, 1e-3)
    slab_thickness = max(step_size / 50.0, 1e-7)

    def _try_at(p_try, t_try, max_dist):
        pieces = cross_section_pieces(
            solid, p_try, t_try,
            in_plane_radius=slab_radius, thickness=slab_thickness)
        return pick_local_piece(pieces, p_try, max_dist=max_dist)

    # Initial probe at the seed point
    c, area, sub = _try_at(p, t, max_dist=10 * step_size)
    if c is None:
        raise RuntimeError("seed plane has no usable cross-section")
    bbox = _piece_bbox_uv(sub, c, t)
    prof = classify_profile(area, bbox)
    points.append(c); tangents.append(t.copy())
    profiles.append(prof); polygons.append(np.array(bbox)); areas.append(area)

    for k in range(1, max_stations):
        # Adaptive step: try full step, on miss halve up to 6 times
        success = False
        s = step_size
        for retry in range(7):
            p_try = points[-1] + s * t
            c, area, sub = _try_at(p_try, t, max_dist=4 * s)
            if c is not None:
                success = True
                break
            s *= 0.5
        if not success:
            if verbose:
                print(f"  step {k}: walk halted (no valid section)")
            break

        bbox = _piece_bbox_uv(sub, c, t)
        prof = classify_profile(area, bbox)
        points.append(c); tangents.append(t.copy())
        profiles.append(prof)
        polygons.append(np.array(bbox)); areas.append(area)

        # Update tangent from last 2 centroids (Frenet first-derivative)
        t_new = points[-1] - points[-2]
        nrm = np.linalg.norm(t_new)
        if nrm > 1e-12:
            t = t_new / nrm

        # Closure check (after a few steps to avoid trivial closure)
        if k > 5 and np.linalg.norm(c - points[0]) < close_tol:
            if verbose:
                print(f"  step {k}: closed loop")
            closed = True
            break

        if verbose and k % 10 == 0:
            print(f"  step {k:3d}: p={c}, area={area:.4e}, used s={s:.4e}")

    if not points:
        raise RuntimeError("no stations extracted")

    pts = np.array(points)
    tans = np.array(tangents)
    arclen = np.zeros(len(pts))
    for i in range(1, len(pts)):
        arclen[i] = arclen[i - 1] + np.linalg.norm(pts[i] - pts[i - 1])

    return CenterlineResult(
        polyline=pts,
        tangents=tans,
        profiles=profiles,
        polygons=polygons,
        arclen=arclen,
        closed=locals().get('closed', False),
    )


# ---------- Sprint 2: polyline -> CoilBuilder segments ----------------------


@dataclass
class SegmentSpec:
    """Geometric specification of one reconstructed segment."""
    kind: str               # 'straight' or 'arc'
    length: float           # arc length [m]
    radius: float = 0.0     # arc radius [m] (0 for straight)
    angle_deg: float = 0.0  # arc angle [deg] (0 for straight)
    profile: Profile = None # cross-section profile
    start_pos: np.ndarray = None      # (3,) world starting point
    start_tangent: np.ndarray = None  # (3,) unit tangent
    start_normal: np.ndarray = None   # (3,) unit cross-section normal
    arc_center: np.ndarray = None     # (3,) arc center (arcs only)
    arc_normal: np.ndarray = None     # (3,) arc plane normal (arcs only)


def _three_point_circle(p0, p1, p2, plane_normal):
    """Return (center, radius) of the circle through three 3D points.

    Uses perpendicular bisectors in the plane defined by ``plane_normal``.
    Returns (None, None) if the three points are collinear or the
    intersection is ill-conditioned.
    """
    n = plane_normal / np.linalg.norm(plane_normal)
    a = p1 - p0
    b = p2 - p1
    # Perpendicular bisectors (in-plane directions)
    da = np.cross(n, a)
    db = np.cross(n, b)
    ma = 0.5 * (p0 + p1)
    mb = 0.5 * (p1 + p2)
    # Solve ma + t*da = mb + s*db for t  -> 2x2 system in the plane.
    # Set up: t*da - s*db = mb - ma
    A = np.column_stack([da, -db])  # (3, 2)
    rhs = mb - ma
    # Least-squares solution (plane-restricted, so rank 2)
    try:
        ts, *_ = np.linalg.lstsq(A, rhs, rcond=None)
    except np.linalg.LinAlgError:
        return None, None
    t = ts[0]
    center = ma + t * da
    radius = float(np.linalg.norm(p0 - center))
    return center, radius


def _smooth_curvature(polyline, closed):
    """Compute discrete curvature kappa_i [1/m] at each polyline station.

    Uses three consecutive points (p_{i-1}, p_i, p_{i+1}) and the
    Menger curvature formula:  kappa = 4 * area / (|a| * |b| * |c|).

    For closed loops the indices wrap; for open chains the endpoints
    use a one-sided estimate (set to neighbor's kappa).
    """
    n = len(polyline)
    if n < 3:
        return np.zeros(n)
    kappa = np.zeros(n)
    for i in range(n):
        if not closed and (i == 0 or i == n - 1):
            continue
        i0 = (i - 1) % n
        i1 = i
        i2 = (i + 1) % n
        a = polyline[i1] - polyline[i0]
        b = polyline[i2] - polyline[i1]
        c = polyline[i2] - polyline[i0]
        an, bn, cn = (np.linalg.norm(a), np.linalg.norm(b),
                      np.linalg.norm(c))
        if an < 1e-12 or bn < 1e-12 or cn < 1e-12:
            continue
        # Triangle area via cross product
        area = 0.5 * np.linalg.norm(np.cross(a, b))
        kappa[i] = 4.0 * area / (an * bn * cn)
    if not closed:
        kappa[0] = kappa[1]
        kappa[-1] = kappa[-2]
    return kappa


def _estimate_plane_normal(polyline):
    """Best-fit plane normal (SVD) for a roughly planar polyline.

    Used to disambiguate the cross-section's normal direction for
    arc segments. For non-planar coils the result is the principal
    out-of-plane direction.
    """
    centroid = polyline.mean(axis=0)
    pts = polyline - centroid
    # SVD: smallest singular value's right-singular vector = plane normal
    _, _, Vt = np.linalg.svd(pts, full_matrices=False)
    return Vt[-1] / np.linalg.norm(Vt[-1])


def polyline_to_segments(result, *,
                         straight_kappa_threshold=2.0,
                         curvature_var_tol=0.20):
    """Group polyline stations into straight + arc segments.

    Args:
        result: CenterlineResult from extract_centerline().
        straight_kappa_threshold: kappa < threshold/avg_radius -> straight
            (default: stations with kappa less than 1/(50 * avg_radius)).
        curvature_var_tol: relative variation of kappa within an arc
            group must be < this fraction.

    Returns:
        list of SegmentSpec.
    """
    pts = result.polyline
    tans = result.tangents
    profiles = result.profiles
    arclen = result.arclen
    closed = result.closed
    n = len(pts)
    if n < 3:
        return []

    plane_normal = _estimate_plane_normal(pts)
    kappa = _smooth_curvature(pts, closed)
    avg_R = np.mean(arclen[1:] - arclen[:-1]) * n / (2 * np.pi) \
        if closed else float(np.linalg.norm(pts[-1] - pts[0]))
    avg_R = max(avg_R, 1e-6)
    kappa_straight = straight_kappa_threshold / (50.0 * avg_R)

    # Classify each station
    kinds = ['straight' if k < kappa_straight else 'arc' for k in kappa]

    # Group consecutive same-kind stations (with similar kappa for arcs)
    groups = []  # list of (kind, [indices])
    cur_kind = kinds[0]
    cur_idx = [0]
    cur_kappa_mean = kappa[0]
    for i in range(1, n):
        if kinds[i] == cur_kind:
            if cur_kind == 'arc':
                # Check kappa stability within group
                grp_mean = np.mean([kappa[j] for j in cur_idx])
                if grp_mean > 0 and abs(kappa[i] - grp_mean) / grp_mean \
                        < curvature_var_tol:
                    cur_idx.append(i)
                    continue
                # Curvature changed too much -> close group and restart
                groups.append((cur_kind, cur_idx))
                cur_kind = 'arc'
                cur_idx = [i]
                cur_kappa_mean = kappa[i]
                continue
            cur_idx.append(i)
        else:
            groups.append((cur_kind, cur_idx))
            cur_kind = kinds[i]
            cur_idx = [i]
            cur_kappa_mean = kappa[i]
    groups.append((cur_kind, cur_idx))

    # Merge small noise groups into a single dominant arc when the
    # polyline is "mostly arc" (e.g., a torus centerline where the
    # first and last few stations have inflated curvature due to the
    # wrap-around gap).  The `_smooth_curvature` endpoint artifacts
    # otherwise split what is geometrically a single circular arc
    # into multiple groups, truncating the reconstructed arc angle.
    n_arc_pts = sum(len(idx) for k, idx in groups if k == 'arc')
    if (n_arc_pts >= 0.7 * n
            and all(k == 'arc' for k, _ in groups)):
        # Find the dominant group (largest) and adopt its kappa as
        # the canonical radius.  Replace everything with one group
        # spanning [0..n-1].
        dom_kind, dom_idx = max(groups, key=lambda g: len(g[1]))
        dom_kappa = float(np.mean([kappa[j] for j in dom_idx]))
        groups = [('arc', list(range(n)))]
        # Override kappa[i] for outliers so the downstream radius
        # reconstruction isn't biased by endpoint noise.
        for i in range(n):
            if dom_kappa > 0 and abs(kappa[i] - dom_kappa) / dom_kappa > curvature_var_tol:
                kappa[i] = dom_kappa

    # Convert groups to SegmentSpec
    segs = []
    for kind, idx in groups:
        if len(idx) < 2:
            continue
        i0, i1 = idx[0], idx[-1]
        # Use full polyline arclength for a closed loop spanning all
        # stations (arclen[-1] is distance of pts[n-1] from pts[0], so
        # we must add the wrap-around chord to close the loop).
        if closed and i0 == 0 and i1 == n - 1:
            length = arclen[i1] + np.linalg.norm(pts[0] - pts[-1])
        elif closed and i1 == n - 1:
            length = arclen[i1] - arclen[i0] + np.linalg.norm(pts[0] - pts[-1])
        else:
            length = arclen[i1] - arclen[i0]
        prof = profiles[i0]
        start_pos = pts[i0].copy()
        start_tan = tans[i0].copy()
        if kind == 'straight':
            segs.append(SegmentSpec(
                kind='straight', length=length, profile=prof,
                start_pos=start_pos, start_tangent=start_tan))
        else:
            kappa_mean = float(np.mean([kappa[j] for j in idx]))
            radius = 1.0 / kappa_mean if kappa_mean > 0 else 0.0
            angle_rad = length * kappa_mean
            # 3-point circle fit over well-spaced polyline samples
            # inside this arc group to recover the arc center.
            # For a closed loop the first and last indices are
            # adjacent in world space, so picking j0=idx[0] and
            # j2=idx[-1] gives a nearly degenerate triangle and an
            # ill-conditioned center.  Use the 1/4, 1/2, 3/4
            # positions instead (well separated around the arc).
            arc_center = None
            if len(idx) >= 4:
                m = len(idx)
                j0 = idx[m // 4]
                j1 = idx[m // 2]
                j2 = idx[(3 * m) // 4]
                c_fit, r_fit = _three_point_circle(
                    pts[j0], pts[j1], pts[j2], plane_normal)
                if c_fit is not None:
                    arc_center = c_fit
                    # Prefer the geometric radius over the
                    # kappa-averaged one — the 3-point fit is exact
                    # for circular arcs, kappa averaging isn't.
                    radius = r_fit
                    # Re-derive angle from the exact radius so the
                    # reconstructed arc spans the correct sweep.
                    angle_rad = length / radius
            segs.append(SegmentSpec(
                kind='arc', length=length, radius=radius,
                angle_deg=math.degrees(angle_rad), profile=prof,
                start_pos=start_pos, start_tangent=start_tan,
                arc_center=arc_center,
                arc_normal=plane_normal.copy()))
    return segs


def to_coil_builder(result, current=1.0):
    """Build a CoilBuilder from an extracted centerline result.

    Currently supports planar coils with straight + circular-arc
    segments (Sprint 2 scope). Lofts (varying cross-section) and
    helical / non-planar arcs are NOT yet handled.
    """
    from radia.coil_builder import CoilBuilder
    segs = polyline_to_segments(result)
    if not segs:
        raise RuntimeError("no segments reconstructed")

    builder = CoilBuilder(current=current)
    first = segs[0]
    # CoilBuilder orientation convention (rows = local X, Y, Z axes):
    #   row 1 (Y) = path tangent (StraightSegment extends in +Y,
    #               ArcSegment tangent at start is +Y)
    #   row 0 (X) = radial outward from arc_center toward start_pos
    #               (ArcSegment: arc_center = start_pos - radius * row0)
    #               For a straight lead-in the X axis is the cross-section
    #               width direction; any choice orthogonal to tangent works.
    #   row 2 (Z) = arc plane normal (cross-section height direction)
    plane_n = _estimate_plane_normal(result.polyline)
    y_axis = first.start_tangent / np.linalg.norm(first.start_tangent)
    if first.kind == 'arc' and first.arc_center is not None:
        # Radial from arc center to start position — this is the
        # authoritative X-axis for the ArcSegment convention.
        x_axis = first.start_pos - first.arc_center
        x_axis = x_axis / np.linalg.norm(x_axis)
        z_axis = np.cross(x_axis, y_axis)
        z_axis = z_axis / np.linalg.norm(z_axis)
        # Re-orthogonalize Y against (X, Z) to absorb tangent estimation
        # noise while keeping X radial-exact (needed for arc_center
        # round-trip correctness).
        y_axis = np.cross(z_axis, x_axis)
    else:
        z_axis = plane_n / np.linalg.norm(plane_n)
        z_axis = z_axis - np.dot(z_axis, y_axis) * y_axis
        z_axis = z_axis / np.linalg.norm(z_axis)
        x_axis = np.cross(y_axis, z_axis)
    orient = np.array([x_axis, y_axis, z_axis])  # row vectors
    builder.set_start(first.start_pos.tolist(), orientation=orient)
    builder.set_profile(first.profile)
    for s in segs:
        if s.profile is not None:
            # Refresh profile in case it changed across segment boundary
            builder.set_profile(s.profile)
        if s.kind == 'straight':
            builder.add_straight(length=s.length)
        else:
            builder.add_arc(radius=s.radius, arc_angle=s.angle_deg)
    return builder, segs


# coil_builder_from_step / coil_builder_from_wire_step were removed
# 2026-04-25.  The EM panel is CoilBuilder-only (.py); the PEEC
# pipeline uses `extract_centerline` + `to_coil_builder` internally
# via `coil_from_cad.filaments_from_step`.  There is no supported
# STEP -> CoilBuilder entry point in this module anymore.
