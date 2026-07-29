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
    # (M,) measured cross-section areas (piece.mass / slab thickness --
    # the volume-based TRUE section area, unlike Profile.total_area()
    # which inherits the world-AABB inflation of _piece_bbox_uv on
    # tilted sections).  None only for hand-built results.
    areas: np.ndarray = None


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
        areas=(None if res.areas is None
               else np.asarray(res.areas, dtype=float) * (s * s)),
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


def _fan_directions(t_in, polar_degs=(45.0, 90.0, 135.0), n_azimuth=8):
    """Unit directions on cones around ``t_in`` (gentlest polar first)."""
    t_in = np.asarray(t_in, dtype=float)
    t_in = t_in / np.linalg.norm(t_in)
    u = _pick_helper(t_in)
    v = np.cross(t_in, u)
    out = []
    for polar_deg in polar_degs:
        ct = math.cos(math.radians(polar_deg))
        st = math.sin(math.radians(polar_deg))
        for j in range(n_azimuth):
            phi = 2.0 * math.pi * j / n_azimuth
            d = ct * t_in + st * (math.cos(phi) * u + math.sin(phi) * v)
            out.append(d / np.linalg.norm(d))
    return out


def _turn_search(solid, p_last, t_in, recent_points, recent_areas,
                 recent_wh_max, *, step_size, slab_radius,
                 slab_thickness, extra_guard_pts=()):
    """Search for the outgoing wire direction at a sharp spine corner.

    Called when the march exhausts its step-halving without finding a
    section: the Frenet tangent cannot bend around a corner sharper
    than the halving can follow (e.g. the 90-deg legs of a rectangular
    racetrack).  Near such a corner, sections through the corner
    region MERGE with the orthogonal member (an L-shaped or full-width
    cut whose centroid is far / whose area is discontinuous), so a
    single coupled probe ``p_last + s * d`` cannot both clear the
    corner block AND cut the outgoing member cleanly.  Two decoupled
    stages instead:

    1. ANCHOR: a point on/near the outgoing member.
       (a) PRIMARY: the centroid of the very merged section that
           halted the march (global slab perpendicular to the incoming
           tangent, one step ahead).  The corner cut is CONNECTED to
           the outgoing member, so its centroid is pulled onto it (for
           a plate frame it lands exactly on the outgoing member's
           mid-plane).  Verified to be material by a local-disk probe.
       (b) FALLBACK GRID: a fan of directions (cones at 45/90/135 deg,
           gentlest first) at reaches 1..3 x max(step, wire size),
           each gated by the area-continuity window so an off-solid
           disk snagging an outer-edge sliver cannot become an anchor.
    2. TANGENT AT THE ANCHOR: minimum-section-area direction fan AT
       the anchor -- the ``_axis_agnostic_seed`` principle: the true
       wire tangent minimises the local section area.  The sign is
       resolved toward progress from ``p_last``.

    A straight outgoing member may be entered mid-member (the (a)
    anchor is the merged-cut centroid): the polyline then represents
    the skipped stretch by its chord, which coincides with the
    centerline for straight members; the following march covers the
    rest.

    Guards (all fail-closed; no survivor = genuine end, walk halts):
      * area continuity: the wire cross-section is continuous through
        a corner -- accepted area must be 0.6..1.8x the recent median
        (rejects slivers and merged cuts);
      * anti-revisit: the new station must not land on the
        already-walked centerline (within 0.7 * step of ANY station of
        this march except the seed region and the last two, plus any
        ``extra_guard_pts`` from a sibling march) -- rejects
        back-turns into the incoming wire at a true open end AND
        wrong-branch turns at a junction onto an already-walked leg
        (which would multi-lap a closed circuit).  The first 3
        stations are exempt so the FINAL corner of a closed loop may
        legitimately turn back toward the seed and let the closure
        check fire;
      * midpoint continuity: local material must exist halfway from
        ``p_last`` to the new station (rejects hops across an air gap
        onto a nearby but disconnected passage, e.g. an adjacent
        turn of a tight winding);
      * minimum motion: |c* - p_last| >= 0.3 * step.

    A ZERO-RADIUS 180-deg hairpin remains out of scope: its re-entry
    lands on the walked centerline and anti-revisit rejects it -- the
    walk halts and the caller's coverage checks fail loud.

    Returns ``(centroid, area, sub_solid, direction)`` or
    ``(None, None, None, None)``.
    """
    t_in = np.asarray(t_in, dtype=float)
    t_in = t_in / np.linalg.norm(t_in)
    p_last = np.asarray(p_last, dtype=float)

    med_area = float(np.median(np.asarray(recent_areas, dtype=float)))
    guard_pts = [np.asarray(q, dtype=float)
                 for q in list(recent_points[3:-2]) + list(extra_guard_pts)]

    reach_unit = max(float(step_size), float(recent_wh_max))
    rho = 2.0 * max(float(recent_wh_max), float(step_size))

    def _local_piece(point, normal, radius):
        pieces = cross_section_pieces(
            solid, point, normal,
            in_plane_radius=radius, thickness=slab_thickness)
        return pick_local_piece(pieces, point, max_dist=0.8 * radius)

    def _anchors():
        # (a) merged-forward-section centroid: the corner cut that
        # halted the march is connected to the outgoing member and its
        # centroid is pulled onto it.  Verify it is material.
        p_fwd = p_last + step_size * t_in
        pieces = cross_section_pieces(
            solid, p_fwd, t_in,
            in_plane_radius=slab_radius, thickness=slab_thickness)
        c_a, _a_a, _s_a = pick_local_piece(pieces, p_fwd, max_dist=None)
        if c_a is not None:
            c_a = np.asarray(c_a, dtype=float)
            c_chk, _chk_a, _chk_s = _local_piece(c_a, t_in, rho)
            if c_chk is not None:
                yield c_a
        # (b) direction-fan grid, gated by the area-continuity window
        # so an off-solid disk snagging an outer-edge sliver cannot
        # become an anchor.
        for d1 in _fan_directions(t_in):
            for reach in (1.0, 2.0, 3.0):
                q = p_last + reach * reach_unit * d1
                c_q, a_q, _sub_q = _local_piece(q, d1, rho)
                if c_q is None or a_q is None or a_q <= 0:
                    continue
                if med_area > 0 and not (0.6 * med_area <= a_q
                                          <= 1.8 * med_area):
                    continue
                yield np.asarray(c_q, dtype=float)
                break  # one anchor per direction

    stage2_budget = 6
    for anchor in _anchors():
        if stage2_budget <= 0:
            break
        stage2_budget -= 1

        # Stage 2: min-area tangent at the anchor.
        best = None  # (area, centroid, sub, direction)
        for d2 in _fan_directions(t_in):
            c2, a2, sub2 = _local_piece(anchor, d2, 1.25 * rho)
            if c2 is None or a2 is None or a2 <= 0:
                continue
            if med_area > 0 and not (0.6 * med_area <= a2
                                      <= 1.8 * med_area):
                continue
            if best is None or a2 < best[0]:
                best = (a2, c2, sub2, d2)
                # Early accept only for a near-perfect section: an
                # oblique cut at 30 deg is already 1.15x, so 1.05x
                # admits only the true perpendicular (a
                # 22.5-deg-misaligned-azimuth worst case falls through
                # to the full-fan minimum).
                if a2 <= 1.05 * med_area:
                    break
        if best is None:
            continue  # material anchor, but no clean section here
        a_star, c_star, sub_star, d_star = best

        # Resolve tangent sign toward progress from p_last.
        if float(np.dot(c_star - p_last, d_star)) < 0:
            d_star = -d_star
        if float(np.linalg.norm(c_star - p_last)) < 0.3 * step_size:
            continue
        if any(np.linalg.norm(c_star - g) < 0.7 * step_size
               for g in guard_pts):
            continue  # lands on the already-walked centerline
        d_mid = c_star - p_last
        c_m, _a_m, _s_m = _local_piece(
            0.5 * (p_last + c_star),
            d_mid / np.linalg.norm(d_mid), rho)
        if c_m is None:
            continue  # air gap mid-bend: disconnected passage
        return c_star, a_star, sub_star, d_star
    return None, None, None, None


def _march(solid, p0, t0, *, step_size, max_stations, close_tol,
           slab_radius, slab_thickness, verbose=False, max_turns=64,
           extra_guard_pts=()):
    """One directional walking-plane march from seed ``(p0, t0)``.

    The first station is the seed section itself.  Sharp spine corners
    are turned via ``_turn_search`` (from step 2 on -- a halt on the
    very first step means the seed tangent points off the wire end,
    which the caller's BACKWARD march covers; a turn there would
    re-walk the same wire and double-cover it).

    Returns ``(points, tangents, profiles, polygons, areas, closed,
    turns)``.  Raises ``RuntimeError`` if the seed plane has no usable
    section.
    """
    def _try_at(p_try, t_try, max_dist):
        pieces = cross_section_pieces(
            solid, p_try, t_try,
            in_plane_radius=slab_radius, thickness=slab_thickness)
        return pick_local_piece(pieces, p_try, max_dist=max_dist)

    p = np.asarray(p0, dtype=float)
    t = np.asarray(t0, dtype=float)
    t = t / np.linalg.norm(t)

    guard_extra = (np.asarray(list(extra_guard_pts), dtype=float)
                   if len(list(extra_guard_pts)) else None)

    def _revisits(c_cand):
        """True if c_cand lands on a NON-recent already-walked station.

        pick_local_piece can snap a normal step across a nearby
        parallel passage (max_dist = 4 * step easily exceeds the gap
        between hairpin legs), silently re-walking a traced leg /
        starting a second lap of a closed circuit.  The seed region
        (first 3) stays exempt so a lap closure can approach the
        start; the last 10 stay exempt for dense halved-step clusters
        at bends.
        """
        own = points[3:-10]
        for grp in (np.asarray(own, dtype=float) if own else None,
                    guard_extra):
            if grp is not None and grp.size:
                if float(np.min(np.linalg.norm(
                        grp - c_cand, axis=1))) < 0.7 * step_size:
                    return True
        return False

    points, tangents, profiles, polygons, areas = [], [], [], [], []
    closed = False
    turns = 0

    # Initial probe at the seed point
    c, area, sub = _try_at(p, t, max_dist=10 * step_size)
    if c is None:
        raise RuntimeError("seed plane has no usable cross-section")
    bbox = _piece_bbox_uv(sub, c, t)
    prof = classify_profile(area, bbox)
    points.append(c); tangents.append(t.copy())
    profiles.append(prof); polygons.append(np.array(bbox)); areas.append(area)

    for k in range(1, max_stations):
        # Adaptive step: try full step, on miss halve up to 6 times.
        # A degenerate sliver cut (< 0.3x the recent median area --
        # e.g. the disk grazing the outer edge of a sharp bend) is
        # treated as a miss, not a station.
        success = False
        turned = False
        s = step_size
        med_recent = float(np.median(np.asarray(areas[-10:], dtype=float)))
        for retry in range(7):
            p_try = points[-1] + s * t
            c, area, sub = _try_at(p_try, t, max_dist=4 * s)
            if (c is not None and area is not None
                    and (med_recent <= 0 or area >= 0.3 * med_recent)
                    and not _revisits(c)):
                success = True
                break
            s *= 0.5
        if not success and k >= 2 and turns < max_turns:
            # Sharp-corner turn (v4.95.x): the halving cannot bend the
            # Frenet tangent around a corner sharper than the local
            # curvature it can follow -- probe a direction fan instead.
            wh_max = max(max(pf.bounding_wh()) for pf in profiles[-5:])
            c, area, sub, d = _turn_search(
                solid, points[-1], t, points, areas, wh_max,
                step_size=step_size, slab_radius=slab_radius,
                slab_thickness=slab_thickness,
                extra_guard_pts=extra_guard_pts)
            if c is not None:
                success = True
                turned = True
                turns += 1
                t = d
                if verbose:
                    print(f"  step {k}: corner turn #{turns}, "
                          f"new tangent {d}")
        if not success:
            if verbose:
                print(f"  step {k}: walk halted (no valid section)")
            break

        bbox = _piece_bbox_uv(sub, c, t)
        prof = classify_profile(area, bbox)
        points.append(c); tangents.append(t.copy())
        profiles.append(prof)
        polygons.append(np.array(bbox)); areas.append(area)

        # Update tangent from last 2 centroids (Frenet first-derivative).
        # Skip right after a corner turn: the centroid difference mixes
        # the incoming and outgoing legs (corner diagonal); the
        # min-area direction from the turn search is the better
        # outgoing-tangent estimate.
        if not turned:
            t_new = points[-1] - points[-2]
            nrm = np.linalg.norm(t_new)
            if nrm > 1e-12:
                t = t_new / nrm

        # Closure check (after a few steps to avoid trivial closure).
        # Position AND direction must both match: a parallel return leg
        # passing within close_tol of the seed (e.g. the two legs of a
        # hairpin 12 mm apart) must NOT close the loop -- at a true
        # closure the walk re-enters the start ALONG the start tangent.
        if (k > 5 and np.linalg.norm(c - points[0]) < close_tol
                and float(np.dot(t, tangents[0])) > 0.5):
            if verbose:
                print(f"  step {k}: closed loop")
            closed = True
            break

        if verbose and k % 10 == 0:
            print(f"  step {k:3d}: p={c}, area={area:.4e}, used s={s:.4e}")

    return points, tangents, profiles, polygons, areas, closed, turns


def extract_centerline(step_path_or_solid, *,
                       start_hint=None,
                       step_size=None,
                       max_stations=2000,
                       close_tol=None,
                       verbose=False):
    """Walk a coil solid and return its centerline + per-station profiles.

    The march turns sharp spine corners (see ``_turn_search``) and is
    BIDIRECTIONAL: when the forward march does not close on itself,
    the other direction is marched from the same seed and the two
    open chains are stitched.  This makes a mid-wire seed (the
    axis-agnostic seed sections through the bbox centre) cover BOTH
    sides of an open coil instead of only the forward side.

    Args:
        step_path_or_solid: STEP filepath, or a netgen.occ Solid.
        start_hint: optional ((px,py,pz), (tx,ty,tz)) seed.
        step_size: walking step [m]. Defaults to the cube root of
            (solid.mass / 100), giving ~100 stations for a typical coil.
        max_stations: hard cap (per direction) to prevent runaway loops.
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

    bb = solid.bounding_box
    diag = math.sqrt(sum((bb[1][i] - bb[0][i]) ** 2 for i in range(3)))
    slab_radius = max(diag, 1e-3)
    slab_thickness = max(step_size / 50.0, 1e-7)

    march_kw = dict(step_size=step_size, max_stations=max_stations,
                    close_tol=close_tol, slab_radius=slab_radius,
                    slab_thickness=slab_thickness, verbose=verbose)
    points, tangents, profiles, polygons, areas, closed, turns = _march(
        solid, p, t, **march_kw)

    if not closed:
        # Bidirectional completion: march the opposite direction from
        # the same seed.  Backward station 0 re-probes the seed section
        # (identical slab, opposite normal sign) -- dropped from the
        # stitch.  If the backward march CLOSES (the forward march was
        # blocked immediately, e.g. seed tangent pointing off a corner
        # the turn search could not resolve from that side), the
        # backward loop IS the complete trace.
        # The forward stations (minus the shared seed region) guard the
        # backward march's steps and turns: at a junction the backward
        # march must not walk onto a leg the forward march already
        # traced (which would double-cover it in the stitched chain).
        (b_points, b_tangents, b_profiles, b_polygons, b_areas,
         b_closed, b_turns) = _march(solid, p, -t,
                                     extra_guard_pts=points[3:],
                                     **march_kw)
        turns += b_turns
        if b_closed:
            points, tangents = b_points, b_tangents
            profiles, polygons = b_profiles, b_polygons
            areas = b_areas
            closed = True
        elif len(b_points) > 1:
            # Backward stations traverse along -t; flip their tangents
            # so the stitched polyline is consistently forward-oriented.
            points = list(reversed(b_points[1:])) + points
            tangents = ([-bt for bt in reversed(b_tangents[1:])]
                        + tangents)
            profiles = list(reversed(b_profiles[1:])) + profiles
            polygons = list(reversed(b_polygons[1:])) + polygons
            areas = list(reversed(b_areas[1:])) + areas

    if not points:
        raise RuntimeError("no stations extracted")

    pts = np.array(points)
    tans = np.array(tangents)
    arclen = np.zeros(len(pts))
    for i in range(1, len(pts)):
        arclen[i] = arclen[i - 1] + np.linalg.norm(pts[i] - pts[i - 1])

    # Self-overlap net (fail-loud): a walk that TURNED at least once
    # may still have multi-lapped or leg-snapped despite the guards
    # (junction geometries are adversarial).  A silently self-
    # overlapping polyline would double-count the circuit in PEEC.
    # A pair counts as an overlap when the stations are CLOSE IN SPACE
    # but FAR ALONG THE WALK (arclength separation > 2 steps) -- dense
    # halved-step clusters at a corner are near in arclength and do
    # not count.  Only checked when turns fired: a smooth tightly-
    # wound multi-turn helix legitimately has spatially-close stations
    # a full turn apart and never turns.
    n_sta = len(pts)
    if turns > 0 and n_sta > 12:
        total = float(arclen[-1])
        if closed:
            total += float(np.linalg.norm(pts[0] - pts[-1]))
        overlap = 0
        for i in range(n_sta):
            d = np.linalg.norm(pts - pts[i], axis=1)
            ds = np.abs(arclen - arclen[i])
            if closed and total > 0:
                ds = np.minimum(ds, total - ds)
            if bool(np.any((ds > 2.0 * step_size)
                           & (d < 0.7 * step_size))):
                overlap += 1
        if overlap > max(4, 0.1 * n_sta):
            raise RuntimeError(
                f"extract_centerline: self-overlapping walk -- "
                f"{overlap}/{n_sta} stations lie within "
                f"{0.7 * step_size:.3g} of a station more than 2 steps "
                f"away along the walk, after {turns} corner turn(s).  "
                "The walk likely multi-lapped a closed circuit or "
                "snapped across parallel legs at a junction.  If this "
                "geometry is a tightly-wound coil with sharp corners, "
                "extract via the n_peri UV path or --coil-solver bem-a "
                "instead.")

    return CenterlineResult(
        polyline=pts,
        tangents=tans,
        profiles=profiles,
        polygons=polygons,
        arclen=arclen,
        closed=closed,
        areas=np.asarray(areas, dtype=float),
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
