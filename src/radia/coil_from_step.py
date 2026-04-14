"""coil_from_step.py

Sprint 1 (CAD-from-STEP): walking-plane centerline + per-station profile
extraction from an arbitrary STEP coil solid.

Pipeline:
    STEP file -> OCCGeometry -> single Solid
              -> walking-plane sweep along the conductor
              -> (centerline polyline, [Profile per station])

The result is consumed directly by ``peec_bundle.build_bundle_solver``,
or it can be promoted to a CoilBuilder reconstruction in Sprint 2.

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

from coil_profile import (Profile, RectProfile, CircleProfile,
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


# ---------- STEP loading ----------------------------------------------------


def load_step_solid(step_path):
    """Load a STEP file and return its single Solid.

    Raises ValueError if the STEP contains 0 or >1 solids.
    """
    from netgen.occ import OCCGeometry
    if not os.path.isfile(step_path):
        raise FileNotFoundError(step_path)
    geo = OCCGeometry(step_path)
    shape = geo.shape
    solids = list(shape.solids)
    if len(solids) != 1:
        raise ValueError(
            f"STEP must contain exactly 1 solid, got {len(solids)} "
            f"in {step_path}")
    return solids[0]


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
