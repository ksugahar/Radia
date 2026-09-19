"""Unified coil-topology classifier for STEP-driven PEEC.

Single source of truth for whether a CAD solid is a CLOSED-loop coil
(full revolution, no source/sink) or an OPEN coil (with two source /
sink cap faces and a small angular gap), plus the rotation axis and
spine arc parameters that all five filament generation Paths in
``coil_from_cad.py`` should consume.

Created 2026-05-02 to replace the scattered, inconsistent OPEN-vs-CLOSED
handling in ``_spine_from_rotation_axis_z`` (full-360deg fallback that
clipped 14deg of conductor on the rect_torus_lofted_united fixture --
visually confirmed in GMSH against the BEM-A surface mesh).

Note: not to be confused with ``coil_geometry.py`` (which builds OCC
shapes for Radia coil objects, a separate concern).

API
===
- ``CoilTopology`` dataclass (``is_open``, ``cap_a/b``, ``theta_a/b``,
  ``sweep_deg``, ``axis``, ``R_spine``).
- ``detect_cap_faces(solid)``: returns ``(cap_a, cap_b)`` if 2 small
  planar end-cap faces exist, else ``None``.
- ``classify_axis(solid, hint)``: returns the rotation-axis unit vector
  (``z`` hard-coded for now; ``hint`` reserved for future PCA-based
  detection).
- ``extract_coil_topology(solid)``: build a fully populated
  ``CoilTopology``.
- ``generate_spine(topo, n_segments)``: ``(n_segments, 3)`` spine
  polyline.  OPEN samples ``[theta_a, theta_b]`` along the arc the
  extractor measured as occupied with ``endpoint=True``; CLOSED samples
  ``[0, 2pi)`` with ``endpoint=False``.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional, Tuple, Any

import numpy as np


_AXIS_Z = np.array([0.0, 0.0, 1.0], dtype=np.float64)


@dataclass
class CoilTopology:
    """Topology + spine parameters of a coil solid.

    Fields:
        is_open: True if 2 cap faces detected (OPEN coil with terminals),
            False if no caps (CLOSED full revolution).
        cap_a, cap_b: build123d Face objects for the two end caps when
            ``is_open``, else ``None``.
        theta_a, theta_b: cap face centroid angles around the rotation
            axis, in **radians** (``atan2``).  Only meaningful for OPEN.
        sweep_deg: arc length spanned by the conductor in degrees,
            MEASURED by probing the material (``360`` for CLOSED).  A
            nearly-closed coil gives ``360 - gap_deg``; a partial-arc
            conductor (e.g. a 120-deg beak-fin segment) gives its own
            short arc, not the complementary one.
        sweep_ccw: True when the occupied arc runs counter-clockwise
            from ``theta_a`` to ``theta_b`` about ``axis``.  Only
            meaningful for OPEN.
        axis: (3,) unit vector for the rotation axis.
        R_spine: scalar spine radius estimate.  Bbox-derived; the
            actual spine radius is recovered by the caller via
            cross-section centroids when sectioning succeeds.
        cross_section_kind: optional hint ('circle' / 'rect' /
            'polygon' / 'bspline' / 'unknown').  Currently unused by
            the spine generator; reserved for future profile selection.
    """
    is_open: bool
    cap_a: Optional[Any] = None
    cap_b: Optional[Any] = None
    theta_a: float = 0.0
    theta_b: float = 0.0
    sweep_deg: float = 360.0
    sweep_ccw: bool = True
    axis: np.ndarray = field(default_factory=lambda: _AXIS_Z.copy())
    R_spine: float = 0.0
    cross_section_kind: str = "unknown"


# ----------------------------------------------------------------------
# Cap detection
# ----------------------------------------------------------------------
def _material_depth_reaches(solid, face, depth: float) -> bool:
    """True if the solid's material extends at least ``depth`` along
    the INWARD normal from the face centroid.

    Centroid and normal are computed directly through OCP
    (``BRepGProp.SurfaceProperties`` + the surface normal at the UV
    midpoint, orientation-corrected) rather than through the
    duck-typed Face wrappers: the real build123d ``Face.center()``
    returned origin / garbage coordinates on an OCCT-swept solid
    (measured 2026-07-29 on a 355-deg JernArc rect sweep) while the
    GProp centroid is exact for both wrapper flavours.

    The inward sign is determined by probing +-0.3 * sqrt(area): the
    side that is inside the solid is inward (guards against OCC face
    orientation quirks).  Returns False when neither or both probes
    are inside (degenerate).
    """
    from OCP.BRep import BRep_Tool
    from OCP.BRepClass3d import BRepClass3d_SolidClassifier
    from OCP.BRepGProp import BRepGProp
    from OCP.BRepTools import BRepTools
    from OCP.GProp import GProp_GProps
    from OCP.GeomLProp import GeomLProp_SLProps
    from OCP.TopAbs import TopAbs_IN, TopAbs_REVERSED
    from OCP.gp import gp_Pnt

    fw = face.wrapped
    props = GProp_GProps()
    BRepGProp.SurfaceProperties_s(fw, props)
    com = props.CentreOfMass()
    c0 = np.array([com.X(), com.Y(), com.Z()], dtype=float)

    surf = BRep_Tool.Surface_s(fw)
    umin, umax, vmin, vmax = BRepTools.UVBounds_s(fw)
    sl = GeomLProp_SLProps(surf, 0.5 * (umin + umax),
                           0.5 * (vmin + vmax), 1, 1e-9)
    if not sl.IsNormalDefined():
        return False
    n = sl.Normal()
    nn = np.array([n.X(), n.Y(), n.Z()], dtype=float)
    ln = float(np.linalg.norm(nn))
    if ln < 1e-12:
        return False
    nn /= ln
    if fw.Orientation() == TopAbs_REVERSED:
        nn = -nn
    probe0 = 0.3 * math.sqrt(max(float(face.area), 1e-30))

    def _inside(p):
        cls = BRepClass3d_SolidClassifier(solid.wrapped)
        cls.Perform(gp_Pnt(float(p[0]), float(p[1]), float(p[2])), 1e-9)
        return cls.State() == TopAbs_IN

    plus_in = _inside(c0 + probe0 * nn)
    minus_in = _inside(c0 - probe0 * nn)
    if plus_in == minus_in:
        return False
    inward = -nn if minus_in else nn
    # Monotone depth: material at half depth AND at full depth.
    return (_inside(c0 + 0.5 * depth * inward)
            and _inside(c0 + depth * inward))


def _face_centroid(face) -> np.ndarray:
    """Exact face centroid (3,) through ``BRepGProp``.

    ``Face.center()`` returned origin / garbage coordinates on an
    OCCT-swept solid (measured 2026-07-29, see
    :func:`_material_depth_reaches`), so the cap angles that decide the
    spine direction are taken from the GProp centroid instead.
    """
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps

    props = GProp_GProps()
    BRepGProp.SurfaceProperties_s(face.wrapped, props)
    com = props.CentreOfMass()
    return np.array([com.X(), com.Y(), com.Z()], dtype=float)


def _lateral_azimuths(solid, caps, n_uv: int = 5) -> np.ndarray:
    """Azimuths (rad, about +z) of UV-grid samples on every non-cap face.

    The conductor's lateral surface covers exactly the arc the conductor
    occupies, so the azimuths of points on it say which of the two arcs
    joining the caps is the real one -- without needing a spine radius.
    Radius-free matters: on a lead-bearing coil both cap centroids sit
    essentially on the rotation axis (measured 2026-09-20 on
    ``keiko_outsideline.step``: cap centroid radius 0.05 mm), so any
    probe riding a cross-section centroid circle misses the conductor.

    Samples are taken on the surface's UV box, so a trimmed face
    contributes a few points outside its trim.  That only blurs the
    azimuth histogram slightly; the caller decides on a clear majority.
    """
    from OCP.BRep import BRep_Tool
    from OCP.BRepTools import BRepTools

    cap_shapes = [cap.wrapped for cap in caps]
    azimuths = []
    for face in solid.faces():
        wrapped = face.wrapped
        if any(wrapped.IsSame(cap) for cap in cap_shapes):
            continue
        surface = BRep_Tool.Surface_s(wrapped)
        umin, umax, vmin, vmax = BRepTools.UVBounds_s(wrapped)
        for iu in range(n_uv):
            u = umin + (umax - umin) * (iu + 0.5) / n_uv
            for iv in range(n_uv):
                v = vmin + (vmax - vmin) * (iv + 0.5) / n_uv
                point = surface.Value(u, v)
                azimuths.append(math.atan2(point.Y(), point.X()))
    return np.asarray(azimuths, dtype=float)


def detect_cap_faces(solid, depth_factor: float = 2.0
                       ) -> Optional[Tuple[Any, Any]]:
    """Return ``(cap_a, cap_b)`` if the solid has 2 source/sink cap faces.

    A true end cap has a decisive GEOMETRIC signature: the wire runs
    AWAY from it, so the solid's material extends several
    cross-section sizes along the cap's inward normal.  Every other
    planar face is shallow in its own normal direction: a lateral
    panel exits after one wire thickness, a fused boss face after the
    boss height.  Classification:

      1. Every PLANE face is depth-tested: material must reach
         ``depth_factor * sqrt(area)`` along the inward normal
         (point-in-solid probes at 0.5x and 1.0x that depth).
      2. EXACTLY 2 survivors -> those are the caps.  The two end
         cross-sections need NOT match: a split or tapered conductor
         has different end areas (measured 2026-07-29 on
         G_equator_split: verified caps of 1.395 and 40.867 mm^2,
         29x apart -- an area-pair requirement here mis-routed it
         CLOSED).  A wire has exactly two places where the material
         runs away from a planar face, so no disambiguation is
         needed.
      3. MORE than 2 survivors -> area disambiguation: the caps must
         be a UNIQUE tight area pair (within 10%) and every other
         survivor must differ from it by more than 30% -- e.g. a deep
         fused-boss top face is allowed as long as its area is
         clearly different.  Two tight pairs (measured on
         B_rect_sweep: 63.6/63.9 mm^2 AND 66.52/66.54 mm^2, a
         butt-jointed multi-segment conductor with four exposed
         wire-end faces) are ambiguous -> None (refuse; the dispatch
         treats the coil as CLOSED and the downstream
         volume-reconciliation guard reports it).

    This replaces the 2026-05-02 area-gap heuristic ("the 2 smallest
    planes, accepted only when the 3rd is > 2x larger"), which could
    not see caps on rectangular-wire sweeps whose lateral panel areas
    are comparable to the cap area (measured 2026-07-29:
    B_rect_sweep caps 63.6/63.9 mm^2 vs a 65.3 mm^2 lateral facet;
    J_boss_fused caps 77.1/77.1 mm^2 vs 133-144 mm^2 boss faces --
    both were refused and mis-routed CLOSED).

    Re-verified fixture behaviour:
      - closed torus (0 PLANE) -> None (CLOSED)
      - gapped torus (2 PLANE + revolution lateral) -> 2 caps
      - rect_torus_lofted_united (28 PLANE) -> the 2 true caps (the
        26 lateral panels are one wire thickness shallow)
      - 3turnCoil (2 PLANE + 617 BSPLINE) -> 2 caps
    """
    from radia._b3d_shim import GeomType
    plane_faces = [f for f in solid.faces() if f.geom_type == GeomType.PLANE]
    if len(plane_faces) < 2:
        return None

    survivors = []
    for f in plane_faces:
        area = float(f.area)
        if area <= 0:
            continue
        if _material_depth_reaches(solid, f,
                                   depth_factor * math.sqrt(area)):
            survivors.append(f)
        if len(survivors) > 8:
            return None  # degenerate: many deep planar faces

    if len(survivors) < 2:
        return None
    survivors.sort(key=lambda f: float(f.area))
    if len(survivors) == 2:
        # Depth signature alone is decisive (docstring point 2); the
        # end areas may legitimately differ (G_equator_split: 29x).
        return (survivors[0], survivors[1])
    areas = [float(f.area) for f in survivors]

    # Tight pairs (within 10%).
    pairs = [(i, j)
             for i in range(len(survivors))
             for j in range(i + 1, len(survivors))
             if areas[j] <= 1.10 * areas[i]]
    if len(pairs) != 1:
        return None  # no pair, or ambiguous multiple pairs
    i, j = pairs[0]
    # Every OTHER survivor must be clearly different from the pair.
    pair_mean = 0.5 * (areas[i] + areas[j])
    for k in range(len(survivors)):
        if k in (i, j):
            continue
        if 0.77 * pair_mean <= areas[k] <= 1.30 * pair_mean:
            return None
    return (survivors[i], survivors[j])


# ----------------------------------------------------------------------
# Rotation axis
# ----------------------------------------------------------------------
def classify_axis(solid, hint: str = "z") -> np.ndarray:
    """Return the rotation-axis unit vector for a coil solid.

    Currently:
      - ``hint == "z"`` (default): return +z.  All current Radia
        fixtures (gapped torus, rect_torus_lofted, 3turnCoil) use the
        z-axis as their rotation axis by convention.
      - Other hints: not yet implemented; raise.

    PCA-based axis detection (e.g. principal direction of a
    rotation-symmetric solid) is reserved for a future iteration.
    """
    if hint == "z":
        return _AXIS_Z.copy()
    raise NotImplementedError(
        f"classify_axis: hint={hint!r} not implemented yet.")


# ----------------------------------------------------------------------
# Spine radius (bbox-based estimate)
# ----------------------------------------------------------------------
def _bbox_spine_radius(solid) -> float:
    """Rough spine radius from solid bbox.

    For a coil swept around the z-axis, the spine sits inside the
    bbox's outer diameter; we use 0.85 * R_outer as a typical
    fraction (cross-section width is usually ~10-20 % of R).
    """
    bbox = solid.bounding_box()
    R_outer = max(abs(bbox.max.X), abs(bbox.min.X),
                   abs(bbox.max.Y), abs(bbox.min.Y))
    return 0.85 * R_outer


# ----------------------------------------------------------------------
# Top-level extractor
# ----------------------------------------------------------------------
def extract_coil_topology(solid, axis_hint: str = "z") -> CoilTopology:
    """Single-source-of-truth classifier for a coil solid.

    Determines OPEN vs CLOSED, cap faces + their angles around the
    rotation axis, the LONG arc sweep, and a bbox-based spine radius.
    The returned ``CoilTopology`` is consumed by all 5 filament
    generation Paths in ``coil_from_cad.py`` so they share a coherent
    view of the geometry.

    Behaviour summary:
      - 2 cap faces detected -> OPEN.  ``theta_a, theta_b`` are the
        cap centroid angles; ``sweep_deg`` / ``sweep_ccw`` describe
        whichever of the two arcs between them was measured to carry
        material.
      - 0 cap faces -> CLOSED.  ``sweep_deg = 360``.
    """
    if solid is None:
        raise ValueError("extract_coil_topology: solid is None")

    R_spine = _bbox_spine_radius(solid)
    if R_spine <= 0.0:
        raise ValueError(
            f"extract_coil_topology: degenerate bbox (R_spine={R_spine})")

    axis = classify_axis(solid, axis_hint)
    caps = detect_cap_faces(solid)

    if caps is None:
        return CoilTopology(
            is_open=False,
            sweep_deg=360.0,
            axis=axis,
            R_spine=R_spine,
        )

    cap_a, cap_b = caps
    ca, cb = _face_centroid(cap_a), _face_centroid(cap_b)
    theta_a = math.atan2(ca[1], ca[0])
    theta_b = math.atan2(cb[1], cb[0])
    delta_ccw = (theta_b - theta_a) % (2.0 * math.pi)

    # Two arcs join the caps: CCW by delta_ccw, or CW by the
    # complement.  MEASURE which one carries the conductor instead of
    # assuming the longer one.  The long-arc assumption is right for a
    # nearly-closed 1-turn coil but sends the spine around the outside
    # of a partial-arc conductor (measured 2026-09-20 on a 120-deg
    # revolved beak fin: the spine took the empty 240-deg arc, so the
    # section planes missed the solid and the downstream spacing guard
    # blamed the CAD for "straight leads").
    azimuths = _lateral_azimuths(solid, caps)
    if azimuths.size == 0:
        raise ValueError(
            "extract_coil_topology: the solid has no lateral face "
            "outside the two caps, so the occupied arc cannot be "
            "measured.")
    relative = (azimuths - theta_a) % (2.0 * math.pi)
    frac_ccw = float(np.count_nonzero(relative < delta_ccw)) / relative.size
    if 0.25 < frac_ccw < 0.75:
        raise ValueError(
            f"extract_coil_topology: the lateral surface straddles both "
            f"arcs between the caps ({frac_ccw:.0%} of "
            f"{relative.size} samples on the CCW side), so the solid is "
            f"not a single arc sweep about {axis.tolist()} "
            f"(multi-turn, self-overlapping, or a different axis).")

    sweep_ccw = frac_ccw >= 0.75
    sweep_deg = math.degrees(
        delta_ccw if sweep_ccw else 2.0 * math.pi - delta_ccw)

    return CoilTopology(
        is_open=True,
        cap_a=cap_a, cap_b=cap_b,
        theta_a=theta_a, theta_b=theta_b,
        sweep_deg=sweep_deg,
        sweep_ccw=sweep_ccw,
        axis=axis,
        R_spine=R_spine,
    )


# ----------------------------------------------------------------------
# Spine generation
# ----------------------------------------------------------------------
def generate_spine(topo: CoilTopology, n_segments: int) -> np.ndarray:
    """Spine polyline along the conductor centerline.

    Args:
        topo: ``CoilTopology`` from ``extract_coil_topology``.
        n_segments: number of stations.  For OPEN, the first station
            sits on cap A and the last station on cap B; for CLOSED,
            ``endpoint=False`` so the last station does NOT wrap back
            to the first.

    Returns:
        ``(n_segments, 3)`` float64 array of (x, y, z).  Currently
        always in the z=0 plane (rotation axis = z assumption).
    """
    if n_segments < 2:
        raise ValueError(
            f"generate_spine: need n_segments >= 2 (got {n_segments})")
    if not np.allclose(topo.axis, _AXIS_Z):
        raise NotImplementedError(
            "generate_spine: non-z rotation axis not yet supported")

    if topo.is_open:
        # Walk the arc the extractor MEASURED as occupied, in its
        # measured direction; theta_a + delta lands on theta_b mod 2pi.
        delta = math.radians(topo.sweep_deg)
        thetas = np.linspace(topo.theta_a,
                              topo.theta_a + (delta if topo.sweep_ccw
                                              else -delta),
                              n_segments)
    else:
        # CLOSED: full revolution, do NOT include the wrap-around point.
        thetas = np.linspace(0.0, 2.0 * math.pi, n_segments,
                              endpoint=False)

    return np.column_stack([
        topo.R_spine * np.cos(thetas),
        topo.R_spine * np.sin(thetas),
        np.zeros_like(thetas),
    ])
