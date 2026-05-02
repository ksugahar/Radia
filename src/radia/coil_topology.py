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
  polyline.  OPEN samples ``[theta_a, theta_b]`` along the LONG arc
  with ``endpoint=True``; CLOSED samples ``[0, 2pi)`` with
  ``endpoint=False``.
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
        sweep_deg: LONG arc length spanned by the conductor in degrees
            (``360 - gap_deg`` for OPEN, ``360`` for CLOSED).
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
    axis: np.ndarray = field(default_factory=lambda: _AXIS_Z.copy())
    R_spine: float = 0.0
    cross_section_kind: str = "unknown"


# ----------------------------------------------------------------------
# Cap detection
# ----------------------------------------------------------------------
def detect_cap_faces(solid, area_ratio_threshold: float = 2.0
                       ) -> Optional[Tuple[Any, Any]]:
    """Return ``(cap_a, cap_b)`` if the solid has 2 source/sink cap faces.

    Heuristic:
      - Caps are the 2 SMALLEST PLANE faces.
      - Caps are accepted iff there are EXACTLY 2 PLANE faces, OR the
        3rd-smallest PLANE face has area > ``area_ratio_threshold`` x
        the 2nd-smallest cap (= the 2 caps stand out as distinctly
        smaller than the rest of the planar lateral panels).
      - 0 PLANE faces -> CLOSED (returns None).

    Verified 2026-05-02 on:
      - ``ih_closed_torus_coil.step``: 4 TORUS / 0 PLANE -> None (CLOSED)
      - ``ih_peec_inductance_coil.step``: 2 PLANE + 4 TORUS -> 2 caps
      - ``rect_torus_lofted_united.step``: 28 PLANE; 2 smallest are
        the 8x6mm caps (~4.8e-5 m^2) vs the 26 lateral panels at
        ~1.4e-4 m^2 (~3x larger) -> 2 caps detected
      - ``3turnCoil_work_coil_m.step``: 2 PLANE + 617 BSPLINE -> 2 caps
    """
    from build123d import GeomType
    plane_faces = [f for f in solid.faces() if f.geom_type == GeomType.PLANE]
    if len(plane_faces) < 2:
        return None
    plane_faces.sort(key=lambda f: f.area)
    cap_a, cap_b = plane_faces[0], plane_faces[1]
    if len(plane_faces) >= 3:
        if plane_faces[2].area <= area_ratio_threshold * cap_b.area:
            # 3rd-smallest plane is comparable to the cap candidates --
            # heuristic doesn't apply (e.g. polygon-cross-section coil
            # with many small lateral facets).  Refuse to classify.
            return None
    return (cap_a, cap_b)


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
        cap centroid angles; ``sweep_deg`` is the LONG arc between
        them (i.e. ``360 - gap_deg``).
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
    ca, cb = cap_a.center(), cap_b.center()
    theta_a = math.atan2(ca.Y, ca.X)
    theta_b = math.atan2(cb.Y, cb.X)
    # LONG arc length: pick whichever direction (CCW from theta_a to
    # theta_b, or the other way) is longer.
    delta_ccw = (theta_b - theta_a) % (2.0 * math.pi)
    if delta_ccw < math.pi:
        # Short arc is CCW; long arc is CW.
        sweep_deg = math.degrees(2.0 * math.pi - delta_ccw)
    else:
        sweep_deg = math.degrees(delta_ccw)

    return CoilTopology(
        is_open=True,
        cap_a=cap_a, cap_b=cap_b,
        theta_a=theta_a, theta_b=theta_b,
        sweep_deg=sweep_deg,
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
        delta_ccw = (topo.theta_b - topo.theta_a) % (2.0 * math.pi)
        if delta_ccw > math.pi:
            # LONG arc is CCW from theta_a to theta_b.
            thetas = np.linspace(topo.theta_a,
                                  topo.theta_a + delta_ccw,
                                  n_segments)
        else:
            # LONG arc is CW from theta_a (i.e. through negative dtheta).
            thetas = np.linspace(topo.theta_a,
                                  topo.theta_a - (2.0 * math.pi - delta_ccw),
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
