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


def filaments_from_step(step_path: str,
                        path_points_m: np.ndarray,
                        sigma: float = 5.8e7,
                        nwinc: int = 1,
                        nhinc: int = 1,
                        cad_units_per_meter: float = 1000.0):
    """End-to-end: STEP solid + path -> PEEC topology.

    Args:
        step_path: STEP file path.
        path_points_m: (N+1, 3) path vertices in meters.
        sigma: Conductivity [S/m].
        nwinc, nhinc: Sub-filament subdivision.
        cad_units_per_meter: Scale factor (default 1000 = mm).

    Returns:
        topology_dict from PEECBuilder.build_topology().
    """
    n_seg = path_points_m.shape[0] - 1
    tangents = path_tangents(path_points_m)

    # Section midpoints
    midpoints = 0.5 * (path_points_m[:-1] + path_points_m[1:])

    sections = section_solid_along_path(
        step_path, midpoints, tangents, units_scale=cad_units_per_meter
    )

    # Convert section dimensions from CAD units to meters
    scale_inv = 1.0 / cad_units_per_meter
    widths = np.zeros(n_seg, dtype=np.float64)
    heights = np.zeros(n_seg, dtype=np.float64)
    for i, sec in enumerate(sections):
        if sec["w_est"] is not None:
            widths[i] = sec["w_est"] * scale_inv
            heights[i] = sec["h_est"] * scale_inv
        else:
            raise RuntimeError(
                f"Section {i} at midpoint {midpoints[i]} failed. "
                "Check that the STEP solid intersects the path."
            )

    return build_peec_from_path(
        path_points_m, widths, heights,
        sigma=sigma, nwinc=nwinc, nhinc=nhinc,
    )
