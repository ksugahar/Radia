"""coil_from_jou.py — Build PEEC filaments from a Cubit .jou explicit centerline.

Complements coil_from_step.py for coils where `extract_centerline` cannot
walk the full path (multi-turn helices, 3turnCoil-class).  The .jou file
already encodes the centerline as explicit `move Surface N x Y y Y z Z`
coordinates, so we can skip walker heuristics entirely.

Canonical input shape (Cubit jou pattern for a lofted coil):

  create surface circle radius 3.15 zplane           # cross-section
  rotate Surface 1 angle ... direction ... include_merged
  move Surface 1 x 129.394895 y -12.500000 z 10.000000 include_merged
  create surface circle radius 3.15 zplane
  ...
  move Surface 383 x 129.428679 y 12.500000 z 10.000000 include_merged
  create volume loft surface 1 2
  create volume loft surface 2 3
  ...

Each `move Surface` anchors one cross-section position → N move commands
define the N-point centerline.  The `create surface circle radius R`
before the first move gives the cross-section radius (assumed constant).

Pipeline:

  parse_jou_centerline(jou)
    -> (path_points_m, radius_m)
    -> _parallel_transport_frame(path_points)
    -> for each of n_peri perimeter angles, offset the centerline
    -> feed n_peri filament polylines to peec_bundle.build_bundle_solver
    -> topology_dict (same shape as filaments_from_step output)

Limitations (current scope):

  - Cross-section: circular only (radius from `create surface circle`).
    Rect / lofted profiles need separate parsers.
  - Units: meters by default (cad_units_per_meter=1.0, CLAUDE.md
    "Unit System Policy: Radia always uses meters").  Override via
    arg if your .jou is in mm or cm.
  - Open path only (no closure).  PEEC solver adds port_plus/port_minus
    at first/last centerline point.

When to use this vs coil_from_step:

  - `filaments_from_step`: clean single-loop coils (torus, simple bend).
    extract_centerline walks the solid boundary.
  - `filaments_from_jou` (this file): multi-turn helices, 3turnCoil-class
    where the walker halts mid-coil.  Requires the original .jou (which
    has the explicit centerline).
"""
from __future__ import annotations

import math
import re
from typing import Optional, Tuple

import numpy as np


RE_CIRCLE = re.compile(r'^\s*create\s+surface\s+circle\s+radius\s+([\d.eE+-]+)')
RE_MOVE = re.compile(
    r'^\s*move\s+Surface\s+\d+\s+x\s+([-\d.eE+]+)'
    r'\s+y\s+([-\d.eE+]+)\s+z\s+([-\d.eE+]+)'
)


def parse_jou_centerline(jou_path: str,
                          cad_units_per_meter: float = 1.0
                          ) -> Tuple[np.ndarray, float]:
    """Return (points_m, radius_m) from a Cubit .jou file.

    Scans for the first `create surface circle radius R` (cross-section)
    and every `move Surface N x Y y Y z Z` (centerline anchor).

    Raises ValueError if the file doesn't match the expected pattern.
    """
    points = []
    radius_raw: Optional[float] = None
    with open(jou_path, encoding='utf-8') as f:
        for line in f:
            m = RE_CIRCLE.match(line)
            if m and radius_raw is None:
                radius_raw = float(m.group(1))
                continue
            m = RE_MOVE.match(line)
            if m:
                points.append([float(m.group(i)) for i in (1, 2, 3)])
    if not points:
        raise ValueError(
            f"no `move Surface N x ... y ... z ...` commands in {jou_path}")
    if radius_raw is None:
        raise ValueError(
            f"no `create surface circle radius R` in {jou_path}")
    scale = 1.0 / cad_units_per_meter
    return np.asarray(points) * scale, radius_raw * scale


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

    # Per-segment unit tangents, then per-vertex (average of adjacents)
    seg_t = np.diff(pts, axis=0)
    seg_len = np.linalg.norm(seg_t, axis=1, keepdims=True)
    # Guard against zero-length segments
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
    # Seed u_0: world axis least parallel to t_0, Gram-Schmidt projected
    t0 = tangent[0]
    pick = int(np.argmin(np.abs(t0)))
    cand = np.zeros(3)
    cand[pick] = 1.0
    u0 = cand - np.dot(cand, t0) * t0
    u0_norm = np.linalg.norm(u0)
    if u0_norm < 1e-12:
        # Very unlikely but guard: fall back to any orthogonal
        cand = np.array([0.0, 1.0, 0.0]) if pick == 0 else np.array([1.0, 0.0, 0.0])
        u0 = cand - np.dot(cand, t0) * t0
        u0_norm = np.linalg.norm(u0)
    u_hat[0] = u0 / u0_norm
    v_hat[0] = np.cross(tangent[0], u_hat[0])

    # Rodrigues rotation forward along the polyline
    for i in range(1, n):
        t_prev = tangent[i - 1]
        t_curr = tangent[i]
        axis = np.cross(t_prev, t_curr)
        sin_a = np.linalg.norm(axis)
        cos_a = float(np.clip(np.dot(t_prev, t_curr), -1.0, 1.0))
        if sin_a < 1e-12:
            # Tangent didn't turn (straight through) → u unchanged
            u_hat[i] = u_hat[i - 1]
        else:
            k = axis / sin_a
            u_prev = u_hat[i - 1]
            # Rodrigues:  u' = u cos + (k × u) sin + k (k·u)(1 - cos)
            u_new = (u_prev * cos_a
                     + np.cross(k, u_prev) * sin_a
                     + k * np.dot(k, u_prev) * (1.0 - cos_a))
            u_hat[i] = u_new
        # Re-orthogonalise against the current tangent (numerical drift)
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

    Core engine shared by ``filaments_from_jou`` (Cubit journal input)
    and the longest-edge STEP fallback in ``coil_from_cad``.

    Args:
        pts_m: (N, 3) float64 centerline points in meters.
        radius_m: circular cross-section radius in meters.
        sigma, n_peri: solver parameters (thin-skin perimeter).
        source_tag: free-form string copied into the result dict's
            ``source`` field (e.g. "jou", "step_longest_edge").

    Returns:
        Same topology_dict shape as ``filaments_from_step``, so
        downstream PEEC (``solve_loop_bundle``, ``PEECCircuitSolver``)
        works unchanged.
    """
    pts = np.asarray(pts_m, dtype=float)
    n_pts = len(pts)
    if n_pts < 2:
        raise ValueError(f"need at least 2 centerline points, got {n_pts}")

    _, u_hat, v_hat = _parallel_transport_frame(pts)

    # CircleProfile perimeter samples (cell-centered, arc-length-spaced)
    theta = 2.0 * np.pi * (np.arange(n_peri) + 0.5) / float(n_peri)
    u_offset = radius_m * np.cos(theta)  # (n_peri,)
    v_offset = radius_m * np.sin(theta)

    # Per-filament cell area: total cross-section / n_peri (isotropic side)
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

    # Deferred import: triggers radia's MKL DLL setup before peec_bundle.
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


def filaments_from_jou(jou_path: str, *,
                       sigma: float = 5.8e7,
                       n_peri: int = 16,
                       cad_units_per_meter: float = 1.0):
    """Build a PEEC topology_dict from a Cubit .jou file.

    Thin wrapper around ``filaments_from_polyline`` that first parses
    the centerline and cross-section radius out of the .jou.
    """
    pts, radius = parse_jou_centerline(jou_path, cad_units_per_meter)
    topo = filaments_from_polyline(pts, radius,
                                    sigma=sigma, n_peri=n_peri,
                                    source_tag="jou")
    return topo
