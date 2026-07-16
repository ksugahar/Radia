"""Unit tests: orient_surface_triangles (global winding consistency).

Background (2026-07-17, Takahashi genus-1 workpiece): the hole extractor's
old per-triangle "centroid-outward" flip actively CREATED an inconsistent
winding on a genus-1 tube (the bore-wall outward normal points TOWARD the
centroid, so the whole inner wall was flipped -- 199 directed-edge
conflicts).  The inconsistent winding corrupts the scalar BIE's
double-layer operator: Takahashi 7 kHz / mu_r=100 gave P_wp = 37.9 kW
instead of 22.5 kW (references 17.0-17.7 kW) -- the dominant share of the
known x2 heating over-estimate.  ``orient_surface_triangles`` replaces the
heuristic with face-BFS flip propagation + per-component signed-volume
outward normalisation.
"""
from __future__ import annotations

import numpy as np

import sys
from pathlib import Path

_PANELS = Path(__file__).resolve().parents[1] / "src" / "radia" / "panels"
if str(_PANELS) not in sys.path:
    sys.path.insert(0, str(_PANELS))

from surface_mesh_extract import orient_surface_triangles


def _octahedron(center=(0.0, 0.0, 0.0), scale=1.0):
    c = np.asarray(center, dtype=float)
    pts = c + scale * np.array([
        [0, 0, 1], [1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0], [0, 0, -1],
    ], dtype=float)
    eq = [1, 2, 3, 4]
    tris = []
    for k in range(4):
        a, b = eq[k], eq[(k + 1) % 4]
        tris.append([0, a, b])          # top fan (outward CCW)
        tris.append([5, b, a])          # bottom fan
    return pts, np.array(tris, dtype=np.int64)


def _signed_volume(pts, tris):
    return sum(float(np.dot(pts[t[0]], np.cross(pts[t[1]], pts[t[2]])))
               for t in tris) / 6.0


def test_consistent_outward_mesh_is_a_noop():
    pts, tris = _octahedron()
    assert _signed_volume(pts, tris) > 0
    out, stats = orient_surface_triangles(pts, tris)
    assert stats["n_flips"] == 0
    assert stats["components_flipped"] == 0
    assert stats["conflicts_after"] == 0
    assert np.array_equal(out, tris)


def test_scrambled_winding_is_repaired():
    pts, tris = _octahedron()
    bad = tris.copy()
    for ti in (1, 3, 6):                      # flip a few triangles
        bad[ti, 1], bad[ti, 2] = bad[ti, 2], bad[ti, 1]
    out, stats = orient_surface_triangles(pts, bad)
    assert stats["conflicts_before"] > 0
    assert stats["conflicts_after"] == 0
    assert _signed_volume(pts, out) > 0


def test_fully_inverted_component_is_flipped_outward():
    pts, tris = _octahedron()
    inward = tris.copy()
    inward[:, [1, 2]] = inward[:, [2, 1]]     # consistent but inward
    out, stats = orient_surface_triangles(pts, inward)
    assert stats["conflicts_after"] == 0
    assert stats["components_flipped"] == 1
    assert _signed_volume(pts, out) > 0


def test_two_components_oriented_independently():
    p1, t1 = _octahedron(center=(0, 0, 0))
    p2, t2 = _octahedron(center=(10, 0, 0))
    t2 = t2.copy()
    t2[:, [1, 2]] = t2[:, [2, 1]]             # second component inward
    pts = np.vstack([p1, p2])
    tris = np.vstack([t1, t2 + len(p1)])
    out, stats = orient_surface_triangles(pts, tris)
    assert stats["n_components"] == 2
    assert stats["components_flipped"] == 1
    assert stats["conflicts_after"] == 0
    # each component outward on its own
    v1 = _signed_volume(pts, out[:len(t1)])
    v2 = _signed_volume(pts, out[len(t1):])
    assert v1 > 0 and v2 > 0


def test_genus1_torus_grid_consistent():
    """A structured torus grid (chi=0) with randomised winding must come
    back conflict-free -- the genus-1 case is exactly where the old
    centroid heuristic broke."""
    R, b = 3.0, 1.0
    nu, nv_ = 12, 8
    pts = []
    for i in range(nu):
        tu = 2 * np.pi * i / nu
        for j in range(nv_):
            tv = 2 * np.pi * j / nv_
            rho = R + b * np.cos(tv)
            pts.append([rho * np.cos(tu), rho * np.sin(tu), b * np.sin(tv)])
    pts = np.array(pts)

    def vid(i, j):
        return (i % nu) * nv_ + (j % nv_)

    tris = []
    for i in range(nu):
        for j in range(nv_):
            a, b2 = vid(i, j), vid(i + 1, j)
            c, d = vid(i + 1, j + 1), vid(i, j + 1)
            tris.append([a, b2, c])
            tris.append([a, c, d])
    tris = np.array(tris, dtype=np.int64)
    rng = np.random.default_rng(7)
    bad = tris.copy()
    flip = rng.random(len(bad)) < 0.5
    bad[flip] = bad[flip][:, [0, 2, 1]]
    out, stats = orient_surface_triangles(pts, bad)
    assert stats["conflicts_after"] == 0
    assert _signed_volume(pts, out) > 0
