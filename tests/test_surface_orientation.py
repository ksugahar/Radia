"""Unit tests: orient_surface_triangles (global winding consistency).

The hole extractor's old per-triangle "centroid-outward" flip actively
created an inconsistent
winding on a genus-1 tube (the bore-wall outward normal points TOWARD the
centroid, so the whole inner wall is flipped).  The inconsistent winding
corrupts the scalar BIE's double-layer operator.
``orient_surface_triangles`` replaces the
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


def _torus_grid(nu=12, nv=8, major_radius=3.0, minor_radius=1.0):
    pts = []
    for i in range(nu):
        tu = 2 * np.pi * i / nu
        for j in range(nv):
            tv = 2 * np.pi * j / nv
            rho = major_radius + minor_radius * np.cos(tv)
            pts.append([rho * np.cos(tu), rho * np.sin(tu), minor_radius * np.sin(tv)])
    pts = np.array(pts)

    def vid(i, j):
        return (i % nu) * nv + (j % nv)

    tris = []
    for i in range(nu):
        for j in range(nv):
            a, b = vid(i, j), vid(i + 1, j)
            c, d = vid(i + 1, j + 1), vid(i, j + 1)
            tris.append([a, b, c])
            tris.append([a, c, d])
    return pts, np.array(tris, dtype=np.int64)


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
    pts, tris = _torus_grid()
    rng = np.random.default_rng(7)
    bad = tris.copy()
    flip = rng.random(len(bad)) < 0.5
    bad[flip] = bad[flip][:, [0, 2, 1]]
    out, stats = orient_surface_triangles(pts, bad)
    assert stats["conflicts_after"] == 0
    assert _signed_volume(pts, out) > 0


def test_p2_extractor_orients_concave_torus_and_preserves_dof_coordinates(
    monkeypatch,
):
    from radia.bem import sibc_hacapk

    major_radius = 3.0
    pts, outward_tris = _torus_grid(major_radius=major_radius)
    rng = np.random.default_rng(17)
    input_tris = outward_tris.copy()
    flip = rng.random(len(input_tris)) < 0.5
    input_tris[flip] = input_tris[flip][:, [0, 2, 1]]

    p2_nodes = np.empty((len(input_tris), 6, 3))
    for index, tri in enumerate(input_tris):
        corners = pts[tri]
        p2_nodes[index, :3] = corners
        p2_nodes[index, 3] = 0.5 * (corners[0] + corners[1])
        p2_nodes[index, 4] = 0.5 * (corners[1] + corners[2])
        p2_nodes[index, 5] = 0.5 * (corners[2] + corners[0])

    original_nodes = p2_nodes.copy()
    monkeypatch.setattr(
        sibc_hacapk,
        "extract_surface_curved",
        lambda mesh, bnd_label=None, geom_order=2: (
            pts.copy(), input_tris.copy(), np.arange(len(pts)), p2_nodes.copy()
        ),
    )

    (
        vertices,
        triangles,
        _vertex_global,
        oriented_nodes,
        dofs_per_tri,
        _n_dof,
        dof_coords,
    ) = sibc_hacapk.extract_surface_p2_lagrange(object())

    for index, (before, after) in enumerate(zip(input_tris, triangles)):
        if np.array_equal(after, before):
            np.testing.assert_allclose(oriented_nodes[index], original_nodes[index])
        else:
            np.testing.assert_array_equal(after, before[[0, 2, 1]])
            np.testing.assert_allclose(
                oriented_nodes[index], original_nodes[index, [0, 2, 1, 5, 4, 3]]
            )

        center = vertices[after].mean(axis=0)
        normal = np.cross(vertices[after[1]] - vertices[after[0]],
                          vertices[after[2]] - vertices[after[0]])
        azimuth = np.arctan2(center[1], center[0])
        tube_center = major_radius * np.array(
            [np.cos(azimuth), np.sin(azimuth), 0.0]
        )
        assert np.dot(normal, center - tube_center) > 0.0

        np.testing.assert_allclose(
            dof_coords[dofs_per_tri[index]], oriented_nodes[index]
        )


def test_p2_extractor_rejects_open_surface(monkeypatch):
    from radia.bem import sibc_hacapk

    pts = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    tris = np.array([[0, 1, 2]], dtype=np.int64)
    nodes = np.array([[
        pts[0], pts[1], pts[2],
        0.5 * (pts[0] + pts[1]),
        0.5 * (pts[1] + pts[2]),
        0.5 * (pts[2] + pts[0]),
    ]])
    monkeypatch.setattr(
        sibc_hacapk,
        "extract_surface_curved",
        lambda mesh, bnd_label=None, geom_order=2: (
            pts, tris, np.arange(3), nodes
        ),
    )

    with np.testing.assert_raises_regex(ValueError, "closed two-manifold"):
        sibc_hacapk.extract_surface_p2_lagrange(object())


def test_p2_extractor_rejects_multiple_closed_shells(monkeypatch):
    from radia.bem import sibc_hacapk

    first_points, first_tris = _octahedron()
    second_points, second_tris = _octahedron(center=(5.0, 0.0, 0.0))
    points = np.vstack([first_points, second_points])
    triangles = np.vstack([first_tris, second_tris + len(first_points)])
    nodes = np.empty((len(triangles), 6, 3))
    for index, tri in enumerate(triangles):
        corners = points[tri]
        nodes[index] = [
            corners[0], corners[1], corners[2],
            0.5 * (corners[0] + corners[1]),
            0.5 * (corners[1] + corners[2]),
            0.5 * (corners[2] + corners[0]),
        ]
    monkeypatch.setattr(
        sibc_hacapk,
        "extract_surface_curved",
        lambda mesh, bnd_label=None, geom_order=2: (
            points, triangles, np.arange(len(points)), nodes
        ),
    )

    with np.testing.assert_raises_regex(ValueError, "exactly one closed component"):
        sibc_hacapk.extract_surface_p2_lagrange(object())
