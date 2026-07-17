"""Goldens: radia.cohomology surface complex (H^1 of triangulated surfaces).

Pure numpy/scipy: the combinatorial Hodge-Laplacian machinery of the
gmsh-free cohomology engine, applied to the triangle-surface complex
(C0 -> C1 -> C2).  Closed orientable surface: b1 = 2 * genus.  This is
the topology engine behind the genus-1 loop DOF of
radia.bem_loop_extension.
"""
from __future__ import annotations

import math

import numpy as np

from radia.cohomology import (surface_chain_complex, surface_cohomology,
                              surface_homology_loops)


def _tetrahedron():
    """Boundary of a tetrahedron: the minimal genus-0 closed surface."""
    pts = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    tris = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])
    return pts, tris


def _torus(nu=14, nw=8, R=3.0, r=1.0):
    """Structured triangulated torus: V=nu*nw, E=3*nu*nw, F=2*nu*nw
    (Euler chi = 0, genus 1)."""
    pts = []
    for i in range(nu):
        u = 2.0 * math.pi * i / nu
        for j in range(nw):
            w = 2.0 * math.pi * j / nw
            pts.append(((R + r * math.cos(w)) * math.cos(u),
                        (R + r * math.cos(w)) * math.sin(u),
                        r * math.sin(w)))

    def idx(i, j):
        return (i % nu) * nw + (j % nw)

    tris = []
    for i in range(nu):
        for j in range(nw):
            a, b = idx(i, j), idx(i + 1, j)
            c, d = idx(i + 1, j + 1), idx(i, j + 1)
            tris += [(a, b, c), (a, c, d)]
    return np.asarray(pts), np.asarray(tris, dtype=np.int64)


def _windings(pts, loop, R):
    """(toroidal, poloidal) winding numbers of a vertex loop on the
    standard torus (implicitly closed: append loop[0])."""
    p = pts[list(loop) + [loop[0]]]
    u = np.unwrap(np.arctan2(p[:, 1], p[:, 0]))
    rho = np.sqrt(p[:, 0] ** 2 + p[:, 1] ** 2)
    w = np.unwrap(np.arctan2(p[:, 2], rho - R))
    return ((u[-1] - u[0]) / (2 * math.pi),
            (w[-1] - w[0]) / (2 * math.pi))


def test_genus0_has_no_loops():
    _, tris = _tetrahedron()
    b1, harm, ctx, _ = surface_cohomology(tris)
    assert b1 == 0 and harm.shape[1] == 0
    assert surface_homology_loops(tris) == []


def test_chain_complex_property_and_euler():
    _, tris = _torus()
    G, Curl, eidx, EI, V, E, F = surface_chain_complex(tris)
    assert abs(Curl @ G).sum() < 1e-12          # d1 . d0 = 0
    assert V - E + F == 0                       # Euler chi of the torus


def test_torus_has_two_independent_generators():
    pts, tris = _torus()
    b1, _, _, _ = surface_cohomology(tris)
    assert b1 == 2                              # 2 * genus

    loops = surface_homology_loops(tris, nv=len(pts))
    assert len(loops) == 2

    # each loop is a CLOSED edge path on the mesh
    edges = {(min(a, b), max(a, b)) for t in tris
             for a, b in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0]))}
    for lp in loops:
        assert len(lp) == len(set(lp)), "loop revisits a vertex"
        ring = list(lp) + [lp[0]]
        for a, b in zip(ring[:-1], ring[1:]):
            assert (min(a, b), max(a, b)) in edges, "non-edge step in loop"

    # the two homology classes are independent: the (toroidal, poloidal)
    # winding matrix is unimodular-ish (nonzero integer determinant)
    W = np.array([_windings(pts, lp, R=3.0) for lp in loops])
    Wi = np.rint(W)
    assert np.allclose(W, Wi, atol=1e-6), f"non-integer windings {W}"
    assert abs(round(np.linalg.det(Wi))) >= 1, f"dependent classes {Wi}"


def test_isolated_vertices_keep_id_space():
    """nv > tris.max()+1 (isolated vertices) must not shift loop ids."""
    pts, tris = _torus(nu=10, nw=6)
    loops_a = surface_homology_loops(tris, nv=len(pts))
    loops_b = surface_homology_loops(tris, nv=len(pts) + 5)
    assert [list(l) for l in loops_a] == [list(l) for l in loops_b]
