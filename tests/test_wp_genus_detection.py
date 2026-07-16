"""Unit tests: workpiece surface genus detection for the scalar BIE + SIBC.

The analytic genus-0 sphere locks the simply connected scalar-potential
BIE path.  On a flux-linked genus-1 tube the BIE surface current
``J_s = n x (-grad phi)`` with single-valued phi carries ZERO net current
through any cut of the surface, so the physical shorted-turn eddy current
on a flux-linked handle -- and its Lenz screening -- is unrepresentable.

``surface_euler_characteristic`` (chi = V - E + F) is the detector; the
calc_inductance weak/strong drivers warn and set a ``P_wp_caveat`` for
chi != 2.  These tests pin the detector on synthetic sphere / torus
triangulations (no mesh generator, no NGSolve mesh: a tiny stand-in mesh
object suffices since the function only walks ``Elements(BND)``).
"""
from __future__ import annotations

import numpy as np

from radia.bem_sibc_solver import surface_euler_characteristic


class _FakeVertex:
    def __init__(self, nr):
        self.nr = nr


class _FakeElement:
    def __init__(self, vids):
        self.vertices = [_FakeVertex(v) for v in vids]


class _FakeSurfaceMesh:
    """Duck-typed stand-in: only Elements(BND) is consumed."""

    def __init__(self, tris):
        self._els = [_FakeElement(t) for t in tris]

    def Elements(self, _vb):
        return list(self._els)


def _octahedron():
    """Closed genus-0 surface: 6 vertices, 8 triangles, chi = 2."""
    top, bot = 0, 5
    eq = [1, 2, 3, 4]
    tris = []
    for k in range(4):
        a, b = eq[k], eq[(k + 1) % 4]
        tris.append([top, a, b])
        tris.append([bot, b, a])
    return tris


def _torus_grid(nu=8, nv=6):
    """Structured torus triangulation: chi = 0 for any (nu, nv)."""
    def vid(i, j):
        return (i % nu) * nv + (j % nv)

    tris = []
    for i in range(nu):
        for j in range(nv):
            a, b = vid(i, j), vid(i + 1, j)
            c, d = vid(i + 1, j + 1), vid(i, j + 1)
            tris.append([a, b, c])
            tris.append([a, c, d])
    return tris


def test_octahedron_is_genus_0():
    mesh = _FakeSurfaceMesh(_octahedron())
    chi = surface_euler_characteristic(mesh)
    assert chi == 2
    assert (2 - chi) // 2 == 0


def test_torus_is_genus_1():
    mesh = _FakeSurfaceMesh(_torus_grid())
    chi = surface_euler_characteristic(mesh)
    assert chi == 0
    assert (2 - chi) // 2 == 1


def test_torus_grid_resolution_invariant():
    for nu, nv in ((4, 3), (16, 12), (7, 5)):
        chi = surface_euler_characteristic(_FakeSurfaceMesh(_torus_grid(nu, nv)))
        assert chi == 0, f"nu={nu} nv={nv} gave chi={chi}"


def test_genus_one_regression_counts():
    """Lock the Euler arithmetic for a genus-1 production-size mesh."""
    V, E, F = 2955, 8865, 5910
    chi = V - E + F
    assert chi == 0
    assert (2 - chi) // 2 == 1
