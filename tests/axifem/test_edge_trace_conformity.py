# -*- coding: utf-8 -*-
r"""Shared-edge trace continuity of the general order-2 Henrotte elements.

Guards: the P2 triangle and the curvedquad Q2 are conforming on edges that are
oblique in (r, z), next to each other and to themselves.  For every global DOF a
unit coefficient is evaluated just inside both elements along the shared edge,
between the interpolation nodes.  A conforming space gives a jump of the order
of the offset times the gradient; the former physical-(s, z)-monomial basis gave
O(1) jumps on oblique edges (quad 0.5, triangle 0.02).
Assumes: straight elements, so the (r, z) and (s, z) element maps share the edge.
"""
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

import numpy as np
import ngsolve as ng
import pytest
from netgen.csg import Pnt
from netgen.meshing import Element2D, FaceDescriptor, Mesh as NGMesh, MeshPoint
from radia.axifem import H1Henrotte

OFFSET = 1.0e-6
JUMP_LIMIT = 1.0e-4     # offset 1e-6 times basis gradients of order 10

CASES = {
    "quad_quad_oblique": ([(1, 0), (2, 0), (2.5, 1), (1, 1), (3, 0), (3.5, 1)],
                          [(0, 1, 2, 3), (1, 4, 5, 2)], (2.0, 0.0), (2.5, 1.0)),
    "quad_trig_oblique": ([(1, 0), (2, 0), (2.5, 1), (1, 1), (3, 0)],
                          [(0, 1, 2, 3), (1, 4, 2)], (2.0, 0.0), (2.5, 1.0)),
    "trig_trig_oblique": ([(1, 0), (2, 0), (2.5, 1), (3, 0)],
                          [(0, 1, 2), (1, 3, 2)], (2.0, 0.0), (2.5, 1.0)),
    "quad_trig_aligned": ([(1, 0), (2, 0), (2, 1), (1, 1), (3, 0)],
                          [(0, 1, 2, 3), (1, 4, 2)], (2.0, 0.0), (2.0, 1.0)),
    "near_axis_oblique": ([(0.0, 0), (0.2, 0), (0.35, 1), (0.0, 1), (0.6, 0)],
                          [(0, 1, 2, 3), (1, 4, 2)], (0.2, 0.0), (0.35, 1.0)),
}


def _mesh(points, elements):
    m = NGMesh(dim=2)
    m.Add(FaceDescriptor(surfnr=1, domin=1, bc=1))
    ids = [m.Add(MeshPoint(Pnt(x, y, 0))) for x, y in points]
    for e in elements:
        m.Add(Element2D(1, [ids[i] for i in e]))
    m.SetMaterial(1, "c")
    return ng.Mesh(m)


@pytest.mark.parametrize("name", sorted(CASES))
def test_shared_edge_trace_is_continuous(name):
    points, elements, a, b = CASES[name]
    mesh = _mesh(points, elements)
    fes = H1Henrotte(mesh, order=2, curvedquad=True)
    gf = ng.GridFunction(fes)
    a, b = np.asarray(a, float), np.asarray(b, float)
    t = b - a
    n = np.array([t[1], -t[0]]) / np.linalg.norm(t)
    worst = (0.0, None, None)
    for dof in range(fes.ndof):
        gf.vec[:] = 0.0
        gf.vec[dof] = 1.0
        for s in np.linspace(0.02, 0.98, 25):
            p = a + s * t
            jump = abs(gf(mesh(*(p - OFFSET * n))) - gf(mesh(*(p + OFFSET * n))))
            if jump > worst[0]:
                worst = (jump, dof, float(s))
    assert worst[0] < JUMP_LIMIT, (
        f"{name}: trace jump {worst[0]:.3e} at DOF {worst[1]}, edge parameter "
        f"{worst[2]}, expected < {JUMP_LIMIT:.0e}")
