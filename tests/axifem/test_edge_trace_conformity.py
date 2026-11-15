# -*- coding: utf-8 -*-
r"""Shared-edge trace continuity of the H1Henrotte elements.

Guards: every element family is conforming on the edges it shares with itself
and with the other families: oblique, horizontal (z = const), vertical and
axis-connected edges, triangles next to quads, the closed-form axis-aligned
quads next to the general elements.  For every global DOF a unit coefficient is
evaluated on the shared edge itself, between the interpolation nodes, from
each of the two elements through its own reference coordinates
(``mesh.GetTrafo``).  No point location is involved, so both values come from
different elements by construction; an offset-and-locate probe can put both
points in one element and pass trivially.
Assumes: straight elements; the physical point of each evaluation is checked
against the edge.
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

JUMP_LIMIT = 1.0e-10
REF_VERTICES = {
    ng.ET.TRIG: [(1.0, 0.0), (0.0, 1.0), (0.0, 0.0)],
    ng.ET.QUAD: [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
}

# name: (points, elements, shared edge as two point indices)
MESHES = {
    "quad_quad_oblique": ([(1, 0), (2, 0), (2.5, 1), (1, 1), (3, 0), (3.5, 1)],
                          [(0, 1, 2, 3), (1, 4, 5, 2)], (1, 2)),
    "quad_trig_oblique": ([(1, 0), (2, 0), (2.5, 1), (1, 1), (3, 0)],
                          [(0, 1, 2, 3), (1, 4, 2)], (1, 2)),
    "trig_trig_oblique": ([(1, 0), (2, 0), (2.5, 1), (3, 0)],
                          [(0, 1, 2), (1, 3, 2)], (1, 2)),
    "quad_trig_vertical": ([(1, 0), (2, 0), (2, 1), (1, 1), (3, 0)],
                           [(0, 1, 2, 3), (1, 4, 2)], (1, 2)),
    "quad_trig_horizontal": ([(1, 0), (2, 0), (2, 1), (1, 1), (1.5, 2)],
                             [(0, 1, 2, 3), (3, 2, 4)], (3, 2)),
    "axis_quad_trig_oblique": ([(0.0, 0), (0.2, 0), (0.35, 1), (0.0, 1), (0.6, 0)],
                               [(0, 1, 2, 3), (1, 4, 2)], (1, 2)),
    "axis_connected_trig": ([(0.0, 0), (0.5, 0.3), (0.0, 1), (0.6, 1.2)],
                            [(0, 1, 2), (2, 1, 3)], (1, 2)),
    "axis_connected_trig_trig": ([(0.0, 0), (0.4, 0.5), (0.0, 1), (0.5, -0.4)],
                                 [(0, 1, 2), (0, 3, 1)], (0, 1)),
}

CASES = [
    # (mesh, order, curvedquad)
    *[(m, 2, True) for m in MESHES],
    ("quad_trig_vertical", 2, False),
    ("quad_trig_horizontal", 2, False),
    ("trig_trig_oblique", 1, False),
    ("quad_trig_vertical", 1, False),
    ("quad_trig_horizontal", 1, False),
    ("axis_connected_trig", 1, False),
    ("axis_connected_trig_trig", 1, False),
]


def _mesh(points, elements):
    m = NGMesh(dim=2)
    m.Add(FaceDescriptor(surfnr=1, domin=1, bc=1))
    ids = [m.Add(MeshPoint(Pnt(x, y, 0))) for x, y in points]
    for e in elements:
        m.Add(Element2D(1, [ids[i] for i in e]))
    m.SetMaterial(1, "c")
    return ng.Mesh(m)


def _edge_evaluators(mesh, a, b):
    """For each element holding vertices a and b: (trafo, ref point of t)."""
    out = []
    for el in mesh.Elements(ng.VOL):
        verts = [v.nr for v in el.vertices]
        if a in verts and b in verts:
            ra = np.array(REF_VERTICES[el.type][verts.index(a)])
            rb = np.array(REF_VERTICES[el.type][verts.index(b)])
            trafo = mesh.GetTrafo(ng.ElementId(ng.VOL, el.nr))
            out.append(lambda t, tr=trafo, ra=ra, rb=rb: tr(*((1 - t) * ra + t * rb)))
    assert len(out) == 2, f"edge ({a},{b}) is shared by {len(out)} elements"
    return out


@pytest.mark.parametrize("name,order,curvedquad", CASES,
                         ids=[f"{m}-p{o}{'-curved' if c else ''}" for m, o, c in CASES])
def test_shared_edge_trace_is_continuous(name, order, curvedquad):
    points, elements, (a, b) = MESHES[name]
    mesh = _mesh(points, elements)
    fes = H1Henrotte(mesh, order=order, curvedquad=curvedquad)
    gf = ng.GridFunction(fes)
    left, right = _edge_evaluators(mesh, a, b)
    pa, pb = np.asarray(points[a], float), np.asarray(points[b], float)
    ts = np.linspace(0.03, 0.97, 21)
    for t in ts:
        for ev in (left, right):
            p = np.array(ev(t).point[:2])
            assert np.linalg.norm(p - (pa + t * (pb - pa))) < 1e-12, (name, t, p)
    worst = (0.0, None, None)
    for dof in range(fes.ndof):
        gf.vec[:] = 0.0
        gf.vec[dof] = 1.0
        for t in ts:
            jump = abs(gf(left(t)) - gf(right(t)))
            if jump > worst[0]:
                worst = (jump, dof, float(t))
    assert worst[0] < JUMP_LIMIT, (
        f"{name} order {order}: trace jump {worst[0]:.3e} at DOF {worst[1]}, edge "
        f"parameter {worst[2]}, expected < {JUMP_LIMIT:.0e}")
