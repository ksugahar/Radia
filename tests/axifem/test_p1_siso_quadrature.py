# -*- coding: utf-8 -*-
r"""Order-1 triangle matrices integrate over the (s = r^2, z) triangle.

Guards: the P1 Henrotte triangle is s-isoparametric: psi_i is linear in (s, z)
and the stiffness and sigma mass are integrated over the triangle whose edges
are straight in (s, z), so the assembled space is conforming on that
triangulation.  The former element integrated the same basis over the (r, z)
triangle, where a basis that is linear in (s, z) is not continuous across an
oblique edge.  The expected matrices come from an independent numpy quadrature;
the (r, z)-triangle integral is computed too and must differ, so the check can
fail.
Assumes: straight elements, mu = sigma = 1, no axis vertex.
"""
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

import numpy as np
import ngsolve as ng
from netgen.csg import Pnt
from netgen.meshing import Element2D, FaceDescriptor, Mesh as NGMesh, MeshPoint
from radia.axifem import AxiHenrotteSigmaMassBFI, AxiHenrotteStiffnessBFI, H1Henrotte

POINTS = [(1.0, 0.0), (2.0, 0.0), (2.5, 1.0), (3.0, 0.0), (1.4, 1.3)]
ELEMENTS = [(0, 1, 2), (1, 3, 2), (0, 2, 4)]
RTOL = 1.0e-6     # the 8x8 element rule on 1/s over an element spanning s = 1..9 (1.5e-7)


def _mesh():
    m = NGMesh(dim=2)
    m.Add(FaceDescriptor(surfnr=1, domin=1, bc=1))
    ids = [m.Add(MeshPoint(Pnt(x, y, 0))) for x, y in POINTS]
    for e in ELEMENTS:
        m.Add(Element2D(1, [ids[i] for i in e]))
    m.SetMaterial(1, "c")
    return ng.Mesh(m)


def _triangle_rule(n=24):
    x, w = np.polynomial.legendre.leggauss(n)
    x, w = 0.5 * (x + 1.0), 0.5 * w
    u, v = np.meshgrid(x, x, indexing="ij")
    wu, wv = np.meshgrid(w, w, indexing="ij")
    # Collapse at the third vertex: barycentrics (l0, l1) = (u (1 - v), u v).
    return (u * (1 - v)).ravel(), (u * v).ravel(), (wu * wv * u).ravel()


def _element(rn, zn, in_s):
    """Element K, M of the (s, z)-linear basis over the (s, z) or (r, z) triangle."""
    rn, zn = np.asarray(rn), np.asarray(zn)
    sn = rn ** 2
    V = np.column_stack([np.ones(3), sn, zn])
    C = np.linalg.inv(V)                         # psi_j = C[0,j] + C[1,j] s + C[2,j] z
    l0, l1, w = _triangle_rule()
    lam = np.stack([l0, l1, 1 - l0 - l1])
    if in_s:
        s = sn @ lam; z = zn @ lam
        jac = abs((sn[0] - sn[2]) * (zn[1] - zn[2]) - (sn[1] - sn[2]) * (zn[0] - zn[2]))
        dA_over_r = w * jac / (2.0 * s)          # dr dz / r = ds dz / (2 s)
    else:
        r = rn @ lam; z = zn @ lam; s = r ** 2
        jac = abs((rn[0] - rn[2]) * (zn[1] - zn[2]) - (rn[1] - rn[2]) * (zn[0] - zn[2]))
        dA_over_r = w * jac / r
    psi = C[0][:, None] + C[1][:, None] * s + C[2][:, None] * z
    dr = 2.0 * np.sqrt(s)[None, :] * C[1][:, None]
    dz = np.broadcast_to(C[2][:, None], psi.shape)
    K = (dr * dA_over_r) @ dr.T + (dz * dA_over_r) @ dz.T
    M = (psi * dA_over_r) @ psi.T
    pre = 2.0 * np.pi * np.outer(rn, rn)
    return pre * K, pre * M


def _global(in_s, nd):
    K = np.zeros((nd, nd)); M = np.zeros((nd, nd))
    for e in ELEMENTS:
        rn = [POINTS[i][0] for i in e]; zn = [POINTS[i][1] for i in e]
        ke, me = _element(rn, zn, in_s)
        K[np.ix_(e, e)] += ke; M[np.ix_(e, e)] += me
    return K, M


def _assembled():
    mesh = _mesh()
    fes = H1Henrotte(mesh, order=1)
    one = ng.CoefficientFunction(1.0)
    a = ng.BilinearForm(fes, symmetric=False); a += AxiHenrotteStiffnessBFI(one); a.Assemble()
    m = ng.BilinearForm(fes, symmetric=False); m += AxiHenrotteSigmaMassBFI(one); m.Assemble()
    vnum = [mesh.vertices[i].nr for i in range(len(POINTS))]
    assert vnum == list(range(len(POINTS))), vnum
    return (np.array(a.mat.ToDense()), np.array(m.mat.ToDense()), fes.ndof)


def _rel(a, b):
    return np.linalg.norm(a - b) / np.linalg.norm(b)


def test_p1_matrices_integrate_over_the_s_triangle():
    K, M, nd = _assembled()
    Ks, Ms = _global(True, nd)
    Kr, Mr = _global(False, nd)
    assert _rel(Ks, Kr) > 1e3 * RTOL and _rel(Ms, Mr) > 1e3 * RTOL, (
        f"(s,z) and (r,z) integrals too close to discriminate: K {_rel(Ks, Kr):.3e}, M {_rel(Ms, Mr):.3e}")
    assert _rel(K, Ks) <= RTOL, f"stiffness: rel {_rel(K, Ks):.3e} vs (s,z) quadrature (limit {RTOL})"
    assert _rel(M, Ms) <= RTOL, f"sigma mass: rel {_rel(M, Ms):.3e} vs (s,z) quadrature (limit {RTOL})"

