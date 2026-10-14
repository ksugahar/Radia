"""A touching HEX pair without a conforming shared entity takes the graded near family.

A 1x1x1 cell against the x=2 face of a 2x2x2 cell shares four Q2 lattice nodes and no whole vertex,
edge or face, so the pair-domain Duffy rule does not apply.  Before 2026-09-28 the pair then fell
to the plain near-band product rule, whose integrand is singular on the contact face; it now takes
the graded near family like every other touching pair.

Reference: the constant-density interaction integral of two boxes, with the inner potential of the
source box in closed form and a Gauss outer rule over the target box.  Two well-separated boxes
calibrate the Gram's charge normalization, so only the touching pair's quadrature is tested
(measured relative error 1.9e-6 fixed, 1.8e-5 before).
"""

import numpy as np
import pytest

import radia._radia_pybind as _rp
from radia.vim import _vim as V

BIG = ((0.0, 0.0, 0.0), (2.0, 2.0, 2.0))
TOUCHING = ((2.0, 0.0, 0.0), (3.0, 1.0, 1.0))
FAR = (((6.0, 0.0, 0.0), (7.0, 1.0, 1.0)), ((9.0, 0.0, 0.0), (10.0, 1.0, 1.0)))


def _q2_nodes(lo, hi):
    return [lo[d] + (hi[d] - lo[d])*i/2.0
            for iz in range(3) for iy in range(3) for ix in range(3)
            for d, i in ((0, ix), (1, iy), (2, iz))]


def _cell_constant_gram(boxes):
    """HEX charge Gram of one constant cell charge per box, with the BDM1 production rules."""
    glo, gwo = V._g01(4)
    gli, gwi = V._g01(5)
    gln, gwn = V._g01(8)
    glp, gwp = V._g01(6)
    tet_p = np.asarray(V._SYM5_TET[0]).ravel().tolist()
    tri_p = np.asarray(V._SYM5_TRI[0]).ravel().tolist()
    n = len(boxes)
    return _rp._ChargeGramHMatrix(
        hex_cell_nodes=sum((_q2_nodes(*box) for box in boxes), []), quad_face_nodes=[],
        n_el=n, n_bf=0, charge_host=list(range(n)), charge_kind=[0]*n, charge_expo=[0]*(3*n),
        sym_tet_pts=tet_p, sym_tet_w=list(V._SYM5_TET[1]),
        sym_tri_pts=tri_p, sym_tri_w=list(V._SYM5_TRI[1]),
        gl_out=list(glo), gw_out=list(gwo), gl_in=list(gli), gw_in=list(gwi),
        far_tet_pts=tet_p, far_tet_w=list(V._SYM5_TET[1]),
        far_tri_pts=tri_p, far_tri_w=list(V._SYM5_TRI[1]),
        near_grade=1.0, far_inner_factor=1.0, eps=1e-12, leaf=64, eta=2.0, build=False,
        gl_near=list(gln), gw_near=list(gwn), near_inner_exact=True,
        gl_in_self=list(gln), gw_in_self=list(gwn), gl_pair=list(glp), gw_pair=list(gwp),
        gl_pair_affine=list(glp), gw_pair_affine=list(gwp))


def _box_antiderivative(u, v, w):
    """Triple antiderivative of 1/r (u, v, w are never zero at the outer Gauss nodes used here)."""
    r = np.sqrt(u*u + v*v + w*w)
    return (v*w*np.log(u + r) + u*w*np.log(v + r) + u*v*np.log(w + r)
            - 0.5*u*u*np.arctan(v*w/(u*r)) - 0.5*v*v*np.arctan(u*w/(v*r))
            - 0.5*w*w*np.arctan(u*v/(w*r)))


def _box_potential(box, x):
    lo, hi = box
    total = 0.0
    for i, a in enumerate((lo[0], hi[0])):
        for j, b in enumerate((lo[1], hi[1])):
            for k, c in enumerate((lo[2], hi[2])):
                total = total + (-1.0)**(i + j + k + 1)*_box_antiderivative(
                    a - x[..., 0], b - x[..., 1], c - x[..., 2])
    return total


def _box_box(source, target, n=32):
    g, w = np.polynomial.legendre.leggauss(n)
    lo, hi = np.asarray(target[0]), np.asarray(target[1])
    axes = [0.5*(lo[d] + hi[d]) + 0.5*(hi[d] - lo[d])*g for d in range(3)]
    points = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1)
    weights = np.einsum("i,j,k->ijk", *[0.5*(hi[d] - lo[d])*w for d in range(3)])
    return float(np.sum(weights*_box_potential(source, points)))


def test_box_potential_reference_is_the_newton_potential():
    far_point = np.array([[30.0, 1.0, 1.0]])
    assert _box_potential(BIG, far_point)[0] == pytest.approx(8.0/np.linalg.norm(far_point[0] - 1.0),
                                                              rel=2e-3)
    g, w = np.polynomial.legendre.leggauss(60)
    grid = np.stack(np.meshgrid(g + 1.0, g + 1.0, g + 1.0, indexing="ij"), axis=-1)
    x0 = np.array([2.3, 0.4, 0.6])
    brute = float(np.sum(np.einsum("i,j,k->ijk", w, w, w)/np.linalg.norm(grid - x0, axis=-1)))
    assert _box_potential(BIG, x0[None])[0] == pytest.approx(brute, rel=1e-10)


def test_nonconforming_touching_hex_pair_matches_the_analytic_interaction():
    gram = _cell_constant_gram([BIG, TOUCHING, *FAR])
    calibration = [gram.entry(0, 2 + k)/_box_box(BIG, box) for k, box in enumerate(FAR)]
    assert calibration[0] == pytest.approx(calibration[1], rel=1e-7)
    reference = calibration[0]*_box_box(BIG, TOUCHING)
    assert gram.entry(0, 1) == pytest.approx(reference, rel=5e-6)
    assert gram.stats()["hex_pair_nonconforming"] > 0
