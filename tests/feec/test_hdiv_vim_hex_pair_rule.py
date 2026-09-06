"""Pair-domain Duffy rule: the affine-pair point count and the exact unit-cube constants.

The rule integrates every touching BDM1 host pair on its product domain.  For AFFINE hosts the
integrand converges exponentially, so affine-affine pairs take 6 points per dimension
(``glpair_affine_n``), distorted pairs 8 (``glpair_n``).  The unit cube has closed-form self
energies: cube-cube 1.88231264438961 and face-face 2.9732095982 (times 1/(4 pi) in the Gram).
"""

import numpy as np
import pytest

pytest.importorskip("ngsolve")

import ngsolve as ng  # noqa: E402
import scipy.sparse as sp  # noqa: E402
from ngsolve.meshes import MakeStructured3DMesh  # noqa: E402

from radia.vim import _vim as V  # noqa: E402

CUBE_SELF = 1.88231264438961
SQUARE_SELF = 2.9732095982
FOUR_PI = 4.0 * np.pi


def _dense(mesh, **kw):
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=1)
        B, G, M, _ = V._build_charge_gram_hex(fes, eps=1e-10, build_hmatrix=False, **kw)
        stats = dict(G.stats())
        cb = V._charge_basis_hex(fes, cob_quad=2, materialize_mass=False)
    n = sp.csr_matrix(B).shape[0]
    Gd = np.array([[G.entry(i, j) for j in range(n)] for i in range(n)])
    return Gd, stats, cb


def _constant_charges(cb):
    kind = np.asarray(cb["kind"], int)
    expo = np.asarray(cb["expo"], int).reshape(-1, 3)
    cell = int(np.flatnonzero((kind == 0) & (expo.sum(axis=1) == 0))[0])
    faces = np.flatnonzero((kind == 1) & (expo.sum(axis=1) == 0))
    return cell, faces


def test_affine_pairs_take_six_points_and_reproduce_the_unit_cube_constants():
    Gd, stats, cb = _dense(MakeStructured3DMesh(hexes=True, nx=1, ny=1, nz=1))
    assert stats["hex_glpair_n"] == 8.0 and stats["hex_glpair_affine_n"] == 6.0
    assert stats["hex_pair_duffy_enabled"] == 1.0
    cell, faces = _constant_charges(cb)
    assert abs(Gd[cell, cell] * FOUR_PI - CUBE_SELF) / CUBE_SELF < 2.0e-8      # 4.9e-9 measured at 6 points
    assert abs(Gd[faces[0], faces[0]] * FOUR_PI - SQUARE_SELF) / SQUARE_SELF < 2.0e-8


def test_eight_points_on_affine_pairs_sharpen_to_the_reference_and_agree_with_six():
    G6, _, cb = _dense(MakeStructured3DMesh(hexes=True, nx=1, ny=1, nz=1))
    G8, stats, _ = _dense(MakeStructured3DMesh(hexes=True, nx=1, ny=1, nz=1), glpair_affine_n=8)
    assert stats["hex_glpair_affine_n"] == 8.0
    cell, faces = _constant_charges(cb)
    assert abs(G8[cell, cell] * FOUR_PI - CUBE_SELF) / CUBE_SELF < 1.0e-11
    scale = np.sqrt(np.outer(np.abs(np.diag(G8)), np.abs(np.diag(G8))))
    assert (np.abs(G6 - G8) / scale).max() < 1.0e-7          # every touching block agrees to the 6-point accuracy


def test_distorted_pairs_keep_the_eight_point_rule():
    """A lattice sheared by bilinear vertex terms: every cell is a non-parallelepiped trilinear hex and
    every face a non-parallelogram (a tapered sector still has rectangular constant-angle faces, and
    a purely quadratic vertex map on a 2-cell lattice still gives parallelepipeds), so
    glpair_affine_n must not change any entry."""
    mesh = lambda: MakeStructured3DMesh(  # noqa: E731
        hexes=True, nx=2, ny=2, nz=1,
        mapping=lambda x, y, z: (x + 0.2 * x * y, y + 0.2 * y * z, z + 0.2 * z * x))
    Ga, _, _ = _dense(mesh())
    Gb, stats, _ = _dense(mesh(), glpair_affine_n=3)
    assert stats["hex_glpair_affine_n"] == 3.0
    assert np.array_equal(Ga, Gb)
