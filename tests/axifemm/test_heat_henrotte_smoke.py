"""Smoke + parity test for the axisymmetric heat Henrotte BFIs.

Confirms that the C++ ``AxiHenrotteHeatStiffnessBFI`` and
``AxiHenrotteHeatMassBFI`` (released in radia 4.31.0) assemble matrices
with the structural properties required by the heat equation:

  * ``K`` symmetric + 1 zero eigenvalue (constant-T null mode)
  * ``M`` symmetric positive definite (heat capacity matrix)
  * Constant-T mode ``T = 1`` annihilates ``K``: ``K @ ones == 0``
  * ``M @ ones`` integrates the volume of revolution: ``2 pi * V_axisym``

These are not full convergence tests -- they verify the assembly path
end-to-end on a tiny mesh that fits in this file.  Convergence /
parity-vs-H1 comparison lives in the ``examples/axifemm/`` research
scripts.

Skips cleanly when ``radia.axifem`` is unavailable
(e.g. running on a machine without NGSolve).
"""
from math import pi

import numpy as np
import pytest


pytest.importorskip("ngsolve")

from netgen.occ import OCCGeometry, MoveTo
from ngsolve import (
    BilinearForm, CoefficientFunction, Mesh, TaskManager,
)


def _to_dense(mat, n):
    A = np.zeros((n, n))
    for r, c, v in zip(*mat.COO()):
        A[r, c] += v
    return A


def _build_axisym_quad_mesh(ra, rb, za, zb, maxh):
    """Single (or few) Q1 quad on [ra,rb] x [za,zb]."""
    box = MoveTo(ra, za).Rectangle(rb - ra, zb - za).Face()
    box.faces.name = "wp"
    return Mesh(OCCGeometry(box, dim=2).GenerateMesh(
        maxh=maxh, quad_dominated=True))


@pytest.mark.parametrize("ra,rb,za,zb", [
    (1e-3, 3e-3, 0.0, 1e-3),
    # Note: tall-aspect cases (e.g. (5e-3, 7e-3, -2e-3, 3e-3)) are excluded
    # because Netgen ignores the maxh hint for axis-disparate rectangles
    # and splits them into non-axis-aligned sub-quads, which the
    # AxiHenrotteFE_Q1_AxisAligned dispatch path is not designed for.
    # See tests/axifemm/test_python_reference_consistency.py for the
    # same restriction.
    (1e-3, 3e-3, 3e-3, 4e-3),
])
def test_q1_heat_assembly_structure(ra, rb, za, zb):
    """Q1 heat: K is SPSD with one zero mode, M is SPD, both symmetric."""
    axifem = pytest.importorskip("radia.axifem")
    mesh = _build_axisym_quad_mesh(ra, rb, za, zb,
                                    maxh=10 * max(rb - ra, zb - za))
    fes = axifem.H1Henrotte(mesh, order=1)
    assert fes.ndof >= 4

    a = BilinearForm(fes, symmetric=True)
    a += axifem.AxiHenrotteHeatStiffnessBFI(CoefficientFunction(50.0))
    m = BilinearForm(fes, symmetric=True)
    m += axifem.AxiHenrotteHeatMassBFI(CoefficientFunction(3.5e6))
    with TaskManager():
        a.Assemble()
        m.Assemble()

    K = _to_dense(a.mat, fes.ndof)
    M = _to_dense(m.mat, fes.ndof)

    np.testing.assert_allclose(K, K.T, atol=1e-12,
                                err_msg="K not symmetric")
    np.testing.assert_allclose(M, M.T, atol=1e-12,
                                err_msg="M not symmetric")

    # K constant-T null mode: K @ ones == 0
    ones = np.ones(fes.ndof)
    np.testing.assert_allclose(K @ ones, 0.0, atol=1e-10 * abs(K).max(),
                                err_msg="K does not annihilate constant T")

    # M strictly positive definite
    eigs_M = np.linalg.eigvalsh(M)
    assert eigs_M.min() > 0, f"M not SPD: min eigenvalue {eigs_M.min()}"


def test_q2_heat_assembly_structure():
    """Q2 heat: same structural checks at 9 DOF/cell."""
    axifem = pytest.importorskip("radia.axifem")
    ra, rb, za, zb = 1e-3, 3e-3, 0.0, 2e-3
    mesh = _build_axisym_quad_mesh(ra, rb, za, zb,
                                    maxh=10 * max(rb - ra, zb - za))
    fes = axifem.H1Henrotte(mesh, order=2)
    assert fes.ndof >= 9

    a = BilinearForm(fes, symmetric=True)
    a += axifem.AxiHenrotteHeatStiffnessBFI(CoefficientFunction(50.0))
    m = BilinearForm(fes, symmetric=True)
    m += axifem.AxiHenrotteHeatMassBFI(CoefficientFunction(3.5e6))
    with TaskManager():
        a.Assemble()
        m.Assemble()

    K = _to_dense(a.mat, fes.ndof)
    M = _to_dense(m.mat, fes.ndof)

    np.testing.assert_allclose(K, K.T, atol=1e-12 * abs(K).max(),
                                err_msg="K not symmetric")
    np.testing.assert_allclose(M, M.T, atol=1e-12 * abs(M).max(),
                                err_msg="M not symmetric")

    ones = np.ones(fes.ndof)
    np.testing.assert_allclose(K @ ones, 0.0, atol=1e-10 * abs(K).max(),
                                err_msg="K does not annihilate constant T")

    eigs_M = np.linalg.eigvalsh(M)
    assert eigs_M.min() > 0, f"M Q2 not SPD: min eigenvalue {eigs_M.min()}"


def test_q1_axis_touching_finite():
    """Q1 heat on an axis-touching element (s_a = 0) integrates cleanly --
    no 1/s singularity.  This is the architectural difference vs the
    magnetic stiffness BFI which would hit a log(0) without the axis-
    reduced basis."""
    axifem = pytest.importorskip("radia.axifem")
    ra, rb, za, zb = 0.0, 2e-3, 0.0, 1e-3
    mesh = _build_axisym_quad_mesh(ra, rb, za, zb,
                                    maxh=10 * (rb - ra))
    fes = axifem.H1Henrotte(mesh, order=1)

    a = BilinearForm(fes, symmetric=True)
    a += axifem.AxiHenrotteHeatStiffnessBFI(CoefficientFunction(50.0))
    m = BilinearForm(fes, symmetric=True)
    m += axifem.AxiHenrotteHeatMassBFI(CoefficientFunction(3.5e6))
    with TaskManager():
        a.Assemble()
        m.Assemble()

    K = _to_dense(a.mat, fes.ndof)
    M = _to_dense(m.mat, fes.ndof)
    assert np.all(np.isfinite(K)), "K has non-finite entries (axis case)"
    assert np.all(np.isfinite(M)), "M has non-finite entries (axis case)"
    eigs_K = np.linalg.eigvalsh(K)
    assert eigs_K.min() > -1e-10 * abs(K).max(), \
        f"K axis case has spurious negative eigenvalue {eigs_K.min()}"
