"""Smoke + parity tests for axisymmetric heat on quadrilateral meshes.

Confirms that the C++ ``AxiHenrotteHeatStiffnessBFI`` and
``AxiHenrotteHeatMassBFI`` (released in radia 4.31.0) assemble matrices
with the structural properties required by the heat equation off axis and at
Q1.  Axis-touching Q2 is deliberately rejected because ``H1Henrotte`` exposes
the magnetic six-function axis basis while the heat matrices require nine
scalar-temperature functions.  The production axis-touching Q2 regression
therefore uses standard NGSolve ``H1(order=2)`` with the ``2 pi r`` weight.

  * ``K`` symmetric + 1 zero eigenvalue (constant-T null mode)
  * ``M`` symmetric positive definite (heat capacity matrix)
  * Constant-T mode ``T = 1`` annihilates ``K``: ``K @ ones == 0``
  * ``M @ ones`` integrates the volume of revolution: ``2 pi * V_axisym``

These are not full convergence tests -- they verify the assembly path
end-to-end on a tiny mesh that fits in this file.  Convergence /
parity-vs-H1 comparison lives in the ``docs/axifem/`` result notebooks and
``validation_test/axifem/`` research checks.

Skips cleanly when ``radia.axifem`` is unavailable
(e.g. running on a machine without NGSolve).
"""
from math import pi

import numpy as np
import pytest


pytest.importorskip("ngsolve")

from netgen.meshing import (
    Element1D, Element2D, FaceDescriptor, Mesh as NgMesh, MeshPoint, Pnt,
)
from ngsolve import (
    BilinearForm, CoefficientFunction, GridFunction, H1, Integrate,
    LinearForm, Mesh, TaskManager, dx, grad, x,
)


def _to_dense(mat, n):
    A = np.zeros((n, n))
    for r, c, v in zip(*mat.COO()):
        A[r, c] += v
    return A


def _build_axisym_quad_mesh(ra, rb, za, zb, maxh):
    """Structured axis-aligned Q1 quads on [ra,rb] x [za,zb].

    Keep the Netgen mesh in memory.  Netgen 2D ``NgMesh.Save()`` currently
    access-violates on Windows, and persistence is irrelevant to these
    element-level regressions.
    """
    nx = max(1, int(np.ceil((rb - ra) / maxh)))
    ny = max(1, int(np.ceil((zb - za) / maxh)))
    ngmesh = NgMesh()
    ngmesh.dim = 2
    ngmesh.SetMaterial(1, "domain")
    for boundary, name in enumerate(("bottom", "right", "top", "left"), 1):
        ngmesh.Add(FaceDescriptor(surfnr=boundary, domin=0, bc=boundary))
        ngmesh.SetBCName(boundary - 1, name)

    pids = [[ngmesh.Add(MeshPoint(Pnt(
        ra + i / nx * (rb - ra),
        za + j / ny * (zb - za),
        0.0,
    ))) for i in range(nx + 1)] for j in range(ny + 1)]

    for j in range(ny):
        for i in range(nx):
            ngmesh.Add(Element2D(1, [
                pids[j][i], pids[j][i + 1],
                pids[j + 1][i + 1], pids[j + 1][i],
            ]))
    for i in range(nx):
        ngmesh.Add(Element1D([pids[0][i], pids[0][i + 1]], index=1))
        ngmesh.Add(Element1D([pids[ny][i], pids[ny][i + 1]], index=3))
    for j in range(ny):
        ngmesh.Add(Element1D([pids[j][nx], pids[j + 1][nx]], index=2))
        ngmesh.Add(Element1D([pids[j][0], pids[j + 1][0]], index=4))
    return Mesh(ngmesh)


@pytest.mark.parametrize("ra,rb,za,zb", [
    (1e-3, 3e-3, 0.0, 1e-3),
    # Note: tall-aspect cases (e.g. (5e-3, 7e-3, -2e-3, 3e-3)) are excluded
    # because Netgen ignores the maxh hint for axis-disparate rectangles
    # and splits them into non-axis-aligned sub-quads, which the
    # AxiHenrotteFE_Q1_AxisAligned dispatch path is not designed for.
    # See tests/axifem/test_python_reference_consistency.py for the
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


def test_q2_axis_touching_henrotte_heat_fails_fast():
    """Both legacy heat BFIs reject the magnetic Q2 axis basis explicitly."""
    axifem = pytest.importorskip("radia.axifem")
    mesh = _build_axisym_quad_mesh(0.0, 1.0, 0.0, 1.0, maxh=2.0)
    fes = axifem.H1Henrotte(mesh, order=2)

    for integrator in (
        axifem.AxiHenrotteHeatStiffnessBFI(CoefficientFunction(2.0)),
        axifem.AxiHenrotteHeatMassBFI(CoefficientFunction(3.0)),
    ):
        form = BilinearForm(fes, symmetric=True)
        form += integrator
        with pytest.raises(Exception, match="axis-touching Q2 AxiHenrotteFE"):
            with TaskManager():
                form.Assemble()


def test_q2_axis_touching_ngsolve_heat_structure_and_manufactured_solve():
    """Standard H1 Q2 is finite on axis and reproduces a quadratic field."""
    radius, height = 1.0, 1.0
    conductivity, capacity, dt = 2.0, 3.0, 0.25
    mesh = _build_axisym_quad_mesh(
        0.0, radius, 0.0, height, maxh=0.5,
    )

    # Structural check without essential boundaries.  H1 order-2 coefficients
    # are hierarchical, so obtain the constant-mode coefficient vector via Set.
    fes_free = H1(mesh, order=2)
    u_free, v_free = fes_free.TnT()
    stiffness = BilinearForm(fes_free, symmetric=True)
    stiffness += conductivity * grad(u_free) * grad(v_free) * 2 * pi * x * dx
    mass = BilinearForm(fes_free, symmetric=True)
    mass += capacity * u_free * v_free * 2 * pi * x * dx
    constant = GridFunction(fes_free)
    with TaskManager():
        constant.Set(CoefficientFunction(1.0))
        stiffness.Assemble()
        mass.Assemble()

    np.testing.assert_allclose(
        [float(constant(mesh(0.0, 0.5))), float(constant(mesh(0.5, 0.5)))],
        1.0,
        atol=1e-13,
    )
    K = _to_dense(stiffness.mat, fes_free.ndof)
    M = _to_dense(mass.mat, fes_free.ndof)
    c = constant.vec.FV().NumPy().copy()
    np.testing.assert_allclose(K, K.T, atol=1e-12 * abs(K).max())
    np.testing.assert_allclose(M, M.T, atol=1e-12 * abs(M).max())
    np.testing.assert_allclose(K @ c, 0.0, atol=1e-11 * abs(K).max())
    assert np.linalg.eigvalsh(M).min() > 0.0
    np.testing.assert_allclose(
        c @ M @ c,
        capacity * pi * radius**2 * height,
        rtol=1e-12,
    )

    # One backward-Euler step with the exact solution T=R^2-r^2.  Its
    # axisymmetric source is q=-k(T_rr+r^-1*T_r)=4*k, and T=0 at r=R.
    fes = H1(mesh, order=2, dirichlet="right")
    u, v = fes.TnT()
    exact = radius**2 - x**2
    source = 4 * conductivity
    system = BilinearForm(fes, symmetric=True)
    system += conductivity * grad(u) * grad(v) * 2 * pi * x * dx
    system += (capacity / dt) * u * v * 2 * pi * x * dx
    rhs = LinearForm(fes)
    rhs += (source + (capacity / dt) * exact) * v * 2 * pi * x * dx
    solution = GridFunction(fes)
    with TaskManager():
        system.Assemble()
        rhs.Assemble()
        solution.vec.data = system.mat.Inverse(
            fes.FreeDofs(), inverse="sparsecholesky",
        ) * rhs.vec

    error = float(Integrate((solution - exact)**2 * 2 * pi * x, mesh, order=10))
    reference = float(Integrate(exact**2 * 2 * pi * x, mesh, order=10))
    assert np.sqrt(error / reference) < 1e-11
