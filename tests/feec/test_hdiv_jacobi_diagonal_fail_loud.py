"""The exact Jacobi diagonal of ``inv_chi*M + B^T G B`` must raise when it is not positive.

``SolveConfiguredLinearMaterialAutoPrec`` (and its multi-RHS twin) used to replace a non-positive or
non-finite diagonal entry by 1.0 and carry on.  A non-positive exact diagonal means the charge Gram lost
positive definiteness on that DOF's own charge support -- the defect the CG breakdown guard reports
later -- so the preconditioner builder must report it (No-Fallbacks, ESRF #6 2026-09-05).
"""

import numpy as np
import pytest

pytest.importorskip("ngsolve")

import ngsolve as ng  # noqa: E402
from ngsolve.meshes import MakeStructured3DMesh  # noqa: E402

from radia.vim import _vim as V  # noqa: E402


def _configured_gram():
    sp = pytest.importorskip("scipy.sparse")
    mesh = MakeStructured3DMesh(hexes=True, nx=2, ny=2, nz=1)
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=1)
        B, G, M, _ = V._build_charge_gram_hex(fes, eps=1e-10)
        Bs = sp.csr_matrix(B)
        n_face = Bs.shape[1]
        G.configure_charge_map(np.ascontiguousarray(Bs.indptr, np.int32),
                               np.ascontiguousarray(Bs.indices, np.int32),
                               np.ascontiguousarray(Bs.data, float), int(n_face))
    return G, sp.csr_matrix(M), n_face


def _configure_mass(G, coo, n_face):
    G.configure_mass_matrix(np.ascontiguousarray(coo.row, np.int32), np.ascontiguousarray(coo.col, np.int32),
                            np.ascontiguousarray(coo.data, float), int(n_face))


def test_negative_jacobi_diagonal_raises():
    sp = pytest.importorskip("scipy.sparse")
    G, Md, n_face = _configured_gram()
    # a large negative diagonal on every DOF: value = inv_chi*(-big) + ndiag < 0
    bad = sp.coo_matrix(Md - 1.0e3 * abs(Md.diagonal()).max() * sp.eye(n_face, format="coo"))
    _configure_mass(G, bad, n_face)
    rhs = np.ones(n_face)
    with ng.TaskManager():
        with pytest.raises(RuntimeError, match="Jacobi diagonal.*not positive"):
            G.solve_configured_linear_material_auto_prec(1.0, rhs, 1e-8, 100)
        with pytest.raises(RuntimeError, match="Jacobi diagonal.*not positive"):
            G.solve_configured_linear_material_auto_prec_many(1.0, np.ones((2, n_face)), 1e-8, 100)


def test_positive_jacobi_diagonal_still_solves():
    G, Md, n_face = _configured_gram()
    _configure_mass(G, sp_coo(Md), n_face)
    rhs = np.asarray(Md @ np.ones(n_face), float)
    with ng.TaskManager():
        res = G.solve_configured_linear_material_auto_prec(1.0e-3, np.ascontiguousarray(rhs), 1e-8, 2000)
    assert res["prec_min"] > 0.0 and res["iters"] < 2000


def sp_coo(matrix):
    sp = pytest.importorskip("scipy.sparse")
    return sp.coo_matrix(matrix)


@pytest.mark.parametrize("value", [np.nan, np.inf])
@pytest.mark.parametrize("target", ["rhs", "x0"])
def test_nonfinite_scalar_pcg_cannot_report_convergence(value, target):
    G, Md, n_face = _configured_gram()
    _configure_mass(G, sp_coo(Md), n_face)
    rhs = np.ones(n_face)
    initial = np.zeros(n_face)
    (rhs if target == "rhs" else initial)[0] = value
    # Zero budget reaches the final true-residual path without a curvature
    # check masking the NaN-to-zero norm regression.
    with ng.TaskManager(), pytest.raises(RuntimeError, match="non-finite .*norm"):
        G.solve_configured_linear_material_auto_prec(1.0, rhs, 1e-8, 0, x0=initial)


def test_indefinite_operator_with_positive_jacobi_diagonal_raises_in_batched_pcg():
    """Batched multi-RHS PCG keeps the scalar SPD contract (2026-09-28 review).

    The Jacobi diagonal stays positive, but a rank-1 shift of the mass matrix moves the smallest
    eigenvalue of inv_chi*M + B^T G B to -1 % of its magnitude.  The batched recurrence must report
    the SPD violation instead of stepping along non-positive curvature.  This locks the contract; the
    pre-fix build reached the same error through its scalar finish on this fixture.
    """
    sp = pytest.importorskip("scipy.sparse")
    mesh = MakeStructured3DMesh(hexes=True, nx=2, ny=2, nz=1)
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=1)
        B, G, M, _ = V._build_charge_gram_hex(fes, eps=1e-10)
        Bs = sp.csr_matrix(B)
        n_face = Bs.shape[1]
        G.configure_charge_map(np.ascontiguousarray(Bs.indptr, np.int32),
                               np.ascontiguousarray(Bs.indices, np.int32),
                               np.ascontiguousarray(Bs.data, float), int(n_face))
    Bd = Bs.toarray()
    N = np.column_stack([Bd.T @ np.asarray(G.matvec_sym(Bd[:, j].tolist()), float)
                         for j in range(n_face)])
    mass = sp.csr_matrix(M).toarray()
    lam, vec = np.linalg.eigh(mass + 0.5*(N + N.T))
    assert lam[0] > 0.0
    shifted = mass - 1.01*lam[0]*np.outer(vec[:, 0], vec[:, 0])
    operator = shifted + 0.5*(N + N.T)
    assert np.min(np.diag(operator)) > 0.0 and np.linalg.eigvalsh(operator)[0] < 0.0
    _configure_mass(G, sp.coo_matrix(shifted), n_face)
    rhs = np.ascontiguousarray(np.random.default_rng(7).standard_normal((2, n_face)))
    with ng.TaskManager(), pytest.raises(RuntimeError, match="not SPD"):
        G.solve_configured_linear_material_auto_prec_many(1.0, rhs, 1e-10, 2000)
