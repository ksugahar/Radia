"""Surrogate tests for PeriodicStepReuseSolver.

Guards: one Radia SparseSolv preconditioner (Compact AMG or Compact AMS)
constructed once and reused by native PCG; exact pitch-step eddy-current state
permutation against an independent dense per-step solve; signed permutations;
per-step true residual; loud PCG non-convergence and indefinite breakdown;
loud rejection of non-integer, non-finite or complex inputs and of a missing or
mismatched preconditioner choice, without counting a failed step.
Assumes: linear SPD A = K + M/dt; a periodic scalar DOF ring for AMG and a
lowest-order HCurl curl-curl + mass system on a small cube mesh for AMS.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse import csc_matrix, csr_matrix, diags

import radia.periodic_step_reuse as periodic_step_reuse
from radia._residual_gate import RELATIVE_LIMIT
from radia.periodic_step_reuse import PeriodicStepReuseSolver


def _periodic_ring_matrices(n: int = 12):
    # Symmetric positive-definite 1D periodic ring K and conductivity mass M
    main = np.full(n, 4.0)
    off = np.full(n - 1, -1.0)
    K = diags([off, main, off], [-1, 0, 1], shape=(n, n), format="lil")
    K[0, n - 1] = -1.0
    K[n - 1, 0] = -1.0
    M = diags([np.full(n, 0.5)], [0], shape=(n, n), format="csc")
    return csc_matrix(K), M


def _hcurl_cube_system():
    ng = pytest.importorskip("ngsolve")
    from netgen.occ import Box, OCCGeometry

    from radia.sparsesolv_ngsolve import LowestOrderGradient

    mesh = ng.Mesh(OCCGeometry(Box((0, 0, 0), (1, 1, 1))).GenerateMesh(maxh=0.4))
    fes = ng.HCurl(mesh, order=0)
    u, v = fes.TnT()
    stiffness = ng.BilinearForm(ng.curl(u) * ng.curl(v) * ng.dx).Assemble()
    mass = ng.BilinearForm(u * v * ng.dx).Assemble()

    def scipy_matrix(native):
        values, columns, pointers = native.CSR()
        return csr_matrix((np.array(values), np.array(columns), np.array(pointers)),
                          shape=(native.height, native.width))

    gradient = scipy_matrix(LowestOrderGradient(fes))
    coordinates = np.asarray(mesh.ngmesh.Coordinates())
    return scipy_matrix(stiffness.mat), scipy_matrix(mass.mat), gradient, coordinates


def _counting_builder(monkeypatch):
    constructions = []
    real_builder = periodic_step_reuse._build_preconditioner

    def counting(kind, native, gradient, coordinates):
        constructions.append(kind)
        return real_builder(kind, native, gradient, coordinates)

    monkeypatch.setattr(periodic_step_reuse, "_build_preconditioner", counting)
    return constructions


def _march_against_dense_reference(solver, K, M, dt, perm, sign, steps, seed):
    # Independent reference: dense matrices, explicit permutation matrix P with
    # P[i, perm[i]] = sign[i], and numpy.linalg.solve per step.
    n = K.shape[0]
    A_dense = K.toarray() + M.toarray() / dt
    P = np.zeros((n, n))
    P[np.arange(n), perm] = 1.0 if sign is None else sign
    rng = np.random.default_rng(seed)
    x_reuse = np.zeros(n)
    x_ref = np.zeros(n)
    errors = []
    for step in range(steps):
        rhs = rng.standard_normal(n)
        x_reuse = solver.step(rhs, x_prev=x_reuse)
        x_ref = np.linalg.solve(A_dense, rhs + M.toarray() @ (P @ x_ref) / dt)
        errors.append(np.linalg.norm(x_reuse - x_ref) / np.linalg.norm(x_ref))
        assert solver.last_relative_residual <= RELATIVE_LIMIT, (
            f"step {step}: actual={solver.last_relative_residual}, expected<={RELATIVE_LIMIT}")
    return errors


@pytest.mark.parametrize("signed", [False, True])
def test_amg_matches_dense_reference_with_single_preconditioner(signed, monkeypatch):
    n = 12
    dt = 0.01
    K, M = _periodic_ring_matrices(n)
    perm = np.roll(np.arange(n), -1)  # 1-pitch rotational shift
    sign = np.where(np.arange(n) % 2 == 0, 1.0, -1.0) if signed else None
    constructions = _counting_builder(monkeypatch)
    solver = PeriodicStepReuseSolver(K, M, dt=dt, perm=perm, sign=sign, preconditioner="amg")
    errors = _march_against_dense_reference(solver, K, M, dt, perm, sign, steps=8, seed=42)
    assert max(errors) <= 1e-6, f"actual errors={errors}, expected<=1e-6"
    assert constructions == ["amg"], f"actual constructions={constructions}, expected one 'amg'"
    assert solver.preconditioner_setups == 1
    assert solver.steps_solved == 8


def test_ams_matches_dense_reference_with_single_preconditioner(monkeypatch):
    K, M, gradient, coordinates = _hcurl_cube_system()
    n = K.shape[0]
    dt = 0.1
    rng = np.random.default_rng(7)
    perm = rng.permutation(n)
    sign = rng.choice([-1.0, 1.0], size=n)
    constructions = _counting_builder(monkeypatch)
    solver = PeriodicStepReuseSolver(K, M, dt=dt, perm=perm, sign=sign, preconditioner="ams",
                                     gradient=gradient, coordinates=coordinates)
    errors = _march_against_dense_reference(solver, K, M, dt, perm, sign, steps=4, seed=3)
    assert max(errors) <= 1e-5, f"actual errors={errors}, expected<=1e-5"
    assert constructions == ["ams"], f"actual constructions={constructions}, expected one 'ams'"
    assert solver.steps_solved == 4

    starved = PeriodicStepReuseSolver(K, M, dt=dt, perm=perm, preconditioner="ams",
                                      gradient=gradient, coordinates=coordinates,
                                      max_iterations=1)
    with pytest.raises(RuntimeError, match="did not converge"):
        starved.step(rng.standard_normal(n))
    assert starved.steps_solved == 0, f"actual={starved.steps_solved}, expected=0"


def test_rejects_invalid_construction():
    K, M = _periodic_ring_matrices(6)
    gradient = csr_matrix(np.ones((6, 2)))
    coordinates = np.zeros((2, 3))
    bad = [
        dict(dt=0.0, perm=np.arange(6)),
        dict(dt=float("nan"), perm=np.arange(6)),
        dict(dt=0.01, perm=[0, 0, 2, 3, 4, 5]),
        dict(dt=0.01, perm=np.arange(6), sign=[1, 1, 0, 1, 1, 1]),
        dict(dt=0.01, perm=np.arange(6), sign=[1, 1, np.nan, 1, 1, 1]),
        dict(dt=0.01, perm=np.arange(6), max_iterations=0),
        dict(dt=0.01, perm=np.arange(6), max_iterations=True),
    ]
    for kwargs in bad:
        with pytest.raises(ValueError):
            PeriodicStepReuseSolver(K, M, preconditioner="amg", **kwargs)
    for kwargs in (
        dict(preconditioner="lu"),
        dict(preconditioner="amg", gradient=gradient, coordinates=coordinates),
        dict(preconditioner="ams"),
        dict(preconditioner="ams", gradient=gradient),
        dict(preconditioner="ams", gradient=csr_matrix(np.ones((5, 2))), coordinates=coordinates),
        dict(preconditioner="ams", gradient=gradient, coordinates=np.zeros((3, 3))),
        dict(preconditioner="ams", gradient=gradient, coordinates=np.full((2, 3), np.nan)),
    ):
        with pytest.raises(ValueError):
            PeriodicStepReuseSolver(K, M, dt=0.01, perm=np.arange(6), **kwargs)
    with pytest.raises(TypeError):
        PeriodicStepReuseSolver(K, M, dt=0.01, perm=np.arange(6))  # no preconditioner choice
    for perm in ([0.0, 1.5, 2.0, 3.0, 4.0, 5.0], [0.0, 1.0, 2.0, 3.0, 4.0, 5.0], [True] * 6):
        with pytest.raises(TypeError):
            PeriodicStepReuseSolver(K, M, dt=0.01, perm=perm, preconditioner="amg")
    with pytest.raises(TypeError):
        PeriodicStepReuseSolver(K.astype(complex), M, dt=0.01, perm=np.arange(6),
                                preconditioner="amg")
    K_nan = K.copy()
    K_nan[0, 0] = np.nan
    with pytest.raises(ValueError):
        PeriodicStepReuseSolver(K_nan, M, dt=0.01, perm=np.arange(6), preconditioner="amg")


def test_rejects_non_spd_before_any_preconditioner(monkeypatch):
    K, M = _periodic_ring_matrices(6)
    constructions = _counting_builder(monkeypatch)
    K_asym = K.tolil()
    K_asym[0, 1] = -2.0
    with pytest.raises(ValueError, match="symmetric positive definite"):
        PeriodicStepReuseSolver(csc_matrix(K_asym), M, dt=0.01, perm=np.arange(6),
                                preconditioner="amg")
    with pytest.raises(ValueError, match="positive definite"):
        PeriodicStepReuseSolver(-K, M * 0.0, dt=0.01, perm=np.arange(6), preconditioner="amg")
    assert constructions == [], f"actual constructions={constructions}, expected none"


def test_rejects_bad_step_inputs_without_counting():
    K, M = _periodic_ring_matrices(6)
    solver = PeriodicStepReuseSolver(K, M, dt=0.01, perm=np.roll(np.arange(6), 1),
                                     preconditioner="amg")
    good = np.ones(6)
    for rhs, x_prev, error in (
        (np.array([np.nan, 1, 1, 1, 1, 1]), None, ValueError),
        (np.array([np.inf, 1, 1, 1, 1, 1]), None, ValueError),
        (np.array([1j, 1, 1, 1, 1, 1]), None, TypeError),
        (good, np.array([np.nan, 1, 1, 1, 1, 1]), ValueError),
        (good, np.ones(5), ValueError),
    ):
        with pytest.raises(error):
            solver.step(rhs, x_prev=x_prev)
    assert solver.steps_solved == 0, f"actual={solver.steps_solved}, expected=0"
    solver.step(good)
    assert solver.steps_solved == 1, f"actual={solver.steps_solved}, expected=1"


def test_raises_on_true_residual_without_counting():
    K, M = _periodic_ring_matrices(6)
    solver = PeriodicStepReuseSolver(K, M, dt=0.01, perm=np.arange(6), preconditioner="amg")
    accurate = solver._solve
    solver._solve = lambda b: accurate(b) * (1.0 + 1e-3)
    with pytest.raises(RuntimeError, match="true relative residual"):
        solver.step(np.ones(6))
    assert solver.steps_solved == 0, f"actual={solver.steps_solved}, expected=0"


def test_zero_load_returns_exact_zero():
    K, M = _periodic_ring_matrices(6)
    solver = PeriodicStepReuseSolver(K, M, dt=0.01, perm=np.arange(6), preconditioner="amg")
    x = solver.step(np.zeros(6))
    assert np.array_equal(x, np.zeros(6)), f"actual={x}, expected exact zeros"
    assert solver.last_relative_residual == 0.0, (
        f"actual={solver.last_relative_residual}, expected=0.0")
    assert solver.steps_solved == 1, f"actual={solver.steps_solved}, expected=1"


def test_indefinite_symmetric_raises_loudly():
    # Symmetric, positive diagonal, indefinite (eigenvalues 3 and -1): passes the
    # cheap SPD screens; PCG must break down loudly rather than return a value.
    K = csc_matrix(np.array([[1.0, 2.0], [2.0, 1.0]]))
    M = csc_matrix((2, 2))
    solver = PeriodicStepReuseSolver(K, M, dt=1.0, perm=np.arange(2), preconditioner="amg")
    with pytest.raises(RuntimeError):
        solver.step(np.array([1.0, 0.0]))
    assert solver.steps_solved == 0, f"actual={solver.steps_solved}, expected=0"
