"""Sequential weak coupling vs the monolithic mixed solve on the coupled
HDiv-MMM x HCurl-VIM system.

Three locks the coupled suite did not have (study 2026-07-28,
C:/temp/vim_topopt/coupled_weak_vs_monolithic.{py,json}):

1. INDEPENDENT block reconstruction: [A_M K; K^H Z(s)] rebuilt from the public
   parts (magnetic_operator / coupling / eddy_impedance / magnetic_rhs /
   eddy_rhs) must reproduce solve_frequency's solution.  This pins the
   semantic contract s = j*2*pi*f AND the SIBC surface term
   Zs = SkinImpedance(s, sigma) inside the eddy block; a sign / transpose /
   term drift in any block breaks it.
2. The production-style SEQUENTIAL weak iteration (block Gauss-Seidel:
   magnetics with currents frozen -> eddy with the fresh magnetization) must
   converge to the SAME fixed point as the monolithic mixed solve.
3. The divergence mechanism: the iteration operator E = A_M^-1 K Z^-1 K^H
   scales exactly as lambda^2 under K -> lambda*K, the weak iteration
   DIVERGES for rho > 1 (lambda = 1.1*lambda_crit) while the monolithic
   solve keeps a machine-precision residual.  Measured physical margin on
   this configuration: rho ~ 1e-7..1e-12 over 1 Hz..100 kHz.
"""
import numpy as np
import pytest

pytest.importorskip("ngsolve")

import ngsolve as ng  # noqa: E402
from netgen import occ  # noqa: E402

from radia import vim  # noqa: E402

SIGMA = 5.8e7
MU_R = 1001.0
FREQ = 100.0


@pytest.fixture(scope="module")
def coupled():
    box = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(1, 1, 1))
    box.mat("cond")
    for face in box.faces:
        face.name = "skin"
    mesh = ng.Mesh(occ.OCCGeometry(box).GenerateMesh(maxh=2.0))

    fes = ng.HCurl(mesh, order=2, nograds=True)
    u, v = fes.TnT()
    stiffness = ng.BilinearForm(fes)
    stiffness += ng.curl(u) * ng.curl(v) * ng.dx + 0.05 * u * v * ng.dx
    mass = ng.BilinearForm(fes)
    mass += u * v * ng.dx
    port = ng.LinearForm(fes)
    port += ng.CoefficientFunction((-ng.y, ng.x, 0.0)) * v * ng.dx
    with ng.TaskManager():
        stiffness.Assemble()
        mass.Assemble()
        port.Assemble()

    def swirl_vector_potential(points):
        points = np.asarray(points)
        return np.column_stack(
            (-points[:, 1], points[:, 0], np.zeros(points.shape[0])))

    with ng.TaskManager():
        mixed = vim.NgsolveBDMEddyBubbleVIM(
            mesh, fes, stiffness, mass, port,
            (ng.CoefficientFunction((1.0, 0.0, 0.0)),),
            hdiv_order=1, mu_r=MU_R,
            external_fields=(ng.CoefficientFunction((1.0, 0.0, 1.0)),),
            hdiv_max_modes=1, magnetic_materials="cond",
            steps=2, sigma=SIGMA, conductive_materials="cond",
            response_backend="dense", intorder=1,
            port_vector_potentials=(swirl_vector_potential,),
            coupling_kernel_epsilon=0.1,
        )
        solution = mixed.solve_frequency(FREQ)
    return mixed, solution


def _blocks_and_rhs(mixed, s):
    A_M = np.asarray(mixed.magnetic_operator)
    K = np.asarray(mixed.coupling)
    Z = np.asarray(mixed.eddy_impedance(s, surface_impedance=vim.SkinImpedance(s, SIGMA)))
    b_M = np.asarray(mixed.magnetic_rhs).reshape(A_M.shape[0], -1)
    b_J = np.asarray(mixed.eddy_rhs).reshape(Z.shape[0], -1)
    return A_M, K, Z, b_M, b_J


def _solution_vectors(mixed, solution):
    x_m = np.asarray(solution.magnetization_coefficients).reshape(
        mixed.n_hdiv_mmm_modes, -1)
    x_j = np.asarray(solution.eddy_coefficients).reshape(
        mixed.n_hcurl_vim_modes, -1)
    return x_m, x_j


def _gauss_seidel(A_M, K, Z, b_M, b_J, maxit=200, tol=1e-12):
    x_m = np.zeros((A_M.shape[0], b_M.shape[1]), complex)
    x_j = np.zeros((Z.shape[0], b_J.shape[1]), complex)
    for it in range(maxit):
        x_m_new = np.linalg.solve(A_M, b_M - K @ x_j)
        x_j_new = np.linalg.solve(Z, b_J - K.conj().T @ x_m_new)
        step = max(np.linalg.norm(x_m_new - x_m), np.linalg.norm(x_j_new - x_j))
        scale = max(np.linalg.norm(x_m_new), np.linalg.norm(x_j_new), 1e-300)
        x_m, x_j = x_m_new, x_j_new
        if step / scale < tol:
            return x_m, x_j, it + 1
        if step / scale > 1e12:
            return x_m, x_j, -(it + 1)
    return x_m, x_j, -maxit


def _rho(A_M, K, Z):
    E = np.linalg.solve(A_M, K @ np.linalg.solve(Z, K.conj().T))
    return float(np.max(np.abs(np.linalg.eigvals(E))))


def test_independent_block_reconstruction_matches_monolithic(coupled):
    mixed, solution = coupled
    assert solution.residual_relative_norm < 1e-10
    s = 2j * np.pi * FREQ
    A_M, K, Z, b_M, b_J = _blocks_and_rhs(mixed, s)
    x_m, x_j = _solution_vectors(mixed, solution)
    scale = max(np.linalg.norm(b_M), np.linalg.norm(b_J))
    r_m = np.linalg.norm(A_M @ x_m + K @ x_j - b_M) / scale
    r_j = np.linalg.norm(K.conj().T @ x_m + Z @ x_j - b_J) / scale
    assert r_m < 1e-6
    assert r_j < 1e-6


def test_sequential_weak_iteration_reaches_the_monolithic_fixed_point(coupled):
    mixed, solution = coupled
    s = 2j * np.pi * FREQ
    A_M, K, Z, b_M, b_J = _blocks_and_rhs(mixed, s)
    x_m_star, x_j_star = _solution_vectors(mixed, solution)
    x_m, x_j, iters = _gauss_seidel(A_M, K, Z, b_M, b_J)
    assert iters > 0
    err = (np.linalg.norm(x_m - x_m_star) + np.linalg.norm(x_j - x_j_star)) \
        / (np.linalg.norm(x_m_star) + np.linalg.norm(x_j_star))
    assert err < 1e-6


def test_weak_iteration_diverges_past_lambda_crit_monolithic_stays_exact(coupled):
    mixed, _ = coupled
    s = 2j * np.pi * FREQ
    A_M, K, Z, b_M, b_J = _blocks_and_rhs(mixed, s)
    rho_base = _rho(A_M, K, Z)
    assert 0.0 < rho_base < 1e-3          # physical margin on this configuration
    lam_crit = 1.0 / np.sqrt(rho_base)

    # exact lambda^2 scaling of the iteration operator
    np.testing.assert_allclose(_rho(A_M, 0.5 * lam_crit * K, Z), 0.25, rtol=1e-8)
    rho_over = _rho(A_M, 1.1 * lam_crit * K, Z)
    np.testing.assert_allclose(rho_over, 1.21, rtol=1e-8)

    # below critical: weak converges; above: weak diverges, monolithic exact
    _, _, it_under = _gauss_seidel(A_M, 0.5 * lam_crit * K, Z, b_M, b_J)
    assert it_under > 0
    Kl = 1.1 * lam_crit * K
    _, _, it_over = _gauss_seidel(A_M, Kl, Z, b_M, b_J)
    assert it_over < 0                    # diverged
    O = np.block([[A_M, Kl], [Kl.conj().T, Z]])
    b = np.vstack([b_M, b_J])
    x = np.linalg.solve(O, b)
    assert np.linalg.norm(O @ x - b) / np.linalg.norm(b) < 1e-12
