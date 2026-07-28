"""Adjoint-gradient and reciprocity locks for radia.isochronous_topopt.

Stage-1 gate of docs/hdiv_vim/TOPOLOGY_OPTIMIZATION_ISOCHRONOUS.md: the
per-element density adjoint gradient must match central finite differences to
the 1e-6 class, and the dipole-array reciprocity load must reproduce the
independent C++ analytic charge evaluator.  Bands are set from measured values
with 2-3 decades of margin (research run 2026-07-28,
C:/temp/vim_topopt/stage1_adjoint_gate.py: directional 8.1e-10, per-element
worst 3.9e-7 at the FD noise floor, reciprocity 1.1e-10 at bonus_intorder=10).
"""
from math import pi
from types import SimpleNamespace

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("radia")

from netgen.occ import OCCGeometry, Pnt, Sphere  # noqa: E402
from ngsolve import HDiv, InnerProduct, Mesh, SetNumThreads, TaskManager  # noqa: E402

from radia.isochronous_topopt import (  # noqa: E402
    CHI_MIN, MU0, DensityAdjointVIM, demag_field_from_solution,
    density_gradient_from_s_gradient, density_to_s, field_functional_load,
    gradient_pair_points, uniform_field_load,
)


@pytest.fixture(scope="module")
def problem():
    SetNumThreads(4)
    with TaskManager():
        mesh = Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), 1.0)).GenerateMesh(maxh=0.45))
        fes = HDiv(mesh, order=1)
        prob = DensityAdjointVIM(fes, eps=1e-7)
        rng = np.random.default_rng(1234)
        s0 = np.exp(rng.uniform(np.log(1e-2), 0.0, prob.n_el))
        theta = np.linspace(0.0, 2.0 * pi, 6, endpoint=False)
        orbit = np.stack([1.6 * np.cos(theta), 1.6 * np.sin(theta),
                          0.3 * np.ones_like(theta)], axis=1)
        w_orbit = np.array([1.0, -0.7, 1.3, 0.5, -1.1, 0.9])
        pts, wts = gradient_pair_points(orbit, w_orbit, delta=0.1, axis=0)
        f_state = uniform_field_load(fes, (0.0, 0.0, 1.0))
        f_adj = field_functional_load(fes, pts, wts, axis=2, scale=MU0,
                                      bonus_intorder=10)
        base = prob.objective_and_gradient(s0, f_state, f_adj)
    return SimpleNamespace(mesh=mesh, fes=fes, prob=prob, s0=s0,
                           dipole_points=pts, dipole_weights=wts,
                           f_state=f_state, f_adj=f_adj, base=base, rng=rng)


def _objective_at(problem, s):
    with TaskManager():
        gf, _ = problem.prob.solve(s, problem.f_state)
    return float(InnerProduct(problem.f_adj.vec, gf.vec))


def test_adjoint_gradient_directional_fd(problem):
    """All-element check: grad . ds == central FD along a random direction.

    h = 1e-4: the FD signal |2 h grad.ds| must sit well above the CG solution
    noise (tol * J ~ 1e-19); at h = 1e-5 this direction's derivative is
    cancellation-small and the comparison hits the solver noise floor.
    """
    d = np.random.default_rng(7).standard_normal(problem.prob.n_el)
    ds = problem.s0 * d              # multiplicative direction keeps s > 0
    predicted = float(problem.base.gradient @ ds)
    h = 1e-4
    fd = (_objective_at(problem, problem.s0 + h * ds)
          - _objective_at(problem, problem.s0 - h * ds)) / (2.0 * h)
    rel = abs(fd - predicted) / abs(predicted)
    assert rel < 1e-6, (fd, predicted, rel)


def test_adjoint_gradient_per_element_fd(problem):
    """Per-element FD on the strongest-gradient element plus random picks."""
    g = problem.base.gradient
    picks = {int(np.abs(g).argmax())}
    rng = np.random.default_rng(21)
    while len(picks) < 4:
        picks.add(int(rng.integers(0, problem.prob.n_el)))
    for e in sorted(picks):
        h = 1e-4 * problem.s0[e]
        sp = problem.s0.copy(); sp[e] += h
        sm = problem.s0.copy(); sm[e] -= h
        fd = (_objective_at(problem, sp) - _objective_at(problem, sm)) / (2.0 * h)
        rel = abs(fd - g[e]) / max(abs(g[e]), abs(fd))
        assert rel < 1e-5, (e, g[e], fd, rel)


def test_reciprocity_matches_cpp_field_evaluator(problem):
    """f_adj^T m == MU0 * sum_i w_i H_d[m]_z(x_i) via the independent C++ path."""
    with TaskManager():
        H_d = demag_field_from_solution(problem.prob.demag, problem.base.gfM,
                                        problem.dipole_points)
    independent = float(MU0 * np.dot(problem.dipole_weights, H_d[:, 2]))
    rel = abs(problem.base.objective - independent) / abs(independent)
    assert rel < 1e-7, (problem.base.objective, independent, rel)


def test_warm_start_reproduces_and_saves_iterations(problem):
    """Warm-started re-solve at the same design matches and converges fast."""
    with TaskManager():
        again = problem.prob.objective_and_gradient(
            problem.s0, problem.f_state, problem.f_adj, warm=problem.base)
    rel = abs(again.objective - problem.base.objective) / abs(problem.base.objective)
    assert rel < 1e-9
    assert again.state_iterations <= max(3, problem.base.state_iterations // 5)


def test_dipole_point_inside_mesh_raises(problem):
    with pytest.raises(ValueError, match="INSIDE"):
        field_functional_load(problem.fes, [[0.0, 0.0, 0.0]], [1.0])


def test_wrong_design_vector_raises(problem):
    with pytest.raises(ValueError, match="entries"):
        problem.prob.solve(np.ones(3), problem.f_state)
    with pytest.raises(ValueError, match="positive"):
        problem.prob.solve(np.zeros(problem.prob.n_el), problem.f_state)


def test_density_mapping_floor_and_chain_rule():
    chi_iron = 100.0
    s = density_to_s(np.array([0.0, 1.0, 0.5]), chi_iron)
    assert s[0] == pytest.approx(1.0 / CHI_MIN)
    assert s[1] == pytest.approx(1.0 / chi_iron)
    rho = np.array([0.1, 0.4, 0.9])
    grad_s = np.array([2.0, -1.0, 0.5])
    h = 1e-7
    fd = (grad_s * (density_to_s(rho + h, chi_iron)
                    - density_to_s(rho - h, chi_iron)) / (2.0 * h))
    chain = density_gradient_from_s_gradient(rho, grad_s, chi_iron)
    assert np.allclose(chain, fd, rtol=1e-5)
    with pytest.raises(ValueError, match=r"\[0, 1\]"):
        density_to_s(np.array([1.5]), chi_iron)


def test_gradient_pair_points_contract():
    pts, wts = gradient_pair_points([[1.0, 2.0, 3.0]], [4.0], delta=0.2, axis=1)
    assert pts.shape == (2, 3) and wts.shape == (2,)
    np.testing.assert_allclose(pts[0], [1.0, 2.1, 3.0])
    np.testing.assert_allclose(pts[1], [1.0, 1.9, 3.0])
    np.testing.assert_allclose(wts, [20.0, -20.0])


def test_element_volumes_sum_to_mesh_volume(problem):
    with TaskManager():
        direct = float(ng.Integrate(ng.CoefficientFunction(1.0), problem.mesh))
    total = float(problem.prob.element_volumes.sum())
    assert abs(total - direct) / direct < 1e-12
    deficit = (4.0 * pi / 3.0 - total) / (4.0 * pi / 3.0)
    assert 0.0 < deficit < 0.12    # measured 7.0% straight-facet deficit at maxh=0.45
