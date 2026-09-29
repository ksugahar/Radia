"""Direct axisymmetric solves: residual, source scaling and constant-BH limit."""
from types import SimpleNamespace

import numpy as np
import pytest

from radia.panels import calc_axisym_volumetric as solver


def test_bh_secant_is_continuous_at_first_knot_and_rejects_extrapolation():
    bh = [(0., 0.), (100., .1), (1000., .5)]
    low = solver._mu_r_from_b_table(np.array([0., .05, .1]), bh)
    np.testing.assert_allclose(low, .001 / solver.MU0)
    near = solver._mu_r_from_b_table(np.array([.1 - 1e-10, .1 + 1e-10]), bh)
    assert near[1] == pytest.approx(near[0], rel=1e-8)
    for invalid in (-.01, .5001, np.nan):
        with pytest.raises(ValueError, match="outside the B-H table"):
            solver._mu_r_from_b_table(np.array([invalid]), bh)


@pytest.mark.parametrize("complex_space", [False, True])
def test_direct_correction_can_reuse_a_converged_solution(complex_space):
    from netgen.geom2d import unit_square
    from ngsolve import Mesh, H1, BilinearForm, LinearForm, GridFunction, TaskManager, dx, grad, x, BND
    from radia.panels.calc_common import apply_fe_inverse

    mesh = Mesh(unit_square.GenerateMesh(maxh=.2))
    fes = H1(mesh, order=2, complex=complex_space, dirichlet=".*")
    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    a += (grad(u) * grad(v) + (1j if complex_space else 1) * u * v) * dx
    f = LinearForm(fes)
    f += v * dx
    solution = GridFunction(fes)
    solution.Set(1 + x, BND)
    with TaskManager():
        a.Assemble()
        f.Assemble()
        inverse = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
        apply_fe_inverse(a.mat, inverse, f.vec, solution.vec, fes.FreeDofs())
        first = solution.vec.FV().NumPy().copy()
        for _ in range(3):
            assert apply_fe_inverse(a.mat, inverse, f.vec, solution.vec, fes.FreeDofs()) < 1e-8
    np.testing.assert_allclose(solution.vec.FV().NumPy(), first, rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize("frequency", [0.0, 1000.0])
@pytest.mark.parametrize("order", [1, 2])
def test_direct_residual_scaling_and_constant_bh(monkeypatch, frequency, order):
    args = SimpleNamespace(
        R_wp=.005, H_wp=.01, R_coil=.02, R_outer=.04,
        maxh_wp=.002, maxh_air=.008, mu_r=10., sigma=2e6,
        frequency=frequency, current=100., order=order, output=None,
    )
    mesh = solver.build_mesh(args.R_wp, args.H_wp, args.R_coil,
                             args.R_outer, args.maxh_wp, args.maxh_air)
    monkeypatch.setattr(solver, "build_mesh", lambda **kwargs: mesh)
    linear = solver.run_axisym_linear(args)
    mu = solver.MU0 * args.mu_r
    nonlinear = solver.run_axisym_nonlinear(
        args, [(0., 0.), (1., mu), (1e8, mu * 1e8)])
    args.current *= 2
    doubled = solver.run_axisym_linear(args)
    for result in (linear, nonlinear, doubled):
        assert np.isfinite(result["linear_relative_residual"])
        assert result["linear_relative_residual"] < 1e-8
    assert nonlinear["picard_convergence"][-1]["max_dmu_r"] < 1e-2
    assert nonlinear["H_t_mean_A_per_m"] == pytest.approx(
        linear["H_t_mean_A_per_m"], rel=1e-8)
    assert doubled["H_t_mean_A_per_m"] == pytest.approx(
        2 * linear["H_t_mean_A_per_m"], rel=1e-8)
    if frequency:
        assert linear["P_wp_W"] > 0
        assert nonlinear["P_wp_W"] == pytest.approx(linear["P_wp_W"], rel=1e-8)
        assert doubled["P_wp_W"] == pytest.approx(4 * linear["P_wp_W"], rel=1e-8)
    else:
        assert linear["P_wp_W"] == nonlinear["P_wp_W"] == doubled["P_wp_W"] == 0


@pytest.mark.parametrize("current", [1000., 5000.])
def test_picard_converges_saturating_case_without_relaxing_acceptance(current):
    args = SimpleNamespace(
        R_wp=.005, H_wp=.01, R_coil=.02, R_outer=.04,
        maxh_wp=.0005, maxh_air=.004, mu_r=1000., sigma=2e6,
        frequency=1000., current=current, order=2, output=None,
    )
    bh = [(0., 0.), (10., .012), (100., .12),
          (1000., .9), (10000., 1.5), (1e6, 2.75)]
    result = solver.run_axisym_nonlinear(args, bh)
    assert result["picard_convergence"][-1]["constitutive_relative_residual"] < 1e-2
    assert result["linear_relative_residual"] < 1e-8
    assert result["P_wp_W"] > 0
    assert result["picard_convergence"][-1]["mu_r_wp_mean"] < args.mu_r
