"""Direct axisymmetric solves: residual, source scaling and constant-BH limit."""
from types import SimpleNamespace

import numpy as np
import pytest

from radia.panels import calc_axisym_volumetric as solver


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


def test_peak_flux_density_of_a_uniform_axial_phasor_is_its_amplitude():
    """A_phi = c B0 r / 2 with |c| = 1 is the uniform field B_z = c B0.

    Phasors are peak amplitudes (loss 0.5 sigma omega^2 |A|^2), so the B-H
    update must read |B| = B0, not sqrt(2) B0.
    """
    from ngsolve import H1, GridFunction, Integrate, x

    mesh = solver.build_mesh(.005, .01, .02, .04, .002, .008)
    fes = H1(mesh, order=2, complex=True)
    gfu = GridFunction(fes)
    B0, phase = 0.7, (1 + 1j) / np.sqrt(2)
    gfu.Set(phase * B0 * x / 2)
    area = Integrate(1, mesh, definedon=mesh.Materials("workpiece"))
    mean_B = Integrate(solver.peak_flux_density_cf(gfu), mesh,
                       definedon=mesh.Materials("workpiece")) / area
    assert mean_B == pytest.approx(B0, rel=1e-8)


def test_secant_permeability_is_continuous_and_starts_at_the_table_slope():
    bh = [(0., 0.), (10., .012), (100., .12), (1000., .9), (10000., 1.5), (1e6, 2.75)]
    mu_r = solver.secant_permeability_from_table(bh)
    initial = .012 / 10. / solver.MU0
    assert mu_r(np.array([0.]))[0] == pytest.approx(initial, rel=1e-6)
    values = mu_r(np.linspace(1e-6, 2.7, 20001))
    assert np.all(np.isfinite(values)) and np.all(values >= 1.0)
    # no jump at the first table row (the old lookup switched to --mu-r there)
    below, above = mu_r(np.array([.012 - 1e-7, .012 + 1e-7]))
    assert above == pytest.approx(below, rel=1e-4)
    for H, Bv in bh[1:]:
        assert mu_r(np.array([Bv]))[0] == pytest.approx(Bv / (solver.MU0 * H), rel=1e-3)
    # beyond the table the law continues with slope mu0
    assert mu_r(np.array([3.0]))[0] == pytest.approx(
        3.0 / (solver.MU0 * (1e6 + (3.0 - 2.75) / solver.MU0)), rel=1e-9)


def test_picard_converges_saturating_case_without_relaxing_acceptance():
    args = SimpleNamespace(
        R_wp=.005, H_wp=.01, R_coil=.02, R_outer=.04,
        maxh_wp=.0005, maxh_air=.004, mu_r=1000., sigma=2e6,
        frequency=1000., current=1000., order=2, output=None,
    )
    bh = [(0., 0.), (10., .012), (100., .12),
          (1000., .9), (10000., 1.5), (1e6, 2.75)]
    result = solver.run_axisym_nonlinear(args, bh)
    assert result["picard_convergence"][-1]["constitutive_relative_residual"] < 1e-2
    assert result["linear_relative_residual"] < 1e-8
    assert result["P_wp_W"] > 0
    assert result["picard_convergence"][-1]["mu_r_wp_mean"] < args.mu_r
