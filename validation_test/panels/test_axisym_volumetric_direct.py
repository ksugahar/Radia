"""Direct axisymmetric solves: residual, source scaling and constant-BH limit."""
from types import SimpleNamespace

import numpy as np
import pytest

from radia.panels import calc_axisym_volumetric as solver


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
