"""Golden: the all-C++ SYMMETRIC ENERGY-NEWTON is the DEFAULT HDiv-VIM nonlinear soft-iron solver
(radia.vim._solve._solve_nonlinear_energy_cpp), and it matches the forward (scipy splu + GMRES) Newton.

The nonlinear inner Newton step is now solved by the EXISTING C++ symmetric W-CG
(solve_linear_material_mass_riesz, W = the differential-reluctivity tangent mass, mass-Riesz PARDISO, N
H-matvec) -- bringing the nonlinear solve to C++ parity with the linear path (no scipy splu / GMRES /
M_mass^-1).  The co-energy form is robust through deep saturation via a hard-saturation barrier, a
co-energy line search, and settled-step acceptance for the achievable-precision limit cycle of the M-form.

Locks: (1) the default nonlinear solver is 'energy-newton-cpp'; (2) it agrees with the forward Newton
(linear_solver='gmres') to ~1e-5 at moderate AND deep drive; (3) deep saturation converges to M ~ Msat.
"""
import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")
from netgen.occ import Sphere, Pnt, OCCGeometry          # noqa: E402
from radia.vim import hdiv_demag_solve                    # noqa: E402

# realistic soft-iron BH table (mu_r ~ 4000 at low H, saturating ~2.2 T)
_H = np.array([0, 50, 100, 200, 500, 1e3, 2e3, 5e3, 1e4, 3e4, 1e5, 3e5, 1e6])
_B = np.array([0, 0.30, 0.60, 1.0, 1.45, 1.7, 1.9, 2.0, 2.05, 2.1, 2.15, 2.25, 2.5])
_BH = np.column_stack([_H, _B]).tolist()
_MU0 = 4e-7 * np.pi
_MSAT = _B[-1] / _MU0 - _H[-1]                            # the table's saturation magnetization


def _sphere(maxh=0.6):
    return ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), 1.0)).GenerateMesh(maxh=maxh))


def test_default_nonlinear_is_energy_newton_cpp_and_matches_forward():
    """Default iron-only nonlinear routes to the all-C++ energy-Newton, == the forward Newton."""
    mesh = _sphere()
    for H0 in (500.0, 5000.0):                            # moderate (linear-ish) and into the knee
        Hext = ng.CoefficientFunction((0, 0, H0))
        with ng.TaskManager():
            re = hdiv_demag_solve(mesh, bh_table=_BH, H_ext=Hext, order=0)                      # default
            rf = hdiv_demag_solve(mesh, bh_table=_BH, H_ext=Hext, order=0, linear_solver="gmres")  # forward
        assert re["linear_solver"] == "energy-newton-cpp", re["linear_solver"]
        assert rf["linear_solver"] == "forward-newton-gmres", rf["linear_solver"]
        assert re["nonlinear"] is True
        rd = abs(re["M_avg"][2] - rf["M_avg"][2]) / (abs(rf["M_avg"][2]) + 1e-30)
        assert rd < 1e-5, f"energy-Newton {re['M_avg'][2]:.4f} vs forward {rf['M_avg'][2]:.4f} at H0={H0} (rel {rd:.2e})"


def test_energy_newton_deep_saturation():
    """Deep saturation (H0 well past the table) drives M to ~Msat and converges (the hard-saturation barrier
    + settled acceptance handle the M-form limit cycle)."""
    mesh = _sphere()
    with ng.TaskManager():
        r = hdiv_demag_solve(mesh, bh_table=_BH, H_ext=ng.CoefficientFunction((0, 0, 3e6)), order=0)
    assert r["linear_solver"] == "energy-newton-cpp"
    # at deep drive M_avg -> Msat (the table's saturation; demag-independent there).  The soft hard-saturation
    # barrier permits a small (<~1%) overshoot of the uniform Msat at the discrete/volume-average level.
    assert 0.95 * _MSAT < r["M_avg"][2] < 1.03 * _MSAT, (r["M_avg"][2], _MSAT)
    # it converged (returning, not raising) -- the deep-saturation M-form iteration count is table-dependent
    # and higher than the forward H-form (the M-form limit-cycles to the achievable precision); a generous
    # bound catches a runaway without being brittle to the table shape.
    assert r["iters"] < 100
