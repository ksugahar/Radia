"""Levenberg-Marquardt nonlinear least squares (lsqnonlin) -- the nonlinear complement of
the linear inverse solvers (TSVD/Tikhonov). Gated on problems with KNOWN optima:
exact parameter recovery from noiseless models, a radia-relevant coil-field fit, the
Rosenbrock least-squares minimum, and agreement with scipy.optimize.least_squares.
"""
import math
import os
import sys

import pytest

np = pytest.importorskip("numpy")

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.topology_optimization.nonlinear_lsq import levenberg_marquardt

MU0 = 4e-7 * math.pi


def test_legacy_import_is_canonical_helper():
    from radia_mcp.optimization.nonlinear_lsq import levenberg_marquardt as canonical
    assert levenberg_marquardt is canonical


def test_no_improvement_is_not_convergence():
    # Deliberately inconsistent Jacobian: every trial leaves the cost unchanged.
    out = levenberg_marquardt(lambda x: np.array([1.0]), [0.0], jac=lambda x: [[1.0]])
    assert not out["converged"]
    assert out["termination_reason"] == "no_improvement"
    assert out["grad_norm"] == 1.0


@pytest.mark.parametrize("options,reason", [
    ({"max_iter": 1}, "iteration_limit"), ({"xtol": 10.0}, "step_tolerance"),
])
def test_final_gradient_and_small_step_are_not_false_success(options, reason):
    out = levenberg_marquardt(lambda x: x - 1.0, [0.0], jac=lambda x: [[1.0]], **options)
    assert not out["converged"]
    assert out["termination_reason"] == reason
    assert out["grad_norm"] == pytest.approx(abs(out["x"][0]-1.0))
    assert out["cost"] == pytest.approx(.5*(out["x"][0]-1.0)**2)


@pytest.mark.parametrize("options", [{"max_iter": 0}, {"max_iter": 1.5},
    {"lam0": -1}, {"eps": 0}, {"gtol": np.inf}, {"xtol": np.nan}])
def test_invalid_controls_fail_loudly(options):
    with pytest.raises(ValueError):
        levenberg_marquardt(lambda x: x, [1.0], **options)


def test_invalid_callback_evidence_fails_loudly():
    with pytest.raises(ValueError, match="residual"):
        levenberg_marquardt(lambda x: [np.nan], [0.0])
    with pytest.raises(ValueError, match="Jacobian"):
        levenberg_marquardt(lambda x: x, [1.0], jac=lambda x: [[np.inf]])
    with pytest.raises(ValueError, match="Jacobian"):
        levenberg_marquardt(lambda x: x, [1.0], jac=lambda x: [[1,2]])


def test_exact_parameter_recovery_exponential():
    # y = a exp(b t) + c, noiseless -> LM recovers (a,b,c) to ~machine precision
    t = np.linspace(0, 2, 25)
    true = np.array([2.5, -1.3, 0.7])
    y = true[0] * np.exp(true[1] * t) + true[2]
    res = levenberg_marquardt(lambda p: p[0] * np.exp(p[1] * t) + p[2] - y, [1.0, -0.5, 0.0])
    assert res["converged"]
    assert np.max(np.abs(res["x"] - true)) < 1e-10
    assert res["cost"] < 1e-20


def test_coil_axial_field_geometry_recovery():
    # radia-relevant inverse fit: circular-loop on-axis field B_z(z)=mu0 I a^2/(2(a^2+z^2)^{3/2})
    # -> recover (current I, radius a) from noiseless samples. Residual NORMALISED to O(1)
    # (the scaling caveat) so recovery reaches machine precision.
    zc = np.linspace(-0.5, 0.5, 21)
    I_true, a_true = 13.0, 0.25
    Bz = MU0 * I_true * a_true**2 / (2.0 * (a_true**2 + zc**2)**1.5)
    scale = Bz.max()

    def res(p):
        I, a = p
        return (MU0 * I * a**2 / (2.0 * (a**2 + zc**2)**1.5) - Bz) / scale

    out = levenberg_marquardt(res, [5.0, 0.1])
    assert out["converged"]
    assert abs(out["x"][0] - I_true) / I_true < 1e-8
    assert abs(out["x"][1] - a_true) / a_true < 1e-8


def test_rosenbrock_least_squares_minimum():
    # the Rosenbrock valley posed as residuals r=[1-x, 10(y-x^2)] -> global min (1,1), cost 0
    def res(p):
        x, y = p
        return np.array([1.0 - x, 10.0 * (y - x * x)])
    out = levenberg_marquardt(res, [-1.2, 1.0])
    assert np.max(np.abs(out["x"] - np.array([1.0, 1.0]))) < 1e-8
    assert out["cost"] < 1e-18


def test_agreement_with_scipy_least_squares():
    sp_opt = pytest.importorskip("scipy.optimize")
    t = np.linspace(0, 2, 25)
    true = np.array([2.5, -1.3, 0.7])
    y = true[0] * np.exp(true[1] * t) + true[2]
    f = lambda p: p[0] * np.exp(p[1] * t) + p[2] - y
    mine = levenberg_marquardt(f, [1.0, -0.5, 0.0])["x"]
    theirs = sp_opt.least_squares(f, [1.0, -0.5, 0.0], method="lm").x
    assert np.max(np.abs(mine - theirs)) < 1e-9


def test_user_supplied_jacobian_matches_numerical():
    # an analytic Jacobian gives the same optimum as forward differences (a linear LS here)
    A = np.array([[1.0, 2.0], [3.0, 1.0], [2.0, 2.0]])
    b = np.array([1.0, 2.0, 1.5])
    res = lambda p: A @ p - b
    jac = lambda p: A
    x_num = levenberg_marquardt(res, [0.0, 0.0])["x"]
    x_jac = levenberg_marquardt(res, [0.0, 0.0], jac=jac)["x"]
    x_exact = np.linalg.lstsq(A, b, rcond=None)[0]      # the normal-equation solution
    assert np.max(np.abs(x_num - x_exact)) < 1e-9
    assert np.max(np.abs(x_jac - x_exact)) < 1e-9
