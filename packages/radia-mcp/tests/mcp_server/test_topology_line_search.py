"""Armijo backtracking line search -- the local optimization stabilizer used by
steepest descent, nonlinear CG, and quasi-Newton methods. Gated on simple
objectives with known sufficient-decrease behaviour.
"""
import os
import sys

import pytest

np = pytest.importorskip("numpy")

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.topology_optimization.line_search import armijo_backtracking


def test_quadratic_accepts_half_step_to_minimum():
    def f(x):
        return float((x[0] - 2.0) ** 2)

    x0 = np.array([0.0])
    g0 = np.array([-4.0])
    out = armijo_backtracking(f, x0, -g0, f0=f(x0), g0=g0)
    assert out["accepted"]
    assert out["n_iter"] == 2
    assert out["alpha"] == pytest.approx(0.5)
    assert out["x_new"][0] == pytest.approx(2.0)
    assert out["f_new"] == pytest.approx(0.0)
    assert out["slope0"] == pytest.approx(-16.0)


def test_rosenbrock_steepest_descent_decreases_objective():
    def f(x):
        return float(100.0 * (x[1] - x[0] * x[0]) ** 2 + (1.0 - x[0]) ** 2)

    def grad(x):
        return np.array([
            -400.0 * x[0] * (x[1] - x[0] * x[0]) - 2.0 * (1.0 - x[0]),
            200.0 * (x[1] - x[0] * x[0]),
        ])

    x0 = np.array([-1.2, 1.0])
    g0 = grad(x0)
    out = armijo_backtracking(f, x0, -g0, grad=grad)
    assert out["accepted"]
    assert out["alpha"] < 1.0
    assert out["f_new"] < f(x0)
    assert out["slope0"] == pytest.approx(-float(g0 @ g0))


def test_rejects_non_descent_direction():
    with pytest.raises(ValueError, match="descent"):
        armijo_backtracking(lambda x: float((x[0] - 2.0) ** 2),
                            np.array([0.0]), np.array([-1.0]),
                            g0=np.array([-4.0]))


def test_requires_gradient_or_cached_gradient():
    with pytest.raises(ValueError, match="g0 or grad"):
        armijo_backtracking(lambda x: float(x[0] ** 2), np.array([1.0]),
                            np.array([-1.0]))


@pytest.mark.parametrize(
    "kwargs",
    [
        {"alpha0": 0.0},
        {"rho": 1.0},
        {"rho": 0.0},
        {"c1": 1.0},
        {"c1": 0.0},
        {"min_alpha": 0.0},
        {"max_iter": 0},
    ],
)
def test_parameter_validation(kwargs):
    with pytest.raises(ValueError):
        armijo_backtracking(lambda x: float(x[0] ** 2), np.array([1.0]),
                            np.array([-1.0]), g0=np.array([2.0]), **kwargs)
