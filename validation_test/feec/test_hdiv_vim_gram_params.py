"""Golden: RT1 charge-Gram default resolution.

HDiv-VIM is now TET / RT1 only.  The public solver rejects RT0, RT2+, the
Gauss point backend, and the H-LU system-A path before assembly.  This pure
helper is therefore locked only for the live RT1 path: ACA eps, far quadrature,
and the high-order near/far threshold.  ``near_factor`` remains an RT0-era knob
and must fail loud if someone tries to pass it into RT1.
"""
from radia.vim._solve import _resolve_gram_params as R   # noqa: E402  (pure logic, no ngsolve needed)
import pytest                                              # noqa: E402

_BASE = dict(gram_backend="analytic", linear_solver="auto", gram_eps=None, far_quad=None)


def test_rt1_defaults():
    assert R(order=1, uniform_linear=False, near_factor=None, ho_far_factor=None, **_BASE) == \
        {"eps": 1e-10, "far_quad": 3, "ho_far_factor": 2.0}


def test_rt1_explicit_values_win():
    g = R(order=1, uniform_linear=False, near_factor=None, ho_far_factor=float("inf"),
          gram_backend="analytic", linear_solver="auto", gram_eps=1e-9, far_quad=5)
    assert g["eps"] == 1e-9
    assert g["far_quad"] == 5
    assert g["ho_far_factor"] == float("inf")


def test_rt1_rejects_rt0_near_factor():
    with pytest.raises(ValueError, match="near_factor is an order=0"):
        R(order=1, uniform_linear=False, near_factor=2.0, ho_far_factor=None, **_BASE)
