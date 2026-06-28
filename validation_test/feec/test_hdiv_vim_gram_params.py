"""Golden: the Gram-build default resolution (radia.vim._solve._resolve_gram_params) is the SINGLE source of
truth for the charge-Gram ACA-eps + near/far-split defaults (Gram-path consolidation step A), and the near/far
knob is ORDER-SPECIFIC (step C): near_factor is RT0 (order=0)-only, ho_far_factor is high-order (order>0)-only,
and a wrong-order knob FAILS LOUD (No-Fallbacks) rather than being silently remapped (the prior overload, where
high-order silently fed near_factor into ho_far_factor, is removed).

Locks:
 (A) the resolved defaults reproduce the prior inline logic exactly -- RT0: eps 1e-12 (uniform-linear) / 1e-10,
     near 2.0 (fast) / 1e30 (H-LU exact), far 4 / 0; high-order: eps 1e-10, far 3, ho_far 2.0;
 (C) RT0 rejects ho_far_factor; high-order rejects near_factor; an explicit value always wins.
"""
from radia.vim._solve import _resolve_gram_params as R   # noqa: E402  (pure logic, no ngsolve needed)
import pytest                                              # noqa: E402

_BASE = dict(gram_backend="analytic", linear_solver="auto", gram_eps=None, far_quad=None)


def test_rt0_defaults_match_prior_inline_logic():
    # uniform-linear auto path: tight eps + fast near/far split
    assert R(order=0, uniform_linear=True, near_factor=None, ho_far_factor=None, **_BASE) == \
        {"eps": 1e-12, "near_factor": 2.0, "far_quad": 4}
    # non-uniform (per-region / nonlinear) keeps 1e-10, still fast build
    assert R(order=0, uniform_linear=False, near_factor=None, ho_far_factor=None, **_BASE) == \
        {"eps": 1e-10, "near_factor": 2.0, "far_quad": 4}


def test_rt0_hlu_keeps_exact_all_analytic():
    g = R(order=0, gram_backend="analytic", linear_solver="hlu", uniform_linear=True,
          gram_eps=None, near_factor=None, far_quad=None, ho_far_factor=None)
    assert g == {"eps": 1e-10, "near_factor": 1e30, "far_quad": 0}


def test_rt0_gauss_backend_not_fast_build():
    # gram_backend != analytic -> fast_build False -> exact near/far (the gauss path ignores these anyway)
    g = R(order=0, gram_backend="gauss", linear_solver="auto", uniform_linear=True,
          gram_eps=None, near_factor=None, far_quad=None, ho_far_factor=None)
    assert g == {"eps": 1e-12, "near_factor": 1e30, "far_quad": 0}


def test_rt0_explicit_near_factor_wins():
    assert R(order=0, uniform_linear=True, near_factor=1e30, ho_far_factor=None, **_BASE)["near_factor"] == 1e30


def test_highorder_defaults():
    assert R(order=2, uniform_linear=False, near_factor=None, ho_far_factor=None, **_BASE) == \
        {"eps": 1e-10, "far_quad": 3, "ho_far_factor": 2.0}


def test_highorder_explicit_ho_far_factor_wins():
    g = R(order=2, uniform_linear=False, near_factor=None, ho_far_factor=float("inf"), **_BASE)
    assert g["ho_far_factor"] == float("inf")


def test_rt0_rejects_ho_far_factor():
    with pytest.raises(ValueError, match="ho_far_factor is an order>0"):
        R(order=0, uniform_linear=True, near_factor=None, ho_far_factor=2.0, **_BASE)


def test_highorder_rejects_near_factor():
    with pytest.raises(ValueError, match="near_factor is an order=0"):
        R(order=2, uniform_linear=False, near_factor=2.0, ho_far_factor=None, **_BASE)
