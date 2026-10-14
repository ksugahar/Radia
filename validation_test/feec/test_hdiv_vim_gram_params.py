"""Golden: HDiv order-1 charge-Gram default resolution."""
from radia.vim._solve import _resolve_gram_params as R   # noqa: E402  (pure logic, no ngsolve needed)

_BASE = dict(gram_eps=None, far_quad=None, ho_far_factor=None)


def test_rt1_defaults():
    assert R(**_BASE) == \
        {"eps": 1e-10, "far_quad": 3, "ho_far_factor": 2.0}


def test_rt1_explicit_values_win():
    g = R(gram_eps=1e-9, far_quad=5, ho_far_factor=float("inf"))
    assert g["eps"] == 1e-9
    assert g["far_quad"] == 5
    assert g["ho_far_factor"] == float("inf")
