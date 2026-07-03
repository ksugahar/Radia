"""Golden: Frontier F1 -- the FOLIATION-CHOICE wall (non-convex outer loop).

Locks the HONEST partly-negative result that choosing the Clebsch foliation mu for
a general (non-axisymmetric) target is a NON-CONVEX outer problem, even though the
inner lambda-solve (mu fixed) is convex.  This is the volume analogue of the plasma
"coil winding surface" choice: fixing the surface gives REGCOIL's convex current-
potential solve, but choosing the surface is non-convex (inverse Biot-Savart
|r-r'|^-3 nonlinearity; arXiv:2408.08267).

The golden asserts the WALL behaviour, NOT a solve:
  * inner solve is convex/feasible -- the field is fit EXACTLY at every foliation
    tilt (residual ~ machine eps);
  * single-axis target -> the coil-complexity outer objective g(t)=||J|| is
    UNIMODAL (one basin): foliation choice is convex when mu can align;
  * two-axis target -> g(t) is MULTI-MODAL (>= 2 local minima) with a real barrier:
    the genuinely non-convex foliation choice that has no closed form / no convex
    reformulation.

See validation_test/feec/vim_legacy/foliation_choice_wall.py and the unification memory file
memory/clebsch_cohomology_streamfunction_unification.md (Honest frontier F1).
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__),
                                "../../validation_test/feec/vim_legacy"))


@pytest.fixture(scope="module")
def result():
    import foliation_choice_wall as fw
    res, ok = fw.main()
    return res, ok


def test_inner_solve_is_convex(result):
    """With mu fixed the lambda-solve is a convex least-norm problem: the target
    field is fit EXACTLY at every foliation tilt."""
    res, _ = result
    assert res["fit_max"] < 1e-6, f"inner fit residual {res['fit_max']:.2e}"


def test_single_axis_choice_is_convex(result):
    """A single-axis (cylinder-representable) target -> the foliation-choice outer
    objective is UNIMODAL (one basin): mu aligned with the loop axis suffices."""
    res, _ = result
    assert res["n_min_single"] == 1, \
        f"single-axis g(t) has {res['n_min_single']} local minima, expected 1"


def test_two_axis_choice_is_non_convex(result):
    """The WALL: a two-axis target -> the SAME outer objective is MULTI-MODAL
    (>= 2 local minima) with a real barrier -- non-convex foliation choice."""
    res, _ = result
    assert res["n_min_two"] >= 2, \
        f"two-axis g(t) has {res['n_min_two']} local minima, expected >= 2 (the wall)"
    assert res["barrier_two"] > 0.05, \
        f"barrier between basins {res['barrier_two']:.2%} too small to be a wall"


def test_overall_pass(result):
    _, ok = result
    assert ok
