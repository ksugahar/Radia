"""Golden: F2 frontier -- Clebsch recovery from a general div-free J is ILL-POSED.

Locks the WALL behavior (an honest partly-negative result), NOT a fake success.
On KNOWN Clebsch currents (ground truth available) the textbook recovery recipe
(solve J.grad(mu)=0 for mu, then lambda) exhibits two unfixable obstructions:

  WALL 1 -- NON-UNIQUENESS: the transport PDE J.grad(mu)=0 is solved to ~1e-9 by
            BOTH the seed mu0=x AND the reparametrization f(mu0)=x^2; the two
            recovered stream surfaces are genuinely distinct (O(1) apart) yet
            related by mu_b = f(mu_a).  Any f(mu0) is a valid second scalar, so
            recovery returns the seed, not a canonical pair.

  WALL 2 -- NULL DEGENERACY: for a Clebsch J=(0,0,x^2+y^2) vanishing on the z-axis,
            the recovery conditioning ~ 1/|J|^2 blows up like 1/r^4 toward the
            null (log-log slope ~ -4); the lambda-solve is singular ON the null.

See examples/vim/clebsch_recovery_wall.py and streamfunction("clebsch_3d").
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__),
                                "../../examples/vim"))


@pytest.fixture(scope="module")
def result():
    import clebsch_recovery_wall as cw
    res, ok = cw.main()
    return res, ok


def test_both_seeds_solve_the_transport_pde(result):
    """Both mu0 and f(mu0) solve J.grad(mu)=0 -- both are valid stream surfaces."""
    res, _ = result
    assert res["wall1"]["res_seed_mu0"] < 1e-5, res["wall1"]["res_seed_mu0"]
    assert res["wall1"]["res_seed_fmu0"] < 1e-5, res["wall1"]["res_seed_fmu0"]


def test_recovered_mu_is_a_reparametrization(result):
    """The f(mu0)=x^2 seed recovers mu_b = f(mu_a) = mu_a^2 (chain-rule gauge)."""
    res, _ = result
    assert res["wall1"]["reparam_match"] < 1e-3, res["wall1"]["reparam_match"]


def test_recovery_is_non_unique(result):
    """...yet mu_b is genuinely distinct from mu_a (O(1) apart): NON-UNIQUE.
    This is the wall: recovery returns whatever the boundary seed dictates."""
    res, _ = result
    assert res["wall1"]["distinct"] > 0.1, res["wall1"]["distinct"]


def test_recovery_is_singular_at_field_null(result):
    """cond(lambda-recovery) ~ 1/|J|^2 ~ 1/r^4 toward the null: log-log slope ~ -4."""
    res, _ = result
    assert res["wall2"]["cond_slope"] < -3.5, res["wall2"]["cond_slope"]
    assert res["wall2"]["cond_ratio"] > 1e3, res["wall2"]["cond_ratio"]


def test_overall_wall_demonstrated(result):
    """The probe PASS asserts the WALL (ill-posedness), not a recovery success."""
    _, ok = result
    assert ok
