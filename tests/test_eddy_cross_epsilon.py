"""Analytic cross-kernel oracle; no full electromagnetic solve is claimed."""

import numpy as np
import pytest

from validation_test.vim_coupled.validate_eddy_cross_epsilon import run, spherical_basis


@pytest.fixture(scope="module")
def evidence():
    return run()


@pytest.mark.parametrize("kind,power,factor", [("volume", 3, 4 / 3), ("surface", 2, 4)])
def test_spherical_rule_has_exact_measure(kind, power, factor):
    basis = spherical_basis(3, kind, np.zeros(3), radius=0.5)
    assert np.sum(basis.weights) == pytest.approx(factor * np.pi * 0.5**power)


def test_native_dense_leaf_matches_independent_sampled_sum(evidence):
    records = [r for r in evidence["records"] if r["compression"] == "dense-leaf"]
    assert len(records) == 18
    assert all(r["lowrank_leaf_count"] == 0 for r in records)
    assert max(abs(r["native_minus_sampled_relative"]) for r in records) < 1e-12


def test_error_budget_closes_without_conflating_epsilon_and_aca(evidence):
    for r in evidence["records"]:
        decomposed = (r["quadrature_relative_error"] + r["epsilon_shift_relative"]
                      + r["native_minus_sampled_relative"])
        assert decomposed == pytest.approx(r["total_relative_error"], abs=1e-14)
        assert r["reciprocity_relative_error"] < 1e-12
    assert evidence["certifies_hdiv_sibc_coupled_solver"] is False


@pytest.mark.parametrize("pair", ["volume-volume", "volume-surface", "surface-surface"])
def test_epsilon_bias_decreases_against_analytic_oracle(evidence, pair):
    records = [r for r in evidence["records"] if r["pair"] == pair
               and r["quadrature_order"] == 5 and r["compression"] == "dense-leaf"]
    assert len(records) == 3
    assert max(abs(r["quadrature_relative_error"]) for r in records) < 1e-8
    shifts = [r["epsilon_shift_relative"] for r in records]
    assert shifts[0] < shifts[1] < shifts[2] < 0
    # Positive epsilon lowers 1/r; fixed quadrature cannot remove that bias.
    assert abs(shifts[-1]) > 1000 * max(abs(r["quadrature_relative_error"]) for r in records)


@pytest.mark.xfail(strict=True, reason=(
    "Known sampled volume-volume ACA error floor: native candidate "
    "966adc2e5 misses the explicit 1e-8 cross-block acceptance budget even "
    "at aca_eps=1e-10; dense-leaf evaluation passes."
))
def test_compressed_volume_cross_meets_explicit_accuracy_budget(evidence):
    record = next(r for r in evidence["records"] if r["pair"] == "volume-volume"
                  and r["quadrature_order"] == 5 and r["epsilon_m"] == 0.3
                  and r["aca_tolerance"] == 1e-10 and r["compression"] == "compressed")
    assert record["lowrank_leaf_count"] > 0
    assert abs(record["native_minus_sampled_relative"]) < 1e-8
