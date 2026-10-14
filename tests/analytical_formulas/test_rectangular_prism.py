"""Tests for rectangular-prism magnetometric demagnetization factors."""

import math

import pytest

from radia.analytical_formulas import (
    linear_prism_average_flux_density,
    rectangular_prism_demag_factors,
)


def test_cube_factors_are_one_third() -> None:
    factors = rectangular_prism_demag_factors(1.0, 1.0, 1.0)
    assert factors == pytest.approx((1.0 / 3.0,) * 3, abs=2.0e-15)


def test_factors_are_scale_invariant_and_sum_to_one() -> None:
    expected = rectangular_prism_demag_factors(0.011, 0.009, 0.027)
    scaled = rectangular_prism_demag_factors(11.0, 9.0, 27.0)
    assert scaled == pytest.approx(expected, abs=3.0e-15)
    assert sum(expected) == pytest.approx(1.0, abs=3.0e-15)
    assert expected == pytest.approx(
        (0.38365893740576007, 0.46352762906500533, 0.15281343352923685),
        abs=3.0e-15,
    )


def test_low_susceptibility_flux_reference() -> None:
    nx = rectangular_prism_demag_factors(0.011, 0.009, 0.027)[0]
    value = linear_prism_average_flux_density(125000.0, 1.04, nx)
    assert value == pytest.approx(0.1608936859629949, rel=2.0e-15)


@pytest.mark.parametrize("sides", [(0.0, 1.0, 1.0), (-1.0, 1.0, 1.0)])
def test_invalid_side_lengths_fail_closed(sides) -> None:
    with pytest.raises(ValueError, match="side lengths"):
        rectangular_prism_demag_factors(*sides)


@pytest.mark.parametrize("demag", [-0.1, 1.1])
def test_invalid_demag_factor_fails_closed(demag) -> None:
    with pytest.raises(ValueError, match="demag_factor"):
        linear_prism_average_flux_density(1.0, 2.0, demag)
