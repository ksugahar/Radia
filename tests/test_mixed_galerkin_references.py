"""Fast regression tests for public mixed-Galerkin reference functions."""

from __future__ import annotations

import math

import pytest

from radia.maglev.mixed_galerkin.references import Y_DC_cylinder, Y_cln_pade


MU0 = 4.0 * math.pi * 1e-7
SIGMA_CU = 5.8e7


@pytest.mark.parametrize("kind", ["L", "R", "l", "r"])
def test_cln_pade_has_exact_dc_anchor(kind: str) -> None:
    radius = 5e-3
    value = Y_cln_pade(0.0j, 4, radius, SIGMA_CU, MU0, kind=kind)
    assert value == pytest.approx(Y_DC_cylinder(radius, SIGMA_CU), rel=2e-14)


@pytest.mark.parametrize("kind", ["L", "R"])
def test_cln_pade_preserves_dimensionless_radius_scaling(kind: str) -> None:
    radius = 5e-3
    scale = 1.7
    s = 2j * math.pi * 12_500.0
    base = Y_cln_pade(s, 4, radius, SIGMA_CU, MU0, kind=kind)
    scaled = Y_cln_pade(
        s / scale**2, 4, radius * scale, SIGMA_CU, MU0, kind=kind)
    assert scaled / scale**2 == pytest.approx(base, rel=2e-12)


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"N": 0}, "N must be at least 1"),
        ({"N": 1.5}, "must be integers"),
        ({"n_modes": 0}, "n_modes must be at least 1"),
        ({"n_modes": 2.5}, "must be integers"),
        ({"n_taylor": -1}, "n_taylor must be non-negative"),
        ({"n_taylor": 4.5}, "must be integers"),
        ({"a": 0.0}, "must be positive"),
        ({"sigma": 0.0}, "must be positive"),
        ({"mu": 0.0}, "must be positive"),
        ({"kind": "X"}, "kind must be 'L' or 'R'"),
        ({"N": 8, "kind": "R", "n_taylor": 15}, "needs 16 Taylor terms"),
    ],
)
def test_cln_pade_rejects_invalid_contract(kwargs: dict, match: str) -> None:
    args = {"s": 1j, "N": 4, "a": 5e-3, "sigma": SIGMA_CU, "mu": MU0}
    args.update(kwargs)
    with pytest.raises(ValueError, match=match):
        Y_cln_pade(**args)
