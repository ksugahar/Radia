"""Fast regression tests for public mixed-Galerkin reference functions."""

from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
# Analytic references do not need the native solver or its DLLs.
spec = importlib.util.spec_from_file_location(
    "mixed_galerkin_references", ROOT / "src/radia/maglev/mixed_galerkin/references.py"
)
if spec is None or spec.loader is None:
    raise ImportError("Cannot load the mixed-Galerkin analytic reference module")
references = importlib.util.module_from_spec(spec)
spec.loader.exec_module(references)
Y_DC_cylinder = references.Y_DC_cylinder
Y_cln_pade = references.Y_cln_pade


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
    scaled = Y_cln_pade(s / scale**2, 4, radius * scale, SIGMA_CU, MU0, kind=kind)
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
        ({"a": float("nan")}, "must be positive"),
        ({"mu": float("inf")}, "must be positive"),
        ({"kind": "X"}, "kind must be 'L' or 'R'"),
        ({"N": 8, "kind": "R", "n_taylor": 15}, "needs 16 Taylor terms"),
    ],
)
def test_cln_pade_rejects_invalid_contract(kwargs: dict, match: str) -> None:
    args = {"s": 1j, "N": 4, "a": 5e-3, "sigma": SIGMA_CU, "mu": MU0}
    args.update(kwargs)
    with pytest.raises(ValueError, match=match):
        Y_cln_pade(**args)


@pytest.mark.parametrize("kind", ["L", "R"])
def test_tenth_order_pade_matches_independent_high_precision_reference(kind):
    evidence = json.loads(
        (ROOT / "validation_test/mixed_galerkin/results/cln_pade_reference.json").read_text(
            encoding="utf-8"
        )
    )
    case = next(
        row
        for row in evidence["cases"]
        if row["N"] == 10 and row["kind"] == kind and row["u_imag"] == 1e6
    )
    expected = complex(*case["reference_ratio"])
    actual = (
        Y_cln_pade(1e6j, 10, 1.0, 1.0, 1.0, kind=kind, n_modes=evidence["n_modes"], n_taylor=20)
        / math.pi
    )
    assert actual == pytest.approx(expected, rel=1e-8, abs=1e-12)


@pytest.mark.parametrize("kind", ["L", "R"])
def test_pade_handles_exhausted_modal_rank(kind):
    from scipy.special import jn_zeros

    u = 30j
    zeros = jn_zeros(0, 2)
    expected = 1 - sum(4 / x**2 * u / (x**2 + u) for x in zeros)
    actual = Y_cln_pade(u, 4, 1.0, 1.0, 1.0, kind=kind, n_modes=2) / math.pi
    assert actual == pytest.approx(expected, rel=2e-13)
