"""Fast regression tests for public mixed-Galerkin reference functions."""

from __future__ import annotations

import importlib.util
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




def test_cylinder_asymptotic_coefficients_follow_bessel_riccati_equation():
    from fractions import Fraction

    coefficients = [Fraction(1)]
    for n in range(1, 7):
        convolution = sum(coefficients[i] * coefficients[n-i] for i in range(1, n))
        coefficients.append(((n - 2) * coefficients[n-1] - convolution) / 2)
    assert tuple(map(float, coefficients[1:])) == references._CYLINDER_BESSEL_RATIO_ASYMPT_COEFFS


@pytest.mark.parametrize("magnitude", [499.9, 500.0, 500.1, 1000.0, 10000.0])
@pytest.mark.parametrize("phase", [-math.pi / 4, 0.0, math.pi / 4])
def test_cylinder_continuation_against_independent_scaled_bessel(magnitude, phase):
    import cmath
    from scipy.special import ive

    z = cmath.rect(magnitude, phase)
    expected = 2 * math.pi * ive(1, z) / (z * ive(0, z))
    actual = references.Y_exact_cylinder(z*z, 1.0, 1.0, 1.0)
    assert actual == pytest.approx(expected, rel=3e-13)
