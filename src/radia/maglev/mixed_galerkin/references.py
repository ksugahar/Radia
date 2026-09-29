"""Analytic reference admittances for Modal-bulk SIBC studies.

These helpers are canonical enough to be imported by docs and lightweight
examples, but they remain analytic references rather than solver machinery.
"""

from __future__ import annotations

import cmath
import math

from scipy.special import iv

__all__ = [
    "Y_DC_cylinder",
    "K_SIBC_cylinder",
    "Y_exact_cylinder",
    "Y_DC_sphere",
    "K_SIBC_sphere",
    "Y_exact_sphere",
]

_CYLINDER_GA_CROSSOVER = 500.0
# I1/I0 = 1 + sum(b_n/z**n). From the modified-Bessel equation and
# I0' = I1 (NIST DLMF 10.25.1, 10.29.3), R' = 1 - R/z - R**2.
# Coefficient matching gives 2*b_n = (n-2)*b_(n-1)
# - sum(b_i*b_(n-i), i=1..n-1), b_0=1. These are Bessel-ratio
# coefficients, not independently attributed to an unspecified Senior paper.
_CYLINDER_BESSEL_RATIO_ASYMPT_COEFFS = (
    -1.0 / 2.0,
    -1.0 / 8.0,
    -1.0 / 8.0,
    -25.0 / 128.0,
    -13.0 / 32.0,
    -1073.0 / 1024.0,
)


def Y_DC_cylinder(a: float, sigma: float) -> float:
    """DC conductance per unit length of a cylinder."""
    return math.pi * a**2 * sigma


def K_SIBC_cylinder(a: float, sigma: float, mu: float) -> float:
    """Leading SIBC coefficient for an infinite cylinder."""
    return 2.0 * math.pi * a * math.sqrt(sigma / mu)


def Y_exact_cylinder(s: complex, a: float, sigma: float, mu: float) -> complex:
    """Cylinder admittance with a six-term Bessel-ratio asymptotic continuation.

    Above the crossover this is an asymptotic approximation, not exact Bessel
    evaluation. See docs/design/notebook_evidence_closeout_20260914.md.
    """
    y_dc = Y_DC_cylinder(a, sigma)
    gamma_a = cmath.sqrt(s * mu * sigma) * a
    if abs(gamma_a) < _CYLINDER_GA_CROSSOVER:
        return y_dc * 2.0 * iv(1, gamma_a) / (gamma_a * iv(0, gamma_a))

    ratio = 1.0 + 0j
    z_pow = 1.0 + 0j
    for coeff in _CYLINDER_BESSEL_RATIO_ASYMPT_COEFFS:
        z_pow *= gamma_a
        ratio += coeff / z_pow
    return y_dc * 2.0 * ratio / gamma_a


def Y_DC_sphere(a: float, sigma: float) -> float:
    """DC admittance of a conducting sphere."""
    return (4.0 / 3.0) * math.pi * a**3 * sigma


def K_SIBC_sphere(a: float, sigma: float, mu: float) -> float:
    """Leading SIBC coefficient for a conducting sphere."""
    return 4.0 * math.pi * a**2 * math.sqrt(sigma / mu)


def _coth_safe(z: complex) -> complex:
    if z.real > 30.0:
        return complex(1.0, 0.0)
    if z.real < -30.0:
        return complex(-1.0, 0.0)
    return 1.0 + 2.0 / (cmath.exp(2.0 * z) - 1.0)


def Y_exact_sphere(s: complex, a: float, sigma: float, mu: float) -> complex:
    """Exact sphere admittance using a guarded coth expression."""
    y_dc = Y_DC_sphere(a, sigma)
    gamma_a = cmath.sqrt(s * mu * sigma) * a
    if abs(gamma_a) < 1e-6:
        return y_dc * (1.0 - gamma_a**2 / 15.0)
    return y_dc * (3.0 / gamma_a) * (_coth_safe(gamma_a) - 1.0 / gamma_a)
