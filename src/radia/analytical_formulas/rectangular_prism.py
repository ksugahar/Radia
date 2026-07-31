"""Magnetometric demagnetization factors of a rectangular prism.

The factors are the closed-form volume averages for a uniformly magnetized
rectangular prism from Aharoni, J. Appl. Phys. 83, 3432 (1998).  They are a
geometry-only reference for low-susceptibility solver validation; a finite-
susceptibility prism does not magnetize uniformly.
"""

from __future__ import annotations

import math


MU_0 = 4.0e-7 * math.pi


def _factor_along_c(a: float, b: float, c: float) -> float:
    abc = math.sqrt(a * a + b * b + c * c)
    ab = math.sqrt(a * a + b * b)
    ac = math.sqrt(a * a + c * c)
    bc = math.sqrt(b * b + c * c)
    value = (
        (b * b - c * c) / (2.0 * b * c)
        * math.log((abc - a) / (abc + a))
        + (a * a - c * c) / (2.0 * a * c)
        * math.log((abc - b) / (abc + b))
        + b / (2.0 * c) * math.log((ab + a) / (ab - a))
        + a / (2.0 * c) * math.log((ab + b) / (ab - b))
        + c / (2.0 * a) * math.log((bc - b) / (bc + b))
        + c / (2.0 * b) * math.log((ac - a) / (ac + a))
        + 2.0 * math.atan2(a * b, c * abc)
        + (a**3 + b**3 - 2.0 * c**3) / (3.0 * a * b * c)
        + (a * a + b * b - 2.0 * c * c) * abc / (3.0 * a * b * c)
        + c * (ac + bc) / (a * b)
        - (ab**3 + bc**3 + ac**3) / (3.0 * a * b * c)
    )
    return value / math.pi


def rectangular_prism_demag_factors(
    a: float,
    b: float,
    c: float,
) -> tuple[float, float, float]:
    """Return magnetometric ``(Nx, Ny, Nz)`` for side lengths ``a,b,c``.

    Only side-length ratios matter.  The factors are dimensionless and sum to
    one to roundoff.
    """
    sides = tuple(float(value) for value in (a, b, c))
    if not all(math.isfinite(value) and value > 0.0 for value in sides):
        raise ValueError("rectangular-prism side lengths must be finite and > 0")
    a, b, c = sides
    return (
        _factor_along_c(b, c, a),
        _factor_along_c(c, a, b),
        _factor_along_c(a, b, c),
    )


def linear_prism_average_flux_density(
    H0: float,
    mu_r: float,
    demag_factor: float,
) -> float:
    """Low-susceptibility magnetometric estimate of average internal ``B``.

    The approximation uses ``H = H0 / (1 + (mu_r - 1) N)`` and
    ``B = mu0 mu_r H``.  It is intended as a convergence reference as
    ``mu_r -> 1``; it is not an exact finite-susceptibility prism solution.
    """
    H0 = float(H0)
    mu_r = float(mu_r)
    demag_factor = float(demag_factor)
    if not math.isfinite(H0):
        raise ValueError("H0 must be finite")
    if not math.isfinite(mu_r) or mu_r <= 0.0:
        raise ValueError("mu_r must be finite and > 0")
    if not math.isfinite(demag_factor) or not 0.0 <= demag_factor <= 1.0:
        raise ValueError("demag_factor must be finite and in [0, 1]")
    return MU_0 * mu_r * H0 / (1.0 + (mu_r - 1.0) * demag_factor)
