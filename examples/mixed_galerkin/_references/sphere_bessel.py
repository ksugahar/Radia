"""
Reference Y_exact for a conducting sphere of radius a.

Used as the "ground truth" against which mixed Galerkin approximations
are measured.  The spherical-Bessel solution involves coth, which is
easy to evaluate at any |gamma a| without overflow concerns (unlike
the modified Bessel functions of the cylinder case).

## Formula

For an axially-uniform applied field driving azimuthal eddy currents
in a sphere of radius a, conductivity sigma, permeability mu:

    gamma  = sqrt(s mu sigma)
    gamma_a = a gamma

    Y(s) / Y_DC = (3 / gamma_a) * (coth(gamma_a) - 1 / gamma_a)

with Y_DC = (4/3) pi a^3 sigma.

## Senior tower of the sphere

The asymptotic expansion of Y/Y_DC in powers of 1/gamma_a is

    Y/Y_DC ~ 3/gamma_a - 3/(gamma_a)^2 + O(exp(-2 gamma_a))

i.e. the Senior tower TERMINATES at gamma_1 = -1/a (Senior 1962;
constructive derivation: project_senior_hoibc_tower_constructive_derivation
in lab memory).  All higher gamma_k vanish identically and the
remainder is exponentially small.  This is why the sphere mixed
Galerkin with the gamma_1 curvature correction achieves machine
precision at deep skin (Phase 7: 0.001% wall band, 0.0% at 10^8 Hz).

## API

    Y_exact_sphere(s, a, sigma, mu)   -> complex Y(s)
    Y_DC_sphere(a, sigma)             -> real Y_DC
    K_SIBC_sphere(a, sigma, mu)       -> real K_SIBC

## Self-test

Run as `python sphere_bessel.py` to verify both DC and deep-skin
limits and the smoothness of the transition.
"""

from __future__ import annotations

import cmath
import math


def Y_DC_sphere(a: float, sigma: float) -> float:
    """DC admittance Y_DC = (4/3) pi a^3 sigma (S*m^2)."""
    return (4.0 / 3.0) * math.pi * a**3 * sigma


def K_SIBC_sphere(a: float, sigma: float, mu: float) -> float:
    """Leading SIBC coefficient K_SIBC = 4 pi a^2 sqrt(sigma/mu).

    Y(s) -> K_SIBC / sqrt(s) as |gamma a| -> infinity.
    """
    return 4.0 * math.pi * a**2 * math.sqrt(sigma / mu)


def _coth_safe(z: complex) -> complex:
    """Numerically stable coth(z) for complex z.

    For |Re(z)| > 30, exp(2z) overflows in float64.  In that regime
    coth(z) -> sign(Re(z)) exponentially fast, with error
    O(exp(-2 |Re(z)|)) which is well below float64 precision.
    """
    if z.real > 30.0:
        return complex(1.0, 0.0)
    if z.real < -30.0:
        return complex(-1.0, 0.0)
    # Stable form: coth(z) = 1 + 2 / (exp(2z) - 1)
    return 1.0 + 2.0 / (cmath.exp(2.0 * z) - 1.0)


def Y_exact_sphere(s, a: float, sigma: float, mu: float):
    """Exact Y(s) for a solid sphere.

    Uses the closed-form spherical-Bessel expression with a guarded
    coth evaluation.

    Parameters
    ----------
    s : complex
        Laplace variable.
    a : float
        Sphere radius (m).
    sigma : float
        Conductivity (S/m).
    mu : float
        Permeability (H/m).

    Returns
    -------
    complex
        Y(s) in (S*m^2).
    """
    Y_DC = Y_DC_sphere(a, sigma)
    gamma_a = cmath.sqrt(s * mu * sigma) * a

    if abs(gamma_a) < 1e-6:
        # Small-argument Taylor expansion of (3/x)(coth x - 1/x).
        return Y_DC * (1.0 - gamma_a**2 / 15.0)

    coth_ga = _coth_safe(gamma_a)
    return Y_DC * (3.0 / gamma_a) * (coth_ga - 1.0 / gamma_a)


def _self_test() -> None:
    sigma = 5.8e7
    mu = 4 * math.pi * 1e-7
    a = 5e-3
    MS = mu * sigma

    print("=== sphere_bessel.py self-test ===")
    print(f"a = {a*1e3} mm, sigma = {sigma:.2e} S/m, mu = {mu:.4e} H/m")
    print(f"Y_DC   = {Y_DC_sphere(a, sigma):.4e} S*m^2")
    print(f"K_SIBC = {K_SIBC_sphere(a, sigma, mu):.4e}")
    print()

    print("DC limit (gamma_a -> 0):  Y/Y_DC should -> 1")
    for f in [1e-3, 1e-2, 1e-1, 1.0]:
        s = 1j * 2 * math.pi * f
        Y = Y_exact_sphere(s, a, sigma, mu)
        gamma_a = cmath.sqrt(s * MS) * a
        print(f"  f = {f:7.0e} Hz, |gamma a| = {abs(gamma_a):.3e}, Y/Y_DC = {abs(Y)/Y_DC_sphere(a, sigma):.6f}")
    print()

    print("Deep-skin limit: Y should -> K_SIBC / sqrt(s) at high f")
    K_SIBC = K_SIBC_sphere(a, sigma, mu)
    for f in [1e6, 1e7, 1e8, 1e9, 1e10]:
        s = 1j * 2 * math.pi * f
        Y = Y_exact_sphere(s, a, sigma, mu)
        Y_asympt = K_SIBC / cmath.sqrt(s)
        gamma_a = cmath.sqrt(s * MS) * a
        rel = abs(Y) / abs(Y_asympt)
        print(f"  f = {f:7.0e} Hz, |gamma a| = {abs(gamma_a):8.2f}, |Y| = {abs(Y):.4e}, "
              f"|Y|/|K_SIBC/sqrt(s)| = {rel:.6f}")
    print()

    print("Round-trip (no overflow at large gamma_a):")
    for f in [1e3, 1e6, 1e9, 1e12]:
        s = 1j * 2 * math.pi * f
        Y = Y_exact_sphere(s, a, sigma, mu)
        gamma_a = cmath.sqrt(s * MS) * a
        print(f"  f = {f:7.0e} Hz, |gamma a| = {abs(gamma_a):10.2f}, |Y| = {abs(Y):.6e}")


if __name__ == "__main__":
    _self_test()
