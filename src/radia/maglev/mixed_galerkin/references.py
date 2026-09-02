"""Analytic reference admittances for mixed-Galerkin examples.

These helpers are canonical enough to be imported by docs and lightweight
examples, but they remain analytic references rather than solver machinery.
"""

from __future__ import annotations

import cmath
import math

from scipy.special import iv, jn_zeros

__all__ = [
    "Y_cln_pade",
    "Y_DC_cylinder",
    "K_SIBC_cylinder",
    "Y_exact_cylinder",
    "Y_DC_sphere",
    "K_SIBC_sphere",
    "Y_exact_sphere",
]

_CYLINDER_GA_CROSSOVER = 500.0
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
    """Exact cylinder admittance with Senior-tower overflow continuation."""
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


def Y_cln_pade(s: complex, N: int, a: float, sigma: float, mu: float, *,
               kind: str = "L", n_modes: int = 200,
               n_taylor: int = 30) -> complex:
    """Admittance of a Cauer ladder truncated to ``N`` stages, for a disk.

    The exact disk admittance is a Foster sum over the Bessel eigenmodes,

        Y(s)/Y_DC = 1 - sum_n (4 / j_n^2) * u / (j_n^2 + u),
        u = s mu sigma a^2

    with ``j_n`` the zeros of J_0 and ``lambda_n = (j_n / a)^2``. A Cauer
    ladder of ``N`` stages is the Pade approximant of that series at ``s = 0``,
    so the truncation is computed here rather than by running the ladder
    recursion: the two agree, and the Pade route does not lose digits to the
    difference-and-divide step the recursion takes.

    That is also why the truncated ladder has the wrong tail. A Pade
    approximant is rational, and a rational function has an integer asymptotic
    slope, whereas the disk decays as ``s^-1/2``. Raising ``N`` moves the poles
    around; it does not change the exponent.

    Args:
        s: Laplace variable.
        N: number of ladder stages.
        a: disk radius, metres.
        sigma: conductivity, S/m.
        mu: permeability, H/m.
        kind: ``"L"`` for the inductance-terminated ladder, Pade [N-1/N];
            ``"R"`` for the resistance-terminated one, Pade [N/N].
        n_modes: Bessel zeros used to build the moments.
        n_taylor: highest Taylor order retained before the Pade solve; must be
            at least ``2N`` for kind ``"R"``.

    Returns:
        The truncated ladder's admittance at ``s``.

    Raises:
        ValueError: on a non-physical geometry/material value, invalid order,
            unknown ``kind``, or insufficient Taylor order.
    """
    import operator

    import numpy as np

    try:
        N = operator.index(N)
        n_modes = operator.index(n_modes)
        n_taylor = operator.index(n_taylor)
    except TypeError as exc:
        raise ValueError("N, n_modes, and n_taylor must be integers") from exc

    if N < 1:
        raise ValueError(f"N must be at least 1, not {N!r}")
    if n_modes < 1:
        raise ValueError(f"n_modes must be at least 1, not {n_modes!r}")
    if n_taylor < 0:
        raise ValueError(f"n_taylor must be non-negative, not {n_taylor!r}")
    if a <= 0.0 or sigma <= 0.0 or mu <= 0.0:
        raise ValueError("a, sigma, and mu must be positive")

    kind = kind.upper()
    if kind == "L":
        m, n = N - 1, N
    elif kind == "R":
        m, n = N, N
    else:
        raise ValueError(f"kind must be 'L' or 'R', not {kind!r}")
    if m + n > n_taylor:
        raise ValueError(
            f"Pade [{m}/{n}] needs {m + n} Taylor terms; n_taylor={n_taylor}")

    zeros = jn_zeros(0, n_modes)
    eigenvalues = zeros**2
    cn = 4.0 / zeros ** 2

    # Build in the dimensionless variable u = s*mu*sigma*a^2. Using the
    # dimensional eigenvalues (j_n/a)^2 makes the Pade system needlessly
    # ill-conditioned and breaks geometric similarity at roundoff level.
    # Y/Y_DC = a_0 + a_1 u + ... with
    # a_0 = 1 and a_k = (-1)^k sum_n cn/(j_n^2)^k.
    coef = np.empty(n_taylor + 1)
    coef[0] = 1.0
    for k in range(1, n_taylor + 1):
        coef[k] = (-1) ** k * float(np.sum(cn / eigenvalues**k))

    # Q from the coefficients of z^(m+1) .. z^(m+n), with q_0 = 1
    mat = np.zeros((n, n))
    rhs = np.zeros(n)
    for i, k in enumerate(range(m + 1, m + n + 1)):
        for j in range(1, n + 1):
            if 0 <= k - j <= n_taylor:
                mat[i, j - 1] = coef[k - j]
        rhs[i] = -coef[k]
    Q = np.concatenate([[1.0], np.linalg.solve(mat, rhs)])

    P = np.array([sum(Q[j] * coef[k - j] for j in range(min(k, n) + 1))
                  for k in range(m + 1)])

    u = s * mu * sigma * a**2
    num = sum(P[k] * u**k for k in range(m + 1))
    den = sum(Q[k] * u**k for k in range(n + 1))
    return Y_DC_cylinder(a, sigma) * num / den
