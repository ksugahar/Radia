"""Complex magnetic polarizability of a solid conducting permeable sphere.

A sphere (radius ``a``, conductivity ``sigma``, relative permeability
``mu_r``) in a uniform harmonic field ``H_ext = H0 e^{j w t} z-hat`` responds
with a magnetic dipole ``m = alpha(w) H0``.  ``alpha(w)`` interpolates between
the magnetostatic sphere and the perfectly shielded (diamagnetic) limit and is
the closed-form reference for coupled eddy-current + magnetization solvers
(``validation_test/vim_coupled/``).

Derivation (assembled and solved per frequency -- no transcribed end formula):
interior ``A_phi = C j1(k r) sin(th)`` with ``k^2 = -j w mu sigma`` (the
``e^{j w t}`` diffusion branch), exterior = uniform field + dipole;
continuity of ``B_r`` and ``H_theta`` at ``r = a`` gives a 2x2 linear system
in ``(C, m)``.  Only the ratio ``G/j1 = j0/j1 - 1/x`` of spherical Bessel
functions enters; it is evaluated with the overflow-safe form
``cot z = j (1 + q)/(1 - q)``, ``q = e^{-2 j z}`` (``|q| <= 1`` on the
diffusion branch), so the strong-skin limit is numerically stable.

Limits (locked by ``tests/analytical_formulas/test_sphere_eddy.py``):

    w -> 0   :  alpha -> 4 pi a^3 (mu_r - 1) / (mu_r + 2)
    w -> oo  :  alpha -> -2 pi a^3          (perfect shielding)
    mu_r = 1, low w :  Im(alpha) prop. +w,  Re(alpha) prop. -w^2

The dipole convention is SI with the ``1/(4 pi)`` in the field expression:
``H_dip = (m / 4 pi) (2 cos(th)/r^3 r-hat + sin(th)/r^3 th-hat)``, i.e.
``m`` in A m^2 for ``H0`` in A/m, so ``m = int M dV + 1/2 int r x J dV`` of
the induced sources.

References
----------
Wait J. R., "A conducting sphere in a time varying magnetic field",
  Geophysics 16 (1951), pp. 666-672.
Landau L. D., Lifshitz E. M., Pitaevskii L. P., "Electrodynamics of
  Continuous Media", 2nd ed., Pergamon (1984), sec. 59 (eddy currents;
  the sphere in an alternating field appears among the worked problems).
"""

from __future__ import annotations

import numpy as np

MU_0 = 4.0e-7 * np.pi


def _cot(z):
    """cot(z), overflow-safe for Im(z) <= 0 (the diffusion branch):
    cot z = j (1+q)/(1-q) with q = e^{-2jz}, |q| = e^{2 Im z} <= 1."""
    q = np.exp(-2j * z)
    return 1j * (1.0 + q) / (1.0 - q)


def _j0_over_j1(z):
    """j0(z)/j1(z) = z / (1 - z cot z), stable at small and large |z|."""
    z = complex(z)
    if abs(z) < 1e-4:
        return 3.0 / z + z / 5.0                   # series of j0/j1
    return z / (1.0 - z * _cot(z))


def sphere_complex_polarizability(freq_hz, radius, sigma, mu_r=1.0):
    """Complex magnetic polarizability ``alpha(w)`` of a solid sphere.

    Parameters
    ----------
    freq_hz : float or array_like
        Frequency in Hz (``>= 0``; 0 returns the magnetostatic limit).
    radius : float
        Sphere radius ``a`` in metres.
    sigma : float
        Electrical conductivity in S/m (``> 0``).
    mu_r : float, optional
        Relative permeability (``>= 1``), default 1.

    Returns
    -------
    complex or ndarray of complex
        ``alpha`` such that the induced dipole is ``m = alpha H0``
        (``m`` in A m^2 per unit ``H0`` in A/m).  Scalar input returns a
        scalar; array input returns an array of the same shape.
    """
    a = float(radius)
    sigma = float(sigma)
    mu_r = float(mu_r)
    if a <= 0.0:
        raise ValueError("radius must be positive")
    if sigma <= 0.0:
        raise ValueError("sigma must be positive")
    if mu_r < 1.0:
        raise ValueError("mu_r must be >= 1")
    mu = mu_r * MU_0
    freq = np.atleast_1d(np.asarray(freq_hz, dtype=float))
    if np.any(freq < 0.0):
        raise ValueError("freq_hz must be non-negative")
    out = np.empty(freq.shape, dtype=complex)
    static = 4.0 * np.pi * a**3 * (mu_r - 1.0) / (mu_r + 2.0)
    for i, f in enumerate(freq.ravel()):
        if f == 0.0:
            out.ravel()[i] = static
            continue
        w = 2.0 * np.pi * f
        k = np.sqrt(-1j * w * mu * sigma)          # principal branch, Im(ka) <= 0
        x = k * a
        ratio = _j0_over_j1(x) - 1.0 / x           # G(x)/j1(x), overflow-safe
        # unknowns (C' = C j1(x), m):
        #   B_r     :  2 C'/a                 = mu0 (H0 + 2 m/(4 pi a^3))
        #   H_theta :  -(C' k/mu) (G/j1)      = -H0 + m/(4 pi a^3)
        A = np.array(
            [[2.0 / a, -2.0 * MU_0 / (4.0 * np.pi * a**3)],
             [-(k / mu) * ratio, -1.0 / (4.0 * np.pi * a**3)]],
            dtype=complex)
        b = np.array([MU_0, -1.0], dtype=complex)  # H0 = 1
        out.ravel()[i] = np.linalg.solve(A, b)[1]
    if np.isscalar(freq_hz) or np.asarray(freq_hz).ndim == 0:
        return complex(out.ravel()[0])
    return out.reshape(np.asarray(freq_hz).shape)
