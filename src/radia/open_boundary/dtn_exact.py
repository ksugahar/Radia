# -*- coding: utf-8 -*-
"""Exact exterior Dirichlet-to-Neumann (DtN) symbols for separable truncations.

For a SEPARABLE truncation surface (a sphere in 3D / a circle in 2D) the exterior
DtN diagonalises in spherical (cylindrical) harmonics; each multipole n has a
scalar frequency symbol that is EXACTLY a reverse-Bessel rational function -- in s
for the wave (Helmholtz) exterior, in q=sqrt(s) for the magneto-quasistatic (eddy /
diffusion) exterior, with the SAME poles roots(theta_n).  This module evaluates
those symbols and reverse-Bessel roots. Wave and diffusion use different
frequency variables: diffusion roots are in q, not physical-time s poles.
The positive-residue approximation of sqrt(s) is a finite-band approximation,
not an exact finite-dimensional realisation of diffusion memory.

SCOPE.  This is the right tool for COMPACT / quasi-spherical, MAGNETO-QUASI-STATIC
problems where an EXACT, DC-well-conditioned open boundary is wanted.  It is not
general: the truncation is locked to a SPHERE (an elongated object wastes mesh on
the spherical shell, where a box CFS-PML hugs better); genuine wave RADIATION is
outside radia's MQS scope (-> PML / NGSolve); for a NON-separable body build the
DtN via Kelvin-FEM / Schur first (`radia.open_boundary.kelvin_dtn`).  See
docs/open_boundary/OPEN_BOUNDARY_MAP.md for the selector.

NOT NOVEL (cite, do not claim).  The exact rational radiation DtN + local
auxiliary-ODE realisation is Grote-Keller (SIAM J. Appl. Math. 1995) /
Hagstrom-Warburton (complete radiation BCs).

The tests compare the frequency symbols and fitted responses on specified
samples. Negative root real parts alone do not prove passivity, global error
bounds, or unconditional stability of a time-discretised coupled problem.

UNITS.  radia is meters / SI.  R0 = truncation radius (m); mu_sigma = mu*sigma
(the MQS diffusion coefficient) so the eddy wavenumber is gamma = sqrt(s*mu_sigma)
and the reverse-Bessel argument is q = R0*gamma.  s is the Laplace variable
(s = i*omega on the imaginary axis).
"""
import numpy as np
from math import factorial
from scipy.special import hankel1, kv

__all__ = [
    "reverse_bessel_theta", "reverse_bessel_roots",
    "eddy_dtn", "eddy_dtn_rational_q", "wave_dtn",
    "companion_poles", "sqrt_s_passive_poles", "eval_sqrt_poles",
]


# ---------------------------------------------------------------------------
# reverse Bessel polynomial theta_n -- the shared structure of BOTH exteriors
# ---------------------------------------------------------------------------
def reverse_bessel_theta(n):
    """Ascending-power coefficients of the reverse Bessel polynomial
    theta_n(x) = sum_{k=0}^{n} (n+k)! / ((n-k)! k! 2^k) * x^{n-k}."""
    c = np.zeros(n + 1)
    for k in range(n + 1):
        c[n - k] = factorial(n + k) / (factorial(n - k) * factorial(k) * 2 ** k)
    return c


def reverse_bessel_roots(n):
    """Roots of theta_n (all Re<0 for n>=1): the shared poles of the wave (in s)
    and the diffusion (in q=sqrt(s)) exterior DtN; these are different variables."""
    if n == 0:
        return np.array([], dtype=complex)
    return np.roots(reverse_bessel_theta(n)[::-1].copy()).astype(complex)


# ---------------------------------------------------------------------------
# wave (Helmholtz) exterior DtN -- rational in s,  z = k R0
# ---------------------------------------------------------------------------
def _sph_h1(n, z):
    z = np.asarray(z, dtype=complex)
    return np.sqrt(np.pi / (2.0 * z)) * hankel1(n + 0.5, z)


def wave_dtn(l, z):
    """Exact wave (Helmholtz) exterior DtN eigenvalue for multipole l at a sphere:
    Lambda_l(z) = z h_l^(1)'(z) / h_l^(1)(z),  z = k R0.  Rational in z with poles
    i*roots(theta_l).  (Provided for the wave<->diffusion unification; genuine wave
    radiation is OUTSIDE radia's MQS scope -- see the module scope note.)"""
    z = np.asarray(z, dtype=complex)
    hp = _sph_h1(l - 1, z) - (l + 1) / z * _sph_h1(l, z)
    return z * hp / _sph_h1(l, z)


# ---------------------------------------------------------------------------
# diffusion / eddy-current exterior DtN -- rational in q = R0 sqrt(s mu_sigma)
# ---------------------------------------------------------------------------
def eddy_dtn(n, s, R0=1.0, mu_sigma=1.0):
    """Exact magneto-quasistatic (eddy / diffusion) exterior DtN eigenvalue for
    multipole n at a sphere of radius R0:
        G_n(s) = -q K_{n-1/2}(q)/K_{n+1/2}(q) - (n+1),   q = R0 sqrt(s*mu_sigma).
    EXACTLY rational in q of degree n (poles = roots(theta_n))."""
    q = R0 * np.sqrt(complex(s) * mu_sigma)
    return -q * kv(n - 0.5, q) / kv(n + 0.5, q) - (n + 1.0)


def eddy_dtn_rational_q(n):
    """The eddy DtN written as the rational A(q)/theta_n(q) in q (ascending-power
    coeffs).  G_n = A(q)/theta_n(q), with A = -q^2 theta_{n-1} - (n+1) theta_n."""
    th_n = reverse_bessel_theta(n)
    th_n1 = reverse_bessel_theta(n - 1) if n >= 1 else np.array([1.0])
    A = np.zeros(max(len(th_n1) + 2, len(th_n)))
    A[2:2 + len(th_n1)] += -th_n1          # -q^2 theta_{n-1}
    A[:len(th_n)] += -(n + 1) * th_n       # -(n+1) theta_n
    return A, th_n


# ---------------------------------------------------------------------------
# transient Robin realisation: companion auxiliary ODEs (Grote-Keller form)
# ---------------------------------------------------------------------------
def companion_poles(n):
    """Return reverse-Bessel roots (dimensionless wave companion rates).

    For diffusion these are q-plane roots, q=R0*sqrt(s*mu_sigma), not
    physical-time ODE rates. Re(root)<0 does not by itself establish
    passivity or unconditional stability of any time-stepping scheme.
    """
    return reverse_bessel_roots(n)


# ---------------------------------------------------------------------------
# finite PASSIVE realisation of the sqrt(s) diffusion-memory element
# ---------------------------------------------------------------------------
def sqrt_s_passive_poles(omega, K):
    """Fit sqrt(s) ~ sum_m g_m * s/(s + p_m) with g_m >= 0 (passive) and p_m
    log-spaced over the band omega -- the time-domain realisation of the diffusion
    memory.  Each term is one first-order ODE; the real poles -p_m < 0 => stable.
    Returns (g, p, nrmse)."""
    from scipy.optimize import nnls
    omega = np.asarray(omega, float)
    p = np.logspace(np.log10(omega[0]) - 0.5, np.log10(omega[-1]) + 0.5, K)
    s = 1j * omega
    Amat = np.column_stack([s / (s + pj) for pj in p])
    target = np.sqrt(s)
    g, _ = nnls(np.vstack([Amat.real, Amat.imag]),
                np.concatenate([target.real, target.imag]))
    fit = Amat @ g
    nrmse = float(np.sqrt(np.mean(np.abs(fit - target) ** 2))
                  / np.sqrt(np.mean(np.abs(target) ** 2)))
    return g, p, nrmse


def eval_sqrt_poles(g, p, s):
    """Evaluate the passive sqrt(s) pole-residue sum sum_m g_m s/(s+p_m) at s."""
    s = complex(s)
    return np.sum(g * s / (s + p))
