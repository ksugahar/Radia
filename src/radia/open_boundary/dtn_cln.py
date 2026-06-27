# -*- coding: utf-8 -*-
"""Exact exterior Dirichlet-to-Neumann (DtN) open boundary as a Cauer Ladder
Network (CLN) -- the "Zs-DtN-CLN" open boundary.

WHAT THIS IS.  For a SEPARABLE truncation surface (a sphere in 3D / a circle in
2D) the exterior DtN diagonalises in spherical (cylindrical) harmonics; each
multipole n has a scalar frequency symbol that is EXACTLY a reverse-Bessel
rational function -- in s for the wave (Helmholtz) exterior, in q=sqrt(s) for the
magneto-quasistatic (eddy / diffusion) exterior, with the SAME poles roots(theta_n).
This module realises that exact symbol as a finite Cauer continued-fraction ladder
(EXACT at n+1 stages, well-conditioned, passive) and exposes the Grote-Keller
companion auxiliary-ODE form for a transient Robin open boundary.

SCOPE -- read this; it is the boundary of where this is the RIGHT tool.
  It WINS over a (CFS-)PML only inside its island: COMPACT / quasi-spherical,
  MAGNETO-QUASI-STATIC (MQS / Laplace-kernel) problems where you want an EXACT,
  DC-well-conditioned open boundary AND a compact passive circuit / ROM.  It does
  NOT win in general:
    * the truncation is locked to a SPHERE (Liouville) -- an elongated / arbitrary
      object wastes mesh on the spherical shell, where a box CFS-PML hugs better;
    * genuine wave RADIATION (finite real kR) is OUTSIDE radia's MQS scope (-> PML
      / NGSolve);
    * for a NON-separable body the per-mode exactness becomes a convergent BAND
      approximation (build the DtN via Kelvin-FEM / Schur first;
      archived act6_06.. sources in docs/kelvin/kelvin_dtn_spectrum_archive_results.json).
  See docs/open_boundary/OPEN_BOUNDARY_MAP.md for the full selector.

NOT NOVEL (cite, do not claim).  The exact rational radiation DtN + local
auxiliary-ODE realisation is Grote-Keller (SIAM J. Appl. Math. 1995) /
Hagstrom-Warburton (complete radiation BCs / continued-fraction ABCs); the
sqrt(s) diffusion (Warburg) impedance as a Cauer ladder is classical network
synthesis (e.g. Phys. Chem. Chem. Phys. 18 (2016) 9498); the Cauer Ladder Network
MOR is Kameari-Ebrahimi-Sugahara-Shindo-Matsuo, IEEE T-Magn 54(3):7201804 (2018).
This module is the VERIFIED, reusable operator -- not a novelty claim.

VERIFIED (tests/open_boundary/test_dtn_cln.py, ported from the research demos
archived `act6_02_cln_dtn_cauer.py` (full source: `docs/kelvin/kelvin_dtn_spectrum_archive_results.json`) +
act6_11_exact_dtn_fetd.py):
  - the eddy DtN is EXACTLY rational in q=sqrt(s); the Cauer-in-q ladder is exact
    at n+1 stages (NRMSE ~1e-15) for n=1..6, well-conditioned (coeff spread <1e3),
    whereas a Foster fit in s floors ~1e-3 and ill-conditions (~1e5);
  - the wave (in s) and diffusion (in sqrt(s)) DtN share the SAME poles roots(theta_n);
  - the companion auxiliary-ODE rates are roots(theta_n), all Re<0 => passive,
    unconditionally-stable transient open boundary (Grote-Keller form).

UNITS.  radia is meters / SI.  R0 = truncation radius (m); mu_sigma = mu*sigma
(the MQS diffusion coefficient) so the eddy wavenumber is gamma = sqrt(s*mu_sigma)
and the reverse-Bessel argument is q = R0*gamma.  s is the Laplace variable
(s = i*omega on the imaginary axis).
"""
import numpy as np
import numpy.polynomial.polynomial as _P
from math import factorial
from scipy.special import hankel1, kv

__all__ = [
    "reverse_bessel_theta", "reverse_bessel_roots",
    "eddy_dtn", "eddy_dtn_rational_q", "wave_dtn",
    "cauer_ladder", "eval_ladder",
    "companion_poles", "sqrt_s_passive_ladder", "eval_sqrt_ladder",
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
    and the diffusion (in q=sqrt(s)) exterior DtN, and the companion-ODE rates."""
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
# Cauer continued fraction in q -- EXACT realisation of the eddy DtN at n+1 stages
# ---------------------------------------------------------------------------
def cauer_ladder(n):
    """Cauer continued-fraction stages of the eddy DtN in q = R0 sqrt(s*mu_sigma).
    EXACT at n+1 stages.  Returns a list of stage polynomials (each a short
    ascending-power ndarray in q); evaluate with eval_ladder."""
    A, den = eddy_dtn_rational_q(n)
    num = np.trim_zeros(np.asarray(A, float), 'b')
    den = np.trim_zeros(np.asarray(den, float), 'b')
    stages = []
    while len(num) and np.any(np.abs(den) > 1e-13) and len(stages) < 40:
        quo, rem = _P.polydiv(num, den)
        stages.append(quo)
        num, den = den, np.trim_zeros(rem, 'b')
        if len(den) == 0:
            break
    return stages


def eval_ladder(stages, s, R0=1.0, mu_sigma=1.0):
    """Evaluate a Cauer ladder (from cauer_ladder) at the Laplace variable s.
    Reproduces eddy_dtn to machine precision (the ladder IS the exact operator)."""
    q = R0 * np.sqrt(complex(s) * mu_sigma)
    val = None
    for stage in reversed(stages):
        sv = sum(stage[k] * q ** k for k in range(len(stage)))
        val = sv if val is None else sv + 1.0 / val
    return val


# ---------------------------------------------------------------------------
# transient Robin realisation: companion auxiliary ODEs (Grote-Keller form)
# ---------------------------------------------------------------------------
def companion_poles(n):
    """Relaxation rates of the auxiliary ODEs for the TIME-DOMAIN Robin realisation
    of the exterior DtN (Grote-Keller form): lambda_j = roots(theta_n), all Re<0
    => passive, unconditionally stable.  For the wave exterior (R0 = c = 1):
        g(t)      = -du/dt - u + sum_j psi_j ,   u = field trace at the truncation,
        dpsi_j/dt =  lambda_j ( psi_j + u ) ,    one first-order ODE per pole.
    (Verified reflectionless in a 1-D radial FETD solve in
    archived `act6_11_exact_dtn_fetd.py` (full source: `docs/kelvin/kelvin_dtn_spectrum_archive_results.json`).)"""
    return reverse_bessel_roots(n)


# ---------------------------------------------------------------------------
# finite PASSIVE realisation of the sqrt(s) diffusion-memory element
# ---------------------------------------------------------------------------
def sqrt_s_passive_ladder(omega, K):
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


def eval_sqrt_ladder(g, p, s):
    """Evaluate the passive sqrt(s) ladder sum_m g_m s/(s+p_m) at s."""
    s = complex(s)
    return np.sum(g * s / (s + p))
