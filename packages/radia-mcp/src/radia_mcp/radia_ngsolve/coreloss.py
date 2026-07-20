# -*- coding: utf-8 -*-
"""Iron / core-loss separation post-processor (Steinmetz + Bertotti) for
rotor-sweep motor analysis.

Marquee cross-learning capability (radia-ngsolve BEATS FEMM here). FEMM computes
iron loss only from a time-HARMONIC solve and ignores the DC (permanent-magnet)
field, forcing the awkward ".ans frequency hand-edit" hack for PM machines.
radia-ngsolve instead takes the per-element flux-density waveform B(theta) from
a rotor-angle sweep (e.g. the committed cogging sweep), FFTs it into harmonics,
and applies the Bertotti three-term loss separation natively -- the same path
JMAG's iron-loss study uses.

Loss separation (specific loss, W/m^3, sinusoidal flux at frequency f, peak Bpk):

    P_hyst   = k_h * f   * Bpk^2          hysteresis (Steinmetz P = k_h f^a Bpk^b)
    P_class  = k_c * (f * Bpk)^2          classical eddy; k_c = pi^2 sigma d^2 / 6
                                          for a lamination of thickness d
    P_excess = k_e * (f * Bpk)^1.5        excess / anomalous
    P_total  = P_hyst + P_class + P_excess

For a non-sinusoidal waveform decomposed into harmonics B_n at n*f1, the
classical and excess terms add per harmonic (they depend on dB/dt):

    P_class  = k_c * sum_n (n f1 B_n)^2 ,   P_excess = k_e * sum_n (n f1 B_n)^1.5

The classical term reduces to the textbook lamination eddy loss
    P_class = pi^2 sigma d^2 f^2 Bpk^2 / 6
when k_c = pi^2 sigma d^2 / 6 -- the exact closed form used for validation
(validation_test/radia_mcp/test_core_loss.py).
"""
from __future__ import annotations

import math

import numpy as np


def lamination_kc(sigma, d):
    """Classical-eddy Bertotti coefficient for a lamination of thickness ``d``
    [m] and conductivity ``sigma`` [S/m]:  ``k_c = pi^2 sigma d^2 / 6``, so that
    ``P_class = k_c (f Bpk)^2 = pi^2 sigma d^2 f^2 Bpk^2 / 6`` [W/m^3]."""
    return math.pi ** 2 * sigma * d ** 2 / 6.0


def classical_eddy_loss_density(Bpk, f, sigma, d):
    """Closed-form classical eddy-current loss density of a lamination [W/m^3]:
    ``P = pi^2 sigma d^2 f^2 Bpk^2 / 6``.  This is the THIN-lamination limit (d << skin depth);
    for d ~ delta multiply by :func:`lamination_eddy_skin_factor` (see
    :func:`classical_eddy_loss_density_skin`)."""
    return math.pi ** 2 * sigma * d ** 2 * f ** 2 * Bpk ** 2 / 6.0


def lamination_eddy_skin_factor(d, delta):
    """Skin-effect correction factor for the classical lamination eddy loss, thickness ``d`` [m]
    and skin depth ``delta`` [m] (``delta = sqrt(2/(omega mu sigma))``):

        F(xi) = (3/xi) (sinh xi - sin xi)/(cosh xi - cos xi),    xi = d/delta,

    the ratio of the ACTUAL eddy loss (imposed surface flux) to the thin-limit classical value.
    F -> 1 as xi -> 0 (thin lamination) and F -> 3/xi as xi -> infinity (skin-limited), so the
    classical formula OVERESTIMATES once d >~ delta (e.g. F(5)=0.61, F(8)=0.37).  This is the
    standard slab/lamination skin-effect result (Stoll, 1974; Lammeraner & Stafl); validated to
    < 1.1 % over xi = 1..8 against an independent eddy-current FE solve (stored regression reference)."""
    xi = d / delta
    if xi < 1e-4:
        return 1.0
    return (3.0 / xi) * (math.sinh(xi) - math.sin(xi)) / (math.cosh(xi) - math.cos(xi))


def classical_eddy_loss_density_skin(Bpk, f, sigma, d, mu):
    """Classical lamination eddy loss density WITH the skin-effect correction [W/m^3]:

        P = (pi^2 sigma d^2 f^2 Bpk^2 / 6) * F(d/delta),   delta = sqrt(2/(omega mu sigma)),

    where ``Bpk`` is the imposed peak SURFACE flux density [T], ``sigma`` [S/m], ``d`` the lamination
    thickness [m], and ``mu = mu0*mu_r`` [H/m] (the skin depth, hence the correction, depends on mu --
    the thin-limit term does not).  Reduces to :func:`classical_eddy_loss_density` for d << delta and
    rolls off as ~1/d at high frequency (flux no longer penetrates).  Validated to < 1.1 % over
    d/delta = 1..8 against an independent eddy-current FE reference (stored regression reference)."""
    omega = 2.0 * math.pi * f
    delta = math.sqrt(2.0 / (omega * mu * sigma))
    return classical_eddy_loss_density(Bpk, f, sigma, d) * lamination_eddy_skin_factor(d, delta)


def steinmetz_loss_density(Bpk, f, k_h, alpha=1.0, beta=2.0):
    """Steinmetz hysteresis loss density [W/m^3]:  ``P = k_h f^alpha Bpk^beta``."""
    return k_h * f ** alpha * Bpk ** beta


def bertotti_loss_density(Bpk, f, k_h, k_c, k_e=0.0):
    """Bertotti 3-term loss separation for SINUSOIDAL flux [W/m^3]:
    ``P = k_h f Bpk^2 + k_c (f Bpk)^2 + k_e (f Bpk)^1.5``."""
    x = f * Bpk
    return k_h * f * Bpk ** 2 + k_c * x ** 2 + k_e * x ** 1.5


def fft_amplitudes(waveform):
    """Single-sided harmonic PEAK amplitudes of a real periodic waveform sampled
    over ONE period (N samples). Returns array ``A`` where ``A[n]`` is the peak
    amplitude of the n-th harmonic and ``A[0]`` is the DC level. (numpy.fft.rfft
    with the standard 2/N scaling; DC unscaled.)"""
    w = np.asarray(waveform, dtype=float)
    N = len(w)
    c = np.fft.rfft(w)
    A = np.abs(c) * (2.0 / N)
    A[0] = abs(c[0]) / N
    return A


def harmonic_loss_density(B_harmonics, f1, k_h, k_c, k_e=0.0):
    """Loss density [W/m^3] of a NON-sinusoidal flux waveform given its harmonic
    peak amplitudes ``B_harmonics[n]`` at frequency ``n*f1`` (index 0 = DC,
    ignored). Classical and excess add per harmonic (dB/dt-based); hysteresis
    uses the fundamental (major-loop) term as the first model:

        P_class  = k_c sum_n (n f1 B_n)^2
        P_excess = k_e sum_n (n f1 B_n)^1.5
        P_hyst   = k_h f1 B_1^2

    Minor-loop hysteresis corrections (the iGSE refinement) are a TODO; for
    sinusoidal flux this reduces exactly to :func:`bertotti_loss_density`.
    Returns a dict with P_hyst / P_class / P_excess / P_total."""
    B = np.asarray(B_harmonics, dtype=float)
    n = np.arange(len(B))
    fh = n * f1
    P_class = k_c * float(np.sum((fh * B) ** 2))
    P_excess = k_e * float(np.sum((fh * B) ** 1.5)) if k_e else 0.0
    B1 = float(B[1]) if len(B) > 1 else 0.0
    P_hyst = k_h * f1 * B1 ** 2
    return {"P_hyst": P_hyst, "P_class": P_class, "P_excess": P_excess,
            "P_total": P_hyst + P_class + P_excess}


def core_loss_in_region(Bpk_cf, mesh, region, f, k_h, k_c, k_e=0.0):
    """Integrate the (sinusoidal) Bertotti loss density over a mesh region given
    a per-point PEAK flux-density magnitude CF ``Bpk_cf``. Returns total loss per
    unit stack length [W/m]. (Per-element Bpk from a rotor sweep is the harness's
    job; this just integrates the density.)"""
    from ngsolve import Integrate
    x = f * Bpk_cf
    dens = k_h * f * Bpk_cf ** 2 + k_c * x ** 2
    if k_e:
        dens = dens + k_e * x ** 1.5
    return Integrate(dens, mesh, definedon=mesh.Materials(region))
