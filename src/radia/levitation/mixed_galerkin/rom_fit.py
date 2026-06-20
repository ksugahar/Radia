"""rom_fit.py -- passive, stable Foster ROM of a sampled polarizability alpha(s).

Turns frequency samples alpha(j omega_i) into a passive + stable LTI

    alpha(s) ~ alpha_inf + sum_k g_k / (1 + s tau_k),   g_k >= 0, tau_k > 0,

i.e. one RC relaxation state per pole plus a feedthrough D = alpha_inf.  This is
the EXTERIOR-MATCHED *physical* polarizability-tensor ROM: it fits the verified
per-frequency 3D HCurl solve (examples/levitation/ellipsoid/
ellipsoid_alpha_tensor_3d.py, which carries the air reaction dipole, the lift /
Re[alpha] part), NOT the interior-PEC bulk eigenmodes
(mixed_galerkin.bulk_foster_via_eigen / bulk_foster_vector_via_eigen).

Why a sample fit and not a Kameari + Kelvin eigen-accumulation: the 3D HCurl
Kameari + Kelvin accumulation structurally BREAKS DOWN on the
isolated-conductor-in-vacuum problem (L_n sign flip at stage 1; see
examples/levitation/research_cln/ngsolve_validation/
cuboid_521_kameari_kelvin_v15_canonical.py).  Sphere axisym Kameari reaches the
Stoll Cauer ladder to 0.000%, but the general 3D body does not.  This module
sidesteps that by building the LTI on the verified per-frequency solve instead.

Recipe (VERIFIED on the analytic sphere Stoll spectrum to < 0.2 %):
  1. scipy.interpolate.AAA discovers the DOMINANT real poles; they land on the
     physical Stoll decay times tau_n = mu0 sigma a^2 / (n pi)^2 to ~0.00 %.
  2. The pole set = those dominant poles UNION a log-spaced filler spanning the
     sample band.  The filler captures the high-order tail AAA buries in
     Froissart pairs and is ALWAYS present, so the fit degrades gracefully if
     AAA discovers nothing (one coherent method, not a fallback chain).
  3. NNLS real residues g_k >= 0 -> passive by construction; the poles are real
     and negative -> stable by construction.

The dominant poles are reported separately (`FosterROM.dominant_tau`) as the
physical Stoll decay-time spectrum; the filler poles are an approximation
basis, not claimed to be individually physical.
"""
from __future__ import annotations

import warnings

import numpy as np
from scipy.interpolate import AAA
from scipy.optimize import nnls

__all__ = ["FosterROM", "passive_foster_fit", "diagonal_tensor_state_space"]


class FosterROM:
    """Passive scalar Foster ROM  alpha(s) = alpha_inf + sum_k g_k/(1+s tau_k).

    Attributes
    ----------
    tau_n, g_n : ndarray
        Relaxation times (s) and non-negative residues of the active poles.
    alpha_inf : float
        Feedthrough D = high-frequency limit alpha(s->inf) (e.g. the
        perfect-conductor flux-exclusion value, negative for a diamagnet).
    dominant_tau : ndarray
        The AAA-discovered dominant poles -- the physical Stoll decay times.
    band_fit_relerr : float
        max|fit - data| / max|data| over the sample band.
    """

    def __init__(self, tau_n, g_n, alpha_inf, *, dominant_tau, band_fit_relerr):
        self.tau_n = np.asarray(tau_n, dtype=float)
        self.g_n = np.asarray(g_n, dtype=float)
        self.alpha_inf = float(alpha_inf)
        self.dominant_tau = np.asarray(dominant_tau, dtype=float)
        self.band_fit_relerr = float(band_fit_relerr)

    @property
    def n_states(self):
        return len(self.tau_n)

    def eval(self, s):
        """Evaluate alpha(s) of the ROM (s complex, array or scalar)."""
        s = np.atleast_1d(np.asarray(s, dtype=complex))
        val = (self.g_n[None, :] / (1.0 + s[:, None] * self.tau_n[None, :])
               ).sum(axis=1) + self.alpha_inf
        return val

    def state_space(self):
        """(A, B, C, D) for this passive, stable scalar Foster ROM.

        x_k' = -x_k/tau_k + u,   y = sum_k (g_k/tau_k) x_k + alpha_inf u.
        Transfer function C (sI - A)^-1 B + D = eval(s) exactly.
        """
        n = self.n_states
        A = np.diag(-1.0 / self.tau_n) if n else np.zeros((0, 0))
        B = np.ones((n, 1))
        C = (self.g_n / self.tau_n)[None, :] if n else np.zeros((1, 0))
        D = np.array([[self.alpha_inf]])
        return A, B, C, D


def passive_foster_fit(s, alpha, *, n_filler=20, max_aaa_terms=24,
                       aaa_rtol=1e-12, prune_rel=1e-4):
    """Fit a passive, stable Foster ROM to samples alpha(s_i).

    Convention: the fit form g_k/(1 + s tau_k) with g_k >= 0, tau_k > 0 has
    Im[alpha(j omega)] < 0 (the causal / passive e^{+j omega t} convention).
    Pass data in THAT convention.  A 3D HCurl eddy solve (e.g.
    ellipsoid_alpha_tensor_3d.alpha_tensor_component) returns the PHYSICS
    convention with Im > 0 -- conjugate it (np.conj(alpha)) before fitting,
    or the NNLS residues collapse and the band fit is ~50-60 %.

    Parameters
    ----------
    s : array_like (complex)
        Sample points on the imaginary axis, s_i = j omega_i (strictly
        increasing |omega| recommended; the last sample sets alpha_inf).
    alpha : array_like (complex)
        Samples alpha(s_i).
    n_filler : int
        Number of log-spaced filler poles spanning the sample band.
    max_aaa_terms, aaa_rtol : int, float
        scipy.interpolate.AAA controls for the dominant-pole discovery.
    prune_rel : float
        Drop active residues below prune_rel * max(g) from the returned ROM.

    Returns
    -------
    FosterROM
    """
    s = np.asarray(s, dtype=complex)
    alpha = np.asarray(alpha, dtype=complex)
    if s.shape != alpha.shape or s.ndim != 1:
        raise ValueError("s and alpha must be 1-D arrays of equal length")
    amax = float(np.max(np.abs(alpha)))
    if amax == 0.0:
        raise ValueError("alpha samples are all zero")

    # --- 1. AAA discovers the dominant real LHP poles ----------------------
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="AAA failed to converge*",
                category=RuntimeWarning,
            )
            rat = AAA(s, alpha, max_terms=max_aaa_terms, rtol=aaa_rtol)
        dominant_tau = np.array(sorted(
            (-1.0 / p.real for p in rat.poles()
             if p.real < 0 and abs(p.imag) < 1e-2 * abs(p.real)),
            reverse=True))
    except Exception:
        dominant_tau = np.empty(0)

    # --- 2. pole set = dominant UNION log-spaced filler over the band -------
    w = np.abs(s.imag)
    w = w[w > 0]
    tau_lo = 1.0 / w.max()          # fastest relaxation resolvable
    tau_hi = 1.0 / w.min()          # slowest
    filler = np.logspace(np.log10(tau_lo), np.log10(tau_hi), int(n_filler))
    tau_set = np.unique(np.concatenate([dominant_tau, filler]))[::-1]

    # --- 3. NNLS passive residues on the fixed real poles ------------------
    alpha_inf = float(alpha[np.argmax(w)].real)   # HF plateau (feedthrough D)
    res = alpha - alpha_inf
    Phi = 1.0 / (1.0 + s[:, None] * tau_set[None, :])
    Areal = np.vstack([Phi.real, Phi.imag])
    breal = np.concatenate([res.real, res.imag])
    g, _ = nnls(Areal, breal)

    # prune negligible residues
    if g.max() > 0:
        keep = g > prune_rel * g.max()
    else:
        keep = np.zeros_like(g, dtype=bool)
    tau_n, g_n = tau_set[keep], g[keep]

    rom = FosterROM(tau_n, g_n, alpha_inf, dominant_tau=dominant_tau,
                    band_fit_relerr=0.0)
    fit = rom.eval(s)
    rom.band_fit_relerr = float(np.max(np.abs(fit - alpha)) / amax)
    return rom


def diagonal_tensor_state_space(roms):
    """Block-diagonal MIMO (A, B, C, D) for a diagonal polarizability tensor.

    Parameters
    ----------
    roms : sequence of FosterROM
        One scalar ROM per principal axis (length P; P=3 for a 3D tensor).

    Returns
    -------
    A, B, C, D : ndarray
        Shapes (n_states, n_states), (n_states, P), (P, n_states), (P, P).
        Off-diagonal coupling is zero (valid for principal-axis-aligned bodies).
    n_states : int
    """
    P = len(roms)
    blocks = [r.state_space() for r in roms]
    ns = [b[0].shape[0] for b in blocks]
    n_states = int(sum(ns))
    A = np.zeros((n_states, n_states))
    B = np.zeros((n_states, P))
    C = np.zeros((P, n_states))
    D = np.zeros((P, P))
    off = 0
    for p, (Ap, Bp, Cp, Dp) in enumerate(blocks):
        n = Ap.shape[0]
        A[off:off + n, off:off + n] = Ap
        B[off:off + n, p:p + 1] = Bp
        C[p:p + 1, off:off + n] = Cp
        D[p, p] = Dp[0, 0]
        off += n
    return A, B, C, D, n_states
