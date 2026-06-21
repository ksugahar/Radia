# -*- coding: utf-8 -*-
r"""Energy-based (vector Play) magnetic hysteresis for the radia-ngsolve core.

This brings the lab's energy-consistent hysteresis lineage -- Henrotte-Nicolet-Hameyer 2006,
Francois-Lavet et al. 2013, Jacques 2018 (see motor.henrotte_lineage_knowledge) -- into the NGSolve
FE core.  The Radia C++ side already has ``MatPlayHysteresis`` / ``MatEnergyHysteresis`` (fed by
``radia.hysteresis_io``); this is the Python/NGSolve-side B-input vector Play operator, usable both
standalone (B-H loop post-processing, iron-loss from the real loop area) and, per Gauss point, inside
a hysteresis-aware FE solve.

Model (B-input vector Play, matching ``radia.hysteresis_io.play_hysteron`` and the unified Radia
Type 5/6 formulation, where the energy-model ``chi`` == the Play threshold ``eta``):

    for each cell k with B-space threshold eta_k and play state p_k (a vector in B-space):
        p_k <- B - eta_k * (B - p_k) / max(eta_k, |B - p_k|)        # the play / backlash update
    H(B) = sum_k f_k(|p_k|) * p_k / |p_k|                            # shape-function H output

The play state lags B by eta_k -- the friction-like pinning (|h_irr| <= eta) that IS the hysteresis.
A cell with eta_k = 0 tracks B exactly (the reversible / anhysteretic part, no loop).

Why energy-based (not Jiles-Atherton / Preisach): the lineage notes J-A does not extend to vector and
Preisach has no energy interpretation; the Play / energy model is intrinsically vector and keeps a
non-negative dissipation (2nd law) by construction -- the friction always dissipates.

VERIFICATION HOOK -- linear cells.  With ``f_k(r) = a_k * r`` (so H_k = a_k * p_k) the major-loop area
of a single dissipative cell under a symmetric cycle of amplitude Bm > eta_k is exactly

    loop_area_k = 4 * a_k * eta_k * (Bm - eta_k)        [A/m * T = J/m^3 per cycle]

(zero if Bm <= eta_k).  The total hysteresis loss is the sum over cells -- a PHYSICAL loop area, not a
fitted Steinmetz k_h.  A distribution of thresholds reproduces the empirical P_hyst ~ Bm^2 (Steinmetz)
and the low-field Rayleigh ~ Bm^3 without any curve fit.  This is the principled replacement for the
fitted ``coreloss.steinmetz_loss_density`` k_h.
"""
from __future__ import annotations

import numpy as np

MU0 = 4e-7 * np.pi

_trapz = getattr(np, "trapezoid", None) or np.trapz   # numpy 2.x renamed trapz -> trapezoid


class PlayHysteresis:
    r"""B-input vector Play (energy-based) hysteresis operator.

    Parameters
    ----------
    eta : array (K,)
        B-space play thresholds [T], ascending; ``eta[0] == 0`` is the reversible cell.
    a : array (K,), optional
        Linear shape-function slopes: ``f_k(r) = a_k * r`` so ``H_k = a_k * p_k`` [(A/m)/T].
        Gives the analytic loop area ``4 a_k eta_k (Bm - eta_k)`` per dissipative cell.
    f_k : list of callables, optional
        General shape functions ``f_k(r) -> H`` (used instead of ``a`` if given, e.g. identified
        from measured loops via ``radia.hysteresis_io.build_shape_functions``).
    """

    def __init__(self, eta, a=None, f_k=None):
        self.eta = np.asarray(eta, dtype=float)
        self.K = len(self.eta)
        self.a = None if a is None else np.asarray(a, dtype=float)
        self.f_k = f_k
        if self.a is None and self.f_k is None:
            raise ValueError("provide either linear slopes `a` or shape functions `f_k`")

    # --- per-step vector update (the FE-per-Gauss-point primitive) ---
    def step(self, Bx, By, px, py):
        """One B-input step: update the play states and return (Hx, Hy, px, py)."""
        eta = self.eta
        Bp = np.sqrt((Bx - px) ** 2 + (By - py) ** 2)
        denom = np.maximum(eta, Bp)
        denom = np.where(denom < 1e-30, 1e-30, denom)
        px = Bx - eta * (Bx - px) / denom
        py = By - eta * (By - py) / denom
        if eta[0] == 0.0:                       # the eta=0 cell tracks B exactly (reversible)
            px[0], py[0] = Bx, By
        pmag = np.sqrt(px * px + py * py)
        fk = (self.a * pmag) if self.a is not None else np.array(
            [self.f_k[k](pmag[k]) for k in range(self.K)])
        safe = np.where(pmag < 1e-30, 1e-30, pmag)
        Hx = float(np.sum(np.where(pmag < 1e-30, 0.0, fk * px / safe)))
        Hy = float(np.sum(np.where(pmag < 1e-30, 0.0, fk * py / safe)))
        return Hx, Hy, px, py

    # --- scalar 1D B-H history (By = 0) ---
    def bh_history(self, B_wave, p0=None):
        """H(t) for a scalar B(t) history.  Returns (H_wave, final_state)."""
        K = self.K
        px = np.zeros(K) if p0 is None else p0[0].copy()
        py = np.zeros(K) if p0 is None else p0[1].copy()
        H = np.empty(len(B_wave))
        for i, B in enumerate(B_wave):
            H[i], _, px, py = self.step(float(B), 0.0, px, py)
        return H, (px, py)

    def steady_loop(self, Bm, n=721, cycles=3):
        """Symmetric major (or minor, Bm<eta_max) loop at amplitude Bm, after `cycles` to settle.
        Returns (B_loop, H_loop) over the final closed cycle (B: +Bm -> -Bm -> +Bm)."""
        th = np.linspace(0.0, 2 * np.pi, n)
        one = Bm * np.cos(th)
        wave = np.concatenate([one[:-1]] * cycles + [one])
        H, _ = self.bh_history(wave)
        return one, H[-n:]

    def loss_per_cycle(self, Bm, n=721, cycles=3):
        """Hysteresis loss density [J/m^3 per cycle] = loop area |∮ H dB| at amplitude Bm."""
        B, H = self.steady_loop(Bm, n=n, cycles=cycles)
        return abs(float(_trapz(H, B)))

    def loss_density(self, Bm, f, n=721, cycles=3):
        """Hysteresis loss density [W/m^3] = f * loop_area -- the physical (loop-area) hysteresis
        loss, the principled replacement for the fitted Steinmetz k_h*f*B^2."""
        return f * self.loss_per_cycle(Bm, n=n, cycles=cycles)

    # --- analytic loop area for the linear-cell model (verification) ---
    def analytic_loss_per_cycle(self, Bm):
        """Closed-form loop area sum_k 4 a_k eta_k (Bm - eta_k) over dissipative cells (eta_k < Bm).
        Linear cells only."""
        if self.a is None:
            raise ValueError("analytic loss only for linear cells (slopes `a`)")
        contrib = 0.0
        for k in range(self.K):
            if 0.0 < self.eta[k] < Bm:
                contrib += 4.0 * self.a[k] * self.eta[k] * (Bm - self.eta[k])
        return contrib


def rayleigh_cells(eta_max, K, a_each):
    """UNIFORM threshold density on (0, eta_max] (plus an eta=0 reversible cell).  The loop-area loss
    is then  loss(Bm) = (2/3) a rho Bm^3  -- the RAYLEIGH cubic law P_hyst ~ Bm^3 -- by the closed
    form 4 a eta (Bm-eta) summed over a flat eta distribution.  Useful as the low-field model."""
    eta = np.linspace(0.0, eta_max, K + 1)
    a = np.full(K + 1, a_each)
    return PlayHysteresis(eta=eta, a=a)


def steinmetz_cells(eta_min, eta_max, K, a_each):
    """GEOMETRIC (log-spaced) thresholds -> density ~ 1/eta -> loop-area loss  loss(Bm) ~ Bm^2, i.e.
    the empirical STEINMETZ hysteresis law P_hyst = k_h f Bm^2, but with k_h DERIVED from the cell
    parameters (a_each, the geometric spacing) -- NOT a fitted coefficient.  This is the principled
    replacement for ``coreloss.steinmetz_loss_density``'s fitted k_h: the same Bm^2 law falls out of
    the physical loop area of an energy-consistent hysteresis model."""
    eta = np.concatenate([[0.0], np.geomspace(eta_min, eta_max, K)])
    a = np.full(K + 1, a_each)
    return PlayHysteresis(eta=eta, a=a)


def dissipation_increments(model, B_wave):
    """Per-step dissipation along a B history.  For the linear-cell Play model the friction
    dissipation is  dD = sum_k a_k eta_k |dp_k|  -- the pinning force (a_k eta_k) times the play
    displacement |dp_k|, NON-NEGATIVE by construction (the 2nd law: friction always dissipates), and
    summing over a closed cycle to the loop area sum_k 4 a_k eta_k (Bm - eta_k).  Returns the array
    of dD (each >= 0)."""
    if model.a is None:
        raise ValueError("dissipation check implemented for linear cells")
    K = model.K
    px = np.zeros(K); py = np.zeros(K)
    dD = []
    for B in B_wave:
        px_old, py_old = px.copy(), py.copy()
        _, _, px, py = model.step(float(B), 0.0, px, py)
        dp = np.sqrt((px - px_old) ** 2 + (py - py_old) ** 2)
        dD.append(float(np.sum(model.a * model.eta * dp)))
    return np.array(dD)
