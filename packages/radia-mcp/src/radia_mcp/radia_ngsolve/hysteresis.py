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

WHY B-INPUT (and not H-input) FOR THE radia-ngsolve CORE.  The lab's own measurements
(Hane-Sugahara-Morishita-Ahagon, "Experimental Verification of Congruency Property ... Using B-input
Play Model", IEEE TMAG 2026) show that real non-oriented silicon steel has CONGRUENCY along the H-axis
(all minor loops the same shape in H for a given B-range) but INCONGRUENCY along the B-axis.  The
B-input play model has H-axis congruency by construction (the shape functions depend only on the play
hysterons, Eq. 4), so it reproduces measured minor loops -- shapes, areas, and losses -- where the
H-input model (which forces B-axis congruency) fails.  And the B-input model is the one that matches an
A-formulation FE solve (B = curl A is primary, H = sum_k f_k(|p_k|) p_k/|p_k| is the constitutive
output) -- exactly the radia-ngsolve AGE / A_z core.  So B-input is both more accurate AND the natural
FE fit here.  The play model is mathematically equivalent to the static Preisach model (Bobbio 1997;
Matsuo-Shimasaki 2003) but much cheaper; the shape-function slope ``a_k`` is a RELUCTIVITY (for B-input).
Standard discretisation: equal-interval thresholds eta_k = (k - 1/2) * Bmax / N (``play_cells``), shape
functions identified from measured loops via the Everett function (next increment).  FE coupling: the
vector play model + Newton-Raphson (Mitsuoka-Mifune-Matsuo 2013) -- the planned hysteresis-aware solve.
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

    # --- differential reluctivity tensor dH/dB: the Newton-Raphson FE ingredient ---
    def play_tangent(self, Bx, By, px, py):
        """Consistent differential-reluctivity tensor dH/dB for the B-input Play with FROZEN states
        (px, py) -- the Newton-Raphson FE-coupled-hysteresis ingredient (Mitsuoka-Mifune-Matsuo 2013).
        Linear cells: dH/dB = a_0 I + sum_{moving} a_k [ I - (eta_k/|s_k|)(I - s^_k (x) s^_k) ], with
        s_k = B - p_k; a PINNED cell (|s_k| <= eta_k) contributes 0, the eta=0 cell contributes a_0 I.
        SPD by construction (eigenvalues a along s^ and a(1-eta/|s|) across), so the Newton tangent
        K_t = INT (dH/dB) curl(dA).curl(w) is positive-definite.  Inputs are arrays of shape (n,) and
        (n, K); returns (nu_xx, nu_xy, nu_yy) each (n,)."""
        Bx = np.atleast_1d(np.asarray(Bx, float)); By = np.atleast_1d(np.asarray(By, float))
        px = np.atleast_2d(px); py = np.atleast_2d(py)
        sx = Bx[:, None] - px; sy = By[:, None] - py
        smag = np.sqrt(sx*sx + sy*sy); smag = np.where(smag < 1e-30, 1e-30, smag)
        moving = (smag > self.eta[None, :]) & (self.eta[None, :] > 0)
        fac = np.where(moving, self.eta[None, :]/smag, 0.0)
        shx = np.where(moving, sx/smag, 0.0); shy = np.where(moving, sy/smag, 0.0)
        a1 = self.a[None, 1:]; m1 = moving[:, 1:]; f1 = fac[:, 1:]; sx1 = shx[:, 1:]; sy1 = shy[:, 1:]
        nu_xx = self.a[0] + np.sum(a1*np.where(m1, 1.0 - f1*(1.0 - sx1**2), 0.0), axis=1)
        nu_yy = self.a[0] + np.sum(a1*np.where(m1, 1.0 - f1*(1.0 - sy1**2), 0.0), axis=1)
        nu_xy = np.sum(a1*np.where(m1, f1*sx1*sy1, 0.0), axis=1)
        return nu_xx, nu_xy, nu_yy

    # --- convex incremental co-energy: the MERIT FUNCTION for the variational (global) Newton ---
    def play_coenergy_density(self, Bx, By, px_prev, py_prev):
        r"""Convex incremental co-energy density psi*(B) whose B-gradient is the model field
        H = a_0 B + sum_{k>=1} a_k p_k, with the FROZEN previous play states (px_prev, py_prev).
        This is the MERIT FUNCTION that makes the FE-coupled-hysteresis Newton GLOBALLY CONVERGENT:
        Pi(A) = INT_iron psi*(B(A)) dx + INT_air (nu0/2)|B|^2 dx - INT J.A dx is CONVEX (psi* is convex
        in B because dH/dB = ``play_tangent`` is SPD, and B = curl A is linear in A), so its unique
        minimiser is the FE solution and a damped Newton with an Armijo line search ON Pi -- not on the
        residual norm -- always finds a decreasing step (the Newton direction d = -K^{-1}R satisfies
        grad(Pi).d = -R^T K^{-1} R < 0).  That cures the kink-overshoot DIVERGENCE of a plain
        residual-norm Newton on the piecewise-linear play H(B): the play model is the exact
        rate-independent-plasticity analogue (p_prev <-> back-stress, eta_k <-> yield radius, the play
        update <-> radial return), and just as in computational plasticity the incremental energy is
        convex and Newton on it with the consistent tangent converges.

        Per dissipative cell k (q = p_prev_k, r = B - q):
            pinned  (|r| <= eta_k):  psi_k = a_k ( q.B - |q|^2 / 2 )            [affine in B]
            slipping(|r| >  eta_k):  psi_k = a_k ( |B|^2/2 - eta_k|r| + eta_k^2/2 )
        (C^1 across |r|=eta_k; the eta_k^2/2 constant enforces value-continuity).  The eta=0 reversible
        cell falls out of the slip branch as a_0|B|^2/2 with no special case.  Linear cells only.
        Inputs (n,) and (n, K); returns psi* (n,)."""
        if self.a is None:
            raise ValueError("co-energy implemented for linear cells (slopes `a`)")
        Bx = np.atleast_1d(np.asarray(Bx, float)); By = np.atleast_1d(np.asarray(By, float))
        qx = np.atleast_2d(px_prev); qy = np.atleast_2d(py_prev)
        rx = Bx[:, None] - qx; ry = By[:, None] - qy
        rmag = np.sqrt(rx * rx + ry * ry)
        eta = self.eta[None, :]
        psi_pin = (qx * Bx[:, None] + qy * By[:, None]) - 0.5 * (qx * qx + qy * qy)
        psi_slip = 0.5 * (Bx * Bx + By * By)[:, None] - eta * rmag + 0.5 * eta * eta
        psi_cell = self.a[None, :] * np.where(rmag <= eta, psi_pin, psi_slip)
        return np.sum(psi_cell, axis=1)

    # --- arbitrary (vector) B-waveform loss: the per-point MOTOR iron-loss kernel ---
    def loss_from_waveform(self, Bx_wave, By_wave=None, settle=2):
        """Hysteresis loss density [J/m^3 per cycle] = |∮ H . dB| for an arbitrary (vector) B(t)
        waveform over one closed cycle (the waveform's first and last samples should coincide).  The
        operator is pre-cycled ``settle`` times to reach the steady minor/major loop.  This is the
        per-point motor iron-loss kernel: feed the per-element B(theta) from a rotor-angle sweep (the
        committed cogging / 3-phase AGE sweeps) and you get the REAL loop-area hysteresis loss --
        minor loops handled natively (every reversal is tracked), no fitted Steinmetz k_h."""
        Bx = np.asarray(Bx_wave, dtype=float)
        By = np.zeros_like(Bx) if By_wave is None else np.asarray(By_wave, dtype=float)
        n = len(Bx)
        px = np.zeros(self.K); py = np.zeros(self.K)
        for _ in range(int(settle)):
            for i in range(n):
                _, _, px, py = self.step(Bx[i], By[i], px, py)
        Hx = np.empty(n); Hy = np.empty(n)
        for i in range(n):
            Hx[i], Hy[i], px, py = self.step(Bx[i], By[i], px, py)
        return abs(float(_trapz(Hx, Bx) + _trapz(Hy, By)))

    # --- vector ROTATIONAL hysteresis (rotating B of constant magnitude) ---
    def rotational_loop(self, Bm, n=721, cycles=3):
        """Steady rotational B-H: B = Bm (cos phi, sin phi) over a full revolution (after `cycles` to
        settle).  Returns (Bx, By, Hx, Hy) over the final revolution."""
        phi = np.linspace(0.0, 2 * np.pi, n)
        bx, by = Bm * np.cos(phi), Bm * np.sin(phi)
        px = np.zeros(self.K); py = np.zeros(self.K)
        Hx = np.empty(n); Hy = np.empty(n)
        for c in range(cycles):
            for i in range(n - 1):
                _, _, px, py = self.step(bx[i], by[i], px, py)
        for i in range(n):
            Hx[i], Hy[i], px, py = self.step(bx[i], by[i], px, py)
        return bx, by, Hx, Hy

    def rotational_loss_per_cycle(self, Bm, n=721, cycles=3):
        """Rotational hysteresis loss density [J/m^3 per revolution] = ∮ H . dB for a rotating B."""
        bx, by, Hx, Hy = self.rotational_loop(Bm, n=n, cycles=cycles)
        return abs(float(_trapz(Hx, bx) + _trapz(Hy, by)))

    def analytic_rotational_loss(self, Bm):
        """Closed-form rotational loss for linear cells:  2*pi * sum_k a_k eta_k * sqrt(Bm^2 - eta_k^2).
        Under a rotating B of radius Bm, the play state p_k traces a circle of radius sqrt(Bm^2-eta_k^2)
        lagging B by sin(delta_k)=eta_k/Bm (|B-p_k|=eta_k), and the work ∮ H.dB = 2*pi*a_k*r_k*Bm*
        sin(delta_k) = 2*pi*a_k*eta_k*sqrt(Bm^2-eta_k^2) per cell.  In the LOW-FIELD limit (eta_k << Bm)
        this -> 2*pi*a_k*eta_k*Bm, so the rotational/alternating loss ratio -> 2*pi/4 = pi/2 (the
        classical Rayleigh-regime result).  At high field each cell's rotational loss falls back toward
        zero as eta_k -> Bm (sqrt -> 0) -- the seed of the rotational-loss peak-and-collapse that a
        SATURATING shape function produces."""
        if self.a is None:
            raise ValueError("analytic rotational loss only for linear cells")
        return float(2.0 * np.pi * np.sum([self.a[k] * self.eta[k] * np.sqrt(Bm ** 2 - self.eta[k] ** 2)
                                           for k in range(self.K) if 0.0 < self.eta[k] < Bm]))

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

    # --- FORC (first-order reversal curves): the Preisach-identification technique (HysterSoft) ---
    def forc_curves(self, Bsat, n=121):
        r"""First-Order Reversal Curves in B-space.  Saturate to +Bsat, reverse DOWN to each reversal
        field Ba on a uniform grid, then sweep B back UP recording H(Ba, Bb).  Returns ``(grid, H)``
        with ``grid = linspace(-Bsat, Bsat, n)`` and ``H[i, j]`` = H at reversal ``Ba=grid[i]``, field
        ``Bb=grid[j]`` (``NaN`` for ``Bb < Ba``).

        FORC (Mayergoyz; Pike-Roberts-Verosub 1999) is the standard experimental identification of the
        Preisach distribution, and the signature feature of HysterSoft (Dimian-Andrei, FAMU-FSU / TU
        Vienna / Cuza University; "Scalar and vector hysteresis simulations using HysterSoft").  The
        FORC distribution :func:`forc_distribution` ``rho = -1/2 d2H/dBa dBb`` of the B-input Play model
        is a set of RIDGES at the coercivity ``(Bb-Ba)/2 = eta_k`` with integrated weight
        ``a_k (Bsat - eta_k)`` (:meth:`analytic_forc_weights`) -- i.e. the play thresholds ARE the
        B-space Preisach / FORC density (Play == static Preisach, Bobbio 1997)."""
        grid = np.linspace(-Bsat, Bsat, n)
        H = np.full((n, n), np.nan)
        for i in range(n):
            px = np.zeros(self.K); py = np.zeros(self.K)
            _, _, px, py = self.step(float(Bsat), 0.0, px, py)        # positive saturation
            _, _, px, py = self.step(float(grid[i]), 0.0, px, py)     # reverse down to Ba
            for j in range(i, n):                                     # sweep up: Bb >= Ba
                Hx, _, px, py = self.step(float(grid[j]), 0.0, px, py)
                H[i, j] = Hx
        return grid, H

    def analytic_forc_weights(self, Bsat):
        r"""Closed-form FORC-distribution ridge weights for the linear-cell model: ``rho`` is
        concentrated on the coercivity lines ``(Bb-Ba)/2 = eta_k`` and the integrated weight of ridge
        ``k`` is ``a_k (Bsat - eta_k)``.  Returns ``[(eta_k, weight_k), ...]`` over the dissipative
        cells (``0 < eta_k < Bsat``)."""
        if self.a is None:
            raise ValueError("analytic FORC weights only for linear cells")
        return [(float(self.eta[k]), float(self.a[k] * (Bsat - self.eta[k])))
                for k in range(self.K) if 0.0 < self.eta[k] < Bsat]


def play_cells(Bmax, N, reluctivities):
    """Standard play discretisation (Matsuo / Hane-Sugahara IEEE TMAG 2026, Eq. 3): equal-interval
    B-space thresholds  eta_k = (k - 1/2) * Bmax / N,  k = 1..N, plus the eta=0 reversible cell.
    ``reluctivities`` is the per-cell linear slope a_k (a reluctivity, for B-input); pass a scalar for
    equal cells or an array of length N+1.  These equal intervals are what the Everett-function shape
    identification assumes."""
    eta = np.concatenate([[0.0], (np.arange(1, N + 1) - 0.5) * Bmax / N])
    a = np.full(N + 1, reluctivities) if np.isscalar(reluctivities) else np.asarray(reluctivities, float)
    return PlayHysteresis(eta=eta, a=a)


def saturating_cells(Bmax, N, a_each, B_sat):
    """Standard equal-interval thresholds but SATURATING shape functions f_k(r) = a*B_sat*tanh(r/B_sat),
    so H saturates at high B.  The ROTATIONAL loss then SATURATES to the closed-form high-field limit
    2*pi * a * B_sat * sum_k eta_k  (each saturated cell contributes 2*pi*a*B_sat*Bm*sin(delta_k) =
    2*pi*a*B_sat*eta_k, independent of Bm), instead of growing without bound like the linear-cell model
    -- a physical, BOUNDED rotational loss.  (Nonlinear cells: the analytic linear-cell loop-area /
    rotational formulas do not apply.)"""
    eta = np.concatenate([[0.0], (np.arange(1, N + 1) - 0.5) * Bmax / N])
    f = [(lambda r, A=a_each, S=B_sat: A * S * np.tanh(r / S)) for _ in range(N + 1)]
    return PlayHysteresis(eta=eta, f_k=f)


def identify_from_loop_areas(Bmax, N, loop_areas, a0=0.0):
    """IDENTIFY the dissipative cell slopes a_k from MEASURED symmetric-loop areas W(m*Delta), m=1..N
    (Delta=Bmax/N, equal-interval thresholds eta_k=(k-1/2)*Delta).  For the linear-cell B-input Play,
    the loop area is the lower-triangular relation

        W(m*Delta) = 4 Delta^2 * sum_{k=1..m} a_k (k-1/2)(m-k+1/2),

    so the cell slopes are recovered by inverting it (a deterministic identification, the linear-cell
    analogue of the Everett-function shape-function identification used for general piecewise-linear
    cells -- see radia.hysteresis_io / the lab Hane-Sugahara TMAG 2026 paper).  ``a0`` is the reversible
    reluctivity (the loop-TIP slope, which the area does not constrain).  Returns a fitted
    :class:`PlayHysteresis`."""
    D = Bmax / N
    M = np.zeros((N, N))
    for mi in range(1, N + 1):
        for k in range(1, mi + 1):
            M[mi - 1, k - 1] = 4.0 * D * D * (k - 0.5) * (mi - k + 0.5)
    a = np.linalg.solve(M, np.asarray(loop_areas, dtype=float))
    eta = np.concatenate([[0.0], (np.arange(1, N + 1) - 0.5) * D])
    return PlayHysteresis(eta=eta, a=np.concatenate([[a0], a]))


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


def forc_distribution(grid, H):
    r"""FORC distribution ``rho(Ba, Bb) = -1/2 d2H / dBa dBb`` from a FORC family ``(grid, H)``
    (:meth:`PlayHysteresis.forc_curves`) -- the central mixed second difference, ``NaN`` wherever the
    3x3 stencil leaves the measured ``Bb >= Ba`` region.  For the linear-cell Play model this is a set
    of ridges at the coercivity ``(Bb-Ba)/2 = eta_k`` (the B-space Preisach density); the per-ridge
    integrated weight ``INT rho dBa dBb`` over a band around ``eta_k`` equals ``a_k (Bsat - eta_k)``
    (:meth:`PlayHysteresis.analytic_forc_weights`).  This is the Play-model side of the FORC /
    Preisach identification that HysterSoft (Dimian-Andrei) performs on measured loops."""
    grid = np.asarray(grid, float)
    n = len(grid)
    d = grid[1] - grid[0]
    rho = np.full((n, n), np.nan)
    for i in range(1, n - 1):
        for j in range(1, n - 1):
            if np.isnan(H[i - 1:i + 2, j - 1:j + 2]).any():
                continue
            rho[i, j] = -0.5 * (H[i + 1, j + 1] - H[i + 1, j - 1]
                                - H[i - 1, j + 1] + H[i - 1, j - 1]) / (4.0 * d * d)
    return rho


def forc_coercivity_weight(grid, rho, Bc, half_band):
    r"""Integrate the FORC distribution ``rho`` over the band ``|(Bb-Ba)/2 - Bc| < half_band`` --
    ``INT_INT rho dBa dBb`` near coercivity ``Bc``.  For the linear-cell Play model this recovers the
    ridge weight ``a_k (Bsat - eta_k)`` when ``Bc = eta_k`` (verified to ~1-3% on a fine grid)."""
    grid = np.asarray(grid, float)
    n = len(grid)
    d = grid[1] - grid[0]
    w = 0.0
    for i in range(n):
        for j in range(n):
            if np.isnan(rho[i, j]):
                continue
            if abs((grid[j] - grid[i]) / 2.0 - Bc) < half_band:
                w += rho[i, j] * d * d
    return w


def identify_from_forc(grid, H, etas, Bsat, a0=0.0, half_band=None):
    r"""INVERSE FORC identification (the HysterSoft workflow): recover the linear-cell B-input Play
    slopes ``a_k`` from a measured FORC family ``(grid, H)`` (:meth:`PlayHysteresis.forc_curves`) by
    reading the RIDGE COMB of its FORC distribution.  For each threshold ``eta_k`` (the known
    equal-interval discretisation), the ridge weight ``INT rho dBa dBb`` over the coercivity band
    ``|(Bb-Ba)/2 - eta_k| < half_band`` equals ``a_k (Bsat - eta_k)``, so
    ``a_k = weight / (Bsat - eta_k)``.  The reversible slope ``a0`` (``eta = 0`` -> no irreversible
    switching -> no FORC ridge) is a separate input (the loop-tip slope).  ``half_band`` defaults to
    0.4 x the smallest threshold spacing.  Returns a fitted :class:`PlayHysteresis` (round-trips a known
    model to ~1-3 %, tighter on a finer ``grid``).

    This is the FORC-based identification HysterSoft (Dimian-Andrei) performs on measured loops, applied
    to a Play material -- the reversal-curve counterpart of the symmetric-loop-area
    :func:`identify_from_loop_areas` (FORC resolves the full Preisach density, not just the diagonal)."""
    etas = np.asarray(etas, dtype=float)
    pos = etas[etas > 0.0]
    if half_band is None:
        knots = np.sort(np.concatenate([[0.0], pos]))
        half_band = 0.4 * float(np.min(np.diff(knots)))
    rho = forc_distribution(grid, H)
    a = [a0]
    for ek in pos:
        a.append(forc_coercivity_weight(grid, rho, float(ek), half_band) / (Bsat - ek))
    return PlayHysteresis(eta=np.concatenate([[0.0], pos]), a=np.array(a, dtype=float))
