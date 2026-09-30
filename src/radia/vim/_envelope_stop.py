"""Coefficient-driven energy-envelope Stop material core.

This module contains only the generic numerical model.  Applications must
provide a finite nonnegative coefficient vector and its supported flux-density
range explicitly; no material measurements or fitted defaults are bundled.

The scalar (``dim=1``) path is the supported production mode.  Multidimensional
use is available only through an explicit ``allow_unvalidated_rotation=True``
opt-in because the isotropic vector lift has structural tests but no rotational
material-accuracy validation.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import BSpline

__all__ = [
    "EnvelopeStopVector3D",
    "envelope_stop_scalar_H",
    "envelope_stop_full_design",
]

MU0 = 4.0e-7 * np.pi

# ---------------------------------------------------------------------------
# I-spline shape bases
# ---------------------------------------------------------------------------


def build_ispline_knots(degree=3, n_internal=8):
    """Clamped knot vector on [0, 1] and the number of I-spline bases."""
    knots_int = np.linspace(0.0, 1.0, n_internal + 2)[1:-1]
    t = np.concatenate([np.zeros(degree + 1), knots_int, np.ones(degree + 1)])
    return t, len(t) - degree - 1


def eval_ispline(z, t_knots, n_basis, degree=3):
    """I-spline basis values ``I_j(z) = int_0^z B_j``; shape ``(len(z), n_basis)``."""
    z = np.clip(np.asarray(z, dtype=float), 0.0, 1.0)
    out = np.empty((z.size, n_basis))
    for j in range(n_basis):
        c = np.zeros(n_basis)
        c[j] = 1.0
        anti = BSpline(t_knots, c, degree, extrapolate=False).antiderivative()
        out[:, j] = anti(z) - anti(0.0)
    return out


# ---------------------------------------------------------------------------
# Model configuration (must stay identical to the reference parameterization)
# ---------------------------------------------------------------------------

FIXED_LINEAR = 0.0          # no separately fixed alpha*B background
ALPHA = FIXED_LINEAR
MARGIN = 2.0                # parameterization tangent lower bound dH/dB >= MARGIN (A/m/T)
Bs = 1.6
Q = 10
M_RAW = 32                  # clamped cubic B-splines on [-Bs, Bs]
M = M_RAW // 2              # 16 even-folded pairs -> H exactly odd
K = 150
eta = np.linspace(0.01, 1.5, K)
_en = (eta - eta.min()) / (eta.max() - eta.min())
_ec = np.linspace(0, 1, Q)
_ew = 1.0 / Q * 1.3
Rband = np.exp(-((_en[None, :] - _ec[:, None]) / _ew) ** 2)
wq = np.gradient(eta)
wq = wq / wq.mean()
degree, n_internal = 4, 8
t_knots, P = build_ispline_knots(degree=degree, n_internal=n_internal)
zg = np.linspace(0, 1, 300)
Ig = eval_ispline(np.minimum(zg, 1 - 1e-9), t_knots, P, degree=degree)
IIg = cumulative_trapezoid(Ig, zg, axis=0, initial=0.0)
Pqm = P * Q * M
DSTATIC = P
DPAR = Pqm + DSTATIC
bc = np.linspace(-Bs, Bs, M_RAW)
_bspl_t = np.concatenate([[bc[0]] * 3, np.linspace(bc[0], bc[-1], M_RAW - 2), [bc[-1]] * 3])
_bspl = [BSpline(_bspl_t, np.eye(M_RAW)[m], 3) for m in range(M_RAW)]
_bspl_d = [s.derivative() for s in _bspl]


def psi(x):
    """Even-folded envelope basis ``psi~_m = psi_m + psi_{M_RAW-1-m}`` and its
    derivative, each of shape (M, N)."""
    xc = np.clip(x, bc[0], bc[-1])
    ps = np.stack([s(xc) for s in _bspl])
    dp = np.stack([s(xc) for s in _bspl_d])
    return ps[:M] + ps[:M - 1 - M_RAW:-1], dp[:M] + dp[:M - 1 - M_RAW:-1]


def _psi_second(x):
    """Second derivative of the even-folded envelope basis, shape (M, N)."""
    xc = np.clip(np.atleast_1d(np.asarray(x, float)), bc[0], bc[-1])
    dd = np.stack([s.derivative(2)(xc) for s in _bspl])
    return dd[:M] + dd[:M - 1 - M_RAW:-1]


# ---------------------------------------------------------------------------
# Scalar design (the identification-side definition of the same model)
# ---------------------------------------------------------------------------


def stop_trajectory(Bw):
    """Scalar Stop states ``s_k(t)`` for the flux-density waveform ``Bw``."""
    Bw = np.asarray(Bw, float)
    s = np.empty((len(Bw), K))
    s[0] = np.clip(Bw[0], -eta, eta)
    for i in range(1, len(Bw)):
        s[i] = np.clip(s[i - 1] + (Bw[i] - Bw[i - 1]), -eta, eta)
    return s


def _SR_JR(Bc, sc):
    N = len(Bc)
    zk = np.clip(sc / eta[None, :], -1, 1)
    az = np.abs(zk)
    sg = np.sign(sc)
    sg[sc == 0] = 0
    az_flat = az.reshape(-1)
    Ik = eval_ispline(np.minimum(az_flat, 1 - 1e-9), t_knots, P, degree=degree).reshape(N, K, P)
    Sp = (sg[:, :, None] * Ik).transpose(0, 2, 1)
    Jp = np.empty((N, P, K))
    for pp in range(P):
        Jp[:, pp, :] = eta[None, :] * np.interp(az_flat, zg, IIg[:, pp]).reshape(N, K)
    SR = np.einsum('npk,qk,k->npq', Sp, Rband, wq)
    JR = np.einsum('npk,qk,k->npq', Jp, Rband, wq)
    return SR, JR


def design_hysteretic(Bw, chunk=2400):
    """(N, 2080) history-dependent design matrix."""
    Bw = np.asarray(Bw, float)
    outs = []
    s_all = stop_trajectory(Bw)
    for i0 in range(0, len(Bw), chunk):
        sl = slice(i0, min(i0 + chunk, len(Bw)))
        Bc = Bw[sl]
        SR, JR = _SR_JR(Bc, s_all[sl])
        ps, dps = psi(Bc)
        outs.append((SR[:, :, :, None] * ps.T[:, None, None, :]
                     + JR[:, :, :, None] * dps.T[:, None, None, :]).reshape(len(Bc), Pqm))
    return np.vstack(outs)


def state_independent_design(Bw):
    """(N, 13) odd B-only field bases with even nonnegative tangents."""
    Bw = np.asarray(Bw, dtype=float)
    z = np.minimum(np.clip(np.abs(Bw) / Bs, 0.0, 1.0), 1.0 - 1e-9)
    return np.sign(Bw)[:, None] * eval_ispline(z, t_knots, P, degree=degree)


def envelope_stop_full_design(Bw, chunk=2400):
    """(N, 2093) total-energy design: history block followed by B-only block."""
    return np.hstack([design_hysteretic(Bw, chunk=chunk), state_independent_design(Bw)])


def envelope_stop_scalar_H(Bw, theta):
    """Scalar response ``H(t)`` for an explicit coefficient vector."""
    return envelope_stop_full_design(Bw) @ np.ravel(np.asarray(theta, float))


# ---------------------------------------------------------------------------
# Cached evaluators used by the vector material
# ---------------------------------------------------------------------------


def _build_ispline_evaluators():
    """I-spline, integrated I-spline and I-spline derivative evaluators.

    The first is bit-checked against :func:`eval_ispline`; the second is its
    exact antiderivative (B-only stored energy); the third its exact
    derivative (analytic tangent).
    """
    antis = []
    for j in range(P):
        c = np.zeros(P)
        c[j] = 1.0
        antis.append(BSpline(t_knots, c, degree, extrapolate=False).antiderivative())
    C = np.stack([a.c[:len(antis[0].c)] for a in antis], axis=1)
    S = BSpline(antis[0].t, C, antis[0].k, extrapolate=False)
    z0 = S(np.array([0.0]))[0]
    SS = S.antiderivative()
    zz0 = SS(np.array([0.0]))[0]
    DS = S.derivative()

    def eval_cached(z):
        z = np.atleast_1d(np.clip(np.asarray(z, float), 0.0, 1.0))
        return S(z) - z0[None, :]

    def eval_integrated(z):
        z = np.atleast_1d(np.clip(np.asarray(z, float), 0.0, 1.0))
        return SS(z) - zz0[None, :]

    def eval_derivative(z):
        z = np.atleast_1d(np.clip(np.asarray(z, float), 0.0, 1.0))
        return DS(z)

    zp = np.linspace(0.0, 1.0 - 1e-9, 257)
    err = float(np.max(np.abs(eval_cached(zp) - eval_ispline(zp, t_knots, P, degree=degree))))
    if err >= 1e-13:
        raise RuntimeError("cached I-spline mismatch: %.3e" % err)
    return eval_cached, eval_integrated, eval_derivative


_EVAL_I, _EVAL_II, _EVAL_DI = _build_ispline_evaluators()
# Piecewise-constant slopes of the II interpolant: the coefficients were
# defined with np.interp'd II columns, so the analytic tangent
# differentiates that interpolant.
_IIG_SLOPES = np.diff(IIg, axis=0) / np.diff(zg)[:, None]


# ---------------------------------------------------------------------------
# Vector material (radia.vim.SolveHysteresis duck-typed protocol)
# ---------------------------------------------------------------------------


class EnvelopeStopVector3D:
    """Isotropic vector lift of the scalar C4-even energy envelope Stop.

    Protocol: ``state0()``, ``forward(B, states)`` (pure), ``commit(B,
    states)``, ``nu_bound()``, ``b_max()``; plus ``energy``, ``tangent``
    (finite difference, audit), ``tangent_analytic`` and ``inverse``.
    State layout per point: ``[s_1..s_K (K*dim), B_last (dim)]``.

    ``theta`` is always explicit.  ``b_max_T`` records the largest flux
    density covered by the caller's coefficient validation and is consumed by
    :func:`radia.vim.SolveHysteresis`'s fail-loud range gate.
    """

    permanent_magnet_model = "b-input-envelope-stop-c4even"
    permanent_magnet_level = 4

    def __init__(self, theta, *, b_max_T, dim=1,
                 allow_unvalidated_rotation=False):
        self.dim = int(dim)
        if self.dim < 1:
            raise ValueError("EnvelopeStopVector3D: dim must be positive")
        if self.dim > 1 and not allow_unvalidated_rotation:
            raise ValueError(
                "EnvelopeStopVector3D: multidimensional/rotational use is not "
                "validated; pass allow_unvalidated_rotation=True only for "
                "explicit experimental or structural studies")
        self.rotational_validation = (
            "unvalidated-explicit-opt-in" if self.dim > 1 else "not-applicable")
        self._b_max_T = float(b_max_T)
        if not np.isfinite(self._b_max_T) or self._b_max_T <= 0.0:
            raise ValueError("EnvelopeStopVector3D: b_max_T must be finite and positive")
        self.coefficient_set = "user"
        th = np.asarray(theta, float)
        th = np.ravel(th)
        if th.size == Pqm:
            history = th
            beta = np.zeros(P)
            self.coefficient_schema = "legacy-history-only-2080"
        elif th.size == DPAR:
            history = th[:Pqm]
            beta = th[Pqm:]
            self.coefficient_schema = "total-energy-2093"
        else:
            raise ValueError(
                "EnvelopeStopVector3D: theta must contain %d legacy history "
                "coefficients or %d total-energy coefficients, got %d" % (Pqm, DPAR, th.size))
        if not np.all(np.isfinite(th)):
            raise ValueError("EnvelopeStopVector3D: theta must be finite")
        if np.any(th < -1.0e-12):
            raise ValueError("EnvelopeStopVector3D: negative coefficients violate the "
                             "passive energy parameterization")
        T3 = history.reshape(P, Q, M)
        self._Wk = np.einsum('pqm,qk,k->kpm', T3, Rband, wq)
        self._beta = beta.copy()
        self._th = th.copy()

    # -- protocol -------------------------------------------------------------
    def state0(self):
        return np.zeros(K * self.dim + self.dim)

    @property
    def state_size(self):
        return K * self.dim + self.dim

    def _advance(self, B, states):
        B = np.asarray(B, float)
        states = np.asarray(states, float)
        if B.ndim != 2 or B.shape[1] != self.dim:
            raise ValueError("EnvelopeStopVector3D: B must have shape (n,%d)" % self.dim)
        if states.ndim != 2 or states.shape != (B.shape[0], self.state_size):
            raise ValueError("EnvelopeStopVector3D: states must have shape (n,%d)"
                             % self.state_size)
        if not np.all(np.isfinite(B)) or not np.all(np.isfinite(states)):
            raise ValueError("EnvelopeStopVector3D: B and states must be finite")
        n = states.shape[0]
        s = states[:, :K * self.dim].reshape(n, K, self.dim)
        B_last = states[:, K * self.dim:]
        s2 = s + (B - B_last)[:, None, :]
        r = np.sqrt((s2 * s2).sum(2))
        scale = np.where(r > eta[None, :], eta[None, :] / np.maximum(r, 1e-300), 1.0)
        s2 = s2 * scale[:, :, None]
        r = np.minimum(r, eta[None, :])
        return s2, r

    def _state_independent_terms(self, Bn):
        """Radial B-only field magnitude and its stored energy."""
        Bn = np.asarray(Bn, float)
        z_limit = 1.0 - 1.0e-9
        z = np.minimum(np.maximum(Bn / Bs, 0.0), z_limit)
        h0 = _EVAL_I(z) @ self._beta
        energy0 = Bs * (_EVAL_II(z) @ self._beta)
        B_cap = Bs * z_limit
        beyond = Bn > B_cap
        if np.any(beyond):
            energy0[beyond] += (Bn[beyond] - B_cap) * h0[beyond]
        return h0, energy0

    def _H_of(self, B, s2, r):
        B = np.asarray(B, float)
        n = r.shape[0]
        z2 = r / eta[None, :]
        z = np.minimum(z2, 1.0 - 1e-9)
        Ik = _EVAL_I(z.ravel()).reshape(n, K, P)
        IIk = np.empty((n, K, P))
        zf = z2.ravel()
        for pp in range(P):
            IIk[:, :, pp] = np.interp(zf, zg, IIg[:, pp]).reshape(n, K)
        Bn = np.linalg.norm(B, axis=1)
        ps, dps = psi(Bn)
        h0, _ = self._state_independent_terms(Bn)
        g = np.einsum('nkp,kpm,mn->nk', Ik, self._Wk, ps)
        c = np.einsum('nkp,k,kpm,mn->n', IIk, eta, self._Wk, dps)
        shat = np.zeros_like(s2)
        nz = r > 1e-12
        shat[nz] = s2[nz] / r[nz, None]
        Bhat = np.zeros_like(B)
        bn = Bn > 1e-12
        Bhat[bn] = B[bn] / Bn[bn, None]
        return (g[:, :, None] * shat).sum(1) + (c + h0)[:, None] * Bhat

    def forward(self, B, states):
        s2, r = self._advance(B, states)
        return self._H_of(B, s2, r)

    def commit(self, B, states):
        B = np.asarray(B, float)
        s2, _ = self._advance(B, states)
        return np.concatenate([s2.reshape(s2.shape[0], -1), B], axis=1)

    def energy(self, B, states):
        """Trial stored energy density ``Psi(B, s(B))`` (J/m^3)."""
        B = np.asarray(B, float)
        s2, r = self._advance(B, states)
        n = r.shape[0]
        IIk = np.empty((n, K, P))
        zf = (r / eta[None, :]).ravel()
        for pp in range(P):
            IIk[:, :, pp] = np.interp(zf, zg, IIg[:, pp]).reshape(n, K)
        Bn = np.linalg.norm(B, axis=1)
        ps, _ = psi(Bn)
        U = np.einsum('nkp,k,kpm,mn->n', IIk, eta, self._Wk, ps)
        _, energy0 = self._state_independent_terms(Bn)
        return U + energy0

    def stored_energy(self, B, states):
        return self.energy(B, states)

    def tangent(self, B, states, fd_eps=1.0e-6):
        """Centered finite-difference tangent ``dH/dB`` (independent audit route)."""
        B = np.asarray(B, float)
        states = np.asarray(states, float)
        self._advance(B, states)
        J = np.empty((B.shape[0], self.dim, self.dim))
        for column in range(self.dim):
            Bp = B.copy()
            Bm = B.copy()
            Bp[:, column] += fd_eps
            Bm[:, column] -= fd_eps
            J[:, :, column] = (self.forward(Bp, states) - self.forward(Bm, states)) / (2.0 * fd_eps)
        return J

    def _state_independent_tangent(self, Bn):
        Bn = np.asarray(Bn, float)
        z_limit = 1.0 - 1.0e-9
        z_raw = np.maximum(Bn / Bs, 0.0)
        z = np.minimum(z_raw, z_limit)
        h0d = (_EVAL_DI(z) @ self._beta) / Bs
        h0d[z_raw > z_limit] = 0.0
        return h0d

    def tangent_analytic(self, B, states, r_tol=1.0e-12, q_tol=1.0e-12):
        """Closed-form history-consistent tangent ``dH/dB`` of :meth:`forward`.

        The projected branch loses its Bhat (x) shat block because projection
        freezes the cell radius; this is the analytic origin of the tangent
        skew under rotation.  At a projection threshold the branch selected by
        the state advance is returned.
        """
        B = np.asarray(B, float)
        states = np.asarray(states, float)
        n = B.shape[0]
        dim = self.dim
        s_prev = states[:, :K * dim].reshape(n, K, dim)
        B_last = states[:, K * dim:]
        t = s_prev + (B - B_last)[:, None, :]
        rt = np.sqrt((t * t).sum(2))
        projected = rt > eta[None, :]
        scale = np.where(projected, eta[None, :] / np.maximum(rt, 1e-300), 1.0)
        s = t * scale[:, :, None]
        r = np.minimum(rt, eta[None, :])

        z2 = r / eta[None, :]
        z = np.minimum(z2, 1.0 - 1e-9)
        Ik = _EVAL_I(z.ravel()).reshape(n, K, P)
        dIk = _EVAL_DI(z.ravel()).reshape(n, K, P)
        IIk = np.empty((n, K, P))
        zf = z2.ravel()
        for pp in range(P):
            IIk[:, :, pp] = np.interp(zf, zg, IIg[:, pp]).reshape(n, K)

        Bn = np.linalg.norm(B, axis=1)
        ps, dps = psi(Bn)
        dds = _psi_second(Bn)
        Gp = np.einsum('kpm,mn->nkp', self._Wk, ps)
        Gq = np.einsum('kpm,mn->nkp', self._Wk, dps)
        g = np.einsum('nkp,nkp->nk', Ik, Gp)
        dgdr = np.einsum('nkp,nkp->nk', dIk, Gp) / eta[None, :]
        dgdq = np.einsum('nkp,nkp->nk', Ik, Gq)
        cell = np.clip(np.searchsorted(zg, zf) - 1, 0, len(zg) - 2)
        IIslope = _IIG_SLOPES[cell].reshape(n, K, P)
        cpr = np.einsum('nkp,nkp->nk', IIslope, Gq)
        c = np.einsum('nkp,k,nkp->n', IIk, eta, Gq)
        dcdq = np.einsum('nkp,k,kpm,mn->n', IIk, eta, self._Wk, dds)
        h0, _ = self._state_independent_terms(Bn)
        h0d = self._state_independent_tangent(Bn)

        act = r > r_tol
        shat = np.zeros_like(s)
        shat[act] = s[act] / r[act, None]
        bq = Bn > q_tol
        Bhat = np.zeros_like(B)
        Bhat[bq] = B[bq] / Bn[bq, None]

        free = ~projected
        a1 = np.where(free & act, dgdr, 0.0)
        r_div = np.where(projected, rt, r)
        a2 = np.where(act, g / np.maximum(r_div, 1e-300), 0.0)
        iso = np.where(free & ~act, dgdr, 0.0)
        b1 = np.where(act, dgdq, 0.0)
        b2 = np.where(free & act, cpr, 0.0)

        I_dim = np.eye(dim)
        J = np.einsum('nk,nka,nkb->nab', a1 - a2, shat, shat)
        J += (a2.sum(1) + iso.sum(1))[:, None, None] * I_dim[None]
        v1 = np.einsum('nk,nka->na', b1, shat)
        v2 = np.einsum('nk,nka->na', b2, shat)
        J += v1[:, :, None] * Bhat[:, None, :]
        J += Bhat[:, :, None] * v2[:, None, :]
        radial = dcdq + h0d
        tangential = np.where(bq, (c + h0) / np.maximum(Bn, 1e-300), radial)
        J += tangential[:, None, None] * I_dim[None]
        J += (radial - tangential)[:, None, None] * Bhat[:, :, None] * Bhat[:, None, :]
        return J

    # -- inverse map B(H) ------------------------------------------------------
    def inverse(self, H_target, states, B0=None, tol=1e-8, maxit=150, fd_eps=1e-6,
                allow_extrapolation=False, require_positive_tangent=True):
        """Batched inverse ``B(H)`` (Levenberg-Marquardt on :meth:`forward`) and
        the differential permeability ``inv(dH/dB)`` at the solution.

        Fails loudly on non-convergence, on a result beyond ``b_max`` (unless
        ``allow_extrapolation``) and, by default, on a non-positive symmetric
        tangent.  ``fd_eps`` is accepted for call-site compatibility.
        """
        H_target = np.asarray(H_target, float)
        n = H_target.shape[0]
        B = np.zeros((n, self.dim)) if B0 is None else np.asarray(B0, float).copy()
        Bn0 = np.linalg.norm(B, axis=1)
        warmstart_radius = 0.95 * self._b_max_T
        far = Bn0 > warmstart_radius
        if far.any():
            B[far] *= (warmstart_radius / Bn0[far])[:, None]
        scale = max(1.0, float(np.max(np.linalg.norm(H_target, axis=1))))
        I3 = np.eye(self.dim)
        lam = np.full(n, 10.0)
        H0 = self.forward(B, states)
        rn = np.linalg.norm(H0 - H_target, axis=1) / scale
        for _it in range(maxit):
            if np.all(rn < tol):
                if (not allow_extrapolation) and float(np.max(np.linalg.norm(B, axis=1))) > self._b_max_T:
                    raise RuntimeError("EnvelopeStopVector3D.inverse: converged outside the "
                                       "validated |B| <= %.3g T range" % self._b_max_T)
                J = self.tangent_analytic(B, states)
                Jsym = 0.5 * (J + np.transpose(J, (0, 2, 1)))
                if require_positive_tangent:
                    mineig = np.linalg.eigvalsh(Jsym)[:, 0]
                    if np.any(mineig <= 0.0):
                        worst = int(np.argmin(mineig))
                        raise RuntimeError(
                            "EnvelopeStopVector3D.inverse: non-positive vector material "
                            "tangent at row %d (lambda_min=%.6g A/m/T)" % (worst, mineig[worst]))
                return B, np.linalg.inv(J)
            J = self.tangent_analytic(B, states)
            better = np.zeros(n, dtype=bool)
            for _ in range(30):
                Jl = J + lam[:, None, None] * I3[None, :, :]
                step = np.linalg.solve(Jl, (H0 - H_target)[..., None])[..., 0]
                Bn = B - step
                Hn = self.forward(Bn, states)
                rnn = np.linalg.norm(Hn - H_target, axis=1) / scale
                better = (rnn <= rn) | (rnn < tol)
                if better.all():
                    break
                lam[~better] = np.minimum(lam[~better] * 10.0, 1.0e8)
            B = np.where(better[:, None], Bn, B)
            H0 = np.where(better[:, None], Hn, H0)
            rn = np.where(better, rnn, rn)
            lam = np.where(better, np.maximum(lam * 0.3, 1e-3), lam)
        raise RuntimeError(
            "EnvelopeStopVector3D.inverse: LM-Newton did not converge in %d iterations "
            "(worst rel residual %.2e); H may be outside the validated/invertible range"
            % (maxit, float(rn.max())))

    # -- Hantila bound ---------------------------------------------------------
    def nu_bound(self, safety=1.5, seed=3):
        """Scanned sup of directional ``|dH|/|dB|`` times ``safety``.

        Empirical scan over validation-style drives and random directions, not a
        closed-form certificate; SolveHysteresis's convergence guard remains
        the backstop.
        """
        if getattr(self, "_nu_cache", None) is not None:
            return self._nu_cache * safety / self._nu_safety
        if self.dim != 3:
            raise NotImplementedError("EnvelopeStopVector3D.nu_bound: the scan drives are "
                                      "defined for dim=3 only")
        rng = np.random.default_rng(seed)
        drives = []
        for a in (0.5, 0.9, 1.3, 1.44):
            bx = np.concatenate([np.linspace(0, a, 60), np.linspace(a, -a, 120),
                                 np.linspace(-a, a, 120)])
            D3 = np.zeros((bx.size, 3))
            D3[:, 2] = bx
            drives.append(D3)
        th = np.linspace(0, 4 * np.pi, 240)
        for a in (0.6, 0.9, 1.2):
            R3 = np.stack([a * np.cos(th), np.zeros_like(th), a * np.sin(th)], 1)
            drives.append(np.vstack([np.linspace(0, 1, 40)[:, None] * R3[0], R3]))
        eps = 1e-6
        dirs = np.vstack([np.eye(3), rng.normal(size=(5, 3))])
        dirs /= np.linalg.norm(dirs, axis=1)[:, None]
        worst = ALPHA
        for D in drives:
            st = self.state0()[None, :]
            for nstep in range(len(D)):
                Bc = D[nstep][None, :]
                H0 = self.forward(Bc, st)
                for d in dirs:
                    H1 = self.forward(Bc + eps * d[None, :], st)
                    worst = max(worst, float(np.linalg.norm(H1 - H0) / eps))
                st = self.commit(Bc, st)
        self._nu_cache = worst * safety
        self._nu_safety = safety
        return self._nu_cache

    def b_max(self):
        return self._b_max_T
