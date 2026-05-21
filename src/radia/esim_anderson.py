"""Anderson-type-II acceleration for the outer Karl iteration.

The damped-Picard "Karl" loop

    x_{k+1} = x_k + alpha * (G(x_k) - x_k)

contracts only when `alpha * L < 1` where `L` is the local Lipschitz
constant of `G`.  For per-element ESIM on a BH-knee-straddling
workpiece, the worst-DOF Lipschitz is ~2-3 and `alpha = 0.3-0.5`
gives a per-DOF `dZ_max` noise floor of ~0.06-0.20 that does NOT
fall below the strict `esim_tol = 1e-3` criterion within reasonable
iteration counts.

Anderson acceleration combines the last `m` iterates via a least-
squares solve to eliminate the slowly-decaying modes that the
damped-Picard cannot suppress.  Reference:

    H. F. Walker and P. Ni, "Anderson Acceleration for Fixed-Point
    Iterations," SIAM J. Numer. Anal. 49 (4), 1715-1735, 2011.

Use this module's :class:`AndersonAccelerator` from any Karl loop
(BIE / FEM-Kelvin / FEM-coilmesh) by replacing

    Z_s_wp = alpha * Z_s_new + (1 - alpha) * Z_s_old

with

    Z_s_wp = anderson.step(Z_s_old, Z_s_new)

where ``Z_s_new = G(Z_s_old)`` is the cell-solver output.

The accelerator works for BOTH scalar and per-DOF Karl: the input
is treated uniformly as a complex 1-D array internally; the public
``step`` method preserves the input shape (scalar -> scalar,
ndarray -> ndarray).
"""
from __future__ import annotations

from collections import deque
from typing import Deque

import numpy as np


class AndersonAccelerator:
    """Anderson-type-II accelerator for a fixed-point map G(x).

    Parameters
    ----------
    m : int
        Memory depth (max number of past iterates retained).  ``m = 0``
        falls back to plain damped Picard (no acceleration).  Common
        choice: 3 - 5 for moderate Lipschitz problems.
    alpha : float
        Damping factor in (0, 1].  ``1.0`` recovers undamped Anderson;
        the lab convention is ``alpha = 0.5 - 0.7`` even when
        acceleration is on, to stabilise the first few iterations.
    """

    def __init__(self, m: int = 5, alpha: float = 0.7,
                 step_clip: float = 2.0,
                 restart_growth: float = 2.0):
        """
        Parameters (continued)
        ----------
        step_clip : float
            Bound on the Anderson correction magnitude **relative** to
            the damped-Picard step.  After computing ``corr = -(dX +
            alpha * dF) @ gamma``, if ``||corr|| > step_clip *
            ||alpha * f_k||``, scale ``corr`` down so the inequality
            holds.  Prevents pathological over-extrapolation when the
            LSQ matrix is ill-conditioned.  Default ``2.0`` allows
            Anderson to make at most 2x the damped-Picard move per
            step.  Set to ``float('inf')`` to disable clipping.
        restart_growth : float
            Residual-monotone restart threshold (Walker-Ni 2011 § 3.4).
            If ``||f_k|| > restart_growth * min(||f_j||, j < k)``, the
            most recent step diverged: clear the history (keep current
            iterate) so the next step is plain damped Picard and
            Anderson rebuilds from scratch.  Default ``2.0`` means
            "restart when residual grew by more than 2x best-so-far".
            Set to ``float('inf')`` to disable restart.
        """
        if m < 0:
            raise ValueError(f"m must be >= 0, got {m}")
        if not (0.0 < alpha <= 1.0):
            raise ValueError(f"alpha must be in (0, 1], got {alpha}")
        if step_clip <= 0:
            raise ValueError(f"step_clip must be > 0, got {step_clip}")
        if restart_growth <= 1.0:
            raise ValueError(
                f"restart_growth must be > 1, got {restart_growth}")
        self.m = int(m)
        self.alpha = float(alpha)
        self.step_clip = float(step_clip)
        self.restart_growth = float(restart_growth)
        self._X: Deque[np.ndarray] = deque(maxlen=m + 1)
        self._F: Deque[np.ndarray] = deque(maxlen=m + 1)
        # Track the previous step's residual norm for restart criterion.
        self._f_norm_prev: float = 0.0
        # Diagnostic counters (exposed via properties below).
        self._n_restarts: int = 0
        self._n_clips: int = 0

    @property
    def n_restarts(self) -> int:
        """Number of residual-monotone restarts triggered so far."""
        return self._n_restarts

    @property
    def n_clips(self) -> int:
        """Number of Anderson steps clipped by ``step_clip``."""
        return self._n_clips

    def reset(self) -> None:
        """Clear iterate / residual history.

        Call this when the underlying map has changed discontinuously
        (e.g. mesh refined, BH curve swapped), or after a restart
        criterion fires.
        """
        self._X.clear()
        self._F.clear()

    def step(self, x_k, g_k):
        """Return the next iterate ``x_{k+1}`` given ``x_k`` and ``G(x_k)``.

        Parameters
        ----------
        x_k : complex or complex ndarray
            Current iterate.
        g_k : complex or complex ndarray
            Fixed-point map evaluated at ``x_k``.  Same shape as ``x_k``.

        Returns
        -------
        x_next : same type as ``x_k``
            The accelerated iterate.  For ``m = 0`` or when no history
            is available yet, this is just the damped-Picard update
            ``x_k + alpha * (g_k - x_k)``.
        """
        scalar_in = np.ndim(x_k) == 0
        x_arr = np.atleast_1d(x_k).astype(complex, copy=True)
        g_arr = np.atleast_1d(g_k).astype(complex, copy=True)
        if x_arr.shape != g_arr.shape:
            raise ValueError(
                f"x_k shape {x_arr.shape} != g_k shape {g_arr.shape}")
        f_arr = g_arr - x_arr

        # Residual-monotone restart (Walker-Ni 2011 § 3.4).  Compare
        # the current RELATIVE inf-norm residual against the
        # IMMEDIATELY PREVIOUS one (not the historical minimum -- per-
        # element ESIM has an intrinsic noise floor that the iteration
        # never falls below, so min-tracking triggers continuous
        # restarts and kills the Anderson history).  Use the inf-norm
        # ratio because per-element ESIM has spatially varying |Z_s|;
        # the 2-norm is dominated by cool-spot DOFs and misses the
        # hot-spot oscillation that matters.
        x_inf = float(np.max(np.abs(x_arr)))
        f_inf = float(np.max(np.abs(f_arr)))
        rel_resid = f_inf / max(x_inf, 1e-30)
        if (rel_resid > self.restart_growth * self._f_norm_prev
                and self._f_norm_prev > 0.0
                and len(self._X) > 0):
            self._X.clear()
            self._F.clear()
            self._n_restarts += 1
        # Update previous-residual tracker.
        self._f_norm_prev = rel_resid

        self._X.append(x_arr.copy())
        self._F.append(f_arr.copy())

        m_k = len(self._X) - 1
        picard_step = self.alpha * f_arr
        if self.m == 0 or m_k == 0:
            x_next = x_arr + picard_step
        else:
            X_list = list(self._X)
            F_list = list(self._F)
            dX = np.column_stack(
                [X_list[i + 1] - X_list[i] for i in range(m_k)])
            dF = np.column_stack(
                [F_list[i + 1] - F_list[i] for i in range(m_k)])
            gamma = self._solve_real_lsq(dF, f_arr)
            anderson_corr = -(dX + self.alpha * dF) @ gamma
            # Step-clip safeguard: limit the Anderson correction to a
            # bounded multiple of the damped-Picard step.
            corr_norm = float(np.linalg.norm(anderson_corr))
            picard_norm = float(np.linalg.norm(picard_step))
            max_corr = self.step_clip * max(picard_norm, 1e-30)
            if corr_norm > max_corr:
                anderson_corr = anderson_corr * (max_corr / corr_norm)
                self._n_clips += 1
            x_next = x_arr + picard_step + anderson_corr

        if scalar_in:
            return complex(x_next.item())
        return x_next

    @staticmethod
    def _solve_real_lsq(dF: np.ndarray, f_k: np.ndarray) -> np.ndarray:
        """Solve ``min_gamma || dF * gamma - f_k ||^2`` over REAL gamma.

        Because Anderson weights are conventionally real-valued for
        complex fixed-point iterations (Walker-Ni 2011 §2.2 remark),
        we stack the real and imaginary rows of ``dF`` and ``f_k`` to
        convert the complex LSQ into a purely real one.

        Numerical conditioning is delegated to ``np.linalg.lstsq``'s
        SVD-based rcond cutoff (default).
        """
        if np.iscomplexobj(dF) or np.iscomplexobj(f_k):
            A = np.vstack([dF.real, dF.imag])
            b = np.concatenate([f_k.real, f_k.imag])
        else:
            A, b = dF, f_k
        gamma, *_ = np.linalg.lstsq(A, b, rcond=None)
        return gamma
