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

    def __init__(self, m: int = 5, alpha: float = 0.7):
        if m < 0:
            raise ValueError(f"m must be >= 0, got {m}")
        if not (0.0 < alpha <= 1.0):
            raise ValueError(f"alpha must be in (0, 1], got {alpha}")
        self.m = int(m)
        self.alpha = float(alpha)
        self._X: Deque[np.ndarray] = deque(maxlen=m + 1)
        self._F: Deque[np.ndarray] = deque(maxlen=m + 1)

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

        self._X.append(x_arr.copy())
        self._F.append(f_arr.copy())

        m_k = len(self._X) - 1
        if self.m == 0 or m_k == 0:
            x_next = x_arr + self.alpha * f_arr
        else:
            X_list = list(self._X)
            F_list = list(self._F)
            dX = np.column_stack(
                [X_list[i + 1] - X_list[i] for i in range(m_k)])
            dF = np.column_stack(
                [F_list[i + 1] - F_list[i] for i in range(m_k)])
            gamma = self._solve_real_lsq(dF, f_arr)
            x_next = x_arr + self.alpha * f_arr - (dX + self.alpha * dF) @ gamma

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
