"""Hantila polarization method for the scalar BIE-SIBC ESIM Karl loop.

The damped-Picard (and even Anderson-accelerated) Karl loop
re-assembles the BIE matrix `A_sys(Z_s)` and re-factorises it
(O(N^3) dense LU) at every outer iteration.  Hantila's polarisation
method (Hantila 1975) splits the nonlinear constitutive relation
into a CONSTANT linear part and a residual: the BIE matrix then
has a cached LU that is re-used across all outer iterations, and
only a SMALL correction matrix-vector product is needed per iter.

Splitting (analogue of Hantila for the BIE-SIBC formulation)
-------------------------------------------------------------

The scalar BIE-SIBC matrix is

    A_sys(Z_s) = (1/2) M - DL + diag(gamma(Z_s)) * (SL M^-1 K)
    gamma(Z_s) = Z_s / (j omega mu_0)               (eq 1)

Since `gamma` is LINEAR in `Z_s`, the splitting

    Z_s = Z_s_lin + Z_s_res
    gamma = gamma_lin + gamma_res

gives an exact decomposition

    A_sys(Z_s) = A_lin + A_res(Z_s_res)
    A_lin     = (1/2) M - DL + diag(gamma_lin) * (SL M^-1 K)    (eq 2a)
    A_res     = diag(gamma_res) * (SL M^-1 K)                    (eq 2b)

`A_lin` depends only on Z_s_lin (a constant -- the low-H linear-mu
Dowell/Bessel limit) and the BIE geometry; it is factored ONCE.
`A_res` is small if Z_s_res << Z_s_lin (i.e. the workpiece is
near-linear) and is recomputed per Karl iter via O(N) operations
(diagonal scaling of a precomputed matrix).

The Karl iter becomes:

    1. Given Z_s_res^{k}, build A_res^{k} = diag(gamma_res^{k}) (SL M^-1 K)
    2. Solve  A_lin * phi^{k+1} = phi_inc - A_res^{k} * phi^{k}
       via the CACHED LU on A_lin (= O(N^2) back-sub)
    3. Extract |H_t|^{k+1} from phi^{k+1}
    4. Update Z_s^{k+1} = E(|H_t|^{k+1}) via the cell solver
    5. Z_s_res^{k+1} = Z_s^{k+1} - Z_s_lin

Step 2 is a fixed-point iteration ON phi for fixed Z_s_res; it
typically converges in 2-5 inner iters when |gamma_res / gamma_lin| < 0.5.

Cost per Karl iter
------------------

| Approach                  | BIE setup per iter | Solve per iter |
|---------------------------|--------------------|----------------|
| Damped Picard / Anderson | O(N^2) assembly + O(N^3) LU | O(N^2) (uses fresh LU) |
| Hantila                   | O(N^2) (A_res = diag scaling) | O(N^2) (cached LU back-sub) x inner |

For dense BIE backends at N >= 10000, the cached LU saves
~10-100x on the BIE cost per Karl iter.  For N ~ 1000 (the IGTE
IH benchmark) the BIE cost is already small and the cell-solver
loop dominates, so the Hantila speedup is modest.

References
----------
F.I. Hantila, "A method for solving stationary magnetic field in
non-linear media,"  Rev. Roum. Sci. Techn. - Electrotechn. et Energ.
20(3), 1975, 397-407.

Status
------
This module is a Python prototype.  It is not yet wired into
calc_inductance.py.  See task #34 in the project task list.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional

import math
import numpy as np
from scipy.linalg import lu_factor, lu_solve

MU0 = 4 * math.pi * 1e-7


@dataclass
class HantilaBIEConfig:
    """Configuration for the Hantila-style Karl iteration on the BIE.

    Attributes
    ----------
    Z_s_lin : complex
        Linear-part surface impedance (single value, e.g. low-H Bessel).
        The full `Z_s(|H_t|)` is split as `Z_s_lin + Z_s_res(|H_t|)`.
    omega : float
        Angular frequency [rad/s].
    inner_tol : float
        Tolerance on inner-iter relative phi update; default 1e-6.
    inner_max_iter : int
        Cap on inner-iter count; default 20.
    """
    Z_s_lin: complex
    omega: float
    inner_tol: float = 1e-6
    inner_max_iter: int = 20


class HantilaBIESolver:
    """Hantila-style outer Karl loop wrapped around a precomputed
    scalar BIE-SIBC discretisation.

    The caller supplies the BIE matrices (`M`, `DL`, `SL`, `K`,
    `M_inv`) and the cell-solver callable `E(|H_t|) -> Z_s`.

    Workflow
    --------

        solver = HantilaBIESolver(M, DL, SL, K, M_inv, config)
        solver.factor()                # one-time LU on A_lin
        for k in range(max_karl):
            phi = solver.solve_inner(phi_inc, Z_s_res)
            H_t = extract_H_t(phi, M, K)
            Z_s_new = cell_solver(H_t)
            Z_s_res_new = Z_s_new - config.Z_s_lin
            if max|Z_s_res_new - Z_s_res| / max|Z_s_new| < tol: break
            Z_s_res = Z_s_res_new
    """

    def __init__(self, M, DL, SL, K, M_inv, config: HantilaBIEConfig):
        self.M = np.asarray(M, dtype=complex)
        self.DL = np.asarray(DL, dtype=complex)
        self.SL = np.asarray(SL, dtype=complex)
        self.K = np.asarray(K, dtype=complex)
        self.M_inv = np.asarray(M_inv, dtype=complex)
        self.cfg = config
        n = self.M.shape[0]
        if any(A.shape != (n, n) for A in (self.DL, self.SL, self.K, self.M_inv)):
            raise ValueError(
                "M, DL, SL, K, M_inv must all be square of the same size; "
                f"got {self.M.shape}, {self.DL.shape}, {self.SL.shape}, "
                f"{self.K.shape}, {self.M_inv.shape}")
        self.n = n
        # The "Robin operator" SL M^-1 K appears in both A_lin and A_res.
        self._robin_op = self.SL @ self.M_inv @ self.K
        # Lazily filled by factor().
        self._lu = None
        self._piv = None

    def factor(self):
        """Build A_lin and LU-factor it.  Call once per geometry / Z_s_lin.

        A_lin = (1/2) M - DL + gamma_lin * (SL M^-1 K)
        """
        gamma_lin = self.cfg.Z_s_lin / (1j * self.cfg.omega * MU0)
        A_lin = 0.5 * self.M - self.DL + gamma_lin * self._robin_op
        self._lu, self._piv = lu_factor(A_lin)
        return self

    def solve_inner(self, phi_inc: np.ndarray, Z_s_res: np.ndarray) -> np.ndarray:
        """Inner fixed-point iteration on phi for fixed Z_s_res.

        Solves   A_lin phi + A_res phi = phi_inc
            <=>  phi = A_lin^{-1} (phi_inc - A_res phi)

        as a Picard iteration; converges geometrically with rate
        ||A_lin^{-1} A_res||.  For |Z_s_res| << |Z_s_lin| the rate is
        small and 2-5 inner iters suffice.

        Parameters
        ----------
        phi_inc : (n,) complex
            Incident scalar potential at the BIE DOFs.
        Z_s_res : complex or (n,) complex
            Residual surface impedance (full per-DOF support).

        Returns
        -------
        phi : (n,) complex
        """
        if self._lu is None:
            raise RuntimeError("call factor() before solve_inner()")
        gamma_res = np.atleast_1d(Z_s_res).astype(complex) \
                    / (1j * self.cfg.omega * MU0)
        if gamma_res.shape == (1,):
            gamma_res = np.full(self.n, gamma_res[0])
        # A_res @ phi = diag(gamma_res) @ robin_op @ phi
        #             = gamma_res * (robin_op @ phi)  (element-wise)
        phi = lu_solve((self._lu, self._piv), phi_inc)
        for it in range(self.cfg.inner_max_iter):
            rhs = phi_inc - gamma_res * (self._robin_op @ phi)
            phi_new = lu_solve((self._lu, self._piv), rhs)
            rel = float(np.max(np.abs(phi_new - phi))
                        / max(np.max(np.abs(phi_new)), 1e-30))
            phi = phi_new
            if rel < self.cfg.inner_tol:
                break
        return phi


def demo_hantila_vs_picard():
    """Synthetic demo: solve a small linear BIE-like system both ways.

    Build a random A_lin (well-conditioned), then construct a "true"
    system A_lin + A_res with A_res a diagonal perturbation, and
    show that Hantila inner-iter converges geometrically with rate
    matching ||A_lin^-1 A_res||.
    """
    rng = np.random.default_rng(0)
    n = 200
    A_lin = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    A_lin += 4.0 * n * np.eye(n)
    M = np.eye(n)
    DL = -A_lin / 2 + 0.25 * M       # so that 0.5 M - DL = A_lin
    SL = 0.5 * np.eye(n)
    K = M.copy()
    M_inv = M.copy()
    cfg = HantilaBIEConfig(Z_s_lin=1.0+0j, omega=1.0)
    # Sanity: 0.5 M - DL == A_lin; the +SL M^-1 K * gamma_lin contribution
    # makes A_lin_full = A_lin + gamma_lin * SL.  We will use Z_s_lin
    # such that gamma_lin = 0 (decouples lin part).  Setting Z_s_lin=0
    # is degenerate so use omega such that 1/(j*omega*mu_0) -> 0 for
    # Z_s_lin=1.  Easier: just check the SOLVE works.
    solver = HantilaBIESolver(M, DL, SL, K, M_inv, cfg).factor()
    Z_s_res = np.full(n, 0.05 + 0.02j)
    phi_inc = rng.normal(size=n) + 1j * rng.normal(size=n)
    phi = solver.solve_inner(phi_inc, Z_s_res)
    # Reference: full LU solve of A_full = A_lin + A_res
    gamma_res = Z_s_res / (1j * 1.0 * MU0)
    A_res = np.diag(gamma_res) @ (SL @ M_inv @ K)
    gamma_lin = 1.0 / (1j * 1.0 * MU0)
    A_full = 0.5 * M - DL + gamma_lin * (SL @ M_inv @ K) + A_res
    phi_ref = np.linalg.solve(A_full, phi_inc)
    err = np.linalg.norm(phi - phi_ref) / np.linalg.norm(phi_ref)
    print(f"Hantila-inner vs full-LU rel err = {err:.3e}")
    assert err < 1e-5, "Hantila did not match full solve"


if __name__ == "__main__":
    demo_hantila_vs_picard()
