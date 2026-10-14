"""
peec_hacapk_solver.py -- Complex BiCGSTAB on top of HACApKPEECManager.

Solves the frequency-domain PEEC filament system

    (diag(R) + j*omega*L + diag(Zs)) @ I = V

where

    R       real [N]        DC resistance per filament
    Zs      complex [N]     surface impedance per filament (skin effect)
    L       real [N,N]      Ruehli mutual inductance, stored as HACApK H-matrix
    I       complex [N]     filament current (unknown)
    V       complex [N]     voltage source per filament (right-hand side)
    omega   real            = 2*pi*freq_hz

HACApK only stores the real symmetric L.  The complex system matrix is
assembled matrix-free by splitting x = x_re + j*x_im and calling the
real HACApK MatVec twice per complex matvec:

    y_re = (R + Re(Zs)) * x_re  -  (omega*L + Im(Zs)) * x_im
    y_im = (omega*L + Im(Zs)) * x_re  +  (R + Re(Zs)) * x_im

i.e. one HACApK MatVec on x_im (for the y_re row) and one on x_re
(for the y_im row).  The diagonal parts are applied elementwise.

Reference: memory/hacapk_peec_prima_plan.md "Option A"
Step 4 stage 2 of the HACApK-PEEC plan.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Optional

import numpy as np
from scipy.sparse.linalg import LinearOperator, bicgstab

from radia._radia_pybind import HACApKPEECManager


@dataclass
class PEECHACApKSolveResult:
    """Return value of :meth:`PEECHACApKSolver.solve`."""

    x: np.ndarray                   # complex [N]
    info: int                       # scipy bicgstab info (0 = converged)
    iterations: int                 # actual callback count
    final_residual: float           # ||b - A x|| / ||b||
    converged: bool
    residual_history: list[float] = field(default_factory=list)
    n_matvec_complex: int = 0       # count of complex matvecs issued
    n_matvec_real: int = 0          # = 2 * n_matvec_complex (HACApK calls)


def _default_R_dc(lengths: np.ndarray,
                  widths: np.ndarray,
                  heights: np.ndarray,
                  sigmas: np.ndarray) -> np.ndarray:
    """DC resistance R = length / (sigma * width * height) per filament."""
    sigmas = np.asarray(sigmas, dtype=np.float64)
    if np.any(sigmas <= 0.0):
        raise ValueError("all sigmas must be positive for default R_dc")
    return np.asarray(lengths, dtype=np.float64) / (
        sigmas * np.asarray(widths, dtype=np.float64)
              * np.asarray(heights, dtype=np.float64)
    )


class PEECHACApKSolver:
    """
    HACApK-accelerated complex BiCGSTAB for PEEC filament systems.

    Build once, solve many frequencies.  The H-matrix of L is independent
    of frequency, so the constructor builds it and each :meth:`solve`
    call reuses it.

    Parameters
    ----------
    centers, directions : (N, 3) float64
    lengths, widths, heights, sigmas : (N,) float64
    aca_eps, leaf_size, eta, max_rank : HACApK build parameters
        Negative values use the PEEC-tuned defaults
        (aca_eps=1e-4, leaf_size=128, eta=3.0, max_rank=400).
    print_level : int
        HACApK build verbosity (0 = silent).
    """

    def __init__(self,
                 centers: np.ndarray,
                 directions: np.ndarray,
                 lengths: np.ndarray,
                 widths: np.ndarray,
                 heights: np.ndarray,
                 sigmas: np.ndarray,
                 *,
                 aca_eps: float = -1.0,
                 leaf_size: int = -1,
                 eta: float = -1.0,
                 max_rank: int = -1,
                 print_level: int = 0):
        centers = np.ascontiguousarray(centers, dtype=np.float64)
        directions = np.ascontiguousarray(directions, dtype=np.float64)
        lengths = np.ascontiguousarray(lengths, dtype=np.float64)
        widths = np.ascontiguousarray(widths, dtype=np.float64)
        heights = np.ascontiguousarray(heights, dtype=np.float64)
        sigmas = np.ascontiguousarray(sigmas, dtype=np.float64)

        self._mgr = HACApKPEECManager(
            centers, directions, lengths, widths, heights, sigmas
        )
        ok = self._mgr.BuildHMatrix(
            aca_eps=aca_eps, leaf_size=leaf_size, eta=eta,
            max_rank=max_rank, print_level=print_level,
        )
        if not ok or not self._mgr.IsValid():
            raise RuntimeError("HACApKPEECManager.BuildHMatrix failed")

        self.N: int = int(self._mgr.GetNDOF())
        if self.N != centers.shape[0]:
            raise RuntimeError(
                f"HACApK NDOF mismatch: mgr={self.N}, geom={centers.shape[0]}"
            )
        self._R_dc_default = _default_R_dc(lengths, widths, heights, sigmas)
        self._lengths = lengths
        self._widths = widths
        self._heights = heights
        self._sigmas = sigmas
        self._L_diag: Optional[np.ndarray] = None   # lazy-extracted

    @property
    def hacapk_stats(self) -> dict:
        return self._mgr.GetStats()

    def apply_L(self, x: np.ndarray) -> np.ndarray:
        """Real matvec y = L @ x via HACApK (exposed for diagnostics)."""
        x = np.ascontiguousarray(x, dtype=np.float64)
        if x.shape != (self.N,):
            raise ValueError(f"x must have shape ({self.N},), got {x.shape}")
        return self._mgr.MatVec(x)

    def get_L_diag(self) -> np.ndarray:
        """Diagonal of L extracted via N HACApK matvecs.  Cached.

        Cost: O(N) H-matrix matvecs, each O(N log N).  Intended for
        building a Jacobi preconditioner of the PEEC system.  For very
        large N (> ~5000) consider replacing this with analytical
        self-inductance (Ruehli 1974) to avoid the O(N^2 log N) cost.
        """
        if self._L_diag is None:
            N = self.N
            diag = np.zeros(N, dtype=np.float64)
            e = np.zeros(N, dtype=np.float64)
            for j in range(N):
                e[:] = 0.0
                e[j] = 1.0
                col = self._mgr.MatVec(e)
                diag[j] = col[j]
            self._L_diag = diag
        return self._L_diag

    def build_jacobi_preconditioner(self,
                                    freq_hz: float,
                                    Zs: Optional[np.ndarray] = None,
                                    R: Optional[np.ndarray] = None
                                    ) -> LinearOperator:
        """Return M = diag(A)^{-1} as a scipy LinearOperator.

        The preconditioner is ``M @ r = r / diag_A`` where
        ``diag_A[i] = R[i] + Zs[i] + j*omega*L[i,i]``.  Using a Jacobi
        preconditioner fixes BiCGSTAB breakdown on the ill-conditioned
        PEEC system at mid-range frequencies.
        """
        N = self.N
        omega = 2.0 * np.pi * float(freq_hz)
        R_diag = self._R_dc_default if R is None else np.asarray(R, dtype=np.float64)
        if Zs is None:
            Zs_diag = np.zeros(N, dtype=np.complex128)
        else:
            Zs_diag = np.asarray(Zs, dtype=np.complex128)
        L_diag = self.get_L_diag()
        diag_A = (R_diag + Zs_diag.real) + 1j * (omega * L_diag + Zs_diag.imag)
        return LinearOperator(
            shape=(N, N),
            matvec=lambda r: np.asarray(r) / diag_A,
            dtype=np.complex128,
        )

    def build_operator(self,
                       freq_hz: float,
                       Zs: Optional[np.ndarray] = None,
                       R: Optional[np.ndarray] = None
                       ) -> tuple[LinearOperator, np.ndarray, dict]:
        """
        Return (A, diag_complex, counter) for the system
        A @ x = (diag(R) + j*omega*L + diag(Zs)) @ x.

        ``counter`` is a dict with keys ``n_complex`` and ``n_real`` that
        gets incremented in place each time ``A.matvec`` is called.
        Useful for benchmarking and for the solve() wrapper.
        """
        omega = 2.0 * np.pi * float(freq_hz)
        N = self.N

        R_diag = self._R_dc_default if R is None else np.asarray(R, dtype=np.float64)
        if R_diag.shape != (N,):
            raise ValueError(f"R must have shape ({N},), got {R_diag.shape}")

        if Zs is None:
            Zs_diag = np.zeros(N, dtype=np.complex128)
        else:
            Zs_diag = np.asarray(Zs, dtype=np.complex128)
            if Zs_diag.shape != (N,):
                raise ValueError(f"Zs must have shape ({N},), got {Zs_diag.shape}")

        diag_complex = (R_diag + Zs_diag.real) + 1j * Zs_diag.imag
        counter = {"n_complex": 0, "n_real": 0}

        def _matvec(x_in):
            x = np.asarray(x_in, dtype=np.complex128).reshape(-1)
            if x.shape != (N,):
                raise ValueError(f"expected x of shape ({N},), got {x.shape}")

            x_re = np.ascontiguousarray(x.real)
            x_im = np.ascontiguousarray(x.imag)

            # Two real HACApK matvecs per complex matvec (Option A):
            Lxr = self._mgr.MatVec(x_re)   # L @ Re(x)
            Lxi = self._mgr.MatVec(x_im)   # L @ Im(x)

            # y = diag_complex * x + j*omega * L * x
            #   diag term:  (R + Re Zs) x + j (Im Zs) x
            #   j*omega*L term:  -omega*Lxi + j*omega*Lxr
            y = diag_complex * x
            y.real -= omega * Lxi
            y.imag += omega * Lxr

            counter["n_complex"] += 1
            counter["n_real"] += 2
            return y

        A = LinearOperator(shape=(N, N), matvec=_matvec, dtype=np.complex128)
        return A, diag_complex, counter

    def solve(self,
              freq_hz: float,
              rhs: np.ndarray,
              *,
              Zs: Optional[np.ndarray] = None,
              R: Optional[np.ndarray] = None,
              rtol: float = 1e-6,
              atol: float = 0.0,
              maxiter: Optional[int] = None,
              x0: Optional[np.ndarray] = None,
              M: Optional[LinearOperator] = None,
              callback: Optional[Callable[[np.ndarray], None]] = None,
              track_residual: bool = True,
              ) -> PEECHACApKSolveResult:
        """
        Solve (R + j*omega*L + Zs) I = rhs via scipy BiCGSTAB.

        ``M`` is a preconditioner LinearOperator (optional).  For PEEC,
        a cheap choice is diagonal of A: ``M.matvec(r) = r / diag_A``.
        If ``track_residual`` is True the relative residual is recomputed
        at every callback tick (one extra complex matvec per iteration).
        """
        N = self.N
        rhs = np.asarray(rhs, dtype=np.complex128).reshape(-1)
        if rhs.shape != (N,):
            raise ValueError(f"rhs must have shape ({N},), got {rhs.shape}")
        b_norm = np.linalg.norm(rhs)
        if b_norm == 0.0:
            return PEECHACApKSolveResult(
                x=np.zeros(N, dtype=np.complex128),
                info=0, iterations=0, final_residual=0.0, converged=True,
            )

        A, _diag, counter = self.build_operator(freq_hz, Zs=Zs, R=R)

        residuals: list[float] = []
        n_cb = [0]

        def _cb(xk):
            n_cb[0] += 1
            if track_residual:
                r = rhs - A.matvec(xk)
                residuals.append(float(np.linalg.norm(r) / b_norm))
            if callback is not None:
                callback(xk)

        x0_arg = None if x0 is None else np.asarray(x0, dtype=np.complex128)

        x, info = bicgstab(
            A, rhs, x0=x0_arg,
            rtol=rtol, atol=atol, maxiter=maxiter, M=M, callback=_cb,
        )

        final_res = float(np.linalg.norm(rhs - A.matvec(x)) / b_norm)
        # scipy BiCGSTAB returns info<0 on breakdown; the final iterate
        # is still usable if the residual already meets rtol.  Treat
        # residual-based convergence as authoritative.
        converged = final_res <= max(rtol, atol / b_norm) + 1e-14

        return PEECHACApKSolveResult(
            x=x,
            info=int(info),
            iterations=int(n_cb[0]),
            final_residual=final_res,
            converged=bool(converged),
            residual_history=residuals,
            n_matvec_complex=int(counter["n_complex"]),
            n_matvec_real=int(counter["n_real"]),
        )


def solve_peec_hacapk(centers: np.ndarray,
                      directions: np.ndarray,
                      lengths: np.ndarray,
                      widths: np.ndarray,
                      heights: np.ndarray,
                      sigmas: np.ndarray,
                      freq_hz: float,
                      rhs: np.ndarray,
                      *,
                      Zs: Optional[np.ndarray] = None,
                      R: Optional[np.ndarray] = None,
                      rtol: float = 1e-6,
                      maxiter: Optional[int] = None,
                      x0: Optional[np.ndarray] = None,
                      hacapk_params: Optional[dict] = None,
                      print_level: int = 0,
                      ) -> PEECHACApKSolveResult:
    """
    Convenience one-shot: build HACApK H-matrix and solve at a single
    frequency.  Use :class:`PEECHACApKSolver` directly for frequency
    sweeps (the H-matrix build is reused).
    """
    hp = hacapk_params or {}
    solver = PEECHACApKSolver(
        centers, directions, lengths, widths, heights, sigmas,
        aca_eps=hp.get("aca_eps", -1.0),
        leaf_size=hp.get("leaf_size", -1),
        eta=hp.get("eta", -1.0),
        max_rank=hp.get("max_rank", -1),
        print_level=print_level,
    )
    return solver.solve(
        freq_hz, rhs, Zs=Zs, R=R, rtol=rtol, maxiter=maxiter, x0=x0,
    )


# ---------------------------------------------------------------------------
# Self-test: helical coil, compare complex BiCGSTAB vs dense LU
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import math
    import time

    def _helical_geometry(turns: int, seg_per_turn: int,
                          radius: float, pitch: float,
                          w: float, h: float, sigma: float):
        n = turns * seg_per_turn
        centers = np.zeros((n, 3), dtype=np.float64)
        directions = np.zeros((n, 3), dtype=np.float64)
        lengths = np.zeros(n, dtype=np.float64)
        for i in range(n):
            t0 = 2.0 * math.pi * i / seg_per_turn
            t1 = 2.0 * math.pi * (i + 1) / seg_per_turn
            z0 = pitch * (i / seg_per_turn)
            z1 = pitch * ((i + 1) / seg_per_turn)
            p0 = np.array([radius * math.cos(t0), radius * math.sin(t0), z0])
            p1 = np.array([radius * math.cos(t1), radius * math.sin(t1), z1])
            d = p1 - p0
            L = float(np.linalg.norm(d))
            centers[i] = 0.5 * (p0 + p1)
            directions[i] = d / L
            lengths[i] = L
        widths = np.full(n, w)
        heights = np.full(n, h)
        sigmas = np.full(n, sigma)
        return centers, directions, lengths, widths, heights, sigmas

    N_TURNS = 8
    SEG_PER_TURN = 16
    N = N_TURNS * SEG_PER_TURN
    centers, directions, lengths, widths, heights, sigmas = _helical_geometry(
        turns=N_TURNS, seg_per_turn=SEG_PER_TURN,
        radius=0.050, pitch=0.006,
        w=0.003, h=0.003, sigma=5.8e7,
    )
    FREQ = 10.0e3
    OMEGA = 2.0 * math.pi * FREQ

    print(f"[self-test] N = {N} filaments, freq = {FREQ*1e-3:.1f} kHz")

    t0 = time.perf_counter()
    solver = PEECHACApKSolver(centers, directions, lengths,
                              widths, heights, sigmas,
                              print_level=0)
    t_build = time.perf_counter() - t0
    print(f"  HACApK build: {t_build*1e3:.1f} ms  stats={solver.hacapk_stats}")

    # --- Extract dense L by N real HACApK matvecs (reference) ---
    t0 = time.perf_counter()
    L_dense = np.zeros((N, N), dtype=np.float64)
    e = np.zeros(N, dtype=np.float64)
    for j in range(N):
        e[:] = 0.0; e[j] = 1.0
        L_dense[:, j] = solver.apply_L(e)
    t_extract = time.perf_counter() - t0
    print(f"  dense-L extraction: {t_extract*1e3:.1f} ms "
          f"(||L - L^T||_F / ||L||_F = "
          f"{np.linalg.norm(L_dense - L_dense.T) / np.linalg.norm(L_dense):.2e})")

    R_dc = _default_R_dc(lengths, widths, heights, sigmas)
    A_dense = (np.diag(R_dc).astype(np.complex128)
               + 1j * OMEGA * L_dense)

    rng = np.random.default_rng(42)
    rhs = rng.standard_normal(N) + 1j * rng.standard_normal(N)

    # --- Dense LU reference ---
    t0 = time.perf_counter()
    x_dense = np.linalg.solve(A_dense, rhs)
    t_lu = time.perf_counter() - t0

    # --- HACApK + complex BiCGSTAB ---
    # Cheap Jacobi preconditioner: M = 1 / diag(A)
    diag_A = R_dc.astype(np.complex128) + 1j * OMEGA * np.diag(L_dense)
    M = LinearOperator((N, N),
                       matvec=lambda r: np.asarray(r) / diag_A,
                       dtype=np.complex128)

    t0 = time.perf_counter()
    res = solver.solve(FREQ, rhs, rtol=1e-10, maxiter=2000, M=M,
                       track_residual=False)
    t_bicg = time.perf_counter() - t0

    rel_err = np.linalg.norm(res.x - x_dense) / np.linalg.norm(x_dense)
    print(f"  dense LU:       {t_lu*1e3:.1f} ms")
    print(f"  HACApK BiCGSTAB:{t_bicg*1e3:.1f} ms "
          f"iters={res.iterations} info={res.info} "
          f"final_res={res.final_residual:.2e} "
          f"n_real_matvec={res.n_matvec_real}")
    print(f"  ||x_hacapk - x_LU|| / ||x_LU|| = {rel_err:.2e}")

    assert res.converged, f"BiCGSTAB did not converge: info={res.info}"
    assert rel_err < 1e-6, (
        f"HACApK solution disagrees with dense LU: rel_err={rel_err:.2e}"
    )
    print("  [self-test PASS]")
