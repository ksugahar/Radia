"""
3D cube mixed Galerkin verification at MID frequency.

At high f (Mellin valid) my closed-form K_ss matches Mellin to 4 decimals.
But at MID f (wall band region, where neither Mellin nor truncated Foster
is a clean ground truth), the framework is NOT yet rigorously verified.

This script uses Foster sum at high N (up to 999 with chunked
vectorisation) + Richardson extrapolation to get a converged
reference at mid f.  Then compares with mixed Galerkin (closed form K_ss).
"""
import cmath
import math
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _references.cube3d_foster import K_SIBC_cube3d, Y_DC_cube3d, Y_mellin_cube3d

SIGMA = 5.8e7
MU = 4 * math.pi * 1e-7
L = 5e-3
MS = MU * SIGMA
Y_DC = Y_DC_cube3d(L, SIGMA)
K_SIBC = K_SIBC_cube3d(L, SIGMA, MU)
V_CUBE = L**3


def Y_foster_3d(s, N_max):
    """3D Foster sum, n-loop with (p,q) vectorized."""
    odd = np.arange(1, N_max + 1, 2, dtype=np.float64)
    pi_L2 = math.pi**2 / L**2
    sMS = s * mu * sigma if False else s * MU * SIGMA  # avoid clash
    sMS = s * MS

    total = 0.0 + 0j
    pp, qq = np.meshgrid(odd, odd, indexing="ij")
    pq_sq = pp**2 * qq**2
    p2_q2 = pp**2 + qq**2
    for n in odd:
        n2 = n * n
        lam = (n2 + p2_q2) * pi_L2
        term = lam / ((lam + sMS) * (n2 * pq_sq))
        total += np.sum(term)
    return Y_DC * (512.0 / math.pi**6) * total


def _tanh_safe(z):
    if abs(z.real) > 50:
        return complex(1.0, 0.0) if z.real > 0 else complex(-1.0, 0.0)
    return cmath.tanh(z)


def _sech2_safe(z):
    if abs(z.real) > 50:
        return complex(0.0, 0.0)
    return 1.0 / cmath.cosh(z)**2


def J_2_exact(t):
    Lt2 = L * t / 2
    return t * _tanh_safe(Lt2) - (t**2 * L / 2) * _sech2_safe(Lt2)


def J_3_exact(t):
    Lt2 = L * t / 2
    return (L / 2) * _sech2_safe(Lt2) - 3 * _tanh_safe(Lt2) / t + L


def F_0_exact(t):
    return (2 / t) * _tanh_safe(L * t / 2) - L


def phi1_overlap(t):
    return 2 * L * math.pi / (math.pi**2 + L**2 * t**2) - 2 * L / math.pi


def Y_mixed_exact_Kss(s):
    t = cmath.sqrt(s * MS)
    J_2 = J_2_exact(t)
    J_3 = J_3_exact(t)
    F_0 = F_0_exact(t)
    phi1 = phi1_overlap(t)

    K_ss = 3 * J_2 * J_3**2 + s * MS * J_3**3
    b_psi = F_0**3
    lambda_111 = 3 * math.pi**2 / L**2
    K_b = (lambda_111 + s * MS) * (L / 2)**3
    b_b = (2 * L / math.pi)**3
    K_bs = (lambda_111 + s * MS) * phi1**3

    K_mat = np.array([[K_b, K_bs], [K_bs, K_ss]], dtype=complex)
    b_vec = np.array([b_b, b_psi], dtype=complex)
    xi = np.linalg.solve(K_mat, -s * MS * b_vec)
    v_avg = (xi @ b_vec) / V_CUBE
    return Y_DC * (1 + v_avg)


def main():
    print("=== 3D cube mid-f verification: Foster N sweep + Richardson ===")
    print(f"L = {L*1e3} mm")
    print()

    # At mid f, compute Foster at high N.
    for f in [1e3, 1e4, 1e5, 1e6, 1e7]:
        s = 1j * 2 * math.pi * f
        gL = abs(cmath.sqrt(s * MS) * L)
        Y_mix = Y_mixed_exact_Kss(s)
        Y_Mellin = Y_mellin_cube3d(s, L, SIGMA, MU)

        print(f"\nf = {f:.2e} Hz, |gamma L| = {gL:.2f}")
        print(f"  |Y_mixed| (closed K_ss) = {abs(Y_mix):.6e}")
        print(f"  |Y_Mellin (3-term)|     = {abs(Y_Mellin):.6e}")

        # Foster N sweep
        Y_F = {}
        for N in [99, 199, 399, 799]:
            t0 = time.time()
            Y_F[N] = Y_foster_3d(s, N)
            dt = time.time() - t0
            print(f"  Foster N={N:4d}: |Y_F|={abs(Y_F[N]):.6e}  (Y_F/Y_mix={abs(Y_F[N])/abs(Y_mix):.5f})  [{dt:.1f}s]")

        # Aitken extrapolation on three highest N
        Y1, Y2, Y3 = Y_F[199], Y_F[399], Y_F[799]
        dY1 = Y2 - Y1
        dY2 = Y3 - Y2
        d2Y = dY2 - dY1
        if abs(d2Y) > abs(Y3) * 1e-13:
            Y_A = Y3 - dY2**2 / d2Y
        else:
            Y_A = Y3
        print(f"  Aitken (199,399,799):   |Y_A|={abs(Y_A):.6e}  (Y_mix/Y_A={abs(Y_mix)/abs(Y_A):.5f})")


if __name__ == "__main__":
    main()
