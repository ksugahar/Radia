"""
Diagnostic for the 3D cube mixed Galerkin SIBC-tail overshoot.

Question: at f >> 1 MHz, |Y_mixed| > |K_SIBC/sqrt(s)| by a factor that
grows with f.  Is this due to:
    (a) Insufficient sine-mode truncation N_PSI in the surface basis?
    (b) Bulk-surface coupling at deep skin?
    (c) Genuine envelope insufficiency (tensor product can't represent
        all 6 faces + 12 edges + 8 corners with correct weights)?

Strategy:
    1. Compute Y_envelope_only (no bulk) at sweep N_PSI = 49, 99, 199.
       If overshoot depends on N_PSI -> (a).
    2. Compare to Y_mixed (rank-1 bulk + 1-DOF surface).
       If different from envelope-only -> (b).
    3. The deep-skin asymptote of envelope-only is computable
       analytically: derive the ratio Y_env/K_SIBC at sqrt(s) -> infty.
"""
import cmath
import math
import sys
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


def f_n(n, t):
    return -4 * L**2 * t**2 / (n * math.pi * (n**2 * math.pi**2 + L**2 * t**2))


def F_0(t):
    return (2.0 / t) * cmath.tanh(L * t / 2) - L


def Y_envelope_only(s, N_psi):
    """Surface envelope alone (no bulk), to isolate the envelope's
    intrinsic K_SIBC^eff."""
    t = cmath.sqrt(s * MS)
    odd_n = np.arange(1, 2 * N_psi, 2)
    f_arr = np.array([f_n(int(n), t) for n in odd_n], dtype=complex)
    nn, pp, qq = np.meshgrid(odd_n, odd_n, odd_n, indexing="ij")
    lam_3d = (nn**2 + pp**2 + qq**2) * math.pi**2 / L**2
    f_3d = np.einsum("i,j,k->ijk", f_arr, f_arr, f_arr)
    K_ss = (L / 2)**3 * np.sum((lam_3d + s * MS) * f_3d**2)
    b_psi = F_0(t)**3
    xi_s = -s * MS * b_psi / K_ss
    v_avg = xi_s * b_psi / V_CUBE
    return Y_DC * (1 + v_avg)


def Y_mixed(s, N_psi):
    """rank-1 bulk + 1-DOF surface."""
    t = cmath.sqrt(s * MS)
    odd_n = np.arange(1, 2 * N_psi, 2)
    f_arr = np.array([f_n(int(n), t) for n in odd_n], dtype=complex)
    nn, pp, qq = np.meshgrid(odd_n, odd_n, odd_n, indexing="ij")
    lam_3d = (nn**2 + pp**2 + qq**2) * math.pi**2 / L**2
    f_3d = np.einsum("i,j,k->ijk", f_arr, f_arr, f_arr)
    K_ss = (L / 2)**3 * np.sum((lam_3d + s * MS) * f_3d**2)
    b_psi = F_0(t)**3
    lambda_111 = 3 * math.pi**2 / L**2
    K_b = (lambda_111 + s * MS) * (L / 2)**3
    K_bs = (lambda_111 + s * MS) * (L / 2)**3 * f_arr[0]**3
    b_b = 8 * L**3 / math.pi**3
    K_mat = np.array([[K_b, K_bs], [K_bs, K_ss]], dtype=complex)
    b_vec = np.array([b_b, b_psi], dtype=complex)
    xi = np.linalg.solve(K_mat, -s * MS * b_vec)
    v_avg = (xi @ b_vec) / V_CUBE
    return Y_DC * (1 + v_avg)


def main():
    print("=== 3D cube SIBC tail diagnostic ===")
    print(f"True K_SIBC = {K_SIBC:.4e}")
    print()

    # (1) N_PSI sensitivity at high f
    print("(1) Envelope-only Y at f=1e8 Hz, sweep N_PSI:")
    print(f"{'N_PSI':>6}  {'|Y_env|':>14}  {'|Y_env|/|K_SIBC/sqrt(s)|':>24}")
    s = 1j * 2 * math.pi * 1e8
    K_sqrts = abs(K_SIBC / cmath.sqrt(s))
    for Np in [25, 49, 99, 199]:
        Y_e = Y_envelope_only(s, Np)
        print(f"{Np:6d}  {abs(Y_e):14.6e}  {abs(Y_e)/K_sqrts:22.4f}")

    # (2) Compare envelope-only vs mixed at sweep f
    print()
    print("(2) Compare envelope-only vs mixed (rank-1 bulk + surface), N_PSI=99:")
    print(f"{'f (Hz)':>10}  {'|Y_env|':>12}  {'|Y_mixed|':>12}  {'|Y_Mellin|':>12}  {'env/Mellin':>10}  {'mix/Mellin':>10}")
    for f in [1e5, 1e6, 1e7, 1e8, 1e9, 1e10]:
        s = 1j * 2 * math.pi * f
        Y_env = Y_envelope_only(s, 99)
        Y_mix = Y_mixed(s, 99)
        Y_M = Y_mellin_cube3d(s, L, SIGMA, MU)
        print(f"{f:10.2e}  {abs(Y_env):12.4e}  {abs(Y_mix):12.4e}  {abs(Y_M):12.4e}  "
              f"{abs(Y_env)/abs(Y_M):8.4f}  {abs(Y_mix)/abs(Y_M):8.4f}")

    # (3) Asymptotic ratio: |Y_env|/(K_SIBC/sqrt(s)) at increasing f
    print()
    print("(3) Asymptotic envelope/K_SIBC ratio at very high f (N_PSI=199):")
    print(f"{'f (Hz)':>10}  {'|gamma L|':>10}  {'|Y_env|/|K_SIBC/sqrt(s)|':>24}")
    for f in [1e6, 1e7, 1e8, 1e9, 1e10, 1e11]:
        s = 1j * 2 * math.pi * f
        gL = abs(cmath.sqrt(s * MS) * L)
        K_sqrts = abs(K_SIBC / cmath.sqrt(s))
        Y_e = Y_envelope_only(s, 199)
        print(f"{f:10.2e}  {gL:10.2f}  {abs(Y_e)/K_sqrts:22.4f}")


if __name__ == "__main__":
    main()
