"""3D cube rank-N bulk + closed-form K_ss surface envelope.

Test if mid-f residual (0.55% at f=1e4) drops as we add more bulk modes.
Cylinder pattern: 3x per DOF. Does cube follow?

Bulk basis: phi_{m,p,q} = sin(m pi x/L) sin(p pi y/L) sin(q pi z/L), all odd
ordered by m^2 + p^2 + q^2.

  N=1:     (1,1,1)                            lam L^2/pi^2 = 3
  N=2-4:   (1,1,3),(1,3,1),(3,1,1)            lam L^2/pi^2 = 11
  N=5-7:   (1,3,3),(3,1,3),(3,3,1)            lam L^2/pi^2 = 19
  N=8-10:  (1,1,5),(1,5,1),(5,1,1)            lam L^2/pi^2 = 27
  N=11:    (3,3,3)                            lam L^2/pi^2 = 27
  ...
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


def phi_m_overlap(m, t):
    """integral_0^L f(x) sin(m pi x / L) dx for odd m.

    For odd m, direct integration of
        f(x) = cosh((x-L/2) t) / cosh(L t / 2) - 1
    against sin(m pi x / L) gives
        phi_m = 2 m pi L / (L^2 t^2 + m^2 pi^2) - 2 L / (m pi)
    """
    return 2 * m * math.pi * L / (L**2 * t**2 + m**2 * math.pi**2) - 2 * L / (m * math.pi)


def enumerate_bulk_modes(N, M_max=15):
    """Return N lowest (m,p,q) bulk modes ordered by m^2+p^2+q^2 (odd m,p,q)."""
    modes = []
    for m in range(1, M_max + 1, 2):
        for p in range(1, M_max + 1, 2):
            for q in range(1, M_max + 1, 2):
                modes.append((m, p, q, m * m + p * p + q * q))
    modes.sort(key=lambda x: x[3])
    return [(m, p, q) for (m, p, q, _) in modes[:N]]


def Y_mixed_rank_N(s, N_bulk):
    """rank-N bulk + 1-DOF surface (closed-form K_ss)."""
    t = cmath.sqrt(s * MS)
    sMS = s * MS

    J_2 = J_2_exact(t)
    J_3 = J_3_exact(t)
    F_0 = F_0_exact(t)

    # Surface
    K_ss = 3 * J_2 * J_3**2 + sMS * J_3**3
    b_psi = F_0**3

    # Bulk modes
    modes = enumerate_bulk_modes(N_bulk)

    K_mat = np.zeros((N_bulk + 1, N_bulk + 1), dtype=complex)
    b_vec = np.zeros(N_bulk + 1, dtype=complex)

    # Bulk-bulk (diagonal; sine basis is orthogonal under both Laplacian and identity)
    for i, (m, p, q) in enumerate(modes):
        lam = (m * m + p * p + q * q) * math.pi**2 / L**2
        K_mat[i, i] = (lam + sMS) * (L / 2)**3
        b_vec[i] = 8 * L**3 / (m * p * q * math.pi**3)

    # Surface
    K_mat[N_bulk, N_bulk] = K_ss
    b_vec[N_bulk] = b_psi

    # Bulk-surface coupling: integral (-Lap + sMS) phi . psi dV
    #                     = (lam + sMS) phi_m_overlap(m,t) phi_m_overlap(p,t) phi_m_overlap(q,t)
    for i, (m, p, q) in enumerate(modes):
        lam = (m * m + p * p + q * q) * math.pi**2 / L**2
        phi_overlap = phi_m_overlap(m, t) * phi_m_overlap(p, t) * phi_m_overlap(q, t)
        K_mat[i, N_bulk] = (lam + sMS) * phi_overlap
        K_mat[N_bulk, i] = K_mat[i, N_bulk]

    # Solve
    xi = np.linalg.solve(K_mat, -sMS * b_vec)
    v_avg = (xi @ b_vec) / V_CUBE
    return Y_DC * (1 + v_avg)


# Reference from cube3d_midf_verification.py (Aitken Foster N=799)
Y_AITKEN = {
    1e3: 6.892354,
    1e4: 3.431919,
    1e5: 1.218962,
    1e6: 0.399611,
    1e7: 0.126191,
}


def main():
    print("=== 3D cube rank-N bulk + closed-form K_ss surface ===")
    print(f"L = {L*1e3} mm, sigma = {SIGMA:.2e}, K_SIBC = {K_SIBC:.4e}")
    print()

    for f in [1e3, 1e4, 1e5, 1e6, 1e7]:
        s = 1j * 2 * math.pi * f
        gL = abs(cmath.sqrt(s * MS) * L)
        Y_A = Y_AITKEN[f]
        print(f"\nf = {f:.2e} Hz, |gamma L| = {gL:.2f}, |Y_Aitken_F799| = {Y_A:.6e}")
        print(f"  {'rank':>4}  {'|Y_mixed|':>14}  {'Y_mix/Y_A':>10}  {'rel err':>12}  {'top mode':>14}")
        for N in [1, 2, 4, 7, 10, 14, 20]:
            Y = Y_mixed_rank_N(s, N)
            modes = enumerate_bulk_modes(N)
            ratio = abs(Y) / Y_A
            err_pct = (ratio - 1) * 100
            print(f"  {N:4d}  {abs(Y):14.6e}  {ratio:10.5f}  {err_pct:+11.4f}%  {str(modes[-1]):>14}")


if __name__ == "__main__":
    main()
