"""
3D cube mixed Galerkin with EXACT closed-form K_ss (no sine truncation).

The deeper question: is the SIBC-tail overshoot a sine-mode truncation
artifact, or is the tensor-product envelope psi = f(x) f(y) f(z)
genuinely incapable of representing the 6 L^2 sqrt(sigma/mu) / sqrt(s)
asymptote of the cube?

## Tensor-product envelope separability

For psi(r) = f(x) f(y) f(z):

    int |psi|^2 dV = J_3^3
    int |grad psi|^2 dV = 3 J_2 J_3^2
    int psi dV = F_0^3
    int psi phi_111 dV = (int f sin(pi x/L) dx)^3 = phi_1^3
                         where phi_1 = sum f_n / n (sin coefficient at n=1)

where:
    f(x) = cosh((x - L/2) t) / cosh(L t / 2) - 1
    J_2 = int_0^L (df/dx)^2 dx = t tanh(L t/2) - (t^2 L) / (2 cosh^2(L t/2))
    J_3 = int_0^L f^2 dx = L/(2 cosh^2(L t/2)) - 3 tanh(L t/2)/t + L
    F_0 = int_0^L f dx = (2/t) tanh(L t/2) - L
    phi_1 = int_0^L f(x) sin(pi x/L) dx
          = 2 pi L / (pi^2 + L^2 t^2) - 2 L / pi

All in closed form, no truncation.  We compare against the sine-mode
truncated version of `01_corner_envelope_uncertain.py` to isolate
which contribution is the truncation.
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
    """tanh(z) without overflow at very large Re(z)."""
    if abs(z.real) > 50:
        return complex(1.0, 0.0) if z.real > 0 else complex(-1.0, 0.0)
    return cmath.tanh(z)


def _sech2_safe(z):
    """1 / cosh^2(z), guarding against overflow.  Result ~ 4 exp(-2 |Re z|)."""
    if abs(z.real) > 50:
        return complex(0.0, 0.0)
    return 1.0 / cmath.cosh(z)**2


def J_2_exact(t):
    """∫_0^L (df/dx)^2 dx."""
    Lt2 = L * t / 2
    return t * _tanh_safe(Lt2) - (t**2 * L / 2) * _sech2_safe(Lt2)


def J_3_exact(t):
    """∫_0^L f^2 dx."""
    Lt2 = L * t / 2
    return (L / 2) * _sech2_safe(Lt2) - 3 * _tanh_safe(Lt2) / t + L


def F_0_exact(t):
    """∫_0^L f dx."""
    return (2 / t) * _tanh_safe(L * t / 2) - L


def phi1_overlap(t):
    """∫_0^L f(x) sin(pi x / L) dx (overlap with bulk phi_{1,1,1} per axis).

    f_n (sine coefficient at n=1) = (2/L) phi1_overlap.
    Direct integration gives: phi1 = 2 L pi / (pi^2 + L^2 t^2) - 2 L / pi
    """
    return 2 * L * math.pi / (math.pi**2 + L**2 * t**2) - 2 * L / math.pi


def Y_mixed_exact_Kss(s):
    """rank-1 bulk + 1-DOF surface, with EXACT closed-form K_ss."""
    t = cmath.sqrt(s * MS)
    J_2 = J_2_exact(t)
    J_3 = J_3_exact(t)
    F_0 = F_0_exact(t)
    phi1 = phi1_overlap(t)

    # Surface block: K_ss = 3 J_2 J_3^2 + sMS J_3^3
    K_ss = 3 * J_2 * J_3**2 + s * MS * J_3**3
    b_psi = F_0**3

    # Bulk: phi_{1,1,1} = sin(pi x/L) sin(pi y/L) sin(pi z/L)
    # K_bb = (lambda_{111} + sMS) (L/2)^3
    # b_b = (2 L / pi)^3
    lambda_111 = 3 * math.pi**2 / L**2
    K_b = (lambda_111 + s * MS) * (L / 2)**3
    b_b = (2 * L / math.pi)**3

    # Cross K_bs: needs ∫(-Lap + sMS) phi_111 . psi dV
    # = (lambda_{111} + sMS) ∫ phi_111 psi dV
    # ∫ phi_111 psi dV = phi1^3 (by separability)
    K_bs = (lambda_111 + s * MS) * phi1**3

    K_mat = np.array([[K_b, K_bs], [K_bs, K_ss]], dtype=complex)
    b_vec = np.array([b_b, b_psi], dtype=complex)
    xi = np.linalg.solve(K_mat, -s * MS * b_vec)
    v_avg = (xi @ b_vec) / V_CUBE
    return Y_DC * (1 + v_avg)


def Y_envelope_only_exact(s):
    """Surface envelope alone (no bulk), exact K_ss."""
    t = cmath.sqrt(s * MS)
    J_2 = J_2_exact(t)
    J_3 = J_3_exact(t)
    F_0 = F_0_exact(t)
    K_ss = 3 * J_2 * J_3**2 + s * MS * J_3**3
    b_psi = F_0**3
    xi_s = -s * MS * b_psi / K_ss
    v_avg = xi_s * b_psi / V_CUBE
    return Y_DC * (1 + v_avg)


def main():
    print("=== 3D cube: EXACT closed-form K_ss (no sine truncation) ===")
    print(f"L = {L*1e3} mm, K_SIBC = {K_SIBC:.4e}")
    print()
    print(f"{'f (Hz)':>10}  {'|gamma L|':>10}  {'|Y_env_exact|':>14}  {'|Y_mixed_exact|':>16}  {'|Y_Mellin|':>14}  {'env/Mellin':>10}  {'mix/Mellin':>10}")
    for f in [1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12]:
        s = 1j * 2 * math.pi * f
        gL = abs(cmath.sqrt(s * MS) * L)
        Y_e = Y_envelope_only_exact(s)
        Y_m = Y_mixed_exact_Kss(s)
        Y_M = Y_mellin_cube3d(s, L, SIGMA, MU)
        print(f"{f:10.2e}  {gL:10.2f}  {abs(Y_e):14.6e}  {abs(Y_m):16.6e}  {abs(Y_M):14.6e}  "
              f"{abs(Y_e)/abs(Y_M):8.4f}  {abs(Y_m)/abs(Y_M):8.4f}")

    print()
    print("If env/Mellin -> 1 at very high f, the BASIS is correct and")
    print("the earlier overshoot was sine-mode truncation. If env/Mellin")
    print("plateaus above 1, the tensor product envelope itself is the limit.")


if __name__ == "__main__":
    main()
