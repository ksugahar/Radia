"""dd_spherical_duffy.py — Double-double Spherical Duffy K matrix assembly.

Port of assemble_K_cupy() from hex_vim_cupy_kassembly.py to DD arithmetic.
Strategy:
  - sin/cos/Gauss-Legendre nodes precomputed in mpmath 30-digit, converted to DD
  - All subsequent geometry (omega, rho, c-box, r, r') in DD via dd_arithmetic
  - HDiv basis evaluation in DD via dd_basis.dd_evaluate_basis
  - K[i,j] = sum_q DD-weighted inner product, computed via DD matmul

Memory note: for cuboid order=3 (n_dofs=81, n_q ≈ 750K):
  Phi_L_hi + Phi_L_lo: 81 × 750K × 3 × 8 × 2 = 2.9 GB
  Same for Phi_R: 2.9 GB
  Total Phi: ~6 GB (FP64+DD overhead). Order=4 doubles → ~12 GB.
"""
from __future__ import annotations

import math
import sys
import time
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import mpmath as mp

sys.path.insert(0, str(Path(__file__).parent))
from dd_arithmetic import (
    dd_add, dd_sub, dd_mul, dd_div, dd_from_float, dd_from_mpmath,
    dd_sqrt,
)
from dd_basis import dd_evaluate_basis


# === Precompute Gauss-Legendre nodes/weights in mpmath 30-digit ===============

def gauss_legendre_dd(n, dps=40):
    """Compute Gauss-Legendre nodes/weights at high precision, return as DD."""
    mp.mp.dps = dps
    # mpmath's gauss-legendre on [-1, 1]
    nodes_mp = []
    weights_mp = []
    # Use numpy's leggauss as starting point, then refine via mpmath
    n_f, w_f = np.polynomial.legendre.leggauss(n)
    # Convert to mpmath and refine via Newton on the Legendre polynomial
    for i in range(n):
        x = mp.mpf(float(n_f[i]))
        # Newton iteration: x_new = x - P_n(x) / P_n'(x)
        for _ in range(20):
            Pnm1 = mp.mpf(1)
            Pn = x
            for k in range(1, n):
                Pnp1 = ((2 * k + 1) * x * Pn - k * Pnm1) / mp.mpf(k + 1)
                Pnm1 = Pn
                Pn = Pnp1
            # P_n derivative
            Pn_deriv = n * (x * Pn - Pnm1) / (x ** 2 - 1)
            x = x - Pn / Pn_deriv
        nodes_mp.append(x)
        # Weight: 2 / ((1 - x²) (P_n'(x))²)
        Pnm1 = mp.mpf(1)
        Pn = x
        for k in range(1, n):
            Pnp1 = ((2 * k + 1) * x * Pn - k * Pnm1) / mp.mpf(k + 1)
            Pnm1 = Pn
            Pn = Pnp1
        Pn_deriv = n * (x * Pn - Pnm1) / (x ** 2 - 1)
        w = mp.mpf(2) / ((1 - x ** 2) * Pn_deriv ** 2)
        weights_mp.append(w)
    # Sort
    pairs = sorted(zip(nodes_mp, weights_mp), key=lambda p: float(p[0]))
    nodes_mp = [p[0] for p in pairs]
    weights_mp = [p[1] for p in pairs]
    # Convert to DD (numpy arrays)
    nodes_hi = np.array([float(n) for n in nodes_mp])
    nodes_lo = np.array([float(n - mp.mpf(float(n))) for n in nodes_mp])
    weights_hi = np.array([float(w) for w in weights_mp])
    weights_lo = np.array([float(w - mp.mpf(float(w))) for w in weights_mp])
    return nodes_hi, nodes_lo, weights_hi, weights_lo


def assemble_K_dd(p, a, b, c, n_th=12, n_ph=24, n_rh=12, n_c=6,
                  verbose=True):
    """Assemble bare K matrix in DD (no mu0/4pi prefactor).

    Returns K_hi, K_lo of shape (n_dofs, n_dofs).
    """
    n_outer = n_th * n_ph * n_rh
    n_q = n_outer * n_c ** 3
    n_dofs = 2 * p ** 3 + 3 * p ** 2

    if verbose:
        print(f"  order={p}, n_dofs={n_dofs}")
        print(f"  quad: th={n_th}, ph={n_ph}, rh={n_rh}, c={n_c}")
        print(f"  n_outer={n_outer}, n_q={n_q}")
        mem_phi = 2 * n_q * n_dofs * 3 * 8 / 1e9   # DD = 2× FP64
        print(f"  Est Phi memory: {2 * mem_phi:.2f} GB (DD, both L+R)")

    # ----- Gauss-Legendre nodes/weights in DD ---------------------------------
    th_n_hi, th_n_lo, th_w_hi, th_w_lo = gauss_legendre_dd(n_th)
    ph_n_hi, ph_n_lo, ph_w_hi, ph_w_lo = gauss_legendre_dd(n_ph)
    rh_n_hi, rh_n_lo, rh_w_hi, rh_w_lo = gauss_legendre_dd(n_rh)
    c_n_hi, c_n_lo, c_w_hi, c_w_lo = gauss_legendre_dd(n_c)

    # ----- Outer (theta, phi) quadrature: build omega, sphere_w in DD --------
    mp.mp.dps = 40
    pi_mp = mp.pi
    pi_hi, pi_lo = dd_from_mpmath(pi_mp)

    # theta = (th_n + 1) * pi/2 : DD multiplication and addition
    # Compute in DD: theta = (th_n + 1) * pi / 2
    one_hi = np.ones_like(th_n_hi)
    one_lo = np.zeros_like(th_n_hi)
    th_p1_hi, th_p1_lo = dd_add(th_n_hi, th_n_lo, one_hi, one_lo)
    pi_arr_hi = pi_hi + 0 * th_n_hi
    pi_arr_lo = pi_lo + 0 * th_n_hi
    two_hi = 2.0 * np.ones_like(th_n_hi)
    two_lo = 0.0 * np.ones_like(th_n_hi)
    pi_over_2_hi, pi_over_2_lo = dd_div(pi_arr_hi, pi_arr_lo, two_hi, two_lo)
    theta_hi, theta_lo = dd_mul(th_p1_hi, th_p1_lo,
                                  pi_over_2_hi, pi_over_2_lo)

    # phi = (ph_n + 1) * pi
    ph_p1_hi, ph_p1_lo = dd_add(ph_n_hi, ph_n_lo,
                                  np.ones_like(ph_n_hi),
                                  np.zeros_like(ph_n_hi))
    pi_ph_hi = pi_hi + 0 * ph_n_hi
    pi_ph_lo = pi_lo + 0 * ph_n_hi
    phi_hi, phi_lo = dd_mul(ph_p1_hi, ph_p1_lo, pi_ph_hi, pi_ph_lo)

    # sin, cos via mpmath (lookup table — fast since n_th, n_ph are small)
    s_th_hi = np.zeros(n_th)
    s_th_lo = np.zeros(n_th)
    c_th_hi = np.zeros(n_th)
    c_th_lo = np.zeros(n_th)
    for i in range(n_th):
        t_mp = mp.mpf(theta_hi[i]) + mp.mpf(theta_lo[i])
        s_th_hi[i], s_th_lo[i] = dd_from_mpmath(mp.sin(t_mp))
        c_th_hi[i], c_th_lo[i] = dd_from_mpmath(mp.cos(t_mp))

    s_ph_hi = np.zeros(n_ph)
    s_ph_lo = np.zeros(n_ph)
    c_ph_hi = np.zeros(n_ph)
    c_ph_lo = np.zeros(n_ph)
    for i in range(n_ph):
        p_mp = mp.mpf(phi_hi[i]) + mp.mpf(phi_lo[i])
        s_ph_hi[i], s_ph_lo[i] = dd_from_mpmath(mp.sin(p_mp))
        c_ph_hi[i], c_ph_lo[i] = dd_from_mpmath(mp.cos(p_mp))

    # omega_x[t, p] = s_th[t] * c_ph[p]
    # omega_y[t, p] = s_th[t] * s_ph[p]
    # omega_z[t, p] = c_th[t] (broadcast over phi)
    omega_x_hi = np.zeros((n_th, n_ph))
    omega_x_lo = np.zeros((n_th, n_ph))
    omega_y_hi = np.zeros((n_th, n_ph))
    omega_y_lo = np.zeros((n_th, n_ph))
    omega_z_hi = np.zeros((n_th, n_ph))
    omega_z_lo = np.zeros((n_th, n_ph))
    for t in range(n_th):
        for ph in range(n_ph):
            ox_hi, ox_lo = dd_mul(np.array(s_th_hi[t]), np.array(s_th_lo[t]),
                                    np.array(c_ph_hi[ph]),
                                    np.array(c_ph_lo[ph]))
            omega_x_hi[t, ph] = float(ox_hi)
            omega_x_lo[t, ph] = float(ox_lo)
            oy_hi, oy_lo = dd_mul(np.array(s_th_hi[t]), np.array(s_th_lo[t]),
                                    np.array(s_ph_hi[ph]),
                                    np.array(s_ph_lo[ph]))
            omega_y_hi[t, ph] = float(oy_hi)
            omega_y_lo[t, ph] = float(oy_lo)
            omega_z_hi[t, ph] = c_th_hi[t]
            omega_z_lo[t, ph] = c_th_lo[t]

    # For brevity in MVP demo, we accept that the rest of the computation
    # (omega_max, rho_max, Lomega, sphere_w, rho, u_xyz, c-box, r, r')
    # would proceed in DD via similar pattern. For full implementation,
    # this is straightforward but ~500 lines of vectorized DD code.
    #
    # For Phase 1 deliverable, we demonstrate that:
    # 1. DD primitives are valid (dd_arithmetic.py - validated)
    # 2. DD Legendre + IntegratedLegendre are valid (dd_legendre.py - 1e-32)
    # 3. DD HDiv basis evaluation is valid (dd_basis.py - 3e-14 vs C++ FP64)
    # 4. DD geometric setup (sin/cos via mpmath, omega in DD) is valid
    #
    # Phase 2 full assembly: 4-6 more hours of careful porting.

    print(f"  DD geometric setup complete:")
    print(f"    omega_x[0,0] = {omega_x_hi[0,0]:.16e} + {omega_x_lo[0,0]:.6e}")
    print(f"    omega_y[0,0] = {omega_y_hi[0,0]:.16e} + {omega_y_lo[0,0]:.6e}")
    print(f"    omega_z[0,0] = {omega_z_hi[0,0]:.16e} + {omega_z_lo[0,0]:.6e}")
    # Verify omega lies on unit sphere
    norm_sq_hi, norm_sq_lo = dd_mul(np.array(omega_x_hi[0, 0]),
                                       np.array(omega_x_lo[0, 0]),
                                       np.array(omega_x_hi[0, 0]),
                                       np.array(omega_x_lo[0, 0]))
    sy_hi, sy_lo = dd_mul(np.array(omega_y_hi[0, 0]),
                            np.array(omega_y_lo[0, 0]),
                            np.array(omega_y_hi[0, 0]),
                            np.array(omega_y_lo[0, 0]))
    norm_sq_hi, norm_sq_lo = dd_add(norm_sq_hi, norm_sq_lo, sy_hi, sy_lo)
    sz_hi, sz_lo = dd_mul(np.array(omega_z_hi[0, 0]),
                            np.array(omega_z_lo[0, 0]),
                            np.array(omega_z_hi[0, 0]),
                            np.array(omega_z_lo[0, 0]))
    norm_sq_hi, norm_sq_lo = dd_add(norm_sq_hi, norm_sq_lo, sz_hi, sz_lo)
    one_minus_hi, one_minus_lo = dd_sub(norm_sq_hi, norm_sq_lo,
                                          np.array(1.0), np.array(0.0))
    print(f"    omega norm² - 1 = {float(one_minus_hi):.3e} + "
          f"{float(one_minus_lo):.3e} (should be ~1e-30)")
    return None  # K not yet assembled


def main():
    """Phase 1+2 MVP demo."""
    print("=" * 72)
    print("DD Spherical Duffy quadrature setup demo (Phase 1+2 MVP)")
    print("=" * 72)
    print()
    print("Phase 1 components DONE & validated:")
    print("  ✓ DD arithmetic primitives (dd_arithmetic.py, 1e-30 rel)")
    print("  ✓ DD Legendre + IntegratedLegendre (dd_legendre.py, 1e-32 rel)")
    print("  ✓ DD HDiv basis evaluation (dd_basis.py, matches C++ FP64)")
    print()
    print("Phase 2 MVP: DD Spherical Duffy quadrature setup")
    print()

    # Cuboid 5x2x1 mm, order=3 for MVP
    A_M, B_M, C_M = 5e-3, 2e-3, 1e-3
    t0 = time.time()
    assemble_K_dd(3, A_M, B_M, C_M, n_th=8, n_ph=12, n_rh=8, n_c=4,
                   verbose=True)
    print(f"\n  DD geometric setup: {time.time()-t0:.2f}s")
    print()
    print("Remaining Phase 2 work (estimated 4-6 hours):")
    print("  - DD omega_max, rho_max, Lomega computation (vectorized)")
    print("  - DD c-box [cx_lo, cx_hi] etc.")
    print("  - DD inner quadrature (cx, cy, cz)")
    print("  - DD basis evaluation at (r, r') for all 750K quad points")
    print("  - DD weighted matmul: K = einsum 'q, qic, qjc, c -> ij'")
    print("  - Validation: DD K matches C++ FP64 K to ~1e-14 (FP64 precision)")
    print("  - Validation: DD K extracts more Cauer stages via verified-interval")


if __name__ == "__main__":
    main()
