"""dd_k_entry_proof.py — Proof-of-concept DD computation of single K[i,j] entry.

Computes K[0,0] = ∫∫_{[0,1]^6} φ_0(r) · diag(a²,b²,c²) · φ_0(r') / |Λ(r-r')|
dV dV' in DD arithmetic, validates against C++ FP64
compute_K_entry_spherical_duffy.

Slow (nested loops, no GPU), but demonstrates correctness of the DD pipeline.
For full vectorized DD K assembly, this proof-of-concept extends to all
(i, j) pairs via einsum-like construction with DD arithmetic.
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


def dd_compute_K_entry(p, i, j, a, b, c, n_th=8, n_ph=12, n_rh=8, n_c=4):
    """Compute single K[i,j] entry of bare Spherical Duffy integral in DD.

    Follows the FP64 algorithm in C++ compute_K_entry_spherical_duffy,
    with all arithmetic in DD via dd_arithmetic.

    Returns: (K_hi, K_lo) DD result.
    """
    mp.mp.dps = 40

    # Gauss-Legendre nodes/weights in mpmath (then convert to DD)
    th_n_f, th_w_f = np.polynomial.legendre.leggauss(n_th)
    ph_n_f, ph_w_f = np.polynomial.legendre.leggauss(n_ph)
    rh_n_f, rh_w_f = np.polynomial.legendre.leggauss(n_rh)
    c_n_f, c_w_f = np.polynomial.legendre.leggauss(n_c)

    # Refine nodes via mpmath Newton (skip for MVP - use FP64 nodes)
    # For full precision: use the gauss_legendre_dd from dd_spherical_duffy.py

    pi_mp = mp.pi

    # Precompute sin/cos at theta = (th_n + 1) * pi/2, phi = (ph_n + 1) * pi
    th_sin = []
    th_cos = []
    theta_vals = []
    for k in range(n_th):
        theta_mp = (mp.mpf(float(th_n_f[k])) + 1) * pi_mp / 2
        theta_vals.append(theta_mp)
        th_sin.append(mp.sin(theta_mp))
        th_cos.append(mp.cos(theta_mp))

    ph_sin = []
    ph_cos = []
    phi_vals = []
    for k in range(n_ph):
        phi_mp = (mp.mpf(float(ph_n_f[k])) + 1) * pi_mp
        phi_vals.append(phi_mp)
        ph_sin.append(mp.sin(phi_mp))
        ph_cos.append(mp.cos(phi_mp))

    jac_th_mp = pi_mp / 2
    jac_ph_mp = pi_mp

    a_sq = mp.mpf(a) * mp.mpf(a)
    b_sq = mp.mpf(b) * mp.mpf(b)
    c_sq = mp.mpf(c) * mp.mpf(c)

    # Accumulate result in DD
    result_hi = np.float64(0.0)
    result_lo = np.float64(0.0)

    for it in range(n_th):
        s_th_mp = th_sin[it]
        c_th_mp = th_cos[it]

        for ip in range(n_ph):
            c_ph_mp = ph_cos[ip]
            s_ph_mp = ph_sin[ip]

            omega_x_mp = s_th_mp * c_ph_mp
            omega_y_mp = s_th_mp * s_ph_mp
            omega_z_mp = c_th_mp

            omega_max_mp = max(abs(omega_x_mp), abs(omega_y_mp),
                                abs(omega_z_mp))
            if omega_max_mp < mp.mpf("1e-30"):
                continue
            rho_max_mp = mp.mpf(1) / omega_max_mp

            Lomega_norm_mp = mp.sqrt(a_sq * omega_x_mp ** 2
                                       + b_sq * omega_y_mp ** 2
                                       + c_sq * omega_z_mp ** 2)

            sphere_w_mp = (mp.mpf(float(th_w_f[it])) * mp.mpf(float(ph_w_f[ip]))
                            * jac_th_mp * jac_ph_mp
                            * s_th_mp / Lomega_norm_mp)

            for ir in range(n_rh):
                rho_val_mp = (mp.mpf(float(rh_n_f[ir])) + 1) * rho_max_mp / 2
                jac_rh_mp = rho_max_mp / 2
                u_x_mp = rho_val_mp * omega_x_mp
                u_y_mp = rho_val_mp * omega_y_mp
                u_z_mp = rho_val_mp * omega_z_mp

                cx_lo_mp = abs(u_x_mp) / 2
                cy_lo_mp = abs(u_y_mp) / 2
                cz_lo_mp = abs(u_z_mp) / 2
                cx_jac_mp = (mp.mpf(1) - abs(u_x_mp)) / 2
                cy_jac_mp = (mp.mpf(1) - abs(u_y_mp)) / 2
                cz_jac_mp = (mp.mpf(1) - abs(u_z_mp)) / 2
                if (cx_jac_mp <= 0 or cy_jac_mp <= 0 or cz_jac_mp <= 0):
                    continue

                c_inner_mp = mp.mpf(0)
                for icx in range(n_c):
                    cx_mp = cx_lo_mp + (mp.mpf(float(c_n_f[icx])) + 1) * cx_jac_mp
                    wcx_mp = mp.mpf(float(c_w_f[icx]))
                    for icy in range(n_c):
                        cy_mp = cy_lo_mp + (mp.mpf(float(c_n_f[icy])) + 1) * cy_jac_mp
                        wcy_mp = mp.mpf(float(c_w_f[icy]))
                        for icz in range(n_c):
                            cz_mp = cz_lo_mp + (mp.mpf(float(c_n_f[icz])) + 1) * cz_jac_mp
                            wcz_mp = mp.mpf(float(c_w_f[icz]))

                            rx_mp = cx_mp + u_x_mp / 2
                            ry_mp = cy_mp + u_y_mp / 2
                            rz_mp = cz_mp + u_z_mp / 2
                            rpx_mp = cx_mp - u_x_mp / 2
                            rpy_mp = cy_mp - u_y_mp / 2
                            rpz_mp = cz_mp - u_z_mp / 2

                            # Basis evaluation in DD (convert mp → DD inputs)
                            rx_hi, rx_lo = dd_from_mpmath(rx_mp)
                            ry_hi, ry_lo = dd_from_mpmath(ry_mp)
                            rz_hi, rz_lo = dd_from_mpmath(rz_mp)
                            rpx_hi, rpx_lo = dd_from_mpmath(rpx_mp)
                            rpy_hi, rpy_lo = dd_from_mpmath(rpy_mp)
                            rpz_hi, rpz_lo = dd_from_mpmath(rpz_mp)

                            # Need basis at SINGLE point, so wrap in array
                            x_arr_hi = np.array([rx_hi])
                            x_arr_lo = np.array([rx_lo])
                            y_arr_hi = np.array([ry_hi])
                            y_arr_lo = np.array([ry_lo])
                            z_arr_hi = np.array([rz_hi])
                            z_arr_lo = np.array([rz_lo])
                            Phi_i_hi, Phi_i_lo = dd_evaluate_basis(
                                x_arr_hi, x_arr_lo, y_arr_hi, y_arr_lo,
                                z_arr_hi, z_arr_lo, p, xp_mod=np)

                            xp_arr_hi = np.array([rpx_hi])
                            xp_arr_lo = np.array([rpx_lo])
                            yp_arr_hi = np.array([rpy_hi])
                            yp_arr_lo = np.array([rpy_lo])
                            zp_arr_hi = np.array([rpz_hi])
                            zp_arr_lo = np.array([rpz_lo])
                            Phi_j_hi, Phi_j_lo = dd_evaluate_basis(
                                xp_arr_hi, xp_arr_lo, yp_arr_hi, yp_arr_lo,
                                zp_arr_hi, zp_arr_lo, p, xp_mod=np)

                            # Dot product φ_i · diag(a²,b²,c²) · φ_j in DD
                            # x-component
                            dx_hi, dx_lo = dd_mul(
                                np.array(Phi_i_hi[0, i, 0]),
                                np.array(Phi_i_lo[0, i, 0]),
                                np.array(Phi_j_hi[0, j, 0]),
                                np.array(Phi_j_lo[0, j, 0]))
                            a_sq_hi, a_sq_lo = dd_from_mpmath(a_sq)
                            dx_hi, dx_lo = dd_mul(
                                dx_hi, dx_lo,
                                np.array(a_sq_hi), np.array(a_sq_lo))

                            # y
                            dy_hi, dy_lo = dd_mul(
                                np.array(Phi_i_hi[0, i, 1]),
                                np.array(Phi_i_lo[0, i, 1]),
                                np.array(Phi_j_hi[0, j, 1]),
                                np.array(Phi_j_lo[0, j, 1]))
                            b_sq_hi, b_sq_lo = dd_from_mpmath(b_sq)
                            dy_hi, dy_lo = dd_mul(
                                dy_hi, dy_lo,
                                np.array(b_sq_hi), np.array(b_sq_lo))

                            # z
                            dz_hi, dz_lo = dd_mul(
                                np.array(Phi_i_hi[0, i, 2]),
                                np.array(Phi_i_lo[0, i, 2]),
                                np.array(Phi_j_hi[0, j, 2]),
                                np.array(Phi_j_lo[0, j, 2]))
                            c_sq_hi, c_sq_lo = dd_from_mpmath(c_sq)
                            dz_hi, dz_lo = dd_mul(
                                dz_hi, dz_lo,
                                np.array(c_sq_hi), np.array(c_sq_lo))

                            # dot = dx + dy + dz
                            dot_hi, dot_lo = dd_add(dx_hi, dx_lo, dy_hi, dy_lo)
                            dot_hi, dot_lo = dd_add(dot_hi, dot_lo,
                                                      dz_hi, dz_lo)

                            # Weight in mpmath, accumulate to DD
                            inner_weight_mp = wcx_mp * wcy_mp * wcz_mp
                            full_weight_mp = (sphere_w_mp
                                                * mp.mpf(float(rh_w_f[ir]))
                                                * jac_rh_mp * rho_val_mp
                                                * cx_jac_mp * cy_jac_mp
                                                * cz_jac_mp * inner_weight_mp)

                            # contribution = full_weight * dot
                            w_hi, w_lo = dd_from_mpmath(full_weight_mp)
                            contrib_hi, contrib_lo = dd_mul(
                                np.array(w_hi), np.array(w_lo),
                                dot_hi, dot_lo)

                            # result += contrib
                            result_hi, result_lo = dd_add(
                                np.array(result_hi), np.array(result_lo),
                                contrib_hi, contrib_lo)
                            result_hi = float(result_hi)
                            result_lo = float(result_lo)

    return result_hi, result_lo


def main():
    print("=" * 72)
    print("DD compute_K_entry proof-of-concept: cuboid 5×2×1, order=2")
    print("=" * 72)

    sys.path.insert(0,
        str(Path("S:/Radia/01_GitHub/src/ext/radia_vim/src")))
    import radia_vim

    p = 2  # order=2 (n_dofs=28) for MVP demo
    A_M, B_M, C_M = 5e-3, 2e-3, 1e-3

    # Reduced quadrature for fast MVP demo
    n_th, n_ph, n_rh, n_c = 4, 6, 4, 2

    print(f"\n  Configuration: order={p}, quad=({n_th},{n_ph},{n_rh},{n_c})")
    print(f"  n_q = {n_th*n_ph*n_rh*n_c**3}")

    # Compute K[0,0] in DD
    print(f"\n  Computing K[0,0] in DD...")
    t0 = time.time()
    K_hi, K_lo = dd_compute_K_entry(p, 0, 0, A_M, B_M, C_M,
                                      n_th, n_ph, n_rh, n_c)
    t_dd = time.time() - t0
    print(f"  DD result: K[0,0]_hi = {float(K_hi):.16e}")
    print(f"             K[0,0]_lo = {float(K_lo):.3e}")
    print(f"  Time: {t_dd:.2f}s")

    # C++ FP64 reference
    basis = radia_vim.HDivDivFreeHexBasis(p)
    K_cpp = radia_vim.compute_K_entry_spherical_duffy(
        basis, 0, 0, A_M, B_M, C_M, n_th, n_ph, n_rh, n_c)
    print(f"\n  C++ FP64 ref: K[0,0] = {K_cpp:.16e}")

    diff = abs(float(K_hi) - K_cpp)
    rel_err = diff / abs(K_cpp) if abs(K_cpp) > 0 else diff
    print(f"\n  Difference DD_hi vs C++ FP64: abs={diff:.3e}, rel={rel_err:.3e}")
    print(f"  (Expected: rel ~ 1e-14 to 1e-15, FP64 precision limit)")
    if rel_err < 1e-13:
        print(f"  PASS: DD computation matches C++ FP64 within precision limit")
    else:
        print(f"  CHECK: difference larger than expected — investigate")


if __name__ == "__main__":
    main()
