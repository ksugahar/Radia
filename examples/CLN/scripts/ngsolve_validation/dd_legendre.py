"""dd_legendre.py — Double-double Legendre polynomials for HDiv basis evaluation.

Implements IntegratedLegendre L_{j+2}(xi) and classical Legendre P_n(xi) in
DD arithmetic, matching the FP64 evaluation in hex_vim_cupy.py.

Validation against mpmath 30+ digit reference at sample points.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
sys.path.insert(0, str(Path(__file__).parent))
from dd_arithmetic import (
    dd_add, dd_sub, dd_mul, dd_div, dd_from_float, dd_from_mpmath,
)


def dd_integrated_legendre_with_deriv(xi_hi, xi_lo, p):
    """DD version of integrated_legendre_with_deriv.

    Mirrors the FP64 recurrence in hex_vim_cupy.py:
      p3=0, p2=-1, p1=xi; p3'=0, p2'=0, p1'=1
      For j=2..p+1:
        new_p1  = ((2j-3)/j) xi p2 - ((j-3)/j) p3
        new_p1p = ((2j-3)/j) (p2 + xi p2') - ((j-3)/j) p3'
        L[j-2]  = new_p1
        Lp[j-2] = new_p1p
        shift down
    """
    # State as DD tuples (hi, lo)
    zero_hi = 0.0 * xi_hi
    zero_lo = 0.0 * xi_hi
    one_hi = 0.0 * xi_hi + 1.0
    one_lo = 0.0 * xi_hi

    p3_hi, p3_lo = zero_hi, zero_lo
    p2_hi, p2_lo = zero_hi - 1.0, zero_lo
    p1_hi, p1_lo = xi_hi.copy() if hasattr(xi_hi, "copy") else xi_hi, \
                   xi_lo.copy() if hasattr(xi_lo, "copy") else xi_lo

    p3p_hi, p3p_lo = zero_hi, zero_lo
    p2p_hi, p2p_lo = zero_hi, zero_lo
    p1p_hi, p1p_lo = one_hi, one_lo

    L_hi = []
    L_lo = []
    Lp_hi = []
    Lp_lo = []

    for j in range(2, p + 2):
        # DD coefficients: a = (2j-3)/j, b = (j-3)/j
        # (numerator and denominator are exact small integers; divide in DD)
        a_num_hi = 0.0 * xi_hi + float(2 * j - 3)
        a_num_lo = 0.0 * xi_lo
        b_num_hi = 0.0 * xi_hi + float(j - 3)
        b_num_lo = 0.0 * xi_lo
        j_hi = 0.0 * xi_hi + float(j)
        j_lo = 0.0 * xi_lo
        a_hi, a_lo = dd_div(a_num_hi, a_num_lo, j_hi, j_lo)
        b_hi, b_lo = dd_div(b_num_hi, b_num_lo, j_hi, j_lo)

        # Save current (new p2 becomes p1 of next iter)
        p3_new_hi, p3_new_lo = p2_hi, p2_lo
        p2_new_hi, p2_new_lo = p1_hi, p1_lo
        p3p_new_hi, p3p_new_lo = p2p_hi, p2p_lo
        p2p_new_hi, p2p_new_lo = p1p_hi, p1p_lo

        # new_p1 = a * xi * p2_new - b * p3_new
        xi_p2_hi, xi_p2_lo = dd_mul(xi_hi, xi_lo, p2_new_hi, p2_new_lo)
        a_xi_p2_hi, a_xi_p2_lo = dd_mul(a_hi, a_lo, xi_p2_hi, xi_p2_lo)
        b_p3_hi, b_p3_lo = dd_mul(b_hi, b_lo, p3_new_hi, p3_new_lo)
        new_p1_hi, new_p1_lo = dd_sub(a_xi_p2_hi, a_xi_p2_lo,
                                       b_p3_hi, b_p3_lo)

        # new_p1p = a (p2_new + xi p2p_new) - b p3p_new
        xi_p2p_hi, xi_p2p_lo = dd_mul(xi_hi, xi_lo, p2p_new_hi, p2p_new_lo)
        sum_hi, sum_lo = dd_add(p2_new_hi, p2_new_lo, xi_p2p_hi, xi_p2p_lo)
        a_sum_hi, a_sum_lo = dd_mul(a_hi, a_lo, sum_hi, sum_lo)
        b_p3p_hi, b_p3p_lo = dd_mul(b_hi, b_lo, p3p_new_hi, p3p_new_lo)
        new_p1p_hi, new_p1p_lo = dd_sub(a_sum_hi, a_sum_lo,
                                          b_p3p_hi, b_p3p_lo)

        L_hi.append(new_p1_hi)
        L_lo.append(new_p1_lo)
        Lp_hi.append(new_p1p_hi)
        Lp_lo.append(new_p1p_lo)

        # Shift for next iteration
        p3_hi, p3_lo = p3_new_hi, p3_new_lo
        p2_hi, p2_lo = p2_new_hi, p2_new_lo
        p1_hi, p1_lo = new_p1_hi, new_p1_lo

        p3p_hi, p3p_lo = p3p_new_hi, p3p_new_lo
        p2p_hi, p2p_lo = p2p_new_hi, p2p_new_lo
        p1p_hi, p1p_lo = new_p1p_hi, new_p1p_lo

    return L_hi, L_lo, Lp_hi, Lp_lo


def dd_legendre(xi_hi, xi_lo, p_max):
    """DD Legendre P_0..P_p_max. Returns lists of DD tuples."""
    P_hi = [0.0 * xi_hi + 1.0]  # P_0 = 1
    P_lo = [0.0 * xi_lo]
    if p_max >= 1:
        P_hi.append(xi_hi)
        P_lo.append(xi_lo)
    for k in range(2, p_max + 1):
        # DD coefficients: inv_k = 1/k, coef1 = (2k-1)/k, coef2 = -(k-1)/k
        one_hi = 0.0 * xi_hi + 1.0
        one_lo = 0.0 * xi_lo
        k_hi = 0.0 * xi_hi + float(k)
        k_lo = 0.0 * xi_lo
        inv_k_hi, inv_k_lo = dd_div(one_hi, one_lo, k_hi, k_lo)
        # coef1 = 2 - inv_k
        two_hi = 0.0 * xi_hi + 2.0
        two_lo = 0.0 * xi_lo
        coef1_hi, coef1_lo = dd_sub(two_hi, two_lo, inv_k_hi, inv_k_lo)
        # coef2 = inv_k - 1
        coef2_hi, coef2_lo = dd_sub(inv_k_hi, inv_k_lo, one_hi, one_lo)

        # P[k] = coef1 * xi * P[k-1] + coef2 * P[k-2]
        xi_Pk1_hi, xi_Pk1_lo = dd_mul(xi_hi, xi_lo, P_hi[k - 1], P_lo[k - 1])
        c1_xi_Pk1_hi, c1_xi_Pk1_lo = dd_mul(coef1_hi, coef1_lo,
                                              xi_Pk1_hi, xi_Pk1_lo)
        c2_Pk2_hi, c2_Pk2_lo = dd_mul(coef2_hi, coef2_lo,
                                        P_hi[k - 2], P_lo[k - 2])
        Pk_hi, Pk_lo = dd_add(c1_xi_Pk1_hi, c1_xi_Pk1_lo,
                                c2_Pk2_hi, c2_Pk2_lo)
        P_hi.append(Pk_hi)
        P_lo.append(Pk_lo)
    return P_hi, P_lo


def main():
    """Validate DD Legendre against mpmath 60-digit reference at sample points."""
    import mpmath as mp
    mp.mp.dps = 60

    print("=" * 72)
    print("DD Legendre polynomial validation vs mpmath 60-digit")
    print("=" * 72)

    rng = np.random.default_rng(45)
    test_pts = [-0.7, -0.3, 0.1, 0.4, 0.85]

    for p in [2, 3, 4, 5]:
        max_err_L = 0.0
        max_err_Lp = 0.0
        max_err_P = 0.0
        for xi in test_pts:
            xi_arr = np.array([xi])
            xi_hi = xi_arr
            xi_lo = np.zeros_like(xi_arr)
            L_hi, L_lo, Lp_hi, Lp_lo = dd_integrated_legendre_with_deriv(
                xi_hi, xi_lo, p)
            P_hi, P_lo = dd_legendre(xi_hi, xi_lo, p + 1)

            # mpmath reference
            xi_mp = mp.mpf(xi)
            p3 = mp.mpf(0); p2 = mp.mpf(-1); p1 = xi_mp
            p3p = mp.mpf(0); p2p = mp.mpf(0); p1p = mp.mpf(1)
            L_mp = []; Lp_mp = []
            for j in range(2, p + 2):
                a = mp.mpf(2 * j - 3) / mp.mpf(j)
                b = mp.mpf(j - 3) / mp.mpf(j)
                p3_new = p2; p2_new = p1
                p3p_new = p2p; p2p_new = p1p
                new_p1 = a * xi_mp * p2_new - b * p3_new
                new_p1p = a * (p2_new + xi_mp * p2p_new) - b * p3p_new
                L_mp.append(new_p1)
                Lp_mp.append(new_p1p)
                p3 = p3_new; p2 = p2_new; p1 = new_p1
                p3p = p3p_new; p2p = p2p_new; p1p = new_p1p

            P_mp = [mp.mpf(1), xi_mp]
            for k in range(2, p + 2):
                inv_k = mp.mpf(1) / mp.mpf(k)
                Pk = (mp.mpf(2) - inv_k) * xi_mp * P_mp[k - 1] \
                     + (inv_k - mp.mpf(1)) * P_mp[k - 2]
                P_mp.append(Pk)

            for j in range(p):
                dd_val = mp.mpf(float(L_hi[j][0])) + mp.mpf(float(L_lo[j][0]))
                err = abs(dd_val - L_mp[j]) / max(abs(L_mp[j]), mp.mpf("1e-300"))
                max_err_L = max(max_err_L, float(err))

                ddp_val = mp.mpf(float(Lp_hi[j][0])) + mp.mpf(float(Lp_lo[j][0]))
                errp = abs(ddp_val - Lp_mp[j]) / max(abs(Lp_mp[j]), mp.mpf("1e-300"))
                max_err_Lp = max(max_err_Lp, float(errp))

            for k in range(p + 2):
                dd_val = mp.mpf(float(P_hi[k][0])) + mp.mpf(float(P_lo[k][0]))
                err = abs(dd_val - P_mp[k]) / max(abs(P_mp[k]), mp.mpf("1e-300"))
                max_err_P = max(max_err_P, float(err))

        print(f"\n  order p={p}:")
        print(f"    max rel error L_{{j+2}}(xi)     = {max_err_L:.3e}")
        print(f"    max rel error L_{{j+2}}'(xi)    = {max_err_Lp:.3e}")
        print(f"    max rel error P_k(xi)         = {max_err_P:.3e}")


if __name__ == "__main__":
    main()
