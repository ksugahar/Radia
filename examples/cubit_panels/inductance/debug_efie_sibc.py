"""
Debug EFIE-SIBC: diagnoses why BEM EFIE-SIBC gives wrong results for finite Z_s.

Findings:
=========

1. FEM BUG (FIXED): The scattered-field FEM formulation was missing the incident
   field contribution. The correct RHS is:
     f(v) = -(jw/Z_s)*<A_inc, v>_sibc  -  <n x H_inc, v>_sibc
   The second term (H_inc boundary) was missing.  Equivalently, using the total
   field decomposition A = A_inc + A_scat:
     f(v) = -a(A_inc, v) = -curl-curl(A_inc, v) - sibc(A_inc, v)
   Fix: add  f_lf += (-1) * nxH_inc * v.Trace() * ds("sibc")  [boundary form]
   or equivalently  f_lf += (-nu0) * curl_A_inc * curl(v) * dx  [volume form]

2. BEM BUG (OPEN): The EFIE-SIBC  Z_s*J + jw*mu0*SL(J) = -jw*A_inc  uses the
   representation A_scat = mu0*SL(J_s), but the Laplace SL eigenvalue for l=1
   on a sphere is R/3 (not R).  This gives denominator (3*Z_s + jw*mu0*R)
   instead of the correct (Z_s + jw*mu0*R).
   The fix requires a different integral equation (MFIE, PMCHWT, or
   Stratton-Chu with both SL and DL).

Usage:
    python debug_efie_sibc.py
"""

import math
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

MU_0 = 4e-7 * np.pi


def analytical_H_t_rms(R, Z_s, B0, omega):
    """Exact analytical SIBC H_t_rms on sphere."""
    beta = 1j * omega * MU_0 * R
    J_max = abs((3.0 / 2.0) * 1j * omega * B0 * R / (beta + Z_s))
    return J_max * math.sqrt(2.0 / 3.0)


def bem_predicted_H_t_rms(R, Z_s, B0, omega):
    """Predicted BEM result from SL eigenvalue = R/3."""
    J_max = abs(1j * omega * B0 * R / (2.0 * (Z_s + 1j * omega * MU_0 * R / 3.0)))
    return J_max * math.sqrt(2.0 / 3.0)


if __name__ == "__main__":
    R = 0.01
    B0 = 0.001
    freq = 1000
    omega = 2 * np.pi * freq
    jw_mu0_R = omega * MU_0 * R

    print("=" * 72)
    print("EFIE-SIBC Diagnostic: BEM formula vs Analytical")
    print("=" * 72)
    print(f"R = {R*1e3:.1f} mm, B0 = {B0*1e3:.1f} mT, f = {freq} Hz")
    print(f"|jw*mu0*R| = {jw_mu0_R:.4e}")
    print()
    print("BEM EFIE: J_0 = jw*B0*R / (2*(Z_s + jw*mu0*R/3))")
    print("Exact:    J_0 = (3/2)*jw*B0*R / (jw*mu0*R + Z_s)")
    print("Ratio:    |jw*mu0*R + Z_s| / |3*Z_s + jw*mu0*R|")
    print()

    from verify_sphere_sibc import bem_sphere_sibc

    print(f"{'ratio':>8s} {'Analytical':>12s} {'BEM':>12s} {'BEM_pred':>12s} "
          f"{'BEM/Ana':>8s} {'Pred/Ana':>8s}")
    print("-" * 72)
    for ratio in [0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]:
        Zs = ratio * jw_mu0_R * (1 + 1j) / np.sqrt(2)
        H_ana = analytical_H_t_rms(R, Zs, B0, omega)
        H_pred = bem_predicted_H_t_rms(R, Zs, B0, omega)
        bem = bem_sphere_sibc(R, Zs, B0, omega)
        H_bem = bem['H_t_rms']
        r_b = H_bem / H_ana if H_ana > 0 else 0
        r_p = H_pred / H_ana if H_ana > 0 else 0
        print(f"{ratio:8.3f} {H_ana:12.4f} {H_bem:12.4f} {H_pred:12.4f} "
              f"{r_b:8.4f} {r_p:8.4f}")
