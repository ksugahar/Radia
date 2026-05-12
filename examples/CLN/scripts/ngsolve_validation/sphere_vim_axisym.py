"""sphere_vim_axisym.py — Axisymmetric VIM for Cu sphere.

Cu sphere R=10 mm, sigma=5.8e7 S/m, uniform B_z=1 T applied.

Formulation: identical to cylinder_vim_axisym.py but on a (rho, theta) grid
covering the sphere interior rho in [0, R], theta in [0, pi].

  - Piecewise constant J_phi(rho, theta) basis (one DoF per cell)
  - Mass matrix M_ii = volume of revolution =
      2 pi * (rho_b^3 - rho_a^3)/3 * (cos(theta_a) - cos(theta_b))
  - Kernel matrix K_ij = (mu0/2) * double integral of axisymmetric
    Green's function G_axisym(r, z; r', z') * 2 pi r over cells i, j
  - G_axisym uses complete elliptic integrals K(m), E(m) (scipy.special)
  - Generalized eigenproblem K v = lambda M v
  - tau_n = sigma * lambda_n
  - Foster amplitudes g_n = v_n . b   with
      b_i = pi B0 * (rho_b^4 - rho_a^4)/4 *
            [(theta_b - theta_a)/2 - (sin(2*theta_b) - sin(2*theta_a))/4]

Output: Foster spectrum -> Hankel-Pade -> Cauer-I rungs (Kameari conv).
Reference: Stoll Bessel analytical spectrum gives leading tau ~ 694 us.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import scipy.special as sp
from scipy.linalg import eigh
import mpmath as mp

mp.mp.dps = 80


# === Geometry & material =====================================================
R_SPHERE = 10e-3
SIGMA = 5.8e7
MU0 = 4 * math.pi * 1e-7
B0 = 1.0


def axisym_kernel(r, z, r_p, z_p):
    """Bare axisymmetric Green's function (without mu0/(4 pi) factor):
    A_phi(r, z) from unit ring source at (r_p, z_p) is
       (mu0/pi) * sqrt(r_p/r) * ((2/k - k) K(m) - (2/k) E(m))
    with m = k^2 = 4 r r_p / ((r + r_p)^2 + (z - z_p)^2).
    Returns the bare bracket without mu0/pi prefactor.
    """
    if r < 1e-15 or r_p < 1e-15:
        return 0.0
    d2 = (r + r_p) ** 2 + (z - z_p) ** 2
    m = 4.0 * r * r_p / d2
    if m >= 1.0 - 1e-12:
        m = 1.0 - 1e-12
    k = math.sqrt(m)
    K = sp.ellipk(m)
    E = sp.ellipe(m)
    coef = math.sqrt(r_p / r) / math.pi
    return coef * ((2.0 / k - k) * K - (2.0 / k) * E)


def assemble_K_M_b(N_rho, N_theta, n_quad=2):
    """Assemble K, M, b for piecewise-constant basis on (N_rho, N_theta) grid
    over sphere interior rho in [0, R], theta in [0, pi]."""
    n = N_rho * N_theta
    drho = R_SPHERE / N_rho
    dtheta = math.pi / N_theta

    cells = []
    for i in range(N_rho):
        for j in range(N_theta):
            rho_a, rho_b = i * drho, (i + 1) * drho
            th_a, th_b = j * dtheta, (j + 1) * dtheta
            cells.append((rho_a, rho_b, th_a, th_b))

    # Mass matrix (diagonal): vol = 2 pi r dr dz  in (rho, theta):
    #   = 2 pi * (rho_b^3 - rho_a^3)/3 * (cos(th_a) - cos(th_b))
    M = np.zeros((n, n))
    for idx, (ra, rb, ta, tb) in enumerate(cells):
        vol = (2.0 * math.pi * (rb ** 3 - ra ** 3) / 3.0
               * (math.cos(ta) - math.cos(tb)))
        M[idx, idx] = vol

    # b vector: b_i = int A_imposed * 2 pi r dr dz over cell  with
    #   A_imposed = (B0/2) * r,  so b = pi B0 * int r^2 dr dz
    #     = pi B0 * (rho_b^4 - rho_a^4)/4 *
    #         [(tb - ta)/2 - (sin(2 tb) - sin(2 ta))/4]
    b_vec = np.zeros(n)
    for idx, (ra, rb, ta, tb) in enumerate(cells):
        rho_int = (rb ** 4 - ra ** 4) / 4.0
        th_int = ((tb - ta) / 2.0
                  - (math.sin(2.0 * tb) - math.sin(2.0 * ta)) / 4.0)
        b_vec[idx] = math.pi * B0 * rho_int * th_int

    # K matrix: (mu0/2) * int int G(r,z;r',z') * 2 pi r dV dV'
    # Gauss-Legendre quadrature in each cell over (rho, theta)
    nodes_m1, w_m1 = np.polynomial.legendre.leggauss(n_quad)
    nodes_01 = 0.5 * (nodes_m1 + 1.0)
    w_01 = 0.5 * w_m1

    K = np.zeros((n, n))
    for i, (ra_i, rb_i, ta_i, tb_i) in enumerate(cells):
        rho_qi = np.array([ra_i + (rb_i - ra_i) * nq for nq in nodes_01])
        th_qi = np.array([ta_i + (tb_i - ta_i) * nq for nq in nodes_01])
        wrho_i = np.array([(rb_i - ra_i) * w for w in w_01])
        wth_i = np.array([(tb_i - ta_i) * w for w in w_01])
        # (r, z) and Jacobian dr dz = rho * drho * dtheta
        r_qi = np.outer(rho_qi, np.sin(th_qi))   # (n_quad, n_quad)
        z_qi = np.outer(rho_qi, np.cos(th_qi))
        jac_i = rho_qi[:, None] * np.ones_like(th_qi)[None, :]

        for j, (ra_j, rb_j, ta_j, tb_j) in enumerate(cells):
            rho_qj = np.array([ra_j + (rb_j - ra_j) * nq for nq in nodes_01])
            th_qj = np.array([ta_j + (tb_j - ta_j) * nq for nq in nodes_01])
            wrho_j = np.array([(rb_j - ra_j) * w for w in w_01])
            wth_j = np.array([(tb_j - ta_j) * w for w in w_01])
            r_qj = np.outer(rho_qj, np.sin(th_qj))
            z_qj = np.outer(rho_qj, np.cos(th_qj))
            jac_j = rho_qj[:, None] * np.ones_like(th_qj)[None, :]

            kij = 0.0
            for ar in range(n_quad):
                for at in range(n_quad):
                    ri = r_qi[ar, at]
                    zi = z_qi[ar, at]
                    if ri < 1e-15:
                        continue
                    wi = (jac_i[ar, at] * wrho_i[ar] * wth_i[at])
                    for br in range(n_quad):
                        for bt in range(n_quad):
                            rj = r_qj[br, bt]
                            zj = z_qj[br, bt]
                            if rj < 1e-15:
                                continue
                            G = axisym_kernel(ri, zi, rj, zj)
                            wj = (jac_j[br, bt] * wrho_j[br] * wth_j[bt])
                            kij += G * (2.0 * math.pi * ri) * wi * wj
            K[i, j] = (MU0 / 2.0) * kij

    K = 0.5 * (K + K.T)
    return K, M, b_vec


def hankel_pade_cauer(alphas, max_stages):
    c = [mp.mpf(a) for a in alphas]
    p = []
    for _ in range(max_stages):
        if len(c) < 2 or abs(c[0]) < mp.mpf(10) ** (-mp.mp.dps + 10):
            break
        n = len(c)
        e = [mp.mpf(0)] * n
        e[0] = 1 / c[0]
        for k in range(1, n):
            e[k] = -sum(c[j + 1] * e[k - j - 1] for j in range(k)
                        if k - j - 1 >= 0) / c[0]
        p.append(e[0])
        c = e[1:]
    return p


def main():
    N_rho, N_theta = 30, 36   # 1080 cells
    print(f"Sphere VIM axisym: N_rho={N_rho}, N_theta={N_theta}, "
          f"n_dofs={N_rho * N_theta}")
    print(f"R={R_SPHERE * 1000} mm, sigma={SIGMA:.2e}")
    print(f"Stoll analytical leading tau = "
          f"{1e6 * MU0 * SIGMA * R_SPHERE**2 / math.pi**2:.2f} us")
    print()

    print("Assembling K, M, b matrices...")
    import time
    t0 = time.time()
    K, M, b = assemble_K_M_b(N_rho, N_theta, n_quad=2)
    print(f"  done in {time.time() - t0:.1f}s")
    print(f"  K symmetry: max|K - K^T| = {np.max(np.abs(K - K.T)):.2e}")

    print("Solving generalized eigenproblem...")
    eigvals, eigvecs = eigh(K, M)
    for k in range(len(eigvals)):
        v = eigvecs[:, k]
        nrm = math.sqrt(float(v @ M @ v))
        if nrm > 1e-30:
            eigvecs[:, k] = v / nrm

    tau_us = np.where(eigvals > 1e-30, 1e6 * SIGMA * eigvals, 0.0)
    g = eigvecs.T @ b
    g2 = g * g

    order = np.argsort(-g2)
    print(f"  Top 10 modes by |g|^2:")
    for rank in range(min(10, len(order))):
        idx = order[rank]
        print(f"    rank={rank}: tau={tau_us[idx]:.4f} us, "
              f"g2={g2[idx]:.4e}")

    n_use = min(80, int(np.sum(g2 > 1e-32)))
    print(f"\nKameari Cauer-I extraction (top {n_use} modes, mpmath 80 dig):")
    sorted_idx = order[:n_use]
    tau_seconds = [mp.mpf(float(tau_us[i])) * mp.mpf("1e-6")
                   for i in sorted_idx]
    g2_mp = [mp.mpf(float(g2[i])) for i in sorted_idx]

    n_moments = min(40, 2 * n_use - 4)
    alphas = []
    for n in range(n_moments):
        a = mp.mpf(0)
        for gv, tv in zip(g2_mp, tau_seconds):
            a += gv * tv ** n
        alphas.append((mp.mpf(-1)) ** n * a)

    p = hankel_pade_cauer(alphas, 12)
    print(f"  k | R_2k       | L_2k+1     | tau_pair us")
    rungs = []
    for k in range(min(6, len(p) // 2)):
        Rv = p[2 * k]
        Linv = p[2 * k + 1]
        if abs(Linv) < mp.mpf(10) ** (-50):
            break
        Lv = 1 / Linv
        tau_p = float(Lv / Rv) * 1e6
        print(f"  {k} | {float(Rv):>10.4e} | {float(Lv):>10.4e} | {tau_p:.4f}")
        rungs.append({"k": k, "R_2k": float(Rv),
                      "L_2k_plus_1": float(Lv), "tau_pair_us": tau_p})

    out = {
        "method": "Sphere axisymmetric VIM (Python, axisym Green's function)",
        "params": {
            "R_mm": R_SPHERE * 1000, "sigma_S_per_m": SIGMA, "B0_T": B0,
            "N_rho": N_rho, "N_theta": N_theta,
            "n_dofs": N_rho * N_theta, "n_quad_inner": 2,
        },
        "stoll_analytical_tau_us":
            1e6 * MU0 * SIGMA * R_SPHERE ** 2 / math.pi ** 2,
        "leading_tau_us": float(tau_us[order[0]]),
        "top_modes_by_g2": [
            {"rank": rank, "idx": int(order[rank]),
             "tau_us": float(tau_us[order[rank]]),
             "g2": float(g2[order[rank]])}
            for rank in range(min(20, len(order)))
        ],
        "Cauer_rungs_Kameari": rungs,
    }
    out_path = (Path(__file__).parent
                / "sphere_vim_axisym_results.json")
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
