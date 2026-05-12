"""cylinder_vim_axisym.py — Axisymmetric Volume Integral Method for Cu disk.

Cu cylinder R=10 mm, t=2 mm, sigma=5.8e7 S/m, uniform B_z=1 T applied.

Formulation:
  - Discretise conductor cross-section [0, R] x [-t/2, t/2] into (Nr, Nz) cells
  - Piecewise constant J_phi(r, z) basis (one DoF per cell)
  - Mass matrix M_ii = volume of cell i with 2 pi r Jacobian
  - Kernel matrix K_ij = mu0/(4 pi) * double integral of axisymmetric Green's
    function G_axisym(r, z; r', z') * 2 pi r * 2 pi r' over cells i and j
  - G_axisym uses complete elliptic integrals K(m), E(m) (scipy.special)
  - Generalized eigenproblem K v = lambda M v
  - tau_n = sigma * lambda_n  (Foster pole time constants)
  - Foster amplitudes g_n = v_n . b where b_i = int_cell_i A_ext . phi_i dV
    For uniform B_z: A_ext_phi = (B0/2) r, so b_i = int (B0/2) r * 2 pi r dr dz

Output: Foster spectrum -> Hankel-Pade -> Cauer-I rungs (Kameari conv)
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


# === Geometry & material =======================================================
R_DISK = 10e-3
T_DISK = 2e-3
SIGMA = 5.8e7
MU0 = 4 * math.pi * 1e-7
B0 = 1.0


def axisym_kernel(r, z, r_p, z_p):
    """Bare axisymmetric Green's function (without mu0/(4 pi) factor):
       A_phi(r, z) from unit ring source at (r_p, z_p) is
       (mu0/pi) * sqrt(r_p/r) * ((2/k - k) K(m) - (2/k) E(m))
       where m = k^2 = 4 r r_p / ((r + r_p)^2 + (z - z_p)^2).
       We return only the bare bracket (no mu0/pi).
    """
    if r < 1e-15 or r_p < 1e-15:
        return 0.0
    d2 = (r + r_p) ** 2 + (z - z_p) ** 2
    m = 4.0 * r * r_p / d2     # scipy uses m parameter
    if m >= 1.0 - 1e-12:
        # Singular limit (collocation near identical) — return finite cap
        # by using m close to 1
        m = 1.0 - 1e-12
    k = math.sqrt(m)
    K = sp.ellipk(m)
    E = sp.ellipe(m)
    coef = math.sqrt(r_p / r) / math.pi
    return coef * ((2.0 / k - k) * K - (2.0 / k) * E)


def assemble_K_M_b(Nr, Nz, n_quad_inner=2):
    """Assemble K, M, b for piecewise-constant basis on (Nr, Nz) grid.

    M_ij = delta_ij * (volume of cell i with 2 pi r weight)
    K_ij = (mu0 / 4 pi) * integral_i integral_j G * 2 pi r * 2 pi r' dV dV'
         (4-point Gauss midpoint approx in each cell for outer/inner)
    b_i  = (B0 / 2) * integral_cell_i r * 2 pi r dr dz  (J_phi expansion of A_ext)
         = (B0 / 2) * 2 pi * (r_b^3 - r_a^3) / 3 * (z_b - z_a)
    """
    n = Nr * Nz
    dr = R_DISK / Nr
    dz = T_DISK / Nz

    # Cell coords
    cells = []
    for i in range(Nr):
        for j in range(Nz):
            r_a, r_b = i * dr, (i + 1) * dr
            z_a, z_b = -T_DISK / 2 + j * dz, -T_DISK / 2 + (j + 1) * dz
            cells.append((r_a, r_b, z_a, z_b, (r_a + r_b) / 2, (z_a + z_b) / 2))

    # Mass matrix (diagonal) — int 2 pi r dr dz over cell
    M = np.zeros((n, n))
    for idx, (r_a, r_b, z_a, z_b, rc, zc) in enumerate(cells):
        M[idx, idx] = math.pi * (r_b ** 2 - r_a ** 2) * (z_b - z_a)

    # b vector — int (B0/2) r * 2 pi r dr dz = (B0/2) * (2 pi/3) * (r_b^3 - r_a^3) * (z_b - z_a)
    b_vec = np.zeros(n)
    for idx, (r_a, r_b, z_a, z_b, _, _) in enumerate(cells):
        b_vec[idx] = (B0 / 2.0) * (2 * math.pi / 3) * (r_b ** 3 - r_a ** 3) * (z_b - z_a)

    # K matrix — double volume integral
    # Use Gauss-Legendre quadrature in each cell (n_quad_inner points per dim)
    nodes_m1, weights_m1 = np.polynomial.legendre.leggauss(n_quad_inner)
    nodes_01 = 0.5 * (nodes_m1 + 1.0)   # [0, 1]
    weights_01 = 0.5 * weights_m1

    K = np.zeros((n, n))
    for i, (r_a, r_b, z_a, z_b, _, _) in enumerate(cells):
        # quadrature points in cell i
        r_qi = np.array([r_a + (r_b - r_a) * nq for nq in nodes_01])
        z_qi = np.array([z_a + (z_b - z_a) * nq for nq in nodes_01])
        wr_i = np.array([(r_b - r_a) * w for w in weights_01])
        wz_i = np.array([(z_b - z_a) * w for w in weights_01])

        for j, (r_a2, r_b2, z_a2, z_b2, _, _) in enumerate(cells):
            r_qj = np.array([r_a2 + (r_b2 - r_a2) * nq for nq in nodes_01])
            z_qj = np.array([z_a2 + (z_b2 - z_a2) * nq for nq in nodes_01])
            wr_j = np.array([(r_b2 - r_a2) * w for w in weights_01])
            wz_j = np.array([(z_b2 - z_a2) * w for w in weights_01])

            # double sum over quadrature points
            # Galerkin: K_ij = (mu0 / 2) * 2pi * int int phi_i G phi_j r_o dr_o dz_o dr' dz'
            # Note: observation side has 2 pi r_o Jacobian, source side has plain dr' dz'
            # (the ring current is J_phi * dr' * dz', no extra 2 pi r').
            kij = 0.0
            for ar, ri in enumerate(r_qi):
                for az, zi in enumerate(z_qi):
                    for br, rj in enumerate(r_qj):
                        for bz, zj in enumerate(z_qj):
                            G = axisym_kernel(ri, zi, rj, zj)
                            kij += (G * (2 * math.pi * ri)
                                    * wr_i[ar] * wz_i[az]
                                    * wr_j[br] * wz_j[bz])
            K[i, j] = (MU0 / 2.0) * kij

    # Symmetrize (numerical)
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
    Nr, Nz = 48, 12   # 576 cells = 576 modes (4x refinement from 144)
    print(f"Cylinder VIM axisym: Nr={Nr}, Nz={Nz}, n_dofs={Nr*Nz}")
    print(f"R={R_DISK*1000} mm, t={T_DISK*1000} mm, sigma={SIGMA:.2e}")
    print()

    print("Assembling K, M, b matrices...")
    import time
    t0 = time.time()
    K, M, b = assemble_K_M_b(Nr, Nz, n_quad_inner=2)
    print(f"  done in {time.time() - t0:.1f}s")
    print(f"  K symmetry: max|K - K^T| = {np.max(np.abs(K - K.T)):.2e}")

    # Generalized eigenvalue: K v = lambda M v
    print("Solving generalized eigenproblem...")
    eigvals, eigvecs = eigh(K, M)
    # M-normalize
    for k in range(len(eigvals)):
        v = eigvecs[:, k]
        nrm = math.sqrt(float(v @ M @ v))
        if nrm > 1e-30:
            eigvecs[:, k] = v / nrm

    # tau_n = sigma * lambda_n
    tau_us = np.where(eigvals > 1e-30, 1e6 * SIGMA * eigvals, 0.0)
    g = eigvecs.T @ b
    g2 = g * g

    # Sort by |g|^2
    order = np.argsort(-g2)
    print(f"  Top 10 modes by |g|^2:")
    for rank in range(min(10, len(order))):
        idx = order[rank]
        print(f"    rank={rank}: tau={tau_us[idx]:.4f} us, "
              f"g2={g2[idx]:.4e}")

    # Foster-to-Kameari Cauer-I extraction
    n_use = min(50, np.sum(g2 > 1e-30))
    print(f"\nKameari Cauer-I extraction (top {n_use} modes, mpmath 80 dig):")
    sorted_idx = order[:n_use]
    tau_seconds = [mp.mpf(float(tau_us[i])) * mp.mpf("1e-6") for i in sorted_idx]
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

    # Save JSON
    out = {
        "method": "Cylinder axisymmetric VIM (Python, axisym Green's function)",
        "params": {
            "R_mm": R_DISK * 1000, "t_mm": T_DISK * 1000,
            "sigma_S_per_m": SIGMA, "B0_T": B0,
            "Nr": Nr, "Nz": Nz, "n_dofs": Nr * Nz,
            "n_quad_inner": 2,
        },
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
                / "cylinder_vim_axisym_results.json")
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
