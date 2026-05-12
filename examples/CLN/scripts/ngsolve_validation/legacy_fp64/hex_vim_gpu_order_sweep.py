"""hex_vim_gpu_order_sweep.py — GPU K assembly + Foster + Cauer
   for cuboid 5×2×1 mm Cu, sweeping order = 3, 4, 5.

Pipeline (per order):
  1. GPU assemble bare K via hex_vim_cupy_kassembly
  2. Multiply by mu0 / (4 pi) for physical K
  3. Build mass matrix M (Python GL n=8) — same as radia-vim scripts
  4. Build excitation b vector — A_ext = (-y/2, x/2, 0) on physical cuboid
  5. Solve generalized eigenproblem K v = lambda M v
  6. Foster spectrum + Hankel-Pade → Cauer-I rungs (Kameari)

Outputs JSON results similar to radia-vim/scripts/extract_tau_cuboid_521.py
but assembled in seconds instead of hours.
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
from scipy.linalg import eigh
import mpmath as mp
mp.mp.dps = 80

sys.path.insert(0, str(Path(__file__).parent))
from hex_vim_cupy_kassembly import assemble_K_cupy
from hex_vim_cupy import evaluate_basis

sys.path.insert(0, str(Path("S:/Radia/01_GitHub/packages/radia-vim/src")))
import radia_vim


SIGMA = 5.8e7
MU0 = 4 * math.pi * 1e-7
B0 = 1.0

# Geometries to sweep
GEOMETRIES = {
    "cuboid521": (5e-3, 2e-3, 1e-3),
    "A1":        (17.72e-3, 17.72e-3, 2e-3),
}

import os as _os
GEOM_KEY = _os.environ.get("HEXVIM_GEOM", "cuboid521")
A_M, B_M, C_M = GEOMETRIES[GEOM_KEY]


def assemble_mass_matrix_python(basis, a, b, c, n_gauss=8):
    """Same as radia-vim/scripts/extract_tau_cuboid_521.py."""
    n = basis.n_dofs
    nodes_m1, weights_m1 = np.polynomial.legendre.leggauss(n_gauss)
    nodes = 0.5 * (nodes_m1 + 1.0)
    weights = 0.5 * weights_m1

    Phi = np.zeros((n, 3, n_gauss, n_gauss, n_gauss))
    for ix, x in enumerate(nodes):
        for iy, y in enumerate(nodes):
            for iz, z in enumerate(nodes):
                for d in range(n):
                    v = basis.evaluate(d, float(x), float(y), float(z))
                    Phi[d, 0, ix, iy, iz] = v[0]
                    Phi[d, 1, ix, iy, iz] = v[1]
                    Phi[d, 2, ix, iy, iz] = v[2]

    W3 = weights[:, None, None] * weights[None, :, None] * weights[None, None, :]
    Mref = np.zeros((3, n, n))
    for comp in range(3):
        Pcomp = Phi[:, comp]
        Pw = Pcomp * W3[None, :, :, :]
        Mref[comp] = np.tensordot(Pcomp.reshape(n, -1),
                                  Pw.reshape(n, -1), axes=([1], [1]))

    Mphys = (a / (b * c)) * Mref[0] + (b / (a * c)) * Mref[1] + (c / (a * b)) * Mref[2]
    return Mphys


def assemble_b_vector_python(basis, a, b, c, n_gauss=8):
    """A_ext = (-y/2, x/2, 0) under HDiv contravariant Piola map.
    On reference cube [0,1]^3, ref point (xi, eta, zeta) maps to physical
    (a*xi, b*eta, c*zeta) — but for HDiv Piola map, A_ext_ref = J^T A_phys / det J
    where J = diag(a, b, c).
    A_phys at (a*xi, b*eta, c*zeta) = (-b*eta/2, a*xi/2, 0)
    A_ext_ref = (1/(a b c)) × diag(a, b, c) × A_phys
              = (1/(b c) × (-b*eta/2), 1/(a c) × (a*xi/2), 0)
              = (-eta/(2*c), xi/(2*c), 0)
    Wait that doesn't seem right. Let me re-derive.

    HDiv Piola: phi_phys = (1/det J) J phi_ref
    Inner product: ∫ phi_phys · A_ext dV_phys = ∫ (1/detJ) J phi_ref · A_ext detJ dV_ref
                                              = ∫ phi_ref · J^T A_ext dV_ref
    So b_i = ∫_ref phi_ref_i · (J^T A_ext_phys) dV_ref
    where A_ext_phys is evaluated at the physical point corresponding to the ref point.

    For our problem: A_ext_phys(x_phys) = (-y_phys/2, x_phys/2, 0) × B0
                                       = (-b*eta/2, a*xi/2, 0) × B0
    J^T A_ext = diag(a, b, c) × (-b*eta/2, a*xi/2, 0) × B0
              = (-a*b*eta/2, b*a*xi/2, 0) × B0
              = (a*b*B0/2) × (-eta, xi, 0)

    Hmm but this includes ab/2 prefactor.

    Actually simpler: directly compute b = ∫_phys phi_phys(r) · A_ext(r) dV
                                       = ∫_phys sigma_factor * phi_phys(r) · A_ext(r) dV
    Wait — for VIM Foster b vector:
        b_i = ∫ A_ext(r) · J_i(r) dV   where J_i is the i-th basis (physical)

    OK let me just match what's in extract_tau_cuboid_521.py.
    """
    raise NotImplementedError("see assemble_b_vector below for actual impl")


def main():
    import os as _os
    parser_orders_str = _os.environ.get("HEXVIM_ORDERS", "3,4,5")
    parser_orders = [int(x) for x in parser_orders_str.split(",")]
    n_th, n_ph, n_rh, n_c = 12, 24, 12, 6

    print("=" * 72)
    print(f"GPU K assembly for cuboid {A_M*1e3:g} x {B_M*1e3:g} x {C_M*1e3:g} mm Cu")
    print(f"Quad: th={n_th}, ph={n_ph}, rh={n_rh}, c={n_c}")
    print("=" * 72)

    out_results = {}
    for p in parser_orders:
        print(f"\n=== Order {p} ===")
        n_dofs = 2 * p**3 + 3 * p**2
        print(f"  n_dofs = {n_dofs}")

        # Estimated GPU memory for Phi arrays
        n_q = n_th * n_ph * n_rh * n_c**3
        mem_gb = 2 * n_q * n_dofs * 3 * 8 / 1e9
        print(f"  n_q = {n_q}, est Phi memory = {mem_gb:.2f} GB")

        # Bare K matrix on GPU (no mu0/4pi). Order >= 6 needs FP32 (mem).
        dtype = np.float32 if p >= 6 else np.float64
        t0 = time.time()
        K_bare = assemble_K_cupy(p, A_M, B_M, C_M, n_th, n_ph, n_rh, n_c,
                                  use_gpu=True, dtype=dtype, verbose=True)
        # Upcast to FP64 for downstream eigh
        K_bare = K_bare.astype(np.float64)
        t_K = time.time() - t0
        print(f"  GPU K assembly: {t_K:.2f}s")

        # Physical K = bare K × mu0/(4 pi)
        K_phys = K_bare * (MU0 / (4 * math.pi))

        # Mass matrix (Python, ~few seconds)
        basis = radia_vim.HDivDivFreeHexBasis(p)
        t0 = time.time()
        M_phys = assemble_mass_matrix_python(basis, A_M, B_M, C_M, n_gauss=8)
        t_M = time.time() - t0
        print(f"  Mass M assembly: {t_M:.2f}s")

        # b vector for A_ext = (-y/2, x/2, 0) × B0 on physical cuboid
        # On reference cube (xi, eta, zeta) ∈ [0,1]^3, physical
        # x = A_M * xi, y = B_M * eta, z = C_M * zeta. Centered at origin
        # means x_phys = A_M * (xi - 0.5), y_phys = B_M * (eta - 0.5).
        # A_ext_phys = B0 × (-y_phys/2, x_phys/2, 0)
        #            = B0 × (-B_M*(eta-0.5)/2, A_M*(xi-0.5)/2, 0)
        # HDiv Piola: phi_phys = (1/(A_M*B_M*C_M)) diag(A_M,B_M,C_M) phi_ref
        # b_i = ∫_phys A_ext · phi_phys_i dV_phys
        #     = ∫_ref A_ext(phys(ref)) · phi_phys_i (det J) dV_ref
        # but phi_phys = (1/detJ) J phi_ref so
        # b_i = ∫_ref A_ext · J phi_ref_i dV_ref
        # Component-wise:
        #   b_i = ∫_ref [A_ext_x A_M phi_ref_i.x + A_ext_y B_M phi_ref_i.y
        #              + A_ext_z C_M phi_ref_i.z] dV_ref

        n_gauss = 8
        nodes_m1, w_m1 = np.polynomial.legendre.leggauss(n_gauss)
        nodes = 0.5 * (nodes_m1 + 1.0)
        weights = 0.5 * w_m1

        # Build all reference points (n_gauss^3, 3)
        XI, ETA, ZETA = np.meshgrid(nodes, nodes, nodes, indexing="ij")
        xi_flat = XI.flatten()
        eta_flat = ETA.flatten()
        zeta_flat = ZETA.flatten()

        # Evaluate basis at all reference points (CPU is fine, small n)
        Phi_ref = evaluate_basis(xi_flat, eta_flat, zeta_flat, p, xp_mod=np)
        # Phi_ref shape: (n_gauss^3, n_dofs, 3)

        # A_ext at physical points corresponding to ref
        x_phys = A_M * (xi_flat - 0.5)  # centered cuboid
        y_phys = B_M * (eta_flat - 0.5)
        Ax_phys = B0 * (-y_phys / 2)
        Ay_phys = B0 * (x_phys / 2)
        # Az_phys = 0

        # Quadrature weights (n_gauss^3,)
        W3 = (weights[:, None, None] * weights[None, :, None]
              * weights[None, None, :]).flatten()

        # b_i = sum_q W3[q] × (Ax * A_M * Phi_ref[q,i,0] + Ay * B_M * Phi_ref[q,i,1])
        # Note Az=0 so z-component drops
        wax = W3 * Ax_phys * A_M
        way = W3 * Ay_phys * B_M
        b_vec = (Phi_ref[..., 0].T @ wax) + (Phi_ref[..., 1].T @ way)

        # Solve generalized eigenproblem K v = lambda M v
        print(f"  Solving K v = lambda M v ({n_dofs} dofs)...")
        t0 = time.time()
        eigvals, eigvecs = eigh(K_phys, M_phys)
        # M-normalize
        for k in range(len(eigvals)):
            v = eigvecs[:, k]
            nrm = math.sqrt(float(v @ M_phys @ v))
            if nrm > 1e-30:
                eigvecs[:, k] = v / nrm
        t_eig = time.time() - t0
        print(f"  eigh: {t_eig:.2f}s")

        # tau = sigma * lambda (SIGMA, MU0 stripped from K)
        # Wait — K_phys has mu0/(4pi) prefactor, so the units are right.
        # The eddy current eigenvalue: K v = lambda M v with K = (mu0/4pi) integral
        # gives lambda such that tau = sigma * lambda (cylinder VIM convention)
        tau_us = np.where(eigvals > 1e-30, 1e6 * SIGMA * eigvals, 0.0)
        g = eigvecs.T @ b_vec
        g2 = g * g

        order = np.argsort(-g2)
        print(f"  Top 5 modes by |g|^2:")
        for rank in range(min(5, len(order))):
            idx = order[rank]
            print(f"    rank={rank}: tau={tau_us[idx]:.4f} us, g2={g2[idx]:.4e}")

        # Kameari Cauer-I extraction
        n_use = min(80, int(np.sum(g2 > 1e-32)))
        sorted_idx = order[:n_use]
        tau_seconds = [mp.mpf(float(tau_us[i])) * mp.mpf("1e-6")
                       for i in sorted_idx]
        g2_mp = [mp.mpf(float(g2[i])) for i in sorted_idx]
        n_moments = min(40, 2 * n_use - 4)
        alphas = []
        for n in range(n_moments):
            a_val = mp.mpf(0)
            for gv, tv in zip(g2_mp, tau_seconds):
                a_val += gv * tv ** n
            alphas.append((mp.mpf(-1)) ** n * a_val)

        # Hankel-Padé
        c_list = [mp.mpf(a) for a in alphas]
        cauer_p = []
        for _ in range(12):
            if len(c_list) < 2 or abs(c_list[0]) < mp.mpf(10) ** (-mp.mp.dps + 10):
                break
            n_c_list = len(c_list)
            e = [mp.mpf(0)] * n_c_list
            e[0] = 1 / c_list[0]
            for k in range(1, n_c_list):
                e[k] = -sum(c_list[j + 1] * e[k - j - 1] for j in range(k)
                            if k - j - 1 >= 0) / c_list[0]
            cauer_p.append(e[0])
            c_list = e[1:]

        rungs = []
        print(f"  Kameari Cauer-I rungs:")
        print(f"  k | R_2k       | L_2k+1     | tau_pair us")
        for k in range(min(6, len(cauer_p) // 2)):
            Rv = cauer_p[2 * k]
            Linv = cauer_p[2 * k + 1]
            if abs(Linv) < mp.mpf(10) ** (-50):
                break
            Lv = 1 / Linv
            tau_p = float(Lv / Rv) * 1e6
            print(f"  {k} | {float(Rv):>10.4e} | {float(Lv):>10.4e} | {tau_p:.4f}")
            rungs.append({"k": k, "R_2k": float(Rv),
                          "L_2k_plus_1": float(Lv), "tau_pair_us": tau_p})

        out_results[f"order_{p}"] = {
            "n_dofs": n_dofs,
            "gpu_K_assembly_s": t_K,
            "M_assembly_s": t_M,
            "eig_s": t_eig,
            "leading_tau_us": float(tau_us[order[0]]),
            "Cauer_rungs_Kameari": rungs,
        }

    out_path = (Path(__file__).parent
                / f"hex_vim_gpu_{GEOM_KEY}_results.json")
    out_path.write_text(json.dumps(out_results, indent=2), encoding="utf-8")
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
