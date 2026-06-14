"""Hex VIM operator assembly; historical CLN extraction retired.

Original method attributions: Hiruma, Sugahara, Kameari.
Original reference retained: Nagamine, FreeFEM++ rectangle/arnoldi.edp.
The reference display below concerns independently stored scalar synthesis.
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import scipy.linalg as la

sys.path.insert(0, str(Path(__file__).parent))
from hex_vim_cupy_kassembly import assemble_K_cupy
from hex_vim_cupy import evaluate_basis

sys.path.insert(0, str(Path("S:/Radia/01_GitHub/packages/radia-vim/src")))
import radia_vim


# Cuboid 5x2x1 mm Cu
A_M = 5e-3
B_M = 2e-3
C_M = 1e-3
SIGMA = 5.8e7
MU0 = 4 * math.pi * 1e-7
B0 = 1.0


def assemble_mass_matrix_python(basis, a, b, c, n_gauss=8):
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


def assemble_b_vector(basis, a, b, c, p, n_gauss=8):
    """b_i = sum_q W3[q] * (Ax_phys * a * Phi_ref[q,i,0] + Ay_phys * b * Phi_ref[q,i,1])"""
    nodes_m1, w_m1 = np.polynomial.legendre.leggauss(n_gauss)
    nodes = 0.5 * (nodes_m1 + 1.0)
    weights = 0.5 * w_m1
    XI, ETA, ZETA = np.meshgrid(nodes, nodes, nodes, indexing="ij")
    xi_flat = XI.flatten()
    eta_flat = ETA.flatten()
    zeta_flat = ZETA.flatten()
    Phi_ref = evaluate_basis(xi_flat, eta_flat, zeta_flat, p, xp_mod=np)
    x_phys = a * (xi_flat - 0.5)
    y_phys = b * (eta_flat - 0.5)
    Ax_phys = B0 * (-y_phys / 2)
    Ay_phys = B0 * (x_phys / 2)
    W3 = (weights[:, None, None] * weights[None, :, None]
          * weights[None, None, :]).flatten()
    wax = W3 * Ax_phys * a
    way = W3 * Ay_phys * b
    return (Phi_ref[..., 0].T @ wax) + (Phi_ref[..., 1].T @ way)




def main():
    p = 4
    n_th, n_ph, n_rh, n_c = 12, 24, 12, 6
    print("=" * 72)
    print(f"Hex VIM cuboid 5x2x1 order={p}: independent operator assembly")
    print("=" * 72)

    # Build K, M, b
    print(f"\nAssembling K, M, b for order={p}...")
    t0 = time.time()
    K_bare = assemble_K_cupy(p, A_M, B_M, C_M, n_th, n_ph, n_rh, n_c,
                              use_gpu=True, verbose=False)
    K_phys = K_bare * (MU0 / (4 * math.pi))
    print(f"  K assembled in {time.time()-t0:.2f}s, shape {K_phys.shape}")

    basis = radia_vim.HDivDivFreeHexBasis(p)
    M = assemble_mass_matrix_python(basis, A_M, B_M, C_M, n_gauss=8)
    b = assemble_b_vector(basis, A_M, B_M, C_M, p, n_gauss=8)
    print(f"  M, b assembled")

    # axifemm convention: K (no sigma), M (with sigma), b (with sigma).
    # tau in seconds directly from lam (no extra sigma factor needed).
    M_sigma = SIGMA * M
    b_sigma = SIGMA * b

    # Reference: QD-Padé Kameari (from existing JSON)
    print(f"\n--- Reference: QD-Padé Kameari τ_pair (existing GPU sweep) ---")
    ref_json = (Path(__file__).parent / "hex_vim_gpu_cuboid521_results.json")
    if ref_json.exists():
        ref = json.loads(ref_json.read_text(encoding="utf-8"))
        for key in ref:
            if "order_4" in key:
                print(f"  {key}:")
                for r in ref[key]["Cauer_rungs_Kameari"][:6]:
                    print(f"    k={r['k']}: tau_pair = {r['tau_pair_us']:.4f}")

    out_path = (Path(__file__).parent
                / "hex_vim_cuboid521_assembly.json")
    out_path.write_text(json.dumps({
        "order": p,
        "n_dofs": basis.n_dofs,
        "K_shape": list(K_phys.shape),
        "M_shape": list(M_sigma.shape),
        "b_size": int(b_sigma.size),
        "status": "assembled; no CLN reduction performed",
    }, indent=2), encoding="utf-8")
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
