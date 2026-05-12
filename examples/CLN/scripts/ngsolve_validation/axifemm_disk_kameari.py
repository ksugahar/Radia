"""axifemm_disk_kameari.py — axifemm Q2 axisym FEM with Kameari Cauer-I extraction.

Existing radia-axifemm tests use Hiruma 3-term Lanczos. This script
takes the same K, M, b from radia-axifemm Q2 disk assembly and applies
the Kameari Foster-to-Kameari Hankel-Pade extraction:
  alpha_n = (-1)^n sum_k g_k^2 tau_k^n
where {tau_k, g_k} is the Foster spectrum from K v = lambda M v.
"""
from __future__ import annotations

import sys
sys.stdout.reconfigure(encoding="utf-8")

import json
import os
from math import pi
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.linalg as la
import mpmath as mp

mp.mp.dps = 80

# Add radia-axifemm tests dir to path so we can reuse the mesh builder
sys.path.insert(0, str(Path("S:/Radia/01_GitHub/packages/radia-axifemm/tests")))

from ngsolve import (
    Mesh, BilinearForm, LinearForm, CoefficientFunction, TaskManager,
    x, dx, ngsglobals,
)
from radia.radia_axifemm import (
    H1Henrotte, AxiHenrotteStiffnessBFI, AxiHenrotteSigmaMassBFI,
)
from test_hiruma_disk_q1 import (  # type: ignore
    make_structured_disk_quad_mesh,
    to_scipy_csr,
    R_DISK, T_DISK, SIGMA_CU, MU0, B0,
)


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
    # Mesh: same as Q2 fine
    NR_disk, Nz_disk, NR_air, Nz_air = 30, 6, 60, 30
    R_air, Z_air = 50e-3, 50e-3

    print(f"axifemm Q2 disk + Kameari extraction")
    print(f"  NR_disk={NR_disk}, Nz_disk={Nz_disk}, NR_air={NR_air}, Nz_air={Nz_air}")
    print(f"  R={R_DISK*1000} mm, t={T_DISK*1000} mm, sigma={SIGMA_CU:.2e}")
    print()

    ngsglobals.msg_level = 0
    mesh = make_structured_disk_quad_mesh(NR_disk, Nz_disk, NR_air, Nz_air,
                                           R_air, Z_air)
    fes = H1Henrotte(mesh, order=2, dirichlet="axis|right|top|bot")
    n_free = sum(1 for f in fes.FreeDofs() if f)
    print(f"  mesh ne={mesh.ne}, fes ndof={fes.ndof}, free={n_free}")

    mu_cf = CoefficientFunction(MU0)
    sigma_cf = mesh.MaterialCF({"conductor": SIGMA_CU}, default=0.0)
    A_imposed = B0 * x / 2  # A_phi imposed for uniform B_z

    a = BilinearForm(fes, symmetric=True)
    a += AxiHenrotteStiffnessBFI(mu_cf)
    with TaskManager():
        a.Assemble()

    m = BilinearForm(fes, symmetric=True)
    m += AxiHenrotteSigmaMassBFI(sigma_cf)
    with TaskManager():
        m.Assemble()

    from ngsolve import LinearForm
    b_form = LinearForm(fes)
    b_form += A_imposed * fes.TestFunction() * sigma_cf * dx
    with TaskManager():
        b_form.Assemble()

    K_csr = to_scipy_csr(a.mat, fes.ndof)
    M_csr = to_scipy_csr(m.mat, fes.ndof)
    b_full = np.array(b_form.vec)

    free = np.array([i for i in range(fes.ndof) if fes.FreeDofs()[i]],
                     dtype=int)
    K = K_csr[free[:, None], free[None, :]].toarray()
    M = M_csr[free[:, None], free[None, :]].toarray()
    b = b_full[free]

    print(f"  K, M, b assembled. K shape = {K.shape}")

    # M is sigma-weighted mass, only nonzero on conductor DoFs
    # Identify conductor DoFs
    M_diag = np.diag(M)
    cond_idx = np.where(np.abs(M_diag) > 1e-30)[0]
    print(f"  conductor DoFs = {len(cond_idx)} of {len(free)}")

    # Kameari moments DIRECT from K, M, b (Hiruma matrix form):
    #   alpha_n_Kameari = u^T M (K^{-1} M)^n u   with u = M^{-1} b
    # Use full K (not just cond) to handle Kelvin / air properly
    # For axifemm without Kelvin, use just cond block (M is zero outside)
    # Solve K_cond w = M_cond u_n  recursively
    K_cond = K[cond_idx[:, None], cond_idx[None, :]]
    M_cond = M[cond_idx[:, None], cond_idx[None, :]]
    b_cond = b[cond_idx]

    # Compute alpha_n via Krylov recurrence (no eigenvalue decomposition needed)
    # u_0 = M^{-1} b, T = K^{-1} M, T^n u_0 = ?
    # alpha_n = u_0^T M T^n u_0
    print(f"  Computing Kameari moments via Krylov on K, M (cond block)...")
    K_lu = la.lu_factor(K_cond)
    M_lu = la.lu_factor(M_cond)
    u = la.lu_solve(M_lu, b_cond)  # u_0 = M^{-1} b
    n_moments = 30
    alphas_float = []
    Tu = u.copy()  # T^0 u = u
    for n in range(n_moments):
        # alpha_n = u^T M (Tu) where Tu = T^n u
        a_val = float(u @ (M_cond @ Tu))
        alphas_float.append(a_val)
        # advance Tu -> T (Tu) = K^{-1} M (Tu)
        MTu = M_cond @ Tu
        Tu = la.lu_solve(K_lu, MTu)

    print(f"  alpha_0 = {alphas_float[0]:.4e}")
    print(f"  alpha_1 = {alphas_float[1]:.4e}")
    print(f"  alpha_2 = {alphas_float[2]:.4e}")
    R0_canonical = 1.0 / alphas_float[0]
    print(f"  R_0_canonical (Kameari) = 1/alpha_0 = {R0_canonical:.4e} (expected 2195 ohm)")

    # alpha_n with sign (-1)^n for Hankel-Padé
    alphas = [(mp.mpf(-1))**n * mp.mpf(alphas_float[n]) for n in range(n_moments)]

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
        "method": "axifemm Q2 axisym FEM + Kameari Hankel-Pade Cauer-I",
        "params": {
            "NR_disk": NR_disk, "Nz_disk": Nz_disk,
            "NR_air": NR_air, "Nz_air": Nz_air,
            "R_air_mm": R_air * 1000, "Z_air_mm": Z_air * 1000,
            "R_disk_mm": R_DISK * 1000, "t_disk_mm": T_DISK * 1000,
            "sigma": SIGMA_CU, "mu0": MU0, "B0": B0,
        },
        "leading_tau_us": float(tau_us[order[0]]),
        "Cauer_rungs_Kameari": rungs,
    }
    out_path = (Path(__file__).parent
                / "axifemm_disk_kameari_results.json")
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\n  Saved: {out_path}")


if __name__ == "__main__":
    main()
