"""axifemm_disk_kameari_v2.py — axifemm cylinder Cauer-I in Kameari convention.

Direct Foster -> Kameari extraction from axifemm K, M, b using Schur
complement onto conductor DOFs. The earlier v1 used K_cond directly which
ignores the air-side coupling (Schur complement) and gave wrong tau ~ 24 us.

Algorithm:
  Partition free DOFs into (cond, air). Build:
    S    = K_cc - K_ca @ K_aa^{-1} @ K_ac     (Schur complement on cond)
    M_cc = M restricted to cond block
    b_c  = b restricted to cond
    u_c  = M_cc^{-1} b_c  (= M^{-1} b restricted to cond, since M_aa = 0)
  Then:
    alpha_n = u_c^T M_cc (S^{-1} M_cc)^n u_c
  For Foster spectrum: solve generalized eigenproblem
    S v = lambda M_cc v   (lambda = 1/tau, in 1/s)
  M_cc-normalize v (v^T M_cc v = 1), so g_k = v_k^T b_c are Foster amplitudes
  matching the Kameari accumulation moments.
"""
from __future__ import annotations

import json
import os
import sys
from math import pi
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import scipy.linalg as la
import mpmath as mp

mp.mp.dps = 80

sys.path.insert(0,
    str(Path("S:/Radia/01_GitHub/examples/axifemm/research/tests")))

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


def schur_foster_kameari(K, M, b, cond_idx, air_idx, label=""):
    """Build Schur complement on cond block, solve Foster eigenproblem,
    compute Kameari moments and Cauer-I rungs.

    K, M, b: full FREE matrices/vector (sparse K, M; numpy b).
    cond_idx, air_idx: indices into FREE DOF range.
    Returns dict with rungs, leading tau, etc.
    """
    print(f"\n  [{label}] cond DOFs={len(cond_idx)}, air DOFs={len(air_idx)}")

    K_dense = K.toarray()
    M_dense = M.toarray()
    K_aa = K_dense[air_idx[:, None], air_idx[None, :]]
    K_ac = K_dense[air_idx[:, None], cond_idx[None, :]]
    K_ca = K_dense[cond_idx[:, None], air_idx[None, :]]
    K_cc = K_dense[cond_idx[:, None], cond_idx[None, :]]
    M_cc = M_dense[cond_idx[:, None], cond_idx[None, :]]
    b_c = b[cond_idx]

    print(f"  Building Schur complement S = K_cc - K_ca K_aa^{-1} K_ac ...")
    K_aa_inv_K_ac = la.solve(K_aa, K_ac, assume_a="pos")
    S = K_cc - K_ca @ K_aa_inv_K_ac
    S = 0.5 * (S + S.T)
    print(f"  S shape = {S.shape}, ||S||_F = {np.linalg.norm(S):.3e}")
    print(f"  M_cc shape = {M_cc.shape}, "
          f"M_cc PSD eig range: [{la.eigh(M_cc, eigvals_only=True)[0]:.3e}, "
          f"{la.eigh(M_cc, eigvals_only=True)[-1]:.3e}]")

    # Generalized eigenproblem S v = lambda M_cc v
    eigvals, eigvecs = la.eigh(S, M_cc)
    # M_cc-normalize
    for k in range(len(eigvals)):
        v = eigvecs[:, k]
        nrm = float(v @ (M_cc @ v))
        if nrm > 1e-30:
            eigvecs[:, k] = v / np.sqrt(nrm)

    # tau_k = 1 / lambda_k  (axifemm convention)
    tau_seconds_all = np.where(eigvals > 1e-30, 1.0 / eigvals, 0.0)
    tau_us = tau_seconds_all * 1e6
    g = eigvecs.T @ b_c
    g2 = g * g

    order = np.argsort(-g2)
    print(f"  Top 10 modes by |g|^2:")
    for rank in range(min(10, len(order))):
        idx = order[rank]
        print(f"    rank={rank}: tau={tau_us[idx]:.4f} us, "
              f"g2={g2[idx]:.4e}, lam={eigvals[idx]:.4e}")

    # Use top n_use modes for Kameari Foster moments
    n_use = min(60, int(np.sum(g2 > 1e-32)))
    sorted_idx = order[:n_use]
    tau_seconds = [mp.mpf(float(tau_seconds_all[i])) for i in sorted_idx]
    g2_mp = [mp.mpf(float(g2[i])) for i in sorted_idx]

    n_moments = min(40, 2 * n_use - 4)
    alphas = []
    for n in range(n_moments):
        a = mp.mpf(0)
        for gv, tv in zip(g2_mp, tau_seconds):
            a += gv * tv ** n
        alphas.append((mp.mpf(-1)) ** n * a)

    p = hankel_pade_cauer(alphas, 12)
    print(f"  Kameari Cauer-I rungs (mpmath 80 dig):")
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

    return {
        "label": label,
        "leading_tau_us": float(tau_us[order[0]]),
        "Cauer_rungs_Kameari": rungs,
        "n_modes_used": n_use,
        "top_modes": [
            {"rank": r, "tau_us": float(tau_us[order[r]]),
             "g2": float(g2[order[r]])}
            for r in range(min(10, len(order)))
        ],
    }


def solve_disk_kameari(NR_disk, Nz_disk, NR_air, Nz_air, R_air, Z_air,
                       label=""):
    print(f"\n=== axifemm Kameari Cu disk: NR_d={NR_disk} Nz_d={Nz_disk} "
          f"NR_a={NR_air} Nz_a={Nz_air} R_air={R_air * 1e3}mm  {label} ===")
    ngsglobals.msg_level = 0
    mesh = make_structured_disk_quad_mesh(NR_disk, Nz_disk, NR_air, Nz_air,
                                          R_air, Z_air)
    fes = H1Henrotte(mesh, dirichlet="axis|right|top|bot")
    n_free = sum(1 for f in fes.FreeDofs() if f)
    print(f"  mesh ne={mesh.ne}, ndof={fes.ndof}, free={n_free}")

    mu_cf = CoefficientFunction(MU0)
    sigma_cf = mesh.MaterialCF({"conductor": SIGMA_CU}, default=0.0)
    A_imposed = B0 * x / 2

    a = BilinearForm(fes, symmetric=True)
    a += AxiHenrotteStiffnessBFI(mu_cf)
    with TaskManager():
        a.Assemble()
    m = BilinearForm(fes, symmetric=True)
    m += AxiHenrotteSigmaMassBFI(sigma_cf)
    with TaskManager():
        m.Assemble()
    b_form = LinearForm(fes)
    v = fes.TestFunction()
    b_form += sigma_cf * A_imposed * v * 2 * pi * x * dx
    with TaskManager():
        b_form.Assemble()

    K_csr = to_scipy_csr(a.mat, fes.ndof)
    M_csr = to_scipy_csr(m.mat, fes.ndof)
    b_full = np.array(b_form.vec)

    free = np.array([i for i in range(fes.ndof) if fes.FreeDofs()[i]],
                    dtype=int)
    K_red = K_csr[free[:, None], free[None, :]]
    M_red = M_csr[free[:, None], free[None, :]]
    b_red = b_full[free]

    # Identify cond vs air based on M_red diagonal (M is sigma-mass,
    # zero on air DOFs)
    M_diag = np.array(M_red.diagonal()).flatten()
    cond_mask = np.abs(M_diag) > 1e-25
    cond_idx = np.where(cond_mask)[0]
    air_idx = np.where(~cond_mask)[0]

    return schur_foster_kameari(K_red, M_red, b_red, cond_idx, air_idx,
                                label=label)


def main():
    cases = [
        (40,  8, 15, 15, 200e-3, 200e-3, "medium"),
        (80, 16, 20, 20, 500e-3, 500e-3, "fine"),
    ]
    all_res = {}
    for NR_d, Nz_d, NR_a, Nz_a, R_a, Z_a, lbl in cases:
        all_res[lbl] = solve_disk_kameari(NR_d, Nz_d, NR_a, Nz_a, R_a, Z_a,
                                          label=lbl)

    print("\n" + "="*78)
    print("Summary: axifemm Q1 Kameari Cauer-I (cylinder Cu disk)")
    print("="*78)
    print(f"  BEM Kameari leading tau (reference) = 209 us")
    print(f"  axifemm Hiruma leading tau (reference) = 218.71 us (Q2 fine)")
    for lbl, r in all_res.items():
        rungs = r["Cauer_rungs_Kameari"]
        if rungs:
            tau0 = rungs[0]["tau_pair_us"]
            print(f"  {lbl:<10} axifemm Q1 Kameari tau_pair[0] = {tau0:.4f} us")

    out_path = Path(__file__).parent / "axifemm_disk_kameari_v2_results.json"
    out_path.write_text(json.dumps(all_res, indent=2), encoding="utf-8")
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
