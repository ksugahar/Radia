"""cln_team28_axisym.py — COMSOL TEAM 28 CLN algorithm in axisymmetric NGSolve.

Implementation entirely in DOF space using axifemm K, M_σ matrices.

Variables (DOF space, restricted to conductor DOFs where M_σ is invertible):
  A_dof_cond: continuous A_phi DOF representation on cond DOFs
  A_acc_cond: accumulated CLN vector potential (cond DOFs)
  A_acc_full: accumulated CLN vector potential (full mesh, for K solves)

Iteration:
  G_n = A_dof_cond^T M_σ_cond A_dof_cond                 (= ∫σA² dV on cond)
  R_n = 1/G_n
  b_n_full = M_σ × A_full (where A_full has A_dof_cond on cond, 0 on air)
  Solve K_full × w_n_full = b_n_full → A_phi response
  A_acc_full += R_n × w_n_full
  τ_n = b_n_cond · A_acc_full[cond_dofs]   (= ∫J·A_acc dV via Galerkin)
  L_n = R_n × τ_n
  A_dof_cond -= A_acc_full[cond_dofs] / L_n   (next iterate)
"""
from __future__ import annotations

import json
import sys
import time
from math import pi
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)
sys.path.insert(0, r"S:/Radia/01_GitHub/examples/axifemm/research/tests")

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from ngsolve import (
    Mesh, BilinearForm, LinearForm, CoefficientFunction, GridFunction,
    TaskManager, x, dx, ngsglobals,
)
from axifem import (
    H1Henrotte, AxiHenrotteStiffnessBFI, AxiHenrotteSigmaMassBFI,
)

from test_hiruma_disk_q1 import (
    make_structured_disk_quad_mesh, to_scipy_csr,
    R_DISK, T_DISK, SIGMA_CU, MU0, B0,
)


def cln_team28_axisym(NR_disk=20, Nz_disk=8, NR_air=12, Nz_air=12,
                       R_air=200e-3, Z_air=200e-3, N_stages=6, order=2,
                       label=""):
    print(f"\n=== TEAM 28 axisym Q{order}: NR_d={NR_disk} Nz_d={Nz_disk} "
          f"NR_a={NR_air} Nz_a={Nz_air} R_air={R_air*1e3}mm  {label} ===",
          flush=True)
    ngsglobals.msg_level = 0
    mesh = make_structured_disk_quad_mesh(NR_disk, Nz_disk, NR_air, Nz_air,
                                          R_air, Z_air)
    fes = H1Henrotte(mesh, order=order, dirichlet="axis|right|top|bot")
    n_free = sum(1 for f in fes.FreeDofs() if f)
    print(f"  mesh ne={mesh.ne} ndof={fes.ndof}  free={n_free}",
          flush=True)

    mu_cf = CoefficientFunction(MU0)
    sigma_cf = mesh.MaterialCF({"conductor": SIGMA_CU}, default=0.0)
    A_imposed = B0 * x / 2
    v = fes.TestFunction()

    # K = AxiHenrotte stiffness with 1/μ
    a = BilinearForm(fes, symmetric=True)
    a += AxiHenrotteStiffnessBFI(mu_cf)
    with TaskManager(): a.Assemble()
    K_csr = to_scipy_csr(a.mat, fes.ndof)

    # M_σ = AxiHenrotte mass with σ (conductor only)
    m = BilinearForm(fes, symmetric=True)
    m += AxiHenrotteSigmaMassBFI(sigma_cf)
    with TaskManager(): m.Assemble()
    M_csr = to_scipy_csr(m.mat, fes.ndof)

    # Initial loading: b = ∫σ A_imposed v · 2πx dx = M_σ × A_imposed_dof
    b_form = LinearForm(fes)
    b_form += sigma_cf * A_imposed * v * 2 * pi * x * dx
    with TaskManager(): b_form.Assemble()
    b_full = np.array(b_form.vec)

    # Free DOFs (Dirichlet excluded)
    free = np.array([i for i in range(fes.ndof) if fes.FreeDofs()[i]],
                     dtype=int)
    K_red = K_csr[free[:, None], free[None, :]]
    M_red = M_csr[free[:, None], free[None, :]]
    b_red = b_full[free]

    # Conductor DOFs (M diag != 0)
    M_diag = M_red.diagonal()
    cond_mask = np.abs(M_diag) > 1e-30
    cond_local_idx = np.where(cond_mask)[0]   # indices into free-DOF space
    n_cond = len(cond_local_idx)
    print(f"  free={len(free)}, cond DOFs={n_cond}", flush=True)

    M_cond = M_red[cond_local_idx[:, None], cond_local_idx[None, :]]
    b_cond = b_red[cond_local_idx]

    K_factor = spla.factorized(K_red.tocsc())
    M_cond_factor = spla.factorized(M_cond.tocsc())

    # Initial A_dof_cond: σ-weighted projection of A_imposed onto cond DOFs
    # M_cond × A_dof_cond = b_cond  →  A_dof_cond = M_cond^{-1} × b_cond
    A_dof_cond = M_cond_factor(b_cond)

    A_acc_full = np.zeros(len(free))     # full free-DOF space accumulator

    Rn, Ln, tau_us = [], [], []
    diag = []

    for n in range(N_stages):
        t_stage = time.time()

        # G_n = A_dof_cond^T M_cond A_dof_cond = ∫σA² dV
        G_n = float(A_dof_cond @ (M_cond @ A_dof_cond))
        if G_n < 1e-50 or not np.isfinite(G_n):
            print(f"  Stage {n}: G_n = {G_n}, breaking", flush=True)
            break
        R_n = 1.0 / G_n
        Rn.append(R_n)

        # Loading for K equation: b_n_full has b_n_cond = M_cond × A_dof_cond,
        # 0 elsewhere
        b_n_cond = M_cond @ A_dof_cond
        b_n_full = np.zeros(len(free))
        b_n_full[cond_local_idx] = b_n_cond

        # Solve K_full × w_n_full = b_n_full
        w_n_full = K_factor(b_n_full)

        # A_acc_full += R_n × w_n_full
        A_acc_full += R_n * w_n_full

        # τ_n = b_n_cond · A_acc_full[cond_local_idx]  (Galerkin: ∫J·A_acc dV)
        A_acc_cond = A_acc_full[cond_local_idx]
        tau_n = float(b_n_cond @ A_acc_cond)
        tau_n_us_val = tau_n * 1e6
        L_n = R_n * tau_n
        Ln.append(L_n)
        tau_us.append(tau_n_us_val)

        # Update A_dof_cond: A_{n+1}_cond = A_n_cond - A_acc_cond / L_n
        A_dof_cond = A_dof_cond - A_acc_cond / L_n

        diag.append({
            "stage": n,
            "G_n": G_n,
            "R_n": R_n,
            "tau_n_us": tau_n_us_val,
            "L_n": L_n,
            "stage_time_s": time.time() - t_stage,
        })
        print(f"  Stage {n}: G={G_n:.4e}, R={R_n:.4e},"
              f" τ={tau_n_us_val:.4f} μs, L={L_n:.4e}"
              f"  ({time.time()-t_stage:.2f} s)",
              flush=True)

    # Reference comparison
    print("\n" + "=" * 78, flush=True)
    print(" Summary (TEAM 28 axisym, COMSOL formula)", flush=True)
    print("=" * 78, flush=True)
    print("\n  Reference (v23 axisym Hiruma 3-term, ORDER=3 fine):",
          flush=True)
    v23_ref = [218.6222, 78.2068, 39.5928, 23.1196, 16.0682, 13.3002]
    for k in range(min(N_stages, len(v23_ref))):
        match = ""
        if k < len(tau_us):
            gap = (tau_us[k] - v23_ref[k]) / v23_ref[k] * 100
            match = (f"  COMSOL τ_{k} = {tau_us[k]:.4f} μs"
                     f"  gap = {gap:+.1f}%")
        print(f"    k={k}: v23 τ = {v23_ref[k]:.4f} us{match}",
              flush=True)

    return {"Rn": Rn, "Ln": Ln, "tau_n_us": tau_us, "stages": diag}


def main():
    # R_air sweep: test if larger air-box approaches v23 = 218.62
    print("\n" + "#" * 78)
    print("# R_air sweep — does enlarging air-box approach v23 218.62?")
    print("#" * 78)

    sweep_cases = [
        (200e-3, "R_air=200mm"),
        (500e-3, "R_air=500mm (= v23 fine)"),
        (1000e-3, "R_air=1000mm"),
        (2000e-3, "R_air=2000mm"),
    ]
    sweep_results = {}
    for R_a, lbl in sweep_cases:
        res = cln_team28_axisym(
            NR_disk=20, Nz_disk=8, NR_air=15, Nz_air=15,
            R_air=R_a, Z_air=R_a, N_stages=6, order=2,
            label=lbl)
        sweep_results[f"R_air_{R_a*1000:.0f}mm"] = res

    # Summary
    print("\n" + "#" * 78)
    print("# SWEEP SUMMARY")
    print("#" * 78)
    print(f"{'R_air':>15} | {'τ_0':>10} {'τ_1':>10} {'τ_2':>10}"
          f" {'τ_3':>10} {'τ_4':>10} {'τ_5':>10}")
    for R_a, lbl in sweep_cases:
        key = f"R_air_{R_a*1000:.0f}mm"
        taus = sweep_results[key]["tau_n_us"]
        row = f"{R_a*1000:>10.0f} mm  |"
        for t in taus[:6]:
            row += f" {t:>10.4f}"
        print(row, flush=True)
    print(f"{'v23 ref (Hiruma)':>15} |  218.6222"
          f"   78.2068   39.5928   23.1196   16.0682   13.3002")

    out_path = (Path(__file__).parent
                 / f"2026-05-10-cln_team28_axisym_Rair_sweep.json")
    out_path.write_text(json.dumps(sweep_results, indent=2), encoding="utf-8")
    print(f"\nResults: {out_path}", flush=True)


if __name__ == "__main__":
    main()
