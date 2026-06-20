"""kameari_mesh_sweep.py — Mesh refinement sweep for Kameari + air-box axisym.

Tests whether Kameari τ_0 ≈ 208 μs (cylinder R=10mm t=2mm Cu, B_0=1T) is
discretization-limited or algorithm-structural.

If algorithm-limited (Padé[1,1] M-average): τ_0 invariant under mesh refinement.
If discretization-limited: τ_0 → 218.7 (Hiruma exact) as h→0.

Reuses Kameari accumulation algorithm in axifemm Henrotte basis (= same
framework as v23 Hiruma 3-term, ensuring apples-to-apples comparison).
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


def run_kameari(NR_disk, Nz_disk, NR_air, Nz_air, R_air, Z_air,
                 N_stages=4, order=2, label=""):
    print(f"\n=== {label}: NR_d={NR_disk} Nz_d={Nz_disk} "
          f"NR_a={NR_air} Nz_a={Nz_air} ===", flush=True)
    ngsglobals.msg_level = 0
    mesh = make_structured_disk_quad_mesh(NR_disk, Nz_disk, NR_air, Nz_air,
                                          R_air, Z_air)
    fes = H1Henrotte(mesh, order=order, dirichlet="axis|right|top|bot")
    n_free = sum(1 for f in fes.FreeDofs() if f)
    print(f"  ne={mesh.ne} ndof={fes.ndof} free={n_free}", flush=True)

    mu_cf = CoefficientFunction(MU0)
    sigma_cf = mesh.MaterialCF({"conductor": SIGMA_CU}, default=0.0)
    A_imposed = B0 * x / 2
    v = fes.TestFunction()

    a = BilinearForm(fes, symmetric=True)
    a += AxiHenrotteStiffnessBFI(mu_cf)
    with TaskManager(): a.Assemble()
    K_csr = to_scipy_csr(a.mat, fes.ndof)

    m = BilinearForm(fes, symmetric=True)
    m += AxiHenrotteSigmaMassBFI(sigma_cf)
    with TaskManager(): m.Assemble()
    M_csr = to_scipy_csr(m.mat, fes.ndof)

    b_form = LinearForm(fes)
    b_form += sigma_cf * A_imposed * v * 2 * pi * x * dx
    with TaskManager(): b_form.Assemble()
    b_full = np.array(b_form.vec)

    free = np.array([i for i in range(fes.ndof) if fes.FreeDofs()[i]],
                     dtype=int)
    K_red = K_csr[free[:, None], free[None, :]]
    M_red = M_csr[free[:, None], free[None, :]]
    b_red = b_full[free]

    M_diag = M_red.diagonal()
    cond_local = np.where(np.abs(M_diag) > 1e-30)[0]
    M_cond = M_red[cond_local[:, None], cond_local[None, :]]
    b_cond = b_red[cond_local]

    K_factor = spla.factorized(K_red.tocsc())
    M_cond_factor = spla.factorized(M_cond.tocsc())

    # Initial u_dof_cond from σ-weighted projection of A_imposed
    u_dof_cond = M_cond_factor(b_cond)
    A_acc_full = np.zeros(len(free))

    Rn, Ln, tau_us = [], [], []
    for n in range(N_stages):
        # G_n = u^T M_cond u
        G_n = float(u_dof_cond @ (M_cond @ u_dof_cond))
        if G_n < 1e-50 or not np.isfinite(G_n):
            print(f"  Stage {n}: G_n bad, breaking", flush=True)
            break
        R_n = 1.0 / G_n
        Rn.append(R_n)

        # Source for K
        b_n_full = np.zeros(len(free))
        b_n_full[cond_local] = M_cond @ u_dof_cond
        w_n = K_factor(b_n_full)

        A_acc_full += R_n * w_n
        A_acc_cond = A_acc_full[cond_local]
        tau_n = float((M_cond @ u_dof_cond) @ A_acc_cond)
        L_n = R_n * tau_n
        Ln.append(L_n)
        tau_us.append(tau_n * 1e6)

        u_dof_cond = u_dof_cond - A_acc_cond / L_n

    return {"Rn": Rn, "Ln": Ln, "tau_n_us": tau_us,
             "ne": mesh.ne, "ndof": fes.ndof}


def main():
    cases = [
        (10, 4,  8,  8, 100e-3, 100e-3, "very-coarse"),
        (20, 8, 12, 12, 200e-3, 200e-3, "coarse-Q2"),
        (40, 16, 15, 15, 200e-3, 200e-3, "medium-Q2"),
        (80, 32, 20, 20, 500e-3, 500e-3, "fine-Q2"),
    ]
    all_res = {}
    for NR_d, Nz_d, NR_a, Nz_a, R_a, Z_a, lbl in cases:
        all_res[lbl] = run_kameari(NR_d, Nz_d, NR_a, Nz_a, R_a, Z_a,
                                     N_stages=4, order=2, label=lbl)

    print("\n" + "=" * 78, flush=True)
    print(" Mesh refinement sweep — Kameari + axifemm Henrotte axisym",
          flush=True)
    print("=" * 78, flush=True)
    print(f"{'case':>15} {'ne':>6} {'ndof':>6} | "
          f"{'τ_0':>10} {'τ_1':>10} {'τ_2':>10} {'τ_3':>10}",
          flush=True)
    for lbl in [c[6] for c in cases]:
        r = all_res[lbl]
        ts = r["tau_n_us"]
        print(f"{lbl:>15} {r['ne']:>6} {r['ndof']:>6} | "
              f"{ts[0]:>10.4f} {ts[1]:>10.4f} {ts[2]:>10.4f}"
              f" {ts[3]:>10.4f}",
              flush=True)
    print(f"{'v23 Hiruma ref':>15} {'-':>6} {'-':>6} | "
          f"{'218.7078':>10} {'78.1340':>10} {'39.5591':>10} {'23.1679':>10}",
          flush=True)
    print(f"{'BEM-Foster':>15} {'-':>6} {'-':>6} | "
          f"{'219.32':>10} {'78.65':>10} {'40.04':>10} {'23.74':>10}",
          flush=True)

    out_path = (Path(__file__).parent
                 / "2026-05-10-kameari_mesh_sweep_results.json")
    out_path.write_text(json.dumps(all_res, indent=2), encoding="utf-8")
    print(f"\nResults: {out_path}", flush=True)


if __name__ == "__main__":
    main()
