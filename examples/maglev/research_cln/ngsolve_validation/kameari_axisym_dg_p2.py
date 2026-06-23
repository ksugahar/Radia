"""kameari_axisym_dg_p2.py — Kameari accumulation in P2 DG (FreeFEM P2DC analog).

Tests user's hypothesis: cond/air interface continuity in H1 P2 causes 5%
Kameari τ_0 gap. P2 Discontinuous (DG) lifts continuity → may eliminate gap.

Implementation:
  - L2(mesh, order=2, dgjumps=True) for u (discontinuous per element)
  - SIPG bilinear form with interior penalty
  - Same Kameari accumulation algorithm as cln_team28_axisym.py
  - cylinder R=10mm t=2mm Cu, B_0=1T axisym, air-box
"""
from __future__ import annotations

import json
import sys
import time
from math import pi
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)
sys.path.insert(0, r"S:/Radia/01_GitHub/examples/axifem/research/tests")

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from ngsolve import (
    Mesh, L2, BilinearForm, LinearForm, CoefficientFunction, GridFunction,
    TaskManager, x, y, dx, IfPos, grad, specialcf,
    Integrate, ngsglobals,
)

from test_hiruma_disk_q1 import (
    make_structured_disk_quad_mesh,
    R_DISK, T_DISK, SIGMA_CU, MU0, B0,
)


def to_csr(bf_mat, n):
    rs, cs, vs = bf_mat.COO()
    K = sp.csr_matrix(
        (np.asarray(vs), (np.asarray(rs), np.asarray(cs))),
        shape=(n, n))
    return (K + K.T) * 0.5


def run_kameari_dg(NR_disk, Nz_disk, NR_air, Nz_air, R_air, Z_air,
                    N_stages=4, order=2, eta=10.0, label=""):
    print(f"\n=== Kameari DG P{order}: {label} ===", flush=True)
    print(f"  NR_d={NR_disk} Nz_d={Nz_disk} NR_a={NR_air} Nz_a={Nz_air}",
          flush=True)
    ngsglobals.msg_level = 0
    mesh = make_structured_disk_quad_mesh(NR_disk, Nz_disk, NR_air, Nz_air,
                                          R_air, Z_air)

    # L2 = P2 discontinuous Galerkin
    fes = L2(mesh, order=order, dgjumps=True)
    print(f"  ne={mesh.ne} ndof={fes.ndof}", flush=True)

    nu_cf = CoefficientFunction(1.0 / MU0)
    sigma_cf = mesh.MaterialCF({"conductor": SIGMA_CU}, default=0.0)
    A_imposed = B0 * x / 2  # = A_θ in NGSolve coords (x = r)

    u, v = fes.TnT()
    n_vec = specialcf.normal(2)
    h = specialcf.mesh_size
    r_weight = IfPos(x - 1e-10, x, 1e-10)

    # u_imposed = r * A_θ_imposed = B_0 r²/2 (continuous representation)
    u_imposed_cf = B0 * x**2 / 2

    # SIPG bilinear form for K: -∇·((nu/r) ∇u) = source
    # K_DG(u, v) = ∫(nu/r) grad(u)·grad(v) dx
    #             - ∫_int {(nu/r)∂_n u} [v] dS
    #             - ∫_int [u] {(nu/r)∂_n v} dS
    #             + (eta/h) ∫_int [u][v] dS
    # where {·} = average across interface, [·] = jump

    # Compute jump and average
    jump_u = u - u.Other()
    jump_v = v - v.Other()
    mean_grad_u = 0.5 * (grad(u) + grad(u.Other()))
    mean_grad_v = 0.5 * (grad(v) + grad(v.Other()))

    a = BilinearForm(fes, symmetric=True)
    a += nu_cf / r_weight * grad(u) * grad(v) * dx
    # Interior penalty (SIPG)
    a += -nu_cf / r_weight * (mean_grad_u * n_vec) * jump_v * dx(skeleton=True)
    a += -nu_cf / r_weight * (mean_grad_v * n_vec) * jump_u * dx(skeleton=True)
    a += eta * (order ** 2) / h * nu_cf / r_weight * jump_u * jump_v * dx(
        skeleton=True)
    # Dirichlet on outer boundary (axis|right|top|bot)
    # For DG: weak Dirichlet via boundary integrals
    a += -nu_cf / r_weight * (grad(u) * n_vec) * v * dx(
        skeleton=True, definedon=mesh.Boundaries("axis|right|top|bot"))
    a += -nu_cf / r_weight * (grad(v) * n_vec) * u * dx(
        skeleton=True, definedon=mesh.Boundaries("axis|right|top|bot"))
    a += eta * (order ** 2) / h * nu_cf / r_weight * u * v * dx(
        skeleton=True, definedon=mesh.Boundaries("axis|right|top|bot"))

    print(" Assembling K (SIPG) ...", flush=True)
    t0 = time.time()
    with TaskManager(): a.Assemble()
    print(f"  K assemble: {time.time()-t0:.1f} s", flush=True)

    # M = ∫(σ/r) u v dx (volume only, no skeleton)
    m = BilinearForm(fes, symmetric=True)
    m += sigma_cf / r_weight * u * v * dx
    print(" Assembling M ...", flush=True)
    with TaskManager(): m.Assemble()

    # b = ∫(σ × u_imposed / r) × v dx  (= M × u_imposed_dof analog)
    b_form = LinearForm(fes)
    b_form += sigma_cf * u_imposed_cf / r_weight * v * dx
    with TaskManager(): b_form.Assemble()
    b_full = np.array(b_form.vec)

    K_csr = to_csr(a.mat, fes.ndof)
    M_csr = to_csr(m.mat, fes.ndof)

    # No Dirichlet DOFs in DG (handled weakly via boundary skeleton)
    free = np.arange(fes.ndof)
    K_red = K_csr
    M_red = M_csr
    b_red = b_full

    M_diag = M_red.diagonal()
    cond_local = np.where(np.abs(M_diag) > 1e-30)[0]
    print(f"  cond DOFs: {len(cond_local)} / {len(free)}", flush=True)
    M_cond = M_red[cond_local[:, None], cond_local[None, :]]
    b_cond = b_red[cond_local]

    K_factor = spla.factorized(K_red.tocsc())
    M_cond_factor = spla.factorized(M_cond.tocsc())

    u_dof_cond = M_cond_factor(b_cond)
    A_acc_full = np.zeros(len(free))

    Rn, Ln, tau_us = [], [], []
    for nstage in range(N_stages):
        G_n = float(u_dof_cond @ (M_cond @ u_dof_cond))
        if G_n < 1e-50 or not np.isfinite(G_n):
            print(f"  Stage {nstage}: G bad", flush=True)
            break
        R_n = 1.0 / G_n
        Rn.append(R_n)

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
             "ne": mesh.ne, "ndof": fes.ndof,
             "n_cond_dof": len(cond_local)}


def main():
    cases = [
        (20, 8, 12, 12, 200e-3, 200e-3, "coarse-DG"),
        (40, 16, 15, 15, 200e-3, 200e-3, "medium-DG"),
    ]
    all_res = {}
    for NR_d, Nz_d, NR_a, Nz_a, R_a, Z_a, lbl in cases:
        all_res[lbl] = run_kameari_dg(NR_d, Nz_d, NR_a, Nz_a, R_a, Z_a,
                                        N_stages=4, order=2, eta=10.0,
                                        label=lbl)

    print("\n" + "=" * 78, flush=True)
    print(" Kameari DG P2 (= FreeFEM P2DC analog) results", flush=True)
    print("=" * 78, flush=True)
    print(f"{'case':>15} {'ne':>6} {'ndof':>6} | "
          f"{'τ_0':>10} {'τ_1':>10} {'τ_2':>10} {'τ_3':>10}",
          flush=True)
    for lbl in [c[6] for c in cases]:
        r = all_res[lbl]
        ts = r["tau_n_us"]
        if len(ts) >= 4:
            print(f"{lbl:>15} {r['ne']:>6} {r['ndof']:>6} | "
                  f"{ts[0]:>10.4f} {ts[1]:>10.4f} {ts[2]:>10.4f}"
                  f" {ts[3]:>10.4f}",
                  flush=True)
        else:
            print(f"{lbl:>15} {r['ne']:>6} {r['ndof']:>6} | "
                  f"{ts}",
                  flush=True)
    print(f"{'v23 Hiruma ref':>15} {'-':>6} {'-':>6} | "
          f"{'218.7078':>10} {'78.1340':>10} {'39.5591':>10} {'23.1679':>10}",
          flush=True)
    print(f"{'BEM-Foster':>15} {'-':>6} {'-':>6} | "
          f"{'219.32':>10} {'78.65':>10} {'40.04':>10} {'23.74':>10}",
          flush=True)
    print(f"{'Kameari H1 P2':>15} {'-':>6} {'-':>6} | "
          f"{'208.45':>10} {'67.54':>10} {'32.58':>10} {'20.41':>10}",
          flush=True)

    out_path = (Path(__file__).parent
                 / "2026-05-10-kameari_dg_p2_results.json")
    out_path.write_text(json.dumps(all_res, indent=2), encoding="utf-8")
    print(f"\nResults: {out_path}", flush=True)


if __name__ == "__main__":
    main()
