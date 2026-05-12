"""Axisymmetric SIPG operator assembly on the original cylinder mesh.

The historical field-recurrence experiment has been retired.
This driver reports assembly dimensions only, not a reduced-model validation.
"""
from __future__ import annotations

import json
import sys
import time
from math import pi
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)
sys.path.insert(0, r"S:/Radia/01_GitHub/packages/radia-axifemm/tests")

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


def assemble_axisymmetric_sipg(NR_disk, Nz_disk, NR_air, Nz_air, R_air, Z_air,
                    order=2, eta=10.0, label=""):
    print(f"\n=== Axisymmetric SIPG P{order}: {label} ===", flush=True)
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

    return {"K": K_red, "M": M_red, "b": b_red, "mesh": mesh, "fes": fes,
            "conductor_dofs": cond_local}


def main():
    cases = [
        (20, 8, 12, 12, 200e-3, 200e-3, "coarse-DG"),
        (40, 16, 15, 15, 200e-3, 200e-3, "medium-DG"),
    ]
    all_res = {}
    for NR_d, Nz_d, NR_a, Nz_a, R_a, Z_a, lbl in cases:
        operators = assemble_axisymmetric_sipg(
            NR_d, Nz_d, NR_a, Nz_a, R_a, Z_a,
            order=2, eta=10.0, label=lbl)
        all_res[lbl] = {
            "ne": operators["mesh"].ne,
            "ndof": operators["fes"].ndof,
            "n_cond_dof": len(operators["conductor_dofs"]),
            "K_nnz": operators["K"].nnz,
            "M_nnz": operators["M"].nnz,
            "status": "assembled; no reduced model extracted",
        }
        print(lbl, all_res[lbl], flush=True)
    out_path = Path(__file__).parent / "axisymmetric_sipg_assembly.json"
    out_path.write_text(json.dumps(all_res, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
