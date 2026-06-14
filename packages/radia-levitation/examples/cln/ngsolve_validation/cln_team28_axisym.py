"""Historical TEAM 28 axisymmetric assembly setup.

The field-reduction recurrence and its dependent reporting are retired.
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
    Mesh, BilinearForm, LinearForm, CoefficientFunction, GridFunction,
    TaskManager, x, dx, ngsglobals,
)
from radia_axifemm import (
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

    raise RuntimeError("Historical field-reduction recurrence retired; no stage coefficients are produced.")


def main():
    raise RuntimeError("Historical field-reduction driver retired; assembly setup is retained for reference.")


if __name__ == "__main__":
    main()
