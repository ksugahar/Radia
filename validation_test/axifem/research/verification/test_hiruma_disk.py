"""Axisymmetric copper-disk operator assembly and generalized eigenmodes.

Uses the original stiffness, conductivity mass and imposed-field load.
The scipy eigsh calculation and stored BEM-Foster reference comparison remain.
The historical Hiruma field-recurrence comparison has been retired."""

import json
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from math import pi

from ngsolve import (
    Mesh, BilinearForm, LinearForm, Integrate, GridFunction,
    CoefficientFunction, TaskManager, x, ngsglobals,
)
from netgen.occ import OCCGeometry, MoveTo, Glue, X, Y
from radia.axifem import (
    H1Henrotte, AxiHenrotteStiffnessBFI, AxiHenrotteSigmaMassBFI,
)


R_DISK = 10e-3
T_DISK = 2e-3
SIGMA_CU = 5.8e7
MU0 = 4 * pi * 1e-7
B0 = 1.0  # imposed axial flux density [T]

BEM_TAU_REF_PATH = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "..", "..", "maglev", "research_cln", "ngsolve_validation",
    "bem_disk_axisym_v3_refined.json",
))


def build_disk_mesh(R_air=200e-3, Z_air=200e-3, maxh_disk=0.2e-3, maxh_air=10e-3):
    air = MoveTo(0, -Z_air).Rectangle(R_air, 2 * Z_air).Face()
    disk = MoveTo(0, -T_DISK / 2).Rectangle(R_DISK, T_DISK).Face()
    disk.faces.name = "conductor"
    air.faces.name = "air"
    disk.maxh = maxh_disk
    air.edges.Min(X).name = "axis"
    air.edges.Max(X).name = "right"
    air.edges.Min(Y).name = "bot"
    air.edges.Max(Y).name = "top"
    return Mesh(OCCGeometry(Glue([air, disk]), dim=2).GenerateMesh(maxh=maxh_air))


def to_scipy_csr(mat, n):
    rs, cs, vs = mat.COO()
    K = sp.csr_matrix(
        (np.asarray(vs, dtype=np.float64),
         (np.asarray(rs, dtype=np.int64), np.asarray(cs, dtype=np.int64))),
        shape=(n, n))
    return (K + K.T) * 0.5






def solve_disk_modes(maxh_disk=0.2e-3, maxh_air=10e-3, R_air=200e-3, Z_air=200e-3,
                      N_stages=6, label=""):
    print(f"\n=== Eigenmodes on Cu disk: maxh_disk={maxh_disk*1e3}mm "
          f"R_air={R_air*1e3}mm  {label} ===")
    ngsglobals.msg_level = 0
    mesh = build_disk_mesh(R_air, Z_air, maxh_disk, maxh_air)
    fes = H1Henrotte(mesh, dirichlet="axis|right|top|bot")
    n_free = sum(1 for f in fes.FreeDofs() if f)
    print(f"  mesh ne={mesh.ne}  ndof={fes.ndof}  free={n_free}")

    mu_cf = CoefficientFunction(MU0)
    sigma_cf = mesh.MaterialCF({"conductor": SIGMA_CU}, default=0.0)
    A_imposed = B0 * x / 2  # uniform B_z = B0 -> A_phi = B0 r / 2

    a = BilinearForm(fes, symmetric=True)
    a += AxiHenrotteStiffnessBFI(mu_cf)
    with TaskManager(): a.Assemble()

    m = BilinearForm(fes, symmetric=True)
    m += AxiHenrotteSigmaMassBFI(sigma_cf)
    with TaskManager(): m.Assemble()

    # b_vec = sigma * A_phi_imposed * v assembled with the same DOF convention.
    # In V-DOF, the linear form has natural integrand sigma * A_imposed * N_A_v
    # but for consistency with M (which is sigma * A * A * 2 pi r dr dz), we use
    # the same 2 pi r weight. NGSolve LinearForm with sigma_cf * A_imposed * v * dx
    # gives Integrate(sigma A v dr dz) — missing 2 pi r weight.
    # We add it explicitly via 2 pi x:
    b_form = LinearForm(fes)
    v = fes.TestFunction()
    b_form += sigma_cf * A_imposed * v * 2 * pi * x * dx_marker(fes)
    with TaskManager(): b_form.Assemble()
    b_vec = np.array(b_form.vec)

    K_csr = to_scipy_csr(a.mat, fes.ndof)
    M_csr = to_scipy_csr(m.mat, fes.ndof)
    free = np.array([i for i in range(fes.ndof) if fes.FreeDofs()[i]], dtype=int)

    # Cross-check: scipy eigsh tau_n
    K_red = K_csr[free[:, None], free[None, :]]
    M_red = M_csr[free[:, None], free[None, :]]
    K_red = (K_red + K_red.T) * 0.5
    M_red = (M_red + M_red.T) * 0.5
    eigs, _ = spla.eigsh(K_red, k=min(N_stages, n_free // 2),
                         M=M_red, sigma=0.0, which="LM",
                         tol=1e-10, maxiter=3000)
    eigsh_taus = sorted((1.0/e) * 1e6 for e in eigs)[::-1]
    print(f"  scipy eigsh tau_n[:{len(eigsh_taus)}] = {[f'{t:.3f}' for t in eigsh_taus]} us")

    res = {}
    res["eigsh_tau_us"] = eigsh_taus
    res["mesh"] = {"ne": mesh.ne, "ndof": fes.ndof, "free": n_free}
    return res


def dx_marker(fes):
    """Helper: NGSolve dx symbol valid in current import scope."""
    from ngsolve import dx
    return dx


def main():
    # Multi-mesh sweep to detect FEM-side convergence
    cases = [
        (0.5e-3, 5e-3, 100e-3, 100e-3, "coarse"),
        (0.3e-3, 5e-3, 200e-3, 200e-3, "medium"),
        (0.2e-3, 10e-3, 200e-3, 200e-3, "fine"),
    ]
    all_results = {}
    for h_d, h_a, R_a, Z_a, lbl in cases:
        all_results[lbl] = solve_disk_modes(
            maxh_disk=h_d, maxh_air=h_a, R_air=R_a, Z_air=Z_a,
            N_stages=6, label=lbl)

    # BEM-Foster reference
    with open(BEM_TAU_REF_PATH) as fp:
        bem = json.load(fp)
    bem_tau = bem["bem_tau_us"]
    print("\n" + "="*70)
    print("Independent scipy eigsh vs BEM-Foster comparison")
    for lbl, res in all_results.items():
        print(lbl, res["mesh"])
        for i, t_eig in enumerate(res["eigsh_tau_us"]):
            t_bem = bem_tau[i] if i < len(bem_tau) else None
            ratio = t_eig / t_bem if t_bem else None
            print(i + 1, t_eig, t_bem, ratio)

    out_path = os.path.join(os.path.dirname(__file__), "test_hiruma_disk_results.json")
    with open(out_path, "w") as fp:
        json.dump({"results": all_results, "bem_tau_us": bem_tau}, fp, indent=2)
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
