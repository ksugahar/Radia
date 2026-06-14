"""Phase A2: Hiruma 3-term on Cu disk with **Q2 quad mesh**.

Goal: close the FEM-side accuracy gap from 0.55% (Q1 quad, very fine mesh,
ne=15170) to << 0.1% with Q2 (9 DOFs / quad: 4 vertex + 4 edge + 1 face),
matching the BEM-Foster Mathematica reference at 224.31 us.

Reuses the structured quad mesh builder from test_hiruma_disk_q1; the only
change is `order=2` on H1Henrotte.
"""

import json
import os
from math import pi

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from ngsolve import (
    Mesh, BilinearForm, LinearForm, CoefficientFunction, TaskManager,
    x, dx, ngsglobals,
)
from radia.axifem import (
    H1Henrotte, AxiHenrotteStiffnessBFI, AxiHenrotteSigmaMassBFI,
)

from test_hiruma_disk_q1 import (  # type: ignore
    make_structured_disk_quad_mesh,
    to_scipy_csr,
    hiruma_3term,
    R_DISK, T_DISK, SIGMA_CU, MU0, B0,
    BEM_TAU, PYTHON_Q1_TAU,
)


def solve_disk_q2(NR_disk, Nz_disk, NR_air, Nz_air, R_air, Z_air, N_stages=6,
                  label=""):
    print(f"\n=== Q2 quad: NR_d={NR_disk} Nz_d={Nz_disk} NR_a={NR_air} "
          f"Nz_a={Nz_air} R_air={R_air*1e3}mm  {label} ===")
    ngsglobals.msg_level = 0
    mesh = make_structured_disk_quad_mesh(NR_disk, Nz_disk, NR_air, Nz_air,
                                          R_air, Z_air)
    fes = H1Henrotte(mesh, order=2, dirichlet="axis|right|top|bot")
    n_free = sum(1 for f in fes.FreeDofs() if f)
    print(f"  mesh ne={mesh.ne} ndof={fes.ndof}  free={n_free}  "
          f"materials={mesh.GetMaterials()}")

    mu_cf = CoefficientFunction(MU0)
    sigma_cf = mesh.MaterialCF({"conductor": SIGMA_CU}, default=0.0)
    A_imposed = B0 * x / 2

    a = BilinearForm(fes, symmetric=True)
    a += AxiHenrotteStiffnessBFI(mu_cf)
    with TaskManager(): a.Assemble()

    m = BilinearForm(fes, symmetric=True)
    m += AxiHenrotteSigmaMassBFI(sigma_cf)
    with TaskManager(): m.Assemble()

    b_form = LinearForm(fes)
    v = fes.TestFunction()
    b_form += sigma_cf * A_imposed * v * 2 * pi * x * dx
    with TaskManager(): b_form.Assemble()
    b_vec = np.array(b_form.vec)

    K_csr = to_scipy_csr(a.mat, fes.ndof)
    M_csr = to_scipy_csr(m.mat, fes.ndof)
    free = np.array([i for i in range(fes.ndof) if fes.FreeDofs()[i]], dtype=int)
    K_red = K_csr[free[:, None], free[None, :]]
    M_red = M_csr[free[:, None], free[None, :]]
    b_red = b_vec[free]

    eigs, _ = spla.eigsh(K_red, k=min(N_stages, len(free) // 2),
                         M=M_red, sigma=0.0, which="LM",
                         tol=1e-10, maxiter=3000)
    eigsh_taus = sorted((1.0/e) * 1e6 for e in eigs)[::-1]
    print(f"  scipy eigsh tau_n[:{len(eigsh_taus)}] = "
          f"{[f'{t:.4f}' for t in eigsh_taus]} us")

    res = hiruma_3term(K_red, M_red, b_red, N_stages, label=label)
    res["eigsh_tau_us"] = eigsh_taus
    res["mesh"] = {"ne": mesh.ne, "ndof": fes.ndof, "free": int(len(free))}
    return res


def main():
    # Q2 has 9 DOFs/quad (vs 4 for Q1) so we start much coarser. Mesh costs
    # scale with ndof, not ne; ne ~ ndof / 9 for Q2 vs ne ~ ndof for Q1, so
    # Q2 with NR_d=20 has roughly the same DOF count as Q1 NR_d=40.
    cases = [
        (10, 4,  8,  8, 100e-3, 100e-3, "coarse"),
        (20, 8, 12, 12, 200e-3, 200e-3, "medium"),
        (40, 16, 15, 15, 500e-3, 500e-3, "fine"),
    ]
    all_res = {}
    for NR_d, Nz_d, NR_a, Nz_a, R_a, Z_a, lbl in cases:
        all_res[lbl] = solve_disk_q2(NR_d, Nz_d, NR_a, Nz_a, R_a, Z_a,
                                     N_stages=6, label=lbl)

    print("\n" + "="*78)
    print("Summary: Q2 quad mesh vs BEM-Foster (Mathematica reference 224.307 us)")
    print("="*78)
    print(f"{'case':<10} {'ne':>8} {'free':>6}  "
          f"{'τ_1 eigsh':>10}  {'/BEM':>7}  {'/PyQ1':>7}  {'gap%':>7}")
    for lbl, r in all_res.items():
        t1 = r["eigsh_tau_us"][0]
        ratio_bem = t1 / BEM_TAU[0]
        ratio_py = t1 / PYTHON_Q1_TAU[0]
        gap_pct = (1.0 - ratio_bem) * 100
        print(f"{lbl:<10} {r['mesh']['ne']:>8} {r['mesh']['free']:>6}  "
              f"{t1:>10.4f}  {ratio_bem:>7.4f}  {ratio_py:>7.4f}  "
              f"{gap_pct:>+7.3f}")

    out_path = os.path.join(os.path.dirname(__file__),
                            "test_hiruma_disk_q2_results.json")
    with open(out_path, "w") as fp:
        json.dump({"results": all_res, "bem_tau_us": BEM_TAU,
                   "python_q1_tau_us": PYTHON_Q1_TAU}, fp, indent=2)
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
