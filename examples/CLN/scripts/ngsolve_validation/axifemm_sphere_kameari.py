"""axifemm_sphere_kameari.py — axifemm Cu sphere Cauer-I in Kameari convention.

Cu sphere R=10 mm, sigma=5.8e7 S/m, uniform B_z=1 T.

Mesh: structured Q1 axis-aligned quad with stair-step approximation of the
sphere boundary (cells classified as 'conductor' or 'air' by checking
whether the cell center is inside r^2 + z^2 < R^2). Q1 is the highest-
precision axifemm path (closed-form integration) so we use stair-step on
a fine grid rather than P1 triangle on an OCC curved boundary.

Foster -> Kameari extraction via Schur complement on conductor block.
"""
from __future__ import annotations

import json
import sys
from math import pi
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import scipy.linalg as la
import mpmath as mp

mp.mp.dps = 80

sys.path.insert(0,
    str(Path("S:/Radia/01_GitHub/packages/radia-axifemm/tests")))

import netgen.meshing as ng_meshing
from netgen.meshing import (
    Mesh as NgMesh, Element1D, Element2D, FaceDescriptor, MeshPoint, Pnt,
)
from ngsolve import (
    Mesh, BilinearForm, LinearForm, CoefficientFunction, TaskManager,
    x, dx, ngsglobals,
)
from radia.radia_axifemm import (
    H1Henrotte, AxiHenrotteStiffnessBFI, AxiHenrotteSigmaMassBFI,
)
from test_hiruma_disk_q1 import to_scipy_csr  # type: ignore


R_SPHERE = 10e-3
SIGMA_CU = 5.8e7
MU0 = 4 * pi * 1e-7
B0 = 1.0


def make_sphere_quad_mesh(NR_sphere=80, Nz_sphere=160, NR_air=15, Nz_air=15,
                          R_air=200e-3, Z_air=200e-3):
    """Structured axis-aligned quad mesh for sphere with stair-step boundary.

    Layout in (r, z):
      - Sphere interior region: r in [0, R], z in [-R, R], split into
        NR_sphere x Nz_sphere uniform cells. Cells with center inside
        r^2 + z^2 < R^2 are 'conductor'; others are 'air'.
      - Air right of sphere: r in [R, R_air], geometric grading.
      - Air above sphere: z in [R, Z_air], geometric grading.
      - Air below sphere: z in [-Z_air, -R], geometric grading.
    """
    r_sphere = np.linspace(0, R_SPHERE, NR_sphere + 1)
    r_air = np.geomspace(R_SPHERE, R_air, NR_air + 1)[1:]
    r_grid = np.concatenate([r_sphere, r_air])

    z_sphere = np.linspace(-R_SPHERE, R_SPHERE, Nz_sphere + 1)
    z_above = np.geomspace(R_SPHERE + 1e-7, Z_air, Nz_air + 1)[1:]
    z_below = -np.geomspace(R_SPHERE + 1e-7, Z_air, Nz_air + 1)[1:][::-1]
    z_grid = np.concatenate([z_below, z_sphere, z_above])

    NR = len(r_grid) - 1
    NZ = len(z_grid) - 1

    ngmesh = NgMesh()
    ngmesh.dim = 2
    ngmesh.SetMaterial(1, "air")
    ngmesh.SetMaterial(2, "conductor")
    ngmesh.Add(FaceDescriptor(surfnr=1, domin=0, bc=1))
    ngmesh.Add(FaceDescriptor(surfnr=2, domin=0, bc=2))
    ngmesh.Add(FaceDescriptor(surfnr=3, domin=0, bc=3))
    ngmesh.Add(FaceDescriptor(surfnr=4, domin=0, bc=4))
    ngmesh.SetBCName(0, "axis")
    ngmesh.SetBCName(1, "right")
    ngmesh.SetBCName(2, "top")
    ngmesh.SetBCName(3, "bot")

    pids = np.empty((NZ + 1, NR + 1), dtype=object)
    for j in range(NZ + 1):
        for i in range(NR + 1):
            pids[j, i] = ngmesh.Add(MeshPoint(Pnt(r_grid[i], z_grid[j], 0)))

    n_cond_cells = 0
    for j in range(NZ):
        for i in range(NR):
            r_c = 0.5 * (r_grid[i] + r_grid[i + 1])
            z_c = 0.5 * (z_grid[j] + z_grid[j + 1])
            in_sphere = (r_c ** 2 + z_c ** 2) < R_SPHERE ** 2
            mat = 2 if in_sphere else 1
            if in_sphere:
                n_cond_cells += 1
            ngmesh.Add(Element2D(mat, [
                pids[j,     i    ],
                pids[j,     i + 1],
                pids[j + 1, i + 1],
                pids[j + 1, i    ],
            ]))

    for j in range(NZ):
        ngmesh.Add(Element1D([pids[j, 0], pids[j + 1, 0]], index=1))
    for j in range(NZ):
        ngmesh.Add(Element1D([pids[j, NR], pids[j + 1, NR]], index=2))
    for i in range(NR):
        ngmesh.Add(Element1D([pids[NZ, i], pids[NZ, i + 1]], index=3))
    for i in range(NR):
        ngmesh.Add(Element1D([pids[0, i], pids[0, i + 1]], index=4))

    print(f"  mesh: NR={NR}, NZ={NZ}, total quads={NR*NZ}, "
          f"cond cells={n_cond_cells}")
    return Mesh(ngmesh)


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


def schur_foster_kameari(K_dense, M_dense, b, cond_idx, air_idx, label=""):
    print(f"\n  [{label}] cond DOFs={len(cond_idx)}, air DOFs={len(air_idx)}")
    K_aa = K_dense[air_idx[:, None], air_idx[None, :]]
    K_ac = K_dense[air_idx[:, None], cond_idx[None, :]]
    K_ca = K_dense[cond_idx[:, None], air_idx[None, :]]
    K_cc = K_dense[cond_idx[:, None], cond_idx[None, :]]
    M_cc = M_dense[cond_idx[:, None], cond_idx[None, :]]
    b_c = b[cond_idx]

    print(f"  Building Schur complement ...")
    K_aa_inv_K_ac = la.solve(K_aa, K_ac, assume_a="pos")
    S = K_cc - K_ca @ K_aa_inv_K_ac
    S = 0.5 * (S + S.T)

    eigvals, eigvecs = la.eigh(S, M_cc)
    for k in range(len(eigvals)):
        v = eigvecs[:, k]
        nrm = float(v @ (M_cc @ v))
        if nrm > 1e-30:
            eigvecs[:, k] = v / np.sqrt(nrm)

    tau_seconds_all = np.where(eigvals > 1e-30, 1.0 / eigvals, 0.0)
    tau_us = tau_seconds_all * 1e6
    g = eigvecs.T @ b_c
    g2 = g * g

    order = np.argsort(-g2)
    print(f"  Top 10 modes by |g|^2:")
    for rank in range(min(10, len(order))):
        idx = order[rank]
        print(f"    rank={rank}: tau={tau_us[idx]:.4f} us, "
              f"g2={g2[idx]:.4e}")

    n_use = min(80, int(np.sum(g2 > 1e-32)))
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
    print(f"  Kameari Cauer-I rungs:")
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
    }


def main():
    cases = [
        (40,  80, 12, 12, 100e-3, 100e-3, "medium"),
        (80, 160, 15, 15, 200e-3, 200e-3, "fine"),
    ]
    all_res = {}
    for NR_s, Nz_s, NR_a, Nz_a, R_a, Z_a, lbl in cases:
        print(f"\n=== axifemm Kameari Cu sphere: NR_s={NR_s} Nz_s={Nz_s} "
              f"NR_a={NR_a} Nz_a={Nz_a} R_air={R_a*1e3}mm  {lbl} ===")
        ngsglobals.msg_level = 0
        mesh = make_sphere_quad_mesh(NR_s, Nz_s, NR_a, Nz_a, R_a, Z_a)
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

        K_dense = K_red.toarray()
        M_dense = M_red.toarray()
        M_diag = np.diag(M_dense)
        cond_mask = np.abs(M_diag) > 1e-25
        cond_idx = np.where(cond_mask)[0]
        air_idx = np.where(~cond_mask)[0]

        all_res[lbl] = schur_foster_kameari(K_dense, M_dense, b_red,
                                            cond_idx, air_idx, label=lbl)

    print("\n" + "="*78)
    print("Summary: axifemm Q1 Kameari Cauer-I (sphere Cu)")
    print("="*78)
    print(f"  Stoll analytical leading tau (Mathematica HP) = 694.14 us")
    for lbl, r in all_res.items():
        rungs = r["Cauer_rungs_Kameari"]
        if rungs:
            tau0 = rungs[0]["tau_pair_us"]
            print(f"  {lbl:<10} axifemm Q1 sphere Kameari "
                  f"tau_pair[0] = {tau0:.4f} us")

    out_path = (Path(__file__).parent
                / "axifemm_sphere_kameari_results.json")
    out_path.write_text(json.dumps(all_res, indent=2), encoding="utf-8")
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
