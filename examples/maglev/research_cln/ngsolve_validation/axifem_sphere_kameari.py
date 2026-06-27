"""axifem_sphere_kameari.py — axifem Cu sphere Cauer-I in Kameari convention.

Cu sphere R=10 mm, sigma=5.8e7 S/m, uniform B_z=1 T.

Mesh: structured Q1 axis-aligned quad with stair-step approximation of the
sphere boundary (cells classified as 'conductor' or 'air' by checking
whether the cell center is inside r^2 + z^2 < R^2). Q1 is the highest-
precision axifem path (closed-form integration) so we use stair-step on
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
import scipy.sparse as sp
import scipy.sparse.linalg as sla
import mpmath as mp

mp.mp.dps = 80

sys.path.insert(0,
    str(Path("S:/Radia/01_GitHub/validation_test/axifem/research/verification")))

import netgen.meshing as ng_meshing
from netgen.meshing import (
    Mesh as NgMesh, Element1D, Element2D, FaceDescriptor, MeshPoint, Pnt,
)
from ngsolve import (
    Mesh, BilinearForm, LinearForm, CoefficientFunction, TaskManager,
    x, dx, ngsglobals,
)
from radia.axifem import (
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


def _generalized_eigh_gpu(S_np, M_np):
    """Generalized eigh via Cholesky transform with GPU/CPU dispatch.

    Solve S v = λ M v with M SPD: decompose M = L L^T, transform
    to standard form C = L^{-1} S L^{-T}, eigh(C) → λ, w; recover
    v = L^{-T} w. Eigenvectors satisfy v^T M v = 1.

    Strategy: do Cholesky + triangular solves on GPU (cheap), then
    eigh of standard-form C on CPU because cuSOLVER eigh needs
    ~5x matrix memory and A100 40GB OOMs at N>~25000. GPU is still
    used for the most expensive non-eigh kernels.
    """
    import cupy as cp
    from cupyx.scipy.linalg import solve_triangular as cp_solve_tri

    S_gpu = cp.asarray(S_np)
    M_gpu = cp.asarray(M_np)
    print(f"  [GPU] Cholesky of M ({M_gpu.shape[0]}x{M_gpu.shape[1]}) ...")
    L_gpu = cp.linalg.cholesky(M_gpu)
    del M_gpu
    print(f"  [GPU] Transform to standard form C = L^-1 S L^-T ...")
    Linv_S = cp_solve_tri(L_gpu, S_gpu, lower=True)
    del S_gpu
    C_gpu = cp_solve_tri(L_gpu, Linv_S.T, lower=True).T
    del Linv_S
    C_gpu = 0.5 * (C_gpu + C_gpu.T)
    # eigh on C: try GPU first, fall back to CPU if cuSOLVER OOM.
    N = C_gpu.shape[0]
    try:
        print(f"  [GPU] eigh of C ({N}x{N}) ...")
        eigvals_gpu, W_gpu = cp.linalg.eigh(C_gpu)
        del C_gpu
        print(f"  [GPU] Back-transform v = L^-T w ...")
        V_gpu = cp_solve_tri(L_gpu.T, W_gpu, lower=False)
        del L_gpu, W_gpu
        eigvals = cp.asnumpy(eigvals_gpu)
        eigvecs = cp.asnumpy(V_gpu)
        del eigvals_gpu, V_gpu
    except Exception as e:
        print(f"  [GPU eigh failed: {e}] falling back to CPU eigh on C ...")
        C_cpu = cp.asnumpy(C_gpu)
        del C_gpu
        cp.get_default_memory_pool().free_all_blocks()
        print(f"  [CPU] eigh of C ({N}x{N}) ... (may take 30-60 min)")
        eigvals, W_cpu = la.eigh(C_cpu, overwrite_a=True)
        del C_cpu
        print(f"  [GPU] Back-transform v = L^-T w (GPU) ...")
        W_gpu = cp.asarray(W_cpu)
        del W_cpu
        V_gpu = cp_solve_tri(L_gpu.T, W_gpu, lower=False)
        del L_gpu, W_gpu
        eigvecs = cp.asnumpy(V_gpu)
        del V_gpu
    cp.get_default_memory_pool().free_all_blocks()
    return eigvals, eigvecs


def schur_foster_kameari(K_sp, M_sp, b, cond_idx, air_idx,
                          label="", use_gpu=True):
    """Foster spectrum via dense Schur complement + Cholesky-transform eigh.

    K_sp, M_sp are sparse matrices on the free DOF set. Block submatrices
    are extracted as sparse then converted to dense one at a time so we
    avoid allocating the full free-DOF dense matrices (saves ~80 GB for
    Q2 fine: 71631² × 8 × 2 = 82 GB).
    """
    print(f"\n  [{label}] cond DOFs={len(cond_idx)}, air DOFs={len(air_idx)}")
    # K_aa is sparse FE stiffness; keep it sparse and use sparse splu
    # for K_aa^{-1} K_ac. Sparse LU fill-in (~N^1.5 entries) is far smaller
    # than the dense 7.8 GB Cholesky workspace that OOMs at Q2 fine.
    print(f"  Extracting K_aa block ({len(air_idx)}^2) sparse ...")
    K_aa_sp = sp.csc_matrix(K_sp[air_idx[:, None], air_idx[None, :]])
    print(f"  K_aa sparse: nnz={K_aa_sp.nnz}, mem={K_aa_sp.data.nbytes/1e6:.1f} MB")
    print(f"  Extracting K_ac block ({len(air_idx)}x{len(cond_idx)}) dense ...")
    K_ac = np.asarray(K_sp[air_idx[:, None], cond_idx[None, :]].todense())
    print(f"  Sparse splu(K_aa) factorization ...")
    K_aa_lu = sla.splu(K_aa_sp)
    print(f"  Sparse solve Y = K_aa^-1 K_ac (full dense RHS) ...")
    Y = K_aa_lu.solve(K_ac)
    del K_aa_lu, K_aa_sp
    print(f"  Extracting K_cc block ({len(cond_idx)}^2) ...")
    K_cc = np.asarray(K_sp[cond_idx[:, None], cond_idx[None, :]].todense())
    print(f"  Building S = K_cc - K_ac.T @ Y ({len(cond_idx)}^2 GEMM) ...")
    S = K_cc - K_ac.T @ Y
    S = 0.5 * (S + S.T)
    del K_ac, Y, K_cc
    print(f"  Extracting M_cc block ({len(cond_idx)}^2) ...")
    M_cc = np.asarray(M_sp[cond_idx[:, None], cond_idx[None, :]].todense())
    b_c = b[cond_idx]

    # Adaptive Tikhonov: try eigh with progressively larger α until Cholesky
    # succeeds. With cell-based cond_idx, this loop typically exits at α=0.
    # In-place diagonal modification avoids np.eye(N) which is 13 GB at N=40k.
    m_max = float(np.max(np.abs(np.diag(M_cc))))
    orig_diag = np.diag(M_cc).copy()
    eigvals = None
    eigvecs = None
    for alpha_rel in [0.0, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4]:
        alpha = alpha_rel * m_max
        if alpha > 0:
            np.fill_diagonal(M_cc, orig_diag + alpha)
        else:
            np.fill_diagonal(M_cc, orig_diag)
        try:
            if use_gpu:
                eigvals, eigvecs = _generalized_eigh_gpu(S, M_cc)
            else:
                print(f"  CPU eigh (α_rel={alpha_rel:.0e}) ...")
                eigvals, eigvecs = la.eigh(S, M_cc)
            print(f"  [eigh OK] α_rel={alpha_rel:.0e}")
            break
        except Exception as e:
            print(f"  [eigh failed] α_rel={alpha_rel:.0e}: {e}")
            continue
    # Restore M_cc for subsequent M-norm normalization
    np.fill_diagonal(M_cc, orig_diag)
    if eigvals is None:
        raise RuntimeError("All Tikhonov α values failed; "
                            "M_cc fundamentally indefinite")

    # M-norm renormalize (eigh from Cholesky transform already satisfies this,
    # but recompute for robustness against GPU FP rounding)
    for k in range(eigvecs.shape[1]):
        v = eigvecs[:, k]
        nrm = float(v @ (M_cc @ v))
        if nrm > 1e-30:
            eigvecs[:, k] = v / np.sqrt(nrm)

    tau_seconds_all = np.where(np.abs(eigvals) > 1e-30, 1.0 / eigvals, 0.0)
    tau_us = tau_seconds_all * 1e6
    g = eigvecs.T @ b_c
    g2 = g * g

    order = np.argsort(-g2)
    g2_total = float(np.sum(g2))
    foster_modes = []
    print(f"  Top 15 modes by |g|^2:")
    for rank in range(min(15, len(order))):
        idx = order[rank]
        print(f"    rank={rank}: tau={tau_us[idx]:.4f} us, "
              f"g2={g2[idx]:.4e}")
        foster_modes.append({
            "rank": rank,
            "tau_us": float(tau_us[idx]),
            "g2": float(g2[idx]),
            "g2_normalized": float(g2[idx]) / g2_total if g2_total > 0 else 0.0,
        })

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
        "top_modes": foster_modes,
        "Cauer_rungs_Kameari": rungs,
    }


def main():
    # Q1 default (order=1): dense Schur + Cholesky-transform GPU eigh.
    # medium = NR=40/Nz=80, fine = NR=80/Nz=160. Q2 (order=2) was attempted
    # as an upgrade but Q2 fine eigh exceeds A100 40GB cuSOLVER memory and
    # Cauer high-k extraction degrades; documented as future work in the
    # supplement paper. Q1 fine remains the production axifem path for
    # sphere Cauer extraction.
    cases = [
        (40,  80, 12, 12, 100e-3, 100e-3, "medium"),
        (80, 160, 15, 15, 200e-3, 200e-3, "fine"),
    ]
    all_res = {}
    for NR_s, Nz_s, NR_a, Nz_a, R_a, Z_a, lbl in cases:
        print(f"\n=== axifem Q2 Kameari Cu sphere: NR_s={NR_s} Nz_s={Nz_s} "
              f"NR_a={NR_a} Nz_a={Nz_a} R_air={R_a*1e3}mm  {lbl} ===")
        ngsglobals.msg_level = 0
        mesh = make_sphere_quad_mesh(NR_s, Nz_s, NR_a, Nz_a, R_a, Z_a)
        fes = H1Henrotte(mesh, order=1, dirichlet="axis|right|top|bot")
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

        # Cell-based cond classification: a DOF is "conductor" iff at
        # least one of its supporting elements is in the conductor
        # material. This is the FEM-correct classification and yields
        # a strictly PD M_cc submatrix on cond_idx (avoids the M_diag>0
        # threshold ambiguity that broke Q2 Cholesky).
        from ngsolve import VOL
        materials = mesh.GetMaterials()
        cond_global = set()
        for el in mesh.Elements(VOL):
            mat_name = materials[el.mat] if isinstance(el.mat, int) else el.mat
            if str(mat_name).startswith("conductor"):
                for d in fes.GetDofNrs(el):
                    if d >= 0:
                        cond_global.add(int(d))
        # Map global DOF index -> position in `free` array
        free_to_local = {int(g): i for i, g in enumerate(free)}
        cond_idx = np.array(sorted(free_to_local[d] for d in cond_global
                                    if d in free_to_local), dtype=int)
        air_mask = np.ones(len(free), dtype=bool)
        air_mask[cond_idx] = False
        air_idx = np.where(air_mask)[0]
        print(f"  Cell-based cond_idx: {len(cond_idx)} DOFs "
              f"(air_idx: {len(air_idx)} DOFs)")

        all_res[lbl] = schur_foster_kameari(
            K_red, M_red, b_red, cond_idx, air_idx,
            label=lbl, use_gpu=True)

    print("\n" + "="*78)
    print("Summary: axifem Q2 Kameari Cauer-I (sphere Cu)")
    print("="*78)
    print(f"  Stoll analytical leading tau (Mathematica HP) = 694.14 us")
    for lbl, r in all_res.items():
        rungs = r["Cauer_rungs_Kameari"]
        if rungs:
            tau0 = rungs[0]["tau_pair_us"]
            print(f"  {lbl:<10} axifem Q2 sphere Kameari "
                  f"tau_pair[0] = {tau0:.4f} us")

    out_path = (Path(__file__).parent
                / "axifem_sphere_kameari_results.json")
    out_path.write_text(json.dumps(all_res, indent=2), encoding="utf-8")
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
