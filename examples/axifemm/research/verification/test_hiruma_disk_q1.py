"""Phase 3 follow-up: Hiruma 3-term on Cu disk with **structured Q1 quad mesh**.

Goal: close the FEM-side accuracy gap from 2% (TRIG mesh) to 0.6% (Q1 quad
mesh, matching axifemm/disk_hiruma_quad.py Python prototype validated against
BEM-Foster Mathematica reference at 224.31 us).

Builds the netgen 2D mesh manually with axis-aligned rectangles so that
AxiHenrotteFE_Q1_AxisAligned is dispatched (instead of AxiHenrotteFE_P1_Triangle
which the OCC mesher produces by default).
"""

import json
import os
from math import pi

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

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


R_DISK = 10e-3
T_DISK = 2e-3
SIGMA_CU = 5.8e7
MU0 = 4 * pi * 1e-7
B0 = 1.0


def make_structured_disk_quad_mesh(NR_disk=40, Nz_disk=8, NR_air=15, Nz_air=15,
                                    R_air=200e-3, Z_air=200e-3):
    """Mirror axifemm/axifemm_quad.py:structured_disk_mesh, output as NGSolve.

    Layout in (r, z):
        - Disk: r in [0, R_disk], z in [-T/2, T/2], NR_disk x Nz_disk uniform cells.
        - Air right of disk: r in [R_disk, R_air], geometric grading.
        - Air above disk: z in [T/2, Z_air], geometric grading.
        - Air below disk: z in [-Z_air, -T/2], geometric grading.

    Boundary names: "axis" (r=0), "right" (r=R_air), "top" (z=Z_air), "bot" (z=-Z_air).
    Material names: "conductor" (in disk), "air" (elsewhere).
    """
    r_disk = np.linspace(0, R_DISK, NR_disk + 1)
    r_air  = np.geomspace(R_DISK, R_air, NR_air + 1)[1:]
    r_grid = np.concatenate([r_disk, r_air])

    z_disk = np.linspace(-T_DISK / 2, T_DISK / 2, Nz_disk + 1)
    z_above = np.geomspace(T_DISK/2 + 1e-7, Z_air, Nz_air + 1)[1:]
    z_below = -np.geomspace(T_DISK/2 + 1e-7, Z_air, Nz_air + 1)[1:][::-1]
    z_grid = np.concatenate([z_below, z_disk, z_above])

    NR = len(r_grid) - 1
    NZ = len(z_grid) - 1

    ngmesh = NgMesh()
    ngmesh.dim = 2

    # Materials: index 1 = air, 2 = conductor.
    ngmesh.SetMaterial(1, "air")
    ngmesh.SetMaterial(2, "conductor")

    # FaceDescriptors for the 4 boundary edges.
    fd_axis  = ngmesh.Add(FaceDescriptor(surfnr=1, domin=0, bc=1))
    fd_right = ngmesh.Add(FaceDescriptor(surfnr=2, domin=0, bc=2))
    fd_top   = ngmesh.Add(FaceDescriptor(surfnr=3, domin=0, bc=3))
    fd_bot   = ngmesh.Add(FaceDescriptor(surfnr=4, domin=0, bc=4))
    ngmesh.SetBCName(0, "axis")
    ngmesh.SetBCName(1, "right")
    ngmesh.SetBCName(2, "top")
    ngmesh.SetBCName(3, "bot")

    # Add nodes.
    pids = np.empty((NZ + 1, NR + 1), dtype=object)
    for j in range(NZ + 1):
        for i in range(NR + 1):
            pids[j, i] = ngmesh.Add(MeshPoint(Pnt(r_grid[i], z_grid[j], 0)))

    # Add quads with material index. Vertex order MUST be CCW in (r, z):
    # (i,j) -> (i+1,j) -> (i+1,j+1) -> (i,j+1).
    for j in range(NZ):
        for i in range(NR):
            r_c = 0.5 * (r_grid[i] + r_grid[i + 1])
            z_c = 0.5 * (z_grid[j] + z_grid[j + 1])
            in_disk = (r_c < R_DISK) and (abs(z_c) < T_DISK / 2)
            mat_index = 2 if in_disk else 1
            ngmesh.Add(Element2D(mat_index, [
                pids[j,     i    ],
                pids[j,     i + 1],
                pids[j + 1, i + 1],
                pids[j + 1, i    ],
            ]))

    # Add boundary edges. Each Element1D needs a face index pointing at a
    # FaceDescriptor (whose bc index gives the BC name).
    # Axis: r=0 -> column i=0
    for j in range(NZ):
        ngmesh.Add(Element1D([pids[j, 0], pids[j + 1, 0]], index=1))
    # Right: r=R_air -> column i=NR
    for j in range(NZ):
        ngmesh.Add(Element1D([pids[j, NR], pids[j + 1, NR]], index=2))
    # Top: z=Z_air -> row j=NZ
    for i in range(NR):
        ngmesh.Add(Element1D([pids[NZ, i], pids[NZ, i + 1]], index=3))
    # Bot: z=-Z_air -> row j=0
    for i in range(NR):
        ngmesh.Add(Element1D([pids[0, i], pids[0, i + 1]], index=4))

    return Mesh(ngmesh)


def to_scipy_csr(mat, n):
    rs, cs, vs = mat.COO()
    K = sp.csr_matrix(
        (np.asarray(vs, dtype=np.float64),
         (np.asarray(rs, dtype=np.int64), np.asarray(cs, dtype=np.int64))),
        shape=(n, n))
    return (K + K.T) * 0.5


def hiruma_3term(K, M, b, N_stages, label=""):
    """Hiruma 3-term Lanczos recursion -> Cauer ladder R_{2k}, L_{2k+1}.

    The Cauer ladder follows Nagamine et al. 2026 Fig. 3 / Eq. (11):
        in --R_0--+--R_2--+--R_4--+-...
                  |        |        |
                  L_1      L_3      L_5  ...
                  |        |        |
                  gnd      gnd      gnd
    -- R_{2k} are the series resistors (even subscripts) and L_{2k+1} are
    the shunt inductors (odd subscripts). Per-pair time constant:
        tau_pair[k] = L_{2k+1} / R_{2k}

    The Hiruma 3-term recursion produces the K^{-1} M Krylov-orthogonal
    basis {w_1, w_2, ...}. The diagonal coefficients lam[i] = w_i^T A w_i
    map to the Cauer ladder via:
        lam[2k+1] = w_{2k+1}^T K w_{2k+1} = 1 / R_{2k}      (conductance)
        lam[2k+2] = w_{2k+2}^T M w_{2k+2} =     L_{2k+1}    (inductance)
    so that  tau_pair[k] = L_{2k+1} / R_{2k} = lam[2k+1] * lam[2k+2].
    """
    factor = spla.factorized(K.tocsc())

    w = [None] * (2 * N_stages + 3)
    lam = [None] * (2 * N_stages + 3)

    w[1] = factor(b)
    lam[1] = float(w[1] @ (K @ w[1]))           # = 1 / R_0
    w[2] = -w[1] / lam[1]

    print(f"\n  [{label}] 1/R_0 = lam_1 = {lam[1]:.6e}  =>  R_0 = {1.0/lam[1]:.6e}")
    print(f"  {'k':>3} {'R_{2k}':>13} {'L_{2k+1}':>13} "
          f"{'tau_pair us':>13} {'orth_drift':>12}")

    out = []
    for n in range(1, N_stages + 1):
        # n is the iteration number; Nagamine k-pair index is k = n-1 the
        # FIRST iteration produces (R_0, L_1) corresponding to k=0; n-th
        # iteration produces (R_{2(n-1)}, L_{2n-1}).
        lam[2*n] = float(w[2*n] @ (M @ w[2*n]))         # = L_{2n-1}
        Cw = M @ w[2*n]
        Aw = -factor(Cw)
        w[2*n+1] = w[2*n-1].copy() - Aw / lam[2*n]
        lam[2*n+1] = float(w[2*n+1] @ (K @ w[2*n+1]))    # = 1 / R_{2n}

        # Cauer values for THIS iteration (Nagamine k-pair index = n-1)
        R_2k = 1.0 / lam[2*n - 1]                # R_{2(n-1)} (already known
                                                 # from prev iter or init)
        L_2k_plus_1 = lam[2*n]                   # L_{2(n-1)+1} = L_{2n-1}
        # tau_pair for the (n-1)-th Nagamine pair:
        tau_pair_us = (L_2k_plus_1 / R_2k * 1e6
                       if R_2k > 0 else None)

        drift = 0.0
        for kk in range(1, n):
            cross = float(w[2*kk+1] @ (K @ w[2*n+1]))
            drift = max(drift, abs(cross) / abs(lam[2*n+1]))

        w[2*n+2] = w[2*n].copy() - w[2*n+1] / lam[2*n+1]

        ts = f"{tau_pair_us:.4f}" if tau_pair_us is not None else "N/A"
        k_idx = n - 1                             # Nagamine pair index
        print(f"  {k_idx:>3} {R_2k:>13.6e} {L_2k_plus_1:>13.6e} "
              f"{ts:>13} {drift:>12.3e}")

        out.append({"k": k_idx,
                    "R_2k": R_2k,
                    "L_2k_plus_1": L_2k_plus_1,
                    "tau_pair_us": tau_pair_us,
                    "orth_drift": drift,
                    # Back-compat aliases (n indexed from 1)
                    "n": n, "L_n": lam[2*n], "inv_R_n_plus_1": lam[2*n+1],
                    "tau_per_stage_us": tau_pair_us})
    return {"lam_1": lam[1], "R_0": 1.0 / lam[1], "stages": out}


def solve_disk(NR_disk, Nz_disk, NR_air, Nz_air, R_air, Z_air, N_stages=6,
               label=""):
    print(f"\n=== Q1 quad: NR_d={NR_disk} Nz_d={Nz_disk} NR_a={NR_air} "
          f"Nz_a={Nz_air} R_air={R_air*1e3}mm  {label} ===")
    ngsglobals.msg_level = 0
    mesh = make_structured_disk_quad_mesh(NR_disk, Nz_disk, NR_air, Nz_air,
                                            R_air, Z_air)
    fes = H1Henrotte(mesh, dirichlet="axis|right|top|bot")
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
          f"{[f'{t:.3f}' for t in eigsh_taus]} us")

    res = hiruma_3term(K_red, M_red, b_red, N_stages, label=label)
    res["eigsh_tau_us"] = eigsh_taus
    res["mesh"] = {"ne": mesh.ne, "ndof": fes.ndof, "free": int(len(free))}
    return res


BEM_TAU = [224.3070587702379, 88.42395618426703, 49.53986759314118,
           32.090204181064834, 23.339598470067966, 22.588437092711622]
PYTHON_Q1_TAU = [223.06, 87.60, 48.86, 31.49]


def main():
    cases = [
        (20,  4,  10, 10, 100e-3, 100e-3, "coarse"),
        (40,  8,  15, 15, 200e-3, 200e-3, "medium"),
        (80, 16,  20, 20, 500e-3, 500e-3, "fine"),
        (160, 32, 25, 25, 500e-3, 500e-3, "very fine"),
    ]
    all_res = {}
    for NR_d, Nz_d, NR_a, Nz_a, R_a, Z_a, lbl in cases:
        all_res[lbl] = solve_disk(NR_d, Nz_d, NR_a, Nz_a, R_a, Z_a,
                                    N_stages=6, label=lbl)

    print("\n" + "="*78)
    print("Summary: Q1 quad mesh vs BEM-Foster (Mathematica reference)")
    print("="*78)
    print(f"{'case':<10} {'ne':>8} {'free':>6}  "
          f"{'τ_1 eigsh':>10}  {'/BEM':>7}  {'/PyQ1':>7}")
    for lbl, r in all_res.items():
        t1 = r["eigsh_tau_us"][0]
        ratio_bem = t1 / BEM_TAU[0]
        ratio_py = t1 / PYTHON_Q1_TAU[0]
        print(f"{lbl:<10} {r['mesh']['ne']:>8} {r['mesh']['free']:>6}  "
              f"{t1:>10.4f}  {ratio_bem:>7.4f}  {ratio_py:>7.4f}")

    out_path = os.path.join(os.path.dirname(__file__),
                            "test_hiruma_disk_q1_results.json")
    with open(out_path, "w") as fp:
        json.dump({"results": all_res, "bem_tau_us": BEM_TAU,
                   "python_q1_tau_us": PYTHON_Q1_TAU}, fp, indent=2)
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
