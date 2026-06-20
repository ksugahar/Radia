"""test_p2_disk_cpp.py — C++ P2 disk benchmark (no air, no Schur).

Compares C++ NGSolve+H1Henrotte order=2 vs Python test_p2_disk_convergence
on the same structured disk mesh. Cu disk only — no air box. This isolates
whether the C++ assembly itself matches Python, independent of Schur logic.
"""
from __future__ import annotations
import sys
sys.stdout.reconfigure(encoding="utf-8")

from math import pi
import time
import numpy as np
import scipy.sparse as sp
import scipy.linalg as sla

import netgen.meshing as ng_meshing
from netgen.meshing import (
    Mesh as NgMesh, Element1D, Element2D, FaceDescriptor, MeshPoint, Pnt,
)
from ngsolve import (
    Mesh, BilinearForm, CoefficientFunction, TaskManager, ngsglobals,
)
from radia.axifem import (
    H1Henrotte, AxiHenrotteStiffnessBFI, AxiHenrotteSigmaMassBFI,
)

R_DISK = 10.0e-3
T_DISK = 2.0e-3
SIGMA = 5.8e7
MU0 = 4 * pi * 1e-7


def make_disk_mesh(N_r, N_z):
    r_pts = np.linspace(0.0, R_DISK, N_r + 1)
    z_pts = np.linspace(-T_DISK / 2.0, T_DISK / 2.0, N_z + 1)

    m = NgMesh()
    m.dim = 2
    m.SetMaterial(1, "conductor")
    # 4 FaceDescriptors with distinct bc, like working Q1 disk test.
    m.Add(FaceDescriptor(surfnr=1, domin=0, bc=1))  # axis
    m.Add(FaceDescriptor(surfnr=2, domin=0, bc=2))  # right
    m.Add(FaceDescriptor(surfnr=3, domin=0, bc=3))  # top
    m.Add(FaceDescriptor(surfnr=4, domin=0, bc=4))  # bot
    m.SetBCName(0, "axis")
    m.SetBCName(1, "right")
    m.SetBCName(2, "top")
    m.SetBCName(3, "bot")

    pids = np.empty((N_z + 1, N_r + 1), dtype=object)
    for j, z in enumerate(z_pts):
        for i, r in enumerate(r_pts):
            pids[j, i] = m.Add(MeshPoint(Pnt(r, z, 0)))

    for j in range(N_z):
        for i in range(N_r):
            m.Add(Element2D(1, [pids[j, i], pids[j, i+1], pids[j+1, i+1]]))
            m.Add(Element2D(1, [pids[j, i], pids[j+1, i+1], pids[j+1, i]]))

    # Boundary 1D edges. Axis (i=0), right (i=N_r), top (j=N_z), bot (j=0).
    for j in range(N_z):
        m.Add(Element1D([pids[j, 0],     pids[j+1, 0]],     index=1))  # axis
        m.Add(Element1D([pids[j, N_r],   pids[j+1, N_r]],   index=2))  # right
    for i in range(N_r):
        m.Add(Element1D([pids[N_z, i],   pids[N_z, i+1]],   index=3))  # top
        m.Add(Element1D([pids[0, i],     pids[0, i+1]],     index=4))  # bot

    return Mesh(m)


def run_one(N_r, N_z):
    ngsglobals.msg_level = 0
    mesh = make_disk_mesh(N_r, N_z)
    fes = H1Henrotte(mesh, order=2, dirichlet="right|top|bot")
    n_free = sum(1 for f in fes.FreeDofs() if f)
    print(f"  Mesh ne={mesh.ne} ndof={fes.ndof} free={n_free}")

    a = BilinearForm(fes, symmetric=True, check_unused=False)
    a += AxiHenrotteStiffnessBFI(CoefficientFunction(MU0))
    t0 = time.time()
    with TaskManager(): a.Assemble()
    t_a = time.time() - t0

    sigma_cf = mesh.MaterialCF({"conductor": SIGMA}, default=0.0)
    m_bf = BilinearForm(fes, symmetric=True, check_unused=False)
    m_bf += AxiHenrotteSigmaMassBFI(sigma_cf)
    t0 = time.time()
    with TaskManager(): m_bf.Assemble()
    t_m = time.time() - t0

    free = np.array([i for i in range(fes.ndof) if fes.FreeDofs()[i]], dtype=int)
    def csr_dense(bfmat, n):
        rs, cs, vs = bfmat.COO()
        return sp.coo_matrix((np.array(vs), (np.array(rs), np.array(cs))),
                             shape=(n, n)).toarray()
    K = csr_dense(a.mat, fes.ndof)[np.ix_(free, free)]
    M = csr_dense(m_bf.mat, fes.ndof)[np.ix_(free, free)]
    K = 0.5 * (K + K.T)
    M = 0.5 * (M + M.T)

    # Drop axis DOFs (M_diag=0 AND K_diag=0).
    diag_M = np.diag(M); diag_K = np.diag(K)
    keep = (diag_M > 1e-30 * (diag_M.max() + 1e-30)) | (diag_K > 1e-30 * (diag_K.max() + 1e-30))
    K = K[np.ix_(keep, keep)]
    M = M[np.ix_(keep, keep)]
    keep_M = np.diag(M) > 1e-30 * np.max(np.diag(M))
    K = K[np.ix_(keep_M, keep_M)]
    M = M[np.ix_(keep_M, keep_M)]
    print(f"  After drop axis & zero-mass: shape {K.shape}")

    t0 = time.time()
    eigvals = sla.eigvalsh(K, M, subset_by_index=[0, min(5, K.shape[0] - 1)])
    t_e = time.time() - t0
    pos = eigvals[eigvals > 0]
    tau_1 = 1.0 / pos[0] if len(pos) > 0 else np.nan
    print(f"  assemble K {t_a:.2f}s, M {t_m:.2f}s, eigh {t_e:.2f}s, "
          f"tau_1 = {tau_1*1e6:.4f} us")
    return tau_1


def main():
    print("=" * 64)
    print("P2 disk C++ assembly (no air, no Schur)")
    print("=" * 64)
    # Same levels as test_p2_disk_convergence so we can compare.
    levels = [(5, 2), (10, 4), (20, 6), (30, 8)]
    for N_r, N_z in levels:
        print(f"\n  (N_r, N_z) = ({N_r}, {N_z})")
        run_one(N_r, N_z)


if __name__ == "__main__":
    main()
