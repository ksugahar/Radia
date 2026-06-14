"""Diagnose where NaN comes from in test_p2_cpp_sphere assembly.

Replicate the (4, 4, 12) mesh, assemble K and M, find rows/cols containing
NaN, and report which DOFs and (if possible) which triangle is responsible.
"""
from __future__ import annotations
import sys
sys.stdout.reconfigure(encoding="utf-8")

from math import pi
import numpy as np
import scipy.sparse as sp

import netgen.meshing as ng_meshing
from netgen.meshing import (
    Mesh as NgMesh, Element1D, Element2D, FaceDescriptor, MeshPoint, Pnt,
)
from ngsolve import (
    Mesh, BilinearForm, CoefficientFunction, TaskManager, ngsglobals,
)
from radia.radia_axifemm import (
    H1Henrotte, AxiHenrotteStiffnessBFI, AxiHenrotteSigmaMassBFI,
)

R_SPHERE = 10.0e-3
R_AIR = 50.0e-3
SIGMA_CU = 5.8e7
MU0 = 4 * pi * 1e-7


def make_sphere_axisym_trig_mesh(N_rho_in, N_rho_out, N_theta):
    rho_in = np.linspace(0.0, R_SPHERE, N_rho_in + 1)
    rho_out = np.linspace(R_SPHERE, R_AIR, N_rho_out + 1)
    rho_grid = np.concatenate([rho_in, rho_out[1:]])
    theta_grid = np.linspace(0.0, pi, N_theta + 1)
    n_rho = len(rho_grid)

    mesh = NgMesh()
    mesh.dim = 2
    mesh.SetMaterial(1, "air")
    mesh.SetMaterial(2, "conductor")
    mesh.Add(FaceDescriptor(surfnr=1, domin=0, bc=1))
    mesh.Add(FaceDescriptor(surfnr=2, domin=0, bc=2))
    mesh.SetBCName(0, "outer")
    mesh.SetBCName(1, "inner")

    pnums = []
    coords_all = []
    for theta in theta_grid:
        for rho in rho_grid:
            r = rho * np.sin(theta)
            z = rho * np.cos(theta)
            if abs(r) < 1e-12: r = 0.0
            pnums.append(mesh.Add(MeshPoint(Pnt(r, z, 0))))
            coords_all.append((r, z))

    def coords_of(pnum_idx):
        j = pnum_idx // n_rho
        i = pnum_idx % n_rho
        theta = theta_grid[j]
        rho = rho_grid[i]
        return rho * np.sin(theta), rho * np.cos(theta)

    def triangle_has_off_axis(node_indices, atol=1e-12):
        for idx in node_indices:
            r, _ = coords_of(idx)
            if abs(r) > atol:
                return True
        return False

    tri_records = []
    for j in range(N_theta):
        for i in range(n_rho - 1):
            i00 = j * n_rho + i
            i10 = j * n_rho + i + 1
            i01 = (j + 1) * n_rho + i
            i11 = (j + 1) * n_rho + i + 1
            n00, n10, n01, n11 = pnums[i00], pnums[i10], pnums[i01], pnums[i11]
            rho_mid = 0.5 * (rho_grid[i] + rho_grid[i + 1])
            mat_idx = 2 if rho_mid < R_SPHERE else 1
            if triangle_has_off_axis([i00, i10, i11]):
                tri_records.append((mat_idx, (i00, i10, i11)))
                mesh.Add(Element2D(mat_idx, [n00, n10, n11]))
            if triangle_has_off_axis([i00, i11, i01]):
                tri_records.append((mat_idx, (i00, i11, i01)))
                mesh.Add(Element2D(mat_idx, [n00, n11, n01]))

    i_outer = n_rho - 1
    for j in range(N_theta):
        p0 = pnums[j * n_rho + i_outer]
        p1 = pnums[(j + 1) * n_rho + i_outer]
        mesh.Add(Element1D([p0, p1], index=1))

    return Mesh(mesh), tri_records, coords_all


def to_csr(mat, n):
    rows, cols, vals = mat.COO()
    return sp.coo_matrix(
        (np.array(vals), (np.array(rows), np.array(cols))),
        shape=(n, n),
    ).tocsr()


def main():
    ngsglobals.msg_level = 0
    mesh, tri_records, coords_all = make_sphere_axisym_trig_mesh(4, 4, 12)
    fes = H1Henrotte(mesh, order=2, dirichlet="outer")
    print(f"ne={mesh.ne} nv={mesh.nv} ndof={fes.ndof}")

    a = BilinearForm(fes, symmetric=True, check_unused=False)
    a += AxiHenrotteStiffnessBFI(CoefficientFunction(MU0))
    with TaskManager(): a.Assemble()

    sigma_cf = mesh.MaterialCF({"conductor": SIGMA_CU}, default=0.0)
    m_bf = BilinearForm(fes, symmetric=True, check_unused=False)
    m_bf += AxiHenrotteSigmaMassBFI(sigma_cf)
    with TaskManager(): m_bf.Assemble()

    K = to_csr(a.mat, fes.ndof).toarray()
    M = to_csr(m_bf.mat, fes.ndof).toarray()

    # Find NaN entries
    nan_K = ~np.isfinite(K)
    nan_M = ~np.isfinite(M)
    print(f"K: total {K.size}, NaN count {nan_K.sum()}, finite max {np.nanmax(np.abs(K)):.3e}")
    print(f"M: total {M.size}, NaN count {nan_M.sum()}, finite max {np.nanmax(np.abs(M)):.3e}")

    nan_rows_K = np.where(nan_K.any(axis=1))[0]
    nan_rows_M = np.where(nan_M.any(axis=1))[0]
    print(f"K rows with NaN: {len(nan_rows_K)} -> {nan_rows_K[:10]}{'...' if len(nan_rows_K)>10 else ''}")
    print(f"M rows with NaN: {len(nan_rows_M)} -> {nan_rows_M[:10]}{'...' if len(nan_rows_M)>10 else ''}")

    # NV = number of vertices = mesh.nv. Edge DOFs start at NV.
    nv = mesh.nv
    print(f"\n  nv={nv}  edge DOFs start at index {nv}")
    print(f"  Are NaN rows vertex DOFs ({nan_rows_K[nan_rows_K < nv]}) or edge DOFs ({nan_rows_K[nan_rows_K >= nv] - nv})?")

    # Check diagonal entries explicitly
    diag_K = np.diag(K)
    diag_M = np.diag(M)
    print(f"\n  K diagonal: {np.isnan(diag_K).sum()} NaN, "
          f"finite max {np.nanmax(np.abs(diag_K)):.3e}")
    print(f"  M diagonal: {np.isnan(diag_M).sum()} NaN, "
          f"finite max {np.nanmax(np.abs(diag_M)):.3e}")

    # Are the NaN DOFs free or Dirichlet?
    free = np.array([fes.FreeDofs()[i] for i in range(fes.ndof)], dtype=bool)
    print(f"\n  NaN-row DOFs free? K: {free[nan_rows_K].sum()}/{len(nan_rows_K)} free")
    print(f"  NaN-row DOFs free? M: {free[nan_rows_M].sum()}/{len(nan_rows_M)} free")


if __name__ == "__main__":
    main()
