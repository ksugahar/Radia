"""test_p2_cpp_sphere.py — Phase B2 sphere benchmark via NGSolve.

Build a netgen triangular axisymmetric mesh of Cu sphere + 50 mm air box,
use H1Henrotte order=2 FESpace which dispatches AxiHenrotteFE_P2_Triangle
for triangles (C++ Phase B2). Assemble K and sigma-mass via NGSolve, solve
generalized eigenvalue, compare with Stoll 738.48 us.

Reference: Python Phase B1c gave -0.51% off Stoll. The C++ Phase B2
implementation must match the Python reference at machine precision (both
use FEMM Henrotte (r_i/r) basis + Duffy-Gauss 8x8 quadrature).
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
    Mesh as NgMesh, Element2D, Element1D, FaceDescriptor, MeshPoint, Pnt,
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

STOLL_TAU_1 = MU0 * SIGMA_CU * R_SPHERE ** 2 / (pi ** 2)


def make_sphere_axisym_trig_mesh(N_rho_in, N_rho_out, N_theta):
    """Build axisymmetric triangulated mesh in (r, z) plane for sphere + air.

    Structured (rho, theta) -> (r, z) = (rho sin theta, rho cos theta).
    Conductor: rho in [0, R_SPHERE]. Air: rho in [R_SPHERE, R_AIR].
    Each (rho, theta) cell split into 2 triangles.
    """
    rho_in = np.linspace(0.0, R_SPHERE, N_rho_in + 1)
    rho_out = np.linspace(R_SPHERE, R_AIR, N_rho_out + 1)
    rho_grid = np.concatenate([rho_in, rho_out[1:]])
    theta_grid = np.linspace(0.0, pi, N_theta + 1)
    n_rho = len(rho_grid)

    # netgen 2D mesh, matching test_hiruma_disk_q1.py pattern:
    # - SetMaterial(N, name) names volume domain N.
    # - One FaceDescriptor must exist for every Element2D index in use;
    #   otherwise Mesh(ngmesh) crashes ("has no facedecoding: ind = K").
    # - Each FaceDescriptor must have a DISTINCT `bc` value, otherwise
    #   GetMaterials() collapses and only reports one domain.
    mesh = NgMesh()
    mesh.dim = 2
    mesh.SetMaterial(1, "air")
    mesh.SetMaterial(2, "conductor")
    mesh.Add(FaceDescriptor(surfnr=1, domin=0, bc=1))   # idx 1 -> air
    mesh.Add(FaceDescriptor(surfnr=2, domin=0, bc=2))   # idx 2 -> conductor
    mesh.SetBCName(0, "outer")
    mesh.SetBCName(1, "inner")  # interior interface -- unused but distinct

    # Add points, deduplicating axis points so each physical (r=0, z) becomes
    # a single mesh vertex. Without this, NGSolve creates degenerate
    # zero-length edges between coincident axis vertices (each theta slice
    # gives a separate vertex at the same (0, z) for rho=0; also theta=0/pi
    # rows have r=rho*sin(0)=0 with z = rho*cos(0)=rho lying on the axis).
    # Both cause P2 Vandermonde singularity → NaN element matrices.
    pnums = [None] * (n_rho * (N_theta + 1))
    axis_z_to_pnum = {}  # (z rounded to 1e-15) -> pnum, for r=0 vertices

    def add_or_get(r, z):
        if abs(r) < 1e-12:
            key = round(z, 15)
            if key in axis_z_to_pnum:
                return axis_z_to_pnum[key]
            p = mesh.Add(MeshPoint(Pnt(0.0, z, 0)))
            axis_z_to_pnum[key] = p
            return p
        return mesh.Add(MeshPoint(Pnt(r, z, 0)))

    for j, theta in enumerate(theta_grid):
        for i, rho in enumerate(rho_grid):
            r = rho * np.sin(theta)
            z = rho * np.cos(theta)
            pnums[j * n_rho + i] = add_or_get(r, z)

    axis_pnums = set(axis_z_to_pnum.values())

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

    def triangle_has_3_distinct_pnums(triple):
        return len(set(triple)) == 3

    for j in range(N_theta):
        for i in range(n_rho - 1):
            i00 = j * n_rho + i
            i10 = j * n_rho + i + 1
            i01 = (j + 1) * n_rho + i
            i11 = (j + 1) * n_rho + i + 1
            n00, n10, n01, n11 = pnums[i00], pnums[i10], pnums[i01], pnums[i11]
            rho_mid = 0.5 * (rho_grid[i] + rho_grid[i + 1])
            mat_idx = 2 if rho_mid < R_SPHERE else 1
            # Skip both: all-axis (zero volume) AND collapsed (two pnums equal
            # after dedup, e.g. (rho=0) axis vertex shared across the 'i00 vs
            # i01' direction).
            if (triangle_has_off_axis([i00, i10, i11])
                    and triangle_has_3_distinct_pnums((n00, n10, n11))):
                mesh.Add(Element2D(mat_idx, [n00, n10, n11]))
            if (triangle_has_off_axis([i00, i11, i01])
                    and triangle_has_3_distinct_pnums((n00, n11, n01))):
                mesh.Add(Element2D(mat_idx, [n00, n11, n01]))

    i_outer = n_rho - 1
    for j in range(N_theta):
        p0 = pnums[j * n_rho + i_outer]
        p1 = pnums[(j + 1) * n_rho + i_outer]
        mesh.Add(Element1D([p0, p1], index=1))

    return Mesh(mesh), axis_pnums


def to_scipy_csr(ngs_mat, ndof):
    """Convert NGSolve BilinearForm matrix to scipy.sparse.csr_matrix."""
    rows, cols, vals = ngs_mat.COO()
    return sp.coo_matrix(
        (np.array(vals), (np.array(rows), np.array(cols))),
        shape=(ndof, ndof),
    ).tocsr()


def run_one(N_rho_in, N_rho_out, N_theta):
    ngsglobals.msg_level = 0
    mesh, axis_pnums = make_sphere_axisym_trig_mesh(N_rho_in, N_rho_out, N_theta)
    print(f"  Materials: {mesh.GetMaterials()}")
    print(f"  Boundaries: {mesh.GetBoundaries()}")
    fes = H1Henrotte(mesh, order=2, dirichlet="outer")
    n_free = sum(1 for f in fes.FreeDofs() if f)
    print(f"  Mesh ne={mesh.ne} ndof={fes.ndof} free={n_free}  "
          f"axis_pnums={len(axis_pnums)}")

    mu_cf = CoefficientFunction(MU0)
    sigma_cf = mesh.MaterialCF({"conductor": SIGMA_CU}, default=0.0)

    a = BilinearForm(fes, symmetric=True, check_unused=False)
    a += AxiHenrotteStiffnessBFI(mu_cf)
    t0 = time.time()
    with TaskManager(): a.Assemble()
    t_a = time.time() - t0

    m = BilinearForm(fes, symmetric=True, check_unused=False)
    m += AxiHenrotteSigmaMassBFI(sigma_cf)
    t0 = time.time()
    with TaskManager(): m.Assemble()
    t_m = time.time() - t0

    K_csr = to_scipy_csr(a.mat, fes.ndof)
    M_csr = to_scipy_csr(m.mat, fes.ndof)
    free = np.array([i for i in range(fes.ndof) if fes.FreeDofs()[i]],
                    dtype=int)
    K_red = K_csr[free[:, None], free[None, :]].toarray()
    M_red = M_csr[free[:, None], free[None, :]].toarray()
    K_red = 0.5 * (K_red + K_red.T)
    M_red = 0.5 * (M_red + M_red.T)
    print(f"  K.max={np.abs(K_red).max():.2e}  M.max={np.abs(M_red).max():.2e}")

    # Schur reduction: eliminate air DOFs from K-side, retain conductor block.
    # cond_idx = DOFs with nonzero M diagonal (conductor mass-mat presence).
    # air_idx  = remaining free DOFs.
    # S = K_cc - K_ca K_aa^{-1} K_ac;  solve S u = lambda M_cc u.
    diag_M = np.diag(M_red)
    diag_K = np.diag(K_red)
    cond_mask = diag_M > 1e-30 * (np.max(diag_M) + 1e-30)
    # Air = free \ conductor \ axis-zero-K (auto-decoupled by r_i factor).
    nonzero_K = diag_K > 1e-30 * (np.max(diag_K) + 1e-30)
    air_mask = (~cond_mask) & nonzero_K
    cond_idx = np.where(cond_mask)[0]
    air_idx = np.where(air_mask)[0]
    print(f"  Schur split: cond={len(cond_idx)}  air={len(air_idx)}  "
          f"axis-zero-dropped={(~nonzero_K).sum()}")

    t0 = time.time()
    if len(air_idx) > 0:
        K_cc = K_red[np.ix_(cond_idx, cond_idx)]
        K_ca = K_red[np.ix_(cond_idx, air_idx)]
        K_aa = K_red[np.ix_(air_idx, air_idx)]
        X = sla.solve(K_aa, K_ca.T, assume_a="sym")
        S = K_cc - K_ca @ X
    else:
        S = K_red[np.ix_(cond_idx, cond_idx)]
    M_cc = M_red[np.ix_(cond_idx, cond_idx)]
    S = 0.5 * (S + S.T)
    M_cc = 0.5 * (M_cc + M_cc.T)
    t_schur = time.time() - t0

    t0 = time.time()
    eigvals = sla.eigvalsh(S, M_cc,
                            subset_by_index=[0, min(5, S.shape[0] - 1)])
    t_eigh = time.time() - t0
    pos = eigvals[eigvals > 0]
    tau_1 = 1.0 / pos[0] if len(pos) > 0 else np.nan
    gap = (tau_1 - STOLL_TAU_1) / STOLL_TAU_1 * 100
    print(f"  assemble K {t_a:.2f}s, M {t_m:.2f}s, schur {t_schur:.2f}s, "
          f"eigh {t_eigh:.2f}s")
    return tau_1, gap


def main():
    print("=" * 64)
    print("Phase B2 sphere benchmark via NGSolve assembly (C++ P2)")
    print("=" * 64)
    print(f"\nStoll tau_1 = {STOLL_TAU_1*1e6:.4f} us  "
          f"(Cu sphere R={R_SPHERE*1e3} mm)")
    print(f"Air box R_air = {R_AIR*1e3} mm (Phase B3 Kelvin removes truncation)\n")

    levels = [
        (4, 4, 12),
        (6, 6, 16),
        (8, 8, 24),
        (10, 10, 30),
    ]
    print(f"  {'(in, out, th)':>15s}  {'tau [us]':>12s}  {'gap':>10s}")
    for N_in, N_out, N_th in levels:
        print(f"\n  ({N_in}, {N_out}, {N_th}):")
        tau, gap = run_one(N_in, N_out, N_th)
        print(f"    tau_1 = {tau*1e6:.4f} us  gap = {gap:+.3f}%")


if __name__ == "__main__":
    main()
