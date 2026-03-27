"""
Verification of ngsolve.bem LaplaceSL/LaplaceDL on a unit sphere.

Analytical solution for PEC sphere (radius a) in uniform field H_0 z_hat:
  B_n = 0 at surface (perfect diamagnet)
  J_pec = (3/2) * H_0 * sin(theta) * phi_hat
        = (3/2) * H_0 * (y, -x, 0)  on unit sphere

Tests:
  1. MFIE: (1/2 M + K) * J_pec == n x H_inc        [LaplaceDL]
     -> DL eigenvalue for l=1 mode: lambda_K = 1/6
  2. EFIE saddle point: PEC sphere (Z_s=0)           [LaplaceSL]
     -> J should recover (3/2) * H_0 * sin(theta)
  3. EFIE+SIBC: screening factor vs Z_s              [LaplaceSL]
     -> Z_s=0: J/J_pec -> 1;  Z_s >> 1: J/J_pec -> 0

Usage:
    python verify_laplace_bem.py
"""

import math
import os
import sys
import time
import numpy as np

MU_0 = 4e-7 * np.pi


def to_dense(mat, n):
    """Extract dense NumPy matrix from NGSolve BaseMatrix."""
    try:
        rows, cols, vals = mat.COO()
        M = np.zeros((n, n), dtype=float)
        for r, c, v in zip(rows, cols, vals):
            M[int(r), int(c)] = float(v)
        return M
    except Exception:
        pass
    D = mat.ToDense()
    return D.NumPy().copy()


def to_dense_rect(mat, nrow, ncol):
    """Extract rectangular dense matrix."""
    try:
        rows, cols, vals = mat.COO()
        M = np.zeros((nrow, ncol), dtype=float)
        for r, c, v in zip(rows, cols, vals):
            M[int(r), int(c)] = float(v)
        return M
    except Exception:
        pass
    D = mat.ToDense()
    return D.NumPy()[:nrow, :ncol].copy()


def make_sphere_mesh(R=1.0, maxh=0.3, order=2):
    """Create sphere surface mesh."""
    from ngsolve import Mesh
    from netgen.occ import Sphere, Pnt, OCCGeometry, Glue
    sph = Sphere(Pnt(0, 0, 0), R)
    for f in sph.faces:
        f.name = "sphere"
    geo = OCCGeometry(Glue(sph.faces))
    ngmesh = geo.GenerateMesh(maxh=maxh)
    mesh = Mesh(ngmesh)
    mesh.Curve(order)
    return mesh


def project_cf_to_hdiv(cf, fes, M_inv):
    """Project a CoefficientFunction onto HDivSurface via M^{-1} * lf."""
    from ngsolve import LinearForm, ds
    _, v = fes.TnT()
    lf = LinearForm(fes)
    lf += cf * v.Trace() * ds
    lf.Assemble()
    rhs = lf.vec.FV().NumPy().copy()
    return np.linalg.solve(M_inv, rhs)


# ============================================================
# Test 1: MFIE on PEC sphere - LaplaceDL eigenvalue
# ============================================================
def test_mfie_pec(mesh, H_0=1.0, verbose=True):
    """Verify (1/2 M + K) * J_pec == n x H_inc on PEC sphere.

    J_pec = (3/2)*H_0*(y, -x, 0),  n x H_inc = H_0*(y, -x, 0).
    Eigenvalue relation: (1/2 + lambda_K) = 2/3.
    """
    from ngsolve import (HDivSurface, BilinearForm, CF, ds, BND, TaskManager, x, y, z)
    from ngsolve.bem import LaplaceDL

    if verbose:
        print("=" * 65)
        print("Test 1: MFIE PEC sphere - LaplaceDL")
        print("=" * 65)

    fes = HDivSurface(mesh, order=0)
    u, v = fes.TnT()
    ndof = fes.ndof

    # Assemble DL and Mass
    t0 = time.perf_counter()
    with TaskManager():
        K_bf = LaplaceDL(u.Trace() * ds, use_fmm=False) * v.Trace() * ds
    K_mat = to_dense(K_bf.mat, ndof)
    t_K = time.perf_counter() - t0

    mass_bf = BilinearForm(fes)
    mass_bf += u.Trace() * v.Trace() * ds
    mass_bf.Assemble()
    M_mat = to_dense(mass_bf.mat, ndof)

    if verbose:
        print(f"  DOFs: {ndof}, DL assembly: {t_K:.1f}s")

    # Project J_pec = (3/2)*H_0*(y, -x, 0) and nxH = H_0*(y, -x, 0)
    J_pec_vec = project_cf_to_hdiv(CF((y, -x, 0)) * 1.5 * H_0, fes, M_mat)
    nxH_vec = project_cf_to_hdiv(CF((y, -x, 0)) * H_0, fes, M_mat)

    # Check that projections are non-zero
    if np.linalg.norm(J_pec_vec) < 1e-12:
        print("  ERROR: J_pec projection is zero!")
        return {'err_mfie': float('inf'), 'lambda_K': float('nan')}

    # LHS: (1/2 M + K) * J_pec
    lhs = (0.5 * M_mat + K_mat) @ J_pec_vec

    # RHS: M * nxH (since nxH_vec = M^{-1} * rhs_lf, we need M * nxH_vec)
    rhs = M_mat @ nxH_vec

    # Relative error
    err = np.linalg.norm(lhs - rhs) / np.linalg.norm(rhs)

    # Eigenvalue: (1/2 M + K)*J = alpha*M*J
    MJ = M_mat @ J_pec_vec
    alpha = (J_pec_vec @ lhs) / (J_pec_vec @ MJ)
    lambda_K = alpha - 0.5

    if verbose:
        print(f"  |(1/2 M + K)*J - M*nxH| / |M*nxH| = {err:.4e}")
        print(f"  Effective (1/2 + lambda_K) = {alpha:.6f} (expected 2/3 = {2/3:.6f})")
        print(f"  lambda_K = {lambda_K:.6f} (analytical 1/6 = {1/6:.6f})")
        print(f"  lambda_K relative error = {abs(lambda_K - 1/6) / (1/6):.2e}")

    return {
        'err_mfie': float(err),
        'lambda_K': float(lambda_K),
    }


# ============================================================
# Test 2: EFIE saddle point on PEC sphere - LaplaceSL
# ============================================================
def test_efie_pec(mesh, H_0=1.0, verbose=True):
    """EFIE saddle point on PEC sphere (Z_s=0).

    System: [mu0*SL, D^T; D, 0] * [J; p] = [-A_inc; 0]
    Expected: J = (3/2)*H_0*(y, -x, 0) = J_pec
    """
    from ngsolve import (HDivSurface, SurfaceL2, BilinearForm, LinearForm,
                         CF, ds, BND, TaskManager, div, x, y, z)
    from ngsolve.bem import LaplaceSL

    if verbose:
        print()
        print("=" * 65)
        print("Test 2: EFIE PEC sphere (Z_s=0) - LaplaceSL")
        print("=" * 65)

    fes_J = HDivSurface(mesh, order=0)
    fes_L2 = SurfaceL2(mesh, order=0)
    u, v = fes_J.TnT()
    ndof = fes_J.ndof
    n_p = fes_L2.ndof

    # SL matrix
    t0 = time.perf_counter()
    with TaskManager():
        SL_bf = LaplaceSL(u.Trace() * ds, use_fmm=False) * v.Trace() * ds
    SL_mat = to_dense(SL_bf.mat, ndof)
    t_SL = time.perf_counter() - t0

    # Mass matrix
    mass_bf = BilinearForm(fes_J)
    mass_bf += u.Trace() * v.Trace() * ds
    mass_bf.Assemble()
    M_mat = to_dense(mass_bf.mat, ndof)

    # Divergence matrix D
    u_J = fes_J.TrialFunction()
    q = fes_L2.TestFunction()
    bf_D = BilinearForm(trialspace=fes_J, testspace=fes_L2)
    bf_D += div(u_J.Trace()) * q * ds
    bf_D.Assemble()
    D_mat = to_dense_rect(bf_D.mat, n_p, ndof)
    rank_D = np.linalg.matrix_rank(D_mat)
    D_red = D_mat[:rank_D, :]

    if verbose:
        print(f"  DOFs: {ndof} (J) + {n_p} (p), SL assembly: {t_SL:.1f}s")
        print(f"  D rank: {rank_D}/{n_p}")

    # RHS: A_inc for uniform B = mu0*H0*z_hat -> A = (mu0*H0/2)*(-y, x, 0)
    lf_A = LinearForm(fes_J)
    lf_A += CF((-y, x, 0)) * (MU_0 * H_0 / 2) * v.Trace() * ds
    lf_A.Assemble()
    f_A = lf_A.vec.FV().NumPy().copy()

    # For PEC: mu0 * SL * J + D^T * p = -A_inc_proj (jw cancels)
    n_total = ndof + rank_D
    A_block = np.zeros((n_total, n_total), dtype=float)
    A_block[:ndof, :ndof] = MU_0 * SL_mat
    A_block[:ndof, ndof:] = D_red.T
    A_block[ndof:, :ndof] = D_red

    rhs_block = np.zeros(n_total, dtype=float)
    rhs_block[:ndof] = -f_A

    J_p = np.linalg.solve(A_block, rhs_block)
    J_vec = J_p[:ndof]

    # Reference: J_pec = M^{-1} * lf(J_pec_cf)
    J_pec_vec = project_cf_to_hdiv(CF((y, -x, 0)) * 1.5 * H_0, fes_J, M_mat)

    # Ratio (least squares)
    ratio = np.dot(J_vec, M_mat @ J_pec_vec) / np.dot(J_pec_vec, M_mat @ J_pec_vec)
    err_shape = np.linalg.norm(J_vec - ratio * J_pec_vec) / np.linalg.norm(ratio * J_pec_vec)

    # Energy norms
    J_M = np.sqrt(J_vec @ M_mat @ J_vec)
    Jpec_M = np.sqrt(J_pec_vec @ M_mat @ J_pec_vec)

    if verbose:
        print(f"  J/J_pec ratio = {ratio:.6f} (expected ~1.0)")
        print(f"  Shape error = {err_shape:.4e}")
        print(f"  ||J||_M = {J_M:.4f}, ||J_pec||_M = {Jpec_M:.4f}")
        print(f"  div(J) residual = {np.linalg.norm(D_mat @ J_vec):.2e}")

    return {
        'ratio': float(ratio),
        'shape_err': float(err_shape),
    }


# ============================================================
# Test 3: EFIE+SIBC screening factor sweep
# ============================================================
def test_screening_sweep(mesh, H_0=1.0, verbose=True):
    """Sweep Z_s: Z_s=0 -> PEC (ratio~1), Z_s >> 1 -> transparent (ratio~0)."""
    from ngsolve import (HDivSurface, SurfaceL2, BilinearForm, LinearForm,
                         CF, ds, BND, TaskManager, div, x, y, z)
    from ngsolve.bem import LaplaceSL

    if verbose:
        print()
        print("=" * 65)
        print("Test 3: Screening factor vs Z_s sweep (EFIE+SIBC)")
        print("=" * 65)

    omega = 2 * np.pi * 1000

    fes_J = HDivSurface(mesh, order=0)
    fes_L2 = SurfaceL2(mesh, order=0)
    u, v = fes_J.TnT()
    ndof = fes_J.ndof
    n_p = fes_L2.ndof

    with TaskManager():
        SL_bf = LaplaceSL(u.Trace() * ds, use_fmm=False) * v.Trace() * ds
    SL_mat = to_dense(SL_bf.mat, ndof)

    mass_bf = BilinearForm(fes_J)
    mass_bf += u.Trace() * v.Trace() * ds
    mass_bf.Assemble()
    M_mat = to_dense(mass_bf.mat, ndof)

    u_J = fes_J.TrialFunction()
    q = fes_L2.TestFunction()
    bf_D = BilinearForm(trialspace=fes_J, testspace=fes_L2)
    bf_D += div(u_J.Trace()) * q * ds
    bf_D.Assemble()
    D_mat = to_dense_rect(bf_D.mat, n_p, ndof)
    rank_D = np.linalg.matrix_rank(D_mat)
    D_red = D_mat[:rank_D, :]

    # RHS: -jw * A_inc
    lf_A = LinearForm(fes_J)
    lf_A += CF((-y, x, 0)) * (MU_0 * H_0 / 2) * v.Trace() * ds
    lf_A.Assemble()
    f_A = -1j * omega * lf_A.vec.FV().NumPy().copy().astype(complex)

    # PEC reference
    J_pec_vec = project_cf_to_hdiv(CF((y, -x, 0)) * 1.5 * H_0, fes_J, M_mat)
    Jpec_M2 = J_pec_vec @ M_mat @ J_pec_vec

    # Sweep
    Z_s_values = np.logspace(-7, 1, 17)
    n_total = ndof + rank_D

    if verbose:
        print(f"  {'Z_s':>10s} {'|J/J_pec|':>10s} {'P [W]':>12s}")

    results = []
    for Z_s_val in Z_s_values:
        Z_s = complex(Z_s_val)
        A_block = np.zeros((n_total, n_total), dtype=complex)
        A_block[:ndof, :ndof] = Z_s * M_mat + 1j * omega * MU_0 * SL_mat
        A_block[:ndof, ndof:] = D_red.T
        A_block[ndof:, :ndof] = D_red
        rhs_block = np.zeros(n_total, dtype=complex)
        rhs_block[:ndof] = f_A

        J_p = np.linalg.solve(A_block, rhs_block)
        J_vec = J_p[:ndof]

        ratio = abs(np.dot(J_vec, M_mat @ J_pec_vec) / Jpec_M2)
        P = 0.5 * Z_s_val * float(abs(np.conj(J_vec) @ M_mat @ J_vec))

        results.append({'Z_s': Z_s_val, 'ratio': ratio, 'P': P})
        if verbose:
            print(f"  {Z_s_val:10.2e} {ratio:10.6f} {P:12.4e}")

    return results


# ============================================================
# Main
# ============================================================
if __name__ == "__main__":
    from ngsolve import BND

    print("Creating unit sphere mesh...")
    mesh = make_sphere_mesh(R=1.0, maxh=0.25, order=2)
    ne = mesh.GetNE(BND)
    print(f"  {ne} surface elements, {mesh.nv} vertices")
    print()

    r1 = test_mfie_pec(mesh)
    r2 = test_efie_pec(mesh)
    r3 = test_screening_sweep(mesh)

    # Summary
    print()
    print("=" * 65)
    print("Summary")
    print("=" * 65)
    # NOTE: LaplaceDL for HDivSurface is the SCALAR double layer (dG/dn_y),
    # NOT the MFIE K operator (n x curl SL). Eigenvalue -1/6 is correct for
    # scalar DL on l=1 mode. This confirms LaplaceDL is NOT suitable for MFIE.
    print(f"  Test 1 (scalar DL): lambda_K = {r1['lambda_K']:.6f} "
          f"(scalar DL exact: -1/6 = {-1/6:.6f})")
    print(f"  Test 2 (EFIE PEC): J/J_pec = {r2['ratio']:.6f}, "
          f"shape err = {r2['shape_err']:.2e}")
    print(f"  Test 3: Z_s=1e-7 -> {r3[0]['ratio']:.4f}, "
          f"Z_s=10 -> {r3[-1]['ratio']:.6f}")

    ok_1 = abs(r1['lambda_K'] - (-1.0/6)) / (1.0/6) < 0.01  # scalar DL: -1/6
    ok_2 = abs(r2['ratio'] - 1.0) < 0.15
    ok_3 = r3[0]['ratio'] > 0.85 and r3[-1]['ratio'] < 0.01
    all_pass = ok_1 and ok_2 and ok_3
    print(f"\n  Test 1 (DL = -1/6): {'PASS' if ok_1 else 'FAIL'}")
    print(f"  Test 2 (EFIE PEC): {'PASS' if ok_2 else 'FAIL'}")
    print(f"  Test 3 (screening): {'PASS' if ok_3 else 'FAIL'}")
    print(f"  Overall: {'PASS' if all_pass else 'FAIL'}")
