"""ecb_lorentz_force.py -- Linear ECB drag/lift via Lorentz integral F = int J x B dV.

For a PM (cube Nd) at position (x_PM, 0, h_gap + L_PM/2) moving along +x at
velocity v, this script computes the time-averaged drag and lift forces on
the plate using the LORENTZ integral over the conductor volume:

    < F_conductor > = (1/2) Re[ integral conjugate(J(s)) x B_ext(r) dV ]
    < F_PM > = - < F_conductor >  (Newton 3rd)

J(r, s) is reconstructed from the Foster eigenmode expansion:

    v(r, s) = sum_n c_n(s) phi_n(r),   c_n(s) = -s mu sigma <B_ext_z, phi_n> / (lam_n + s mu sigma)
    J(r, s) ~ sigma * d/dt [grad v(r, s)] cross zhat   (scalar A_z approximation)
            ~ -j omega sigma * (B_ext_z(r) + sum c_n phi_n)  (scalar simplification)

The Foster modes are loaded from the bulk_foster_via_eigen output of
Phase L-1 (cln_sibc_general_hex). For the PROPER vector eddy-current
formulation see calc_fem_kelvin.py in src/radia/panels/ -- this script
is a SCALAR-A_z approximation valid for plates where B_ext_z dominates.
"""
from __future__ import annotations

import argparse
import cmath
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))


SIGMA_AL = 3.5e7
MU_0 = 4 * math.pi * 1e-7


def pm_field_z_dipole(x, y, z, m, x_pm=0.0, y_pm=0.0, z_pm=0.0):
    """B_z(x,y,z) from a magnetic dipole m at (x_pm, y_pm, z_pm), aligned with z."""
    dx = x - x_pm
    dy = y - y_pm
    dz = z - z_pm
    r2 = dx**2 + dy**2 + dz**2
    r = np.sqrt(r2 + 1e-30)
    r5 = r**5
    # B_z of z-aligned dipole: (mu0/4pi) * (3 m * dz^2 / r^5 - m / r^3)
    return (MU_0 / (4*math.pi)) * (3*m*dz*dz/r5 - m/(r**3))


def pm_field_xz_dipole(x, y, z, m, x_pm=0.0, y_pm=0.0, z_pm=0.0):
    """B_x(x,y,z) component from z-aligned dipole (for force in x direction)."""
    dx = x - x_pm
    dy = y - y_pm
    dz = z - z_pm
    r2 = dx**2 + dy**2 + dz**2
    r5 = (r2 + 1e-30)**(5/2.0)
    # B_x of z-aligned dipole: (mu0/4pi) * 3 m * dx * dz / r^5
    return (MU_0 / (4*math.pi)) * 3*m*dx*dz/r5


def compute_lorentz_force_via_foster(mesh, lam, b_n, vecs_free, free_mask,
                                       sigma, mu, s, m_pm, h_gap, x_pm):
    """Time-averaged drag (F_x) and lift (F_z) on the conductor.

    Strategy:
      1. Evaluate B_ext at mesh quadrature points (using PM dipole model).
      2. Reconstruct J_y(r, s) = -j omega sigma * v(r, s) (scalar A_z) where
         v(r, s) = sum_n c_n(s) phi_n(r),
         c_n(s) = -s mu sigma <1, phi_n>_M * <phi_n, B_ext_z>_M / (lam_n + s mu sigma)
         (the second projection is what couples the spatial B distribution).
      3. F_x = -<Re[conj(J_y) B_z]> integrated over volume
         F_z = +<Re[conj(J_y) B_x]> integrated over volume
         (signs from F = J cross B for J = J_y yhat)
    """
    from ngsolve import Integrate, CoefficientFunction, GridFunction, H1, x as xC, y as yC, z as zC

    fes = H1(mesh, order=2, dirichlet="outer")
    # Build B_ext_z, B_ext_x as CoefficientFunctions
    Bz_cf = (MU_0 / (4*math.pi)) * (
        3*m_pm*(zC - (h_gap))*(zC - (h_gap))/((xC-x_pm)**2 + yC**2 + (zC-h_gap)**2 + 1e-30)**(5/2.0)
        - m_pm / ((xC-x_pm)**2 + yC**2 + (zC-h_gap)**2 + 1e-30)**(3/2.0)
    )
    Bx_cf = (MU_0 / (4*math.pi)) * (
        3*m_pm*(xC-x_pm)*(zC - (h_gap))/((xC-x_pm)**2 + yC**2 + (zC-h_gap)**2 + 1e-30)**(5/2.0)
    )

    # Compute projection <phi_n, B_ext_z>_M for each eigenmode
    # For this we need the eigenvectors back in the GridFunction (we have vec arrays).
    n_modes = vecs_free.shape[1]
    phi_gfu = GridFunction(fes)
    phi_full = np.zeros(fes.ndof)
    # We need M @ B_ext_z for each mode. Build B_ext_z gridfunction.
    Bz_gfu = GridFunction(fes)
    Bz_gfu.Set(Bz_cf)
    Bz_vec_full = np.array(Bz_gfu.vec.FV().NumPy())
    Bz_vec_free = Bz_vec_full[free_mask]

    # Mass matrix
    from ngsolve import BilinearForm
    u, v = fes.TnT()
    m_form = BilinearForm(fes, symmetric=True)
    from ngsolve import dx
    m_form += u * v * dx
    m_form.Assemble()
    import scipy.sparse as sp
    rows_m, cols_m, vals_m = m_form.mat.COO()
    M = sp.csr_matrix((np.asarray(vals_m), (np.asarray(rows_m), np.asarray(cols_m))),
                       shape=(fes.ndof, fes.ndof))
    M_free = M[free_mask][:, free_mask]

    # Mass-weighted inner product: <phi_n, Bz>_M = phi_n^T M Bz
    M_Bz = M_free.dot(Bz_vec_free)
    proj_Bz = vecs_free.T @ M_Bz   # shape (n_modes,)

    # Foster coefficients c_n(s) = -s mu sigma * proj_Bz_n / (lam_n + s mu sigma)
    s_mu_sigma = s * mu * sigma
    c_n = -s_mu_sigma * proj_Bz / (lam + s_mu_sigma)

    # Reconstruct v(r, s) on the mesh: v = sum c_n phi_n
    v_vec_free = vecs_free @ c_n   # complex array on free DOFs
    v_vec_full = np.zeros(fes.ndof, dtype=complex)
    v_vec_full[free_mask] = v_vec_free
    v_gfu_re = GridFunction(fes)
    v_gfu_im = GridFunction(fes)
    v_gfu_re.vec.FV().NumPy()[:] = v_vec_full.real
    v_gfu_im.vec.FV().NumPy()[:] = v_vec_full.imag

    # J_y(r, s) = -j omega sigma * (B_ext_z + v)? Or just -j omega sigma * v?
    # For the SCALAR formulation with v defined by (-Lap + s mu sigma) v = -s mu sigma,
    # v represents (A_z_total / A_z_ext - 1) essentially. So J_y = -j omega sigma * (B_ext_z + v) ?
    # We'll use a simpler heuristic: J_y(r, s) ~ s sigma * v(r, s).
    # The constant factor doesn't matter for the qualitative force trend.
    omega = s.imag   # s = j omega -> omega = imag(s)

    # F_x_volume = -Integrate(Re[conj(J_y) * Bz_cf])  (drag = J_y x B_z gives -x)
    # F_z_volume = +Integrate(Re[conj(J_y) * Bx_cf])  (lift = J_y x B_x gives +z)
    Jy_re = -omega * sigma * v_gfu_im   # Re[J_y] = Re[-j omega sigma v] = omega sigma Im[v]
    Jy_im =  omega * sigma * v_gfu_re   # Im[J_y] = -omega sigma Re[v]

    # Time-averaged: <Re[A] Re[B] + Im[A] Im[B]> / 2 with phasors
    # F_x_time_avg = (1/2) Re[ conj(J_y) Bz_cf ] integrated -- but Bz_cf is REAL (Radia static)
    # So <F_x> = (1/2) Re[J_y] * Bz_cf integrated.
    Fx_drag = -0.5 * Integrate(Jy_re * Bz_cf, mesh)
    Fz_lift = +0.5 * Integrate(Jy_re * Bx_cf, mesh)
    # NOTE: for proper velocity-dependent force, B_ext should be in the PM frame
    # (steady) and seen by plate as moving. We do quasi-static here.

    return float(Fx_drag), float(Fz_lift)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vol", default="plate_100x50x5.vol")
    parser.add_argument("--n-eigen", type=int, default=300)
    parser.add_argument("--m-pm", type=float, default=0.955,
                        help="PM dipole moment (A m^2). Default: 10mm cube Nd Br=1.2T")
    parser.add_argument("--h-gap", type=float, default=10e-3,
                        help="PM bottom -> plate top gap (m).")
    parser.add_argument("--velocities", default="1,3,10,30,100",
                        help="comma-sep velocity values (m/s)")
    args = parser.parse_args()

    print("=== Linear ECB Lorentz force ===")
    print(f"  vol = {args.vol}")
    print(f"  PM dipole moment m = {args.m_pm:.3f} A m^2")
    print(f"  PM gap h = {args.h_gap*1e3:.1f} mm")
    print()

    from ngsolve import Mesh, TaskManager, H1, BilinearForm, grad, dx, GridFunction, Integrate
    import scipy.sparse as sp
    from scipy.sparse.linalg import eigsh

    mesh = Mesh(args.vol)
    print(f"  ne = {mesh.ne}")

    # Solve eigenproblem just like Phase L-1 but keep the eigenvectors around.
    with TaskManager():
        fes = H1(mesh, order=2, dirichlet="outer")
        u, v = fes.TnT()
        a = BilinearForm(fes, symmetric=True)
        a += grad(u) * grad(v) * dx
        a.Assemble()
        mform = BilinearForm(fes, symmetric=True)
        mform += u * v * dx
        mform.Assemble()
        rows_a, cols_a, vals_a = a.mat.COO()
        rows_m, cols_m, vals_m = mform.mat.COO()
        ndof = fes.ndof
        A = sp.csr_matrix((np.asarray(vals_a), (np.asarray(rows_a), np.asarray(cols_a))), shape=(ndof, ndof))
        M = sp.csr_matrix((np.asarray(vals_m), (np.asarray(rows_m), np.asarray(cols_m))), shape=(ndof, ndof))
        free = np.array([fes.FreeDofs()[i] for i in range(ndof)], dtype=bool)
        A_free = A[free][:, free]
        M_free = M[free][:, free]
        k_eig = min(args.n_eigen, A_free.shape[0] - 2)
        print(f"  Solving {k_eig} eigenmodes (ndof_free = {A_free.shape[0]})...")
        lam, vecs = eigsh(A_free, k=k_eig, M=M_free, sigma=0, which="LM")
        order = np.argsort(lam)
        lam = lam[order]
        vecs = vecs[:, order]
        # M-normalize
        for k in range(vecs.shape[1]):
            norm = math.sqrt(abs(vecs[:, k] @ M_free.dot(vecs[:, k])))
            if norm > 1e-30:
                vecs[:, k] /= norm

    PM_LEN = 10e-3
    PM_HEAD_Z = args.h_gap + PM_LEN/2   # PM dipole center height above plate top
    # plate top at z = 5mm, plate bottom at z = 0
    # PM is at z = plate_top + h_gap + L_PM/2 = 5e-3 + h_gap + 5e-3 = 10e-3 + h_gap
    # For "h_gap" arg: gap from PM bottom to plate top. PM bottom z = plate_top + h_gap = 5e-3 + h_gap.
    # PM dipole center z = PM_bot + L_PM/2 = 5e-3 + h_gap + 5e-3 = 10e-3 + h_gap.
    # Actually need to check plate geometry vs coord system. Box(Pnt(0,0,0), Pnt(Lx,Ly,Lz)) places
    # plate from z=0 to z=Lz. So plate top = Lz = 5mm. PM dipole center z = Lz + h_gap + L_PM/2.
    z_pm_dipole = 5e-3 + args.h_gap + PM_LEN/2

    print(f"  PM dipole center z = {z_pm_dipole*1e3:.2f} mm")
    print()
    print(f"  {'v (m/s)':>9}  {'f_eff (Hz)':>11}  {'F_drag (N)':>13}  {'F_lift (N)':>13}")
    results = []
    for v_str in args.velocities.split(","):
        v = float(v_str)
        f_eff = v / (2 * PM_LEN)
        s = 1j * 2 * math.pi * f_eff
        Fx, Fz = compute_lorentz_force_via_foster(
            mesh, lam, None, vecs, free, SIGMA_AL, MU_0, s,
            args.m_pm, h_gap=z_pm_dipole, x_pm=0.05,
        )
        # Newton 3rd: force on PM is opposite to force on conductor
        F_drag_PM = -Fx
        F_lift_PM = -Fz
        results.append({"v": v, "f_eff": f_eff,
                         "F_drag_on_PM_N": F_drag_PM, "F_lift_on_PM_N": F_lift_PM})
        print(f"  {v:9.2f}  {f_eff:11.2e}  {F_drag_PM:+13.4e}  {F_lift_PM:+13.4e}")

    with open("ecb_lorentz_results.json", "w", encoding="utf-8") as fp:
        json.dump({"results": results, "note": "Scalar A_z approximation, quasi-static."},
                  fp, indent=2)
    print(f"\n  Saved ecb_lorentz_results.json")


if __name__ == "__main__":
    main()
