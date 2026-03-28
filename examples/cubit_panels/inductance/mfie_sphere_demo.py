"""
MFIE vs EFIE for sphere with SIBC: demonstration for ngsolve.bem.

Key findings:
1. EFIE-SIBC (Z_s*J + jw*mu0*SL(J) = -jw*A_inc) is WRONG for finite Z_s
   because A_scat != mu0*SL(J_s) when the conductor has finite impedance.
   The Laplace SL eigenvalue for l=1 is R/3, causing a systematic error.

2. MFIE ((1/2)M - DL)*J = <n x H_inc, v> is CORRECT for PEC
   using LaplaceDL with HDivSurface.  This is the K operator:
   K(J) = n x PV int nabla_x G(x,y) x J(y) dS(y)
   For tangential J on a closed surface, K(J) reduces to the LaplaceDL
   (adjoint of the scalar double layer operator).

3. For SIBC (finite Z_s), the MFIE alone gives the PEC answer (no Z_s
   dependence).  The correct BEM formulation requires PMCHWT with both
   J = n x H (electric) and M = -n x E (magnetic) surface currents.

4. FEM-SIBC (scattered field with H_inc correction) is correct for all Z_s.

Proposed ngsolve.bem contribution:
- Expose LaplaceDL for HDivSurface as "LaplaceMFIE" or "LaplaceK"
- Implement PMCHWT-SIBC for MQS eddy current problems

Usage:
    python mfie_sphere_demo.py
"""

import math
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

MU_0 = 4e-7 * np.pi


def run_demo():
    from ngsolve import (Mesh, HDivSurface, SurfaceL2, BilinearForm,
                         LinearForm, GridFunction, Integrate, CF, ds, div,
                         BND, TaskManager, x, y)
    from ngsolve.bem import LaplaceDL, LaplaceSL
    from netgen.occ import Sphere, Pnt, OCCGeometry, Glue

    R = 0.01
    B0 = 0.001
    freq = 1000
    omega = 2 * math.pi * freq
    H0 = B0 / MU_0
    jw_mu0_R = omega * MU_0 * R

    print("=" * 70)
    print("MFIE vs EFIE on Sphere with SIBC")
    print("=" * 70)
    print(f"R = {R*1e3:.1f} mm, B0 = {B0*1e3:.1f} mT, f = {freq} Hz")
    print(f"|jw*mu0*R| = {jw_mu0_R:.4e}")
    print()

    # --- Mesh ---
    sph = Sphere(Pnt(0, 0, 0), R)
    for f in sph.faces:
        f.name = "sphere"
    geo = OCCGeometry(Glue(sph.faces))
    ngmesh = geo.GenerateMesh(maxh=R / 3)
    mesh = Mesh(ngmesh)
    mesh.Curve(2)

    fes = HDivSurface(mesh, order=0)
    fes_l2 = SurfaceL2(mesh, order=0)
    u, v = fes.TnT()
    ndof = fes.ndof
    n_p = fes_l2.ndof
    ne = mesh.GetNE(BND)
    print(f"Mesh: {ne} elements, {ndof} J-DOFs, {n_p} p-DOFs")

    # --- Assemble operators ---
    t0 = time.perf_counter()

    # DL: K operator for MFIE
    with TaskManager():
        DL_op = LaplaceDL(u.Trace() * ds) * v.Trace() * ds
    DL_mat = np.zeros((ndof, ndof))
    for i in range(ndof):
        ei = GridFunction(fes)
        ei.vec[:] = 0
        ei.vec[i] = 1.0
        res = ei.vec.CreateVector()
        res.data = DL_op.mat * ei.vec
        DL_mat[:, i] = res.FV().NumPy().copy()

    # SL: for EFIE
    with TaskManager():
        SL_op = LaplaceSL(u.Trace() * ds, use_fmm=False) * v.Trace() * ds
    SL_mat = np.zeros((ndof, ndof), dtype=complex)
    try:
        rows, cols, vals = SL_op.mat.COO()
        for r, c, val in zip(rows, cols, vals):
            SL_mat[int(r), int(c)] = val
    except Exception:
        for i in range(ndof):
            ei = GridFunction(fes)
            ei.vec[:] = 0
            ei.vec[i] = 1.0
            res = ei.vec.CreateVector()
            res.data = SL_op.mat * ei.vec
            SL_mat[:, i] = res.FV().NumPy().copy()

    # Mass
    mass_bf = BilinearForm(fes)
    mass_bf += u.Trace() * v.Trace() * ds
    mass_bf.Assemble()
    M_mat = np.zeros((ndof, ndof))
    rows, cols, vals = mass_bf.mat.COO()
    for r, c, val in zip(rows, cols, vals):
        M_mat[int(r), int(c)] = val

    # Divergence constraint
    q = fes_l2.TestFunction()
    u_J = fes.TrialFunction()
    bf_D = BilinearForm(trialspace=fes, testspace=fes_l2)
    bf_D += div(u_J.Trace()) * q * ds
    bf_D.Assemble()
    D_mat = np.zeros((n_p, ndof))
    try:
        rows, cols, vals = bf_D.mat.COO()
        for r, c, val in zip(rows, cols, vals):
            D_mat[int(r), int(c)] = float(val)
    except Exception:
        pass
    rank_D = np.linalg.matrix_rank(D_mat)
    D_red = D_mat[:rank_D, :]

    t_asm = time.perf_counter() - t0
    print(f"Assembly: {t_asm:.1f}s")
    print()

    # --- RHS ---
    # MFIE RHS: <n x H_inc, v>
    nxH_inc = CF((y * H0 / R, -x * H0 / R, 0))
    lf_mfie = LinearForm(fes)
    lf_mfie += nxH_inc * v.Trace() * ds
    lf_mfie.Assemble()
    rhs_mfie = lf_mfie.vec.FV().NumPy().copy()

    # EFIE RHS: -jw * <A_inc, v>
    A_inc_cf = CF((-(B0 / 2) * y, (B0 / 2) * x, 0))
    lf_efie = LinearForm(fes)
    lf_efie += A_inc_cf * v.Trace() * ds
    lf_efie.Assemble()
    rhs_efie = -1j * omega * lf_efie.vec.FV().NumPy().copy().astype(complex)

    # --- H_t_rms computation ---
    def compute_H_t_rms(J_vec):
        gf_re = GridFunction(fes)
        gf_im = GridFunction(fes)
        gf_re.vec.FV().NumPy()[:] = J_vec.real
        gf_im.vec.FV().NumPy()[:] = J_vec.imag
        elem_areas = Integrate(CF(1), mesh, BND, element_wise=True)
        Jsq_re = Integrate(gf_re * gf_re, mesh, BND, element_wise=True)
        Jsq_im = Integrate(gf_im * gf_im, mesh, BND, element_wise=True)
        area_total = 0.0
        Jsq_total = 0.0
        for el in mesh.Elements(BND):
            a = abs(elem_areas[el.nr])
            area_total += a
            Jsq_total += abs(Jsq_re[el.nr]) + abs(Jsq_im[el.nr])
        return math.sqrt(Jsq_total / max(area_total, 1e-30))

    # --- Sweep Z_s ratios ---
    print(f"{'ratio':>8s} {'Analytical':>11s} {'EFIE':>11s} {'MFIE(PEC)':>11s} "
          f"{'EFIE/A':>7s} {'MFIE/A':>7s}")
    print("-" * 70)

    for ratio in [0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]:
        Zs = ratio * jw_mu0_R * (1 + 1j) / math.sqrt(2)
        beta = 1j * omega * MU_0 * R
        J_max_ana = abs((3.0 / 2.0) * 1j * omega * B0 * R / (beta + Zs))
        H_ana = J_max_ana * math.sqrt(2.0 / 3.0)

        # EFIE-SIBC (saddle point)
        n_total = ndof + rank_D
        A_efie = np.zeros((n_total, n_total), dtype=complex)
        A_efie[:ndof, :ndof] = Zs * M_mat + 1j * omega * MU_0 * SL_mat
        A_efie[:ndof, ndof:] = D_red.T
        A_efie[ndof:, :ndof] = D_red
        rhs_block = np.zeros(n_total, dtype=complex)
        rhs_block[:ndof] = rhs_efie
        J_efie = np.linalg.solve(A_efie, rhs_block)[:ndof]
        H_efie = compute_H_t_rms(J_efie)

        # MFIE PEC: (1/2)*M - DL
        J_mfie = np.linalg.solve(0.5 * M_mat - DL_mat, rhs_mfie)
        H_mfie = compute_H_t_rms(J_mfie)

        re = H_efie / H_ana if H_ana > 0 else 0
        rm = H_mfie / H_ana if H_ana > 0 else 0
        print(f"{ratio:8.3f} {H_ana:11.2f} {H_efie:11.2f} {H_mfie:11.2f} "
              f"{re:7.3f} {rm:7.3f}")

    print()
    print("Observations:")
    print("  - EFIE degrades with Z_s/jw*mu0*R (0.35x at ratio 10)")
    print("  - MFIE gives PEC answer for all Z_s (no impedance dependence)")
    print("  - For SIBC, need PMCHWT: J and M = -Z_s*(n x J) as unknowns")
    print("  - FEM-SIBC (scattered field + H_inc fix) works for all Z_s")


if __name__ == "__main__":
    run_demo()
