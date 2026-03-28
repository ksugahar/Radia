"""
Sphere + SIBC + uniform field: cross-validation of FEM (Kelvin) vs BEM vs analytical.

Analytical solution for a sphere of radius R with SIBC Z_s in uniform B0 z-hat:
  J_s(theta) = (3/2) * jw*B0*R / (jw*mu0*R + Z_s) * sin(theta) * phi-hat
  H_t_rms = sqrt(3/2) * w*B0*R / |jw*mu0*R + Z_s|

FEM: scattered field formulation with 3D Kelvin open boundary.
  curl(nu curl A_scat) = 0 in air
  (jw/Z_s) * A_scat_t = -(jw/Z_s) * A_inc_t  on sphere surface
  A_scat -> 0 at infinity (Kelvin GND)

BEM: EFIE-SIBC on sphere surface (LaplaceSL + Z_s mass + div constraint).

Usage:
    python verify_sphere_sibc.py
    python verify_sphere_sibc.py --Zs-ratio 0.1    # Z_s/(jw*mu0*R) ratio
    python verify_sphere_sibc.py --Zs-ratio 10.0   # PMC-like
"""

import math
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

MU_0 = 4e-7 * np.pi


# ============================================================
# Analytical solution
# ============================================================
def analytical_sphere_sibc(R, Z_s, B0, omega):
    """Analytical H_t for sphere with SIBC in uniform field B0.

    Returns:
        J_s_max: peak surface current at equator [A/m]
        H_t_rms: RMS surface current over sphere [A/m]
        A_t_max: peak tangential A at equator [T*m]
    """
    beta = 1j * omega * MU_0 * R  # jw*mu0*R
    # J_s(theta) = (3/2) * jw*B0*R / (jw*mu0*R + Z_s) * sin(theta)
    J_s_max = abs((3.0 / 2.0) * 1j * omega * B0 * R / (beta + Z_s))
    # <sin^2(theta)> over sphere = 2/3
    H_t_rms = J_s_max * math.sqrt(2.0 / 3.0)
    # A_t(theta) = -(3/2) * B0*R * Z_s / (jw*mu0*R + Z_s) * sin(theta)
    A_t_max = abs((3.0 / 2.0) * B0 * R * Z_s / (beta + Z_s))
    return {
        'J_s_max': J_s_max,
        'H_t_rms': H_t_rms,
        'A_t_max': A_t_max,
        'beta': beta,
    }


# ============================================================
# FEM: scattered field + Kelvin
# ============================================================
def fem_sphere_sibc(R, Z_s, B0, omega, a_kelvin=None, maxh=None, order=1):
    """FEM with 3D Kelvin: scattered field formulation."""
    from ngsolve import (Mesh, HCurl, BilinearForm, LinearForm, GridFunction,
                         Integrate, curl, dx, ds, CF, BND, TaskManager,
                         InnerProduct, Periodic)
    from netgen.occ import (Sphere, Pnt, OCCGeometry, Glue, Vertex)
    from netgen.meshing import IdentificationType

    if a_kelvin is None:
        a_kelvin = 5 * R
    if maxh is None:
        maxh = R / 3

    nu0 = 1.0 / MU_0
    t0 = time.perf_counter()

    # --- Geometry: sphere shell + Kelvin ---
    # Physical domain: R < r < a_kelvin (air shell around conductor sphere)
    inner_sphere = Sphere(Pnt(0, 0, 0), a_kelvin)
    conductor = Sphere(Pnt(0, 0, 0), R)
    for f in conductor.faces:
        f.name = "sibc"
        f.maxh = R / 4

    air = inner_sphere - conductor
    air.name = "air"

    # Kelvin domain: 0 < r < a_kelvin (maps exterior r > a_kelvin)
    kelvin_ball = Sphere(Pnt(0, 0, 0), a_kelvin)
    kelvin_ball.name = "kelvin"

    # Identify inner sphere surface (physical side) with Kelvin surface
    # In NGSolve: the faces at r = a_kelvin
    for f_air in air.faces:
        if f_air.name != "sibc":
            f_air.name = "kelvin_int"
    for f_kel in kelvin_ball.faces:
        f_kel.name = "kelvin_ext"

    # Periodic identification
    for f_int in [f for f in air.faces if f.name == "kelvin_int"]:
        for f_ext in kelvin_ball.faces:
            try:
                f_int.Identify(f_ext, "kelvin",
                               IdentificationType.CLOSESURFACES)
                break
            except Exception:
                pass

    # GND at center of Kelvin domain (A_scat = 0 at infinity)
    gnd = Vertex(Pnt(0, 0, 0))
    gnd.name = "GND"

    shape = Glue([air, kelvin_ball, gnd])
    geo = OCCGeometry(shape)
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=maxh)
    mesh = Mesh(ngmesh)
    mesh.Curve(order + 1)
    t_mesh = time.perf_counter() - t0

    ne = mesh.ne
    print(f"  FEM mesh: {ne} tets, {mesh.nv} verts ({t_mesh:.1f}s)")
    print(f"  Materials: {mesh.GetMaterials()}")
    print(f"  Boundaries: {mesh.GetBoundaries()}")

    # --- Kelvin reluctivity ---
    from ngsolve import x, y, z, sqrt, IfPos
    r_cf = sqrt(x * x + y * y + z * z)
    r_safe = IfPos(r_cf - 1e-12, r_cf, 1e-12)
    # nu_kelvin = nu0 * (r/a)^4  (3D spherical Kelvin)
    nu_kelvin = nu0 * (r_safe / a_kelvin)**4
    nu_cf = mesh.MaterialCF({"air": nu0, "kelvin": nu_kelvin}, default=nu0)

    # --- FE space ---
    fes = HCurl(mesh, order=order, dirichlet="GND", complex=True)
    u, v = fes.TnT()
    print(f"  DOFs: {fes.ndof}")

    # --- Incident field A_inc for uniform B0 z-hat ---
    # A_inc = (B0/2) * (-y, x, 0)  [gives curl A = B0 z-hat]
    A_inc = CF((-(B0 / 2) * y, (B0 / 2) * x, 0))

    # --- SIBC: scattered field formulation ---
    sibc_coeff = 1j * omega / Z_s

    a_bf = BilinearForm(fes)
    a_bf += nu_cf * curl(u) * curl(v) * dx
    a_bf += sibc_coeff * u.Trace() * v.Trace() * ds("sibc")
    a_bf.Assemble()

    # RHS: -a(A_inc, v) from total field decomposition A = A_inc + A_scat
    # f(v) = -(jw/Z_s)*<A_inc, v>_sibc - <n x H_inc, v>_sibc
    # The second term (H_inc boundary) is equivalent to -integral nu*curl(A_inc)*curl(v) dx
    # in a single-domain mesh.  For the Kelvin formulation, use the boundary form
    # directly to avoid incorrect contributions from the Kelvin domain.
    # n x H_inc = r_hat x (H0*z) = H0*(r_hat x z_hat)
    # In Cartesian: cross([x,y,z]/R, [0,0,H0]) = [y*H0/R, -x*H0/R, 0]
    H0 = B0 / MU_0
    nxH_inc = CF((y * H0 / R, -x * H0 / R, 0))
    f_lf = LinearForm(fes)
    f_lf += (-sibc_coeff) * A_inc * v.Trace() * ds("sibc")
    f_lf += (-1.0) * nxH_inc * v.Trace() * ds("sibc")
    f_lf.Assemble()

    t0 = time.perf_counter()
    gfu = GridFunction(fes)
    with TaskManager():
        gfu.vec.data = a_bf.mat.Inverse(
            fes.FreeDofs(), inverse="pardiso") * f_lf.vec
    t_solve = time.perf_counter() - t0
    print(f"  Solve: {t_solve:.1f}s")

    # --- Extract H_t from SIBC ---
    # Total A on sibc surface: A_total = A_inc + A_scat
    A_total = A_inc + gfu
    At_sq = sum(A_total[i].real * A_total[i].real +
                A_total[i].imag * A_total[i].imag for i in range(3))
    sibc_region = mesh.Boundaries("sibc")
    A_sphere = Integrate(CF(1), mesh, BND, definedon=sibc_region).real
    At_int = Integrate(At_sq, mesh, BND, definedon=sibc_region).real
    At_rms = math.sqrt(max(At_int, 0) / max(A_sphere, 1e-30))
    H_t_rms = abs(1j * omega / Z_s) * At_rms

    # Also compute H_t from J_s = (jw/Z_s) * A_total_t
    # and peak A_t at equator (y=0, x=R, z=0)
    try:
        mip = mesh(R * 1.001, 0, 0)
        At_pt = [A_total[i](mip) for i in range(3)]
        At_pt_mag = math.sqrt(sum(abs(a)**2 for a in At_pt))
    except Exception:
        At_pt_mag = float('nan')

    return {
        'H_t_rms': float(H_t_rms),
        'At_rms': float(At_rms),
        'At_pt_equator': float(At_pt_mag),
        'A_sphere': float(A_sphere),
        'ndof': fes.ndof,
        'ne': ne,
        't_mesh': t_mesh,
        't_solve': t_solve,
    }


# ============================================================
# BEM: EFIE-SIBC on sphere surface
# ============================================================
def bem_sphere_sibc(R, Z_s, B0, omega, maxh=None):
    """BEM with LaplaceSL + Z_s mass on sphere surface."""
    from ngsolve import (Mesh, HDivSurface, SurfaceL2, BilinearForm,
                         LinearForm, GridFunction, Integrate, CF, ds, div,
                         BND, TaskManager)
    from ngsolve.bem import LaplaceSL
    from netgen.occ import Sphere, Pnt, OCCGeometry, Glue

    if maxh is None:
        maxh = R / 3

    t0 = time.perf_counter()

    # Sphere surface mesh
    sph = Sphere(Pnt(0, 0, 0), R)
    for f in sph.faces:
        f.name = "sphere"
    geo = OCCGeometry(Glue(sph.faces))
    ngmesh = geo.GenerateMesh(maxh=maxh)
    mesh = Mesh(ngmesh)
    mesh.Curve(2)
    t_mesh = time.perf_counter() - t0

    fes = HDivSurface(mesh, order=0)
    fes_l2 = SurfaceL2(mesh, order=0)
    u, v = fes.TnT()
    ndof = fes.ndof
    n_p = fes_l2.ndof
    ne = mesh.GetNE(BND)
    print(f"  BEM mesh: {ne} surface elements, {ndof} J-DOFs + {n_p} p-DOFs ({t_mesh:.1f}s)")

    # --- Assemble SL ---
    t0 = time.perf_counter()
    with TaskManager():
        SL_bf = LaplaceSL(u.Trace() * ds, use_fmm=False) * v.Trace() * ds
    # Extract dense matrix (column-by-column via mat-vec)
    SL_mat = np.zeros((ndof, ndof), dtype=complex)
    mat = SL_bf if not hasattr(SL_bf, 'mat') else SL_bf.mat
    try:
        rows, cols, vals = mat.COO()
        for r, c, val in zip(rows, cols, vals):
            SL_mat[int(r), int(c)] = val
    except Exception:
        for i in range(ndof):
            ei = GridFunction(fes)
            ei.vec[:] = 0; ei.vec[i] = 1.0
            res = ei.vec.CreateVector()
            res.data = mat * ei.vec
            SL_mat[:, i] = res.FV().NumPy().copy()
    t_sl = time.perf_counter() - t0
    print(f"  SL assembled: {ndof}x{ndof} ({t_sl:.1f}s)")

    # --- Mass matrix ---
    mass_bf = BilinearForm(fes)
    mass_bf += u.Trace() * v.Trace() * ds
    mass_bf.Assemble()
    Mass_mat = np.zeros((ndof, ndof), dtype=complex)
    rows, cols, vals = mass_bf.mat.COO()
    for r, c, val in zip(rows, cols, vals):
        Mass_mat[int(r), int(c)] = val

    # --- Divergence constraint ---
    q = fes_l2.TestFunction()
    u_J = fes.TrialFunction()
    bf_D = BilinearForm(trialspace=fes, testspace=fes_l2)
    bf_D += div(u_J.Trace()) * q * ds
    bf_D.Assemble()
    D_mat = np.zeros((n_p, ndof), dtype=float)
    try:
        rows, cols, vals = bf_D.mat.COO()
        for r, c, val in zip(rows, cols, vals):
            D_mat[int(r), int(c)] = float(val)
    except Exception:
        for i in range(ndof):
            ei = GridFunction(fes)
            ei.vec[:] = 0; ei.vec[i] = 1.0
            res = ei.vec.CreateVector()
            res.data = bf_D.mat * ei.vec
            D_mat[:, i] = res.FV().NumPy()[:n_p].copy()
    rank_D = np.linalg.matrix_rank(D_mat)
    D_red = D_mat[:rank_D, :]

    # --- RHS: -jw * A_inc on sphere ---
    from ngsolve import x, y
    A_inc_cf = CF((-(B0 / 2) * y, (B0 / 2) * x, 0))
    lf_rhs = LinearForm(fes)
    lf_rhs += A_inc_cf * v.Trace() * ds
    lf_rhs.Assemble()
    f_rhs = -1j * omega * lf_rhs.vec.FV().NumPy().copy().astype(complex)

    # --- Saddle point solve ---
    n_total = ndof + rank_D
    A_block = np.zeros((n_total, n_total), dtype=complex)
    A_block[:ndof, :ndof] = Z_s * Mass_mat + 1j * omega * MU_0 * SL_mat
    A_block[:ndof, ndof:] = D_red.T
    A_block[ndof:, :ndof] = D_red

    rhs_block = np.zeros(n_total, dtype=complex)
    rhs_block[:ndof] = f_rhs

    t0 = time.perf_counter()
    J_p = np.linalg.solve(A_block, rhs_block)
    t_solve = time.perf_counter() - t0
    J_vec = J_p[:ndof]

    # --- H_t from |J| ---
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
    H_t_rms = math.sqrt(Jsq_total / max(area_total, 1e-30))

    print(f"  BEM solve: {t_solve:.3f}s")

    return {
        'H_t_rms': float(H_t_rms),
        'J_vec': J_vec,
        'ndof': ndof,
        'ne': ne,
        't_solve': t_solve,
    }


# ============================================================
# Main
# ============================================================
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Sphere SIBC validation")
    parser.add_argument("--R", type=float, default=0.01, help="Sphere radius [m]")
    parser.add_argument("--B0", type=float, default=0.001, help="Incident B0 [T]")
    parser.add_argument("--freq", type=float, default=1000, help="Frequency [Hz]")
    parser.add_argument("--Zs-ratio", type=float, default=None,
                        help="Z_s / (jw*mu0*R) ratio. Overrides --Zs")
    parser.add_argument("--Zs", type=float, default=None, help="Z_s magnitude [Ohm]")
    parser.add_argument("--maxh", type=float, default=None, help="Mesh size")
    args = parser.parse_args()

    R = args.R
    B0 = args.B0
    omega = 2 * np.pi * args.freq
    H0 = B0 / MU_0

    # Determine Z_s
    jw_mu0_R = omega * MU_0 * R
    if args.Zs_ratio is not None:
        Z_s = complex(args.Zs_ratio * jw_mu0_R, args.Zs_ratio * jw_mu0_R)
        Z_s = args.Zs_ratio * jw_mu0_R * (1 + 1j) / math.sqrt(2)
    elif args.Zs is not None:
        Z_s = args.Zs * (1 + 1j) / math.sqrt(2)
    else:
        # Default: steel-like Z_s/(jw*mu0*R) ~ 10
        Z_s = 10.0 * jw_mu0_R * (1 + 1j) / math.sqrt(2)

    print("=" * 65)
    print("Sphere + SIBC + Uniform Field: Cross-Validation")
    print("=" * 65)
    print(f"R = {R*1e3:.1f} mm, B0 = {B0*1e3:.1f} mT, H0 = {H0:.1f} A/m")
    print(f"f = {args.freq} Hz, omega = {omega:.1f}")
    print(f"|Z_s| = {abs(Z_s):.4e}, jw*mu0*R = {jw_mu0_R:.4e}")
    print(f"Z_s/(jw*mu0*R) = {abs(Z_s)/jw_mu0_R:.2f}")
    print()

    # --- Analytical ---
    print("[1/3] Analytical...")
    ana = analytical_sphere_sibc(R, Z_s, B0, omega)
    print(f"  H_t_rms = {ana['H_t_rms']:.4f} A/m")
    print(f"  J_s_max = {ana['J_s_max']:.4f} A/m")
    print(f"  A_t_max = {ana['A_t_max']:.4e} T*m")
    print()

    # --- BEM ---
    print("[2/3] BEM (EFIE-SIBC)...")
    bem = bem_sphere_sibc(R, Z_s, B0, omega, maxh=args.maxh)
    print(f"  H_t_rms = {bem['H_t_rms']:.4f} A/m")
    print()

    # --- FEM ---
    print("[3/3] FEM (Kelvin + scattered field)...")
    try:
        fem = fem_sphere_sibc(R, Z_s, B0, omega, maxh=args.maxh)
        print(f"  H_t_rms = {fem['H_t_rms']:.4f} A/m")
        print(f"  At_rms  = {fem['At_rms']:.4e} T*m")
    except Exception as e:
        print(f"  FEM failed: {e}")
        fem = {'H_t_rms': float('nan')}
    print()

    # --- Comparison ---
    print("=" * 65)
    print("Comparison")
    print("=" * 65)
    print(f"{'':>15s} {'Analytical':>12s} {'BEM':>12s} {'FEM':>12s}")
    ha, hb, hf = ana['H_t_rms'], bem['H_t_rms'], fem['H_t_rms']
    print(f"{'H_t_rms [A/m]':>15s} {ha:12.4f} {hb:12.4f} {hf:12.4f}")
    if ha > 1e-10:
        print(f"{'vs analytical':>15s} {'ref':>12s} "
              f"{(hb-ha)/ha*100:+12.1f}% {(hf-ha)/ha*100:+12.1f}%")
    print()

    # --- Sweep Z_s ratio ---
    if args.Zs_ratio is None and args.Zs is None:
        print("Z_s ratio sweep:")
        print(f"{'ratio':>8s} {'Analytical':>12s} {'BEM':>12s} {'BEM/Ana':>8s}")
        for ratio in [0.01, 0.1, 1.0, 5.0, 10.0, 50.0]:
            Zs_i = ratio * jw_mu0_R * (1 + 1j) / math.sqrt(2)
            ana_i = analytical_sphere_sibc(R, Zs_i, B0, omega)
            bem_i = bem_sphere_sibc(R, Zs_i, B0, omega, maxh=args.maxh)
            r_ba = bem_i['H_t_rms'] / ana_i['H_t_rms'] if ana_i['H_t_rms'] > 0 else 0
            print(f"{ratio:8.2f} {ana_i['H_t_rms']:12.4f} "
                  f"{bem_i['H_t_rms']:12.4f} {r_ba:8.3f}")
