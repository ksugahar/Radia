"""Conducting sphere in uniform AC field: analytical vs BEM vs FEM.

Analytical solution (Smythe):
  J_s(theta) = (3/2) * jw * B0 * a / (jw*mu0*a + Z_s) * sin(theta)
  H_t_rms   = |J_s_max| * sqrt(2/3)

Three numerical solvers:
  1. BEM:  EFIE-SIBC on sphere surface (LaplaceSL + Z_s mass)
  2. FEM scattered: HCurl Kelvin, A = A_inc + A_scat
  3. FEM total:     HCurl Kelvin, solve for total A directly

Usage:
    python sphere_uniform_field.py
    python sphere_uniform_field.py --material copper --freq 7000
    python sphere_uniform_field.py --sweep
"""

import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))
from em_material import EMMaterial, MU_0


# ============================================================
# Analytical
# ============================================================

def analytical_H_t(R, Z_s, B0, omega):
    """Smythe analytical H_t_rms for sphere with SIBC in uniform B0."""
    beta = 1j * omega * MU_0 * R
    J_s_max = abs((3.0 / 2.0) * 1j * omega * B0 * R / (beta + Z_s))
    H_t_rms = J_s_max * math.sqrt(2.0 / 3.0)
    return H_t_rms


# ============================================================
# BEM: EFIE-SIBC on sphere surface
# ============================================================

def bem_solve(R, Z_s, B0, omega, maxh=None):
    """BEM EFIE-SIBC on sphere surface. Returns H_t_rms."""
    from ngsolve import (Mesh, HDivSurface, SurfaceL2, BilinearForm,
                         LinearForm, GridFunction, Integrate, CF, ds,
                         div, BND, TaskManager, x, y)
    from ngsolve.bem import LaplaceSL
    from netgen.occ import Sphere, Pnt, OCCGeometry, Glue

    if maxh is None:
        maxh = R / 3

    sph = Sphere(Pnt(0, 0, 0), R)
    for f in sph.faces:
        f.name = "sphere"
    geo = OCCGeometry(Glue(sph.faces))
    ngmesh = geo.GenerateMesh(maxh=maxh)
    mesh = Mesh(ngmesh)
    mesh.Curve(2)

    fes = HDivSurface(mesh, order=0)
    fes_l2 = SurfaceL2(mesh, order=0)
    u, v = fes.TnT()
    ndof = fes.ndof
    n_p = fes_l2.ndof

    # SL matrix
    with TaskManager():
        SL_bf = LaplaceSL(u.Trace() * ds, use_fmm=False) * v.Trace() * ds
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

    # Mass matrix
    mass_bf = BilinearForm(fes)
    mass_bf += u.Trace() * v.Trace() * ds
    mass_bf.Assemble()
    Mass_mat = np.zeros((ndof, ndof), dtype=complex)
    rows, cols, vals = mass_bf.mat.COO()
    for r, c, val in zip(rows, cols, vals):
        Mass_mat[int(r), int(c)] = val

    # Divergence constraint
    q = fes_l2.TestFunction()
    u_J = fes.TrialFunction()
    bf_D = BilinearForm(trialspace=fes, testspace=fes_l2)
    bf_D += div(u_J.Trace()) * q * ds
    bf_D.Assemble()
    D_mat = np.zeros((n_p, ndof), dtype=float)
    rows, cols, vals = bf_D.mat.COO()
    for r, c, val in zip(rows, cols, vals):
        D_mat[int(r), int(c)] = float(val)
    rank_D = np.linalg.matrix_rank(D_mat)
    D_red = D_mat[:rank_D, :]

    # RHS: -jw * A_inc
    A_inc_cf = CF((-(B0 / 2) * y, (B0 / 2) * x, 0))
    lf_rhs = LinearForm(fes)
    lf_rhs += A_inc_cf * v.Trace() * ds
    lf_rhs.Assemble()
    f_rhs = -1j * omega * lf_rhs.vec.FV().NumPy().copy().astype(complex)

    # Saddle point solve
    n_total = ndof + rank_D
    A_block = np.zeros((n_total, n_total), dtype=complex)
    A_block[:ndof, :ndof] = Z_s * Mass_mat + 1j * omega * MU_0 * SL_mat
    A_block[:ndof, ndof:] = D_red.T
    A_block[ndof:, :ndof] = D_red
    rhs_block = np.zeros(n_total, dtype=complex)
    rhs_block[:ndof] = f_rhs

    J_p = np.linalg.solve(A_block, rhs_block)
    J_vec = J_p[:ndof]

    # H_t from |J|
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
    return H_t_rms


# ============================================================
# FEM: scattered-field + Kelvin
# ============================================================

def fem_solve(R, Z_s, B0, omega, formulation="scattered",
              a_kelvin=None, maxh=None, order=1):
    """FEM HCurl + Kelvin. Returns H_t_rms."""
    from ngsolve import (Mesh, HCurl, BilinearForm, LinearForm, GridFunction,
                         Integrate, curl, dx, ds, CF, BND, TaskManager,
                         Periodic)
    from netgen.occ import Sphere, Pnt, OCCGeometry, Glue, Vertex
    from netgen.meshing import IdentificationType
    from ngsolve import x, y, z, sqrt, IfPos

    if a_kelvin is None:
        a_kelvin = 5 * R
    if maxh is None:
        maxh = R / 3

    nu0 = 1.0 / MU_0

    # Geometry: air shell (R..a_kelvin) + Kelvin ball (0..a_kelvin)
    inner_sphere = Sphere(Pnt(0, 0, 0), a_kelvin)
    conductor = Sphere(Pnt(0, 0, 0), R)
    for f in conductor.faces:
        f.name = "sibc"
        f.maxh = R / 4

    air = inner_sphere - conductor
    air.name = "air"

    kelvin_ball = Sphere(Pnt(0, 0, 0), a_kelvin)
    kelvin_ball.name = "kelvin"

    for f_air in air.faces:
        if f_air.name != "sibc":
            f_air.name = "kelvin_int"
    for f_kel in kelvin_ball.faces:
        f_kel.name = "kelvin_ext"

    for f_int in [f for f in air.faces if f.name == "kelvin_int"]:
        for f_ext in kelvin_ball.faces:
            try:
                f_int.Identify(f_ext, "kelvin",
                               IdentificationType.CLOSESURFACES)
                break
            except Exception:
                pass

    gnd = Vertex(Pnt(0, 0, 0))
    gnd.name = "GND"

    shape = Glue([air, kelvin_ball, gnd])
    geo = OCCGeometry(shape)
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=maxh)
    mesh = Mesh(ngmesh)
    mesh.Curve(order + 1)

    # Kelvin reluctivity: nu = nu0 * (r/a)^4
    r_cf = sqrt(x * x + y * y + z * z)
    r_safe = IfPos(r_cf - 1e-12, r_cf, 1e-12)
    nu_kelvin = nu0 * (r_safe / a_kelvin)**4
    nu_cf = mesh.MaterialCF({"air": nu0, "kelvin": nu_kelvin}, default=nu0)

    fes = HCurl(mesh, order=order, dirichlet="GND", complex=True)
    u, v = fes.TnT()

    sibc_coeff = 1j * omega / Z_s

    # A_inc for uniform B0 z-hat: A = (B0/2)(-y, x, 0)
    A_inc = CF((-(B0 / 2) * y, (B0 / 2) * x, 0))
    H0 = B0 / MU_0

    if formulation == "scattered":
        # Bilinear form: curl-curl + Robin
        a_bf = BilinearForm(fes)
        a_bf += nu_cf * curl(u) * curl(v) * dx
        a_bf += sibc_coeff * u.Trace() * v.Trace() * ds("sibc")
        a_bf.Assemble()

        # RHS: surface terms only (volume terms cancel for exact A_inc)
        # f(v) = -(jw/Z_s)*<A_inc, v>_sibc - <n x H_inc, v>_sibc
        nxH_inc = CF((y * H0 / R, -x * H0 / R, 0))
        f_lf = LinearForm(fes)
        f_lf += (-sibc_coeff) * A_inc * v.Trace() * ds("sibc")
        f_lf += (-1.0) * nxH_inc * v.Trace() * ds("sibc")
        f_lf.Assemble()

        gfu_scat = GridFunction(fes)
        with TaskManager():
            gfu_scat.vec.data = a_bf.mat.Inverse(
                fes.FreeDofs(), inverse="pardiso") * f_lf.vec

        A_total = A_inc + gfu_scat

    else:
        # Total-field: solve for A directly
        # The "source" is embedded in the boundary condition:
        #   curl(nu curl A) = 0 in air (no volume source)
        #   (jw/Z_s) * A_t = n x H_inc on SIBC
        # Rearranged: a(A, v) = <n x H_inc, v>_sibc
        a_bf = BilinearForm(fes)
        a_bf += nu_cf * curl(u) * curl(v) * dx
        a_bf += sibc_coeff * u.Trace() * v.Trace() * ds("sibc")
        a_bf.Assemble()

        # For uniform field problem without a coil, the RHS is the
        # incident field boundary condition.
        # Total field A satisfies: (jw/Z_s)*A_t + n x H = 0 on SIBC
        # where H = H_inc + H_scat. This gives:
        #   a(A, v) = integral of source terms
        # For uniform field with no volume source, the total-field
        # formulation needs A_inc as Dirichlet-like driving.
        # Use: f(v) = <n x H_inc, v>_sibc (tangential H_inc drives A)
        nxH_inc = CF((y * H0 / R, -x * H0 / R, 0))
        f_lf = LinearForm(fes)
        f_lf += nxH_inc * v.Trace() * ds("sibc")
        f_lf.Assemble()

        gfu = GridFunction(fes)
        with TaskManager():
            gfu.vec.data = a_bf.mat.Inverse(
                fes.FreeDofs(), inverse="pardiso") * f_lf.vec

        A_total = gfu

    # H_t extraction from A_total on SIBC
    At_sq = sum(A_total[i].real * A_total[i].real +
                A_total[i].imag * A_total[i].imag for i in range(3))
    sibc_region = mesh.Boundaries("sibc")
    A_sphere = Integrate(CF(1), mesh, BND, definedon=sibc_region).real
    At_int = Integrate(At_sq, mesh, BND, definedon=sibc_region).real
    At_rms = math.sqrt(max(At_int, 0) / max(A_sphere, 1e-30))
    H_t_rms = abs(1j * omega / Z_s) * At_rms

    return H_t_rms


# ============================================================
# Main
# ============================================================

def run_single(mat, freq, R=0.01, B0=0.001, maxh=None):
    """Run analytical + BEM + FEM for one (material, frequency) point."""
    omega = 2 * math.pi * freq
    Z_s = mat.dowell_Zs(freq, R)
    delta = mat.skin_depth(freq)
    ratio = R / delta if delta < 1e10 else 0

    H_ana = analytical_H_t(R, Z_s, B0, omega)

    t0 = time.perf_counter()
    H_bem = bem_solve(R, Z_s, B0, omega, maxh=maxh)
    t_bem = time.perf_counter() - t0

    t0 = time.perf_counter()
    H_fem_scat = fem_solve(R, Z_s, B0, omega, formulation="scattered",
                           maxh=maxh)
    t_fem_s = time.perf_counter() - t0

    t0 = time.perf_counter()
    H_fem_total = fem_solve(R, Z_s, B0, omega, formulation="total",
                            maxh=maxh)
    t_fem_t = time.perf_counter() - t0

    return {
        'material': mat.name, 'freq': freq, 'R': R,
        'delta': delta, 'R_over_delta': ratio,
        'Z_s': Z_s, 'B0': B0,
        'H_ana': H_ana,
        'H_bem': H_bem, 'err_bem': (H_bem - H_ana) / H_ana * 100,
        'H_fem_scat': H_fem_scat,
        'err_fem_scat': (H_fem_scat - H_ana) / H_ana * 100,
        'H_fem_total': H_fem_total,
        'err_fem_total': (H_fem_total - H_ana) / H_ana * 100,
        't_bem': t_bem, 't_fem_scat': t_fem_s, 't_fem_total': t_fem_t,
    }


def print_result(r):
    """Print one result row."""
    print(f"  {r['material']:>8s} {r['freq']:>8.0f} Hz  "
          f"R/d={r['R_over_delta']:>6.1f}  "
          f"|Z_s|={abs(r['Z_s']):.2e}  "
          f"H_ana={r['H_ana']:.4f}  "
          f"BEM={r['err_bem']:+.1f}%  "
          f"FEM-scat={r['err_fem_scat']:+.1f}%  "
          f"FEM-total={r['err_fem_total']:+.1f}%")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Sphere in uniform field: eddy current validation")
    parser.add_argument("--material", default="copper",
                        choices=["copper", "steel", "aluminum"])
    parser.add_argument("--freq", type=float, default=7000,
                        help="Frequency [Hz]")
    parser.add_argument("--R", type=float, default=0.01,
                        help="Sphere radius [m]")
    parser.add_argument("--B0", type=float, default=0.001,
                        help="Incident field [T]")
    parser.add_argument("--maxh", type=float, default=None,
                        help="Mesh size [m]")
    parser.add_argument("--sweep", action="store_true",
                        help="Run frequency sweep for both materials")
    args = parser.parse_args()

    print("=" * 78)
    print("Conducting Sphere in Uniform AC Field: Analytical Validation")
    print("=" * 78)

    if args.sweep:
        freqs = [100, 500, 1000, 5000, 7000, 10000, 50000, 100000]
        materials = ["copper", "steel"]
        print(f"\n{'Material':>8s} {'f [Hz]':>10s}  {'R/d':>6s}  "
              f"{'|Z_s|':>10s}  {'H_ana':>8s}  "
              f"{'BEM':>8s}  {'FEM-scat':>10s}  {'FEM-total':>10s}")
        print("-" * 78)
        for mname in materials:
            mat = EMMaterial.from_name(mname)
            for f in freqs:
                try:
                    r = run_single(mat, f, R=args.R, B0=args.B0,
                                   maxh=args.maxh)
                    print_result(r)
                except Exception as e:
                    print(f"  {mname:>8s} {f:>8.0f} Hz  FAILED: {e}")
    else:
        mat = EMMaterial.from_name(args.material)
        print(f"\nMaterial: {mat.name} (sigma={mat.sigma:.2e}, mu_r={mat.mu_r})")
        print(f"R = {args.R*1e3:.1f} mm, B0 = {args.B0*1e3:.1f} mT, "
              f"f = {args.freq:.0f} Hz")
        delta = mat.skin_depth(args.freq)
        print(f"Skin depth = {delta*1e3:.3f} mm, R/delta = {args.R/delta:.1f}")
        print()

        r = run_single(mat, args.freq, R=args.R, B0=args.B0, maxh=args.maxh)
        print(f"Analytical:  H_t_rms = {r['H_ana']:.6f} A/m")
        print(f"BEM:         H_t_rms = {r['H_bem']:.6f} A/m  "
              f"({r['err_bem']:+.2f}%)  [{r['t_bem']:.1f}s]")
        print(f"FEM-scat:    H_t_rms = {r['H_fem_scat']:.6f} A/m  "
              f"({r['err_fem_scat']:+.2f}%)  [{r['t_fem_scat']:.1f}s]")
        print(f"FEM-total:   H_t_rms = {r['H_fem_total']:.6f} A/m  "
              f"({r['err_fem_total']:+.2f}%)  [{r['t_fem_total']:.1f}s]")
