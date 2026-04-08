"""
Diagnostic: compare RHS assembly methods for scattered-field SIBC.

Sphere + uniform B0 case. Tests three ways to compute the n x H_inc boundary term:
1. Explicit analytical n x H (known geometry-dependent formula) [REFERENCE]
2. specialcf.normal(3) cross product with H_inc CF
3. Volume form: -nu0 * curl(A_inc) * curl(v) * dx (replaces boundary term)

Also tests BiotSavartCF approach for H_inc.
"""

import math
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

MU_0 = 4e-7 * np.pi


def build_sphere_mesh(R=0.01, a_kelvin=0.05, maxh=None, order=1):
    """Build sphere + Kelvin mesh (same as verify_sphere_sibc.py)."""
    from ngsolve import Mesh, TaskManager
    from netgen.occ import Sphere, Pnt, OCCGeometry, Glue, Vertex
    from netgen.meshing import IdentificationType

    if maxh is None:
        maxh = R / 3

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

    return mesh


def test_rhs_methods():
    R = 0.01
    B0 = 0.001
    freq = 1000
    omega = 2 * np.pi * freq
    H0 = B0 / MU_0
    nu0 = 1.0 / MU_0
    a_kelvin = 5 * R

    # Z_s ratio = 10 (same as verify_sphere_sibc.py default)
    jw_mu0_R = omega * MU_0 * R
    Z_s = 10.0 * jw_mu0_R * (1 + 1j) / math.sqrt(2)
    sibc_coeff = 1j * omega / Z_s

    print("=" * 65)
    print("Diagnostic: n x H_inc RHS assembly methods")
    print("=" * 65)
    print(f"R = {R*1e3:.1f} mm, B0 = {B0*1e3:.1f} mT, H0 = {H0:.1f} A/m")
    print(f"|Z_s| = {abs(Z_s):.4e}, |sibc_coeff| = {abs(sibc_coeff):.4e}")
    print()

    from ngsolve import (HCurl, LinearForm, CF, ds, dx, curl, Integrate, BND,
                         specialcf, InnerProduct)
    from ngsolve import x, y, z, sqrt, IfPos

    mesh = build_sphere_mesh(R=R, a_kelvin=a_kelvin, order=1)
    fes = HCurl(mesh, order=1, dirichlet="GND", complex=True)
    v = fes.TestFunction()
    ndof = fes.ndof
    print(f"Mesh: {mesh.ne} tets, DOFs: {ndof}")
    print(f"Materials: {mesh.GetMaterials()}")
    print(f"Boundaries: {mesh.GetBoundaries()}")

    # --- Incident field ---
    A_inc = CF((-(B0 / 2) * y, (B0 / 2) * x, 0))
    # curl(A_inc) = B0 * z_hat (constant)
    B_inc = CF((0, 0, B0))
    H_inc_cf = CF((0, 0, H0))

    # --- Kelvin reluctivity ---
    r_cf = sqrt(x * x + y * y + z * z)
    r_safe = IfPos(r_cf - 1e-12, r_cf, 1e-12)
    nu_kelvin = nu0 * (r_safe / a_kelvin)**4
    nu_cf = mesh.MaterialCF({"air": nu0, "kelvin": nu_kelvin}, default=nu0)

    # ================================================================
    # Method 0: REFERENCE - explicit analytical n x H (geometry-dependent)
    # ================================================================
    nxH_explicit = CF((y * H0 / R, -x * H0 / R, 0))
    f0 = LinearForm(fes)
    f0 += (-sibc_coeff) * A_inc * v.Trace() * ds("sibc")
    f0 += (-1.0) * nxH_explicit * v.Trace() * ds("sibc")
    f0.Assemble()
    norm0 = f0.vec.Norm()

    # ================================================================
    # Method 1: specialcf.normal(3) cross product (WRONG SIGN)
    # specialcf.normal points OUTWARD from mesh domain (= air outward)
    # For a hole (conductor removed), this points INTO the conductor.
    # SIBC convention: n = outward from conductor = -specialcf.normal
    # ================================================================
    n_mesh = specialcf.normal(3)  # outward from air = inward into conductor
    n_cond = -n_mesh  # outward from conductor (SIBC convention)

    # Method 1a: WRONG sign (n_mesh)
    nxH_wrong = CF((n_mesh[1] * H_inc_cf[2] - n_mesh[2] * H_inc_cf[1],
                    n_mesh[2] * H_inc_cf[0] - n_mesh[0] * H_inc_cf[2],
                    n_mesh[0] * H_inc_cf[1] - n_mesh[1] * H_inc_cf[0]))
    f1a = LinearForm(fes)
    f1a += (-sibc_coeff) * A_inc * v.Trace() * ds("sibc")
    f1a += (-1.0) * nxH_wrong * v.Trace() * ds("sibc")
    f1a.Assemble()
    norm1a = f1a.vec.Norm()

    # Method 1b: CORRECT sign (n_cond = -n_mesh)
    nxH_correct = CF((n_cond[1] * H_inc_cf[2] - n_cond[2] * H_inc_cf[1],
                      n_cond[2] * H_inc_cf[0] - n_cond[0] * H_inc_cf[2],
                      n_cond[0] * H_inc_cf[1] - n_cond[1] * H_inc_cf[0]))
    f1 = LinearForm(fes)
    f1 += (-sibc_coeff) * A_inc * v.Trace() * ds("sibc")
    f1 += (-1.0) * nxH_correct * v.Trace() * ds("sibc")
    f1.Assemble()
    norm1 = f1.vec.Norm()

    # ================================================================
    # Method 2: Cross() helper with corrected sign
    # ================================================================
    try:
        from ngsolve import Cross
        f2 = LinearForm(fes)
        f2 += (-sibc_coeff) * A_inc * v.Trace() * ds("sibc")
        # n_cond x H = -(n_mesh x H) = Cross(H, n_mesh) [since AxB = -BxA]
        nxH_cross = Cross(H_inc_cf, n_mesh)  # = H x n_mesh = -(n_mesh x H) = n_cond x H
        f2 += (-1.0) * nxH_cross * v.Trace() * ds("sibc")
        f2.Assemble()
        norm2 = f2.vec.Norm()
    except Exception as e:
        print(f"  Method 2 (Cross): failed - {e}")
        norm2 = 0.0

    # ================================================================
    # Method 3: Volume form (air only)
    # -int_air nu0 * B_inc . curl(v) dx = -int_air H_inc . curl(v) dx
    # Since nu0 * B_inc = nu0 * mu0 * H_inc = H_inc
    # ================================================================
    f3 = LinearForm(fes)
    f3 += (-sibc_coeff) * A_inc * v.Trace() * ds("sibc")
    f3 += (-H_inc_cf) * curl(v) * dx("air")
    f3.Assemble()
    norm3 = f3.vec.Norm()

    # ================================================================
    # Method 4: Volume form (all domains, with nu_cf)
    # -int nu_cf * B_inc . curl(v) dx
    # ================================================================
    f4 = LinearForm(fes)
    f4 += (-sibc_coeff) * A_inc * v.Trace() * ds("sibc")
    f4 += (-nu_cf) * B_inc * curl(v) * dx
    f4.Assemble()
    norm4 = f4.vec.Norm()

    # ================================================================
    # Method 5: HCurl GridFunction interpolation of A_inc, then curl(gf)
    # ================================================================
    from ngsolve import GridFunction as GF_cls
    gf_A_inc = GF_cls(fes)
    gf_A_inc.Set(A_inc)
    f5 = LinearForm(fes)
    f5 += (-nu_cf) * curl(gf_A_inc) * curl(v) * dx
    f5 += (-sibc_coeff) * gf_A_inc * v.Trace() * ds("sibc")
    f5.Assemble()
    norm5 = f5.vec.Norm()

    # ================================================================
    # Method 6: BiotSavartCF approach (simulate coil as distant solenoid)
    # For uniform B0 z-hat, use a large loop at z = +/- large_dist
    # This tests the full pipeline: BiotSavartCF -> n x H on boundary
    # ================================================================
    try:
        from biot_savart import biot_savart_loop
        # Large loop to approximate uniform field B0 z-hat
        # B_center = mu0 * I / (2*R_loop) -> I = B0 * 2 * R_loop / mu0
        R_loop = 1.0  # 1 meter (much larger than sphere)
        I_loop = B0 * 2 * R_loop / MU_0
        H_biot = biot_savart_loop(center=(0, 0, 0), radius=R_loop,
                                   normal=(0, 0, 1), current=I_loop,
                                   n_segments=200, order=20)
        B_biot = MU_0 * H_biot
        nxH_biot = CF((n_cond[1] * H_biot[2] - n_cond[2] * H_biot[1],
                       n_cond[2] * H_biot[0] - n_cond[0] * H_biot[2],
                       n_cond[0] * H_biot[1] - n_cond[1] * H_biot[0]))
        f6 = LinearForm(fes)
        f6 += (-sibc_coeff) * A_inc * v.Trace() * ds("sibc")
        f6 += (-1.0) * nxH_biot * v.Trace() * ds("sibc")
        f6.Assemble()
        norm6 = f6.vec.Norm()
    except Exception as e:
        print(f"  Method 6 (BiotSavartCF): failed - {e}")
        norm6 = 0.0

    # ================================================================
    # Compare
    # ================================================================
    print()
    print("RHS vector norms:")
    print("-" * 65)
    print(f"  Method 0 (explicit n x H, REFERENCE):  {norm0:.6e}")
    print(f"  Method 1a(n_mesh x H, WRONG sign):     {norm1a:.6e}  "
          f"ratio={norm1a/norm0:.4f}")
    print(f"  Method 1b(-n_mesh x H, CORRECT sign):  {norm1:.6e}  "
          f"ratio={norm1/norm0:.4f}")
    print(f"  Method 2 (Cross(H, n_mesh)):            {norm2:.6e}  "
          f"ratio={norm2/norm0:.4f}" if norm2 > 0 else
          f"  Method 2: FAILED")
    print(f"  Method 3 (volume air only):             {norm3:.6e}  "
          f"ratio={norm3/norm0:.4f}")
    print(f"  Method 4 (volume all, nu_cf):           {norm4:.6e}  "
          f"ratio={norm4/norm0:.4f}")
    print(f"  Method 5 (HCurl GF curl(gf)):           {norm5:.6e}  "
          f"ratio={norm5/norm0:.4f}")
    print(f"  Method 6 (BiotSavartCF):                {norm6:.6e}  "
          f"ratio={norm6/norm0:.4f}" if norm6 > 0 else
          f"  Method 6: FAILED")

    # Difference vectors
    print()
    print("Difference from reference (relative):")
    print("-" * 65)
    for name, fv in [("Method 1a (n_mesh, WRONG)", f1a),
                     ("Method 1b (-n_mesh, CORRECT)", f1),
                     ("Method 3 (vol air)", f3),
                     ("Method 4 (vol all)", f4)]:
        diff = fv.vec.CreateVector()
        diff.data = fv.vec - f0.vec
        rel = diff.Norm() / norm0 if norm0 > 0 else float('inf')
        print(f"  {name:40s}: ||f-f0||/||f0|| = {rel:.4e}")

    # ================================================================
    # Now solve and compare H_t_rms (the actual metric)
    # ================================================================
    from ngsolve import BilinearForm, GridFunction, TaskManager

    u = fes.TrialFunction()
    a_bf = BilinearForm(fes)
    a_bf += nu_cf * curl(u) * curl(v) * dx
    a_bf += sibc_coeff * u.Trace() * v.Trace() * ds("sibc")
    a_bf.Assemble()

    # Analytical reference
    beta = 1j * omega * MU_0 * R
    J_s_max = abs((3.0 / 2.0) * 1j * omega * B0 * R / (beta + Z_s))
    H_t_rms_analytical = J_s_max * math.sqrt(2.0 / 3.0)

    print()
    print(f"Analytical H_t_rms = {H_t_rms_analytical:.4f} A/m")
    print()
    print("Solve results:")
    print("-" * 65)

    sibc_region = mesh.Boundaries("sibc")
    A_sphere = Integrate(CF(1), mesh, BND, definedon=sibc_region).real

    methods = [("Method 0 (explicit)", f0),
               ("Method 1b (-n_mesh, CORRECT)", f1)]
    if norm6 > 0:
        methods.append(("Method 6 (BiotSavartCF)", f6))
    for name, fv in methods:
        gfu = GridFunction(fes)
        with TaskManager():
            gfu.vec.data = a_bf.mat.Inverse(
                fes.FreeDofs(), inverse="pardiso") * fv.vec
        A_total = A_inc + gfu
        At_sq = sum(A_total[i].real * A_total[i].real +
                    A_total[i].imag * A_total[i].imag for i in range(3))
        At_int = Integrate(At_sq, mesh, BND, definedon=sibc_region).real
        At_rms = math.sqrt(max(At_int, 0) / max(A_sphere, 1e-30))
        H_t_rms = abs(1j * omega / Z_s) * At_rms
        err = (H_t_rms - H_t_rms_analytical) / H_t_rms_analytical * 100
        print(f"  {name:40s}: H_t_rms = {H_t_rms:.4f}  ({err:+.1f}%)")


if __name__ == "__main__":
    test_rhs_methods()
