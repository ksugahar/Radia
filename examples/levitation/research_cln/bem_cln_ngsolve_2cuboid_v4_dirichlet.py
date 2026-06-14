#!/usr/bin/env python3
"""
bem_cln_ngsolve_2cuboid_v4_dirichlet.py

Phase 5b: NGSolve cross-check of the 3D 2-cuboid BEM-CLN polarizability
tensor (A1) using truncation Dirichlet BC at a large outer sphere
(instead of Periodic Kelvin which fails for full-sphere auto-mesh).

Why this approach:
  - Periodic Kelvin (v3) requires copy-meshed congruent surface
    tessellations on inner/outer sphere boundaries; netgen.occ does
    not generate these automatically for closed full spheres
    (slaved = 0 even with explicit gp_Trsf translation).
  - Truncation Dirichlet with R_outer = 10 R_cuboid_max gives
    truncation error << 1% for the dipole tail (Bz ~ 1/R^3).

Production-validated pattern (from src/radia/panels/calc_fem_kelvin.py
when has_kelvin_periodic=False):
  - HCurl(..., nograds=True, dirichlet="outer")
  - curl-curl + AC mass + reg * NU_0 * u * v * dx everywhere (no
    Kelvin restriction needed because no Kelvin domain)
  - BDDC preconditioner + CG
  - bonus_intorder=4 (not needed here since nu is constant, but kept
    for consistency with production)

Expected: m_y per cuboid ~ -7e-3 A m^2 at f=1e5 Hz (per A1 with Aharoni
N_x=0.117, N_y=0.302, N_z=0.581 for 5x2x1 mm Cu).
"""

import math
import sys


def build_truncation_2cuboid_geometry(test_freq=1e5, sigma=5.8e7,
                                       mu0=4*math.pi*1e-7, R_outer=0.05):
    """
    Build 2-cuboid + single large outer sphere geometry.

    Layout:
      Cuboids 1, 2 at +/- 7.5 mm along x, each 5x2x1 mm Cu.
      Outer sphere R_outer (~50 mm) -- large enough that 1/R^3 dipole
      tail truncation error is << 1%.

    Returns:
      geo: OCCGeometry
    """
    from netgen.occ import Box, Pnt, Sphere, OCCGeometry, Glue

    a, b, c = 5e-3, 2e-3, 1e-3
    D = 15e-3
    omega = 2 * math.pi * test_freq
    delta = math.sqrt(2 / (omega * mu0 * sigma))
    face_maxh = delta / 3
    print(f"At f = {test_freq:.2e} Hz: skin depth = {delta*1e3:.4f} mm, "
          f"face maxh = {face_maxh*1e3:.4f} mm")
    print(f"Outer truncation sphere radius = {R_outer*1e3:.1f} mm "
          f"(= {R_outer / max(a, b, c, D/2):.1f} x max cuboid extent)")

    cu1 = Box(Pnt(-D/2 - a/2, -b/2, -c/2),
              Pnt(-D/2 + a/2,  b/2,  c/2))
    cu1.mat("cu")
    cu1.maxh = 0.5e-3
    for f in cu1.faces:
        f.maxh = face_maxh

    cu2 = Box(Pnt(D/2 - a/2, -b/2, -c/2),
              Pnt(D/2 + a/2,  b/2,  c/2))
    cu2.mat("cu")
    cu2.maxh = 0.5e-3
    for f in cu2.faces:
        f.maxh = face_maxh

    # Outer truncation sphere -- tag boundary face BEFORE subtraction
    outer_sphere = Sphere(Pnt(0, 0, 0), R_outer)
    outer_sphere.mat("air")
    outer_sphere.maxh = 5e-3
    for face in outer_sphere.faces:
        face.name = "outer"
    air = outer_sphere - cu1 - cu2

    return OCCGeometry(Glue([cu1, cu2, air]))


def solve_eddy_truncation(geo, freq, sigma_cu=5.8e7, mu0=4*math.pi*1e-7,
                            B0_y=1.0, order=2, reg=1e-6):
    """
    Solve eddy current with Dirichlet truncation BC + BDDC iterative.

    Production pattern (no Kelvin path):
      - HCurl(..., nograds=True, dirichlet="outer", complex=True)
      - curl-curl + AC mass + reg gauge regularization (everywhere)
      - BDDC preconditioner
    """
    from ngsolve import (Mesh, HCurl, BilinearForm, LinearForm,
                          GridFunction, CoefficientFunction, Integrate,
                          curl, x, y, z, dx, ds, InnerProduct,
                          Preconditioner, CGSolver, IfPos)
    from ngsolve import TaskManager

    NU_0 = 1.0 / mu0

    mesh = Mesh(geo.GenerateMesh(maxh=10e-3))
    print(f"Mesh: {mesh.ne} elements, {mesh.nv} vertices")

    omega = 2 * math.pi * freq
    print(f"  freq = {freq:.2e} Hz")

    materials = mesh.GetMaterials()
    print(f"  materials = {materials}")
    boundaries = mesh.GetBoundaries()
    has_outer = any("outer" in b for b in boundaries)
    print(f"  boundaries has 'outer' = {has_outer}")
    if not has_outer:
        raise RuntimeError(f"outer boundary missing; have: {boundaries}")

    nu_cf_list = []
    sigma_cf_list = []
    for mat in materials:
        if mat == "cu":
            nu_cf_list.append(NU_0)
            sigma_cf_list.append(sigma_cu)
        else:  # air
            nu_cf_list.append(NU_0)
            sigma_cf_list.append(0)
    nu_cf = CoefficientFunction(nu_cf_list)
    sigma_cf = CoefficientFunction(sigma_cf_list)

    # HCurl with Dirichlet on outer truncation face and nograds=True
    fes = HCurl(mesh, order=order, complex=True,
                 nograds=True, dirichlet="outer")
    print(f"  HCurl order={order} nograds=True dirichlet='outer' DOF: {fes.ndof}")

    # Applied A: A_z = -B0 x for uniform B_y
    A_ext = CoefficientFunction((0, 0, -B0_y * x))

    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    a += nu_cf * InnerProduct(curl(u), curl(v)) * dx(bonus_intorder=4)
    a += 1j * omega * sigma_cf * InnerProduct(u, v) * dx("cu")
    # Gauge regularization: everywhere (no Kelvin restriction needed)
    a += reg * NU_0 * InnerProduct(u, v) * dx

    f = LinearForm(fes)
    f += -1j * omega * sigma_cf * InnerProduct(A_ext, v) * dx("cu")

    prec = Preconditioner(a, "bddc")

    print("  Assembling...")
    with TaskManager():
        a.Assemble()
        f.Assemble()

        gfu = GridFunction(fes)
        print("  Solving (BDDC + CG)...")
        inv = CGSolver(a.mat, prec.mat, complex=True, maxsteps=2000,
                        precision=1e-8, printrates=True)
        gfu.vec.data = inv * f.vec

        A_total = gfu + A_ext
        J = -1j * omega * sigma_cf * A_total

        r1_x = -7.5e-3; r2_x = 7.5e-3
        # m_y = (1/2) integral_V (r x J)_y dV = (1/2) integral (z J_x - x J_z) dV
        # for cuboid 1, shift x by -r1_x = +7.5e-3 (center at -7.5)
        m1_y_cf = 0.5 * (z * J[0] - (x - r1_x) * J[2])
        m2_y_cf = 0.5 * (z * J[0] - (x - r2_x) * J[2])

        in_cu1 = IfPos(2.5e-3 - (x - r1_x),
                       IfPos(2.5e-3 + (x - r1_x), 1, 0), 0)
        in_cu2 = IfPos(2.5e-3 - (x - r2_x),
                       IfPos(2.5e-3 + (x - r2_x), 1, 0), 0)

        m1_y = Integrate(m1_y_cf * in_cu1, mesh, definedon=mesh.Materials("cu"))
        m2_y = Integrate(m2_y_cf * in_cu2, mesh, definedon=mesh.Materials("cu"))

        return m1_y, m2_y


def main():
    print("=== bem_cln_ngsolve_2cuboid_v4_dirichlet.py ===")
    print("Phase 5b: Dirichlet truncation + BDDC + order=2 HCurl")
    print()

    test_freq = 1e5
    geo = build_truncation_2cuboid_geometry(test_freq=test_freq, R_outer=0.05)

    print(f"\n=== Solving at f = {test_freq:.2e} Hz ===")
    try:
        m1, m2 = solve_eddy_truncation(geo, test_freq, order=2)
        print(f"\n  m1_y = {m1:.4e}")
        print(f"  m2_y = {m2:.4e}")
        print()
        # A1 prediction at f=1e5 Hz: m_y ~ -7e-3 A m^2 (per cuboid)
        print(f"  A1 prediction (per cuboid)    ~ -7e-3 A m^2")
        if m1 != 0:
            print(f"  |m1| / |A1|                    = {abs(m1)/7e-3:.3f}")
            print(f"  (target ~ 1.0 with proper gauge fix)")
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()

    print()
    print("=== DONE ===")


if __name__ == "__main__":
    main()
