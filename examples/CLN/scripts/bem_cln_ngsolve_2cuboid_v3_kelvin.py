#!/usr/bin/env python3
"""
bem_cln_ngsolve_2cuboid_v3_kelvin.py

Phase 5: PROPER Kelvin transformation + iterative-solver NGSolve
cross-check of the 3D 2-cuboid BEM-CLN polarizability tensor (A1).

Improvements over v2:
  (a) HCurl order=2 with BDDC preconditioned iterative solver
      (avoids sparse-direct memory limit at large DOF count).
  (b) Proper Kelvin inversion transformation via Radia's
      kelvin_geometry + kelvin_material infrastructure
      (Sugahara 2022 IEEE TransMag, Nagamine CEFC 2026).
      Inner sphere = physical domain containing cuboids.
      Outer sphere = Kelvin image of physical infinity.
      Periodic identification between sphere outer faces.

Expected: with proper Kelvin BC and order=2 + BDDC, the FEM solve
should converge to the A1 polarizability tensor prediction within
a few percent.

Status: foundation script; production run requires careful mesh
sizing + iterative-solve tuning.  Compute time estimate: 30-90 min
per frequency on LAB hardware.
"""

import math
import os
import sys


def build_kelvin_2cuboid_geometry(test_freq=1e6, sigma=5.8e7,
                                   mu0=4*math.pi*1e-7, R_K=0.05, offset=(0, 0, 0)):
    """
    Build 2-cuboid + Kelvin-sphere-pair geometry using Radia helpers.

    Layout:
      Cuboids 1, 2 at +/- 7.5 mm along x, each 5x2x1 mm Cu.
      Inner sphere R_K=50 mm contains cuboids + surrounding air.
      Outer sphere R_K=50 mm at same offset = Kelvin image of physical
      exterior (Periodic identification with inner sphere outer face).
    """
    from netgen.occ import Box, Pnt, Sphere, OCCGeometry, Glue
    # Note: radia.kelvin_geometry is in Radia package; for portability
    # we set up the Kelvin pair manually following its convention.

    a, b, c = 5e-3, 2e-3, 1e-3
    D = 15e-3
    omega = 2 * math.pi * test_freq
    delta = math.sqrt(2 / (omega * mu0 * sigma))
    face_maxh = delta / 3
    print(f"At f = {test_freq:.2e} Hz: skin depth = {delta*1e3:.4f} mm, "
          f"face maxh = {face_maxh*1e3:.4f} mm")

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

    # Inner sphere: physical air around cuboids
    inner_sphere = Sphere(Pnt(*offset), R_K)
    inner_sphere.mat("air")
    inner_sphere.maxh = 5e-3
    for f in inner_sphere.faces:
        f.name = "kelvin_int"
    inner_air = inner_sphere - cu1 - cu2

    # Outer sphere: Kelvin image
    outer_sphere = Sphere(Pnt(*offset), R_K)
    outer_sphere.mat("kelvin")
    outer_sphere.maxh = 5e-3
    for f in outer_sphere.faces:
        f.name = "kelvin_ext"

    # Periodic identification between inner/outer spheres
    inner_face = None
    for f in inner_sphere.faces:
        if f.name == "kelvin_int":
            inner_face = f; break
    outer_face = None
    for f in outer_sphere.faces:
        if f.name == "kelvin_ext":
            outer_face = f; break

    from netgen.occ import IdentificationType
    inner_face.Identify(outer_face, "kelvin_periodic",
                        IdentificationType.PERIODIC)

    return OCCGeometry(Glue([cu1, cu2, inner_air, outer_sphere]))


def solve_eddy_kelvin(geo, freq, sigma_cu=5.8e7, mu0=4*math.pi*1e-7,
                       B0_y=1.0, order=2, R_K=0.05):
    """
    Solve eddy current with Kelvin-transformed exterior + BDDC iterative.

    nu in Kelvin region: (R_K / |r - offset|)^2 * (1/mu_0)  (Sugahara 2022)
    A_ext background: uniform B_y in inner region only; zero (or
    pull-back) in Kelvin region.
    """
    from ngsolve import (Mesh, HCurl, Periodic, BilinearForm, LinearForm,
                          GridFunction, CoefficientFunction, Integrate,
                          curl, x, y, z, dx, ds, InnerProduct,
                          Preconditioner, CGSolver, IfPos)

    mesh = Mesh(geo.GenerateMesh(maxh=10e-3))
    print(f"Mesh: {mesh.ne} elements, {mesh.nv} vertices")

    omega = 2 * math.pi * freq
    print(f"  freq = {freq:.2e} Hz")

    # Kelvin-modulated nu CoefficientFunction
    r_from_offset = ((x)**2 + (y)**2 + (z)**2)**0.5
    nu_inner = 1/mu0
    # In Kelvin region, nu' = (R_K / r')^2 * nu_0 (spherical 3D conformal)
    nu_kelvin = (R_K / r_from_offset)**2 * (1/mu0)

    # Per-material switch
    materials = mesh.GetMaterials()
    nu_cf_list = []
    sigma_cf_list = []
    for mat in materials:
        if mat == "cu":
            nu_cf_list.append(1/mu0)
            sigma_cf_list.append(sigma_cu)
        elif mat == "kelvin":
            nu_cf_list.append(nu_kelvin)
            sigma_cf_list.append(0)
        else:  # air or default
            nu_cf_list.append(1/mu0)
            sigma_cf_list.append(0)
    nu_cf = CoefficientFunction(nu_cf_list)
    sigma_cf = CoefficientFunction(sigma_cf_list)

    # HCurl with periodic identification on kelvin_int/ext faces
    fes = Periodic(HCurl(mesh, order=order, complex=True))
    print(f"  Periodic HCurl order={order} DOF: {fes.ndof}")

    # Applied A: A_z = -B0 x for uniform B_y (only in physical/Cu region)
    A_ext_inner = CoefficientFunction((0, 0, -B0_y * x))
    # In Kelvin region, no background field (will be carried by periodic BC)
    # For simplicity, set to zero in Kelvin region:
    A_ext_kelvin = CoefficientFunction((0, 0, 0))

    A_ext_cflist = []
    for mat in materials:
        if mat == "kelvin":
            A_ext_cflist.append(A_ext_kelvin)
        else:
            A_ext_cflist.append(A_ext_inner)
    A_ext = CoefficientFunction(A_ext_cflist)

    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    a += nu_cf * InnerProduct(curl(u), curl(v)) * dx
    a += 1j * omega * sigma_cf * InnerProduct(u, v) * dx("cu")

    f = LinearForm(fes)
    f += -1j * omega * sigma_cf * InnerProduct(A_ext, v) * dx("cu")

    # Use BDDC preconditioner (works at order >= 1)
    prec = Preconditioner(a, "bddc")

    a.Assemble()
    f.Assemble()

    gfu = GridFunction(fes)

    # CG solve with BDDC preconditioner (NGSolve API: complex=True, maxsteps, precision)
    inv = CGSolver(a.mat, prec.mat, complex=True, maxsteps=2000, precision=1e-8,
                   printrates=False)
    gfu.vec.data = inv * f.vec

    A_total = gfu + A_ext
    J = -1j * omega * sigma_cf * A_total

    r1_x = -7.5e-3; r2_x = 7.5e-3
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
    print("=== bem_cln_ngsolve_2cuboid_v3_kelvin.py ===")
    print("Phase 5: order=2 HCurl + BDDC + Kelvin transformation")
    print()
    try:
        from ngsolve import Mesh
    except ImportError:
        print("NGSolve not available")
        sys.exit(1)

    test_freq = 1e5
    geo = build_kelvin_2cuboid_geometry(test_freq=test_freq)

    print(f"\n=== Solving at f = {test_freq:.2e} Hz ===")
    try:
        m1, m2 = solve_eddy_kelvin(geo, test_freq, order=2)
        print(f"  m1_y = {m1:.4e}")
        print(f"  m2_y = {m2:.4e}")
        print()
        # A1 prediction at f=1e5 Hz: m_y ~ -7e-3 A m^2 (per cuboid)
        print(f"  A1 prediction (per cuboid)    ~ -7e-3 A m^2")
        print(f"  Agreement ratio NGSolve/A1     ~ {abs(m1)/7e-3:.3f}")
        print(f"  (target ~ 1.0 with proper Kelvin + order=2)")
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()

    print()
    print("=== DONE ===")


if __name__ == "__main__":
    main()
