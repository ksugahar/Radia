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
                                   mu0=4*math.pi*1e-7, R_K=0.05):
    """
    Build 2-cuboid + Kelvin-sphere-pair geometry using the canonical
    offset-Kelvin convention (cf.
    examples/kelvin_transformation/AdaptiveMesh/A-formulation/
    CircularCoil_A_formulation_with_Kelvin.py).

    Layout:
      Cuboids 1, 2 at +/- 7.5 mm along x in inner sphere.
      Inner sphere R_K=50 mm at origin contains cuboids + air.
      Outer sphere R_K=50 mm at (offset_x, 0, 0) with offset_x = 2*R_K
      is the Kelvin image of physical exterior; its CENTER maps to
      physical infinity (r' = 0 there).
      Periodic identification: inner sphere outer face <-> outer sphere
      outer face.

    Returns:
      geo: OCCGeometry
      kelvin_center: (offset_x, 0, 0) -- center of outer (Kelvin) sphere
    """
    from netgen.occ import (Box, Pnt, Sphere, OCCGeometry, Glue,
                              IdentificationType, gp_Trsf, gp_Vec)

    a, b, c = 5e-3, 2e-3, 1e-3
    D = 15e-3
    offset_x = 2 * R_K
    omega = 2 * math.pi * test_freq
    delta = math.sqrt(2 / (omega * mu0 * sigma))
    face_maxh = delta / 3
    print(f"At f = {test_freq:.2e} Hz: skin depth = {delta*1e3:.4f} mm, "
          f"face maxh = {face_maxh*1e3:.4f} mm")
    print(f"Kelvin pair: inner sphere @ (0,0,0), outer sphere @ ({offset_x*1e3:.1f}, 0, 0) mm, R_K = {R_K*1e3:.1f} mm")

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

    # Inner sphere at origin: tag the sphere face BEFORE the subtraction
    # so the name propagates through cu1/cu2 subtraction to inner_air.
    inner_sphere = Sphere(Pnt(0, 0, 0), R_K)
    inner_sphere.mat("air")
    inner_sphere.maxh = 5e-3
    for face in inner_sphere.faces:
        face.name = "kelvin_int"
    inner_air = inner_sphere - cu1 - cu2

    outer_sphere = Sphere(Pnt(offset_x, 0, 0), R_K)
    outer_sphere.mat("kelvin")
    outer_sphere.maxh = 5e-3
    for face in outer_sphere.faces:
        face.name = "kelvin_ext"

    kelvin_int_face = None
    for face in inner_air.faces:
        if face.name == "kelvin_int":
            kelvin_int_face = face; break
    kelvin_ext_face = None
    for face in outer_sphere.faces:
        if face.name == "kelvin_ext":
            kelvin_ext_face = face; break

    if kelvin_int_face is None or kelvin_ext_face is None:
        raise RuntimeError(f"Failed to identify Kelvin boundary faces: "
                           f"int={kelvin_int_face}, ext={kelvin_ext_face}")

    # Closed sphere has no vertex anchor; auto-detect fails (slaved=0).
    # Pass explicit translation: kelvin_ext (slave) -> kelvin_int (master).
    trsf = gp_Trsf.Translation(gp_Vec(-offset_x, 0, 0))
    kelvin_ext_face.Identify(kelvin_int_face, "kelvin_periodic",
                              IdentificationType.PERIODIC, trsf)

    geo = OCCGeometry(Glue([cu1, cu2, inner_air, outer_sphere]))
    return geo, (offset_x, 0.0, 0.0)


def solve_eddy_kelvin(geo, freq, kelvin_center=(0.1, 0, 0),
                       sigma_cu=5.8e7, mu0=4*math.pi*1e-7,
                       B0_y=1.0, order=2, R_K=0.05, reg=1e-6):
    """
    Solve eddy current with Kelvin-transformed exterior + BDDC iterative.

    Production pattern (from src/radia/panels/calc_fem_kelvin.py):
      - HCurl(..., nograds=True): removes pure-gradient basis subspace
      - Periodic(base_fes): far-field closed by Kelvin identification
      - Gauge regularization reg * NU_0 * u * v * dx(non_kelvin):
        nograds removes H1-gradient subspace but residual harmonic
        cohomology still leaves curl-curl singular.  Restrict mass
        regularization to non-Kelvin materials because nu_cf -> 0 in
        Kelvin region would let constant mass term dominate.
      - bonus_intorder=4 on curl-curl integral for curved boundary
        accuracy at order >= 2.

    Reference: radia-mcp kelvin_transformation(topic="verified_recipe").

    nu in Kelvin region: (R_K / |r - offset|)^2 * (1/mu_0)  (Sugahara 2022)
    A_ext background: uniform B_y in inner region only; zero in Kelvin region.
    """
    from ngsolve import (Mesh, HCurl, Periodic, BilinearForm, LinearForm,
                          GridFunction, CoefficientFunction, Integrate,
                          curl, x, y, z, dx, ds, InnerProduct,
                          Preconditioner, CGSolver, IfPos)

    NU_0 = 1.0 / mu0

    mesh = Mesh(geo.GenerateMesh(maxh=10e-3))
    print(f"Mesh: {mesh.ne} elements, {mesh.nv} vertices")

    omega = 2 * math.pi * freq
    print(f"  freq = {freq:.2e} Hz")

    # Kelvin-modulated nu CoefficientFunction.  Distance r' is measured
    # from the Kelvin sphere CENTER (image of physical infinity).
    kx, ky, kz = kelvin_center
    r_prime_sq = (x - kx)**2 + (y - ky)**2 + (z - kz)**2 + 1e-20
    nu_inner = NU_0
    # Kelvin-domain nu: nu' = (r'/R_K)^2 * nu_0 (Nagamine CEFC 2026 eq. 9)
    nu_kelvin = r_prime_sq / R_K**2 * NU_0

    # Per-material switch
    materials = mesh.GetMaterials()
    print(f"  materials = {materials}")
    nu_cf_list = []
    sigma_cf_list = []
    for mat in materials:
        if mat == "cu":
            nu_cf_list.append(NU_0)
            sigma_cf_list.append(sigma_cu)
        elif mat == "kelvin":
            nu_cf_list.append(nu_kelvin)
            sigma_cf_list.append(0)
        else:  # air or default
            nu_cf_list.append(NU_0)
            sigma_cf_list.append(0)
    nu_cf = CoefficientFunction(nu_cf_list)
    sigma_cf = CoefficientFunction(sigma_cf_list)

    non_kelvin_mats = [m for m in materials if "kelvin" not in m.lower()]
    print(f"  non_kelvin_mats = {non_kelvin_mats}")

    # HCurl with periodic identification on kelvin_int/ext faces.
    # nograds=True is REQUIRED for curl-curl gauge fixing.  Without it
    # BDDC factorization hits ZGETRF::info > 0 and direct solvers
    # return gradient-polluted results.
    base_fes = HCurl(mesh, order=order, complex=True, nograds=True)
    fes = Periodic(base_fes)
    print(f"  Periodic HCurl order={order} nograds=True DOF: {fes.ndof}")

    # Applied A: A_z = -B0 x for uniform B_y (only in physical/Cu region)
    A_ext_inner = CoefficientFunction((0, 0, -B0_y * x))
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
    a += nu_cf * InnerProduct(curl(u), curl(v)) * dx(bonus_intorder=4)
    a += 1j * omega * sigma_cf * InnerProduct(u, v) * dx("cu")
    # Gauge regularization restricted to non-Kelvin materials
    # (nu_kelvin -> 0 at r=0; constant mass term would dominate there).
    if reg > 0 and non_kelvin_mats:
        a += reg * NU_0 * InnerProduct(u, v) * dx("|".join(non_kelvin_mats))

    f = LinearForm(fes)
    f += -1j * omega * sigma_cf * InnerProduct(A_ext, v) * dx("cu")

    prec = Preconditioner(a, "bddc")

    with TaskManager():
        a.Assemble()
        f.Assemble()

        gfu = GridFunction(fes)

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
        from ngsolve import TaskManager
    except ImportError:
        print("NGSolve not available")
        sys.exit(1)

    test_freq = 1e5
    geo, kelvin_center = build_kelvin_2cuboid_geometry(test_freq=test_freq)

    print(f"\n=== Solving at f = {test_freq:.2e} Hz ===")
    try:
        m1, m2 = solve_eddy_kelvin(geo, test_freq,
                                    kelvin_center=kelvin_center, order=2)
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
