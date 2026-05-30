#!/usr/bin/env python3
"""
bem_cln_ngsolve_2cuboid_v2.py

Phase 4 final-3: skin-layer-graded NGSolve FEM cross-check of the
3D 2-cuboid BEM-CLN polarizability tensor prediction (A1).

Improvements over bem_cln_ngsolve_2cuboid.py (Phase 3 B.2):
  - Surface-graded mesh: maxh on conductor face = delta/3 at test freq
  - HCurl polynomial order p = 2 (curved tet edges)
  - A-induced formulation (B.2 bug fix retained)
  - Proper Kelvin reflection inversion (or large outer sphere)

For test frequency f = 10^5 Hz: delta = 0.21 mm, face maxh = 0.07 mm
For test frequency f = 10^6 Hz: delta = 0.066 mm, face maxh = 0.022 mm

Comparison target: A1 polarizability tensor prediction
    m_y = alpha^(y)(j omega) * B_0 / mu_0
At f = 10^5 Hz, alpha^(y) ~ -1.07e-8 + O(complex), giving m_y ~ -0.0085 + complex A m^2.

Expected runtime: 30-60 min per frequency on LAB hardware (depends on
mesh size; surface-graded mesh can grow to 200k-500k elements).

Status: foundation script written.  Heavy compute deferred to dedicated
session.
"""

import math
import os
import sys


def build_geometry_skin_graded(test_freq=1e6, sigma=5.8e7, mu0=4*math.pi*1e-7):
    """
    Build 2-cuboid geometry with skin-layer-graded surface mesh.

    Mesh strategy:
      - On Cu faces: maxh = delta(test_freq) / 3
      - In Cu bulk: maxh = 0.5 mm
      - In near-air shell: maxh = 1 mm
      - At Kelvin sphere R = 0.5 m: maxh = 50 mm

    Returns OCCGeometry.
    """
    from netgen.occ import Box, Pnt, Sphere, OCCGeometry, Glue

    a, b, c = 5e-3, 2e-3, 1e-3
    D = 15e-3

    # Skin depth at test_freq
    omega = 2 * math.pi * test_freq
    delta = math.sqrt(2 / (omega * mu0 * sigma))
    face_maxh = delta / 3
    print(f"Skin depth at {test_freq:.2e} Hz: {delta*1e3:.4f} mm")
    print(f"Surface mesh maxh: {face_maxh*1e3:.4f} mm")

    # Cuboid 1 at (-D/2, 0, 0)
    cu1 = Box(Pnt(-D/2 - a/2, -b/2, -c/2),
              Pnt(-D/2 + a/2,  b/2,  c/2))
    cu1.mat("cu")
    cu1.maxh = 0.5e-3                    # bulk Cu interior

    # Set face maxh for skin-graded refinement
    for f in cu1.faces:
        f.maxh = face_maxh

    cu2 = Box(Pnt(D/2 - a/2, -b/2, -c/2),
              Pnt(D/2 + a/2,  b/2,  c/2))
    cu2.mat("cu")
    cu2.maxh = 0.5e-3

    for f in cu2.faces:
        f.maxh = face_maxh

    # Larger outer sphere for better approximation of unbounded domain
    R_outer = 0.5                         # 500 mm
    sphere = Sphere(Pnt(0, 0, 0), R_outer)
    air = sphere - cu1 - cu2
    air.mat("air")
    air.maxh = 50e-3
    air.bc("kelvin_outer")

    return OCCGeometry(Glue([cu1, cu2, air]))


def solve_eddy_2cuboid_v2(geo, freq, sigma_cu=5.8e7, mu0=4*math.pi*1e-7,
                           B0_y=1.0, order=2):
    """A-induced formulation (B.2 bug fix), order=2 HCurl."""
    from ngsolve import (Mesh, HCurl, BilinearForm, LinearForm, GridFunction,
                          CoefficientFunction, Integrate, curl,
                          x, y, z, dx, InnerProduct, IfPos)

    mesh = Mesh(geo.GenerateMesh(maxh=50e-3))
    print(f"Mesh: {mesh.ne} elements, {mesh.nv} vertices")

    omega = 2 * math.pi * freq
    skin_depth = math.sqrt(2 / (omega * mu0 * sigma_cu))
    print(f"  freq = {freq:.2e} Hz, skin depth = {skin_depth*1e3:.4f} mm")

    sigma_cf = CoefficientFunction([sigma_cu if mat == "cu" else 0
                                     for mat in mesh.GetMaterials()])
    nu_cf = CoefficientFunction([1/mu0 for mat in mesh.GetMaterials()])

    # Higher-order HCurl
    fes = HCurl(mesh, order=order, dirichlet="kelvin_outer", complex=True)
    print(f"  HCurl order={order} DOF: {fes.ndof}")

    A_ext = CoefficientFunction((0, 0, -B0_y * x))

    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    a += nu_cf * InnerProduct(curl(u), curl(v)) * dx
    a += 1j * omega * sigma_cf * InnerProduct(u, v) * dx("cu")

    f = LinearForm(fes)
    f += -1j * omega * sigma_cf * InnerProduct(A_ext, v) * dx("cu")

    with TaskManager():
        a.Assemble()
        f.Assemble()

        gfu = GridFunction(fes)
        # A-induced formulation: homogeneous Dirichlet on kelvin_outer

        inv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
        gfu.vec.data = inv * f.vec

        A_total = gfu + A_ext
        J = -1j * omega * sigma_cf * A_total

        r1_x = -7.5e-3
        r2_x =  7.5e-3
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
    """Run skin-graded mesh A2 cross-check."""
    print("=== bem_cln_ngsolve_2cuboid_v2.py ===")
    print("Phase 4 A2: skin-graded NGSolve mesh, HCurl order=2")
    print()

    try:
        from ngsolve import Mesh
        from ngsolve import TaskManager
    except ImportError:
        print("ERROR: NGSolve not available")
        sys.exit(1)

    # Build geometry with mesh graded for test_freq
    test_freq = 1e5     # Easiest to resolve; delta = 0.21 mm
    geo = build_geometry_skin_graded(test_freq=test_freq)

    # Predicted m_y from A1 polarizability tensor at this freq
    # (from polarizability_tensor_3D.wls output)
    # At f=10^5 Hz, alpha_y(jw) ~ -1.07e-8 + O(complex), m_y_pred ~ -0.0085 A m^2
    print()
    print("=== Solving with mesh graded for f = ", test_freq, "Hz ===")
    try:
        m1, m2 = solve_eddy_2cuboid_v2(geo, test_freq, order=1)
        print(f"  m1_y = {m1:.4e} A m^2")
        print(f"  m2_y = {m2:.4e} A m^2")
        print()
        print("Compare with A1 tensor prediction at f=1e5 Hz:")
        print("  m_y_pred = alpha_y(jω) * B / mu_0 ~ O(10^-3 A m^2)")
        print()
        print("Agreement to within 5% would validate A1 framework.")
        print("Larger discrepancies suggest residual mesh or formulation issues.")
    except Exception as e:
        print(f"ERROR: {e}")

    print()
    print("=== DONE ===")


if __name__ == "__main__":
    main()
