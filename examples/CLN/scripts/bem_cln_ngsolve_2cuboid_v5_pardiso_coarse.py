#!/usr/bin/env python3
"""
v5: Fast sanity check at coarse mesh + PARDISO direct solver.

Goal: prove the gauge fix (nograds=True + reg) gives results in the
right ballpark, before investing in fine-mesh order=2 BDDC.

Coarse mesh: maxh=20mm outer, maxh=1mm on Cu surface.  Skin depth
0.2mm at 1e5 Hz, so mesh is under-resolved for skin; expect maybe
30-50% magnitude error vs A1 PEC, but should be O(1) magnitude --
proves the formulation is correct.
"""
import math
import sys


def build_coarse_2cuboid(R_outer=0.05):
    from netgen.occ import Box, Pnt, Sphere, OCCGeometry, Glue

    a, b, c = 5e-3, 2e-3, 1e-3
    D = 15e-3
    print(f"Coarse mesh: maxh=20mm outer, maxh=1mm cuboid")

    cu1 = Box(Pnt(-D/2 - a/2, -b/2, -c/2),
              Pnt(-D/2 + a/2,  b/2,  c/2))
    cu1.mat("cu")
    cu1.maxh = 1e-3

    cu2 = Box(Pnt(D/2 - a/2, -b/2, -c/2),
              Pnt(D/2 + a/2,  b/2,  c/2))
    cu2.mat("cu")
    cu2.maxh = 1e-3

    outer = Sphere(Pnt(0, 0, 0), R_outer)
    outer.mat("air")
    outer.maxh = 20e-3
    for face in outer.faces:
        face.name = "outer"
    air = outer - cu1 - cu2

    return OCCGeometry(Glue([cu1, cu2, air]))


def solve_pardiso(geo, freq, sigma_cu=5.8e7, mu0=4*math.pi*1e-7,
                   B0_y=1.0, order=1, reg=1e-6):
    from ngsolve import (Mesh, HCurl, BilinearForm, LinearForm,
                          GridFunction, CoefficientFunction, Integrate,
                          curl, x, y, z, dx, InnerProduct, IfPos)
    from ngsolve import TaskManager

    NU_0 = 1.0 / mu0

    mesh = Mesh(geo.GenerateMesh(maxh=20e-3))
    print(f"Mesh: {mesh.ne} elements, {mesh.nv} vertices")

    omega = 2 * math.pi * freq
    print(f"  freq = {freq:.2e} Hz")
    materials = mesh.GetMaterials()
    print(f"  materials = {materials}")

    nu_cf_list = []
    sigma_cf_list = []
    for mat in materials:
        if mat == "cu":
            nu_cf_list.append(NU_0)
            sigma_cf_list.append(sigma_cu)
        else:
            nu_cf_list.append(NU_0)
            sigma_cf_list.append(0)
    nu_cf = CoefficientFunction(nu_cf_list)
    sigma_cf = CoefficientFunction(sigma_cf_list)

    fes = HCurl(mesh, order=order, complex=True,
                 nograds=True, dirichlet="outer")
    print(f"  HCurl order={order} nograds=True dirichlet='outer' DOF: {fes.ndof}")

    A_ext = CoefficientFunction((0, 0, -B0_y * x))

    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    a += nu_cf * InnerProduct(curl(u), curl(v)) * dx(bonus_intorder=4)
    a += 1j * omega * sigma_cf * InnerProduct(u, v) * dx("cu")
    a += reg * NU_0 * InnerProduct(u, v) * dx

    f = LinearForm(fes)
    f += -1j * omega * sigma_cf * InnerProduct(A_ext, v) * dx("cu")

    print("  Assembling...")
    with TaskManager():
        a.Assemble()
        f.Assemble()

        print("  Solving (PARDISO direct)...")
        gfu = GridFunction(fes)
        gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="pardiso") * f.vec

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
    print("=== bem_cln_ngsolve_2cuboid_v5_pardiso_coarse.py ===")
    print("Fast gauge-fix sanity: coarse mesh + PARDISO direct")
    print()

    test_freq = 1e5
    geo = build_coarse_2cuboid(R_outer=0.05)

    print(f"\n=== Solving at f = {test_freq:.2e} Hz ===")
    try:
        m1, m2 = solve_pardiso(geo, test_freq, order=1)
        print(f"\n  m1_y = {m1:.4e}")
        print(f"  m2_y = {m2:.4e}")
        print()
        print(f"  A1 PEC bound per cuboid:    -1.14e-2 A m^2")
        print(f"  Expected at f=1e5 (c/delta=5): ~50% of PEC, so ~-5e-3 A m^2")
        if m1 != 0:
            print(f"  Ratio |m1| / |A1_PEC|       = {abs(m1)/1.14e-2:.3f}")
            print(f"  Ratio |m1| / |~5e-3|        = {abs(m1)/5e-3:.3f}")
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback; traceback.print_exc()

    print()
    print("=== DONE ===")


if __name__ == "__main__":
    main()
