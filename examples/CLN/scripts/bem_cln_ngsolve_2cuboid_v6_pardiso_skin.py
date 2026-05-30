#!/usr/bin/env python3
"""
v6: PARDISO + skin-resolved mesh.  Builds on v5 by refining the
conductor surface to delta/2 (= 0.1 mm at f=1e5).

Expected: m1, m2 magnitudes should drop toward A1 PEC bound
(-1.14e-2 A m^2 per cuboid) as the surface mesh resolves skin.
"""
import math
import sys


def build_skin_2cuboid(test_freq=1e5, sigma=5.8e7, mu0=4*math.pi*1e-7,
                       R_outer=0.05):
    from netgen.occ import Box, Pnt, Sphere, OCCGeometry, Glue

    a, b, c = 5e-3, 2e-3, 1e-3
    D = 15e-3
    omega = 2 * math.pi * test_freq
    delta = math.sqrt(2 / (omega * mu0 * sigma))
    face_maxh = delta / 2
    print(f"At f = {test_freq:.2e} Hz: skin depth = {delta*1e3:.4f} mm, "
          f"face maxh = {face_maxh*1e3:.4f} mm")

    cu1 = Box(Pnt(-D/2 - a/2, -b/2, -c/2),
              Pnt(-D/2 + a/2,  b/2,  c/2))
    cu1.mat("cu")
    cu1.maxh = 0.3e-3
    for f in cu1.faces:
        f.maxh = face_maxh

    cu2 = Box(Pnt(D/2 - a/2, -b/2, -c/2),
              Pnt(D/2 + a/2,  b/2,  c/2))
    cu2.mat("cu")
    cu2.maxh = 0.3e-3
    for f in cu2.faces:
        f.maxh = face_maxh

    outer = Sphere(Pnt(0, 0, 0), R_outer)
    outer.mat("air")
    outer.maxh = 10e-3
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
    mesh = Mesh(geo.GenerateMesh(maxh=10e-3))
    print(f"Mesh: {mesh.ne} elements, {mesh.nv} vertices")

    omega = 2 * math.pi * freq
    materials = mesh.GetMaterials()
    print(f"  materials = {materials}")

    nu_cf_list = [NU_0 for _ in materials]
    sigma_cf_list = [sigma_cu if mat == "cu" else 0 for mat in materials]
    nu_cf = CoefficientFunction(nu_cf_list)
    sigma_cf = CoefficientFunction(sigma_cf_list)

    fes = HCurl(mesh, order=order, complex=True,
                 nograds=True, dirichlet="outer")
    print(f"  HCurl order={order} nograds=True DOF: {fes.ndof}")

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
    print("=== bem_cln_ngsolve_2cuboid_v6_pardiso_skin.py ===")
    print("Skin-resolved mesh + PARDISO direct + nograds gauge fix")
    print()

    test_freq = 1e5
    geo = build_skin_2cuboid(test_freq=test_freq, R_outer=0.05)

    print(f"\n=== Solving at f = {test_freq:.2e} Hz ===")
    try:
        m1, m2 = solve_pardiso(geo, test_freq, order=1)
        print(f"\n  m1_y = {m1:.4e}")
        print(f"  m2_y = {m2:.4e}")
        print()
        print(f"  A1 PEC bound per cuboid:    -1.14e-2 A m^2")
        print(f"  |m1|                        = {abs(m1):.4e} A m^2")
        if m1 != 0:
            print(f"  Ratio |m1| / |A1_PEC|       = {abs(m1)/1.14e-2:.3f}")
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback; traceback.print_exc()

    print()
    print("=== DONE ===")


if __name__ == "__main__":
    main()
