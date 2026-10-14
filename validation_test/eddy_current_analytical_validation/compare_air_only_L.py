"""Air-only L comparison: 2D axisym vs 3D.

No workpiece, no eddy currents, no SIBC. Just a circular coil + Kelvin
open boundary. Checks whether 2D axisym H1(A_phi) and 3D HCurl(A) agree
on the pure self-inductance of a coil in free space.

If this test shows < 1% gap, the 2D vs 3D SIBC gap is from SIBC physics
or mesh-at-WP-boundary. If this test shows ~5% gap, the gap is
discretization error in the base formulation.

Usage:
    python compare_air_only_L.py
    python compare_air_only_L.py --order 2     # HCurl order 2 (slow)
"""

import argparse
import math
import sys
import os
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

import numpy as np
from ngsolve import (Mesh, H1, Periodic, HCurl, BilinearForm, LinearForm,
                     GridFunction, Integrate, Conj, grad, curl,
                     dx, ds, CF, CoefficientFunction, BND, VOL,
                     TaskManager, IfPos, specialcf)
from ngsolve import x, y, z, sqrt as ngsqrt
from netgen.occ import (WorkPlane, OCCGeometry, Glue, Axes, Pnt, Dir, Axis,
                        Vertex, MoveTo, X, Y, Z, Sphere)
from netgen.meshing import IdentificationType

MU_0 = 4e-7 * math.pi
NU_0 = 1.0 / MU_0


def solve_2d_axisym_air(R_coil, a_coil, a_kelvin, z_offset,
                         maxh_air, maxh_coil, order, I_total=1.0):
    """2D axisym, no workpiece, circular coil, Kelvin open BC.
    Returns L_air [H] = 2*W/I^2 from volume energy integral.
    """
    # Interior half-disk
    wp_in = WorkPlane()
    inner_full = wp_in.Circle(a_kelvin).Face()

    # Circular coil cross-section (matches 3D torus)
    wp_c = WorkPlane()
    wp_c.MoveTo(R_coil + a_coil, 0).Direction(0, 1)
    wp_c.Arc(a_coil, 180).Arc(a_coil, 180)
    coil = wp_c.Face()
    coil.name = "coil"
    coil.maxh = maxh_coil

    # Cut to half x >= 0
    cut_in = MoveTo(-a_kelvin - 0.1, -a_kelvin - 0.1).Rectangle(
        a_kelvin + 0.1, 2 * a_kelvin + 0.2).Face()
    inner_half = inner_full - cut_in
    coil_half = coil - cut_in
    air_in = inner_half - coil_half
    air_in.name = "air"

    # Exterior half-disk (Kelvin)
    wp_ex = WorkPlane(Axes((0, z_offset, 0), n=Z, h=X))
    outer_full = wp_ex.Circle(a_kelvin).Face()
    cut_ex = MoveTo(-a_kelvin - 0.1, z_offset - a_kelvin - 0.1).Rectangle(
        a_kelvin + 0.1, 2 * a_kelvin + 0.2).Face()
    outer_half = outer_full - cut_ex
    outer_half.name = "air_kelvin"

    # Name axis/kelvin edges + Periodic pairing
    inner_kelvin = []
    outer_kelvin = []
    for e in air_in.edges:
        cx = e.center.x
        try:
            v0, v1 = e.vertices
            d0 = math.sqrt(v0.p.x**2 + v0.p.y**2)
            d1 = math.sqrt(v1.p.x**2 + v1.p.y**2)
            if abs(d0 - a_kelvin) < 0.01 and abs(d1 - a_kelvin) < 0.01 and cx > 0.01:
                e.name = "kelvin_int"
                inner_kelvin.append(e)
            elif cx < 0.01:
                e.name = "axis"
            else:
                e.name = "default"
        except:
            pass
    for e in outer_half.edges:
        cx = e.center.x
        try:
            v0, v1 = e.vertices
            d0 = math.sqrt(v0.p.x**2 + (v0.p.y - z_offset)**2)
            d1 = math.sqrt(v1.p.x**2 + (v1.p.y - z_offset)**2)
            if abs(d0 - a_kelvin) < 0.01 and abs(d1 - a_kelvin) < 0.01 and cx > 0.01:
                e.name = "kelvin_ext"
                outer_kelvin.append(e)
            elif cx < 0.01:
                e.name = "axis_ext"
        except:
            pass
    for e in coil_half.edges:
        if e.center.x < 0.01:
            e.name = "axis"
    for ie in inner_kelvin:
        for oe in outer_kelvin:
            iy = ie.center.y
            oy = oe.center.y - z_offset
            if (iy > 0 and oy > 0) or (iy < 0 and oy < 0):
                ie.Identify(oe, "kelvin", IdentificationType.PERIODIC)
                break

    gnd = Vertex(Pnt(0, z_offset, 0))
    gnd.name = "GND"

    shape = Glue([air_in, coil_half, outer_half, gnd])
    geo = OCCGeometry(shape, dim=2)
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=maxh_air, grading=0.3)
    mesh = Mesh(ngmesh)
    mesh.Curve(order)

    fes_pre = H1(mesh, order=order, complex=False,
                 dirichlet="axis|axis_ext", dirichlet_bbnd="GND")
    fes = Periodic(fes_pre)

    # Kelvin reluctivity factor (rho'/a_kelvin)^2 in exterior
    y_loc = y - z_offset
    rho_p = ngsqrt(x * x + y_loc * y_loc)
    rho_safe = IfPos(rho_p - 1e-10, rho_p, 1e-10)
    kelvin_fac = (rho_safe / a_kelvin) ** 2

    mat_names = list(mesh.GetMaterials())
    nu_vals = []
    for m in mat_names:
        if "kelvin" in m:
            nu_vals.append(NU_0 * kelvin_fac)
        else:
            nu_vals.append(NU_0)
    nu_cf = CoefficientFunction(nu_vals)

    r_safe = IfPos(x - 1e-10, x, 1e-10)
    J0 = I_total / (math.pi * a_coil ** 2)

    uh, vh = fes.TnT()
    a_bf = BilinearForm(fes, symmetric=True)
    a_bf += nu_cf / r_safe * grad(uh) * grad(vh) * dx
    f_lf = LinearForm(fes)
    f_lf += J0 * vh * dx("coil")

    with TaskManager():
        a_bf.Assemble()
        f_lf.Assemble()
        gfu = GridFunction(fes)
        gfu.vec.data = a_bf.mat.Inverse(fes.FreeDofs(), inverse="pardiso") * f_lf.vec

    # Energy: W = 0.5 * int nu/r |grad u|^2 dV_axisym = 2*pi* int (1/2) nu/r |grad u|^2 r dr dz
    # With 2D dx = dr*dz (no 2*pi), 3D energy = 2*pi * 2D integral of (1/2) nu/r |grad u|^2
    W_2D = Integrate(0.5 * nu_cf / r_safe * grad(gfu) * grad(gfu), mesh).real
    W_3D = 2 * math.pi * W_2D
    L_air = 2 * W_3D / I_total ** 2
    return {
        'L_air': L_air,
        'ndof': fes.ndof,
        'ne': mesh.ne,
        'W_3D': W_3D,
    }


def solve_3d_air(R_coil, a_coil, R_air, R_kelvin, gap_deg,
                 maxh_air, maxh_coil, order, I_total=1.0):
    """3D, no workpiece, full (or gapped) torus, spherical Kelvin exterior domain.
    Returns L_air [H] from volume energy integral.
    """
    from netgen.occ import Sphere, Cylinder

    wp = WorkPlane(Axes(p=Pnt(R_coil, 0, 0),
                        n=Dir(0, 1, 0), h=Dir(0, 0, 1)))
    circle = wp.Circle(a_coil).Face()
    torus = circle.Revolve(Axis(Pnt(0, 0, 0), Dir(0, 0, 1)), 360 - gap_deg)
    torus.name = "coil"
    torus.maxh = maxh_coil

    air_sphere = Sphere(Pnt(0, 0, 0), R_air)
    air_sphere.name = "air"

    kelvin_outer = Sphere(Pnt(0, 0, 0), R_kelvin)
    kelvin_ext = kelvin_outer - Sphere(Pnt(0, 0, 0), R_air)
    kelvin_ext.name = "kelvin"
    kelvin_ext.maxh = maxh_air * 2
    for f in kelvin_outer.faces:
        f.name = "outer"

    air = air_sphere - torus
    air.name = "air"
    shape = Glue([air, torus, kelvin_ext])

    geo = OCCGeometry(shape)
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=maxh_air, grading=0.3)
    mesh = Mesh(ngmesh)
    mesh.Curve(max(order + 1, 3))

    fes = HCurl(mesh, order=order, dirichlet="outer", complex=False)

    r_xy = ngsqrt(x * x + y * y)
    r_safe = IfPos(r_xy - 1e-10, r_xy, 1e-10)
    J_source = CF((-y / r_safe, x / r_safe, 0)) * (I_total / (math.pi * a_coil**2))

    # Kelvin reluctivity
    r_sq = x * x + y * y + z * z
    kelvin_fac = R_air ** 2 / (r_sq + 1e-20)
    mats = mesh.GetMaterials()
    nu_vals = []
    for m in mats:
        if m == "kelvin":
            nu_vals.append(NU_0 * kelvin_fac)
        else:
            nu_vals.append(NU_0)
    nu_cf = CoefficientFunction(nu_vals)

    u, v = fes.TnT()
    a_bf = BilinearForm(fes, symmetric=True)
    a_bf += nu_cf * curl(u) * curl(v) * dx
    # Regularization (no mass) — HCurl curl-curl has kernel; tiny shift helps
    # pardiso but doesn't change physics.
    a_bf += 1e-10 * NU_0 * u * v * dx

    f_lf = LinearForm(fes)
    f_lf += J_source * v * dx("coil")

    with TaskManager():
        a_bf.Assemble()
        f_lf.Assemble()
        gfu = GridFunction(fes)
        gfu.vec.data = a_bf.mat.Inverse(fes.FreeDofs(), inverse="pardiso") * f_lf.vec

    # L = 2*W/I^2, W = 0.5 * int nu |curl A|^2 dV = 0.5 * int nu |B|^2 dV / (mu0)
    # Wait: B = curl A, nu = 1/mu. Energy = 0.5 * int (1/mu) |B|^2 dV = 0.5 * int nu |curl A|^2 dV.
    B = curl(gfu)
    W = 0.5 * Integrate(nu_cf * B * B, mesh).real
    L_air = 2 * W / I_total ** 2

    return {
        'L_air': L_air,
        'ndof': fes.ndof,
        'ne': mesh.ne,
        'W': W,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--R-coil", type=float, default=0.030)
    ap.add_argument("--a-coil", type=float, default=0.003)
    ap.add_argument("--R-air", type=float, default=0.120)
    ap.add_argument("--R-kelvin", type=float, default=0.240)
    ap.add_argument("--gap-deg", type=float, default=0.0)
    ap.add_argument("--maxh-air", type=float, default=0.012)
    ap.add_argument("--maxh-coil", type=float, default=0.002)
    ap.add_argument("--order", type=int, default=1,
                    help="3D HCurl order; 2D H1 fixed at 3")
    ap.add_argument("--skip-2d", action="store_true")
    ap.add_argument("--skip-3d", action="store_true")
    args = ap.parse_args()

    # Analytical: filamentary loop self-inductance
    R = args.R_coil
    a = args.a_coil
    L_thin_loop = MU_0 * R * (math.log(8 * R / a) - 2)
    # With finite cross-section, slightly different. Use thin-wire as guide.
    L_wheeler = MU_0 * R * (math.log(8 * R / a) - 1.75)  # with internal self

    print()
    print("=" * 92)
    print(f"  Air-only L:  2D axisym (H1 order=3)  vs  3D (HCurl order={args.order})")
    print("=" * 92)
    print(f"  Coil: R={R*1e3:.1f}mm, a={a*1e3:.1f}mm, gap_deg={args.gap_deg} (3D)")
    print(f"  Kelvin: 2D a_kelvin={args.R_air*1e3:.0f}mm (z-offset)")
    print(f"          3D R_air={args.R_air*1e3:.0f}mm, R_kelvin={args.R_kelvin*1e3:.0f}mm")
    print(f"  Mesh: maxh_air={args.maxh_air*1e3:.1f}mm, "
          f"maxh_coil={args.maxh_coil*1e3:.2f}mm")
    print(f"  Analytical (thin-wire, external): L = {L_thin_loop*1e9:.3f} nH")
    print(f"  Analytical (Wheeler,   ext+int):  L = {L_wheeler*1e9:.3f} nH")
    print("=" * 92)

    a_kelvin_2d = args.R_air
    z_offset_2d = a_kelvin_2d * 2.5

    out = {}
    if not args.skip_2d:
        print("\n[a] 2D axisym (no WP)")
        t0 = time.perf_counter()
        try:
            r2d = solve_2d_axisym_air(
                R_coil=R, a_coil=a,
                a_kelvin=a_kelvin_2d, z_offset=z_offset_2d,
                maxh_air=args.maxh_air, maxh_coil=args.maxh_coil,
                order=3,
            )
            dt = time.perf_counter() - t0
            out['2d'] = r2d
            print(f"  L_air = {r2d['L_air']*1e9:.4f} nH   "
                  f"ndof={r2d['ndof']}  ne={r2d['ne']}  t={dt:.1f}s")
        except Exception as e:
            import traceback; traceback.print_exc()
            out['2d'] = None

    if not args.skip_3d:
        print(f"\n[b] 3D (no WP)")
        t0 = time.perf_counter()
        try:
            r3d = solve_3d_air(
                R_coil=R, a_coil=a,
                R_air=args.R_air, R_kelvin=args.R_kelvin,
                gap_deg=args.gap_deg,
                maxh_air=args.maxh_air, maxh_coil=args.maxh_coil,
                order=args.order,
            )
            dt = time.perf_counter() - t0
            out['3d'] = r3d
            print(f"  L_air = {r3d['L_air']*1e9:.4f} nH   "
                  f"ndof={r3d['ndof']}  ne={r3d['ne']}  t={dt:.1f}s")
        except Exception as e:
            import traceback; traceback.print_exc()
            out['3d'] = None

    print()
    print("-" * 80)
    print("  SUMMARY (air only)")
    print("-" * 80)
    if out.get('2d') and out.get('3d'):
        L2 = out['2d']['L_air'] * 1e9
        L3 = out['3d']['L_air'] * 1e9
        print(f"  2D axisym L_air = {L2:.4f} nH")
        print(f"  3D         L_air = {L3:.4f} nH")
        print(f"  dL (3D - 2D)/2D = {(L3-L2)/L2*100:+.3f}%")
        print(f"  Thin-wire  L    = {L_thin_loop*1e9:.4f} nH")
        print(f"  Wheeler    L    = {L_wheeler*1e9:.4f} nH")
        print(f"  2D - Wheeler    = {(L2 - L_wheeler*1e9)/(L_wheeler*1e9)*100:+.3f}%")
        print(f"  3D - Wheeler    = {(L3 - L_wheeler*1e9)/(L_wheeler*1e9)*100:+.3f}%")


if __name__ == "__main__":
    main()
